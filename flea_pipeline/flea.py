#!/usr/bin/env python
"""
Runs the complete pipeline on a set of fasta files.

Input is a file with the following space-seperated information on each
line for each time point:

<file/id> <visit id> <date>

If an alignment is not provided, the first entry is the filename of
the raw sequences for a timepoint. Otherwise it is the base id for a
timepoint.

Usage:
  flea_pipeline.py [options] <file>
  flea_pipeline.py -h | --help

"""

# FIXME: for now, filenames cannot have spaces. Make this more
# robust. For instance use tab-seperated values.

# FIXME: HYPHY does not like dots in sequence names.

import os
import csv
import json
from collections import namedtuple, defaultdict

from configparser import ConfigParser, ExtendedInterpolation

from ruffus import cmdline

from Bio import SeqIO

from flea_pipeline import __version__
from flea_pipeline.util import name_key_to_label
from flea_pipeline.util import str_to_type
import flea_pipeline.pipeline_globals as globals_


parser = cmdline.get_argparse(description='Run complete env pipeline',
                              ignored_args=['jobs', 'use_threads'])
parser.add_argument('file')
parser.add_argument('--config', type=str,
                    help='Configuration file.')
parser.add_argument('--align', type=str,
                    help='Fasta file.')
parser.add_argument('--analyze', type=str,
                    help='Codon-aligned fasta file.')
parser.add_argument('--copynumbers', type=str,
                    help='Copynumber file.')
parser.add_argument('--local', action='store_true',
                    help='Run locally, with processes instead of threads.')
parser.add_argument('--jobs', type=int, default=-1,
                    help=('Number of threads/processes.'
                          ' If negative, use config value'))
parser.add_argument('--ppn', type=int, default=-1,
                    help=('Processors to use on cluster nodes.'
                          ' If negative, use config value.'))


options = parser.parse_args()

# useful directories
data_dir = os.path.dirname(os.path.abspath(options.file))
script_dir = os.path.abspath(os.path.split(__file__)[0])

# standard python logger which can be synchronised across concurrent
# Ruffus tasks
if options.log_file is None:
    options.log_file = os.path.join(data_dir, 'flea.log')
if not options.verbose:
    options.verbose = 1
logger, logger_mutex = cmdline.setup_logging(__name__,
                                             options.log_file,
                                             options.verbose)

do_full = options.align is None and options.analyze is None

# read configuration
if options.config is None:
    configfile = os.path.join(data_dir, 'flea.config')
else:
    configfile = options.config
if not os.path.exists(configfile):
    raise Exception('pipeline configuration file not found')

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(configfile)

# write a copy of the configuration file
with open(os.path.join(data_dir, 'run_parameters.config'), 'w') as handle:
    config.write(handle)


################################################################################
# TODO: encapsulate this timepoint business
Timepoint = namedtuple('Timepoint', ['key', 'label', 'date'])

with open(options.file, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    keymod = lambda k: k
    if do_full:
        keymod = lambda f: os.path.abspath(os.path.join(data_dir, f))
    timepoints = list(Timepoint(keymod(k), a, d) for k, a, d in reader)

key_to_label = dict((t.key, t.label) for t in timepoints)
label_to_date = dict((t.label, t.date) for t in timepoints)

if len(set(t.label for t in timepoints)) != len(timepoints):
    raise Exception('non-unique timepoint labels')


import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()

# set globals
globals_.config = config
globals_.options = options
globals_.script_dir = script_dir
globals_.data_dir = data_dir
globals_.job_script_dir = os.path.join(data_dir, 'drmaa_job_scripts')
globals_.logger = logger
globals_.logger_mutex = logger_mutex
globals_.drmaa_session = drmaa_session
globals_.run_locally = options.local
if options.jobs <= 0:
    options.jobs = config.getint('Jobs', 'jobs')
if options.ppn > 0:
    globals_.ppn = options.ppn
else:
    globals_.ppn = config.getint('Jobs', 'ppn')

globals_.timepoints = timepoints
globals_.key_to_label = key_to_label
globals_.label_to_date = label_to_date

if not os.path.exists(globals_.job_script_dir):
    os.mkdir(globals_.job_script_dir)

#################################################################################

# check alignment options

if do_full:
    if options.copynumbers is not None:
        raise Exception('--copynumber option is invalid without --alignment')
else:
    # check every timepoint has at least three sequences, and
    # every sequence belongs to a timepoint
    if options.alignment is not None:
        records = list(SeqIO.parse(options.alignment, 'fasta'))
    else:
        records = list(SeqIO.parse(options.analyze, 'fasta'))
    if len(set(r.id for r in records)) != len(records):
        raise Exception('non-unique ids')
    timepoint_counts = defaultdict(int)
    for r in records:
        label = name_key_to_label(r.id)
        timepoint_counts[label] += 1
    for k, v in timepoint_counts.items():
        if v < 3:
            raise Exception('timepoint {} has only {} sequences'.format(k, v))

# check reference database
references = SeqIO.parse(globals_.config.get('Parameters', 'reference_db'), 'fasta')
for r in references:
    if "-" in str(r.seq):
        raise Exception('reference sequence {} contains gaps')
    if "*" in r.seq.translate()[-1]:
        raise Exception('reference sequence {} has stop codon'.format(r.id))

# check reference sequence and coordinates
ref_seq = next(SeqIO.parse(globals_.config.get('Parameters', 'reference_sequence'), 'fasta'))
ref_coords = open(globals_.config.get('Parameters', 'reference_coordinates')).read().strip().split()
if '*' in str(ref_seq.seq):
    raise Exception('translated reference sequence has stop codon')
if '-' in str(ref_seq.seq):
    raise Exception('translated reference sequence has a gap character')
if len(ref_seq) != len(ref_coords):
    raise Exception('reference sequence does not match reference coordinates')

# make pipelines
from flea_pipeline.quality_pipeline import make_quality_pipeline
from flea_pipeline.consensus_pipeline import make_consensus_pipeline
from flea_pipeline.alignment_pipeline import make_alignment_pipeline
from flea_pipeline.analysis_pipeline import make_analysis_pipeline
from flea_pipeline.diagnosis_pipeline import make_diagnosis_pipeline

from flea_pipeline.preanalysis_pipeline import make_preanalysis_pipeline


if do_full:
    inputs = list(t.key for t in timepoints)
    for f in inputs:
        if not os.path.exists(f):
            raise Exception('file does not exist: "{}"'.format(f))
    if globals_.config.getboolean('Tasks', 'quality_pipeline'):
        p_qual = make_quality_pipeline()
        p_qual.set_input(input=inputs)
        high_qual_inputs = p_qual
    else:
        high_qual_inputs = inputs

    p_cons = make_consensus_pipeline()
    p_cons['make_input'].set_input(input=high_qual_inputs)

    p_aln = make_alignment_pipeline()
    p_aln.set_input(input=p_cons['cat_all_hqcs'])

    if globals_.config.getboolean('Tasks', 'analysis'):
        p_anl = make_analysis_pipeline()
        p_anl.set_input(input=[p_aln['backtranslate_alignment'],
                               p_cons['cat_copynumbers']])

    if globals_.config.getboolean('Tasks', 'align_ccs'):
        p_diag = make_diagnosis_pipeline()
        p_diag.set_input(input=[p_aln['copy_protein_alignment'],
                                p_aln['backtranslate_alignment'],
                                p_cons['cat_copynumbers'],
                                high_qual_inputs])
else:
    p_pre = make_preanalysis_pipeline()
    if config.align is not None:
        p_pre.set_input(input=options.align)
        p_aln = make_alignment_pipeline()
        p_aln.set_input(input=p_pre['rename_records'])

        if globals_.config.getboolean('Tasks', 'analysis'):
            p_anl = make_analysis_pipeline()
            p_anl.set_input(input=[p_aln['backtranslate_alignment'],
                                   p_pre['make_copynumbers']])

    else:
        pre.set_input(input=options.analyze)
        p_anl = make_analysis_pipeline()
        p_anl.set_input(input=p_pre)


def config_to_dict(config):
    config_dict = {}
    for section in config:
        if section == "Paths":
            continue
        d = dict(config[section])
        config_dict[section] = dict((k, str_to_type(v))
                                    for k, v in d.items())
    return config_dict


if __name__ == '__main__':
    # write run info file
    run_info = {}
    run_info['version'] = __version__
    run_info['configuration'] = config_to_dict(config)
    with open(os.path.join(data_dir, 'run_info.json'), 'w') as handle:
        json.dump(run_info, handle, separators=(",\n", ":"))

    checksum_level = config.getint('Misc', 'checksum_level')
    kwargs = {}
    if options.jobs > 1 and not options.local:
        kwargs['multithread'] = options.jobs
    try:
        logger.info("pipeline start")
        cmdline.run(options, checksum_level=checksum_level, **kwargs)
        logger.info("pipeline stop")
    except Exception as e:
        with logger_mutex:
            logger.error(str(e))
        raise
