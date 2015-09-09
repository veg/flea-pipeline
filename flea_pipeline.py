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

Options:
  --alignment [STR]  Fasta file containing codon-aligned sequences.
  --copynumbers [STR]  tsv file containing "<id>\t<number>" lines.
  --config [STR]   Configuration file.

"""

__version__ = "0.2.0"

# FIXME: for now, filenames cannot have spaces. Make this more
# robust. For instance use tab-seperated values.

# FIXME: HYPHY does not like dots in sequence names.

import os
import csv
import json
from collections import namedtuple, defaultdict

from configparser import ConfigParser, ExtendedInterpolation

from ruffus import cmdline, Pipeline

from Bio import SeqIO

from util import n_jobs, name_key_to_label
import pipeline_globals as globals_


parser = cmdline.get_argparse(description='Run complete env pipeline',
                              ignored_args=['jobs'])
parser.add_argument('file')
parser.add_argument('--config', type=str,
                    help='Configuration file.')
parser.add_argument('--alignment', type=str,
                    help='Codon-aligned fasta file.')
parser.add_argument('--copynumbers', type=str,
                    help='Copynumber file.')


options = parser.parse_args()

# useful directories
data_dir = os.path.dirname(os.path.abspath(options.file))
script_dir = os.path.abspath(os.path.split(__file__)[0])

# standard python logger which can be synchronised across concurrent
# Ruffus tasks
if options.log_file is None:
    options.log_file = os.path.join(data_dir, 'flea.log')
logger, logger_mutex = cmdline.setup_logging(__name__,
                                             options.log_file,
                                             options.verbose)

do_alignment = options.alignment is None

# read configuration
if options.config is None:
    configfile = os.path.join(data_dir, 'flea_pipeline.config')
    if not os.path.exists(configfile):
        configfile = os.path.join(script_dir, 'flea_pipeline.config')
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
    if do_alignment:
        keymod = lambda f: os.path.abspath(os.path.join(data_dir, f))
    timepoints = list(Timepoint(keymod(k), a, d) for k, a, d in reader)

key_to_label = dict((t.key, t.label) for t in timepoints)
label_to_date = dict((t.label, t.date) for t in timepoints)

if len(set(t.label for t in timepoints)) != len(timepoints):
    raise Exception('non-unique timepoint labels')


# set globals
globals_.config = config
globals_.options = options
globals_.script_dir = script_dir
globals_.data_dir = data_dir
globals_.qsub_dir = os.path.join(data_dir, 'qsub_files')
globals_.logger = logger
globals_.logger_mutex = logger_mutex

globals_.timepoints = timepoints
globals_.key_to_label = key_to_label
globals_.label_to_date = label_to_date

if not os.path.exists(globals_.qsub_dir):
    os.mkdir(globals_.qsub_dir)

# job handling
n_local_jobs, n_remote_jobs = n_jobs()
options.jobs = max(n_local_jobs, n_remote_jobs)


#################################################################################

# check alignment options

if do_alignment:
    if options.copynumbers is not None:
        raise Exception('--copynumber option is invalid without --alignment')
else:
    # check every timepoint has at least three sequences, and
    # every sequence belongs to a timepoint
    records = list(SeqIO.parse(options.alignment, 'fasta'))
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
from alignment_pipeline import make_alignment_pipeline
from analysis_pipeline import make_analysis_pipeline
from preanalysis_pipeline import make_preanalysis_pipeline


if do_alignment:
    inputs = list(t.key for t in timepoints)
    for f in inputs:
        if not os.path.exists(f):
            raise Exception('file does not exist: "{}"'.format(f))
    p1 = make_alignment_pipeline()
    p1.set_input(input=inputs)
    if globals_.config.getboolean('Tasks', 'analysis'):
        p2 = make_analysis_pipeline(do_hyphy=globals_.config.getboolean('Tasks', 'hyphy_analysis'))
        p2.set_input(input=p1)
        for task in p2.head_tasks:
            task.set_input(input=p1.tail_tasks)
else:
    if globals_.config.getboolean('Tasks', 'analysis'):
        p1 = make_preanalysis_pipeline(options.copynumbers)
        p1.set_input(input=options.alignment)
        p2 = make_analysis_pipeline(do_hyphy=globals_.config.getboolean('Tasks', 'hyphy_analysis'))
        for task in p2.head_tasks:
            task.set_input(input=p1.tail_tasks)


if __name__ == '__main__':
    # write run info file
    run_info = {}
    run_info['version'] = __version__
    config_dict = {}
    for section in config:
        config_dict[section] = dict(config[section])
    run_info['configuration'] = config_dict
    with open(os.path.join(data_dir, 'run_info.json'), 'w') as handle:
        json.dump(run_info, handle, separators=(",\n", ":"))

    checksum_level = config.getint('Misc', 'checksum_level')
    try:
        logger.info("pipeline start")
        cmdline.run(options, checksum_level=checksum_level)
        logger.info("pipeline stop")
    except Exception as e:
        with logger_mutex:
            logger.error(str(e))
        raise
