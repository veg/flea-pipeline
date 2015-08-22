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
  env_pipeline.py [options] <file>
  env_pipeline.py -h | --help

Options:
  --alignment [STR]  Fasta file containing codon-aligned sequences.
  --config [STR]   Configuration file.

"""

# TODO:
# - run as many tasks as possible on the cluster
# - order of input files not guaranteed; make more robust
# - check everywhere usearch is used; summarize # seqs lost
# - are identities re-used correctly?
# - rename files; flea and test subdirectories
# - logging
# - update to latest versions of dependencies
# - replace unnecessary dependencies
# - update README

# FIXME: for now, filenames cannot have spaces. Make this more
# robust. For instance use tab-seperated values.

# FIXME: HYPHY does not like dots in sequence names.

import os
import csv
from collections import namedtuple, defaultdict

from configparser import ConfigParser, ExtendedInterpolation

from ruffus import cmdline, Pipeline

from Bio import SeqIO

from util import n_jobs
import pipeline_globals as globals_


parser = cmdline.get_argparse(description='Run complete env pipeline',
                              ignored_args=['jobs'])
parser.add_argument('file')
parser.add_argument('--config', type=str,
                    help='Configuration file.')
parser.add_argument('--alignment', type=str,
                    help='Codon-aligned fasta file.')


options = parser.parse_args()

if not options.verbose:
    options.verbose = 1

# standard python logger which can be synchronised across concurrent
# Ruffus tasks
logger, logger_mutex = cmdline.setup_logging(__name__,
                                             options.log_file,
                                             options.verbose)


################################################################################
# TODO: encapsulate this timepoint business
Timepoint = namedtuple('Timepoint', ['file', 'id', 'date'])

with open(options.file, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    timepoints = list(Timepoint(f, i, d) for f, i, d in reader)

timepoint_ids = {t.file: t.id for t in timepoints}

do_alignment = options.alignment is None

if do_alignment:
    start_files = list(t.file for t in timepoints)
    for f in start_files:
        if not os.path.exists(f):
            raise Exception('file does not exist: "{}"'.format(f))
else:
    # check every timepoint has at least three sequences, and
    # every sequence belongs to a timepoint
    records = list(SeqIO.parse(options.alignment, 'fasta'))
    timepoint_counts = defaultdict(int)

    if len(set(r.id for r in records)) != len(records):
        raise Exception('non-unique ids')

    # TODO: code duplication
    for r in records:
        found = False
        # check every possible key, in case they are different lengths.
        for k, v in timepoint_ids.items():
            if r.id.startswith(k):
                rest = r.id[len(k):]
                rest = rest.strip('_')
                if not rest:
                    raise Exception('sequence id has no unique part')
                found = True
                timepoint_counts[k] += 1
                break
        if not found:
            raise Exception('record id not found in input file')
    for k, v in timepoint_counts.items():
        if v < 3:
            raise Exception('timepoint {} has only {} sequences'.format(k, v))


if len(set(t.id for t in timepoints)) != len(timepoints):
    raise Exception('non-unique sequence ids')

################################################################################

# useful directories
data_dir = os.path.dirname(os.path.abspath(options.file))
script_dir = os.path.abspath(os.path.split(__file__)[0])

# read configuration
if options.config is None:
    configfile = os.path.join(data_dir, 'env_pipeline.config')
    if not os.path.exists(configfile):
        configfile = os.path.join(script_dir, 'env_pipeline.config')
else:
    configfile = options.config
if not os.path.exists(configfile):
    raise Exception('pipeline configuration file not found')

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(configfile)

# write a copy of the configuration file
with open(os.path.join(data_dir, 'run_parameters.config'), 'w') as handle:
    config.write(handle)

# set globals
globals_.config = config
globals_.options = options
globals_.script_dir = script_dir
globals_.data_dir = data_dir
globals_.qsub_dir = os.path.join(data_dir, 'qsub_files')
globals_.timepoints = timepoints
globals_.timepoint_ids = timepoint_ids
globals_.logger = logger
globals_.logger_mutex = logger_mutex


if not os.path.exists(globals_.qsub_dir):
    os.mkdir(globals_.qsub_dir)

# job handling
n_local_jobs, n_remote_jobs = n_jobs()
options.jobs = max(n_local_jobs, n_remote_jobs)

from alignment_pipeline import make_alignment_pipeline
from hyphy_pipeline import make_hyphy_pipeline
from fasttree_pipeline import make_fasttree_pipeline

if do_alignment:
    inputs = list(os.path.join(data_dir, t.file) for t in timepoints)
    p1 = make_alignment_pipeline()
    p1.set_input(input=inputs)
    if globals_.config.getboolean('Tasks', 'hyphy'):
        p2 = make_hyphy_pipeline(standalone=False)
        # p2.set_input(input=p1)
        for task in p2.head_tasks:
            task.set_input(input=p1.tail_tasks)
    if globals_.config.getboolean('Tasks', 'fasttree'):
        p3 = make_fasttree_pipeline()
        p3.set_input(input=p1)
        for task in p3.head_tasks:
            task.set_input(input=p1.tail_tasks)
else:
    if globals_.config.getboolean('Tasks', 'hyphy'):
        p1 = make_hyphy_pipeline(standalone=True)
        p1.set_input(input=options.alignment)
    if globals_.config.getboolean('Tasks', 'fasttree'):
        p2 = make_fasttree_pipeline()
        p2.set_input(input=options.alignment)


if __name__ == '__main__':
    checksum_level = config.getint('Misc', 'checksum_level')
    cmdline.run(options, checksum_level=checksum_level)
