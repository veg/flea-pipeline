#!/usr/bin/env python

import sys
import re

import click

from Bio import SeqIO

from flea.util import run_regexp, iter_sample, is_in_frame
from flea.util import parse_copynumbers
from flea.util import handle_codon, replace_stop_codons, replace_gapped_codons
from flea.util import update_record_seq
from flea.util import replace_id, id_with_copynumber


caln_pattern = r'[0-9]*[IDM]'


def n_del(part):
    """

    >>> n_del("D")
    1

    >>> n_del("34D")
    34

    >>> n_del("I")
    0

    >>> n_del("3I")
    0

    """
    if not re.match(caln_pattern, part):
        raise Exception('argument does not match regexp')
    if len(part) == 1:
        return 1 if part == 'D' else 0
    return int(part[:-1]) if part[-1] == 'D' else 0


def n_trim(caln):
    """Get amount to trim from start and end.

    >>> n_trim("3D4M2D")
    (3, 2)

    >>> n_trim("3I4M2D")
    (0, 2)

    >>> n_trim("3I4M")
    (0, 0)

    """
    parts = re.findall(caln_pattern, caln)
    return n_del(parts[0]), n_del(parts[-1])


def handle_record(record, userfields):
    query_rcomp, db_rcomp, caln = userfields
    if query_rcomp != db_rcomp:
        rcomp = record.reverse_complement()
        rcomp.id = record.id
        rcomp.description = record.description
        record = rcomp
    a, b = n_trim(caln)
    start = a
    stop = len(record) - b
    return record[start:stop]


def parse_userfile(userfile):
    with open(userfile) as handle:
        lines = list(line.split('\t') for line in handle.read().strip().split('\n'))
    result = dict((name, (qstrand, tstrand, caln)) for name, qstrand, tstrand, caln in lines)
    return result


@click.command()
@click.argument('mode')
@click.argument('informat')
@click.argument('outformat')
@click.argument('params', nargs=-1)
def main(mode, informat, outformat, params):
    records = SeqIO.parse(sys.stdin, informat)
    result = ()
    if mode == "convert":
        # convert file format
        result = records
    elif mode == "runs":
        # filter seqs with homopolymer runs
        run_length = int(params[0])
        cregexp = run_regexp(run_length)
        result = (r for r in records if cregexp.search(str(r.seq)) is None)
    elif mode == 'sample':
        # random sample of inputs
        min_n, max_n = map(int, params)        
        records = list(records)
        n = len(records)
        if min_n <= n <= max_n:
            result = records
        elif n > max_n:
            records = iter_sample(records, max_n)
        else:
            records = ()
    elif mode == "length":
        # filter out sequences that fall outside length limits
        min_length, max_length = map(int, params)
        result = (r for r in records if min_length <= len(r.seq) <= max_length)
    elif mode == "inframe":
        # filter out-of-frame sequences
        allow_stops = True if params[0] == 'true' else False
        result = (r for r in records if is_in_frame(r.seq, allow_stops))
    elif mode == "copynumber":
        # filter out 0-copynumber sequences
        abfile = params[0]
        abdict = parse_copynumbers(abfile)
        result = (r for r in records if r.id in abdict and abdict[r.id] > 0)
    elif mode == "stop_codons":
        # replace stop codons
        result = (replace_stop_codons(record) for record in records)
    elif mode == "gap_codons":
        # replace codons containing gaps
        result = (replace_gapped_codons(record) for record in records)
    elif mode == "degap":
        # remove gaps
        result = (update_record_seq(r, r.seq.ungap('-')) for r in records)
    elif mode == 'prefix':
        # keep sequences with id that matches prefix
        prefix = params[0]
        result = (r for r in records if r.id.startswith(prefix))
    elif mode == 'userout':
        # reverse complement and trim terminal gaps
        userfile = params[0]
        ud = parse_userfile(userfile)
        result = (handle_record(r, ud[r.id]) for r in records if r.id in ud)
    elif mode == "add_copynumber":
        abfile = params[0]
        abdict = parse_copynumbers(abfile)
        result = (replace_id(r, id_with_copynumber(r.id, abdict[r.id])) for r in records)
    else:
        raise Exception('unknown mode: {}'.format(mode))
    SeqIO.write(result, sys.stdout, outformat)
    

if __name__ == '__main__':
    main()
