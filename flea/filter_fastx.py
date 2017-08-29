#!/usr/bin/env python

import sys

from Bio import SeqIO

from flea.util import run_regexp, iter_sample, is_in_frame
from flea.util import parse_copynumbers


def main(mode, informat, outformat, params):
    records = SeqIO.parse(sys.stdin, informat)
    result = ()
    if mode == "runs":
        run_length = int(params[0])
        cregexp = run_regexp(run_length)
        result = (r for r in records if cregexp.search(str(r.seq)) is None)
    elif mode == 'sample':
        min_n, max_n = map(int, params)        
        records = list(records)
        n = len(records)
        if min_n <= n <= max_n:
            result = records
        if n > max_n:
            records = iter_sample(records, maxsize)
    elif mode == "length":
        min_length, max_length = map(int, params)
        result = (r for r in records if min_length <= len(r.seq) <= max_length)
    elif mode == "inframe":
        allow_stops = True if params[0] == 'true' else False
        result = (r for r in records if is_in_frame(r.seq, allow_stops))
    elif mode == "copynumber":
        cnfile = params[0]
        cndict = parse_copynumbers(cnfile)
        result = (r for r in records if r.id in cndict and cndict[r.id] > 0)
    SeqIO.write(result, sys.stdout, outformat)
    

if __name__ == '__main__':
    mode = sys.argv[1]
    informat = sys.argv[2]
    outformat = sys.argv[3]
    params = sys.argv[4:]
    main(mode, informat, outformat, params)
