import sys
from Bio import SeqIO
from flea.util import parse_copynumbers, replace_id


def id_with_cn(id_, cn):
    return "{}_cn_{}".format(id_, cn)


def main(msafile, cnfile, outfile):
    id_to_cn = parse_copynumbers(cnfile)
    records = SeqIO.parse(msafile, 'fasta')
    result = (replace_id(r, id_with_cn(r.id, id_to_cn[r.id]))
              for r in records)
    SeqIO.write(result, outfile, "fasta")


if __name__ == "__main__":
    msafile, cnfile, outfile = sys.argv[1:]
    main(msafile, cnfile, outfile)
