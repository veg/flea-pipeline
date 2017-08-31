import sys

from Bio import SeqIO

from flea.util import new_record_seq_str, insert_gaps


def main(bealign_file, msafile, outfile):
    ref, *seqs = list(SeqIO.parse(bealign_file, 'fasta'))
    ref_gapped = next(r for r in SeqIO.parse(msafile, 'fasta')
                      if r.id == ref.id)
    seqs_gapped = (new_record_seq_str(s, insert_gaps(str(ref_gapped.seq),
                                                     str(s.seq),
                                                     '-', '-'))
                   for s in seqs)
    SeqIO.write(seqs_gapped, outfile, "fasta")


if __name__ == "__main__":
    main(*(sys.argv[1:]))
