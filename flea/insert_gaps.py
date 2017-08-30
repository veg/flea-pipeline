from Bio import SeqIO

from flea.util import new_record_seq_str, insert_gaps


def main(infile, msafile, outfile):
    ref, *seqs = list(SeqIO.parse(infile, 'fasta'))
    ref_gapped = next(r for r in SeqIO.parse(msafile, 'fasta')
                      if r.id == ref.id)
    seqs_gapped = (new_record_seq_str(r, insert_gaps(str(ref_gapped.seq),
                                                     str(r.seq),
                                                     '-', '-'))
                   for r in seqs)
    SeqIO.write(seqs_gapped, outfile, "fasta")


if __name__ == "__main__":
    infile, msafile, outfile = sys.argv[1:]
    main(infile, msafile, outfile)
