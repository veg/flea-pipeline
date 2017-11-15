import sys
import json

import click

from Bio import SeqIO

from flea.util import extend_coordinates
from flea.util import read_single_record


@click.command()
@click.argument('mrca')
@click.argument('aligned')
@click.argument('ref_coords')
@click.argument('outfile')
def main(mrca, aligned, ref_coords, outfile):
    pairwise_mrca, pairwise_ref = list(SeqIO.parse(aligned, 'fasta'))
    ref_coords = open(ref_coords).read().strip().split()

    # transfer coordinates to MRCA
    pairwise_coords = extend_coordinates(ref_coords, str(pairwise_ref.seq))
    mrca_coords = list(coord for char, coord in zip(str(pairwise_mrca.seq),
                                                    pairwise_coords)
                       if char != "-")

    # extend MRCA coords to full MSA
    mrca = read_single_record(mrca, 'fasta', True)
    result = extend_coordinates(mrca_coords, str(mrca.seq))

    rdict = {'coordinates': result}
    with open(outfile, 'w') as handle:
        json.dump(rdict, handle, separators=(",\n", ":"))


if __name__ == "__main__":
    main()
