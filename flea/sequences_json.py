import sys
import json

from Bio import SeqIO

from flea.util import id_to_label
from flea.util import read_single_record


def main(aln_file, mrca_file, coords_file, ref_protein, ref_coords, outfile):

    result = {}
    observed = {}

    # add sequences
    # TODO: check if MRCA differs from our inferred MRCA
    for r in SeqIO.parse(aln_file, "fasta"):
        if r.name == "Node0":
            r.name = 'MRCA'
        if 'ancestor' in r.name or r.name == 'MRCA':
            if 'Ancestors' not in result:
                result['Ancestors'] = {}
            result['Ancestors'][r.id] = str(r.seq)
        else:
            label = id_to_label(r.id)
            if label not in observed:
                observed[label] = {}
            observed[label][r.id] = str(r.seq)
    result['Observed'] = observed

    # add MRCA
    mrca = read_single_record(mrca_file, 'fasta', True)
    result['MRCA'] = str(mrca.seq)

    # add reference
    reference = read_single_record(ref_protein, 'fasta', True)
    rstr = str(reference.seq)
    with open(coords_file) as handle:
        alignment_coords = json.load(handle)['coordinates']
    ref_coords = open(ref_coords).read().strip().split()
    cmap = dict((c, i) for i, c in enumerate(ref_coords))
    ref_result = ''.join(list(rstr[cmap[c]] for c in alignment_coords))
    result['Reference'] = ref_result

    with open(outfile, 'w') as handle:
        json.dump(result, handle, separators=(",\n", ":"))


if __name__ == "__main__":
    main(*(sys.argv[1:]))
