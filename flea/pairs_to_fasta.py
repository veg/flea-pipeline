from Bio import SeqIO


def main(qcsfile, hqcsfile, pairfile, outdir):
    with open(pairfile) as handle:
        pairs = list(line.strip().split("\t") for line in handle.readlines())
    match_dict = defaultdict(list)
    for pair in pairs:
        if len(pair) != 2:
            warnings.warn('QCS {} did not match any HQCS'.format(pair))
            continue
        qcs, hqcs = pair
        match_dict[hqcs].append(qcs)
    hqcs_records = list(SeqIO.parse(hqcsfile, "fasta"))
    qcs_records = list(SeqIO.parse(qcsfile, "fasta"))
    hqcs_dict = {r.id : r for r in hqcs_records}
    qcs_dict = {r.id : r for r in qcs_records}

    for hqcs_id, qcs_ids in match_dict.items():
        outfile = os.path.join(outdir, '{}.unaligned.fasta'.format(hqcs_id))
        records = [hqcs_dict[hqcs_id]]
        records.extend(list(qcs_dict[i] for i in qcs_ids))
        SeqIO.write(records, outfile, "fasta")


if __name__ == "__main__":
    qcsfile, hqcsfile, pairfile, outdir = sys.argv[1:]
    main(qcsfile, hqcsfile, pairfile, outdir)
