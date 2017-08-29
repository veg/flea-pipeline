import sys
import fileinput
from collections import defaultdict


def main():
    pairs = list(line.strip().split("\t") for line in fileinput.input())
    hqcs_counts = defaultdict(lambda: 0)
    for pair in pairs:
        if len(pair) != 2:
            warnings.warn('CCS {} did not match any HQCS'.format(pair))
            continue
        ccs_id, hqcs_id = pair
        hqcs_counts[hqcs_id] += 1
    for id_, count in hqcs_counts.items():
        sys.stdout.write('{}\t{}\n'.format(id_, count))


if __name__ == "__main__":
    main()
