[Paths]
mafft: /usr/local/bin/mafft
usearch: /usr/local/bin/usearch8
prinseq: prinseq
bealign: /opt/talign/env/bin/bealign.py
bam2msa: /opt/talign/env/bin/bam2msa.py
hyphy: HYPHYMP
FastTree: FastTree

[Jobs]
# maximum number of concurrent jobs
local_jobs: 1

# maximum jobs submitted to the cluster
remote_jobs: 64

use_cluster: True
queue: fast
nodes: 1
ppn: 16
walltime: 86400

[Parameters]
reference_dir: /path/to/references
reference_sequence: ${reference_dir}/reference.translated.fasta
reference_coordinates: ${reference_dir}/reference.numbers.csv
# TODO: generate database files on the fly
reference_db: ${reference_dir}/reference.database.udb
contaminants_db: ${reference_dir}/contaminants.database.udb
min_sequence_length: .8
max_err_rate: 0.01
qmax: 55
max_accepts: 300
max_rejects: 600

# obtainEvolutionaryHistory.bf fails if there are fewer than 3
# consensus sequences per time point.
min_n_clusters: 3
contaminant_identity: 0.98
reference_identity: 0.8
cluster_identity: 0.99
raw_to_consensus_identity: 0.95
min_cluster_size: 5
max_cluster_size: 30

[Tasks]
pause_after_codon_alignment: False
align_full: False
# if False, stop after generating alignment
analysis: True
# if True, run evolutionary history and FUBAR scripts
hyphy_analysis: True

[Misc]
# when a task gets re-run
# level 0: file timestamps
# level 1: + history
# level 2: + function bodies
# level 3: + function parameters
checksum_level: 1
