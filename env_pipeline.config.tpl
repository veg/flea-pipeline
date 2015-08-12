[Paths]
mafft: /usr/local/bin/mafft
usearch: /usr/local/bin/usearch8
prinseq: prinseq
bealign: /opt/talign/env/bin/bealign.py
bam2msa: /opt/talign/env/bin/bam2msa.py
hyphy: HYPHYMP

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
contaminants_db: ${reference_dir}/ContaminantRef.udb
reference_db: ${reference_dir}/LANLsubtypeRef.udb
min_sequence_length: 2200
max_sequence_length: 2900
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
min_orf_length: 750

[Tasks]
pause_after_codon_alignment: False
align_full: False
hyphy: True
fasttree: True

[Misc]
# when a task gets re-run
# level 0: file timestamps
# level 1: + history
# level 2: + function bodies
# level 3: + function parameters
checksum_level: 1

