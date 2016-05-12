[Paths]
bealign: bealign.py
python: python
julia: julia
consensus_script: /path/to/script.jl

# the version of bam2msa in tristan's fork of BioExt is outdated.
bam2msa: /opt/python-3.4.3/bin/bam2msa

mafft: mafft
usearch: usearch8
hyphy: HYPHYMP
FastTree: FastTree

[Jobs]
# maximum number of concurrent jobs
local_jobs: 1

# maximum jobs submitted to the cluster
remote_jobs: 4

use_cluster: True
queue: fast
nodes: 1
ppn: 16
walltime: 3600
hyphy_walltime: 86400

[Parameters]
reference_dir: /path/to/references
reference_sequence: ${reference_dir}/reference.translated.fasta
reference_coordinates: ${reference_dir}/reference.numbers.csv
reference_db: ${reference_dir}/reference.database.fasta
contaminants_db: ${reference_dir}/contaminants.database.fasta
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
ccs_to_hqcs_identity: 0.95
min_cluster_size: 5

# for clustering
min_length_ratio: 0.995

consensus_batch_size: 30
consensus_p_ins: 0.01
consensus_p_del: 0.01

# for copy number
max_query_target_length_ratio: 1.005

# for run filtering
run_length: 16

[Tasks]
# if True, use all in-frame HQCS sequences as reference DB for shift
# correction of HQCSs. Also pause to allow user to edit database.
use_inframe_hqcs: False
pause_for_inframe_hqcs: False

pause_after_codon_alignment: False
align_ccs: False
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

