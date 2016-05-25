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
queue: fast
jobs: 4
ppn: 4
walltime: 3600
hyphy_walltime: 86400

[Parameters]

# references
reference_dir: /path/to/references
reference_dna: ${reference_dir}/reference.dna.fasta
reference_protein: ${reference_dir}/reference.aa.fasta
reference_coordinates: ${reference_dir}/reference.numbers.csv
reference_db: ${reference_dir}/reference.database.fasta
contaminants_db: ${reference_dir}/contaminants.database.fasta

# quality pipeline
min_sequence_length: .8
max_err_rate: 0.01
qmax: 55
max_accepts: 300
max_rejects: 600
contaminant_identity: 0.98
reference_identity: 0.8
run_length: 16

# consensus pipeline
cluster_identity: 0.99
min_cluster_size: 3
max_cluster_size: 100
min_length_ratio: 0.995
# obtainEvolutionaryHistory.bf fails if there are fewer than 3
# consensus sequences per time point.
min_n_clusters: 3

consensus_batch_size: 10
consensus_p_ins: 0.01
consensus_p_del: 0.01

copynumber_identity: 0.95
cn_max_length_ratio: 1.005

# diagnosis pipeline
ccs_to_hqcs_identity: 0.95
max_ccs_length_ratio: 1.005

[Tasks]
quality_pipeline: True
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

