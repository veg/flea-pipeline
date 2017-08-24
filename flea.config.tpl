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
fasttree: FastTree

[Jobs]
queue: fast
jobs: 4
ppn: 4
walltime: 3600
long_walltime: 86400

[Parameters]

# references
reference_dir: /path/to/references
reference_protein: ${reference_dir}/reference.aa.fasta
reference_coordinates: ${reference_dir}/reference.numbers.csv
reference_db: ${reference_dir}/reference.database.fasta
contaminants_db: ${reference_dir}/contaminants.database.fasta

# quality pipeline
min_sequence_length: .7
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

use_rifraf: True

# only used if use_rifraf is false
do_shift_correction: True

# if using rifraf
consensus_multiprocess: False
consensus_max_iters: 10
seq_errors: "1,3,3,0,0"
ref_errors: "8,0.1,0.1,1,1"
phred_cap: 30
hqcs_max_err_rate: 1e-6
hqcs_max_base_err_rate: 1e-3

copynumber_identity: 0.95
cn_max_length_ratio: 1.005

# diagnosis pipeline
ccs_to_hqcs_identity: 0.95
max_ccs_length_ratio: 1.005

[Tasks]
quality_pipeline: True

# if True, use all in-frame HQCS sequences as reference DB for shift
# correction or RIFRAF references.
use_inframe_db: True

# pause to edit the database
pause_for_inframe_db: False

# pause to correct codon alignment
pause_after_codon_alignment: False

# run analysis subpipeline
analysis: True

# hyphy scripts
hyphy_analysis: True
evo_hisotry: True
fubar: True

# full CCS alignment, SMD metrics
diagnose: False


[Misc]
# when a task gets re-run
# level 0: file timestamps
# level 1: + history
# level 2: + function bodies
# level 3: + function parameters
checksum_level: 1

