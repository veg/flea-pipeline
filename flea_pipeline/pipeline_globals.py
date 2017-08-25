# TODO: Sharing globals like this tightly couples modules and is bad
# design. However, multiprocessing cannot pickle nested functions, and
# this was the easiest way to un-nest them. Come up with a better way
# in the future.

config = None
options = None
script_dir = None
data_dir = None
output_dir = None
results_dir = None
job_script_dir = None
timepoints = None
key_to_label = None
label_to_date = None
logger = None
logger_mutex = None
drmaa_session = None
run_locally = False
ppn = 1
ppn_large = 1
