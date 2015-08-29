# TODO: Sharing globals like this tightly couples modules and is bad
# design. However, multiprocessing cannot pickle nested functions, and
# this was the easiest way to un-nest them. Come up with a better way
# in the future.

config = None
options = None
script_dir = None
data_dir = None
qsub_dir = None
timepoints = None
key_to_label = None
logger = None
logger_mutex = None
