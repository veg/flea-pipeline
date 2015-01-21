#!/usr/bin/env python

"""
Publish a run to test.datamonkey.org.

Usage:
  publish [options] <directory> <name>
  publish -h | --help

Options:
  --host=<string>         Hostname [default: test.datamonkey.org]
  -u --username=<string>  Username on target host
  -h --help               Show this screen.

"""

import os
import sys
import fnmatch
import socket
import getpass
import subprocess
import shlex

from docopt import docopt

REMOTE_DIR = "/var/www/html/veg/FLEA/"

RSYNC_CMD = "rsync {src} {remote_user}@{remote_host}:{dest_path}"

HTML_CMD = ("ssh {remote_user}@{remote_host} \"cd {remote_dir} &&"
            " cp test.html {name}.html &&"
            " sed -i 's/test-data/{dest}/g' {name}.html\"")


SRC_FILE = "{local_user}@{local_host}:{path}"


def recursive_glob(directory, glob):
    """Find all files in `directory` that match `glob`."""
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, glob):
            matches.append(os.path.join(root, filename))
    return matches


def call(cmd):
    ssh = subprocess.Popen(shlex.split(cmd),
                           shell=False,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE)
    result = ssh.stdout.readlines()
    if result == []:
        error = ssh.stderr.readlines()
        sys.stderr.write("\n".join(s.decode("utf-8") for s in error))
    else:
        print("\n".join(s.decode("utf-8") for s in result))


if __name__ == "__main__":
    args = docopt(__doc__)
    directory = args["<directory>"]
    name = args["<name>"]
    remote_user = args["--username"]
    remote_host = args["--host"]

    local_user = getpass.getuser()
    local_host = socket.gethostname()

    if remote_user is None:
        remote_user = local_user
    
    files = recursive_glob(directory, "*json")
    files.extend(recursive_glob(directory, "*tsv"))

    src = " ".join(files)
    dest = "{name}-data".format(name=name)
    dest_path = os.path.join(REMOTE_DIR, dest)

    rsync_cmd = RSYNC_CMD.format(src=src, remote_user=remote_user,
                                 remote_host=remote_host,
                                 dest_path=dest_path)
    html_cmd = HTML_CMD.format(remote_user=remote_user,
                               remote_host=remote_host,
                               remote_dir=REMOTE_DIR,
                               name=name, dest=dest)
    r1 = call(rsync_cmd)
    r2 = call(html_cmd)

    print("view results at http://test.datamonkey.org/veg/FLEA/{}.html".format(name))
