#!/usr/bin/env python

"""
Publish a run to test.datamonkey.org.

Usage:
  publish [options] <directory> [<target>]
  publish -h | --help

Options:
  --dry-run  Print command and exit.
  -h --help  Show this screen.

"""

import os
import sys
import fnmatch
import socket
import getpass
import subprocess
import shlex

from docopt import docopt

DEFAULT_HOST = "test.datamonkey.org"
DEFAULT_PATH = "/var/www/html/veg/FLEA3"

RSYNC_CMD = "rsync {src} {remote_user}@{remote_host}:{dest_path}"


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


def strip_trailing_slash(directory):
    if directory[-1] == "/":
        directory = directory[:-1]
    return directory


if __name__ == "__main__":
    args = docopt(__doc__)
    directory = args["<directory>"]
    directory = strip_trailing_slash(directory)

    local_user = getpass.getuser()
    local_host = socket.gethostname()

    target = args["<target>"]
    if target is None:
        target = "{}@{}:{}".format(local_user, DEFAULT_HOST, DEFAULT_PATH)
    if "@" in target:
        remote_user, target = target.split("@")
    else:
        remote_user= local_user
    if ":" in target:
        remote_host, target = target.split(":")
    else:
        remote_host = DEFAULT_HOST

    if target:
        remote_path = target
    else:
        remote_path = DEFAULT_PATH

    if not os.path.exists(directory):
        raise Exception('source directory does not exist: {}'.format(directory))
    files = recursive_glob(directory, "*json")
    files.extend(recursive_glob(directory, "copynumbers.tsv"))
    if not files:
        raise Exception('no files found')
    files = list(os.path.abspath(f) for f in files)
    src = " ".join(files)

    _, dir_name = os.path.split(os.path.abspath(directory))
    dir_name = strip_trailing_slash(dir_name)
    dest_path = "{}/".format(os.path.join(remote_path, dir_name))

    rsync_cmd = RSYNC_CMD.format(src=src, remote_user=remote_user,
                                 remote_host=remote_host,
                                 dest_path=dest_path)
    if args["--dry-run"]:
        print(rsync_cmd)
    else:
        result = call(rsync_cmd)
        print("view results at http://test.datamonkey.org:5062/{}/".format(dir_name))
