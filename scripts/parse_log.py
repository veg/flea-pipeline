#!/usr/bin/env python
"""
Parse a flea log file. Extracts pipeline start events and task
completion events.

Usage:
  parse_log [options] <logfile>
  parse_log -h | --help

Options:
  -h --help  Show this screen.

"""

import re
from datetime import datetime
from dateutil.parser import parse

from docopt import docopt


uuid_regex = r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}'


class Task:
    def __init__(self, name):
        self.name = name
        self.stop = None
        self.elapsed = None
        self.commands = ''
        self.inputs = []
        self.outputs = []

    @property
    def start(self):
        return self.stop - self.elapsed


def first(line):
    return line.split()[0].strip()


def last(line):
    return line.split()[-1].strip()


def after(line, phrase):
    idx = line.find(phrase)
    return line[idx + len(phrase)]


def parse_log(logfile):
    with open(logfile) as f:
        lines = f.read().strip().split('\n')
    starts = []
    tasks = []
    task = None
    for i in range(len(lines)):
        line = lines[i]
        if 'INFO' not in line:
            continue
        if 'pipeline start' in line:
            date = parse(first(lines[i - 1]))
            time = parse(first(line))
            d = datetime(date.year, date.month, date.day, time.hour, time.minute, time.second)
            starts.append(d)
            print('start: {}'.format(d))
        elif 'Finish task' in line:
            if task is not None:
                tasks.append(task)
                print('task: {} {}'.format(task.name, task.elapsed))
            name = last(line)
            if re.search(uuid_regex, name):
                name = name[:-37]
            task = Task(name)
        elif 'Time:' in line:
            t = parse(last(line))
            task.stop = time
        elif 'Elapsed:' in line:
            task.elapsed = parse(last(line))
        elif 'Commands:' in line:
            task.commands = after(line, 'Commands: ')
        elif 'Input:' in line:
            task.inputs.append(after(line, 'Input: '))
        elif 'Output: ' in line:
            task.inputs.append(after(line, 'Output: '))
    if task is not None:
        tasks.append(task)
        print('task: {}'.format(task.name))
    return starts, tasks


if __name__ == '__main__':
    args = docopt(__doc__)
    starts, tasks = parse_log(args['<logfile>'])

