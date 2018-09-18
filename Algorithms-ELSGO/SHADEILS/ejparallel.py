#!/usr/bin/python
from __future__ import print_function
import sys
import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Prepare the config file for run in server")

parser.add_argument("-n", required=True, type=str, dest="name", help='Name')
parser.add_argument("-e", required=True, type=str, dest="filename",
                    help='Executable name')
parser.add_argument("-o", dest="output", help='Output name')
 
args = parser.parse_args()

runfile = args.filename

if not os.path.isfile(runfile):
    parser.print_help()
    sys.exit(1)

if args.output:
    file_output = open(args.output, 'w')
else:
    file_output = sys.stdout

funs = range(1, 15+1)
total = len(funs)
head = """\
#!/bin/bash
#$ -N {name}
#$ -q media
#$ -o out_{name}.txt
#$ -e err_{name}.txt
#$ -cwd
#$ -t 1-{total}\
""".format(name=args.name, total=total)

print(head, file=file_output)
count = 1

for fun in funs:
    params = "PARAMS[{}]='-f {} -s 2 -r 25'".format(count, fun)
    print(params, file=file_output)
    count += 1

print('./', runfile, ' ${PARAMS[$SGE_TASK_ID]}', sep='', file=file_output)
file_output.close()
