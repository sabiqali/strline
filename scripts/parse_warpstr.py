#! /usr/bin/env python

import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from straglr', required=True)

args = parser.parse_args()

input_file = open(args.input)

print("read_name\tstrand\tcount")

f = open(args.input)

header = f.readline()

for record in f:
    read_name,sample,run_id,reverse,sam_dist = record.rstrip().split("\t")
    strand = '-' if reverse == "True" else '+'
    print("%s\t%s\t%s" % (read_name,strand,sam_dist))
