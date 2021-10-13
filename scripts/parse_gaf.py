#! /usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from graphaligner', required=True)
args = parser.parse_args()

input_file = open(args.input)
print("\t".join(["read_name","strand", "spanned", "count","align_score"]))

with open(args.input) as f:
    for record in f:
        fields = record.rstrip().split("\t")
        read_id = fields[0].split(' ')[0] # remove FASTQ metadata that graphaligner emits
        path_str = fields[5]
        align_score = fields[13].split(":")[2]
        path_dir = path_str[0]
        segments = path_str.split(path_dir)
        has_prefix = False
        has_suffix = False

        count = 0
        for s in segments:
            if s.find("prefix") >= 0:
                has_prefix = True
            if s.find("suffix") >= 0:
                has_suffix = True
            count += s.find("repeat") >= 0

        valid = has_prefix and has_suffix
        strand = fields[4]
        assert(strand == "+")
        if path_dir == "<":
            strand = "-"
        print("%s\t%s\t%d\t%d\t%s" % (read_id, strand, valid, count, align_score))
