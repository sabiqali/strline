#! /usr/bin/env python

import argparse
import os

class GraphAlignment:
    def __init__(self, read_name, strand, spanned, count, alignment_score, identity, aligned_fraction):
        self.read_name = read_name
        self.strand = strand
        self.spanned = spanned
        self.count = count
        self.alignment_score = alignment_score
        self.identity = identity
        self.aligned_fraction = aligned_fraction

parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from graphaligner', required=True)
parser.add_argument('--min-identity', type=float, default=0.50, help='only use reads with identity greater than this', required=False)
parser.add_argument('--min-aligned-fraction', type=float, default=0.8, help='require alignments cover this proportion of the query sequence', required=False)
parser.add_argument('--write-non-spanned', action='store_true', default=False, help='do not require the reads to span the prefix/suffix region', required=False)
args = parser.parse_args()

input_file = open(args.input)

alignments = dict()
with open(args.input) as f:
    for record in f:
        fields = record.rstrip().split("\t")
        read_id = fields[0].split(' ')[0] # remove FASTQ metadata that graphaligner emits

        tags = dict()
        for t in fields[20:]:
            try:
                (key, data_type, value) = t.split(":")
                tags[key] = value
            except ValueError:
                key = t.split(":")[0]
                tags[key] = 0

        align_score = float(tags["AS"])
        identity = float(tags["id"])

        query_len = int(fields[9])
        query_start = int(fields[10])
        query_end = int(fields[11])
        query_af = float(query_end - query_start) / float(query_len)

        path_str = fields[13]
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
        strand = fields[12]
        assert(strand == "+")
        if path_dir == "<":
            strand = "-"
        ga = GraphAlignment(read_id, strand, valid, count, align_score, identity, query_af)
        if (valid or args.write_non_spanned) and ga.identity > args.min_identity and ga.aligned_fraction > args.min_aligned_fraction:
            if ga.read_name not in alignments or alignments[ga.read_name].alignment_score < ga.alignment_score:
                alignments[ga.read_name] = ga

print("\t".join(["read_name","strand", "spanned", "count", "align_score", "identity", "query_aligned_fraction"]))
for ga in alignments.values():
    print("%s\t%s\t%d\t%d\t%.1f\t%.3f\t%.3f" % (ga.read_name, ga.strand, ga.spanned, ga.count, ga.alignment_score, ga.identity, ga.aligned_fraction))
