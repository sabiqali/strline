#! /usr/bin/env python

import functools
import argparse
import numpy as np
import sys
import csv

class STRCall:
    def __init__(self, strand, count):
        self.strand = strand
        self.count = count

def load_strique(filename):
    out = dict()
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            out[row['ID']] = STRCall(row['strand'], row['count'])
    return out

def load_graphaligner(filename):
    out = dict()
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            read_id = row['read_name']
            if int(row['spanned']) == 1 and read_id not in out:
                out[read_id] = STRCall(row['strand'], row['count'])
    return out

def load_strscore(filename):
    out = dict()
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            out[row['read_name']] = STRCall(row['strand'], row['count'])
    return out

parser = argparse.ArgumentParser()
parser.add_argument('--strique', required=False)
parser.add_argument('--strscore', required=False)
parser.add_argument('--graphaligner', required=False)
parser.add_argument('--read-ids', required=False)
args = parser.parse_args()

data = dict()
if args.strique is not None:
    data['strique'] = load_strique(args.strique)

if args.graphaligner is not None:
    data['graphaligner'] = load_graphaligner(args.graphaligner)

if args.strscore is not None:
    data['strscore'] = load_strscore(args.strscore)

all_reads = list()
if args.read_ids is not None:
    with open(args.read_ids) as f:
        for row in f:
            all_reads.append(row.rstrip())
else:
    all_reads = functools.reduce(set.union, (set(d.keys()) for d in data.values()))
read_set = set(all_reads)

methods = list(data.keys())
print("\t".join(["read_id", "strand"] + methods))
for read_id in all_reads:
    if read_id not in read_set:
        continue

    out = [ read_id ]

    # determine the consensus strand
    strand_count = { "+":0, "-":0 }
    for m in methods:
        if read_id in data[m]:
            strand_count[data[m][read_id].strand] += 1

    consensus_strand = "?"
    if strand_count["+"] > 0 and strand_count["-"] == 0:
        consensus_strand = "+"
    elif strand_count["-"] > 0 and strand_count["+"] == 0:
        consensus_strand = "-"
    
    out.append(consensus_strand)

    for m in methods:
        v = "NA"
        if read_id in data[m]:
            v = data[m][read_id].count
        out.append(v)

    print("\t".join(out))
