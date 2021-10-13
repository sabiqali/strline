#! /usr/bin/env python

import functools
import argparse
import numpy as np
import sys
import csv

def load_strique(filename):
    out = dict()
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            out[row['ID']] = row['count']
    return out

def load_graphaligner(filename):
    out = dict()
    with open(filename) as f:
        for row in f:
            (read_id, count, strand, spanned) = row.rstrip().split()
            if read_id not in out:
                out[read_id] = count
    return out

def load_strscore(filename):
    out = dict()
    with open(filename) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            out[row['read_name']] = row['count']
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
print("\t".join(["read_id"] + methods))
for read_id in all_reads:
    if read_id not in read_set:
        continue

    out = [ read_id ]
    for m in methods:
        v = "NA"
        if read_id in data[m]:
            v = data[m][read_id]
        out.append(v)

    print("\t".join(out))

sys.stderr.write("method\tnum_reads\tnum_called\tnum_not_zero\tmean\tstdv\n")
for m in methods:
    called = [ int(data[m][i]) for i in data[m] if i in read_set ]
    non_zero = [ c for c in called if c > 0 ]
    mean = np.mean(non_zero)
    stdv = np.std(non_zero)
    sys.stderr.write("%s\t%d\t%d\t%d\t%.2lf\t%.2lf\n" % (m, len(all_reads), len(called), len(non_zero), mean, stdv))
