#! /usr/bin/env python
import functools
import argparse
import numpy as np
import sys
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
args = parser.parse_args()

methods = None
data = dict()

with open(args.input) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        if methods is None:

            # initialize
            methods = list(row.keys())[2:]
            for m in methods:
                data[m] = list()
        
        for m in methods:
            data[m].append(row[m])        

sys.stderr.write("method\tnum_reads\tnum_called\tnum_not_zero\tmean\tstdv\n")
for m in methods:
    called = [ int(x) for x in data[m] if x != "NA" ]
    non_zero = [ c for c in called if c > 0 ]
    mean = np.mean(non_zero)
    stdv = np.std(non_zero)
    sys.stderr.write("%s\t%d\t%d\t%.2lf\t%.2lf\n" % (m, len(called), len(non_zero), mean, stdv))

