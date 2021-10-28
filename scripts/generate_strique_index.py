#! /usr/bin/env python

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument('--fast5_dir', type=str, required=True)
parser.add_argument('--summary', type=str, required=True)
args = parser.parse_args()

with open(args.summary) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for record in reader:
        print("%s/%s/read_%s\t%s" % (args.fast5_dir, record['filename'], record['read_id'], record['read_id']))
