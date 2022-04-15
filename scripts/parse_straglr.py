#! /usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--input', help='the input file which was generated from straglr', required=True)

args = parser.parse_args()

input_file = open(args.input)

print("read_name\tstrand\tcount")

f = open(args.input)
header1 = f.readline()
header2 = f.readline()
header3 = f.readline()
#with open(args.input) as f:
for record in f:
    chromosome,start,end,repeat_unit,genotype,read,cn,size,read_start,allele = record.rstrip().split("\t")
    #print(read)
    print("%s\t+\t%s" % (read, int(float(cn))))
