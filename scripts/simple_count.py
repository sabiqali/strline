#! /usr/bin/env python

import argparse
import pysam
import sys
import csv

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def count_motif(sequence, motif):
    
    motif_len = len(motif)
    match_vector = list()
    count = 0
    for i in range(0, len(sequence) - motif_len + 1):
        match = int(read.sequence[i:i+motif_len] == motif)
        match_vector.append(match)
        count += match
    
    # cluster matches together
    clusters = list()
    
    d = args.max_distance
    current_cluster = list()
    for (index, is_match) in enumerate(match_vector):
        if is_match:
            if len(current_cluster) == 0 or index - current_cluster[-1] > d:
                # end cluster, start new one
                if len(current_cluster) > 0:
                    clusters.append(current_cluster)
                current_cluster = [ index ]
            else:
                current_cluster.append(index)
    if len(current_cluster) > 0:
        clusters.append(current_cluster)

    max_size = 0
    for c in clusters:
        #print(",".join([str(x) for x in c]))
        if len(c) > max_size:
            max_size = len(c)
    return max_size

parser = argparse.ArgumentParser()
parser.add_argument('--motif', type=str, required=False)
parser.add_argument('--config', type=str, required=False)
parser.add_argument('--max-distance', type=int, required=False, default=25)
args, input_file = parser.parse_known_args()

motif_fwd = None
if args.config is not None:
    with open(args.config) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for record in reader:
            motif_fwd = record['repeat']
            break
elif args.motif is not None:
    motif_fwd = args.motif
else:
    sys.stderr.write("One of --config or --motif must be provided")
    sys.exit(1)

motif_rev = reverse_complement(motif_fwd)
motif_len = len(motif_fwd)

print("read_name\tstrand\tcount\tfwd_count\trev_count")
for read in pysam.FastxFile(input_file[0]):

    fwd_count = count_motif(read.sequence, motif_fwd)
    rev_count = count_motif(read.sequence, motif_rev)

    count = None
    strand = '?'
    if fwd_count > rev_count:
        count = fwd_count
        strand = '+'
    else:
        count = rev_count
        strand = '-'

    print("%s\t%s\t%d\t%d\t%d" % (read.name, strand, count, fwd_count, rev_count))
