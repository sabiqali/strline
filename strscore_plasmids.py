#! /usr/env/python

import parasail
import sys
import pysam
import argparse
import math

class ReadAlignment:
    def __init__(self, name, strand):
        self.read_name: name
        self.has_prefix_match = False
        self.has_suffix_match = False
        self.strand = strand

    count = 0
    bad_mapping = 0
    prefix_start = 0
    suffix_end = 0
    prefix_length = 0
    suffix_length = 0

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def percentage_identity(cigar_exp):
    m = 0
    nm = 0
    for chr in cigar_exp:
        if chr == '|':
            m = m + 1
        else:
            nm = nm + 1
    pi = float(m/(nm + m))
    return pi

def roundup(x):
    return int(math.ceil(x / 100000.0)) * 100000

def alignment_contains_str_prefix(alignment,start):
    pair_out = alignment.get_aligned_pairs(True)
    c = 0
    for tmp_pairs in  pair_out:
        if tmp_pairs[1] < int(start) and tmp_pairs[1] >= (int(start) - 100):
            c = c + 1
    #result = parasail.sw_trace_scan_32(flank, alignment.query_sequence, 5, 4, scoring_matrix)
    if((c/100) > 0.6):
        return True
    else: 
        return False

def alignment_contains_str_suffix(alignment,end):
    pair_out = alignment.get_aligned_pairs(True)
    c = 0
    for tmp_pairs in  pair_out:
        if tmp_pairs[1] > int(end) and tmp_pairs[1] <= (int(end) + 100):
            c = c + 1
    #result = parasail.sw_trace_scan_32(flank, alignment.query_sequence, 5, 4, scoring_matrix)
    if((c/100) > 0.6):
        return True
    else: 
        return False

def get_alignment_points_prefix(alignment,start,flength):
    pair_out = alignment.get_aligned_pairs(True)
    c = 0
    last_index = int(start)
    for i in range(int(start) - flength, int(start)):
        for tmp_pairs in pair_out:
            if tmp_pairs[1] == i:
                last_index = tmp_pairs[0]
                c= c + 1
    start_index = last_index - c
    return (start_index, last_index)

def get_alignment_points_suffix(alignment,end,flength):
    pair_out = alignment.get_aligned_pairs(True)
    c = 0
    last_index = int(end)
    for i in range(int(end), int(end) + flength):
        for tmp_pairs in pair_out:
            if tmp_pairs[1] == i:
                last_index = tmp_pairs[0]
                c= c + 1
    start_index = last_index - c
    return (start_index, last_index)

def get_alignment_points(alignment,start,end):
    pair_out = alignment.get_aligned_pairs(True)
    prefix_indexes = []
    suffix_indexes = []
    c = 0
    for tmp_pairs in pair_out:
        if tmp_pairs[1] < int(start) and tmp_pairs[1] >= (int(start) - 100):
            prefix_indexes.append(tmp_pairs[0])
        if tmp_pairs[1] > int(end) and tmp_pairs[1] <= (int(end) + 100):
            suffix_indexes.append(tmp_pairs[0])
        c = c + 1
    #return (prefix_indexes[0],suffix_indexes[(len(suffix_indexes) - 1) if len(suffix_indexes) != 0 else 0])
    return (prefix_indexes,suffix_indexes)

def get_alignment_points_parasail(read,flank,matrix):
    result = parasail.ssw(flank,read,5,4,matrix)
    return (result.ref_begin1, result.ref_end1)

parser = argparse.ArgumentParser()
parser.add_argument('--bam', help='the bam file', required=False)
parser.add_argument('--read', help='the read file', required=False)
parser.add_argument('--ref', help='the ref file', required=True)
parser.add_argument('--config', help='the config file', required=True)
parser.add_argument('--output', help='the output file', required=False)
parser.add_argument('--verbose', help='display the alignments of the different regions', type=int, required=False, default=0)

args = parser.parse_args()

reads_file = args.read
reference_file = args.ref
config = args.config
in_bam = args.bam

# read in the configs and store them in the respective variables
configs = list()
config_fh = open(config)
header = config_fh.readline()

for line in config_fh:
    configs.append(line.rstrip().split()[0:7])

chromosome,begin,end,name,repeat,prefix,suffix = configs[0]
    
for read in pysam.FastxFile(reference_file):
    if(read.name == chromosome):
        ref_seq = read.sequence

max_repeats = 250

#First pass through the alignment file to determine what are the matches.
bamfile = pysam.AlignmentFile(in_bam)
upper_limit = int(end) 
lower_limit = int(begin)
idx = 0
scoring_matrix = parasail.matrix_create("ACGT", 5, -1)
reads = dict()
for alignment in bamfile.fetch(chromosome,lower_limit,upper_limit):
    strand = '-' if alignment.is_reverse else '+'
    if alignment.qname not in reads:  
        reads[alignment.qname] = ReadAlignment(alignment.qname,strand)
        #reads[alignment.qname].strand == '-' if alignment.is_reverse else '+'
    if reads[alignment.qname].strand != '' and reads[alignment.qname].strand != strand:
        reads[alignment.qname].bad_mapping = 1
        #print("test")
    if alignment_contains_str_prefix( alignment , begin):
        #print("test1")
        reads[alignment.qname].has_prefix_match = True
    if alignment_contains_str_suffix( alignment, end):
        #print("test2")
        reads[alignment.qname].has_suffix_match = True
    #(prefix_start_tmp,suffix_end_tmp) = get_alignment_points(alignment,begin,end)  #This is the usual way, trying something new with the new function
    #if prefix_start_tmp:
    #    reads[alignment.qname].prefix_start = prefix_start_tmp[0]
    #if suffix_end_tmp:
    #    reads[alignment.qname].suffix_end = suffix_end_tmp[-1]
    (prefix_start_tmp, prefix_end_tmp) = get_alignment_points_prefix(alignment,begin, len(prefix))
    reads[alignment.qname].prefix_length = prefix_end_tmp - prefix_start_tmp
    (suffix_start_tmp, suffix_end_tmp) = get_alignment_points_suffix(alignment, end, len(suffix))
    reads[alignment.qname].suffix_length = suffix_end_tmp - suffix_start_tmp
    reads[alignment.qname].prefix_start = prefix_start_tmp
    reads[alignment.qname].suffix_end = suffix_end_tmp
 
#Once we have all the matches, we can iterate through them to get the count
print("\t".join(["read_name","chromosome","repeat_name","count","strand","aligned_query","aligned_ref"]))
fh = pysam.FastaFile(reads_file)
for read_name , alignment in reads.items():
    ideal_read = ""
    #print("test2")
    if not reads[read_name].has_prefix_match or not reads[read_name].has_suffix_match:
        continue
    if reads[read_name].bad_mapping == 1:
        continue
    #if reads[read_name].strand == '-':
    #    repeat_unit = reverse_complement(repeat)
    #    suffix_unit = reverse_complement(prefix)
    #    prefix_unit = reverse_complement(suffix)
    #else:
    #print("test")
    repeat_unit = repeat
    prefix_unit = prefix
    suffix_unit = suffix
    read_seq = fh.fetch(read_name)

    max_length = len(ref_seq) + (max_repeats * len(repeat_unit))
    if len(read_seq) > max_length:
        continue

    prefix_start_tmp = reads[read_name].prefix_start
    suffix_end_tmp = reads[read_name].suffix_end

    #(prefix_start_tmp,prefix_end_tmp) = get_alignment_points_parasail(read_seq,prefix,scoring_matrix)
    #(suffix_start_tmp,suffix_end_tmp) = get_alignment_points_parasail(read_seq,suffix,scoring_matrix)

    if(prefix_start_tmp < suffix_end_tmp):
        repeat_region_with_flanks = read_seq[prefix_start_tmp:suffix_end_tmp]
    else:
        continue

    if(not len(repeat_region_with_flanks)):
        continue

    #print("test1")
    prev_score = 0 
    c = 0
    ideal_read = prefix_unit[100 - reads[read_name].prefix_length :  ] + ( repeat_unit * c ) + suffix_unit[: reads[read_name].suffix_length - 1]
    #print("test4")
    #print(ideal_read)
    #print(repeat_region_with_flanks)
    #try:
    result = parasail.nw_trace_scan_32(repeat_region_with_flanks, ideal_read, 5, 4, scoring_matrix)
    #except Exception:
    #    continue
    #print("test5")
    prev_result_ref = result.traceback.ref
    result_ref = result.traceback.ref
    result_comp = result.traceback.comp
    prev_result_query = result.traceback.query
    result_query = result.traceback.query
    score = result.score
    prev_score = score
    #print("test2")
    while (score >= prev_score):
        c = c + 1
        prev_score = score
        prev_result_ref = result_ref
        prev_result_query = result_query
        ideal_read = prefix_unit + ( repeat_unit * c ) + suffix_unit
        result = parasail.nw_trace_scan_32(repeat_region_with_flanks, ideal_read, 5, 4, scoring_matrix)
        score = result.score
        result_comp = result.traceback.comp
        result_ref = result.traceback.ref
        result_query = result.traceback.query
    max_score = prev_score
    reads[read_name].count = c - 1
    #print("test3")
    #if read_name == "2341b29f-7f5d-4496-8185-90bbe539251a" or read_name == "98b0e60a-0b39-4e44-9a8b-8fc78d53bc72":
    #    print(read_name)
    #    print(read_seq)
    #    print(c-1)
    #    print(len(repeat))
    #    print(len(prefix))
    #    print(len(suffix))
    #    print(len(read_seq))
    
    print( "\t".join([read_name,chromosome,name,str(reads[read_name].count),reads[read_name].strand,prev_result_query,prev_result_ref,str(reads[read_name].prefix_start),str(reads[read_name].suffix_end)]))
