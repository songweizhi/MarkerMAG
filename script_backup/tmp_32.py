import os
import glob
import shutil
import argparse
import pandas as pd
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
from distutils.spawn import find_executable


def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc


def sam_flag_to_rc(flag_value):

    binary_flag = "{0:b}".format(int(flag_value))
    binary_flag_len = len(str(binary_flag))
    binary_flag_polished = '0' * (12 - binary_flag_len) + str(binary_flag)

    read_rced = False
    if binary_flag_polished[7] == '1':
        read_rced = True

    return read_rced


class MappingRecord:

    #  sequences store in r1_seq and r2_seq should NOT been reverse complemented

    def __init__(self):

        self.r1_seq = ''
        self.r2_seq = ''
        self.r1_seq_qual = '*'
        self.r2_seq_qual = '*'

        self.r1_refs = dict()
        self.r2_refs = dict()

        self.r1_cigar_to_flag = dict()
        self.r2_cigar_to_flag = dict()

        self.r1_longest_clp_cigar = ''
        self.r1_longest_clp_falg  = ''
        self.r2_longest_clp_cigar = ''
        self.r2_longest_clp_falg  = ''

        self.qualified_reads           = False
        self.consider_round_2          = False
        self.consider_r1_unmapped_mate = False
        self.consider_r1_clipping_part = False
        self.consider_r2_unmapped_mate = False
        self.consider_r2_clipping_part = False

        self.r1_filtered_refs = set()
        self.r2_filtered_refs = set()

        self.r1_clipping_seq = ''
        self.r2_clipping_seq = ''

        self.r1_clipping_seq_qual = '*'
        self.r2_clipping_seq_qual = '*'

        self.unmapped_r1_refs = set()
        self.unmapped_r2_refs = set()

        self.clipping_r1_refs = set()
        self.clipping_r2_refs = set()


def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if (each_element.isalpha() is True) or (each_element == '='):
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


def get_cigar_stats(cigar_splitted):

    # aligned_len: M I X =
    # clipping_len: S
    # mismatch_len: X I D
    # mismatch_pct = mismatch_len / aligned_len
    # aligned_pct  = aligned_len  / (aligned_len + clipping_len)
    # clipping_pct = clipping_len / (aligned_len + clipping_len)

    aligned_len = 0
    clipping_len = 0
    mismatch_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]

        # get aligned_len
        if each_part_cate in {'M', 'm', 'I', 'i', 'X', 'x', '='}:
            aligned_len += each_part_len

        # get clipping_len
        if each_part_cate in ['S', 's']:
            clipping_len += each_part_len

        # get mismatch_len
        if each_part_cate in {'I', 'i', 'X', 'x', 'D', 'd'}:
            mismatch_len += each_part_len

    aligned_pct  = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
    clipping_pct = float("{0:.2f}".format(clipping_len * 100 / (aligned_len + clipping_len)))
    mismatch_pct = float("{0:.2f}".format(mismatch_len * 100 / (aligned_len)))

    return aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct


########################################################################################################################

# file in
pwd_samfile                         = '/Users/songweizhi/Desktop/new_algorithm_Kelp/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.sam'
min_M_pct                           = 20
min_M_len                           = 30
min_clp_len                         = 30
min_clp_M_len                       = 25
max_mis_pct                         = 3

# file out
pwd_samfile_reads                   = '/Users/songweizhi/Desktop/new_algorithm_Kelp/Kelp_SILVA138_id99_assembled_16S_uclust_0.999_reads.fa'
clipping_part_seq                   = '/Users/songweizhi/Desktop/new_algorithm_Kelp/clipping_parts.fa'
unmapped_mates                      = '/Users/songweizhi/Desktop/new_algorithm_Kelp/unmapped_mates.fa'
clipping_part_seq_sam               = '/Users/songweizhi/Desktop/new_algorithm_Kelp/clipping_parts.sam'
unmapped_mates_sam                  = '/Users/songweizhi/Desktop/new_algorithm_Kelp/unmapped_mates.sam'
stats_combined                      = '/Users/songweizhi/Desktop/new_algorithm_Kelp/stats_combined.txt'
stats_combined_filtered             = '/Users/songweizhi/Desktop/new_algorithm_Kelp/stats_combined_filtered.txt'
blast_results_all_vs_all_16s        = '/Users/songweizhi/Desktop/new_algorithm_Kelp/16S_all_vs_all_blastn.tab'
mock_final_op                       = '/Users/songweizhi/Desktop/new_algorithm_Kelp/mock_final_op.txt'
mock_final_op_ctg_level             = '/Users/songweizhi/Desktop/new_algorithm_Kelp/mock_final_op_ctg_level.txt'

marker_to_ctg_Key_connector_str = '___M___'
marker_to_gnm_Key_connector_str = '___M___'
gnm_ctg_connector = '___'


########################################################################################################################

MappingRecord_dict = {}
for each_read in open(pwd_samfile):
    if not each_read.startswith('@'):
        store_read_seq = False
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        read_flag = int(each_read_split[1])
        read_seq = each_read_split[9]
        read_seq_qual = each_read_split[10]
        cigar = each_read_split[5]

        if cigar != '*':
            ref_id = each_read_split[2]
            ref_pos = each_read_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

            if (aligned_len >= min_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= (max_mis_pct)):
                if read_id_base not in MappingRecord_dict:
                    MappingRecord_dict[read_id_base] = MappingRecord()
                if read_strand == '1':
                    MappingRecord_dict[read_id_base].r1_refs[ref_id_with_pos] = cigar
                    MappingRecord_dict[read_id_base].r1_cigar_to_flag[cigar] = read_flag
                    store_read_seq = True
                if read_strand == '2':
                    MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar
                    MappingRecord_dict[read_id_base].r2_cigar_to_flag[cigar] = read_flag
                    store_read_seq = True




            else:
                store_read_seq = True
        else:
            store_read_seq = True

        # store_read_seq into dict
        if store_read_seq is True:
            read_rc = sam_flag_to_rc(read_flag)

            if read_id_base not in MappingRecord_dict:
                MappingRecord_dict[read_id_base] = MappingRecord()

            # turn back if read reverse complemented
            read_seq_to_store = read_seq
            read_seq_qual_to_store = read_seq_qual
            if read_rc is True:
                read_seq_to_store = get_rc(read_seq)
                read_seq_qual_to_store = read_seq_qual[::-1]

            if read_strand == '1':
                MappingRecord_dict[read_id_base].r1_seq = read_seq_to_store
                MappingRecord_dict[read_id_base].r1_seq_qual = read_seq_qual_to_store
            if read_strand == '2':
                MappingRecord_dict[read_id_base].r2_seq = read_seq_to_store
                MappingRecord_dict[read_id_base].r2_seq_qual = read_seq_qual_to_store


        print(each_read_split)
        print(MappingRecord_dict[read_id_base].r1_cigar_to_flag)
        print(MappingRecord_dict[read_id_base].r2_cigar_to_flag)
        print()

