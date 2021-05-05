import os
import glob
import shutil
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
import multiprocessing as mp
from datetime import datetime
import plotly.graph_objects as go
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from distutils.spawn import find_executable


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


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


class MappingRecord:

    #  sequences store in r1_seq and r2_seq should NOT been reverse complemented

    def __init__(self):

        self.r1_seq = ''
        self.r2_seq = ''

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

        self.unmapped_r1_refs = set()
        self.unmapped_r2_refs = set()

        self.clipping_r1_refs = set()
        self.clipping_r2_refs = set()

        self.unmapped_r1_refs_with_pos = set()
        self.unmapped_r2_refs_with_pos = set()

        self.clipping_r1_refs_with_pos = set()
        self.clipping_r2_refs_with_pos = set()


def sam_flag_to_rc(flag_value):

    read_rced = 'na'
    if flag_value != '':
        binary_flag = "{0:b}".format(int(flag_value))
        binary_flag_len = len(str(binary_flag))
        binary_flag_polished = '0' * (12 - binary_flag_len) + str(binary_flag)

        if binary_flag_polished[7] == '0':
            read_rced = False
        if binary_flag_polished[7] == '1':
            read_rced = True

    return read_rced


def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc


def get_max_clp_and_index(r1_cigar_list, r2_cigar_list):

    r1_cigar_list_split = [cigar_splitter(i) for i in r1_cigar_list]
    r2_cigar_list_split = [cigar_splitter(i) for i in r2_cigar_list]

    r1_cigar_list_split_only_clp = []
    for each_r1_cigar_split in r1_cigar_list_split:
        clp_len_l = 0
        if each_r1_cigar_split[0][-1] in ['S', 's']:
            clp_len_l = int(each_r1_cigar_split[0][:-1])
        clp_len_r = 0
        if each_r1_cigar_split[-1][-1] in ['S', 's']:
            clp_len_r = int(each_r1_cigar_split[-1][:-1])
        r1_cigar_list_split_only_clp.append([clp_len_l, clp_len_r])

    r2_cigar_list_split_only_clp = []
    for each_r2_cigar_split in r2_cigar_list_split:
        clp_len_l = 0
        if each_r2_cigar_split[0][-1] in ['S', 's']:
            clp_len_l = int(each_r2_cigar_split[0][:-1])
        clp_len_r = 0
        if each_r2_cigar_split[-1][-1] in ['S', 's']:
            clp_len_r = int(each_r2_cigar_split[-1][:-1])
        r2_cigar_list_split_only_clp.append([clp_len_l, clp_len_r])

    cigar_list_split_only_clp_r1_r2 = [r1_cigar_list_split_only_clp, r2_cigar_list_split_only_clp]

    max_value = 0
    max_value_index = ''
    for num_list_1 in cigar_list_split_only_clp_r1_r2[0]:
        if num_list_1[0] > max_value:
            max_value = num_list_1[0]
            max_value_index = 'r1_l'
        if num_list_1[1] > max_value:
            max_value = num_list_1[1]
            max_value_index = 'r1_r'
    for num_list_2 in cigar_list_split_only_clp_r1_r2[1]:
        if num_list_2[0] > max_value:
            max_value = num_list_2[0]
            max_value_index = 'r2_l'
        if num_list_2[1] > max_value:
            max_value = num_list_2[1]
            max_value_index = 'r2_r'

    # get the best cigar
    best_cigar = ''
    if max_value_index == 'r1_l':
        if best_cigar == '':
            for each_cigar in r1_cigar_list:
                if (each_cigar.startswith('%sS' % max_value)) or (each_cigar.startswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r1_r':
        if best_cigar == '':
            for each_cigar in r1_cigar_list:
                if (each_cigar.endswith('%sS' % max_value)) or (each_cigar.endswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r2_l':
        if best_cigar == '':
            for each_cigar in r2_cigar_list:
                if (each_cigar.startswith('%sS' % max_value)) or (each_cigar.startswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r2_r':
        if best_cigar == '':
            for each_cigar in r2_cigar_list:
                if (each_cigar.endswith('%sS' % max_value)) or (each_cigar.endswith('%ss' % max_value)):
                    best_cigar = each_cigar

    return best_cigar, max_value, max_value_index


wd = '/Users/songweizhi/Desktop/tunning'
#input_reads_to_16s_sam_best_match       = '%s/GI_0503_mis1_75_45_input_reads_to_16S_best_match.sam' % wd
input_reads_to_16s_sam_best_match       = '%s/sub.sam'                                                             % wd
sam_best_match_unmapped_mates_seq_file  = '%s/GI_0504_mis2_75_45_input_reads_to_16S_best_match_unmapped_mates.fa'   % wd
unmapped_mates_seq_file                 = '%s/unmapped_mates.fa'                                                    % wd
clipping_parts_seq_file                 = '%s/clipping_parts.fa'                                                    % wd
min_M_len_16s                           = 75
min_M_len_ctg                           = 45
min_M_pct                               = 35
mismatch_cutoff                         = 2


MappingRecord_dict = {}
marker_len_dict = {}
for each_read in open(input_reads_to_16s_sam_best_match):
    each_read_split = each_read.strip().split('\t')

    if each_read.startswith('@'):
        marker_id = ''
        marker_len = 0
        for each_element in each_read_split:
            if each_element.startswith('SN:'):
                marker_id = each_element[3:]
            if each_element.startswith('LN:'):
                marker_len = int(each_element[3:])
        marker_len_dict[marker_id] = marker_len
    else:
        store_read_seq = False
        read_id = each_read_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        read_flag = int(each_read_split[1])
        read_seq = each_read_split[9]
        cigar = each_read_split[5]

        if cigar != '*':
            ref_id = each_read_split[2]
            ref_pos = each_read_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
            both_ends_clp = check_both_ends_clipping(cigar_splitted)

            if (aligned_pct >= min_M_pct) and (mismatch_pct <= mismatch_cutoff) and (both_ends_clp is False):
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

            # turn back if read reverse complemented
            read_rc = sam_flag_to_rc(read_flag)
            read_seq_to_store = read_seq
            if read_rc is True:
                read_seq_to_store = get_rc(read_seq)

            if read_id_base not in MappingRecord_dict:
                MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if MappingRecord_dict[read_id_base].r1_seq == '':
                    MappingRecord_dict[read_id_base].r1_seq = read_seq_to_store
            if read_strand == '2':
                if MappingRecord_dict[read_id_base].r2_seq == '':
                    MappingRecord_dict[read_id_base].r2_seq = read_seq_to_store

# add sequences of unmapped mates to mp dict
for each_read in SeqIO.parse(sam_best_match_unmapped_mates_seq_file, 'fasta'):
    read_id = str(each_read.id)
    read_basename = '.'.join(read_id.split('.')[:-1])
    read_strand = read_id.split('.')[-1]

    if read_basename in MappingRecord_dict:
        if read_strand == '1':
            MappingRecord_dict[read_basename].r1_seq = str(each_read.seq)
        if read_strand == '2':
            MappingRecord_dict[read_basename].r2_seq = str(each_read.seq)


##################################################### parse MappingRecord_dict ####################################################

for each_mp in MappingRecord_dict:
    current_mp_record = MappingRecord_dict[each_mp]
    print('%s\tcurrent_mp_r1_refs\t%s' % (each_mp, current_mp_record.r1_refs))
    print('%s\tcurrent_mp_r2_refs\t%s' % (each_mp, current_mp_record.r2_refs))

    # only r1 mapped
    if (current_mp_record.r1_refs != {}) and (current_mp_record.r2_refs == {}):

        #print('%s\tcurrent_mp_r1_refs\t%s' % (each_mp, current_mp_record.r1_refs))
        #print('%s\tcurrent_mp_r2_refs\t%s' % (each_mp, current_mp_r2_refs))

        r1_refs_filtered_by_M_len = {}
        for r1_ref in current_mp_record.r1_refs:
            r1_ref_cigar = current_mp_record.r1_refs[r1_ref]
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r1_ref_cigar))
            if aligned_len >= min_M_len_16s:
                r1_refs_filtered_by_M_len[r1_ref] = current_mp_record.r1_refs[r1_ref]
        #print('%s\tr1_refs_filtered\t%s' % (each_mp, r1_refs_filtered_by_M_len))
        #print()
        if len(r1_refs_filtered_by_M_len) == 0:
            current_mp_record.r1_refs = {}
        else:
            # consider unmapped mate
            current_mp_record.qualified_reads = True
            current_mp_record.consider_r1_unmapped_mate = True
            current_mp_record.r1_refs = r1_refs_filtered_by_M_len
            current_mp_record.r1_filtered_refs = {i.split('_pos_')[0] for i in r1_refs_filtered_by_M_len}

            # consider as clipping mapped
            r1_cigar_list = list(r1_refs_filtered_by_M_len.values())
            best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, [])
            best_cigar_flag = current_mp_record.r1_cigar_to_flag.get(best_cigar, '')
            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

            if max_value >= min_M_len_ctg:
                current_mp_record.consider_r1_clipping_part = True

                current_clipping_seq = ''
                if best_cigar_rc is False:
                    if max_value_index == 'r1_l':
                        current_clipping_seq = current_mp_record.r1_seq[:max_value]
                    if max_value_index == 'r1_r':
                        current_clipping_seq = current_mp_record.r1_seq[-max_value:]
                else:
                    r1_seq_rc = get_rc(current_mp_record.r1_seq)
                    if max_value_index == 'r1_l':
                        current_clipping_seq_rc = r1_seq_rc[:max_value]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)
                    if max_value_index == 'r1_r':
                        current_clipping_seq_rc = r1_seq_rc[-max_value:]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)
                current_mp_record.r1_clipping_seq = current_clipping_seq

    # only r2 mapped
    elif (current_mp_record.r1_refs == {}) and (current_mp_record.r2_refs != {}):

        r2_refs_filtered_by_M_len = {}
        for r2_ref in current_mp_record.r2_refs:
            r2_ref_cigar = current_mp_record.r2_refs[r2_ref]
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r2_ref_cigar))
            if aligned_len >= min_M_len_16s:
                r2_refs_filtered_by_M_len[r2_ref] = current_mp_record.r2_refs[r2_ref]

        if len(r2_refs_filtered_by_M_len) == 0:
            current_mp_record.r2_refs = {}
        else:
            # consider as paired
            current_mp_record.qualified_reads = True
            current_mp_record.consider_r2_unmapped_mate = True
            current_mp_record.r2_refs = r2_refs_filtered_by_M_len
            current_mp_record.r2_filtered_refs =  {i.split('_pos_')[0] for i in r2_refs_filtered_by_M_len}

            # consider as clipping mapped
            r2_cigar_list = list(r2_refs_filtered_by_M_len.values())
            best_cigar, max_value, max_value_index = get_max_clp_and_index([], r2_cigar_list)
            best_cigar_flag = current_mp_record.r2_cigar_to_flag.get(best_cigar, '')
            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

            if max_value >= min_M_len_ctg:
                current_mp_record.consider_r2_clipping_part = True

                current_clipping_seq = ''
                if best_cigar_rc is False:
                    if max_value_index == 'r2_l':
                        current_clipping_seq = current_mp_record.r2_seq[:max_value]
                    if max_value_index == 'r2_r':
                        current_clipping_seq = current_mp_record.r2_seq[-max_value:]
                else:
                    r2_seq_rc = get_rc(current_mp_record.r2_seq)
                    if max_value_index == 'r2_l':
                        current_clipping_seq_rc = r2_seq_rc[:max_value]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)
                    if max_value_index == 'r2_r':
                        current_clipping_seq_rc = r2_seq_rc[-max_value:]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                current_mp_record.r2_clipping_seq = current_clipping_seq

    # both of r1 and r2 mapped
    elif (current_mp_record.r1_refs != {}) and (current_mp_record.r2_refs != {}):

        r1_cigar_list = list(current_mp_record.r1_refs.values())
        r2_cigar_list = list(current_mp_record.r2_refs.values())

        print(r1_cigar_list)
        print(r2_cigar_list)

        r1_cigar_all_clp = True
        for each_cigar_r1 in r1_cigar_list:
            if ('S' not in each_cigar_r1) and ('s' not in each_cigar_r1):
                r1_cigar_all_clp = False

        r2_cigar_all_clp = True
        for each_cigar_r2 in r2_cigar_list:
            if ('S' not in each_cigar_r2) and ('s' not in each_cigar_r2):
                r2_cigar_all_clp = False

        best_cigar = ''
        max_value = 0
        max_value_index = ''
        if (r1_cigar_all_clp is True) and (r2_cigar_all_clp is False):
            best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, [])
        if (r1_cigar_all_clp is False) and (r2_cigar_all_clp is True):
            best_cigar, max_value, max_value_index = get_max_clp_and_index([], r2_cigar_list)
        if (r1_cigar_all_clp is True) and (r2_cigar_all_clp is True):
            best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, r2_cigar_list)

        if max_value >= min_M_len_ctg:

            best_cigar_flag = ''
            if max_value_index in ['r1_l', 'r1_r']:
                best_cigar_flag = current_mp_record.r1_cigar_to_flag.get(best_cigar, '')
            if max_value_index in ['r2_l', 'r2_r']:
                best_cigar_flag = current_mp_record.r2_cigar_to_flag.get(best_cigar, '')
            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

            current_mp_r1_refs_no_pos = {i.split('_pos_')[0] for i in current_mp_record.r1_refs}
            current_mp_r2_refs_no_pos = {i.split('_pos_')[0] for i in current_mp_record.r2_refs}
            current_mp_shared_refs_no_pos = current_mp_r1_refs_no_pos.intersection(current_mp_r2_refs_no_pos)

            if len(current_mp_shared_refs_no_pos) > 0:
                current_mp_record.qualified_reads = True

                # for clipping part, only consider shared refs
                current_clipping_seq = ''
                if best_cigar_rc is False:

                    if max_value_index == 'r1_l':
                        current_clipping_seq = current_mp_record.r1_seq[:max_value]

                    if max_value_index == 'r1_r':
                        current_clipping_seq = current_mp_record.r1_seq[-max_value:]

                    if max_value_index == 'r2_l':
                        current_clipping_seq = current_mp_record.r2_seq[:max_value]

                    if max_value_index == 'r2_r':
                        current_clipping_seq = current_mp_record.r2_seq[-max_value:]
                else:
                    r1_seq_rc = get_rc(current_mp_record.r1_seq)
                    r2_seq_rc = get_rc(current_mp_record.r2_seq)

                    if max_value_index == 'r1_l':
                        current_clipping_seq_rc = r1_seq_rc[:max_value]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if max_value_index == 'r1_r':
                        current_clipping_seq_rc = r1_seq_rc[-max_value:]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if max_value_index == 'r2_l':
                        current_clipping_seq_rc = r2_seq_rc[:max_value]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if max_value_index == 'r2_r':
                        current_clipping_seq_rc = r2_seq_rc[-max_value:]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                if max_value_index in ['r1_l', 'r1_r']:
                    current_mp_record.consider_r1_clipping_part = True
                    current_mp_record.r1_clipping_seq = current_clipping_seq
                    current_mp_record.r1_filtered_refs = current_mp_shared_refs_no_pos
                if max_value_index in ['r2_l', 'r2_r']:
                    current_mp_record.consider_r2_clipping_part = True
                    current_mp_record.r2_clipping_seq = current_clipping_seq
                    current_mp_record.r2_filtered_refs = current_mp_shared_refs_no_pos



# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    if MappingRecord_dict[each_mp].qualified_reads is False:
        MappingRecord_dict.pop(each_mp)

# write out sequences
unmapped_mates_handle = open(unmapped_mates_seq_file, 'w')
clipping_part_seq_handle = open(clipping_parts_seq_file, 'w')
for qualified_read in MappingRecord_dict:

    r1_name = '%s.1' % qualified_read
    r2_name = '%s.2' % qualified_read
    read_mr = MappingRecord_dict[qualified_read]

    if read_mr.consider_r1_unmapped_mate is True:
        unmapped_mates_handle.write('>%s\n' % r2_name)
        unmapped_mates_handle.write('%s\n' % read_mr.r2_seq)

    if read_mr.consider_r2_unmapped_mate is True:
        unmapped_mates_handle.write('>%s\n' % r1_name)
        unmapped_mates_handle.write('%s\n' % read_mr.r1_seq)

    if read_mr.consider_r1_clipping_part is True:
        clipping_part_seq_handle.write('>%s\n' % r1_name)
        clipping_part_seq_handle.write('%s\n' % read_mr.r1_clipping_seq)

    if read_mr.consider_r2_clipping_part is True:
        clipping_part_seq_handle.write('>%s\n' % r2_name)
        clipping_part_seq_handle.write('%s\n' % read_mr.r2_clipping_seq)

unmapped_mates_handle.close()
clipping_part_seq_handle.close()
