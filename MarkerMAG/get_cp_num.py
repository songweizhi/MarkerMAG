#!/usr/bin/env python3

# Copyright (C) 2020, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com or t.thomas@unsw.edu.au

# MarkerMAG is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MarkerMAG is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import glob
import shutil
import random
import argparse
import numpy as np
from Bio import SeqIO
from time import sleep
import multiprocessing as mp
from itertools import groupby
from datetime import datetime
from operator import itemgetter
from distutils.spawn import find_executable
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


get_cp_num_usage = '''
=======================================  Example commands =======================================

MarkerMAG get_cp_num -p prefix -r1 R1.fa -r2 R2.fa -linkages linkages_by_genome.txt -marker 16S.fa -mag mag_folder -x fa -t 12 -force

=================================================================================================
'''


def force_create_folder(folder_to_create):

    rm_rd = 0
    while os.path.isdir(folder_to_create) is True:
        shutil.rmtree(folder_to_create, ignore_errors=True)

        if rm_rd >= 10:
            print('Failed in removing %s, program exited!' % folder_to_create)
            exit()

        rm_rd += 1
        sleep(1)

    os.mkdir(folder_to_create)


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


def group_consecutive_numbers(num_list):

    group_list = []
    for k, g in groupby(enumerate(num_list), lambda ix: ix[0] - ix[1]):
        current_group = list(map(itemgetter(1), g))
        current_group_str = '%s-%s' % (current_group[0], current_group[-1])
        group_list.append(current_group_str)
    group_list_str = ','.join(group_list)

    return group_list_str


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def plot_gc_bias(GC_bias_txt, gc_content_mag, gc_content_16s_list, mag_global_mean_coverage, MS_GC_bias_plot):
    gc_content_list = []
    cov_list = []
    norm_cov_list = []
    count_num_list = []
    for each_line in open(GC_bias_txt):
        each_line_split = each_line.strip().split()
        if not each_line.startswith('GC'):
            gc = int(each_line_split[0])
            cov = float(each_line_split[1])
            count_num = int(each_line_split[2])
            gc_content_list.append(gc)
            cov_list.append(cov)
            count_num_list.append(count_num)

    gc_content_array = np.array(gc_content_list)
    cov_array = np.array(cov_list)
    count_num_array = np.array(count_num_list)

    fig = plt.figure(figsize=(12, 6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    # left one
    ax1.margins(0.05)
    ax1.plot(gc_content_array, cov_array)
    ax1.axvline(x=gc_content_mag, c='skyblue', linewidth=1, linestyle='dashed')
    for each_gc_content_16s in gc_content_16s_list:
        ax1.axvline(x=each_gc_content_16s, c='orange', linewidth=1, linestyle='dashed')
    ax1.axhline(y=mag_global_mean_coverage, c='black', linewidth=1, linestyle='dashed')
    ax1.set_xlabel('GC content (%)')
    ax1.set_ylabel('Coverage (X)')

    # right one
    ax2.margins(0.05)
    ax2.plot(gc_content_array, count_num_array)
    ax2.axvline(x=gc_content_mag, c='skyblue', linewidth=1, linestyle='dashed')
    for each_gc_content_16s in gc_content_16s_list:
        ax2.axvline(x=each_gc_content_16s, c='orange', linewidth=1, linestyle='dashed')
    ax2.axhline(y=0, c='black', linewidth=1, linestyle='dashed')
    ax2.set_xlabel('GC content (%)')
    ax2.set_ylabel('Number of Kmers (100 bp)')

    # save plot
    plt.tight_layout()
    plt.savefig(MS_GC_bias_plot)
    plt.close()


def get_scatter_plot(num_list_1, num_list_2, png_file):

    if (len(num_list_1) > 0) and (len(num_list_2) > 0):

        num_arrary_1 = np.array(num_list_1)
        num_arrary_2 = np.array(num_list_2)

        fig = plt.figure(figsize=(6, 6))
        plt.margins(0)

        plt.scatter(num_arrary_1, num_arrary_2)
        plt.xlabel("Estimated copy number", fontsize=15)
        plt.ylabel("User provided copy number", fontsize=15)

        # set axis range
        max_value = max([max(num_list_1), max(num_list_2)])
        plt.xlim(0, round(max_value + 1))
        plt.ylim(0, round(max_value + 1))

        # add fit line
        coeffs = np.polyfit(num_arrary_1, num_arrary_2, 1)
        slope = coeffs[0]
        intercept = coeffs[1]
        plt.plot(num_arrary_1, slope * num_arrary_1 + intercept)

        # get R-squared value
        p = np.poly1d(coeffs)
        yhat = p(num_arrary_1)
        ybar = np.sum(num_arrary_2)/len(num_arrary_2)
        ssreg = np.sum((yhat-ybar)**2)
        sstot = np.sum((num_arrary_2 - ybar)**2)
        r_squared = ssreg/sstot

        if intercept >= 0:
            title_str = 'y = %sx + %s, R2 = %s' % (float("{0:.4f}".format(slope)), float("{0:.4f}".format(intercept)), float("{0:.4f}".format(r_squared)))
        else:
            title_str = 'y = %sx - %s, R2 = %s' % (float("{0:.4f}".format(slope)), abs(float("{0:.4f}".format(intercept))), float("{0:.4f}".format(r_squared)))

        plt.title(title_str, fontsize=16)

        # save plot
        plt.tight_layout()
        plt.savefig(png_file, format='svg')
        plt.close()


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


class MappingRecord:

    def __init__(self):

        #################### overall ####################

        self.qualified_reads = False

        #################### round 1 16s ####################

        self.consider_r1_unmapped_mate = False
        self.consider_r2_unmapped_mate = False

        self.r1_16s_ref_dict = dict()
        self.r2_16s_ref_dict = dict()

        self.r1_16s_refs_lowest_mismatch = None
        self.r2_16s_refs_lowest_mismatch = None

        self.r1_16s_refs_no_ignored = dict()
        self.r2_16s_refs_no_ignored = dict()
        self.shared_16s_refs_no_ignored = dict()

        self.both_mapped_to_16s = False

        #################### round 1 ctg ####################

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()

        self.r1_ctg_refs_lowest_mismatch = None
        self.r2_ctg_refs_lowest_mismatch = None

        self.r1_ctg_refs_no_ignored = dict()
        self.r2_ctg_refs_no_ignored = dict()
        self.shared_ctg_refs_no_ignored = dict()

        self.matched_to_ctg = False

        #################### round 2 ####################

        self.qualified_reads_rd2 = False

        self.r1_ctg_ref_dict_rd2 = dict()
        self.r2_ctg_ref_dict_rd2 = dict()

        self.r1_ctg_refs_lowest_mismatch_rd2 = None
        self.r2_ctg_refs_lowest_mismatch_rd2 = None

        #################### round 2 mini_assembly ####################

        self.r1_mini_ref_dict = dict()
        self.r2_mini_ref_dict = dict()

        self.r1_mini_refs_lowest_mismatch = None
        self.r2_mini_refs_lowest_mismatch = None

        self.r1_mini_refs_no_ignored = dict()
        self.r2_mini_refs_no_ignored = dict()
        self.shared_mini_refs_no_ignored = dict()


def check_cigar_quality(cigar_str, consider_ref_min_mismatch, ref_min_mismatch, mismatch_cutoff, min_M_len, ref_pos, marker_len):

    # check the following:
    # 1. both end clip  2. mismatch     3. aligned length   4. clp in the middle    5. ref_min_mismatch

    cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        aln_len, aln_pct, clp_len, clp_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

        # check mismatch
        passed_mismatch_check = False
        if consider_ref_min_mismatch is False:
            if mismatch_pct <= mismatch_cutoff:
                passed_mismatch_check = True
        else:
            if ref_min_mismatch != 'NA':
                if (mismatch_pct <= ref_min_mismatch) and (mismatch_pct <= mismatch_cutoff):
                    passed_mismatch_check = True

        if passed_mismatch_check is True:
            # check aligned length
            if aln_len >= min_M_len:
                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar_str) or ('s' in cigar_str):
                    clip_in_middle = True
                    if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + aln_len - 1) == marker_len:
                            clip_in_middle = False

                if clip_in_middle is False:
                    qualified_cigar = True

    return qualified_cigar


def get_cigar_stats(cigar_splitted):

    # aligned_len: M X =
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
        if each_part_cate in {'M', 'm', 'X', 'x', '='}:
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


def get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len):

    mismatch_set_all_cigar = set()
    mismatch_set_long_M_cigars = set()
    for each_cigar in r1_ref_cigar_set:
        aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
        mismatch_set_all_cigar.add(mismatch_pct)
        if aligned_len >= min_M_len:
            mismatch_set_long_M_cigars.add(mismatch_pct)

    min_mismatch = 'NA'
    if len(mismatch_set_all_cigar) > 0:
        min_mismatch = min(mismatch_set_all_cigar)
        if len(mismatch_set_long_M_cigars) > 0:
            min_mismatch = min(mismatch_set_long_M_cigars)

    return min_mismatch


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


def r12_16s_ref_dict_to_str(r2_16s_ref_dict):

    top_dict_str = 'na'

    if r2_16s_ref_dict != {}:
        top_dict_list = []
        for sub_dict_key in r2_16s_ref_dict:
            sub_dict_value = r2_16s_ref_dict[sub_dict_key]
            bottom_dict_list = []
            for bottom_dict_key in sub_dict_value:
                bottom_dict_value = sub_dict_value[bottom_dict_key]
                bottom_dict_list.append('%s:%s' % (bottom_dict_key, bottom_dict_value))
            bottom_dict_str = '%s' % (';'.join(bottom_dict_list))
            top_dict_list.append('%s:::%s' % (sub_dict_key, bottom_dict_str))

        top_dict_str = ';;;'.join(top_dict_list)

    return top_dict_str


def get_r12_16s_ref_dict_from_str(r2_16s_ref_dict_str):

    ref_dict_from_str = {}
    if r2_16s_ref_dict_str != 'na':
        r2_16s_ref_dict_str_split = r2_16s_ref_dict_str.split(';;;')

        ref_dict_from_str = {}
        for each_sub_dict_str in r2_16s_ref_dict_str_split:
            each_sub_dict_str_split = each_sub_dict_str.split(':::')
            sub_dict_key = each_sub_dict_str_split[0]
            sub_dict_str = each_sub_dict_str_split[1]
            sub_dict_str_split = sub_dict_str.split(';')
            bottom_dict = dict()
            for each_bottom_dict_str in sub_dict_str_split:
                bottom_dict_str_split = each_bottom_dict_str.split(':')
                bottom_dict_key = int(bottom_dict_str_split[0])
                bottom_dict_value = bottom_dict_str_split[1]
                bottom_dict[bottom_dict_key] = bottom_dict_value
            ref_dict_from_str[sub_dict_key] = bottom_dict

    return ref_dict_from_str


def no_ignored_dict_to_str(r2_16s_refs_no_ignored_dict):
    top_dict_str = 'na'

    if r2_16s_refs_no_ignored_dict != {}:
        top_dict_list = []
        for each_ref in r2_16s_refs_no_ignored_dict:
            each_ref_cigar_list = r2_16s_refs_no_ignored_dict[each_ref]
            each_ref_cigar_str = ','.join(each_ref_cigar_list)
            top_dict_list.append('%s:%s' % (each_ref, each_ref_cigar_str))
        top_dict_str = ';'.join(top_dict_list)

    return top_dict_str


def get_no_ignored_dict_from_str(r2_16s_refs_no_ignored_dict_str):

    no_ignored_dict_from_str = {}
    if r2_16s_refs_no_ignored_dict_str != 'na':
        r2_16s_refs_no_ignored_dict_str_split = r2_16s_refs_no_ignored_dict_str.split(';')
        for each_sub_dict_str in r2_16s_refs_no_ignored_dict_str_split:
            each_sub_dict_str_split = each_sub_dict_str.split(':')
            ref_id = each_sub_dict_str_split[0]
            cigar_list_str = each_sub_dict_str_split[1]
            cigar_list = cigar_list_str.split(',')
            no_ignored_dict_from_str[ref_id] = cigar_list

    return no_ignored_dict_from_str


def get_mag_gc_bias_worker(argument_list):

    mag                                 = argument_list[0]
    current_mag_ctg_set                 = argument_list[1]
    current_mag_pos_depth_dict          = argument_list[2]
    ctg_seq_dict                        = argument_list[3]
    ctg_len_dict                        = argument_list[4]
    ctg_len_cutoff                      = argument_list[5]
    current_mag_region_to_ignore_dict   = argument_list[6]
    ignore_ends_len_mag                 = argument_list[7]
    window_len                          = argument_list[8]
    mag_gc_content_dict                 = argument_list[9]
    gnm_to_linked_16s_gc_content_dict   = argument_list[10]
    current_mag_gc_bias_txt             = argument_list[11]
    current_mag_gc_bias_png             = argument_list[12]
    current_mag_depth_GC_content_file   = argument_list[13]

    current_mag_total_mapped_read_len = 0
    current_mag_total_counted_ctg_len = 0
    current_mag_gc_content_to_depth_list_dict = {}
    for each_ctg in current_mag_ctg_set:
        current_ctg_pos_depth_dict = current_mag_pos_depth_dict.get(each_ctg, {})
        current_ctg_len = ctg_len_dict.get(each_ctg, 0)
        if (len(current_ctg_pos_depth_dict) > 0) and (current_ctg_len >= ctg_len_cutoff):
            current_ctg_seq = ctg_seq_dict[each_ctg]
            current_ctg_regions_to_ignore = current_mag_region_to_ignore_dict.get(each_ctg, set())
            window_start_pos = ignore_ends_len_mag + 1
            while window_start_pos <= (current_ctg_len - (ignore_ends_len_mag + window_len) + 1):
                window_end_pos = window_start_pos + window_len - 1

                if (str(window_start_pos) not in current_ctg_regions_to_ignore) and (str(window_end_pos) not in current_ctg_regions_to_ignore):
                    if (str(window_start_pos) in current_ctg_pos_depth_dict) and (str(window_end_pos) in current_ctg_pos_depth_dict):

                        # get sequence of current window
                        if window_start_pos < (len(current_ctg_seq) - window_len + 1):
                            window_seq = current_ctg_seq[(window_start_pos - 1):window_end_pos]
                        else:
                            window_seq = current_ctg_seq[(window_start_pos - 1):]
                        window_seq_upper = window_seq.upper()

                        masked_base_num = window_seq_upper.count('N')
                        if masked_base_num <= 4:
                            window_pos_list = list(range(window_start_pos, (window_end_pos + 1)))
                            window_pos_depth_list = [current_ctg_pos_depth_dict[str(pos)] for pos in window_pos_list]
                            window_mean_depth = sum(window_pos_depth_list) / len(window_pos_depth_list)
                            window_gc_content = (window_seq_upper.count('G') + window_seq_upper.count('C')) * 100 / len(window_seq)

                            current_mag_total_mapped_read_len += window_mean_depth
                            current_mag_total_counted_ctg_len += 1

                            # when the length of sliding window is 100 bp
                            window_gc_content = int(window_gc_content)

                            if window_gc_content not in current_mag_gc_content_to_depth_list_dict:
                                current_mag_gc_content_to_depth_list_dict[window_gc_content] = [window_mean_depth]
                            else:
                                current_mag_gc_content_to_depth_list_dict[window_gc_content].append(window_mean_depth)

                window_start_pos += 1

    # get current_mag_average_global_coverage
    if current_mag_total_counted_ctg_len == 0:
        current_mag_average_global_coverage = 0
    else:
        current_mag_average_global_coverage = current_mag_total_mapped_read_len / current_mag_total_counted_ctg_len
    current_mag_GG_content = mag_gc_content_dict.get(mag, 'NA')

    # write out current_mag_depth_GC_content_file
    current_mag_depth_GC_content_file_handle = open(current_mag_depth_GC_content_file, 'w')
    current_mag_depth_GC_content_file_handle.write('%s\t%s\t%s\n' % (mag, float("{0:.2f}".format(current_mag_average_global_coverage)), current_mag_GG_content))
    current_mag_depth_GC_content_file_handle.close()

    # write out current_mag_gc_bias_txt
    current_mag_gc_bias_txt_handle = open(current_mag_gc_bias_txt, 'w')
    current_mag_gc_bias_txt_handle.write('GC\tCoverage\tNumber\tProportion\n')
    for each_gc_group in sorted([i for i in current_mag_gc_content_to_depth_list_dict]):
        current_gc_group_depth_list = current_mag_gc_content_to_depth_list_dict[each_gc_group]
        current_gc_group_proportion = len(current_gc_group_depth_list) * 100 / current_mag_total_counted_ctg_len
        current_gc_group_mean_depth = sum(current_gc_group_depth_list) / len(current_gc_group_depth_list)

        current_mag_gc_bias_txt_handle.write('%s\t%s\t%s\t%s\n' % (each_gc_group,
                                                                   float("{0:.2f}".format(current_gc_group_mean_depth)),
                                                                   len(current_gc_group_depth_list),
                                                                   float("{0:.4f}".format(current_gc_group_proportion))))
    current_mag_gc_bias_txt_handle.close()

    # plot GC bias
    linked_16s_gc_content_list = gnm_to_linked_16s_gc_content_dict[mag]
    plot_gc_bias(current_mag_gc_bias_txt, current_mag_GG_content, linked_16s_gc_content_list,
                 current_mag_average_global_coverage, current_mag_gc_bias_png)


def parse_sam16s_worker(argument_list):

    sorted_sam              = argument_list[0]
    MappingRecord_file      = argument_list[1]
    min_M_len_16s           = argument_list[2]
    mismatch_cutoff         = argument_list[3]
    marker_len_dict         = argument_list[4]
    unlinked_refs_to_ignore = argument_list[5]

    MappingRecord_file_handle = open(MappingRecord_file, 'w')
    current_read_base = ''
    current_read_base_r1_16s_ref_dict = dict()
    current_read_base_r2_16s_ref_dict = dict()
    with open(sorted_sam) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if not each_read.startswith('@'):
                each_read_split = each_read.strip().split('\t')
                ref_id = each_read_split[2]
                if ref_id not in unlinked_refs_to_ignore:
                    cigar = each_read_split[5]
                    read_id = each_read_split[0]
                    read_id_base = '.'.join(read_id.split('.')[:-1])
                    read_strand = read_id.split('.')[-1]
                    ref_pos = int(each_read_split[3])

                    if current_read_base == '':
                        current_read_base = read_id_base

                        if cigar != '*':
                            if read_strand == '1':
                                current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                            if read_strand == '2':
                                current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}

                    elif read_id_base == current_read_base:

                        if cigar != '*':
                            if read_strand == '1':
                                if ref_id not in current_read_base_r1_16s_ref_dict:
                                    current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                                else:
                                    current_read_base_r1_16s_ref_dict[ref_id][ref_pos] = cigar
                            if read_strand == '2':
                                if ref_id not in current_read_base_r2_16s_ref_dict:
                                    current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}
                                else:
                                    current_read_base_r2_16s_ref_dict[ref_id][ref_pos] = cigar
                    else:
                        ################################### analysis previous read refs ####################################

                        current_read_base___qualified_reads             = 0
                        current_read_base___consider_r1_unmapped_mate   = 0
                        current_read_base___consider_r2_unmapped_mate   = 0
                        current_read_base___both_mapped_to_16s          = 0
                        current_read_base___r1_16s_ref_dict             = dict()
                        current_read_base___r2_16s_ref_dict             = dict()
                        current_read_base___r1_16s_refs_no_ignored      = dict()
                        current_read_base___r2_16s_refs_no_ignored      = dict()
                        current_read_base___shared_16s_refs_no_ignored  = dict()

                        ########## get lowest mismatch for r1/r2 16s refs ##########

                        # get r1_ref_cigar_set
                        r1_ref_cigar_set = set()
                        for each_pos_dict in current_read_base_r1_16s_ref_dict.values():
                            each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                            r1_ref_cigar_set.update(each_pos_dict_values)

                        # get r2_ref_cigar_set
                        r2_ref_cigar_set = set()
                        for each_pos_dict in current_read_base_r2_16s_ref_dict.values():
                            each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                            r2_ref_cigar_set.update(each_pos_dict_values)

                        r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_16s)
                        r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_16s)
                        refs_to_ignore = set()

                        ########## filter r1 16s refs ##########

                        r1_16s_refs_passed_qc = {}
                        r1_16s_refs_passed_qc_with_pos = {}
                        for r1_16s_ref in current_read_base_r1_16s_ref_dict:
                            r1_matched_pos_dict = current_read_base_r1_16s_ref_dict[r1_16s_ref]

                            # one read need to mapped to one 16S only for one time
                            if len(r1_matched_pos_dict) > 1:
                                refs_to_ignore.add(r1_16s_ref)
                            else:
                                r1_16s_ref_pos = list(r1_matched_pos_dict.keys())[0]
                                r1_16s_ref_cigar = r1_matched_pos_dict[r1_16s_ref_pos]

                                r1_16s_ref_qualified_cigar = check_cigar_quality(r1_16s_ref_cigar, True, r1_ref_min_mismatch,
                                                                                 mismatch_cutoff, min_M_len_16s, r1_16s_ref_pos,
                                                                                 marker_len_dict[r1_16s_ref])

                                # exclude the ref if cigar not qualified
                                if r1_16s_ref_qualified_cigar is False:
                                    refs_to_ignore.add(r1_16s_ref)
                                else:
                                    r1_16s_refs_passed_qc[r1_16s_ref] = [r1_16s_ref_cigar]
                                    r1_16s_refs_passed_qc_with_pos[r1_16s_ref] = {r1_16s_ref_pos: r1_16s_ref_cigar}

                        ########## filter r2 16s refs ##########

                        r2_16s_refs_passed_qc = {}
                        r2_16s_refs_passed_qc_with_pos = {}
                        for r2_16s_ref in current_read_base_r2_16s_ref_dict:
                            r2_matched_pos_dict = current_read_base_r2_16s_ref_dict[r2_16s_ref]

                            # one read need to mapped to one 16S only once
                            if len(r2_matched_pos_dict) > 1:
                                refs_to_ignore.add(r2_16s_ref)
                            else:
                                r2_16s_ref_pos = list(r2_matched_pos_dict.keys())[0]
                                r2_16s_ref_cigar = r2_matched_pos_dict[r2_16s_ref_pos]

                                r2_16s_ref_qualified_cigar = check_cigar_quality(r2_16s_ref_cigar, True, r2_ref_min_mismatch,
                                                                                 mismatch_cutoff, min_M_len_16s, r2_16s_ref_pos,
                                                                                 marker_len_dict[r2_16s_ref])

                                # exclude the ref if cigar not qualified
                                if r2_16s_ref_qualified_cigar is False:
                                    refs_to_ignore.add(r2_16s_ref)
                                else:
                                    r2_16s_refs_passed_qc[r2_16s_ref] = [r2_16s_ref_cigar]
                                    r2_16s_refs_passed_qc_with_pos[r2_16s_ref] = {r2_16s_ref_pos: r2_16s_ref_cigar}

                        #################################

                        r1_16s_refs_no_ignored          = {key: value for key, value in r1_16s_refs_passed_qc.items() if key not in refs_to_ignore}
                        r2_16s_refs_no_ignored          = {key: value for key, value in r2_16s_refs_passed_qc.items() if key not in refs_to_ignore}
                        r1_16s_refs_no_ignored_with_pos = {key: value for key, value in r1_16s_refs_passed_qc_with_pos.items() if key not in refs_to_ignore}
                        r2_16s_refs_no_ignored_with_pos = {key: value for key, value in r2_16s_refs_passed_qc_with_pos.items() if key not in refs_to_ignore}

                        # no mate has no_ignored alignments
                        if (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) == 0):
                            pass

                        # only r1 has no_ignored alignments
                        elif (len(r1_16s_refs_no_ignored) > 0) and (len(r2_16s_refs_no_ignored) == 0):

                            current_read_base___qualified_reads = 1
                            current_read_base___consider_r1_unmapped_mate = 1
                            current_read_base___r1_16s_refs_no_ignored = r1_16s_refs_no_ignored
                            current_read_base___r1_16s_ref_dict = current_read_base_r1_16s_ref_dict

                        # only r2 has no_ignored alignments
                        elif (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) > 0):

                            current_read_base___qualified_reads = 1
                            current_read_base___consider_r2_unmapped_mate = 1
                            current_read_base___r2_16s_refs_no_ignored = r2_16s_refs_no_ignored
                            current_read_base___r2_16s_ref_dict = current_read_base_r2_16s_ref_dict

                        # both r1 and r2 have no_ignored alignments
                        else:
                            shared_16s_ref_dict = {key: [r1_16s_refs_no_ignored[key][0], r2_16s_refs_no_ignored[key][0]] for key
                                                   in set(r1_16s_refs_no_ignored).intersection(set(r2_16s_refs_no_ignored))}
                            if len(shared_16s_ref_dict) > 0:

                                current_read_base___qualified_reads = 1
                                current_read_base___both_mapped_to_16s = 1
                                current_read_base___shared_16s_refs_no_ignored = shared_16s_ref_dict
                                current_read_base___r1_16s_ref_dict = current_read_base_r1_16s_ref_dict
                                current_read_base___r2_16s_ref_dict = current_read_base_r2_16s_ref_dict

                        if current_read_base___qualified_reads == 1:
                            current_read_base___r1_16s_ref_dict_str            = r12_16s_ref_dict_to_str(current_read_base___r1_16s_ref_dict)
                            current_read_base___r2_16s_ref_dict_str            = r12_16s_ref_dict_to_str(current_read_base___r2_16s_ref_dict)
                            current_read_base___r1_16s_refs_no_ignored_str     = no_ignored_dict_to_str(current_read_base___r1_16s_refs_no_ignored)
                            current_read_base___r2_16s_refs_no_ignored_str     = no_ignored_dict_to_str(current_read_base___r2_16s_refs_no_ignored)
                            current_read_base___shared_16s_refs_no_ignored_str = no_ignored_dict_to_str(current_read_base___shared_16s_refs_no_ignored)
                            MappingRecord_file_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (current_read_base,
                                                                                                          current_read_base___qualified_reads,
                                                                                                          current_read_base___consider_r1_unmapped_mate,
                                                                                                          current_read_base___consider_r2_unmapped_mate,
                                                                                                          current_read_base___both_mapped_to_16s,
                                                                                                          current_read_base___r1_16s_ref_dict_str,
                                                                                                          current_read_base___r2_16s_ref_dict_str,
                                                                                                          current_read_base___r1_16s_refs_no_ignored_str,
                                                                                                          current_read_base___r2_16s_refs_no_ignored_str,
                                                                                                          current_read_base___shared_16s_refs_no_ignored_str))


                        ########################################### reset values ###########################################

                        current_read_base = read_id_base
                        current_read_base_r1_16s_ref_dict = dict()
                        current_read_base_r2_16s_ref_dict = dict()

                        if cigar != '*':
                            if read_strand == '1':
                                current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                            if read_strand == '2':
                                current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}

    MappingRecord_file_handle.close()


def remove_high_mismatch(sam_in, aln_len_cutoff, mismatch_cutoff, sam_out):

    sam_out_handle = open(sam_out, 'w')
    ref_len_dict = {}
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)

            marker_id = ''
            marker_len = 0
            for each_element in each_read_split:
                if each_element.startswith('SN:'):
                    marker_id = each_element[3:]
                if each_element.startswith('LN:'):
                    marker_len = int(each_element[3:])
            ref_len_dict[marker_id] = marker_len

        else:
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                cigar_splitted = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitted)
                if both_ends_clp is False:
                    r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)
                    if r1_mismatch_pct <= mismatch_cutoff:

                        if r1_aligned_len >= aln_len_cutoff:

                            # check if clp in middle
                            if ('S' not in cigar) and ('s' not in cigar):
                                sam_out_handle.write(each_read)
                            else:
                                clip_in_middle = True
                                if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                                    clip_in_middle = False
                                if (cigar_splitted[-1][-1] in ['S', 's']):
                                    if (ref_pos + r1_aligned_len - 1) == ref_len_dict[ref_id]:
                                        clip_in_middle = False

                                if clip_in_middle is False:
                                    sam_out_handle.write(each_read)
    sam_out_handle.close()


def parse_sam_file_with_multiprocessing(sam_sorted_by_read, sam_sorted_line_num_txt, sam_sorted_by_read_split, MappingRecord_folder, min_M_len_16s, mismatch_cutoff, marker_len_dict, num_threads):

    # get the number of lines per file
    os.system('wc -l %s > %s' % (sam_sorted_by_read, sam_sorted_line_num_txt))
    sam16s_line_num = int(open(sam_sorted_line_num_txt).readline().strip().split(' ')[0])
    os.remove(sam_sorted_line_num_txt)
    line_num_per_file = int(round(sam16s_line_num / (num_threads * 10))) + 10

    # split sam file
    os.mkdir(sam_sorted_by_read_split)
    split_sam_cmd = 'split -l %s %s %s/splitted_sam_ ' % (line_num_per_file, sam_sorted_by_read, sam_sorted_by_read_split)
    os.system(split_sam_cmd)

    # get splitted sam file list
    splitted_sam_file_list = [os.path.basename(file_name) for file_name in glob.glob('%s/splitted_sam_*' % sam_sorted_by_read_split)]

    # prepare lol for mp worker
    list_for_parse_sam16s_worker = []
    splitted_sam_mp_file_set = set()
    for splitted_sam_file in splitted_sam_file_list:
        pwd_splitted_sam_file = '%s/%s' % (sam_sorted_by_read_split, splitted_sam_file)
        pwd_splitted_sam_mp_file = '%s/%s_mp.txt' % (MappingRecord_folder, splitted_sam_file)
        splitted_sam_mp_file_set.add(pwd_splitted_sam_mp_file)
        list_for_parse_sam16s_worker.append([pwd_splitted_sam_file,
                                             pwd_splitted_sam_mp_file,
                                             min_M_len_16s,
                                             mismatch_cutoff,
                                             marker_len_dict])

    # parsing mappping results with multiprocessing
    pool_parse_sam16s = mp.Pool(processes=num_threads)
    pool_parse_sam16s.map(parse_sam16s_worker, list_for_parse_sam16s_worker)
    pool_parse_sam16s.close()
    pool_parse_sam16s.join()

    #report_and_log(('Round 1: removing splitted subsets from disk'), pwd_log_file, keep_quiet)
    os.system('rm -r %s' % sam_sorted_by_read_split)

    return splitted_sam_mp_file_set


def get_cigar_aln_len(cigar_splitted):
    # aligned_len: M I X =
    aligned_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]
        # get aligned_len
        if each_part_cate in {'M', 'm', 'I', 'i', 'X', 'x', '='}:
            aligned_len += each_part_len
    return aligned_len


def sep_paired_and_singleton_reads(seq_in, seq_out_r1, seq_out_r2, seq_out_singleton):

    input_seq_ext = os.path.splitext(seq_in)[1]
    input_seq_fmt = 'fasta'
    if ('q' in input_seq_ext) or ('Q' in input_seq_ext):
        input_seq_fmt = 'fastq'

    fmt_to_provide_for_write = 'fasta-2line'
    if input_seq_fmt == 'fastq':
        fmt_to_provide_for_write = 'fastq'

    reads_pair_dict = {}
    for read_record in SeqIO.parse(seq_in, input_seq_fmt):
        read_id_base = '.'.join(read_record.id.split('.')[:-1])
        read_strand = read_record.id.split('.')[-1]
        if read_id_base not in reads_pair_dict:
            reads_pair_dict[read_id_base] = {read_strand}
        else:
            reads_pair_dict[read_id_base].add(read_strand)

    read_list_paired = set()
    read_list_singleton = set()
    for read_base in reads_pair_dict:
        if len(reads_pair_dict[read_base]) == 1:
            read_list_singleton.add(read_base)
        if len(reads_pair_dict[read_base]) == 2:
            read_list_paired.add(read_base)

    fasta_out_r1_handle = open(seq_out_r1, 'w')
    fasta_out_r2_handle = open(seq_out_r2, 'w')
    fasta_out_singleton_handle = open(seq_out_singleton, 'w')
    for read_record in SeqIO.parse(seq_in, input_seq_fmt):
        read_id_base = '.'.join(read_record.id.split('.')[:-1])
        read_strand = read_record.id.split('.')[-1]
        if read_id_base in read_list_singleton:
            SeqIO.write(read_record, fasta_out_singleton_handle, fmt_to_provide_for_write)
        if read_id_base in read_list_paired:
            if read_strand == '1':
                SeqIO.write(read_record, fasta_out_r1_handle, fmt_to_provide_for_write)
            if read_strand == '2':
                SeqIO.write(read_record, fasta_out_r2_handle, fmt_to_provide_for_write)
    fasta_out_r1_handle.close()
    fasta_out_r2_handle.close()
    fasta_out_singleton_handle.close()


def keep_only_both_mapped_reads(sam_in, sam_out):
    mapping_dic = {}
    with open(sam_in) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if not each_read.startswith('@'):
                each_read_split = each_read.strip().split('\t')
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]

                if read_id_base not in mapping_dic:
                    mapping_dic[read_id_base] = {read_strand}
                else:
                    mapping_dic[read_id_base].add(read_strand)

    mapping_dic_both_mapped = {}
    for read_base in mapping_dic:
        if len(mapping_dic[read_base]) == 2:
            mapping_dic_both_mapped[read_base] = mapping_dic[read_base]

    sam_out_handle = open(sam_out, 'w')
    with open(sam_in) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if each_read.startswith('@'):
                sam_out_handle.write(each_read)
            else:
                each_read_split = each_read.strip().split('\t')
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                if read_id_base in mapping_dic_both_mapped:
                    sam_out_handle.write(each_read)
    sam_out_handle.close()


def get_well_aligned_reads_pct_dict(mapped_to_linked_16s_sam, mapped_to_all_16s_sam, consider_as_mis_assignment_pct):

    marker_to_check = ''

    marker_id_list = []
    ref_to_read_dict_linked = {}
    with open(mapped_to_linked_16s_sam) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            each_read_split = each_read.strip().split('\t')

            if each_read.startswith('@'):
                mini_assembly_id = ''
                for each_element in each_read_split:
                    if each_element.startswith('SN:'):
                        mini_assembly_id = each_element[3:]
                if mini_assembly_id != '':
                    marker_id_list.append(mini_assembly_id)
            else:
                read_id = each_read_split[0]
                ref_id = each_read_split[2]
                if ref_id not in ref_to_read_dict_linked:
                    ref_to_read_dict_linked[ref_id] = {read_id}
                else:
                    ref_to_read_dict_linked[ref_id].add(read_id)

    read_to_ref_dict_all = {}
    ref_to_read_dict_all = {}
    with open(mapped_to_all_16s_sam) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if not each_read.startswith('@'):
                each_read_split = each_read.strip().split('\t')
                read_id = each_read_split[0]
                ref_id = each_read_split[2]
                if read_id not in read_to_ref_dict_all:
                    read_to_ref_dict_all[read_id] = {ref_id}
                else:
                    read_to_ref_dict_all[read_id].add(ref_id)

                if ref_id not in ref_to_read_dict_all:
                    ref_to_read_dict_all[ref_id] = {read_id}
                else:
                    ref_to_read_dict_all[ref_id].add(read_id)

    if marker_to_check != '':
        print('Total reads mapped to %s in %s: %s' % (marker_to_check, mapped_to_linked_16s_sam, len(ref_to_read_dict_linked[marker_to_check])))
        print('Total reads mapped to %s in %s: %s' % (marker_to_check, mapped_to_all_16s_sam, len(ref_to_read_dict_all[marker_to_check])))
        print('The %s reads were aligned to the following references in %s:' % (len(ref_to_read_dict_linked[marker_to_check]), mapped_to_all_16s_sam))

    well_assigned_reads_pct_dict = {}
    for marker_id in marker_id_list:

        count_dict = {}
        for each_read in ref_to_read_dict_linked[marker_id]:
            read_ref_in_all_map = read_to_ref_dict_all.get(each_read, set())
            for each_ref in read_ref_in_all_map:
                if each_ref not in count_dict:
                    count_dict[each_ref] = 1
                else:
                    count_dict[each_ref] += 1

        current_marker_assigned_reads_when_map_to_all = count_dict[marker_id]

        read_num_well_assigned = 0
        read_num_total_assigned = 0
        for each_ref in count_dict:
            current_ref_num = count_dict[each_ref]
            current_ref_pct = current_ref_num * 100 / current_marker_assigned_reads_when_map_to_all
            current_ref_pct = float("{0:.2f}".format(current_ref_pct))
            if current_ref_pct < consider_as_mis_assignment_pct:
                read_num_total_assigned += current_ref_num
                #print('%s\t%s\t%s\tbad' % (each_ref, current_ref_num, current_ref_pct))
            else:
                read_num_total_assigned += current_ref_num
                read_num_well_assigned += current_ref_num
                #print('%s\t%s\t%s\tgood' % (each_ref, current_ref_num, current_ref_pct))

        well_aligned_reads_pct = 0
        if read_num_total_assigned != 0:
            well_aligned_reads_pct = read_num_well_assigned * 100 / read_num_total_assigned
            well_aligned_reads_pct = float("{0:.2f}".format(well_aligned_reads_pct))
        well_assigned_reads_pct_dict[marker_id] = well_aligned_reads_pct

    return well_assigned_reads_pct_dict


def filter_sam_file(splitted_sam_mp_file_set, min_insert_size_16s, both_pair_mapped, keep_uniq,
                    sam_file_sorted, sam_file_filtered, sam_file_filtered_random):

    # read filter results into dict
    large_insert_size_paired_read_num = 0
    small_insert_size_paired_read_num = 0
    MappingRecord_dict = {}
    for each_mp_file in splitted_sam_mp_file_set:
        first_line_base = True
        with open(each_mp_file) as each_mp_file_opened:
            for each_read_base in each_mp_file_opened:
                each_read_base_split = each_read_base.strip().split('\t')
                current_read_base___id = each_read_base_split[0]
                current_read_base___qualified_reads = each_read_base_split[1]
                current_read_base___consider_r1_unmapped_mate = each_read_base_split[2]
                current_read_base___consider_r2_unmapped_mate = each_read_base_split[3]
                current_read_base___both_mapped_to_16s = each_read_base_split[4]
                current_read_base___r1_16s_ref_dict = get_r12_16s_ref_dict_from_str(each_read_base_split[5])
                current_read_base___r2_16s_ref_dict = get_r12_16s_ref_dict_from_str(each_read_base_split[6])
                current_read_base___r1_16s_refs_no_ignored = get_no_ignored_dict_from_str(
                    each_read_base_split[7])
                current_read_base___r2_16s_refs_no_ignored = get_no_ignored_dict_from_str(
                    each_read_base_split[8])
                current_read_base___shared_16s_refs_no_ignored = get_no_ignored_dict_from_str(
                    each_read_base_split[9])

                if first_line_base is True:
                    first_line_base = False
                else:
                    MappingRecord_dict[current_read_base___id] = MappingRecord()
                    MappingRecord_dict[current_read_base___id].qualified_reads = True

                    if current_read_base___consider_r1_unmapped_mate == '1':
                        MappingRecord_dict[current_read_base___id].consider_r1_unmapped_mate = True
                    if current_read_base___consider_r2_unmapped_mate == '1':
                        MappingRecord_dict[current_read_base___id].consider_r2_unmapped_mate = True
                    if current_read_base___both_mapped_to_16s == '1':
                        MappingRecord_dict[current_read_base___id].both_mapped_to_16s = True

                    MappingRecord_dict[
                        current_read_base___id].r1_16s_ref_dict = current_read_base___r1_16s_ref_dict
                    MappingRecord_dict[
                        current_read_base___id].r2_16s_ref_dict = current_read_base___r2_16s_ref_dict
                    MappingRecord_dict[
                        current_read_base___id].r1_16s_refs_no_ignored = current_read_base___r1_16s_refs_no_ignored
                    MappingRecord_dict[
                        current_read_base___id].r2_16s_refs_no_ignored = current_read_base___r2_16s_refs_no_ignored

                    current_read_base___shared_16s_refs_no_ignored_with_large_insert_size = {}
                    if len(current_read_base___shared_16s_refs_no_ignored) > 0:
                        for each_shared_ref in current_read_base___shared_16s_refs_no_ignored:
                            r1_ref_pos = list(current_read_base___r1_16s_ref_dict[each_shared_ref].keys())[0]
                            r1_ref_cigar = current_read_base___r1_16s_ref_dict[each_shared_ref][r1_ref_pos]
                            r1_ref_cigar_aln_len = get_cigar_aln_len(cigar_splitter(r1_ref_cigar))
                            r2_ref_pos = list(current_read_base___r2_16s_ref_dict[each_shared_ref].keys())[0]
                            r2_ref_cigar = current_read_base___r2_16s_ref_dict[each_shared_ref][r2_ref_pos]
                            r2_ref_cigar_aln_len = get_cigar_aln_len(cigar_splitter(r2_ref_cigar))

                            left_cigar_right_end = 0
                            right_cigar_left_end = 0
                            if r1_ref_pos == min([r1_ref_pos, r2_ref_pos]):
                                left_cigar_right_end = r1_ref_pos + r1_ref_cigar_aln_len
                                right_cigar_left_end = r2_ref_pos
                            if r2_ref_pos == min([r1_ref_pos, r2_ref_pos]):
                                left_cigar_right_end = r2_ref_pos + r2_ref_cigar_aln_len
                                right_cigar_left_end = r1_ref_pos
                            current_insert_size = right_cigar_left_end - left_cigar_right_end

                            if current_insert_size >= min_insert_size_16s:
                                current_read_base___shared_16s_refs_no_ignored_with_large_insert_size[
                                    each_shared_ref] = \
                                    current_read_base___shared_16s_refs_no_ignored[each_shared_ref]
                                large_insert_size_paired_read_num += 1
                            else:
                                small_insert_size_paired_read_num += 1

                    MappingRecord_dict[
                        current_read_base___id].shared_16s_refs_no_ignored = current_read_base___shared_16s_refs_no_ignored_with_large_insert_size
    # print('large_insert_size_paired_read_num\t%s' % large_insert_size_paired_read_num)
    # print('small_insert_size_paired_read_num\t%s' % small_insert_size_paired_read_num)

    # filter sam file
    pwd_sam_file_filtered_handle = open(sam_file_filtered, 'w')
    for each_line in open(sam_file_sorted):
        if each_line.startswith('@'):
            pwd_sam_file_filtered_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            read_base = read_id[:-2]
            read_strand = read_id[-1]
            ref_id = each_line_split[2]
            if read_base in MappingRecord_dict:

                to_write = False
                if (ref_id in MappingRecord_dict[read_base].shared_16s_refs_no_ignored):
                    if keep_uniq is True:
                        if len(MappingRecord_dict[read_base].shared_16s_refs_no_ignored) == 1:
                            to_write = True
                    else:
                        to_write = True

                if both_pair_mapped is False:
                    if (ref_id not in MappingRecord_dict[read_base].shared_16s_refs_no_ignored):
                        if read_strand == '1':
                            if (ref_id in MappingRecord_dict[read_base].r1_16s_refs_no_ignored):
                                if keep_uniq is True:
                                    if len(MappingRecord_dict[read_base].r1_16s_refs_no_ignored) == 1:
                                        to_write = True
                                else:
                                    to_write = True
                        if read_strand == '2':
                            if (ref_id in MappingRecord_dict[read_base].r2_16s_refs_no_ignored):
                                if keep_uniq is True:
                                    if len(MappingRecord_dict[read_base].r2_16s_refs_no_ignored) == 1:
                                        to_write = True
                                else:
                                    to_write = True
                if to_write is True:
                    pwd_sam_file_filtered_handle.write(each_line)
    pwd_sam_file_filtered_handle.close()

    # read filtered sam into dict
    read_to_ref_dict = {}
    ref_len_dict = {}
    read_len_dict = {}
    for each_line in open(sam_file_filtered):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('@'):
            mini_assembly_id = ''
            mini_assembly_len = 0
            for each_element in each_line_split:
                if each_element.startswith('SN:'):
                    mini_assembly_id = each_element[3:]
                if each_element.startswith('LN:'):
                    mini_assembly_len = int(each_element[3:])
            ref_len_dict[mini_assembly_id] = mini_assembly_len
        else:
            read_id = each_line_split[0]
            ref_id = each_line_split[2]
            read_len = len(each_line_split[9])

            # if ref_id in linked_16s_set:
            if read_id not in read_to_ref_dict:
                read_to_ref_dict[read_id] = {ref_id}
            else:
                read_to_ref_dict[read_id].add(ref_id)
            read_len_dict[read_id] = read_len

    # randomly select one ref for reads mapped to multiple 16S rRNA gene sequences
    ref_to_read_dict_overall = {}
    ref_to_read_dict_random_one = {}
    read_to_ref_dict_random_one = {}
    for each_read in read_to_ref_dict:
        matched_ref_list = [i for i in read_to_ref_dict[each_read]]
        random_ref = random.sample(matched_ref_list, 1)[0]
        read_to_ref_dict_random_one[each_read] = random_ref
        if random_ref not in ref_to_read_dict_random_one:
            ref_to_read_dict_random_one[random_ref] = {each_read}
        else:
            ref_to_read_dict_random_one[random_ref].add(each_read)

        for each_matched_ref in matched_ref_list:
            if each_matched_ref not in ref_to_read_dict_overall:
                ref_to_read_dict_overall[each_matched_ref] = {each_read}
            else:
                ref_to_read_dict_overall[each_matched_ref].add(each_read)

    # get random one sam file
    pwd_sam_file_filtered_random_handle = open(sam_file_filtered_random, 'w')
    for each_line in open(sam_file_filtered):
        if each_line.startswith('@'):
            pwd_sam_file_filtered_random_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            ref_id = each_line_split[2]
            random_ref = read_to_ref_dict_random_one.get(read_id, 'NA')
            if ref_id == random_ref:
                pwd_sam_file_filtered_random_handle.write(each_line)
    pwd_sam_file_filtered_random_handle.close()


def get_depth_from_sam(sam_to_use, linked_16s_to_gnm_dict, marker_seq_dict, marker_gc_content_dict, ignore_ends_len_16s,
                       pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest,
                       ref_depth_by_bp_file, assessment_txt_pos, all_16s_depth_GC_content_file):

    marker_pos_depth_dict = {}
    with open(sam_to_use) as pwd_sam_file_filtered_random_opened:
        for each_read in pwd_sam_file_filtered_random_opened:
            each_read_split = each_read.strip().split('\t')

            # initialize marker_pos_depth_dict dict
            if each_read.startswith('@'):
                marker_id = ''
                marker_len = 0
                for each_element in each_read_split:
                    if each_element.startswith('SN:'):
                        marker_id = each_element[3:]
                    if each_element.startswith('LN:'):
                        marker_len = int(each_element[3:])
                pos_list = list(range(1, (marker_len + 1)))
                if marker_id != '':
                    if marker_id in linked_16s_to_gnm_dict:
                        current_marker_pos_dict = {}
                        for each_pos in pos_list:
                            current_marker_pos_dict[str(each_pos)] = 0
                        marker_pos_depth_dict[marker_id] = current_marker_pos_dict
            else:
                cigar = each_read_split[5]
                ref_id = each_read_split[2]
                if ref_id in linked_16s_to_gnm_dict:
                    ref_pos = int(each_read_split[3])
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(
                        cigar_splitter(cigar))
                    ref_pos_end = ref_pos + aligned_len
                    current_read_covered_region = list(range(ref_pos, ref_pos_end))
                    for each_covered_pos in current_read_covered_region:
                        covered_pos_str = str(each_covered_pos)
                        marker_pos_depth_dict[ref_id][covered_pos_str] += 1

    # write out marker_pos_depth_dict
    ref_depth_by_bp_file_handle = open(ref_depth_by_bp_file, 'w')
    for each_ref in marker_pos_depth_dict:
        ref_depth_by_bp_dict = marker_pos_depth_dict[each_ref]
        for each_bp in list(range(1, (len(ref_depth_by_bp_dict) + 1))):
            each_bp_depth = ref_depth_by_bp_dict[str(each_bp)]
            ref_depth_by_bp_file_handle.write('%s\t%s\t%s\n' % (each_ref, each_bp, each_bp_depth))
    ref_depth_by_bp_file_handle.close()

    # get outliers
    outlier_depth_dict_16s = {}
    for each_16s in marker_pos_depth_dict:
        len_16s = len(marker_pos_depth_dict[each_16s])
        current_16s_pos_depth_dict = marker_pos_depth_dict[each_16s]
        current_16s_pos_list = sorted(list([int(i) for i in current_16s_pos_depth_dict.keys()]))

        non_end_pos_depth_list = []
        for each_pos in current_16s_pos_list:
            if (each_pos > ignore_ends_len_16s) and (each_pos <= (len_16s - ignore_ends_len_16s)):
                current_pos_depth = marker_pos_depth_dict[each_16s][str(each_pos)]
                non_end_pos_depth_list.append(current_pos_depth)
        non_end_pos_depth_list_sorted = sorted(non_end_pos_depth_list)

        value_num_to_ignore_lowest = round(len(non_end_pos_depth_list) * (pos_pct_to_ignore_lowest / 100))
        value_num_to_ignore_highest = round(len(non_end_pos_depth_list) * (pos_pct_to_ignore_highest / 100))

        smallest_values = non_end_pos_depth_list_sorted[:value_num_to_ignore_lowest]
        largest_values = non_end_pos_depth_list_sorted[(len(non_end_pos_depth_list) - value_num_to_ignore_highest):]

        if pos_pct_to_ignore_lowest == 0:
            small_value_cutoff = 0
        else:
            small_value_cutoff = smallest_values[-1]

        if pos_pct_to_ignore_highest == 0:
            large_value_cutoff = max(non_end_pos_depth_list) + 1
        else:
            large_value_cutoff = largest_values[0]

        outlier_depth_dict_16s[each_16s] = [value_num_to_ignore_lowest, value_num_to_ignore_highest, small_value_cutoff, large_value_cutoff]

    # write 16S depth and GC content
    assessment_txt_pos_handle = open(assessment_txt_pos, 'w')
    assessment_txt_pos_handle.write('16S\tCounted_depth_range(X)\tIgnored_low_depth_bps\tIgnored_high_depth_bps\tCounted_bps/Total_length\n')
    depth_GC_content_file_16s_handle = open(all_16s_depth_GC_content_file, 'w')
    depth_GC_content_file_16s_handle.write('16S\tDepth(x)\tGC(%)\n')
    for each_linked_16s in linked_16s_to_gnm_dict:
        linked_16s_seq = marker_seq_dict[each_linked_16s]
        linked_16s_len = len(linked_16s_seq)
        linked_16s_pos_depth_dict = marker_pos_depth_dict[each_linked_16s]
        current_16s_total_mapped_read_len = 0
        current_16s_total_counted_ctg_len = 0
        currrent_16s_depth_num_to_ignore_lowest = outlier_depth_dict_16s[each_linked_16s][0]
        currrent_16s_depth_num_to_ignore_highest = outlier_depth_dict_16s[each_linked_16s][1]
        currrent_16s_low_depth_cutoff = outlier_depth_dict_16s[each_linked_16s][2]
        currrent_16s_high_depth_cutoff = outlier_depth_dict_16s[each_linked_16s][3]

        ignored_low_depth_pos_num = 0
        ignored_high_depth_pos_num = 0
        pos_index = ignore_ends_len_16s + 1
        current_counted_pos_list = []
        while pos_index <= (linked_16s_len - ignore_ends_len_16s):
            current_pos_depth = linked_16s_pos_depth_dict[str(pos_index)]

            if ((current_pos_depth <= currrent_16s_low_depth_cutoff) and (
                    ignored_low_depth_pos_num < currrent_16s_depth_num_to_ignore_lowest)) or (
                    (current_pos_depth >= currrent_16s_high_depth_cutoff) and (
                    ignored_high_depth_pos_num < currrent_16s_depth_num_to_ignore_highest)):
                if current_pos_depth <= currrent_16s_low_depth_cutoff:
                    ignored_low_depth_pos_num += 1
                if current_pos_depth >= currrent_16s_high_depth_cutoff:
                    ignored_high_depth_pos_num += 1
            else:
                current_counted_pos_list.append(pos_index)
                current_16s_total_mapped_read_len += current_pos_depth
                current_16s_total_counted_ctg_len += 1
            pos_index += 1

        current_counted_pos_str = group_consecutive_numbers(current_counted_pos_list)

        assessment_txt_pos_handle.write('%s\t%s-%s\t%s\t%s\t%s/%s\t%s\n' % (each_linked_16s,
                                                                            currrent_16s_low_depth_cutoff,
                                                                            currrent_16s_high_depth_cutoff,
                                                                            currrent_16s_depth_num_to_ignore_lowest,
                                                                            currrent_16s_depth_num_to_ignore_highest,
                                                                            len(current_counted_pos_list),
                                                                            linked_16s_len,
                                                                            current_counted_pos_str))

        current_16s_average_coverage = current_16s_total_mapped_read_len / current_16s_total_counted_ctg_len
        current_16s_average_coverage = float("{0:.2f}".format(current_16s_average_coverage))
        current_16s_gc_content = marker_gc_content_dict[each_linked_16s]
        depth_GC_content_file_16s_handle.write(
            '%s\t%s\t%s\n' % (each_linked_16s, current_16s_average_coverage, current_16s_gc_content))
    depth_GC_content_file_16s_handle.close()
    assessment_txt_pos_handle.close()


def calculate_cp_num(gnm_to_linked_16s_dict, mag_gc_bias_op_folder, depth_dict_mag, depth_dict_16s,
                     gc_content_dict_mag, gc_content_dict_16s, estimation_with_depth_txt):

    depth_txt_handle = open(estimation_with_depth_txt, 'w')
    depth_txt_handle.write('MAG\t16S\tMAG_depth\t16S_depth\t16S_depth_norm\t16S_cp_num\t16S_cp_num_norm\n')
    for each_mag in gnm_to_linked_16s_dict:
        current_mag_linked_16s = gnm_to_linked_16s_dict[each_mag]
        gc_bias_txt = '%s/%s_GC_bias.txt' % (mag_gc_bias_op_folder, each_mag)

        if os.path.isfile(gc_bias_txt):
            # read in gc bias
            gc_bias_dict = {}
            for each_gc in open(gc_bias_txt):
                if not each_gc.startswith('GC\tCoverage\tNumber\tProportion'):
                    each_gc_split = each_gc.strip().split('\t')
                    gc_value = int(each_gc_split[0])
                    gc_depth = float(each_gc_split[1])
                    gc_bias_dict[gc_value] = gc_depth

            for each_linked_16s in current_mag_linked_16s:
                depth_mag = depth_dict_mag.get(each_mag, 'NA')
                depth_16s = depth_dict_16s.get(each_linked_16s, 'NA')
                gc_content_mag = gc_content_dict_mag.get(each_mag, 'NA')
                gc_content_16s = gc_content_dict_16s.get(each_linked_16s, 'NA')
                gc_content_mag_rounded = round(float(gc_content_mag))
                gc_content_16s_rounded = round(float(gc_content_16s))
                mean_gc_depth_mag = gc_bias_dict[gc_content_mag_rounded]
                mean_gc_depth_16s = gc_bias_dict[gc_content_16s_rounded]
                depth_ratio = mean_gc_depth_16s / mean_gc_depth_mag
                depth_16s_normalized_by_gc_bias = depth_16s / depth_ratio
                depth_16s_normalized_by_gc_bias = float("{0:.2f}".format(depth_16s_normalized_by_gc_bias))

                depth_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_mag, each_linked_16s,
                                                                 float("{0:.2f}".format(depth_mag)),
                                                                 depth_16s, depth_16s_normalized_by_gc_bias))
        else:
            print('Ignored %s, as %s not found, ' % (each_mag, gc_bias_txt))
    depth_txt_handle.close()


def vxtractor_output_to_dict(vxtractor_op, seq_len_dict, end_len_to_ignore, ignore_region_end_pct, min_region_len):

    seq_to_check = ''

    # initialize var con region id list
    var_region_id_list     = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9']
    con_region_id_list     = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10']
    var_con_region_id_list = ['C1', 'V1', 'C2', 'V2', 'C3', 'V3', 'C4', 'V4', 'C5', 'V5', 'C6', 'V6', 'C7', 'V7', 'C8', 'V8', 'C9', 'V9', 'C10']

    # read in variable regions
    seq_var_con_region_dict = {}
    for each_line in open(vxtractor_op):
        if not each_line.startswith(('Command:', 'Options:', 'Sequence,')):
            each_line_no_quote_symbol = ''
            for each_char in each_line.strip():
                if each_char not in ['"', '\'', '']:
                    each_line_no_quote_symbol += each_char

            each_line_split = each_line_no_quote_symbol.strip().split(',')
            seq_id = each_line_split[0]
            seq_var_con_region_dict[seq_id] = {}
            seq_var_regions = each_line_split[37:46]

            for (var_id, var_range) in zip(var_region_id_list, seq_var_regions):
                var_range_to_store = 'na'
                if var_range not in ['', 'notfound']:
                    if 'wrongorder' not in var_range:
                        if 'HMM=' not in var_range:
                            var_range_to_store = sorted([int(i) for i in var_range.split('-')])
                seq_var_con_region_dict[seq_id][var_id] = var_range_to_store

    # get conserved regions
    for each_seq in seq_var_con_region_dict:
        current_var_con_dict = seq_var_con_region_dict.get(each_seq, dict())
        for each_con_index in range(1, 11):

            ##### get C1 #####
            if each_con_index == 1:
                if current_var_con_dict['V1'] == 'na':
                    seq_var_con_region_dict[each_seq]['C1'] = 'na'
                else:
                    c1_left = 1
                    c1_right = current_var_con_dict['V1'][0] - 1
                    if c1_right >= 1:
                        seq_var_con_region_dict[each_seq]['C1'] = [c1_left, c1_right]
                    else:
                        seq_var_con_region_dict[each_seq]['C1'] = 'na'

            ##### get C2 - C9 #####
            elif 1 < each_con_index < 10:
                con_id = 'C%s' % each_con_index
                var_id_left = 'V%s' % (each_con_index - 1)
                var_id_right = 'V%s' % each_con_index
                if (current_var_con_dict[var_id_left] == 'na') or (current_var_con_dict[var_id_right] == 'na'):
                    seq_var_con_region_dict[each_seq][con_id] = 'na'
                else:
                    con_left = current_var_con_dict[var_id_left][1] + 1
                    con_right = current_var_con_dict[var_id_right][0] - 1
                    if con_right - con_left >= 1:
                        seq_var_con_region_dict[each_seq][con_id] = [con_left, con_right]
                    else:
                        seq_var_con_region_dict[each_seq][con_id] = 'na'

            ##### get C10 #####
            else:
                if current_var_con_dict['V9'] == 'na':
                    seq_var_con_region_dict[each_seq]['C10'] = 'na'
                else:
                    c10_left = current_var_con_dict['V9'][1] + 1
                    c10_right = seq_len_dict.get(each_seq, 0)
                    if c10_right - c10_left >= 1:
                        seq_var_con_region_dict[each_seq]['C10'] = [c10_left, c10_right]
                    else:
                        seq_var_con_region_dict[each_seq]['C10'] = 'na'

    # filter seq_var_con_region_dict
    seq_var_con_region_dict_filtered = {}
    for each_seq in seq_var_con_region_dict:
        current_seq_len = seq_len_dict.get(each_seq, 0)
        seq_var_con_region_dict_filtered[each_seq] = {}
        current_var_con_region_dict = seq_var_con_region_dict[each_seq]
        for each_con_var in var_con_region_id_list:
            con_var_range = current_var_con_region_dict[each_con_var]

            if con_var_range == 'na':
                seq_var_con_region_dict_filtered[each_seq][each_con_var] = 'na'
            else:
                pos_l = con_var_range[0]
                pos_r = con_var_range[1]
                region_len = pos_r - pos_l + 1
                region_end_len_to_ignore = round(region_len * ignore_region_end_pct / 100)
                pos_l_no_ignored = pos_l + region_end_len_to_ignore
                pos_r_no_ignored = pos_r - region_end_len_to_ignore

                region_to_consider = False
                pos_l_no_ignored_no_seq_end = pos_l_no_ignored
                pos_r_no_ignored_no_seq_end = pos_r_no_ignored
                if pos_l_no_ignored < end_len_to_ignore:
                    pos_l_no_ignored_no_seq_end = end_len_to_ignore + 1
                    if (pos_r_no_ignored - pos_l_no_ignored_no_seq_end) >= min_region_len:
                        region_to_consider = True
                elif pos_r_no_ignored > (current_seq_len - end_len_to_ignore):
                    pos_r_no_ignored_no_seq_end = current_seq_len - end_len_to_ignore
                    if (pos_r_no_ignored_no_seq_end - pos_l_no_ignored) >= min_region_len:
                        region_to_consider = True
                else:
                    if (pos_r_no_ignored - pos_l_no_ignored + 1) >= min_region_len:
                        region_to_consider = True

                if region_to_consider is False:
                    if each_seq == seq_to_check:
                        print('%s\t%s\t%s-%s\t%s\t%s-%s\t%s\t%s' % (each_seq, each_con_var, end_len_to_ignore, (current_seq_len - end_len_to_ignore), con_var_range, pos_l_no_ignored_no_seq_end, pos_r_no_ignored_no_seq_end, (pos_r_no_ignored_no_seq_end - pos_l_no_ignored_no_seq_end + 1), 'bad'))
                    seq_var_con_region_dict_filtered[each_seq][each_con_var] = 'na'
                else:
                    if each_seq == seq_to_check:
                        print('%s\t%s\t%s-%s\t%s\t%s-%s\t%s\t%s' % (each_seq, each_con_var, end_len_to_ignore, (current_seq_len - end_len_to_ignore), con_var_range, pos_l_no_ignored_no_seq_end, pos_r_no_ignored_no_seq_end, (pos_r_no_ignored_no_seq_end - pos_l_no_ignored_no_seq_end + 1), 'good'))
                    seq_var_con_region_dict_filtered[each_seq][each_con_var] = [pos_l_no_ignored_no_seq_end, pos_r_no_ignored_no_seq_end]

    # ignore region if its two flankings are na
    seq_var_con_region_dict_filtered2 = {}
    for each_seq in seq_var_con_region_dict_filtered:
        current_var_con_region_dict = seq_var_con_region_dict_filtered[each_seq]
        current_var_con_region_dict2 = {}

        # check vars
        for each_var in var_region_id_list:
            var_index = int(each_var[1:])
            flk_con_l_id = 'C%s' % var_index
            flk_con_r_id = 'C%s' % (var_index + 1)
            flk_con_l = current_var_con_region_dict.get(flk_con_l_id, 'na')
            flk_con_r = current_var_con_region_dict.get(flk_con_r_id, 'na')
            if (flk_con_l == 'na') and (flk_con_r == 'na'):
                current_var_con_region_dict2[each_var] = 'na'
            else:
                current_var_con_region_dict2[each_var] = current_var_con_region_dict[each_var]

        # check cons
        for each_con in con_region_id_list:
            con_index = int(each_con[1:])
            flk_var_l_id = 'V%s' % (con_index - 1)
            flk_var_r_id = 'V%s' % (con_index)
            flk_var_l = current_var_con_region_dict.get(flk_var_l_id, 'na')
            flk_var_r = current_var_con_region_dict.get(flk_var_r_id, 'na')
            if (flk_var_l == 'na') and (flk_var_r == 'na'):
                current_var_con_region_dict2[each_con] = 'na'
            else:
                current_var_con_region_dict2[each_con] = current_var_con_region_dict[each_con]

        seq_var_con_region_dict_filtered2[each_seq] = current_var_con_region_dict2

    seq_var_con_combined_pos_dict = {}
    for each_seq in seq_var_con_region_dict_filtered2:
        current_var_con_region_dict = seq_var_con_region_dict_filtered2[each_seq]
        seq_var_con_combined_pos_dict[each_seq] = {}
        seq_var_con_combined_pos_dict[each_seq]['Con'] = set()
        seq_var_con_combined_pos_dict[each_seq]['Var'] = set()
        for each_con_var in current_var_con_region_dict:
            current_var_con_range = current_var_con_region_dict[each_con_var]
            if current_var_con_range != 'na':

                # get var pos list
                if each_con_var in var_region_id_list:
                    current_var_pos_list = list(range((current_var_con_range[0] - 1), current_var_con_range[1]))
                    for each_var_pos in current_var_pos_list:
                        seq_var_con_combined_pos_dict[each_seq]['Var'].add(each_var_pos)

                # get con pos list
                if each_con_var in con_region_id_list:
                    current_con_pos_list = list(range((current_var_con_range[0] - 1), current_var_con_range[1]))
                    for each_con_pos in current_con_pos_list:
                        seq_var_con_combined_pos_dict[each_seq]['Con'].add(each_con_pos)

    return seq_var_con_combined_pos_dict


def opt_estimation(depth_file_all, depth_file_linked, pos_cov_file_all, pos_cov_file_linked, vxtractor_op, end_len_to_ignore, ignore_region_end_pct, min_region_len_no_ignored, seq_len_dict, depth_file_opt):

    # read in mean depth (all 16S)
    marker_to_mag_dict = {}
    mean_depth_mag_dict = {}
    mean_depth_16s_all = {}
    mean_depth_16s_all_norm = {}
    for each_16s in open(depth_file_all):
        if not each_16s.startswith('MAG	16S	MAG_depth	16S_depth'):
            each_16s_split = each_16s.strip().split('\t')
            id_mag = each_16s_split[0]
            id_16s = each_16s_split[1]
            mag_depth = float(each_16s_split[2])
            mean_depth = float(each_16s_split[3])
            mean_depth_norm = float(each_16s_split[4])
            mean_depth_mag_dict[id_mag] = mag_depth
            mean_depth_16s_all[id_16s] = mean_depth
            mean_depth_16s_all_norm[id_16s] = mean_depth_norm
            marker_to_mag_dict[id_16s] = id_mag

    # read in mean depth (linked 16S)
    mean_depth_16s_linked = {}
    mean_depth_16s_linked_norm = {}
    for each_16s in open(depth_file_linked):
        if not each_16s.startswith('MAG	16S	MAG_depth	16S_depth'):
            each_16s_split = each_16s.strip().split('\t')
            id_16s = each_16s_split[1]
            mean_depth = float(each_16s_split[3])
            mean_depth_norm = float(each_16s_split[4])
            mean_depth_16s_linked[id_16s] = mean_depth
            mean_depth_16s_linked_norm[id_16s] = mean_depth_norm

    # read in depth by bp (all 16S)
    depth_by_bp_dict_all = {}
    for each_bp_all in open(pos_cov_file_all):
        each_bp_all_split = each_bp_all.strip().split('\t')
        ref_if = each_bp_all_split[0]
        pos = each_bp_all_split[1]
        depth = int(each_bp_all_split[2])
        if ref_if not in depth_by_bp_dict_all:
            depth_by_bp_dict_all[ref_if] = {}
        depth_by_bp_dict_all[ref_if][pos] = depth

    # read in depth by bp (linked 16S)
    depth_by_bp_dict_linked = {}
    for each_bp_all in open(pos_cov_file_linked):
        each_bp_all_split = each_bp_all.strip().split('\t')
        ref_if = each_bp_all_split[0]
        pos = each_bp_all_split[1]
        depth = int(each_bp_all_split[2])
        if ref_if not in depth_by_bp_dict_linked:
            depth_by_bp_dict_linked[ref_if] = {}
        depth_by_bp_dict_linked[ref_if][pos] = depth

    seq_var_con_combined_pos_dict = vxtractor_output_to_dict(vxtractor_op, seq_len_dict, end_len_to_ignore, ignore_region_end_pct, min_region_len_no_ignored)

    estimation_with_depth_opt_handle = open(depth_file_opt, 'w')
    estimation_with_depth_opt_handle.write('MAG\t16S\tCoverage(cov_variable/cov_conserved)(linked)\tCoverage(cov_variable/cov_conserved)(all)\tCoverage_optimal\tCopies_optimal\tCopies_optimal_normalized\n')
    for each_seq in seq_var_con_combined_pos_dict:
        seq_var_pos_list = sorted([i for i in seq_var_con_combined_pos_dict[each_seq]['Var']])
        seq_con_pos_list = sorted([i for i in seq_var_con_combined_pos_dict[each_seq]['Con']])

        seq_var_pos_total_depth_all = 0
        for each_var_pos in seq_var_pos_list:
            var_pos_depth_all = depth_by_bp_dict_all.get(each_seq, dict()).get(str(each_var_pos), 0)
            seq_var_pos_total_depth_all += var_pos_depth_all

        seq_con_pos_total_depth_all = 0
        for each_con_pos in seq_con_pos_list:
            con_pos_depth_all = depth_by_bp_dict_all.get(each_seq, dict()).get(str(each_con_pos), 0)
            seq_con_pos_total_depth_all += con_pos_depth_all

        seq_var_pos_total_depth_linked = 0
        for each_var_pos in seq_var_pos_list:
            var_pos_depth_linked = depth_by_bp_dict_linked.get(each_seq, dict()).get(str(each_var_pos), 0)
            seq_var_pos_total_depth_linked += var_pos_depth_linked

        seq_con_pos_total_depth_linked = 0
        for each_con_pos in seq_con_pos_list:
            con_pos_depth_linked = depth_by_bp_dict_linked.get(each_seq, dict()).get(str(each_con_pos), 0)
            seq_con_pos_total_depth_linked += con_pos_depth_linked

        var_mean_depth_all = 0
        var_mean_depth_linked = 0
        if len(seq_var_pos_list) != 0:
            var_mean_depth_all = seq_var_pos_total_depth_all/len(seq_var_pos_list)
            var_mean_depth_linked = seq_var_pos_total_depth_linked/len(seq_var_pos_list)

        con_mean_depth_all = 0
        con_mean_depth_linked = 0
        if len(seq_con_pos_list) != 0:
            con_mean_depth_all = seq_con_pos_total_depth_all/len(seq_con_pos_list)
            con_mean_depth_linked = seq_con_pos_total_depth_linked/len(seq_con_pos_list)

        var_to_con_cov_ratio_all = 0
        if con_mean_depth_all > 0:
            var_to_con_cov_ratio_all = var_mean_depth_all/con_mean_depth_all
            var_to_con_cov_ratio_all = float("{0:.2f}".format(var_to_con_cov_ratio_all))

        var_to_con_cov_ratio_linked = 0
        if con_mean_depth_linked > 0:
            var_to_con_cov_ratio_linked = var_mean_depth_linked/con_mean_depth_linked
            var_to_con_cov_ratio_linked = float("{0:.2f}".format(var_to_con_cov_ratio_linked))

        if (var_to_con_cov_ratio_linked != 0) and (var_to_con_cov_ratio_all != 0):

            seq_cov_by_linked      = mean_depth_16s_linked.get(each_seq, 0)
            seq_cov_by_linked_norm = mean_depth_16s_linked_norm.get(each_seq, 0)
            seq_cov_by_all         = mean_depth_16s_all.get(each_seq, 0)
            seq_cov_by_all_norm    = mean_depth_16s_all_norm.get(each_seq, 0)
            linked_mag             = marker_to_mag_dict[each_seq]
            mag_cov                = mean_depth_mag_dict[linked_mag]

            # chose the higher one
            if (0.85 <= var_to_con_cov_ratio_linked <= 1.15) and (0.85 <= var_to_con_cov_ratio_all <= 1.15):
                seq_cov_opt = max([seq_cov_by_linked, seq_cov_by_all])
                seq_cov_opt_norm = max([seq_cov_by_linked_norm, seq_cov_by_all_norm])

            # reads from other genome mapped to conserved regions
            elif (var_to_con_cov_ratio_linked < 0.85) and (var_to_con_cov_ratio_all < 0.85):
                seq_cov_opt = seq_cov_by_all*var_to_con_cov_ratio_all
                seq_cov_opt_norm = seq_cov_by_all_norm*var_to_con_cov_ratio_all

            # reads from conserved regions mapped to other genomes
            elif (var_to_con_cov_ratio_linked > 1.15) and (var_to_con_cov_ratio_all > 1.15):
                seq_cov_opt = seq_cov_by_all/var_to_con_cov_ratio_all
                seq_cov_opt_norm = seq_cov_by_all_norm/var_to_con_cov_ratio_all
            else:
                distance_to_1_linked = abs(var_to_con_cov_ratio_linked - 1)
                distance_to_1_all = abs(var_to_con_cov_ratio_all - 1)
                seq_cov_opt = seq_cov_by_linked/var_to_con_cov_ratio_linked
                seq_cov_opt_norm = seq_cov_by_linked_norm/var_to_con_cov_ratio_linked
                if distance_to_1_linked > distance_to_1_all:
                    seq_cov_opt = seq_cov_by_all/var_to_con_cov_ratio_all
                    seq_cov_opt_norm = seq_cov_by_all_norm/var_to_con_cov_ratio_all
                seq_cov_opt = float("{0:.2f}".format(seq_cov_opt))
                seq_cov_opt_norm = float("{0:.2f}".format(seq_cov_opt_norm))

            seq_cp_num_opt = float("{0:.2f}".format(seq_cov_opt/mag_cov))
            seq_cp_num_opt_norm = float("{0:.2f}".format(seq_cov_opt_norm/mag_cov))
            estimation_with_depth_opt_handle.write('%s\t%s\t%s(%s)\t%s(%s)\t%s\t%s\t%s\n' % (linked_mag, each_seq, seq_cov_by_linked, var_to_con_cov_ratio_linked, seq_cov_by_all, var_to_con_cov_ratio_all, seq_cov_opt, seq_cp_num_opt, seq_cp_num_opt_norm))
    estimation_with_depth_opt_handle.close()


def keep_best_blast_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def rename_seq(ctg_file_in, ctg_file_out, prefix, str_connector):

    ctg_file_out_handle = open(ctg_file_out, 'w')
    for Seq_record in SeqIO.parse(ctg_file_in, 'fasta'):
        Seq_record.id = '%s%s%s' % (prefix, str_connector, Seq_record.id)
        SeqIO.write(Seq_record, ctg_file_out_handle, 'fasta')
    ctg_file_out_handle.close()


def run_vxtractor_on_input_16s(input_16s, silva_ref_seq, get_16s_domain_wd, num_threads, vxtractor_perl, hmm_folder_bac,
                               hmm_folder_arc, vxtractor_op_combined):
    pwd_silva_16s = '%s/SILVA_16S.fa' % get_16s_domain_wd
    makeblastdb_log = '%s/SILVA_16S_makeblastdb.log' % get_16s_domain_wd
    blast_op = '%s/Input_16S_vs_SILVA.tab' % get_16s_domain_wd
    blast_op_best_hit = '%s/Input_16S_vs_SILVA_best_hit.tab' % get_16s_domain_wd
    input_16s_bac = '%s/Input_16S_bac.fa' % get_16s_domain_wd
    input_16s_arc = '%s/Input_16S_arc.fa' % get_16s_domain_wd
    vxtractor_op_csv_bac = '%s/vxtractor_op_bac.csv' % get_16s_domain_wd
    vxtractor_op_csv_arc = '%s/vxtractor_op_arc.csv' % get_16s_domain_wd

    os.system('cp %s %s' % (silva_ref_seq, pwd_silva_16s))

    # read in ref taxon
    silva_seq_domain_dict = {}
    for each_seq in SeqIO.parse(pwd_silva_16s, 'fasta'):
        seq_description = each_seq.description
        seq_des_split = seq_description.split(' ')
        seq_id = seq_des_split[0]
        taxon_list = ' '.join(seq_description.split(' ')[1:]).split(';')
        silva_seq_domain_dict[seq_id] = taxon_list[0]

    make_blastdb_cmd = 'makeblastdb -in %s -dbtype nucl -parse_seqids -logfile %s' % (pwd_silva_16s, makeblastdb_log)
    os.system(make_blastdb_cmd)
    blastn_cmd = 'blastn -query %s -db %s -out %s -outfmt 6 -num_threads %s' % (input_16s, pwd_silva_16s, blast_op, num_threads)
    os.system(blastn_cmd)

    keep_best_blast_hit(blast_op, blast_op_best_hit)

    query_domain_dict = {}
    for each_hit in open(blast_op_best_hit):
        each_hit_split = each_hit.strip().split('\t')
        query_id = each_hit_split[0]
        ref_id = each_hit_split[1]
        ref_domain = silva_seq_domain_dict[ref_id]
        query_domain_dict[query_id] = ref_domain

    bac_16s_num = 0
    arc_16s_num = 0
    input_16s_bac_handle = open(input_16s_bac, 'w')
    input_16s_arc_handle = open(input_16s_arc, 'w')
    for each_16s in SeqIO.parse(input_16s, 'fasta'):
        seq_domain = query_domain_dict.get(each_16s.id, 'Bacteria')
        if seq_domain == 'Bacteria':
            input_16s_bac_handle.write('>%s\n' % each_16s.id)
            input_16s_bac_handle.write('%s\n' % str(each_16s.seq))
            bac_16s_num += 1
        if seq_domain == 'Archaea':
            input_16s_arc_handle.write('>%s\n' % each_16s.id)
            input_16s_arc_handle.write('%s\n' % str(each_16s.seq))
            arc_16s_num += 1
    input_16s_bac_handle.close()
    input_16s_arc_handle.close()

    vxtractor_cmd_bac = 'perl %s -a -h %s -c %s -o /dev/null %s &> /dev/null' % (vxtractor_perl, hmm_folder_bac, vxtractor_op_csv_bac, input_16s_bac)
    vxtractor_cmd_arc = 'perl %s -a -h %s -c %s -o /dev/null %s &> /dev/null' % (vxtractor_perl, hmm_folder_arc, vxtractor_op_csv_arc, input_16s_arc)
    if bac_16s_num > 0:
        os.system(vxtractor_cmd_bac)
    if arc_16s_num > 0:
        os.system(vxtractor_cmd_arc)

    cat_cmd = ''
    if (bac_16s_num > 0) and (arc_16s_num > 0):
        cat_cmd = 'cat %s %s > %s' % (vxtractor_op_csv_bac, vxtractor_op_csv_arc, vxtractor_op_combined)
    elif bac_16s_num > 0:
        cat_cmd = 'cat %s > %s' % (vxtractor_op_csv_bac, vxtractor_op_combined)
    elif arc_16s_num > 0:
        cat_cmd = 'cat %s > %s' % (vxtractor_op_csv_arc, vxtractor_op_combined)
    os.system(cat_cmd)

    # check if output empty
    processed_seq_num = 0
    for each_line in open(vxtractor_op_combined):
        if not each_line.startswith(('Command:', 'Options:', 'Sequence,')):
            processed_seq_num += 1
    if processed_seq_num == 0:
        print('vxtractor failed, program exited!')
        exit()


def subsample_paired_reads(seq_in_r1, seq_in_r2, subsample_pct, seq_out_r1, seq_out_r2, num_threads):

    # get the number of paired reads
    input_r1_line_num_file = '%s_line_num.txt' % (seq_out_r1)
    input_r2_line_num_file = '%s_line_num.txt' % (seq_out_r2)

    os.system('wc -l %s > %s' % (seq_in_r1, input_r1_line_num_file))
    os.system('wc -l %s > %s' % (seq_in_r2, input_r2_line_num_file))
    paired_r1_num = int(int(open(input_r1_line_num_file).readline().strip().split(' ')[0]) / 2)
    paired_r2_num = int(int(open(input_r2_line_num_file).readline().strip().split(' ')[0]) / 2)
    if seq_in_r1[-1] in ['Q', 'q']:
        paired_r1_num = int(int(open(input_r1_line_num_file).readline().strip().split(' ')[0]) / 4)
        paired_r2_num = int(int(open(input_r2_line_num_file).readline().strip().split(' ')[0]) / 4)

    if paired_r1_num != paired_r2_num:
        print('Inconsistent number of reads found in r1 and r2, program exited!')
        exit()

    # remove tmp files
    os.remove(input_r1_line_num_file)
    os.remove(input_r2_line_num_file)

    # get the number of reads paired to subset
    to_extract_reads_num = round(paired_r1_num*subsample_pct/100)

    # remember to use the same random seed to keep pairing
    subsample_r1_cmd = 'seqtk sample -s100 %s %s > %s' % (seq_in_r1, to_extract_reads_num, seq_out_r1)
    subsample_r2_cmd = 'seqtk sample -s100 %s %s > %s' % (seq_in_r2, to_extract_reads_num, seq_out_r2)

    # subsample with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, [subsample_r1_cmd, subsample_r2_cmd])
    pool.close()
    pool.join()


def get_16s_copy_num(args):

    ############################################### file and parameter in ##############################################

    op_folder                           = args['o']
    op_prefix                           = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    matam_16s_reads                     = args['r16s']
    mag_folder                          = args['mag']
    mag_file_ext                        = args['x']
    combined_prefixed_mags              = args['cp_mags']
    combined_prefixed_mags_gff          = args['cp_mag_gff']
    sorted_read_to_16s_sam_file         = args['sam_16s_sorted_by_read']
    identified_linkage_gnm_level        = args['linkages']
    marker_seq_file                     = args['marker']
    min_M_len                           = args['aln_len']
    mismatch_cutoff                     = args['mismatch']
    provided_cp_num_txt                 = args['ref_16s_cp_num']
    both_pair_mapped                    = args['both_pair_mapped']
    ignore_ends_len_16s                 = args['ignore_ends_len_16s']
    pos_pct_to_ignore_lowest            = args['ignore_lowest_pct']
    pos_pct_to_ignore_highest           = args['ignore_highest_pct']
    ignore_gc_bias                      = args['ignore_gc_bias']
    min_insert_size_16s                 = int(args['min_insert_size_16s'])
    mag_depth_gc_txt                    = args['mag_cov_gc']
    mag_gc_bias_folder                  = args['mag_gc_bias']
    num_threads                         = args['t']
    force_overwrite                     = args['force']
    keep_quiet                          = args['quiet']
    vxtractor_perl                      = args['vxtractor']
    silva_ref_seq                       = args['silva_order_refs']
    hmm_folder_bac                      = args['hmm_bac']
    hmm_folder_arc                      = args['hmm_arc']
    reads_subsample_pct_for_mag_cov_gc  = args['subsample_pct']

    keep_uniq                           = False
    window_len                          = 100
    ignore_ends_len_mag                 = 500
    ctg_len_cutoff                      = 10000
    gnm_to_ctg_connector                = '___C___'
    end_len_to_ignore                   = 200
    ignore_region_end_pct               = 25
    min_region_len_no_ignored           = 25

    ################################################ check dependencies ################################################

    # check whether executables exist
    program_list = ['makeblastdb', 'blastn', 'bowtie2-build', 'bowtie2', 'samtools', 'seqtk', 'hmmscan']
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not found, program exited!' % ','.join(not_detected_programs))
        exit()

    ################################################### check output folder ###################################################

    op_folder_to_generate = '%s_get_16S_cp_num_wd' % op_prefix

    if os.path.isdir(op_folder) == '.':
        pwd_op_folder = op_folder_to_generate
    else:
        pwd_op_folder = '%s/%s' % (op_folder, op_folder_to_generate)

    if (os.path.isdir(pwd_op_folder) is True) and (force_overwrite is False):
        print('Output folder detected, please specify a different location or provide a different prefix')
        print('Program exited!')
        exit()
    else:
        force_create_folder(pwd_op_folder)

    ############################################## prefix and combine MAGs #############################################

    if (combined_prefixed_mags is None) or (combined_prefixed_mags_gff is None):

        if mag_folder is None:
            print('No MAG files or combined prefixed MAG file provided, program exited')
            exit()

        if find_executable('barrnap') is None:
            print('barrnap not found, program exited!' )
            exit()

        if mag_folder[-1] == '/':
            mag_folder = mag_folder[:-1]

        # get input mag file list
        mag_file_re = '%s/*%s' % (mag_folder, mag_file_ext)
        mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
        if len(mag_file_list) == 0:
            print('No MAG found, program exited!')
            exit()

        mag_folder_in_get_cn_wd             = '%s/input_MAGs'                           % pwd_op_folder
        prefixed_mag_folder                 = '%s/input_MAGs_prefixed'                  % mag_folder_in_get_cn_wd
        combined_input_gnms                 = '%s/input_MAGs_combined.fa'               % mag_folder_in_get_cn_wd
        barrnap_wd                          = '%s/input_MAGs/barrnap_wd'                % pwd_op_folder
        combined_barrnap_gff                = '%s/input_MAGs/combined_barrnap.gff'      % pwd_op_folder

        # create folder
        os.mkdir(mag_folder_in_get_cn_wd)
        os.mkdir(prefixed_mag_folder)

        # add mag id to its sequences
        argument_list_for_barrnap = []
        for mag_in in mag_file_list:
            mag_basename = '.'.join(mag_in.split('.')[:-1])
            pwd_mag_in = '%s/%s' % (mag_folder, mag_in)
            pwd_mag_renamed = '%s/%s' % (prefixed_mag_folder, mag_in)
            pwd_barrnap_ffn = '%s/%s.ffn' % (barrnap_wd, mag_basename)
            pwd_barrnap_gff = '%s/%s.gff' % (barrnap_wd, mag_basename)
            pwd_barrnap_log = '%s/%s.log' % (barrnap_wd, mag_basename)
            barrnap_cmd = 'barrnap --quiet -o %s %s > %s 2> %s' % (pwd_barrnap_ffn, pwd_mag_renamed, pwd_barrnap_gff, pwd_barrnap_log)
            argument_list_for_barrnap.append(barrnap_cmd)
            rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename, gnm_to_ctg_connector)

        # combine prefixed MAGs
        os.system('cat %s/*%s > %s' % (prefixed_mag_folder, mag_file_ext, combined_input_gnms))

        ########## run barrnap on prefixed MAGs ##########

        os.mkdir(barrnap_wd)

        # run barrnap with multiprocessing
        pool = mp.Pool(processes=num_threads)
        pool.map(os.system, argument_list_for_barrnap)
        pool.close()
        pool.join()

        # remove index files
        os.system('rm %s/*.fai' % prefixed_mag_folder)
        os.system('cat %s/*.gff > %s' % (barrnap_wd, combined_barrnap_gff))

        combined_prefixed_mags = combined_input_gnms
        combined_prefixed_mags_gff = combined_barrnap_gff


    ################################################### define file name ###################################################

    pwd_log_file = args.get('log_file', None)
    if pwd_log_file is None:
        pwd_log_file = '%s/%s.log' % (pwd_op_folder, op_prefix)

    file_path_mag, combined_prefixed_mags_no_ext, ext_mag = sep_path_basename_ext(combined_prefixed_mags)
    marker_seq_path, marker_seq_no_ext, marker_seq_ext = sep_path_basename_ext(marker_seq_file)

    reads_file_r1_subset                            = '%s/%s_R1_subset.fa'                          % (pwd_op_folder, op_prefix)
    reads_file_r2_subset                            = '%s/%s_R2_subset.fa'                          % (pwd_op_folder, op_prefix)
    get_16s_domain_folder                           = '%s/%s_get_16s_domain'                        % (pwd_op_folder, op_prefix)
    index_folder                                    = '%s/%s_index'                                 % (pwd_op_folder, op_prefix)
    marker_seq_index                                = '%s/%s'                                       % (index_folder, marker_seq_no_ext)
    sam_file_reads_to_mag                           = '%s/%s_%s.sam'                                % (pwd_op_folder, op_prefix, combined_prefixed_mags_no_ext)
    sam_file_reads_to_mag_filtered                  = '%s/%s_%s_%s_%s.sam'                          % (pwd_op_folder, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    sam_file_reads_to_mag_filtered_both_mapped      = '%s/%s_%s_%s_%s_both_mapped.sam'              % (pwd_op_folder, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    sam_file_reads_to_mag_filtered_sorted           = '%s/%s_%s_%s_%s_both_mapped_sorted.sam'       % (pwd_op_folder, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    sam_file_reads_to_mag_filtered_sorted_depth     = '%s/%s_%s_%s_%s_both_mapped_sorted_depth.txt' % (pwd_op_folder, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    mag_gc_bias_op_folder                           = '%s/%s_GC_bias'                               % (pwd_op_folder, op_prefix)
    mag_depth_GC_content_folder                     = '%s/%s_MAG_depth_GC_content'                  % (pwd_op_folder, op_prefix)
    mag_depth_GC_content_file                       = '%s/%s_MAG_depth_GC_content.txt'              % (pwd_op_folder, op_prefix)
    vxtractor_combined_output                       = '%s/%s_16S_vxtractor.csv'                     % (pwd_op_folder, op_prefix)
    all_16s_sam_file                                = '%s/%s_all16S.sam'                            % (pwd_op_folder, op_prefix)
    all_16s_sam_file_sorted                         = '%s/%s_all16S_sorted_by_read.sam'             % (pwd_op_folder, op_prefix)
    all_16s_sam_file_sorted_line_num                = '%s/%s_all16S_sorted_by_read_line_num.txt'    % (pwd_op_folder, op_prefix)
    all_16s_sam_file_sorted_split                   = '%s/%s_all16S_sorted_by_read_split'           % (pwd_op_folder, op_prefix)
    all_16s_sam_mp_folder                           = '%s/%s_all16S_mr'                             % (pwd_op_folder, op_prefix)
    linked_16s_sam_mp_folder                        = '%s/%s_linked_16S_mr'                         % (pwd_op_folder, op_prefix)
    all_16s_sam_file_filtered                       = '%s/%s_all16S.sam'                            % (pwd_op_folder, op_prefix)
    all_16s_sam_file_filtered_random                = '%s/%s_all16S_random.sam'                     % (pwd_op_folder, op_prefix)
    all_16s_sam_file_filtered_random_sorted         = '%s/%s_all16S_random_sorted.sam'              % (pwd_op_folder, op_prefix)
    linked_16s_sam_file_filtered                    = '%s/%s_linked16S.sam'                         % (pwd_op_folder, op_prefix)
    linked_16s_sam_file_filtered_random             = '%s/%s_linked16S_random.sam'                  % (pwd_op_folder, op_prefix)
    linked_16s_sam_file_filtered_random_sorted      = '%s/%s_linked16S_random_sorted.sam'           % (pwd_op_folder, op_prefix)
    depth_GC_content_file_all                       = '%s/%s_16S_cov_GC_content_all.txt'            % (pwd_op_folder, op_prefix)
    depth_GC_content_file_linked                    = '%s/%s_16S_cov_GC_content_linked.txt'         % (pwd_op_folder, op_prefix)
    depth_by_bp_file_all                            = '%s/%s_16S_pos_cov_all.txt'                   % (pwd_op_folder, op_prefix)
    depth_by_bp_file_linked                         = '%s/%s_16S_pos_cov_linked.txt'                % (pwd_op_folder, op_prefix)
    counted_pos_txt_all                             = '%s/%s_16S_pos_counted_all.txt'               % (pwd_op_folder, op_prefix)
    counted_pos_txt_linked                          = '%s/%s_16S_pos_counted_linked.txt'            % (pwd_op_folder, op_prefix)
    depth_16s_all                                   = '%s/%s_MAG_16S_cov_all.txt'                   % (pwd_op_folder, op_prefix)
    depth_16s_linked                                = '%s/%s_MAG_16S_cov_linked.txt'                % (pwd_op_folder, op_prefix)
    estimated_cp_num_txt_by_16s                     = '%s/%s_copy_num_by_16S.txt'                   % (pwd_op_folder, op_prefix)
    estimated_cp_num_txt_by_mag                     = '%s/%s_copy_num_by_MAG.txt'                   % (pwd_op_folder, op_prefix)
    scatter_plot                                    = '%s/%s_copy_num_by_MAG.svg'                   % (pwd_op_folder, op_prefix)
    scatter_plot_norm                               = '%s/%s_copy_num_by_MAG_norm.svg'              % (pwd_op_folder, op_prefix)

    os.mkdir(index_folder)

    ############################################## run_vxtractor_on_input_16s ##############################################

    report_and_log(('Get_cn: running V-Xtractor on 16S rRNA genes'), pwd_log_file, keep_quiet)

    os.mkdir(get_16s_domain_folder)
    run_vxtractor_on_input_16s(marker_seq_file, silva_ref_seq, get_16s_domain_folder, num_threads, vxtractor_perl, hmm_folder_bac,
                               hmm_folder_arc, vxtractor_combined_output)

    ########################################################################################################################
    ################################################### estimate GC bias ###################################################
    ########################################################################################################################

    ########################################## map reads to combined prefixed MAGs #########################################

    if mag_depth_gc_txt is None:

        if reads_subsample_pct_for_mag_cov_gc != 100:
            report_and_log(('Get_cn: subsampling reads for MAG coverage estimation'), pwd_log_file, keep_quiet)
            subsample_paired_reads(reads_file_r1, reads_file_r2, reads_subsample_pct_for_mag_cov_gc, reads_file_r1_subset, reads_file_r2_subset, num_threads)

        report_and_log(('Get_cn: mapping reads to combined prefixed MAGs'), pwd_log_file, keep_quiet)
        bowtie_build_cmd = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, combined_prefixed_mags, index_folder, combined_prefixed_mags_no_ext)
        report_and_log((bowtie_build_cmd), pwd_log_file, True)
        os.system(bowtie_build_cmd)

        if reads_subsample_pct_for_mag_cov_gc == 100:
            #mapping_cmd_mag = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s -f --xeq --no-mixed --no-discordant --no-1mm-upfront --no-unal &> /dev/null' % (index_folder, combined_prefixed_mags_no_ext, reads_file_r1, reads_file_r2, sam_file_reads_to_mag, num_threads)
            mapping_cmd_mag = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s -f --xeq --no-unal &> /dev/null' % (index_folder, combined_prefixed_mags_no_ext, reads_file_r1, reads_file_r2, sam_file_reads_to_mag, num_threads)
        else:
            #mapping_cmd_mag = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s -f --xeq --no-mixed --no-discordant --no-1mm-upfront --no-unal &> /dev/null' % (index_folder, combined_prefixed_mags_no_ext, reads_file_r1_subset, reads_file_r2_subset, sam_file_reads_to_mag, num_threads)
            mapping_cmd_mag = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s -f --xeq --no-unal &> /dev/null' % (index_folder, combined_prefixed_mags_no_ext, reads_file_r1_subset, reads_file_r2_subset, sam_file_reads_to_mag, num_threads)
        report_and_log((mapping_cmd_mag), pwd_log_file, True)
        os.system(mapping_cmd_mag)

        # remove reads subsets
        if reads_subsample_pct_for_mag_cov_gc != 100:
            os.remove(reads_file_r1_subset)
            os.remove(reads_file_r2_subset)

        # filter sam file
        report_and_log(('Get_cn: removing high mismatch alignments'), pwd_log_file, keep_quiet)
        remove_high_mismatch(sam_file_reads_to_mag, min_M_len, mismatch_cutoff, sam_file_reads_to_mag_filtered)
        os.remove(sam_file_reads_to_mag)

        #keep_only_both_mapped_reads
        report_and_log(('Get_cn: removing singletons'), pwd_log_file, keep_quiet)
        keep_only_both_mapped_reads(sam_file_reads_to_mag_filtered, sam_file_reads_to_mag_filtered_both_mapped)
        os.remove(sam_file_reads_to_mag_filtered)

        # sort filtered_sam
        report_and_log(('Get_cn: sorting filtered sam file'), pwd_log_file, keep_quiet)
        #sort_cmd_mag = 'samtools sort -O sam --threads %s -o %s %s &> /dev/null' % (num_threads, sam_file_reads_to_mag_filtered_sorted, sam_file_reads_to_mag)
        sort_cmd_mag = 'samtools sort -O sam --threads %s -o %s %s &> /dev/null' % (num_threads, sam_file_reads_to_mag_filtered_sorted, sam_file_reads_to_mag_filtered_both_mapped)
        report_and_log((sort_cmd_mag), pwd_log_file, True)
        os.system(sort_cmd_mag)
        #os.remove(sam_file_reads_to_mag)
        os.remove(sam_file_reads_to_mag_filtered_both_mapped)

        # get depth file
        report_and_log(('Get_cn: getting depth with samtools'), pwd_log_file, keep_quiet)
        get_depth_cmd_mag = 'samtools depth -a %s > %s' % (sam_file_reads_to_mag_filtered_sorted, sam_file_reads_to_mag_filtered_sorted_depth)
        report_and_log((get_depth_cmd_mag), pwd_log_file, True)
        os.system(get_depth_cmd_mag)
        #os.remove(sam_file_reads_to_mag_filtered_sorted)

    ############################################ read in 16S sequences into dict ###########################################

    marker_seq_dict = {}
    marker_len_dict = {}
    marker_gc_content_dict = {}
    for each_marker in SeqIO.parse(marker_seq_file, 'fasta'):
        marker_id = each_marker.id
        marker_seq = str(each_marker.seq).upper()
        marker_GC_content = ((marker_seq.count('G')) + (marker_seq.count('C'))) * 100 / len(marker_seq)
        marker_GC_content = float("{0:.2f}".format(marker_GC_content))
        marker_seq_dict[marker_id] = marker_seq
        marker_len_dict[marker_id] = len(marker_seq)
        marker_gc_content_dict[marker_id] = marker_GC_content

    ############################################### read in MarkerMAG output ###############################################

    gnm_to_linked_16s_gc_content_dict = {}
    gnm_to_linked_16s_dict = {}
    linked_16s_to_gnm_dict = {}
    for each_linkage in open(identified_linkage_gnm_level):
        each_linkage_split = each_linkage.strip().split('\t')
        if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Round'):
            id_16s = each_linkage_split[0]
            id_gnm = each_linkage_split[1]
            gc_content_16s = marker_gc_content_dict[id_16s]
            if id_gnm not in gnm_to_linked_16s_dict:
                gnm_to_linked_16s_dict[id_gnm] = {id_16s}
                gnm_to_linked_16s_gc_content_dict[id_gnm] = [gc_content_16s]
            else:
                gnm_to_linked_16s_dict[id_gnm].add(id_16s)
                gnm_to_linked_16s_gc_content_dict[id_gnm].append(gc_content_16s)
            linked_16s_to_gnm_dict[id_16s] = id_gnm

    unlinked_16s_set = set()
    for each_16s in marker_len_dict:
        if each_16s not in linked_16s_to_gnm_dict:
            unlinked_16s_set.add(each_16s)

    ############################################ read in MAG sequences into dict ###########################################

    report_and_log(('Get_cn: reading in MAG sequences'), pwd_log_file, keep_quiet)

    ctg_seq_dict = {}
    ctg_len_dict = {}
    mag_to_ctg_dict = {}
    mag_total_len_and_gc_count_dict = {}
    for each_ctg in SeqIO.parse(combined_prefixed_mags, 'fasta'):
        ctg_id = each_ctg.id
        mag_id = ctg_id.split(gnm_to_ctg_connector)[0]
        ctg_seq = str(each_ctg.seq)
        ctg_len = len(ctg_seq)
        ctg_gc_num = ctg_seq.upper().count('G') + ctg_seq.upper().count('C')

        # get ctg_seq_dict and ctg_len_dict
        ctg_seq_dict[ctg_id] = ctg_seq
        ctg_len_dict[ctg_id] = ctg_len

        # get mag_to_ctg_dict
        if mag_id not in mag_to_ctg_dict:
            mag_to_ctg_dict[mag_id] = {ctg_id}
        else:
            mag_to_ctg_dict[mag_id].add(ctg_id)

        # get mag_gc_count_dict
        if mag_id not in mag_total_len_and_gc_count_dict:
            mag_total_len_and_gc_count_dict[mag_id] = [ctg_len, ctg_gc_num]
        else:
            mag_total_len_and_gc_count_dict[mag_id][0] += ctg_len
            mag_total_len_and_gc_count_dict[mag_id][1] += ctg_gc_num

    # get mag_gc_content_dict
    mag_gc_content_dict = {}
    for each_mag in mag_total_len_and_gc_count_dict:
        mag_total_len = mag_total_len_and_gc_count_dict[each_mag][0]
        mag_gc_count  = mag_total_len_and_gc_count_dict[each_mag][1]
        mag_gc_content = mag_gc_count*100/mag_total_len
        mag_gc_content = float("{0:.2f}".format(mag_gc_content))
        mag_gc_content_dict[each_mag] = mag_gc_content

    ########################################## get regions (5S, 16S, 23S) to ignore ########################################

    report_and_log(('Get_cn: reading in Barrnap outputs'), pwd_log_file, keep_quiet)

    # get region_to_ignore_dict
    region_to_ignore_dict = {}
    for each_gene in open(combined_prefixed_mags_gff):
        if not each_gene.startswith('#'):
            each_gene_split = each_gene.strip().split('\t')
            ctg_id = each_gene_split[0]
            mag_id = ctg_id.split(gnm_to_ctg_connector)[0]
            hit_start = int(each_gene_split[3])

            # get to_ignore_region_start
            to_ignore_region_start = hit_start - window_len
            if to_ignore_region_start < 1:
                to_ignore_region_start = 1

            # get to_ignore_region_end
            hit_end = int(each_gene_split[4])
            to_ignore_region_end = hit_end

            # initialize dict
            if mag_id not in region_to_ignore_dict:
                region_to_ignore_dict[mag_id] = dict()
            if ctg_id not in region_to_ignore_dict[mag_id]:
                region_to_ignore_dict[mag_id][ctg_id] = set()

            for pos in list(range(to_ignore_region_start, (to_ignore_region_end + 1))):
                region_to_ignore_dict[mag_id][ctg_id].add(pos)

    ################################################## read in depth file ##################################################

    if (mag_depth_gc_txt is None) or (mag_gc_bias_folder is None):

        if os.path.isdir(mag_gc_bias_op_folder) is True:
            os.system('rm -r %s' % mag_gc_bias_op_folder)
        if os.path.isdir(mag_depth_GC_content_folder) is True:
            os.system('rm -r %s' % mag_depth_GC_content_folder)
        os.mkdir(mag_gc_bias_op_folder)
        os.mkdir(mag_depth_GC_content_folder)

        report_and_log(('Get_cn: reading in depth file'), pwd_log_file, keep_quiet)
        seq_pos_depth_dict_gnm_ctg = {}
        for each_bp in open(sam_file_reads_to_mag_filtered_sorted_depth):
            each_bp_split = each_bp.strip().split('\t')
            seq_id = each_bp_split[0]
            gnm_id = seq_id.split(gnm_to_ctg_connector)[0]
            seq_pos = each_bp_split[1]
            pos_depth = int(each_bp_split[2])
            if reads_subsample_pct_for_mag_cov_gc != 100:
                pos_depth = pos_depth * round(100/reads_subsample_pct_for_mag_cov_gc)

            # get seq_pos_depth_dict_gnm_ctg
            if gnm_id not in seq_pos_depth_dict_gnm_ctg:
                seq_pos_depth_dict_gnm_ctg[gnm_id] = dict()
            if seq_id not in seq_pos_depth_dict_gnm_ctg[gnm_id]:
                seq_pos_depth_dict_gnm_ctg[gnm_id][seq_id] = dict()
            seq_pos_depth_dict_gnm_ctg[gnm_id][seq_id][seq_pos] = pos_depth

        # make this multiprocessing
        get_gc_bias_argument_lol = []
        for mag in gnm_to_linked_16s_dict:
            current_mag_depth_GC_content_file   = '%s/%s_depth_GC_content_MAG.txt' % (mag_depth_GC_content_folder, mag)
            current_mag_gc_bias_txt             = '%s/%s_GC_bias.txt' % (mag_gc_bias_op_folder, mag)
            current_mag_gc_bias_png             = '%s/%s_GC_bias.png' % (mag_gc_bias_op_folder, mag)
            current_mag_pos_depth_dict          = seq_pos_depth_dict_gnm_ctg.get(mag, dict())
            current_mag_region_to_ignore_dict   = region_to_ignore_dict.get(mag, dict())
            current_mag_ctg_set                 = mag_to_ctg_dict[mag]

            current_mag_argument_list = [mag, current_mag_ctg_set, current_mag_pos_depth_dict,
                                         ctg_seq_dict, ctg_len_dict, ctg_len_cutoff,
                                         current_mag_region_to_ignore_dict,
                                         ignore_ends_len_mag, window_len, mag_gc_content_dict,
                                         gnm_to_linked_16s_gc_content_dict,
                                         current_mag_gc_bias_txt, current_mag_gc_bias_png, current_mag_depth_GC_content_file]

            get_gc_bias_argument_lol.append(current_mag_argument_list)

        # get GC bias with multi
        report_and_log(('Get_cn: getting Coverage and GC bias for %s Genomes with %s cores' % (len(get_gc_bias_argument_lol), num_threads)), pwd_log_file, keep_quiet)
        pool = mp.Pool(processes=num_threads)
        pool.map(get_mag_gc_bias_worker, get_gc_bias_argument_lol)
        pool.close()
        pool.join()
        report_and_log('Get_cn: get Coverage and GC bias done', pwd_log_file, keep_quiet)

        # combine depth_GC_content_file_mag
        # header: 'MAG\tDepth(x)\tGC(%)\n'
        os.system('cat %s/*.txt > %s' % (mag_depth_GC_content_folder, mag_depth_GC_content_file))
        os.system('rm -r %s' % mag_depth_GC_content_folder)
    else:
        mag_depth_GC_content_file = mag_depth_gc_txt
        mag_gc_bias_op_folder = mag_gc_bias_folder

    ########################################################################################################################
    ############################################# get 16S depth and GC content #############################################
    ########################################################################################################################

    ################################################### map reads to 16S ###################################################

    if sorted_read_to_16s_sam_file is None:

        # mapping 16S reads to all 16S
        report_and_log(('Get_cn: mapping 16S reads to all 16S rRNA gene sequences'), pwd_log_file, keep_quiet)

        index_all_16s_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, marker_seq_file, marker_seq_index)
        report_and_log((index_all_16s_cmd), pwd_log_file, True)
        os.system(index_all_16s_cmd)

        if matam_16s_reads is None:
            map_16s_reads_to_linked_16s_cmd = 'bowtie2 -x %s -U %s,%s -S %s -p %s -f --local --xeq --all --no-unal -N 1 -L 30 &> /dev/null' % (marker_seq_index, reads_file_r1, reads_file_r2, all_16s_sam_file, num_threads)
        else:
            map_16s_reads_to_linked_16s_cmd = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --xeq --all --no-unal -N 1 -L 30 &> /dev/null' % (marker_seq_index, matam_16s_reads, all_16s_sam_file, num_threads)
        report_and_log((map_16s_reads_to_linked_16s_cmd), pwd_log_file, True)
        os.system(map_16s_reads_to_linked_16s_cmd)

        if os.stat(all_16s_sam_file).st_size == 0:
            print('Sam file is empty, program exited!')
            print(all_16s_sam_file)
            exit()

        sort_by_read_cmd = 'samtools sort -n -O sam --threads %s -o %s %s &> /dev/null' % (num_threads, all_16s_sam_file_sorted, all_16s_sam_file)
        report_and_log(('Get_cn: sorting sam file by read id'), pwd_log_file, keep_quiet)
        report_and_log((sort_by_read_cmd), pwd_log_file, True)
        os.system(sort_by_read_cmd)
    else:
        all_16s_sam_file_sorted = sorted_read_to_16s_sam_file

    ############################################# split sorted sam file ############################################

    report_and_log(('Get_cn: splitting sorted sam file '), pwd_log_file, keep_quiet)

    # get the number of lines per file
    os.system('wc -l %s > %s' % (all_16s_sam_file_sorted, all_16s_sam_file_sorted_line_num))
    sam16s_line_num = int(open(all_16s_sam_file_sorted_line_num).readline().strip().split(' ')[0])
    os.remove(all_16s_sam_file_sorted_line_num)
    line_num_per_file = int(round(sam16s_line_num / (num_threads * 10))) + 10

    # split sam file
    os.mkdir(all_16s_sam_file_sorted_split)
    split_sam_cmd = 'split -l %s %s %s/splitted_sam_ ' % (line_num_per_file, all_16s_sam_file_sorted, all_16s_sam_file_sorted_split)
    os.system(split_sam_cmd)

    # get splitted sam file list
    splitted_sam_file_list = [os.path.basename(file_name) for file_name in glob.glob('%s/splitted_sam_*' % all_16s_sam_file_sorted_split)]

    ############################################# parse sorted sam file (all 16S) ############################################

    os.mkdir(all_16s_sam_mp_folder)
    os.mkdir(linked_16s_sam_mp_folder)

    # prepare lol for mp worker
    list_for_parse_sam16s_worker_all = []
    list_for_parse_sam16s_worker_linked = []
    splitted_sam_mp_file_set_all = set()
    splitted_sam_mp_file_set_linked = set()
    for splitted_sam_file in splitted_sam_file_list:
        pwd_splitted_sam_file = '%s/%s' % (all_16s_sam_file_sorted_split, splitted_sam_file)
        pwd_splitted_sam_mp_file_all    = '%s/%s_mp.txt' % (all_16s_sam_mp_folder, splitted_sam_file)
        pwd_splitted_sam_mp_file_linked = '%s/%s_mp.txt' % (linked_16s_sam_mp_folder, splitted_sam_file)
        splitted_sam_mp_file_set_all.add(pwd_splitted_sam_mp_file_all)
        splitted_sam_mp_file_set_linked.add(pwd_splitted_sam_mp_file_linked)
        list_for_parse_sam16s_worker_all.append([pwd_splitted_sam_file,    pwd_splitted_sam_mp_file_all,    min_M_len, mismatch_cutoff, marker_len_dict, []])
        list_for_parse_sam16s_worker_linked.append([pwd_splitted_sam_file, pwd_splitted_sam_mp_file_linked, min_M_len, mismatch_cutoff, marker_len_dict, unlinked_16s_set])

    # parsing mappping results for all 16S
    report_and_log(('Get_cn: parsing alignments of reads to all 16S with multiple threads'), pwd_log_file, keep_quiet)
    pool_parse_sam16s = mp.Pool(processes=num_threads)
    pool_parse_sam16s.map(parse_sam16s_worker, list_for_parse_sam16s_worker_all)
    pool_parse_sam16s.close()
    pool_parse_sam16s.join()

    # parsing mappping results for linked 16S
    report_and_log(('Get_cn: parsing alignments of reads to linked 16S with multiple threads'), pwd_log_file, keep_quiet)
    pool_parse_sam16s = mp.Pool(processes=num_threads)
    pool_parse_sam16s.map(parse_sam16s_worker, list_for_parse_sam16s_worker_linked)
    pool_parse_sam16s.close()
    pool_parse_sam16s.join()

    # report_and_log(('Round 1: removing splitted subsets from disk'), pwd_log_file, keep_quiet)
    os.system('rm -r %s' % all_16s_sam_file_sorted_split)

    ############################################# filter sam file ############################################

    report_and_log(('Get_cn: filtering sam file for all 16S'), pwd_log_file, keep_quiet)
    filter_sam_file(splitted_sam_mp_file_set_all, min_insert_size_16s, both_pair_mapped, keep_uniq,
                    all_16s_sam_file_sorted, all_16s_sam_file_filtered, all_16s_sam_file_filtered_random)

    report_and_log(('Get_cn: filtering sam file for linked 16S'), pwd_log_file, keep_quiet)
    filter_sam_file(splitted_sam_mp_file_set_linked, min_insert_size_16s, both_pair_mapped, keep_uniq,
                    all_16s_sam_file_sorted, linked_16s_sam_file_filtered, linked_16s_sam_file_filtered_random)

    # sort sam file
    report_and_log(('Get_cn: sorting filtered sam file by contig id'), pwd_log_file, keep_quiet)

    sort_by_ctg_cmd        = 'samtools sort -O sam --threads %s -o %s %s &> /dev/null' % (num_threads, all_16s_sam_file_filtered_random_sorted, all_16s_sam_file_filtered_random)
    report_and_log((sort_by_ctg_cmd), pwd_log_file, True)
    os.system(sort_by_ctg_cmd)
    os.remove(all_16s_sam_file_filtered_random)

    sort_by_ctg_cmd_linked = 'samtools sort -O sam --threads %s -o %s %s &> /dev/null' % (num_threads, linked_16s_sam_file_filtered_random_sorted, linked_16s_sam_file_filtered_random)
    report_and_log((sort_by_ctg_cmd_linked), pwd_log_file, True)
    os.system(sort_by_ctg_cmd_linked)
    os.remove(linked_16s_sam_file_filtered_random)

    ############################################# get depth from sam ############################################

    report_and_log(('Get_cn: getting depth directly from sam file (all 16S)'), pwd_log_file, keep_quiet)
    get_depth_from_sam(all_16s_sam_file_filtered_random_sorted, linked_16s_to_gnm_dict, marker_seq_dict, marker_gc_content_dict, ignore_ends_len_16s,
                       pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest,
                       depth_by_bp_file_all, counted_pos_txt_all, depth_GC_content_file_all)

    report_and_log(('Get_cn: getting depth directly from sam file (linked 16S)'), pwd_log_file, keep_quiet)
    get_depth_from_sam(linked_16s_sam_file_filtered_random_sorted, linked_16s_to_gnm_dict, marker_seq_dict, marker_gc_content_dict, ignore_ends_len_16s,
                       pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest,
                       depth_by_bp_file_linked, counted_pos_txt_linked, depth_GC_content_file_linked)

    ########################################################################################################################
    ################################################# get 16S copy number ##################################################
    ########################################################################################################################

    ######################################## read in MAG depth and GC content ########################################

    depth_dict_mag = {}
    gc_content_dict_mag = {}
    for each_mag in open(mag_depth_GC_content_file):
        if not each_mag.startswith('MAG\tDepth(x)\tGC(%)'):
            each_mag_split = each_mag.strip().split('\t')
            mag_id = each_mag_split[0]
            mag_depth = float(each_mag_split[1])
            mag_gc_content = float(each_mag_split[2])
            depth_dict_mag[mag_id] = mag_depth
            gc_content_dict_mag[mag_id] = mag_gc_content

    ######################################## read in 16S depth and GC content ########################################

    depth_dict_16s_all = {}
    gc_content_dict_16s = {}
    for each_16s in open(depth_GC_content_file_all):
        if not each_16s.startswith('16S\tDepth(x)\tGC(%)'):
            each_16s_split = each_16s.strip().split('\t')
            s16_id = each_16s_split[0]
            s16_depth = float(each_16s_split[1])
            s16_gc_content = float(each_16s_split[2])
            depth_dict_16s_all[s16_id] = s16_depth
            gc_content_dict_16s[s16_id] = s16_gc_content

    depth_dict_16s_linked = {}
    for each_16s in open(depth_GC_content_file_linked):
        if not each_16s.startswith('16S\tDepth(x)\tGC(%)'):
            each_16s_split = each_16s.strip().split('\t')
            s16_id = each_16s_split[0]
            s16_depth = float(each_16s_split[1])
            depth_dict_16s_linked[s16_id] = s16_depth

    ######################################## get 16S copy number ########################################

    calculate_cp_num(gnm_to_linked_16s_dict, mag_gc_bias_op_folder, depth_dict_mag, depth_dict_16s_all, gc_content_dict_mag,
                     gc_content_dict_16s, depth_16s_all)

    calculate_cp_num(gnm_to_linked_16s_dict, mag_gc_bias_op_folder, depth_dict_mag, depth_dict_16s_linked, gc_content_dict_mag,
                     gc_content_dict_16s, depth_16s_linked)

    opt_estimation(depth_16s_all, depth_16s_linked, depth_by_bp_file_all, depth_by_bp_file_linked, vxtractor_combined_output,
                   end_len_to_ignore, ignore_region_end_pct, min_region_len_no_ignored,
                   marker_len_dict, estimated_cp_num_txt_by_16s)

    ########################################################################################################################
    ################################################## summarize by genome #################################################
    ########################################################################################################################

    # get reported_ref_16s_copy_num_dict
    reported_ref_16s_copy_num_dict = {}
    if provided_cp_num_txt is not None:
        for each_ref in open(provided_cp_num_txt):
            if not each_ref.startswith('Genome\t'):
                each_ref_split = each_ref.strip().split('\t')
                ref_id = each_ref_split[0]
                reported_copy_num = float(each_ref_split[1])
                reported_ref_16s_copy_num_dict[ref_id] = reported_copy_num

    mag_cp_num_dict = {}
    mag_cp_num_dict_norm = {}
    for each_16s in open(estimated_cp_num_txt_by_16s):
        if not each_16s.startswith('MAG	16S	Coverage'):
            each_16s_split = each_16s.strip().split('\t')
            mag_id = each_16s_split[0]
            cp_num = float(each_16s_split[5])
            cp_num_norm = float(each_16s_split[6])
            if mag_id not in mag_cp_num_dict:
                mag_cp_num_dict[mag_id] = cp_num
                mag_cp_num_dict_norm[mag_id] = cp_num_norm
            else:
                mag_cp_num_dict[mag_id] += cp_num
                mag_cp_num_dict_norm[mag_id] += cp_num_norm

    estimated_16s_cp_num_list = []
    estimated_16s_cp_num_list_norm = []
    reported_16s_cp_num_list = []
    estimation_by_gnm_handle = open(estimated_cp_num_txt_by_mag, 'w')
    if provided_cp_num_txt is None:
        estimation_by_gnm_handle.write('MAG\tCopies\n')
    else:
        estimation_by_gnm_handle.write('MAG\tCopies\tProvided\n')
    for each_mag in mag_cp_num_dict:
        mag_cp_num = mag_cp_num_dict.get(each_mag, 0)
        mag_cp_num_norm = mag_cp_num_dict_norm.get(each_mag, 0)
        mag_cp_num_provided = reported_ref_16s_copy_num_dict.get(each_mag, 0)

        # add to cp num list
        estimated_16s_cp_num_list.append(mag_cp_num)
        estimated_16s_cp_num_list_norm.append(mag_cp_num_norm)
        reported_16s_cp_num_list.append(mag_cp_num_provided)

        # write out estimations
        if provided_cp_num_txt is None:
            estimation_by_gnm_handle.write('%s\t%s\n' % (each_mag, mag_cp_num_norm))
        else:
            estimation_by_gnm_handle.write('%s\t%s\t%s\n' % (each_mag, mag_cp_num_norm, mag_cp_num_provided))

    estimation_by_gnm_handle.close()

    # get scatter plot
    if provided_cp_num_txt is not None:
        get_scatter_plot(estimated_16s_cp_num_list, reported_16s_cp_num_list, scatter_plot)
        get_scatter_plot(estimated_16s_cp_num_list_norm, reported_16s_cp_num_list, scatter_plot_norm)

    # Final report
    report_and_log(('Get_cn: estimated copy number exported to %s_copy_num_by_MAG.txt' % op_prefix), pwd_log_file, keep_quiet)

    report_and_log(('Get_cn: removing tmp files'), pwd_log_file, keep_quiet)
    # os.system('rm %s' % all_16s_sam_file_filtered_random_sorted)
    # os.system('rm %s' % all_16s_sam_file_sorted)
    # os.system('rm %s' % all_16s_sam_file_filtered)
    # os.system('rm %s' % counted_pos_txt_all)
    # os.system('rm %s' % counted_pos_txt_linked)
    # os.system('rm %s' % depth_16s_all)
    # os.system('rm %s' % depth_16s_linked)
    # os.system('rm %s' % depth_GC_content_file_all)
    # os.system('rm %s' % depth_GC_content_file_linked)
    # os.system('rm %s' % depth_by_bp_file_all)
    # os.system('rm %s' % depth_by_bp_file_linked)
    # os.system('rm %s' % linked_16s_sam_file_filtered_random_sorted)
    # os.system('rm %s' % linked_16s_sam_file_filtered)
    # os.system('rm %s' % sam_file_reads_to_mag_filtered_sorted_depth)
    # os.system('rm -r %s' % index_folder)
    # os.system('rm -r %s' % linked_16s_sam_mp_folder)
    # os.system('rm -r %s' % all_16s_sam_mp_folder)

    report_and_log(('Get_cn: Done!'), pwd_log_file, keep_quiet)


    ########################################################################################################################

if __name__ == '__main__':

    get_16s_copy_num_parser = argparse.ArgumentParser()
    get_16s_copy_num_parser.add_argument('-o',                      required=False, metavar='',                 default='.',            help='output folder, default current folder')
    get_16s_copy_num_parser.add_argument('-p',                      required=True,                                                      help='output prefix')
    get_16s_copy_num_parser.add_argument('-r1',                     required=False, metavar='',                                         help='paired reads r1 (fasta format)')
    get_16s_copy_num_parser.add_argument('-r2',                     required=False, metavar='',                                         help='paired reads r2 (fasta format)')
    get_16s_copy_num_parser.add_argument('-r16s',                   required=False, metavar='',                 default=None,           help='matam_16s_reads')
    get_16s_copy_num_parser.add_argument('-linkages',               required=False, metavar='',                                         help='identified_linkage_gnm_level')
    get_16s_copy_num_parser.add_argument('-marker',                 required=False, metavar='',                                         help='markers clustered at 99')
    get_16s_copy_num_parser.add_argument('-mag',                    required=False, metavar='',                 default=None,           help='MAG folder')
    get_16s_copy_num_parser.add_argument('-x',                      required=False, metavar='',                 default=None,           help='MAG extension')
    get_16s_copy_num_parser.add_argument('-cp_mags',                required=False, metavar='',                 default=None,           help='combined_prefixed_mags')
    get_16s_copy_num_parser.add_argument('-cp_mag_gff',             required=False, metavar='',                 default=None,           help='combined_prefixed_mags_gff')
    get_16s_copy_num_parser.add_argument('-mismatch',               required=False, metavar='', type=float,     default=2,              help='mismatch_cutoff')
    get_16s_copy_num_parser.add_argument('-aln_len',                required=False, metavar='', type=int,       default=50,             help='min_M_len')
    get_16s_copy_num_parser.add_argument('-subsample_pct',          required=False, metavar='', type=float,     default=25,             help='used a fraction of reads for MAG coverage estimation (in percentage), default:25')
    get_16s_copy_num_parser.add_argument('-min_insert_size_16s',    required=False, metavar='',                 default=-1000,          help='min_insert_size_16s')
    get_16s_copy_num_parser.add_argument('-ignore_gc_bias',         required=False, action="store_true",                                help='ignore GC bias')
    get_16s_copy_num_parser.add_argument('-ignore_ends_len_16s',    required=False, metavar='', type=int,       default=150,            help='ignore_ends_len_16s')
    get_16s_copy_num_parser.add_argument('-ignore_lowest_pct',      required=False, metavar='', type=float,     default=25,             help='pos_pct_to_ignore_lowest')
    get_16s_copy_num_parser.add_argument('-ignore_highest_pct',     required=False, metavar='', type=float,     default=25,             help='pos_pct_to_ignore_highest')
    get_16s_copy_num_parser.add_argument('-both_pair_mapped',       required=False, action='store_true',                                help='both_pair_mapped')
    get_16s_copy_num_parser.add_argument('-mag_cov_gc',             required=False,                             default=None,           help='mag_depth_gc txt file')
    get_16s_copy_num_parser.add_argument('-mag_gc_bias',            required=False,                             default=None,           help='mag_gc_bias folder')
    get_16s_copy_num_parser.add_argument('-sam_16s_sorted_by_read', required=False,                             default=None,           help='read_to_16s_sam')
    get_16s_copy_num_parser.add_argument('-vxtractor',              required=False,                             default='vxtractor.pl', help='path to vxtractor.pl')
    get_16s_copy_num_parser.add_argument('-ref_16s_cp_num',         required=False,                             default=None,           help='provided_cp_num_txt')
    get_16s_copy_num_parser.add_argument('-t',                      required=False, type=int,                   default=1,              help='number of threads')
    get_16s_copy_num_parser.add_argument('-force',                  required=False, action="store_true",                                help='Force overwrite existing results')
    get_16s_copy_num_parser.add_argument('-quiet',                  required=False, action="store_true",                                help='Not report progress')
    get_16s_copy_num_parser.add_argument('-silva_order_refs',       required=False,                                                     help='subsampled SILVA 16S sequences (order level)')
    get_16s_copy_num_parser.add_argument('-hmm_bac',                required=False,                                                     help='vxtractor hmm_bac folder')
    get_16s_copy_num_parser.add_argument('-hmm_arc',                required=False,                                                     help='vxtractor hmm_arc folder')

    args = vars(get_16s_copy_num_parser.parse_args())
    get_16s_copy_num(args)


'''
check if reads are in pair
'''