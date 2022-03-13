import os
import glob
import shutil
import random
import argparse
import numpy as np
from Bio import SeqIO
import multiprocessing as mp
from itertools import groupby
from operator import itemgetter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
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

    num_arrary_1 = np.array(num_list_1)
    num_arrary_2 = np.array(num_list_2)

    fig = plt.figure(figsize=(6, 6))
    plt.margins(0)

    plt.scatter(num_arrary_1, num_arrary_2)
    plt.xlabel("Estimated copy number", fontsize=12)
    plt.ylabel("User provided copy number", fontsize=12)

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

    title_str = 'y = %sx + %s, R-squared = %s' % (float("{0:.4f}".format(slope)), float("{0:.4f}".format(intercept)), float("{0:.4f}".format(r_squared)))
    plt.title(title_str)

    # save plot
    plt.tight_layout()
    plt.savefig(png_file)
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

    print('Processing %s' % mag)
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
                            window_gc_content = (window_seq_upper.count('G') + window_seq_upper.count('C')) * 100 / len(
                                window_seq)
                            # print('%s-%s\t%sx\t%s%s' % (window_start_pos, window_end_pos, window_mean_depth, window_gc_content, '%'))

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

    sorted_sam          = argument_list[0]
    MappingRecord_file  = argument_list[1]
    min_M_len_16s       = argument_list[2]
    mismatch_cutoff     = argument_list[3]
    marker_len_dict     = argument_list[4]

    MappingRecord_file_handle = open(MappingRecord_file, 'w')
    current_read_base = ''
    current_read_base_r1_16s_ref_dict = dict()
    current_read_base_r2_16s_ref_dict = dict()
    with open(sorted_sam) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if not each_read.startswith('@'):
                each_read_split = each_read.strip().split('\t')
                cigar = each_read_split[5]
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_read_split[2]
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


def get_16s_copy_num(args):

    ############################################### file and parameter in ##############################################

    op_folder                           = args['o']
    op_prefix                           = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    matam_16s_reads                     = args['r16s']
    combined_prefixed_mags              = args['cp_mags']
    combined_prefixed_mags_gff          = args['cp_mag_gff']
    sorted_read_to_16s_sam_file         = args['sam_16s_sorted_by_read']
    identified_linkage_gnm_level        = args['linkages']
    marker_seq_file                     = args['marker']
    marker_seq_file_unclustered         = args['unclustered_marker']
    min_M_len                           = args['aln_len']
    mismatch_cutoff                     = args['mismatch']
    provided_cp_num_txt                 = args['ref']
    both_pair_mapped                    = args['both_pair_mapped']
    ignore_ends_len_16s                 = args['ignore_ends_len_16s']
    pos_pct_to_ignore_lowest            = args['ignore_lowest_pct']
    pos_pct_to_ignore_highest           = args['ignore_highest_pct']
    ignore_gc_bias                      = args['ignore_gc_bias']
    map_to_all_16s                      = args['map_to_all_16s']
    sam_for_16s_depth                   = args['sam_for_16s_depth']
    min_insert_size_16s                 = int(args['min_insert_size_16s'])
    keep_uniq                           = args['by_uniq_aln']
    mag_depth_gc_txt                    = args['mag_cov_gc']
    mag_gc_bias_folder                  = args['mag_gc_bias']
    num_threads                         = args['t']
    force_overwrite                     = args['force']
    keep_quiet                          = args['quiet']
    consider_as_mis_assignment_pct_cutoff = args['mapc']

    window_len                          = 100
    ignore_ends_len_mag                 = 500
    ctg_len_cutoff                      = 10000
    gnm_to_ctg_connector                = '___C___'

    rm_mis_alignments = False
    if marker_seq_file_unclustered is not None:
        rm_mis_alignments = True

    ################################################### check output folder ###################################################

    op_folder_to_generate = '%s_get_16S_cp_num_wd' % op_prefix

    if os.path.isdir(op_folder) == '.':
        pwd_op_folder_to_generate = op_folder_to_generate
    else:
        pwd_op_folder_to_generate = '%s/%s' % (op_folder, op_folder_to_generate)

    if (os.path.isdir(pwd_op_folder_to_generate) is True) and (force_overwrite is False):
        print('Output folder detected, please specify a different location or provide a different prefix')
        print('Program exited!')
        exit()
    else:
        force_create_folder(pwd_op_folder_to_generate)

    ################################################### define file name ###################################################

    str_in_op = 'all'
    if both_pair_mapped is True:
        str_in_op = 'both'

    file_path_mag, combined_prefixed_mags_no_ext, ext_mag = sep_path_basename_ext(combined_prefixed_mags)
    marker_seq_path, marker_seq_no_ext, marker_seq_ext = sep_path_basename_ext(marker_seq_file)
    marker_seq_unclustered_path, marker_seq_unclustered_no_ext, marker_seq_unclustered_ext = sep_path_basename_ext(marker_seq_file_unclustered)

    index_folder                                = '%s/%s_index'                                         % (pwd_op_folder_to_generate, op_prefix)

    # get_well_assigned_reads_pct
    get_well_assigned_reads_pct_folder          = '%s/%s_get_well_assigned_reads_pct'                   % (pwd_op_folder_to_generate, op_prefix)
    marker_seq_index                            = '%s/%s'                                               % (index_folder, marker_seq_no_ext)
    marker_seq_unclustered_index                = '%s/%s'                                               % (index_folder, marker_seq_unclustered_no_ext)
    r16s_r1                                     = '%s/r16s_R1.fasta'                                    % get_well_assigned_reads_pct_folder
    r16s_r2                                     = '%s/r16s_R2.fasta'                                    % get_well_assigned_reads_pct_folder
    r16s_singleton                              = '%s/r16s_UP.fasta'                                    % get_well_assigned_reads_pct_folder
    clustered_16s_sam                           = '%s/clustered_16s.sam'                                % get_well_assigned_reads_pct_folder
    clustered_16s_sam_log                       = '%s/clustered_16s.log'                                % get_well_assigned_reads_pct_folder
    clustered_16s_sam_filtered                  = '%s/clustered_16s_filtered.sam'                       % get_well_assigned_reads_pct_folder
    clustered_16s_sam_filtered_both_mapped      = '%s/clustered_16s_filtered_both_mapped.sam'           % get_well_assigned_reads_pct_folder
    unclustered_16s_sam                         = '%s/unclustered_16s.sam'                              % get_well_assigned_reads_pct_folder
    unclustered_16s_sam_log                     = '%s/unclustered_16s.log'                              % get_well_assigned_reads_pct_folder
    unclustered_16s_sam_filtered                = '%s/unclustered_16s_filtered.sam'                     % get_well_assigned_reads_pct_folder
    unclustered_16s_sam_filtered_both_mapped    = '%s/unclustered_16s_filtered_both_mapped.sam'         % get_well_assigned_reads_pct_folder
    well_assigned_reads_pct_txt                 = '%s/well_assigned_reads_pct.txt'                      % get_well_assigned_reads_pct_folder

    pwd_log_file                                = '%s/%s.log'                                                                                   % (pwd_op_folder_to_generate, op_prefix)
    sam_file_reads_to_mag                       = '%s/%s_%s.sam'                                                                                % (pwd_op_folder_to_generate, op_prefix, combined_prefixed_mags_no_ext)
    sam_file_reads_to_mag_filtered              = '%s/%s_%s_%s_%s.sam'                                                                          % (pwd_op_folder_to_generate, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    sam_file_reads_to_mag_filtered_both_mapped  = '%s/%s_%s_%s_%s_both_mapped.sam'                                                              % (pwd_op_folder_to_generate, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    sam_file_reads_to_mag_filtered_sorted       = '%s/%s_%s_%s_%s_both_mapped_sorted.sam'                                                       % (pwd_op_folder_to_generate, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    sam_file_reads_to_mag_filtered_sorted_depth = '%s/%s_%s_%s_%s_both_mapped_sorted_depth.txt'                                                 % (pwd_op_folder_to_generate, op_prefix, combined_prefixed_mags_no_ext, mismatch_cutoff, min_M_len)
    mag_gc_bias_op_folder                       = '%s/%s_MAGs_%s_%s_GC_bias'                                                                    % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len)
    mag_depth_GC_content_folder                 = '%s/%s_MAGs_%s_%s_MAG_depth_GC_content'                                                       % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len)
    mag_depth_GC_content_file                   = '%s/%s_MAGs_%s_%s_MAG_depth_GC_content.txt'                                                   % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len)
    linked_16s_seq_file                         = '%s/%s_linked_16s.fa'                                                                         % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_seq_file_no_ext                  = '%s/%s_linked_16s'                                                                            % (index_folder, op_prefix)
    linked_16s_sam_file                         = '%s/%s_linked_16s.sam'                                                                        % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_sam_file_log                     = '%s/%s_linked_16s.log'                                                                        % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_sam_file_sorted                  = '%s/%s_linked_16s_sorted_by_read.sam'                                                         % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_sam_file_sorted_log              = '%s/%s_linked_16s_sorted_by_read.log'                                                         % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_sam_file_sorted_line_num         = '%s/%s_linked_16s_sorted_by_read_line_num.txt'                                                % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_sam_file_sorted_split            = '%s/%s_linked_16s_sorted_by_read_split'                                                       % (pwd_op_folder_to_generate, op_prefix)
    linked_16s_sam_mp_folder                    = '%s/%s_linked_16S_%s_%s_mr'                                                                   % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len)
    linked_16s_sam_file_filtered                    = '%s/%s_linked_16S_%s_%s_%s.sam'                                                           % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, str_in_op)
    linked_16s_sam_file_filtered_random             = '%s/%s_linked_16S_%s_%s_%s_random.sam'                                                    % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, str_in_op)
    linked_16s_sam_file_filtered_random_sorted      = '%s/%s_linked_16S_%s_%s_%s_random_sorted.sam'                                             % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, str_in_op)
    linked_16s_sam_file_filtered_random_sorted_log  = '%s/%s_linked_16S_%s_%s_%s_random_sorted.log'                                             % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, str_in_op)
    linked_16s_depth_GC_content_file            = '%s/%s_linked_16S_%s_%s_%s_depth_GC_content.txt'                                              % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, str_in_op)
    assessment_txt_pos                          = '%s/%s_MAGs_%s_%s__16S_%s_%s__ignore_ed%s_ld%s_hd%s__%s_counted_pos.txt'                      % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, mismatch_cutoff, min_M_len, ignore_ends_len_16s, pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest, str_in_op)
    assessment_txt_depth                        = '%s/%s_MAGs_%s_%s__16S_%s_%s__ignore_ed%s_ld%s_hd%s__%s_depth.txt'                            % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, mismatch_cutoff, min_M_len, ignore_ends_len_16s, pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest, str_in_op)
    assessment_txt                              = '%s/%s_MAGs_%s_%s__16S_%s_%s__ignore_ed%s_ld%s_hd%s__%s_estimation.txt'                       % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, mismatch_cutoff, min_M_len, ignore_ends_len_16s, pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest, str_in_op)
    scatter_plot                                = '%s/%s_MAGs_%s_%s__16S_%s_%s__ignore_ed%s_ld%s_hd%s__%s_assessment.png'                       % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, mismatch_cutoff, min_M_len, ignore_ends_len_16s, pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest, str_in_op)
    scatter_plot_normalized_by_gc_bias          = '%s/%s_MAGs_%s_%s__16S_%s_%s__ignore_ed%s_ld%s_hd%s__%s_assessment_normalized_by_gc_bias.png' % (pwd_op_folder_to_generate, op_prefix, mismatch_cutoff, min_M_len, mismatch_cutoff, min_M_len, ignore_ends_len_16s, pos_pct_to_ignore_lowest, pos_pct_to_ignore_highest, str_in_op)

    ########################################################################################################################
    ################################################### estimate GC bias ###################################################
    ########################################################################################################################

    ########################################## map reads to combined prefixed MAGs #########################################

    if mag_depth_gc_txt is None:
        report_and_log(('mapping reads to combined prefixed MAGs'), pwd_log_file, keep_quiet)

        bowtie_build_cmd = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, combined_prefixed_mags, index_folder, combined_prefixed_mags_no_ext)
        report_and_log((bowtie_build_cmd), pwd_log_file, True)
        os.system(bowtie_build_cmd)

        # mapping_cmd_mag = 'bowtie2 -x %s/%s -U %s,%s -S %s -p %s -f --xeq --no-unal -N 1 -L 30' % (pwd_op_folder_to_generate, combined_prefixed_mags_no_ext, reads_file_r1, reads_file_r2, sam_file_reads_to_mag, num_threads)
        mapping_cmd_mag = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s -f --xeq --no-unal' % (index_folder, combined_prefixed_mags_no_ext, reads_file_r1, reads_file_r2, sam_file_reads_to_mag, num_threads)
        report_and_log((mapping_cmd_mag), pwd_log_file, True)
        os.system(mapping_cmd_mag)

        # filter sam file
        report_and_log(('removing high mismatch alignments'), pwd_log_file, keep_quiet)
        remove_high_mismatch(sam_file_reads_to_mag, min_M_len, mismatch_cutoff, sam_file_reads_to_mag_filtered)
        os.remove(sam_file_reads_to_mag)

        #keep_only_both_mapped_reads
        report_and_log(('removing singletons'), pwd_log_file, keep_quiet)
        keep_only_both_mapped_reads(sam_file_reads_to_mag_filtered, sam_file_reads_to_mag_filtered_both_mapped)
        os.remove(sam_file_reads_to_mag_filtered)

        # sort filtered_sam
        report_and_log(('sorting filtered sam file'), pwd_log_file, keep_quiet)
        sort_cmd_mag = 'samtools sort -O sam --threads %s -o %s %s' % (num_threads, sam_file_reads_to_mag_filtered_sorted, sam_file_reads_to_mag_filtered_both_mapped)
        report_and_log((sort_cmd_mag), pwd_log_file, True)
        os.system(sort_cmd_mag)
        os.remove(sam_file_reads_to_mag_filtered_both_mapped)

        # get depth file
        report_and_log(('getting depth with samtools'), pwd_log_file, keep_quiet)
        get_depth_cmd_mag = 'samtools depth -a %s > %s' % (sam_file_reads_to_mag_filtered_sorted, sam_file_reads_to_mag_filtered_sorted_depth)
        report_and_log((get_depth_cmd_mag), pwd_log_file, True)
        os.system(get_depth_cmd_mag)
        os.remove(sam_file_reads_to_mag_filtered_sorted)

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

    ############################################ read in MAG sequences into dict ###########################################

    report_and_log(('reading in MAG sequences'), pwd_log_file, keep_quiet)

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

    report_and_log(('reading in Barrnap outputs'), pwd_log_file, keep_quiet)

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


        report_and_log(('reading in depth file'), pwd_log_file, keep_quiet)
        seq_pos_depth_dict_gnm_ctg = {}
        for each_bp in open(sam_file_reads_to_mag_filtered_sorted_depth):
            each_bp_split = each_bp.strip().split('\t')
            seq_id = each_bp_split[0]
            gnm_id = seq_id.split(gnm_to_ctg_connector)[0]
            seq_pos = each_bp_split[1]
            pos_depth = int(each_bp_split[2])

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
        report_and_log(('getting Coverage and GC bias for %s Genomes with %s cores' % (len(get_gc_bias_argument_lol), num_threads)), pwd_log_file, keep_quiet)
        pool = mp.Pool(processes=num_threads)
        pool.map(get_mag_gc_bias_worker, get_gc_bias_argument_lol)
        pool.close()
        pool.join()

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

    if map_to_all_16s is False:
        report_and_log(('getting sequences of linked 16S'), pwd_log_file, keep_quiet)
        linked_16s_seq_file_handle = open(linked_16s_seq_file, 'w')
        for each_16s in marker_seq_dict:
            if each_16s in linked_16s_to_gnm_dict:
                linked_16s_seq_file_handle.write('>%s\n' % each_16s)
                linked_16s_seq_file_handle.write('%s\n' % marker_seq_dict[each_16s])
        linked_16s_seq_file_handle.close()

    os.mkdir(index_folder)

    ############################################# get_well_assigned_reads_pct ##############################################

    if map_to_all_16s is False:
        index_linked_16s_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, linked_16s_seq_file, linked_16s_seq_file_no_ext)
        os.system(index_linked_16s_cmd)

    index_clustered_16s_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, marker_seq_file, marker_seq_index)
    os.system(index_clustered_16s_cmd)

    if rm_mis_alignments is True:

        os.mkdir(get_well_assigned_reads_pct_folder)

        report_and_log(('Assessing mapping accuracy'), pwd_log_file, keep_quiet)

        # sep_paired_and_singleton_reads
        sep_paired_and_singleton_reads(matam_16s_reads, r16s_r1, r16s_r2, r16s_singleton)

        # index references
        index_unclustered_16s_cmd   = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, marker_seq_file_unclustered, marker_seq_unclustered_index)
        os.system(index_unclustered_16s_cmd)

        # mapping
        report_and_log(('Mapping reads to clustered 16S rRNA gene sequences'), pwd_log_file, keep_quiet)
        if map_to_all_16s is False:
            map_to_clustered_16s_cmd    = 'bowtie2 -x %s -1 %s -2 %s -S %s -p %s -f --local --xeq --no-unal 2> %s' % (linked_16s_seq_file_no_ext, r16s_r1, r16s_r2, clustered_16s_sam, num_threads, clustered_16s_sam_log)
            os.system(map_to_clustered_16s_cmd)
        else:
            map_to_clustered_16s_cmd    = 'bowtie2 -x %s -1 %s -2 %s -S %s -p %s -f --local --xeq --no-unal 2> %s' % (marker_seq_index, r16s_r1, r16s_r2, clustered_16s_sam, num_threads, clustered_16s_sam_log)
            os.system(map_to_clustered_16s_cmd)

        report_and_log(('Mapping reads to unclustered 16S rRNA gene sequences'), pwd_log_file, keep_quiet)
        map_to_unclustered_16s_cmd  = 'bowtie2 -x %s -1 %s -2 %s -S %s -p %s -f --local --xeq --no-unal 2> %s' % (marker_seq_unclustered_index, r16s_r1, r16s_r2, unclustered_16s_sam, num_threads, unclustered_16s_sam_log)
        os.system(map_to_unclustered_16s_cmd)

        # remove_high_mismatch
        report_and_log(('Filtering mapping results'), pwd_log_file, keep_quiet)
        remove_high_mismatch(clustered_16s_sam, min_M_len, mismatch_cutoff, clustered_16s_sam_filtered)
        remove_high_mismatch(unclustered_16s_sam, min_M_len, mismatch_cutoff, unclustered_16s_sam_filtered)
        os.remove(clustered_16s_sam)
        os.remove(unclustered_16s_sam)

        #keep_only_both_mapped_reads
        keep_only_both_mapped_reads(clustered_16s_sam_filtered, clustered_16s_sam_filtered_both_mapped)
        keep_only_both_mapped_reads(unclustered_16s_sam_filtered, unclustered_16s_sam_filtered_both_mapped)
        os.remove(clustered_16s_sam_filtered)
        os.remove(unclustered_16s_sam_filtered)

        # get well_assigned_reads_pct_dict
        well_assigned_reads_pct_dict = get_well_aligned_reads_pct_dict(clustered_16s_sam_filtered_both_mapped, unclustered_16s_sam_filtered_both_mapped, consider_as_mis_assignment_pct_cutoff)

        # write out assessment
        well_assigned_reads_pct_txt_handle = open(well_assigned_reads_pct_txt, 'w')
        well_assigned_reads_pct_txt_handle.write('Marker\tAccuracy\n')
        for each_marker in well_assigned_reads_pct_dict:
            well_assigned_reads_pct_txt_handle.write('%s\t%s\n' % (each_marker, well_assigned_reads_pct_dict[each_marker]))
        well_assigned_reads_pct_txt_handle.close()

        report_and_log(('Assess mapping accuracy done'), pwd_log_file, keep_quiet)

    ########################################################################################################################

    if sam_for_16s_depth is None:
        if sorted_read_to_16s_sam_file is None:

            # pretend the overall one are the linked sequences
            if map_to_all_16s is True:
                linked_16s_seq_file = marker_seq_file
                linked_16s_seq_file_no_ext = marker_seq_index

            # mapping 16S reads to linked 16S
            if map_to_all_16s is False:
                report_and_log(('mapping 16S reads to linked 16S rRNA gene sequences'), pwd_log_file, keep_quiet)
            else:
                report_and_log(('mapping 16S reads to clustered 16S rRNA gene sequences'), pwd_log_file, keep_quiet)

            index_16s_cmd = 'bowtie2-build --quiet --threads 16 -f %s %s' % (linked_16s_seq_file, linked_16s_seq_file_no_ext)
            report_and_log((index_16s_cmd), pwd_log_file, True)
            os.system(index_16s_cmd)

            map_16s_reads_to_linked_16s_cmd = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --xeq --all --no-unal -N 1 -L 30 2> %s' % (linked_16s_seq_file_no_ext, matam_16s_reads, linked_16s_sam_file, num_threads, linked_16s_sam_file_log)
            report_and_log((map_16s_reads_to_linked_16s_cmd), pwd_log_file, True)
            os.system(map_16s_reads_to_linked_16s_cmd)

            sort_by_read_cmd = 'samtools sort -n -O sam --threads %s -o %s %s 2> %s' % (num_threads, linked_16s_sam_file_sorted, linked_16s_sam_file, linked_16s_sam_file_sorted_log)
            report_and_log(('sorting sam file by read id'), pwd_log_file, keep_quiet)
            report_and_log((sort_by_read_cmd), pwd_log_file, True)
            os.system(sort_by_read_cmd)
        else:
            linked_16s_sam_file_sorted = sorted_read_to_16s_sam_file

        # parse sorted sam file with multiple threads
        report_and_log(('parsing sorted sam file with multiple threads'), pwd_log_file, keep_quiet)

        os.mkdir(linked_16s_sam_mp_folder)
        splitted_sam_mp_file_set = parse_sam_file_with_multiprocessing(linked_16s_sam_file_sorted,
                                                                       linked_16s_sam_file_sorted_line_num,
                                                                       linked_16s_sam_file_sorted_split,
                                                                       linked_16s_sam_mp_folder,
                                                                       min_M_len, mismatch_cutoff,
                                                                       marker_len_dict, num_threads)
        report_and_log(('parse sam done'), pwd_log_file, keep_quiet)

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
                    current_read_base___r1_16s_refs_no_ignored = get_no_ignored_dict_from_str(each_read_base_split[7])
                    current_read_base___r2_16s_refs_no_ignored = get_no_ignored_dict_from_str(each_read_base_split[8])
                    current_read_base___shared_16s_refs_no_ignored = get_no_ignored_dict_from_str(each_read_base_split[9])

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

                        MappingRecord_dict[current_read_base___id].r1_16s_ref_dict = current_read_base___r1_16s_ref_dict
                        MappingRecord_dict[current_read_base___id].r2_16s_ref_dict = current_read_base___r2_16s_ref_dict
                        MappingRecord_dict[current_read_base___id].r1_16s_refs_no_ignored = current_read_base___r1_16s_refs_no_ignored
                        MappingRecord_dict[current_read_base___id].r2_16s_refs_no_ignored = current_read_base___r2_16s_refs_no_ignored

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
                                    current_read_base___shared_16s_refs_no_ignored_with_large_insert_size[each_shared_ref] = current_read_base___shared_16s_refs_no_ignored[each_shared_ref]
                                    large_insert_size_paired_read_num += 1
                                else:
                                    small_insert_size_paired_read_num += 1

                        MappingRecord_dict[current_read_base___id].shared_16s_refs_no_ignored = current_read_base___shared_16s_refs_no_ignored_with_large_insert_size

        # print('large_insert_size_paired_read_num\t%s' % large_insert_size_paired_read_num)
        # print('small_insert_size_paired_read_num\t%s' % small_insert_size_paired_read_num)

        # filter sam file
        pwd_sam_file_filtered_handle = open(linked_16s_sam_file_filtered, 'w')
        for each_line in open(linked_16s_sam_file_sorted):
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
        for each_line in open(linked_16s_sam_file_filtered):
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
        pwd_sam_file_filtered_random_handle = open(linked_16s_sam_file_filtered_random, 'w')
        for each_line in open(linked_16s_sam_file_filtered):
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

        # sort sam file
        sort_by_ctg_cmd = 'samtools sort -O sam --threads %s -o %s %s 2> %s' % (num_threads, linked_16s_sam_file_filtered_random_sorted, linked_16s_sam_file_filtered_random, linked_16s_sam_file_filtered_random_sorted_log)
        report_and_log(('sorting sam file by contig id'), pwd_log_file, keep_quiet)
        report_and_log((sort_by_ctg_cmd), pwd_log_file, True)
        os.system(sort_by_ctg_cmd)
        os.remove(linked_16s_sam_file_filtered_random)

        # sam file to use
        sam_to_use = linked_16s_sam_file_filtered_random_sorted
    else:
        sam_to_use = sam_for_16s_depth


    report_and_log(('getting depth directly from sam file'), pwd_log_file, keep_quiet)
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
                    current_marker_pos_dict = {}
                    for each_pos in pos_list:
                        current_marker_pos_dict[str(each_pos)] = 0
                    marker_pos_depth_dict[marker_id] = current_marker_pos_dict
            else:
                cigar = each_read_split[5]
                ref_id = each_read_split[2]
                if ref_id in linked_16s_to_gnm_dict:
                    ref_pos = int(each_read_split[3])
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(cigar))
                    ref_pos_end = ref_pos + aligned_len
                    current_read_covered_region = list(range(ref_pos, ref_pos_end))
                    for each_covered_pos in current_read_covered_region:
                        covered_pos_str = str(each_covered_pos)
                        marker_pos_depth_dict[ref_id][covered_pos_str] += 1

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

        value_num_to_ignore_lowest = round(len(non_end_pos_depth_list)*(pos_pct_to_ignore_lowest/100))
        value_num_to_ignore_highest = round(len(non_end_pos_depth_list)*(pos_pct_to_ignore_highest/100))

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
    report_and_log(('writing out 16S rRNA gene depth and GC content'), pwd_log_file, keep_quiet)
    assessment_txt_pos_handle = open(assessment_txt_pos, 'w')
    assessment_txt_pos_handle.write('16S\tCounted_depth_range(X)\tIgnored_low_depth_bps\tIgnored_high_depth_bps\tCounted_bps/Total_length\n')
    depth_GC_content_file_16s_handle = open(linked_16s_depth_GC_content_file, 'w')
    depth_GC_content_file_16s_handle.write('16S\tDepth(x)\tGC(%)\n')
    for each_linked_16s in linked_16s_to_gnm_dict:
        linked_16s_seq = marker_seq_dict[each_linked_16s]
        linked_16s_len = len(linked_16s_seq)
        linked_16s_pos_depth_dict = marker_pos_depth_dict[each_linked_16s]
        current_16s_total_mapped_read_len = 0
        current_16s_total_counted_ctg_len = 0
        currrent_16s_depth_num_to_ignore_lowest     = outlier_depth_dict_16s[each_linked_16s][0]
        currrent_16s_depth_num_to_ignore_highest    = outlier_depth_dict_16s[each_linked_16s][1]

        currrent_16s_low_depth_cutoff = outlier_depth_dict_16s[each_linked_16s][2]
        currrent_16s_high_depth_cutoff = outlier_depth_dict_16s[each_linked_16s][3]

        ignored_low_depth_pos_num = 0
        ignored_high_depth_pos_num = 0
        pos_index = ignore_ends_len_16s + 1
        current_counted_pos_list = []
        while pos_index <= (linked_16s_len - ignore_ends_len_16s):
            current_pos_depth = linked_16s_pos_depth_dict[str(pos_index)]

            if ((current_pos_depth <= currrent_16s_low_depth_cutoff) and (ignored_low_depth_pos_num < currrent_16s_depth_num_to_ignore_lowest)) or ( (current_pos_depth >= currrent_16s_high_depth_cutoff) and (ignored_high_depth_pos_num < currrent_16s_depth_num_to_ignore_highest)):
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
        depth_GC_content_file_16s_handle.write('%s\t%s\t%s\n' % (each_linked_16s, current_16s_average_coverage, current_16s_gc_content))
    depth_GC_content_file_16s_handle.close()
    assessment_txt_pos_handle.close()

    ########################################################################################################################
    ################################################# get 16S copy number ##################################################
    ########################################################################################################################

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

    depth_dict_16s = {}
    gc_content_dict_16s = {}
    for each_16s in open(linked_16s_depth_GC_content_file):
        if not each_16s.startswith('16S\tDepth(x)\tGC(%)'):
            each_16s_split = each_16s.strip().split('\t')
            s16_id = each_16s_split[0]
            s16_depth = float(each_16s_split[1])
            s16_gc_content = float(each_16s_split[2])
            depth_dict_16s[s16_id] = s16_depth
            gc_content_dict_16s[s16_id] = s16_gc_content

    mag_16s_cp_num_dict = {}
    mag_16s_cp_num_dict_normalized_by_gc_bias = {}
    assessment_txt_depth_handle = open(assessment_txt_depth, 'w')
    assessment_txt_depth_handle.write('MAG\t16S\tMAG_depth\t16S_depth\t16S_depth_norm\n')
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
                depth_ratio = mean_gc_depth_16s/mean_gc_depth_mag
                depth_16s_normalized_by_gc_bias = depth_16s/depth_ratio
                depth_16s_normalized_by_gc_bias = float("{0:.2f}".format(depth_16s_normalized_by_gc_bias))
                assessment_txt_depth_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_mag, each_linked_16s, float("{0:.2f}".format(depth_mag)), depth_16s, depth_16s_normalized_by_gc_bias))


                estimate_cp_num                                 = float("{0:.2f}".format(depth_16s/depth_mag))
                estimate_cp_num_normalized_by_gc_bias           = float("{0:.2f}".format(depth_16s_normalized_by_gc_bias/depth_mag))
                if rm_mis_alignments is True:
                    current_16s_well_assigned_reads_pct         = well_assigned_reads_pct_dict.get(each_linked_16s, 100)
                    estimate_cp_num                             = float("{0:.2f}".format(depth_16s * (current_16s_well_assigned_reads_pct/100) / depth_mag))
                    estimate_cp_num_normalized_by_gc_bias       = float("{0:.2f}".format(depth_16s_normalized_by_gc_bias * (current_16s_well_assigned_reads_pct/100) / depth_mag))

                if each_mag not in mag_16s_cp_num_dict:
                    mag_16s_cp_num_dict[each_mag] = estimate_cp_num
                else:
                    mag_16s_cp_num_dict[each_mag] += estimate_cp_num

                if each_mag not in mag_16s_cp_num_dict_normalized_by_gc_bias:
                    mag_16s_cp_num_dict_normalized_by_gc_bias[each_mag] = estimate_cp_num_normalized_by_gc_bias
                else:
                    mag_16s_cp_num_dict_normalized_by_gc_bias[each_mag] += estimate_cp_num_normalized_by_gc_bias
        else:
            print('Ignored %s, as %s not found, ' % (each_mag, gc_bias_txt))

    assessment_txt_depth_handle.close()

    ########################################################################################################################
    ################################################### assess estimation ##################################################
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

    # estimated_16s_cp_num_list = []
    estimated_16s_cp_num_list = []
    estimated_16s_cp_num_list_normalized_by_gc_bias = []
    reported_16s_cp_num_list = []
    assessment_txt_handle = open(assessment_txt, 'w')
    if provided_cp_num_txt is not None:
        assessment_txt_handle.write('MAG\tCopy_number_raw\tCopy_number_normalized\tUser_provided\t\n')
    else:
        assessment_txt_handle.write('MAG\tCopy_number_raw\tCopy_number_normalized\n')
    for each_mag in mag_16s_cp_num_dict_normalized_by_gc_bias:
        estimate_cp_num = float("{0:.2f}".format(mag_16s_cp_num_dict.get(each_mag, 0)))
        estimate_cp_num_normalized_by_gc_bias = float("{0:.2f}".format(mag_16s_cp_num_dict_normalized_by_gc_bias.get(each_mag, 0)))
        if provided_cp_num_txt is not None:
            assessment_txt_handle.write('%s\t%s\t%s\t%s\n' % (each_mag, estimate_cp_num, estimate_cp_num_normalized_by_gc_bias, reported_ref_16s_copy_num_dict.get(each_mag, 0)))
            estimated_16s_cp_num_list.append(mag_16s_cp_num_dict.get(each_mag, 0))
            estimated_16s_cp_num_list_normalized_by_gc_bias.append(mag_16s_cp_num_dict_normalized_by_gc_bias.get(each_mag, 0))
            reported_16s_cp_num_list.append(reported_ref_16s_copy_num_dict.get(each_mag, 0))
        else:
            assessment_txt_handle.write('%s\t%s\t%s\n' % (each_mag, estimate_cp_num, estimate_cp_num_normalized_by_gc_bias))
    assessment_txt_handle.close()

    # get scatter plot
    if provided_cp_num_txt is not None:
        get_scatter_plot(estimated_16s_cp_num_list, reported_16s_cp_num_list, scatter_plot)
        get_scatter_plot(estimated_16s_cp_num_list_normalized_by_gc_bias, reported_16s_cp_num_list, scatter_plot_normalized_by_gc_bias)

    ########################################################################################################################

if __name__ == '__main__':

    get_16s_copy_num_parser = argparse.ArgumentParser()
    get_16s_copy_num_parser.add_argument('-o',                      required=False,  metavar='', default='.',               help='output folder, default current folder')
    get_16s_copy_num_parser.add_argument('-p',                      required=True,                                          help='output prefix')
    get_16s_copy_num_parser.add_argument('-r1',                     required=False,  metavar='',                            help='paired reads r1 (fasta format)')
    get_16s_copy_num_parser.add_argument('-r2',                     required=False,  metavar='',                            help='paired reads r2 (fasta format)')
    get_16s_copy_num_parser.add_argument('-r16s',                   required=False,                                         help='matam_16s_reads')
    get_16s_copy_num_parser.add_argument('-linkages',               required=False,                                         help='identified_linkage_gnm_level')
    get_16s_copy_num_parser.add_argument('-marker',                 required=False,                                         help='markers clustered at 99')
    get_16s_copy_num_parser.add_argument('-unclustered_marker',     required=False, default=None,                           help='unclustered_marker, or markers clustered at 99.9')
    get_16s_copy_num_parser.add_argument('-cp_mags',                required=False,                                         help='combined_prefixed_mags')
    get_16s_copy_num_parser.add_argument('-cp_mag_gff',             required=False,                                         help='combined_prefixed_mags_gff')
    get_16s_copy_num_parser.add_argument('-mismatch',               required=False, metavar='', type=float, default=2,      help='mismatch_cutoff')
    get_16s_copy_num_parser.add_argument('-aln_len',                required=False, metavar='', type=int,   default=50,     help='min_M_len')
    get_16s_copy_num_parser.add_argument('-min_insert_size_16s',    required=False, metavar='',             default=-1000,  help='min_insert_size_16s')
    get_16s_copy_num_parser.add_argument('-ignore_gc_bias',         required=False, metavar='', type=int,   default=150,    help='ignore GC bias')
    get_16s_copy_num_parser.add_argument('-ignore_ends_len_16s',    required=False, metavar='', type=int,   default=150,    help='ignore_ends_len_16s')
    get_16s_copy_num_parser.add_argument('-ignore_lowest_pct',      required=False, metavar='', type=float, default=45,     help='pos_pct_to_ignore_lowest')
    get_16s_copy_num_parser.add_argument('-ignore_highest_pct',     required=False, metavar='', type=float, default=5,      help='pos_pct_to_ignore_highest')
    get_16s_copy_num_parser.add_argument('-mapc',                   required=False, metavar='', type=float, default=15,     help='mapping accuracy percentage cutoff')
    get_16s_copy_num_parser.add_argument('-by_uniq_aln',            required=False, action='store_true',                    help='use uniquely aligned reads for 16S depth estimation')
    get_16s_copy_num_parser.add_argument('-both_pair_mapped',       required=False, action='store_true',                    help='both_pair_mapped')
    get_16s_copy_num_parser.add_argument('-map_to_all_16s',         required=False, action='store_true',                    help='map_to_all_16s, not the linked')
    get_16s_copy_num_parser.add_argument('-mag_cov_gc',             required=False, default=None,                           help='mag_depth_gc txt file')
    get_16s_copy_num_parser.add_argument('-mag_gc_bias',            required=False, default=None,                           help='mag_gc_bias folder')
    get_16s_copy_num_parser.add_argument('-sam_16s_sorted_by_read', required=False, default=None,                           help='read_to_16s_sam')
    get_16s_copy_num_parser.add_argument('-sam_for_16s_depth',      required=False, default=None,                           help='sam_for_16s_depth')
    get_16s_copy_num_parser.add_argument('-ref',                    required=False, default=None,                           help='provided_cp_num_txt')
    get_16s_copy_num_parser.add_argument('-t',                      required=False, type=int, default=1,                    help='number of threads')
    get_16s_copy_num_parser.add_argument('-force',                  required=False, action="store_true",                    help='Force overwrite existing results')
    get_16s_copy_num_parser.add_argument('-quiet',                  required=False, action="store_true",                    help='Not report progress')

    args = vars(get_16s_copy_num_parser.parse_args())
    get_16s_copy_num(args)
