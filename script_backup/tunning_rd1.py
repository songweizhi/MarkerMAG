import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
from datetime import datetime
from Bio.SeqRecord import SeqRecord


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

    aligned_pct = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
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

        self.r1_16s_refs_no_ignored_with_pos = dict()
        self.r2_16s_refs_no_ignored_with_pos = dict()
        self.shared_16s_refs_no_ignored_with_pos = dict()

        self.both_mapped_to_16s = False

        #################### round 1 ctg ####################

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()

        self.r1_ctg_refs_lowest_mismatch = None
        self.r2_ctg_refs_lowest_mismatch = None

        self.r1_ctg_refs_no_ignored = dict()
        self.r2_ctg_refs_no_ignored = dict()
        self.shared_ctg_refs_no_ignored = dict()

        self.r1_ctg_refs_no_ignored_with_pos = dict()
        self.r2_ctg_refs_no_ignored_with_pos = dict()
        self.shared_ctg_refs_no_ignored_with_pos = dict()

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

        self.r1_mini_refs_no_ignored_with_pos = dict()
        self.r2_mini_refs_no_ignored_with_pos = dict()
        self.shared_mini_refs_no_ignored_with_pos = dict()


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


def blast_results_to_pairwise_16s_iden_dict(blastn_output, align_len_cutoff, cov_cutoff):
    pairwise_iden_dict = {}
    for match in open(blastn_output):
        match_split = match.strip().split('\t')
        query = match_split[0]
        subject = match_split[1]
        iden = float(match_split[2])
        align_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        coverage_q = float(align_len) * 100 / float(query_len)
        coverage_s = float(align_len) * 100 / float(subject_len)

        if (align_len >= align_len_cutoff) and (query != subject) and (coverage_q >= cov_cutoff) and (
                coverage_s >= cov_cutoff):
            query_to_subject_key = '__|__'.join(sorted([query, subject]))
            if query_to_subject_key not in pairwise_iden_dict:
                pairwise_iden_dict[query_to_subject_key] = iden
            else:
                if iden > pairwise_iden_dict[query_to_subject_key]:
                    pairwise_iden_dict[query_to_subject_key] = iden

    return pairwise_iden_dict


def check_cigar_quality(cigar_str, mismatch_cutoff, min_M_len, ref_pos, ref_len):
    r2_ctg_ref_cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        r2_aligned_len_ctg, r2_aligned_pct_ctg, r2_clipping_len_ctg, r2_clipping_pct_ctg, r2_mismatch_pct_ctg = get_cigar_stats(
            r2_ctg_ref_cigar_splitted)
        if r2_mismatch_pct_ctg <= mismatch_cutoff:
            # check aligned length
            if r2_aligned_len_ctg >= min_M_len:
                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar_str) or ('s' in cigar_str):
                    clip_in_middle = True
                    if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + r2_aligned_len_ctg - 1) == ref_len:
                            clip_in_middle = False

                if clip_in_middle is False:
                    qualified_cigar = True

    return qualified_cigar


def sep_path_basename_ext(file_in):
    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = ''

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def sort_csv_by_col(file_in, file_out, col_header):
    df_in = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


def get_mean_iden_list(linked_16s_list, pairwise_16s_iden_dict):

    mean_iden_list = []
    for each_16s_1 in linked_16s_list:
        current_16s_iden_list = []
        for each_16s_2 in linked_16s_list:
            if each_16s_1 != each_16s_2:
                key = '__|__'.join(sorted([each_16s_1, each_16s_2]))
                key_iden = pairwise_16s_iden_dict.get(key, 0)
                current_16s_iden_list.append(key_iden)
        mean_iden = sum(current_16s_iden_list) / len(current_16s_iden_list)
        mean_iden = float("{0:.3f}".format(mean_iden))
        mean_iden_list.append(mean_iden)

    return mean_iden_list


# def filter_linkages_iteratively_new(sorted_file_in, pairwise_16s_iden_dict, within_genome_16s_divergence_cutoff, marker_len_dict,
#                                     min_linkages, within_gnm_linkage_num_diff, file_out,
#                                     marker_to_gnm_linking_cigar_dict_16s_side,
#                                     marker_to_gnm_linking_cigar_dict_ctg_side,
#                                     marker_to_ctg_gnm_Key_connector):
#
#     # filter linkage
#     gnm_max_link_num_dict = {}
#     file_out_handle = open(file_out, 'w')
#     MarkerGene_with_assignment = set()
#     MarkerGene_to_be_ignored = set()
#     current_gnm = ''
#     current_gnm_best_16s_list = []
#     current_gnm_highest_link_num = 0
#     gnm_with_assignment = set()
#     gnm_to_assignmed_16s_dict = dict()
#     best_16s_processed_gnm_list = set()
#     for each_match in open(sorted_file_in):
#         if each_match.startswith('MarkerGene,GenomicSeq,Number'):
#             file_out_handle.write(each_match)
#         else:
#             match_split = each_match.strip().split(',')
#             MarkerGene = match_split[0][12:]
#             MarkerGene_len = marker_len_dict[MarkerGene]
#             GenomicSeq = match_split[1][12:]
#             linkage_num = int(match_split[2])
#             MarkerGene_to_GenomicSeq_key = '%s%s%s' % (MarkerGene, marker_to_ctg_gnm_Key_connector, GenomicSeq)
#
#             if linkage_num >= min_linkages:
#
#                 ########### quality pre-check (I guess) ###########
#                 # failed linkage will be added to MarkerGene_to_be_ignored
#
#                 # first check if linked to conserved regions ???
#                 already_assigned_16s_list = []
#                 iden_with_already_assigned_16s_list = []
#                 for already_assigned_16s in MarkerGene_with_assignment:
#                     current_key = '__|__'.join(sorted([already_assigned_16s, MarkerGene]))
#                     current_key_value = pairwise_16s_iden_dict.get(current_key, 0)
#                     already_assigned_16s_list.append(already_assigned_16s)
#                     iden_with_already_assigned_16s_list.append(current_key_value)
#
#                 if len(already_assigned_16s_list) > 0:
#                     sorted_best_matched_16s_list = [[seq_id, mean_iden] for mean_iden, seq_id in sorted(zip(iden_with_already_assigned_16s_list, already_assigned_16s_list), reverse=True)]
#                     best_matched_marker      = sorted_best_matched_16s_list[0][0]
#                     best_matched_marker_iden = sorted_best_matched_16s_list[0][1]
#                     best_matched_marker_len  = marker_len_dict[best_matched_marker]
#                     if ((best_matched_marker_len - MarkerGene_len) >= 200) and (best_matched_marker_iden >= 99):
#
#                         # get clp pct at gnm level
#                         linking_cigar_16s_side = marker_to_gnm_linking_cigar_dict_16s_side[MarkerGene_to_GenomicSeq_key]
#                         linking_cigar_ctg_side = marker_to_gnm_linking_cigar_dict_ctg_side[MarkerGene_to_GenomicSeq_key]
#                         linking_cigar_16s_side_clp = [i for i in linking_cigar_16s_side if (('S' in i) or ('s' in i))]
#                         linking_cigar_ctg_side_clp = [i for i in linking_cigar_ctg_side if (('S' in i) or ('s' in i))]
#                         linking_cigar_16s_side_clp_pct = len(linking_cigar_16s_side_clp) * 100 / len(linking_cigar_16s_side)
#                         linking_cigar_ctg_side_clp_pct = len(linking_cigar_ctg_side_clp) * 100 / len(linking_cigar_ctg_side)
#
#                         if (linking_cigar_16s_side_clp_pct >= 60) and (linking_cigar_ctg_side_clp_pct >= 60):
#                             MarkerGene_to_be_ignored.add(MarkerGene)
#
#                 ###################### process here ######################
#
#                 if (MarkerGene not in MarkerGene_with_assignment) and (MarkerGene not in MarkerGene_to_be_ignored):
#
#                     if current_gnm == '':
#                         current_gnm = GenomicSeq
#                         current_gnm_best_16s_list.append(MarkerGene)
#                         current_gnm_highest_link_num = linkage_num
#
#                     elif current_gnm == GenomicSeq:
#
#                         if linkage_num == current_gnm_highest_link_num:
#                             current_gnm_best_16s_list.append(MarkerGene)
#                         else:
#                             # process markers with highest number of linking reads
#                             if current_gnm not in best_16s_processed_gnm_list:
#
#                                 gnm_max_link_num_dict[current_gnm] = current_gnm_highest_link_num
#
#                                 if len(current_gnm_best_16s_list) == 1:
#                                     file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (current_gnm_best_16s_list[0], current_gnm, current_gnm_highest_link_num))
#                                     gnm_with_assignment.add(current_gnm)
#                                     MarkerGene_with_assignment.add(current_gnm_best_16s_list[0])
#                                     gnm_to_assignmed_16s_dict[current_gnm] = {current_gnm_best_16s_list[0]}
#                                     best_16s_processed_gnm_list.add(current_gnm)
#                                 else:
#                                     # process here
#                                     current_gnm_best_16s_mean_iden_list = get_mean_iden_list(current_gnm_best_16s_list, pairwise_16s_iden_dict)
#                                     sorted_best_16s_list = [seq_id for mean_iden, seq_id in sorted(zip(current_gnm_best_16s_mean_iden_list, current_gnm_best_16s_list), reverse=True)]
#
#                                     add_index = 0
#                                     for linked_16s in sorted_best_16s_list:
#                                         if add_index == 0:
#                                             file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
#                                             gnm_with_assignment.add(current_gnm)
#                                             MarkerGene_with_assignment.add(linked_16s)
#                                             if current_gnm not in gnm_to_assignmed_16s_dict:
#                                                 gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
#                                             else:
#                                                 gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
#                                         else:
#                                             # get identity list with already assigned 16s
#                                             already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(current_gnm)
#                                             iden_with_already_assigned_16s_list = []
#                                             for each_assigned_16s in already_assigned_16s_list:
#                                                 marker_key = '__|__'.join(sorted([each_assigned_16s, linked_16s]))
#                                                 marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
#                                                 iden_with_already_assigned_16s_list.append(marker_key_iden)
#                                             if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:
#                                                 file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
#                                                 gnm_with_assignment.add(current_gnm)
#                                                 MarkerGene_with_assignment.add(linked_16s)
#                                                 if current_gnm not in gnm_to_assignmed_16s_dict:
#                                                     gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
#                                                 else:
#                                                     gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
#                                         add_index += 1
#                                     best_16s_processed_gnm_list.add(current_gnm)
#
#                             # process current one here
#                             # get identity list with already assigned 16s
#                             already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(GenomicSeq)
#                             iden_with_already_assigned_16s_list = []
#                             for each_assigned_16s in already_assigned_16s_list:
#                                 marker_key = '__|__'.join(sorted([each_assigned_16s, MarkerGene]))
#                                 marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
#                                 iden_with_already_assigned_16s_list.append(marker_key_iden)
#                             if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:
#
#                                 gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
#                                 if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
#                                     file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (MarkerGene, GenomicSeq, linkage_num))
#                                     gnm_with_assignment.add(GenomicSeq)
#                                     MarkerGene_with_assignment.add(MarkerGene)
#                                     if GenomicSeq not in gnm_to_assignmed_16s_dict:
#                                         gnm_to_assignmed_16s_dict[GenomicSeq] = {MarkerGene}
#                                     else:
#                                         gnm_to_assignmed_16s_dict[GenomicSeq].add(MarkerGene)
#                                 else:
#                                     # check length, if shorter than assigned, ignore this one by adding it to MarkerGene_to_be_ignored
#                                     MarkerGene_len = marker_len_dict[MarkerGene]
#                                     already_assigned_16s_len_list = [marker_len_dict[s16] for s16 in already_assigned_16s_list]
#                                     median_assign_16s_len = np.median(already_assigned_16s_len_list)
#                                     if (median_assign_16s_len - MarkerGene_len) > 50:
#                                         MarkerGene_to_be_ignored.add(MarkerGene)
#                     else:
#                         # first check if previous one processed
#                         if current_gnm not in best_16s_processed_gnm_list:
#
#                             gnm_max_link_num_dict[current_gnm] = current_gnm_highest_link_num
#
#                             # process unprocessed
#                             if len(current_gnm_best_16s_list) == 1:
#                                 file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (current_gnm_best_16s_list[0], current_gnm, current_gnm_highest_link_num))
#                                 gnm_with_assignment.add(current_gnm)
#                                 MarkerGene_with_assignment.add(current_gnm_best_16s_list[0])
#                                 gnm_to_assignmed_16s_dict[current_gnm] = {current_gnm_best_16s_list[0]}
#                                 best_16s_processed_gnm_list.add(current_gnm)
#                             else:
#                                 current_gnm_best_16s_mean_iden_list = get_mean_iden_list(current_gnm_best_16s_list, pairwise_16s_iden_dict)
#                                 sorted_best_16s_list = [seq_id for mean_iden, seq_id in sorted(zip(current_gnm_best_16s_mean_iden_list, current_gnm_best_16s_list), reverse=True)]
#
#                                 add_index = 0
#                                 for linked_16s in sorted_best_16s_list:
#                                     if add_index == 0:
#                                         file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
#                                         gnm_with_assignment.add(current_gnm)
#                                         MarkerGene_with_assignment.add(linked_16s)
#                                         if current_gnm not in gnm_to_assignmed_16s_dict:
#                                             gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
#                                         else:
#                                             gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
#                                     else:
#                                         # get identity list with already assigned 16s
#                                         already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(current_gnm)
#                                         iden_with_already_assigned_16s_list = []
#                                         for each_assigned_16s in already_assigned_16s_list:
#                                             marker_key = '__|__'.join(sorted([each_assigned_16s, linked_16s]))
#                                             marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
#                                             iden_with_already_assigned_16s_list.append(marker_key_iden)
#
#                                         if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:
#                                             file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
#                                             gnm_with_assignment.add(current_gnm)
#                                             MarkerGene_with_assignment.add(linked_16s)
#                                             if current_gnm not in gnm_to_assignmed_16s_dict:
#                                                 gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
#                                             else:
#                                                 gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
#                                     add_index += 1
#                                 best_16s_processed_gnm_list.add(current_gnm)
#
#                         # then process current one
#                         if GenomicSeq in best_16s_processed_gnm_list:
#                             # get identity list with already assigned 16s
#                             already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(GenomicSeq)
#
#
#                             iden_with_already_assigned_16s_list = []
#                             for each_assigned_16s in already_assigned_16s_list:
#                                 marker_key = '__|__'.join(sorted([each_assigned_16s, MarkerGene]))
#                                 marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
#                                 iden_with_already_assigned_16s_list.append(marker_key_iden)
#
#                             if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:
#
#                                 gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
#                                 if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
#                                     file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (MarkerGene, GenomicSeq, linkage_num))
#                                     gnm_with_assignment.add(GenomicSeq)
#                                     MarkerGene_with_assignment.add(MarkerGene)
#                                     if GenomicSeq not in gnm_to_assignmed_16s_dict:
#                                         gnm_to_assignmed_16s_dict[GenomicSeq] = {MarkerGene}
#                                     else:
#                                         gnm_to_assignmed_16s_dict[GenomicSeq].add(MarkerGene)
#                                 else:
#                                     # check length, if shorter than assigned, ignore this one by adding it to MarkerGene_with_assignment
#                                     already_assigned_16s_len_list = [marker_len_dict[s16] for s16 in already_assigned_16s_list]
#                                     median_assign_16s_len = np.median(already_assigned_16s_len_list)
#                                     if (median_assign_16s_len - MarkerGene_len) > 50:
#                                         MarkerGene_to_be_ignored.add(MarkerGene)
#                         else:
#                             current_gnm = GenomicSeq
#                             current_gnm_best_16s_list = [MarkerGene]
#                             current_gnm_highest_link_num = linkage_num
#     file_out_handle.close()
#

def get_gnm_16s_to_ignore(link_stats_combined_sorted, pairwise_16s_iden_dict, min_iden_16s, min_linking_num):

    # get gnm_best_16s_dict
    gnm_best_16s_link_num_dict = {}
    gnm_best_16s_dict = {}
    for each_match in open(link_stats_combined_sorted):
        if not each_match.endswith(',Number\n'):
            match_split = each_match.strip().split(',')
            MarkerGene = match_split[0][12:]
            GenomicSeq = match_split[1][12:]
            linkage_num = int(match_split[2])
            if linkage_num >= min_linking_num:
                if GenomicSeq not in gnm_best_16s_link_num_dict:
                    gnm_best_16s_link_num_dict[GenomicSeq] = linkage_num
                    gnm_best_16s_dict[GenomicSeq] = {MarkerGene}
                else:
                    if linkage_num == gnm_best_16s_link_num_dict[GenomicSeq]:
                        gnm_best_16s_dict[GenomicSeq].add(MarkerGene)

    # get gnm to ignore and best markers to ignore
    gnm_to_ignore = set()
    markers_to_ignore_best = set()
    for each_gnm in gnm_best_16s_dict:
        gnm_best_16s_list = gnm_best_16s_dict[each_gnm]
        if len(gnm_best_16s_list) > 1:
            all_combination_set = set()
            for left_16s in gnm_best_16s_list:
                for right_16s in gnm_best_16s_list:
                    if left_16s != right_16s:
                        marker_to_marker_key = '__|__'.join(sorted([left_16s, right_16s]))
                        all_combination_set.add(marker_to_marker_key)

            best_16s_pairwise_iden_list = []
            for each_combination in all_combination_set:
                marker_to_marker_iden = pairwise_16s_iden_dict.get(each_combination, 0)
                best_16s_pairwise_iden_list.append(marker_to_marker_iden)

            min_best_16s_iden = min(best_16s_pairwise_iden_list)
            if min_best_16s_iden <= min_iden_16s:
                gnm_to_ignore.add(each_gnm)
                for j in gnm_best_16s_list:
                    markers_to_ignore_best.add(j)

    # get all markers to ignore
    markers_to_ignore_all = set()
    for each_pair_16s in pairwise_16s_iden_dict:
        each_pair_16s_split = each_pair_16s.split('__|__')
        l_16s = each_pair_16s_split[0]
        r_16s = each_pair_16s_split[1]
        iden_16s = pairwise_16s_iden_dict[each_pair_16s]
        if iden_16s >= min_iden_16s:
            if l_16s in markers_to_ignore_best:
                markers_to_ignore_all.add(r_16s)
            if r_16s in markers_to_ignore_best:
                markers_to_ignore_all.add(l_16s)

    return gnm_to_ignore, markers_to_ignore_best, markers_to_ignore_all


def filter_linkages_iteratively_new2(sorted_file_in, pairwise_16s_iden_dict, within_genome_16s_divergence_cutoff,
                                    marker_len_dict,
                                    min_linkages, within_gnm_linkage_num_diff, file_out,
                                    marker_to_gnm_linking_cigar_dict_16s_side,
                                    marker_to_gnm_linking_cigar_dict_ctg_side,
                                    marker_to_ctg_gnm_Key_connector):

    # filter linkage
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    MarkerGene_to_be_ignored = set()
    gnm_with_assignment = set()
    gnm_best_16s_link_num_dict = {}
    gnm_assigned_16s_dict = {}
    for each_match in open(sorted_file_in):
        if each_match.startswith('MarkerGene,GenomicSeq,Number'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            MarkerGene = match_split[0][12:]
            MarkerGene_len = marker_len_dict[MarkerGene]
            GenomicSeq = match_split[1][12:]
            linkage_num = int(match_split[2])
            MarkerGene_to_GenomicSeq_key = '%s%s%s' % (MarkerGene, marker_to_ctg_gnm_Key_connector, GenomicSeq)

            if linkage_num >= min_linkages:

                ########### quality pre-check (I guess) ###########
                # failed linkage will be added to MarkerGene_to_be_ignored

                # first check if linked to conserved regions ???
                already_assigned_16s_list = []
                iden_with_already_assigned_16s_list = []
                for already_assigned_16s in MarkerGene_with_assignment:
                    current_key = '__|__'.join(sorted([already_assigned_16s, MarkerGene]))
                    current_key_value = pairwise_16s_iden_dict.get(current_key, 0)
                    already_assigned_16s_list.append(already_assigned_16s)
                    iden_with_already_assigned_16s_list.append(current_key_value)

                if len(already_assigned_16s_list) > 0:
                    sorted_best_matched_16s_list = [[seq_id, mean_iden] for mean_iden, seq_id in sorted(
                        zip(iden_with_already_assigned_16s_list, already_assigned_16s_list), reverse=True)]
                    best_matched_marker = sorted_best_matched_16s_list[0][0]
                    best_matched_marker_iden = sorted_best_matched_16s_list[0][1]
                    best_matched_marker_len = marker_len_dict[best_matched_marker]
                    if ((best_matched_marker_len - MarkerGene_len) >= 200) and (best_matched_marker_iden >= 99):

                        # get clp pct at gnm level
                        linking_cigar_16s_side = marker_to_gnm_linking_cigar_dict_16s_side[MarkerGene_to_GenomicSeq_key]
                        linking_cigar_ctg_side = marker_to_gnm_linking_cigar_dict_ctg_side[MarkerGene_to_GenomicSeq_key]
                        linking_cigar_16s_side_clp = [i for i in linking_cigar_16s_side if (('S' in i) or ('s' in i))]
                        linking_cigar_ctg_side_clp = [i for i in linking_cigar_ctg_side if (('S' in i) or ('s' in i))]
                        linking_cigar_16s_side_clp_pct = len(linking_cigar_16s_side_clp) * 100 / len(linking_cigar_16s_side)
                        linking_cigar_ctg_side_clp_pct = len(linking_cigar_ctg_side_clp) * 100 / len(linking_cigar_ctg_side)

                        if (linking_cigar_16s_side_clp_pct >= 60) and (linking_cigar_ctg_side_clp_pct >= 60):
                            MarkerGene_to_be_ignored.add(MarkerGene)

                ###################### process here ######################

                # get gnm_to_ignore, markers_to_ignore_best, markers_to_ignore_all
                gnm_to_ignore_multi_best, markers_to_ignore_multi_best, markers_to_ignore_multi_best_all = get_gnm_16s_to_ignore(sorted_file_in,
                                                                                                     pairwise_16s_iden_dict,
                                                                                                     within_genome_16s_divergence_cutoff,
                                                                                                     linkage_num)

                if (MarkerGene not in markers_to_ignore_multi_best_all) and (GenomicSeq not in gnm_to_ignore_multi_best):
                    if (MarkerGene not in MarkerGene_with_assignment) and (MarkerGene not in MarkerGene_to_be_ignored):

                        if GenomicSeq not in gnm_with_assignment:
                            file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (MarkerGene, GenomicSeq, linkage_num))
                            gnm_with_assignment.add(GenomicSeq)
                            MarkerGene_with_assignment.add(MarkerGene)
                            gnm_assigned_16s_dict[GenomicSeq] = [MarkerGene]
                            gnm_best_16s_link_num_dict[GenomicSeq] = linkage_num
                        else:
                            gnm_assigned_16s_list = gnm_assigned_16s_dict[GenomicSeq]
                            gnm_best_assigned_16s_link_num = gnm_best_16s_link_num_dict[GenomicSeq]

                            # get iden list with already assigned 16s
                            iden_with_assigned_16s_list = []
                            for each_assigned_16s in gnm_assigned_16s_list:
                                marker_to_maker_key = '__|__'.join(sorted([MarkerGene, each_assigned_16s]))
                                marker_key_iden = pairwise_16s_iden_dict.get(marker_to_maker_key, 0)
                                iden_with_assigned_16s_list.append(marker_key_iden)

                            # get min iden with already assigned 16s
                            min_iden_with_assigned_16s = min(iden_with_assigned_16s_list)

                            # if min iden higher than within_genome_16s_divergence_cutoff
                            if min_iden_with_assigned_16s >= within_genome_16s_divergence_cutoff:

                                # filter here
                                if (linkage_num * 100 / gnm_best_assigned_16s_link_num) >= within_gnm_linkage_num_diff:
                                    file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (MarkerGene, GenomicSeq, linkage_num))
                                    MarkerGene_with_assignment.add(MarkerGene)
                                    gnm_assigned_16s_dict[GenomicSeq].append(MarkerGene)
                                else:
                                    MarkerGene_to_be_ignored.add(MarkerGene)
    file_out_handle.close()


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
                    for r1_16s_ref in current_read_base_r1_16s_ref_dict:
                        r1_matched_pos_dict = current_read_base_r1_16s_ref_dict[r1_16s_ref]

                        # one read need to mapped to one 16S only for one time
                        if len(r1_matched_pos_dict) > 1:
                            refs_to_ignore.add(r1_16s_ref)
                        else:
                            r1_16s_ref_pos = list(r1_matched_pos_dict.keys())[0]
                            r1_16s_ref_cigar = r1_matched_pos_dict[r1_16s_ref_pos]
                            r1_16s_ref_cigar_splitted = cigar_splitter(r1_16s_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r1_16s_ref_cigar_splitted)
                            if both_end_clp is True:
                                refs_to_ignore.add(r1_16s_ref)
                            else:
                                # check mismatch
                                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                    r1_16s_ref_cigar_splitted)
                                if r1_ref_min_mismatch == 'NA':
                                    refs_to_ignore.add(r1_16s_ref)
                                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                                    refs_to_ignore.add(r1_16s_ref)
                                else:
                                    # check aligned length
                                    if r1_aligned_len < min_M_len_16s:
                                        refs_to_ignore.add(r1_16s_ref)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r1_16s_ref_cigar) or ('s' in r1_16s_ref_cigar):
                                            clip_in_middle = True
                                            if (r1_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r1_16s_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r1_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r1_16s_ref_pos + r1_aligned_len - 1) == marker_len_dict[r1_16s_ref]:
                                                    clip_in_middle = False

                                        # exclude the ref if clp in the middle is True
                                        if clip_in_middle is True:
                                            refs_to_ignore.add(r1_16s_ref)
                                        else:
                                            r1_16s_refs_passed_qc[r1_16s_ref] = [r1_16s_ref_cigar]

                    ########## filter r2 16s refs ##########

                    r2_16s_refs_passed_qc = {}
                    for r2_16s_ref in current_read_base_r2_16s_ref_dict:
                        r2_matched_pos_dict = current_read_base_r2_16s_ref_dict[r2_16s_ref]

                        # one read need to mapped to one 16S only once
                        if len(r2_matched_pos_dict) > 1:
                            refs_to_ignore.add(r2_16s_ref)
                        else:
                            r2_16s_ref_pos = list(r2_matched_pos_dict.keys())[0]
                            r2_16s_ref_cigar = r2_matched_pos_dict[r2_16s_ref_pos]
                            r2_16s_ref_cigar_splitted = cigar_splitter(r2_16s_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r2_16s_ref_cigar_splitted)
                            if both_end_clp is True:
                                refs_to_ignore.add(r2_16s_ref)
                            else:
                                # check mismatch
                                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                    r2_16s_ref_cigar_splitted)
                                if r2_ref_min_mismatch == 'NA':
                                    refs_to_ignore.add(r2_16s_ref)
                                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                                    refs_to_ignore.add(r2_16s_ref)
                                else:
                                    # check aligned length
                                    if r2_aligned_len < min_M_len_16s:
                                        refs_to_ignore.add(r2_16s_ref)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r2_16s_ref_cigar) or ('s' in r2_16s_ref_cigar):
                                            clip_in_middle = True
                                            if (r2_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r2_16s_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r2_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r2_16s_ref_pos + r2_aligned_len - 1) == marker_len_dict[r2_16s_ref]:
                                                    clip_in_middle = False

                                        # exclude the ref if clp in the middle is True
                                        if clip_in_middle is True:
                                            refs_to_ignore.add(r2_16s_ref)
                                        else:
                                            r2_16s_refs_passed_qc[r2_16s_ref] = [r2_16s_ref_cigar]

                    #################################

                    r1_16s_refs_no_ignored = {key: value for key, value in r1_16s_refs_passed_qc.items() if
                                              key not in refs_to_ignore}
                    r2_16s_refs_no_ignored = {key: value for key, value in r2_16s_refs_passed_qc.items() if
                                              key not in refs_to_ignore}

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

                    # if current_read_base in MappingRecord_dict:
                    #     print('%s\tr1_16s_ref_dict\t%s' % (current_read_base, MappingRecord_dict[current_read_base].r1_16s_ref_dict))
                    #     print('%s\tr1_16s_ref_dict\t%s' % (current_read_base, r12_16s_ref_dict_to_str(MappingRecord_dict[current_read_base].r1_16s_ref_dict)))
                    #
                    #     print('%s\tr2_16s_ref_dict\t%s' % (current_read_base, MappingRecord_dict[current_read_base].r2_16s_ref_dict))
                    #     print('%s\tr2_16s_ref_dict\t%s' % (current_read_base, r12_16s_ref_dict_to_str(MappingRecord_dict[current_read_base].r2_16s_ref_dict)))
                    #
                    #     print('%s\tr1_16s_refs_no_ignored\t%s' % (current_read_base, MappingRecord_dict[current_read_base].r1_16s_refs_no_ignored))
                    #     print('%s\tr1_16s_refs_no_ignored\t%s' % (current_read_base, no_ignored_dict_to_str(MappingRecord_dict[current_read_base].r1_16s_refs_no_ignored)))
                    #
                    #     print('%s\tr2_16s_refs_no_ignored\t%s' % (current_read_base, MappingRecord_dict[current_read_base].r2_16s_refs_no_ignored))
                    #     print('%s\tr2_16s_refs_no_ignored\t%s' % (current_read_base, no_ignored_dict_to_str(MappingRecord_dict[current_read_base].r2_16s_refs_no_ignored)))
                    #
                    #     print('%s\tshared_16s_refs_no_ignored\t%s' % (current_read_base, MappingRecord_dict[current_read_base].shared_16s_refs_no_ignored))
                    #     print('%s\tshared_16s_refs_no_ignored\t%s' % (current_read_base, no_ignored_dict_to_str(MappingRecord_dict[current_read_base].shared_16s_refs_no_ignored)))
                    #
                    #     print()

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


def get_regions_to_ignore_from_barrnap_output(combined_barrnap_gff, ctg_len_dict):

    ctg_ignore_region_dict = {}

    for each_line in open(combined_barrnap_gff):
        if (not each_line.startswith('#')) and ('16S_rRNA' in each_line):
            each_line_split = each_line.strip().split('\t')
            ctg_id = each_line_split[0]
            ctg_len = ctg_len_dict[ctg_id]
            start_pos = int(each_line_split[3])
            end_pos = int(each_line_split[4])
            len_16s = end_pos - start_pos + 1
            left_gap = start_pos - 1
            right_gap = ctg_len - end_pos - 1

            if left_gap <= 50:
                if ctg_id not in ctg_ignore_region_dict:
                    ctg_ignore_region_dict[ctg_id] = {'left_end'}
                else:
                    ctg_ignore_region_dict[ctg_id].add('left_end')

            if right_gap <= 50:
                if ctg_id not in ctg_ignore_region_dict:
                    ctg_ignore_region_dict[ctg_id] = {'right_end'}
                else:
                    ctg_ignore_region_dict[ctg_id].add('right_end')

    return ctg_ignore_region_dict


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


def get_min_dist_to_ref_end(cigar_str, cigar_pos, ref_len):
    cigar_aln_len = get_cigar_aln_len(cigar_splitter(cigar_str))
    cigar_dist_to_left = cigar_pos - 1
    cigar_dist_to_right = ref_len - cigar_pos - cigar_aln_len + 1
    min_dist_to_mini_end = min(cigar_dist_to_left, cigar_dist_to_right)
    return min_dist_to_mini_end


def get_short_cigar_pct(cigar_list, short_M_len):
    short_cigar_num = 0
    for each_cigar in cigar_list:
        aln_len, aln_pct, clp_len, clp_pct, mis_pct = get_cigar_stats(cigar_splitter(each_cigar))
        if aln_len <= short_M_len:
            short_cigar_num += 1
    short_cigar_pct = short_cigar_num * 100 / len(cigar_list)
    return short_cigar_pct


class LinkingRecord:

    def __init__(self):

        self.linked_seq_l = ''
        self.linked_seq_r = ''

        self.linked_seq_len_l = 0
        self.linked_seq_len_r = 0

        self.linking_reads_base = []

        self.linking_reads_l = []
        self.linking_reads_r = []

        self.linking_cigar_l = []
        self.linking_cigar_r = []

        self.linking_pos_l = []
        self.linking_pos_r = []

        self.min_dist_to_end_l = []
        self.min_dist_to_end_r = []

        self.all_matched_pos_l = set()
        self.all_matched_pos_r = set()

def get_linked_to_ref_end_cigar_num_and_pct(min_dist_to_ref_end_list):
    linked_to_ref_end_cigar_num = 0
    for each_min_dist in min_dist_to_ref_end_list:
        if each_min_dist <= 10:
            linked_to_ref_end_cigar_num += 1
    linked_to_ref_end_cigar_pct = linked_to_ref_end_cigar_num * 100 / len(min_dist_to_ref_end_list)
    return linked_to_ref_end_cigar_num, linked_to_ref_end_cigar_pct


def get_all_matched_pos(linking_cigar_list, linking_pos_list):

    all_matched_pos = set()
    for (each_linking_cigar, each_linking_pos) in zip(linking_cigar_list, linking_pos_list):
        current_cigar_aligned_len = get_cigar_aln_len(cigar_splitter(each_linking_cigar))
        matched_seq_end = each_linking_pos + current_cigar_aligned_len
        pos_list = list(range(each_linking_pos, matched_seq_end))
        all_matched_pos.update(pos_list)
    return all_matched_pos


if __name__ ==  '__main__':

    ##################################################### file in ####################################################

    # step_1_wd                                   = '/Users/songweizhi/Desktop/tunning_rd1_GI'
    # input_16s_polished                          = '%s/file_in/GI_128_16S_0.999.QC.fasta'           % step_1_wd
    # combined_input_gnms                         = '%s/file_in/GI_refined_bins_combined_no_ending_16s.fa'            % step_1_wd
    # input_reads_to_16s_sam_MappingRecord_folder = '%s/file_in/GI_0719_128_45_45_1200_diff10_input_reads_to_16S_MappingRecord' % step_1_wd
    # blast_results_all_vs_all_16s                = '%s/file_in/GI_0719_128_45_45_1200_diff10_16S_all_vs_all_blastn.tab'        % step_1_wd
    # rd1_extracted_to_gnm_sam_reformatted_sorted = '%s/file_in/rd1_extracted_to_gnm_sorted.sam'                          % step_1_wd
    # combined_barrnap_gff                        = '%s/file_in/combined_barrnap.gff'                                     % step_1_wd

    step_1_wd                                   = '/Users/songweizhi/Desktop/tunning_rd1_Kelp'
    input_16s_polished                          = '%s/file_in/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.QC.fasta'   % step_1_wd
    combined_input_gnms                         = '%s/file_in/BH_ER_050417_refined_bins_combined_no_ending_16s.fa'      % step_1_wd
    input_reads_to_16s_sam_MappingRecord_folder = '%s/file_in/Kelp_0720_60_60_900_input_reads_to_16S_MappingRecord'     % step_1_wd
    blast_results_all_vs_all_16s                = '%s/file_in/Kelp_0720_60_60_900_16S_all_vs_all_blastn.tab'            % step_1_wd
    rd1_extracted_to_gnm_sam_reformatted_sorted = '%s/file_in/rd1_extracted_to_gnm_sorted.sam'                          % step_1_wd
    combined_barrnap_gff                        = '%s/file_in/combined_barrnap.gff'                                     % step_1_wd

    # step_1_wd                                   = '/Users/songweizhi/Desktop/tunning_rd1_Oral'
    # input_16s_polished                          = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_min1200.QC.fa'           % step_1_wd
    # combined_input_gnms                         = '%s/file_in/3_Oral_refined_MAGs_combined_no_ending_16s.fa'            % step_1_wd
    # input_reads_to_16s_sam_MappingRecord_folder = '%s/file_in/Oral_0713_45_45_min1200_input_reads_to_16S_MappingRecord' % step_1_wd
    # blast_results_all_vs_all_16s                = '%s/file_in/Oral_0713_45_45_min1200_16S_all_vs_all_blastn.tab'        % step_1_wd
    # rd1_extracted_to_gnm_sam_reformatted_sorted = '%s/file_in/rd1_extracted_to_gnm_sorted.sam'                          % step_1_wd
    # combined_barrnap_gff                        = '%s/file_in/combined_barrnap.gff'                                     % step_1_wd


    # step_1_wd                                   = '/Users/songweizhi/Desktop/tunning_rd1'
    # input_16s_polished                          = '%s/file_in/MBARC26_SILVA138_polished.QC.fasta'                             % step_1_wd
    # combined_input_gnms                         = '%s/file_in/Refined_refined_bins_renamed_combined_no_ending_16s.fa'         % step_1_wd
    # input_reads_to_16s_sam_MappingRecord_folder = '%s/file_in/MBARC26_0708_45_45_min1200_input_reads_to_16S_MappingRecord'    % step_1_wd
    # blast_results_all_vs_all_16s                = '%s/file_in/MBARC26_0708_45_45_min1200_16S_all_vs_all_blastn.tab'           % step_1_wd
    # rd1_extracted_to_gnm_sam_reformatted_sorted = '%s/file_in/rd1_extracted_to_gnm_reformatted_sorted.sam'                    % step_1_wd
    # combined_barrnap_gff                        = '%s/file_in/combined_barrnap.gff'                                           % step_1_wd

    min_M_len_16s = 45
    min_M_len_ctg = 45
    min_M_pct = 35
    short_M_len_16s = 75
    short_M_len_ctg = 75

    mismatch_cutoff = 2
    marker_to_ctg_gnm_Key_connector = '___M___'
    read_to_marker_connector = '___r___'
    gnm_to_ctg_connector = '___C___'
    min_aln_16s = 500
    min_cov_16s = 30
    mean_depth_dict_gnm = {}
    mean_depth_dict_16s = {}
    min_16s_gnm_multiple = 0
    min_iden_16s = 98.5
    min_link_num = 8
    within_gnm_linkage_num_diff = 10
    time_format = '[%Y-%m-%d %H:%M:%S]'
    num_threads = 4
    min_M_len_mini                                  = 75
    short_M_len_16s                                 = 75
    short_M_len_ctg                                 = 75
    short_M_len_mini                                = 75


    ################################################# define file name #################################################

    rd1_clp_pct_diff_txt            = '%s/rd1_clp_pct_diff.txt'             % step_1_wd
    rd1_clp_pct_diff_txt_to_ignore  = '%s/rd1_clp_pct_diff_to_ignore.txt'   % step_1_wd
    linkages_QC_txt                 = '%s/linkages_QC.txt'                  % step_1_wd
    link_stats_combined             = '%s/stats_combined.txt'               % step_1_wd
    link_stats_combined_sorted      = '%s/stats_combined_sorted.txt'        % step_1_wd
    link_stats_combined_filtered_s1 = '%s/stats_combined_filtered.txt'      % step_1_wd
    linking_reads_rd1               = '%s/linking_reads_rd1.txt'            % step_1_wd
    linked_by_clp_pct               = '%s/rd1_linked_by_clp_pct.txt'        % step_1_wd


    ##################################################### read in mp file ####################################################

    print('%s %s' % ((datetime.now().strftime(time_format)), 'read in mp file'))

    # get marker len dict
    marker_len_dict = {}
    for each_marker_record in SeqIO.parse(input_16s_polished, 'fasta'):
        marker_len_dict[each_marker_record.id] = len(each_marker_record.seq)

    splitted_sam_mp_file_re = '%s/*_mp.txt' % input_reads_to_16s_sam_MappingRecord_folder
    splitted_sam_mp_file_set = glob.glob(splitted_sam_mp_file_re)

    processing_index = 1
    MappingRecord_dict = {}
    read_base_to_pop = set()
    for each_mp_file in splitted_sam_mp_file_set:
        if processing_index % 30 == 0:
            print('Processing %s/%s' % (processing_index, len(splitted_sam_mp_file_set)))
        processing_index += 1
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
                    read_base_to_pop.add(current_read_base___id)
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
                    MappingRecord_dict[current_read_base___id].shared_16s_refs_no_ignored = current_read_base___shared_16s_refs_no_ignored

    for each_mp in MappingRecord_dict.copy():
        if each_mp in read_base_to_pop:
            MappingRecord_dict.pop(each_mp)


    ######################################### read mapping results of rd1 extracted mates into mp dict  #########################################

    print('%s %s' % ((datetime.now().strftime(time_format)), 'read in sam file'))

    ctg_len_dict = {}
    for each_ctg_record in SeqIO.parse(combined_input_gnms, 'fasta'):
        ctg_len_dict[each_ctg_record.id] = len(each_ctg_record.seq)

    ctg_ignore_region_dict = get_regions_to_ignore_from_barrnap_output(combined_barrnap_gff, ctg_len_dict)

    current_read_base = ''
    current_read_base_r1_ctg_ref_dict = dict()
    current_read_base_r2_ctg_ref_dict = dict()
    with open(rd1_extracted_to_gnm_sam_reformatted_sorted) as rd1_extracted_to_gnm_sam_reformatted_opened:
        for each_read in rd1_extracted_to_gnm_sam_reformatted_opened:
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
                            current_read_base_r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}

                elif read_id_base == current_read_base:
                    if cigar != '*':
                        if read_strand == '1':
                            if ref_id not in current_read_base_r1_ctg_ref_dict:
                                current_read_base_r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r1_ctg_ref_dict[ref_id][ref_pos] = cigar

                        if read_strand == '2':
                            if ref_id not in current_read_base_r2_ctg_ref_dict:
                                current_read_base_r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r2_ctg_ref_dict[ref_id][ref_pos] = cigar
                else:
                    ################################### analysis previous read refs ####################################

                    refs_to_ignore_ctg = set()

                    ########## filter r1 ctg refs ##########

                    r1_ctg_refs_passed_qc = {}
                    r1_ctg_refs_passed_qc_with_pos = {}
                    for r1_ctg_ref in current_read_base_r1_ctg_ref_dict:
                        r1_ctg_ref_matched_pos_dict = current_read_base_r1_ctg_ref_dict[r1_ctg_ref]

                        # one read need to mapped to one ctg only for one time
                        if len(r1_ctg_ref_matched_pos_dict) > 1:
                            refs_to_ignore_ctg.add(r1_ctg_ref)
                        else:
                            r1_ctg_ref_pos = list(r1_ctg_ref_matched_pos_dict.keys())[0]
                            r1_ctg_ref_cigar = r1_ctg_ref_matched_pos_dict[r1_ctg_ref_pos]
                            r1_ctg_ref_len = ctg_len_dict[r1_ctg_ref]
                            qualified_cigar = check_cigar_quality(r1_ctg_ref_cigar, mismatch_cutoff, min_M_len_ctg, r1_ctg_ref_pos, r1_ctg_ref_len)

                            if qualified_cigar is True:

                                # check if matched to regions need to be ignored
                                matched_to_r1_ref_ignored_region = False
                                if r1_ctg_ref in ctg_ignore_region_dict:
                                    r1_ctg_ref_ends_to_ignore = ctg_ignore_region_dict[r1_ctg_ref]
                                    for to_ignore_region in r1_ctg_ref_ends_to_ignore:
                                        if to_ignore_region == 'left_end':
                                            if r1_ctg_ref_pos <= 50:
                                                matched_to_r1_ref_ignored_region = True
                                        if to_ignore_region == 'right_end':
                                            aln_len, aln_pct, clp_len, clp_pct, mis_pct = get_cigar_stats(cigar_splitter(r1_ctg_ref_cigar))
                                            if (ctg_len_dict[r1_ctg_ref] - r1_ctg_ref_pos - aln_len) <= 50:
                                                matched_to_r1_ref_ignored_region = True

                                if matched_to_r1_ref_ignored_region is False:
                                    r1_ctg_refs_passed_qc[r1_ctg_ref] = [r1_ctg_ref_cigar]
                                    r1_ctg_refs_passed_qc_with_pos[r1_ctg_ref] = {r1_ctg_ref_pos: r1_ctg_ref_cigar}
                            else:
                                refs_to_ignore_ctg.add(r1_ctg_ref)

                    ########## filter r2 ctg refs ##########

                    r2_ctg_refs_passed_qc = {}
                    r2_ctg_refs_passed_qc_with_pos = {}
                    for r2_ctg_ref in current_read_base_r2_ctg_ref_dict:
                        r2_ctg_ref_matched_pos_dict = current_read_base_r2_ctg_ref_dict[r2_ctg_ref]

                        # one read need to mapped to one ctg only for one time
                        if len(r2_ctg_ref_matched_pos_dict) > 1:
                            refs_to_ignore_ctg.add(r2_ctg_ref)
                        else:
                            r2_ctg_ref_pos = list(r2_ctg_ref_matched_pos_dict.keys())[0]
                            r2_ctg_ref_cigar = r2_ctg_ref_matched_pos_dict[r2_ctg_ref_pos]
                            r2_ctg_ref_len = ctg_len_dict.get(r2_ctg_ref, 0)
                            qualified_cigar = check_cigar_quality(r2_ctg_ref_cigar, mismatch_cutoff, min_M_len_ctg, r2_ctg_ref_pos, r2_ctg_ref_len)

                            if qualified_cigar is True:

                                # check if matched to regions need to be ignored
                                matched_to_r2_ref_ignored_region = False
                                if r2_ctg_ref in ctg_ignore_region_dict:
                                    r2_ctg_ref_ends_to_ignore = ctg_ignore_region_dict[r2_ctg_ref]
                                    for to_ignore_region in r2_ctg_ref_ends_to_ignore:
                                        if to_ignore_region == 'left_end':
                                            if r2_ctg_ref_pos <= 50:
                                                matched_to_r2_ref_ignored_region = True
                                        if to_ignore_region == 'right_end':
                                            aln_len, aln_pct, clp_len, clp_pct, mis_pct = get_cigar_stats(cigar_splitter(r2_ctg_ref_cigar))
                                            if (ctg_len_dict[r2_ctg_ref] - r2_ctg_ref_pos - aln_len) <= 50:
                                                matched_to_r2_ref_ignored_region = True

                                if matched_to_r2_ref_ignored_region is False:
                                    r2_ctg_refs_passed_qc[r2_ctg_ref] = [r2_ctg_ref_cigar]
                                    r2_ctg_refs_passed_qc_with_pos[r2_ctg_ref] = {r2_ctg_ref_pos: r2_ctg_ref_cigar}
                            else:
                                refs_to_ignore_ctg.add(r2_ctg_ref)

                    ####################################################################################################

                    r1_ctg_refs_no_ignored          = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}
                    r2_ctg_refs_no_ignored          = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}
                    r1_ctg_refs_no_ignored_with_pos = {key: value for key, value in r1_ctg_refs_passed_qc_with_pos.items() if key not in refs_to_ignore_ctg}
                    r2_ctg_refs_no_ignored_with_pos = {key: value for key, value in r2_ctg_refs_passed_qc_with_pos.items() if key not in refs_to_ignore_ctg}

                    # no mate matched to ctg
                    if (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) == 0):
                        pass

                    # only r1 matched to ctg
                    elif (len(r1_ctg_refs_no_ignored) > 0) and (len(r2_ctg_refs_no_ignored) == 0):

                        if current_read_base not in MappingRecord_dict:
                            MappingRecord_dict[current_read_base] = MappingRecord()

                        MappingRecord_dict[current_read_base].matched_to_ctg = True
                        MappingRecord_dict[current_read_base].r1_ctg_refs_no_ignored = r1_ctg_refs_no_ignored
                        MappingRecord_dict[current_read_base].r1_ctg_ref_dict = current_read_base_r1_ctg_ref_dict

                    # only r2 matched to ctg
                    elif (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) > 0):

                        if current_read_base not in MappingRecord_dict:
                            MappingRecord_dict[current_read_base] = MappingRecord()

                        MappingRecord_dict[current_read_base].matched_to_ctg = True
                        MappingRecord_dict[current_read_base].r2_ctg_refs_no_ignored = r2_ctg_refs_no_ignored
                        MappingRecord_dict[current_read_base].r2_ctg_ref_dict = current_read_base_r2_ctg_ref_dict

                    # both r1 and r2 matched to ctg
                    else:
                        shared_ctg_ref_set = {key: [r1_ctg_refs_no_ignored[key][0], r2_ctg_refs_no_ignored[key][0]] for key
                                              in set(r1_ctg_refs_no_ignored).intersection(set(r2_ctg_refs_no_ignored))}
                        if len(shared_ctg_ref_set) > 0:

                            if current_read_base not in MappingRecord_dict:
                                MappingRecord_dict[current_read_base] = MappingRecord()

                            MappingRecord_dict[current_read_base].matched_to_ctg = True
                            MappingRecord_dict[current_read_base].shared_ctg_refs_no_ignored = shared_ctg_ref_set
                            MappingRecord_dict[current_read_base].r1_ctg_ref_dict = current_read_base_r1_ctg_ref_dict
                            MappingRecord_dict[current_read_base].r2_ctg_ref_dict = current_read_base_r2_ctg_ref_dict

                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_ctg_ref_dict = dict()
                    current_read_base_r2_ctg_ref_dict = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}

    ##################################################### get linkages from MappingRecord_dict #####################################################


    def add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
                                          ref_dict_16s, ref_dict_16s_from,
                                          ref_dict_ctg, ref_dict_ctg_from,
                                          marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
                                          marker_to_ctg_gnm_Key_connector):

        for each_16s_ref in ref_dict_16s:
            for each_ctg_ref in ref_dict_ctg:
                marker_to_ctg_key = '%s%s%s' % (each_16s_ref, marker_to_ctg_gnm_Key_connector, each_ctg_ref)

                if marker_to_ctg_key not in marker_to_ctg_LinkingRecord_dict:
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key] = LinkingRecord()
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linked_seq_l = each_16s_ref
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linked_seq_r = each_ctg_ref
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linked_seq_len_l = marker_len_dict[each_16s_ref]
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linked_seq_len_r = ctg_len_dict[each_ctg_ref]

                marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_reads_base.append(qualified_read)

                if ref_dict_16s_from in ['r1', 'both']:
                    marker_side_cigar_r1 = \
                    list(MappingRecord_dict[qualified_read].r1_16s_ref_dict[each_16s_ref].values())[0]
                    marker_side_cigar_pos_r1 = \
                    list(MappingRecord_dict[qualified_read].r1_16s_ref_dict[each_16s_ref].keys())[0]
                    marker_side_cigar_r1_min_end_dist = get_min_dist_to_ref_end(marker_side_cigar_r1,
                                                                                marker_side_cigar_pos_r1,
                                                                                marker_len_dict[each_16s_ref])
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_pos_l.append(marker_side_cigar_pos_r1)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_cigar_l.append(marker_side_cigar_r1)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].min_dist_to_end_l.append(
                        marker_side_cigar_r1_min_end_dist)

                if ref_dict_16s_from in ['r2', 'both']:
                    marker_side_cigar_r2 = \
                    list(MappingRecord_dict[qualified_read].r2_16s_ref_dict[each_16s_ref].values())[0]
                    marker_side_cigar_pos_r2 = \
                    list(MappingRecord_dict[qualified_read].r2_16s_ref_dict[each_16s_ref].keys())[0]
                    marker_side_cigar_r2_min_end_dist = get_min_dist_to_ref_end(marker_side_cigar_r2,
                                                                                marker_side_cigar_pos_r2,
                                                                                marker_len_dict[each_16s_ref])
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_pos_l.append(marker_side_cigar_pos_r2)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_cigar_l.append(marker_side_cigar_r2)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].min_dist_to_end_l.append(
                        marker_side_cigar_r2_min_end_dist)

                if ref_dict_ctg_from in ['r1', 'both']:
                    ctg_side_cigar_r1 = list(MappingRecord_dict[qualified_read].r1_ctg_ref_dict[each_ctg_ref].values())[
                        0]
                    ctg_side_cigar_pos_r1 = \
                    list(MappingRecord_dict[qualified_read].r1_ctg_ref_dict[each_ctg_ref].keys())[0]
                    ctg_side_cigar_r1_min_end_dist = get_min_dist_to_ref_end(ctg_side_cigar_r1, ctg_side_cigar_pos_r1,
                                                                             ctg_len_dict[each_ctg_ref])
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_pos_r.append(ctg_side_cigar_pos_r1)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_cigar_r.append(ctg_side_cigar_r1)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].min_dist_to_end_r.append(
                        ctg_side_cigar_r1_min_end_dist)

                if ref_dict_ctg_from in ['r2', 'both']:
                    ctg_side_cigar_r2 = list(MappingRecord_dict[qualified_read].r2_ctg_ref_dict[each_ctg_ref].values())[
                        0]
                    ctg_side_cigar_pos_r2 = \
                    list(MappingRecord_dict[qualified_read].r2_ctg_ref_dict[each_ctg_ref].keys())[0]
                    ctg_side_cigar_r2_min_end_dist = get_min_dist_to_ref_end(ctg_side_cigar_r2, ctg_side_cigar_pos_r2,
                                                                             ctg_len_dict[each_ctg_ref])
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_pos_r.append(ctg_side_cigar_pos_r2)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].linking_cigar_r.append(ctg_side_cigar_r2)
                    marker_to_ctg_LinkingRecord_dict[marker_to_ctg_key].min_dist_to_end_r.append(
                        ctg_side_cigar_r2_min_end_dist)


    print('%s %s' % ((datetime.now().strftime(time_format)), 'get linkages from MappingRecord_dict'))


    marker_to_ctg_LinkingRecord_dict = {}
    for qualified_read in MappingRecord_dict:

        r1_16s_refs     = MappingRecord_dict[qualified_read].r1_16s_refs_no_ignored
        r2_16s_refs     = MappingRecord_dict[qualified_read].r2_16s_refs_no_ignored
        shared_16s_refs = MappingRecord_dict[qualified_read].shared_16s_refs_no_ignored
        r1_ctg_refs     = MappingRecord_dict[qualified_read].r1_ctg_refs_no_ignored
        r2_ctg_refs     = MappingRecord_dict[qualified_read].r2_ctg_refs_no_ignored
        shared_ctg_refs = MappingRecord_dict[qualified_read].shared_ctg_refs_no_ignored

        if (len(shared_16s_refs) > 0) or (len(shared_ctg_refs) > 0):
            if (len(shared_16s_refs) > 0) and (len(shared_ctg_refs) > 0):
                add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
                                                  shared_16s_refs, 'both', shared_ctg_refs, 'both',
                                                  marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
                                                  marker_to_ctg_gnm_Key_connector)

            # elif (len(shared_16s_refs) > 0) and (len(shared_ctg_refs) == 0):
            #
            #     if len(r1_ctg_refs) > 0:
            #         add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
            #                                           shared_16s_refs, 'both', r1_ctg_refs, 'r1',
            #                                           marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
            #                                           marker_to_ctg_gnm_Key_connector)
            #
            #     if len(r2_ctg_refs) > 0:
            #         add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
            #                                           shared_16s_refs, 'both', r2_ctg_refs, 'r2',
            #                                           marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
            #                                           marker_to_ctg_gnm_Key_connector)
            #
            # elif (len(shared_16s_refs) == 0) and (len(shared_ctg_refs) > 0):
            #     if len(r1_16s_refs) > 0:
            #         add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
            #                                           r1_16s_refs, 'r1', shared_ctg_refs, 'both',
            #                                           marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
            #                                           marker_to_ctg_gnm_Key_connector)
            #
            #     if len(r2_16s_refs) > 0:
            #         add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
            #                                           r2_16s_refs, 'r2', shared_ctg_refs, 'both',
            #                                           marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
            #                                           marker_to_ctg_gnm_Key_connector)

        else:
            if (len(r1_16s_refs) > 0) and (len(r2_ctg_refs) > 0):
                r2_gnm_refs = {ctg.split(gnm_to_ctg_connector)[0] for ctg in r2_ctg_refs}
                if len(r2_gnm_refs) == 1:
                    add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
                                                      r1_16s_refs, 'r1', r2_ctg_refs, 'r2',
                                                      marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
                                                      marker_to_ctg_gnm_Key_connector)

            if (len(r2_16s_refs) > 0) and (len(r1_ctg_refs) > 0):
                r1_gnm_refs = {ctg.split(gnm_to_ctg_connector)[0] for ctg in r1_ctg_refs}
                if len(r1_gnm_refs) == 1:
                    add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
                                                      r2_16s_refs, 'r2', r1_ctg_refs, 'r1',
                                                      marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
                                                      marker_to_ctg_gnm_Key_connector)

            if (len(r1_16s_refs) > 0) and (len(r1_ctg_refs) > 0):
                r1_gnm_refs = {ctg.split(gnm_to_ctg_connector)[0] for ctg in r1_ctg_refs}
                if len(r1_gnm_refs) == 1:
                    add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
                                                      r1_16s_refs, 'r1', r1_ctg_refs, 'r1',
                                                      marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
                                                      marker_to_ctg_gnm_Key_connector)

            if (len(r2_16s_refs) > 0) and (len(r2_ctg_refs) > 0):
                r2_gnm_refs = {ctg.split(gnm_to_ctg_connector)[0] for ctg in r2_ctg_refs}
                if len(r2_gnm_refs) == 1:
                    add_linkage_to_LinkingRecord_dict(MappingRecord_dict, qualified_read,
                                                      r2_16s_refs, 'r2', r2_ctg_refs, 'r2',
                                                      marker_to_ctg_LinkingRecord_dict, marker_len_dict, ctg_len_dict,
                                                      marker_to_ctg_gnm_Key_connector)

    # filter with link num, short cigar pct and linked to ref end pct
    linkages_QC_txt_handle = open(linkages_QC_txt, 'w')
    linkages_QC_txt_handle.write('Marker\tContig\tLinkages\tLinked_to_16s_end(num pct)\tLinked_to_ctg_end(num pct)\tshort_cigar_pct_16s\tshort_cigar_pct_ctg\tmatched_len_marker\tmatched_len_contig\n')
    marker_to_ctg_linkage_num_dict_after_qc = {}
    for each_link in marker_to_ctg_LinkingRecord_dict:
        linking_reads = marker_to_ctg_LinkingRecord_dict[each_link].linking_reads_base
        if len(linking_reads) >= 3:
            each_link_split = each_link.split(marker_to_ctg_gnm_Key_connector)
            marker_id = each_link_split[0]
            gnm_id = each_link_split[1].split(gnm_to_ctg_connector)[0]
            ctg_id = each_link_split[1].split(gnm_to_ctg_connector)[1]

            linkage_cigar_16s_side_all = marker_to_ctg_LinkingRecord_dict[each_link].linking_cigar_l
            linkage_cigar_ctg_side_all = marker_to_ctg_LinkingRecord_dict[each_link].linking_cigar_r

            # get pct of short aligned cigar
            short_cigar_pct_16s = get_short_cigar_pct(linkage_cigar_16s_side_all, short_M_len_16s)
            short_cigar_pct_ctg = get_short_cigar_pct(linkage_cigar_ctg_side_all, short_M_len_ctg)

            # check num/pct of reads linked to 16s and ctg end
            min_dist_list_to_16s_end = marker_to_ctg_LinkingRecord_dict[each_link].min_dist_to_end_l
            min_dist_list_to_ctg_end = marker_to_ctg_LinkingRecord_dict[each_link].min_dist_to_end_r
            linked_to_16s_end_cigar_num, linked_to_16s_end_cigar_pct = get_linked_to_ref_end_cigar_num_and_pct(min_dist_list_to_16s_end)
            linked_to_ctg_end_cigar_num, linked_to_ctg_end_cigar_pct = get_linked_to_ref_end_cigar_num_and_pct(min_dist_list_to_ctg_end)

            # get all matched positions
            all_matched_pos_l = get_all_matched_pos(marker_to_ctg_LinkingRecord_dict[each_link].linking_cigar_l, marker_to_ctg_LinkingRecord_dict[each_link].linking_pos_l)
            all_matched_pos_r = get_all_matched_pos(marker_to_ctg_LinkingRecord_dict[each_link].linking_cigar_r, marker_to_ctg_LinkingRecord_dict[each_link].linking_pos_r)

            linked_to_16s_end = False
            if (linked_to_16s_end_cigar_num > 0) and (linked_to_16s_end_cigar_pct >= 5):
                linked_to_16s_end = True

            linked_to_ctg_end = False
            if (linked_to_ctg_end_cigar_num > 0) and (linked_to_ctg_end_cigar_pct >= 5):
                linked_to_ctg_end = True

            linkages_QC_txt_handle.write('%s\t%s\t%s(%s %s)\t%s(%s %s)\t%s\t%s\t%sbp\t%sbp\n' % ('\t'.join(each_link.split('___M___')), len(linking_reads),
                                                                  linked_to_16s_end, linked_to_16s_end_cigar_num, float("{0:.2f}".format(linked_to_16s_end_cigar_pct)),
                                                                  linked_to_ctg_end, linked_to_ctg_end_cigar_num, float("{0:.2f}".format(linked_to_ctg_end_cigar_pct)),
                                                                  float("{0:.2f}".format(short_cigar_pct_16s)),
                                                                  float("{0:.2f}".format(short_cigar_pct_ctg)),
                                                                  len(all_matched_pos_l),
                                                                  len(all_matched_pos_r)))

            if (short_cigar_pct_16s < 85) and (short_cigar_pct_ctg < 85):
                if (linked_to_16s_end is True) and (linked_to_ctg_end is True):
                    marker_to_ctg_linkage_num_dict_after_qc[each_link] = len(linking_reads)
                else:
                    if (len(all_matched_pos_l) >= 500) and (len(all_matched_pos_r) >= 500):
                        marker_to_ctg_linkage_num_dict_after_qc[each_link] = len(linking_reads)
    linkages_QC_txt_handle.close()

    # get number of linkages at genome level
    marker_to_gnm_linkage_cigar_dict_16s_side = {}
    marker_to_gnm_linkage_cigar_dict_ctg_side = {}
    marker_to_gnm_link_num = {}
    for each_marker_to_ctg_key in marker_to_ctg_linkage_num_dict_after_qc:
        marker_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[0]
        ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[1]
        gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
        marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
        marker_to_ctg_linkage_cigar_16s_side = marker_to_ctg_LinkingRecord_dict[each_marker_to_ctg_key].linking_cigar_l
        marker_to_ctg_linkage_cigar_ctg_side = marker_to_ctg_LinkingRecord_dict[each_marker_to_ctg_key].linking_cigar_r

        if marker_to_gnm_key not in marker_to_gnm_link_num:
            marker_to_gnm_link_num[marker_to_gnm_key] = marker_to_ctg_linkage_num_dict_after_qc[each_marker_to_ctg_key]
            marker_to_gnm_linkage_cigar_dict_16s_side[marker_to_gnm_key] = marker_to_ctg_linkage_cigar_16s_side
            marker_to_gnm_linkage_cigar_dict_ctg_side[marker_to_gnm_key] = marker_to_ctg_linkage_cigar_ctg_side
        else:
            marker_to_gnm_link_num[marker_to_gnm_key] += marker_to_ctg_linkage_num_dict_after_qc[each_marker_to_ctg_key]
            for each in marker_to_ctg_linkage_cigar_16s_side:
                marker_to_gnm_linkage_cigar_dict_16s_side[marker_to_gnm_key].append(each)
            for each in marker_to_ctg_linkage_cigar_ctg_side:
                marker_to_gnm_linkage_cigar_dict_ctg_side[marker_to_gnm_key].append(each)

    # write out linkages at genome level
    sankey_file_in_handle = open(link_stats_combined, 'w')
    sankey_file_in_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_linkage in marker_to_gnm_link_num:
        sankey_file_in_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
            each_linkage.split(marker_to_ctg_gnm_Key_connector)[0],
            each_linkage.split(marker_to_ctg_gnm_Key_connector)[1],
            marker_to_gnm_link_num[each_linkage]))
    sankey_file_in_handle.close()

    # write out linking reads
    linking_reads_rd1_handle = open(linking_reads_rd1, 'w')
    for each_link in marker_to_ctg_LinkingRecord_dict:
        current_linking_reads = marker_to_ctg_LinkingRecord_dict[each_link].linking_reads_base
        linking_reads_rd1_handle.write('%s\t%s\n' % (each_link, ','.join(current_linking_reads)))
    linking_reads_rd1_handle.close()


    ####################################### filter_linkages_iteratively ########################################

    print('%s %s' % ((datetime.now().strftime(time_format)), 'filter_linkages_iteratively'))

    pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

    # sort file in
    sort_csv_by_col(link_stats_combined, link_stats_combined_sorted, 'Number')

    filter_linkages_iteratively_new2(link_stats_combined_sorted, pairwise_16s_iden_dict, min_iden_16s, marker_len_dict,
                                    min_link_num, within_gnm_linkage_num_diff, link_stats_combined_filtered_s1,
                                    marker_to_gnm_linkage_cigar_dict_16s_side,
                                    marker_to_gnm_linkage_cigar_dict_ctg_side,
                                    marker_to_ctg_gnm_Key_connector)

    #os.system('python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/d7_Oral/assess_linkages_Oral.py')
    print('Done!')
