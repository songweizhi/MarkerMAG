import os
import glob
import shutil
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
from datetime import datetime
import plotly.graph_objects as go
from Bio.SeqRecord import SeqRecord
from distutils.spawn import find_executable


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

        #################### round 2 mini_assembly ####################


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


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


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


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, within_gnm_linkage_num_diff, file_out):

    # get MarkerGene_to_GenomicSeq_dict
    MarkerGene_to_GenomicSeq_dict = {}
    for each_linkage in open(file_in):
        if not each_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            each_linkage_split = each_linkage.strip().split(',')
            MarkerGene_id = each_linkage_split[0][12:]
            GenomicSeq_id = each_linkage_split[1][12:]
            linkage_num = int(each_linkage_split[2])
            if linkage_num > 1:
                if MarkerGene_id not in MarkerGene_to_GenomicSeq_dict:
                    MarkerGene_to_GenomicSeq_dict[MarkerGene_id] = {GenomicSeq_id}
                else:
                    MarkerGene_to_GenomicSeq_dict[MarkerGene_id].add(GenomicSeq_id)

    file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
    file_in_sorted = '%s/%s_sorted%s' % (file_in_path, file_in_basename, file_in_extension)

    # sort file in
    sort_csv_by_col(file_in, file_in_sorted, sort_by_col_header)

    # fileter linkage
    gnm_max_link_num_dict = {}
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    GenomicSeq_best_marker_dict = {}
    for each_match in open(file_in_sorted):
        if each_match.startswith('MarkerGene,GenomicSeq,Number'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            MarkerGene = match_split[0][12:]
            GenomicSeq = match_split[1][12:]
            linkage_num = int(match_split[2])

            current_min_linkage = min_linkages_for_uniq_linked_16s
            if MarkerGene in MarkerGene_to_GenomicSeq_dict:
                if len(MarkerGene_to_GenomicSeq_dict[MarkerGene]) > 1:
                    current_min_linkage = min_linkages

            if linkage_num >= current_min_linkage:
                if MarkerGene not in MarkerGene_with_assignment:
                    if GenomicSeq not in GenomicSeq_best_marker_dict:
                        GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                        gnm_max_link_num_dict[GenomicSeq] = linkage_num
                        file_out_handle.write(each_match)
                        MarkerGene_with_assignment.add(MarkerGene)
                    else:
                        # get identity with best marker
                        current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                        key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))
                        iden_with_best_marker = pairwise_16s_iden_dict.get(key_str, 0)
                        if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                            gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
                            if (linkage_num*100/gnm_max_link_num) >= within_gnm_linkage_num_diff:
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)
                            else:
                                MarkerGene_with_assignment.add(MarkerGene)
    file_out_handle.close()


def filter_linkages_iteratively_mini_assembly_to_ctg(file_in_sorted, min_linkages, file_out):

    # do mini-assemblies assigned to the same mag need to have roughly the same number of linkages? think about this later
    mag_ctg_max_link_num_dict = {}
    mini_assembly_to_mag_dict = {}
    mini_assembly_to_ctg_dict = {}
    file_out_handle = open(file_out, 'w')
    mini_assembly_with_assignment = set()
    mag_ctg_set_with_linked_mini_assembly = set()
    for each_match in open(file_in_sorted):
        if each_match.startswith('MiniAssembly,GenomicSeq,Number'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            mini_assembly = match_split[0]
            mag_ctg_id = match_split[1]
            mag_id = mag_ctg_id.split('___C___')[0]

            linkage_num = int(match_split[2])
            if linkage_num >= min_linkages:
                if mini_assembly not in mini_assembly_with_assignment:
                    if mag_ctg_id not in mag_ctg_max_link_num_dict:
                        mag_ctg_max_link_num_dict[mag_ctg_id] = linkage_num
                        file_out_handle.write(each_match)
                        mini_assembly_to_ctg_dict[mini_assembly] = mag_ctg_id
                        mini_assembly_to_mag_dict[mini_assembly] = mag_id
                        mini_assembly_with_assignment.add(mini_assembly)
                        mag_ctg_set_with_linked_mini_assembly.add(mag_ctg_id)
                    else:
                        ratio_with_best_assignment = linkage_num/(mag_ctg_max_link_num_dict[mag_ctg_id])
                        if ratio_with_best_assignment >= 0.8:
                            file_out_handle.write(each_match)
                            mini_assembly_to_ctg_dict[mini_assembly] = mag_ctg_id
                            mini_assembly_to_mag_dict[mini_assembly] = mag_id
                            mini_assembly_with_assignment.add(mini_assembly)
                            mag_ctg_set_with_linked_mini_assembly.add(mag_ctg_id)
                        else:
                            mini_assembly_with_assignment.add(mini_assembly)
    file_out_handle.close()

    return mini_assembly_to_ctg_dict, mini_assembly_to_mag_dict, mag_ctg_set_with_linked_mini_assembly


def filter_linkages_iteratively_16s_to_mini_assembly(file_in_sorted, min_linkages, file_out):

    # do mini-assemblies assigned to the same mag need to have roughly the same number of linkages? think about this later
    mag_ctg_max_link_num_dict = {}
    mini_assembly_to_mag_dict = {}
    mini_assembly_to_ctg_dict = {}
    file_out_handle = open(file_out, 'w')
    mini_assembly_with_assignment = set()
    mag_ctg_set_with_linked_mini_assembly = set()
    for each_match in open(file_in_sorted):
        if each_match.endswith(',Number\n'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            mini_assembly = match_split[0]
            mag_ctg_id = match_split[1]

            linkage_num = int(match_split[2])
            if linkage_num >= min_linkages:
                if mini_assembly not in mini_assembly_with_assignment:
                    if mag_ctg_id not in mag_ctg_max_link_num_dict:
                        mag_ctg_max_link_num_dict[mag_ctg_id] = linkage_num
                        file_out_handle.write(each_match)
                        mini_assembly_to_ctg_dict[mini_assembly] = mag_ctg_id
                        mini_assembly_with_assignment.add(mini_assembly)
                        mag_ctg_set_with_linked_mini_assembly.add(mag_ctg_id)
                    else:
                        ratio_with_best_assignment = linkage_num / (mag_ctg_max_link_num_dict[mag_ctg_id])
                        if ratio_with_best_assignment >= 0.8:
                            file_out_handle.write(each_match)
                            mini_assembly_to_ctg_dict[mini_assembly] = mag_ctg_id
                            mini_assembly_with_assignment.add(mini_assembly)
                            mag_ctg_set_with_linked_mini_assembly.add(mag_ctg_id)
                        else:
                            mini_assembly_with_assignment.add(mini_assembly)
    file_out_handle.close()

    return mini_assembly_to_ctg_dict


def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc


def linkage_vis_worker(arguments_list):
    reads_file_base = arguments_list[0]
    mafft_seq_folder = arguments_list[1]
    marker_seq = arguments_list[2]
    contig_seq = arguments_list[3]
    end_ctg_len_for_mafft = arguments_list[4]
    gap_N_num = arguments_list[5]
    bowtie_parameter = arguments_list[6]
    marker_pos_list = arguments_list[7]
    contig_pos_list = arguments_list[8]
    pwd_bowtie2_build_exe = arguments_list[9]
    pwd_bowtie2_exe = arguments_list[10]
    marker_seq_name = arguments_list[11]
    contig_seq_name = arguments_list[12]

    marker_len = len(marker_seq)
    contig_len = len(contig_seq)

    # get marker linked end
    marker_pos_median = np.median(marker_pos_list)
    linked_end_marker = 'middle'
    if marker_pos_median <= (marker_len / 3):
        linked_end_marker = 'left'
    elif marker_pos_median >= (marker_len * 2 / 3):
        linked_end_marker = 'right'

    # get contig linked end
    contig_pos_median = np.median(contig_pos_list)
    contig_pos_middle = int(round(float(contig_pos_median)))
    linked_end_contig = 'middle(%s)' % contig_pos_middle
    if contig_pos_median <= 200:
        linked_end_contig = 'left'
    elif (contig_len - contig_pos_median) <= 200:
        linked_end_contig = 'right'

    # get contig end sequence
    if linked_end_contig == 'left':
        contig_seq_for_mafft = contig_seq[:end_ctg_len_for_mafft]
    elif linked_end_contig == 'right':
        contig_seq_for_mafft = contig_seq[-end_ctg_len_for_mafft:]
    else:
        left_end_pos = 0
        if (contig_pos_middle - end_ctg_len_for_mafft) > 0:
            left_end_pos = contig_pos_middle - end_ctg_len_for_mafft
        right_end_pos = contig_pos_middle + end_ctg_len_for_mafft
        if right_end_pos > contig_len:
            right_end_pos = contig_len
        contig_seq_for_mafft = contig_seq[left_end_pos:(right_end_pos - 1)]

    # concatenate 16s and contig sequences
    to_concatenate = False
    concatenated_seq_id = ''
    concatenated_seq = ''
    concatenate_pos = 0
    if linked_end_contig in ['left', 'right']:
        if (linked_end_marker == 'right') and (linked_end_contig == 'left'):
            concatenated_seq_id = '%s_NNN_%s' % (marker_seq_name, contig_seq_name)
            concatenated_seq = '%s%s%s' % (marker_seq, 'N' * gap_N_num, contig_seq_for_mafft)
            to_concatenate = True
            concatenate_pos = len(marker_seq) + round(gap_N_num / 2)

        if (linked_end_marker == 'left') and (linked_end_contig == 'right'):
            concatenated_seq_id = '%s_NNN_%s' % (contig_seq_name, marker_seq_name)
            concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N' * gap_N_num, marker_seq)
            to_concatenate = True
            concatenate_pos = len(contig_seq_for_mafft) + round(gap_N_num / 2)

        if (linked_end_marker == 'left') and (linked_end_contig == 'left'):
            marker_seq_rc = get_rc(marker_seq)
            concatenated_seq_id = '%s_RC_NNN_%s' % (marker_seq_name, contig_seq_name)
            concatenated_seq = '%s%s%s' % (marker_seq_rc, 'N' * gap_N_num, contig_seq_for_mafft)
            to_concatenate = True
            concatenate_pos = len(marker_seq_rc) + round(gap_N_num / 2)

        if (linked_end_marker == 'right') and (linked_end_contig == 'right'):
            marker_seq_rc = get_rc(marker_seq)
            concatenated_seq_id = '%s_NNN_%s_RC' % (contig_seq_name, marker_seq_name)
            concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N' * gap_N_num, marker_seq_rc)
            to_concatenate = True
            concatenate_pos = len(contig_seq_for_mafft) + round(gap_N_num / 2)

    # write out sequences
    pwd_seq_file_cbd = '%s/%s/%s_cbd.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s = '%s/%s/%s_16s.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg = '%s/%s/%s_ctg.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    if to_concatenate is True:
        pwd_seq_file_cbd_handle = open(pwd_seq_file_cbd, 'w')
        pwd_seq_file_cbd_handle.write('>%s\n' % concatenated_seq_id)
        pwd_seq_file_cbd_handle.write('%s\n' % concatenated_seq)
        pwd_seq_file_cbd_handle.close()
    else:
        # write out 16s sequence
        pwd_seq_file_16s_handle = open(pwd_seq_file_16s, 'w')
        pwd_seq_file_16s_handle.write('>%s\n' % marker_seq_name)
        pwd_seq_file_16s_handle.write('%s\n' % marker_seq)
        pwd_seq_file_16s_handle.close()
        # write out ctg sequence
        pwd_seq_file_ctg_handle = open(pwd_seq_file_ctg, 'w')
        pwd_seq_file_ctg_handle.write('>%s\n' % contig_seq_name)
        pwd_seq_file_ctg_handle.write('%s\n' % contig_seq_for_mafft)
        pwd_seq_file_ctg_handle.close()

    ########## mapping ##########

    pwd_seq_file_cbd = '%s/%s/%s_cbd.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s = '%s/%s/%s_16s.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg = '%s/%s/%s_ctg.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_reads_r1 = '%s/%s/%s_R1.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_reads_r2 = '%s/%s/%s_R2.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_index = '%s/%s/%s_cbd' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_index = '%s/%s/%s_16s' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_index = '%s/%s/%s_ctg' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_sam = '%s/%s/%s_cbd.sam' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_sam = '%s/%s/%s_16s.sam' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_sam = '%s/%s/%s_ctg.sam' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_sam_log = '%s/%s/%s_cbd.log' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_sam_log = '%s/%s/%s_16s.log' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_sam_log = '%s/%s/%s_ctg.log' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_Tablet_xml = '%s/%s/%s_cbd.tablet' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_Tablet_xml = '%s/%s/%s_16s.tablet' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_Tablet_xml = '%s/%s/%s_ctg.tablet' % (mafft_seq_folder, reads_file_base, reads_file_base)

    if to_concatenate is True:
        index_ref_cmd = '%s --quiet -f %s %s' % (pwd_bowtie2_build_exe, pwd_seq_file_cbd, pwd_seq_file_cbd_index)
        bowtie2_cmd = '%s -x %s -U %s,%s -S %s -p 1 -f %s 2> %s' % (
            pwd_bowtie2_exe, pwd_seq_file_cbd_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_cbd_sam,
            bowtie_parameter, pwd_seq_file_cbd_sam_log)
        os.system(index_ref_cmd)
        os.system(bowtie2_cmd)

        # write out Tablet xml file
        pwd_seq_file_cbd_Tablet_xml_handle = open(pwd_seq_file_cbd_Tablet_xml, 'w')
        pwd_seq_file_cbd_Tablet_xml_handle.write('<tablet>\n')
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <assembly>%s_cbd.sam</assembly>\n' % reads_file_base)
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <reference>%s_cbd.fa</reference>\n' % reads_file_base)
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <contig>%s</contig>\n' % concatenated_seq_id)
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <position>%s</position>\n' % concatenate_pos)
        pwd_seq_file_cbd_Tablet_xml_handle.write('</tablet>\n')
        pwd_seq_file_cbd_Tablet_xml_handle.close()
    else:
        index_ref_cmd_16s = '%s --quiet -f %s %s' % (pwd_bowtie2_build_exe, pwd_seq_file_16s, pwd_seq_file_16s_index)
        index_ref_cmd_ctg = '%s --quiet -f %s %s' % (pwd_bowtie2_build_exe, pwd_seq_file_ctg, pwd_seq_file_ctg_index)
        os.system(index_ref_cmd_16s)
        os.system(index_ref_cmd_ctg)
        bowtie2_cmd_16s = '%s -x %s -U %s,%s -S %s -p 6 -f %s 2> %s' % (pwd_bowtie2_exe, pwd_seq_file_16s_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_16s_sam,
            bowtie_parameter, pwd_seq_file_16s_sam_log)
        bowtie2_cmd_ctg = '%s -x %s -U %s,%s -S %s -p 6 -f %s 2> %s' % (pwd_bowtie2_exe, pwd_seq_file_ctg_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_ctg_sam,
            bowtie_parameter, pwd_seq_file_ctg_sam_log)
        os.system(bowtie2_cmd_16s)
        os.system(bowtie2_cmd_ctg)

        # write out Tablet xml file
        pwd_seq_file_16s_Tablet_xml_handle = open(pwd_seq_file_16s_Tablet_xml, 'w')
        pwd_seq_file_16s_Tablet_xml_handle.write('<tablet>\n')
        pwd_seq_file_16s_Tablet_xml_handle.write('        <assembly>%s_16s.sam</assembly>\n' % reads_file_base)
        pwd_seq_file_16s_Tablet_xml_handle.write('        <reference>%s_16s.fa</reference>\n' % reads_file_base)
        pwd_seq_file_16s_Tablet_xml_handle.write('        <contig>Marker</contig>\n')
        pwd_seq_file_16s_Tablet_xml_handle.write('</tablet>\n')
        pwd_seq_file_16s_Tablet_xml_handle.close()

        pwd_seq_file_ctg_Tablet_xml_handle = open(pwd_seq_file_ctg_Tablet_xml, 'w')
        pwd_seq_file_ctg_Tablet_xml_handle.write('<tablet>\n')
        pwd_seq_file_ctg_Tablet_xml_handle.write('        <assembly>%s_ctg.sam</assembly>\n' % reads_file_base)
        pwd_seq_file_ctg_Tablet_xml_handle.write('        <reference>%s_ctg.fa</reference>\n' % reads_file_base)
        pwd_seq_file_ctg_Tablet_xml_handle.write('        <contig>Contig</contig>\n')
        pwd_seq_file_ctg_Tablet_xml_handle.write('</tablet>\n')
        pwd_seq_file_ctg_Tablet_xml_handle.close()

    # remove tmp files
    os.system('rm %s/%s/%s*.bt2' % (mafft_seq_folder, reads_file_base, reads_file_base))


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


def get_qualified_ref_dict(ref_dict, ref_len_dict, ref_min_mismatch, min_M_len, mismatch_cutoff, refs_to_ignore):

    refs_passed_qc = {}
    refs_passed_qc_with_pos = {}

    for each_ref in ref_dict:

        # one read can only mapped to on ref once
        matched_pos_dict = ref_dict[each_ref]
        if len(matched_pos_dict) > 1:
            refs_to_ignore.add(each_ref)
        else:
            ref_pos = list(matched_pos_dict.keys())[0]
            ref_cigar = matched_pos_dict[ref_pos]
            qualified_cigar = check_cigar_quality(ref_cigar, True, ref_min_mismatch, mismatch_cutoff, min_M_len, ref_pos, ref_len_dict[each_ref])

            if qualified_cigar is False:
                refs_to_ignore.add(each_ref)
            else:
                refs_passed_qc[each_ref] = [ref_cigar]
                refs_passed_qc_with_pos[each_ref] = {ref_pos: ref_cigar}

    return refs_passed_qc, refs_passed_qc_with_pos


def get_short_cigar_pct(cigar_list, short_M_len):
    short_cigar_num = 0
    for each_cigar in cigar_list:
        aln_len, aln_pct, clp_len, clp_pct, mis_pct = get_cigar_stats(cigar_splitter(each_cigar))
        if aln_len <= short_M_len:
            short_cigar_num += 1
    short_cigar_pct = short_cigar_num * 100 / len(cigar_list)
    short_cigar_pct = float("{0:.2f}".format(short_cigar_pct))
    return short_cigar_pct


def get_min_dist_to_ref_end(cigar_str, cigar_pos, ref_len):
    cigar_aln_len = get_cigar_aln_len(cigar_splitter(cigar_str))
    cigar_dist_to_left = cigar_pos - 1
    cigar_dist_to_right = ref_len - cigar_pos - cigar_aln_len + 1
    min_dist_to_mini_end = min(cigar_dist_to_left, cigar_dist_to_right)
    return min_dist_to_mini_end


def get_all_matched_pos(linking_cigar_list, linking_pos_list):

    all_matched_pos = set()
    for (each_linking_cigar, each_linking_pos) in zip(linking_cigar_list, linking_pos_list):
        current_cigar_aligned_len = get_cigar_aln_len(cigar_splitter(each_linking_cigar))
        matched_seq_end = each_linking_pos + current_cigar_aligned_len
        pos_list = list(range(each_linking_pos, matched_seq_end))
        all_matched_pos.update(pos_list)
    return all_matched_pos


def get_matched_pos_num_pct(linking_cigar_list, linking_pos_list, ref_len):

    matched_pos = set()
    for (each_linking_cigar, each_linking_pos) in zip(linking_cigar_list, linking_pos_list):
        current_cigar_aligned_len = get_cigar_aln_len(cigar_splitter(each_linking_cigar))
        matched_seq_end = each_linking_pos + current_cigar_aligned_len
        pos_list = list(range(each_linking_pos, matched_seq_end))
        matched_pos.update(pos_list)

    matched_pos_num = len(matched_pos)
    matched_pos_pct = matched_pos_num * 100 / ref_len
    matched_pos_pct = float("{0:.2f}".format(matched_pos_pct))

    return matched_pos_num, matched_pos_pct


####################################################### file in ########################################################

if __name__ == '__main__':

    # step_2_wd                                       = '/Users/songweizhi/Desktop/tunning_rd2_Oral'
    # blast_results_all_vs_all_16s                    = '%s/file_in/Oral_0717_45_45_min1200_16S_all_vs_all_blastn.tab'            % step_2_wd
    # splitted_sam_mp_file_folder                     = '%s/file_in/Oral_0717_45_45_min1200_input_reads_to_16S_MappingRecord'     % step_2_wd
    # marker_gene_seqs                                = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_min1200.QC.fa'               % step_2_wd
    # mini_assemblies                                 = '%s/file_in/scaffolds.fasta'                                              % step_2_wd
    # sam_file_mini_assembly_reformatted              = '%s/file_in/scaffolds_bowtie_16s.sam'                                     % step_2_wd
    # combined_1st_round_unlinked_mag_end_seq         = '%s/file_in/round_1_unlinked_gnm_end_500bp.fa'                            % step_2_wd
    # rd2_to_extract_flking_both_r1_fa                = '%s/file_in/rd2_to_extract_flking_both_r1.fa'                             % step_2_wd
    # rd2_to_extract_flking_both_r2_fa                = '%s/file_in/rd2_to_extract_flking_both_r2.fa'                             % step_2_wd
    # sam_file_mini_assembly_reformatted_sorted       = '%s/file_in/scaffolds_bowtie_sorted.sam'                                  % step_2_wd
    # free_living_ctg_ref_file_with_pos_cigar         = '%s/file_in/round2_free_living_ctg_refs_with_pos_cigar.txt'               % step_2_wd
    # free_living_16s_ref_file_no_linked              = '%s/file_in/round2_free_living_16s_refs_no_linked.txt'                        % step_2_wd
    # free_living_ctg_ref_file_no_linked              = '%s/file_in/round2_free_living_ctg_refs_no_linked.txt'                        % step_2_wd

    step_2_wd                                       = '/Users/songweizhi/Desktop/tunning_rd2_MBARC26'
    marker_gene_seqs                                = '%s/file_in/MBARC26_SILVA138_polished.QC.fasta'                                               % step_2_wd
    blast_results_all_vs_all_16s                    = '%s/file_in/MBARC26_0719_60_60_min1200_diff10_mismatch0_3_16S_all_vs_all_blastn.tab'          % step_2_wd
    splitted_sam_mp_file_folder                     = '%s/file_in/MBARC26_0719_60_60_min1200_diff10_mismatch0_3_input_reads_to_16S_MappingRecord'   % step_2_wd
    link_stats_combined_filtered_s1                 = '%s/file_in/MBARC26_0719_60_60_min1200_diff10_mismatch0_3_stats_combined_filtered.txt'        % step_2_wd
    mini_assemblies                                 = '%s/file_in/scaffolds.fasta'                                                                  % step_2_wd
    combined_1st_round_unlinked_mag_end_seq         = '%s/file_in/round_1_unlinked_gnm_end_500bp.fa'                                                % step_2_wd
    rd2_to_extract_flking_both_r1_fa                = '%s/file_in/rd2_to_extract_flking_both_r1.fa'                                                 % step_2_wd
    rd2_to_extract_flking_both_r2_fa                = '%s/file_in/rd2_to_extract_flking_both_r2.fa'                                                 % step_2_wd
    sam_file_mini_assembly_reformatted_sorted       = '%s/file_in/scaffolds_bowtie_sorted.sam'                                                      % step_2_wd
    free_living_ctg_ref_file_with_pos_cigar         = '%s/file_in/round2_free_living_ctg_refs_with_pos_cigar.txt'                                   % step_2_wd
    free_living_16s_ref_file_no_linked              = '%s/file_in/round2_free_living_16s_refs_no_linked.txt'                                        % step_2_wd
    free_living_ctg_ref_file_no_linked              = '%s/file_in/round2_free_living_ctg_refs_no_linked.txt'                                        % step_2_wd

    within_gnm_linkage_num_diff                     = 10
    max_short_cigar_pct_cutoff                      = 75
    linked_to_ctg_end_cigar_num_cutoff              = 1
    linked_to_ctg_end_cigar_pct_cutoff              = 20
    min_iden_16s                                    = 98.5
    min_M_len_mini                                  = 75
    mismatch_cutoff                                 = 0
    ctg_level_min_link                              = 3
    min_link_num                                    = 8
    max_mini_assembly_link_num_diff_between_ctg_16s = 10
    min_aln_16s                                     = 500
    min_cov_16s                                     = 30
    marker_to_ctg_gnm_Key_connector                 = '___M___'
    gnm_to_ctg_connector                            = '___C___'
    mini_assembly_to_16s_ctg_connector              = '___Mini___'
    read_to_marker_connector                        = '___r___'
    gap_N_num                                       = 50
    end_ctg_len_for_mafft                           = 1000
    num_threads                                     = 4
    rd2_with_both_mates                             = False
    short_M_len                                     = 75

    pwd_bowtie2_build_exe                           = '/Users/songweizhi/Softwares/bowtie2/bowtie2-build'
    pwd_bowtie2_exe                                 = '/Users/songweizhi/Softwares/bowtie2/bowtie2'
    bowtie_parameter                                = '--local --all --no-unal -N 1 -L 30'
    seqtk_exe                                       = 'seqtk'

    ##################################################### file out #####################################################

    stats_mini_assembly_to_ctg_qc_report            = '%s/stats_mini_assembly_to_ctg_qc_report.txt'                         % step_2_wd
    stats_mini_assembly_to_16s_qc_report            = '%s/stats_mini_assembly_to_16s_qc_report.txt'                         % step_2_wd
    stats_mini_assembly_to_ctg                      = '%s/stats_mini_assembly_to_ctg.txt'                                   % step_2_wd
    stats_mini_assembly_to_ctg_sorted               = '%s/stats_mini_assembly_to_ctg_sorted.txt'                            % step_2_wd
    stats_mini_assembly_to_ctg_filtered             = '%s/stats_mini_assembly_to_ctg_filtered.txt'                          % step_2_wd

    # stats_mini_assembly_to_16s                      = '%s/stats_mini_assembly_to_16s.txt'                                   % step_2_wd
    # stats_mini_assembly_to_16s_sorted               = '%s/stats_mini_assembly_to_16s_sorted.txt'                            % step_2_wd
    # stats_mini_assembly_to_16s_filtered             = '%s/stats_mini_assembly_to_16s_filtered.txt'                          % step_2_wd

    stats_GapFilling_ctg                            = '%s/stats_GapFilling_ctg.txt'                                         % step_2_wd
    stats_GapFilling_file                           = '%s/stats_GapFilling_gnm.txt'                                         % step_2_wd
    stats_GapFilling_file_filtered                  = '%s/stats_GapFilling_gnm_filtered.txt'                                % step_2_wd

    mafft_seq_folder_mini_to_16s                    = '%s/vis_folder_mini_to_16S'                                           % step_2_wd
    mafft_seq_folder_mini_to_ctg                    = '%s/vis_folder_mini_to_ctg'                                           % step_2_wd
    linking_reads_tab_rd2                           = '%s/linking_reads_rd2.txt'                                            % step_2_wd
    linking_reads_r1_txt                            = '%s/linking_reads_r1.txt'                                             % step_2_wd
    linking_reads_r2_txt                            = '%s/linking_reads_r2.txt'                                             % step_2_wd
    linking_reads_r1_fasta                          = '%s/linking_reads_r1.fa'                                              % step_2_wd
    linking_reads_r2_fasta                          = '%s/linking_reads_r2.fa'                                              % step_2_wd
    vis_folder_rd2                                  = '%s/vis_folder_rd2'                                                   % step_2_wd
    sam_file_mini_assembly                          = '%s/scaffolds_bowtie.sam'                                             % step_2_wd
    sam_file_mini_assembly_log                      = '%s/scaffolds_bowtie.log'                                             % step_2_wd
    rd2_read_extracted_flanking_both_r12_seq        = '%s/rd2_read_to_extract_flanking_both_R12.fa'                         % step_2_wd
    rd2_read_extracted_flanking_both_r12_up_seq     = '%s/rd2_read_to_extract_flanking_both_R12_up.fa'                      % step_2_wd
    rd2_mini_to_16s_short_M_pct                     = '%s/rd2_mini_to_16s_short_M_pct.txt'                                  % step_2_wd

    ####################################################################################################################

    pairwise_16s_iden_dict      = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

    print('Round 2: read sequences into dict')

    # read sequence of flanking reads into dict
    flk_both_read_seq_dict = {}
    for linking_r1 in SeqIO.parse(rd2_to_extract_flking_both_r1_fa, 'fasta'):
        flk_both_read_seq_dict[linking_r1.id] = str(linking_r1.seq)
    for linking_r2 in SeqIO.parse(rd2_to_extract_flking_both_r2_fa, 'fasta'):
        flk_both_read_seq_dict[linking_r2.id] = str(linking_r2.seq)

    # read sequence of mini_assembly into dict
    mini_assembly_seq_dict = {}
    for linked_ctg in SeqIO.parse(mini_assemblies, 'fasta'):
        mini_assembly_seq_dict[linked_ctg.id] = str(linked_ctg.seq)

    # read sequence of unlinked mag end seq into dict
    unlinked_mag_end_seq_dict = {}
    for linked_ctg in SeqIO.parse(combined_1st_round_unlinked_mag_end_seq, 'fasta'):
        unlinked_mag_end_seq_dict[linked_ctg.id] = str(linked_ctg.seq)

    # read marker sequences into dict
    marker_seq_dict = {}
    for each_16s in SeqIO.parse(marker_gene_seqs, 'fasta'):
        marker_seq_dict[each_16s.id] = str(each_16s.seq)

    ##################################################### get reads flanking 16s and ctg ####################################################

    flanking_16s_reads = set()
    flanking_16s_reads_base = set()
    for free_living_read_16s in open(free_living_16s_ref_file_no_linked):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_16s_id = free_living_read_16s_split[0]
            read_16s_id_base = '.'.join(read_16s_id.split('.')[:-1])
            flanking_16s_reads.add(read_16s_id)
            flanking_16s_reads_base.add(read_16s_id_base)

    flanking_ctg_reads = set()
    flanking_ctg_reads_base = set()
    for each_read_to_ctg_ref in open(free_living_ctg_ref_file_no_linked):
        read_id = each_read_to_ctg_ref.strip().split('\t')[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        flanking_ctg_reads.add(read_id)
        flanking_ctg_reads_base.add(read_id_base)

    ##################################################### combine reads flanking 16s and contig ####################################################

    rd2_read_to_extract_flanking_both = set.union(flanking_16s_reads, flanking_ctg_reads)
    rd2_read_to_extract_flanking_both_base = set.union(flanking_16s_reads_base, flanking_ctg_reads_base)

    # this step produce the input file for spades
    rd2_read_extracted_flanking_both_r12_up_seq_handle = open(rd2_read_extracted_flanking_both_r12_up_seq, 'w')
    for each_flk_read in flk_both_read_seq_dict:
        if each_flk_read in rd2_read_to_extract_flanking_both:
            rd2_read_extracted_flanking_both_r12_up_seq_handle.write('>%s\n' % each_flk_read)
            rd2_read_extracted_flanking_both_r12_up_seq_handle.write('%s\n' % flk_both_read_seq_dict[each_flk_read])
    rd2_read_extracted_flanking_both_r12_up_seq_handle.close()

    # file produced here will be mapped to mini-assemblies
    rd2_read_extracted_flanking_both_r12_seq_handle = open(rd2_read_extracted_flanking_both_r12_seq, 'w')
    for flanking_both_base in rd2_read_to_extract_flanking_both_base:
        flanking_both_r1 = '%s.1' % flanking_both_base
        flanking_both_r2 = '%s.2' % flanking_both_base
        flanking_both_r1_seq = flk_both_read_seq_dict[flanking_both_r1]
        flanking_both_r2_seq = flk_both_read_seq_dict[flanking_both_r2]
        rd2_read_extracted_flanking_both_r12_seq_handle.write('>%s\n' % flanking_both_r1)
        rd2_read_extracted_flanking_both_r12_seq_handle.write('%s\n'  % flanking_both_r1_seq)
        rd2_read_extracted_flanking_both_r12_seq_handle.write('>%s\n' % flanking_both_r2)
        rd2_read_extracted_flanking_both_r12_seq_handle.write('%s\n'  % flanking_both_r2_seq)
    rd2_read_extracted_flanking_both_r12_seq_handle.close()

    ######################################### second round linking by assembly #########################################

    # running spades, skipped here

    ##################################################### read in mini-assembly sam file ####################################################

    print('Round 2: read in mini-assembly sam file')

    mini_assembly_len_dict = {}
    MappingRecord_dict_mini = {}
    current_read_base = ''
    current_read_base_r1_mini_ref_dict = dict()
    current_read_base_r2_mini_ref_dict = dict()
    with open(sam_file_mini_assembly_reformatted_sorted) as sam_file_mini_assembly_reformatted_sorted_opened:
        for each_line in sam_file_mini_assembly_reformatted_sorted_opened:
            each_line_split = each_line.strip().split('\t')
            if each_line.startswith('@'):
                mini_assembly_id = ''
                mini_assembly_len = 0
                for each_element in each_line_split:
                    if each_element.startswith('SN:'):
                        mini_assembly_id = each_element[3:]
                    if each_element.startswith('LN:'):
                        mini_assembly_len = int(each_element[3:])
                mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
            else:
                cigar = each_line_split[5]
                read_id = each_line_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_line_split[2]
                ref_pos = int(each_line_split[3])

                if current_read_base == '':
                    current_read_base = read_id_base

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_mini_ref_dict[ref_id] = {ref_pos: cigar}

                elif read_id_base == current_read_base:

                    if cigar != '*':
                        if read_strand == '1':
                            if ref_id not in current_read_base_r1_mini_ref_dict:
                                current_read_base_r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r1_mini_ref_dict[ref_id][ref_pos] = cigar
                        if read_strand == '2':
                            if ref_id not in current_read_base_r2_mini_ref_dict:
                                current_read_base_r2_mini_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r2_mini_ref_dict[ref_id][ref_pos] = cigar
                else:
                    ################################### analysis previous read refs ####################################

                    mini_refs_to_ignore = set()

                    ########## get lowest mismatch for r1/r2 mini refs ##########

                    # get r1_ref_cigar_set
                    r1_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r1_mini_ref_dict.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r1_ref_cigar_set.update(each_pos_dict_values)

                    # get r2_ref_cigar_set
                    r2_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r2_mini_ref_dict.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r2_ref_cigar_set.update(each_pos_dict_values)

                    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_mini)
                    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_mini)

                    ########## filter mini refs for r1 and r2 ##########

                    r1_mini_refs_passed_qc, r1_mini_refs_passed_qc_with_pos = get_qualified_ref_dict(current_read_base_r1_mini_ref_dict, mini_assembly_len_dict, r1_ref_min_mismatch, min_M_len_mini, mismatch_cutoff, mini_refs_to_ignore)
                    r2_mini_refs_passed_qc, r2_mini_refs_passed_qc_with_pos = get_qualified_ref_dict(current_read_base_r2_mini_ref_dict, mini_assembly_len_dict, r2_ref_min_mismatch, min_M_len_mini, mismatch_cutoff, mini_refs_to_ignore)

                    ####################################################################################################

                    r1_mini_refs_no_ignored = {key: value for key, value in r1_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
                    r2_mini_refs_no_ignored = {key: value for key, value in r2_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}

                    # no mate has no_ignored alignments
                    if (len(r1_mini_refs_no_ignored) == 0) and (len(r2_mini_refs_no_ignored) == 0):
                        pass

                    # only r1 has no_ignored alignments
                    elif (len(r1_mini_refs_no_ignored) > 0) and (len(r2_mini_refs_no_ignored) == 0):
                        if current_read_base not in MappingRecord_dict_mini:
                            MappingRecord_dict_mini[current_read_base] = MappingRecord()
                        MappingRecord_dict_mini[current_read_base].r1_mini_ref_dict = current_read_base_r1_mini_ref_dict
                        MappingRecord_dict_mini[current_read_base].r1_mini_refs_no_ignored = r1_mini_refs_no_ignored

                    # only r2 has no_ignored alignments
                    elif (len(r1_mini_refs_no_ignored) == 0) and (len(r2_mini_refs_no_ignored) > 0):
                        if current_read_base not in MappingRecord_dict_mini:
                            MappingRecord_dict_mini[current_read_base] = MappingRecord()
                        MappingRecord_dict_mini[current_read_base].r2_mini_ref_dict = current_read_base_r2_mini_ref_dict
                        MappingRecord_dict_mini[current_read_base].r2_mini_refs_no_ignored = r2_mini_refs_no_ignored

                    # both r1 and r2 have no_ignored alignments
                    else:
                        shared_mini_refs_no_ignored = {key: [r1_mini_refs_no_ignored[key][0], r2_mini_refs_no_ignored[key][0]] for key in
                            set(r1_mini_refs_no_ignored).intersection(set(r2_mini_refs_no_ignored))}
                        if len(shared_mini_refs_no_ignored) > 0:
                            if current_read_base not in MappingRecord_dict_mini:
                                MappingRecord_dict_mini[current_read_base] = MappingRecord()
                            MappingRecord_dict_mini[current_read_base].r1_mini_ref_dict = current_read_base_r1_mini_ref_dict
                            MappingRecord_dict_mini[current_read_base].r2_mini_ref_dict = current_read_base_r2_mini_ref_dict
                            MappingRecord_dict_mini[current_read_base].shared_mini_refs_no_ignored = shared_mini_refs_no_ignored

                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_mini_ref_dict = dict()
                    current_read_base_r2_mini_ref_dict = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_mini_ref_dict[ref_id] = {ref_pos: cigar}

    ####################################################################################################################
    #################################################### rd2 linking ###################################################
    ####################################################################################################################

    splitted_sam_mp_file_re  = '%s/*.txt' % splitted_sam_mp_file_folder
    splitted_sam_mp_file_set = glob.glob(splitted_sam_mp_file_re)

    print('Round 2: read in %s 16S MP text file' % len(splitted_sam_mp_file_set))

    MappingRecord_dict = {}
    for each_mp_file in splitted_sam_mp_file_set:
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


    ############################################ link mini-assembly to MAGs ############################################

    # read ctg side aln_pos and cigar into dict
    mini_ctg_side_pos_dict = dict()
    mini_ctg_side_cigar_dict = dict()
    for each_read_to_ctg_ref in open(free_living_ctg_ref_file_with_pos_cigar):
        each_read_to_ctg_ref_split = each_read_to_ctg_ref.strip().split('\t')
        read_id = each_read_to_ctg_ref_split[0]
        ctg_refs = each_read_to_ctg_ref_split[1].split(',')
        for each_ctg_ref in ctg_refs:
            each_ctg_ref_split = each_ctg_ref.split('__pc__')
            ctg_ref_id = each_ctg_ref_split[0]
            ctg_ref_pos = int(each_ctg_ref_split[1])
            ctg_ref_cigar = each_ctg_ref_split[2]
            read_to_ctg_ref_key = '%s__ctg__%s' % (read_id, ctg_ref_id)
            mini_ctg_side_pos_dict[read_to_ctg_ref_key] = ctg_ref_pos
            mini_ctg_side_cigar_dict[read_to_ctg_ref_key] = ctg_ref_cigar

    mini_to_mag_LinkingRecord_dict = {}
    for free_living_read_ctg in open(free_living_ctg_ref_file_no_linked):
        free_living_read_ctg_split = free_living_read_ctg.strip().split('\t')
        read_id = free_living_read_ctg_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        read_ctg_refs = free_living_read_ctg_split[1].split(',')
        read_mini_mp = MappingRecord_dict_mini.get(read_id_base, None)
        if read_mini_mp is not None:

            read_mini_refs = set()
            for mini_ref in read_mini_mp.shared_mini_refs_no_ignored:
                read_mini_refs.add(mini_ref)
            if read_strand == '1':
                for mini_ref in read_mini_mp.r1_mini_refs_no_ignored:
                    read_mini_refs.add(mini_ref)
            if read_strand == '2':
                for mini_ref in read_mini_mp.r2_mini_refs_no_ignored:
                    read_mini_refs.add(mini_ref)

            for each_mini_ref in read_mini_refs:
                for each_ctg_ref in read_ctg_refs:
                    mini_to_ctg_key_with_l_r = '%s%s%s' % (each_mini_ref, mini_assembly_to_16s_ctg_connector, each_ctg_ref)
                    mini_ref_len = mini_assembly_len_dict[each_mini_ref]
                    ctg_ref_len = len(unlinked_mag_end_seq_dict[each_ctg_ref])
                    if mini_to_ctg_key_with_l_r not in mini_to_mag_LinkingRecord_dict:
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r] = LinkingRecord()
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_l = each_mini_ref
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_r = each_ctg_ref
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_len_l = mini_ref_len
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_len_r = ctg_ref_len

                    # get mini side cigar here
                    current_mini_mp = MappingRecord_dict_mini[read_id_base]
                    current_linking_cigar = ''
                    current_linking_pos = ''
                    if read_strand == '1':
                        current_linking_cigar = list(current_mini_mp.r1_mini_ref_dict[each_mini_ref].values())[0]
                        current_linking_pos = list(current_mini_mp.r1_mini_ref_dict[each_mini_ref].keys())[0]
                    if read_strand == '2':
                        current_linking_cigar = list(current_mini_mp.r2_mini_ref_dict[each_mini_ref].values())[0]
                        current_linking_pos = list(current_mini_mp.r2_mini_ref_dict[each_mini_ref].keys())[0]

                    # get min dist to mini end
                    current_linking_cigar_aln_len = get_cigar_aln_len(cigar_splitter(current_linking_cigar))
                    current_linking_cigar_dist_to_left = current_linking_pos - 1
                    current_linking_cigar_dist_to_right = mini_ref_len - current_linking_pos - current_linking_cigar_aln_len + 1
                    min_dist_to_mini_end = min(current_linking_cigar_dist_to_left, current_linking_cigar_dist_to_right)

                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_reads_base.append(read_id_base)
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_cigar_l.append(current_linking_cigar)
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_pos_l.append(current_linking_pos)
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].min_dist_to_end_l.append(min_dist_to_mini_end)

                    # get ctg side cigar here
                    read_to_ctg_ref_key = '%s__ctg__%s' % (read_id, each_ctg_ref)
                    current_linking_cigar_ctg_side = mini_ctg_side_cigar_dict.get(read_to_ctg_ref_key, None)
                    current_linking_pos_ctg_side = mini_ctg_side_pos_dict.get(read_to_ctg_ref_key, None)
                    if current_linking_cigar_ctg_side is not None:
                        # get min dist to ctg end
                        current_linking_cigar_ctg_side_aln_len = get_cigar_aln_len(cigar_splitter(current_linking_cigar_ctg_side))
                        current_linking_cigar_ctg_side_dist_to_left = current_linking_pos_ctg_side - 1
                        current_linking_cigar_ctg_side_dist_to_right = ctg_ref_len - current_linking_pos_ctg_side - current_linking_cigar_ctg_side_aln_len + 1
                        min_dist_to_mini_end_ctg_side = min(current_linking_cigar_ctg_side_dist_to_left, current_linking_cigar_ctg_side_dist_to_right)

                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_cigar_r.append(current_linking_cigar_ctg_side)
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_pos_r.append(current_linking_pos_ctg_side)
                        mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].min_dist_to_end_r.append(min_dist_to_mini_end_ctg_side)

    all_linking_reads_base_set_rd2_mini_to_ctg = set()
    stats_mini_assembly_to_ctg_qc_report_handle = open(stats_mini_assembly_to_ctg_qc_report, 'w')
    stats_mini_assembly_to_ctg_qc_report_handle.write('Mini\tContig\tLinkages\tLinked_to_mini_end(num pct)\tLinked_to_ctg_end(num pct)\tshort_cigar_pct_mini\tshort_cigar_pct_ctg\n')
    stats_mini_assembly_to_ctg_handle = open(stats_mini_assembly_to_ctg, 'w')
    stats_mini_assembly_to_ctg_handle.write('MiniAssembly,GenomicSeq,Number\n')
    mini_assembly_to_ctg_dict_with_l_r_read_base_min3 = {}
    for each_link in mini_to_mag_LinkingRecord_dict.copy():
        current_linkage_linking_read_base = mini_to_mag_LinkingRecord_dict[each_link].linking_reads_base
        if len(current_linkage_linking_read_base) < 3:
            mini_to_mag_LinkingRecord_dict.pop(each_link)
        else:
            id_mini_assembly = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_l
            id_ctg = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_r
            current_linkage_cigar_mini_side = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_l
            current_linkage_cigar_ctg_side = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_r

            # get pct of short aligned cigar
            short_cigar_pct_mini = get_short_cigar_pct(current_linkage_cigar_mini_side, short_M_len)
            short_cigar_pct_ctg = get_short_cigar_pct(current_linkage_cigar_ctg_side, short_M_len)

            # check num/pct of reads linked to mini end
            min_dist_list_to_mini_end  = mini_to_mag_LinkingRecord_dict[each_link].min_dist_to_end_l
            linked_to_mini_end_cigar_num = 0
            for each_min_dist in min_dist_list_to_mini_end:
                if each_min_dist <= 10:
                    linked_to_mini_end_cigar_num += 1
            linked_to_mini_end_cigar_pct = linked_to_mini_end_cigar_num * 100 / len(min_dist_list_to_mini_end)
            linked_to_mini_end_cigar_pct = float("{0:.2f}".format(linked_to_mini_end_cigar_pct))

            # check num/pct of reads linked to ctg end
            min_dist_list_to_ctg_end = mini_to_mag_LinkingRecord_dict[each_link].min_dist_to_end_r
            linked_to_ctg_end_cigar_num = 0
            for each_min_dist in min_dist_list_to_ctg_end:
                if each_min_dist <= 10:
                    linked_to_ctg_end_cigar_num += 1
            linked_to_ctg_end_cigar_pct = linked_to_ctg_end_cigar_num * 100 / len(min_dist_list_to_ctg_end)
            linked_to_ctg_end_cigar_pct = float("{0:.2f}".format(linked_to_ctg_end_cigar_pct))

            linked_to_mini_end = False
            if (linked_to_mini_end_cigar_num > 0) and (linked_to_mini_end_cigar_pct >= 5):
                linked_to_mini_end = True

            linked_to_ctg_end = False
            if (linked_to_ctg_end_cigar_num >= linked_to_ctg_end_cigar_num_cutoff) and (linked_to_ctg_end_cigar_pct >= linked_to_ctg_end_cigar_pct_cutoff):
                linked_to_ctg_end = True

            short_cigar_pct_mini_bool = False
            if short_cigar_pct_mini < max_short_cigar_pct_cutoff:
                short_cigar_pct_mini_bool = True

            short_cigar_pct_ctg_bool = False
            if short_cigar_pct_ctg < max_short_cigar_pct_cutoff:
                short_cigar_pct_ctg_bool = True

            # get bp and pct of matched positions
            linking_cigar_list_mini = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_l
            linking_pos_list_mini   = mini_to_mag_LinkingRecord_dict[each_link].linking_pos_l
            ref_len_mini            = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_len_l
            matched_pos_mini_num, matched_pos_mini_pct = get_matched_pos_num_pct(linking_cigar_list_mini, linking_pos_list_mini, ref_len_mini)

            linking_cigar_list_ctg  = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_r
            linking_pos_list_ctg    = mini_to_mag_LinkingRecord_dict[each_link].linking_pos_r
            ref_len_ctg             = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_len_r
            matched_pos_ctg_num, matched_pos_ctg_pct   = get_matched_pos_num_pct(linking_cigar_list_ctg, linking_pos_list_ctg, ref_len_ctg)

            matched_region_passed_qc_mini = False
            if (matched_pos_mini_num >= 300) or (matched_pos_mini_pct >= 50):
                matched_region_passed_qc_mini = True

            matched_region_passed_qc_ctg = False
            if matched_pos_ctg_num >= 300:
                matched_region_passed_qc_ctg = True

            # write out qc
            stats_mini_assembly_to_ctg_qc_report_handle.write('%s\t%s\t%s\t%s(%s %s)__v__%s(%sbp %s)\t%s(%s %s)__v__%s(%sbp)\t%s(%s)\t%s(%s)\n' % (id_mini_assembly, id_ctg, len(current_linkage_linking_read_base),
                                                                  linked_to_mini_end, linked_to_mini_end_cigar_num, linked_to_mini_end_cigar_pct,
                                                                  matched_region_passed_qc_mini, matched_pos_mini_num, matched_pos_mini_pct,
                                                                  linked_to_ctg_end, linked_to_ctg_end_cigar_num, linked_to_ctg_end_cigar_pct,
                                                                  matched_region_passed_qc_ctg, matched_pos_ctg_num,
                                                                  short_cigar_pct_mini_bool, short_cigar_pct_mini,
                                                                  short_cigar_pct_ctg_bool, short_cigar_pct_ctg))

            if (short_cigar_pct_mini_bool is True) and (short_cigar_pct_ctg_bool is True):
                if ((linked_to_mini_end is True) or (matched_region_passed_qc_mini is True)) and ((linked_to_ctg_end is True)):
                    mini_assembly_to_ctg_dict_with_l_r_read_base_min3[each_link] = current_linkage_linking_read_base
                    all_linking_reads_base_set_rd2_mini_to_ctg.update(current_linkage_linking_read_base)
                    stats_mini_assembly_to_ctg_handle.write('%s,%s,%s\n' % (id_mini_assembly, id_ctg, len(current_linkage_linking_read_base)))
    stats_mini_assembly_to_ctg_handle.close()
    stats_mini_assembly_to_ctg_qc_report_handle.close()

    # sort and filter
    sort_csv_by_col(stats_mini_assembly_to_ctg, stats_mini_assembly_to_ctg_sorted, 'Number')
    mini_assembly_to_ctg_dict, mini_assembly_to_mag_dict, mag_ctg_set_with_linked_mini_assembly = filter_linkages_iteratively_mini_assembly_to_ctg(stats_mini_assembly_to_ctg_sorted, 3,
                                                                                 stats_mini_assembly_to_ctg_filtered)

    ################################## link 16S to mini-assemblies with MAG assignment #################################

    mini_to_16s_LinkingRecord_dict = {}  # left: mini, right: 16S
    for free_living_read_16s in open(free_living_16s_ref_file_no_linked):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_id = free_living_read_16s_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            read_16s_refs = free_living_read_16s_split[1].split(',')
            read_mini_mp = MappingRecord_dict_mini.get(read_id_base, None)
            read_16s_mp  = MappingRecord_dict[read_id_base]
            if read_mini_mp is not None:

                read_mini_refs = set()
                for mini_ref in read_mini_mp.shared_mini_refs_no_ignored:
                    read_mini_refs.add(mini_ref)
                if read_strand == '1':
                    for mini_ref in read_mini_mp.r1_mini_refs_no_ignored:
                        read_mini_refs.add(mini_ref)
                if read_strand == '2':
                    for mini_ref in read_mini_mp.r2_mini_refs_no_ignored:
                        read_mini_refs.add(mini_ref)

                for each_mini_ref in read_mini_refs:
                    for each_16s_ref in read_16s_refs:
                        mini_assembly_to_16s_key = '%s%s%s' % (each_mini_ref, mini_assembly_to_16s_ctg_connector, each_16s_ref)
                        mini_ref_len = mini_assembly_len_dict[each_mini_ref]
                        s16_ref_len = len(marker_seq_dict[each_16s_ref])

                        if mini_assembly_to_16s_key not in mini_to_16s_LinkingRecord_dict:
                            mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key] = LinkingRecord()
                            mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linked_seq_l     = each_mini_ref
                            mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linked_seq_r     = each_16s_ref
                            mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linked_seq_len_l = mini_ref_len
                            mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linked_seq_len_r = s16_ref_len

                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linking_reads_base.append(read_id_base)

                        # get mini side cigar here
                        current_linking_cigar_mini_side = ''
                        current_linking_pos_mini_side = ''
                        if read_strand == '1':
                            current_linking_cigar_mini_side = list(read_mini_mp.r1_mini_ref_dict[each_mini_ref].values())[0]
                            current_linking_pos_mini_side   = list(read_mini_mp.r1_mini_ref_dict[each_mini_ref].keys())[0]
                        if read_strand == '2':
                            current_linking_cigar_mini_side = list(read_mini_mp.r2_mini_ref_dict[each_mini_ref].values())[0]
                            current_linking_pos_mini_side   = list(read_mini_mp.r2_mini_ref_dict[each_mini_ref].keys())[0]
                        min_dist_to_mini_end = get_min_dist_to_ref_end(current_linking_cigar_mini_side, current_linking_pos_mini_side, mini_ref_len)
                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linking_cigar_l.append(current_linking_cigar_mini_side)
                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linking_pos_l.append(current_linking_pos_mini_side)
                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].min_dist_to_end_l.append(min_dist_to_mini_end)

                        # get 16s side cigar here
                        current_linking_cigar_16s_side = ''
                        current_linking_pos_16s_side = ''
                        if read_strand == '1':
                            current_linking_cigar_16s_side = list(read_16s_mp.r2_16s_ref_dict[each_16s_ref].values())[0]
                            current_linking_pos_16s_side   = list(read_16s_mp.r2_16s_ref_dict[each_16s_ref].keys())[0]
                        if read_strand == '2':
                            current_linking_cigar_16s_side = list(read_16s_mp.r1_16s_ref_dict[each_16s_ref].values())[0]
                            current_linking_pos_16s_side   = list(read_16s_mp.r1_16s_ref_dict[each_16s_ref].keys())[0]
                        min_dist_to_16s_end = get_min_dist_to_ref_end(current_linking_cigar_16s_side, current_linking_pos_16s_side, s16_ref_len)
                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linking_cigar_r.append(current_linking_cigar_16s_side)
                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].linking_pos_r.append(current_linking_pos_16s_side)
                        mini_to_16s_LinkingRecord_dict[mini_assembly_to_16s_key].min_dist_to_end_r.append(min_dist_to_16s_end)

    # filter with link num
    stats_mini_assembly_to_16s_qc_report_handle = open(stats_mini_assembly_to_16s_qc_report, 'w')
    stats_mini_assembly_to_16s_qc_report_handle.write('Mini\t16S\tLinkages\tLinked_to_mini_end(num pct)\tLinked_to_16S_end(num pct)\tshort_cigar_pct_mini\tshort_cigar_pct_16S\n')
    ctgs_to_extract_mini = set()
    all_linking_reads_base_set_rd2_mini_to_16s = set()
    marker_to_gnm_link_num_dict_rd2 = {}
    linking_reads_tab_rd2_handle = open(linking_reads_tab_rd2, 'w')
    stats_GapFilling_ctg_handle = open(stats_GapFilling_ctg, 'w')
    stats_GapFilling_ctg_handle.write('MarkerGene,MiniAssembly,Number\n')
    rd2_mini_to_16s_short_M_pct_handle= open(rd2_mini_to_16s_short_M_pct, 'w')
    rd2_mini_to_16s_short_M_pct_handle.write('linkage\tMini-assembly\tMarker\n')
    mini_assembly_to_16s_passed_qc = set()
    for each_link in mini_to_16s_LinkingRecord_dict.copy():
        linking_reads = mini_to_16s_LinkingRecord_dict[each_link].linking_reads_base
        if len(linking_reads) < ctg_level_min_link:
            mini_to_16s_LinkingRecord_dict.pop(each_link)
        else:
            id_mini = mini_to_16s_LinkingRecord_dict[each_link].linked_seq_l
            id_16s = mini_to_16s_LinkingRecord_dict[each_link].linked_seq_r
            mini_assembly_mag = mini_assembly_to_mag_dict.get(id_mini, None)
            current_linkage_cigar_mini_side = mini_to_16s_LinkingRecord_dict[each_link].linking_cigar_l
            current_linkage_cigar_16s_side = mini_to_16s_LinkingRecord_dict[each_link].linking_cigar_r

            # get pct of short aligned cigar
            short_cigar_pct_mini = get_short_cigar_pct(current_linkage_cigar_mini_side, short_M_len)
            short_cigar_pct_16s = get_short_cigar_pct(current_linkage_cigar_16s_side, short_M_len)

            # check num/pct of reads linked to mini end
            min_dist_list_to_mini_end = mini_to_16s_LinkingRecord_dict[each_link].min_dist_to_end_l
            linked_to_mini_end_cigar_num = 0
            for each_min_dist in min_dist_list_to_mini_end:
                if each_min_dist <= 10:
                    linked_to_mini_end_cigar_num += 1
            linked_to_mini_end_cigar_pct = linked_to_mini_end_cigar_num * 100 / len(min_dist_list_to_mini_end)
            linked_to_mini_end_cigar_pct = float("{0:.2f}".format(linked_to_mini_end_cigar_pct))

            # check num/pct of reads linked to 16s end
            min_dist_list_to_16s_end = mini_to_16s_LinkingRecord_dict[each_link].min_dist_to_end_r
            linked_to_16s_end_cigar_num = 0
            for each_min_dist in min_dist_list_to_16s_end:
                if each_min_dist <= 10:
                    linked_to_16s_end_cigar_num += 1
            linked_to_16s_end_cigar_pct = linked_to_16s_end_cigar_num * 100 / len(min_dist_list_to_16s_end)
            linked_to_16s_end_cigar_pct = float("{0:.2f}".format(linked_to_16s_end_cigar_pct))

            linked_to_mini_end = False
            if (linked_to_mini_end_cigar_num > 0) and (linked_to_mini_end_cigar_pct >= 5):
                linked_to_mini_end = True

            linked_to_16s_end = False
            if (linked_to_16s_end_cigar_num > 0) and (linked_to_16s_end_cigar_pct >= 5):
                linked_to_16s_end = True

            short_cigar_pct_mini_bool = False
            if short_cigar_pct_mini < max_short_cigar_pct_cutoff:
                short_cigar_pct_mini_bool = True

            short_cigar_pct_16s_bool = False
            if short_cigar_pct_16s < max_short_cigar_pct_cutoff:
                short_cigar_pct_16s_bool = True

            # get bp and pct of matched position
            linking_cigar_list_mini = mini_to_16s_LinkingRecord_dict[each_link].linking_cigar_l
            linking_pos_list_mini   = mini_to_16s_LinkingRecord_dict[each_link].linking_pos_l
            ref_len_mini            = mini_to_16s_LinkingRecord_dict[each_link].linked_seq_len_l
            matched_pos_mini_num, matched_pos_mini_pct = get_matched_pos_num_pct(linking_cigar_list_mini, linking_pos_list_mini, ref_len_mini)

            linking_cigar_list_16s = mini_to_16s_LinkingRecord_dict[each_link].linking_cigar_r
            linking_pos_list_16s   = mini_to_16s_LinkingRecord_dict[each_link].linking_pos_r
            ref_len_16s            = mini_to_16s_LinkingRecord_dict[each_link].linked_seq_len_r
            matched_pos_16s_num, matched_pos_16s_pct = get_matched_pos_num_pct(linking_cigar_list_16s, linking_pos_list_16s, ref_len_16s)

            matched_region_passed_qc_mini = False
            if (matched_pos_mini_num >= 300) or (matched_pos_mini_pct >= 50):
                matched_region_passed_qc_mini = True

            matched_region_passed_qc_16s = False
            if (matched_pos_16s_num >= 300) or (matched_pos_16s_pct >= 50):
                matched_region_passed_qc_16s = True

            # write out qc
            stats_mini_assembly_to_16s_qc_report_handle.write('%s\t%s\t%s\t%s(%s %s)__v__%s(%s %s)\t%s(%s %s)__v__%s(%s %s)\t%s(%s)\t%s(%s)\n' % (id_mini, id_16s, len(linking_reads),
                                                                  linked_to_mini_end, linked_to_mini_end_cigar_num, linked_to_mini_end_cigar_pct,
                                                                  matched_region_passed_qc_mini, matched_pos_mini_num, matched_pos_mini_pct,
                                                                  linked_to_16s_end, linked_to_16s_end_cigar_num, linked_to_16s_end_cigar_pct,
                                                                  matched_region_passed_qc_16s, matched_pos_16s_num, matched_pos_16s_pct,
                                                                  short_cigar_pct_mini_bool, short_cigar_pct_mini,
                                                                  short_cigar_pct_16s_bool, short_cigar_pct_16s))

            if (short_cigar_pct_mini_bool is True) and (short_cigar_pct_16s_bool is True):
                if ((linked_to_mini_end is True) or (matched_region_passed_qc_mini is True)) and ((linked_to_16s_end is True) or (matched_region_passed_qc_16s is True)):
                    rd2_mini_to_16s_short_M_pct_handle.write('%s\t%s\t%s\n' % (each_link, short_cigar_pct_mini, short_cigar_pct_16s))
                    mini_assembly_to_16s_passed_qc.add(each_link)

                    if mini_assembly_mag != None:
                        stats_GapFilling_ctg_handle.write('%s,%s%s%s,%s\n' % (id_16s, mini_assembly_mag, gnm_to_ctg_connector, id_mini, len(linking_reads)))
                        linking_reads_tab_rd2_handle.write('%s\t%s%s%s\t%s\n' % (id_16s, mini_assembly_mag, gnm_to_ctg_connector, id_mini, ','.join(linking_reads)))
                        all_linking_reads_base_set_rd2_mini_to_16s.update(linking_reads)
                        ctgs_to_extract_mini.add(id_mini)
                        marker_to_gnm_key_rd2 = '%s%s%s' % (id_16s, marker_to_ctg_gnm_Key_connector, mini_assembly_mag)
                        if marker_to_gnm_key_rd2 not in marker_to_gnm_link_num_dict_rd2:
                            marker_to_gnm_link_num_dict_rd2[marker_to_gnm_key_rd2] = len(linking_reads)
                        else:
                            marker_to_gnm_link_num_dict_rd2[marker_to_gnm_key_rd2] += len(linking_reads)
    rd2_mini_to_16s_short_M_pct_handle.close()
    stats_GapFilling_ctg_handle.close()
    linking_reads_tab_rd2_handle.close()
    stats_mini_assembly_to_16s_qc_report_handle.close()

    # write out linkages at genome level
    stats_GapFilling_gnm_handle = open(stats_GapFilling_file, 'w')
    stats_GapFilling_gnm_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_linkage in marker_to_gnm_link_num_dict_rd2:
        stats_GapFilling_gnm_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (each_linkage.split(marker_to_ctg_gnm_Key_connector)[0], each_linkage.split(marker_to_ctg_gnm_Key_connector)[1], marker_to_gnm_link_num_dict_rd2[each_linkage]))
    stats_GapFilling_gnm_handle.close()

    # filter_linkages_iteratively
    filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, min_iden_16s, min_link_num, min_link_num, within_gnm_linkage_num_diff, stats_GapFilling_file_filtered)

    ################################################# rd2 visualizing ##################################################

    #os.system('python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/d7_Oral/assess_linkages_Oral.py')

    ####################################################################################################################
    ################################################# rd2 visualizing ##################################################
    ####################################################################################################################

    if os.path.isdir(mafft_seq_folder_mini_to_ctg) is True:
        os.system('rm -r %s' % mafft_seq_folder_mini_to_ctg)
    os.mkdir(mafft_seq_folder_mini_to_ctg)

    if os.path.isdir(mafft_seq_folder_mini_to_16s) is True:
        os.system('rm -r %s' % mafft_seq_folder_mini_to_16s)
    os.mkdir(mafft_seq_folder_mini_to_16s)

    #################### extract linking reads ####################

    # combine linking reads
    all_linking_reads_base_set_rd2 = set.union(all_linking_reads_base_set_rd2_mini_to_ctg, all_linking_reads_base_set_rd2_mini_to_16s)

    #########################################################################################################

    argument_lol_for_linkage_vis_worker_mini_to_ctg = []
    for mini_to_ctg in mini_to_mag_LinkingRecord_dict:
        reads_file_base_tmp = mini_to_ctg.replace(mini_assembly_to_16s_ctg_connector, '___')
        reads_file_base     = reads_file_base_tmp.replace(gnm_to_ctg_connector, '___')
        mini_id             = mini_to_ctg.split(mini_assembly_to_16s_ctg_connector)[0]
        ctg_id              = mini_to_ctg.split(mini_assembly_to_16s_ctg_connector)[1]
        mini_linked_ctg = mini_assembly_to_ctg_dict.get(mini_id, None)

        if mini_linked_ctg == ctg_id:
            mini_seq = mini_assembly_seq_dict[mini_id]
            ctg_seq  = unlinked_mag_end_seq_dict[ctg_id]
            linking_reads       = mini_assembly_to_ctg_dict_with_l_r_read_base_min3[mini_to_ctg]

            # create sub folders
            os.mkdir('%s/%s' % (mafft_seq_folder_mini_to_ctg, reads_file_base))
            vis_reads_file_r1 = '%s/%s/%s_R1.fa' % (mafft_seq_folder_mini_to_ctg, reads_file_base, reads_file_base)
            vis_reads_file_r2 = '%s/%s/%s_R2.fa' % (mafft_seq_folder_mini_to_ctg, reads_file_base, reads_file_base)

            # get marker_pos_list and contig_pos_list and write out sequences of linking reads
            marker_pos_list = []
            contig_pos_list = []
            vis_reads_file_r1_handle = open(vis_reads_file_r1, 'w')
            vis_reads_file_r2_handle = open(vis_reads_file_r2, 'w')
            for each_linking_read in linking_reads:
                linking_r1_id = '%s.1' % each_linking_read
                linking_r1_seq = flk_both_read_seq_dict[linking_r1_id]
                linking_r2_id = '%s.2' % each_linking_read
                linking_r2_seq = flk_both_read_seq_dict[linking_r2_id]
                vis_reads_file_r1_handle.write('>%s\n' % linking_r1_id)
                vis_reads_file_r1_handle.write('%s\n'  % linking_r1_seq)
                vis_reads_file_r2_handle.write('>%s\n' % linking_r2_id)
                vis_reads_file_r2_handle.write('%s\n'  % linking_r2_seq)
            vis_reads_file_r1_handle.close()
            vis_reads_file_r2_handle.close()
            fake_pos_list = [1, 1, 1, 1, 1, 1, 1, 1, 1]
            argument_lol_for_linkage_vis_worker_mini_to_ctg.append([reads_file_base, mafft_seq_folder_mini_to_ctg,
                                                                    mini_seq, ctg_seq,
                                                                    end_ctg_len_for_mafft, gap_N_num, bowtie_parameter,
                                                                    fake_pos_list, fake_pos_list,
                                                                    pwd_bowtie2_build_exe, pwd_bowtie2_exe,
                                                                    'Mini', 'MAG'])

    # visualize linkages
    print('Round 2: visualizing %s linkages with %s threads' % (len(argument_lol_for_linkage_vis_worker_mini_to_ctg), num_threads))
    vis_linkages_pool = mp.Pool(processes=num_threads)
    vis_linkages_pool.map(linkage_vis_worker, argument_lol_for_linkage_vis_worker_mini_to_ctg)
    vis_linkages_pool.close()
    vis_linkages_pool.join()


    ####################################### get linking reads for visualization ########################################

    argument_lol_for_linkage_vis_worker = []
    for marker_to_mini in mini_assembly_to_16s_passed_qc:

        reads_file_base_tmp = marker_to_mini.replace(mini_assembly_to_16s_ctg_connector, '___')
        reads_file_base = reads_file_base_tmp.replace(gnm_to_ctg_connector, '___')
        linking_reads = mini_to_16s_LinkingRecord_dict[marker_to_mini].linking_reads_base

        marker_id  = marker_to_mini.split(mini_assembly_to_16s_ctg_connector)[1]
        mini_id    = marker_to_mini.split(mini_assembly_to_16s_ctg_connector)[0]
        marker_seq = marker_seq_dict[marker_id]
        contig_seq = mini_assembly_seq_dict[mini_id]
        mini_assembly_mag = mini_assembly_to_mag_dict.get(mini_id, None)

        if mini_assembly_mag != None:

            # create sub folders
            os.mkdir('%s/%s' % (mafft_seq_folder_mini_to_16s, reads_file_base))
            vis_reads_file_r1 = '%s/%s/%s_R1.fa' % (mafft_seq_folder_mini_to_16s, reads_file_base, reads_file_base)
            vis_reads_file_r2 = '%s/%s/%s_R2.fa' % (mafft_seq_folder_mini_to_16s, reads_file_base, reads_file_base)

            # get marker_pos_list and contig_pos_list and write out sequences of linking reads
            marker_pos_list = []
            contig_pos_list = []
            vis_reads_file_r1_handle = open(vis_reads_file_r1, 'w')
            vis_reads_file_r2_handle = open(vis_reads_file_r2, 'w')
            for each_linking_read in linking_reads:
                linking_r1_id = '%s.1' % each_linking_read
                linking_r1_seq = flk_both_read_seq_dict[linking_r1_id]
                linking_r2_id = '%s.2' % each_linking_read
                linking_r2_seq = flk_both_read_seq_dict[linking_r2_id]
                vis_reads_file_r1_handle.write('>%s\n' % linking_r1_id)
                vis_reads_file_r1_handle.write('%s\n' % linking_r1_seq)
                vis_reads_file_r2_handle.write('>%s\n' % linking_r2_id)
                vis_reads_file_r2_handle.write('%s\n' % linking_r2_seq)

                # get matched position on makrer and contig
                marker_pos_r1 = list(MappingRecord_dict[each_linking_read].r1_16s_ref_dict.get(marker_id, dict()).keys())
                marker_pos_r2 = list(MappingRecord_dict[each_linking_read].r2_16s_ref_dict.get(marker_id, dict()).keys())
                mini_pos_r1   = list(MappingRecord_dict_mini[each_linking_read].r1_mini_ref_dict.get(mini_id, dict()).keys())
                mini_pos_r2   = list(MappingRecord_dict_mini[each_linking_read].r2_mini_ref_dict.get(mini_id, dict()).keys())

                if len(marker_pos_r1) == 1:
                    marker_pos_list.append(marker_pos_r1[0])
                if len(marker_pos_r2) == 1:
                    marker_pos_list.append(marker_pos_r2[0])
                if len(mini_pos_r1) == 1:
                    contig_pos_list.append(mini_pos_r1[0])
                if len(mini_pos_r2) == 1:
                    contig_pos_list.append(mini_pos_r2[0])
            vis_reads_file_r1_handle.close()
            vis_reads_file_r2_handle.close()

            argument_lol_for_linkage_vis_worker.append([reads_file_base, mafft_seq_folder_mini_to_16s,
                                                        marker_seq, contig_seq,
                                                        end_ctg_len_for_mafft, gap_N_num, bowtie_parameter,
                                                        marker_pos_list, contig_pos_list,
                                                        pwd_bowtie2_build_exe, pwd_bowtie2_exe,
                                                        'Marker', 'Mini'])

    ##########################################################################################

    # visualize linkages
    print('Round 2: visualizing %s linkages with %s threads' % (len(argument_lol_for_linkage_vis_worker), num_threads))
    vis_linkages_pool = mp.Pool(processes=num_threads)
    vis_linkages_pool.map(linkage_vis_worker, argument_lol_for_linkage_vis_worker)
    vis_linkages_pool.close()
    vis_linkages_pool.join()

    ####################################################################################################################

    if os.path.isdir(vis_folder_rd2) is True:
        os.system('rm -r %s' % vis_folder_rd2)
    os.mkdir(vis_folder_rd2)

    mag_to_mini_assembly_dict = {}
    for each_mini_assembly in mini_assembly_to_mag_dict:
        assigned_mag = mini_assembly_to_mag_dict[each_mini_assembly]
        if assigned_mag not in mag_to_mini_assembly_dict:
            mag_to_mini_assembly_dict[assigned_mag] = {each_mini_assembly}
        else:
            mag_to_mini_assembly_dict[assigned_mag].add(each_mini_assembly)

    # read in final linkages
    matam_16s_to_mag_rd2 = {}
    for each_rd2_link in open(stats_GapFilling_file_filtered):
        if ',Number\n' not in each_rd2_link:
            each_rd2_link_split = each_rd2_link.strip().split(',')
            id_16s = each_rd2_link_split[0][12:]
            id_mag = each_rd2_link_split[1][12:]
            matam_16s_to_mag_rd2[id_16s] = id_mag

    for each_linked_16s in matam_16s_to_mag_rd2:
        linked_mag = matam_16s_to_mag_rd2[each_linked_16s]
        linked_mini_set = mag_to_mini_assembly_dict[linked_mag]
        matam_to_mag_folder = '%s/%s___%s' % (vis_folder_rd2, each_linked_16s, linked_mag)
        os.mkdir(matam_to_mag_folder)

        for each_linked_mini in linked_mini_set:
            vis_mini_to_16s_folder_name = '%s/%s___%s' % (mafft_seq_folder_mini_to_16s, each_linked_mini, each_linked_16s)
            if os.path.isdir(vis_mini_to_16s_folder_name) is True:
                os.system('cp -r %s %s/' % (vis_mini_to_16s_folder_name, matam_to_mag_folder))

                vis_mini_to_ctg_folder_re = '%s/%s___%s___*' % (mafft_seq_folder_mini_to_ctg, each_linked_mini, linked_mag)
                vis_mini_to_ctg_folder_list = glob.glob(vis_mini_to_ctg_folder_re)
                for each_vis_mini_to_ctg_folder in vis_mini_to_ctg_folder_list:
                    if os.path.isdir(each_vis_mini_to_ctg_folder) is True:
                        os.system('cp -r %s %s/' % (each_vis_mini_to_ctg_folder, matam_to_mag_folder))

    ####################################################################################################################

