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


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, within_gnm_linkage_num_diff, file_out):

    # get MarkerGene_to_GenomicSeq_dict
    MarkerGene_to_GenomicSeq_dict = {}
    for each_linkage in open(file_in):
        if not ((each_linkage.startswith('MarkerGene,GenomicSeq,Number')) or (each_linkage.startswith('MarkerGene,MiniAssembly,Number'))):
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
        if (each_match.startswith('MarkerGene,GenomicSeq,Number')) or (each_match.startswith('MarkerGene,MiniAssembly,Number')):
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

                    # consider depth
                    if min_16s_gnm_multiple > 0:

                        # get marker and genome depth
                        MarkerGene_depth = marker_gene_depth_dict.get(MarkerGene, 'na')
                        GenomicSeq_depth = genomic_seq_depth_dict.get(GenomicSeq, 'na')
                        marker_genome_depth_ratio = 'na'
                        if (MarkerGene_depth != 'na') and (GenomicSeq_depth != 'na'):
                            if GenomicSeq_depth > 0:
                                marker_genome_depth_ratio = MarkerGene_depth / GenomicSeq_depth

                        if marker_genome_depth_ratio >= min_16s_gnm_multiple:
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
                                    if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
                                        file_out_handle.write(each_match)
                                        MarkerGene_with_assignment.add(MarkerGene)
                                    else:
                                        MarkerGene_with_assignment.add(MarkerGene)
                    # ignore depth
                    else:
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


####################################################### file in ########################################################

if __name__ == '__main__':

    step_2_wd                                       = '/Users/songweizhi/Desktop/tunning_rd2'
    reads_file_r1_fasta                             = 'Oral_5_25_R1.fasta'
    reads_file_r2_fasta                             = 'Oral_5_25_R2.fasta'
    free_living_16s_ref_file                        = '%s/file_in/round2_free_living_16s_refs.txt'                              % step_2_wd
    free_living_ctg_ref_file                        = '%s/file_in/round2_free_living_ctg_refs.txt'                              % step_2_wd
    mini_assembly_to_16s_reads                      = '%s/file_in/mini_assembly_to_16s_reads.txt'                               % step_2_wd
    mini_assembly_to_ctg_reads                      = '%s/file_in/mini_assembly_to_ctg_reads.txt'                               % step_2_wd
    blast_results_all_vs_all_16s                    = '%s/file_in/Oral_0705_60_60_16S_all_vs_all_blastn.tab'                    % step_2_wd
    mini_assemblies                                 = '%s/file_in/scaffolds.fasta'                                              % step_2_wd
    marker_gene_seqs                                = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_min1200.QC.fa'               % step_2_wd
    splitted_sam_mp_file_folder                     = '%s/file_in/Oral_0705_60_60_input_reads_to_16S_MappingRecord'             % step_2_wd
    sam_file_mini_assembly_reformatted              = '%s/file_in/scaffolds_bowtie_reformatted_16s.sam'                         % step_2_wd
    combined_1st_round_unlinked_mag_end_seq         = '%s/file_in/round_1_unlinked_gnm_end_500bp.fa'                            % step_2_wd
    rd1_unlinked_mags_sam_bowtie_reformat_sorted    = '%s/file_in/round_1_unlinked_gnm_bowtie_reformatted_sorted.sam'           % step_2_wd
    rd2_to_extract_flking_both_r1_fa                = '%s/file_in/rd2_to_extract_flking_both_r1.fa'                             % step_2_wd
    rd2_to_extract_flking_both_r2_fa                = '%s/file_in/rd2_to_extract_flking_both_r2.fa'                             % step_2_wd

    ctg_level_min_link                              = 3
    within_gnm_linkage_num_diff                     = 80
    max_mini_assembly_link_num_diff_between_ctg_16s = 10
    min_aln_16s                                     = 500
    min_cov_16s                                     = 30
    min_iden_16s                                    = 98
    min_link_num                                    = 8
    marker_to_ctg_gnm_Key_connector                 = '___M___'
    gnm_to_ctg_connector                            = '___C___'
    mini_assembly_to_16s_ctg_connector              = '___Mini___'
    read_to_marker_connector                        = '___r___'
    mean_depth_dict_gnm                             = {}
    mean_depth_dict_16s                             = {}
    min_16s_gnm_multiple                            = 0
    gap_N_num                                       = 50
    end_ctg_len_for_mafft                           = 1000
    num_threads                                     = 4
    pwd_bowtie2_build_exe                           = '/Users/songweizhi/Softwares/bowtie2/bowtie2-build'
    pwd_bowtie2_exe                                 = '/Users/songweizhi/Softwares/bowtie2/bowtie2'
    bowtie_parameter                                = '--local --all --no-unal -N 1 -L 30'
    bowtie_parameter_mini_assembly                  = '--local --all --no-unal -N 1 -L 30'
    seqtk_exe                                       = 'seqtk'

    ##################################################### file out #####################################################

    stats_mini_assembly_to_ctg                      = '%s/stats_mini_assembly_to_ctg.txt'                                   % step_2_wd
    stats_mini_assembly_to_ctg_sorted               = '%s/stats_mini_assembly_to_ctg_sorted.txt'                            % step_2_wd
    stats_mini_assembly_to_ctg_filtered             = '%s/stats_mini_assembly_to_ctg_filtered.txt'                          % step_2_wd
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


    ####################################################################################################################

    splitted_sam_mp_file_re     = '%s/*.txt' % splitted_sam_mp_file_folder
    splitted_sam_mp_file_set    = glob.glob(splitted_sam_mp_file_re)

    pairwise_16s_iden_dict      = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

    ####################################################################################################################

    # read sequence of flanking reads into dict
    flk_both_read_seq_dict = {}
    for linking_r1 in SeqIO.parse(rd2_to_extract_flking_both_r1_fa, 'fasta'):
        flk_both_read_seq_dict[linking_r1.id] = str(linking_r1.seq)
    for linking_r2 in SeqIO.parse(rd2_to_extract_flking_both_r2_fa, 'fasta'):
        flk_both_read_seq_dict[linking_r2.id] = str(linking_r2.seq)

    ##################################################### get reads flanking 16s and ctg ####################################################

    flanking_16s_reads = set()
    flanking_16s_reads_base = set()
    for free_living_read_16s in open(free_living_16s_ref_file):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_16s_id = free_living_read_16s_split[0]
            read_16s_id_base = '.'.join(read_16s_id.split('.')[:-1])
            flanking_16s_reads.add(read_16s_id)
            flanking_16s_reads_base.add(read_16s_id_base)

    flanking_ctg_reads = set()
    flanking_ctg_reads_base = set()
    for each_read_to_ctg_ref in open(free_living_ctg_ref_file):
        read_id = each_read_to_ctg_ref.strip().split('\t')[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        flanking_ctg_reads.add(read_id)
        flanking_ctg_reads_base.add(read_id_base)

    ##################################################### combine reads flanking 16s and contig ####################################################

    rd2_read_extracted_flanking_both_r1_seq     = '%s/rd2_read_to_extract_flanking_both_R1.fa'     % step_2_wd
    rd2_read_extracted_flanking_both_r2_seq     = '%s/rd2_read_to_extract_flanking_both_R2.fa'     % step_2_wd
    rd2_read_extracted_flanking_both_r12_seq    = '%s/rd2_read_to_extract_flanking_both_R12.fa'    % step_2_wd
    rd2_read_extracted_flanking_both_r12_up_seq = '%s/rd2_read_to_extract_flanking_both_R12_up.fa' % step_2_wd

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
    rd2_read_extracted_flanking_both_r1_seq_handle = open(rd2_read_extracted_flanking_both_r1_seq, 'w')
    rd2_read_extracted_flanking_both_r2_seq_handle = open(rd2_read_extracted_flanking_both_r2_seq, 'w')
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

        rd2_read_extracted_flanking_both_r1_seq_handle.write('>%s\n' % flanking_both_r1)
        rd2_read_extracted_flanking_both_r1_seq_handle.write('%s\n'  % flanking_both_r1_seq)
        rd2_read_extracted_flanking_both_r2_seq_handle.write('>%s\n' % flanking_both_r2)
        rd2_read_extracted_flanking_both_r2_seq_handle.write('%s\n'  % flanking_both_r2_seq)

    rd2_read_extracted_flanking_both_r1_seq_handle.close()
    rd2_read_extracted_flanking_both_r2_seq_handle.close()
    rd2_read_extracted_flanking_both_r12_seq_handle.close()

    ######################################### second round linking by assembly #########################################

    sam_file_mini_assembly                          = '%s/scaffolds_bowtie.sam'                         % step_2_wd
    sam_file_mini_assembly_log                      = '%s/scaffolds_bowtie.log'                         % step_2_wd
    sam_file_mini_assembly_reformatted              = '%s/scaffolds_bowtie_reformatted.sam'             % step_2_wd
    sam_file_mini_assembly_reformatted_log          = '%s/scaffolds_bowtie_reformatted.log'             % step_2_wd
    sam_file_mini_assembly_reformatted_sorted       = '%s/file_in/scaffolds_bowtie_reformatted_sorted.sam'      % step_2_wd
    spades_wd                                       = '%s/file_in/mini_assembly_SPAdes_wd'              % step_2_wd
    mini_assemblies                                 = '%s/scaffolds.fasta'                              % spades_wd
    mini_assemblies_no_ext = '.'.join(mini_assemblies.split('.')[:-1])
    rd2_with_both_mates = False

    min_M_len_mini = 45
    mismatch_cutoff = 2


    # mapping with bowtie
    if rd2_with_both_mates is False:
        bowtie_cmd_miniassembly = '%s -x %s -U %s -S %s -p %s -f %s 2> %s'    % (pwd_bowtie2_exe, mini_assemblies_no_ext, rd2_read_extracted_flanking_both_r12_seq, sam_file_mini_assembly, num_threads, bowtie_parameter_mini_assembly, sam_file_mini_assembly_log)
        #bowtie_cmd_miniassembly = '%s -x %s -U %s,%s -S %s -p %s -f %s 2> %s' % (pwd_bowtie2_exe, mini_assemblies_no_ext, rd2_read_extracted_flanking_both_r1_seq, rd2_read_extracted_flanking_both_r2_seq, sam_file_mini_assembly, num_threads, bowtie_parameter_mini_assembly, sam_file_mini_assembly_log)
        print(bowtie_cmd_miniassembly)
        #os.system(bowtie_cmd_miniassembly)

        # reformat cigar string
        mini_assembly_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly, sam_file_mini_assembly_reformatted, sam_file_mini_assembly_reformatted_log)
        print(mini_assembly_reformat_cmd)
        #os.system(mini_assembly_reformat_cmd)

        # sort sam files
        sort_by_read_cmd_mini_assembly = 'samtools sort -n -O sam --threads %s -o %s %s' % (num_threads, sam_file_mini_assembly_reformatted_sorted, sam_file_mini_assembly_reformatted)
        print(sort_by_read_cmd_mini_assembly)
        #os.system(sort_by_read_cmd_mini_assembly)

    ##################################################### read in sam file ####################################################

    print('Round 2: read in sam file')

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

                    ########## filter r1 mini refs ##########
                    r1_mini_refs_passed_qc = {}
                    for r1_ctg_ref_rd2 in current_read_base_r1_mini_ref_dict:
                        r1_matched_pos_dict = current_read_base_r1_mini_ref_dict[r1_ctg_ref_rd2]
                        if len(r1_matched_pos_dict) > 1:
                            mini_refs_to_ignore.add(r1_ctg_ref_rd2)
                        else:
                            r1_ctg_ref_pos = list(r1_matched_pos_dict.keys())[0]
                            r1_ctg_ref_cigar = r1_matched_pos_dict[r1_ctg_ref_pos]
                            r1_ctg_ref_cigar_splitted = cigar_splitter(r1_ctg_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r1_ctg_ref_cigar_splitted)
                            if both_end_clp is True:
                                mini_refs_to_ignore.add(r1_ctg_ref_rd2)
                            else:
                                # check mismatch
                                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                    r1_ctg_ref_cigar_splitted)
                                if r1_ref_min_mismatch == 'NA':
                                    mini_refs_to_ignore.add(r1_ctg_ref_rd2)
                                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                                    mini_refs_to_ignore.add(r1_ctg_ref_rd2)
                                else:
                                    # check if clp in the middle
                                    clip_in_middle = False
                                    if ('S' in r1_ctg_ref_cigar) or ('s' in r1_ctg_ref_cigar):
                                        clip_in_middle = True
                                        if (r1_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                r1_ctg_ref_pos == 1):
                                            clip_in_middle = False
                                        if (r1_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                            if (r1_ctg_ref_pos + r1_aligned_len - 1) == mini_assembly_len_dict[r1_ctg_ref_rd2]:
                                                clip_in_middle = False

                                    if clip_in_middle is True:
                                        mini_refs_to_ignore.add(r1_ctg_ref_rd2)
                                    else:
                                        # check aligned length
                                        if r1_aligned_len < min_M_len_mini:
                                            # mini_refs_to_ignore.add(r1_ctg_ref_rd2)
                                            pass
                                        else:
                                            r1_mini_refs_passed_qc[r1_ctg_ref_rd2] = [r1_ctg_ref_cigar]

                    ########## filter r2 ctg refs ##########
                    r2_mini_refs_passed_qc = {}
                    for r2_ctg_ref_rd2 in current_read_base_r2_mini_ref_dict:
                        r2_matched_pos_dict = current_read_base_r2_mini_ref_dict[r2_ctg_ref_rd2]
                        if len(r2_matched_pos_dict) > 1:
                            mini_refs_to_ignore.add(r2_ctg_ref_rd2)
                        else:
                            r2_ctg_ref_pos = list(r2_matched_pos_dict.keys())[0]
                            r2_ctg_ref_cigar = r2_matched_pos_dict[r2_ctg_ref_pos]
                            r2_ctg_ref_cigar_splitted = cigar_splitter(r2_ctg_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
                            if both_end_clp is True:
                                mini_refs_to_ignore.add(r2_ctg_ref_rd2)
                            else:
                                # check mismatch
                                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                    r2_ctg_ref_cigar_splitted)
                                if r2_ref_min_mismatch == 'NA':
                                    mini_refs_to_ignore.add(r2_ctg_ref_rd2)
                                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                                    mini_refs_to_ignore.add(r2_ctg_ref_rd2)
                                else:
                                    # check if clp in the middle
                                    clip_in_middle = False
                                    if ('S' in r2_ctg_ref_cigar) or ('s' in r2_ctg_ref_cigar):
                                        clip_in_middle = True
                                        if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                r2_ctg_ref_pos == 1):
                                            clip_in_middle = False
                                        if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                            if (r2_ctg_ref_pos + r2_aligned_len - 1) == mini_assembly_len_dict[r2_ctg_ref_rd2]:
                                                clip_in_middle = False

                                    if clip_in_middle is True:
                                        mini_refs_to_ignore.add(r2_ctg_ref_rd2)
                                    else:
                                        # check aligned length
                                        if r2_aligned_len < min_M_len_mini:
                                            # mini_refs_to_ignore.add(r2_ctg_ref_rd2)
                                            pass
                                        else:
                                            r2_mini_refs_passed_qc[r2_ctg_ref_rd2] = [r2_ctg_ref_cigar]

                    ####################################################################################################

                    r1_mini_refs_no_ignored = {key: value for key, value in r1_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
                    r2_mini_refs_no_ignored = {key: value for key, value in r2_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}

                    # if current_read_base == 'S8_242469':
                    #     print(current_read_base_r1_mini_ref_dict)
                    #     print(current_read_base_r2_mini_ref_dict)
                    #     print(r1_mini_refs_no_ignored)
                    #     print(r2_mini_refs_no_ignored)

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
                        shared_mini_refs_no_ignored = {key: [r1_mini_refs_no_ignored[key][0], r2_mini_refs_no_ignored[key][0]] for key in set(r1_mini_refs_no_ignored).intersection(set(r2_mini_refs_no_ignored))}
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

    #read_base = 'S8_242469'
    #print('%s\tr1_mini_ref_dict\t%s'            % (read_base, MappingRecord_dict_mini[read_base].r1_mini_ref_dict))
    #print('%s\tr2_mini_ref_dict\t%s'            % (read_base, MappingRecord_dict_mini[read_base].r2_mini_ref_dict))
    #print('%s\tr1_mini_refs_no_ignored\t%s'     % (read_base, MappingRecord_dict_mini[read_base].r1_mini_refs_no_ignored))
    #print('%s\tr2_mini_refs_no_ignored\t%s'     % (read_base, MappingRecord_dict_mini[read_base].r2_mini_refs_no_ignored))
    #print('%s\tshared_mini_refs_no_ignored\t%s' % (read_base, MappingRecord_dict_mini[read_base].shared_mini_refs_no_ignored))

    ############################################ link mini-assembly to MAGs ############################################

    mini_assembly_to_ctg_dict_with_l_r_read_base = {}
    for free_living_read_ctg in open(free_living_ctg_ref_file):
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
                    mini_assembly_to_ctg_key_with_l_r = '%s%s%s' % (each_mini_ref, mini_assembly_to_16s_ctg_connector, each_ctg_ref)
                    if mini_assembly_to_ctg_key_with_l_r not in mini_assembly_to_ctg_dict_with_l_r_read_base:
                        mini_assembly_to_ctg_dict_with_l_r_read_base[mini_assembly_to_ctg_key_with_l_r] = {read_id_base}
                    else:
                        mini_assembly_to_ctg_dict_with_l_r_read_base[mini_assembly_to_ctg_key_with_l_r].add(read_id_base)


    ############################################ link mini-assembly to MAGs ############################################

    # round2_free_living_ctg_ref_dict = {}
    # round2_free_living_ctg_ref_dict_with_l_r = {}
    # for free_living_read_ctg in open(free_living_ctg_ref_file):
    #     free_living_read_ctg_split = free_living_read_ctg.strip().split('\t')
    #     read_ctg_id = free_living_read_ctg_split[0]
    #     read_ctg_refs = free_living_read_ctg_split[1].split(',')
    #     read_ctg_refs_no_suffix = []
    #     read_ctg_refs_with_suffix = []
    #     for each_read_ctg_ref in read_ctg_refs:
    #         if each_read_ctg_ref[-2:] in ['_l', '_r']:
    #             each_read_ctg_ref_no_suffix = each_read_ctg_ref[:-2]
    #             read_ctg_refs_no_suffix.append(each_read_ctg_ref_no_suffix)
    #             read_ctg_refs_with_suffix.append(each_read_ctg_ref)
    #     round2_free_living_ctg_ref_dict[read_ctg_id] = read_ctg_refs_no_suffix
    #     round2_free_living_ctg_ref_dict_with_l_r[read_ctg_id] = read_ctg_refs_with_suffix
    #
    # mini_assembly_to_ctg_dict = {}
    # mini_assembly_to_ctg_dict_read_base = {}
    # mini_assembly_to_ctg_dict_with_l_r = {}
    # mini_assembly_to_ctg_dict_with_l_r_read_base = {}
    # for each_mini_assembly in open(mini_assembly_to_ctg_reads):
    #     mini_assembly_split = each_mini_assembly.strip().split('\t')
    #     mini_assembly_id = mini_assembly_split[0]
    #     mini_assembly_mapped_reads = mini_assembly_split[1].split(',')
    #     for each_mapped_read in mini_assembly_mapped_reads:
    #         each_mapped_read_base = '.'.join(each_mapped_read.split('.')[:-1])
    #         mapped_read_ctg_refs = round2_free_living_ctg_ref_dict.get(each_mapped_read, [])
    #         mapped_read_ctg_refs_with_l_r = round2_free_living_ctg_ref_dict_with_l_r.get(each_mapped_read, [])
    #         for each_mapped_read_ctg_ref_with_l_r in mapped_read_ctg_refs_with_l_r:
    #             mini_assembly_to_ctg_key_with_l_r = '%s%s%s' % (mini_assembly_id, mini_assembly_to_16s_ctg_connector, each_mapped_read_ctg_ref_with_l_r)
    #             if mini_assembly_to_ctg_key_with_l_r not in mini_assembly_to_ctg_dict_with_l_r:
    #                 mini_assembly_to_ctg_dict_with_l_r[mini_assembly_to_ctg_key_with_l_r] = 1
    #                 mini_assembly_to_ctg_dict_with_l_r_read_base[mini_assembly_to_ctg_key_with_l_r] = {each_mapped_read_base}
    #             else:
    #                 mini_assembly_to_ctg_dict_with_l_r[mini_assembly_to_ctg_key_with_l_r] += 1
    #                 mini_assembly_to_ctg_dict_with_l_r_read_base[mini_assembly_to_ctg_key_with_l_r].add(each_mapped_read_base)

    mini_assembly_to_ctg_dict_with_l_r_read_base_min3 = {}
    for each_link in mini_assembly_to_ctg_dict_with_l_r_read_base:
        if len(mini_assembly_to_ctg_dict_with_l_r_read_base[each_link]) >= 3:
            mini_assembly_to_ctg_dict_with_l_r_read_base_min3[each_link] = mini_assembly_to_ctg_dict_with_l_r_read_base[each_link]

    all_linking_reads_base_set_rd2_mini_to_ctg = set()
    stats_mini_assembly_to_ctg_handle = open(stats_mini_assembly_to_ctg, 'w')
    stats_mini_assembly_to_ctg_handle.write('MiniAssembly,GenomicSeq,Number\n')
    for each_mini_assembly_to_ctg in mini_assembly_to_ctg_dict_with_l_r_read_base_min3:
        linking_reads_mini_to_ctg = mini_assembly_to_ctg_dict_with_l_r_read_base_min3[each_mini_assembly_to_ctg]
        id_mini_assembly = each_mini_assembly_to_ctg.split(mini_assembly_to_16s_ctg_connector)[0]
        id_ctg = each_mini_assembly_to_ctg.split(mini_assembly_to_16s_ctg_connector)[1]
        all_linking_reads_base_set_rd2_mini_to_ctg.update(linking_reads_mini_to_ctg)
        stats_mini_assembly_to_ctg_handle.write('%s,%s,%s\n' % (id_mini_assembly, id_ctg, len(linking_reads_mini_to_ctg)))
    stats_mini_assembly_to_ctg_handle.close()

    # sort and filter
    sort_csv_by_col(stats_mini_assembly_to_ctg, stats_mini_assembly_to_ctg_sorted, 'Number')
    mini_assembly_to_ctg_dict, mini_assembly_to_mag_dict, mag_ctg_set_with_linked_mini_assembly = filter_linkages_iteratively_mini_assembly_to_ctg(stats_mini_assembly_to_ctg_sorted, 3,
                                                                                 stats_mini_assembly_to_ctg_filtered)

    ################################## link 16S to mini-assemblies with MAG assignment #################################

    mini_assembly_to_16s_dict_read_base = {}
    for free_living_read_16s in open(free_living_16s_ref_file):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_id = free_living_read_16s_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            read_16s_refs = free_living_read_16s_split[1].split(',')
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
                    for each_16s_ref in read_16s_refs:
                        mini_assembly_to_16s_key = '%s%s%s' % (each_mini_ref, mini_assembly_to_16s_ctg_connector, each_16s_ref)
                        if mini_assembly_to_16s_key not in mini_assembly_to_16s_dict_read_base:
                            mini_assembly_to_16s_dict_read_base[mini_assembly_to_16s_key] = {read_id_base}
                        else:
                            mini_assembly_to_16s_dict_read_base[mini_assembly_to_16s_key].add(read_id_base)

    # filter with link num
    mini_assembly_to_16s_dict_read_base_min3 = {}
    for each_link in mini_assembly_to_16s_dict_read_base:
        if len(mini_assembly_to_16s_dict_read_base[each_link]) >= 3:
            mini_assembly_to_16s_dict_read_base_min3[each_link] = mini_assembly_to_16s_dict_read_base[each_link]

    ctgs_to_extract_mini = set()
    all_linking_reads_base_set_rd2_mini_to_16s = set()
    marker_to_gnm_link_num_dict_rd2 = {}
    linking_reads_tab_rd2_handle = open(linking_reads_tab_rd2, 'w')
    stats_mini_assembly_to_16s_handle = open(stats_GapFilling_ctg, 'w')
    stats_mini_assembly_to_16s_handle.write('MarkerGene,MiniAssembly,Number\n')
    for each_mini_assembly_to_16s in mini_assembly_to_16s_dict_read_base_min3:
        linking_reads = mini_assembly_to_16s_dict_read_base_min3[each_mini_assembly_to_16s]
        id_mini = each_mini_assembly_to_16s.split(mini_assembly_to_16s_ctg_connector)[0]
        id_16s = each_mini_assembly_to_16s.split(mini_assembly_to_16s_ctg_connector)[1]
        mini_assembly_mag = mini_assembly_to_mag_dict.get(id_mini, None)
        if len(linking_reads) >= ctg_level_min_link:
            if mini_assembly_mag != None:
                stats_mini_assembly_to_16s_handle.write('%s,%s%s%s,%s\n' % (id_16s, mini_assembly_mag, gnm_to_ctg_connector, id_mini, len(linking_reads)))
                linking_reads_tab_rd2_handle.write('%s\t%s%s%s\t%s\n' % (id_16s, mini_assembly_mag, gnm_to_ctg_connector, id_mini, ','.join(linking_reads)))
                all_linking_reads_base_set_rd2_mini_to_16s.update(linking_reads)
                ctgs_to_extract_mini.add(id_mini)
                marker_to_gnm_key_rd2 = '%s%s%s' % (id_16s, marker_to_ctg_gnm_Key_connector, mini_assembly_mag)
                if marker_to_gnm_key_rd2 not in marker_to_gnm_link_num_dict_rd2:
                    marker_to_gnm_link_num_dict_rd2[marker_to_gnm_key_rd2] = len(linking_reads)
                else:
                    marker_to_gnm_link_num_dict_rd2[marker_to_gnm_key_rd2] += len(linking_reads)
    stats_mini_assembly_to_16s_handle.close()
    linking_reads_tab_rd2_handle.close()

    # write out linkages at genome level
    stats_GapFilling_gnm_handle = open(stats_GapFilling_file, 'w')
    stats_GapFilling_gnm_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_linkage in marker_to_gnm_link_num_dict_rd2:
        stats_GapFilling_gnm_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
            each_linkage.split(marker_to_ctg_gnm_Key_connector)[0],
            each_linkage.split(marker_to_ctg_gnm_Key_connector)[1],
            marker_to_gnm_link_num_dict_rd2[each_linkage]))
    stats_GapFilling_gnm_handle.close()

    # filter_linkages_iteratively
    filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm,
                                mean_depth_dict_16s, min_16s_gnm_multiple, min_iden_16s, min_link_num, min_link_num,
                                within_gnm_linkage_num_diff, stats_GapFilling_file_filtered)

    os.system('python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/d7_Oral/assess_linkages_Oral.py')

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

    #################### read sequences into dict ####################

    print('Round 2: read sequences into dict')

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

    #################### prepare seq with multi-processing ####################

    MappingRecord_dict = {}
    for each_mp_file in splitted_sam_mp_file_set:
        with open(each_mp_file) as each_mp_file_opened:
            for each_read_base in each_mp_file_opened:
                each_read_base_split = each_read_base.strip().split('\t')
                current_read_base___id = each_read_base_split[0]
                if current_read_base___id in all_linking_reads_base_set_rd2:
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


    ##################################################### read in sam file ####################################################

    # mini_assembly_len_dict_rd2_vis = {}
    # mini_assembly_mp_dict_rd2_vis = {}
    # with open(sam_file_mini_assembly_reformatted) as sam_file_mini_assembly_reformatted_opened:
    #     for each_read in sam_file_mini_assembly_reformatted_opened:
    #         each_read_split = each_read.strip().split('\t')
    #         if each_read.startswith('@'):
    #             mini_assembly_id = ''
    #             mini_assembly_len = 0
    #             for each_element in each_read_split:
    #                 if each_element.startswith('SN:'):
    #                     mini_assembly_id = each_element[3:]
    #                 if each_element.startswith('LN:'):
    #                     mini_assembly_len = int(each_element[3:])
    #             mini_assembly_len_dict_rd2_vis[mini_assembly_id] = mini_assembly_len
    #         else:
    #             cigar = each_read_split[5]
    #             if cigar != '*':
    #                 read_id = each_read_split[0]
    #                 read_id_base = '.'.join(read_id.split('.')[:-1])
    #                 read_strand = read_id.split('.')[-1]
    #                 ref_id = each_read_split[2]
    #                 ref_pos = int(each_read_split[3])
    #
    #                 if read_id_base not in mini_assembly_mp_dict_rd2_vis:
    #                     mini_assembly_mp_dict_rd2_vis[read_id_base] = MappingRecord()
    #
    #                 if read_strand == '1':
    #                     if ref_id not in mini_assembly_mp_dict_rd2_vis[read_id_base].r1_mini_ref_dict:
    #                         mini_assembly_mp_dict_rd2_vis[read_id_base].r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
    #                     else:
    #                         mini_assembly_mp_dict_rd2_vis[read_id_base].r1_mini_ref_dict[ref_id][ref_pos] = cigar
    #
    #                 if read_strand == '2':
    #                     if ref_id not in mini_assembly_mp_dict_rd2_vis[read_id_base].r2_mini_ref_dict:
    #                         mini_assembly_mp_dict_rd2_vis[read_id_base].r2_mini_ref_dict[ref_id] = {ref_pos: cigar}
    #                     else:
    #                         mini_assembly_mp_dict_rd2_vis[read_id_base].r2_mini_ref_dict[ref_id][ref_pos] = cigar

    #########################################################################################################

    argument_lol_for_linkage_vis_worker_mini_to_ctg = []
    for mini_to_ctg in mini_assembly_to_ctg_dict_with_l_r_read_base_min3:
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
                                                        end_ctg_len_for_mafft, gap_N_num, bowtie_parameter_mini_assembly,
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
    for marker_to_mini in mini_assembly_to_16s_dict_read_base_min3:

        reads_file_base_tmp = marker_to_mini.replace(mini_assembly_to_16s_ctg_connector, '___')
        reads_file_base = reads_file_base_tmp.replace(gnm_to_ctg_connector, '___')
        linking_reads = mini_assembly_to_16s_dict_read_base_min3[marker_to_mini]
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
                                                        end_ctg_len_for_mafft, gap_N_num, bowtie_parameter_mini_assembly,
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
    ####################################################################################################################
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


''' 



bowtie2 -x mini_assembly_SPAdes_wd/scaffolds -U rd2_read_to_extract_flanking_both_R1.fa,rd2_read_to_extract_flanking_both_R2.fa -S bbb.sam -p 4 -f --local --all --no-unal -N 1 -L 30 2> scaffolds_bowtie.log
bowtie2 -x mini_assembly_SPAdes_wd/scaffolds -U rd2_read_to_extract_flanking_both_R1.fa,rd2_read_to_extract_flanking_both_R2.fa -S aaa.sam -p 4 -f --local --all --no-unal -N 1 -L 30 2> scaffolds_bowtie.log

reformat.sh in=scaffolds_bowtie.sam out=scaffolds_bowtie_reformatted.sam sam=1.4 2> scaffolds_bowtie_reformatted.log
samtools sort -n -O sam --threads 4 -o scaffolds_bowtie_reformatted_sorted.sam scaffolds_bowtie_reformatted.sam


bowtie2 -x /Users/songweizhi/Desktop/tunning_rd2/file_in/mini_assembly_SPAdes_wd/scaffolds -U /Users/songweizhi/Desktop/tunning_rd2/rd2_read_to_extract_flanking_both_R1.fa,/Users/songweizhi/Desktop/tunning_rd2/rd2_read_to_extract_flanking_both_R2.fa -S /Users/songweizhi/Desktop/tunning_rd2/scaffolds_bowtie.sam -p 4 -f --local --all --no-unal -N 1 -L 30 2> /Users/songweizhi/Desktop/tunning_rd2/scaffolds_bowtie.log
reformat.sh in=/Users/songweizhi/Desktop/tunning_rd2/scaffolds_bowtie.sam out=/Users/songweizhi/Desktop/tunning_rd2/scaffolds_bowtie_reformatted.sam sam=1.4 2> /Users/songweizhi/Desktop/tunning_rd2/scaffolds_bowtie_reformatted.log
samtools sort -n -O sam --threads 4 -o /Users/songweizhi/Desktop/tunning_rd2/file_in/scaffolds_bowtie_reformatted_sorted.sam /Users/songweizhi/Desktop/tunning_rd2/scaffolds_bowtie_reformatted.sam



# why S7_9224725.2 did not mapped to mini NODE_4

bowtie2 -x mini_assembly_SPAdes_wd/scaffolds -U rd2_read_to_extract_flanking_both_R12.fa -S scaffolds_bowtie.sam -p 12 -f --local --all --no-unal -N 1 -L 30 2> scaffolds_bowtie.log
bowtie2 -x mini_assembly_SPAdes_wd/scaffolds -U rd2_read_to_extract_flanking_both_R1.fa,rd2_read_to_extract_flanking_both_R2.fa -S scaffolds_bowtie.sam -p 4 -f --local --all --no-unal -N 1 -L 30 2> scaffolds_bowtie.log

bowtie2 -x mini_assembly_SPAdes_wd/scaffolds -U rd2_read_to_extract_flanking_both_R1.fa,rd2_read_to_extract_flanking_both_R2.fa -S scaffolds_bowtie.sam -p 4 -f --all --no-unal -N 1 -L 30 2> scaffolds_bowtie.log
reformat.sh in=scaffolds_bowtie.sam out=scaffolds_bowtie_reformatted.sam sam=1.4 2> scaffolds_bowtie_reformatted.log
samtools sort -n -O sam --threads 4 -o scaffolds_bowtie_reformatted_sorted.sam scaffolds_bowtie_reformatted.sam


bowtie2 -x mini_assembly_SPAdes_wd/scaffolds -U S7_9224725.fa -S S7_9224725.sam -p 4 -f --all --no-unal -N 1 -L 30 2> S7_9224725_scaffolds_bowtie.log



bowtie2-build --quiet -f NODE4.fa NODE4
bowtie2 -x NODE4 -U S7_9224725.fa -S S7_9224725_vs_NODE4.sam -p 4 -f --all --no-unal -N 1 -L 30 


# on Katana
bowtie2-build --quiet -f NODE26.fa NODE26
bowtie2 -x NODE26 -U S15_645227.fa -S S15_645227_vs_NODE26_Katana.sam -p 4 -f --all --no-unal -N 1 -L 30 

# on MAC
cd /Users/songweizhi/Desktop
/Users/songweizhi/Softwares/bowtie2/bowtie2-build --quiet -f NODE26.fa NODE26
/Users/songweizhi/Softwares/bowtie2/bowtie2 -x NODE26 -U S15_645227.fa -S S15_645227_vs_NODE26_Mac.sam -p 4 -f --all --no-unal -N 1 -L 30 

index_ref_cmd = '%s --quiet -f %s %s' % (pwd_bowtie2_build_exe, pwd_seq_file_cbd, pwd_seq_file_cbd_index)
bowtie2_cmd = '%s -x %s -U %s,%s -S %s -p 1 -f %s 2> %s' % (
cd /Users/songweizhi/Desktop/tunning_rd2/vis_folder_rd2/CAMI_Oral_subsample_50_571___Oral_57/NODE_26_length_405_cov_0.571942___CAMI_Oral_subsample_50_571
/Users/songweizhi/Softwares/bowtie2/bowtie2-build --quiet -f cbd.fa cbd
/Users/songweizhi/Softwares/bowtie2/bowtie2 -x cbd -U R1.fa,R2.fa -S test.sam -p 1 -f --local --all --no-unal -N 1 -L 30
/Users/songweizhi/Softwares/bowtie2/bowtie2 -x cbd -U S15_645227_R1.fa,S15_645227_R2.fa -S S15_645227_test.sam -p 1 -f --local --all --no-unal -N 1 -L 30
/Users/songweizhi/Softwares/bowtie2/bowtie2-build --quiet -f cbd_2.fa cbd_2
/Users/songweizhi/Softwares/bowtie2/bowtie2 -x cbd_2 -U S15_645227_R1.fa,S15_645227_R2.fa -S S15_645227_test_2.sam -p 1 -f --local --all --no-unal -N 1 -L 30



'''
