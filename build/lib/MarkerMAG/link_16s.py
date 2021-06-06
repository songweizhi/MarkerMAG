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

link_Marker_MAG_usage = '''
=================================== MarkerMAG example commands ===================================

# example commands
MarkerMAG link -p Test -r1 R1.fasta -r2 R2.fasta -marker 16S_seqs.fa -mag MAG_folder -x fa -t 6

# For more details: https://github.com/songweizhi/MarkerMAG

==================================================================================================
'''


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


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
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


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_read_num_and_length(reads_file, tmp_file_location, seqtk_exe):

    reads_file_line_num = '%s/R1_line_num.txt'  % (tmp_file_location)
    reads_file_sub1000  = '%s/R1_sub1000.fasta' % (tmp_file_location)

    # get the number of paired reads
    os.system('wc -l %s > %s' % (reads_file, reads_file_line_num))
    paired_reads_num = int(int(open(reads_file_line_num).readline().strip().split(' ')[0]) / 2)
    if reads_file[-1] in ['Q', 'q']:
        paired_reads_num = int(int(open(reads_file_line_num).readline().strip().split(' ')[0])/4)

    # subsample 1000 reads
    os.system('%s sample -s100 %s 1000 > %s' % (seqtk_exe, reads_file, reads_file_sub1000))

    read_len_list = []
    for each_seq in open(reads_file_sub1000):
        if each_seq[0] not in ['>', '@', '+']:
            read_len_list.append(len(each_seq.strip()))

    read_len_median = np.median(read_len_list)
    read_len_max    = np.max(read_len_list)

    os.system('rm %s' % reads_file_line_num)
    os.system('rm %s' % reads_file_sub1000)

    return paired_reads_num, read_len_median, read_len_max


def sep_paired_and_singleton_reads(fasta_in, fasta_out_r1, fasta_out_r2, fasta_out_singleton):
    reads_pair_dict = {}
    for read_record in SeqIO.parse(fasta_in, 'fasta'):
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

    fasta_out_r1_handle = open(fasta_out_r1, 'w')
    fasta_out_r2_handle = open(fasta_out_r2, 'w')
    fasta_out_singleton_handle = open(fasta_out_singleton, 'w')

    for read_record in SeqIO.parse(fasta_in, 'fasta'):

        read_id_base = '.'.join(read_record.id.split('.')[:-1])
        read_strand = read_record.id.split('.')[-1]

        if read_id_base in read_list_singleton:
            fasta_out_singleton_handle.write('>%s\n' % read_record.id)
            fasta_out_singleton_handle.write('%s\n' % str(read_record.seq))

        if read_id_base in read_list_paired:

            if read_strand == '1':
                fasta_out_r1_handle.write('>%s\n' % read_record.id)
                fasta_out_r1_handle.write('%s\n' % str(read_record.seq))

            if read_strand == '2':
                fasta_out_r2_handle.write('>%s\n' % read_record.id)
                fasta_out_r2_handle.write('%s\n' % str(read_record.seq))

    fasta_out_r1_handle.close()
    fasta_out_r2_handle.close()
    fasta_out_singleton_handle.close()


def remove_clp_in_middle(sam_in, sam_out):

    sam_out_handle = open(sam_out, 'w')

    marker_len_dict = {}
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
            marker_len_dict[marker_id] = marker_len
        else:
            cigar = each_read_split[5]

            # check if clp in the middle
            if ('S' not in cigar) and ('s' not in cigar):
                sam_out_handle.write(each_read)
            else:
                ref_id = each_read_split[2]
                ref_pos = int(each_read_split[3])
                cigar_splitted = cigar_splitter(cigar)
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)

                clip_in_middle = True
                if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                    clip_in_middle = False
                if (cigar_splitted[-1][-1] in ['S', 's']):
                    if (ref_pos + r1_aligned_len - 1) == marker_len_dict[ref_id]:
                        clip_in_middle = False

                if clip_in_middle is False:
                    sam_out_handle.write(each_read)

    sam_out_handle.close()


def remove_high_mismatch(sam_in, mismatch_cutoff, sam_out):
    sam_out_handle = open(sam_out, 'w')
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)
        else:
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                cigar_splitted = cigar_splitter(cigar)
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)
                if r1_mismatch_pct <= mismatch_cutoff:
                    sam_out_handle.write(each_read)
    sam_out_handle.close()


def get_ctg_mean_depth_by_samtools_coverage(index_ref, ref_seq, reads_r1, reads_r2, reads_unpaired, subsample_rate, num_threads):

    ref_seq_file_path, ref_seq_file_basename, ref_seq_file_extension = sep_path_basename_ext(ref_seq)

    sam_file                                                            = '%s/%s.sam'                                                           % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_one_end_clp                                                = '%s/%s_one_end_clp.sam'                                               % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_one_end_clp_no_middle                                      = '%s/%s_one_end_clp_no_middle.sam'                                     % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_one_end_clp_no_middle_reformatted                          = '%s/%s_one_end_clp_no_middle_reformatted.sam'                         % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_one_end_clp_no_middle_reformatted_log                      = '%s/%s_one_end_clp_no_middle_reformatted.log'                         % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_one_end_clp_no_middle_reformatted_best_match               = '%s/%s_one_end_clp_no_middle_reformatted_best_match.sam'              % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_one_end_clp_no_middle_reformatted_best_match_low_mismatch  = '%s/%s_one_end_clp_no_middle_reformatted_best_match_low_mismatch.sam' % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_sorted                                                     = '%s/%s_sorted.sam'                                                    % (ref_seq_file_path, ref_seq_file_basename)
    coverage_file                                                       = '%s/%s_cov.txt'                                                       % (ref_seq_file_path, ref_seq_file_basename)

    # build reference index
    cmd_bowtie2_build   = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, ref_seq, ref_seq_file_path, ref_seq_file_basename)
    if index_ref is True:
        os.system(cmd_bowtie2_build)

    # mapping
    # if reads_unpaired == '':
    #     cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s --all -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
    # else:
    #     cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -1 %s -2 %s -U %s -S %s -p %s --all -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, reads_unpaired, sam_file, num_threads)

    # mapping
    if reads_unpaired == '':
        cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -U %s,%s -S %s -p %s --local --all --no-unal -N 1 -L 30 -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
    else:
        cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -U %s,%s,%s -S %s -p %s --local --all --no-unal -N 1 -L 30 -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, reads_unpaired, sam_file, num_threads)

    os.system(cmd_bowtie2_mapping)

    # filter mapping
    remove_both_ends_clp(sam_file, sam_file_one_end_clp)
    remove_clp_in_middle(sam_file_one_end_clp, sam_file_one_end_clp_no_middle)
    bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_one_end_clp_no_middle, sam_file_one_end_clp_no_middle_reformatted, sam_file_one_end_clp_no_middle_reformatted_log)
    os.system(bbmap_reformat_cmd)
    keep_best_matches_in_sam_keep_short_M(sam_file_one_end_clp_no_middle_reformatted, 35, sam_file_one_end_clp_no_middle_reformatted_best_match)
    remove_high_mismatch(sam_file_one_end_clp_no_middle_reformatted_best_match, 2, sam_file_one_end_clp_no_middle_reformatted_best_match_low_mismatch)

    # sort mapping
    cmd_samtools_sort = 'samtools sort %s -o %s' % (sam_file_one_end_clp_no_middle_reformatted_best_match_low_mismatch, sam_file_sorted)
    os.system(cmd_samtools_sort)

    # get mean depth
    cmd_samtools_coverage = 'samtools coverage --ff 4 %s -o %s' % (sam_file_sorted, coverage_file)
    os.system(cmd_samtools_coverage)

    # remove sam files
    os.system('rm %s' % sam_file)
    # os.system('rm %s' % sam_file_sorted)

    # store mean depth into dict
    mean_depth_dict_ctg = {}
    ctg_len_dict = {}
    for each_ctg_depth in open(coverage_file):
        if not each_ctg_depth.startswith('#'):
            ctg_depth_split = each_ctg_depth.strip().split('\t')
            ctg_id = ctg_depth_split[0]
            ctg_len = int(ctg_depth_split[2])
            ctg_depth = float(ctg_depth_split[6]) * (1 / subsample_rate)
            mean_depth_dict_ctg[ctg_id] = ctg_depth
            ctg_len_dict[ctg_id] = ctg_len

    return mean_depth_dict_ctg, ctg_len_dict


def remove_reads_with_multi_best_aln(sam_in, sam_out):

    sam_out_tmp = '%s.tmp' % sam_out

    multi_aligned_reads = set()
    best_hit_cigar_dict = {}
    sam_out_best_hits_handle = open(sam_out_tmp, 'w')
    for each_line in open(sam_in):
        if each_line.startswith('@'):
            sam_out_best_hits_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            cigar   = each_line_split[5]
            if read_id not in best_hit_cigar_dict:
                best_hit_cigar_dict[read_id] = cigar
                sam_out_best_hits_handle.write(each_line)
            else:
                if cigar == best_hit_cigar_dict[read_id]:
                    sam_out_best_hits_handle.write(each_line)
                    multi_aligned_reads.add(read_id)
    sam_out_best_hits_handle.close()

    sam_out_no_ambiguous_handle = open(sam_out, 'w')
    for best_aln in open(sam_out_tmp):
        if best_aln.startswith('@'):
            sam_out_no_ambiguous_handle.write(best_aln)
        else:
            read_id = best_aln.strip().split('\t')[0]
            if read_id not in multi_aligned_reads:
                sam_out_no_ambiguous_handle.write(best_aln)
    sam_out_no_ambiguous_handle.close()


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


def split_list(list_in, subset_num):

    list_in_formatted = [i for i in list_in]

    # get the number of element per subset
    file_num_per_folder = round(len(list_in_formatted) / subset_num)

    n = 1
    lol_out = []
    while n <= subset_num:

        if n < subset_num:
            current_subset_elements = {i for i in list_in_formatted[(file_num_per_folder * (n - 1)):(file_num_per_folder * n)]}
            lol_out.append(current_subset_elements)
        else:
            current_subset_elements = {i for i in list_in_formatted[(file_num_per_folder * (n - 1)):]}
            lol_out.append(current_subset_elements)

        n += 1

    return lol_out


def stats_dict_to_sankey_file_in(clipping_stats_dict, paired_stats_dict, sankey_file_in_clipping, sankey_file_in_paired):

    # prepare input file for plot of clipping mapped reads
    sankey_file_in_clipping_handle = open(sankey_file_in_clipping, 'w')
    sankey_file_in_clipping_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_clipping in clipping_stats_dict:
        sankey_file_in_clipping_handle.write('%s,%s\n' % (','.join(each_clipping.split('_|_')), clipping_stats_dict[each_clipping]))
    sankey_file_in_clipping_handle.close()

    # prepare input file for plot of paired reads
    sankey_file_in_paired_handle = open(sankey_file_in_paired, 'w')
    sankey_file_in_paired_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_paired in paired_stats_dict:
        sankey_file_in_paired_handle.write('%s,%s\n' % (','.join(each_paired.split('_|_')), paired_stats_dict[each_paired]))
    sankey_file_in_paired_handle.close()


def sort_csv_by_col(file_in, file_out, col_header):
    df_in        = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


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

        if (align_len >= align_len_cutoff) and (query != subject) and (coverage_q >= cov_cutoff) and (coverage_s >= cov_cutoff):
            query_to_subject_key = '__|__'.join(sorted([query, subject]))
            if query_to_subject_key not in pairwise_iden_dict:
                pairwise_iden_dict[query_to_subject_key] = iden
            else:
                if iden > pairwise_iden_dict[query_to_subject_key]:
                    pairwise_iden_dict[query_to_subject_key] = iden

    return pairwise_iden_dict


def filter_linkages_iteratively_backup(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, file_out):

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

                # consider depth
                if min_16s_gnm_multiple > 0:
                    MarkerGene_depth = marker_gene_depth_dict[MarkerGene]
                    GenomicSeq_depth = genomic_seq_depth_dict[GenomicSeq]
                    if (MarkerGene_depth/GenomicSeq_depth) >= min_16s_gnm_multiple:
                        if MarkerGene not in MarkerGene_with_assignment:

                            if GenomicSeq not in GenomicSeq_best_marker_dict:
                                GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)
                            else:
                                current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                                key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))

                                iden_with_best_marker = 0
                                if key_str in pairwise_16s_iden_dict:
                                    iden_with_best_marker = pairwise_16s_iden_dict[key_str]

                                if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                                    file_out_handle.write(each_match)
                                    MarkerGene_with_assignment.add(MarkerGene)
                # ignore depth
                else:
                    if MarkerGene not in MarkerGene_with_assignment:
                        if GenomicSeq not in GenomicSeq_best_marker_dict:
                            GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                            file_out_handle.write(each_match)
                            MarkerGene_with_assignment.add(MarkerGene)
                        else:
                            current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                            key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))

                            iden_with_best_marker = 0
                            if key_str in pairwise_16s_iden_dict:
                                iden_with_best_marker = pairwise_16s_iden_dict[key_str]

                            if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)

    file_out_handle.close()

    # remove tmp file
    # os.remove(file_in_sorted)


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, within_gnm_linkage_num_diff, file_out):

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


def filter_linkages_iteratively_new(sorted_file_in, pairwise_16s_iden_dict, within_genome_16s_divergence_cutoff, marker_len_dict,
                                    min_linkages, within_gnm_linkage_num_diff, file_out,
                                    marker_to_gnm_linking_cigar_dict_16s_side,
                                    marker_to_gnm_linking_cigar_dict_ctg_side,
                                    marker_to_ctg_gnm_Key_connector):

    # filter linkage
    gnm_max_link_num_dict = {}
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    MarkerGene_to_be_ignored = set()
    current_gnm = ''
    current_gnm_best_16s_list = []
    current_gnm_highest_link_num = 0
    gnm_with_assignment = set()
    gnm_to_assignmed_16s_dict = dict()
    best_16s_processed_gnm_list = set()
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

                # first check if linked to conserved regions
                already_assigned_16s_list = []
                iden_with_already_assigned_16s_list = []
                for already_assigned_16s in MarkerGene_with_assignment:
                    current_key = '__|__'.join(sorted([already_assigned_16s, MarkerGene]))
                    current_key_value = pairwise_16s_iden_dict.get(current_key, 0)
                    already_assigned_16s_list.append(already_assigned_16s)
                    iden_with_already_assigned_16s_list.append(current_key_value)

                if len(already_assigned_16s_list) > 0:
                    sorted_best_matched_16s_list = [[seq_id, mean_iden] for mean_iden, seq_id in sorted(zip(iden_with_already_assigned_16s_list, already_assigned_16s_list), reverse=True)]
                    best_matched_marker      = sorted_best_matched_16s_list[0][0]
                    best_matched_marker_iden = sorted_best_matched_16s_list[0][1]
                    best_matched_marker_len  = marker_len_dict[best_matched_marker]
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

                if (MarkerGene in MarkerGene_with_assignment) or (MarkerGene in MarkerGene_to_be_ignored):
                    pass
                else:
                    if current_gnm == '':
                        current_gnm = GenomicSeq
                        current_gnm_best_16s_list.append(MarkerGene)
                        current_gnm_highest_link_num = linkage_num

                    elif current_gnm == GenomicSeq:

                        if linkage_num == current_gnm_highest_link_num:
                            current_gnm_best_16s_list.append(MarkerGene)
                        else:
                            # process markers with highest number of linking reads
                            if current_gnm not in best_16s_processed_gnm_list:

                                gnm_max_link_num_dict[current_gnm] = current_gnm_highest_link_num

                                if len(current_gnm_best_16s_list) == 1:
                                    file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (current_gnm_best_16s_list[0], current_gnm, current_gnm_highest_link_num))
                                    gnm_with_assignment.add(current_gnm)
                                    MarkerGene_with_assignment.add(current_gnm_best_16s_list[0])
                                    gnm_to_assignmed_16s_dict[current_gnm] = {current_gnm_best_16s_list[0]}
                                    best_16s_processed_gnm_list.add(current_gnm)

                                else:
                                    # process here
                                    current_gnm_best_16s_mean_iden_list = get_mean_iden_list(current_gnm_best_16s_list, pairwise_16s_iden_dict)
                                    sorted_best_16s_list = [seq_id for mean_iden, seq_id in sorted(zip(current_gnm_best_16s_mean_iden_list, current_gnm_best_16s_list), reverse=True)]

                                    add_index = 0
                                    for linked_16s in sorted_best_16s_list:
                                        if add_index == 0:
                                            file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
                                            gnm_with_assignment.add(current_gnm)
                                            MarkerGene_with_assignment.add(linked_16s)
                                            if current_gnm not in gnm_to_assignmed_16s_dict:
                                                gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
                                            else:
                                                gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
                                        else:
                                            # get identity list with already assigned 16s
                                            already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(current_gnm)
                                            iden_with_already_assigned_16s_list = []
                                            for each_assigned_16s in already_assigned_16s_list:
                                                marker_key = '__|__'.join(sorted([each_assigned_16s, linked_16s]))
                                                marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
                                                iden_with_already_assigned_16s_list.append(marker_key_iden)
                                            if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:
                                                file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
                                                gnm_with_assignment.add(current_gnm)
                                                MarkerGene_with_assignment.add(linked_16s)
                                                if current_gnm not in gnm_to_assignmed_16s_dict:
                                                    gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
                                                else:
                                                    gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
                                        add_index += 1
                                    best_16s_processed_gnm_list.add(current_gnm)

                            # process current one here
                            # get identity list with already assigned 16s
                            already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(GenomicSeq)
                            iden_with_already_assigned_16s_list = []
                            for each_assigned_16s in already_assigned_16s_list:
                                marker_key = '__|__'.join(sorted([each_assigned_16s, MarkerGene]))
                                marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
                                iden_with_already_assigned_16s_list.append(marker_key_iden)
                            if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:

                                gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
                                if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
                                    file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (MarkerGene, GenomicSeq, linkage_num))
                                    gnm_with_assignment.add(GenomicSeq)
                                    MarkerGene_with_assignment.add(MarkerGene)
                                    if GenomicSeq not in gnm_to_assignmed_16s_dict:
                                        gnm_to_assignmed_16s_dict[GenomicSeq] = {MarkerGene}
                                    else:
                                        gnm_to_assignmed_16s_dict[GenomicSeq].add(MarkerGene)
                                else:
                                    # check length, if shorter than assigned, ignore this one by adding it to MarkerGene_to_be_ignored
                                    MarkerGene_len = marker_len_dict[MarkerGene]
                                    already_assigned_16s_len_list = [marker_len_dict[s16] for s16 in already_assigned_16s_list]
                                    median_assign_16s_len = np.median(already_assigned_16s_len_list)
                                    if (median_assign_16s_len - MarkerGene_len) > 50:
                                        MarkerGene_to_be_ignored.add(MarkerGene)
                    else:
                        # first check if previous one processed
                        if current_gnm not in best_16s_processed_gnm_list:

                            gnm_max_link_num_dict[current_gnm] = current_gnm_highest_link_num

                            # process unprocessed
                            if len(current_gnm_best_16s_list) == 1:
                                file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (current_gnm_best_16s_list[0], current_gnm, current_gnm_highest_link_num))
                                gnm_with_assignment.add(current_gnm)
                                MarkerGene_with_assignment.add(current_gnm_best_16s_list[0])
                                gnm_to_assignmed_16s_dict[current_gnm] = {current_gnm_best_16s_list[0]}
                                best_16s_processed_gnm_list.add(current_gnm)
                            else:
                                current_gnm_best_16s_mean_iden_list = get_mean_iden_list(current_gnm_best_16s_list, pairwise_16s_iden_dict)
                                sorted_best_16s_list = [seq_id for mean_iden, seq_id in sorted(zip(current_gnm_best_16s_mean_iden_list, current_gnm_best_16s_list), reverse=True)]

                                add_index = 0
                                for linked_16s in sorted_best_16s_list:
                                    if add_index == 0:
                                        file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
                                        gnm_with_assignment.add(current_gnm)
                                        MarkerGene_with_assignment.add(linked_16s)
                                        if current_gnm not in gnm_to_assignmed_16s_dict:
                                            gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
                                        else:
                                            gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
                                    else:
                                        # get identity list with already assigned 16s
                                        already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(current_gnm)
                                        iden_with_already_assigned_16s_list = []
                                        for each_assigned_16s in already_assigned_16s_list:
                                            marker_key = '__|__'.join(sorted([each_assigned_16s, linked_16s]))
                                            marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
                                            iden_with_already_assigned_16s_list.append(marker_key_iden)

                                        if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:
                                            file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (linked_16s, current_gnm, current_gnm_highest_link_num))
                                            gnm_with_assignment.add(current_gnm)
                                            MarkerGene_with_assignment.add(linked_16s)
                                            if current_gnm not in gnm_to_assignmed_16s_dict:
                                                gnm_to_assignmed_16s_dict[current_gnm] = {linked_16s}
                                            else:
                                                gnm_to_assignmed_16s_dict[current_gnm].add(linked_16s)
                                    add_index += 1
                                best_16s_processed_gnm_list.add(current_gnm)

                        # then process current one
                        if GenomicSeq in best_16s_processed_gnm_list:
                            # get identity list with already assigned 16s
                            already_assigned_16s_list = gnm_to_assignmed_16s_dict.get(GenomicSeq)


                            iden_with_already_assigned_16s_list = []
                            for each_assigned_16s in already_assigned_16s_list:
                                marker_key = '__|__'.join(sorted([each_assigned_16s, MarkerGene]))
                                marker_key_iden = pairwise_16s_iden_dict.get(marker_key, 0)
                                iden_with_already_assigned_16s_list.append(marker_key_iden)

                            if min(iden_with_already_assigned_16s_list) >= within_genome_16s_divergence_cutoff:

                                gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
                                if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
                                    file_out_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (MarkerGene, GenomicSeq, linkage_num))
                                    gnm_with_assignment.add(GenomicSeq)
                                    MarkerGene_with_assignment.add(MarkerGene)
                                    if GenomicSeq not in gnm_to_assignmed_16s_dict:
                                        gnm_to_assignmed_16s_dict[GenomicSeq] = {MarkerGene}
                                    else:
                                        gnm_to_assignmed_16s_dict[GenomicSeq].add(MarkerGene)
                                else:
                                    # check length, if shorter than assigned, ignore this one by adding it to MarkerGene_with_assignment
                                    already_assigned_16s_len_list = [marker_len_dict[s16] for s16 in already_assigned_16s_list]
                                    median_assign_16s_len = np.median(already_assigned_16s_len_list)
                                    if (median_assign_16s_len - MarkerGene_len) > 50:
                                        MarkerGene_to_be_ignored.add(MarkerGene)
                        else:
                            current_gnm = GenomicSeq
                            current_gnm_best_16s_list = [MarkerGene]
                            current_gnm_highest_link_num = linkage_num

    #print('current_gnm: %s' % current_gnm)
    #print('current_gnm_highest_link_num: %s' % current_gnm_highest_link_num)
    #print('current_gnm_best_16s_list: %s' % current_gnm_best_16s_list)

    file_out_handle.close()


def combine_paired_and_clipping_linkages(paired_linkages, clipping_linkages, file_out_summary, file_out_intersect_linkages):

    # file in:   file_in_paired    and  file_in_clipping
    # file out:  file_out_summary  and  file_out_intersection

    combined_paired_and_clipping_keys = set()

    # read in paired linkages
    paired_linkages_dict = {}
    for paired_linkage in open(paired_linkages):
        if not paired_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            paired_linkage_split = paired_linkage.strip().split(',')
            paired_key = '%s__|__%s' % (paired_linkage_split[0], paired_linkage_split[1])
            paired_value = int(paired_linkage_split[2])
            paired_linkages_dict[paired_key] = paired_value
            combined_paired_and_clipping_keys.add(paired_key)

    # read in clipping linkages
    clipping_linkages_dict = {}
    for clipping_linkage in open(clipping_linkages):
        if not clipping_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            clipping_linkage_split = clipping_linkage.strip().split(',')
            clipping_key = '%s__|__%s' % (clipping_linkage_split[0], clipping_linkage_split[1])
            clipping_value = int(clipping_linkage_split[2])
            clipping_linkages_dict[clipping_key] = clipping_value
            combined_paired_and_clipping_keys.add(clipping_key)

    combined_paired_and_clipping_keys_sorted = sorted([i for i in combined_paired_and_clipping_keys])

    # combine paired and clipping linkages
    file_out_summary_handle = open(file_out_summary, 'w')
    file_out_intersect_linkages_handle = open(file_out_intersect_linkages, 'w')
    file_out_summary_handle.write('MarkerGene\tGenomicSeq\tPaired\tClipping\n')
    file_out_intersect_linkages_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_key in combined_paired_and_clipping_keys_sorted:

        current_key_paired_value = 0
        if each_key in paired_linkages_dict:
            current_key_paired_value = paired_linkages_dict[each_key]

        current_key_clipping_value = 0
        if each_key in clipping_linkages_dict:
            current_key_clipping_value = clipping_linkages_dict[each_key]

        if current_key_paired_value > 0:

            current_key_combined = current_key_paired_value + current_key_clipping_value

            # write out
            file_out_summary_handle.write('%s\t%s\t%s\n' % ('\t'.join([i[12:] for i in each_key.split('__|__')]), current_key_paired_value, current_key_clipping_value))
            file_out_intersect_linkages_handle.write('%s,%s\n' % (','.join(each_key.split('__|__')), current_key_combined))

    file_out_summary_handle.close()
    file_out_intersect_linkages_handle.close()


def get_unlinked_mag_end_seq(ref_in, ref_in_end_seq, end_seq_len):

    # get ref seqs subset
    ref_subset_handle = open(ref_in_end_seq, 'w')
    for ref_seq in SeqIO.parse(ref_in, 'fasta'):

        ref_seq_id = ref_seq.id
        ref_seq_len = len(ref_seq.seq)

        if ref_seq_len < end_seq_len * 2:
            ref_subset_handle.write('>%s\n' % ref_seq_id)
            ref_subset_handle.write('%s\n' % ref_seq.seq)
        else:
            ref_seq_left_end_id = '%s_l' % ref_seq_id
            ref_seq_right_end_id = '%s_r' % ref_seq_id
            ref_seq_left_end = ref_seq.seq[:end_seq_len]
            ref_seq_right_end = ref_seq.seq[-end_seq_len:]

            # write out left end
            ref_subset_handle.write('>%s\n' % ref_seq_left_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_left_end)

            # write out right end
            ref_subset_handle.write('>%s\n' % ref_seq_right_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_right_end)
    ref_subset_handle.close()


def get_free_living_mate(ref_in, sam_file, reads_r1, reads_r2, end_seq_len, num_threads, pwd_bbmap_exe, bbmap_memory):

    ref_in_path, ref_in_basename, ref_in_ext = sep_path_basename_ext(ref_in)

    ref_subset      = '%s/%s_ends_%sbp%s'                % (ref_in_path, ref_in_basename, end_seq_len, ref_in_ext)
    bbmap_stderr    = '%s/%s_ends_%sbp_bbmap_stderr.txt' % (ref_in_path, ref_in_basename, end_seq_len)

    # get ref seqs subset
    ref_subset_handle = open(ref_subset, 'w')
    for ref_seq in SeqIO.parse(ref_in, 'fasta'):

        ref_seq_id = ref_seq.id
        ref_seq_len = len(ref_seq.seq)

        if ref_seq_len < end_seq_len * 2:
            ref_subset_handle.write('>%s\n' % ref_seq_id)
            ref_subset_handle.write('%s\n' % ref_seq.seq)
        else:
            ref_seq_left_end_id = '%s_l' % ref_seq_id
            ref_seq_right_end_id = '%s_r' % ref_seq_id
            ref_seq_left_end = ref_seq.seq[:end_seq_len]
            ref_seq_right_end = ref_seq.seq[-end_seq_len:]

            # write out left end
            ref_subset_handle.write('>%s\n' % ref_seq_left_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_left_end)

            # write out right end
            ref_subset_handle.write('>%s\n' % ref_seq_right_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_right_end)
    ref_subset_handle.close()

    # mapping with bbmap
    bbmap_parameter_round2 = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % (num_threads, bbmap_memory)
    bbmap_cmd_round2 = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, ref_subset, reads_r1, reads_r2, sam_file, bbmap_parameter_round2, bbmap_stderr)
    os.system(bbmap_cmd_round2)

    # mapping with bowtie


def get_best_ctg_or_16s_for_gap_seq_iteratively(file_in, sort_by_col_header, min_linkages, file_out):

    file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
    file_in_sorted = '%s/%s_sorted%s' % (file_in_path, file_in_basename, file_in_extension)

    # sort file in
    sort_csv_by_col(file_in, file_in_sorted, sort_by_col_header)

    # fileter linkage
    file_out_handle = open(file_out, 'w')
    gap_seq_with_assignment = set()
    for each_match in open(file_in_sorted):
        if each_match.startswith('Gap_seq,'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            gap_seq_id = match_split[0]
            linkage_num = int(match_split[2])
            if (linkage_num >= min_linkages) and (gap_seq_id not in gap_seq_with_assignment):
                file_out_handle.write(each_match)
                gap_seq_with_assignment.add(gap_seq_id)
    file_out_handle.close()

    # remove tmp file
    # os.remove(file_in_sorted)


def get_accuracy(file_in, marker_num):

    linkage_num_total = 0
    linkage_num_correct = 0
    recovered_markers = set()
    for each_match in open(file_in):
        if not each_match.startswith('MarkerGene\tGenomicSeq\tLinkage'):
            match_split = each_match.strip().split('\t')
            linkage_num = int(match_split[2])
            MarkerGene_genome = match_split[0][:2]
            GenomicSeq_genome = match_split[1]

            linkage_num_total += linkage_num
            if MarkerGene_genome == GenomicSeq_genome:
                linkage_num_correct += linkage_num
                recovered_markers.add(match_split[0])

    marker_recovery = float("{0:.2f}".format(len(recovered_markers)*100/marker_num))

    link_accuracy = 0
    if linkage_num_total > 0:
        link_accuracy = float("{0:.2f}".format(linkage_num_correct*100/linkage_num_total))

    marker_recovery = '%s/%s(%s)' % (len(recovered_markers), marker_num, marker_recovery)

    return marker_recovery, link_accuracy, recovered_markers


def get_accuracy_by_genome(file_in, mag_folder, mag_file_extension):

    # get MAG file list
    mag_file_re             = '%s/*%s' % (mag_folder, mag_file_extension)
    mag_file_list           = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
    mag_file_list_no_ext    = {'.'.join(i.split('.')[:-1]) for i in mag_file_list}

    genome_with_right_16s_assignment_tmp = set()
    genome_with_wrong_16s_assignment = set()
    for each_match in open(file_in):
        if not each_match.startswith('MarkerGene\tGenomicSeq\tLinkage'):
            match_split = each_match.strip().split('\t')
            MarkerGene_genome = match_split[0][:2]
            GenomicSeq_genome = match_split[1]

            if GenomicSeq_genome == MarkerGene_genome:
                genome_with_right_16s_assignment_tmp.add(GenomicSeq_genome)
            else:
                genome_with_wrong_16s_assignment.add(GenomicSeq_genome)

    genome_with_right_16s_assignment_always = []
    genome_without_right_16s_assignment = []
    for input_genome in mag_file_list_no_ext:
        if (input_genome in genome_with_right_16s_assignment_tmp) and (input_genome not in genome_with_wrong_16s_assignment):
            genome_with_right_16s_assignment_always.append(input_genome)
        else:
            genome_without_right_16s_assignment.append(input_genome)


    marker_gene_assignment_rate = float("{0:.2f}".format(len(genome_with_right_16s_assignment_always)*100/len(mag_file_list_no_ext)))

    marker_gene_assignment_accuracy = 0
    if (len(genome_with_right_16s_assignment_always) + len(genome_with_wrong_16s_assignment)) > 0:
        marker_gene_assignment_accuracy = float("{0:.2f}".format(len(genome_with_right_16s_assignment_always)*100/(len(genome_with_right_16s_assignment_always) + len(genome_with_wrong_16s_assignment))))
    marker_gene_assignment_rate = '%s/%s(%s)' % (len(genome_with_right_16s_assignment_always), len(mag_file_list_no_ext), marker_gene_assignment_rate)

    return marker_gene_assignment_rate, marker_gene_assignment_accuracy, genome_with_right_16s_assignment_always, genome_without_right_16s_assignment


def rename_seq(ctg_file_in, ctg_file_out, prefix, str_connector):

    ctg_file_out_handle = open(ctg_file_out, 'w')
    for Seq_record in SeqIO.parse(ctg_file_in, 'fasta'):
        Seq_record.id = '%s%s%s' % (prefix, str_connector, Seq_record.id)
        SeqIO.write(Seq_record, ctg_file_out_handle, 'fasta')
    ctg_file_out_handle.close()


def SeqIO_convert_worker(argument_list):

    file_in         = argument_list[0]
    file_in_fmt     = argument_list[1]
    file_out        = argument_list[2]
    file_out_fmt    = argument_list[3]
    SeqIO.convert(file_in, file_in_fmt, file_out, file_out_fmt)


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


def get_sankey_plot(node_list, source_list, target_list, value_list, color_list, plot_title, plot_height, output_html):

    node_index_dict = {y: x for x, y in enumerate(node_list)}
    source_index = [node_index_dict[x] for x in source_list]
    target_index = [node_index_dict[x] for x in target_list]

    # https://anvil.works/docs/api/plotly.graph_objs.sankey
    fig = go.Figure(data=[go.Sankey(node=dict(label=node_list,  # line=0,
                                              pad=5,  # space between node
                                              thickness=12,  # node width
                                              line=dict(width=0)),  # set width of node border to 0
                                    link=dict(source=source_index,
                                              target=target_index,
                                              value=value_list,
                                              color=color_list))])

    fig.update_layout(autosize=False, width=1200, height=plot_height, margin=dict(l=50, r=50, b=50, t=125), paper_bgcolor="white", title=plot_title)
    fig.update_traces(textfont_size=11)
    fig.write_html(output_html)


def sankey_linkages(combined_linkage_file_ctg_level, linkage_plot_rd1_html, linkage_plot_rd2_html):

    dict_for_sankey_key_connector   = '___X___'

    linkage_num_dict_rd1 = {}
    linkage_num_dict_rd2 = {}
    ctg_to_gnm_dict_rd1 = {}
    ctg_to_gnm_dict_rd2 = {}
    node_set_rd1 = set()
    node_set_rd2 = set()
    genome_set_rd1 = set()
    genome_set_rd2 = set()
    contig_set_rd1 = set()
    contig_set_rd2 = set()
    marker_gene_set_rd1 = set()
    marker_gene_set_rd2 = set()
    for each_linkage in open(combined_linkage_file_ctg_level):
        if not each_linkage.startswith('Marker___Genome(total)\tContig\tRd1\tRd2'):
            each_linkage_split = each_linkage.strip().split('\t')

            marker_id = each_linkage_split[0].split('___')[0]
            gnm_id = each_linkage_split[0].split('___')[1].split('(')[0]
            ctg_id = each_linkage_split[1]
            total_link_num = int(each_linkage_split[2]) + int(each_linkage_split[3])
            marker_to_ctg_key = '%s%s%s' % (marker_id, dict_for_sankey_key_connector, ctg_id)
            ctg_to_gnm_key = '%s%s%s' % (ctg_id, dict_for_sankey_key_connector, gnm_id)

            if int(each_linkage_split[3]) == 0:
                genome_set_rd1.add(gnm_id)
                contig_set_rd1.add(ctg_id)
                marker_gene_set_rd1.add(marker_id)
                node_set_rd1.add(marker_id)
                node_set_rd1.add(ctg_id)
                node_set_rd1.add(gnm_id)

                if ctg_id not in ctg_to_gnm_dict_rd1:
                    ctg_to_gnm_dict_rd1[ctg_id] = gnm_id

                if marker_to_ctg_key not in linkage_num_dict_rd1:
                    linkage_num_dict_rd1[marker_to_ctg_key] = total_link_num
                else:
                    linkage_num_dict_rd1[marker_to_ctg_key] += total_link_num

                if ctg_to_gnm_key not in linkage_num_dict_rd1:
                    linkage_num_dict_rd1[ctg_to_gnm_key] = total_link_num
                else:
                    linkage_num_dict_rd1[ctg_to_gnm_key] += total_link_num

            if int(each_linkage_split[2]) == 0:
                genome_set_rd2.add(gnm_id)
                contig_set_rd2.add(ctg_id)
                marker_gene_set_rd2.add(marker_id)
                node_set_rd2.add(marker_id)
                node_set_rd2.add(ctg_id)
                node_set_rd2.add(gnm_id)
                if ctg_id not in ctg_to_gnm_dict_rd2:
                    ctg_to_gnm_dict_rd2[ctg_id] = gnm_id

                if marker_to_ctg_key not in linkage_num_dict_rd2:
                    linkage_num_dict_rd2[marker_to_ctg_key] = total_link_num
                else:
                    linkage_num_dict_rd2[marker_to_ctg_key] += total_link_num

                if ctg_to_gnm_key not in linkage_num_dict_rd2:
                    linkage_num_dict_rd2[ctg_to_gnm_key] = total_link_num
                else:
                    linkage_num_dict_rd2[ctg_to_gnm_key] += total_link_num

    source_list_rd1 = []
    target_list_rd1 = []
    value_list_rd1 = []
    for each_rd1_linkage in linkage_num_dict_rd1:
        each_rd1_linkage_split = each_rd1_linkage.split(dict_for_sankey_key_connector)
        source_list_rd1.append(each_rd1_linkage_split[0])
        target_list_rd1.append(each_rd1_linkage_split[1])
        value_list_rd1.append(linkage_num_dict_rd1[each_rd1_linkage])

    source_list_rd2 = []
    target_list_rd2 = []
    value_list_rd2 = []
    for each_rd2_linkage in linkage_num_dict_rd2:
        each_rd2_linkage_split = each_rd2_linkage.split(dict_for_sankey_key_connector)
        source_list_rd2.append(each_rd2_linkage_split[0])
        target_list_rd2.append(each_rd2_linkage_split[1])
        value_list_rd2.append(linkage_num_dict_rd2[each_rd2_linkage])

    gnm_color_list_rd1 = sns.color_palette('tab20', len(genome_set_rd1)).as_hex()
    gnm_color_list_rd2 = sns.color_palette('tab20', len(genome_set_rd2)).as_hex()

    genome_to_color_dict_rd1 = {gnm: color for gnm, color in zip(genome_set_rd1, gnm_color_list_rd1)}
    genome_to_color_dict_rd2 = {gnm: color for gnm, color in zip(genome_set_rd2, gnm_color_list_rd2)}

    color_list_rd1 = []
    for each_target in target_list_rd1:
        if each_target in genome_to_color_dict_rd1:
            color_list_rd1.append(genome_to_color_dict_rd1[each_target])
        else:
            target_genome = ctg_to_gnm_dict_rd1[each_target]
            color_list_rd1.append(genome_to_color_dict_rd1[target_genome])

    color_list_rd2 = []
    for each_target in target_list_rd2:
        if each_target in genome_to_color_dict_rd2:
            color_list_rd2.append(genome_to_color_dict_rd2[each_target])
        else:
            target_genome = ctg_to_gnm_dict_rd2[each_target]
            color_list_rd2.append(genome_to_color_dict_rd2[target_genome])

    node_list_rd1 = sorted([i for i in node_set_rd1])
    node_list_rd2 = sorted([i for i in node_set_rd2])

    plot_title_text_rd1 = 'MarkerMAG detected linkages (round 1)<br>Number of linked genomes: %s<br>Number of linked markers: %s' % (len(genome_set_rd1), len(marker_gene_set_rd1))
    plot_title_text_rd2 = 'MarkerMAG detected linkages (round 2)<br>Number of linked genomes: %s<br>Number of linked markers: %s' % (len(genome_set_rd2), len(marker_gene_set_rd2))

    plot_height_rd1 = 900 if max([len(contig_set_rd1), len(marker_gene_set_rd1)]) <= 25 else max([len(contig_set_rd1), len(marker_gene_set_rd1)]) * 32
    plot_height_rd2 = 900 if max([len(contig_set_rd2), len(marker_gene_set_rd2)]) <= 25 else max([len(contig_set_rd2), len(marker_gene_set_rd2)]) * 32

    plot_title_dict_rd1 = dict(text=plot_title_text_rd1, x=0.05, y=(1-(50/plot_height_rd1)))
    plot_title_dict_rd2 = dict(text=plot_title_text_rd2, x=0.05, y=(1-(50/plot_height_rd2)))

    get_sankey_plot(node_list_rd1, source_list_rd1, target_list_rd1, value_list_rd1, color_list_rd1, plot_title_dict_rd1, plot_height_rd1, linkage_plot_rd1_html)
    get_sankey_plot(node_list_rd2, source_list_rd2, target_list_rd2, value_list_rd2, color_list_rd2, plot_title_dict_rd2, plot_height_rd2, linkage_plot_rd2_html)


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


def remove_both_ends_clp(sam_in, sam_out):
    sam_out_handle = open(sam_out, 'w')
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)
        else:
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    sam_out_handle.write(each_read)
    sam_out_handle.close()


def keep_best_matches_in_sam(sam_in, sam_out):
    # get read_to_cigar_dict
    read_to_cigar_dict = {}
    for each_line in open(sam_in):
        each_line_split = each_line.strip().split('\t')
        if not each_line.startswith('@'):
            read_id = each_line_split[0]
            cigar = each_line_split[5]
            if cigar != '*':
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    if read_id not in read_to_cigar_dict:
                        read_to_cigar_dict[read_id] = {cigar}
                    else:
                        read_to_cigar_dict[read_id].add(cigar)

    # get min_mismatch for each read
    read_min_mismatch_dict = {}
    for each_read in read_to_cigar_dict:
        read_mismatch_set = set()
        for each_cigar in read_to_cigar_dict[each_read]:
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
            read_mismatch_set.add(mismatch_pct)
        read_min_mismatch = min(read_mismatch_set)
        read_min_mismatch_dict[each_read] = read_min_mismatch

    sam_file_best_match_handle = open(sam_out, 'w')
    for each_line in open(sam_in):
        if each_line.startswith('@'):
            sam_file_best_match_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            cigar = each_line_split[5]
            if cigar == '*':
                sam_file_best_match_handle.write(each_line)
            else:
                cigar_split = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_split)
                    if mismatch_pct <= (read_min_mismatch_dict[read_id] * 1.5):

                        sam_file_best_match_handle.write(each_line)
    sam_file_best_match_handle.close()


def keep_best_matches_in_sam_keep_short_M(sam_in, min_M_len, sam_out):

    # get read_to_cigar_dict
    read_to_cigar_dict = {}
    for each_line in open(sam_in):
        each_line_split = each_line.strip().split('\t')
        if not each_line.startswith('@'):
            read_id = each_line_split[0]
            cigar = each_line_split[5]
            if cigar != '*':
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    if read_id not in read_to_cigar_dict:
                        read_to_cigar_dict[read_id] = {cigar}
                    else:
                        read_to_cigar_dict[read_id].add(cigar)

    # get min_mismatch for each read
    read_min_mismatch_dict = {}
    for each_read in read_to_cigar_dict:
        read_mismatch_set_all_M = set()
        read_mismatch_set_long_M = set()
        for each_cigar in read_to_cigar_dict[each_read]:
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
            read_mismatch_set_all_M.add(mismatch_pct)
            if aligned_len >= min_M_len:
                read_mismatch_set_long_M.add(mismatch_pct)
        read_min_mismatch = min(read_mismatch_set_all_M)
        if len(read_mismatch_set_long_M) > 0:
            read_min_mismatch = min(read_mismatch_set_long_M)
        read_min_mismatch_dict[each_read] = read_min_mismatch

    sam_file_best_match_handle = open(sam_out, 'w')
    for each_line in open(sam_in):
        if each_line.startswith('@'):
            sam_file_best_match_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            cigar = each_line_split[5]
            if cigar == '*':
                sam_file_best_match_handle.write(each_line)
            else:
                read_id = each_line_split[0]
                cigar_split = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_split)
                    if mismatch_pct <= (read_min_mismatch_dict[read_id] * 1.5):
                        sam_file_best_match_handle.write(each_line)
    sam_file_best_match_handle.close()


def run_mira5(output_prefix, mira_tmp_dir, step_2_wd, mira_manifest, unpaired_fastq, mira_stdout, force_overwrite):

    # prepare manifest file
    mira_manifest_handle = open(mira_manifest, 'w')
    mira_manifest_handle.write('project = %s_mira_est_no_chimera\n' % output_prefix)
    mira_manifest_handle.write('job=est,denovo,accurate\n')
    mira_manifest_handle.write('parameters = -CL:ascdc\n')
    mira_manifest_handle.write('readgroup = SomeUnpairedIlluminaReadsIGotFromTheLab\n')
    mira_manifest_handle.write('data = %s\n' % os.path.abspath(unpaired_fastq))
    mira_manifest_handle.write('technology = solexa\n')
    mira_manifest_handle.close()

    if os.path.isdir(mira_tmp_dir) is False:
        os.mkdir(mira_tmp_dir)

    # run Mira
    mira_cmd = 'mira -c %s %s > %s' % (step_2_wd, os.path.abspath(mira_manifest), mira_stdout)
    if mira_tmp_dir is not None:
        mira_cmd = 'mira -c %s %s > %s' % (mira_tmp_dir, os.path.abspath(mira_manifest), mira_stdout)
    os.system(mira_cmd)

    # parse mira output
    if (mira_tmp_dir is not None) and (mira_tmp_dir != step_2_wd):
        os.system('cp -r %s/%s_mira_est_no_chimera_assembly %s/' % (mira_tmp_dir, output_prefix, step_2_wd))


def extract_reads_worker(argument_list):

    reads_file_in = argument_list[0]
    reads_fmt = argument_list[1]
    reads_to_extract = argument_list[2]
    reads_file_out = argument_list[3]

    reads_file_out_handle = open(reads_file_out, 'w')
    for read_record in SeqIO.parse(reads_file_in, reads_fmt):
        if read_record.id in reads_to_extract:
            if reads_fmt == 'fasta':
                reads_file_out_handle.write('>%s\n' % read_record.id)
                reads_file_out_handle.write('%s\n' % read_record.seq)
            if reads_fmt == 'fastq':
                SeqIO.write(read_record, reads_file_out_handle, 'fastq')
    reads_file_out_handle.close()


def get_GapFilling_stats_by_assembly(free_living_16s_ref_file,
                                     free_living_ctg_ref_file,
                                     mini_assembly_to_16s_reads,
                                     mini_assembly_to_ctg_reads,
                                     ctg_level_min_link,
                                     mini_assembly_to_16s_ctg_connector,
                                     gnm_to_ctg_connector,
                                     marker_to_ctg_gnm_Key_connector,
                                     max_within_cate_diff_pct,
                                     max_between_cate_diff_pct,
                                     stats_GapFilling_ctg,
                                     stats_GapFilling_gnm):

    round2_free_living_16s_ref_dict = {}
    for free_living_read_16s in open(free_living_16s_ref_file):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_16s_id = free_living_read_16s_split[0]
            read_16s_refs = free_living_read_16s_split[1].split(',')
            round2_free_living_16s_ref_dict[read_16s_id] = read_16s_refs

    round2_free_living_ctg_ref_dict = {}
    for free_living_read_ctg in open(free_living_ctg_ref_file):
        free_living_read_ctg_split = free_living_read_ctg.strip().split('\t')
        read_ctg_id = free_living_read_ctg_split[0]
        read_ctg_refs = free_living_read_ctg_split[1].split(',')
        read_ctg_refs_no_suffix = []
        for each_read_ctg_ref in read_ctg_refs:
            if each_read_ctg_ref[-2:] in ['_l', '_r']:
                each_read_ctg_ref_no_suffix = each_read_ctg_ref[:-2]
                read_ctg_refs_no_suffix.append(each_read_ctg_ref_no_suffix)
        round2_free_living_ctg_ref_dict[read_ctg_id] = read_ctg_refs_no_suffix

    mini_assembly_to_16s_dict = {}
    for each_mini_assembly in open(mini_assembly_to_16s_reads):
        mini_assembly_split = each_mini_assembly.strip().split('\t')
        mini_assembly_id = mini_assembly_split[0]
        mini_assembly_mapped_reads = mini_assembly_split[1].split(',')
        for each_mapped_read in mini_assembly_mapped_reads:
            mapped_read_16s_refs = round2_free_living_16s_ref_dict.get(each_mapped_read, [])
            for each_mapped_read_16s_ref in mapped_read_16s_refs:
                mini_assembly_to_16s_key = '%s%s%s' % (mini_assembly_id, mini_assembly_to_16s_ctg_connector, each_mapped_read_16s_ref)
                if mini_assembly_to_16s_key not in mini_assembly_to_16s_dict:
                    mini_assembly_to_16s_dict[mini_assembly_to_16s_key] = 1
                else:
                    mini_assembly_to_16s_dict[mini_assembly_to_16s_key] += 1

    mini_assembly_to_ctg_dict = {}
    for each_mini_assembly in open(mini_assembly_to_ctg_reads):
        mini_assembly_split = each_mini_assembly.strip().split('\t')
        mini_assembly_id = mini_assembly_split[0]
        mini_assembly_mapped_reads = mini_assembly_split[1].split(',')
        for each_mapped_read in mini_assembly_mapped_reads:
            mapped_read_ctg_refs = round2_free_living_ctg_ref_dict.get(each_mapped_read, [])
            for each_mapped_read_ctg_ref in mapped_read_ctg_refs:
                mini_assembly_to_ctg_key = '%s%s%s' % (mini_assembly_id, mini_assembly_to_16s_ctg_connector, each_mapped_read_ctg_ref)
                if mini_assembly_to_ctg_key not in mini_assembly_to_ctg_dict:
                    mini_assembly_to_ctg_dict[mini_assembly_to_ctg_key] = 1
                else:
                    mini_assembly_to_ctg_dict[mini_assembly_to_ctg_key] += 1

    mini_assembly_to_16s_dict_reformatted = {}
    max_link_nun_dict_16s = {}
    for each in mini_assembly_to_16s_dict:
        mini_assembly_id = each.split(mini_assembly_to_16s_ctg_connector)[0]
        seq_16s_id = each.split(mini_assembly_to_16s_ctg_connector)[1]
        linkage_num = mini_assembly_to_16s_dict[each]
        seq_16s_with_num = '%s__num__%s' % (seq_16s_id, linkage_num)
        if linkage_num >= ctg_level_min_link:

            # add to mini_assembly_to_16s_dict_reformatted
            if mini_assembly_id not in mini_assembly_to_16s_dict_reformatted:
                mini_assembly_to_16s_dict_reformatted[mini_assembly_id] = {seq_16s_with_num}
            else:
                mini_assembly_to_16s_dict_reformatted[mini_assembly_id].add(seq_16s_with_num)

            # add to max_link_nun_dict_16s
            if seq_16s_id not in max_link_nun_dict_16s:
                max_link_nun_dict_16s[seq_16s_id] = linkage_num
            else:
                if linkage_num > max_link_nun_dict_16s[seq_16s_id]:
                    max_link_nun_dict_16s[seq_16s_id] = linkage_num

    mini_assembly_to_ctg_dict_reformatted = {}
    max_link_nun_dict_ctg = {}
    for each in mini_assembly_to_ctg_dict:
        mini_assembly_id = each.split(mini_assembly_to_16s_ctg_connector)[0]
        ctg_id = each.split(mini_assembly_to_16s_ctg_connector)[1]
        linkage_num = mini_assembly_to_ctg_dict[each]
        ctg_with_num = '%s__num__%s' % (ctg_id, linkage_num)
        if linkage_num >= ctg_level_min_link:

            # add to mini_assembly_to_ctg_dict_reformatted
            if mini_assembly_id not in mini_assembly_to_ctg_dict_reformatted:
                mini_assembly_to_ctg_dict_reformatted[mini_assembly_id] = {ctg_with_num}
            else:
                mini_assembly_to_ctg_dict_reformatted[mini_assembly_id].add(ctg_with_num)

            # add to max_link_nun_dict_ctg
            if ctg_id not in max_link_nun_dict_ctg:
                max_link_nun_dict_ctg[ctg_id] = linkage_num
            else:
                if linkage_num > max_link_nun_dict_ctg[ctg_id]:
                    max_link_nun_dict_ctg[ctg_id] = linkage_num

    mini_assembly_linked_both = set(mini_assembly_to_16s_dict_reformatted).intersection(mini_assembly_to_ctg_dict_reformatted)

    stats_GapFilling_ctg_handle = open(stats_GapFilling_ctg, 'w')
    stats_GapFilling_gnm_dict = {}
    for each_mini_assembly in mini_assembly_linked_both:
        linked_16s = mini_assembly_to_16s_dict_reformatted[each_mini_assembly]
        linked_ctg = mini_assembly_to_ctg_dict_reformatted[each_mini_assembly]
        linked_16s_num_list = [int(i.split('__num__')[1]) for i in linked_16s]
        linked_ctg_num_list = [int(i.split('__num__')[1]) for i in linked_ctg]
        linked_16s_num_max = max(linked_16s_num_list)
        linked_ctg_num_max = max(linked_ctg_num_list)

        if (min(linked_16s_num_max, linked_ctg_num_max) * 100 / max(linked_16s_num_max, linked_ctg_num_max)) >= max_between_cate_diff_pct:

            linked_16s_filtered = [i for i in linked_16s if int(i.split('__num__')[1])*100/linked_16s_num_max >= max_within_cate_diff_pct]
            linked_ctg_filtered = [i for i in linked_ctg if int(i.split('__num__')[1])*100/linked_ctg_num_max >= max_within_cate_diff_pct]

            for each_linked_16s in linked_16s_filtered:
                linked_16s_id = each_linked_16s.split('__num__')[0]
                linked_16s_num = int(each_linked_16s.split('__num__')[1])
                linked_16s_num_pct_by_max = linked_16s_num * 100 / max_link_nun_dict_16s[linked_16s_id]

                for each_linked_ctg in linked_ctg_filtered:
                    linked_ctg_id = each_linked_ctg.split('__num__')[0]
                    linked_gnm_id = linked_ctg_id.split(gnm_to_ctg_connector)[0]
                    linked_ctg_num = int(each_linked_ctg.split('__num__')[1])
                    linked_ctg_num_pct_by_max = linked_ctg_num*100/max_link_nun_dict_ctg[linked_ctg_id]

                    if (linked_16s_num_pct_by_max >= 50) and (linked_ctg_num_pct_by_max >= 50):
                        stats_GapFilling_ctg_handle.write('%s\t%s\t%s\n' % (linked_16s_id, linked_ctg_id, (linked_16s_num + linked_ctg_num)))
                        marker_to_gnm_key = '%s%s%s' % (linked_16s_id, marker_to_ctg_gnm_Key_connector, linked_gnm_id)
                        if marker_to_gnm_key not in stats_GapFilling_gnm_dict:
                            stats_GapFilling_gnm_dict[marker_to_gnm_key] = (linked_16s_num + linked_ctg_num)
                        else:
                            stats_GapFilling_gnm_dict[marker_to_gnm_key] += (linked_16s_num + linked_ctg_num)
    stats_GapFilling_ctg_handle.close()

    stats_GapFilling_gnm_handle = open(stats_GapFilling_gnm, 'w')
    stats_GapFilling_gnm_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_16s_to_gnm in stats_GapFilling_gnm_dict:
        each_16s_to_gnm_split = each_16s_to_gnm.split(marker_to_ctg_gnm_Key_connector)
        id_16s = each_16s_to_gnm_split[0]
        id_gnm = each_16s_to_gnm_split[1]
        linkage_num = stats_GapFilling_gnm_dict[each_16s_to_gnm]
        stats_GapFilling_gnm_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (id_16s, id_gnm, linkage_num))
    stats_GapFilling_gnm_handle.close()


def get_unmapped_mates_seq(sam_file, input_r1_fasta, input_r2_fasta, extracted_seq_file):

    sam_path, sam_basename, sam_ext = sep_path_basename_ext(sam_file)
    reads_to_extract_r1_txt = '%s/%s_unmapped_mates_R1.txt' % (sam_path, sam_basename)
    reads_to_extract_r2_txt = '%s/%s_unmapped_mates_R2.txt' % (sam_path, sam_basename)
    reads_to_extract_r1_fa  = '%s/%s_unmapped_mates_R1.fa'  % (sam_path, sam_basename)
    reads_to_extract_r2_fa  = '%s/%s_unmapped_mates_R2.fa'  % (sam_path, sam_basename)

    prescreening_qualified_reads_set = set()
    prescreening_qualified_reads_base_set = set()
    for each_read in open(sam_file):
        each_read_split = each_read.strip().split('\t')
        if not each_read.startswith('@'):
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            prescreening_qualified_reads_set.add(read_id)
            prescreening_qualified_reads_base_set.add(read_id_base)

    # get id files
    reads_to_extract_r1 = set()
    reads_to_extract_r2 = set()
    for each_read_base in prescreening_qualified_reads_base_set:
        read_r1 = '%s.1' % each_read_base
        read_r2 = '%s.2' % each_read_base
        if (read_r1 not in prescreening_qualified_reads_set) and (read_r2 in prescreening_qualified_reads_set):
            reads_to_extract_r1.add(read_r1)
        if (read_r1 in prescreening_qualified_reads_set) and (read_r2 not in prescreening_qualified_reads_set):
            reads_to_extract_r2.add(read_r2)

    reads_to_extract_r1_txt_handle = open(reads_to_extract_r1_txt, 'w')
    reads_to_extract_r1_txt_handle.write('%s\n' % '\n'.join(reads_to_extract_r1))
    reads_to_extract_r1_txt_handle.close()

    reads_to_extract_r2_txt_handle = open(reads_to_extract_r2_txt, 'w')
    reads_to_extract_r2_txt_handle.write('%s\n' % '\n'.join(reads_to_extract_r2))
    reads_to_extract_r2_txt_handle.close()

    seqtk_extract_r1_read_cmd = 'seqtk subseq %s %s > %s' % (
    input_r1_fasta, reads_to_extract_r1_txt, reads_to_extract_r1_fa)
    seqtk_extract_r2_read_cmd = 'seqtk subseq %s %s > %s' % (
    input_r2_fasta, reads_to_extract_r2_txt, reads_to_extract_r2_fa)
    os.system(seqtk_extract_r1_read_cmd)
    os.system(seqtk_extract_r2_read_cmd)

    os.system('cat %s %s > %s' % (reads_to_extract_r1_fa, reads_to_extract_r2_fa, extracted_seq_file))

    # rm tmp files
    os.system('rm %s' % reads_to_extract_r1_txt)
    os.system('rm %s' % reads_to_extract_r2_txt)
    os.system('rm %s' % reads_to_extract_r1_fa)
    os.system('rm %s' % reads_to_extract_r2_fa)


def polish_16s(file_in, file_out_ffn):

    file_out_path, file_out_base, file_out_ext = sep_path_basename_ext(file_out_ffn)

    barrnap_stdout   = '%s/%s.log'    % (file_out_path, file_out_base)
    file_out_gff     = '%s/%s.gff'    % (file_out_path, file_out_base)
    file_out_ffn_tmp = '%s/%s_tmp%s' % (file_out_path, file_out_base, file_out_ext)

    barrnap_cmd = 'barrnap --quiet -o %s %s 2> %s > %s' % (file_out_ffn_tmp, file_in, barrnap_stdout, file_out_gff)
    os.system(barrnap_cmd)

    wrote_id = []
    file_out_ffn_handle = open(file_out_ffn, 'w')
    for each_16s in SeqIO.parse(file_out_ffn_tmp, 'fasta'):
        seq_id = each_16s.id
        if seq_id.startswith('16S_rRNA::'):
            seq_id_polished = seq_id[10:].split(':')[0]

            if seq_id_polished not in wrote_id:
                file_out_ffn_handle.write('>%s\n' % seq_id_polished)
                file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
                wrote_id.append(seq_id_polished)
            else:
                file_out_ffn_handle.write('>%s_%s\n' % (seq_id_polished, (wrote_id.count(seq_id_polished) + 1)))
                file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
                wrote_id.append(seq_id_polished)

    file_out_ffn_handle.close()

    #os.system('rm %s' % file_out_ffn_tmp)
    #os.system('rm %s.fai' % file_in)


def get_min_max_cigar_S(cigar_list):

    # get cigar_S_set
    cigar_S_set = set()
    for each_cigar in cigar_list:
        if ('S' not in each_cigar) and ('s' not in each_cigar):
            cigar_S_set.add(0)
        else:
            cigar_splitted = cigar_splitter(each_cigar)

            if cigar_splitted[0][-1] in ['S', 's']:
                cigar_S_set.add(int(cigar_splitted[0][:-1]))

            if cigar_splitted[-1][-1] in ['S', 's']:
                cigar_S_set.add(int(cigar_splitted[1][:-1]))

    # get min_cigar_S and max_cigar_S
    min_cigar_S = 0
    max_cigar_S = 0
    if len(cigar_S_set) > 0:
        min_cigar_S = min(cigar_S_set)
        max_cigar_S = max(cigar_S_set)

    return min_cigar_S, max_cigar_S


def get_ctg_mean_depth_by_samtools_coverage_global(index_ref, ref_seq, reads_r1, reads_r2, reads_unpaired, subsample_rate, num_threads):

    ref_seq_file_path, ref_seq_file_basename, ref_seq_file_extension = sep_path_basename_ext(ref_seq)

    sam_file                                      = '%s/%s.sam'                                     % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_reformatted                          = '%s/%s_reformatted.sam'                         % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_reformatted_log                      = '%s/%s_reformatted.log'                         % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_reformatted_best_match               = '%s/%s_reformatted_best_match.sam'              % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_reformatted_best_match_low_mismatch  = '%s/%s_reformatted_best_match_low_mismatch.sam' % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_sorted                               = '%s/%s_sorted.sam'                              % (ref_seq_file_path, ref_seq_file_basename)
    coverage_file                                 = '%s/%s_cov.txt'                                 % (ref_seq_file_path, ref_seq_file_basename)

    # build reference index
    cmd_bowtie2_build   = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, ref_seq, ref_seq_file_path, ref_seq_file_basename)
    if index_ref is True:
        os.system(cmd_bowtie2_build)

    # mapping
    cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -U %s,%s -S %s -p %s --all --no-unal -N 1 -L 30 -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
    os.system(cmd_bowtie2_mapping)

    # filter mapping
    bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file, sam_file_reformatted, sam_file_reformatted_log)
    os.system(bbmap_reformat_cmd)
    keep_best_matches_in_sam_keep_short_M(sam_file_reformatted, 35, sam_file_reformatted_best_match)
    remove_high_mismatch(sam_file_reformatted_best_match, 2, sam_file_reformatted_best_match_low_mismatch)

    # sort mapping
    cmd_samtools_sort = 'samtools sort %s -o %s' % (sam_file_reformatted_best_match_low_mismatch, sam_file_sorted)
    os.system(cmd_samtools_sort)

    # get mean depth
    cmd_samtools_coverage = 'samtools coverage --ff 4 %s -o %s' % (sam_file_sorted, coverage_file)
    os.system(cmd_samtools_coverage)

    # remove sam files
    os.system('rm %s' % sam_file)
    # os.system('rm %s' % sam_file_sorted)

    # store mean depth into dict
    mean_depth_dict_ctg = {}
    ctg_len_dict = {}
    for each_ctg_depth in open(coverage_file):
        if not each_ctg_depth.startswith('#'):
            ctg_depth_split = each_ctg_depth.strip().split('\t')
            ctg_id = ctg_depth_split[0]
            ctg_len = int(ctg_depth_split[2])
            ctg_depth = float(ctg_depth_split[6]) * (1 / subsample_rate)
            mean_depth_dict_ctg[ctg_id] = ctg_depth
            ctg_len_dict[ctg_id] = ctg_len

    return mean_depth_dict_ctg, ctg_len_dict


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


def linkage_vis_worker(arguments_list):

    reads_file_base             = arguments_list[0]
    mafft_seq_folder            = arguments_list[1]
    marker_seq                  = arguments_list[2]
    contig_seq                  = arguments_list[3]
    end_ctg_len_for_mafft       = arguments_list[4]
    gap_N_num                   = arguments_list[5]
    bowtie_parameter            = arguments_list[6]
    marker_pos_list             = arguments_list[7]
    contig_pos_list             = arguments_list[8]

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
    if contig_pos_median <= 500:
        linked_end_contig = 'left'
    elif (contig_len - contig_pos_median) <= 500:
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
            concatenated_seq_id = 'Marker_NNN_Contig'
            concatenated_seq = '%s%s%s' % (marker_seq, 'N' * gap_N_num, contig_seq_for_mafft)
            to_concatenate = True
            concatenate_pos = len(marker_seq) + round(gap_N_num / 2)

        if (linked_end_marker == 'left') and (linked_end_contig == 'right'):
            concatenated_seq_id = 'Contig_NNN_Marker'
            concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N' * gap_N_num, marker_seq)
            to_concatenate = True
            concatenate_pos = len(contig_seq_for_mafft) + round(gap_N_num / 2)

        if (linked_end_marker == 'left') and (linked_end_contig == 'left'):
            marker_seq_rc = get_rc(marker_seq)
            concatenated_seq_id = 'Marker_RC_NNN_Contig'
            concatenated_seq = '%s%s%s' % (marker_seq_rc, 'N' * gap_N_num, contig_seq_for_mafft)
            to_concatenate = True
            concatenate_pos = len(marker_seq_rc) + round(gap_N_num / 2)

        if (linked_end_marker == 'right') and (linked_end_contig == 'right'):
            marker_seq_rc = get_rc(marker_seq)
            concatenated_seq_id = 'Contig_NNN_Marker_RC'
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
        pwd_seq_file_16s_handle.write('>Marker\n')
        pwd_seq_file_16s_handle.write('%s\n' % marker_seq)
        pwd_seq_file_16s_handle.close()
        # write out ctg sequence
        pwd_seq_file_ctg_handle = open(pwd_seq_file_ctg, 'w')
        pwd_seq_file_ctg_handle.write('>Contig\n')
        pwd_seq_file_ctg_handle.write('%s\n' % contig_seq_for_mafft)
        pwd_seq_file_ctg_handle.close()

    ########## mapping ##########

    pwd_seq_file_cbd            = '%s/%s/%s_cbd.fa'         % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s            = '%s/%s/%s_16s.fa'         % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg            = '%s/%s/%s_ctg.fa'         % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_reads_r1       = '%s/%s/%s_R1.fa'          % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_reads_r2       = '%s/%s/%s_R2.fa'          % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_index      = '%s/%s/%s_cbd'            % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_index      = '%s/%s/%s_16s'            % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_index      = '%s/%s/%s_ctg'            % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_sam        = '%s/%s/%s_cbd.sam'        % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_sam        = '%s/%s/%s_16s.sam'        % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_sam        = '%s/%s/%s_ctg.sam'        % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_sam_log    = '%s/%s/%s_cbd.log'        % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_sam_log    = '%s/%s/%s_16s.log'        % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_sam_log    = '%s/%s/%s_ctg.log'        % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_cbd_Tablet_xml = '%s/%s/%s_cbd.tablet'     % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s_Tablet_xml = '%s/%s/%s_16s.tablet'     % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg_Tablet_xml = '%s/%s/%s_ctg.tablet'     % (mafft_seq_folder, reads_file_base, reads_file_base)

    if to_concatenate is True:
        index_ref_cmd = 'bowtie2-build --quiet -f %s %s' % (pwd_seq_file_cbd, pwd_seq_file_cbd_index)
        bowtie2_cmd = 'bowtie2 -x %s -U %s,%s -S %s -p 1 -f %s 2> %s' % (pwd_seq_file_cbd_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_cbd_sam, bowtie_parameter, pwd_seq_file_cbd_sam_log)
        os.system(index_ref_cmd)
        os.system(bowtie2_cmd)

        # write out Tablet xml file
        pwd_seq_file_cbd_Tablet_xml_handle = open(pwd_seq_file_cbd_Tablet_xml, 'w')
        pwd_seq_file_cbd_Tablet_xml_handle.write('<tablet>\n')
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <assembly>%s_cbd.sam</assembly>\n'    % reads_file_base)
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <reference>%s_cbd.fa</reference>\n'   % reads_file_base)
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <contig>%s</contig>\n'                % concatenated_seq_id)
        pwd_seq_file_cbd_Tablet_xml_handle.write('        <position>%s</position>\n'            % concatenate_pos)
        pwd_seq_file_cbd_Tablet_xml_handle.write('</tablet>\n')
        pwd_seq_file_cbd_Tablet_xml_handle.close()
    else:
        index_ref_cmd_16s = 'bowtie2-build --quiet -f %s %s' % (pwd_seq_file_16s, pwd_seq_file_16s_index)
        index_ref_cmd_ctg = 'bowtie2-build --quiet -f %s %s' % (pwd_seq_file_ctg, pwd_seq_file_ctg_index)
        os.system(index_ref_cmd_16s)
        os.system(index_ref_cmd_ctg)
        bowtie2_cmd_16s = 'bowtie2 -x %s -U %s,%s -S %s -p 6 -f %s 2> %s' % (pwd_seq_file_16s_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_16s_sam, bowtie_parameter, pwd_seq_file_16s_sam_log)
        bowtie2_cmd_ctg = 'bowtie2 -x %s -U %s,%s -S %s -p 6 -f %s 2> %s' % (pwd_seq_file_ctg_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_ctg_sam, bowtie_parameter, pwd_seq_file_ctg_sam_log)
        os.system(bowtie2_cmd_16s)
        os.system(bowtie2_cmd_ctg)

        # write out Tablet xml file
        pwd_seq_file_16s_Tablet_xml_handle = open(pwd_seq_file_16s_Tablet_xml, 'w')
        pwd_seq_file_16s_Tablet_xml_handle.write('<tablet>\n')
        pwd_seq_file_16s_Tablet_xml_handle.write('        <assembly>%s_16s.sam</assembly>\n'  % reads_file_base)
        pwd_seq_file_16s_Tablet_xml_handle.write('        <reference>%s_16s.fa</reference>\n' % reads_file_base)
        pwd_seq_file_16s_Tablet_xml_handle.write('        <contig>Marker</contig>\n')
        pwd_seq_file_16s_Tablet_xml_handle.write('</tablet>\n')
        pwd_seq_file_16s_Tablet_xml_handle.close()

        pwd_seq_file_ctg_Tablet_xml_handle = open(pwd_seq_file_ctg_Tablet_xml, 'w')
        pwd_seq_file_ctg_Tablet_xml_handle.write('<tablet>\n')
        pwd_seq_file_ctg_Tablet_xml_handle.write('        <assembly>%s_ctg.sam</assembly>\n'  % reads_file_base)
        pwd_seq_file_ctg_Tablet_xml_handle.write('        <reference>%s_ctg.fa</reference>\n' % reads_file_base)
        pwd_seq_file_ctg_Tablet_xml_handle.write('        <contig>Contig</contig>\n')
        pwd_seq_file_ctg_Tablet_xml_handle.write('</tablet>\n')
        pwd_seq_file_ctg_Tablet_xml_handle.close()

    # remove tmp files
    os.system('rm %s/%s/%s*.bt2' % (mafft_seq_folder, reads_file_base, reads_file_base))


def parse_sam_gnm_worker(arguments_list):

    rd1_unlinked_mags_sam_bowtie_reformat_sorted    = arguments_list[0]
    free_living_ctg_ref_file                        = arguments_list[1]
    min_M_len_ctg                                   = arguments_list[2]
    mismatch_cutoff                                 = arguments_list[3]
    round_2_ctg_end_seq_len_dict                    = arguments_list[4]
    rd2_with_both_mates                             = arguments_list[5]

    free_living_ctg_ref_file_handle = open(free_living_ctg_ref_file, 'w')
    current_read_base = ''
    current_read_base_r1_ctg_ref_dict_rd2 = dict()
    current_read_base_r2_ctg_ref_dict_rd2 = dict()
    with open(rd1_unlinked_mags_sam_bowtie_reformat_sorted) as rd1_unlinked_mags_sam_bowtie_reformat_sorted_opened:
        for each_line in rd1_unlinked_mags_sam_bowtie_reformat_sorted_opened:
            if not each_line.startswith('@'):
                each_line_split = each_line.strip().split('\t')
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
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}

                elif read_id_base == current_read_base:

                    if cigar != '*':
                        if read_strand == '1':
                            if ref_id not in current_read_base_r1_ctg_ref_dict_rd2:
                                current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r1_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar
                        if read_strand == '2':
                            if ref_id not in current_read_base_r2_ctg_ref_dict_rd2:
                                current_read_base_r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r2_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar
                else:
                    ################################### analysis previous read refs ####################################

                    ctg_refs_to_ignore_rd2 = set()

                    ########## get lowest mismatch for r1/r2 ctg refs ##########

                    # get r1_ref_cigar_set
                    r1_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r1_ctg_ref_dict_rd2.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r1_ref_cigar_set.update(each_pos_dict_values)

                    # get r2_ref_cigar_set
                    r2_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r2_ctg_ref_dict_rd2.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r2_ref_cigar_set.update(each_pos_dict_values)

                    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_ctg)
                    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_ctg)

                    ########## filter r1 ctg refs ##########
                    r1_ctg_refs_passed_qc = {}
                    for r1_ctg_ref_rd2 in current_read_base_r1_ctg_ref_dict_rd2:
                        r1_matched_pos_dict = current_read_base_r1_ctg_ref_dict_rd2[r1_ctg_ref_rd2]
                        if len(r1_matched_pos_dict) > 1:
                            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                        else:
                            r1_ctg_ref_pos = list(r1_matched_pos_dict.keys())[0]
                            r1_ctg_ref_cigar = r1_matched_pos_dict[r1_ctg_ref_pos]
                            r1_ctg_ref_cigar_splitted = cigar_splitter(r1_ctg_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r1_ctg_ref_cigar_splitted)
                            if both_end_clp is True:
                                ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                            else:
                                # check mismatch
                                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                    r1_ctg_ref_cigar_splitted)
                                if r1_ref_min_mismatch == 'NA':
                                    ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                                    ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                else:
                                    # check aligned length
                                    if r1_aligned_len < min_M_len_ctg:
                                        ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r1_ctg_ref_cigar) or ('s' in r1_ctg_ref_cigar):
                                            clip_in_middle = True
                                            if (r1_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r1_ctg_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r1_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r1_ctg_ref_pos + r1_aligned_len - 1) == \
                                                        round_2_ctg_end_seq_len_dict[
                                                            r1_ctg_ref_rd2]:
                                                    clip_in_middle = False

                                        if clip_in_middle is True:
                                            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                        else:
                                            r1_ctg_refs_passed_qc[r1_ctg_ref_rd2] = [r1_ctg_ref_cigar]

                    ########## filter r2 ctg refs ##########
                    r2_ctg_refs_passed_qc = {}
                    for r2_ctg_ref_rd2 in current_read_base_r2_ctg_ref_dict_rd2:
                        r2_matched_pos_dict = current_read_base_r2_ctg_ref_dict_rd2[r2_ctg_ref_rd2]
                        if len(r2_matched_pos_dict) > 1:
                            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                        else:
                            r2_ctg_ref_pos = list(r2_matched_pos_dict.keys())[0]
                            r2_ctg_ref_cigar = r2_matched_pos_dict[r2_ctg_ref_pos]
                            r2_ctg_ref_cigar_splitted = cigar_splitter(r2_ctg_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
                            if both_end_clp is True:
                                ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                            else:
                                # check mismatch
                                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                    r2_ctg_ref_cigar_splitted)
                                if r2_ref_min_mismatch == 'NA':
                                    ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                                    ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                else:
                                    # check aligned length
                                    if r2_aligned_len < min_M_len_ctg:
                                        ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r2_ctg_ref_cigar) or ('s' in r2_ctg_ref_cigar):
                                            clip_in_middle = True
                                            if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r2_ctg_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r2_ctg_ref_pos + r2_aligned_len - 1) == \
                                                        round_2_ctg_end_seq_len_dict[r2_ctg_ref_rd2]:
                                                    clip_in_middle = False

                                        if clip_in_middle is True:
                                            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                        else:
                                            r2_ctg_refs_passed_qc[r2_ctg_ref_rd2] = [r2_ctg_ref_cigar]

                    ####################################################################################################

                    r1_ctg_refs_rd2_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in ctg_refs_to_ignore_rd2}
                    r2_ctg_refs_rd2_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in ctg_refs_to_ignore_rd2}

                    # only r1 has no_ignored alignments
                    if (len(r1_ctg_refs_rd2_no_ignored) > 0) and (len(r2_ctg_refs_rd2_no_ignored) == 0):
                        if rd2_with_both_mates is True:
                            free_living_ctg_ref_file_handle.write('%s\t%s\n' % (current_read_base, ','.join(r1_ctg_refs_rd2_no_ignored)))
                        else:
                            free_living_ctg_ref_file_handle.write('%s.2\t%s\n' % (current_read_base, ','.join(r1_ctg_refs_rd2_no_ignored)))

                    # only r2 has no_ignored alignments
                    if (len(r1_ctg_refs_rd2_no_ignored) == 0) and (len(r2_ctg_refs_rd2_no_ignored) > 0):
                        if rd2_with_both_mates is True:
                            free_living_ctg_ref_file_handle.write('%s\t%s\n' % (current_read_base, ','.join(r2_ctg_refs_rd2_no_ignored)))
                        else:
                            free_living_ctg_ref_file_handle.write('%s.1\t%s\n' % (current_read_base, ','.join(r2_ctg_refs_rd2_no_ignored)))


                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_ctg_ref_dict_rd2 = dict()
                    current_read_base_r2_ctg_ref_dict_rd2 = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}

    free_living_ctg_ref_file_handle.close()


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

            # print(each_line_split)
            # print('%s\t%s\t%s\t%s\t%s' % (ctg_id, ctg_len, len_16s, left_gap, right_gap))
            # print()

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


def link_16s(args):

    ###################################################### file in/out #####################################################

    # file in
    output_prefix                       = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    mag_folder                          = args['mag']
    mag_file_extension                  = args['x']
    marker_gene_seqs                    = args['marker']
    min_iden_16s                        = args['min_iden_16s']
    min_cov_16s                         = args['min_cov_16s']
    min_aln_16s                         = args['min_aln_16s']
    min_link_num                        = args['min_link']
    num_threads                         = args['t']
    keep_quiet                          = args['quiet']
    force_overwrite                     = args['force']
    keep_temp                           = args['tmp']
    test_mode                           = args['test_mode']
    bbmap_memory                        = args['bbmap_mem']
    mismatch_cutoff                     = args['mismatch']
    min_M_pct                           = args['min_M_pct']
    within_gnm_linkage_num_diff         = args['link_num_diff']
    min_M_len_16s                       = args['min_M_len_16s']
    min_M_len_ctg                       = args['min_M_len_ctg']
    reads_vs_16s_sam                    = args['sorted_sam16s']
    no_polish                           = args['no_polish']

    # depth related
    min_16s_gnm_multiple                = args['depth_ratio']
    depth_file_16s                      = args['depth_16s']
    depth_file_mag                      = args['depth_mag']

    # by assembly
    round_2_mira                        = args['mira']
    mira_tmp_dir                        = args['mira_tmp']
    clp_read_for_assembly               = args['assemble_clp']
    max_mini_assembly_link_num_diff_between_ctg_16s = args['link_bias_rd2']

    pwd_makeblastdb_exe                 = 'makeblastdb'
    pwd_blastn_exe                      = 'blastn'
    pwd_bowtie2_build_exe               = 'bowtie2-build'
    pwd_bowtie2_exe                     = 'bowtie2'
    pwd_samtools_exe                    = 'samtools'
    pwd_bbmap_exe                       = 'bbmap.sh'
    pwd_spades_exe                      = 'spades.py'
    seqtk_exe                           = 'seqtk'

    marker_to_ctg_gnm_Key_connector                 = '___M___'
    gnm_to_ctg_connector                            = '___C___'
    mini_assembly_to_16s_ctg_connector              = '___Mini___'
    read_to_marker_connector                        = '___r___'
    end_seq_len                                     = 500
    ctg_level_min_link                              = 3
    end_ctg_len_for_mafft                           = 1000
    keep_short_M                                    = True
    gap_N_num                                       = 50
    report_interval                                 = 25000
    clp_pct_ctg_side_max_num                        = 65
    clp_pct_ratio_cutoff                            = 3.5
    # mismatch cutoff for filtering matches between unmapped mates and clipping reads against contig end
    mismatch_ctg_ends = 1
    subsample_rate_for_depth_estimation = 0.1  # between 0 and 1
    min_M_len_mini = min_M_len_ctg


    ################################################ check dependencies ################################################

    # check whether executables exist
    program_list = [pwd_makeblastdb_exe, pwd_blastn_exe, pwd_bowtie2_build_exe, pwd_bowtie2_exe, pwd_samtools_exe]
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()

    ################################################# check input files ################################################

    if os.path.isfile(marker_gene_seqs) is False:
        print('%s not found, program exited!' % os.path.basename(marker_gene_seqs))
        exit()

    # get input mag file list
    mag_file_re = '%s/*%s' % (mag_folder, mag_file_extension)
    mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
    if len(mag_file_list) == 0:
        print('No MAG detected, program exited!')
        exit()

    ############################################# create working directory #############################################

    # create working directory
    working_directory = '%s_MarkerMAG_wd' % output_prefix
    pwd_log_file      = '%s/%s.log'       % (working_directory, output_prefix)

    if (os.path.isdir(working_directory) is True) and (force_overwrite is False):
        print('Working directory detected, program exited!')
        exit()
    else:
        force_create_folder(working_directory)

    step_1_wd = '%s/%s_step_1_wd' % (working_directory, output_prefix)
    step_2_wd = '%s/%s_step_2_wd' % (working_directory, output_prefix)
    os.mkdir(step_1_wd)

    ############################################## check input reads format ############################################

    r1_path, r1_basename, r1_ext = sep_path_basename_ext(reads_file_r1)
    r2_path, r2_basename, r2_ext = sep_path_basename_ext(reads_file_r2)

    reads_file_r1_fasta = reads_file_r1
    reads_file_r2_fasta = reads_file_r2
    if ('q' in r1_ext) and ('q' in r2_ext):

        reads_file_r1_fasta_to_check = '%s/%s.fasta' % (r1_path, r1_basename)
        reads_file_r2_fasta_to_check = '%s/%s.fasta' % (r2_path, r2_basename)

        if (os.path.isfile(reads_file_r1_fasta_to_check) is True) and (os.path.isfile(reads_file_r2_fasta_to_check) is True):
            reads_file_r1_fasta = reads_file_r1_fasta_to_check
            reads_file_r2_fasta = reads_file_r2_fasta_to_check

        else:
            reads_file_r1_fasta = '%s/%s.fasta' % (step_1_wd, r1_basename)
            reads_file_r2_fasta = '%s/%s.fasta' % (step_1_wd, r2_basename)

            if num_threads >= 2:
                num_threads_SeqIO_convert_worker = 2
            else:
                num_threads_SeqIO_convert_worker = 1

            pool = mp.Pool(processes=num_threads_SeqIO_convert_worker)
            pool.map(SeqIO_convert_worker, [[reads_file_r1, 'fastq', reads_file_r1_fasta, 'fasta-2line'], [reads_file_r2, 'fastq', reads_file_r2_fasta, 'fasta-2line']])
            pool.close()
            pool.join()


    ################################################ prepare preset parameters to use ################################################

    # get reads_num, read_len and total len
    # paired_reads_num, read_len_median, read_len_max = get_read_num_and_length(reads_file_r1_fasta, working_directory, seqtk_exe)
    # estimated_total_read_len_gbp = (read_len_median*paired_reads_num*2)/(1024*1024*1024)
    # estimated_total_read_len_gbp = float("{0:.1f}".format(estimated_total_read_len_gbp))

    report_and_log(('mismatch_cutoff:\t%s%s'                % (mismatch_cutoff, '%')), pwd_log_file, keep_quiet)
    report_and_log(('min_M_len_16s:\t%sbp'                  % min_M_len_16s), pwd_log_file, keep_quiet)
    report_and_log(('min_M_len_ctg:\t%sbp'                  % min_M_len_ctg), pwd_log_file, keep_quiet)
    report_and_log(('min_M_pct:\t%s%s'                      % (min_M_pct, '%')), pwd_log_file, keep_quiet)
    report_and_log(('min_link_num:\t%s'                     % min_link_num), pwd_log_file, keep_quiet)
    report_and_log(('ctg_level_min_link:\t%s'               % (ctg_level_min_link)), pwd_log_file, keep_quiet)
    report_and_log(('end_seq_len:\t%sbp'                    % (end_seq_len)), pwd_log_file, keep_quiet)
    report_and_log(('max_mini_assembly_link_num_diff_between_ctg_16s:\t%s%s'   % (max_mini_assembly_link_num_diff_between_ctg_16s, '%')), pwd_log_file, keep_quiet)
    #report_and_log(('Number of paired reads: %s'            % paired_reads_num), pwd_log_file, keep_quiet)
    #report_and_log(('Read length median: %sbp'              % read_len_median), pwd_log_file, keep_quiet)
    #report_and_log(('Read length max: %sbp'                 % read_len_max), pwd_log_file, keep_quiet)
    #report_and_log(('Estimated total length: %sGbp'         % estimated_total_read_len_gbp), pwd_log_file, keep_quiet)


    ######################## check genomic sequence type and prepare files for making blast db #########################

    input_mag_folder_no_path        = mag_folder.split('/')[-1]
    mag_folder_in_wd                = '%s/input_MAGs'                       % step_1_wd
    prefixed_mag_folder             = '%s/%s_prefixed'                      % (mag_folder_in_wd, input_mag_folder_no_path)
    combined_input_gnms             = '%s/%s_combined.fa'                   % (mag_folder_in_wd, input_mag_folder_no_path)
    barrnap_wd                      = '%s/input_MAGs/barrnap_wd'            % step_1_wd
    combined_barrnap_gff            = '%s/input_MAGs/combined_barrnap.gff'  % step_1_wd

    # create folder
    os.mkdir(mag_folder_in_wd)
    os.mkdir(prefixed_mag_folder)

    # add mag id to its sequences
    argument_list_for_barrnap = []
    for mag_in in mag_file_list:
        mag_basename    = '.'.join(mag_in.split('.')[:-1])
        pwd_mag_in      = '%s/%s' % (mag_folder, mag_in)
        pwd_mag_renamed = '%s/%s' % (prefixed_mag_folder, mag_in)
        pwd_barrnap_ffn = '%s/%s.ffn' % (barrnap_wd, mag_basename)
        pwd_barrnap_gff = '%s/%s.gff' % (barrnap_wd, mag_basename)
        pwd_barrnap_log = '%s/%s.log' % (barrnap_wd, mag_basename)
        barrnap_cmd = 'barrnap --quiet -o %s %s > %s 2> %s' % (pwd_barrnap_ffn, pwd_mag_renamed, pwd_barrnap_gff, pwd_barrnap_log)
        argument_list_for_barrnap.append(barrnap_cmd)
        rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename, gnm_to_ctg_connector)

    # combine prefixed MAGs
    os.system('cat %s/*%s > %s' % (prefixed_mag_folder, mag_file_extension, combined_input_gnms))

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


    ########################################### define folder and file name ############################################

    input_16s_folder_in_wd                      = '%s/input_16S'                                           % step_1_wd
    input_reads_to_16s_sam_bowtie               = '%s/%s_input_reads_to_16S_bowtie.sam'                    % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_bowtie_log           = '%s/%s_input_reads_to_16S_bowtie.log'                    % (step_1_wd, output_prefix)
    input_reads_to_16s_sam                      = '%s/%s_input_reads_to_16S_reformatted.sam'               % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_reformat_log         = '%s/%s_input_reads_to_16S_reformat.log'                  % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_sorted               = '%s/%s_input_reads_to_16S_reformatted_sorted.sam'        % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_sorted_line_num      = '%s/%s_input_reads_to_16S_reformatted_sorted_lines.txt'  % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_sorted_split_folder  = '%s/%s_input_reads_to_16S_reformatted_sorted_split'      % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_MappingRecord_folder = '%s/%s_input_reads_to_16S_MappingRecord'                 % (step_1_wd, output_prefix)
    blast_results_all_vs_all_16s                = '%s/%s_16S_all_vs_all_blastn.tab'                        % (step_1_wd, output_prefix)
    pairwise_marker_similarity                  = '%s/%s_pairwise_marker_similarity.txt'                   % (step_1_wd, output_prefix)
    depth_file_ctg                              = '%s/%s_mean_depth_ctg.txt'                               % (step_1_wd, output_prefix)
    depth_file_gnm                              = '%s/%s_mean_depth_gnm.txt'                               % (step_1_wd, output_prefix)
    depth_file_16s_calculated                   = '%s/%s_mean_depth_16s.txt'                               % (step_1_wd, output_prefix)
    link_stats_combined                         = '%s/%s_stats_combined.txt'                               % (step_1_wd, output_prefix)
    link_stats_combined_sorted                  = '%s/%s_stats_combined_sorted.txt'                        % (step_1_wd, output_prefix)
    link_stats_combined_filtered_s1             = '%s/%s_stats_combined_filtered.txt'                      % (step_1_wd, output_prefix)
    linking_reads_rd1                           = '%s/%s_linking_reads_rd1.txt'                            % (step_1_wd, output_prefix)
    mafft_seq_folder                            = '%s/%s_linkage_visualization_rd1'                        % (step_1_wd, output_prefix)
    rd1_r1_to_extract                           = '%s/rd1_r1_to_extract.txt'                               % step_1_wd
    rd1_r2_to_extract                           = '%s/rd1_r2_to_extract.txt'                               % step_1_wd
    rd1_extracted_all_r1                        = '%s/rd1_extracted_all_r1.fasta'                          % step_1_wd
    rd1_extracted_all_r2                        = '%s/rd1_extracted_all_r2.fasta'                          % step_1_wd
    rd1_extracted_p_r1                          = '%s/rd1_extracted_R1.fasta'                              % step_1_wd
    rd1_extracted_p_r2                          = '%s/rd1_extracted_R2.fasta'                              % step_1_wd
    rd1_extracted_up                            = '%s/rd1_extracted_UP.fasta'                              % step_1_wd
    rd1_extracted_to_gnm_sam                    = '%s/rd1_extracted_to_gnm.sam'                            % step_1_wd
    rd1_extracted_to_gnm_sam_log                = '%s/rd1_extracted_to_gnm.sam.log'                        % step_1_wd
    rd1_extracted_to_gnm_sam_reformatted        = '%s/rd1_extracted_to_gnm_reformatted.sam'                % step_1_wd
    rd1_extracted_to_gnm_sam_reformat_log       = '%s/rd1_extracted_to_gnm_reformat.log'                   % step_1_wd
    rd1_extracted_to_gnm_sam_reformatted_sorted = '%s/rd1_extracted_to_gnm_reformatted_sorted.sam'         % step_1_wd
    linking_reads_tab                           = '%s/linking_reads.txt'                                   % step_1_wd
    linking_reads_r1_txt                        = '%s/rd1_linking_reads_R1.txt'                            % step_1_wd
    linking_reads_r2_txt                        = '%s/rd1_linking_reads_R2.txt'                            % step_1_wd
    linking_reads_r1_fasta                      = '%s/rd1_linking_reads_R1.fasta'                          % step_1_wd
    linking_reads_r2_fasta                      = '%s/rd1_linking_reads_R2.fasta'                          % step_1_wd
    linked_contigs_txt                          = '%s/linked_contigs_rd1.txt'                              % step_1_wd
    linked_contigs_fasta                        = '%s/linked_contigs_rd1.fasta'                            % step_1_wd
    rd1_clp_pct_diff_txt                        = '%s/rd1_clp_pct_diff.txt'                                % step_1_wd
    rd1_clp_pct_diff_txt_to_ignore              = '%s/rd1_clp_pct_diff_to_ignore.txt'                      % step_1_wd
    linked_by_clp_pct                           = '%s/rd1_linked_by_clp_pct.txt'                      % step_1_wd

    blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % num_threads
    bbmap_parameter  = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % (num_threads, bbmap_memory)
    # --very-sensitive-local    -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
    bowtie_parameter                = '--local --all --no-unal -N 1 -L 30'
    bowtie_parameter_mini_assembly  = ' --all --no-unal -N 1 -L 30'
    #bowtie_parameter = '--local --all --no-unal --very-sensitive-local'
    #bowtie_parameter = '--local --all --no-unal -D 30 -R 5 -N 0 -L 20 -i S,1,0.50'


    ####################################################################################################################
    ################################################### depth related ##################################################
    ####################################################################################################################

    ###################################### check input file for depth calculation ######################################

    #################################### calculate mean depth for genome/assemblies ####################################

    reads_file_r1_subset = '%s/input_R1_subset.fa' % mag_folder_in_wd
    reads_file_r2_subset = '%s/input_R2_subset.fa' % mag_folder_in_wd

    mean_depth_dict_gnm = {}
    if min_16s_gnm_multiple > 0:

        if depth_file_mag is None:

            report_and_log(('Round 1: depth info will be considered (but not provided!) for linking MAGs and 16S rRNA genes'), pwd_log_file, keep_quiet)
            report_and_log(('Round 1: depth estimation is time consuming, if you already have them, specify with -depth_16s and -depth_mag'), pwd_log_file, keep_quiet)
            report_and_log(('Round 1: calculating MAG depth now, be patient!'), pwd_log_file, keep_quiet)

            # get the number of paired reads
            input_r1_line_num_file = '%s/R1_line_num.txt' % (mag_folder_in_wd)
            input_r2_line_num_file = '%s/R2_line_num.txt' % (mag_folder_in_wd)
            os.system('wc -l %s > %s' % (reads_file_r1, input_r1_line_num_file))
            os.system('wc -l %s > %s' % (reads_file_r2, input_r2_line_num_file))
            paired_r1_num = int(int(open(input_r1_line_num_file).readline().strip().split(' ')[0]) / 2)
            paired_r2_num = int(int(open(input_r2_line_num_file).readline().strip().split(' ')[0]) / 2)
            if reads_file_r1[-1] in ['Q', 'q']:
                paired_r1_num = int(int(open(input_r1_line_num_file).readline().strip().split(' ')[0]) / 2)
                paired_r2_num = int(int(open(input_r2_line_num_file).readline().strip().split(' ')[0]) / 2)

            if paired_r1_num != paired_r2_num:
                print('Inconsistent number of reads found in r1 and r2, program exited!')
                exit()

            # get the number of reads paired to subset
            to_extract_reads_num = round(paired_r1_num * subsample_rate_for_depth_estimation)

            # remember to use the same random seed to keep pairing
            subsample_r1_cmd = 'seqtk sample -s100 %s %s > %s' % (reads_file_r1, to_extract_reads_num, reads_file_r1_subset)
            subsample_r2_cmd = 'seqtk sample -s100 %s %s > %s' % (reads_file_r2, to_extract_reads_num, reads_file_r2_subset)

            # subsample with multiprocessing
            pool = mp.Pool(processes=2)
            pool.map(os.system, [subsample_r1_cmd, subsample_r2_cmd])
            pool.close()
            pool.join()

            # get mean depth for contig
            mean_depth_dict_ctg, ctg_len_dict = get_ctg_mean_depth_by_samtools_coverage_global(True, combined_input_gnms, reads_file_r1_subset, reads_file_r2_subset, '', subsample_rate_for_depth_estimation, num_threads)
            #os.system('rm %s' % reads_file_r1_subset)
            #os.system('rm %s' % reads_file_r2_subset)

            # write out ctg depth
            depth_file_ctg_handle = open(depth_file_ctg, 'w')
            for ctg in mean_depth_dict_ctg:
                depth_file_ctg_handle.write('%s\t%s\n' % (ctg, mean_depth_dict_ctg[ctg]))
            depth_file_ctg_handle.close()

            # get mean_depth_dict_gnm
            gnm_len_total_depth_dict = {}
            for ctg in mean_depth_dict_ctg:
                ctg_genome = ctg.split(gnm_to_ctg_connector)[0]
                ctg_len = ctg_len_dict[ctg]
                ctg_depth = mean_depth_dict_ctg[ctg]
                ctg_total_depth = ctg_depth * ctg_len
                if ctg_genome not in gnm_len_total_depth_dict:
                    gnm_len_total_depth_dict[ctg_genome] = [ctg_len, ctg_total_depth]
                else:
                    gnm_len_total_depth_dict[ctg_genome][0] += ctg_len
                    gnm_len_total_depth_dict[ctg_genome][1] += ctg_total_depth

            for each_gnm in gnm_len_total_depth_dict:
                gnm_len = gnm_len_total_depth_dict[each_gnm][0]
                gnm_total_depth = gnm_len_total_depth_dict[each_gnm][1]
                gnm_mean_depth = float("{0:.6f}".format(gnm_total_depth / gnm_len))
                mean_depth_dict_gnm[each_gnm] = gnm_mean_depth

            # write out gnm depth
            depth_file_gnm_handle = open(depth_file_gnm, 'w')
            for gnm in mean_depth_dict_gnm:
                depth_file_gnm_handle.write('%s\t%s\n' % (gnm, mean_depth_dict_gnm[gnm]))
            depth_file_gnm_handle.close()
        else:
            report_and_log(('Round 1: read in provided MAG depth'), pwd_log_file, keep_quiet)

            # read in depth and store in mean_depth_dict_gnm here
            for each_mag_depth in open(depth_file_mag):
                each_mag_depth_split = each_mag_depth.strip().split('\t')
                mag_id = each_mag_depth_split[0]
                mag_depth = float(each_mag_depth_split[1])
                mean_depth_dict_gnm[mag_id] = mag_depth


    ###################################### calculate mean depth for 16S sequences ######################################

    # copy input 16S into wd
    os.mkdir(input_16s_folder_in_wd)
    os.system('cp %s %s/' % (marker_gene_seqs, input_16s_folder_in_wd))

    # polish input 16S by default
    marker_path, marker_base, marker_ext = sep_path_basename_ext(marker_gene_seqs)
    input_16s_in_wd              = '%s/%s%s'             % (input_16s_folder_in_wd, marker_base, marker_ext)
    input_16s_in_wd_no_ext       = '%s/%s'               % (input_16s_folder_in_wd, marker_base)
    input_16s_polished           = '%s/%s.polished%s'    % (input_16s_folder_in_wd, marker_base, marker_ext)
    input_16s_polished_gff       = '%s/%s.polished.gff'  % (input_16s_folder_in_wd, marker_base)
    input_16s_polished_no_ext    = '%s/%s.polished'      % (input_16s_folder_in_wd, marker_base)

    if no_polish is False:
        report_and_log(('Round 1: polishing input 16S rRNA genes'), pwd_log_file, keep_quiet)
        polish_16s(input_16s_in_wd, input_16s_polished)
        os.system('cp %s %s/' % (input_16s_polished, working_directory))
        os.system('cp %s %s/' % (input_16s_polished_gff, working_directory))
    else:
        input_16s_polished        = input_16s_in_wd
        input_16s_polished_no_ext = input_16s_in_wd_no_ext

    # calculate depth for input 16S
    mean_depth_dict_16s = {}
    if min_16s_gnm_multiple > 0:

        if depth_file_16s is None:

            report_and_log(('Round 1: calculating depth for 16S'), pwd_log_file, keep_quiet)

            # get mean depth for 16S sequences
            mean_depth_dict_16s, s16_len_dict = get_ctg_mean_depth_by_samtools_coverage_global(True, input_16s_polished, reads_file_r1_subset, reads_file_r2_subset, '', subsample_rate_for_depth_estimation, num_threads)

            # write out 16s depth
            depth_file_16s_handle = open(depth_file_16s_calculated, 'w')
            for s16 in mean_depth_dict_16s:
                depth_file_16s_handle.write('%s\t%s\n' % (s16, mean_depth_dict_16s[s16]))
            depth_file_16s_handle.close()
        else:
            report_and_log(('Round 1: read in provided 16S depth'), pwd_log_file, keep_quiet)

            # read in depth and store in mean_depth_dict_gnm here
            for each_16s_depth in open(depth_file_16s):
                each_16s_depth_split = each_16s_depth.strip().split('\t')
                s16_id = each_16s_depth_split[0]
                s16_depth = float(each_16s_depth_split[1])
                mean_depth_dict_16s[s16_id] = s16_depth

    # made 16s depth and genome depth comparable


    ####################################################################################################################
    ############################################### first round linking ################################################
    ####################################################################################################################

    ######################################## map reads to marker gene sequences ########################################

    #bbmap_index_and_mapping_cmd = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, input_16s_in_wd, reads_file_r1, reads_file_r2, input_reads_to_16s_sam, bbmap_parameter, input_reads_to_16s_sam_bbmap_stderr)
    #report_and_log((bbmap_index_and_mapping_cmd), pwd_log_file, True)
    #os.system(bbmap_index_and_mapping_cmd)

    bowtie_build_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, input_16s_polished, input_16s_polished_no_ext)
    os.system(bowtie_build_cmd)
    report_and_log((bowtie_build_cmd), pwd_log_file, True)

    if reads_vs_16s_sam is None:
        report_and_log(('Round 1: Mapping input reads to marker genes'), pwd_log_file, keep_quiet)

        bowtie_read_to_16s_cmd = 'bowtie2 -x %s -U %s,%s -S %s -p %s -f %s 2> %s' % (input_16s_polished_no_ext, reads_file_r1_fasta, reads_file_r2_fasta, input_reads_to_16s_sam_bowtie, num_threads, bowtie_parameter, input_reads_to_16s_sam_bowtie_log)
        report_and_log((bowtie_read_to_16s_cmd), pwd_log_file, True)
        os.system(bowtie_read_to_16s_cmd)

        report_and_log(('Round 1: transforming cigar format from 1.3 to 1.4'), pwd_log_file, keep_quiet)
        bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (input_reads_to_16s_sam_bowtie, input_reads_to_16s_sam, input_reads_to_16s_sam_reformat_log)
        report_and_log((bbmap_reformat_cmd), pwd_log_file, True)
        os.system(bbmap_reformat_cmd)
        os.remove(input_reads_to_16s_sam_bowtie)

        # sort sam file first
        report_and_log(('Round 1: sorting mappping results'), pwd_log_file, keep_quiet)
        sort_by_read_cmd = 'samtools sort -n -O sam --threads %s -o %s %s' % (num_threads, input_reads_to_16s_sam_sorted, input_reads_to_16s_sam)
        os.system(sort_by_read_cmd)
        os.remove(input_reads_to_16s_sam)
    else:
        report_and_log(('Round 1: Sorted reformatted sam file detected'), pwd_log_file, keep_quiet)
        input_reads_to_16s_sam_sorted = reads_vs_16s_sam


    ##################################################### read in sam file ####################################################

    # get marker len dict
    marker_len_dict = {}
    for each_marker_record in SeqIO.parse(input_16s_polished, 'fasta'):
        marker_len_dict[each_marker_record.id] = len(each_marker_record.seq)

    # get the number of lines per file
    report_and_log(('Round 1: calculating the number of lines per subset'), pwd_log_file, keep_quiet)
    os.system('wc -l %s > %s' % (input_reads_to_16s_sam_sorted, input_reads_to_16s_sam_sorted_line_num))
    sam16s_line_num = int(open(input_reads_to_16s_sam_sorted_line_num).readline().strip().split(' ')[0])
    os.remove(input_reads_to_16s_sam_sorted_line_num)
    line_num_per_file = int(round(sam16s_line_num/(num_threads*10))) + 10

    report_and_log(('Round 1: splitting sam file'), pwd_log_file, keep_quiet)
    os.mkdir(input_reads_to_16s_sam_sorted_split_folder)
    split_sam_cmd = 'split -l %s %s %s/splitted_sam_ ' % (line_num_per_file, input_reads_to_16s_sam_sorted, input_reads_to_16s_sam_sorted_split_folder)
    os.system(split_sam_cmd)

    report_and_log(('Round 1: analysing mappping results with %s threads' % num_threads), pwd_log_file, keep_quiet)
    os.mkdir(input_reads_to_16s_sam_MappingRecord_folder)

    # get splitted sam file list
    splitted_sam_file_list = [os.path.basename(file_name) for file_name in glob.glob('%s/splitted_sam_*' % input_reads_to_16s_sam_sorted_split_folder)]

    # prepare lol for mp worker
    list_for_parse_sam16s_worker = []
    splitted_sam_mp_file_set = set()
    for splitted_sam_file in splitted_sam_file_list:
        pwd_splitted_sam_file    = '%s/%s'        % (input_reads_to_16s_sam_sorted_split_folder, splitted_sam_file)
        pwd_splitted_sam_mp_file = '%s/%s_mp.txt' % (input_reads_to_16s_sam_MappingRecord_folder, splitted_sam_file)
        splitted_sam_mp_file_set.add(pwd_splitted_sam_mp_file)
        list_for_parse_sam16s_worker.append([pwd_splitted_sam_file,
                                             pwd_splitted_sam_mp_file,
                                             min_M_len_16s,
                                             mismatch_cutoff,
                                             marker_len_dict])

    pool_parse_sam16s = mp.Pool(processes=num_threads)
    pool_parse_sam16s.map(parse_sam16s_worker, list_for_parse_sam16s_worker)
    pool_parse_sam16s.close()
    pool_parse_sam16s.join()

    report_and_log(('Round 1: removing splitted subsets from disk'), pwd_log_file, keep_quiet)
    os.system('rm -r %s' % input_reads_to_16s_sam_sorted_split_folder)

    # reads filter results into dict
    report_and_log(('Round 1: reading filtered alignments into dict'), pwd_log_file, keep_quiet)
    MappingRecord_dict = {}
    read_base_to_pop = set()
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

    # print('MappingRecord_dict length: %s' % len(MappingRecord_dict))
    # for each_mp in MappingRecord_dict.copy():
    #     if each_mp in read_base_to_pop:
    #         MappingRecord_dict.pop(each_mp)
    # print('MappingRecord_dict length: %s' % len(MappingRecord_dict))


    ############################# map extracted reads to combined input genomes #############################

    report_and_log(('Round 1: extracting sequences of reads matched to 16S'), pwd_log_file, keep_quiet)

    # write out id of reads to extract
    to_extract_read_base_list = sorted(list(MappingRecord_dict.keys()))
    with open(rd1_r1_to_extract, 'w') as rd1_r1_to_extract_handle:
        rd1_r1_to_extract_handle.write('%s\n' % '\n'.join([('%s.1' % i) for i in to_extract_read_base_list]))
    with open(rd1_r2_to_extract, 'w') as rd1_r2_to_extract_handle:
        rd1_r2_to_extract_handle.write('%s\n' % '\n'.join([('%s.2' % i) for i in to_extract_read_base_list]))

    # extract reads with seqtk
    seqtk_extract_cmd_rd1_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd1_r1_to_extract, rd1_extracted_all_r1)
    seqtk_extract_cmd_rd1_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd1_r2_to_extract, rd1_extracted_all_r2)
    report_and_log((seqtk_extract_cmd_rd1_r1), pwd_log_file, True)
    report_and_log((seqtk_extract_cmd_rd1_r2), pwd_log_file, True)
    os.system(seqtk_extract_cmd_rd1_r1)
    os.system(seqtk_extract_cmd_rd1_r2)

    # read extracted read sequences into dict
    extract_read_seq_dict = {}
    for extracted_r1 in SeqIO.parse(rd1_extracted_all_r1, 'fasta'):
        extract_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
    for extracted_r2 in SeqIO.parse(rd1_extracted_all_r2, 'fasta'):
        extract_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)

    # write out paired in the same order
    rd1_extracted_p_r1_handle = open(rd1_extracted_p_r1, 'w')
    rd1_extracted_p_r2_handle = open(rd1_extracted_p_r2, 'w')
    for each_read_base in to_extract_read_base_list:
        current_r1 = '%s.1' % each_read_base
        current_r2 = '%s.2' % each_read_base
        current_r1_seq = extract_read_seq_dict.get(current_r1, '')
        current_r2_seq = extract_read_seq_dict.get(current_r2, '')
        rd1_extracted_p_r1_handle.write('>%s\n' % current_r1)
        rd1_extracted_p_r1_handle.write('%s\n' % current_r1_seq)
        rd1_extracted_p_r2_handle.write('>%s\n' % current_r2)
        rd1_extracted_p_r2_handle.write('%s\n' % current_r2_seq)
    rd1_extracted_p_r1_handle.close()
    rd1_extracted_p_r2_handle.close()

    os.remove(rd1_extracted_all_r1)
    os.remove(rd1_extracted_all_r2)

    ############################# map extracted reads to combined input genomes #############################

    report_and_log(('Round 1: mapping extracted reads to input genomes'), pwd_log_file, keep_quiet)

    combined_input_gnms_no_ext = '.'.join(combined_input_gnms.split('.')[:-1])
    bowtie_build_input_gnm_cmd      = 'bowtie2-build --quiet --threads %s -f %s %s'       % (num_threads, combined_input_gnms, combined_input_gnms_no_ext)
    bowtie_cmd_rd1_extracted_to_mag = 'bowtie2 -x %s -U %s,%s -S %s -p %s -f %s 2> %s' % (combined_input_gnms_no_ext, rd1_extracted_p_r1, rd1_extracted_p_r2, rd1_extracted_to_gnm_sam, num_threads, bowtie_parameter, rd1_extracted_to_gnm_sam_log)

    # index ref
    os.system(bowtie_build_input_gnm_cmd)
    report_and_log((bowtie_cmd_rd1_extracted_to_mag), pwd_log_file, True)
    os.system(bowtie_cmd_rd1_extracted_to_mag)

    bbmap_reformat_cmd_rd1_extracted_to_mag = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (rd1_extracted_to_gnm_sam, rd1_extracted_to_gnm_sam_reformatted, rd1_extracted_to_gnm_sam_reformat_log)
    report_and_log((bbmap_reformat_cmd_rd1_extracted_to_mag), pwd_log_file, True)
    os.system(bbmap_reformat_cmd_rd1_extracted_to_mag)

    # sort sam file first
    sort_by_read_cmd = 'samtools sort -n -O sam --threads %s -o %s %s' % (num_threads, rd1_extracted_to_gnm_sam_reformatted_sorted, rd1_extracted_to_gnm_sam_reformatted)
    os.system(sort_by_read_cmd)


    ######################################### read mapping results of rd1 extracted mates into mp dict  #########################################

    report_and_log(('Round 1: analysing mappping results'), pwd_log_file, keep_quiet)

    ctg_len_dict = {}
    for each_ctg_record in SeqIO.parse(combined_input_gnms, 'fasta'):
        ctg_len_dict[each_ctg_record.id] = len(each_ctg_record.seq)

    ctg_ignore_region_dict = get_regions_to_ignore_from_barrnap_output(combined_barrnap_gff, ctg_len_dict)

    processed_num = 0
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
                    for r1_ctg_ref in current_read_base_r1_ctg_ref_dict:
                        r1_ctg_ref_matched_pos_dict = current_read_base_r1_ctg_ref_dict[r1_ctg_ref]

                        # one read need to mapped to one ctg only for one time
                        if len(r1_ctg_ref_matched_pos_dict) > 1:
                            refs_to_ignore_ctg.add(r1_ctg_ref)
                        else:
                            r1_ctg_ref_pos = list(r1_ctg_ref_matched_pos_dict.keys())[0]
                            r1_ctg_ref_cigar = r1_ctg_ref_matched_pos_dict[r1_ctg_ref_pos]
                            r1_ctg_ref_len = ctg_len_dict[r1_ctg_ref]
                            qualified_cigar = check_cigar_quality(r1_ctg_ref_cigar, mismatch_cutoff, min_M_len_ctg,
                                                                  r1_ctg_ref_pos, r1_ctg_ref_len)
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
                                            if (ctg_len_dict[r1_ctg_ref] - r1_ctg_ref_pos) <= 50:
                                                matched_to_r1_ref_ignored_region = True

                                if matched_to_r1_ref_ignored_region is False:
                                    r1_ctg_refs_passed_qc[r1_ctg_ref] = [r1_ctg_ref_cigar]
                            else:
                                refs_to_ignore_ctg.add(r1_ctg_ref)

                    ########## filter r2 ctg refs ##########

                    r2_ctg_refs_passed_qc = {}
                    for r2_ctg_ref in current_read_base_r2_ctg_ref_dict:
                        r2_ctg_ref_matched_pos_dict = current_read_base_r2_ctg_ref_dict[r2_ctg_ref]

                        # one read need to mapped to one ctg only for one time
                        if len(r2_ctg_ref_matched_pos_dict) > 1:
                            refs_to_ignore_ctg.add(r2_ctg_ref)
                        else:
                            r2_ctg_ref_pos = list(r2_ctg_ref_matched_pos_dict.keys())[0]
                            r2_ctg_ref_cigar = r2_ctg_ref_matched_pos_dict[r2_ctg_ref_pos]
                            r2_ctg_ref_len = ctg_len_dict.get(r2_ctg_ref, 0)
                            qualified_cigar = check_cigar_quality(r2_ctg_ref_cigar, mismatch_cutoff, min_M_len_ctg,
                                                                  r2_ctg_ref_pos, r2_ctg_ref_len)
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
                                            if (ctg_len_dict[r2_ctg_ref] - r2_ctg_ref_pos) <= 50:
                                                matched_to_r2_ref_ignored_region = True

                                if matched_to_r2_ref_ignored_region is False:
                                    r2_ctg_refs_passed_qc[r2_ctg_ref] = [r2_ctg_ref_cigar]
                            else:
                                refs_to_ignore_ctg.add(r2_ctg_ref)

                    ####################################################################################################

                    r1_ctg_refs_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}
                    r2_ctg_refs_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}

                    # no mate matched to ctg
                    if (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) == 0):
                        pass

                    # only r1 matched to ctg
                    elif (len(r1_ctg_refs_no_ignored) > 0) and (len(r2_ctg_refs_no_ignored) == 0):
                        MappingRecord_dict[current_read_base].matched_to_ctg = True
                        MappingRecord_dict[current_read_base].r1_ctg_refs_no_ignored = r1_ctg_refs_no_ignored
                        MappingRecord_dict[current_read_base].r1_ctg_ref_dict = current_read_base_r1_ctg_ref_dict
                    # only r2 matched to ctg
                    elif (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) > 0):
                        MappingRecord_dict[current_read_base].matched_to_ctg = True
                        MappingRecord_dict[current_read_base].r2_ctg_refs_no_ignored = r2_ctg_refs_no_ignored
                        MappingRecord_dict[current_read_base].r2_ctg_ref_dict = current_read_base_r2_ctg_ref_dict

                    # both r1 and r2 matched to ctg
                    else:
                        shared_ctg_ref_set = {key: [r1_ctg_refs_no_ignored[key][0], r2_ctg_refs_no_ignored[key][0]] for key in set(r1_ctg_refs_no_ignored).intersection(set(r2_ctg_refs_no_ignored))}
                        if len(shared_ctg_ref_set) > 0:
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

                    ########################################### report ###########################################

                    processed_num += 1
                    if (processed_num % 50000 == 0):
                        report_and_log(('Round 1: processed %sk pairs of reads' % int(processed_num / 1000)), pwd_log_file, keep_quiet)

        report_and_log(('Round 1: processed %sk' % float("{0:.2f}".format(processed_num / 1000))), pwd_log_file, keep_quiet)

    # remove sam files from disk
    os.remove(rd1_extracted_to_gnm_sam)
    os.remove(rd1_extracted_to_gnm_sam_reformatted)


    ##################################################### get linkages from MappingRecord_dict #####################################################

    # def add_linkage_from_dict(qualified_read, refs_dict_16s, refs_dict_ctg,
    #                           link_num_dict, linking_reads_dict,
    #                           cigar_dict_16s_side, cigar_dict_ctg_side,
    #                           linked_by, linked_by_clip_num_dict, linked_by_pair_num_dict,
    #                           key_connector):
    #
    #     for each_16s_ref in refs_dict_16s:
    #         current_cigar_list_16s_side = refs_dict_16s[each_16s_ref]
    #         for each_ctg_ref in refs_dict_ctg:
    #             current_cigar_list_ctg_side = refs_dict_ctg[each_ctg_ref]
    #
    #             marker_to_ctg_key = '%s%s%s' % (each_16s_ref, key_connector, each_ctg_ref)
    #
    #             if linked_by == 'pair':
    #                 if marker_to_ctg_key not in linked_by_pair_num_dict:
    #                     linked_by_pair_num_dict[marker_to_ctg_key] = 1
    #                 else:
    #                     linked_by_pair_num_dict[marker_to_ctg_key] += 1
    #
    #             if linked_by == 'clip':
    #                 if marker_to_ctg_key not in linked_by_clip_num_dict:
    #                     linked_by_clip_num_dict[marker_to_ctg_key] = 1
    #                 else:
    #                     linked_by_clip_num_dict[marker_to_ctg_key] += 1
    #
    #             # initialize key values if not in the dict
    #             if marker_to_ctg_key not in link_num_dict:
    #                 link_num_dict[marker_to_ctg_key] = 0
    #                 linking_reads_dict[marker_to_ctg_key] = set()
    #                 cigar_dict_16s_side[marker_to_ctg_key] = []
    #                 cigar_dict_ctg_side[marker_to_ctg_key] = []
    #
    #             # add values to dict
    #             link_num_dict[marker_to_ctg_key] += 1
    #             linking_reads_dict[marker_to_ctg_key].add(qualified_read)
    #             for each in current_cigar_list_16s_side:
    #                 cigar_dict_16s_side[marker_to_ctg_key].append(each)
    #             for each in current_cigar_list_ctg_side:
    #                 cigar_dict_ctg_side[marker_to_ctg_key].append(each)

    def add_linkage_from_dict(qualified_read, refs_dict_16s, refs_dict_ctg,
                              link_num_dict, linking_reads_dict,
                              cigar_dict_16s_side, cigar_dict_ctg_side,
                              linked_by, linked_by_clip_num_dict, linked_by_pair_num_dict,
                              key_connector):

        for each_16s_ref in refs_dict_16s:
            current_cigar_list_16s_side = refs_dict_16s[each_16s_ref]
            for each_ctg_ref in refs_dict_ctg:
                current_cigar_list_ctg_side = refs_dict_ctg[each_ctg_ref]
                marker_to_ctg_key = '%s%s%s' % (each_16s_ref, key_connector, each_ctg_ref)
                read_base_to_count = False

                if linked_by == 'both':
                    read_base_to_count = True

                if linked_by == 'pair':
                    read_base_to_count = True
                    if marker_to_ctg_key not in linked_by_pair_num_dict:
                        linked_by_pair_num_dict[marker_to_ctg_key] = 1
                    else:
                        linked_by_pair_num_dict[marker_to_ctg_key] += 1

                if linked_by == 'clip':
                    all_cigar_clipping = True
                    for each in current_cigar_list_16s_side:
                        if ('S' not in each) and ('s' not in each):
                            all_cigar_clipping = False
                    for each in current_cigar_list_ctg_side:
                        if ('S' not in each) and ('s' not in each):
                            all_cigar_clipping = False

                    if all_cigar_clipping is True:
                        read_base_to_count = True

                    if read_base_to_count is True:
                        if marker_to_ctg_key not in linked_by_clip_num_dict:
                            linked_by_clip_num_dict[marker_to_ctg_key] = 1
                        else:
                            linked_by_clip_num_dict[marker_to_ctg_key] += 1

                if read_base_to_count is True:
                    # initialize key values if not in the dict
                    if marker_to_ctg_key not in link_num_dict:
                        link_num_dict[marker_to_ctg_key] = 0
                        linking_reads_dict[marker_to_ctg_key] = set()
                        cigar_dict_16s_side[marker_to_ctg_key] = []
                        cigar_dict_ctg_side[marker_to_ctg_key] = []

                    # add values to dict
                    link_num_dict[marker_to_ctg_key] += 1
                    linking_reads_dict[marker_to_ctg_key].add(qualified_read)
                    for each in current_cigar_list_16s_side:
                        cigar_dict_16s_side[marker_to_ctg_key].append(each)
                    for each in current_cigar_list_ctg_side:
                        cigar_dict_ctg_side[marker_to_ctg_key].append(each)

    report_and_log(('Round 1: Parsing MappingRecord dictionary to get linkages'), pwd_log_file, keep_quiet)

    marker_to_ctg_linkage_cigar_dict_16s_side = {}
    marker_to_ctg_linkage_cigar_dict_ctg_side = {}
    marker_to_ctg_linkage_num_dict = {}
    marker_to_ctg_linking_reads_dict = {}
    linked_by_clip_num_dict = {}
    linked_by_pair_num_dict = {}
    for qualified_read in MappingRecord_dict:

        r1_16s_refs     = MappingRecord_dict[qualified_read].r1_16s_refs_no_ignored
        r2_16s_refs     = MappingRecord_dict[qualified_read].r2_16s_refs_no_ignored
        shared_16s_refs = MappingRecord_dict[qualified_read].shared_16s_refs_no_ignored
        r1_ctg_refs     = MappingRecord_dict[qualified_read].r1_ctg_refs_no_ignored
        r2_ctg_refs     = MappingRecord_dict[qualified_read].r2_ctg_refs_no_ignored
        shared_ctg_refs = MappingRecord_dict[qualified_read].shared_ctg_refs_no_ignored

        if (len(shared_16s_refs) > 0) or (len(shared_ctg_refs) > 0):
            if (len(shared_16s_refs) > 0) and (len(shared_ctg_refs) > 0):
                add_linkage_from_dict(qualified_read, shared_16s_refs, shared_ctg_refs,
                                      marker_to_ctg_linkage_num_dict, marker_to_ctg_linking_reads_dict,
                                      marker_to_ctg_linkage_cigar_dict_16s_side,
                                      marker_to_ctg_linkage_cigar_dict_ctg_side,
                                      'both', linked_by_clip_num_dict, linked_by_pair_num_dict,
                                      marker_to_ctg_gnm_Key_connector)
        else:
            if (len(r1_16s_refs) > 0) and (len(r2_ctg_refs) > 0):
                add_linkage_from_dict(qualified_read, r1_16s_refs, r2_ctg_refs,
                                      marker_to_ctg_linkage_num_dict, marker_to_ctg_linking_reads_dict,
                                      marker_to_ctg_linkage_cigar_dict_16s_side,
                                      marker_to_ctg_linkage_cigar_dict_ctg_side,
                                      'pair', linked_by_clip_num_dict, linked_by_pair_num_dict,
                                      marker_to_ctg_gnm_Key_connector)

            if (len(r2_16s_refs) > 0) and (len(r1_ctg_refs) > 0):
                add_linkage_from_dict(qualified_read, r2_16s_refs, r1_ctg_refs,
                                      marker_to_ctg_linkage_num_dict, marker_to_ctg_linking_reads_dict,
                                      marker_to_ctg_linkage_cigar_dict_16s_side,
                                      marker_to_ctg_linkage_cigar_dict_ctg_side,
                                      'pair', linked_by_clip_num_dict, linked_by_pair_num_dict,
                                      marker_to_ctg_gnm_Key_connector)

            if (len(r1_16s_refs) > 0) and (len(r1_ctg_refs) > 0):
                add_linkage_from_dict(qualified_read, r1_16s_refs, r1_ctg_refs,
                                      marker_to_ctg_linkage_num_dict, marker_to_ctg_linking_reads_dict,
                                      marker_to_ctg_linkage_cigar_dict_16s_side,
                                      marker_to_ctg_linkage_cigar_dict_ctg_side,
                                      'clip', linked_by_clip_num_dict, linked_by_pair_num_dict,
                                      marker_to_ctg_gnm_Key_connector)

            if (len(r2_16s_refs) > 0) and (len(r2_ctg_refs) > 0):
                add_linkage_from_dict(qualified_read, r2_16s_refs, r2_ctg_refs,
                                      marker_to_ctg_linkage_num_dict, marker_to_ctg_linking_reads_dict,
                                      marker_to_ctg_linkage_cigar_dict_16s_side,
                                      marker_to_ctg_linkage_cigar_dict_ctg_side,
                                      'clip', linked_by_clip_num_dict, linked_by_pair_num_dict,
                                      marker_to_ctg_gnm_Key_connector)

    # filter with link num
    marker_to_ctg_linkage_num_dict_min3 = {}
    for each_link in marker_to_ctg_linkage_num_dict:
        if marker_to_ctg_linkage_num_dict[each_link] >= 3:
            marker_to_ctg_linkage_num_dict_min3[each_link] = marker_to_ctg_linkage_num_dict[each_link]

    # calculate clipping reads ratio
    rd1_clp_pct_diff_txt_handle = open(rd1_clp_pct_diff_txt, 'w')
    rd1_clp_pct_diff_txt_to_ignore_handle = open(rd1_clp_pct_diff_txt_to_ignore, 'w')
    rd1_clp_pct_diff_txt_handle.write('Marker\tGenome\tContig\tclp_pct_ctg\tclp_pct_16s\tRatio\n')
    rd1_clp_pct_diff_txt_to_ignore_handle.write('Marker\tGenome\tContig\tclp_pct_ctg\tclp_pct_16s\tRatio\n')
    linkages_to_ignore = set()
    for each_link in marker_to_ctg_linkage_num_dict_min3:
        each_link_split = each_link.split(marker_to_ctg_gnm_Key_connector)
        marker_id = each_link_split[0]
        gnm_id = each_link_split[1].split(gnm_to_ctg_connector)[0]
        ctg_id = each_link_split[1].split(gnm_to_ctg_connector)[1]

        linkage_cigar_16s_side_all = marker_to_ctg_linkage_cigar_dict_16s_side[each_link]
        linkage_cigar_16s_side_clp = [i for i in linkage_cigar_16s_side_all if (('S' in i) or ('s' in i))]
        linkage_cigar_16s_side_clp_pct = float("{0:.2f}".format(len(linkage_cigar_16s_side_clp) * 100 / len(linkage_cigar_16s_side_all)))

        linkage_cigar_ctg_side_all = marker_to_ctg_linkage_cigar_dict_ctg_side[each_link]
        linkage_cigar_ctg_side_clp = [i for i in linkage_cigar_ctg_side_all if (('S' in i) or ('s' in i))]
        linkage_cigar_ctg_side_clp_pct = float("{0:.2f}".format(len(linkage_cigar_ctg_side_clp) * 100 / len(linkage_cigar_ctg_side_all)))

        to_ignore = False
        if (len(linkage_cigar_16s_side_all) >= 50) or (len(linkage_cigar_ctg_side_all) >= 50):
            if (linkage_cigar_16s_side_clp_pct == 0) and (linkage_cigar_ctg_side_clp_pct >= 25):
                to_ignore = True
            elif (linkage_cigar_16s_side_clp_pct >= 25) and (linkage_cigar_ctg_side_clp_pct == 0):
                to_ignore = True
            elif (linkage_cigar_16s_side_clp_pct > 0) and (linkage_cigar_ctg_side_clp_pct > 0):
                if ((linkage_cigar_16s_side_clp_pct/linkage_cigar_ctg_side_clp_pct) >= 3) or ((linkage_cigar_ctg_side_clp_pct/linkage_cigar_16s_side_clp_pct) >= 3):
                    to_ignore = True

        # mean_clp_num = (linkage_cigar_16s_side_clp_num + linkage_cigar_ctg_side_clp_num)/2
        # if (linkage_cigar_16s_side_clp_num == 0) and (linkage_cigar_ctg_side_clp_num >= 10):
        #     to_ignore = True
        # elif (linkage_cigar_16s_side_clp_num >= 10) and (linkage_cigar_ctg_side_clp_num == 0):
        #     to_ignore = True
        # elif (linkage_cigar_16s_side_clp_num > 0) and (linkage_cigar_ctg_side_clp_num > 0):
        #     if mean_clp_num <= 20:
        #         if not (0.33 <= (linkage_cigar_16s_side_clp_num/linkage_cigar_ctg_side_clp_num) <= 3):
        #             to_ignore = True
        #     # elif 20 < mean_clp_num <= 50:
        #     #     if not (0.5<= (linkage_cigar_16s_side_clp_num/linkage_cigar_ctg_side_clp_num) <= 2):
        #     #         to_ignore = True
        #     elif mean_clp_num > 20:
        #         if not (0.5 <= (linkage_cigar_16s_side_clp_num/linkage_cigar_ctg_side_clp_num) <= 2):
        #             to_ignore = True

        # clp_pct_16s_side = 'na'
        # if len(linkage_cigar_16s_side_all) > 0:
        #     clp_pct_16s_side = len(linkage_cigar_16s_side_clp) * 100 / len(linkage_cigar_16s_side_all)
        #     clp_pct_16s_side = float("{0:.2f}".format(clp_pct_16s_side))
        #
        # clp_pct_ctg_side = 'na'
        # if len(linkage_cigar_ctg_side_all) > 0:
        #     clp_pct_ctg_side = len(linkage_cigar_ctg_side_clp) * 100 / len(linkage_cigar_ctg_side_all)
        #     clp_pct_ctg_side = float("{0:.2f}".format(clp_pct_ctg_side))
        #
        # clp_pct_ratio = 'na'
        # if (clp_pct_16s_side != 'na') and (clp_pct_ctg_side != 'na'):
        #     if clp_pct_16s_side > 0:
        #         clp_pct_ratio = float("{0:.2f}".format(clp_pct_ctg_side / clp_pct_16s_side))
        #
        # to_ignore = False
        # if clp_pct_ratio != 'na':
        #     if (clp_pct_ctg_side >= clp_pct_ctg_side_max_num) and (clp_pct_ratio >= clp_pct_ratio_cutoff):
        #         to_ignore = True
        # else:
        #     if (clp_pct_ctg_side != 'na') and (clp_pct_16s_side != 'na'):
        #         if clp_pct_ctg_side >= clp_pct_ctg_side_max_num:
        #             to_ignore = True

        if to_ignore is True:
            linkages_to_ignore.add(each_link)
            rd1_clp_pct_diff_txt_to_ignore_handle.write('%s\t%s\t%s\t(%s/%s)\t(%s/%s)\n' % (marker_id, gnm_id, ctg_id, len(linkage_cigar_ctg_side_clp), len(linkage_cigar_ctg_side_all), len(linkage_cigar_16s_side_clp),len(linkage_cigar_16s_side_all)))
        else:
            rd1_clp_pct_diff_txt_handle.write('%s\t%s\t%s\t(%s/%s)\t(%s/%s)\n' % (marker_id, gnm_id, ctg_id, len(linkage_cigar_ctg_side_clp), len(linkage_cigar_ctg_side_all), len(linkage_cigar_16s_side_clp), len(linkage_cigar_16s_side_all)))

    rd1_clp_pct_diff_txt_handle.close()
    rd1_clp_pct_diff_txt_to_ignore_handle.close()

    # ignore linkages if did not pass the clp pct ratio test
    linked_by_clp_pct_handle = open(linked_by_clp_pct, 'w')
    marker_to_ctg_linkage_num_dict_min3_passed_ratio_check = {}
    linkages_to_ignore = set()
    for each_link in marker_to_ctg_linkage_num_dict_min3:
        linked_by_clip_num = linked_by_clip_num_dict.get(each_link, 0)
        linked_by_clp_pct = linked_by_clip_num * 100 / marker_to_ctg_linkage_num_dict_min3[each_link]
        linked_by_clp_pct = float("{0:.2f}".format(linked_by_clp_pct))
        linked_by_clp_pct_handle.write('%s\t%s(%s/%s)\n' % (each_link, linked_by_clp_pct, linked_by_clip_num, marker_to_ctg_linkage_num_dict_min3[each_link]))
        if each_link not in linkages_to_ignore:
            if linked_by_clp_pct <= 70:
                #print('%s\t%s\t%s\t%s' % (each_link, marker_to_ctg_linkage_num_dict[each_link], linked_by_clip_num, linked_by_clp_pct))
                marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_link] = marker_to_ctg_linkage_num_dict_min3[each_link]
    linked_by_clp_pct_handle.close()

    # get number of linkages at genome level
    marker_to_gnm_linkage_cigar_dict_16s_side = {}
    marker_to_gnm_linkage_cigar_dict_ctg_side = {}
    marker_to_gnm_link_num = {}
    for each_marker_to_ctg_key in marker_to_ctg_linkage_num_dict_min3_passed_ratio_check:
        marker_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[0]
        ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[1]
        gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
        marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
        marker_to_ctg_linkage_cigar_16s_side = marker_to_ctg_linkage_cigar_dict_16s_side[each_marker_to_ctg_key]
        marker_to_ctg_linkage_cigar_ctg_side = marker_to_ctg_linkage_cigar_dict_ctg_side[each_marker_to_ctg_key]

        if marker_to_gnm_key not in marker_to_gnm_link_num:
            marker_to_gnm_link_num[marker_to_gnm_key] = marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_marker_to_ctg_key]
            marker_to_gnm_linkage_cigar_dict_16s_side[marker_to_gnm_key] = marker_to_ctg_linkage_cigar_16s_side
            marker_to_gnm_linkage_cigar_dict_ctg_side[marker_to_gnm_key] = marker_to_ctg_linkage_cigar_ctg_side
        else:
            marker_to_gnm_link_num[marker_to_gnm_key] += marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_marker_to_ctg_key]
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
    for each_link in marker_to_ctg_linkage_num_dict:
        current_linking_reads = marker_to_ctg_linking_reads_dict[each_link]
        linking_reads_rd1_handle.write('%s\t%s\n' % (each_link, ','.join(current_linking_reads)))
    linking_reads_rd1_handle.close()


    ############################################## get pairwise_16s_iden_dict ##############################################

    report_and_log(('Round 1: Get pairwise 16S rRNA gene identities'), pwd_log_file, keep_quiet)

    # makeblastdn with marker gene sequences
    makeblastdb_16s_cmd = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, input_16s_polished)
    os.system(makeblastdb_16s_cmd)

    all_vs_all_16s_blastn_cmd = '%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, input_16s_polished, input_16s_polished, blast_results_all_vs_all_16s, blast_parameters)
    os.system(all_vs_all_16s_blastn_cmd)

    pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

    # write out to file
    pairwise_marker_similarity_handle = open(pairwise_marker_similarity, 'w')
    pairwise_marker_similarity_handle.write('Marker1\tMarker2\tSimilarity\n')
    for marker_pair in pairwise_16s_iden_dict:
        pairwise_marker_similarity_handle.write('%s\t%s\n' % ('\t'.join(marker_pair.split('__|__')), pairwise_16s_iden_dict[marker_pair]))
    pairwise_marker_similarity_handle.close()


    ####################################### filter_linkages_iteratively ########################################

    report_and_log(('Round 1: filtering linkages iteratively'), pwd_log_file, keep_quiet)

    # sort file in
    sort_csv_by_col(link_stats_combined, link_stats_combined_sorted, 'Number')

    filter_linkages_iteratively_new(link_stats_combined_sorted, pairwise_16s_iden_dict, min_iden_16s, marker_len_dict,
                                    min_link_num, within_gnm_linkage_num_diff, link_stats_combined_filtered_s1,
                                    marker_to_gnm_linkage_cigar_dict_16s_side,
                                    marker_to_gnm_linkage_cigar_dict_ctg_side,
                                    marker_to_ctg_gnm_Key_connector)


    ####################################### get linking reads for visualization ########################################

    report_and_log(('Round 1: Extracting linking reads for visualization'), pwd_log_file, keep_quiet)

    os.mkdir(mafft_seq_folder)

    ctgs_to_extract = set()
    all_linking_reads_base_set = set()
    linking_reads_txt_handle = open(linking_reads_tab, 'w')
    for marker_to_ctg in marker_to_ctg_linkage_num_dict_min3:
        marker_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[0]
        ctg_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[1]
        linking_reads = marker_to_ctg_linking_reads_dict[marker_to_ctg]
        all_linking_reads_base_set.update(linking_reads)
        ctgs_to_extract.add(ctg_id)
        linking_reads_txt_handle.write('%s\t%s\t%s\n' % (marker_id, ctg_id, ','.join(linking_reads)))
    linking_reads_txt_handle.close()

    # write out id of linking reads for extraction
    with open(linking_reads_r1_txt, 'w') as linking_reads_r1_txt_handle:
        linking_reads_r1_txt_handle.write('%s\n' % '\n'.join(sorted([('%s.1' % i) for i in all_linking_reads_base_set])))
    with open(linking_reads_r2_txt, 'w') as linking_reads_r2_txt_handle:
        linking_reads_r2_txt_handle.write('%s\n' % '\n'.join(sorted([('%s.2' % i) for i in all_linking_reads_base_set])))

    # extract linking reads with seqtk
    seqtk_extract_cmd_rd1_linking_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, linking_reads_r1_txt, linking_reads_r1_fasta)
    seqtk_extract_cmd_rd1_linking_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, linking_reads_r2_txt, linking_reads_r2_fasta)
    os.system(seqtk_extract_cmd_rd1_linking_r1)
    os.system(seqtk_extract_cmd_rd1_linking_r2)

    # subset combined genome file
    linked_contigs_txt_handle = open(linked_contigs_txt, 'w')
    linked_contigs_txt_handle.write('\n'.join(ctgs_to_extract) + '\n')
    linked_contigs_txt_handle.close()
    subset_linked_ctgs_cmd = 'seqtk subseq %s %s > %s' % (combined_input_gnms, linked_contigs_txt, linked_contigs_fasta)
    os.system(subset_linked_ctgs_cmd)

    # read sequence of linked contigs into dict
    linked_ctg_seq_dict = {}
    for linked_ctg in SeqIO.parse(linked_contigs_fasta, 'fasta'):
        linked_ctg_seq_dict[linked_ctg.id] = str(linked_ctg.seq)

    # read sequence of linking reads into dict
    linking_read_seq_dict = {}
    for linking_r1 in SeqIO.parse(linking_reads_r1_fasta, 'fasta'):
        linking_read_seq_dict[linking_r1.id] = str(linking_r1.seq)
    for linking_r2 in SeqIO.parse(linking_reads_r2_fasta, 'fasta'):
        linking_read_seq_dict[linking_r2.id] = str(linking_r2.seq)

    # read marker sequences into dict
    marker_seq_dict = {}
    for each_16s in SeqIO.parse(marker_gene_seqs, 'fasta'):
        marker_seq_dict[each_16s.id] = str(each_16s.seq)

    ########## prepare seq with multi-processing ##########

    argument_lol_for_linkage_vis_worker = []
    for marker_to_ctg in marker_to_ctg_linkage_num_dict_min3:

        reads_file_base_tmp = marker_to_ctg.replace(marker_to_ctg_gnm_Key_connector, '___')
        reads_file_base = reads_file_base_tmp.replace(gnm_to_ctg_connector, '___')

        linking_reads = marker_to_ctg_linking_reads_dict[marker_to_ctg]
        marker_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[0]
        marker_seq = marker_seq_dict[marker_id]
        ctg_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[1]
        contig_seq = linked_ctg_seq_dict[ctg_id]

        # create sub folders
        os.mkdir('%s/%s' % (mafft_seq_folder, reads_file_base))
        vis_reads_file_r1 = '%s/%s/%s_R1.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
        vis_reads_file_r2 = '%s/%s/%s_R2.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)

        # get marker_pos_list and contig_pos_list and write out sequences of linking reads
        marker_pos_list = []
        contig_pos_list = []
        vis_reads_file_r1_handle = open(vis_reads_file_r1, 'w')
        vis_reads_file_r2_handle = open(vis_reads_file_r2, 'w')
        for each_linking_read in linking_reads:
            linking_r1_id = '%s.1' % each_linking_read
            linking_r1_seq = linking_read_seq_dict[linking_r1_id]
            linking_r2_id = '%s.2' % each_linking_read
            linking_r2_seq = linking_read_seq_dict[linking_r2_id]
            vis_reads_file_r1_handle.write('>%s\n' % linking_r1_id)
            vis_reads_file_r1_handle.write('%s\n' % linking_r1_seq)
            vis_reads_file_r2_handle.write('>%s\n' % linking_r2_id)
            vis_reads_file_r2_handle.write('%s\n' % linking_r2_seq)
            # get matched position on makrer and contig
            marker_pos_r1 = list(MappingRecord_dict[each_linking_read].r1_16s_ref_dict.get(marker_id, dict()).keys())
            marker_pos_r2 = list(MappingRecord_dict[each_linking_read].r2_16s_ref_dict.get(marker_id, dict()).keys())
            contig_pos_r1 = list(MappingRecord_dict[each_linking_read].r1_ctg_ref_dict.get(ctg_id, dict()).keys())
            contig_pos_r2 = list(MappingRecord_dict[each_linking_read].r2_ctg_ref_dict.get(ctg_id, dict()).keys())
            if len(marker_pos_r1) == 1:
                marker_pos_list.append(marker_pos_r1[0])
            if len(marker_pos_r2) == 1:
                marker_pos_list.append(marker_pos_r2[0])
            if len(contig_pos_r1) == 1:
                contig_pos_list.append(contig_pos_r1[0])
            if len(contig_pos_r2) == 1:
                contig_pos_list.append(contig_pos_r2[0])
        vis_reads_file_r1_handle.close()
        vis_reads_file_r2_handle.close()

        argument_lol_for_linkage_vis_worker.append([reads_file_base, mafft_seq_folder,
                                                    marker_seq, contig_seq,
                                                    end_ctg_len_for_mafft, gap_N_num, bowtie_parameter,
                                                    marker_pos_list, contig_pos_list])

    # visualize linkages
    report_and_log(('Round 1: visualizing %s linkages with %s threads' % (len(argument_lol_for_linkage_vis_worker), num_threads)), pwd_log_file, keep_quiet)
    vis_linkages_pool = mp.Pool(processes=num_threads)
    vis_linkages_pool.map(linkage_vis_worker, argument_lol_for_linkage_vis_worker)
    vis_linkages_pool.close()
    vis_linkages_pool.join()


    ####################################################################################################################
    ############################################### second round linking ###############################################
    ####################################################################################################################

    ################################################# define file name #################################################

    combined_1st_round_unlinked_mags                = '%s/round_1_unlinked_gnm.fa'                      % step_2_wd
    combined_1st_round_unlinked_mag_end_seq         = '%s/round_1_unlinked_gnm_end_%sbp.fa'             % (step_2_wd, end_seq_len)
    rd1_unlinked_mag_end_seq_no_ext                 = '%s/round_1_unlinked_gnm_end_%sbp'                % (step_2_wd, end_seq_len)
    rd1_unlinked_mags_sam_bowtie_log                = '%s/round_1_unlinked_gnm_bowtie.log'              % step_2_wd
    rd1_unlinked_mags_sam_bowtie                    = '%s/round_1_unlinked_gnm_bowtie.sam'              % step_2_wd
    rd1_unlinked_mags_sam_bowtie_reformat           = '%s/round_1_unlinked_gnm_bowtie_reformatted.sam'  % step_2_wd
    rd1_unlinked_mags_sam_bowtie_reformat_log       = '%s/round_1_unlinked_gnm_bowtie_reformat.log'     % step_2_wd
    rd1_unlinked_mags_sam_bowtie_reformat_sorted    = '%s/round_1_unlinked_gnm_bowtie_reformatted_sorted.sam' % step_2_wd
    rd1_unlinked_mags_sam_line_num                  = '%s/round_1_unlinked_gnm_line_num.txt'            % step_2_wd
    rd1_unlinked_mags_sam_split_folder              = '%s/round_1_unlinked_gnm_split'                   % step_2_wd
    rd1_unlinked_mags_sam_MappingRecord_folder      = '%s/round_1_unlinked_gnm_MappingRecord'           % step_2_wd
    stats_GapFilling_file                           = '%s/stats_GapFilling_gnm.txt'                     % step_2_wd
    stats_GapFilling_file_filtered                  = '%s/stats_GapFilling_gnm_filtered.txt'            % step_2_wd
    free_living_16s_ref_file                        = '%s/round2_free_living_16s_refs.txt'              % step_2_wd
    free_living_ctg_ref_file                        = '%s/round2_free_living_ctg_refs.txt'              % step_2_wd
    free_living_all                                 = '%s/round2_free_living_all.fa'                    % step_2_wd
    free_living_all_fq_r1                           = '%s/round2_free_living_all_R1.fastq'              % step_2_wd
    free_living_all_fq_r2                           = '%s/round2_free_living_all_R2.fastq'              % step_2_wd
    free_living_all_fq                              = '%s/round2_free_living_all.fastq'                 % step_2_wd
    spades_wd                                       = '%s/mini_assembly_SPAdes_wd'                      % step_2_wd
    spades_log                                      = '%s/SPAdes_stdout.txt'                            % step_2_wd
    mira_manifest                                   = '%s/mira_manifest.txt'                            % step_2_wd
    mira_stdout                                     = '%s/mira_stdout.txt'                              % step_2_wd
    mini_assembly_to_16s_reads                      = '%s/mini_assembly_to_16s_reads.txt'               % step_2_wd
    mini_assembly_to_ctg_reads                      = '%s/mini_assembly_to_ctg_reads.txt'               % step_2_wd
    stats_GapFilling_ctg                            = '%s/stats_GapFilling_ctg.txt'                     % step_2_wd
    sam_file_mini_assembly                          = '%s/scaffolds_bowtie.sam'                         % step_2_wd
    sam_file_mini_assembly_log                      = '%s/scaffolds_bowtie.log'                         % step_2_wd
    sam_file_mini_assembly_reformatted              = '%s/scaffolds_bowtie_reformatted.sam'             % step_2_wd
    sam_file_mini_assembly_reformatted_log          = '%s/scaffolds_bowtie_reformatted.log'             % step_2_wd
    rd2_to_extract_flking_16s_r1_id                 = '%s/rd2_to_extract_flking_16s_r1_id.txt'          % step_2_wd
    rd2_to_extract_flking_16s_r2_id                 = '%s/rd2_to_extract_flking_16s_r2_id.txt'          % step_2_wd
    rd2_extracted_flking_16s_r1_seq_tmp             = '%s/rd2_extracted_flking_16s_r1_tmp.fa'           % step_2_wd
    rd2_extracted_flking_16s_r2_seq_tmp             = '%s/rd2_extracted_flking_16s_r2_tmp.fa'           % step_2_wd
    rd2_to_extract_flking_ctg_r1_id                 = '%s/rd2_to_extract_flking_ctg_r1_id.txt'          % step_2_wd
    rd2_to_extract_flking_ctg_r2_id                 = '%s/rd2_to_extract_flking_ctg_r2_id.txt'          % step_2_wd
    rd2_extracted_flking_ctg_r1_seq_tmp             = '%s/rd2_extracted_flking_ctg_r1_tmp.fa'           % step_2_wd
    rd2_extracted_flking_ctg_r2_seq_tmp             = '%s/rd2_extracted_flking_ctg_r2_tmp.fa'           % step_2_wd
    rd2_extracted_r1_combined                       = '%s/rd2_extracted_r1_combined.fa'                 % step_2_wd
    rd2_extracted_r2_combined                       = '%s/rd2_extracted_r2_combined.fa'                 % step_2_wd
    rd2_to_extract_up_r1_id                         = '%s/rd2_to_extract_UP_R1_id.txt'                  % step_2_wd
    rd2_to_extract_up_r2_id                         = '%s/rd2_to_extract_UP_R2_id.txt'                  % step_2_wd
    rd2_extracted_r1_combined_up                    = '%s/rd2_extracted_r1_combined_up.fa'              % step_2_wd
    rd2_extracted_r2_combined_up                    = '%s/rd2_extracted_r2_combined_up.fa'              % step_2_wd
    rd2_extracted_combined_up                       = '%s/rd2_extracted_combined_up.fa'                 % step_2_wd

    rd2_read_to_extract_flanking_both_r1_up_id_file = '%s/rd2_read_to_extract_flanking_both_R1_up_id.txt'   % step_2_wd
    rd2_read_to_extract_flanking_both_r2_up_id_file = '%s/rd2_read_to_extract_flanking_both_R2_up_id.txt'   % step_2_wd
    rd2_read_extracted_flanking_both_r1_up_seq      = '%s/rd2_read_to_extract_flanking_both_R1_up.fa'       % step_2_wd
    rd2_read_extracted_flanking_both_r2_up_seq      = '%s/rd2_read_to_extract_flanking_both_R2_up.fa'       % step_2_wd
    rd2_read_extracted_flanking_both_r12_up_seq     = '%s/rd2_read_to_extract_flanking_both_R12_up.fa'      % step_2_wd
    rd2_read_extracted_flanking_16s_r12_up_seq      = '%s/rd2_read_extracted_flanking_16s_r12_up.fa'        % step_2_wd
    rd2_read_extracted_flanking_ctg_r12_up_seq      = '%s/rd2_read_extracted_flanking_ctg_r12_up.fa'        % step_2_wd

    sam_file_mini_assembly_16s                      = '%s/scaffolds_bowtie_16s.sam'                         % step_2_wd
    sam_file_mini_assembly_ctg                      = '%s/scaffolds_bowtie_ctg.sam'                         % step_2_wd
    sam_file_mini_assembly_log_16s                  = '%s/scaffolds_bowtie_16s.log'                         % step_2_wd
    sam_file_mini_assembly_log_ctg                  = '%s/scaffolds_bowtie_ctg.log'                         % step_2_wd
    sam_file_mini_assembly_reformatted_16s          = '%s/scaffolds_bowtie_reformatted_16s.sam'             % step_2_wd
    sam_file_mini_assembly_reformatted_ctg          = '%s/scaffolds_bowtie_reformatted_ctg.sam'             % step_2_wd
    sam_file_mini_assembly_reformatted_log_16s      = '%s/scaffolds_bowtie_reformatted_16s.log'             % step_2_wd
    sam_file_mini_assembly_reformatted_log_ctg      = '%s/scaffolds_bowtie_reformatted_ctg.log'             % step_2_wd
    sam_file_mini_assembly_reformatted_sorted_16s   = '%s/scaffolds_bowtie_reformatted_sorted_16s.sam'      % step_2_wd
    sam_file_mini_assembly_reformatted_sorted_ctg   = '%s/scaffolds_bowtie_reformatted_sorted_ctg.sam'      % step_2_wd

    rd2_with_both_mates = False
    filter_mini_assembly_separately = False

    #################### get the sequences of 1st round unlinked marker genes and genomic sequences ####################

    os.mkdir(step_2_wd)
    report_and_log(('Round 2: get unlinked marker genes and genomes'), pwd_log_file, keep_quiet)

    # get linked marker genes and genomic sequences in step 1
    linked_marker_gene_set = set()
    linked_genomic_seq_set = set()
    linked_16s_to_gnm_set = set()
    for each_link in open(link_stats_combined_filtered_s1):
        if not each_link.startswith('MarkerGene,GenomicSeq,Number'):
            each_link_split = each_link.strip().split(',')
            id_16s = each_link_split[0][12:]
            id_gnm = each_link_split[1][12:]
            linked_marker_gene_set.add(id_16s)
            linked_genomic_seq_set.add(id_gnm)
            linked_16s_to_gnm_set.add('%s%s%s' % (id_16s, marker_to_ctg_gnm_Key_connector, id_gnm))

    # get the sequence of unlinked genomic seqs
    renamed_gnm_re = '%s/*.%s' % (prefixed_mag_folder, mag_file_extension)
    renamed_gnm_list = [os.path.basename(file_name) for file_name in glob.glob(renamed_gnm_re)]
    renamed_gnm_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in renamed_gnm_list]

    # keep only unlinked mags
    unlinked_mag_list_with_pwd = []
    for renamed_mag in renamed_gnm_list_no_ext:
        if renamed_mag not in linked_genomic_seq_set:
            pwd_renamed_mag = '%s/%s.%s' % (prefixed_mag_folder, renamed_mag, mag_file_extension)
            unlinked_mag_list_with_pwd.append(pwd_renamed_mag)

    # combine unlinked mags
    cat_cmd = 'cat %s > %s' % (' '.join(unlinked_mag_list_with_pwd), combined_1st_round_unlinked_mags)
    os.system(cat_cmd)


    ######################################## extract reads matched to unlinked 16S #######################################

    # get all rd1 linking reads
    all_linking_reads_for_rd1_linkages = set()
    for each_link in marker_to_ctg_linkage_num_dict_min3_passed_ratio_check:
        current_linking_reads = marker_to_ctg_linking_reads_dict[each_link]
        marker_to_gnm_key = each_link.split(gnm_to_ctg_connector)[0]
        if marker_to_gnm_key in linked_16s_to_gnm_set:
            all_linking_reads_for_rd1_linkages.update(current_linking_reads)

    # get the id of reads to extract
    free_living_16s_ref_file_handle = open(free_living_16s_ref_file, 'w')
    rd2_read_to_extract_flanking_16s = set()
    rd2_read_to_extract_flanking_16s_r1_up = set()
    rd2_read_to_extract_flanking_16s_r2_up = set()
    for each_mp in MappingRecord_dict.copy():
        if each_mp in all_linking_reads_for_rd1_linkages:
            MappingRecord_dict.pop(each_mp)
        else:
            r1_16s_refs = MappingRecord_dict[each_mp].r1_16s_refs_no_ignored
            r2_16s_refs = MappingRecord_dict[each_mp].r2_16s_refs_no_ignored
            if len(MappingRecord_dict[each_mp].shared_16s_refs_no_ignored) > 0:
                MappingRecord_dict.pop(each_mp)
            elif (len(r1_16s_refs) > 0) or (len(r2_16s_refs) == 0):
                rd2_read_to_extract_flanking_16s.add(each_mp)
                rd2_read_to_extract_flanking_16s_r2_up.add('%s.2' % each_mp)
                if rd2_with_both_mates is True:
                    free_living_16s_ref_file_handle.write('%s\t%s\n' % (each_mp, ','.join(r1_16s_refs)))
                else:
                    free_living_16s_ref_file_handle.write('%s.2\t%s\n' % (each_mp, ','.join(r1_16s_refs)))

            elif (len(r1_16s_refs) == 0) or (len(r2_16s_refs) > 0):
                rd2_read_to_extract_flanking_16s.add(each_mp)
                rd2_read_to_extract_flanking_16s_r1_up.add('%s.1' % each_mp)
                if rd2_with_both_mates is True:
                    free_living_16s_ref_file_handle.write('%s\t%s\n' % (each_mp, ','.join(r2_16s_refs)))
                else:
                    free_living_16s_ref_file_handle.write('%s.1\t%s\n' % (each_mp, ','.join(r2_16s_refs)))
    free_living_16s_ref_file_handle.close()

    # write out id of linking reads for extraction
    if rd2_with_both_mates is True:
        with open(rd2_to_extract_flking_16s_r1_id, 'w') as rd2_to_extract_flking_16s_r1_id_handle:
            rd2_to_extract_flking_16s_r1_id_handle.write('%s\n' % '\n'.join(sorted([('%s.1' % i) for i in rd2_read_to_extract_flanking_16s])))
        with open(rd2_to_extract_flking_16s_r2_id, 'w') as rd2_to_extract_flking_16s_r2_id_handle:
            rd2_to_extract_flking_16s_r2_id_handle.write('%s\n' % '\n'.join(sorted([('%s.2' % i) for i in rd2_read_to_extract_flanking_16s])))

        # extract reads with seqtk
        seqtk_extract_cmd_rd2_flk_16s_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_to_extract_flking_16s_r1_id, rd2_extracted_flking_16s_r1_seq_tmp)
        seqtk_extract_cmd_rd2_flk_16s_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_to_extract_flking_16s_r2_id, rd2_extracted_flking_16s_r2_seq_tmp)
        os.system(seqtk_extract_cmd_rd2_flk_16s_r1)
        os.system(seqtk_extract_cmd_rd2_flk_16s_r2)


    ######################################## extract sequences flanking ctg ends #######################################

    report_and_log(('Round 2: Mapping input reads to the ends of contigs in round 1 unlinked genomes'), pwd_log_file, keep_quiet)

    get_unlinked_mag_end_seq(combined_1st_round_unlinked_mags, combined_1st_round_unlinked_mag_end_seq, end_seq_len)

    # index reference
    bowtie_build_unlinked_ctg_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, combined_1st_round_unlinked_mag_end_seq, rd1_unlinked_mag_end_seq_no_ext)
    os.system(bowtie_build_unlinked_ctg_cmd)

    # mapping with bowtie
    bowtie_cmd_unlinked_ctg = 'bowtie2 -x %s -U %s,%s -S %s -p %s -f %s 2> %s' % (rd1_unlinked_mag_end_seq_no_ext, reads_file_r1_fasta, reads_file_r2_fasta, rd1_unlinked_mags_sam_bowtie, num_threads, bowtie_parameter, rd1_unlinked_mags_sam_bowtie_log)
    os.system(bowtie_cmd_unlinked_ctg)

    # convert cigar format
    report_and_log(('Round 2: transforming cigar format from 1.3 to 1.4'), pwd_log_file, keep_quiet)
    bbmap_reformat_rd1_unlinked_mags_sam  = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (rd1_unlinked_mags_sam_bowtie, rd1_unlinked_mags_sam_bowtie_reformat, rd1_unlinked_mags_sam_bowtie_reformat_log)
    os.system(bbmap_reformat_rd1_unlinked_mags_sam)
    os.remove(rd1_unlinked_mags_sam_bowtie)

    # sort sam file first
    report_and_log(('Round 2: sorting mappping results'), pwd_log_file, keep_quiet)
    sort_by_read_cmd = 'samtools sort -n -O sam --threads %s -o %s %s' % (num_threads, rd1_unlinked_mags_sam_bowtie_reformat_sorted, rd1_unlinked_mags_sam_bowtie_reformat)
    os.system(sort_by_read_cmd)
    os.remove(rd1_unlinked_mags_sam_bowtie_reformat)


    ##################################################### read in sam file ####################################################

    round_2_ctg_end_seq_len_dict = {}
    for each_ctg_end_record in SeqIO.parse(combined_1st_round_unlinked_mag_end_seq, 'fasta'):
        round_2_ctg_end_seq_len_dict[each_ctg_end_record.id] = len(each_ctg_end_record.seq)

    # get the number of lines per file
    report_and_log(('Round 2: calculating the number of lines per subset'), pwd_log_file, keep_quiet)
    os.system('wc -l %s > %s' % (rd1_unlinked_mags_sam_bowtie_reformat_sorted, rd1_unlinked_mags_sam_line_num))
    rd1_unlinked_mag_sam_line_num = int(open(rd1_unlinked_mags_sam_line_num).readline().strip().split(' ')[0])
    os.remove(rd1_unlinked_mags_sam_line_num)
    line_num_per_file = int(round(rd1_unlinked_mag_sam_line_num/(num_threads*10))) + 10

    report_and_log(('Round 2: splitting sam file'), pwd_log_file, keep_quiet)
    os.mkdir(rd1_unlinked_mags_sam_split_folder)
    split_gnm_sam_cmd = 'split -l %s %s %s/splitted_sam_ ' % (line_num_per_file, rd1_unlinked_mags_sam_bowtie_reformat_sorted, rd1_unlinked_mags_sam_split_folder)
    os.system(split_gnm_sam_cmd)

    report_and_log(('Round 2: analysing mappping results with %s threads' % num_threads), pwd_log_file, keep_quiet)
    os.mkdir(rd1_unlinked_mags_sam_MappingRecord_folder)

    # get splitted sam file list
    splitted_gnm_sam_file_list = [os.path.basename(file_name) for file_name in glob.glob('%s/splitted_sam_*' % rd1_unlinked_mags_sam_split_folder)]

    # prepare lol for mp worker
    list_for_parse_sam_gnm_worker = []
    splitted_sam_mp_file_set = set()
    for splitted_gnm_sam_file in splitted_gnm_sam_file_list:
        pwd_splitted_gnm_sam_file                       = '%s/%s'                           % (rd1_unlinked_mags_sam_split_folder, splitted_gnm_sam_file)
        pwd_splitted_gnm_sam_free_living_ctg_ref_file   = '%s/%s_free_living_ctg_refs.txt'  % (rd1_unlinked_mags_sam_MappingRecord_folder, splitted_gnm_sam_file)
        splitted_sam_mp_file_set.add(pwd_splitted_gnm_sam_free_living_ctg_ref_file)
        list_for_parse_sam_gnm_worker.append([pwd_splitted_gnm_sam_file,
                                              pwd_splitted_gnm_sam_free_living_ctg_ref_file,
                                              min_M_len_ctg,
                                              mismatch_cutoff,
                                              round_2_ctg_end_seq_len_dict,
                                              rd2_with_both_mates])

    pool_parse_sam_gnm = mp.Pool(processes=num_threads)
    pool_parse_sam_gnm.map(parse_sam_gnm_worker, list_for_parse_sam_gnm_worker)
    pool_parse_sam_gnm.close()
    pool_parse_sam_gnm.join()

    report_and_log(('Round 2: removing splitted subsets from disk'), pwd_log_file, keep_quiet)
    os.system('rm -r %s' % rd1_unlinked_mags_sam_split_folder)

    # combine free_living_ctg_ref_files
    os.system('cat %s > %s' % (' '.join(splitted_sam_mp_file_set), free_living_ctg_ref_file))

    ##################################################### extract reads flanking contig ends ####################################################

    if rd2_with_both_mates is True:
        to_extract_read_base_rd2_ctg = set()
        for each_read_base in open(free_living_ctg_ref_file):
            to_extract_read_base_rd2_ctg.add(each_read_base.strip().split('\t')[0])

        # write out id of linking reads for extraction
        with open(rd2_to_extract_flking_ctg_r1_id, 'w') as rd2_to_extract_flking_ctg_r1_id_handle:
            rd2_to_extract_flking_ctg_r1_id_handle.write('%s\n' % '\n'.join(sorted([('%s.1' % i) for i in to_extract_read_base_rd2_ctg])))
        with open(rd2_to_extract_flking_ctg_r2_id, 'w') as rd2_to_extract_flking_ctg_r2_id_handle:
            rd2_to_extract_flking_ctg_r2_id_handle.write('%s\n' % '\n'.join(sorted([('%s.2' % i) for i in to_extract_read_base_rd2_ctg])))

        # extract reads with seqtk
        seqtk_extract_cmd_rd2_flk_ctg_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_to_extract_flking_ctg_r1_id, rd2_extracted_flking_ctg_r1_seq_tmp)
        seqtk_extract_cmd_rd2_flk_ctg_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_to_extract_flking_ctg_r2_id, rd2_extracted_flking_ctg_r2_seq_tmp)
        os.system(seqtk_extract_cmd_rd2_flk_ctg_r1)
        os.system(seqtk_extract_cmd_rd2_flk_ctg_r2)

        # read extracted read sequences into dict
        extract_rd2_read_seq_dict = {}
        for extracted_r1 in SeqIO.parse(rd2_extracted_flking_16s_r1_seq_tmp, 'fasta'):
            extract_rd2_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
        for extracted_r2 in SeqIO.parse(rd2_extracted_flking_16s_r2_seq_tmp, 'fasta'):
            extract_rd2_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)
        for extracted_r1 in SeqIO.parse(rd2_extracted_flking_ctg_r1_seq_tmp, 'fasta'):
            extract_rd2_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
        for extracted_r2 in SeqIO.parse(rd2_extracted_flking_ctg_r2_seq_tmp, 'fasta'):
            extract_rd2_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)

        # write out paired in the same order
        rd2_extracted_r1_combined_handle = open(rd2_extracted_r1_combined, 'w')
        rd2_extracted_r2_combined_handle = open(rd2_extracted_r2_combined, 'w')
        for each_read_base in set.union(rd2_read_to_extract_flanking_16s, to_extract_read_base_rd2_ctg):
            current_r1 = '%s.1' % each_read_base
            current_r2 = '%s.2' % each_read_base
            current_r1_seq = extract_rd2_read_seq_dict.get(current_r1, '')
            current_r2_seq = extract_rd2_read_seq_dict.get(current_r2, '')
            rd2_extracted_r1_combined_handle.write('>%s\n' % current_r1)
            rd2_extracted_r1_combined_handle.write('%s\n' % current_r1_seq)
            rd2_extracted_r2_combined_handle.write('>%s\n' % current_r2)
            rd2_extracted_r2_combined_handle.write('%s\n' % current_r2_seq)
        rd2_extracted_r1_combined_handle.close()
        rd2_extracted_r2_combined_handle.close()

    else:
        rd2_read_to_extract_flanking_ctg_r1_up = set()
        rd2_read_to_extract_flanking_ctg_r2_up = set()
        for each_read_to_ctg_ref in open(free_living_ctg_ref_file):
            read_id = each_read_to_ctg_ref.strip().split('\t')[0]
            if read_id.endswith('.1') is True:
                rd2_read_to_extract_flanking_ctg_r1_up.add(read_id)
            if read_id.endswith('.2') is True:
                rd2_read_to_extract_flanking_ctg_r2_up.add(read_id)

        # write out id of linking reads for extraction
        rd2_read_to_extract_flanking_both_r1_up = set.union(rd2_read_to_extract_flanking_16s_r1_up, rd2_read_to_extract_flanking_ctg_r1_up)
        rd2_read_to_extract_flanking_both_r2_up = set.union(rd2_read_to_extract_flanking_16s_r2_up, rd2_read_to_extract_flanking_ctg_r2_up)
        with open(rd2_read_to_extract_flanking_both_r1_up_id_file, 'w') as rd2_read_to_extract_flanking_both_r1_up_id_file_handle:
            rd2_read_to_extract_flanking_both_r1_up_id_file_handle.write('%s\n' % '\n'.join(sorted([i for i in rd2_read_to_extract_flanking_both_r1_up])))
        with open(rd2_read_to_extract_flanking_both_r2_up_id_file, 'w') as rd2_read_to_extract_flanking_both_r2_up_id_file_handle:
            rd2_read_to_extract_flanking_both_r2_up_id_file_handle.write('%s\n' % '\n'.join(sorted([i for i in rd2_read_to_extract_flanking_both_r2_up])))

        # extract reads with seqtk
        seqtk_extract_cmd_rd2_flk_both_up_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_read_to_extract_flanking_both_r1_up_id_file, rd2_read_extracted_flanking_both_r1_up_seq)
        seqtk_extract_cmd_rd2_flk_both_up_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_read_to_extract_flanking_both_r2_up_id_file, rd2_read_extracted_flanking_both_r2_up_seq)
        os.system(seqtk_extract_cmd_rd2_flk_both_up_r1)
        os.system(seqtk_extract_cmd_rd2_flk_both_up_r2)

        # combine unpaired r1 and r2 for assembly
        os.system('cat %s %s > %s' % (rd2_read_extracted_flanking_both_r1_up_seq, rd2_read_extracted_flanking_both_r2_up_seq, rd2_read_extracted_flanking_both_r12_up_seq))

        # read extracted reads into dict
        rd2_read_seq_dict = dict()
        for each_read in SeqIO.parse(rd2_read_extracted_flanking_both_r12_up_seq, 'fasta'):
            rd2_read_seq_dict[each_read.id] = str(each_read.seq)

        # write out flanking_16s and flanking_ctg reads separately
        rd2_read_extracted_flanking_16s_r12_up_seq_handle = open(rd2_read_extracted_flanking_16s_r12_up_seq, 'w')
        rd2_read_extracted_flanking_ctg_r12_up_seq_handle = open(rd2_read_extracted_flanking_ctg_r12_up_seq, 'w')
        for each_read in rd2_read_seq_dict:
            read_seq = rd2_read_seq_dict[each_read]
            if (each_read in rd2_read_to_extract_flanking_16s_r1_up) or (each_read in rd2_read_to_extract_flanking_16s_r2_up):
                rd2_read_extracted_flanking_16s_r12_up_seq_handle.write('>%s\n' % each_read)
                rd2_read_extracted_flanking_16s_r12_up_seq_handle.write('%s\n' % read_seq)
            if (each_read in rd2_read_to_extract_flanking_ctg_r1_up) or (each_read in rd2_read_to_extract_flanking_ctg_r2_up):
                rd2_read_extracted_flanking_ctg_r12_up_seq_handle.write('>%s\n' % each_read)
                rd2_read_extracted_flanking_ctg_r12_up_seq_handle.write('%s\n' % read_seq)
        rd2_read_extracted_flanking_16s_r12_up_seq_handle.close()
        rd2_read_extracted_flanking_ctg_r12_up_seq_handle.close()

        # clean from memory
        rd2_read_seq_dict = dict()

    ######################################### second round linking by assembly #########################################

    # assemble
    if round_2_mira is True:

        free_living_all_id_r1 = set()
        free_living_all_id_r2 = set()
        for each_read in open(free_living_all):
            if each_read.startswith('>'):
                read_id = each_read.strip()[1:].split(' ')[0]
                if read_id[-1] == '1':
                    free_living_all_id_r1.add(read_id)
                if read_id[-1] == '2':
                    free_living_all_id_r2.add(read_id)

        argument_list_r1 = [reads_file_r1, 'fastq', free_living_all_id_r1, free_living_all_fq_r1]
        argument_list_r2 = [reads_file_r2, 'fastq', free_living_all_id_r2, free_living_all_fq_r2]

        # extract reads with multiprocessing
        pool = mp.Pool(processes=2)
        pool.map(extract_reads_worker, [argument_list_r1, argument_list_r2])
        pool.close()
        pool.join()

        os.system('cat %s %s > %s' % (free_living_all_fq_r1, free_living_all_fq_r2, free_living_all_fq))
        report_and_log(('Round 2: running Mira on extracted reads'), pwd_log_file, keep_quiet)
        run_mira5(output_prefix, mira_tmp_dir, step_2_wd, mira_manifest, free_living_all_fq, mira_stdout, force_overwrite)
        mini_assemblies = '%s/%s_mira_est_no_chimera_assembly/%s_mira_est_no_chimera_d_results/%s_mira_est_no_chimera_out.unpadded.fasta' % (step_2_wd, output_prefix, output_prefix, output_prefix)
    else:
        report_and_log(('Round 2: running SPAdes on extracted reads'), pwd_log_file, keep_quiet)
        if rd2_with_both_mates is True:
            #spades_cmd = '%s --only-assembler --meta -1 %s -2 %s -o %s -t %s -k 55,75,99,127 > %s' % (pwd_spades_exe, rd2_extracted_r1_combined, rd2_extracted_r2_combined, spades_wd, num_threads, spades_log)
            spades_cmd = '%s --only-assembler --meta -1 %s -2 %s -o %s -t %s -k 55,75,99 > %s' % (pwd_spades_exe, rd2_extracted_r1_combined, rd2_extracted_r2_combined, spades_wd, num_threads, spades_log)
        else:
            spades_cmd = '%s --only-assembler -s %s -o %s -t %s -k 55,75,99,127 > %s' % (pwd_spades_exe, rd2_read_extracted_flanking_both_r12_up_seq, spades_wd, num_threads, spades_log)

        report_and_log((spades_cmd), pwd_log_file, True)
        os.system(spades_cmd)
        mini_assemblies = '%s/scaffolds.fasta' % spades_wd

    ##################################################### mapping reads to mini_assemblies ####################################################

    # index miniassembly
    mini_assemblies_no_ext = '.'.join(mini_assemblies.split('.')[:-1])
    bowtie_build_mini_assemblies_cmd = 'bowtie2-build --quiet --threads %s -f %s %s'   % (num_threads, mini_assemblies, mini_assemblies_no_ext)
    os.system(bowtie_build_mini_assemblies_cmd)

    # mapping with bowtie
    if rd2_with_both_mates is True:
        bowtie_cmd_miniassembly = 'bowtie2 -x %s -U %s,%s -S %s -p %s -f %s 2> %s'     % (mini_assemblies_no_ext, rd2_extracted_r1_combined, rd2_extracted_r2_combined, sam_file_mini_assembly, num_threads, bowtie_parameter_mini_assembly, sam_file_mini_assembly_log)
        os.system(bowtie_cmd_miniassembly)

        # reformat cigar string
        mini_assembly_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly, sam_file_mini_assembly_reformatted, sam_file_mini_assembly_reformatted_log)
        os.system(mini_assembly_reformat_cmd)
    else:
        bowtie_cmd_miniassembly_16s = 'bowtie2 -x %s -U %s -S %s -p %s -f %s 2> %s'    % (mini_assemblies_no_ext, rd2_read_extracted_flanking_16s_r12_up_seq, sam_file_mini_assembly_16s, num_threads, bowtie_parameter_mini_assembly, sam_file_mini_assembly_log_16s)
        bowtie_cmd_miniassembly_ctg = 'bowtie2 -x %s -U %s -S %s -p %s -f %s 2> %s'    % (mini_assemblies_no_ext, rd2_read_extracted_flanking_ctg_r12_up_seq, sam_file_mini_assembly_ctg, num_threads, bowtie_parameter_mini_assembly, sam_file_mini_assembly_log_ctg)
        os.system(bowtie_cmd_miniassembly_16s)
        os.system(bowtie_cmd_miniassembly_ctg)

        # reformat cigar string
        mini_assembly_reformat_cmd_16s = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly_16s, sam_file_mini_assembly_reformatted_16s, sam_file_mini_assembly_reformatted_log_16s)
        mini_assembly_reformat_cmd_ctg = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly_ctg, sam_file_mini_assembly_reformatted_ctg, sam_file_mini_assembly_reformatted_log_ctg)
        os.system(mini_assembly_reformat_cmd_16s)
        os.system(mini_assembly_reformat_cmd_ctg)

        # sort sam files
        sort_by_read_cmd_mini_assembly_16s = 'samtools sort -n -O sam --threads %s -o %s %s' % (num_threads, sam_file_mini_assembly_reformatted_sorted_16s, sam_file_mini_assembly_reformatted_16s)
        sort_by_read_cmd_mini_assembly_ctg = 'samtools sort -n -O sam --threads %s -o %s %s' % (num_threads, sam_file_mini_assembly_reformatted_sorted_ctg, sam_file_mini_assembly_reformatted_ctg)
        os.system(sort_by_read_cmd_mini_assembly_16s)
        os.system(sort_by_read_cmd_mini_assembly_ctg)


    ##################################################### read in sam file ####################################################

    report_and_log(('Round 2: read in sam file'), pwd_log_file, keep_quiet)

    if rd2_with_both_mates is True:

        mini_assembly_len_dict = {}
        mini_assembly_mp_dict = {}
        with open(sam_file_mini_assembly_reformatted) as sam_file_mini_assembly_reformatted_opened:
            for each_read in sam_file_mini_assembly_reformatted_opened:
                each_read_split = each_read.strip().split('\t')
                if each_read.startswith('@'):
                    mini_assembly_id = ''
                    mini_assembly_len = 0
                    for each_element in each_read_split:
                        if each_element.startswith('SN:'):
                            mini_assembly_id = each_element[3:]
                        if each_element.startswith('LN:'):
                            mini_assembly_len = int(each_element[3:])
                    mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
                else:
                    cigar = each_read_split[5]
                    if cigar != '*':
                        read_id = each_read_split[0]
                        read_id_base = '.'.join(read_id.split('.')[:-1])
                        read_strand = read_id.split('.')[-1]
                        ref_id = each_read_split[2]
                        ref_pos = int(each_read_split[3])

                        if read_id_base not in mini_assembly_mp_dict:
                            mini_assembly_mp_dict[read_id_base] = MappingRecord()

                        if read_strand == '1':
                            if ref_id not in mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict:
                                mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict[ref_id][ref_pos] = cigar

                        if read_strand == '2':
                            if ref_id not in mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict:
                                mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict[ref_id][ref_pos] = cigar


        #################################################### parse sam file ####################################################

        report_and_log(('Round 2: parsing sam file'), pwd_log_file, keep_quiet)

        # get mini assembly to 16s reads dict
        mini_assembly_to_16s_reads_dict = {}
        mini_assembly_to_ctg_reads_dict = {}
        for each_mp in mini_assembly_mp_dict:

            current_mp_record = mini_assembly_mp_dict[each_mp]
            mini_refs_to_ignore = set()

            ########## get lowest mismatch for r1/r2 16s refs ##########

            # get r1_ref_cigar_set
            r1_ref_cigar_set = set()
            for each_pos_dict in current_mp_record.r1_mini_ref_dict.values():
                each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                r1_ref_cigar_set.update(each_pos_dict_values)

            # get r2_ref_cigar_set
            r2_ref_cigar_set = set()
            for each_pos_dict in current_mp_record.r2_mini_ref_dict.values():
                each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                r2_ref_cigar_set.update(each_pos_dict_values)

            r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_mini)
            r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_mini)
            mini_assembly_mp_dict[each_mp].r1_mini_refs_lowest_mismatch = r1_ref_min_mismatch
            mini_assembly_mp_dict[each_mp].r2_mini_refs_lowest_mismatch = r2_ref_min_mismatch

            ########## filter r1 mini refs ##########

            r1_mini_refs_passed_qc = {}
            for r1_mini_ref in current_mp_record.r1_mini_ref_dict:
                r1_matched_pos_dict = current_mp_record.r1_mini_ref_dict[r1_mini_ref]

                # one read need to mapped to one 16S only for one time
                if len(r1_matched_pos_dict) > 1:
                    mini_refs_to_ignore.add(r1_mini_ref)
                else:
                    r1_mini_ref_pos = list(r1_matched_pos_dict.keys())[0]
                    r1_mini_ref_cigar = r1_matched_pos_dict[r1_mini_ref_pos]
                    r1_mini_ref_cigar_splitted = cigar_splitter(r1_mini_ref_cigar)

                    # check both end clip
                    both_end_clp = check_both_ends_clipping(r1_mini_ref_cigar_splitted)
                    if both_end_clp is True:
                        mini_refs_to_ignore.add(r1_mini_ref)
                    else:
                        # check mismatch
                        r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                            r1_mini_ref_cigar_splitted)
                        if r1_ref_min_mismatch == 'NA':
                            mini_refs_to_ignore.add(r1_mini_ref)
                        elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                            mini_refs_to_ignore.add(r1_mini_ref)
                        else:
                            # check aligned length
                            if r1_aligned_len < min_M_len_mini:
                                mini_refs_to_ignore.add(r1_mini_ref)
                            else:
                                # check if clp in the middle
                                clip_in_middle = False
                                if ('S' in r1_mini_ref_cigar) or ('s' in r1_mini_ref_cigar):
                                    clip_in_middle = True
                                    if (r1_mini_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_mini_ref_pos == 1):
                                        clip_in_middle = False
                                    if (r1_mini_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                        if (r1_mini_ref_pos + r1_aligned_len - 1) == mini_assembly_len_dict[r1_mini_ref]:
                                            clip_in_middle = False

                                # exclude the ref if clp in the middle is True
                                if clip_in_middle is True:
                                    mini_refs_to_ignore.add(r1_mini_ref)
                                else:
                                    r1_mini_refs_passed_qc[r1_mini_ref] = [r1_mini_ref_cigar]

            ########## filter r2 mini refs ##########

            r2_mini_refs_passed_qc = {}
            for r2_mini_ref in current_mp_record.r2_mini_ref_dict:
                r2_matched_pos_dict = current_mp_record.r2_mini_ref_dict[r2_mini_ref]

                # one read need to mapped to one 16S only once
                if len(r2_matched_pos_dict) > 1:
                    mini_refs_to_ignore.add(r2_mini_ref)
                else:
                    r2_mini_ref_pos = list(r2_matched_pos_dict.keys())[0]
                    r2_mini_ref_cigar = r2_matched_pos_dict[r2_mini_ref_pos]
                    r2_mini_ref_cigar_splitted = cigar_splitter(r2_mini_ref_cigar)

                    # check both end clip
                    both_end_clp = check_both_ends_clipping(r2_mini_ref_cigar_splitted)
                    if both_end_clp is True:
                        mini_refs_to_ignore.add(r2_mini_ref)
                    else:
                        # check mismatch
                        r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                            r2_mini_ref_cigar_splitted)
                        if r2_ref_min_mismatch == 'NA':
                            mini_refs_to_ignore.add(r2_mini_ref)
                        elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                            mini_refs_to_ignore.add(r2_mini_ref)
                        else:
                            # check aligned length
                            if r2_aligned_len < min_M_len_mini:
                                mini_refs_to_ignore.add(r2_mini_ref)
                            else:
                                # check if clp in the middle
                                clip_in_middle = False
                                if ('S' in r2_mini_ref_cigar) or ('s' in r2_mini_ref_cigar):
                                    clip_in_middle = True
                                    if (r2_mini_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_mini_ref_pos == 1):
                                        clip_in_middle = False
                                    if (r2_mini_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                        if (r2_mini_ref_pos + r2_aligned_len - 1) == mini_assembly_len_dict[r2_mini_ref]:
                                            clip_in_middle = False

                                # exclude the ref if clp in the middle is True
                                if clip_in_middle is True:
                                    mini_refs_to_ignore.add(r2_mini_ref)
                                else:
                                    r2_mini_refs_passed_qc[r2_mini_ref] = [r2_mini_ref_cigar]

            ####################################################################################################

            r1_mini_refs_no_ignored = {key: value for key, value in r1_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
            r2_mini_refs_no_ignored = {key: value for key, value in r2_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
            shared_mini_ref_dict = {key: [r1_mini_refs_no_ignored[key][0], r2_mini_refs_no_ignored[key][0]] for key in set(r1_mini_refs_no_ignored).intersection(set(r2_mini_refs_no_ignored))}

            # both r1 and r2 mapped to the same mini assembly
            if len(shared_mini_ref_dict) > 0:

                if (each_mp in rd2_read_to_extract_flanking_16s) and (each_mp not in to_extract_read_base_rd2_ctg):
                    for each_mini_assembly in shared_mini_ref_dict:
                        if each_mini_assembly not in mini_assembly_to_16s_reads_dict:
                            mini_assembly_to_16s_reads_dict[each_mini_assembly] = {each_mp}
                        else:
                            mini_assembly_to_16s_reads_dict[each_mini_assembly].add(each_mp)

                elif (each_mp not in rd2_read_to_extract_flanking_16s) and (each_mp in to_extract_read_base_rd2_ctg):
                    for each_mini_assembly in shared_mini_ref_dict:
                        if each_mini_assembly not in mini_assembly_to_ctg_reads_dict:
                            mini_assembly_to_ctg_reads_dict[each_mini_assembly] = {each_mp}
                        else:
                            mini_assembly_to_ctg_reads_dict[each_mini_assembly].add(each_mp)
    else:
        mini_assembly_to_16s_reads_dict = {}
        mini_assembly_to_ctg_reads_dict = {}

        ############################################# get mini_assembly_to_16s_reads_dict #############################################

        mini_assembly_to_16s_reads_dict = {}
        mini_assembly_len_dict = {}
        current_read = ''
        current_read_ref_dict = dict()
        with open(sam_file_mini_assembly_reformatted_sorted_16s) as sam_file_mini_assembly_reformatted_sorted_16s_opened:
            for each_read in sam_file_mini_assembly_reformatted_sorted_16s_opened:
                each_read_split = each_read.strip().split('\t')
                if each_read.startswith('@'):
                    mini_assembly_id = ''
                    mini_assembly_len = 0
                    for each_element in each_read_split:
                        if each_element.startswith('SN:'):
                            mini_assembly_id = each_element[3:]
                        if each_element.startswith('LN:'):
                            mini_assembly_len = int(each_element[3:])
                    mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
                else:
                    cigar = each_read_split[5]
                    read_id = each_read_split[0]
                    read_id_base = '.'.join(read_id.split('.')[:-1])
                    read_strand = read_id.split('.')[-1]
                    ref_id = each_read_split[2]
                    ref_pos = int(each_read_split[3])

                    if current_read == '':
                        current_read = read_id
                        if cigar != '*':
                            current_read_ref_dict[ref_id] = {ref_pos: cigar}

                    elif read_id == current_read:
                        if cigar != '*':
                            if ref_id not in current_read_ref_dict:
                                current_read_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_ref_dict[ref_id][ref_pos] = cigar
                    else:
                        ################################### analysis previous read refs ####################################

                        ########## get lowest mismatch for current read refs ##########

                        # get r1_ref_cigar_set
                        ref_cigar_set = set()
                        for each_pos_dict in current_read_ref_dict.values():
                            each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                            ref_cigar_set.update(each_pos_dict_values)

                        ref_min_mismatch = get_min_mismatch_from_cigar_list(ref_cigar_set, min_M_len_ctg)

                        ########## filter refs ##########

                        r1_16s_refs_passed_qc = {}
                        for each_ref in current_read_ref_dict:
                            matched_pos_dict = current_read_ref_dict[each_ref]

                            # one read need to mapped to one 16S only for one time
                            if len(matched_pos_dict) == 1:
                                ref_pos = list(matched_pos_dict.keys())[0]
                                ref_cigar = matched_pos_dict[ref_pos]
                                ref_cigar_splitted = cigar_splitter(ref_cigar)

                                # check both end clip
                                both_end_clp = check_both_ends_clipping(ref_cigar_splitted)
                                if both_end_clp is False:
                                    # check mismatch
                                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(ref_cigar_splitted)
                                    if (ref_min_mismatch != 'NA'):
                                        if (mismatch_pct <= ref_min_mismatch) and (mismatch_pct <= mismatch_cutoff):
                                            # check aligned length
                                            if aligned_len >= min_M_len_ctg:
                                                # check if clp in the middle
                                                clip_in_middle = False
                                                if ('S' in ref_cigar) or ('s' in ref_cigar):
                                                    clip_in_middle = True
                                                    if (ref_cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                                                        clip_in_middle = False
                                                    if (ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                        if (ref_pos + aligned_len - 1) == mini_assembly_len_dict[each_ref]:
                                                            clip_in_middle = False

                                                # exclude the ref if clp in the middle is True
                                                if clip_in_middle is False:
                                                    if each_ref not in mini_assembly_to_16s_reads_dict:
                                                        mini_assembly_to_16s_reads_dict[each_ref] = {current_read}
                                                    else:
                                                        mini_assembly_to_16s_reads_dict[each_ref].add(current_read)

                        ########################################### reset values ###########################################

                        current_read = read_id
                        current_read_ref_dict = dict()
                        if cigar != '*':
                            current_read_ref_dict[ref_id] = {ref_pos: cigar}


        ############################################# get mini_assembly_to_ctg_reads_dict #############################################

        mini_assembly_to_ctg_reads_dict = {}
        mini_assembly_len_dict = {}
        current_read = ''
        current_read_ref_dict = dict()
        with open(sam_file_mini_assembly_reformatted_sorted_ctg) as sam_file_mini_assembly_reformatted_sorted_ctg_opened:
            for each_read in sam_file_mini_assembly_reformatted_sorted_ctg_opened:
                each_read_split = each_read.strip().split('\t')
                if each_read.startswith('@'):
                    mini_assembly_id = ''
                    mini_assembly_len = 0
                    for each_element in each_read_split:
                        if each_element.startswith('SN:'):
                            mini_assembly_id = each_element[3:]
                        if each_element.startswith('LN:'):
                            mini_assembly_len = int(each_element[3:])
                    mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
                else:
                    cigar = each_read_split[5]
                    read_id = each_read_split[0]
                    read_id_base = '.'.join(read_id.split('.')[:-1])
                    read_strand = read_id.split('.')[-1]
                    ref_id = each_read_split[2]
                    ref_pos = int(each_read_split[3])

                    if current_read == '':
                        current_read = read_id
                        if cigar != '*':
                            current_read_ref_dict[ref_id] = {ref_pos: cigar}

                    elif read_id == current_read:
                        if cigar != '*':
                            if ref_id not in current_read_ref_dict:
                                current_read_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_ref_dict[ref_id][ref_pos] = cigar
                    else:
                        ################################### analysis previous read refs ####################################

                        ########## get lowest mismatch for current read refs ##########

                        # get r1_ref_cigar_set
                        ref_cigar_set = set()
                        for each_pos_dict in current_read_ref_dict.values():
                            each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                            ref_cigar_set.update(each_pos_dict_values)

                        ref_min_mismatch = get_min_mismatch_from_cigar_list(ref_cigar_set, min_M_len_ctg)

                        ########## filter refs ##########

                        r1_16s_refs_passed_qc = {}
                        for each_ref in current_read_ref_dict:
                            matched_pos_dict = current_read_ref_dict[each_ref]

                            # one read need to mapped to one 16S only for one time
                            if len(matched_pos_dict) == 1:
                                ref_pos = list(matched_pos_dict.keys())[0]
                                ref_cigar = matched_pos_dict[ref_pos]
                                ref_cigar_splitted = cigar_splitter(ref_cigar)

                                # check both end clip
                                both_end_clp = check_both_ends_clipping(ref_cigar_splitted)
                                if both_end_clp is False:
                                    # check mismatch
                                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(ref_cigar_splitted)
                                    if (ref_min_mismatch != 'NA'):
                                        if (mismatch_pct <= ref_min_mismatch) and (mismatch_pct <= mismatch_cutoff):
                                            # check aligned length
                                            if aligned_len >= min_M_len_ctg:

                                                # check if clp in the middle
                                                clip_in_middle = False
                                                if ('S' in ref_cigar) or ('s' in ref_cigar):
                                                    clip_in_middle = True
                                                    if (ref_cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                                                        clip_in_middle = False
                                                    if (ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                        if (ref_pos + aligned_len - 1) == mini_assembly_len_dict[each_ref]:
                                                            clip_in_middle = False

                                                # exclude the ref if clp in the middle is True
                                                if clip_in_middle is False:
                                                    if each_ref not in mini_assembly_to_ctg_reads_dict:
                                                        mini_assembly_to_ctg_reads_dict[each_ref] = {current_read}
                                                    else:
                                                        mini_assembly_to_ctg_reads_dict[each_ref].add(current_read)

                        ########################################### reset values ###########################################

                        current_read = read_id
                        current_read_ref_dict = dict()
                        if cigar != '*':
                            current_read_ref_dict[ref_id] = {ref_pos: cigar}

    # write out mini_assembly_to_16s_reads
    mini_assembly_to_16s_reads_handle = open(mini_assembly_to_16s_reads, 'w')
    for each_mini_assembly in mini_assembly_to_16s_reads_dict:
        mini_assembly_to_16s_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join(mini_assembly_to_16s_reads_dict[each_mini_assembly])))
    mini_assembly_to_16s_reads_handle.close()

    # write out mini_assembly_to_ctg_reads
    mini_assembly_to_ctg_reads_handle = open(mini_assembly_to_ctg_reads, 'w')
    for each_mini_assembly in mini_assembly_to_ctg_reads_dict:
        mini_assembly_to_ctg_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join(mini_assembly_to_ctg_reads_dict[each_mini_assembly])))
    mini_assembly_to_ctg_reads_handle.close()


    ############################################# get_GapFilling_stats #############################################

    if filter_mini_assembly_separately is False:

        get_GapFilling_stats_by_assembly(free_living_16s_ref_file,
                                         free_living_ctg_ref_file,
                                         mini_assembly_to_16s_reads,
                                         mini_assembly_to_ctg_reads,
                                         ctg_level_min_link,
                                         mini_assembly_to_16s_ctg_connector,
                                         gnm_to_ctg_connector,
                                         marker_to_ctg_gnm_Key_connector,
                                         within_gnm_linkage_num_diff,
                                         max_mini_assembly_link_num_diff_between_ctg_16s,
                                         stats_GapFilling_ctg,
                                         stats_GapFilling_file)

        filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, min_iden_16s, min_link_num, min_link_num, within_gnm_linkage_num_diff, stats_GapFilling_file_filtered)

        free_living_16s_to_ctg_linkage_dict_to_use = {}
        for each_ctg_level_link in open(stats_GapFilling_ctg):
            each_ctg_level_link_split = each_ctg_level_link.split('\t')
            marker_id = each_ctg_level_link_split[0]
            ctg_id = each_ctg_level_link_split[1]
            link_num = int(each_ctg_level_link_split[2])
            current_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, ctg_id)
            free_living_16s_to_ctg_linkage_dict_to_use[current_key] = link_num

    else:
        # get_GapFilling_stats_by_assembly_separately()
        pass

    ####################################################################################################################
    ####################################### combine linkages from step 1 and 2  ########################################
    ####################################################################################################################

    ################################################# define file name #################################################

    combined_linkage_file_by_gnm = '%s/%s_linkages_by_genome.txt'    % (working_directory, output_prefix)
    combined_linkage_file_by_ctg = '%s/%s_linkages_by_contig.txt'    % (working_directory, output_prefix)
    linkage_plot_rd1_html        = '%s/%s_linkages_plot_round1.html' % (working_directory, output_prefix)
    linkage_plot_rd2_html        = '%s/%s_linkages_plot_round2.html' % (working_directory, output_prefix)

    report_and_log(('Combining linkages from step 1 and 2'), pwd_log_file, keep_quiet)

    combined_linkage_file_handle     = open(combined_linkage_file_by_gnm, 'w')
    combined_linkage_file_handle.write('MarkerGene\tGenomicSeq\tLinkage\tRound\n')
    for step_1_link in open(link_stats_combined_filtered_s1):
        if not step_1_link.startswith('MarkerGene,GenomicSeq,Number'):
            marker_id = step_1_link.strip().split(',')[0][12:]
            genome_id = step_1_link.strip().split(',')[1][12:]
            link_num  = step_1_link.strip().split(',')[2]
            combined_linkage_file_handle.write('%s\t%s\t%s\tRd1\n' % (marker_id, genome_id, link_num))
    for step_2_link in open(stats_GapFilling_file_filtered):
        if not step_2_link.startswith('MarkerGene,GenomicSeq,Number'):
            marker_id = step_2_link.strip().split(',')[0][12:]
            genome_id = step_2_link.strip().split(',')[1][12:]
            link_num  = step_2_link.strip().split(',')[2]
            combined_linkage_file_handle.write('%s\t%s\t%s\tRd2\n' % (marker_id, genome_id, link_num))
    combined_linkage_file_handle.close()

    #################### summarize linkages at contig level ####################

    combined_linkage_file_ctg_level_handle = open(combined_linkage_file_by_ctg, 'w')
    combined_linkage_file_ctg_level_handle.write('Marker___Genome(total)\tContig\tRd1\tRd2\n')
    for each_linkage in open(combined_linkage_file_by_gnm):
        if not each_linkage.startswith('MarkerGene\tGenomicSeq\tLinkage\tRound'):
            each_linkage_split = each_linkage.strip().split('\t')
            marker_id = each_linkage_split[0]
            mag_id = each_linkage_split[1]
            total_link_num = int(each_linkage_split[2])
            link_step = each_linkage_split[3]

            if link_step == 'Rd1':
                for each_link in marker_to_ctg_linkage_num_dict_min3_passed_ratio_check:
                    link_16s_id = each_link.split(marker_to_ctg_gnm_Key_connector)[0]
                    link_ctg_id = each_link.split(marker_to_ctg_gnm_Key_connector)[1]
                    link_ctg_id_no_gnm = link_ctg_id.split(gnm_to_ctg_connector)[1]
                    link_gnm_id = link_ctg_id.split(gnm_to_ctg_connector)[0]
                    if (link_16s_id == marker_id) and (link_gnm_id == mag_id):
                        current_link_num = marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_link]
                        combined_linkage_file_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t0\n' % (link_16s_id, link_gnm_id, total_link_num, link_ctg_id_no_gnm, current_link_num))

            if link_step == 'Rd2':
                for each_rd2_linkage in free_living_16s_to_ctg_linkage_dict_to_use:
                    rd2_link_16s_id = each_rd2_linkage.split(marker_to_ctg_gnm_Key_connector)[0]
                    rd2_link_ctg_id = each_rd2_linkage.split(marker_to_ctg_gnm_Key_connector)[1]
                    rd2_link_ctg_id_no_gnm = rd2_link_ctg_id.split(gnm_to_ctg_connector)[1]
                    rd2_link_gnm_id = rd2_link_ctg_id.split(gnm_to_ctg_connector)[0]
                    if (rd2_link_16s_id == marker_id) and (rd2_link_gnm_id == mag_id):
                        current_ctg_link_num = free_living_16s_to_ctg_linkage_dict_to_use[each_rd2_linkage]
                        combined_linkage_file_ctg_level_handle.write('%s___%s(%s)\t%s\t0\t%s\n' % (rd2_link_16s_id, rd2_link_gnm_id, total_link_num, rd2_link_ctg_id_no_gnm, current_ctg_link_num))

    combined_linkage_file_ctg_level_handle.close()

    # plot
    report_and_log(('Visualising linkages'), pwd_log_file, keep_quiet)
    sankey_linkages(combined_linkage_file_by_ctg, linkage_plot_rd1_html, linkage_plot_rd2_html)


    ######################################## report assessment under test mode #########################################

    def get_unrecovered_markers(marker_all, marker_recovered):
        unrecovered_markers = []
        for each_marker in marker_all:
            if each_marker not in marker_recovered:
                unrecovered_markers.append(each_marker)
        return sorted(unrecovered_markers)


    if test_mode is True:

        report_and_log(('Test mode on, assessing linkages'), pwd_log_file, keep_quiet)

        marker_id_set = set()
        for marker_seq_record in SeqIO.parse(marker_gene_seqs, 'fasta'):
            marker_id_set.add(marker_seq_record.id)

        # get recovery and accuracy
        recovery_combined, accuracy_combined, recovered_combined = get_accuracy(combined_linkage_file_by_gnm, len(marker_id_set))

        # get unrecovered markers
        unrecovered_markers_paired = get_unrecovered_markers(marker_id_set, recovered_combined)
        unrecovered_markers_paired_str = 'Unrecovered(%s):%s' % (len(unrecovered_markers_paired), ','.join(sorted([i for i in unrecovered_markers_paired])))

        # assessment by genome
        assign_rate, assign_accuracy, right_assign, wrong_assign = get_accuracy_by_genome(combined_linkage_file_by_gnm, mag_folder, mag_file_extension)
        unrecovered_paired_report_str = 'Unrecovered(%s):%s' % (len(wrong_assign), ','.join(sorted([i for i in wrong_assign])))

        # report
        report_and_log(('Prefix\tBy\tRecovery\tAccuracy\tUnrecovered'), pwd_log_file, keep_quiet)
        report_and_log(('%s\tMarker\t%s\t%s\t%s' % (output_prefix, recovery_combined, accuracy_combined, unrecovered_markers_paired_str)), pwd_log_file, keep_quiet)
        report_and_log(('%s\tGenome\t%s\t%s\t%s' % (output_prefix, assign_rate, assign_accuracy, unrecovered_paired_report_str)), pwd_log_file, keep_quiet)


    ################################################### remove tmp files ###################################################

    report_and_log(('Removing temporary files'), pwd_log_file, keep_quiet)

    if keep_temp is False:
        os.remove(input_reads_to_16s_sam)

    # Final report
    report_and_log(('Done!'), pwd_log_file, keep_quiet)


######################################################### main #########################################################

if __name__ == '__main__':

    default_prefix = 'MyRun_%s' % datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
    link_16s_parser = argparse.ArgumentParser(description='Linking MAGs with marker genes', usage=link_Marker_MAG_usage)

    # specify argument group
    link_16s_parser_input_files = link_16s_parser.add_argument_group("input files")
    link_16s_parser_16s         = link_16s_parser.add_argument_group("16S rRNA gene related parameters")
    link_16s_parser_both_rds    = link_16s_parser.add_argument_group("parameters for linking (round 1 and 2)")
    link_16s_parser_rd1         = link_16s_parser.add_argument_group("parameters for linking (round 1)")
    link_16s_parser_rd2         = link_16s_parser.add_argument_group("parameters for linking (round 2)")
    link_16s_parser_preset      = link_16s_parser.add_argument_group("preset parameters, decide automatically if not specified")
    link_16s_parser_dependency  = link_16s_parser.add_argument_group("provide if dependencies are not in your system path")
    link_16s_parser_others      = link_16s_parser.add_argument_group("program settings")
    link_16s_parser_debug       = link_16s_parser.add_argument_group("for debugging, do NOT specify")

    # input files
    link_16s_parser_input_files.add_argument('-p',          required=False, metavar='',             default=default_prefix, help='output prefix, (default: MyRun_SystemTime)')
    link_16s_parser_input_files.add_argument('-r1',         required=True,  metavar='',                                     help='paired reads r1 (fasta format)')
    link_16s_parser_input_files.add_argument('-r2',         required=True,  metavar='',                                     help='paired reads r2 (fasta format)')
    link_16s_parser_input_files.add_argument('-marker',     required=True,  metavar='',                                     help='marker gene sequences')
    link_16s_parser_input_files.add_argument('-mag',        required=True, metavar='',                                      help='metagenome-assembled-genome (MAG) folder')
    link_16s_parser_input_files.add_argument('-x',          required=False, metavar='',             default='fasta',        help='MAG file extension, (default: %(default)s)')

    # depth related
    link_16s_parser_input_files.add_argument('-depth_ratio',required=False, metavar='', type=float, default=0,              help='minimum depth multiple between 16S and  genomic sequences, a value of no higher than 0.2 is recommended, (default: %(default)s)')
    link_16s_parser_input_files.add_argument('-depth_16s',  required=False, metavar='',             default=None,           help='depth info for 16S rRNA genes')
    link_16s_parser_input_files.add_argument('-depth_mag',  required=False, metavar='',             default=None,           help='depth info for MAGs')

    # 16S rRNA gene related parameters
    link_16s_parser_16s.add_argument('-no_polish',          required=False, action="store_true",                            help='skip polishing 16S before linking')
    link_16s_parser_16s.add_argument('-min_iden_16s',       required=False, metavar='', type=float, default=98,             help='minimum similarity for 16S sequences to be assigned to the same genome, (default: %(default)s)')
    link_16s_parser_16s.add_argument('-min_cov_16s',        required=False, metavar='', type=float, default=30,             help='coverage cutoff for calculating pairwise 16S similarity, (default: %(default)s)')
    link_16s_parser_16s.add_argument('-min_aln_16s',        required=False, metavar='', type=int,   default=350,            help='alignment length cutoff for calculating pairwise 16S similarity, (default: %(default)s)')

    # parameters for both rounds linking
    link_16s_parser_both_rds.add_argument('-mismatch',      required=False, metavar='', type=float, default=2,              help='maximum mismatch percentage, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_len_16s', required=False, metavar='', type=int,   default=45,             help='minimum length aligned to 16S, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_len_ctg', required=False, metavar='', type=int,   default=45,             help='minimum length aligned to ctg, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_pct',     required=False, metavar='', type=float, default=35,             help='minimum aligned percentage, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_link',      required=False, metavar='', type=int,   default=8,              help='minimum number of linkages to report, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-link_num_diff', required=False, metavar='', type=float, default=80,             help='within_gnm_linkage_num_diff, (default: %(default)s)')

    # parameters for 1st round linking
    #link_16s_parser_rd1.add_argument('-min_clp_len',        required=False, metavar='', type=int,   default=45,             help='minimum clipping sequence length (bp), (default: %(default)s)')
    #link_16s_parser_rd1.add_argument('-min_clp_M_len',      required=False, metavar='', type=int,   default=35,             help='minimum aligned clipping sequence length (bp), (default: %(default)s)')

    # parameters for 2nd round linking
    #link_16s_parser_rd2.add_argument('-min_overlap_iden',   required=False, metavar='', type=float, default=100,            help='min_overlap_iden, (default: %(default)s)')
    #link_16s_parser_rd2.add_argument('-min_overlap_cov',    required=False, metavar='', type=float, default=50,             help='min_overlap_cov, (default: %(default)s)')
    #link_16s_parser_rd2.add_argument('-min_overlap_len',    required=False, metavar='', type=int,   default=50,             help='min_overlap_len, (default: %(default)s)')
    #link_16s_parser_rd2.add_argument('-min_overlap_num',    required=False, metavar='', type=int,   default=10,             help='minimum number of overlapping reads for a linkages to be reported, (default: %(default)s)')
    link_16s_parser_rd2.add_argument('-assemble_clp',       required=False, action="store_true",                            help='use clipping mapped reads for mini-assembly')
    link_16s_parser_rd2.add_argument('-mira',               required=False, action="store_true",                            help='run Mira, instead of Spades')
    link_16s_parser_rd2.add_argument('-mira_tmp',           required=False, default=None,                                   help='tmp dir for mira')
    link_16s_parser_rd2.add_argument('-link_bias_rd2',      required=False, metavar='', type=float, default=10,             help='max_mini_assembly_link_num_diff_between_ctg_16s, (default: %(default)s)')

    # # preset parameters
    # link_16s_parser_preset.add_argument('-very_sensitive',  required=False, action="store_true",                            help='for greater sensitivity, shortcut for  "min_overlap_iden 99.5 min_overlap_cov 25 min_overlap_len 50 min_overlap_num 3"')
    # link_16s_parser_preset.add_argument('-sensitive',       required=False, action="store_true",                            help='for better sensitivity, shortcut for   "min_overlap_iden 99.5 min_overlap_cov 35 min_overlap_len 50 min_overlap_num 5"')
    # link_16s_parser_preset.add_argument('-specific',        required=False, action="store_true",                            help='for better specificity, shortcut for   "min_overlap_iden 100  min_overlap_cov 55 min_overlap_len 50 min_overlap_num 8"')
    # link_16s_parser_preset.add_argument('-very_specific',   required=False, action="store_true",                            help='for greater specificity, shortcut for  "min_overlap_iden 100  min_overlap_cov 75 min_overlap_len 50 min_overlap_num 10"')

    # program settings
    link_16s_parser_others.add_argument('-bbmap_mem',       required=False, metavar='', type=int,   default=10,             help='bbmap memory allocation (in gigabyte), (default: %(default)s)')
    link_16s_parser_others.add_argument('-t',               required=False, metavar='', type=int,   default=1,              help='number of threads, (default: %(default)s)')
    link_16s_parser_others.add_argument('-tmp',             required=False, action="store_true",                            help='keep temporary files')
    link_16s_parser_others.add_argument('-quiet',           required=False, action="store_true",                            help='not report progress')
    link_16s_parser_others.add_argument('-force',           required=False, action="store_true",                            help='force overwrite existing results')

    # dependency related

    # for debugging
    link_16s_parser_debug.add_argument('-test_mode',        required=False, action="store_true",                            help='only for debugging, do not provide')
    link_16s_parser_debug.add_argument('-rd2_only',         required=False, action="store_true",                            help='run round 2 only')
    link_16s_parser_debug.add_argument('-sorted_sam16s',    required=False, default=None,                                   help='mapping of input reads to 16S')

    args = vars(link_16s_parser.parse_args())
    link_16s(args)


'''

1. the depth of 16S sequences always not lower than the genome they come from
2. with no_ambiguous option, 16S rRNA gene sequences need to be dereplicated. (include dereplication step? with identity and coverage cutoffs?)
3. (doesn't work)!!! to ignore list even without assignment (to handle situations like DM_m4, meanwhile capicable of not assign very diverde 16S (e.g. <98% identity) to the same genome)
4. for clipping mapped reads, the mismatch of clipping part must be 0
5. insert size is important
6. check duplicate sequences in input files
7. a genome depth of 0 will triggle error !!!!!!
8. no header in mag depth and 16s depth files !!!!!!
9. remove short mini-assemblies?

split -l 50000 rd1_extracted_to_gnm_reformatted.sam sam_split_ --additional-suffix=.sam

split -n 36 GI_0522_input_reads_to_16S_reformatted_sorted.sam GI_0522_input_reads_to_16S_reformatted_sorted_ --additional-suffix=.sam


'''
