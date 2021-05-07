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
from Bio import AlignIO
from Bio.Seq import Seq
import multiprocessing as mp
from datetime import datetime
import plotly.graph_objects as go
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from distutils.spawn import find_executable

link_Marker_MAG_usage = '''
=================================== MarkerMAG example commands ===================================

# example commands
MarkerMAG link -p Test -r1 R1.fasta -r2 R2.fasta -m 16S_seqs.fa -mag MAG_folder -x fa -t 6

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


def get_ctg_mean_depth_by_samtools_coverage(index_ref, ref_seq, reads_r1, reads_r2, reads_unpaired, num_threads):

    ref_seq_file_path, ref_seq_file_basename, ref_seq_file_extension = sep_path_basename_ext(ref_seq)

    sam_file        = '%s/%s.sam'          % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_sorted = '%s/%s_sorted.sam'   % (ref_seq_file_path, ref_seq_file_basename)
    coverage_file   = '%s/%s_cov.txt'      % (ref_seq_file_path, ref_seq_file_basename)

    # build reference index
    cmd_bowtie2_build   = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, ref_seq, ref_seq_file_path, ref_seq_file_basename)
    if index_ref is True:
        os.system(cmd_bowtie2_build)

    # mapping
    if reads_unpaired == '':
        cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
    else:
        cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -1 %s -2 %s -U %s -S %s -p %s -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, reads_unpaired, sam_file, num_threads)
    os.system(cmd_bowtie2_mapping)

    # sort mapping
    cmd_samtools_sort = 'samtools sort %s -o %s' % (sam_file, sam_file_sorted)
    os.system(cmd_samtools_sort)

    # get mean depth
    cmd_samtools_coverage = 'samtools coverage %s -o %s' % (sam_file_sorted, coverage_file)
    os.system(cmd_samtools_coverage)

    # remove sam files
    os.system('rm %s' % sam_file)
    os.system('rm %s' % sam_file_sorted)

    # store mean depth into dict
    mean_depth_dict_ctg = {}
    ctg_len_dict = {}
    for each_ctg_depth in open(coverage_file):
        if not each_ctg_depth.startswith('#'):
            ctg_depth_split = each_ctg_depth.strip().split('\t')
            ctg_id = ctg_depth_split[0]
            ctg_len = int(ctg_depth_split[2])
            ctg_depth = float(ctg_depth_split[6])
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
        if not each_match.startswith('MarkerGene,GenomicSeq,Number'):
            match_split = each_match.strip().split(',')
            linkage_num = int(match_split[2])
            MarkerGene_genome = match_split[0][12:][:2]
            GenomicSeq_genome = match_split[1][12:]

            linkage_num_total += linkage_num
            if MarkerGene_genome == GenomicSeq_genome:
                linkage_num_correct += linkage_num
                recovered_markers.add(match_split[0][12:])

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
        if not each_match.startswith('MarkerGene,GenomicSeq,Number'):
            match_split = each_match.strip().split(',')
            MarkerGene_genome = match_split[0][12:][:2]
            GenomicSeq_genome = match_split[1][12:]

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
        if not each_linkage.startswith('Marker___Genome(total)	Contig	Paired	Clipping	Overlapped	Step'):
            each_linkage_split = each_linkage.strip().split('\t')

            marker_id = each_linkage_split[0].split('___')[0]
            gnm_id = each_linkage_split[0].split('___')[1].split('(')[0]
            ctg_id = each_linkage_split[1]
            total_link_num = int(each_linkage_split[2]) + int(each_linkage_split[3]) + int(each_linkage_split[4])
            marker_to_ctg_key = '%s%s%s' % (marker_id, dict_for_sankey_key_connector, ctg_id)
            ctg_to_gnm_key = '%s%s%s' % (ctg_id, dict_for_sankey_key_connector, gnm_id)

            if each_linkage_split[5] == 'S1':
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

            if each_linkage_split[5] == 'S2':
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


def keep_best_matches_in_sam_drop_short_M(sam_in, min_M_len, sam_out):

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
                    if (mismatch_pct <= (read_min_mismatch_dict[read_id] * 1.5)) and (aligned_len >= min_M_len):
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


def mview_linkage(seq_file, plot_title, mafft_exe, mview_exe):

    # define file name
    seq_path, seq_basename, seq_ext = sep_path_basename_ext(seq_file)
    msa_file                    = '%s/%s.aln'               % (seq_path, seq_basename)
    msa_file_mview              = '%s/%s_MView.aln'         % (seq_path, seq_basename)
    msa_file_mviewd_html_tmp    = '%s/%s_MView_tmp.html'    % (seq_path, seq_basename)
    msa_file_mviewd_html        = '%s/%s_MView.html'        % (seq_path, seq_basename)

    gap_char                = '-'
    break_line_char         = '='

    # align sequences
    mafft_cmd = '%s --quiet --retree 1 %s > %s' % (mafft_exe, seq_file, msa_file)
    os.system(mafft_cmd)

    mapped_reads_dict = {}
    ref_id = ''
    ref_seq = ''
    current_line = 0
    for each_seq in AlignIO.read(msa_file, "fasta"):
        seq_id = each_seq.id
        seq_seq = str(each_seq.seq).upper()
        if current_line == 0:
            ref_id = seq_id
            ref_seq = seq_seq
        else:
            seq_id_base = '.'.join(seq_id.split('.')[:-1])
            seq_id_strand = seq_id.split('.')[-1]
            if seq_id_base not in mapped_reads_dict:
                if seq_id_strand == '1':
                    mapped_reads_dict[seq_id_base] = [seq_seq, '']
                if seq_id_strand == '2':
                    mapped_reads_dict[seq_id_base] = ['', seq_seq]
            else:
                if seq_id_strand == '1':
                    mapped_reads_dict[seq_id_base][0] = seq_seq
                if seq_id_strand == '2':
                    mapped_reads_dict[seq_id_base][1] = seq_seq

        current_line += 1

    # create an empty list to hold all sequences in a msa
    seq_record_list = []

    # add ref_seq to seq_record_list
    align_record_mview = MultipleSeqAlignment([])
    ref_seq_split_by_n = ref_seq.split('N')
    ref_seq_split_by_n_updated = []
    for segment in ref_seq_split_by_n:
        if ('-' in segment) and (segment == '-' * len(segment)):
            segment = 'N' * len(segment)
        ref_seq_split_by_n_updated.append(segment)
    ref_seq_updated = 'N'.join(ref_seq_split_by_n_updated)
    seq_record_list.append(SeqRecord(Seq(ref_seq_updated), id=ref_id, description='Reference'))

    # add break line
    seq_record_list.append(SeqRecord(Seq(break_line_char*len(ref_seq_updated)), id='#', description=''))

    # add paired reads to seq_record_list
    singleton_dict = {}
    overlapping_reads_dict = {}
    for each_read_base in mapped_reads_dict:

        r1_seq = mapped_reads_dict[each_read_base][0]
        r2_seq = mapped_reads_dict[each_read_base][1]

        if (r1_seq == '') or (r2_seq == ''):
            if not ((r1_seq == '') and (r2_seq == '')):
                singleton_dict[each_read_base] = mapped_reads_dict[each_read_base]
        else:
            overlapping_bps = 0
            for bp1, bp2 in zip(r1_seq, r2_seq):
                if (bp1 != '-') and (bp2 != '-'):
                    overlapping_bps += 1

            if overlapping_bps > 0:
                overlapping_reads_dict[each_read_base] = mapped_reads_dict[each_read_base]
            else:
                merge_r1_r2 = ''
                for bp1, bp2 in zip(r1_seq, r2_seq):
                    if (bp1 == '-') and (bp2 == '-'):
                        merge_r1_r2 += '-'
                    if (bp1 != '-') and (bp2 == '-'):
                        merge_r1_r2 += bp1
                    if (bp1 == '-') and (bp2 != '-'):
                        merge_r1_r2 += bp2
                seq_record_list.append(SeqRecord(Seq(merge_r1_r2), id=each_read_base, description='In_Pair'))

    # add break line
    seq_record_list.append(SeqRecord(Seq(break_line_char*len(ref_seq_updated)), id='#', description=''))

    # add overlapped reads to seq_record_list
    for each_overlapping_reads in overlapping_reads_dict:
        r1_id  = '%s.1' % each_overlapping_reads
        r2_id  = '%s.2' % each_overlapping_reads
        r1_seq = overlapping_reads_dict[each_overlapping_reads][0]
        r2_seq = overlapping_reads_dict[each_overlapping_reads][1]
        seq_record_list.append(SeqRecord(Seq(r1_seq), id=r1_id, description='Overlapped'))
        seq_record_list.append(SeqRecord(Seq(r2_seq), id=r2_id, description='Overlapped'))

    # add break line
    seq_record_list.append(SeqRecord(Seq(break_line_char*len(ref_seq_updated)), id='#', description=''))

    # add singleton to seq_record_list
    for each_singleton in singleton_dict:
        r1_id  = '%s.1' % each_singleton
        r2_id  = '%s.2' % each_singleton
        r1_seq = singleton_dict[each_singleton][0]
        r2_seq = singleton_dict[each_singleton][1]
        if (r1_seq != '') and (r2_seq == ''):
            seq_record_list.append(SeqRecord(Seq(r1_seq), id=r1_id, description='Unpaired'))
        if (r1_seq == '') and (r2_seq != ''):
            seq_record_list.append(SeqRecord(Seq(r2_seq), id=r2_id, description='Unpaired'))

    # add break line
    seq_record_list.append(SeqRecord(Seq(break_line_char*len(ref_seq_updated)), id='#', description=''))

    # get updated msa record
    align_record_mview = MultipleSeqAlignment(seq_record_list)

    # write out updated msa
    msa_file_mview_handle = open(msa_file_mview, 'w')
    AlignIO.write(align_record_mview, msa_file_mview_handle, 'fasta')
    msa_file_mview_handle.close()

    # visualize update msa
    # coloring: any,identity,mismatch,consensus,group
    mview_parameter_str = '-in fasta -moltype dna -colormap CLUSTAL_NUC -coloring any -css on -html head -ruler off -label0 -label4 -label5 -gap "%s"' % gap_char
    mview_cmd = '%s %s -title %s %s > %s' % (mview_exe, mview_parameter_str, plot_title, msa_file_mview, msa_file_mviewd_html_tmp)
    os.system(mview_cmd)

    msa_file_mviewd_html_handle = open(msa_file_mviewd_html, 'w')
    for each_line in open(msa_file_mviewd_html_tmp):
        if not (('Reference sequence' in each_line) or ('Colored by' in each_line)):
            msa_file_mviewd_html_handle.write(each_line)
    msa_file_mviewd_html_handle.close()


def get_unmapped_mates_seq(sam_file, input_r1_fasta, input_r2_fasta, extracted_seq_file):
    sam_path, sam_basename, sam_ext = sep_path_basename_ext(sam_file)
    reads_to_extract_r1_txt = '%s/%s_unmapped_mates_R1.txt' % (sam_path, sam_basename)
    reads_to_extract_r2_txt = '%s/%s_unmapped_mates_R2.txt' % (sam_path, sam_basename)
    reads_to_extract_r1_fa = '%s/%s_unmapped_mates_R1.fa' % (sam_path, sam_basename)
    reads_to_extract_r2_fa = '%s/%s_unmapped_mates_R2.fa' % (sam_path, sam_basename)

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


def link_16s(args):

    ###################################################### file in/out #####################################################

    # file in
    output_prefix                       = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    reads_file_16s                      = args['r16s']
    #genomic_assemblies                  = args['g']
    mag_folder                          = args['mag']
    mag_file_extension                  = args['x']
    marker_gene_seqs                    = args['marker']
    min_16s_gnm_multiple                = args['depth']
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

    # by assembly
    round_2_spades                      = args['spades']
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
    mafft_exe                           = 'mafft'
    mview_exe                           = 'mview'

    marker_to_ctg_gnm_Key_connector                 = '___M___'
    gnm_to_ctg_connector                            = '___C___'
    mini_assembly_to_16s_ctg_connector              = '___Mini___'
    end_seq_len                                     = 500
    ctg_level_min_link                              = 3
    end_ctg_len_for_mafft                           = 1000
    vis_linkages                                    = False


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

    ############################################# create working directory #############################################

    # create working directory
    working_directory = '%s_MarkerMAG_wd' % output_prefix

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
    paired_reads_num, read_len_median, read_len_max = get_read_num_and_length(reads_file_r1_fasta, working_directory, seqtk_exe)
    estimated_total_read_len_gbp = (read_len_median*paired_reads_num*2)/(1024*1024*1024)
    estimated_total_read_len_gbp = float("{0:.1f}".format(estimated_total_read_len_gbp))

    parameter_list = []
    parameter_list.append('mismatch:%s'             % mismatch_cutoff)
    parameter_list.append('min_M_len_16s:%s'        % min_M_len_16s)
    parameter_list.append('min_M_len_ctg:%s'        % min_M_len_ctg)
    parameter_list.append('min_M_pct:%s'            % min_M_pct)
    parameter_list.append('min_link_num:%s'         % min_link_num)
    parameter_list.append('ctg_level_min_link:%s'   % (ctg_level_min_link))
    parameter_list.append('end_seq_len:%s'          % (end_seq_len))
    parameter_list.append('link_bias_rd2:%s'        % (max_mini_assembly_link_num_diff_between_ctg_16s))
    parameter_str = ';'.join(sorted(parameter_list))
    print(parameter_str)


    ######################## check genomic sequence type and prepare files for making blast db #########################

    combined_input_gnms = ''
    # genomic_seq_type    = ''  # ctg or mag
    renamed_mag_folder  = ''

    # # check the type of input genomic sequences
    # if (genomic_assemblies is not None) and (mag_folder is None):
    #     genomic_seq_type = 'ctg'
    #     metagenomic_assemblies_file_path, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension = sep_path_basename_ext(genomic_assemblies)
    #     blast_db_dir = '%s/%s_%s_db' % (step_1_wd, output_prefix, metagenomic_assemblies_file_basename)
    #     combined_input_gnms     = '%s/%s%s'     % (blast_db_dir, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension)
    #
    #     os.mkdir(blast_db_dir)
    #     os.system('cp %s %s/' % (genomic_assemblies, blast_db_dir))
    # elif (genomic_assemblies is None) and (mag_folder is not None):
    if mag_folder is not None:

        # get input mag file list
        mag_file_re = '%s/*%s' % (mag_folder, mag_file_extension)
        mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
        if len(mag_file_list) == 0:
            print('No MAG detected, program exited!')
            exit()

        # genomic_seq_type    = 'mag'
        mag_folder_name     = mag_folder.split('/')[-1]
        blast_db_dir        = '%s/%s_db'            % (step_1_wd, mag_folder_name)
        renamed_mag_folder  = '%s/%s_db/%s_renamed' % (step_1_wd, mag_folder_name, mag_folder_name)
        os.mkdir(blast_db_dir)
        os.mkdir(renamed_mag_folder)

        # add mag id to its sequences
        for mag_in in mag_file_list:
            pwd_mag_in      = '%s/%s' % (mag_folder, mag_in)
            pwd_mag_renamed = '%s/%s' % (renamed_mag_folder, mag_in)
            mag_basename    = '.'.join(mag_in.split('.')[:-1])
            rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename, gnm_to_ctg_connector)

        # combine renamed MAGs
        combined_input_gnms = '%s/%s_combined.fa' % (blast_db_dir, mag_folder_name)
        os.system('cat %s/*%s > %s' % (renamed_mag_folder, mag_file_extension, combined_input_gnms))
    else:
        print('MAGs (-mag) not provided')
        exit()


    ########################################### define folder and file name ############################################

    marker_gene_seqs_file_path, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension = sep_path_basename_ext(marker_gene_seqs)

    pwd_log_file                                = '%s/%s.log'                                              % (working_directory, output_prefix)
    bowtie_index_dir                            = '%s/%s_index'                                            % (step_1_wd, marker_gene_seqs_file_basename)
    input_reads_to_16s_sam_bowtie               = '%s/%s_input_reads_to_16S_bowtie.sam'                    % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_bowtie_log           = '%s/%s_input_reads_to_16S_bowtie.log'                    % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_bowtie_single_clp    = '%s/%s_input_reads_to_16S_bowtie_single_clp.sam'         % (step_1_wd, output_prefix)
    input_reads_to_16s_sam                      = '%s/%s_input_reads_to_16S.sam'                           % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_reformat_log         = '%s/%s_input_reads_to_16S_reformat.log'                  % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_best_match           = '%s/%s_input_reads_to_16S_best_match.sam'                % (step_1_wd, output_prefix)
    input_reads_to_16s_sam_bbmap_stderr         = '%s/%s_input_reads_to_16S_bbmap_stderr.txt'              % (step_1_wd, output_prefix)
    sam_best_match_unmapped_mates_seq_file      = '%s/%s_input_reads_to_16S_best_match_unmapped_mates.fa'  % (step_1_wd, output_prefix)
    unmapped_mates_seq_file                     = '%s/%s_unmapped_mates.fa'                                % (step_1_wd, output_prefix)
    clipping_parts_seq_file                     = '%s/%s_clipping_parts.fa'                                % (step_1_wd, output_prefix)
    unmapped_to_gnm_sam_bowtie                  = '%s/%s_unmapped_mates_bowtie.sam'                        % (step_1_wd, output_prefix)
    unmapped_to_gnm_sam_bowtie_log              = '%s/%s_unmapped_mates_bowtie.log'                        % (step_1_wd, output_prefix)
    unmapped_to_gnm_sam_bowtie_single_clp       = '%s/%s_unmapped_mates_bowtie_single_clp.sam'             % (step_1_wd, output_prefix)
    unmapped_to_gnm_sam                         = '%s/%s_unmapped_mates.sam'                               % (step_1_wd, output_prefix)
    unmapped_to_gnm_sam_reformat_log            = '%s/%s_unmapped_mates_reformat.log'                      % (step_1_wd, output_prefix)
    unmapped_to_gnm_sam_best_match              = '%s/%s_unmapped_mates_best_match.sam'                    % (step_1_wd, output_prefix)
    clipping_to_gnm_sam_bowtie                  = '%s/%s_clipping_parts_bowtie.sam'                        % (step_1_wd, output_prefix)
    clipping_to_gnm_sam_bowtie_log              = '%s/%s_clipping_parts_bowtie.log'                        % (step_1_wd, output_prefix)
    clipping_to_gnm_sam_bowtie_single_clp       = '%s/%s_clipping_parts_bowtie_single_clp.sam'             % (step_1_wd, output_prefix)
    clipping_to_gnm_sam                         = '%s/%s_clipping_parts.sam'                               % (step_1_wd, output_prefix)
    clipping_to_gnm_sam_reformat_log            = '%s/%s_clipping_parts_reformat.log'                      % (step_1_wd, output_prefix)
    clipping_to_gnm_sam_best_match              = '%s/%s_clipping_parts_best_match.sam'                    % (step_1_wd, output_prefix)
    unmapped_to_gnm_bbmap_stderr                = '%s/%s_unmapped_mates_bbmap_stderr.txt'                  % (step_1_wd, output_prefix)
    clipping_to_gnm_bbmap_stderr                = '%s/%s_clipping_parts_bbmap_stderr.txt'                  % (step_1_wd, output_prefix)
    blast_results_all_vs_all_16s                = '%s/%s_16S_all_vs_all_blastn.tab'                        % (step_1_wd, output_prefix)
    pairwise_marker_similarity                  = '%s/%s_pairwise_marker_similarity.txt'                   % (step_1_wd, output_prefix)
    depth_file_ctg                              = '%s/%s_mean_depth_ctg.txt'                               % (step_1_wd, output_prefix)
    depth_file_gnm                              = '%s/%s_mean_depth_gnm.txt'                               % (step_1_wd, output_prefix)
    depth_file_16s                              = '%s/%s_mean_depth_16s.txt'                               % (step_1_wd, output_prefix)
    link_stats_combined                         = '%s/%s_stats_combined.txt'                               % (step_1_wd, output_prefix)
    link_stats_combined_filtered_s1             = '%s/%s_stats_combined_filtered.txt'                      % (step_1_wd, output_prefix)
    linking_reads_pair                          = '%s/%s_linking_read_rd1_paired.txt'                      % (step_1_wd, output_prefix)
    linking_reads_clip                          = '%s/%s_linking_read_rd1_clipping.txt'                    % (step_1_wd, output_prefix)
    linked_ends_rd1_txt                         = '%s/%s_linked_ends_rd1.txt'                              % (step_1_wd, output_prefix)
    linked_contigs_txt                          = '%s/%s_linked_contigs.txt'                               % (step_1_wd, output_prefix)
    linked_contigs_fasta                        = '%s/%s_linked_contigs.fasta'                             % (step_1_wd, output_prefix)
    mafft_seq_folder                            = '%s/%s_mafft_seq_folder'                                 % (step_1_wd, output_prefix)

    blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % num_threads
    bbmap_parameter  = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % (num_threads, bbmap_memory)

    ################################################# step 2 #################################################

    marker_gene_seqs_1st_round_unlinked                             = '%s/round_1_unlinked_16s.fa'                              % step_2_wd
    combined_1st_round_unlinked_mags                                = '%s/round_1_unlinked_gnm.fa'                              % step_2_wd
    combined_1st_round_unlinked_mag_end_seq                         = '%s/round_1_unlinked_gnm_end_%sbp.fa'                     % (step_2_wd, end_seq_len)
    rd1_unlinked_mag_end_seq_no_ext                                 = '%s/round_1_unlinked_gnm_end_%sbp'                        % (step_2_wd, end_seq_len)
    combined_1st_round_unlinked_mags_sam                            = '%s/round_1_unlinked_gnm.sam'                             % step_2_wd
    rd1_unlinked_mags_sam_bowtie_log                                = '%s/round_1_unlinked_gnm_bowtie.log'                      % step_2_wd
    rd1_unlinked_mags_sam_bowtie                                    = '%s/round_1_unlinked_gnm_bowtie.sam'                      % step_2_wd
    rd1_unlinked_mags_sam_bowtie_single_clp                         = '%s/round_1_unlinked_gnm_bowtie_single_clp.sam'           % step_2_wd
    rd1_unlinked_mags_sam_bowtie_single_clp_reformat                = '%s/round_1_unlinked_gnm_bowtie_single_clp_reformat.sam'  % step_2_wd
    rd1_unlinked_mags_sam_bowtie_single_clp_reformat_log            = '%s/round_1_unlinked_gnm_bowtie_single_clp_reformat.log'  % step_2_wd
    combined_1st_round_unlinked_mags_sam_best_match                 = '%s/round_1_unlinked_gnm_best_match.sam'                  % step_2_wd
    combined_1st_round_unlinked_mags_sam_best_match_unmapped_mates  = '%s/round_1_unlinked_gnm_best_match_unmapped_mates.fa'    % step_2_wd

    combined_1st_round_unlinked_ctgs                = '%s/round_1_unlinked_ctg.fa'                      % step_2_wd
    stats_GapFilling_file                           = '%s/stats_GapFilling_gnm.txt'                     % step_2_wd
    stats_GapFilling_file_filtered                  = '%s/stats_GapFilling_gnm_filtered.txt'            % step_2_wd
    free_living_16s_R1                              = '%s/round2_free_living_16s_R1.fa'                 % step_2_wd
    free_living_16s_R2                              = '%s/round2_free_living_16s_R2.fa'                 % step_2_wd
    free_living_16s_UP                              = '%s/round2_free_living_16s_UP.fa'                 % step_2_wd
    free_living_ctg_R1                              = '%s/round2_free_living_ctg_R1.fa'                 % step_2_wd
    free_living_ctg_R2                              = '%s/round2_free_living_ctg_R2.fa'                 % step_2_wd
    free_living_ctg_UP                              = '%s/round2_free_living_ctg_UP.fa'                 % step_2_wd
    free_living_16s                                 = '%s/round2_free_living_16s.fa'                    % step_2_wd
    free_living_ctg                                 = '%s/round2_free_living_ctg.fa'                    % step_2_wd
    free_living_16s_ref_file                        = '%s/round2_free_living_16s_refs.txt'              % step_2_wd
    free_living_ctg_ref_file                        = '%s/round2_free_living_ctg_refs.txt'              % step_2_wd
    free_living_blast_result                        = '%s/free_living_reads_blastn.tab'                 % step_2_wd
    free_living_all                                 = '%s/round2_free_living_all.fa'                    % step_2_wd
    free_living_all_fq_r1                           = '%s/round2_free_living_all_R1.fastq'              % step_2_wd
    free_living_all_fq_r2                           = '%s/round2_free_living_all_R2.fastq'              % step_2_wd
    free_living_all_fq                              = '%s/round2_free_living_all.fastq'                 % step_2_wd
    free_living_R1                                  = '%s/round2_free_living_R1.fa'                     % step_2_wd
    free_living_R2                                  = '%s/round2_free_living_R2.fa'                     % step_2_wd
    free_living_UP                                  = '%s/round2_free_living_UP.fa'                     % step_2_wd
    spades_wd                                       = '%s/mini_assembly_SPAdes_wd'                      % step_2_wd
    spades_log                                      = '%s/SPAdes_stdout.txt'                            % step_2_wd
    mira_manifest                                   = '%s/mira_manifest.txt'                            % step_2_wd
    mira_stdout                                     = '%s/mira_stdout.txt'                              % step_2_wd
    sam_file_mini_assembly_16s_P                    = '%s/scaffolds_16s_P.sam'                          % step_2_wd
    sam_file_mini_assembly_16s_UP                   = '%s/scaffolds_16s_UP.sam'                         % step_2_wd
    sam_file_mini_assembly_ctg_P                    = '%s/scaffolds_ctg_P.sam'                          % step_2_wd
    sam_file_mini_assembly_ctg_UP                   = '%s/scaffolds_ctg_UP.sam'                         % step_2_wd
    sam_file_mini_assembly_16s_P_stderr             = '%s/scaffolds_16s_P_bbmap_stderr.txt'             % step_2_wd
    sam_file_mini_assembly_16s_UP_stderr            = '%s/scaffolds_16s_UP_bbmap_stderr.txt'            % step_2_wd
    sam_file_mini_assembly_ctg_P_stderr             = '%s/scaffolds_ctg_P_bbmap_stderr.txt'             % step_2_wd
    sam_file_mini_assembly_ctg_UP_stderr            = '%s/scaffolds_ctg_UP_bbmap_stderr.txt'            % step_2_wd
    sam_file_mini_assembly_16s_P_best_match         = '%s/scaffolds_16s_P_best_match.sam'               % step_2_wd
    sam_file_mini_assembly_16s_UP_best_match        = '%s/scaffolds_16s_UP_best_match.sam'              % step_2_wd
    sam_file_mini_assembly_ctg_P_best_match         = '%s/scaffolds_ctg_P_best_match.sam'               % step_2_wd
    sam_file_mini_assembly_ctg_UP_best_match        = '%s/scaffolds_ctg_UP_best_match.sam'              % step_2_wd
    sam_file_mini_assembly_P                        = '%s/scaffolds_P.sam'                              % step_2_wd
    sam_file_mini_assembly_UP                       = '%s/scaffolds_UP.sam'                             % step_2_wd
    sam_file_mini_assembly_stderr_P                 = '%s/scaffolds_bbmap_stderr_P.txt'                 % step_2_wd
    sam_file_mini_assembly_stderr_UP                = '%s/scaffolds_bbmap_stderr_UP.txt'                % step_2_wd
    sam_file_mini_assembly_combined                 = '%s/scaffolds_combined.sam'                       % step_2_wd
    sam_file_mini_assembly_16s_bowtie_log           = '%s/scaffolds_16s_bowtie.log'                     % step_2_wd
    sam_file_mini_assembly_16s                      = '%s/scaffolds_16s_bowtie.sam'                     % step_2_wd
    sam_file_mini_assembly_16s_single_clp           = '%s/scaffolds_16s_bowtie_single_clp.sam'          % step_2_wd
    sam_file_mini_assembly_16s_single_clp_fmt       = '%s/scaffolds_16s_bowtie_single_clp_reformat.sam' % step_2_wd
    sam_file_mini_assembly_16s_single_clp_fmt_log   = '%s/scaffolds_16s_bowtie_single_clp_reformat.log' % step_2_wd
    sam_file_mini_assembly_ctg_bowtie_log           = '%s/scaffolds_ctg_bowtie.log'                     % step_2_wd
    sam_file_mini_assembly_ctg                      = '%s/scaffolds_ctg_bowtie.sam'                     % step_2_wd
    sam_file_mini_assembly_ctg_single_clp           = '%s/scaffolds_ctg_bowtie_single_clp.sam'          % step_2_wd
    sam_file_mini_assembly_ctg_single_clp_fmt       = '%s/scaffolds_ctg_bowtie_single_clp_reformat.sam' % step_2_wd
    sam_file_mini_assembly_ctg_single_clp_fmt_log   = '%s/scaffolds_ctg_bowtie_single_clp_reformat.log' % step_2_wd
    sam_file_mini_assembly_16s_stderr               = '%s/scaffolds_16s_stderr.txt'                     % step_2_wd
    sam_file_mini_assembly_ctg_stderr               = '%s/scaffolds_ctg_stderr.txt'                     % step_2_wd
    sam_file_mini_assembly_16s_best_match           = '%s/scaffolds_16s_best_match.sam'                 % step_2_wd
    sam_file_mini_assembly_ctg_best_match           = '%s/scaffolds_ctg_best_match.sam'                 % step_2_wd
    sam_file_mini_assembly_both_best_match          = '%s/scaffolds_both_best_match.sam'                % step_2_wd
    sam_file_mini_assembly_both_best_match_fa       = '%s/scaffolds_both_best_match.fa'                 % step_2_wd
    mini_assembly_to_16s_reads                      = '%s/mini_assembly_to_16s_reads.txt'               % step_2_wd
    mini_assembly_to_ctg_reads                      = '%s/mini_assembly_to_ctg_reads.txt'               % step_2_wd
    stats_GapFilling_ctg                            = '%s/stats_GapFilling_ctg.txt'                     % step_2_wd


    ################################################# combine linkages from two steps #################################################

    combined_linkage_file_tmp                       = '%s/combined_linkages_tmp.txt'                    % step_2_wd
    combined_linkage_file                           = '%s/%s_identified_linkages_genome_level.txt'      % (working_directory, output_prefix)
    combined_linkage_file_ctg_level                 = '%s/%s_identified_linkages_contig_level.txt'      % (working_directory, output_prefix)
    linkage_plot_rd1_html                           = '%s/%s_identified_linkages_round1.html'           % (working_directory, output_prefix)
    linkage_plot_rd2_html                           = '%s/%s_identified_linkages_round2.html'           % (working_directory, output_prefix)

    parameter_str = ';'.join(sorted(parameter_list))
    report_and_log(('Specified parameters: %s'          % parameter_str), pwd_log_file, keep_quiet)
    report_and_log(('Number of paired reads: %s'        % paired_reads_num), pwd_log_file, keep_quiet)
    report_and_log(('Read length median: %s'            % read_len_median), pwd_log_file, keep_quiet)
    report_and_log(('Read length max: %s'               % read_len_max), pwd_log_file, keep_quiet)
    report_and_log(('Estimated total length: %s Gbp'    % estimated_total_read_len_gbp), pwd_log_file, keep_quiet)


    #################################### calculate mean depth for genome/assemblies ####################################

    # mean_depth_dict_ctg = {}
    mean_depth_dict_gnm = {}
    if min_16s_gnm_multiple > 0:
        # if genomic_seq_type == 'ctg':
        #     report_and_log(('Round 1: calculating depth for %s' % genomic_assemblies), pwd_log_file, keep_quiet)
        # if genomic_seq_type == 'mag':
        report_and_log(('Round 1: calculating depth for genomes in %s' % mag_folder), pwd_log_file, keep_quiet)

        # get mean depth for contig
        mean_depth_dict_ctg, ctg_len_dict = get_ctg_mean_depth_by_samtools_coverage(True, combined_input_gnms, reads_file_r1, reads_file_r2, '', num_threads)

        # write out ctg depth
        depth_file_ctg_handle = open(depth_file_ctg, 'w')
        for ctg in mean_depth_dict_ctg:
            depth_file_ctg_handle.write('%s\t%s\n' % (ctg, mean_depth_dict_ctg[ctg]))
        depth_file_ctg_handle.close()

        # get mean_depth_dict_gnm
        # if genomic_seq_type == 'mag':
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


    ###################################### calculate mean depth for 16S sequences ######################################

    os.mkdir(bowtie_index_dir)
    os.system('cp %s %s/' % (marker_gene_seqs, bowtie_index_dir))

    mean_depth_dict_16s = {}
    if min_16s_gnm_multiple > 0:

        report_and_log(('Round 1: calculating depth for %s' % marker_gene_seqs), pwd_log_file, keep_quiet)

        marker_gene_seqs_path, marker_gene_seqs_basename, marker_gene_seqs_extension = sep_path_basename_ext(marker_gene_seqs)
        pwd_16s = '%s/%s%s' % (bowtie_index_dir, marker_gene_seqs_basename, marker_gene_seqs_extension)

        sortmerna_exe = 'sortmerna'
        if reads_file_16s is None:
            print('Run SortMeRNA, be patient')
            sortmerna_cmd = '%s --ref /srv/scratch/z5039045/DB/Matam/SILVA_128_SSURef_NR95.clustered.fasta,/srv/scratch/z5039045/DB/Matam/SILVA_128_SSURef_NR95.clustered --reads MBARC26_R1_R2.fasta --aligned MBARC26_Matam16S_wd/MBARC26 --fastx --sam --blast "1" --log --best 10 --min_lis 10 -e 1.00e-05 -a 16 -v > MBARC26_Matam16S_wd/MBARC26.SortMeRNA_stdout.txt' % (sortmerna_exe)
            # os.system(sortmerna_cmd)
            print(sortmerna_cmd)

        # separate paired and singleton reads
        reads_file_16s_path, reads_file_16s_basename, reads_file_16s_extension = sep_path_basename_ext(reads_file_16s)
        reads_file_16s_r1           = '%s/%s_r1.fasta'          % (working_directory, reads_file_16s_basename)
        reads_file_16s_r2           = '%s/%s_r2.fasta'          % (working_directory, reads_file_16s_basename)
        reads_file_16s_singleton    = '%s/%s_singleton.fasta'   % (working_directory, reads_file_16s_basename)
        sep_paired_and_singleton_reads(reads_file_16s, reads_file_16s_r1, reads_file_16s_r2, reads_file_16s_singleton)

        # get mean depth for 16S sequences
        mean_depth_dict_16s, s16_len_dict = get_ctg_mean_depth_by_samtools_coverage(True, pwd_16s, reads_file_16s_r1, reads_file_16s_r2, reads_file_16s_singleton, num_threads)

        # write out 16s depth
        depth_file_16s_handle = open(depth_file_16s, 'w')
        for s16 in mean_depth_dict_16s:
            depth_file_16s_handle.write('%s\t%s\n' % (s16, mean_depth_dict_16s[s16]))
        depth_file_16s_handle.close()


    ####################################################################################################################
    ############################################### first round linking ################################################
    ####################################################################################################################

    ######################################## map reads to marker gene sequences ########################################

    marker_gene_seqs_in_wd        = '%s/%s%s' % (bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension)
    marker_gene_seqs_in_wd_no_ext = '%s/%s'   % (bowtie_index_dir, marker_gene_seqs_file_basename)

    report_and_log(('Round 1: Mapping input reads to marker genes'), pwd_log_file, keep_quiet)
    #bbmap_index_and_mapping_cmd = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, marker_gene_seqs_in_wd, reads_file_r1, reads_file_r2, input_reads_to_16s_sam, bbmap_parameter, input_reads_to_16s_sam_bbmap_stderr)
    #report_and_log((bbmap_index_and_mapping_cmd), pwd_log_file, True)
    #os.system(bbmap_index_and_mapping_cmd)

    bowtie_build_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, marker_gene_seqs_in_wd, marker_gene_seqs_in_wd_no_ext)
    os.system(bowtie_build_cmd)

    bowtie_read_to_16s_cmd = 'bowtie2 -x %s -1 %s -2 %s -S %s -p %s -f --local --all --no-unal --very-fast-local 2> %s' % (marker_gene_seqs_in_wd_no_ext, reads_file_r1_fasta, reads_file_r2_fasta, input_reads_to_16s_sam_bowtie, num_threads, input_reads_to_16s_sam_bowtie_log)
    os.system(bowtie_read_to_16s_cmd)

    report_and_log(('Round 1: removing both ends clipping alignments'), pwd_log_file, keep_quiet)
    remove_both_ends_clp(input_reads_to_16s_sam_bowtie, input_reads_to_16s_sam_bowtie_single_clp)
    os.system('rm %s' % input_reads_to_16s_sam_bowtie)

    report_and_log(('Round 1: transforming cigar format from 1.3 to 1.4'), pwd_log_file, keep_quiet)
    bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (input_reads_to_16s_sam_bowtie_single_clp, input_reads_to_16s_sam, input_reads_to_16s_sam_reformat_log)
    os.system(bbmap_reformat_cmd)
    os.system('rm %s' % input_reads_to_16s_sam_bowtie_single_clp )


    ##################################################### parse sam file ####################################################

    report_and_log(('Round 1: keep only best alignments (lowest mismatch) for each read'), pwd_log_file, keep_quiet)
    #keep_best_matches_in_sam(input_reads_to_16s_sam, input_reads_to_16s_sam_best_match)
    #keep_best_matches_in_sam_keep_short_M(input_reads_to_16s_sam, min_M_len_16s, input_reads_to_16s_sam_best_match)
    keep_best_matches_in_sam_drop_short_M(input_reads_to_16s_sam, min_M_len_16s, input_reads_to_16s_sam_best_match)


    #os.system('rm %s' % input_reads_to_16s_sam)

    report_and_log(('Round 1: read in filtered sam file'), pwd_log_file, keep_quiet)
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
    report_and_log(('Round 1: extracting sequences of unmapped mates'), pwd_log_file, keep_quiet)
    get_unmapped_mates_seq(input_reads_to_16s_sam_best_match, reads_file_r1_fasta, reads_file_r2_fasta, sam_best_match_unmapped_mates_seq_file)

    report_and_log(('Round 1: adding sequences of unmapped mates to MappingRecord dictionary'), pwd_log_file, keep_quiet)
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

    report_and_log(('Round 1: analysing mappping results stored in MappingRecord dictionary'), pwd_log_file, keep_quiet)

    for each_mp in MappingRecord_dict:
        current_mp_record = MappingRecord_dict[each_mp]

        # only r1 mapped
        if (current_mp_record.r1_refs != {}) and (current_mp_record.r2_refs == {}):

            # print('%s\tcurrent_mp_r1_refs\t%s' % (each_mp, current_mp_record.r1_refs))
            # print('%s\tcurrent_mp_r2_refs\t%s' % (each_mp, current_mp_r2_refs))

            r1_refs_filtered_by_M_len = {}
            for r1_ref in current_mp_record.r1_refs:
                r1_ref_cigar = current_mp_record.r1_refs[r1_ref]
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(
                    cigar_splitter(r1_ref_cigar))
                if aligned_len >= min_M_len_16s:
                    r1_refs_filtered_by_M_len[r1_ref] = current_mp_record.r1_refs[r1_ref]

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
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(
                    cigar_splitter(r2_ref_cigar))
                if aligned_len >= min_M_len_16s:
                    r2_refs_filtered_by_M_len[r2_ref] = current_mp_record.r2_refs[r2_ref]

            if len(r2_refs_filtered_by_M_len) == 0:
                current_mp_record.r2_refs = {}
            else:
                # consider as paired
                current_mp_record.qualified_reads = True
                current_mp_record.consider_r2_unmapped_mate = True
                current_mp_record.r2_refs = r2_refs_filtered_by_M_len
                current_mp_record.r2_filtered_refs = {i.split('_pos_')[0] for i in r2_refs_filtered_by_M_len}

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
    report_and_log(('Round 1: writing out sequences of unmapped mates and clipping parts'), pwd_log_file, keep_quiet)

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
            #clipping_part_seq_handle.write('%s\n' % read_mr.r1_clipping_seq)
            clipping_part_seq_handle.write('%s\n' % read_mr.r1_seq)

        if read_mr.consider_r2_clipping_part is True:
            clipping_part_seq_handle.write('>%s\n' % r2_name)
            #clipping_part_seq_handle.write('%s\n' % read_mr.r2_clipping_seq)
            clipping_part_seq_handle.write('%s\n' % read_mr.r2_seq)

    unmapped_mates_handle.close()
    clipping_part_seq_handle.close()

    ############################# map clipping sequences and unmapped mates to combined input genomes #############################

    combined_input_gnms_no_ext = '.'.join(combined_input_gnms.split('.')[:-1])
    bowtie_build_input_gnm_cmd          = 'bowtie2-build --quiet --threads %s -f %s %s'                         % (num_threads, combined_input_gnms, combined_input_gnms_no_ext)
    bowtie_cmd_unmapped_to_mag          = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --all --no-unal --very-fast-local 2> %s'    % (combined_input_gnms_no_ext, unmapped_mates_seq_file, unmapped_to_gnm_sam_bowtie, num_threads, unmapped_to_gnm_sam_bowtie_log)
    bowtie_cmd_clipping_to_mag          = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --all --no-unal --very-fast-local 2> %s'    % (combined_input_gnms_no_ext, clipping_parts_seq_file, clipping_to_gnm_sam_bowtie, num_threads, clipping_to_gnm_sam_bowtie_log)
    bbmap_reformat_cmd_unmapped_to_mag  = 'reformat.sh in=%s out=%s sam=1.4 2> %s'                                    % (unmapped_to_gnm_sam_bowtie_single_clp, unmapped_to_gnm_sam, unmapped_to_gnm_sam_reformat_log)
    bbmap_reformat_cmd_clipping_to_mag  = 'reformat.sh in=%s out=%s sam=1.4 2> %s'                                    % (clipping_to_gnm_sam_bowtie_single_clp, clipping_to_gnm_sam, clipping_to_gnm_sam_reformat_log)
    #bbmap_cmd_unmapped_to_mag = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, combined_input_gnms, unmapped_mates_seq_file, unmapped_to_gnm_sam, bbmap_parameter, unmapped_to_gnm_bbmap_stderr)
    #bbmap_cmd_clipping_to_mag = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, combined_input_gnms, clipping_parts_seq_file, clipping_to_gnm_sam, bbmap_parameter, clipping_to_gnm_bbmap_stderr)

    # index ref
    os.system(bowtie_build_input_gnm_cmd)

    # map unmapped mates
    report_and_log(('Round 1: Mapping unmapped mates to genomic sequences'), pwd_log_file, keep_quiet)
    #report_and_log((bbmap_cmd_unmapped_to_mag), pwd_log_file, True)
    #os.system(bbmap_cmd_unmapped_to_mag)
    os.system(bowtie_cmd_unmapped_to_mag)
    remove_both_ends_clp(unmapped_to_gnm_sam_bowtie, unmapped_to_gnm_sam_bowtie_single_clp)
    os.system(bbmap_reformat_cmd_unmapped_to_mag)
    #keep_best_matches_in_sam(unmapped_to_gnm_sam, unmapped_to_gnm_sam_best_match)
    #keep_best_matches_in_sam_keep_short_M(unmapped_to_gnm_sam, min_M_len_ctg, unmapped_to_gnm_sam_best_match)
    keep_best_matches_in_sam_drop_short_M(unmapped_to_gnm_sam, min_M_len_ctg, unmapped_to_gnm_sam_best_match)

    # map clipping sequences
    report_and_log(('Round 1: Mapping clipping sequences to genomic sequences'), pwd_log_file, keep_quiet)
    #report_and_log((bbmap_cmd_clipping_to_mag), pwd_log_file, True)
    #os.system(bbmap_cmd_clipping_to_mag)
    os.system(bowtie_cmd_clipping_to_mag)
    remove_both_ends_clp(clipping_to_gnm_sam_bowtie, clipping_to_gnm_sam_bowtie_single_clp)
    os.system(bbmap_reformat_cmd_clipping_to_mag)
    #keep_best_matches_in_sam(clipping_to_gnm_sam, clipping_to_gnm_sam_best_match)
    #keep_best_matches_in_sam_keep_short_M(clipping_to_gnm_sam, min_M_len_ctg, clipping_to_gnm_sam_best_match)
    keep_best_matches_in_sam_drop_short_M(clipping_to_gnm_sam, min_M_len_ctg, clipping_to_gnm_sam_best_match)

    # rm tmp files
    os.system('rm %s' % unmapped_to_gnm_sam_bowtie)
    os.system('rm %s' % unmapped_to_gnm_sam_bowtie_single_clp )
    os.system('rm %s' % unmapped_to_gnm_sam)
    os.system('rm %s' % clipping_to_gnm_sam_bowtie)
    os.system('rm %s' % clipping_to_gnm_sam_bowtie_single_clp )
    os.system('rm %s' % clipping_to_gnm_sam)


    ######################################### parse mapping results for unmapped mates #########################################

    # full length alignment, or partially mapped to ctg end


    report_and_log(('Round 1: Parsing mapping results for unmapped mates'), pwd_log_file, keep_quiet)

    ctg_len_dict = {}
    for each_read in open(unmapped_to_gnm_sam_best_match):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            ctg_id = ''
            ctg_len = 0
            for each_element in each_read_split:
                if each_element.startswith('SN:'):
                    ctg_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ctg_len = int(each_element[3:])
            ctg_len_dict[ctg_id] = ctg_len
        else:
            read_id = each_read_split[0]
            cigar = each_read_split[5]
            if cigar != '*':
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_read_split[2]
                ref_pos = int(each_read_split[3])
                ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
                cigar_splitted = cigar_splitter(cigar)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                if (aligned_len >= min_M_len_ctg) and (aligned_pct >= min_M_pct) and (mismatch_pct <= mismatch_cutoff):

                    # full length alignment
                    if ('S' not in cigar) and ('s' not in cigar):
                        if read_strand == '1':
                            MappingRecord_dict[read_id_base].unmapped_r1_refs.add(ref_id)
                            MappingRecord_dict[read_id_base].unmapped_r1_refs_with_pos.add(ref_id_with_pos)
                        if read_strand == '2':
                            MappingRecord_dict[read_id_base].unmapped_r2_refs.add(ref_id)
                            MappingRecord_dict[read_id_base].unmapped_r2_refs_with_pos.add(ref_id_with_pos)
                    else:
                        # partially mapped to ctg end
                        ref_len = ctg_len_dict.get(ref_id, 0)
                        if (ref_pos <= 1) or ((ref_pos + aligned_len - 1) == ref_len):
                            if read_strand == '1':
                                MappingRecord_dict[read_id_base].unmapped_r1_refs.add(ref_id)
                                MappingRecord_dict[read_id_base].unmapped_r1_refs_with_pos.add(ref_id_with_pos)
                            if read_strand == '2':
                                MappingRecord_dict[read_id_base].unmapped_r2_refs.add(ref_id)
                                MappingRecord_dict[read_id_base].unmapped_r2_refs_with_pos.add(ref_id_with_pos)


    #################################### parse mapping results for clipping mapped reads #####################################

    # full length alignment, or partially mapped to ctg end

    report_and_log(('Round 1: Parsing mapping results for clipping sequences'), pwd_log_file, keep_quiet)

    for each_read in open(clipping_to_gnm_sam_best_match):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            ctg_id = ''
            ctg_len = 0
            for each_element in each_read_split:
                if each_element.startswith('SN:'):
                    ctg_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ctg_len = int(each_element[3:])
            ctg_len_dict[ctg_id] = ctg_len
        else:
            read_id = each_read_split[0]
            cigar = each_read_split[5]
            if cigar != '*':
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_read_split[2]
                ref_pos = int(each_read_split[3])
                ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
                cigar_splitted = cigar_splitter(cigar)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                if (aligned_len >= min_M_len_ctg) and (aligned_pct >= min_M_pct) and (mismatch_pct <= mismatch_cutoff):

                    # full length alignment
                    if ('S' not in cigar) and ('s' not in cigar):
                        if read_strand == '1':
                            MappingRecord_dict[read_id_base].clipping_r1_refs.add(ref_id)
                            MappingRecord_dict[read_id_base].clipping_r1_refs_with_pos.add(ref_id_with_pos)
                        if read_strand == '2':
                            MappingRecord_dict[read_id_base].clipping_r2_refs.add(ref_id)
                            MappingRecord_dict[read_id_base].clipping_r2_refs_with_pos.add(ref_id_with_pos)
                    else:
                        # partially mapped to ctg end
                        ref_len = ctg_len_dict.get(ref_id, 0)
                        if (ref_pos <= 1) or ((ref_pos + aligned_len - 1) == ref_len):
                            if read_strand == '1':
                                MappingRecord_dict[read_id_base].clipping_r1_refs.add(ref_id)
                                MappingRecord_dict[read_id_base].clipping_r1_refs_with_pos.add(ref_id_with_pos)
                            if read_strand == '2':
                                MappingRecord_dict[read_id_base].clipping_r2_refs.add(ref_id)
                                MappingRecord_dict[read_id_base].clipping_r2_refs_with_pos.add(ref_id_with_pos)


    ############################################## get pairwise_16s_iden_dict ##############################################

    report_and_log(('Round 1: Get pairwise 16S rRNA gene identities'), pwd_log_file, keep_quiet)

    # makeblastdn with marker gene sequences
    blastdb_16s         = '%s/%s%s' % (bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension)
    makeblastdb_16s_cmd = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blastdb_16s)
    os.system(makeblastdb_16s_cmd)

    all_vs_all_16s_blastn_cmd = '%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, blastdb_16s, blastdb_16s, blast_results_all_vs_all_16s, blast_parameters)
    os.system(all_vs_all_16s_blastn_cmd)

    pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

    # write out to file
    pairwise_marker_similarity_handle = open(pairwise_marker_similarity, 'w')
    pairwise_marker_similarity_handle.write('Marker1\tMarker2\tSimilarity\n')
    for marker_pair in pairwise_16s_iden_dict:
        pairwise_marker_similarity_handle.write('%s\t%s\n' % ('\t'.join(marker_pair.split('__|__')), pairwise_16s_iden_dict[marker_pair]))
    pairwise_marker_similarity_handle.close()


    ##################################################### get linkages from MappingRecord_dict #####################################################

    report_and_log(('Round 1: Parsing MappingRecord dictionary to get linkages'), pwd_log_file, keep_quiet)

    marker_to_ctg_linking_reads_dict_pair = {}
    marker_to_ctg_linking_reads_dict_clip = {}
    marker_to_ctg_link_num_dict_pair = {}
    marker_to_ctg_link_num_dict_clip = {}
    for qualified_read in MappingRecord_dict:

        read_mr = MappingRecord_dict[qualified_read]

        if (len(read_mr.unmapped_r2_refs) > 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) == 0):

            for r1_ref in read_mr.r1_filtered_refs:
                for unmapped_r2_ref in read_mr.unmapped_r2_refs:
                    marker_to_ctg_key = '%s%s%s' % (r1_ref, marker_to_ctg_gnm_Key_connector, unmapped_r2_ref)
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key] = {qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key].add(qualified_read)

        elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) > 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) == 0):

            for r2_ref in read_mr.r2_filtered_refs:
                for unmapped_r1_ref in read_mr.unmapped_r1_refs:
                    marker_to_ctg_key = '%s%s%s' % (r2_ref, marker_to_ctg_gnm_Key_connector, unmapped_r1_ref)
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key] = {qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key].add(qualified_read)


        elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) > 0) and (len(read_mr.clipping_r2_refs) == 0):

            for r1_ref in read_mr.r1_filtered_refs:
                for r1_clipping_ref in read_mr.clipping_r1_refs:
                    marker_to_ctg_key = '%s%s%s' % (r1_ref, marker_to_ctg_gnm_Key_connector, r1_clipping_ref)
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key] = {'%s.1' % qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key].add('%s.1' % qualified_read)

        elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) > 0):

            for r2_ref in read_mr.r2_filtered_refs:
                for r2_clipping_ref in read_mr.clipping_r2_refs:
                    marker_to_ctg_key = '%s%s%s' % (r2_ref, marker_to_ctg_gnm_Key_connector, r2_clipping_ref)
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key] = {'%s.2' % qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key].add('%s.2' % qualified_read)


        elif (len(read_mr.unmapped_r2_refs) > 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) > 0) and (len(read_mr.clipping_r2_refs) == 0):

            shared_refs = set(read_mr.unmapped_r2_refs).intersection(read_mr.clipping_r1_refs)

            for r1_ref in read_mr.r1_filtered_refs:
                for shared_ref in shared_refs:
                    marker_to_ctg_key = '%s%s%s' % (r1_ref, marker_to_ctg_gnm_Key_connector, shared_ref)

                    # add to marker_to_ctg_link_num_dict_pair
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key] = {qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key].add(qualified_read)

                    # add to marker_to_ctg_link_num_dict_clip
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key] = {'%s.1' % qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key].add('%s.1' % qualified_read)

        elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) > 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) > 0):

            shared_refs = set(read_mr.unmapped_r1_refs).intersection(read_mr.clipping_r2_refs)

            for r2_ref in read_mr.r2_filtered_refs:
                for shared_ref in shared_refs:
                    marker_to_ctg_key = '%s%s%s' % (r2_ref, marker_to_ctg_gnm_Key_connector, shared_ref)

                    # add to marker_to_ctg_link_num_dict_pair
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key] = {qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_pair[marker_to_ctg_key].add(qualified_read)

                    # add to marker_to_ctg_link_num_dict_clip
                    if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key] = {'%s.2' % qualified_read}
                    else:
                        marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1
                        marker_to_ctg_linking_reads_dict_clip[marker_to_ctg_key].add('%s.2' % qualified_read)

    # combine clip and pair dict
    marker_to_ctg_link_num_combined = {}
    for each_marker_to_ctg_key in marker_to_ctg_link_num_dict_pair:
        marker_to_ctg_link_num_combined[each_marker_to_ctg_key] = marker_to_ctg_link_num_dict_pair[each_marker_to_ctg_key]
    for each_marker_to_ctg_key in marker_to_ctg_link_num_dict_clip:
        if each_marker_to_ctg_key not in marker_to_ctg_link_num_combined:
            marker_to_ctg_link_num_combined[each_marker_to_ctg_key] = marker_to_ctg_link_num_dict_clip[each_marker_to_ctg_key]
        else:
            marker_to_ctg_link_num_combined[each_marker_to_ctg_key] += marker_to_ctg_link_num_dict_clip[each_marker_to_ctg_key]

    # remove linkages less than 3
    marker_to_ctg_link_num_combined_min_3 = {}
    for each_key in marker_to_ctg_link_num_combined:
        if marker_to_ctg_link_num_combined[each_key] >= 3:
            marker_to_ctg_link_num_combined_min_3[each_key] = marker_to_ctg_link_num_combined[each_key]

    marker_to_ctg_link_num_combined_to_use = marker_to_ctg_link_num_combined_min_3

    # get number of linkages at genome level
    marker_to_gnm_link_num_combined = {}
    for each_marker_to_ctg_key in marker_to_ctg_link_num_combined_to_use:
        marker_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[0]
        ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[1]
        gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
        marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
        if marker_to_gnm_key not in marker_to_gnm_link_num_combined:
            marker_to_gnm_link_num_combined[marker_to_gnm_key] = marker_to_ctg_link_num_combined_to_use[each_marker_to_ctg_key]
        else:
            marker_to_gnm_link_num_combined[marker_to_gnm_key] += marker_to_ctg_link_num_combined_to_use[each_marker_to_ctg_key]

    sankey_file_in_clipping_handle = open(link_stats_combined, 'w')
    sankey_file_in_clipping_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_clipping in marker_to_gnm_link_num_combined:
        sankey_file_in_clipping_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
        each_clipping.split(marker_to_ctg_gnm_Key_connector)[0],
        each_clipping.split(marker_to_ctg_gnm_Key_connector)[1],
        marker_to_gnm_link_num_combined[each_clipping]))
    sankey_file_in_clipping_handle.close()

    filter_linkages_iteratively(link_stats_combined, 'Number', pairwise_16s_iden_dict,
                                mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple,
                                min_iden_16s, min_link_num,
                                min_link_num, within_gnm_linkage_num_diff, link_stats_combined_filtered_s1)


    ####################################### get linking reads for visualization ########################################

    report_and_log(('Round 1: Extracting linking reads for visualization'), pwd_log_file, keep_quiet)

    linked_marker_to_gnm_key_set = set()
    for each_filtered_rd1_linkage in open(link_stats_combined_filtered_s1):
        if not each_filtered_rd1_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            marker_id = each_filtered_rd1_linkage.strip().split(',')[0][12:]
            genome_id = each_filtered_rd1_linkage.strip().split(',')[1][12:]
            marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, genome_id)
            linked_marker_to_gnm_key_set.add(marker_to_gnm_key)

    # write out linking reads
    reads_to_ref_dict = {}
    ctgs_to_extract = set()
    linked_ends_dict = {}
    reads_to_extract_all = set()
    reads_to_extract_pair = set()
    reads_to_extract_clip = set()
    linking_reads_pair_handle = open(linking_reads_pair, 'w')
    linking_reads_clip_handle = open(linking_reads_clip, 'w')
    linked_ends_rd1_txt_handle = open(linked_ends_rd1_txt, 'w')
    linked_ends_rd1_txt_handle.write('Marker\tContig\tMarker_end\tContig_end\tMarker_pos\tContig_pos\n')
    for each_linking_read in marker_to_ctg_link_num_combined_to_use:
        marker_id = each_linking_read.split(marker_to_ctg_gnm_Key_connector)[0]
        ctg_id = each_linking_read.split(marker_to_ctg_gnm_Key_connector)[1]
        gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
        marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
        if marker_to_gnm_key in linked_marker_to_gnm_key_set:
            current_linking_reads_pair = marker_to_ctg_linking_reads_dict_pair.get(each_linking_read, [])
            current_linking_reads_clip = marker_to_ctg_linking_reads_dict_clip.get(each_linking_read, [])
            marker_pos_list = []
            contig_pos_list = []
            if current_linking_reads_pair != []:
                linking_reads_pair_handle.write('%s\t%s\t%s\n' % (marker_id, ctg_id, ','.join(current_linking_reads_pair)))
                reads_to_extract_pair.update(current_linking_reads_pair)
                reads_to_extract_all.update(current_linking_reads_pair)

                for each_paired_read in current_linking_reads_pair:

                    # get reads_to_ref_dict
                    paired_read_r1_id = '%s.1' % each_paired_read
                    paired_read_r2_id = '%s.2' % each_paired_read
                    if paired_read_r1_id not in reads_to_ref_dict:
                        reads_to_ref_dict[paired_read_r1_id] = {each_linking_read}
                    else:
                        reads_to_ref_dict[paired_read_r1_id].add(each_linking_read)
                    if paired_read_r2_id not in reads_to_ref_dict:
                        reads_to_ref_dict[paired_read_r2_id] = {each_linking_read}
                    else:
                        reads_to_ref_dict[paired_read_r2_id].add(each_linking_read)

                    # get linked pos list
                    if MappingRecord_dict[each_paired_read].consider_r1_unmapped_mate is True:
                        for r1_ref in MappingRecord_dict[each_paired_read].r1_refs:
                            if (marker_id + '_pos_') in r1_ref:
                                marker_pos_list.append(int(r1_ref.split('_pos_')[1]))
                        for unmapped_r2_ref in MappingRecord_dict[each_paired_read].unmapped_r2_refs_with_pos:
                            if (ctg_id + '_pos_') in unmapped_r2_ref:
                                contig_pos_list.append(int(unmapped_r2_ref.split('_pos_')[1]))
                    if MappingRecord_dict[each_paired_read].consider_r2_unmapped_mate is True:
                        for r2_ref in MappingRecord_dict[each_paired_read].r2_refs:
                            if (marker_id + '_pos_') in r2_ref:
                                marker_pos_list.append(int(r2_ref.split('_pos_')[1]))
                        for unmapped_r1_ref in MappingRecord_dict[each_paired_read].unmapped_r1_refs_with_pos:
                            if (ctg_id + '_pos_') in unmapped_r1_ref:
                                contig_pos_list.append(int(unmapped_r1_ref.split('_pos_')[1]))

            if current_linking_reads_clip != []:
                linking_reads_clip_handle.write('%s\t%s\t%s\n' % (marker_id, ctg_id, ','.join(current_linking_reads_clip)))
                reads_to_extract_clip.update(current_linking_reads_clip)
                reads_to_extract_all.update({i.split('.')[0] for i in current_linking_reads_clip})

                for each_clipping_read in current_linking_reads_clip:

                    # get reads_to_ref_dict
                    if each_clipping_read not in reads_to_ref_dict:
                        reads_to_ref_dict[each_clipping_read] = {each_linking_read}
                    else:
                        reads_to_ref_dict[each_clipping_read].add(each_linking_read)

                    # get linked pos list
                    read_id_no_strand = '.'.join(each_clipping_read.split('.')[:-1])
                    if MappingRecord_dict[read_id_no_strand].consider_r1_clipping_part is True:
                        for r1_ref in MappingRecord_dict[read_id_no_strand].r1_refs:
                            if (marker_id + '_pos_') in r1_ref:
                                marker_pos_list.append(int(r1_ref.split('_pos_')[1]))
                        for clipping_r1_ref in MappingRecord_dict[read_id_no_strand].clipping_r1_refs_with_pos:
                            if (ctg_id + '_pos_') in clipping_r1_ref:
                                contig_pos_list.append(int(clipping_r1_ref.split('_pos_')[1]))
                    if MappingRecord_dict[read_id_no_strand].consider_r2_clipping_part is True:
                        for r2_ref in MappingRecord_dict[read_id_no_strand].r2_refs:
                            if (marker_id + '_pos_') in r2_ref:
                                marker_pos_list.append(int(r2_ref.split('_pos_')[1]))
                        for clipping_r2_ref in MappingRecord_dict[read_id_no_strand].clipping_r2_refs_with_pos:
                            if (ctg_id + '_pos_') in clipping_r2_ref:
                                contig_pos_list.append(int(clipping_r2_ref.split('_pos_')[1]))

            marker_pos_median = np.median(marker_pos_list)
            linked_end_marker = 'left'
            if marker_pos_median > (marker_len_dict[marker_id]/2):
                linked_end_marker = 'right'

            contig_pos_median = np.median(contig_pos_list)
            if contig_pos_median <= 500:
                linked_end_contig = 'left'
            elif (ctg_len_dict[ctg_id] - contig_pos_median) <= 500:
                linked_end_contig = 'right'
            else:
                linked_end_contig = 'middle(%s)' % int(round(contig_pos_median))

            ctgs_to_extract.add(ctg_id)
            linked_ends_rd1_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (marker_id, ctg_id, linked_end_marker, linked_end_contig, ','.join([str(i) for i in marker_pos_list]), ','.join([str(i) for i in contig_pos_list])))
            linked_ends_dict[each_linking_read] = {'marker': linked_end_marker, 'contig': linked_end_contig}

    linking_reads_pair_handle.close()
    linking_reads_clip_handle.close()
    linked_ends_rd1_txt_handle.close()

    ########## subset combined genome file ##########

    linked_contigs_txt_handle = open(linked_contigs_txt, 'w')
    linked_contigs_txt_handle.write('\n'.join(ctgs_to_extract) + '\n')
    linked_contigs_txt_handle.close()
    subset_cmd = 'seqtk subseq %s %s > %s' % (combined_input_gnms, linked_contigs_txt, linked_contigs_fasta)
    os.system(subset_cmd)

    ########## extract sequences for visualization ##########

    os.mkdir(mafft_seq_folder)

    mafft_seq_files = set()
    for each_link in linked_ends_dict:
        #print('%s\t%s' % (each_link, linked_ends_dict[each_link]))

        marker_id = each_link.split(marker_to_ctg_gnm_Key_connector)[0]
        contig_id = each_link.split(marker_to_ctg_gnm_Key_connector)[1]

        seq_file_name_tmp = each_link.replace(marker_to_ctg_gnm_Key_connector, '___')
        pwd_seq_file = '%s/%s.fasta' % (mafft_seq_folder, seq_file_name_tmp.replace(gnm_to_ctg_connector, '___'))
        mafft_seq_files.add(pwd_seq_file)
        pwd_seq_file_handle = open(pwd_seq_file, 'w')

        # get 16S sequences
        marker_seq = ''
        for each_16s in SeqIO.parse(marker_gene_seqs, 'fasta'):
            if each_16s.id == marker_id:
                marker_seq = str(each_16s.seq)

        # extract contig sequence
        contig_seq = ''
        for each_ctg in SeqIO.parse(linked_contigs_fasta, 'fasta'):
            if each_ctg.id == contig_id:
                contig_seq = str(each_ctg.seq)

        # get contig end sequence
        contig_seq_for_mafft = ''
        if linked_ends_dict[each_link]['contig'] == 'left':
            contig_seq_for_mafft = contig_seq[:end_ctg_len_for_mafft]
        elif linked_ends_dict[each_link]['contig'] == 'right':
            contig_seq_for_mafft = contig_seq[-end_ctg_len_for_mafft:]
        else:
            median_value = int(linked_ends_dict[each_link]['contig'][7:-1])

            left_end_pos = 0
            if (median_value - end_ctg_len_for_mafft) > 0:
                left_end_pos = median_value - end_ctg_len_for_mafft
            right_end_pos = median_value + end_ctg_len_for_mafft
            if right_end_pos > ctg_len_dict[contig_id]:
                right_end_pos = ctg_len_dict[contig_id]

            contig_seq_for_mafft = contig_seq[left_end_pos:(right_end_pos - 1)]

        # concatenate 16s and contig sequences
        if linked_ends_dict[each_link]['contig'] in ['left', 'right']:
            concatenated_seq = ''
            if (linked_ends_dict[each_link]['marker'] == 'right') and (linked_ends_dict[each_link]['contig'] == 'left'):
                concatenated_seq = '%s%s%s' % (marker_seq, 'N'*50, contig_seq_for_mafft)
            if (linked_ends_dict[each_link]['marker'] == 'left') and (linked_ends_dict[each_link]['contig'] == 'right'):
                concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N'*50, marker_seq)
            if (linked_ends_dict[each_link]['marker'] == 'left') and (linked_ends_dict[each_link]['contig'] == 'left'):
                marker_seq_rc = get_rc(marker_seq)
                concatenated_seq = '%s%s%s' % (marker_seq_rc, 'N'*50, contig_seq_for_mafft)
            if (linked_ends_dict[each_link]['marker'] == 'right') and (linked_ends_dict[each_link]['contig'] == 'right'):
                marker_seq_rc = get_rc(marker_seq)
                concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N'*50, marker_seq_rc)
            pwd_seq_file_handle.write('>Reference\n')
            pwd_seq_file_handle.write('%s\n' % concatenated_seq)
        else:
            pwd_seq_file_handle.write('>Contig\n')
            pwd_seq_file_handle.write('%s\n' % contig_seq_for_mafft)
            pwd_seq_file_handle.write('>Marker\n')
            pwd_seq_file_handle.write('%s\n' % marker_seq)

        pwd_seq_file_handle.close()

    ########## extract reads for visualization ##########

    for each_read in MappingRecord_dict:

        read_id_r1 = '%s.1' % each_read
        read_id_r2 = '%s.2' % each_read
        r1_refs = reads_to_ref_dict.get(read_id_r1, [])
        r2_refs = reads_to_ref_dict.get(read_id_r2, [])

        for r1_ref in r1_refs:
            r1_ref_refmt_tmp = r1_ref.replace(marker_to_ctg_gnm_Key_connector, '___')
            r1_ref_refmt = r1_ref_refmt_tmp.replace(gnm_to_ctg_connector, '___')
            pwd_seq_file = '%s/%s.fasta' % (mafft_seq_folder, r1_ref_refmt)
            with open(pwd_seq_file, 'a') as r1_mafft_handle:
                r1_mafft_handle.write('>%s\n' % read_id_r1)
                r1_mafft_handle.write('%s\n' % MappingRecord_dict[each_read].r1_seq)

        for r2_ref in r2_refs:
            r2_ref_refmt_tmp = r2_ref.replace(marker_to_ctg_gnm_Key_connector, '___')
            r2_ref_refmt = r2_ref_refmt_tmp.replace(gnm_to_ctg_connector, '___')
            pwd_seq_file = '%s/%s.fasta' % (mafft_seq_folder, r2_ref_refmt)
            with open(pwd_seq_file, 'a') as r2_mafft_handle:
                r2_mafft_handle.write('>%s\n' % read_id_r2)
                r2_mafft_handle.write('%s\n' % MappingRecord_dict[each_read].r2_seq)

    ########## align and visualize ##########

    # if vis_linkages is True:
    #     report_and_log(('Round 1: Generating linking plots'), pwd_log_file, keep_quiet)
    #
    #     for mafft_seq_file in mafft_seq_files:
    #         html_page_title = mafft_seq_file.split('/')[-1]
    #         mview_linkage(mafft_seq_file, html_page_title, mafft_exe, mview_exe)


    ####################################################################################################################
    ############################################### second round linking ###############################################
    ####################################################################################################################

    #################### get the sequences of 1st round unlinked marker genes and genomic sequences ####################

    os.mkdir(step_2_wd)
    report_and_log(('Round 2: get unlinked marker genes and genomes'), pwd_log_file, keep_quiet)

    # get linked marker genes and genomic sequences in step 1
    linked_marker_gene_set = set()
    linked_genomic_seq_set = set()
    for each_link in open(link_stats_combined_filtered_s1):
        if not each_link.startswith('MarkerGene,GenomicSeq,Number'):
            each_link_split = each_link.strip().split(',')
            linked_marker_gene_set.add(each_link_split[0][12:])
            linked_genomic_seq_set.add(each_link_split[1][12:])

    # get the sequence of unlinked genomic seqs
    # if genomic_seq_type == 'mag':

    # put all renamed mag into list
    renamed_gnm_re = '%s/*.%s' % (renamed_mag_folder, mag_file_extension)
    renamed_gnm_list = [os.path.basename(file_name) for file_name in glob.glob(renamed_gnm_re)]
    renamed_gnm_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in renamed_gnm_list]

    # keep only unlinked mags
    unlinked_mag_list_with_pwd = []
    for renamed_mag in renamed_gnm_list_no_ext:
        if renamed_mag not in linked_genomic_seq_set:
            pwd_renamed_mag = '%s/%s.%s' % (renamed_mag_folder, renamed_mag, mag_file_extension)
            unlinked_mag_list_with_pwd.append(pwd_renamed_mag)

    # combine unlinked mags
    cat_cmd = 'cat %s > %s' % (' '.join(unlinked_mag_list_with_pwd), combined_1st_round_unlinked_mags)
    os.system(cat_cmd)

    # # get the sequence of unlinked metagenomic assemblies
    # if genomic_seq_type == 'ctg':
    #     combined_1st_round_unlinked_ctgs_handle = open(combined_1st_round_unlinked_ctgs, 'w')
    #     for ctg_record in SeqIO.parse(genomic_assemblies, 'fasta'):
    #         if ctg_record.id not in linked_genomic_seq_set:
    #             SeqIO.write(ctg_record, combined_1st_round_unlinked_ctgs_handle, 'fasta')
    #     combined_1st_round_unlinked_ctgs_handle.close()
    #

    ######################################## extract sequences flanking 16S ends #######################################

    free_living_16s_refs_file_handle = open(free_living_16s_ref_file, 'w')
    free_living_16s_R1_handle = open(free_living_16s_R1, 'w')
    free_living_16s_R2_handle = open(free_living_16s_R2, 'w')
    free_living_16s_UP_handle = open(free_living_16s_UP, 'w')
    for qualified_read in MappingRecord_dict:
        read_mr = MappingRecord_dict[qualified_read]

        if (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) == 0):

            ref_been_linked = False
            for r1_ref in read_mr.r1_filtered_refs:
                if r1_ref in linked_marker_gene_set:
                    ref_been_linked = True
            for r2_ref in read_mr.r2_filtered_refs:
                if r2_ref in linked_marker_gene_set:
                    ref_been_linked = True

            if ref_been_linked is False:
                read_mr.consider_round_2 = True

                # filter r1 refs according to round 2 mismatch cutoff
                r1_filtered_refs_rd2 = set()
                for each_r1_ref in read_mr.r1_refs:
                    r1_ref_no_pos = each_r1_ref.split('_pos_')[0]
                    r1_ref_cigar = read_mr.r1_refs[each_r1_ref]
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r1_ref_cigar))
                    if mismatch_pct <= mismatch_cutoff:
                        r1_filtered_refs_rd2.add(r1_ref_no_pos)

                # filter r2 refs according to round 2 mismatch cutoff
                r2_filtered_refs_rd2 = set()
                current_r2_refs = read_mr.r2_refs
                for each_r2_ref in read_mr.r2_refs:
                    r2_ref_no_pos = each_r2_ref.split('_pos_')[0]
                    r2_ref_cigar = read_mr.r2_refs[each_r2_ref]
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r2_ref_cigar))
                    if mismatch_pct <= mismatch_cutoff:
                        r2_filtered_refs_rd2.add(r2_ref_no_pos)

                if (read_mr.consider_r1_unmapped_mate is True) and (read_mr.consider_r1_clipping_part is True):

                    if clp_read_for_assembly is True:
                        free_living_16s_refs_file_handle.write('%s.1\t%s\n' % (qualified_read, ','.join(read_mr.r1_filtered_refs)))
                        free_living_16s_R1_handle.write('>%s.1\n' % qualified_read)
                        free_living_16s_R1_handle.write('%s\n' % read_mr.r1_seq)
                        free_living_16s_refs_file_handle.write('%s.2\t%s\n' % (qualified_read, ','.join(r1_filtered_refs_rd2)))
                        free_living_16s_R2_handle.write('>%s.2\n' % qualified_read)
                        free_living_16s_R2_handle.write('%s\n' % read_mr.r2_seq)
                    else:
                        free_living_16s_refs_file_handle.write('%s.2\t%s\n' % (qualified_read, ','.join(r1_filtered_refs_rd2)))
                        free_living_16s_UP_handle.write('>%s.2\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r2_seq)

                elif (read_mr.consider_r2_unmapped_mate is True) and (read_mr.consider_r2_clipping_part is True):

                    if clp_read_for_assembly is True:
                        free_living_16s_refs_file_handle.write('%s.1\t%s\n' % (qualified_read, ','.join(r2_filtered_refs_rd2)))
                        free_living_16s_R1_handle.write('>%s.1\n' % qualified_read)
                        free_living_16s_R1_handle.write('%s\n' % read_mr.r1_seq)
                        free_living_16s_refs_file_handle.write('%s.2\t%s\n' % (qualified_read, ','.join(read_mr.r2_filtered_refs)))
                        free_living_16s_R2_handle.write('>%s.2\n' % qualified_read)
                        free_living_16s_R2_handle.write('%s\n' % read_mr.r2_seq)
                    else:
                        free_living_16s_refs_file_handle.write('%s.1\t%s\n' % (qualified_read, ','.join(r2_filtered_refs_rd2)))
                        free_living_16s_UP_handle.write('>%s.1\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_seq)

                else:
                    if read_mr.consider_r1_unmapped_mate is True:
                        free_living_16s_refs_file_handle.write('%s.2\t%s\n' % (qualified_read, ','.join(r1_filtered_refs_rd2)))
                        free_living_16s_UP_handle.write('>%s.2\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r2_seq)

                    if read_mr.consider_r2_unmapped_mate is True:
                        free_living_16s_refs_file_handle.write('%s.1\t%s\n' % (qualified_read, ','.join(r2_filtered_refs_rd2)))
                        free_living_16s_UP_handle.write('>%s.1\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_seq)

                    if clp_read_for_assembly is True:

                        if read_mr.consider_r1_clipping_part is True:
                            free_living_16s_refs_file_handle.write('%s.1\t%s\n' % (qualified_read, ','.join(read_mr.r1_filtered_refs)))
                            free_living_16s_UP_handle.write('>%s.1\n' % qualified_read)
                            free_living_16s_UP_handle.write('%s\n' % read_mr.r1_seq)

                        if read_mr.consider_r2_clipping_part is True:
                            free_living_16s_refs_file_handle.write('%s.2\t%s\n' % (qualified_read, ','.join(read_mr.r2_filtered_refs)))
                            free_living_16s_UP_handle.write('>%s.2\n' % qualified_read)
                            free_living_16s_UP_handle.write('%s\n' % read_mr.r2_seq)

    free_living_16s_refs_file_handle.close()
    free_living_16s_R1_handle.close()
    free_living_16s_R2_handle.close()
    free_living_16s_UP_handle.close()


    ######################################## extract sequences flanking ctg ends #######################################

    get_unlinked_mag_end_seq(combined_1st_round_unlinked_mags, combined_1st_round_unlinked_mag_end_seq, end_seq_len)

    # index reference
    bowtie_build_unlinked_ctg_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, combined_1st_round_unlinked_mag_end_seq, rd1_unlinked_mag_end_seq_no_ext)
    os.system(bowtie_build_unlinked_ctg_cmd)

    # mapping with bowtie
    bowtie_cmd_unlinked_ctg = 'bowtie2 -x %s -1 %s -2 %s -S %s -p %s -f --local --all --no-unal --very-fast-local 2> %s' % (rd1_unlinked_mag_end_seq_no_ext, reads_file_r1_fasta, reads_file_r2_fasta, rd1_unlinked_mags_sam_bowtie, num_threads, rd1_unlinked_mags_sam_bowtie_log)
    os.system(bowtie_cmd_unlinked_ctg)

    # remove both_ends_clp
    remove_both_ends_clp(rd1_unlinked_mags_sam_bowtie, rd1_unlinked_mags_sam_bowtie_single_clp)
    bbmap_reformat_rd1_unlinked_mags_sam  = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (rd1_unlinked_mags_sam_bowtie_single_clp, rd1_unlinked_mags_sam_bowtie_single_clp_reformat, rd1_unlinked_mags_sam_bowtie_single_clp_reformat_log)
    os.system(bbmap_reformat_rd1_unlinked_mags_sam)
    #keep_best_matches_in_sam(rd1_unlinked_mags_sam_bowtie_single_clp_reformat, combined_1st_round_unlinked_mags_sam_best_match)
    #keep_best_matches_in_sam_keep_short_M(rd1_unlinked_mags_sam_bowtie_single_clp_reformat, min_M_len_ctg, combined_1st_round_unlinked_mags_sam_best_match)
    keep_best_matches_in_sam_drop_short_M(rd1_unlinked_mags_sam_bowtie_single_clp_reformat, min_M_len_ctg, combined_1st_round_unlinked_mags_sam_best_match)


    get_unmapped_mates_seq(combined_1st_round_unlinked_mags_sam_best_match, reads_file_r1_fasta, reads_file_r2_fasta, combined_1st_round_unlinked_mags_sam_best_match_unmapped_mates)

    # rm tmp files
    os.system('rm %s' % rd1_unlinked_mags_sam_bowtie)
    os.system('rm %s' % rd1_unlinked_mags_sam_bowtie_single_clp )
    os.system('rm %s' % rd1_unlinked_mags_sam_bowtie_single_clp_reformat)

    report_and_log(('Round 2: Analysing mapping results'), pwd_log_file, keep_quiet)
    # parse sam file
    round_2_ctg_end_seq_len_dict = {}
    round_2_MappingRecord_dict = {}
    for each_line in open(combined_1st_round_unlinked_mags_sam_best_match):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('@'):
            rd2_ref_id = ''
            rd2_ref_len = 0
            for each_element in each_line_split:
                if each_element.startswith('SN:'):
                    rd2_ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    rd2_ref_len = int(each_element[3:])
            round_2_ctg_end_seq_len_dict[rd2_ref_id] = rd2_ref_len
        else:
            store_read_seq = False
            read_id = each_line_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            read_flag = int(each_line_split[1])
            cigar = each_line_split[5]
            read_seq = each_line_split[9]
            if cigar != '*':
                ref_id = each_line_split[2]
                ref_pos = each_line_split[3]
                ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
                cigar_splitted = cigar_splitter(cigar)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                if (mismatch_pct <= mismatch_cutoff) and (aligned_len >= min_M_len_ctg):

                    if read_id_base not in round_2_MappingRecord_dict:
                        round_2_MappingRecord_dict[read_id_base] = MappingRecord()

                    if read_strand == '1':
                        round_2_MappingRecord_dict[read_id_base].r1_refs[ref_id_with_pos] = cigar
                        round_2_MappingRecord_dict[read_id_base].r1_cigar_to_flag[cigar] = read_flag
                        round_2_MappingRecord_dict[read_id_base].r1_filtered_refs.add(ref_id)

                    if read_strand == '2':
                        round_2_MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar
                        round_2_MappingRecord_dict[read_id_base].r2_cigar_to_flag[cigar] = read_flag
                        round_2_MappingRecord_dict[read_id_base].r2_filtered_refs.add(ref_id)

                    if clipping_len >= min_M_len_ctg:
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

                if read_id_base not in round_2_MappingRecord_dict:
                    round_2_MappingRecord_dict[read_id_base] = MappingRecord()
                if read_strand == '1':
                    if round_2_MappingRecord_dict[read_id_base].r1_seq == '':
                        round_2_MappingRecord_dict[read_id_base].r1_seq = read_seq_to_store
                if read_strand == '2':
                    if round_2_MappingRecord_dict[read_id_base].r2_seq == '':
                        round_2_MappingRecord_dict[read_id_base].r2_seq = read_seq_to_store

    # add sequences of unmapped mates to mp dict
    for each_read in SeqIO.parse(combined_1st_round_unlinked_mags_sam_best_match_unmapped_mates, 'fasta'):
        read_id = str(each_read.id)
        read_basename = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        if read_basename in round_2_MappingRecord_dict:
            if read_strand == '1':
                round_2_MappingRecord_dict[read_basename].r1_seq = str(each_read.seq)
            if read_strand == '2':
                round_2_MappingRecord_dict[read_basename].r2_seq = str(each_read.seq)

    # parse round_2_MappingRecord_dict
    free_living_ctg_refs_file_handle = open(free_living_ctg_ref_file, 'w')
    free_living_ctg_R1_handle = open(free_living_ctg_R1, 'w')
    free_living_ctg_R2_handle = open(free_living_ctg_R2, 'w')
    free_living_ctg_UP_handle = open(free_living_ctg_UP, 'w')
    for read_basename in round_2_MappingRecord_dict.copy():
        read_mr = round_2_MappingRecord_dict[read_basename]
        r1_ref_cigar_list = list(read_mr.r1_refs.values())
        r2_ref_cigar_list = list(read_mr.r2_refs.values())
        best_cigar, max_clp, max_clp_location = get_max_clp_and_index(r1_ref_cigar_list, r2_ref_cigar_list)

        if (read_mr.r1_refs == {}) and (read_mr.r2_refs != {}):

            if len(read_mr.r2_refs) == 1:

                r2_ref_cigar = r2_ref_cigar_list[0]
                r2_ref_flag = read_mr.r2_cigar_to_flag[r2_ref_cigar]
                r2_ref_cigar_rc = sam_flag_to_rc(r2_ref_flag)

                # consider the unmapped mate only
                if max_clp < min_M_len_ctg:
                    read_mr.consider_round_2 = True
                    read_mr.consider_r2_unmapped_mate = True

                    # write out refs
                    free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))

                    # write out sequence
                    free_living_ctg_UP_handle.write('>%s.1\n' % read_basename)
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_seq)

                # consider both of unmapped mate and clipping part
                else:
                    r2_ref_no_pos = list(read_mr.r2_refs.keys())[0].split('_pos_')[0]
                    r2_ref_pos = int(list(read_mr.r2_refs.keys())[0].split('_pos_')[1])

                    if (r2_ref_no_pos[-1]) == (max_clp_location[-1]) == 'l':
                        read_mr.consider_round_2 = True
                        read_mr.consider_r2_unmapped_mate = True
                        if clp_read_for_assembly is True:
                            read_mr.consider_r2_clipping_part = True

                        # write out refs
                        free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))
                        if clp_read_for_assembly is True:
                            free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))

                        # if r2_ref_cigar_rc is False:
                        #     read_mr.r2_clipping_seq = read_mr.r2_seq[:max_clp]
                        # if r2_ref_cigar_rc is True:
                        #     r2_seq_rc = get_rc(read_mr.r2_seq)
                        #     r2_clipping_seq_rc = r2_seq_rc[:max_clp]
                        #     read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)

                    elif (r2_ref_no_pos[-1]) == (max_clp_location[-1]) == 'r':
                        read_mr.consider_round_2 = True
                        read_mr.consider_r2_unmapped_mate = True
                        if clp_read_for_assembly is True:
                            read_mr.consider_r2_clipping_part = True

                        # write out refs
                        free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))
                        if clp_read_for_assembly is True:
                            free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))

                        # if r2_ref_cigar_rc is False:
                        #     read_mr.r2_clipping_seq = read_mr.r2_seq[-max_clp:]
                        # if r2_ref_cigar_rc is True:
                        #     r2_seq_rc = get_rc(read_mr.r2_seq)
                        #     r2_clipping_seq_rc = r2_seq_rc[-max_clp:]
                        #     read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)

                    else:  # mapped to unwanted end, ignore
                        round_2_MappingRecord_dict.pop(read_basename)

                    # write out sequence
                    if read_mr.consider_round_2 is True:

                        if clp_read_for_assembly is True:
                            free_living_ctg_R1_handle.write('>%s.1\n' % read_basename)
                            free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_seq)
                            free_living_ctg_R2_handle.write('>%s.2\n' % read_basename)
                            free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_seq)
                        else:
                            free_living_ctg_UP_handle.write('>%s.1\n' % read_basename)
                            free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_seq)

            else:  # r2 mapped to multiple refs, ignore
                round_2_MappingRecord_dict.pop(read_basename)

        elif (read_mr.r1_refs != {}) and (read_mr.r2_refs == {}):

            if len(read_mr.r1_refs) == 1:

                r1_ref_cigar = r1_ref_cigar_list[0]
                r1_ref_flag = read_mr.r1_cigar_to_flag[r1_ref_cigar]
                r1_ref_cigar_rc = sam_flag_to_rc(r1_ref_flag)

                # consider the unmapped mate only
                if max_clp < min_M_len_ctg:
                    read_mr.consider_round_2 = True
                    read_mr.consider_r1_unmapped_mate = True

                    # write out refs
                    free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))

                    # write out sequence
                    free_living_ctg_UP_handle.write('>%s.2\n' % read_basename)
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_seq)

                # consider both of unmapped mate and clipping part
                else:
                    r1_ref_no_pos = list(read_mr.r1_refs.keys())[0].split('_pos_')[0]
                    r1_ref_pos = int(list(read_mr.r1_refs.keys())[0].split('_pos_')[1])

                    if (r1_ref_no_pos[-1]) == (max_clp_location[-1]) == 'l':
                        read_mr.consider_round_2 = True
                        read_mr.consider_r1_unmapped_mate = True
                        if clp_read_for_assembly is True:
                            read_mr.consider_r1_clipping_part = True

                        # write out refs
                        if clp_read_for_assembly is True:
                            free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))
                        free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))

                        # if r1_ref_cigar_rc is False:
                        #     read_mr.r1_clipping_seq = read_mr.r1_seq[:max_clp]
                        # if r1_ref_cigar_rc is True:
                        #     r1_seq_rc = get_rc(read_mr.r1_seq)
                        #     r1_clipping_seq_rc = r1_seq_rc[:max_clp]
                        #     read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)

                    elif (r1_ref_no_pos[-1]) == (max_clp_location[-1]) == 'r':
                        read_mr.consider_round_2 = True
                        read_mr.consider_r1_unmapped_mate = True
                        if clp_read_for_assembly is True:
                            read_mr.consider_r1_clipping_part = True

                        # write out refs
                        if clp_read_for_assembly is True:
                            free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))
                        free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))

                        # if r1_ref_cigar_rc is False:
                        #     read_mr.r1_clipping_seq = read_mr.r1_seq[-max_clp:]
                        # if r1_ref_cigar_rc is True:
                        #     r1_seq_rc = get_rc(read_mr.r1_seq)
                        #     r1_clipping_seq_rc = r1_seq_rc[-max_clp:]
                        #     read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)

                    else:  # mapped to unwanted end, ignore
                        round_2_MappingRecord_dict.pop(read_basename)

                    # write out sequence
                    if read_mr.consider_round_2 is True:
                        # write out R1 fa
                        if clp_read_for_assembly is True:
                            free_living_ctg_R1_handle.write('>%s.1\n' % read_basename)
                            free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_seq)
                            free_living_ctg_R2_handle.write('>%s.2\n' % read_basename)
                            free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_seq)
                        else:
                            free_living_ctg_UP_handle.write('>%s.2\n' % read_basename)
                            free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_seq)

            else:  # r1 mapped to multiple refs, ignore
                round_2_MappingRecord_dict.pop(read_basename)

        elif (read_mr.r1_refs != {}) and (read_mr.r2_refs != {}):
            if max_clp >= min_M_len_ctg:
                if clp_read_for_assembly is True:

                    if (len(read_mr.r1_refs) == 1) and (len(read_mr.r2_refs) == 1):
                        r1_ref_no_pos   = list(read_mr.r1_refs.keys())[0].split('_pos_')[0]
                        r2_ref_no_pos   = list(read_mr.r2_refs.keys())[0].split('_pos_')[0]
                        r1_ref_pos      = int(list(read_mr.r1_refs.keys())[0].split('_pos_')[1])
                        r2_ref_pos      = int(list(read_mr.r2_refs.keys())[0].split('_pos_')[1])

                        if r1_ref_no_pos == r2_ref_no_pos:
                            if (r1_ref_no_pos[-1]) == (r2_ref_no_pos[-1]) == (max_clp_location[-1]):

                                if (max_clp_location == 'r1_l') and (r1_ref_pos == 1):
                                    read_mr.consider_round_2 = True
                                    read_mr.consider_r1_clipping_part = True

                                    # write out refs
                                    free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))

                                    # best_cigar_flag = read_mr.r1_cigar_to_flag[best_cigar]
                                    # best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                                    # if best_cigar_rc is False:
                                    #     read_mr.r1_clipping_seq = read_mr.r1_seq[:max_clp]
                                    # if best_cigar_rc is True:
                                    #     r1_seq_rc = get_rc(read_mr.r1_seq)
                                    #     r1_clipping_seq_rc = r1_seq_rc[:max_clp]
                                    #     read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)

                                elif (max_clp_location == 'r2_l') and (r2_ref_pos == 1):
                                    read_mr.consider_round_2 = True
                                    read_mr.consider_r2_clipping_part = True

                                    # write out refs
                                    free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))

                                    # best_cigar_flag = read_mr.r2_cigar_to_flag[best_cigar]
                                    # best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                                    # if best_cigar_rc is False:
                                    #     read_mr.r2_clipping_seq = read_mr.r2_seq[:max_clp]
                                    # if best_cigar_rc is True:
                                    #     r2_seq_rc = get_rc(read_mr.r2_seq)
                                    #     r2_clipping_seq_rc = r2_seq_rc[:max_clp]
                                    #     read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)

                                elif (max_clp_location == 'r1_r') and (r1_ref_pos >= (end_seq_len / 5)):
                                    read_mr.consider_round_2 = True
                                    read_mr.consider_r1_clipping_part = True

                                    # write out refs
                                    free_living_ctg_refs_file_handle.write('%s.1\t%s\n' % (read_basename, ','.join(read_mr.r1_filtered_refs)))

                                    # best_cigar_flag = read_mr.r1_cigar_to_flag[best_cigar]
                                    # best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                                    # if best_cigar_rc is False:
                                    #     read_mr.r1_clipping_seq = read_mr.r1_seq[-max_clp:]
                                    # if best_cigar_rc is True:
                                    #     r1_seq_rc = get_rc(read_mr.r1_seq)
                                    #     r1_clipping_seq_rc = r1_seq_rc[-max_clp:]
                                    #     read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)

                                elif (max_clp_location == 'r2_r') and (r2_ref_pos >= (end_seq_len / 5)):
                                    read_mr.consider_round_2 = True
                                    read_mr.consider_r2_clipping_part = True

                                    # write out refs
                                    free_living_ctg_refs_file_handle.write('%s.2\t%s\n' % (read_basename, ','.join(read_mr.r2_filtered_refs)))

                                    # best_cigar_flag = read_mr.r2_cigar_to_flag[best_cigar]
                                    # best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                                    # if best_cigar_rc is False:
                                    #     read_mr.r2_clipping_seq = read_mr.r2_seq[-max_clp:]
                                    # if best_cigar_rc is True:
                                    #     r2_seq_rc = get_rc(read_mr.r2_seq)
                                    #     r2_clipping_seq_rc = r2_seq_rc[-max_clp:]
                                    #     read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)

                                else:  # not too many of them, ignore now, maybe worth check later,
                                    round_2_MappingRecord_dict.pop(read_basename)
                                    # print('%s\tr1_refs:\t%s\t%s' % (read_basename, read_mr.r1_refs, read_mr.r1_seq))
                                    # print('%s\tr2_refs:\t%s\t%s' % (read_basename, read_mr.r2_refs, read_mr.r2_seq))

                                # write out sequence
                                if read_mr.consider_round_2 is True:

                                    if read_mr.consider_r1_clipping_part is True:
                                        free_living_ctg_UP_handle.write('>%s.1\n' % read_basename)
                                        free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_seq)

                                    if read_mr.consider_r2_clipping_part is True:
                                        free_living_ctg_UP_handle.write('>%s.2\n' % read_basename)
                                        free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_seq)

                            else:  # mapped to unwanted end, ignore
                                round_2_MappingRecord_dict.pop(read_basename)

                        else:  # r1 and r2 mapped to different refs
                            # not too many of them, ignore now
                            round_2_MappingRecord_dict.pop(read_basename)

                    else:  # r1 or r2 mapped to multiple refs
                        # not too many of them, ignore now
                        round_2_MappingRecord_dict.pop(read_basename)

            else:  # ignore and remove element from dict
                round_2_MappingRecord_dict.pop(read_basename)

        else:  # ignore and remove element from dict
            round_2_MappingRecord_dict.pop(read_basename)

    free_living_ctg_R1_handle.close()
    free_living_ctg_R2_handle.close()
    free_living_ctg_UP_handle.close()
    free_living_ctg_refs_file_handle.close()

    os.system('cat %s %s %s > %s' % (free_living_16s_R1, free_living_16s_R2, free_living_16s_UP, free_living_16s))
    os.system('cat %s %s %s > %s' % (free_living_ctg_R1, free_living_ctg_R2, free_living_ctg_UP, free_living_ctg))

    os.system('cat %s %s > %s' % (free_living_16s_R1, free_living_ctg_R1, free_living_R1))
    os.system('cat %s %s > %s' % (free_living_16s_R2, free_living_ctg_R2, free_living_R2))
    os.system('cat %s %s > %s' % (free_living_16s_UP, free_living_ctg_UP, free_living_UP))


    ####################################################################################################################
    ######################################### second round linking by assembly #########################################
    ####################################################################################################################

    os.system('cat %s %s > %s' % (free_living_16s, free_living_ctg, free_living_all))

    # assemble
    if round_2_spades is False:

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
        #spades_cmd = '%s --only-assembler -s %s -o %s -t %s -k 59,75,99,127 > %s' % (pwd_spades_exe, free_living_all, spades_wd, num_threads, spades_log)
        #spades_cmd = '%s --only-assembler -s %s -o %s -t %s -k 75,99,127 > %s' % (pwd_spades_exe, free_living_all, spades_wd, num_threads, spades_log)
        #spades_cmd = '%s --only-assembler --careful -s %s -o %s -t %s -k 59,75,99,127 > %s' % (pwd_spades_exe, free_living_all, spades_wd, num_threads, spades_log)
        if clp_read_for_assembly is True:
            spades_cmd = '%s --only-assembler --meta -1 %s -2 %s -s %s -o %s -t %s -k 49,75,99,127 > %s' % (pwd_spades_exe, free_living_R1, free_living_R2, free_living_UP, spades_wd, num_threads, spades_log)
        else:
            spades_cmd = '%s --only-assembler --careful -s %s -o %s -t %s -k 49,75,99,127 > %s' % (pwd_spades_exe, free_living_all, spades_wd, num_threads, spades_log)

        report_and_log((spades_cmd), pwd_log_file, True)

        os.system(spades_cmd)
        mini_assemblies = '%s/scaffolds.fasta' % spades_wd

    # mapping extracted reads to mini assemblies with bbmap
    report_and_log(('Round 2: mapping extracted reads to mini assemblies with bbmap'), pwd_log_file, keep_quiet)
    # bbmap_cmd_miniassembly_paired   = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, mini_assemblies, free_living_R1, free_living_R2, sam_file_mini_assembly_P, bbmap_parameter, sam_file_mini_assembly_stderr_P)
    # bbmap_cmd_miniassembly_unpaired = '%s ref=%s in=%s outm=%s %s 2> %s'        % (pwd_bbmap_exe, mini_assemblies, free_living_UP, sam_file_mini_assembly_UP, bbmap_parameter, sam_file_mini_assembly_stderr_UP)
    # os.system(bbmap_cmd_miniassembly_paired)
    # os.system(bbmap_cmd_miniassembly_unpaired)

    # # mapping with bbmap
    # if clp_read_for_assembly is True:
    #     bbmap_cmd_miniassembly_free_living_16s_paired   = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, mini_assemblies, free_living_16s_R1, free_living_16s_R2,  sam_file_mini_assembly_16s_P,  bbmap_parameter, sam_file_mini_assembly_16s_P_stderr)
    #     bbmap_cmd_miniassembly_free_living_16s_unpaired = '%s ref=%s in=%s outm=%s %s 2> %s'        % (pwd_bbmap_exe, mini_assemblies, free_living_16s_UP,                      sam_file_mini_assembly_16s_UP, bbmap_parameter, sam_file_mini_assembly_16s_UP_stderr)
    #     bbmap_cmd_miniassembly_free_living_ctg_paired   = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, mini_assemblies, free_living_ctg_R1, free_living_ctg_R2,  sam_file_mini_assembly_ctg_P,  bbmap_parameter, sam_file_mini_assembly_ctg_P_stderr)
    #     bbmap_cmd_miniassembly_free_living_ctg_unpaired = '%s ref=%s in=%s outm=%s %s 2> %s'        % (pwd_bbmap_exe, mini_assemblies, free_living_ctg_UP,                      sam_file_mini_assembly_ctg_UP, bbmap_parameter, sam_file_mini_assembly_ctg_UP_stderr)
    #     os.system(bbmap_cmd_miniassembly_free_living_16s_paired)
    #     os.system(bbmap_cmd_miniassembly_free_living_16s_unpaired)
    #     os.system(bbmap_cmd_miniassembly_free_living_ctg_paired)
    #     os.system(bbmap_cmd_miniassembly_free_living_ctg_unpaired)
    #     keep_best_matches_in_sam(sam_file_mini_assembly_16s_P, sam_file_mini_assembly_16s_P_best_match)
    #     keep_best_matches_in_sam(sam_file_mini_assembly_16s_UP, sam_file_mini_assembly_16s_UP_best_match)
    #     keep_best_matches_in_sam(sam_file_mini_assembly_ctg_P, sam_file_mini_assembly_ctg_P_best_match)
    #     keep_best_matches_in_sam(sam_file_mini_assembly_ctg_UP, sam_file_mini_assembly_ctg_UP_best_match)
    #     os.system('cat %s %s > %s' % (sam_file_mini_assembly_16s_P_best_match, sam_file_mini_assembly_16s_UP_best_match, sam_file_mini_assembly_16s_best_match))
    #     os.system('cat %s %s > %s' % (sam_file_mini_assembly_ctg_P_best_match, sam_file_mini_assembly_ctg_UP_best_match, sam_file_mini_assembly_ctg_best_match))
    # else:
    #     bbmap_cmd_miniassembly_free_living_16s = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, mini_assemblies, free_living_16s, sam_file_mini_assembly_16s, bbmap_parameter, sam_file_mini_assembly_16s_stderr)
    #     bbmap_cmd_miniassembly_free_living_ctg = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, mini_assemblies, free_living_ctg, sam_file_mini_assembly_ctg, bbmap_parameter, sam_file_mini_assembly_ctg_stderr)
    #     os.system(bbmap_cmd_miniassembly_free_living_16s)
    #     os.system(bbmap_cmd_miniassembly_free_living_ctg)
    #     keep_best_matches_in_sam(sam_file_mini_assembly_16s, sam_file_mini_assembly_16s_best_match)
    #     keep_best_matches_in_sam(sam_file_mini_assembly_ctg, sam_file_mini_assembly_ctg_best_match)

    # index miniassembly
    mini_assemblies_no_ext = '.'.join(mini_assemblies.split('.')[:-1])
    bowtie_build_mini_assemblies_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, mini_assemblies, mini_assemblies_no_ext)
    os.system(bowtie_build_mini_assemblies_cmd)

    # mapping with bowtie
    if clp_read_for_assembly is True:
        bowtie_cmd_miniassembly_free_living_16s = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p %s -f --local --all --no-unal 2> %s'    % (mini_assemblies_no_ext, free_living_16s_R1, free_living_16s_R2, free_living_16s_UP, sam_file_mini_assembly_16s, num_threads, sam_file_mini_assembly_16s_bowtie_log)
        bowtie_cmd_miniassembly_free_living_ctg = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p %s -f --local --all --no-unal 2> %s'    % (mini_assemblies_no_ext, free_living_ctg_R1, free_living_ctg_R2, free_living_ctg_UP, sam_file_mini_assembly_ctg, num_threads, sam_file_mini_assembly_ctg_bowtie_log)
    else:
        bowtie_cmd_miniassembly_free_living_16s = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --all --no-unal 2> %s'                % (mini_assemblies_no_ext, free_living_16s_UP, sam_file_mini_assembly_16s, num_threads, sam_file_mini_assembly_16s_bowtie_log)
        bowtie_cmd_miniassembly_free_living_ctg = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --all --no-unal 2> %s'                % (mini_assemblies_no_ext, free_living_ctg_UP, sam_file_mini_assembly_ctg, num_threads, sam_file_mini_assembly_ctg_bowtie_log)
    os.system(bowtie_cmd_miniassembly_free_living_16s)
    os.system(bowtie_cmd_miniassembly_free_living_ctg)

    # remove_both_ends_clp
    remove_both_ends_clp(sam_file_mini_assembly_16s, sam_file_mini_assembly_16s_single_clp)
    remove_both_ends_clp(sam_file_mini_assembly_ctg, sam_file_mini_assembly_ctg_single_clp)

    # reformat cigar string
    mini_assembly_16s_bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly_16s_single_clp, sam_file_mini_assembly_16s_single_clp_fmt, sam_file_mini_assembly_16s_single_clp_fmt_log)
    mini_assembly_ctg_bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly_ctg_single_clp, sam_file_mini_assembly_ctg_single_clp_fmt, sam_file_mini_assembly_ctg_single_clp_fmt_log)
    os.system(mini_assembly_16s_bbmap_reformat_cmd)
    os.system(mini_assembly_ctg_bbmap_reformat_cmd)

    # keep only best matches
    #keep_best_matches_in_sam(sam_file_mini_assembly_16s_single_clp_fmt, sam_file_mini_assembly_16s_best_match)
    #keep_best_matches_in_sam(sam_file_mini_assembly_ctg_single_clp_fmt, sam_file_mini_assembly_ctg_best_match)

    #keep_best_matches_in_sam_keep_short_M(sam_file_mini_assembly_16s_single_clp_fmt, min_M_len_ctg, sam_file_mini_assembly_16s_best_match)
    #keep_best_matches_in_sam_keep_short_M(sam_file_mini_assembly_ctg_single_clp_fmt, min_M_len_ctg, sam_file_mini_assembly_ctg_best_match)

    keep_best_matches_in_sam_drop_short_M(sam_file_mini_assembly_16s_single_clp_fmt, min_M_len_ctg, sam_file_mini_assembly_16s_best_match)
    keep_best_matches_in_sam_drop_short_M(sam_file_mini_assembly_ctg_single_clp_fmt, min_M_len_ctg, sam_file_mini_assembly_ctg_best_match)


    #################################################### parse sam file ####################################################

    # combine sam files
    # os.system('cat %s %s > %s' % (sam_file_mini_assembly_P, sam_file_mini_assembly_UP, sam_file_mini_assembly_combined))

    report_and_log(('Round 2: parsing sam file'), pwd_log_file, keep_quiet)

    # get gap seq to 16s reads dict
    gap_seq_to_16s_reads_dict = {}
    for each_read in open(sam_file_mini_assembly_16s_best_match):
        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            cigar = each_read_split[5]
            if cigar != '*':
                read_id = each_read_split[0]
                ref_id = each_read_split[2]
                cigar_splitted = cigar_splitter(cigar)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                if (aligned_len >= min_M_len_ctg) and (aligned_pct >= min_M_pct) and (mismatch_pct <= mismatch_cutoff):
                    if ref_id not in gap_seq_to_16s_reads_dict:
                        gap_seq_to_16s_reads_dict[ref_id] = [read_id]
                    else:
                        gap_seq_to_16s_reads_dict[ref_id].append(read_id)

    # get gap seq to ctg reads dict
    gap_seq_to_ctg_reads_dict = {}
    for each_read in open(sam_file_mini_assembly_ctg_best_match):
        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            cigar = each_read_split[5]
            if cigar != '*':
                read_id = each_read_split[0]
                ref_id = each_read_split[2]
                cigar_splitted = cigar_splitter(cigar)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                if (aligned_len >= min_M_len_ctg) and (aligned_pct >= min_M_pct) and (mismatch_pct <= mismatch_cutoff):
                    if ref_id not in gap_seq_to_ctg_reads_dict:
                        gap_seq_to_ctg_reads_dict[ref_id] = [read_id]
                    else:
                        gap_seq_to_ctg_reads_dict[ref_id].append(read_id)


    mini_assembly_to_16s_reads_handle = open(mini_assembly_to_16s_reads, 'w')
    for each_mini_assembly in gap_seq_to_16s_reads_dict:
        mini_assembly_to_16s_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join(gap_seq_to_16s_reads_dict[each_mini_assembly])))
    mini_assembly_to_16s_reads_handle.close()

    mini_assembly_to_ctg_reads_handle = open(mini_assembly_to_ctg_reads, 'w')
    for each_mini_assembly in gap_seq_to_ctg_reads_dict:
        mini_assembly_to_ctg_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join(gap_seq_to_ctg_reads_dict[each_mini_assembly])))
    mini_assembly_to_ctg_reads_handle.close()


    ############################################# get_GapFilling_stats #############################################

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


    ####################################################################################################################
    ####################################### combine linkages from step 1 and 2  ########################################
    ####################################################################################################################

    report_and_log(('Combining linkages from step 1 and 2'), pwd_log_file, keep_quiet)

    combined_linkage_file_handle     = open(combined_linkage_file, 'w')
    combined_linkage_file_tmp_handle = open(combined_linkage_file_tmp, 'w')
    combined_linkage_file_handle.write('MarkerGene\tGenomicSeq\tLinkage\tStep\n')
    combined_linkage_file_tmp_handle.write('MarkerGene,GenomicSeq,Number\n')
    for step_1_link in open(link_stats_combined_filtered_s1):
        if not step_1_link.startswith('MarkerGene,GenomicSeq,Number'):
            marker_id = step_1_link.strip().split(',')[0][12:]
            genome_id = step_1_link.strip().split(',')[1][12:]
            link_num  = step_1_link.strip().split(',')[2]
            combined_linkage_file_handle.write('%s\t%s\t%s\tS1\n' % (marker_id, genome_id, link_num))
            combined_linkage_file_tmp_handle.write(step_1_link)
    for step_2_link in open(stats_GapFilling_file_filtered):
        if not step_2_link.startswith('MarkerGene,GenomicSeq,Number'):
            marker_id = step_2_link.strip().split(',')[0][12:]
            genome_id = step_2_link.strip().split(',')[1][12:]
            link_num  = step_2_link.strip().split(',')[2]
            combined_linkage_file_handle.write('%s\t%s\t%s\tS2\n' % (marker_id, genome_id, link_num))
            combined_linkage_file_tmp_handle.write(step_2_link)
    combined_linkage_file_handle.close()
    combined_linkage_file_tmp_handle.close()

    #################### summarize linkages at contig level ####################

    combined_linkage_file_ctg_level_handle = open(combined_linkage_file_ctg_level, 'w')
    combined_linkage_file_ctg_level_handle.write('Marker___Genome(total)\tContig\tPaired\tClipping\tOverlapped\tStep\n')
    for each_linkage in open(combined_linkage_file):
        if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
            each_linkage_split = each_linkage.strip().split('\t')
            marker_id = each_linkage_split[0]
            mag_id = each_linkage_split[1]
            total_link_num = int(each_linkage_split[2])
            link_step = each_linkage_split[3]

            if link_step == 'S1':
                # first go through link num dict by paired reads
                counted_16s_to_ctg_key = set()
                for each_paired_link in marker_to_ctg_link_num_dict_pair:
                    if each_paired_link in marker_to_ctg_link_num_combined_to_use:
                        paired_link_16s_id = each_paired_link.split(marker_to_ctg_gnm_Key_connector)[0]
                        paired_link_ctg_id = each_paired_link.split(marker_to_ctg_gnm_Key_connector)[1]
                        paired_link_ctg_id_no_gnm = paired_link_ctg_id.split(gnm_to_ctg_connector)[1]
                        paired_link_gnm_id = paired_link_ctg_id.split(gnm_to_ctg_connector)[0]
                        if (paired_link_16s_id == marker_id) and (paired_link_gnm_id == mag_id):
                            current_pair_link_num = marker_to_ctg_link_num_dict_pair[each_paired_link]
                            current_clip_link_num = marker_to_ctg_link_num_dict_clip.get(each_paired_link, 0)
                            combined_linkage_file_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t%s\t0\tS1\n' % (paired_link_16s_id, paired_link_gnm_id, total_link_num, paired_link_ctg_id_no_gnm, current_pair_link_num, current_clip_link_num))
                            counted_16s_to_ctg_key.add(each_paired_link)

                # then go through link num dict by clipping reads
                for each_clip_link in marker_to_ctg_link_num_dict_clip:
                    if each_clip_link in marker_to_ctg_link_num_combined_to_use:
                        if each_clip_link not in counted_16s_to_ctg_key:
                            clip_link_16s_id = each_clip_link.split(marker_to_ctg_gnm_Key_connector)[0]
                            clip_link_ctg_id = each_clip_link.split(marker_to_ctg_gnm_Key_connector)[1]
                            clip_link_ctg_id_no_gnm = clip_link_ctg_id.split(gnm_to_ctg_connector)[1]
                            clip_link_gnm_id = clip_link_ctg_id.split(gnm_to_ctg_connector)[0]
                            if (clip_link_16s_id == marker_id) and (clip_link_gnm_id == mag_id):
                                current_pair_link_num = marker_to_ctg_link_num_dict_pair.get(each_clip_link, 0)
                                current_clip_link_num = marker_to_ctg_link_num_dict_clip[each_clip_link]
                                combined_linkage_file_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t%s\t0\tS1\n' % (clip_link_16s_id, clip_link_gnm_id, total_link_num, clip_link_ctg_id_no_gnm, current_pair_link_num, current_clip_link_num))
                                counted_16s_to_ctg_key.add(each_clip_link)

            if link_step == 'S2':
                for each_rd2_linkage in free_living_16s_to_ctg_linkage_dict_to_use:
                    rd2_link_16s_id = each_rd2_linkage.split(marker_to_ctg_gnm_Key_connector)[0]
                    rd2_link_ctg_id = each_rd2_linkage.split(marker_to_ctg_gnm_Key_connector)[1]
                    rd2_link_ctg_id_no_gnm = rd2_link_ctg_id.split(gnm_to_ctg_connector)[1]
                    rd2_link_gnm_id = rd2_link_ctg_id.split(gnm_to_ctg_connector)[0]
                    if (rd2_link_16s_id == marker_id) and (rd2_link_gnm_id == mag_id):
                        current_ctg_link_num = free_living_16s_to_ctg_linkage_dict_to_use[each_rd2_linkage]
                        combined_linkage_file_ctg_level_handle.write('%s___%s(%s)\t%s\t0\t0\t%s\tS2\n' % (rd2_link_16s_id, rd2_link_gnm_id, total_link_num, rd2_link_ctg_id_no_gnm, current_ctg_link_num))
    combined_linkage_file_ctg_level_handle.close()


    ####################################################################################################################
    ####################################################### plot #######################################################
    ####################################################################################################################

    report_and_log(('Visualising linkages'), pwd_log_file, keep_quiet)
    sankey_linkages(combined_linkage_file_ctg_level, linkage_plot_rd1_html, linkage_plot_rd2_html)


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
        recovery_combined, accuracy_combined, recovered_combined = get_accuracy(combined_linkage_file_tmp, len(marker_id_set))

        # get unrecovered markers
        unrecovered_markers_paired = get_unrecovered_markers(marker_id_set, recovered_combined)
        unrecovered_markers_paired_str = 'Unrecovered(%s):%s' % (len(unrecovered_markers_paired), ','.join(sorted([i for i in unrecovered_markers_paired])))

        # assessment by genome
        assign_rate, assign_accuracy, right_assign, wrong_assign = get_accuracy_by_genome(combined_linkage_file_tmp, mag_folder, mag_file_extension)
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
    #link_16s_parser_input_files.add_argument('-g',          required=False, metavar='',             default=None,           help='genomic sequences')
    link_16s_parser_input_files.add_argument('-mag',        required=False, metavar='',             default=None,           help='metagenome-assembled-genome (MAG) folder')
    link_16s_parser_input_files.add_argument('-x',          required=False, metavar='',             default='fasta',        help='MAG file extension, (default: %(default)s)')
    link_16s_parser_input_files.add_argument('-depth',      required=False, metavar='', type=float, default=0,              help='minimum depth multiple between 16S and  genomic sequences, a value of no higher than 0.2 is recommended, (default: %(default)s)')
    link_16s_parser_input_files.add_argument('-r16s',       required=False, metavar='',                                     help='16S reads')

    # 16S rRNA gene related parameters
    link_16s_parser_16s.add_argument('-min_iden_16s',       required=False, metavar='', type=float, default=98,             help='minimum similarity for 16S sequences to be assigned to the same genome, (default: %(default)s)')
    link_16s_parser_16s.add_argument('-min_cov_16s',        required=False, metavar='', type=float, default=30,             help='coverage cutoff for calculating pairwise 16S similarity, (default: %(default)s)')
    link_16s_parser_16s.add_argument('-min_aln_16s',        required=False, metavar='', type=int,   default=500,            help='alignment length cutoff for calculating pairwise 16S similarity, (default: %(default)s)')

    # parameters for both rounds linking
    link_16s_parser_both_rds.add_argument('-mismatch',      required=False, metavar='', type=float, default=1,              help='maximum mismatch percentage, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_len_16s', required=False, metavar='', type=int,   default=75,             help='minimum length aligned to 16S, (default: %(default)s)')
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
    link_16s_parser_rd2.add_argument('-mira_tmp',           required=False, default=None,                                   help='tmp dir for mira')
    link_16s_parser_rd2.add_argument('-spades',             required=False, action="store_true",                            help='run spades, instead of Mira')
    link_16s_parser_rd2.add_argument('-link_bias_rd2',      required=False, metavar='', type=float, default=40,             help='max_mini_assembly_link_num_diff_between_ctg_16s, (default: %(default)s)')

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

    args = vars(link_16s_parser.parse_args())
    link_16s(args)

'''
1. how to incorporate the taxonomy of MAGs and 16S sequences
2. the depth of 16S sequences always not lower than the genome they come from
3. with no_ambiguous option, 16S rRNA gene sequences need to be dereplicated. (include dereplication step? with identity and coverage cutoffs?)
4. (doesn't work)!!! to ignore list even without assignment (to handle situations like DM_m4, meanwhile capicable of not assign very diverde 16S (e.g. <98% identity) to the same genome)
5. check the structure of assembled sequences? v1, 2, 3 or v4, 5, 6? how?
6. add "16S_reads" to SortMeRNA's output prefix
7. estimate cutoffs to use based sensitive or specific
8. add --careful to spades?
9. for clipping mapped reads, the mismatch of clipping part must be 0
10. insert size is important
'''
