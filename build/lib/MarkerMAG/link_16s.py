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
import pandas as pd
from Bio import SeqIO
from time import sleep
import multiprocessing as mp
from datetime import datetime
from distutils.spawn import find_executable


config_dict = {'config_file_path'   : '/'.join(os.path.realpath(__file__).split('/')[:-1]),
               'bowtie2'            : 'bowtie2',
               'bowtie2_build'      : 'bowtie2-build',
               'samtools'           : 'samtools',
               'blastn'             : 'blastn',
               'makeblastdb'        : 'makeblastdb',
               'spades'             : 'spades.py',
               'get_sankey_plot_R'  : '%s/get_sankey_plot.R' % '/'.join(os.path.realpath(__file__).split('/')[:-1])}


link_Marker_MAG_usage = '''
=================================== MarkerMAG example commands ===================================

MarkerMAG -p Test -r1 R1.fasta -r2 R2.fasta -m 16S_seqs.fa -g contig.fasta -t 4
MarkerMAG -p Test -r1 R1.fasta -r2 R2.fasta -m 16S_seqs.fa -mag MAGs -x fa -t 4

Note:
MarkerMAG assumes the id of paired reads in a format of XXXX.1 and XXXX.2. The only difference 
is the last character. You can rename your reads with "Rename_reads".

==================================================================================================
'''


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
    cmd_bowtie2_build   = 'bowtie2-build -f %s %s/%s --quiet --threads %s' % (ref_seq, ref_seq_file_path, ref_seq_file_basename, num_threads)
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
        if each_element.isalpha() is True:
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


def extract_reads_worker(argument_list):

    reads_file_in    = argument_list[0]
    reads_to_extract = argument_list[1]
    reads_file_out   = argument_list[2]

    reads_file_out_handle = open(reads_file_out, 'w')
    for read_record in SeqIO.parse(reads_file_in, 'fasta'):
        if read_record.id in reads_to_extract:
            reads_file_out_handle.write('>%s\n' % read_record.id)
            reads_file_out_handle.write('%s\n' % read_record.seq)
    reads_file_out_handle.close()


def extracted_reads_with_multiprocessing(reads_r1, reads_r2, r1_to_extract, r2_to_extract, output_folder, num_threads):
    solely_perfectly_mapped_reads_r1_splitted = split_list(r1_to_extract, num_threads // 2)
    solely_perfectly_mapped_reads_r2_splitted = split_list(r2_to_extract, num_threads // 2)

    argument_list_for_extract_reads_worker = []
    extract_reads_file_index_r1 = 1
    for reads_subset_r1 in solely_perfectly_mapped_reads_r1_splitted:
        current_output_file = '%s/extract_r1_subset_%s.fasta' % (output_folder, extract_reads_file_index_r1)
        argument_list_for_extract_reads_worker.append([reads_r1, reads_subset_r1, current_output_file])
        extract_reads_file_index_r1 += 1

    extract_reads_file_index_r2 = 1
    for reads_subset_r2 in solely_perfectly_mapped_reads_r2_splitted:
        current_output_file = '%s/extract_r2_subset_%s.fasta' % (output_folder, extract_reads_file_index_r2)
        argument_list_for_extract_reads_worker.append([reads_r2, reads_subset_r2, current_output_file])
        extract_reads_file_index_r2 += 1

    # extract reads with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(extract_reads_worker, argument_list_for_extract_reads_worker)
    pool.close()
    pool.join()


def blast_results_to_dict(blastn_results, iden_cutoff, query_cov_cutoff):
    query_to_subject_list_dict = {}

    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        subject = blast_hit_split[1]
        subject_with_prefix = 'GenomicSeq__%s' % subject
        iden = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_len = int(blast_hit_split[12])
        coverage_q = float(align_len) * 100 / float(query_len)
        if (iden >= iden_cutoff) and (coverage_q >= query_cov_cutoff):
            if query not in query_to_subject_list_dict:
                query_to_subject_list_dict[query] = [subject_with_prefix]
            else:
                query_to_subject_list_dict[query].append(subject_with_prefix)

    return query_to_subject_list_dict


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


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, file_out):

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

            if linkage_num >= min_linkages:

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


def get_free_living_mate(ref_in, reads_r1, reads_r2, end_seq_len, minCigarM, max_gap_to_end, bowtie_build_exe, bowtie2_exe, num_threads):

    ref_in_path, ref_in_basename, ref_in_ext = sep_path_basename_ext(ref_in)

    ref_subset = '%s/%s_ends_%sbp%s' % (ref_in_path, ref_in_basename, end_seq_len, ref_in_ext)
    sam_file = '%s/%s_ends_%sbp.sam' % (ref_in_path, ref_in_basename, end_seq_len)

    # get ref seqs subset
    ref_subset_len_dict = {}
    ref_subset_handle = open(ref_subset, 'w')
    for ref_seq in SeqIO.parse(ref_in, 'fasta'):
        ref_seq_id = ref_seq.id
        ref_seq_len = len(ref_seq.seq)
        if ref_seq_len < end_seq_len * 3:
            ref_subset_handle.write('>%s\n' % ref_seq_id)
            ref_subset_handle.write('%s\n' % ref_seq.seq)
            ref_subset_len_dict[ref_seq_id] = ref_seq_len
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

    # index ref seq subset
    ref_subset_no_ext = '.'.join(ref_subset.split('.')[:-1])
    bowtie2_index_ref_cmd = '%s -f %s %s --quiet --threads %s' % (
    bowtie_build_exe, ref_subset, ref_subset_no_ext, num_threads)
    os.system(bowtie2_index_ref_cmd)

    # mapping
    bowtie2_mapping_cmd = '%s -x %s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s' % (
    bowtie2_exe, ref_subset_no_ext, reads_r1, reads_r2, sam_file, num_threads)
    os.system(bowtie2_mapping_cmd)

    # parse sam file
    all_mapped_reads = set()
    qualified_reads_dict = {}
    qualified_reads_to_ref_dict = {}
    for each_line in open(sam_file):
        if not each_line.startswith('@'):
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_line_split[2]
            ref_pos = int(each_line_split[3])
            cigar = each_line_split[5]
            cigar_splitted = cigar_splitter(cigar)

            all_mapped_reads.add(read_id)

            qualified_mapping = False
            if ref_id[-2:] == '_l':
                if ref_pos <= max_gap_to_end:
                    if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                        qualified_mapping = True

                    elif (len(cigar_splitted) == 2) and (cigar_splitted[-1][-1] == 'M') and (
                            cigar_splitted[0][-1] == 'S') and (int(cigar_splitted[-1][:-1]) >= minCigarM) and (
                            ref_pos == 1):
                        qualified_mapping = True

            elif ref_id[-2:] == '_r':

                if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                    ref_pos_end = ref_pos + int(cigar[:-1])

                    if (end_seq_len - ref_pos_end) <= max_gap_to_end:
                        qualified_mapping = True

                elif (len(cigar_splitted) == 2) and (cigar_splitted[0][-1] == 'M') and (
                        cigar_splitted[1][-1] == 'S') and (int(cigar_splitted[0][:-1]) >= minCigarM):
                    if (ref_pos + int(cigar_splitted[0][:-1]) - 1) == end_seq_len:
                        qualified_mapping = True

            else:
                ref_len = ref_subset_len_dict[ref_id]
                if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                    ref_pos_end = ref_pos + int(cigar[:-1])

                    # left side
                    if ref_pos <= max_gap_to_end:
                        qualified_mapping = True

                    # right side
                    elif (ref_len - ref_pos_end) <= max_gap_to_end:
                        qualified_mapping = True

                if len(cigar_splitted) == 2:

                    # left side
                    if (cigar_splitted[-1][-1] == 'M') and (cigar_splitted[0][-1] == 'S') and (
                            int(cigar_splitted[-1][:-1]) >= minCigarM) and (ref_pos == 1):
                        qualified_mapping = True

                    # right side
                    elif (cigar_splitted[-0][-1] == 'M') and (cigar_splitted[1][-1] == 'S') and (
                            int(cigar_splitted[0][:-1]) >= minCigarM):
                        if (ref_pos + int(cigar_splitted[0][:-1]) - 1) == ref_len:
                            qualified_mapping = True

            if qualified_mapping is True:
                qualified_reads_to_ref_dict[read_id] = ref_id

                if read_id_base not in qualified_reads_dict:
                    qualified_reads_dict[read_id_base] = [read_strand]
                else:
                    qualified_reads_dict[read_id_base].append(read_strand)

    reads_to_extract_to_ref_dict = {}
    for qualified_read in qualified_reads_dict:
        read_strand = qualified_reads_dict[qualified_read]
        if len(read_strand) == 1:

            mapped_mate = ''
            mate_to_extract = ''
            if read_strand == ['1']:
                mapped_mate = '%s.1' % (qualified_read)
                mate_to_extract = '%s.2' % (qualified_read)
            if read_strand == ['2']:
                mapped_mate = '%s.2' % (qualified_read)
                mate_to_extract = '%s.1' % (qualified_read)

            if mate_to_extract not in all_mapped_reads:
                reads_to_extract_to_ref_dict[mate_to_extract] = qualified_reads_to_ref_dict[mapped_mate]

    return reads_to_extract_to_ref_dict


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


def link_16s(args, config_dict):

    ###################################################### file in/out #####################################################

    # file in
    output_prefix                       = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    reads_file_16s                      = args['r16s']
    genomic_assemblies                  = args['g']
    mag_folder                          = args['mag']
    mag_file_extension                  = args['x']
    marker_gene_seqs                    = args['m']
    min_16s_gnm_multiple                = args['depth']
    min_cigar_M                         = args['s1_cigarM']
    min_cigar_S                         = args['s1_cigarS']
    reads_iden_cutoff                   = args['s1_ri']
    reads_cov_cutoff                    = args['s1_rc']
    within_genome_minimum_iden16s       = args['s1_mi']
    cov16s                              = args['s1_mc']
    aln16s                              = args['s1_ma']
    min_paired_linkages                 = args['s1_mpl']
    end_seq_len                         = args['s2_e']
    minCigarM                           = args['s2_m']
    max_gap_to_end                      = args['s2_g']
    min_read_num                        = args['s2_r']
    num_threads                         = args['t']
    keep_quiet                          = args['quiet']
    force_overwrite                     = args['force']
    keep_temp                           = args['tmp']
    test_mode                           = args['test_mode']

    pwd_plot_sankey_R                   = config_dict['get_sankey_plot_R']
    pwd_makeblastdb_exe                 = config_dict['makeblastdb']
    pwd_blastn_exe                      = config_dict['blastn']
    pwd_bowtie2_build_exe               = config_dict['bowtie2_build']
    pwd_bowtie2_exe                     = config_dict['bowtie2']
    pwd_samtools_exe                    = config_dict['samtools']
    pwd_spades_exe                      = config_dict['spades']

    str_connector = '___'


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


    ######################## check genomic sequence type and prepare files for making blast db #########################

    blast_db           = ''
    genomic_seq_type   = ''  # ctg or mag
    renamed_mag_folder = ''

    # check the type of input genomic sequences
    if (genomic_assemblies is not None) and (mag_folder is None):
        genomic_seq_type = 'ctg'
        metagenomic_assemblies_file_path, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension = sep_path_basename_ext(genomic_assemblies)
        blast_db_dir = '%s/%s_%s_db' % (step_1_wd, output_prefix, metagenomic_assemblies_file_basename)
        blast_db     = '%s/%s%s'     % (blast_db_dir, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension)
        os.mkdir(blast_db_dir)
        os.system('cp %s %s/' % (genomic_assemblies, blast_db_dir))

    elif (genomic_assemblies is None) and (mag_folder is not None):
        genomic_seq_type = 'mag'
        mag_folder_name     = mag_folder.split('/')[-1]
        blast_db_dir        = '%s/%s_db'            % (step_1_wd, mag_folder_name)
        renamed_mag_folder  = '%s/%s_db/%s_renamed' % (step_1_wd, mag_folder_name, mag_folder_name)
        os.mkdir(blast_db_dir)
        os.mkdir(renamed_mag_folder)

        # get input mag file list
        mag_file_re = '%s/*%s' % (mag_folder, mag_file_extension)
        mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
        if len(mag_file_list) == 0:
            print('No MAG detected, program exited!')
            exit()

        # add mag id to its sequences
        for mag_in in mag_file_list:
            pwd_mag_in      = '%s/%s' % (mag_folder, mag_in)
            pwd_mag_renamed = '%s/%s' % (renamed_mag_folder, mag_in)
            mag_basename    = '.'.join(mag_in.split('.')[:-1])
            rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename, str_connector)

        # combine renamed MAGs
        blast_db = '%s/%s_combined.fa' % (blast_db_dir, mag_folder_name)
        os.system('cat %s/*%s > %s' % (renamed_mag_folder, mag_file_extension, blast_db))

    else:
        print('Please provide genomic sequences either as raw assemblies (-g) or as MAGs (-mag)')
        exit()


    ########################################### define folder and file name ############################################

    marker_gene_seqs_file_path, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension = sep_path_basename_ext(marker_gene_seqs)

    pwd_log_file                                = '%s/%s.log'                                    % (working_directory, output_prefix)

    bowtie_index_dir                            = '%s/%s_index'                                  % (step_1_wd, marker_gene_seqs_file_basename)
    pwd_samfile                                 = '%s/%s.sam'                                    % (step_1_wd, marker_gene_seqs_file_basename)
    clipping_reads_matched_part                 = '%s/clipping_matched_part.txt'                 % step_1_wd
    clipping_reads_not_matched_part_seq         = '%s/clipping_not_matched_part_seq.fasta'       % step_1_wd
    clipping_reads_not_matched_part_seq_blastn  = '%s/clipping_not_matched_part_seq_blast.txt'   % step_1_wd
    clipping_reads_match_profile                = '%s/match_profile_clipping.txt'                % step_1_wd
    unmapped_paired_reads_folder                = '%s/unmapped_paired_reads'                     % step_1_wd
    unmapped_paired_reads_file                  = '%s/unmapped_paired_reads.fasta'               % step_1_wd
    unmapped_paired_reads_blastn                = '%s/unmapped_paired_reads_blast.txt'           % step_1_wd
    paired_reads_match_profile                  = '%s/match_profile_paired.txt'                  % step_1_wd
    blast_results_all_vs_all_16s                = '%s/16S_all_vs_all_blastn.tab'                 % step_1_wd
    link_stats_clipping                         = '%s/stats_clipping.txt'                        % step_1_wd
    link_stats_clipping_filtered                = '%s/stats_clipping_filtered.txt'               % step_1_wd
    link_stats_paired                           = '%s/stats_paired.txt'                          % step_1_wd
    link_stats_paired_filtered                  = '%s/stats_paired_filtered.txt'                 % step_1_wd
    link_stats_combined_table                   = '%s/stats_combined_table.txt'                  % step_1_wd
    link_stats_combined                         = '%s/stats_combined.txt'                        % step_1_wd
    pairwise_marker_similarity                  = '%s/pairwise_marker_similarity.txt'            % step_1_wd
    depth_file_ctg                              = '%s/mean_depth_ctg.txt'                        % step_1_wd
    depth_file_gnm                              = '%s/mean_depth_gnm.txt'                        % step_1_wd
    depth_file_16s                              = '%s/mean_depth_16s.txt'                        % step_1_wd

    blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % num_threads


    ################################################# step 2 #################################################

    marker_gene_seqs_1st_round_unlinked         = '%s/step_1_unlinked_marker_genes.fasta'        % step_2_wd
    combined_1st_round_unlinked_mags            = '%s/step_1_unlinked_combined_gnms.fasta'       % step_2_wd
    combined_1st_round_unlinked_ctgs            = '%s/step_1_unlinked_combined_ctgs.fasta'       % step_2_wd
    free_living_mate_gnm                        = '%s/free_living_mate_ctg.txt'                  % step_2_wd
    free_living_mate_16s                        = '%s/free_living_mate_16s.txt'                  % step_2_wd
    extracted_reads_folder                      = '%s/free_living_reads'                         % step_2_wd
    extracted_reads_cbd                         = '%s/free_living_read_combined.fasta'           % step_2_wd
    spades_wd                                   = '%s/combined_free_living_reads_SPAdes_wd'      % step_2_wd
    spades_assemblies                           = '%s/scaffolds.fasta'                           % spades_wd
    sam_file                                    = '%s/scaffolds.sam'                             % step_2_wd
    stats_GapFilling_file_16s                   = '%s/stats_GapFilling_16s.txt'                  % step_2_wd
    stats_GapFilling_file_ctg                   = '%s/stats_GapFilling_ctg.txt'                  % step_2_wd
    stats_GapFilling_file_filtered_16s          = '%s/stats_GapFilling_16s_filtered.txt'         % step_2_wd
    stats_GapFilling_file_filtered_ctg          = '%s/stats_GapFilling_ctg_filtered.txt'         % step_2_wd
    stats_GapFilling_file                       = '%s/stats_GapFilling.txt'                      % step_2_wd
    stats_GapFilling_file_filtered              = '%s/stats_GapFilling_filtered.txt'             % step_2_wd
    spades_log                                  = '%s/SPAdes.log'                                % step_2_wd
    combined_linkage_file_tmp                   = '%s/combined_linkages_tmp.txt'                 % step_2_wd
    combined_linkage_file_tmp_html              = '%s/combined_linkages_tmp.html'                % step_2_wd
    combined_linkage_file                       = '%s/%s_combined_linkages.txt'                  % (working_directory, output_prefix)
    combined_linkage_file_html                  = '%s/%s_combined_linkages.html'                 % (working_directory, output_prefix)


    #################################### calculate mean depth for genome/assemblies ####################################

    mean_depth_dict_ctg = {}
    mean_depth_dict_gnm = {}
    if min_16s_gnm_multiple > 0:

        if genomic_seq_type == 'ctg':
            report_and_log(('Step 1: calculating depth for %s' % genomic_assemblies), pwd_log_file, keep_quiet)
        if genomic_seq_type == 'mag':
            report_and_log(('Step 1: calculating depth for genomes in %s' % mag_folder), pwd_log_file, keep_quiet)

        # get mean depth for contig
        mean_depth_dict_ctg, ctg_len_dict = get_ctg_mean_depth_by_samtools_coverage(True, blast_db, reads_file_r1, reads_file_r2, '', num_threads)

        # write out ctg depth
        depth_file_ctg_handle = open(depth_file_ctg, 'w')
        for ctg in mean_depth_dict_ctg:
            depth_file_ctg_handle.write('%s\t%s\n' % (ctg, mean_depth_dict_ctg[ctg]))
        depth_file_ctg_handle.close()

        # get mean_depth_dict_gnm
        if genomic_seq_type == 'mag':

            gnm_ctg_connector = '___'

            gnm_len_total_depth_dict = {}
            for ctg in mean_depth_dict_ctg:
                ctg_genome = ctg.split(gnm_ctg_connector)[0]
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

        report_and_log(('Step 1: calculating depth for %s' % marker_gene_seqs), pwd_log_file, keep_quiet)

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

    # copy marker gene sequence file to index folder
    report_and_log(('Step 1: Indexing marker gene sequences for mapping'), pwd_log_file, keep_quiet)
    bowtie2_index_ref_cmd = '%s -f %s/%s%s %s/%s --quiet --threads %s' % (pwd_bowtie2_build_exe, bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension, bowtie_index_dir, marker_gene_seqs_file_basename, num_threads)
    os.system(bowtie2_index_ref_cmd)

    # mapping
    report_and_log(('Step 1: Mapping reads to marker gene sequences'), pwd_log_file, keep_quiet)
    sleep(1)
    report_and_log(('Step 1: Please ignore warnings starting with "Use of uninitialized value" during Bowtie mapping.'), pwd_log_file, keep_quiet)
    sleep(1)
    bowtie2_mapping_cmd = '%s -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s' % (pwd_bowtie2_exe, bowtie_index_dir, marker_gene_seqs_file_basename, reads_file_r1, reads_file_r2, pwd_samfile, num_threads)
    os.system(bowtie2_mapping_cmd)
    sleep(1)
    report_and_log(('Step 1: Please ignore warnings starting with "Use of uninitialized value" during Bowtie mapping.'), pwd_log_file, keep_quiet)
    sleep(1)


    ##################################################### extract reads ####################################################

    report_and_log(('Step 1: Extracting unmapped part of clipping mapped reads from sam file'), pwd_log_file, keep_quiet)

    # export clipping mapped reads and perfectly mapped reads
    clipping_mapped_reads_list = set()
    clipping_reads_mapped_part_dict = {}
    perfectly_mapped_reads_dict = {}
    clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
    for each_read in open(pwd_samfile):

        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
            ref_pos = int(each_read_split[3])
            cigar = each_read_split[5]
            read_seq = each_read_split[9]
            cigar_splitted  = cigar_splitter(cigar)
            read_id_with_ref_pos = '%s__x__%s__x__%s' % (read_id, ref_id, ref_pos)

            # for perfectly mapped reads
            if ('M' in cigar) and (len(cigar_splitted) == 1):
                if read_id_base not in perfectly_mapped_reads_dict:
                    perfectly_mapped_reads_dict[read_id_base] = {read_strand:[ref_id_with_prefix]}
                else:
                    if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                        perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                    else:
                        perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

            # for clipping mapped reads
            if ('S' in cigar) and (len(cigar_splitted) == 2):
                cigar_M_len = 0
                cigar_S_len = 0
                split_pos = 0
                if cigar_splitted[0][-1] == 'M':
                    cigar_M_len = int(cigar_splitted[0][:-1])
                    cigar_S_len = int(cigar_splitted[1][:-1])
                    split_pos = ref_pos + cigar_M_len
                if cigar_splitted[1][-1] == 'M':
                    cigar_M_len = int(cigar_splitted[1][:-1])
                    cigar_S_len = int(cigar_splitted[0][:-1])
                    split_pos = ref_pos

                if (cigar_M_len >= min_cigar_M) and (cigar_S_len >= min_cigar_S):
                    read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                    read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                    if cigar_splitted[0][-1] == 'M':

                        # write out the sequence of unmapped part
                        clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id_with_ref_pos)
                        clipping_reads_not_matched_part_seq_handle.write(read_seq_right + '\n')

                        # store the match info of mapped part
                        if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                            clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                        else:
                            clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                    if cigar_splitted[1][-1] == 'M':

                        # write out the sequence of unmapped part
                        clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id_with_ref_pos)
                        clipping_reads_not_matched_part_seq_handle.write(read_seq_left + '\n')

                        # store the match info of mapped part
                        if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                            clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                        else:
                            clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                    clipping_mapped_reads_list.add(read_id)
    clipping_reads_not_matched_part_seq_handle.close()


    ########################################## extract reads with multiprocessing ##########################################

    perfectly_mapped_read_singleton_dict = {}
    for perfectly_mapped_read in perfectly_mapped_reads_dict:
        current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
        if len(current_value) == 1:
            perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value

    # get the id of paired reads to extract
    solely_perfectly_mapped_reads_r1 = set()
    solely_perfectly_mapped_reads_r2 = set()
    for perfectly_mapped_read in perfectly_mapped_read_singleton_dict:
        current_value = perfectly_mapped_read_singleton_dict[perfectly_mapped_read]
        strand = list(current_value.keys())[0]
        if strand == '1':
            r2_to_extract = '%s.2' % perfectly_mapped_read
            if r2_to_extract not in clipping_mapped_reads_list:
                solely_perfectly_mapped_reads_r2.add(r2_to_extract)
        if strand == '2':
            r1_to_extract = '%s.1' % perfectly_mapped_read
            if r1_to_extract not in clipping_mapped_reads_list:
                solely_perfectly_mapped_reads_r1.add(r1_to_extract)

    # extract reads
    report_and_log(('Step 1: Extracting unmapped paired reads with %s cores' % num_threads), pwd_log_file, keep_quiet)

    os.mkdir(unmapped_paired_reads_folder)

    # extract reads with multiprocessing
    extracted_reads_with_multiprocessing(reads_file_r1, reads_file_r2, solely_perfectly_mapped_reads_r1, solely_perfectly_mapped_reads_r2, unmapped_paired_reads_folder, num_threads)

    # combine extracted reads
    os.system('cat %s/*.fasta > %s' % (unmapped_paired_reads_folder, unmapped_paired_reads_file))


    ############################# run blast between extracted reads and metagenomic assemblies #############################

    # run blastn
    makeblastdb_cmd     = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blast_db)
    blastn_cmd_paired   = '%s -query %s -db %s -out %s %s'                          % (pwd_blastn_exe, unmapped_paired_reads_file, blast_db, unmapped_paired_reads_blastn, blast_parameters)
    blastn_cmd_clipping = '%s -query %s -db %s -out %s %s'                          % (pwd_blastn_exe, clipping_reads_not_matched_part_seq, blast_db, clipping_reads_not_matched_part_seq_blastn, blast_parameters)

    report_and_log(('Step 1: Making blastn database'), pwd_log_file, keep_quiet)
    os.system(makeblastdb_cmd)

    report_and_log(('Step 1: Running blastn for unmapped paired reads'), pwd_log_file, keep_quiet)
    os.system(blastn_cmd_paired)

    report_and_log(('Step 1: Running blastn for unmapped parts of clipping mapped reads'), pwd_log_file, keep_quiet)
    os.system(blastn_cmd_clipping)


    ######################################### parse blast results for paired reads #########################################

    report_and_log(('Step 1: Parsing blast results for paired reads'), pwd_log_file, keep_quiet)

    # filter blast results for paired reads
    unmapped_paired_reads_to_ctg_dict = blast_results_to_dict(unmapped_paired_reads_blastn, reads_iden_cutoff, reads_cov_cutoff)

    paired_stats_dict_num             = {}
    paired_reads_match_profile_handle = open(paired_reads_match_profile, 'w')
    paired_reads_match_profile_handle.write('ID\tR1\tR2\n')
    for unmapped_read in unmapped_paired_reads_to_ctg_dict:

        unmapped_read_base = '.'.join(unmapped_read.split('.')[:-1])
        unmapped_read_strand = unmapped_read.split('.')[-1]

        unmapped_read_base_r1_matched = []
        unmapped_read_base_r2_matched = []
        if unmapped_read_strand == '1':
            unmapped_read_base_r1_matched = unmapped_paired_reads_to_ctg_dict[unmapped_read]
            if unmapped_read_base in perfectly_mapped_read_singleton_dict:
                unmapped_read_base_r2_matched = perfectly_mapped_read_singleton_dict[unmapped_read_base]['2']
        if unmapped_read_strand == '2':
            unmapped_read_base_r2_matched = unmapped_paired_reads_to_ctg_dict[unmapped_read]
            if unmapped_read_base in perfectly_mapped_read_singleton_dict:
                unmapped_read_base_r1_matched = perfectly_mapped_read_singleton_dict[unmapped_read_base]['1']

        if (unmapped_read_base_r1_matched != []) and (unmapped_read_base_r2_matched != []):
            for r1 in unmapped_read_base_r1_matched:
                for r2 in unmapped_read_base_r2_matched:

                    # write out to file
                    paired_reads_match_profile_handle.write('%s\t%s\t%s\n' % (unmapped_read_base, r1, r2))

                    # store in dict
                    paired_key = '_|_'.join(sorted([r1, r2])[::-1])
                    if genomic_seq_type == 'mag':
                        paired_key = '___'.join(paired_key.split('___')[:-1])
                    if paired_key not in paired_stats_dict_num:
                        paired_stats_dict_num[paired_key] = 1
                    else:
                        paired_stats_dict_num[paired_key] += 1

    paired_reads_match_profile_handle.close()


    #################################### parse blast results for clipping mapped reads #####################################

    report_and_log(('Step 1: Parsing blast results for clipping mapped reads'), pwd_log_file, keep_quiet)

    # filter blast results for clipping mapped reads
    clipping_parts_to_ctg_dict  = blast_results_to_dict(clipping_reads_not_matched_part_seq_blastn, reads_iden_cutoff, reads_cov_cutoff)
    clipping_reads_match_profile_handle = open(clipping_reads_match_profile, 'w')
    clipping_reads_match_profile_handle.write('ID\tLeft\tRight\n')
    clipping_stats_dict_num = {}
    for clipping_mapped_read in clipping_reads_mapped_part_dict:

        #clipping_mapped_read_id = clipping_mapped_read[:-2]
        clipping_mapped_read_id = clipping_mapped_read.split('__x__')[0]
        mapped_part = clipping_mapped_read[-1]

        clipping_mapped_read_matches_l = []
        clipping_mapped_read_matches_r = []
        if mapped_part == 'l':
            clipping_mapped_read_matches_l = clipping_reads_mapped_part_dict[clipping_mapped_read]
            if ('%s_r' % clipping_mapped_read[:-2]) in clipping_parts_to_ctg_dict:
                clipping_mapped_read_matches_r = clipping_parts_to_ctg_dict[('%s_r' % clipping_mapped_read[:-2])]
        if mapped_part == 'r':
            clipping_mapped_read_matches_r = clipping_reads_mapped_part_dict[clipping_mapped_read]
            if ('%s_l' % clipping_mapped_read[:-2]) in clipping_parts_to_ctg_dict:
                clipping_mapped_read_matches_l = clipping_parts_to_ctg_dict[('%s_l' % clipping_mapped_read[:-2])]

        if (clipping_mapped_read_matches_l != []) and (clipping_mapped_read_matches_r != []):

            for l in clipping_mapped_read_matches_l:
                for r in clipping_mapped_read_matches_r:

                    # write out to file
                    clipping_reads_match_profile_handle.write('%s\t%s\t%s\n' % (clipping_mapped_read_id, l, r))

                    # store in dict
                    clipping_key = '_|_'.join(sorted([l, r])[::-1])
                    if genomic_seq_type == 'mag':
                        clipping_key = '___'.join(clipping_key.split('___')[:-1])

                    if clipping_key not in clipping_stats_dict_num:
                        clipping_stats_dict_num[clipping_key] = 1
                    else:
                        clipping_stats_dict_num[clipping_key] += 1

    clipping_reads_match_profile_handle.close()


    ############################################## get pairwise_16s_iden_dict ##############################################

    # makeblastdn with marker gene sequences
    blastdb_16s         = '%s/%s%s' % (bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension)
    makeblastdb_16s_cmd = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blastdb_16s)
    os.system(makeblastdb_16s_cmd)

    all_vs_all_16s_blastn_cmd = '%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, blastdb_16s, blastdb_16s, blast_results_all_vs_all_16s, blast_parameters)
    os.system(all_vs_all_16s_blastn_cmd)

    pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, aln16s, cov16s)

    # write out to file
    pairwise_marker_similarity_handle = open(pairwise_marker_similarity, 'w')
    pairwise_marker_similarity_handle.write('Marker1\tMarker2\tSimilarity\n')
    for marker_pair in pairwise_16s_iden_dict:
        pairwise_marker_similarity_handle.write('%s\t%s\n' % ('\t'.join(marker_pair.split('__|__')), pairwise_16s_iden_dict[marker_pair]))
    pairwise_marker_similarity_handle.close()


    ##################################################### get linkages #####################################################

    report_and_log(('Step 1: Parsing linkages'), pwd_log_file, keep_quiet)

    # prepare input file for sankey plot
    stats_dict_to_sankey_file_in(clipping_stats_dict_num, paired_stats_dict_num, link_stats_clipping, link_stats_paired)

    # filter paired and clipping linkages
    if genomic_seq_type == 'mag':
        filter_linkages_iteratively(link_stats_paired, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, min_paired_linkages, link_stats_paired_filtered)
        filter_linkages_iteratively(link_stats_clipping, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, 0, link_stats_clipping_filtered)
    if genomic_seq_type == 'ctg':
        filter_linkages_iteratively(link_stats_paired, 'Number', pairwise_16s_iden_dict, mean_depth_dict_ctg, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, min_paired_linkages, link_stats_paired_filtered)
        filter_linkages_iteratively(link_stats_clipping, 'Number', pairwise_16s_iden_dict, mean_depth_dict_ctg, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, 0, link_stats_clipping_filtered)

    # combine_paired_and_clipping_linkages and get summary table
    combine_paired_and_clipping_linkages(link_stats_paired_filtered, link_stats_clipping_filtered, link_stats_combined_table, link_stats_combined)


    ####################################################################################################################
    ############################################### second round linking ###############################################
    ####################################################################################################################

    gnm_ctg_connector  = '___'
    dict_key_connector = '__|__'
    os.mkdir(step_2_wd)


    #################### get the sequences of 1st round unlinked marker genes and genomic sequences ####################

    report_and_log(('Step 2: get unlinked marker genes and genomes'), pwd_log_file, keep_quiet)

    # get linked marker genes and genomic sequences in step 1
    linked_marker_gene_set = set()
    linked_genomic_seq_set = set()
    for each_link in open(link_stats_combined_table):
        if not each_link.startswith('MarkerGene	GenomicSeq	Paired	Clipping'):
            each_link_split = each_link.strip().split('\t')
            linked_marker_gene_set.add(each_link_split[0])
            linked_genomic_seq_set.add(each_link_split[1])

    # get the sequence of unlinked marker genes
    marker_gene_seqs_1st_round_unlinked_handle = open(marker_gene_seqs_1st_round_unlinked, 'w')
    for marker_gene_record in SeqIO.parse(marker_gene_seqs, 'fasta'):
        if marker_gene_record.id not in linked_marker_gene_set:
            marker_gene_seqs_1st_round_unlinked_handle.write('>%s\n' % marker_gene_record.id)
            marker_gene_seqs_1st_round_unlinked_handle.write('%s\n' % marker_gene_record.seq)
    marker_gene_seqs_1st_round_unlinked_handle.close()

    # get the sequence of unlinked genomic seqs
    if genomic_seq_type == 'mag':

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

    # get the sequence of unlinked metagenomic assemblies
    if genomic_seq_type == 'ctg':
        combined_1st_round_unlinked_ctgs_handle = open(combined_1st_round_unlinked_ctgs, 'w')
        for ctg_record in SeqIO.parse(genomic_assemblies, 'fasta'):
            if ctg_record.id not in linked_genomic_seq_set:
                SeqIO.write(ctg_record, combined_1st_round_unlinked_ctgs_handle, 'fasta')
        combined_1st_round_unlinked_ctgs_handle.close()


    ############################################### get reads to extract ###############################################

    report_and_log(('Step 2: get unmapped reads with mates mapped to contig ends'), pwd_log_file, keep_quiet)
    reads_to_extract_to_ref_dict_gnm = get_free_living_mate(combined_1st_round_unlinked_mags, reads_file_r1, reads_file_r2, end_seq_len, minCigarM, max_gap_to_end, pwd_bowtie2_build_exe, pwd_bowtie2_exe, num_threads)

    report_and_log(('Step 2: get unmapped reads with mates mapped to 16S ends'), pwd_log_file, keep_quiet)
    reads_to_extract_to_ref_dict_16s = get_free_living_mate(marker_gene_seqs_1st_round_unlinked, reads_file_r1, reads_file_r2, end_seq_len, minCigarM, max_gap_to_end, pwd_bowtie2_build_exe, pwd_bowtie2_exe, num_threads)

    # write out free_living_mate_gnm
    free_living_mate_gnm_handle = open(free_living_mate_gnm, 'w')
    for fl_read_gnm in reads_to_extract_to_ref_dict_gnm:
        free_living_mate_gnm_handle.write('%s\t%s\n' % (fl_read_gnm, reads_to_extract_to_ref_dict_gnm[fl_read_gnm]))
    free_living_mate_gnm_handle.close()

    # write out free_living_mate_16s
    free_living_mate_16s_handle = open(free_living_mate_16s, 'w')
    for fl_read_16s in reads_to_extract_to_ref_dict_16s:
        free_living_mate_16s_handle.write('%s\t%s\n' % (fl_read_16s, reads_to_extract_to_ref_dict_16s[fl_read_16s]))
    free_living_mate_16s_handle.close()


    ##################################################  extract reads ##################################################

    report_and_log(('Step 2: extracting unmapped reads with mates mapped to contig/16S ends'), pwd_log_file, keep_quiet)

    os.mkdir(extracted_reads_folder)

    extract_list_gnm = set()
    extract_list_16s = set()
    extract_list_combined_r1 = set()
    extract_list_combined_r2 = set()
    for fl_mate_gnm in reads_to_extract_to_ref_dict_gnm:
        extract_list_gnm.add(fl_mate_gnm)
        if fl_mate_gnm[-2:] == '.1':
            extract_list_combined_r1.add(fl_mate_gnm)
        if fl_mate_gnm[-2:] == '.2':
            extract_list_combined_r2.add(fl_mate_gnm)
    for fl_mate_16s in reads_to_extract_to_ref_dict_16s:
        extract_list_16s.add(fl_mate_16s)
        if fl_mate_16s[-2:] == '.1':
            extract_list_combined_r1.add(fl_mate_16s)
        if fl_mate_16s[-2:] == '.2':
            extract_list_combined_r2.add(fl_mate_16s)

    # remove overlap
    extract_list_combined_r1_no_overlap = set()
    for r1 in extract_list_combined_r1:
        if (r1 in extract_list_gnm) and (r1 in extract_list_16s):
            pass
        else:
            extract_list_combined_r1_no_overlap.add(r1)

    extract_list_combined_r2_no_overlap = set()
    for r2 in extract_list_combined_r2:
        if (r2 in extract_list_gnm) and (r2 in extract_list_16s):
            pass
        else:
            extract_list_combined_r2_no_overlap.add(r2)

    # extract reads with multiprocessing
    extracted_reads_with_multiprocessing(reads_file_r1, reads_file_r2, extract_list_combined_r1_no_overlap, extract_list_combined_r2_no_overlap, extracted_reads_folder, num_threads)

    # combine extracted reads
    os.system('cat %s/*.fasta > %s' % (extracted_reads_folder, extracted_reads_cbd))
    os.system('rm -r %s' % extracted_reads_folder)


    ############################################### assemble and mapping ###############################################

    report_and_log(('Step 2: running SPAdes with extracted reads'), pwd_log_file, keep_quiet)

    #spades_cmd = '%s --meta --careful -s %s -o %s -t %s -k 21,33,55,75,99,127 --only-assembler > %s' % (pwd_spades_exe, extracted_reads_cbd, spades_wd, num_threads, spades_log)
    #spades_cmd = '%s --meta --careful -s %s -o %s -t %s -k 21,33,55,75,99,127 > %s' % (pwd_spades_exe, extracted_reads_cbd, spades_wd, num_threads, spades_log)
    spades_cmd = '%s --meta -s %s -o %s -t %s -k 21,33,55,75,99,127 --only-assembler > %s' % (pwd_spades_exe, extracted_reads_cbd, spades_wd, num_threads, spades_log)
    print(spades_cmd)
    os.system(spades_cmd)

    report_and_log(('Step 2: mapping extracted reads to SPAdes assemblies'), pwd_log_file, keep_quiet)

    spades_assemblies_no_ext = '.'.join(spades_assemblies.split('.')[:-1])
    index_ref_cmd = '%s -f %s %s --quiet --threads %s' % (pwd_bowtie2_build_exe, spades_assemblies, spades_assemblies_no_ext, num_threads)
    mapping_cmd = '%s -x %s -U %s -S %s --threads %s -f --quiet' % (pwd_bowtie2_exe, spades_assemblies_no_ext, extracted_reads_cbd, sam_file, num_threads)
    os.system(index_ref_cmd)
    os.system(mapping_cmd)


    #################################################### parse sam file ####################################################

    report_and_log(('Step 2: parsing sam file'), pwd_log_file, keep_quiet)

    gap_seq_to_reads_dict = {}
    for each_line in open(sam_file):
        if not each_line.startswith('@'):
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            ref_id = each_line_split[2]
            cigar = each_line_split[5]
            cigar_splitted = cigar_splitter(cigar)
            if (len(cigar_splitted) == 1) and (cigar[-1] == 'M'):
                if ref_id not in gap_seq_to_reads_dict:
                    gap_seq_to_reads_dict[ref_id] = [read_id]
                else:
                    gap_seq_to_reads_dict[ref_id].append(read_id)


    report_and_log(('Step 2: linking genomes/16Ss to Spades assemblies'), pwd_log_file, keep_quiet)

    stats_GapFilling_file_16s_handle = open(stats_GapFilling_file_16s, 'w')
    stats_GapFilling_file_ctg_handle = open(stats_GapFilling_file_ctg, 'w')
    stats_GapFilling_file_16s_handle.write('Gap_seq,s16,Number\n')
    stats_GapFilling_file_ctg_handle.write('Gap_seq,ctg,Number\n')
    for gap_seq in gap_seq_to_reads_dict:
        gap_seq_mapped_reads = gap_seq_to_reads_dict[gap_seq]
        gap_seq_mapped_reads_linkages_16s = {}
        gap_seq_mapped_reads_linkages_ctg = {}
        gap_seq_mapped_reads_linkages_gnm = {}
        for mapped_read in gap_seq_mapped_reads:

            # get mate ref
            mapped_read_mate_ref = ''
            mapped_read_mate_ref_gnm = ''
            mate_ref_is_16s = False
            mate_ref_is_ctg = False
            if mapped_read in reads_to_extract_to_ref_dict_gnm:
                mapped_read_mate_ref = reads_to_extract_to_ref_dict_gnm[mapped_read]
                mapped_read_mate_ref_gnm = mapped_read_mate_ref.split(gnm_ctg_connector)[0]
                mate_ref_is_ctg = True
            if mapped_read in reads_to_extract_to_ref_dict_16s:
                mapped_read_mate_ref = reads_to_extract_to_ref_dict_16s[mapped_read]
                mate_ref_is_16s = True

            # store 16s mate in dict
            if (mate_ref_is_16s is True) and (mate_ref_is_ctg is False):
                if mapped_read_mate_ref not in gap_seq_mapped_reads_linkages_16s:
                    gap_seq_mapped_reads_linkages_16s[mapped_read_mate_ref] = 1
                else:
                    gap_seq_mapped_reads_linkages_16s[mapped_read_mate_ref] += 1

            # store ctg mate in dict
            if (mate_ref_is_16s is False) and (mate_ref_is_ctg is True):

                if mapped_read_mate_ref not in gap_seq_mapped_reads_linkages_ctg:
                    gap_seq_mapped_reads_linkages_ctg[mapped_read_mate_ref] = 1
                else:
                    gap_seq_mapped_reads_linkages_ctg[mapped_read_mate_ref] += 1

                if mapped_read_mate_ref_gnm not in gap_seq_mapped_reads_linkages_gnm:
                    gap_seq_mapped_reads_linkages_gnm[mapped_read_mate_ref_gnm] = 1
                else:
                    gap_seq_mapped_reads_linkages_gnm[mapped_read_mate_ref_gnm] += 1

        for each_16s in gap_seq_mapped_reads_linkages_16s:
            stats_GapFilling_file_16s_handle.write('%s,%s,%s\n' % (gap_seq, each_16s, gap_seq_mapped_reads_linkages_16s[each_16s]))

        for each_ctg in gap_seq_mapped_reads_linkages_ctg:
            stats_GapFilling_file_ctg_handle.write('%s,%s,%s\n' % (gap_seq, each_ctg, gap_seq_mapped_reads_linkages_ctg[each_ctg]))

    stats_GapFilling_file_16s_handle.close()
    stats_GapFilling_file_ctg_handle.close()

    get_best_ctg_or_16s_for_gap_seq_iteratively(stats_GapFilling_file_16s, 'Number', min_read_num, stats_GapFilling_file_filtered_16s)
    get_best_ctg_or_16s_for_gap_seq_iteratively(stats_GapFilling_file_ctg, 'Number', min_read_num, stats_GapFilling_file_filtered_ctg)

    gap_seq_to_16s_dict = {}
    for each_line in open(stats_GapFilling_file_filtered_16s):
        if not each_line.startswith('Gap_seq,'):
            each_line_split = each_line.strip().split(',')
            gap_seq_id = each_line_split[0]
            s16_id = each_line_split[1]
            link_num = int(each_line_split[2])
            gap_seq_to_16s_dict[gap_seq_id] = {s16_id: link_num}

    gap_seq_to_gnm_dict = {}
    for each_match in open(stats_GapFilling_file_filtered_ctg):
        if not each_match.startswith('Gap_seq,'):
            match_split = each_match.strip().split(',')
            gap_seq_id = match_split[0]
            ctg_id = match_split[1]
            ctg_gnm = ctg_id.split(gnm_ctg_connector)[0]
            linkage_num = int(match_split[2])
            gap_seq_to_gnm_dict[gap_seq_id] = {ctg_gnm: linkage_num}

    gnm_to_16s_linkage_dict = {}
    for gap_seq in gap_seq_to_gnm_dict:
        if gap_seq in gap_seq_to_16s_dict:
            current_gap_seq_matched_gnm = list(gap_seq_to_gnm_dict[gap_seq].items())[0][0]
            current_gap_seq_matched_gnm_link_num = list(gap_seq_to_gnm_dict[gap_seq].items())[0][1]
            current_gap_seq_matched_16s = list(gap_seq_to_16s_dict[gap_seq].items())[0][0]
            current_gap_seq_matched_16s_link_num = list(gap_seq_to_16s_dict[gap_seq].items())[0][1]
            gnm_to_16s_key = '%s%s%s' % (current_gap_seq_matched_gnm, dict_key_connector, current_gap_seq_matched_16s)
            gnm_to_16s_link_num = current_gap_seq_matched_gnm_link_num + current_gap_seq_matched_16s_link_num
            if gnm_to_16s_key not in gnm_to_16s_linkage_dict:
                gnm_to_16s_linkage_dict[gnm_to_16s_key] = gnm_to_16s_link_num
            else:
                gnm_to_16s_linkage_dict[gnm_to_16s_key] += gnm_to_16s_link_num

    stats_GapFilling_file_handle = open(stats_GapFilling_file, 'w')
    stats_GapFilling_file_handle.write('MarkerGene,GenomicSeq,Number\n')
    for gnm_to_16s in gnm_to_16s_linkage_dict:
        id_gnm = gnm_to_16s.split(dict_key_connector)[0]
        id_16s = gnm_to_16s.split(dict_key_connector)[1]
        stats_GapFilling_file_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (id_16s, id_gnm, gnm_to_16s_linkage_dict[gnm_to_16s]))
    stats_GapFilling_file_handle.close()

    filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, min_read_num, stats_GapFilling_file_filtered)


    ####################################################################################################################
    ####################################### combine linkages from the two steps ########################################
    ####################################################################################################################

    report_and_log(('Combining linkages from step 1 and 2'), pwd_log_file, keep_quiet)

    combined_linkage_file_handle     = open(combined_linkage_file, 'w')
    combined_linkage_file_tmp_handle = open(combined_linkage_file_tmp, 'w')
    combined_linkage_file_handle.write('MarkerGene\tGenomicSeq\tLinkage\tStep\n')
    combined_linkage_file_tmp_handle.write('MarkerGene,GenomicSeq,Number\n')
    for step_1_link in open(link_stats_paired_filtered):
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


    ####################################################################################################################
    ####################################################### plot #######################################################
    ####################################################################################################################

    report_and_log(('Visualising linkages'), pwd_log_file, keep_quiet)

    # get plot height
    MarkerGenes_paired = set()
    GenomicSeqs_paired = set()
    for filtered_paired in open(combined_linkage_file_tmp):
        if not filtered_paired.startswith('MarkerGene,GenomicSeq,Number'):
            filtered_paired_split = filtered_paired.strip().split(',')
            MarkerGenes_paired.add(filtered_paired_split[0])
            GenomicSeqs_paired.add(filtered_paired_split[1])

    # calculate plot height
    plot_height_paired = 500 if max([len(MarkerGenes_paired), len(GenomicSeqs_paired)]) <= 25 else max([len(MarkerGenes_paired), len(GenomicSeqs_paired)]) * 20

    # prepare commands
    cmd_sankey_paired = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, combined_linkage_file_tmp,   600, plot_height_paired)

    # plot
    if len(open(link_stats_paired_filtered).readlines()) == 1:
        report_and_log(('No data in %s, plotting skipped' % (combined_linkage_file_tmp)), pwd_log_file, keep_quiet)
    else:
        os.system(cmd_sankey_paired)

    os.system('mv %s %s' % (combined_linkage_file_tmp_html, combined_linkage_file_html))


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

    if keep_temp is False:

        report_and_log(('Removing temporary files'), pwd_log_file, keep_quiet)
        os.remove(pwd_samfile)
        os.remove(clipping_reads_matched_part)
        os.remove(clipping_reads_not_matched_part_seq)
        os.remove(clipping_reads_not_matched_part_seq_blastn)
        os.remove(unmapped_paired_reads_file)
        os.remove(unmapped_paired_reads_blastn)
        os.remove(link_stats_clipping)
        os.remove(link_stats_paired)


    # Final report
    report_and_log(('Done!'), pwd_log_file, keep_quiet)


######################################################### main #########################################################

if __name__ == '__main__':

    link_16s_parser = argparse.ArgumentParser(description='Link MAGs with marker genes', usage=link_Marker_MAG_usage)

    link_16s_parser.add_argument('-p',               required=True,                                     help='output prefix')
    link_16s_parser.add_argument('-r1',              required=True,                                     help='paired reads r1')
    link_16s_parser.add_argument('-r2',              required=True,                                     help='paired reads r2')
    link_16s_parser.add_argument('-r16s',            required=False,                                    help='16S reads')
    link_16s_parser.add_argument('-m',               required=True,                                     help='marker gene sequences')
    link_16s_parser.add_argument('-g',               required=False,                default=None,       help='genomic sequences')
    link_16s_parser.add_argument('-mag',             required=False,                default=None,       help='metagenome-assembled-genome (MAG) folder')
    link_16s_parser.add_argument('-x',               required=False,                default='fasta',    help='MAG file extension, default: fasta')
    link_16s_parser.add_argument('-depth',           required=False, type=float,    default=0,          help='minimum depth multiple between 16S and  genomic sequences, a value of no higher than 0.2 is recommended, default: 0')
    link_16s_parser.add_argument('-s1_cigarM',       required=False, type=int,      default=30,         help='cigarM cutoff, default: 30')
    link_16s_parser.add_argument('-s1_cigarS',       required=False, type=int,      default=30,         help='cigarS cutoff, default: 30')
    link_16s_parser.add_argument('-s1_ri',           required=False, type=float,    default=100,        help='identity cutoff, default: 100')
    link_16s_parser.add_argument('-s1_rc',           required=False, type=float,    default=100,        help='coverage cutoff, default: 100')
    link_16s_parser.add_argument('-s1_mi',           required=False, type=float,    default=99.5,       help='within genome 16S identity cutoff, default: 99.5')
    link_16s_parser.add_argument('-s1_mc',           required=False, type=float,    default=90,         help='alignment coverage cutoff for calculating 16S identity, default: 90')
    link_16s_parser.add_argument('-s1_ma',           required=False, type=int,      default=500,        help='alignment length cutoff for calculating 16S identity, default: 500')
    link_16s_parser.add_argument('-s1_mpl',          required=False, type=int,      default=10,         help='minimum number of paired reads provided linkages to report, default: 10')
    link_16s_parser.add_argument('-s2_e',            required=False, type=int,      default=3000,       help='end length for mapping, default: 3000')
    link_16s_parser.add_argument('-s2_m',            required=False, type=int,      default=50,         help='minCigarM, default: 50')
    link_16s_parser.add_argument('-s2_g',            required=False, type=int,      default=200,        help='max_gap_to_end, default: 200')
    link_16s_parser.add_argument('-s2_r',            required=False, type=int,      default=3,          help='min_read_num, default: 3')
    link_16s_parser.add_argument('-t',               required=False, type=int,      default=1,          help='number of threads, default: 1')
    link_16s_parser.add_argument('-quiet',           required=False, action="store_true",               help='not report progress')
    link_16s_parser.add_argument('-force',           required=False, action="store_true",               help='force overwrite existing results')
    link_16s_parser.add_argument('-tmp',             required=False, action="store_true",               help='keep temporary files')
    link_16s_parser.add_argument('-test_mode',       required=False, action="store_true",               help='only for debugging, do not provide')
    link_16s_parser.add_argument('-bbmap',           required=False, action="store_true",               help='run bbmap, instead of bowtie')

    args = vars(link_16s_parser.parse_args())

    link_16s(args, config_dict)


To_do = '''

2. where does the paired read of the clipping mapped read mapped to? (should take into consideration!!!)
4. how to incorporate the taxonomy of MAGs and 16S sequences
5. the effect of sequencing depth, insert size and read length
6. the depth of 16S sequences always not lower than the genome they come from
7. split sam file
8. with no_ambiguous option, 16S rRNA gene sequences need to be dereplicated. (include dereplication step? with identity and coverage cutoffs?)
9. check whether input file exist!
10. minimum number of linkages to report?
11. (doesn't work)!!! to ignore list even without assignment (to handle situations like DM_m4, meanwhile capicable of not assign very diverde 16S (e.g. <98% identity) to the same genome)
12. check the structure of assembled sequences? v1, 2, 3 or v4, 5, 6? how?
13. add "16S_reads" to SortMeRNA's output prefix
14. which depth to use, genome level or contig level
15. also at ctg level
16. check if spades failed 

# on Mac
export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH
export PATH=/Users/songweizhi/Softwares/usearch:$PATH
export PATH=/Users/songweizhi/Softwares/samtools-1.11/bin:$PATH
cd /Users/songweizhi/Desktop/MarkerMAG_wd/Mac_test
~/PycharmProjects/MarkerMAG/MarkerMAG/link_16s.py -p Test_0.1 -r1 combined_5x_R1_0.1.fasta -r2 combined_5x_R2_0.1.fasta -r16s ISS_2x_even_16S_reads.fasta -m combined_16S.ffn -mag ref_genomes -x fna -t 4 -tmp -force -test_mode
~/PycharmProjects/MarkerMAG/MarkerMAG/link_16s.py -p Test -r1 combined_5x_R1.fasta -r2 combined_5x_R2.fasta -r16s ISS_5x_even.fasta -m combined_16S.ffn -mag ref_genomes -x fna -t 4 -tmp -force -test_mode

# on Katana
module unload python
module load python/3.7.3
source ~/mypython3env/bin/activate
module unload R
module load R/4.0.2
module load blast+/2.9.0
module load bowtie/2.3.5.1
module load samtools/1.10
module load spades/3.14.0
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26
python3 link_16s.py -p Test_s2g100 -mi 98 -s2g 100 -mc 30 -mpl 10 -r1 MBARC26_R1.fasta -r2 MBARC26_R2.fasta -r16s MBARC26_Matam16S_wd/MBARC26.fasta -m MBARC26_Matam16S_wd/MBARC26_all_depth_assemblies.dereplicated_99.5_500bp.fasta -mag /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/link_mag_ref_wd_spade/Refined_refined_bins_renamed -x fna -t 16 -tmp -test_mode -force -depth 0
python3 link_16s.py -p Test_s2g200 -mi 98 -s2g 200 -mc 30 -mpl 10 -r1 MBARC26_R1.fasta -r2 MBARC26_R2.fasta -r16s MBARC26_Matam16S_wd/MBARC26.fasta -m MBARC26_Matam16S_wd/MBARC26_all_depth_assemblies.dereplicated_99.5_500bp.fasta -mag /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/link_mag_ref_wd_spade/Refined_refined_bins_renamed -x fna -t 16 -tmp -test_mode -force -depth 0
python3 link_16s.py -p Two_Steps_0.03 -mi 98 -mc 30 -mpl 10 -r1 MBARC26_R1_0.05.fasta -r2 MBARC26_R2_0.05.fasta -r16s MBARC26_Matam16S_wd/MBARC26.fasta -m MBARC26_Matam16S_wd/MBARC26_all_depth_assemblies.dereplicated_99.5_500bp.fasta -mag /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/link_mag_ref_wd_spade/Refined_refined_bins_renamed -x fna -t 16 -tmp -test_mode -force -depth 0

'''
