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
               'get_sankey_plot_R'  : '%s/get_sankey_plot.R' % '/'.join(os.path.realpath(__file__).split('/')[:-1])}


link_Marker_MAG_usage = '''
=================================== MarkerMAG example commands ===================================

MarkerMAG link -p Test -r1 R1.fasta -r2 R2.fasta -m 16S_seqs.fa -g contig.fasta -t 4
MarkerMAG link -p Test -r1 R1.fasta -r2 R2.fasta -m 16S_seqs.fa -mag MAGs -x fa -t 4

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


def get_cigar_matched_and_mismatch_pct(cigar):

    cigar_splitted = cigar_splitter(cigar)

    total_len_with_s = 0
    total_len_without_s = 0
    matched_seq_len = 0
    mismatched_seq_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]

        # get total len with/without s
        if each_part_cate in {'M', 'I', '=', 'X'}:
            total_len_with_s += each_part_len
            total_len_without_s += each_part_len

        # get total len with s
        if each_part_cate == 'S':
            total_len_with_s += each_part_len

        # get matched part len
        if each_part_cate == '=':
            matched_seq_len += each_part_len

        # get mismatched part len
        if each_part_cate in {'I', 'X', 'D'}:
            mismatched_seq_len += each_part_len

    matched_pct    = float("{0:.2f}".format(matched_seq_len*100/total_len_with_s))
    mismatch_pct = float("{0:.2f}".format(mismatched_seq_len*100/total_len_without_s))

    return matched_pct, mismatch_pct


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
    reads_fmt        = argument_list[1]
    reads_to_extract = argument_list[2]
    reads_file_out   = argument_list[3]

    reads_file_out_handle = open(reads_file_out, 'w')
    for read_record in SeqIO.parse(reads_file_in, reads_fmt):
        if read_record.id in reads_to_extract:
            if reads_fmt == 'fasta':
                reads_file_out_handle.write('>%s\n' % read_record.id)
                reads_file_out_handle.write('%s\n' % read_record.seq)
            if reads_fmt == 'fastq':
                SeqIO.write(read_record, reads_file_out_handle, 'fastq')
    reads_file_out_handle.close()


def extracted_reads_with_multiprocessing(reads_r1, reads_r2, reads_fmt, r1_to_extract, r2_to_extract, output_folder, num_threads):

    solely_perfectly_mapped_reads_r1_splitted = split_list(r1_to_extract, num_threads // 2)
    solely_perfectly_mapped_reads_r2_splitted = split_list(r2_to_extract, num_threads // 2)

    argument_list_for_extract_reads_worker = []
    extract_reads_file_index_r1 = 1
    for reads_subset_r1 in solely_perfectly_mapped_reads_r1_splitted:
        current_output_file = '%s/extract_r1_subset_%s.%s' % (output_folder, extract_reads_file_index_r1, reads_fmt)
        argument_list_for_extract_reads_worker.append([reads_r1, reads_fmt, reads_subset_r1, current_output_file])
        extract_reads_file_index_r1 += 1

    extract_reads_file_index_r2 = 1
    for reads_subset_r2 in solely_perfectly_mapped_reads_r2_splitted:
        current_output_file = '%s/extract_r2_subset_%s.%s' % (output_folder, extract_reads_file_index_r2, reads_fmt)
        argument_list_for_extract_reads_worker.append([reads_r2, reads_fmt, reads_subset_r2, current_output_file])
        extract_reads_file_index_r2 += 1

    # extract reads with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(extract_reads_worker, argument_list_for_extract_reads_worker)
    pool.close()
    pool.join()


def paired_blast_results_to_dict(blastn_results, iden_cutoff, query_cov_cutoff):

    query_to_subject_list_dict = {}
    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        subject = blast_hit_split[1]
        subject_with_prefix = 'GenomicSeq__%s' % subject
        iden = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_len = int(blast_hit_split[12])
        subject_len = int(blast_hit_split[13])
        coverage_q = float(align_len) * 100 / float(query_len)

        # for perfect hits
        if (iden >= iden_cutoff) and (coverage_q == 100):
            if query not in query_to_subject_list_dict:
                query_to_subject_list_dict[query] = [subject_with_prefix]
            else:
                query_to_subject_list_dict[query].append(subject_with_prefix)

        # for nearly perfect hits
        elif (iden >= iden_cutoff) and (query_cov_cutoff <= coverage_q < 100):
            s_l = sorted([int(blast_hit_split[8]), int(blast_hit_split[9])])[0]
            s_r = sorted([int(blast_hit_split[8]), int(blast_hit_split[9])])[1]
            subject_min_gap = min([s_l, (subject_len - s_r)])
            if subject_min_gap <= 5:
                if query not in query_to_subject_list_dict:
                    query_to_subject_list_dict[query] = [subject_with_prefix]
                else:
                    query_to_subject_list_dict[query].append(subject_with_prefix)

    return query_to_subject_list_dict


def paired_blast_results_to_dict_by_mapping(unmapped_paired_reads_mapping_results):

    cigar_M_pct_min_value = 50
    mismatch_pct_max_value = 1

    unmapped_paired_reads_to_ctg_dict_by_mapping = {}

    ref_len_dict = {}
    for unmapped_read in open(unmapped_paired_reads_mapping_results):
        unmapped_read_split = unmapped_read.strip().split('\t')

        # get ref len dict
        if unmapped_read.startswith('@'):
            ref_id = ''
            ref_len = 0
            for each_element in unmapped_read_split:
                if each_element.startswith('SN:'):
                    ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_len = int(each_element[3:])
            ref_len_dict[ref_id] = ref_len

        else:
            cigar = unmapped_read_split[5]
            if cigar != '*':
                unmapped_read_split = unmapped_read.strip().split('\t')
                read_id = unmapped_read_split[0]
                ref_id = unmapped_read_split[2]
                ref_len = ref_len_dict[ref_id]
                ref_id_with_prefix = 'GenomicSeq__%s' % unmapped_read_split[2]
                ref_pos = int(unmapped_read_split[3])
                read_seq = unmapped_read_split[9]
                read_len = len(read_seq)
                cigar_splitted = cigar_splitter(cigar)
                qualified_unmapped_read = False

                # e.g. 189=
                if ('=' in cigar) and (len(cigar_splitted) == 1):
                    qualified_unmapped_read = True

                elif len(cigar_splitted) == 2:

                    # e.g. ['44S', '176=']
                    if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == '='):
                        cigar_M_pct = int(cigar_splitted[1][:-1]) * 100 / read_len
                        if (ref_pos == 1) and (cigar_M_pct >= cigar_M_pct_min_value):
                            qualified_unmapped_read = True

                    # e.g. ['154=', '66S']
                    if (cigar_splitted[0][-1] == '=') and (cigar_splitted[1][-1] == 'S'):
                        cigar_M_pct = int(cigar_splitted[0][:-1]) * 100 / read_len
                        matched_to_bp = ref_pos + int(cigar_splitted[0][:-1]) - 1
                        if (matched_to_bp == ref_len) and (cigar_M_pct >= cigar_M_pct_min_value):
                            qualified_unmapped_read = True


                elif len(cigar_splitted) == 3:

                    # e.g. ['134=', '1X', '22='], ['215=', '1X', '4=']
                    if (cigar_splitted[0][-1] == '=') and (cigar_splitted[2][-1] == '='):
                        mismatch_pct = int(cigar_splitted[1][:-1])*100/read_len
                        if mismatch_pct <= mismatch_pct_max_value:
                            qualified_unmapped_read = True

                elif len(cigar_splitted) == 4:

                    # e.g. ['118=', '1X', '99=', '1S'], ['8=', '1D', '127=', '84S']
                    if (cigar_splitted[0][-1] == '=') and (cigar_splitted[2][-1] == '=') and (cigar_splitted[3][-1] == 'S'):
                        cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1]))*100/read_len
                        if ((ref_pos + int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1])) == ref_len) and (int(cigar_splitted[1][:-1]) <= 2):
                            if cigar_M_pct >= cigar_M_pct_min_value:
                                qualified_unmapped_read = True

                    # e.g. ['51S', '158=', '1X', '10='] and ref_pos = 1
                    elif (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == '=') and (cigar_splitted[3][-1] == '='):
                        cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1]))*100/read_len
                        if (ref_pos == 1) and (int(cigar_splitted[2][:-1]) <= 1) and (cigar_M_pct >= cigar_M_pct_min_value):
                            qualified_unmapped_read = True

                # add to dict
                if qualified_unmapped_read is True:
                    if read_id not in unmapped_paired_reads_to_ctg_dict_by_mapping:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id] = [ref_id_with_prefix]
                    else:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id].append(ref_id_with_prefix)

    return unmapped_paired_reads_to_ctg_dict_by_mapping


def paired_blast_results_to_dict_by_mapping_new(unmapped_paired_reads_mapping_results, perfect_match_min_cigar_M_pct, global_max_mismatch_pct):

    unmapped_paired_reads_to_ctg_dict_by_mapping = {}
    for unmapped_read in open(unmapped_paired_reads_mapping_results):

        if not unmapped_read.startswith('@'):
            unmapped_read_split = unmapped_read.strip().split('\t')
            cigar = unmapped_read_split[5]
            if cigar != '*':
                read_id = unmapped_read_split[0]
                ref_id_with_prefix = 'GenomicSeq__%s' % unmapped_read_split[2]
                cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)

                if (cigar_match_pct >= perfect_match_min_cigar_M_pct) and (cigar_mismatch_pct <= global_max_mismatch_pct):

                    if read_id not in unmapped_paired_reads_to_ctg_dict_by_mapping:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id] = [ref_id_with_prefix]
                    else:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id].append(ref_id_with_prefix)

    return unmapped_paired_reads_to_ctg_dict_by_mapping


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


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, file_out):

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


def get_free_living_mate(ref_in, reads_r1, reads_r2, end_seq_len, max_gap_to_end, num_threads, pwd_bbmap_exe, bbmap_memory, perfect_match_min_cigar_M_pct, global_max_mismatch_pct, pwd_log_file, keep_quiet):

    ref_in_path, ref_in_basename, ref_in_ext = sep_path_basename_ext(ref_in)

    ref_subset      = '%s/%s_ends_%sbp%s'                % (ref_in_path, ref_in_basename, end_seq_len, ref_in_ext)
    sam_file        = '%s/%s_ends_%sbp.sam'              % (ref_in_path, ref_in_basename, end_seq_len)
    bbmap_stderr    = '%s/%s_ends_%sbp_bbmap_stderr.txt' % (ref_in_path, ref_in_basename, end_seq_len)

    # get ref seqs subset
    report_and_log(('Round 2: map reads to %s, get ends sequences' % ref_in_basename), pwd_log_file, keep_quiet)
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
    report_and_log(('Round 2: map reads to %s, mapping' % ref_in_basename), pwd_log_file, keep_quiet)
    bbmap_parameter_round2 = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % (num_threads, bbmap_memory)
    bbmap_cmd_round2 = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, ref_subset, reads_r1, reads_r2, sam_file, bbmap_parameter_round2, bbmap_stderr)
    os.system(bbmap_cmd_round2)

    # parse sam file
    report_and_log(('Round 2: map reads to %s, filter sam' % ref_in_basename), pwd_log_file, keep_quiet)
    all_mapped_reads = set()
    qualified_reads_dict = {}
    qualified_reads_to_ref_dict = {}
    ref_in_sam_len_dict = {}
    for each_line in open(sam_file):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('@'):
            ref_id = ''
            ref_len = 0
            for each_element in each_line_split:
                if each_element.startswith('SN:'):
                    ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_len = int(each_element[3:])
            ref_in_sam_len_dict[ref_id] = ref_len
        else:
            cigar = each_line_split[5]
            if cigar != '*':
                read_id = each_line_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_line_split[2]
                ref_len = ref_in_sam_len_dict[ref_id]
                ref_pos = int(each_line_split[3])
                cigar_splitted = cigar_splitter(cigar)
                cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)
                all_mapped_reads.add(read_id)

                # get aligned length
                aligned_len = 0
                for each_section in cigar_splitted:
                    each_section_len = int(each_section[:-1])
                    each_section_cate = each_section[-1]
                    if each_section_cate in {'D', '=', 'X'}:
                        aligned_len += each_section_len

                # filter mapped reads
                qualified_mapping = False
                if (cigar_match_pct >= perfect_match_min_cigar_M_pct) and (cigar_mismatch_pct <= global_max_mismatch_pct):

                    # check left end for contig's left end sequence
                    if ref_id[-2:] == '_l':
                        if ref_pos <= max_gap_to_end:
                            qualified_mapping = True

                    # check right end for contig's right end sequence
                    elif ref_id[-2:] == '_r':
                        if cigar_splitted[0][-1] != 'S':
                            ref_pos_end = ref_pos + aligned_len
                            if (ref_len - ref_pos_end) <= max_gap_to_end:
                                qualified_mapping = True

                    # check both ends for contigs without subset
                    else:
                        if (ref_pos <= max_gap_to_end) or ((ref_len - ref_pos - aligned_len) <= max_gap_to_end):
                            qualified_mapping = True

                # proceed with qualified mappings
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


def SeqIO_convert_worker(argument_list):

    file_in         = argument_list[0]
    file_in_fmt     = argument_list[1]
    file_out        = argument_list[2]
    file_out_fmt    = argument_list[3]
    SeqIO.convert(file_in, file_in_fmt, file_out, file_out_fmt)


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
    marker_gene_seqs                    = args['marker']
    min_16s_gnm_multiple                = args['depth']
    within_genome_minimum_iden16s       = args['s1_mi']
    cov16s                              = args['s1_mc']
    aln16s                              = args['s1_ma']
    min_paired_linkages                 = args['s1_mpl']
    min_paired_linkages_for_uniq_linked_16s  = args['s1_mplu']
    max_gap_to_end                      = args['s2_g']
    min_read_num                        = args['s2_r']
    num_threads                         = args['t']
    keep_quiet                          = args['quiet']
    force_overwrite                     = args['force']
    keep_temp                           = args['tmp']
    test_mode                           = args['test_mode']
    round_2_spades                      = args['spades']
    mira_tmp_dir                        = args['mira_tmp']
    bbmap_memory                        = args['bbmap_mem']
    global_max_mismatch_pct             = args['mismatch']
    perfect_match_min_cigar_M_pct       = args['min_M_pct']

    pwd_plot_sankey_R                   = config_dict['get_sankey_plot_R']
    pwd_makeblastdb_exe                 = 'makeblastdb'
    pwd_blastn_exe                      = 'blastn'
    pwd_bowtie2_build_exe               = 'bowtie2-build'
    pwd_bowtie2_exe                     = 'bowtie2'
    pwd_samtools_exe                    = 'samtools'
    pwd_spades_exe                      = 'spades.py'
    pwd_bbmap_exe                       = 'bbmap.sh'

    reads_iden_cutoff   = 100             # %
    reads_cov_cutoff    = 90              # %
    end_seq_len         = 1000            # bp

    # for clipping reads
    gnm_ctg_connector  = '___'
    dict_key_connector = '__|__'


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
    if 'q' in r1_ext:

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


    ######################## check genomic sequence type and prepare files for making blast db #########################

    blast_db           = ''
    blast_db_no_ext    = ''
    genomic_seq_type   = ''  # ctg or mag
    renamed_mag_folder = ''

    # check the type of input genomic sequences
    if (genomic_assemblies is not None) and (mag_folder is None):
        genomic_seq_type = 'ctg'
        metagenomic_assemblies_file_path, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension = sep_path_basename_ext(genomic_assemblies)
        blast_db_dir = '%s/%s_%s_db' % (step_1_wd, output_prefix, metagenomic_assemblies_file_basename)
        blast_db     = '%s/%s%s'     % (blast_db_dir, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension)
        blast_db_no_ext = '%s/%s'     % (blast_db_dir, metagenomic_assemblies_file_basename)

        os.mkdir(blast_db_dir)
        os.system('cp %s %s/' % (genomic_assemblies, blast_db_dir))

    elif (genomic_assemblies is None) and (mag_folder is not None):
        genomic_seq_type    = 'mag'
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
            rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename, gnm_ctg_connector)

        # combine renamed MAGs
        blast_db = '%s/%s_combined.fa' % (blast_db_dir, mag_folder_name)
        blast_db_no_ext = '%s/%s_combined' % (blast_db_dir, mag_folder_name)
        os.system('cat %s/*%s > %s' % (renamed_mag_folder, mag_file_extension, blast_db))

    else:
        print('Please provide genomic sequences either as raw assemblies (-g) or as MAGs (-mag)')
        exit()


    ########################################### define folder and file name ############################################

    marker_gene_seqs_file_path, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension = sep_path_basename_ext(marker_gene_seqs)

    pwd_log_file                                = '%s/%s.log'                                    % (working_directory, output_prefix)
    bowtie_index_dir                            = '%s/%s_index'                                  % (step_1_wd, marker_gene_seqs_file_basename)
    pwd_samfile                                 = '%s/%s.sam'                                    % (step_1_wd, marker_gene_seqs_file_basename)
    pwd_samfile_stderr                          = '%s/%s_bbmap_stderr.txt'                       % (step_1_wd, marker_gene_seqs_file_basename)
    pwd_samfile_to_mag                          = '%s/%s_combined.sam'                           % (step_1_wd, mag_folder_name)
    pwd_samfile_to_mag_stderr                   = '%s/%s_combined_bbmap_stderr.txt'              % (step_1_wd, mag_folder_name)
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
    bbmap_parameter  = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % (num_threads, bbmap_memory)

    ################################################# step 2 #################################################

    marker_gene_seqs_1st_round_unlinked         = '%s/step_1_unlinked_marker_genes.fasta'        % step_2_wd
    combined_1st_round_unlinked_mags            = '%s/step_1_unlinked_combined_gnms.fasta'       % step_2_wd
    combined_1st_round_unlinked_ctgs            = '%s/step_1_unlinked_combined_ctgs.fasta'       % step_2_wd
    free_living_mate_gnm                        = '%s/free_living_mate_ctg.txt'                  % step_2_wd
    free_living_mate_16s                        = '%s/free_living_mate_16s.txt'                  % step_2_wd
    extracted_reads_folder                      = '%s/free_living_reads'                         % step_2_wd
    extracted_reads_cbd                         = '%s/free_living_read_combined.fastq'           % step_2_wd
    extracted_reads_cbd_fasta                   = '%s/free_living_read_combined.fasta'           % step_2_wd
    spades_wd                                   = '%s/combined_free_living_reads_SPAdes_wd'      % step_2_wd
    mira_manifest                               = '%s/mira_manifest.txt'                         % step_2_wd
    mira_stdout                                 = '%s/mira_stdout.txt'                           % step_2_wd
    sam_file_mini_assembly                      = '%s/scaffolds.sam'                             % step_2_wd
    sam_file_mini_assembly_stderr               = '%s/scaffolds_bbmap_stderr.txt'                % step_2_wd
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
            report_and_log(('Round 1: calculating depth for %s' % genomic_assemblies), pwd_log_file, keep_quiet)
        if genomic_seq_type == 'mag':
            report_and_log(('Round 1: calculating depth for genomes in %s' % mag_folder), pwd_log_file, keep_quiet)

        # get mean depth for contig
        mean_depth_dict_ctg, ctg_len_dict = get_ctg_mean_depth_by_samtools_coverage(True, blast_db, reads_file_r1_fasta, reads_file_r2_fasta, '', num_threads)

        # write out ctg depth
        depth_file_ctg_handle = open(depth_file_ctg, 'w')
        for ctg in mean_depth_dict_ctg:
            depth_file_ctg_handle.write('%s\t%s\n' % (ctg, mean_depth_dict_ctg[ctg]))
        depth_file_ctg_handle.close()

        # get mean_depth_dict_gnm
        if genomic_seq_type == 'mag':

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

    marker_gene_seqs_in_wd  = '%s/%s%s' % (bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension)

    report_and_log(('Round 1: Mapping reads to marker gene sequences with bbmap'), pwd_log_file, keep_quiet)
    bbmap_index_and_mapping_cmd = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, marker_gene_seqs_in_wd, reads_file_r1, reads_file_r2, pwd_samfile, bbmap_parameter, pwd_samfile_stderr)
    report_and_log(('Round 1: Command for running bbmap exported to log file'), pwd_log_file, keep_quiet)
    report_and_log((bbmap_index_and_mapping_cmd), pwd_log_file, True)
    os.system(bbmap_index_and_mapping_cmd)


    ##################################################### extract reads ####################################################

    report_and_log(('Round 1: Extracting unmapped part of clipping mapped reads from sam file'), pwd_log_file, keep_quiet)

    # export clipping mapped reads and perfectly mapped reads
    all_mapped_reads_set = set()
    clipping_mapped_reads_list = set()
    clipping_reads_mapped_part_dict = {}
    perfectly_mapped_reads_dict = {}
    clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
    for each_read in open(pwd_samfile):
        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            cigar = each_read_split[5]
            if cigar != '*':
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_read_split[2]
                ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
                ref_pos = int(each_read_split[3])
                read_seq = each_read_split[9]
                cigar_splitted = cigar_splitter(cigar)
                cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)
                read_id_with_ref_pos = '%s__x__%s__x__%s' % (read_id, ref_id, ref_pos)
                all_mapped_reads_set.add(read_id)

                # treat_as_full_match and store into dict
                if (cigar_match_pct >= perfect_match_min_cigar_M_pct) and (
                        cigar_mismatch_pct <= global_max_mismatch_pct):

                    if read_id_base not in perfectly_mapped_reads_dict:
                        perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
                    else:
                        if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                            perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                        else:
                            perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

                # treat as clipping matched
                if (cigar_match_pct < perfect_match_min_cigar_M_pct) and (
                        cigar_mismatch_pct <= global_max_mismatch_pct):

                    # only one end of cigar is S
                    if ((cigar_splitted[0][-1] == 'S') and (cigar_splitted[-1][-1] != 'S')) or (
                            (cigar_splitted[0][-1] != 'S') and (cigar_splitted[-1][-1] == 'S')):

                        # if clipped at left
                        if cigar_splitted[0][-1] == 'S':
                            read_seq_clipped = read_seq[:int(cigar_splitted[0][:-1])]

                            # write out the clipped part
                            clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id_with_ref_pos)
                            clipping_reads_not_matched_part_seq_handle.write(read_seq_clipped + '\n')

                            # store the matching info of aligned part
                            if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                                clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                            else:
                                clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                        # if clipped at right
                        if cigar_splitted[-1][-1] == 'S':
                            read_seq_clipped = read_seq[-int(cigar_splitted[-1][:-1]):]

                            # write out the clipped part
                            clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id_with_ref_pos)
                            clipping_reads_not_matched_part_seq_handle.write(read_seq_clipped + '\n')

                            # store the matching info of aligned part
                            if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                                clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                            else:
                                clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

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
            if r2_to_extract not in all_mapped_reads_set:
                solely_perfectly_mapped_reads_r2.add(r2_to_extract)
        if strand == '2':
            r1_to_extract = '%s.1' % perfectly_mapped_read
            if r1_to_extract not in all_mapped_reads_set:
                solely_perfectly_mapped_reads_r1.add(r1_to_extract)

    # extract reads
    report_and_log(('Round 1: Extracting unmapped paired reads with %s cores' % num_threads), pwd_log_file, keep_quiet)
    os.mkdir(unmapped_paired_reads_folder)

    # extract reads with multiprocessing
    extracted_reads_with_multiprocessing(reads_file_r1_fasta, reads_file_r2_fasta, 'fasta', solely_perfectly_mapped_reads_r1, solely_perfectly_mapped_reads_r2, unmapped_paired_reads_folder, num_threads)

    # combine extracted reads
    os.system('cat %s/*.fasta > %s' % (unmapped_paired_reads_folder, unmapped_paired_reads_file))


    ############################# run blast between extracted reads and metagenomic assemblies #############################

    # run blastn
    report_and_log(('Round 1: Running blastn for unmapped parts of clipping mapped reads'), pwd_log_file, keep_quiet)
    makeblastdb_cmd     = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blast_db)
    blastn_cmd_clipping = '%s -query %s -db %s -out %s %s'                          % (pwd_blastn_exe, clipping_reads_not_matched_part_seq, blast_db, clipping_reads_not_matched_part_seq_blastn, blast_parameters)
    os.system(makeblastdb_cmd)
    os.system(blastn_cmd_clipping)

    # mapping with bbmap
    report_and_log(('Round 1: Mapping extracted reads to genomic sequences with bbmap'), pwd_log_file, keep_quiet)
    bbmap_index_and_mapping_cmd_to_mag = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, blast_db, unmapped_paired_reads_file, pwd_samfile_to_mag, bbmap_parameter, pwd_samfile_to_mag_stderr)
    report_and_log(('Round 1: Command for running bbmap exported to log file'), pwd_log_file, keep_quiet)
    report_and_log((bbmap_index_and_mapping_cmd_to_mag), pwd_log_file, True)
    os.system(bbmap_index_and_mapping_cmd_to_mag)


    ######################################### parse blast results for paired reads #########################################

    report_and_log(('Round 1: Parsing mapping results for paired reads'), pwd_log_file, keep_quiet)

    # filter mapping results
    unmapped_paired_reads_to_ctg_dict = paired_blast_results_to_dict_by_mapping_new(pwd_samfile_to_mag, perfect_match_min_cigar_M_pct, global_max_mismatch_pct)

    paired_stats_dict_num = {}
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

    report_and_log(('Round 1: Parsing blast results for clipping mapped reads'), pwd_log_file, keep_quiet)

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

    report_and_log(('Round 1: Parsing linkages'), pwd_log_file, keep_quiet)

    # prepare input file for sankey plot
    stats_dict_to_sankey_file_in(clipping_stats_dict_num, paired_stats_dict_num, link_stats_clipping, link_stats_paired)

    # filter paired and clipping linkages
    if genomic_seq_type == 'mag':
        filter_linkages_iteratively(link_stats_paired, 'Number', pairwise_16s_iden_dict,
                                    mean_depth_dict_gnm, mean_depth_dict_16s,
                                    min_16s_gnm_multiple,
                                    within_genome_minimum_iden16s, min_paired_linkages, min_paired_linkages_for_uniq_linked_16s, link_stats_paired_filtered)


        filter_linkages_iteratively(link_stats_clipping, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, 0, 0, link_stats_clipping_filtered)
    if genomic_seq_type == 'ctg':
        filter_linkages_iteratively(link_stats_paired, 'Number', pairwise_16s_iden_dict, mean_depth_dict_ctg, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, min_paired_linkages, min_paired_linkages_for_uniq_linked_16s, link_stats_paired_filtered)
        filter_linkages_iteratively(link_stats_clipping, 'Number', pairwise_16s_iden_dict, mean_depth_dict_ctg, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, 0, 0, link_stats_clipping_filtered)

    # combine_paired_and_clipping_linkages and get summary table
    combine_paired_and_clipping_linkages(link_stats_paired_filtered, link_stats_clipping_filtered, link_stats_combined_table, link_stats_combined)


    ####################################################################################################################
    ############################################### second round linking ###############################################
    ####################################################################################################################

    #################### get the sequences of 1st round unlinked marker genes and genomic sequences ####################

    os.mkdir(step_2_wd)

    report_and_log(('Round 2: get unlinked marker genes and genomes'), pwd_log_file, keep_quiet)

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

    report_and_log(('Round 2: get unmapped reads with mates mapped to contig ends'), pwd_log_file, keep_quiet)
    reads_to_extract_to_ref_dict_gnm = get_free_living_mate(combined_1st_round_unlinked_mags, reads_file_r1_fasta, reads_file_r2_fasta, end_seq_len, max_gap_to_end, num_threads, pwd_bbmap_exe, bbmap_memory, perfect_match_min_cigar_M_pct, global_max_mismatch_pct, pwd_log_file, keep_quiet)

    report_and_log(('Round 2: get unmapped reads with mates mapped to 16S ends'), pwd_log_file, keep_quiet)
    reads_to_extract_to_ref_dict_16s = get_free_living_mate(marker_gene_seqs_1st_round_unlinked, reads_file_r1_fasta, reads_file_r2_fasta, end_seq_len, max_gap_to_end, num_threads, pwd_bbmap_exe, bbmap_memory, perfect_match_min_cigar_M_pct, global_max_mismatch_pct, pwd_log_file, keep_quiet)

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

    report_and_log(('Round 2: extracting unmapped reads with mates mapped to contig/16S ends'), pwd_log_file, keep_quiet)

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
    if round_2_spades is False:
        extracted_reads_with_multiprocessing(reads_file_r1, reads_file_r2, 'fastq', extract_list_combined_r1_no_overlap, extract_list_combined_r2_no_overlap, extracted_reads_folder, num_threads)
        # combine extracted reads
        os.system('cat %s/*.fastq > %s' % (extracted_reads_folder, extracted_reads_cbd))
        os.system('rm -r %s' % extracted_reads_folder)
        SeqIO.convert(extracted_reads_cbd, 'fastq', extracted_reads_cbd_fasta, 'fasta-2line')
    else:
        extracted_reads_with_multiprocessing(reads_file_r1_fasta, reads_file_r2_fasta, 'fasta', extract_list_combined_r1_no_overlap, extract_list_combined_r2_no_overlap, extracted_reads_folder, num_threads)
        os.system('cat %s/*.fasta > %s' % (extracted_reads_folder, extracted_reads_cbd_fasta))


    ############################################### assemble and mapping ###############################################

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


    if round_2_spades is False:
        report_and_log(('Round 2: running Mira on extracted reads'), pwd_log_file, keep_quiet)
        run_mira5(output_prefix, mira_tmp_dir, step_2_wd, mira_manifest, extracted_reads_cbd, mira_stdout, force_overwrite)
        mini_assemblies = '%s/%s_mira_est_no_chimera_assembly/%s_mira_est_no_chimera_d_results/%s_mira_est_no_chimera_out.unpadded.fasta' % (step_2_wd, output_prefix, output_prefix, output_prefix)
    else:
        report_and_log(('Round 2: running SPAdes on extracted reads'), pwd_log_file, keep_quiet)
        spades_cmd = '%s -s %s -o %s -t %s -k 55,75,99,127 --only-assembler > %s' % (pwd_spades_exe, extracted_reads_cbd_fasta, spades_wd, num_threads, spades_log)
        os.system(spades_cmd)
        mini_assemblies = '%s/scaffolds.fasta' % spades_wd


    # mapping extracted reads to mini assemblies with bbmap
    report_and_log(('Round 2: mapping extracted reads to mini assemblies with bbmap'), pwd_log_file, keep_quiet)
    bbmap_cmd_miniassembly = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, mini_assemblies, extracted_reads_cbd_fasta, sam_file_mini_assembly, bbmap_parameter, sam_file_mini_assembly_stderr)
    os.system(bbmap_cmd_miniassembly)


    #################################################### parse sam file ####################################################

    report_and_log(('Round 2: parsing sam file'), pwd_log_file, keep_quiet)

    gap_seq_to_reads_dict = {}
    for each_line in open(sam_file_mini_assembly):
        if not each_line.startswith('@'):
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            ref_id = each_line_split[2]
            cigar = each_line_split[5]
            cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)

            qualified_mapping = False
            if (cigar_match_pct >= perfect_match_min_cigar_M_pct) and (cigar_mismatch_pct <= global_max_mismatch_pct):
                qualified_mapping = True

            if qualified_mapping is True:
                if ref_id not in gap_seq_to_reads_dict:
                    gap_seq_to_reads_dict[ref_id] = [read_id]
                else:
                    gap_seq_to_reads_dict[ref_id].append(read_id)


    report_and_log(('Round 2: linking genomes/16Ss to Spades assemblies'), pwd_log_file, keep_quiet)

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

            # check num difference
            link_num_max = max(int(current_gap_seq_matched_gnm_link_num), int(current_gap_seq_matched_16s_link_num))
            link_num_min = min(int(current_gap_seq_matched_gnm_link_num), int(current_gap_seq_matched_16s_link_num))
            if (link_num_min*100/link_num_max) >= 33.3:
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

    filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, min_read_num, min_read_num, stats_GapFilling_file_filtered)


    ####################################################################################################################
    ####################################### combine linkages from step 1 and 2  ########################################
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
    cmd_sankey_paired = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, combined_linkage_file_tmp, 600, plot_height_paired)

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
    link_16s_parser.add_argument('-marker',          required=True,                                     help='marker gene sequences')
    link_16s_parser.add_argument('-g',               required=False,                default=None,       help='genomic sequences')
    link_16s_parser.add_argument('-mag',             required=False,                default=None,       help='metagenome-assembled-genome (MAG) folder')
    link_16s_parser.add_argument('-x',               required=False,                default='fasta',    help='MAG file extension, default: fasta')
    link_16s_parser.add_argument('-depth',           required=False, type=float,    default=0,          help='minimum depth multiple between 16S and  genomic sequences, a value of no higher than 0.2 is recommended, default: 0')
    link_16s_parser.add_argument('-s1_mi',           required=False, type=float,    default=98,         help='within genome 16S identity cutoff, default: 98')
    link_16s_parser.add_argument('-s1_mc',           required=False, type=float,    default=30,         help='alignment coverage cutoff for calculating 16S identity, default: 30')
    link_16s_parser.add_argument('-s1_ma',           required=False, type=int,      default=500,        help='alignment length cutoff for calculating 16S identity, default: 500')
    link_16s_parser.add_argument('-s1_mpl',          required=False, type=int,      default=10,          help='minimum number of paired reads provided linkages to report, default: 10')
    link_16s_parser.add_argument('-s1_mplu',         required=False, type=int,      default=5,          help='minimum number of paired reads provided linkages to report (for uniq linked 16S, default: 5')
    link_16s_parser.add_argument('-s2_g',            required=False, type=int,      default=300,        help='max_gap_to_end, default: 300')
    link_16s_parser.add_argument('-s2_r',            required=False, type=int,      default=3,          help='min_read_num, default: 3')
    link_16s_parser.add_argument('-t',               required=False, type=int,      default=1,          help='number of threads, default: 1')
    link_16s_parser.add_argument('-quiet',           required=False, action="store_true",               help='not report progress')
    link_16s_parser.add_argument('-force',           required=False, action="store_true",               help='force overwrite existing results')
    link_16s_parser.add_argument('-tmp',             required=False, action="store_true",               help='keep temporary files')
    link_16s_parser.add_argument('-test_mode',       required=False, action="store_true",               help='only for debugging, do not provide')
    link_16s_parser.add_argument('-bbmap_mem',       required=False, type=int,      default=10,         help='bbmap memory allocation (in gigabyte), default: 10')
    link_16s_parser.add_argument('-mira_tmp',        required=False, default=None,                      help='tmp dir for mira')
    link_16s_parser.add_argument('-spades',          required=False, action="store_true",               help='run spades, instead of Mira')
    link_16s_parser.add_argument('-min_M_pct',       required=False, type=float,    default=70,         help='perfect_match_min_cigar_M_pct, default: 70')
    link_16s_parser.add_argument('-mismatch',        required=False, type=float,    default=3,          help='maximum mismatch percentage, default: 3')
    args = vars(link_16s_parser.parse_args())

    link_16s(args, config_dict)


To_do = '''

1. where do the mates of clipping mapped read mapped to? (should take into consideration!!!)
2. how to incorporate the taxonomy of MAGs and 16S sequences
4. the depth of 16S sequences always not lower than the genome they come from
5. split sam file
6. with no_ambiguous option, 16S rRNA gene sequences need to be dereplicated. (include dereplication step? with identity and coverage cutoffs?)
9. (doesn't work)!!! to ignore list even without assignment (to handle situations like DM_m4, meanwhile capicable of not assign very diverde 16S (e.g. <98% identity) to the same genome)
10. check the structure of assembled sequences? v1, 2, 3 or v4, 5, 6? how?
11. add "16S_reads" to SortMeRNA's output prefix
12. which depth to use, genome level or contig level
14. check if spades failed 
16. check Mira tmp_dir at the beginning
17. BH_ER_050417_refined_bins_combined.sam: >Kelp_659345.1__x__Refined_9___NODE_5430_length_12683_cov_6.017999 NODE_5430_length_12683_cov_6.017999__x__12552_r (space in ref id !!!)
18. !!!!!! if isfile(mira_assembly) is False: report and only export first round linkages !!!!!!
19. use wc -l to check if fastq and fasta match
20. faster way to rename and extract reads: seqtk (test it)?



# on Katana
module load python/3.7.3
source ~/mypython3env/bin/activate
module load R/4.0.2
module load blast+/2.11.0
module load bowtie/2.3.5.1
module load samtools/1.10
module load spades/3.14.0
module load gcc/8.4.0
module load boost/1.73.0-gcc8   
module load mira/v5rc2
module load java/8u201-jdk
module load bbmap/38.51

cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp
python3 link_16s.py -p Kelp_0.999_bbmap -r1 Kelp_R1.fastq -r2 Kelp_R2.fastq -marker BH_ER_050417_Matam16S_wd/BH_ER_050417_assembled_16S_uclust_0.999.fasta -mag BH_ER_050417_refined_bins -x fasta -t 6 -tmp -force -mira_tmp $TMPDIR


~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:136: MarkerGene__PS_m1,GenomicSeq__DA,18
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:151: MarkerGene__PS_m1,GenomicSeq__FP,16
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:210: MarkerGene__PS_m1,GenomicSeq__SS,12
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:293: MarkerGene__PS_m1,GenomicSeq__DM,7
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:368: MarkerGene__PS_m1,GenomicSeq__TC,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:402: MarkerGene__PS_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:431: MarkerGene__PS_m1,GenomicSeq__FA,4
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:439: MarkerGene__PS_m1,GenomicSeq__DG,4
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:494: MarkerGene__PS_m1,GenomicSeq__MS,4
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:511: MarkerGene__PS_m1,GenomicSeq__OU,3
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:595: MarkerGene__PS_m1,GenomicSeq__PS,2
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:707: MarkerGene__PS_m1,GenomicSeq__NG,2
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:783: MarkerGene__PS_m1,GenomicSeq__HR,1
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:863: MarkerGene__PS_m1,GenomicSeq__NO,1
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:961: MarkerGene__PS_m1,GenomicSeq__TR,1

~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:142: MarkerGene__DA_m3,GenomicSeq__CA,17
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:155: MarkerGene__DA_m1,GenomicSeq__CA,16
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:174: MarkerGene__DA_m2,GenomicSeq__CA,14
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:194: MarkerGene__DA_m7,GenomicSeq__CA,13
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:196: MarkerGene__DA_m6,GenomicSeq__CA,12
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:200: MarkerGene__DM_m1,GenomicSeq__CA,12
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:221: MarkerGene__DM_m2,GenomicSeq__CA,11
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:237: MarkerGene__DM_m3,GenomicSeq__CA,10
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:247: MarkerGene__DG_m8,GenomicSeq__CA,9
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:262: MarkerGene__DM_m9,GenomicSeq__CA,9
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:269: MarkerGene__DM_m5,GenomicSeq__CA,8
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:273: MarkerGene__DG_m7,GenomicSeq__CA,8
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:288: MarkerGene__TC_m1,GenomicSeq__CA,7
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:324: MarkerGene__TC_m2,GenomicSeq__CA,6
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:337: MarkerGene__DG_m1,GenomicSeq__CA,6
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:356: MarkerGene__DG_m10,GenomicSeq__CA,6
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:371: MarkerGene__SS_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:383: MarkerGene__NO_m3,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:384: MarkerGene__FP_m2,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:385: MarkerGene__NG_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:387: MarkerGene__NO_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:402: MarkerGene__PS_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:413: MarkerGene__FA_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:414: MarkerGene__FA_m2,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:419: MarkerGene__MS_m1,GenomicSeq__CA,5
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:460: MarkerGene__DG_m5,GenomicSeq__CA,4
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:519: MarkerGene__DA_m5,GenomicSeq__CA,3
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:536: MarkerGene__HB_m1,GenomicSeq__CA,3
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:561: MarkerGene__NO_m2,GenomicSeq__CA,3
~/Library/Caches/Transmit/E4749E2E-40C9-4EE0-9B33-475CC0D160FD/kdm.restech.unsw.edu.au/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_mis_5_MarkerMAG_wd/MBARC26_mis_5_step_1_wd/stats_paired_sorted.txt:582: MarkerGene__TC_m4,GenomicSeq__CA,3

MarkerGene__CA_m1,GenomicSeq__CA,2



'''
