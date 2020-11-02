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
               'blastn'             : 'blastn',
               'makeblastdb'        : 'makeblastdb',
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


def remove_unmapped_reads(sam_file, sam_file_only_mapped):

    sam_file_only_mapped_handle = open(sam_file_only_mapped, 'w')

    for each_read in open(sam_file):
        if not each_read.startswith('@'):
            ref_id = each_read.strip().split('\t')[2]
            cigar = each_read.strip().split('\t')[5]

            if (ref_id != '*') and (cigar != '*'):
                sam_file_only_mapped_handle.write(each_read)

    sam_file_only_mapped_handle.close()


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


def blast_results_to_dict_by_bitscore(blastn_results, iden_cutoff):
    # get query_to_best_bitscore_dict
    query_to_best_bitscore_dict = {}
    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        bitscore = float(blast_hit_split[11])

        if query not in query_to_best_bitscore_dict:
            query_to_best_bitscore_dict[query] = bitscore
        else:
            if bitscore > query_to_best_bitscore_dict[query]:
                query_to_best_bitscore_dict[query] = bitscore

    # filter balst results
    query_to_subject_list_dict = {}
    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        subject = blast_hit_split[1]
        subject_with_prefix = 'GenomicSeq__%s' % subject
        iden = float(blast_hit_split[2])
        bitscore = float(blast_hit_split[11])
        current_query_best_bitscore = query_to_best_bitscore_dict[query]
        if (iden >= iden_cutoff) and (bitscore == current_query_best_bitscore):
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


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, within_genome_16s_divergence_cutoff, file_out):

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


def rename_seq(ctg_file_in, ctg_file_out, prefix):

    ctg_file_out_handle = open(ctg_file_out, 'w')
    for Seq_record in SeqIO.parse(ctg_file_in, 'fasta'):
        Seq_record.id = '%s___%s' % (prefix, Seq_record.id)
        SeqIO.write(Seq_record, ctg_file_out_handle, 'fasta')
    ctg_file_out_handle.close()


def link_Marker_MAG(args, config_dict):

    ###################################################### file in/out #####################################################

    # file in
    output_prefix                       = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    genomic_assemblies                  = args['g']
    mag_folder                          = args['mag']
    mag_file_extension                  = args['x']
    marker_gene_seqs                    = args['m']
    min_cigar_M                         = args['cigarM']
    min_cigar_S                         = args['cigarS']
    reads_iden_cutoff                   = args['ri']
    reads_cov_cutoff                    = args['rc']
    iden16s                             = args['mi']
    cov16s                              = args['mc']
    aln16s                              = args['ma']
    force_overwrite                     = args['force']
    keep_quiet                          = args['quiet']
    num_threads                         = args['t']
    keep_temp                           = args['tmp']
    test_mode                           = args['test_mode']
    run_bbmap                           = args['bbmap']

    pwd_plot_sankey_R                   = config_dict['get_sankey_plot_R']
    pwd_bowtie2_exe                     = config_dict['bowtie2']
    pwd_bowtie2_build_exe               = config_dict['bowtie2_build']
    pwd_makeblastdb_exe                 = config_dict['makeblastdb']
    pwd_blastn_exe                      = config_dict['blastn']


    ################################################ check dependencies ################################################

    # check whether executables exist
    program_list = [pwd_makeblastdb_exe, pwd_blastn_exe, pwd_bowtie2_build_exe, pwd_bowtie2_exe]
    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


    ############################################# create working directory #############################################

    # create working directory
    working_directory = '%s_MarkerMAG_wd_M%s_S%s_ri%s_rc%s_mi%s_mc%s_ma%sbp' % (output_prefix, min_cigar_M, min_cigar_S, reads_iden_cutoff, reads_cov_cutoff, iden16s, cov16s, aln16s)
    if (os.path.isdir(working_directory) is True) and (force_overwrite is False):
        print('Working directory detected, program exited!')
        exit()
    else:
        force_create_folder(working_directory)


    ######################## check genomic sequence type and prepare files for making blast db #########################

    blast_db         = ''
    genomic_seq_type = ''  # ctg or mag

    # check the type of input genomic sequences
    if (genomic_assemblies is not None) and (mag_folder is None):
        genomic_seq_type = 'ctg'
        metagenomic_assemblies_file_path, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension = sep_path_basename_ext(genomic_assemblies)
        blast_db_dir = '%s/%s_%s_blast_db' % (working_directory, output_prefix, metagenomic_assemblies_file_basename)
        blast_db     = '%s/%s%s'           % (blast_db_dir, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension)
        os.mkdir(blast_db_dir)
        os.system('cp %s %s/' % (genomic_assemblies, blast_db_dir))

    elif (genomic_assemblies is None) and (mag_folder is not None):
        genomic_seq_type = 'mag'
        mag_folder_name     = mag_folder.split('/')[-1]
        blast_db_dir        = '%s/%s_blast_db'            % (working_directory, mag_folder_name)
        renamed_mag_folder  = '%s/%s_blast_db/%s_renamed' % (working_directory, mag_folder_name, mag_folder_name)
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
            rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename)

        # combine renamed MAGs
        blast_db = '%s/%s_combined.fa' % (blast_db_dir, mag_folder_name)
        os.system('cat %s/*%s > %s' % (renamed_mag_folder, mag_file_extension, blast_db))

    else:
        print('Please provide genomic sequences either as raw assemblies (-g) or as MAGs (-mag)')
        exit()


    ########################################### define folder and file name ############################################

    marker_gene_seqs_file_path, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension = sep_path_basename_ext(marker_gene_seqs)

    pwd_log_file                                = '%s/%s_%s.log'                                    % (working_directory, output_prefix, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))
    bowtie_index_dir                            = '%s/%s_%s_index'                                  % (working_directory, output_prefix, marker_gene_seqs_file_basename)
    pwd_samfile                                 = '%s/%s.sam'                                       % (working_directory, marker_gene_seqs_file_basename)
    clipping_reads_matched_part                 = '%s/%s_clipping_matched_part.txt'                 % (working_directory, output_prefix)
    clipping_reads_not_matched_part_seq         = '%s/%s_clipping_not_matched_part_seq.fasta'       % (working_directory, output_prefix)
    clipping_reads_not_matched_part_seq_blastn  = '%s/%s_clipping_not_matched_part_seq_blast.txt'   % (working_directory, output_prefix)
    clipping_reads_match_profile                = '%s/%s_match_profile_clipping.txt'                % (working_directory, output_prefix)
    unmapped_paired_reads_folder                = '%s/%s_unmapped_paired_reads'                     % (working_directory, output_prefix)
    unmapped_paired_reads_file                  = '%s/%s_unmapped_paired_reads.fasta'               % (working_directory, output_prefix)
    unmapped_paired_reads_blastn                = '%s/%s_unmapped_paired_reads_blast.txt'           % (working_directory, output_prefix)
    paired_reads_match_profile                  = '%s/%s_match_profile_paired.txt'                  % (working_directory, output_prefix)
    blast_results_all_vs_all_16s                = '%s/%s_16S_all_vs_all_blastn.tab'                 % (working_directory, output_prefix)
    link_stats_clipping                         = '%s/%s_stats_clipping.txt'                        % (working_directory, output_prefix)
    link_stats_clipping_filtered                = '%s/%s_stats_clipping_filtered.txt'               % (working_directory, output_prefix)
    link_stats_paired                           = '%s/%s_stats_paired.txt'                          % (working_directory, output_prefix)
    link_stats_paired_filtered                  = '%s/%s_stats_paired_filtered.txt'                 % (working_directory, output_prefix)
    link_stats_combined_table                   = '%s/%s_stats_combined_table.txt'                  % (working_directory, output_prefix)
    link_stats_combined                         = '%s/%s_stats_combined.txt'                        % (working_directory, output_prefix)
    pairwise_marker_similarity                  = '%s/%s_pairwise_marker_similarity.txt'            % (working_directory, output_prefix)

    blast_parameters    = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % num_threads


    ######################################## map reads to marker gene sequences ########################################

    # copy marker gene sequence file to index folder
    report_and_log(('Indexing marker gene sequences for mapping'), pwd_log_file, keep_quiet)
    os.mkdir(bowtie_index_dir)
    os.system('cp %s %s/' % (marker_gene_seqs, bowtie_index_dir))

    # run mapping
    report_and_log(('Mapping reads to marker gene sequences'), pwd_log_file, keep_quiet)
    sleep(1)
    report_and_log(('Please ignore warnings starting with "Use of uninitialized value" during Bowtie mapping.'), pwd_log_file, keep_quiet)
    sleep(1)

    bowtie2_index_ref_cmd = '%s -f %s/%s%s %s/%s --quiet --threads %s' % (pwd_bowtie2_build_exe, bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension, bowtie_index_dir, marker_gene_seqs_file_basename, num_threads)

    bowtie2_mapping_cmd = '%s -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s'       % (pwd_bowtie2_exe,  bowtie_index_dir, marker_gene_seqs_file_basename, reads_file_r1, reads_file_r2, pwd_samfile, num_threads)

    # if multiple_placement == '0':
    #     bowtie2_mapping_cmd = '%s -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s'       % (pwd_bowtie2_exe,  bowtie_index_dir, marker_gene_seqs_file_basename, reads_file_r1, reads_file_r2, pwd_samfile, num_threads)
    # elif multiple_placement == 'a':
    #     bowtie2_mapping_cmd = '%s -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal -a --quiet --threads %s'    % (pwd_bowtie2_exe,  bowtie_index_dir, marker_gene_seqs_file_basename, reads_file_r1, reads_file_r2, pwd_samfile, num_threads)
    # else:
    #     bowtie2_mapping_cmd = '%s -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal -k %s --quiet --threads %s' % (pwd_bowtie2_exe,  bowtie_index_dir, marker_gene_seqs_file_basename, reads_file_r1, reads_file_r2, pwd_samfile, multiple_placement, num_threads)

    bbmap_index_and_mapping_cmd = 'bbmap.sh ref=%s/%s%s in=%s in2=%s outm=%s local=t ambig=toss' % (bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension, reads_file_r1, reads_file_r2, pwd_samfile)


    if run_bbmap is True:
        os.system(bbmap_index_and_mapping_cmd)
    else:
        os.system(bowtie2_index_ref_cmd)
        os.system(bowtie2_mapping_cmd)

    sleep(1)
    report_and_log(('Please ignore warnings starting with "Use of uninitialized value" during Bowtie mapping.'), pwd_log_file, keep_quiet)
    sleep(1)


    ##################################################### extract reads ####################################################

    report_and_log(('Extracting unmapped part of clipping mapped reads from sam file'), pwd_log_file, keep_quiet)

    sam_file_to_parse = pwd_samfile


    # export clipping mapped reads and perfectly mapped reads
    clipping_mapped_reads_list = set()
    clipping_reads_mapped_part_dict = {}
    perfectly_mapped_reads_dict = {}
    clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
    for each_read in open(sam_file_to_parse):

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
    report_and_log(('Extracting unmapped paired reads with %s cores' % num_threads), pwd_log_file, keep_quiet)

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

    report_and_log(('Making blastn database'), pwd_log_file, keep_quiet)
    os.system(makeblastdb_cmd)

    report_and_log(('Running blastn for unmapped paired reads'), pwd_log_file, keep_quiet)
    os.system(blastn_cmd_paired)

    report_and_log(('Running blastn for unmapped parts of clipping mapped reads'), pwd_log_file, keep_quiet)
    os.system(blastn_cmd_clipping)


    ######################################### parse blast results for paired reads #########################################

    report_and_log(('Parsing blast results for paired reads'), pwd_log_file, keep_quiet)

    # filter blast results for paired reads
    unmapped_paired_reads_to_ctg_dict = blast_results_to_dict(unmapped_paired_reads_blastn, reads_iden_cutoff, reads_cov_cutoff)
    #unmapped_paired_reads_to_ctg_dict = blast_results_to_dict_by_bitscore(unmapped_paired_reads_blastn, reads_iden_cutoff)

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

    report_and_log(('Parsing blast results for clipping mapped reads'), pwd_log_file, keep_quiet)

    # filter blast results for clipping mapped reads
    clipping_parts_to_ctg_dict        = blast_results_to_dict(clipping_reads_not_matched_part_seq_blastn, reads_iden_cutoff, reads_cov_cutoff)
    #clipping_parts_to_ctg_dict        = blast_results_to_dict_by_bitscore(clipping_reads_not_matched_part_seq_blastn, reads_iden_cutoff)
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

    all_vs_all_16s_blastn_cmd = '%s -query %s -db %s -out %s %s'          % (pwd_blastn_exe, blastdb_16s, blastdb_16s, blast_results_all_vs_all_16s, blast_parameters)
    os.system(all_vs_all_16s_blastn_cmd)

    pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, aln16s, cov16s)

    # write out to file
    pairwise_marker_similarity_handle = open(pairwise_marker_similarity, 'w')
    pairwise_marker_similarity_handle.write('Marker1\tMarker2\tSimilarity\n')
    for marker_pair in pairwise_16s_iden_dict:
        pairwise_marker_similarity_handle.write('%s\t%s\n' % ('\t'.join(marker_pair.split('__|__')), pairwise_16s_iden_dict[marker_pair]))
    pairwise_marker_similarity_handle.close()


    ##################################################### get linkages #####################################################

    report_and_log(('Parsing linkages'), pwd_log_file, keep_quiet)

    # prepare input file for sankey plot
    stats_dict_to_sankey_file_in(clipping_stats_dict_num, paired_stats_dict_num, link_stats_clipping, link_stats_paired)

    # filter paired and clipping linkages
    filter_linkages_iteratively(link_stats_paired, 'Number', pairwise_16s_iden_dict, iden16s, link_stats_paired_filtered)
    filter_linkages_iteratively(link_stats_clipping, 'Number', pairwise_16s_iden_dict, iden16s, link_stats_clipping_filtered)

    # combine_paired_and_clipping_linkages and get summary table
    combine_paired_and_clipping_linkages(link_stats_paired_filtered, link_stats_clipping_filtered, link_stats_combined_table, link_stats_combined)


    ######################################################### plot #########################################################

    report_and_log(('Visualising linkages'), pwd_log_file, keep_quiet)

    # get plot height
    MarkerGenes_clipping = set()
    GenomicSeqs_clipping = set()
    for filtered_clipping in open(link_stats_clipping_filtered):
        if not filtered_clipping.startswith('MarkerGene,GenomicSeq,Number'):
            filtered_clipping_split = filtered_clipping.strip().split(',')
            MarkerGenes_clipping.add(filtered_clipping_split[0])
            GenomicSeqs_clipping.add(filtered_clipping_split[1])

    MarkerGenes_paired = set()
    GenomicSeqs_paired = set()
    for filtered_paired in open(link_stats_paired_filtered):
        if not filtered_paired.startswith('MarkerGene,GenomicSeq,Number'):
            filtered_paired_split = filtered_paired.strip().split(',')
            MarkerGenes_paired.add(filtered_paired_split[0])
            GenomicSeqs_paired.add(filtered_paired_split[1])

    MarkerGenes_combined = set()
    GenomicSeqs_combined = set()
    for filtered_combined in open(link_stats_combined):
        if not filtered_combined.startswith('MarkerGene,GenomicSeq,Number'):
            filtered_combined_split = filtered_combined.strip().split(',')
            MarkerGenes_combined.add(filtered_combined_split[0])
            GenomicSeqs_combined.add(filtered_combined_split[1])

    # calculate plot height
    plot_height_clipping     = 500 if max([len(MarkerGenes_clipping), len(GenomicSeqs_clipping)]) <= 25 else max([len(MarkerGenes_clipping), len(GenomicSeqs_clipping)]) * 20
    plot_height_paired       = 500 if max([len(MarkerGenes_paired), len(GenomicSeqs_paired)]) <= 25 else max([len(MarkerGenes_paired), len(GenomicSeqs_paired)]) * 20
    plot_height_intersection = 500 if max([len(MarkerGenes_combined), len(GenomicSeqs_combined)]) <= 25 else max([len(MarkerGenes_combined), len(GenomicSeqs_combined)]) * 20

    # prepare commands
    cmd_sankey_clipping     = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, link_stats_clipping_filtered, 600, plot_height_clipping)
    cmd_sankey_paired       = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, link_stats_paired_filtered,   600, plot_height_paired)
    cmd_sankey_intersection = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, link_stats_combined,          600, plot_height_intersection)

    # plot
    if len(open(link_stats_clipping_filtered).readlines()) == 1:
        report_and_log(('No data in %s, plotting skipped' % (link_stats_clipping_filtered)), pwd_log_file, keep_quiet)
    else:
        os.system(cmd_sankey_clipping)

    if len(open(link_stats_paired_filtered).readlines()) == 1:
        report_and_log(('No data in %s, plotting skipped' % (link_stats_paired_filtered)), pwd_log_file, keep_quiet)
    else:
        os.system(cmd_sankey_paired)

    if len(open(link_stats_combined).readlines()) == 1:
        report_and_log(('No data in %s, plotting skipped' % (link_stats_combined)), pwd_log_file, keep_quiet)
    else:
        os.system(cmd_sankey_intersection)


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
        recovery_paired, accuracy_paired, recovered_paired       = get_accuracy(link_stats_paired_filtered, len(marker_id_set))
        recovery_clipping, accuracy_clipping, recovered_clipping = get_accuracy(link_stats_clipping_filtered, len(marker_id_set))

        # get unrecovered markers
        unrecovered_markers_paired       = get_unrecovered_markers(marker_id_set, recovered_paired)
        unrecovered_markers_clipping     = get_unrecovered_markers(marker_id_set, recovered_clipping)
        unrecovered_markers_paired_str   = 'Unrecovered_Paired(%s):%s'   % (len(unrecovered_markers_paired), ','.join(sorted([i for i in unrecovered_markers_paired])))
        unrecovered_markers_clipping_str = 'Unrecovered_Clipping(%s):%s' % (len(unrecovered_markers_clipping), ','.join(sorted([i for i in unrecovered_markers_clipping])))

        # assessment by genome
        assign_rate_paired, assign_accuracy_paired, right_assign_paired, wrong_assign_paired         = get_accuracy_by_genome(link_stats_paired_filtered, mag_folder, mag_file_extension)
        assign_rate_clipping, assign_accuracy_clipping, right_assign_clipping, wrong_assign_clipping = get_accuracy_by_genome(link_stats_clipping_filtered, mag_folder, mag_file_extension)
        unrecovered_paired_report_str   = 'Unrecovered_Paired(%s):%s'   % (len(wrong_assign_paired), ','.join(sorted([i for i in wrong_assign_paired])))
        unrecovered_clipping_report_str = 'Unrecovered_Clipping(%s):%s' % (len(wrong_assign_clipping), ','.join(sorted([i for i in wrong_assign_clipping])))

        # report
        report_and_log(('Prefix\tBy\tPaired_r\tPaired_a\tClip_r\tClip_a\tUnrecovered_Paired\tUnrecovered_Clipping'), pwd_log_file, keep_quiet)
        report_and_log(('%s\tMarker\t%s\t%s\t%s\t%s\t%s\t%s' % (output_prefix, recovery_paired, accuracy_paired, recovery_clipping, accuracy_clipping, unrecovered_markers_paired_str, unrecovered_markers_clipping_str)), pwd_log_file, keep_quiet)
        report_and_log(('%s\tGenome\t%s\t%s\t%s\t%s\t%s\t%s' % (output_prefix, assign_rate_paired, assign_accuracy_paired, assign_rate_clipping, assign_accuracy_clipping, unrecovered_paired_report_str, unrecovered_clipping_report_str)), pwd_log_file, keep_quiet)


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

    link_Marker_MAG_parser = argparse.ArgumentParser(description='Link MAGs with marker genes', usage=link_Marker_MAG_usage)

    link_Marker_MAG_parser.add_argument('-p',               required=True,                              help='output prefix')
    link_Marker_MAG_parser.add_argument('-r1',              required=True,                              help='paired reads r1')
    link_Marker_MAG_parser.add_argument('-r2',              required=True,                              help='paired reads r2')
    link_Marker_MAG_parser.add_argument('-m',               required=True,                              help='marker gene sequences')
    link_Marker_MAG_parser.add_argument('-g',               required=False, default=None,               help='genomic sequences')
    link_Marker_MAG_parser.add_argument('-mag',             required=False, default=None,               help='metagenome-assembled-genome (MAG) folder')
    link_Marker_MAG_parser.add_argument('-x',               required=False, default='fasta',            help='MAG file extension, default: fasta')
    link_Marker_MAG_parser.add_argument('-mp',              required=False, type=str, default=0,        help='multiple placement during bowtie mapping, default: 0. choose from a,1,2,3...')
    link_Marker_MAG_parser.add_argument('-cigarM',          required=False, type=int, default=30,       help='cigarM cutoff, default: 30')
    link_Marker_MAG_parser.add_argument('-cigarS',          required=False, type=int, default=30,       help='cigarS cutoff, default: 30')
    link_Marker_MAG_parser.add_argument('-ri',              required=False, type=float, default=100,    help='identity cutoff, default: 100')
    link_Marker_MAG_parser.add_argument('-rc',              required=False, type=float, default=100,    help='coverage cutoff, default: 100')
    link_Marker_MAG_parser.add_argument('-mi',              required=False, type=float, default=99.5,   help='within genome 16S identity cutoff, default: 99.5')
    link_Marker_MAG_parser.add_argument('-mc',              required=False, type=float, default=90,     help='alignment coverage cutoff for calculating 16S identity, default: 90')
    link_Marker_MAG_parser.add_argument('-ma',              required=False, type=int, default=500,      help='alignment length cutoff for calculating 16S identity, default: 500')
    link_Marker_MAG_parser.add_argument('-t',               required=False, type=int, default=1,        help='number of threads, default: 1')
    link_Marker_MAG_parser.add_argument('-quiet',           required=False, action="store_true",        help='not report progress')
    link_Marker_MAG_parser.add_argument('-force',           required=False, action="store_true",        help='force overwrite existing results')
    link_Marker_MAG_parser.add_argument('-tmp',             required=False, action="store_true",        help='keep temporary files')
    link_Marker_MAG_parser.add_argument('-test_mode',       required=False, action="store_true",        help='only for debugging, do not provide')
    link_Marker_MAG_parser.add_argument('-bbmap',           required=False, action="store_true",        help='run bbmap, instead of bowtie')

    args = vars(link_Marker_MAG_parser.parse_args())

    link_Marker_MAG(args, config_dict)


Notes = '''

alias twine='/Users/songweizhi/Library/Python/3.7/bin/twine'
cd /Users/songweizhi/PycharmProjects/MarkerMAG
rm -r build
rm -r dist
rm -r MarkerMAG.egg-info
python setup.py sdist bdist_wheel

twine upload dist/*
songweizhi

twine upload --repository-url https://test.pypi.org/legacy/ dist/*
songweizhi

shan88

pip3 install --upgrade MarkerMAG
pip3 install --upgrade -i https://test.pypi.org/simple/ MarkerMAG

'''

To_do = '''

2. where does the paired read of the clipping mapped read mapped to? (should take into consideration!!!)
4. how to incorporate the taxonomy of MAGs and 16S sequences
5. the effect of sequencing depth, insert size and read length
6. the depth of 16S sequences always not lower than the genome they come from
7. split sam file
8. with no_ambiguous option, 16S rRNA gene sequences need to be dereplicated. (include dereplication step? with identity and coverage cutoffs?)
9.check whether input file exist!

Notes:
linkages only supported by clipping reads were ignored !!!


export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH
cd /Users/songweizhi/Desktop/MarkerMAG_wd/Mac_test
~/PycharmProjects/MarkerMAG/bin/MarkerMAG -p EvenDepth_5x -r1 combined_5x_R1.fasta -r2 combined_5x_R2.fasta -m combined_16S.ffn -mag selected_genomes_renamed_no_plasmid -x fna -t 4 -tmp -force -test_mode


module load python/3.7.3
source ~/mypython3env/bin/activate
module load R/4.0.2
module load java/8u121
module load bbmap/38.51
module load blast+/2.9.0
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/Mac_test
./MarkerMAG -p EvenDepth_5x_bbmap -r1 combined_5x_R1.fasta -r2 combined_5x_R2.fasta -m combined_16S.ffn -mag selected_genomes_renamed_no_plasmid -x fna -t 12 -tmp -force -test_mode


s1_00644	s2_00537	99.282
p5_00515	p5_03129	99.934
p5_00515	p5_02332	99.213
p5_00515	p5_01648	99.213
p5_00515	p5_01238	99.213
p5_00515	p5_00656	99.082
p5_00515	p5_00944	99.016

'''
