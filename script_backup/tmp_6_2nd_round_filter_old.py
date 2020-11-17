
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


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_max_value_in_dict(dict_in):
    max_value = 0
    max_value_key = ''
    for each_key in dict_in:
        if dict_in[each_key] > max_value:
            max_value = dict_in[each_key]
            max_value_key = each_key
    return max_value, max_value_key


def get_total_value_in_dict(dict_in):
    total_value = 0
    for each_key in dict_in:
        total_value += dict_in[each_key]
    return total_value


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


def filter_linkages_iteratively_simple(file_in, sort_by_col_header, min_linkages, file_out):

    file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
    file_in_sorted = '%s/%s_sorted%s' % (file_in_path, file_in_basename, file_in_extension)

    # sort file in
    sort_csv_by_col(file_in, file_in_sorted, sort_by_col_header)

    # fileter linkage
    file_out_handle = open(file_out, 'w')
    gnm_ctg_16s_with_assignment = set()
    for each_match in open(file_in_sorted):
        if each_match.startswith('Gap_seq,'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            gnm_ctg_16s_id = match_split[1]
            linkage_num = int(match_split[2])
            if (linkage_num >= min_linkages) and (gnm_ctg_16s_id not in gnm_ctg_16s_with_assignment):
                file_out_handle.write(each_match)
                gnm_ctg_16s_with_assignment.add(gnm_ctg_16s_id)
    file_out_handle.close()

    # remove tmp file
    # os.remove(file_in_sorted)


################################################### files/parameters ###################################################

wd = '/Users/songweizhi/Desktop/GapFilling_all'
wd = '/Users/songweizhi/Desktop/GapFilling_sep'

free_living_mate_gnm = '%s/free_living_mate_ctg.txt'                                        % wd
free_living_mate_16s = '%s/free_living_mate_16s.txt'                                        % wd
sam_file             = '%s/scaffolds.sam'                                                   % wd
gnm_ctg_connector    = '___'
min_read_num         = 10

if wd.endswith('_all'):
    blastdb_16s = '%s/MBARC26_all_depth_assemblies.dereplicated_99.5_500bp.fasta'           % wd
else:
    blastdb_16s = '%s/rest_16S.fasta'                                                       % wd


########################################################################################################################

# define output file name
blast_results_all_vs_all_16s   = '%s/16s_all_vs_all.txt'                                    % wd
stats_GapFilling_file          = '%s/stats_GapFilling.txt'                                  % wd
stats_GapFilling_file_filtered = '%s/stats_GapFilling_filtered.txt'                         % wd

stats_GapFilling_file_16s          = '%s/stats_GapFilling_16s.txt'                              % wd
stats_GapFilling_file_ctg          = '%s/stats_GapFilling_ctg.txt'                              % wd
stats_GapFilling_file_gnm          = '%s/stats_GapFilling_gnm.txt'                              % wd

stats_GapFilling_file_filtered_16s = '%s/stats_GapFilling_16s_filtered.txt'                         % wd
stats_GapFilling_file_filtered_ctg = '%s/stats_GapFilling_ctg_filtered.txt'                         % wd
stats_GapFilling_file_filtered_gnm = '%s/stats_GapFilling_gnm_filtered.txt'                         % wd


########################################################################################################################

reads_to_extract_to_ref_dict_gnm = {}
for i in open(free_living_mate_gnm):
    i_split = i.strip().split('\t')
    reads_to_extract_to_ref_dict_gnm[i_split[0]] = i_split[1]

reads_to_extract_to_ref_dict_16s = {}
for j in open(free_living_mate_16s):
    j_split = j.strip().split('\t')
    reads_to_extract_to_ref_dict_16s[j_split[0]] = j_split[1]

gap_seq_to_reads_dict = {}
for each_line in open(sam_file):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id         = each_line_split[0]
        read_id_base    = '.'.join(read_id.split('.')[:-1])
        read_strand     = read_id.split('.')[-1]
        ref_id          = each_line_split[2]
        ref_pos         = int(each_line_split[3])
        cigar           = each_line_split[5]
        cigar_splitted  = cigar_splitter(cigar)
        if (len(cigar_splitted) == 1) and (cigar[-1] == 'M'):
            if ref_id not in gap_seq_to_reads_dict:
                gap_seq_to_reads_dict[ref_id] = [read_id]
            else:
                gap_seq_to_reads_dict[ref_id].append(read_id)


stats_GapFilling_file_handle = open(stats_GapFilling_file, 'w')
stats_GapFilling_file_handle.write('MarkerGene,GenomicSeq,Number\n')
for short_ctg in gap_seq_to_reads_dict:
    short_ctg_mapped_reads = gap_seq_to_reads_dict[short_ctg]
    short_ctg_mapped_reads_mate_linkage_dict_16s = {}
    short_ctg_mapped_reads_mate_linkage_dict_ctg = {}
    marker_mate = False
    contig_mate = False
    for mapped_read in short_ctg_mapped_reads:

        mapped_read_mate_ref = ''
        if mapped_read in reads_to_extract_to_ref_dict_gnm:
            mapped_read_mate_ref = reads_to_extract_to_ref_dict_gnm[mapped_read]
        if mapped_read in reads_to_extract_to_ref_dict_16s:
            mapped_read_mate_ref = reads_to_extract_to_ref_dict_16s[mapped_read]

        # check mate mapped to marker gene or genomic sequence
        if '_m' in mapped_read_mate_ref:
            marker_mate = True

            if mapped_read_mate_ref not in short_ctg_mapped_reads_mate_linkage_dict_16s:
                short_ctg_mapped_reads_mate_linkage_dict_16s[mapped_read_mate_ref] = 1
            else:
                short_ctg_mapped_reads_mate_linkage_dict_16s[mapped_read_mate_ref] += 1

        mapped_read_mate_ref_genome = ''
        if 'NODE' in mapped_read_mate_ref:
            contig_mate = True
            mapped_read_mate_ref_genome = mapped_read_mate_ref.split(gnm_ctg_connector)[0]
            mapped_read_mate_ref = mapped_read_mate_ref_genome

            if mapped_read_mate_ref not in short_ctg_mapped_reads_mate_linkage_dict_ctg:
                short_ctg_mapped_reads_mate_linkage_dict_ctg[mapped_read_mate_ref] = 1
            else:
                short_ctg_mapped_reads_mate_linkage_dict_ctg[mapped_read_mate_ref] += 1

    if (len(short_ctg_mapped_reads_mate_linkage_dict_ctg) > 0) and (len(short_ctg_mapped_reads_mate_linkage_dict_16s) > 0):
        print(short_ctg)
        print(short_ctg_mapped_reads_mate_linkage_dict_ctg)
        print(short_ctg_mapped_reads_mate_linkage_dict_16s)
        print()
        for each_16s in short_ctg_mapped_reads_mate_linkage_dict_16s:
            num_16s = short_ctg_mapped_reads_mate_linkage_dict_16s[each_16s]
            for each_ctg in short_ctg_mapped_reads_mate_linkage_dict_ctg:
                num_ctg = short_ctg_mapped_reads_mate_linkage_dict_ctg[each_ctg]

                if (num_16s >= min_read_num) and (num_ctg >= min_read_num):
                    stats_GapFilling_file_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (each_16s, each_ctg, (num_16s + num_ctg)))

stats_GapFilling_file_handle.close()


############################################# get pairwise_16s_iden_dict ##############################################

pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
num_threads = 4
blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % num_threads

# makeblastdn with marker gene sequences
makeblastdb_16s_cmd = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blastdb_16s)
os.system(makeblastdb_16s_cmd)

all_vs_all_16s_blastn_cmd = '%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, blastdb_16s, blastdb_16s, blast_results_all_vs_all_16s, blast_parameters)
os.system(all_vs_all_16s_blastn_cmd)

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, 500, 90)


#######################################################################################################################

filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, {}, {}, 0, 98, 0, stats_GapFilling_file_filtered)
