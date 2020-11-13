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
            MarkerGene_depth = marker_gene_depth_dict[MarkerGene]
            GenomicSeq_depth = genomic_seq_depth_dict[GenomicSeq]


            linkage_num = int(match_split[2])
            if linkage_num >= min_linkages:
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
                                print('%s(%s)\t%s(%s)' % (GenomicSeq, GenomicSeq_depth, MarkerGene, MarkerGene_depth))

    file_out_handle.close()

    # remove tmp file
    # os.remove(file_in_sorted)


wd = '/Users/songweizhi/Desktop/test'
blast_results_all_vs_all_16s    = '%s/Test_16S_all_vs_all_blastn.tab' % wd
link_stats_paired               = '%s/Test_stats_paired.txt' % wd
link_stats_paired_filtered      = '%s/MBARC26_link_Matam_FakeBins_stats_paired_filtered.txt' % wd
within_genome_minimum_iden16s   = 98
aln16s                          = 500
cov16s                          = 30
mean_depth_dict_ctg             = {'c1___CP018800.1': 5.031, 'c2___CP020538.1': 5.031, 'c3___CP000471.1': 5.036, 'f1___CP007152.1': 5.032, 'f2___CP007577.1': 5.037, 'f3___CP019343.1': 5.033, 'f4___AM286690.1': 5.036, 'f5___FN869568.2': 5.033, 'g1___CM001020.1': 5.01, 'g2___CCSF01000001.1': 5.031, 'g3___AP014862.1': 5.034, 'g4___CP012358.1': 5.037, 'g5___CP005094.1': 5.032, 'o1___CP005990.1': 5.035, 'o2___CP020046.1': 5.032, 'o3___LN614827.1': 5.032, 'o4___CP011849.2': 5.037, 'o5___LR133904.1': 5.035, 'p1___HG937516.1': 5.025, 'p2___AE015925.1': 5.032, 'p3___LS483446.1': 5.031, 'p4___LT607413.1': 5.034, 'p5___LR590484.1': 5.034, 's1___CP015638.1': 5.033, 's2___CP020369.1': 5.028, 's3___CP003190.1': 5.033, 's4___CP009533.1': 5.033, 's5___CT573326.1': 5.036}
mean_depth_dict_gnm             = {'c1': 5.031, 'c2': 5.031, 'c3': 5.036, 'f1': 5.032, 'f2': 5.037, 'f3': 5.033, 'f4': 5.036, 'f5': 5.033, 'g1': 5.01, 'g2': 5.031, 'g3': 5.034, 'g4': 5.037, 'g5': 5.032, 'o1': 5.035, 'o2': 5.032, 'o3': 5.032, 'o4': 5.037, 'o5': 5.035, 'p1': 5.025, 'p2': 5.032, 'p3': 5.031, 'p4': 5.034, 'p5': 5.034, 's1': 5.033, 's2': 5.028, 's3': 5.033, 's4': 5.033, 's5': 5.036}
mean_depth_dict_16s             = {'c1_01889': 5.559, 'c2_01869': 5.365, 'c3_00259': 5.702, 'c3_02816': 4.581, 'c3_03081': 3.14, 'f1_05369': 4.339, 'f2_03131': 7.324, 'f3_00209': 5.765, 'f3_02041': 5.1, 'f4_00381': 4.894, 'f4_00512': 5.673, 'f4_02069': 5.216, 'f5_00225': 4.438, 'f5_00355': 4.755, 'f5_00795': 4.896, 'f5_03420': 6.314, 'g1_00650': 4.586, 'g2_03471': 4.81, 'g3_00108': 4.298, 'g3_00113': 5.179, 'g3_00537': 4.558, 'g3_02127': 6.094, 'g3_04871': 4.962, 'g4_00414': 4.647, 'g4_01091': 4.379, 'g4_01774': 5.02, 'g4_01998': 4.27, 'g4_02419': 5.728, 'g5_00174': 3.98, 'g5_01291': 4.308, 'g5_01784': 6.012, 'g5_02575': 4.9, 'g5_03966': 4.321, 'g5_04304': 4.54, 'o1_01731': 4.095, 'o2_01394': 5.094, 'o3_00282': 4.952, 'o3_00426': 4.953, 'o3_02626': 4.888, 'o3_03185': 5.707, 'o4_00099': 4.168, 'o4_00249': 4.583, 'o4_01748': 4.685, 'o4_02347': 4.636, 'o4_02542': 4.985, 'o4_03086': 4.834, 'o5_00048': 1.59, 'o5_00192': 1.259, 'o5_00255': 1.0, 'o5_03318': 3.657, 'o5_03773': 1.431, 'o5_04011': 2.802, 'o5_04162': 2.015, 'p1_01191': 4.07, 'p2_00901': 4.42, 'p3_00101': 4.435, 'p3_00509': 5.044, 'p4_02304': 3.907, 'p4_03395': 4.732, 'p4_06318': 4.663, 'p5_00515': 5.596, 'p5_00656': 4.281, 'p5_00944': 4.996, 'p5_01238': 6.61, 'p5_01648': 5.721, 'p5_02332': 5.411, 'p5_03129': 6.113, 's1_00644': 6.008, 's2_00537': 3.243, 's3_00124': 5.04, 's3_00752': 4.497, 's3_04232': 4.854, 's3_05244': 4.176, 's3_05583': 4.274, 's4_00580': 4.968, 's4_02579': 4.221, 's4_02679': 5.444, 's4_03140': 4.815, 's4_03687': 5.628, 's4_04216': 3.589, 's5_00111': 5.295, 's5_00430': 5.011, 's5_00636': 4.627, 's5_01224': 5.174, 's5_02435': 6.221, 's5_03277': 5.956, 's5_04377': 5.391}
min_16s_gnm_multiple            = 1

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, aln16s, cov16s)

filter_linkages_iteratively(link_stats_paired, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, within_genome_minimum_iden16s, 2, link_stats_paired_filtered)


