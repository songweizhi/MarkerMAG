import os
import glob
import numpy as np
import pandas as pd
from Bio import SeqIO


def sep_path_basename_ext(file_in):
    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def sort_csv_by_col(file_in, file_out, col_header):
    df_in = pd.read_csv(file_in)
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

        if (align_len >= align_len_cutoff) and (query != subject) and (coverage_q >= cov_cutoff) and (
                coverage_s >= cov_cutoff):
            query_to_subject_key = '__|__'.join(sorted([query, subject]))
            if query_to_subject_key not in pairwise_iden_dict:
                pairwise_iden_dict[query_to_subject_key] = iden
            else:
                if iden > pairwise_iden_dict[query_to_subject_key]:
                    pairwise_iden_dict[query_to_subject_key] = iden

    return pairwise_iden_dict


# def filter_linkages_iteratively_backup(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict,
#                                 marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff,
#                                 min_linkages, min_linkages_for_uniq_linked_16s, within_gnm_linkage_num_diff, file_out):
#
#     # get MarkerGene_to_GenomicSeq_dict
#     MarkerGene_to_GenomicSeq_dict = {}
#     for each_linkage in open(file_in):
#         if not each_linkage.startswith('MarkerGene,GenomicSeq,Number'):
#             each_linkage_split = each_linkage.strip().split(',')
#             MarkerGene_id = each_linkage_split[0][12:]
#             GenomicSeq_id = each_linkage_split[1][12:]
#             linkage_num = int(each_linkage_split[2])
#             if linkage_num > 1:
#                 if MarkerGene_id not in MarkerGene_to_GenomicSeq_dict:
#                     MarkerGene_to_GenomicSeq_dict[MarkerGene_id] = {GenomicSeq_id}
#                 else:
#                     MarkerGene_to_GenomicSeq_dict[MarkerGene_id].add(GenomicSeq_id)
#
#     print('MarkerGene_to_GenomicSeq_dict')
#     print(MarkerGene_to_GenomicSeq_dict)
#
#     file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
#     file_in_sorted = '%s/%s_sorted%s' % (file_in_path, file_in_basename, file_in_extension)
#
#     # sort file in
#     sort_csv_by_col(file_in, file_in_sorted, sort_by_col_header)
#
#     # fileter linkage
#     gnm_max_link_num_dict = {}
#     file_out_handle = open(file_out, 'w')
#     MarkerGene_with_assignment = set()
#     GenomicSeq_best_marker_dict = {}
#     for each_match in open(file_in_sorted):
#         if each_match.startswith('MarkerGene,GenomicSeq,Number'):
#             file_out_handle.write(each_match)
#         else:
#             match_split = each_match.strip().split(',')
#             MarkerGene = match_split[0][12:]
#             GenomicSeq = match_split[1][12:]
#             linkage_num = int(match_split[2])
#
#             current_min_linkage = min_linkages_for_uniq_linked_16s
#             if MarkerGene in MarkerGene_to_GenomicSeq_dict:
#                 if len(MarkerGene_to_GenomicSeq_dict[MarkerGene]) > 1:
#                     current_min_linkage = min_linkages
#
#             if linkage_num >= current_min_linkage:
#                 if MarkerGene not in MarkerGene_with_assignment:
#
#                     # consider depth
#                     if min_16s_gnm_multiple > 0:
#
#                         # get marker and genome depth
#                         MarkerGene_depth = marker_gene_depth_dict.get(MarkerGene, 'na')
#                         GenomicSeq_depth = genomic_seq_depth_dict.get(GenomicSeq, 'na')
#                         marker_genome_depth_ratio = 'na'
#                         if (MarkerGene_depth != 'na') and (GenomicSeq_depth != 'na'):
#                             if GenomicSeq_depth > 0:
#                                 marker_genome_depth_ratio = MarkerGene_depth / GenomicSeq_depth
#
#                         if marker_genome_depth_ratio >= min_16s_gnm_multiple:
#                             if GenomicSeq not in GenomicSeq_best_marker_dict:
#                                 GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
#                                 gnm_max_link_num_dict[GenomicSeq] = linkage_num
#                                 file_out_handle.write(each_match)
#                                 MarkerGene_with_assignment.add(MarkerGene)
#                             else:
#                                 # get identity with best marker
#                                 current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
#                                 key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))
#                                 iden_with_best_marker = pairwise_16s_iden_dict.get(key_str, 0)
#                                 if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
#                                     gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
#                                     if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
#                                         file_out_handle.write(each_match)
#                                         MarkerGene_with_assignment.add(MarkerGene)
#                                     else:
#                                         MarkerGene_with_assignment.add(MarkerGene)
#                     # ignore depth
#                     else:
#                         if GenomicSeq not in GenomicSeq_best_marker_dict:
#                             GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
#                             gnm_max_link_num_dict[GenomicSeq] = linkage_num
#                             file_out_handle.write(each_match)
#                             MarkerGene_with_assignment.add(MarkerGene)
#                         else:
#                             # get identity with best marker
#                             current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
#                             key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))
#                             iden_with_best_marker = pairwise_16s_iden_dict.get(key_str, 0)
#                             if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
#                                 gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
#                                 if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
#                                     file_out_handle.write(each_match)
#                                     MarkerGene_with_assignment.add(MarkerGene)
#                                 else:
#                                     MarkerGene_with_assignment.add(MarkerGene)
#     file_out_handle.close()

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


def filter_linkages_iteratively(sorted_file_in, pairwise_16s_iden_dict, within_genome_16s_divergence_cutoff, marker_len_dict,
                                min_linkages, within_gnm_linkage_num_diff, file_out):

    # fileter linkage
    gnm_max_link_num_dict = {}
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
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
            GenomicSeq = match_split[1][12:]
            linkage_num = int(match_split[2])

            if linkage_num >= min_linkages:
                if MarkerGene in MarkerGene_with_assignment:
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
                                    # check length, if shorter than assigned, ignore this one by adding it to MarkerGene_with_assignment
                                    MarkerGene_len = marker_len_dict[MarkerGene]
                                    already_assigned_16s_len_list = [marker_len_dict[s16] for s16 in already_assigned_16s_list]
                                    median_assign_16s_len = np.median(already_assigned_16s_len_list)
                                    if (median_assign_16s_len - MarkerGene_len) > 50:
                                        MarkerGene_with_assignment.add(MarkerGene)
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
                                    MarkerGene_len = marker_len_dict[MarkerGene]
                                    already_assigned_16s_len_list = [marker_len_dict[s16] for s16 in already_assigned_16s_list]
                                    median_assign_16s_len = np.median(already_assigned_16s_len_list)
                                    if (median_assign_16s_len - MarkerGene_len) > 50:
                                        MarkerGene_with_assignment.add(MarkerGene)
                        else:
                            current_gnm = GenomicSeq
                            current_gnm_best_16s_list = [MarkerGene]
                            current_gnm_highest_link_num = linkage_num

    #print('current_gnm: %s' % current_gnm)
    #print('current_gnm_highest_link_num: %s' % current_gnm_highest_link_num)
    #print('current_gnm_best_16s_list: %s' % current_gnm_best_16s_list)

    file_out_handle.close()


file_id = '45'
input_16s_polished                  = '/Users/songweizhi/Desktop/test_filter/GI_128_16S_0.999.fasta'
blast_results_all_vs_all_16s        = '/Users/songweizhi/Desktop/test_filter/GI_0524_128_60_60_16S_all_vs_all_blastn.tab'
min_aln_16s                         = 350
min_cov_16s                         = 30
mean_depth_dict_gnm                 = {}
mean_depth_dict_16s                 = {}
min_16s_gnm_multiple                = 0
min_iden_16s                        = 98.5
min_link_num                        = 8
within_gnm_linkage_num_diff         = 80
link_stats_combined                 = '/Users/songweizhi/Desktop/test_filter/GI_0524_128_%s_%s_stats_combined.txt'          % (file_id, file_id)
link_stats_combined_sorted          = '/Users/songweizhi/Desktop/test_filter/GI_0524_128_%s_%s_stats_combined_sorted.txt'   % (file_id, file_id)
link_stats_combined_sorted_filtered = '/Users/songweizhi/Desktop/test_filter/GI_0524_128_%s_%s_stats_combined_filtered.txt' % (file_id, file_id)


# get marker len dict
marker_len_dict = {}
for each_marker_record in SeqIO.parse(input_16s_polished, 'fasta'):
    marker_len_dict[each_marker_record.id] = len(each_marker_record.seq)

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

# sort file in
sort_csv_by_col(link_stats_combined, link_stats_combined_sorted, 'Number')

filter_linkages_iteratively(link_stats_combined_sorted, pairwise_16s_iden_dict, min_iden_16s, marker_len_dict,
                            min_link_num, within_gnm_linkage_num_diff, link_stats_combined_sorted_filtered)




#print(pairwise_16s_iden_dict['3_GI_subsample_100_2289__|__3_GI_subsample_100_2289'])