import os
import pandas as pd


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


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


def sort_csv_by_col(file_in, file_out, col_header):
    df_in        = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


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

                            # store max link num
                            gnm_max_link_num_dict[GenomicSeq] = linkage_num

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


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict,
                                marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff,
                                min_linkages, min_linkages_for_uniq_linked_16s, file_out):

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
                            if (linkage_num*100/gnm_max_link_num) >= 80:
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)
                            else:
                                MarkerGene_with_assignment.add(MarkerGene)

    file_out_handle.close()




'''
	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	Accuracy
raw	Rd_1	|	140	119	9	12	119/131(90.84)	|	32	29	0	0	3	29/97(29.9)	    29/32(90.62)
33	Rd_1	|	86	82	2	2	82/84(97.62)	|	32	30	0	0	2	30/97(30.93)	30/32(93.75)
50	Rd_1	|	73	69	1	3	69/72(95.83)	|	34	31	0	2	1	31/97(31.96)	31/34(91.18)
75  Rd_1	|	55	53	1	1	53/54(98.15)	|	34	33	0	1	0	33/97(34.02)	33/34(97.06)
85	Rd_1	|	46	44	1	1	44/45(97.78)	|	34	33	0	1	0	33/97(34.02)	33/34(97.06)
90	Rd_1	|	43	41	1	1	41/42(97.62)	|	34	33	0	1	0	33/97(34.02)	33/34(97.06)


33  Rd_1	|	83	79	2	2	79/81(97.53)	|	32	30	0	0	2	30/97(30.93)	30/32(93.75)
50  Rd_1	|	71	69	1	1	69/70(98.57)	|	32	31	0	0	1	31/97(31.96)	31/32(96.88)
75  Rd_1	|	51	50	1	0	50/50(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
80	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
85	Rd_1	|	44	43	1	0	43/43(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
90  Rd_1	|	41	40	1	0	40/40(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)

'''

step_1_wd                       = '/Users/songweizhi/Desktop/step_1_wd'
link_stats_combined             = '%s/stats_combined.txt'               % step_1_wd
link_stats_combined_filtered_s1 = '%s/stats_combined_filtered.txt'      % step_1_wd
blast_results_all_vs_all_16s    = '%s/16S_all_vs_all_blastn.tab' % step_1_wd
min_aln_16s                     = 500
min_cov_16s                     = 30
mean_depth_dict_gnm             = {}
mean_depth_dict_16s             = {}
min_16s_gnm_multiple            = 0
min_iden_16s                    = 98
min_link_num_rd1                = 10
min_uniq_link_num_rd1           = 7

combined_linkage_file = '%s/Test_identified_linkages_genome_level.txt' % step_1_wd

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)


filter_linkages_iteratively(link_stats_combined, 'Number', pairwise_16s_iden_dict,
                            mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple,
                            min_iden_16s, min_link_num_rd1,
                            min_uniq_link_num_rd1, link_stats_combined_filtered_s1)

combined_linkage_file_handle = open(combined_linkage_file, 'w')
combined_linkage_file_handle.write('MarkerGene\tGenomicSeq\tLinkage\tStep\n')
for step_1_link in open(link_stats_combined_filtered_s1):
    if not step_1_link.startswith('MarkerGene,GenomicSeq,Number'):
        marker_id = step_1_link.strip().split(',')[0][12:]
        genome_id = step_1_link.strip().split(',')[1][12:]
        link_num = step_1_link.strip().split(',')[2]
        combined_linkage_file_handle.write('%s\t%s\t%s\tS1\n' % (marker_id, genome_id, link_num))
combined_linkage_file_handle.close()

os.system('python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/GI/assess_linkages_GI.py')

