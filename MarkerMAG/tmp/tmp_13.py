
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


blast_results_all_vs_all_16s = '/Users/songweizhi/Desktop/GI_0524_128_60_60_16S_all_vs_all_blastn.tab'
min_aln_16s = 500
min_cov_16s = 30

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)
highest_link_16s_list = ['3_GI_subsample_100_1336', '3_GI_subsample_75_2502', '3_GI_subsample_25_763', '3_GI_subsample_100_2066', '3_GI_subsample_75_69', '3_GI_subsample_100_3364', '3_GI_subsample_75_1008', '3_GI_subsample_25_565',
                         '3_GI_subsample_25_345', '3_GI_subsample_10_173', '3_GI_subsample_75_988',
                         '3_GI_subsample_75_1644']


def get_mean_iden_list(linked_16s_list):

    mean_iden_list = []
    for each_16s_1 in linked_16s_list:
        current_16s_iden_list = []
        for each_16s_2 in linked_16s_list:
            if each_16s_1 != each_16s_2:
                key = '__|__'.join(sorted([each_16s_1, each_16s_2]))
                key_iden = pairwise_16s_iden_dict.get(key, 'na')
                if key_iden != 'na':
                    current_16s_iden_list.append(key_iden)
        mean_iden = sum(current_16s_iden_list) / len(current_16s_iden_list)
        mean_iden = float("{0:.3f}".format(mean_iden))
        # print('%s\t%s\t%s' % (each_16s_1, mean_iden, current_16s_iden_list))
        mean_iden_list.append(mean_iden)

    return mean_iden_list


mean_iden_list = get_mean_iden_list(highest_link_16s_list)

sorted_seq_id_list = [seq_id for mean_iden, seq_id in sorted(zip(mean_iden_list, highest_link_16s_list), reverse=True)]

print(highest_link_16s_list)

print(sorted(mean_iden_list, reverse=True))

print(sorted_seq_id_list)