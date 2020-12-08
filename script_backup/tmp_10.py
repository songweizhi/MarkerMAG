
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


unmapped_paired_reads_blastn    = '/Users/songweizhi/Desktop/unmapped_paired_reads_blast.txt'
reads_iden_cutoff               = 100  # s1_ri
reads_cov_cutoff                = 90   # s1_rc


# filter blast results for paired reads
unmapped_paired_reads_to_ctg_dict = paired_blast_results_to_dict(unmapped_paired_reads_blastn, reads_iden_cutoff, reads_cov_cutoff)
print(unmapped_paired_reads_to_ctg_dict)
print(len(unmapped_paired_reads_to_ctg_dict))




