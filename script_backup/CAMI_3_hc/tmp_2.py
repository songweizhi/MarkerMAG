
# wd = '/Users/songweizhi/Desktop/assess_linkages_hc'
# blastn_bin_vs_ref           = '%s/file_in/cami_hc_bin_vs_ref.tab'                                   % wd
# iden_cutoff                 = 99
# aln_len_cutoff              = 1000
# cov_q_cutoff                = 90
#
# for match in open(blastn_bin_vs_ref):
#     match_split = match.strip().split('\t')
#     query = match_split[0]
#     query_genome = '_'.join(query.split('_')[:3])
#     subject = match_split[1]
#     subject_genome = '_'.join(subject.split('_')[:2])
#     iden = float(match_split[2])
#     aln_len = int(match_split[3])
#     query_len = int(match_split[12])
#     subject_len = int(match_split[13])
#     coverage_q = float(aln_len) * 100 / float(query_len)
#     if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff) and (coverage_q >= cov_q_cutoff):
#         if query_genome == 'cami_hc_3':
#             print(match_split)



wd = '/Users/songweizhi/Desktop/assess_linkages_hc'
matam_16s_blastn            = '%s/file_in/matam_16S_uclust_0.999_vs_ref.tab'                        % wd
iden_cutoff_16s             = 99.8  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70

for match in open(matam_16s_blastn):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_gnm = subject.split('_16S_')[0]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    if (iden >= iden_cutoff_16s) and (aln_len > aln_len_cutoff_16s) and (coverage_q >= cov_q_cutoff_16s):
        if subject_gnm == 'tax_134676':
            print(match_split)



