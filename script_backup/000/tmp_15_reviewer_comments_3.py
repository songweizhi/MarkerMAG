from Bio import SeqIO

'''
unclustered
1	1/10	0.1
5	2/37	0.05405405405405406
10	5/37	0.13513513513513514
25	4/23	0.17391304347826086
50	3/22	0.13636363636363635
75	5/27	0.18518518518518517
100	6/23	0.2608695652173913

clustered (0.99)
1	1/6	    0.16666666666666666
5	2/22	0.09090909090909091
10	5/27	0.18518518518518517
25	3/14	0.21428571428571427
50	2/12	0.16666666666666666
75	4/14	0.2857142857142857
100	6/14	0.42857142857142855
Combined	20/60	0.3333333333333333

'''

wd              = '/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/submission_Bioinformatics/spurious_Matam_assemblies_sep'

subsample_rate  = '1'

seq_file        = '%s/MBARC26_SILVA138_id99_subsample_%s_0.99.fasta'        % (wd, subsample_rate)
blast_op        = '%s/MBARC26_SILVA138_id99_subsample_%s_0.99_vs_refs.txt'  % (wd, subsample_rate)

seq_file        = '%s/MBARC26_SILVA138_id99_subsample_%s.fasta'             % (wd, subsample_rate)
blast_op        = '%s/MBARC26_SILVA138_id99_subsample_%s_vs_refs.txt'       % (wd, subsample_rate)


subsample_rate  = 'All'
seq_file        = '%s/MBARC26_SILVA138_id99_assembled_16S_0.99.fasta'       % (wd)
blast_op        = '%s/MBARC26_SILVA138_id99_assembled_16S_0.99_vs_refs.txt' % (wd)



matam_16s_len_dict = {}
for each_16s in SeqIO.parse(seq_file, 'fasta'):
    matam_16s_len_dict[each_16s.id] = len(each_16s.seq)


seq_matched_to_refs = set()
seq_matched_to_refs_bad_match = set()
for each in open(blast_op):
    each_split = each.strip().split('\t')
    matam_16s_id = each_split[0]
    ref_id = each_split[1]
    matam_16s_len = matam_16s_len_dict[matam_16s_id]
    aln_len = int(each_split[3])
    query_cov = aln_len*100/matam_16s_len
    query_cov = float("{0:.2f}".format(query_cov))
    iden = float(each_split[2])

    if (iden >= 98) and (query_cov >= 85):
        seq_matched_to_refs.add(matam_16s_id)
    else:
        if matam_16s_id not in seq_matched_to_refs:
            print('%s\t%s\t%s\t%s/%s\t%s' % (matam_16s_id, ref_id, iden, aln_len, matam_16s_len, query_cov))
            seq_matched_to_refs_bad_match.add(matam_16s_id)


spurious_seq_list = []
for each in matam_16s_len_dict:
    if each not in seq_matched_to_refs:
        spurious_seq_list.append(each)


print('%s\t%s/%s\t%s' % (subsample_rate, len(spurious_seq_list), len(matam_16s_len_dict), (len(spurious_seq_list)/len(matam_16s_len_dict))))


