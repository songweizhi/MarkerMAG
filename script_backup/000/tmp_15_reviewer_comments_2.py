from Bio import SeqIO


seq_id = 'MBARC26_SILVA138_id99_subsample_75_165'


matam_16s_len_dict = {}
for each_16s in SeqIO.parse('/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/submission_Bioinformatics/spurious_Matam_assemblies_sep/MBARC26_SILVA138_id99_assembled_16S_unclustered.fasta', 'fasta'):
    matam_16s_len_dict[each_16s.id] = len(each_16s.seq)


seq_matched_to_refs = set()
seq_matched_to_refs_bad_match = set()
for each in open('/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/submission_Bioinformatics/MarkerMAG_Matam/unclustered_vs_refs_BestHit.tab'):
    each_split = each.strip().split('\t')
    matam_16s_id = each_split[0]
    ref_id = each_split[1]
    matam_16s_len = matam_16s_len_dict[matam_16s_id]
    aln_len = int(each_split[3])
    query_cov = aln_len*100/matam_16s_len
    query_cov = float("{0:.2f}".format(query_cov))
    iden = float(each_split[2])

    if (iden >= 98) and (query_cov >= 95):
        seq_matched_to_refs.add(matam_16s_id)
    else:
        if matam_16s_id not in seq_matched_to_refs:
            print('%s\t%s\t%s\t%s/%s\t%s' % (matam_16s_id, ref_id, iden, aln_len, matam_16s_len, query_cov))
            seq_matched_to_refs_bad_match.add(matam_16s_id)



for each in matam_16s_len_dict:
    if (each not in seq_matched_to_refs) and (each not in seq_matched_to_refs_bad_match):
        print(each)




# seq_id = 'MBARC26_SILVA138_id99_subsample_75_165'
# for each in open('/Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/submission_Bioinformatics/MarkerMAG_Matam/unclustered_vs_refs.tab'):
#     each_split = each.strip().split('\t')
#     matam_16s_id = each_split[0]
#     matam_16s_len = matam_16s_len_dict[matam_16s_id]
#     aln_len = int(each_split[3])
#     if matam_16s_id == seq_id:
#         print(each_split)
# print('%s\t%s bp' % (seq_id, matam_16s_len_dict[seq_id]))



'''

cd /Users/songweizhi/Documents/Research/MarkerMAG/0_Manuscript/submission_Bioinformatics/MarkerMAG_Matam
/Users/songweizhi/Software/usearch/usearch -makeudb_usearch ref_16S.fasta -output ref_16S.fasta.udb 

/Users/songweizhi/Software/usearch/usearch -makeudb_usearch ref_16S_uclust_0.99.fasta -output ref_16S_uclust_0.99.fasta.udb 

/Users/songweizhi/Software/usearch/usearch -uchime2_ref MBARC26_SILVA138_id99_assembled_16S_unclustered.fasta -db ref_16S.fasta.udb -uchimeout out_uchime.txt -strand plus -mode sensitive
00:02 52Mb    100.0% Chimeras 10/179 (5.6%), in db 117 (65.4%), not matched 52 (29.1%)

/Users/songweizhi/Software/usearch/usearch -uchime2_ref MBARC26_SILVA138_id99_assembled_16S_unclustered.fasta -db ref_16S_uclust_0.99.fasta.udb -uchimeout out_uchime_uclust.txt -strand plus -mode sensitive


1   1   10
5   2   37
10  5   37
25  4   23
50  3   22
75  5   27
100 6   23
179

MBARC26_SILVA138_id99_subsample_1_13

MBARC26_SILVA138_id99_subsample_5_106	PS_16S_1	100.0	252/518	48.65
MBARC26_SILVA138_id99_subsample_5_35

MBARC26_SILVA138_id99_subsample_10_105	CP_16S_8	100.0	681/900	75.67
MBARC26_SILVA138_id99_subsample_10_114	CP_16S_8	100.0	681/831	81.95
MBARC26_SILVA138_id99_subsample_10_248
MBARC26_SILVA138_id99_subsample_10_45
MBARC26_SILVA138_id99_subsample_10_94	PS_16S_1	100.0	271/578	46.89

MBARC26_SILVA138_id99_subsample_25_105	CP_16S_2	99.844	642/960	66.88
MBARC26_SILVA138_id99_subsample_25_125	CP_16S_2	99.844	642/863	74.39
MBARC26_SILVA138_id99_subsample_25_223
MBARC26_SILVA138_id99_subsample_25_40

MBARC26_SILVA138_id99_subsample_50_111	CP_16S_2	99.814	539/771	69.91
MBARC26_SILVA138_id99_subsample_50_27
MBARC26_SILVA138_id99_subsample_50_96	CP_16S_2	99.814	539/860	62.67

MBARC26_SILVA138_id99_subsample_75_105	CP_16S_5	99.738	381/613	62.15
MBARC26_SILVA138_id99_subsample_75_25
MBARC26_SILVA138_id99_subsample_75_261	HB_16S_2	86.496	822/821	100.12
MBARC26_SILVA138_id99_subsample_75_92	CP_16S_5	99.738	381/702	54.27
MBARC26_SILVA138_id99_subsample_75_99	CP_16S_5	99.738	381/703	54.2

MBARC26_SILVA138_id99_subsample_100_102	CP_16S_3	100.0	338/662	51.06
MBARC26_SILVA138_id99_subsample_100_20
MBARC26_SILVA138_id99_subsample_100_21
MBARC26_SILVA138_id99_subsample_100_25
MBARC26_SILVA138_id99_subsample_100_266	HB_16S_2	86.596	1134/1136	99.82
MBARC26_SILVA138_id99_subsample_100_94	CP_16S_5	99.704	338/661	51.13

'''