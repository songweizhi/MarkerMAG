import os


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def overlap_between_list(list_1, list_2):

    overlap = False
    for list_1_element in list_1:
        if list_1_element in list_2:
            overlap = True

    return overlap


def gnm_level_asessment(linked_mag_dict_rd1):

    rd1_correct_num = 0
    rd1_wrong_gnm = set()
    rd1_unknown_gnm = set()
    rd1_ambiguous_gnm = set()
    for each_rd_1_mag in linked_mag_dict_rd1:
        mag_assess = linked_mag_dict_rd1[each_rd_1_mag]
        if len(mag_assess) == 1:
            if mag_assess == {'Correct'}:
                rd1_correct_num += 1
            elif mag_assess == {'Wrong'}:
                rd1_wrong_gnm.add(each_rd_1_mag)
            elif mag_assess == {'Unknown'}:
                rd1_unknown_gnm.add(each_rd_1_mag)
        elif len(mag_assess) == 2:
            if ('Correct' in mag_assess) and ('Unknown' in mag_assess):
                rd1_correct_num += 1
            if ('Correct' in mag_assess) and ('Wrong' in mag_assess):
                rd1_ambiguous_gnm.add(each_rd_1_mag)
            if ('Wrong' in mag_assess) and ('Unknown' in mag_assess):
                rd1_wrong_gnm.add(each_rd_1_mag)
        elif len(mag_assess) == 3:
            rd1_ambiguous_gnm.add(each_rd_1_mag)

    return rd1_correct_num, rd1_unknown_gnm, rd1_wrong_gnm, rd1_ambiguous_gnm


########################################################################################################################

wd = '/Users/songweizhi/Desktop/assess_linkages_high'

########## reference to cluster ##########

drep_ani_cutoff             = 97
drep_cdb_file               = '%s/file_in/Cdb_%s.csv'                                       % (wd, drep_ani_cutoff)
#ref_to_strain_file          = '%s/file_in/ref_to_strain.txt'                                % wd

########## bin to reference ##########

parse_blastn_bin_vs_ref     = False  # True or False
blastn_bin_vs_ref           = '%s/file_in/cami_hc_bin_vs_ref.tab'                                   % wd
iden_cutoff                 = 99.5
aln_len_cutoff              = 1500
cov_q_cutoff                = 90
min_match_length            = 102400  # 100Kbp:102400, 500Kbp:524288,
#mag_metadata                = '%s/file_in/MAG_metadata.txt'                                 % wd

########## 16S to reference ##########

perform_blastn_16s_vs_refs  = False  # True or False
combined_GI_ref_16S         = '%s/file_in/combined_hc_ref_16S.ffn'                                  % wd
#matam_16s_seqs              = '%s/file_in/cami_hc_SILVA138_id99_assembled_16S_uclust_0.999.fasta'   % wd
matam_16s_seqs              = '%s/file_in/cami_hc_SILVA138_uclust_0.999_polished.fa'                % wd
matam_16s_blastn            = '%s/file_in/matam_16S_uclust_0.999_vs_ref.tab'                        % wd
iden_cutoff_16s             = 99.5  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 97

########## assessment results ##########

MarkerMAG_linkages          = '%s/hc_0514_stats_combined_filtered.txt'                  % wd
linkages_from_rd1           = True


'''
------------------------------------------------------------------------------------------------------------

# cami_hc_SILVA138_id99_50_subsample_100_3569
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_tmp_MarkerMAG_wd_up/hc_0512_tmp_step_1_wd
grep -E '@|cami_hc_SILVA138_id99_50_subsample_100_3569' /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0511_very_sensitive_MarkerMAG_wd/hc_0511_very_sensitive_step_1_wd/hc_0511_very_sensitive_input_reads_to_16S_reformatted.sam > input_reads_to_16S_subset.sam

# cami_hc_14___NODE_13073_length_6967_cov_0.639487
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_tmp_MarkerMAG_wd_up/hc_0512_tmp_step_1_wd
grep -E '@|cami_hc_14___C___NODE_13073_length_6967_cov_0.639487' rd1_extracted_to_gnm_reformatted.sam > rd1_extracted_to_gnm_subset.sam

------------------------------------------------------------------------------------------------------------


cami_hc_S1_24842618.2

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_MarkerMAG_wd/hc_0512_step_1_wd
grep 'cami_hc_S1_24842618.' hc_0512_input_reads_to_16S_reformatted.sam > reads_to_16s_subset.sam
grep 'cami_hc_S1_24842618.' rd1_extracted_to_gnm_reformatted.sam > extracted_to_gnm_subset.sam

bowtie2-build --quiet --threads 12 -f cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_cbd.fa cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_cbd
bowtie2 -x cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_cbd -1 cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R1.fa -2 cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R2.fa -S test_p.sam -p 12 -f --local --all --no-unal

bowtie2 -x cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_cbd -U cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R1.fa,cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R2.fa -S test_up.sam -p 12 -f --local --all --no-unal

bowtie2 -x cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_cbd -U cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R1.fa,cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R2.fa -S test_up2.sam -p 12 -f --local --all --no-unal -N 1
bowtie2 -x cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_cbd -U cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R1.fa,cami_hc_SILVA138_id99_50_subsample_100_3770___cami_hc_14___NODE_13073_length_6967_cov_0.639487_R2.fa -S test_up_L30.sam -p 12 -f --local --all --no-unal -N 1 -L 30


------------------------------------------------------------------------------------------------------------

cami_hc_SILVA138_id99_50_subsample_100_2419	cami_hc_66	126	S1	Correct
cami_hc_SILVA138_id99_50_subsample_100_1159	cami_hc_66	69
cami_hc_SILVA138_id99_50_subsample_100_1198	cami_hc_66	69

cami_hc_SILVA138_id99_50_subsample_100_1159	cami_hc_66	41	S1	Wrong
cami_hc_SILVA138_id99_50_subsample_100_1198	cami_hc_66	41	S1	Wrong
cami_hc_SILVA138_id99_50_subsample_100_2419	cami_hc_66	29                  why much less this time?


why lost the following three linkages in the new run?
cami_hc_SILVA138_id99_50_subsample_100_2419___cami_hc_66___NODE_1112_length_54304_cov_1.219168
cami_hc_S2_9871673.1
cami_hc_S2_9871673.2

cami_hc_SILVA138_id99_50_subsample_100_2419___cami_hc_66___NODE_1393_length_45743_cov_1.217991
cami_hc_SILVA138_id99_50_subsample_100_2419___cami_hc_66___NODE_5764_length_14276_cov_1.210482

# check cami_hc_S2_9871673
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_MarkerMAG_wd
grep 'cami_hc_S2_9871673.' ../hc_0511_very_sensitive_MarkerMAG_wd/hc_0511_very_sensitive_step_1_wd/hc_0511_very_sensitive_input_reads_to_16S_reformatted.sam > linked_to_16s.sam
grep 'cami_hc_S2_9871673.' hc_0512_step_1_wd/rd1_extracted_to_gnm_reformatted.sam > linked_to_gnm.sam

>cami_hc_S2_9871673.1
AATCAATAAAAACGTCAACGTTTTAGTTTTAAATCTTATCTGGCAGATGCCCGTCATCTGCAGATAAACAATTGAGCAAGTCAAACACTTTTTGGAGAGT
>cami_hc_S2_9871673.2
CTCACCCGTCCGCCGCTAACTTGAACGGAGCAAGCTCCGTCAAGTCCGCTCGACTTGCATGTATTAGGCACGCCGCCAGCGTTCGTCCTGAGCCAGGATC

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_MarkerMAG_wd/hc_0512_step_1_wd/cami_hc_refined_bins_db

# got it
bowtie2-build --quiet --threads 12 -f cami_hc_66___NODE_1112_length_54304_cov_1.219168.fa cami_hc_66___NODE_1112_length_54304_cov_1.219168
bowtie2 -x cami_hc_66___NODE_1112_length_54304_cov_1.219168 -1 cami_hc_S2_9871673.1.fa -2 cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673.sam -p 12 -f --local --all --no-unal
# cami_hc_S2_9871673.1	89	cami_hc_66___C___NODE_1112_length_54304_cov_1.219168	1	255	50S50M	=	1	0	ACTCTCCAAAAAGTGTTTGACTTGCTCAATTGTTTATCTGCAGATGACGGGCATCTGCCAGATAAGATTTAAAACTAAAACGTTGACGTTTTTATTGATT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	AS:i:100	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:50	YT:Z:UP

# not here !!!
bowtie2 -x cami_hc_refined_bins_combined -1 cami_hc_S2_9871673.1.fa -2 cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673_vs_all.sam -p 12 -f --local --all --no-unal

# not here !!!
bowtie2-build --quiet --threads 12 -f cami_hc_66.fasta cami_hc_66
bowtie2 -x cami_hc_66 -1 cami_hc_S2_9871673.1.fa -2 cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673_vs_hc_66_as_p.sam -p 12 -f --local --all --no-unal
bowtie2 -x cami_hc_66 -U cami_hc_S2_9871673.1.fa,cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673_vs_hc_66_as_up.sam -p 12 -f --local --all --no-unal
bowtie2 -x cami_hc_66 -1 cami_hc_S2_9871673.1.fa -2 cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673_vs_hc_66_as_p_all.sam -p 12 -f --local --all

bowtie2 -x cami_hc_refined_bins_renamed/cami_hc_66 -1 cami_hc_S2_9871673.1.fa -2 cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673_vs_hc_66_as_p_all.sam -p 12 -f --local --all


# got it
bowtie2-build --quiet --threads 12 -f Node_1112.fa Node_1112
bowtie2 -x Node_1112 -1 cami_hc_S2_9871673.1.fa -2 cami_hc_S2_9871673.2.fa -S cami_hc_S2_9871673_vs_1112.sam -p 12 -f --local --all --no-unal


# why not here
bowtie2 -x cami_hc_refined_bins_db/cami_hc_refined_bins_combined -1 rd1_extracted_R1.fasta -2 rd1_extracted_R2.fasta -S double_check.sam -p 12 -f --local --all --no-unal

# also not here
bowtie2 -x cami_hc_refined_bins_db/cami_hc_refined_bins_combined -1 rd1_extracted_R1.fasta -2 rd1_extracted_R2.fasta -S double_check_very_sensitive.sam -p 12 -f --local --all --no-unal --very-sensitive-local

bowtie2 -x cami_hc_refined_bins_db/cami_hc_refined_bins_combined -1 rd1_extracted_R1.fasta -2 rd1_extracted_R2.fasta -S double_check_super_sensitive.sam -p 12 -f --local --all --no-unal --very-sensitive-local



------------------------------------------------------------------------------------------------------------

cami_hc_S4_52228456.2       middle clip on      cami_hc_SILVA138_id99_75_subsample_75_3121      why?  A: not mapped, need to try very sensitive parameters
 
 
 
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_05010_MarkerMAG_wd/hc_05010_step_1_wd/cami_hc_SILVA138_uclust_0.999_polished_index




cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0511_super_sensitive_MarkerMAG_wd/hc_0511_super_sensitive_step_1_wd
grep cami_hc_S4_52228456.2 hc_0511_super_sensitive_input_reads_to_16S_best_match.sam > super_sensitive_subset.sam

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0511_very_sensitive_MarkerMAG_wd/hc_0511_very_sensitive_step_1_wd
grep cami_hc_S4_52228456.2 hc_0511_very_sensitive_input_reads_to_16S_best_match.sam > very_sensitive_subset.sam

grep cami_hc_S4_52228456.2 ../hc_05010_input_reads_to_16S.sam > vs_all_vs_all.sam


module load bowtie
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_05010_MarkerMAG_wd/hc_05010_step_1_wd/cami_hc_SILVA138_uclust_0.999_polished_index
bowtie2-build --quiet --threads 12 -f cami_hc_SILVA138_id99_75_subsample_75_3121.fa cami_hc_SILVA138_id99_75_subsample_75_3121
bowtie2-build --quiet --threads 12 -f cami_hc_SILVA138_uclust_0.999_polished.fa cami_hc_SILVA138_uclust_0.999_polished
bowtie2-build --quiet --threads 12 -f Five_16s.fa Five_16s

BioSAK select_seq -seq cami_hc_SILVA138_uclust_0.999_polished.fa -id Ten_16s.txt -out Ten_16s.fa -option 1 -oneline
BioSAK select_seq -seq cami_hc_SILVA138_uclust_0.999_polished.fa -id Twenty_16s.txt -out Twenty_16s.fa -option 1 -oneline
bowtie2-build --quiet --threads 12 -f Ten_16s.fa Ten_16s
bowtie2-build --quiet --threads 12 -f Twenty_16s.fa Twenty_16s

bowtie2 -x cami_hc_SILVA138_uclust_0.999_polished -1 cami_hc_S4_52228456.1.fas -2 cami_hc_S4_52228456.2.fas -S vs_all_default.sam -p 12 -f --local --all --no-unal
bowtie2 -x cami_hc_SILVA138_id99_75_subsample_75_3121 -1 cami_hc_S4_52228456.1.fas -2 cami_hc_S4_52228456.2.fas -S vs_single_default.sam -p 12 -f --local --all --no-unal
bowtie2 -x Five_16s -1 cami_hc_S4_52228456.1.fas -2 cami_hc_S4_52228456.2.fas -S vs_Five_16s_default.sam -p 12 -f --local --all --no-unal

bowtie2 -x Ten_16s -1 cami_hc_S4_52228456.1.fas -2 cami_hc_S4_52228456.2.fas -S vs_Ten_16s_default.sam -p 12 -f --local --all --no-unal
bowtie2 -x Twenty_16s -1 cami_hc_S4_52228456.1.fas -2 cami_hc_S4_52228456.2.fas -S vs_Twenty_16s_default.sam -p 12 -f --local --all --no-unal


bowtie2 -x cami_hc_SILVA138_id99_75_subsample_75_3121 -1 cami_hc_S4_52228456.1.fas -2 cami_hc_S4_52228456.2.fas -S cami_hc_S4_52228456.sam -p 12 -f --local --all --no-unal



cami_hc_SILVA138_id99_50_subsample_100_2419	cami_hc_66	126	S1	Correct
cami_hc_SILVA138_id99_50_subsample_100_1159	cami_hc_66	69
cami_hc_SILVA138_id99_50_subsample_100_1198	cami_hc_66	69

cami_hc_SILVA138_id99_50_subsample_100_1159	cami_hc_66	41	S1	Wrong
cami_hc_SILVA138_id99_50_subsample_100_1198	cami_hc_66	41	S1	Wrong
cami_hc_SILVA138_id99_50_subsample_100_2419	cami_hc_66	29                  why much less this time?



hc_0514_stats_combined_filtered.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	Accuracy
hc_0514_stats_combined_filtered.txt	Rd_1	|	93	87	0	6	87/93(93.55)	|	42	38	0	4	0	38/97(39.18)	38/42(90.48)
hc_0514_stats_combined_filtered.txt	Rd_2	|	0	0	0	0	0/0(0)	|	0	0	0	0	0	0/97(0.0)	0/0(0)
hc_0514_stats_combined_filtered.txt	Both	|	93	87	0	6	87/93(93.55)	|	42	38	0	4	0	38/97(39.18)	38/42(90.48)


hc_05010_stats_combined_filtered.txt	Rd_1	|	84	80	0	4	80/84(95.24)	|	40	36	0	4	0	36/97(37.11)	36/40(90.0)
hc_05010_stats_combined_filtered.txt	Rd_2	|	0	0	0	0	0/0(0)	        |	0	0	0	0	0	0/97(0.0)	    0/0(0)
hc_05010_stats_combined_filtered.txt	Both	|	84	80	0	4	80/84(95.24)	|	40	36	0	4	0	36/97(37.11)	36/40(90.0)

hc_0507_mis2_75_45_stats_combined_filtered.txt	Rd_1	|	39	33	1	5	33/38(86.84)	|	24	22	1	1	0	22/97(22.68)	22/23(95.65)
hc_0507_mis2_80_45_stats_combined_filtered.txt	Rd_1	|	45	36	1	8	36/44(81.82)	|	28	24	1	3	0	24/97(24.74)	24/27(88.89)
hc_0507_mis2_85_45_stats_combined_filtered.txt	Rd_1	|	45	36	1	8	36/44(81.82)	|	28	24	1	3	0	24/97(24.74)	24/27(88.89)

'''

'''
# cami_hc_S3_49165521.1

cami_hc_S5_42397362.1	{'cami_hc_SILVA138_id99_50_subsample_50_1792'}
cami_hc_S5_42397362.1	345	cami_hc_SILVA138_id99_75_subsample_75_2741	1	11	24S13=1X62=	=	1	0	CAGTTTTTCAACAGATTTTGTTGGAGAGTTTGATCCTCGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGAAAGGCCCTTCGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	AS:i:144	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:13G62	YT:Z:UP
cami_hc_S5_42397362.1	345	cami_hc_SILVA138_id99_50_subsample_10_876	1	11	24S13=1X62=	=	1	0	CAGTTTTTCAACAGATTTTGTTGGAGAGTTTGATCCTCGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGAAAGGCCCTTCGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	AS:i:144	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:13G62	YT:Z:UP
cami_hc_S5_42397362.1	89	cami_hc_SILVA138_id99_50_subsample_50_1792	1	11	16S21=1X62=	=	1	0	CAGTTTTTCAACAGATTTTGTTGGAGAGTTTGATCCTCGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAGCGGAAAGGCCCTTCGGG	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	AS:i:160	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:21G62	YT:Z:UP

correct
cami_hc_SILVA138_id99_75_subsample_75_2741	cami_hc_3
cami_hc_SILVA138_id99_50_subsample_10_876	cami_hc_3
cami_hc_SILVA138_id99_50_subsample_50_1792  cami_hc_13

wrong
MarkerGene__cami_hc_SILVA138_id99_50_subsample_50_1792,GenomicSeq__cami_hc_3,534
MarkerGene__cami_hc_SILVA138_id99_75_subsample_75_2741,GenomicSeq__cami_hc_3,294
MarkerGene__cami_hc_SILVA138_id99_50_subsample_10_876,GenomicSeq__cami_hc_3,294
MarkerGene__cami_hc_SILVA138_id99_50_subsample_50_1792,GenomicSeq__cami_hc_13,72
MarkerGene__cami_hc_SILVA138_id99_50_subsample_10_876,GenomicSeq__cami_hc_13,43
MarkerGene__cami_hc_SILVA138_id99_75_subsample_75_2741,GenomicSeq__cami_hc_13,38

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0505_mis2_75_45_keep_short_M_MarkerMAG_wd/hc_0505_mis2_75_45_keep_short_M_step_1_wd
grep cami_hc_SILVA138_id99_75_subsample_75_2741 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_1.sam
grep cami_hc_SILVA138_id99_50_subsample_10_876 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_2.sam
grep cami_hc_SILVA138_id99_50_subsample_50_1792 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_3.sam
cat sub_1.sam sub_2.sam sub_3.sam > sub.sam

grep cami_hc_SILVA138_id99_75_subsample_75_3138 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_4.sam
cat sub_1.sam sub_2.sam sub_3.sam sub_4.sam > sub.sam

grep cami_hc_SILVA138_id99_75_subsample_75_2741 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_1.sam
grep cami_hc_SILVA138_id99_50_subsample_10_876 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_2.sam
grep cami_hc_SILVA138_id99_50_subsample_50_1792 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_3.sam
cat sub_1.sam sub_2.sam sub_3.sam > sub.sam


grep cami_hc_SILVA138_id99_75_subsample_75_2741 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S.sam > sub_1all.sam
grep cami_hc_SILVA138_id99_50_subsample_10_876 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S.sam > sub_2all.sam
grep cami_hc_SILVA138_id99_50_subsample_50_1792 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S.sam > sub_3all.sam
cat sub_1all.sam sub_2all.sam sub_3all.sam > sub_all.sam


wrong
cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3(534)	NODE_780_length_68394_cov_3.028655	115	0	0	S1
cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3(534)	NODE_208_length_155603_cov_3.027781	136	0	0	S1
cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3(534)	NODE_701_length_72282_cov_3.046022	185	0	0	S1
cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3(534)	NODE_1353_length_46895_cov_3.028037	98	0	0	S1
MarkerGene__cami_hc_SILVA138_id99_50_subsample_50_1792,GenomicSeq__cami_hc_3,534
MarkerGene__cami_hc_SILVA138_id99_75_subsample_75_2741,GenomicSeq__cami_hc_3,294
MarkerGene__cami_hc_SILVA138_id99_50_subsample_10_876,GenomicSeq__cami_hc_3,294


hc_0505_mis2_75_45_drop_short_M	Rd_1	|	56	49	2	5	49/54(90.74)	|	40	35	2	3	0	35/97(36.08)	35/38(92.11)
hc_0505_mis2_75_45_drop_short_M	Rd_2	|	5	4	0	1	4/5(80.0)	    |	5	4	0	1	0	4/97(4.12)	    4/5(80.0)
hc_0505_mis2_75_45_drop_short_M	Both	|	61	53	2	6	53/59(89.83)	|	45	39	2	4	0	39/97(40.21)	39/43(90.7)

hc_0505_mis2_75_45_keep_short_M	Rd_1	|	53	46	2	5	46/51(90.2)	    |	40	35	2	3	0	35/97(36.08)	35/38(92.11)
hc_0505_mis2_75_45_keep_short_M	Rd_2	|	8	4	2	2	4/6(66.67)	    |	6	3	1	2	0	3/97(3.09)	    3/5(60.0)
hc_0505_mis2_75_45_keep_short_M	Both	|	61	50	4	7	50/57(87.72)	|	46	38	3	5	0	38/97(39.18)	38/43(88.37)

hc_0505_mis2_75_45	    Rd_1	|	60	41	2	17	41/58(70.69)	|	37	30	2	5	0	30/97(30.93)	30/35(85.71)
hc_0505_mis2_75_45	    Rd_2	|	12	9	2	1	9/10(90.0)	    |	8	6	1	1	0	6/97(6.19)	    6/7(85.71)
hc_0505_mis2_75_45	    Both	|	72	50	4	18	50/68(73.53)	|	45	36	3	6	0	36/97(37.11)	36/42(85.71)

hc_0505_mis1.9_75_45	Rd_1	|	72	42	2	28	42/70(60.0)	    |	37	27	2	8	0	27/97(27.84)	27/35(77.14)
hc_0505_mis1.9_75_45	Rd_2	|	13	7	2	4	7/11(63.64)	    |	9	6	1	2	0	6/97(6.19)	    6/8(75.0)
hc_0505_mis1.9_75_45	Both	|	85	49	4	32	49/81(60.49)	|	46	33	3	10	0	33/97(34.02)	33/43(76.74)

hc_0505_mis1.8_75_45	Rd_1	|	73	42	3	28	42/70(60.0)	    |	37	27	2	8	0	27/97(27.84)	27/35(77.14)
hc_0505_mis1.8_75_45	Rd_2	|	15	7	2	6	7/13(53.85)	    |	10	6	1	3	0	6/97(6.19)	    6/9(66.67)
hc_0505_mis1.8_75_45	Both	|	88	49	5	34	49/83(59.04)	|	47	33	3	11	0	33/97(34.02)	33/44(75.0)


stats_combined_filtered.txt	Rd_1	|	57	55	2	0	55/55(100.0)	|	35	33	2	0	0	33/97(34.02)	33/33(100.0)
stats_combined_filtered.txt	Rd_1	|	55	53	2	0	53/53(100.0)	|	35	33	2	0	0	33/97(34.02)	33/33(100.0)


	                                                    Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
CAMI_hc_0428_stats_combined_filtered.txt	            Rd_1	|	49	41	4	4	41/45(91.11)	|	41	33	4	4	0	33/97(34.02)	33/37(89.19)
CAMI_hc_0428_specific_stats_combined_filtered.txt	    Rd_1	|	48	41	3	4	41/45(91.11)	|	40	33	3	4	0	33/97(34.02)	33/37(89.19)
CAMI_hc_0428_very_specific_stats_combined_filtered.txt	Rd_1	|	42	34	3	5	34/39(87.18)	|	38	30	3	5	0	30/97(30.93)	30/35(85.71)


1	cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3	16S	C246_1	tax_1137264	NA
1	cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3	16S	C246_1	tax_1169165	NA
1	cami_hc_SILVA138_id99_50_subsample_50_1792___cami_hc_3	MAG	C248_0	tax_134676	NA

2	cami_hc_SILVA138_id99_75_subsample_75_1629___cami_hc_13	16S	C229_0	tax_1051629	NA
2	cami_hc_SILVA138_id99_75_subsample_75_1629___cami_hc_13	MAG	C246_1	tax_1137264	NA
2	cami_hc_SILVA138_id99_75_subsample_75_1629___cami_hc_13	MAG	C246_1	tax_1169165	NA

3	cami_hc_SILVA138_id99_75_subsample_75_5112___cami_hc_18	16S	C103_1	tax_1124469	NA
3	cami_hc_SILVA138_id99_75_subsample_75_5112___cami_hc_18	16S	C103_1	tax_991945	NA
3	cami_hc_SILVA138_id99_75_subsample_75_5112___cami_hc_18	16S	C103_1	tax_991968	NA
3	cami_hc_SILVA138_id99_75_subsample_75_5112___cami_hc_18	MAG	C100_0	tax_349966	NA

4	cami_hc_SILVA138_id99_50_subsample_50_2324___cami_hc_11	16S	C259_1	tax_1160164	NA
4	cami_hc_SILVA138_id99_50_subsample_50_2324___cami_hc_11	16S	C259_1	tax_1328339	NA
4	cami_hc_SILVA138_id99_50_subsample_50_2324___cami_hc_11	16S	C259_1	tax_520453	NA
4	cami_hc_SILVA138_id99_50_subsample_50_2324___cami_hc_11	MAG	C277_0	tax_1287182	NA

'''

########## script ##########

pwd_plot_sankey_R = '%s/file_in/get_sankey_plot.R' % wd

############################################### define file/folder name ################################################

# bin to reference
bin_vs_ref_txt              = '%s/bin_vs_ref_imag%s.txt'                        % (wd, iden_cutoff)
stats_bin_to_ref_txt        = '%s/stats_bin_to_ref_imag%s.txt'                  % (wd, iden_cutoff)
stats_ref_to_bin_txt        = '%s/stats_ref_to_bin_imag%s.txt'                  % (wd, iden_cutoff)
stats_bin_to_cluster_txt    = '%s/stats_bin_to_cluster_ani%s_imag%s.txt'        % (wd, drep_ani_cutoff, iden_cutoff)
stats_cluster_to_bin_txt    = '%s/stats_cluster_to_bin_ani%s_imag%s.txt'        % (wd, drep_ani_cutoff, iden_cutoff)

# assessment results
linkage_file_path, linkage_file_basename, linkage_file_extension = sep_path_basename_ext(MarkerMAG_linkages)
MarkerMAG_linkages_assessed = '%s/%s_with_assessment_ani%s_imag%s_i16S%s.txt'   % (linkage_file_path, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)
wrong_linkages_txt          = '%s/%s_wrong_ani%s_imag%s_i16S%s.txt'             % (linkage_file_path, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)
unknown_linkages_txt        = '%s/%s_unknown_ani%s_imag%s_i16S%s.txt'           % (linkage_file_path, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)


################################################# reference to cluster #################################################

cluster_to_ref_dict = {}
ref_to_cluster_dict = {}
for each_ref in open(drep_cdb_file):
    if not each_ref.startswith('genome,secondary_cluster'):
        each_ref_split = each_ref.strip().split(',')
        ref_file_name = each_ref_split[0]
        ref_file_name_no_ext = '.'.join(ref_file_name.split('.')[:-1])
        ref_cluster = 'C' + each_ref_split[1]
        ref_to_cluster_dict[ref_file_name_no_ext] = ref_cluster
        if ref_cluster not in cluster_to_ref_dict:
            cluster_to_ref_dict[ref_cluster] = [ref_file_name_no_ext]
        else:
            cluster_to_ref_dict[ref_cluster].append(ref_file_name_no_ext)


################################################### bin to reference ###################################################

bin_ref_connector           = '__|__'

# get ref_to_strain_dict
ref_to_strain_dict = {}
# for ref in open(ref_to_strain_file):
#     ref_split = ref.strip().split('\t')
#     ref_to_strain_dict[ref_split[0]] = ref_split[1]

if parse_blastn_bin_vs_ref is True:
    bin_vs_ref_dict = {}
    for match in open(blastn_bin_vs_ref):
        match_split = match.strip().split('\t')
        query = match_split[0]
        query_genome = '_'.join(query.split('_')[:3])
        subject = match_split[1]
        subject_genome = '_'.join(subject.split('_')[:2])
        iden = float(match_split[2])
        aln_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        coverage_q = float(aln_len) * 100 / float(query_len)
        if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff) and (coverage_q >= cov_q_cutoff):
            key_bin_ref = '%s%s%s' % (query_genome, bin_ref_connector, subject_genome)
            if key_bin_ref not in bin_vs_ref_dict:
                bin_vs_ref_dict[key_bin_ref] = aln_len
            else:
                bin_vs_ref_dict[key_bin_ref] += aln_len

    bin_vs_ref_dict_filtered = {}
    for each_key in bin_vs_ref_dict:
        if bin_vs_ref_dict[each_key] >= min_match_length:
            bin_vs_ref_dict_filtered[each_key] = bin_vs_ref_dict[each_key]

    # write out linkages
    linkage_txt_handle = open(bin_vs_ref_txt, 'w')
    linkage_txt_handle.write('Bin,Ref,Length\n')
    bin_to_ref_dict = {}
    ref_to_bin_dict = {}
    bin_to_cluster_dict = {}
    cluster_to_bin_dict = {}
    for each_linkage in bin_vs_ref_dict_filtered:
        linkage_split = each_linkage.split(bin_ref_connector)
        bin_id = linkage_split[0]
        ref_id = linkage_split[1]
        ref_cluster = ref_to_cluster_dict[ref_id]

        if bin_id not in bin_to_ref_dict:
            bin_to_ref_dict[bin_id] = {ref_id}
        else:
            bin_to_ref_dict[bin_id].add(ref_id)

        if bin_id not in bin_to_cluster_dict:
            bin_to_cluster_dict[bin_id] = {ref_cluster}
        else:
            bin_to_cluster_dict[bin_id].add(ref_cluster)

        if ref_id not in ref_to_bin_dict:
            ref_to_bin_dict[ref_id] = {bin_id}
        else:
            ref_to_bin_dict[ref_id].add(bin_id)

        if ref_cluster not in cluster_to_bin_dict:
            cluster_to_bin_dict[ref_cluster] = {bin_id}
        else:
            cluster_to_bin_dict[ref_cluster].add(bin_id)

        linkage_txt_handle.write('%s,%s,%s\n' % (bin_id, ref_id, bin_vs_ref_dict_filtered[each_linkage]))
    linkage_txt_handle.close()

    # visualize
    cmd_sankey_bin_to_ref = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, bin_vs_ref_txt, 600, 1800)
    os.system(cmd_sankey_bin_to_ref)

    bin_to_ref_txt_handle = open(stats_bin_to_ref_txt, 'w')
    for each_bin in bin_to_ref_dict:
        bin_to_ref_txt_handle.write('%s\t%s\n' % (each_bin, ','.join(bin_to_ref_dict[each_bin])))
    bin_to_ref_txt_handle.close()

    ref_to_bin_txt_handle = open(stats_ref_to_bin_txt, 'w')
    for each_ref in ref_to_bin_dict:
        ref_to_bin_txt_handle.write('%s\t%s\t%s\n' % (each_ref, ','.join(ref_to_bin_dict[each_ref]), ref_to_strain_dict.get(each_ref, 'NA')))
    ref_to_bin_txt_handle.close()

    stats_bin_to_cluster_txt_handle = open(stats_bin_to_cluster_txt, 'w')
    for each_bin in bin_to_cluster_dict:
        stats_bin_to_cluster_txt_handle.write('%s\t%s\n' % (each_bin, ','.join(bin_to_cluster_dict[each_bin])))
    stats_bin_to_cluster_txt_handle.close()

    stats_cluster_to_bin_txt_handle = open(stats_cluster_to_bin_txt, 'w')
    for each_ref in cluster_to_bin_dict:
        stats_cluster_to_bin_txt_handle.write('%s\t%s\n' % (each_ref, ','.join(cluster_to_bin_dict[each_ref])))
    stats_cluster_to_bin_txt_handle.close()

bin_to_cluster_dict = {}
for each_match in open(stats_bin_to_cluster_txt):
    each_match_split = each_match.strip().split('\t')
    bin_to_cluster_dict[each_match_split[0]] = {i for i in each_match_split[1].split(',')}


cluster_to_bin_dict = {}
for each_bin in bin_to_cluster_dict:
    matched_clusters =  bin_to_cluster_dict[each_bin]
    for matched_cluster in matched_clusters:
        if matched_cluster not in cluster_to_bin_dict:
            cluster_to_bin_dict[matched_cluster] = {each_bin}
        else:
            cluster_to_bin_dict[matched_cluster].add(each_bin)


################################################### 16S to reference ###################################################

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 1'
blast_cmd = 'blastn -query %s -subject %s -out %s %s' % (matam_16s_seqs, combined_GI_ref_16S, matam_16s_blastn, blast_parameters)
if perform_blastn_16s_vs_refs is True:
    os.system(blast_cmd)

# get matam_16s_to_cluster_dict
matam_16s_to_cluster_dict = {}
for match in open(matam_16s_blastn):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_gnm = subject.split('_16S_')[0]
    subject_cluster = ref_to_cluster_dict[subject_gnm]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    if (iden >= iden_cutoff_16s) and (aln_len > aln_len_cutoff_16s) and (coverage_q >= cov_q_cutoff_16s):
        if query not in matam_16s_to_cluster_dict:
            matam_16s_to_cluster_dict[query] = {subject_cluster}
        else:
            matam_16s_to_cluster_dict[query].add(subject_cluster)

cluster_to_matam_16s_dict = {}
for matam_16s in matam_16s_to_cluster_dict:
    matched_clusters =  matam_16s_to_cluster_dict[matam_16s]
    for matched_cluster in matched_clusters:
        if matched_cluster not in cluster_to_matam_16s_dict:
            cluster_to_matam_16s_dict[matched_cluster] = {matam_16s}
        else:
            cluster_to_matam_16s_dict[matched_cluster].add(matam_16s)


###################################################### assessment ######################################################

MarkerMAG_linkages_assessed_handle = open(MarkerMAG_linkages_assessed, 'w')
wrong_linkages_txt_handle = open(wrong_linkages_txt, 'w')
unknown_linkages_txt_handle = open(unknown_linkages_txt, 'w')
linkage_num_right = 0
linkage_num_wrong = 0
linkage_num_unknown = 0
linkage_assessment_dict = {}
for each_linkage in open(MarkerMAG_linkages):
    #if each_linkage.startswith('MarkerGene\tGenomicSeq\tLinkage\tStep'):
    if ('MarkerGene\tGenomicSeq\tLinkage\tStep' in each_linkage) or ('MarkerGene,GenomicSeq,Number' in each_linkage):
        MarkerMAG_linkages_assessed_handle.write('MarkerGene\tGenomicSeq\tLinkage\tStep\tAssessment\n')
    else:
        id_16s = ''
        id_mag = ''
        link_num = ''
        if linkages_from_rd1 is False:
            each_linkage_split = each_linkage.strip().split('\t')
            id_16s = each_linkage_split[0]
            id_mag = each_linkage_split[1]
            link_num = each_linkage_split[2]
        else:
            each_linkage_split = each_linkage.strip().split(',')
            id_16s = each_linkage_split[0][12:]
            id_mag = each_linkage_split[1][12:]
            link_num = each_linkage_split[2]

        #print('id_mag: %s' % id_mag)
        key_16s_mag = '%s___%s' % (id_16s, id_mag)

        matched_cluster_16s = {}
        if id_16s in matam_16s_to_cluster_dict:
            matched_cluster_16s = matam_16s_to_cluster_dict[id_16s]

        matched_cluster_mag = {}
        if id_mag in bin_to_cluster_dict:
            matched_cluster_mag = bin_to_cluster_dict[id_mag]

        if (matched_cluster_16s != {}) and (matched_cluster_mag != {}):

            if overlap_between_list(matched_cluster_mag, matched_cluster_16s) is True:
                linkage_num_right += 1

                if linkages_from_rd1 is False:
                    MarkerMAG_linkages_assessed_handle.write('%s\tCorrect\n' % each_linkage.strip())
                else:
                    MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tS1\tCorrect\n' % (id_16s, id_mag, link_num))

                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Correct']
                else:
                    linkage_assessment_dict[id_mag].append('Correct')

            else:
                print('matched_cluster_16s: %s' % matched_cluster_16s)
                print('matched_cluster_mag: %s' % matched_cluster_mag)
                print()
                linkage_num_wrong += 1

                if linkages_from_rd1 is False:
                    MarkerMAG_linkages_assessed_handle.write('%s\tWrong\n' % each_linkage.strip())
                else:
                    MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tS1\tWrong\n' % (id_16s, id_mag, link_num))

                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Wrong']
                else:
                    linkage_assessment_dict[id_mag].append('Wrong')

                for each_16s_cluster in matched_cluster_16s:
                    current_16s_cluster_ref = cluster_to_ref_dict[each_16s_cluster]
                    for each_ref_1 in current_16s_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict.get(each_ref_1, 'NA')))

                for each_mag_cluster in matched_cluster_mag:
                    current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                    for each_ref_2 in current_mag_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict.get(each_ref_2, 'NA')))
                wrong_linkages_txt_handle.write('\n')
        else:
            linkage_num_unknown += 1

            if linkages_from_rd1 is False:
                MarkerMAG_linkages_assessed_handle.write('%s\tUnknown\n' % each_linkage.strip())
            else:
                MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tS1\tUnknown\n' % (id_16s, id_mag, link_num))

            if id_mag not in linkage_assessment_dict:
                linkage_assessment_dict[id_mag] = ['Unknown']
            else:
                linkage_assessment_dict[id_mag].append('Unknown')

            for each_16s_cluster in matched_cluster_16s:
                current_16s_cluster_ref = cluster_to_ref_dict[each_16s_cluster]
                for each_ref_1 in current_16s_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict.get(each_ref_1, 'NA')))
            for each_mag_cluster in matched_cluster_mag:
                current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                for each_ref_2 in current_mag_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict.get(each_ref_2, 'NA')))
            unknown_linkages_txt_handle.write('\n')
wrong_linkages_txt_handle.close()
unknown_linkages_txt_handle.close()
MarkerMAG_linkages_assessed_handle.close()


# get assess stats at mag level
linked_mag_dict_rd1 = {}
linked_mag_dict_rd2 = {}
total_linkage_num_rd1 = 0
total_linkage_num_rd2 = 0
correct_link_rd1 = 0
correct_link_rd2 = 0
wrong_link_rd1 = 0
wrong_link_rd2 = 0
unknown_link_rd1 = 0
unknown_link_rd2 = 0
for each_link in open(MarkerMAG_linkages_assessed):
    if not each_link.startswith('MarkerGene	GenomicSeq	Linkage	Step	Assessment'):
        each_link_split = each_link.strip().split('\t')
        gnm_id = each_link_split[1]
        linked_rd = each_link_split[3]
        assessment = each_link_split[4]
        if linked_rd == 'S1':
            total_linkage_num_rd1 += 1

            if assessment == 'Correct':
                correct_link_rd1 += 1
            if assessment == 'Wrong':
                wrong_link_rd1 += 1
            if assessment == 'Unknown':
                unknown_link_rd1 += 1

            if gnm_id not in linked_mag_dict_rd1:
                linked_mag_dict_rd1[gnm_id] = {assessment}
            else:
                linked_mag_dict_rd1[gnm_id].add(assessment)
        if linked_rd == 'S2':
            total_linkage_num_rd2 += 1

            if assessment == 'Correct':
                correct_link_rd2 += 1
            if assessment == 'Wrong':
                wrong_link_rd2 += 1
            if assessment == 'Unknown':
                unknown_link_rd2 += 1

            if gnm_id not in linked_mag_dict_rd1:
                linked_mag_dict_rd2[gnm_id] = {assessment}
            else:
                linked_mag_dict_rd2[gnm_id].add(assessment)


total_linkage_num_both = total_linkage_num_rd1 + total_linkage_num_rd2
total_linked_mag_num = len(linked_mag_dict_rd1) + len(linked_mag_dict_rd2)

rd1_correct_num, rd1_unknown_gnm, rd1_wrong_gnm, rd1_ambiguous_gnm = gnm_level_asessment(linked_mag_dict_rd1)
rd2_correct_num, rd2_unknown_gnm, rd2_wrong_gnm, rd2_ambiguous_gnm = gnm_level_asessment(linked_mag_dict_rd2)

correct_num_both        = rd1_correct_num + rd2_correct_num
unknown_gnm_both_num    = len(rd1_unknown_gnm) + len(rd2_unknown_gnm)
wrong_gnm_both_num      = len(rd1_wrong_gnm) + len(rd2_wrong_gnm)
ambiguous_gnm_both_num  = len(rd1_ambiguous_gnm) + len(rd2_ambiguous_gnm)

recovery_str_rd1   = '%s/%s(%s)' % (rd1_correct_num,  total_query_mag_num, float("{0:.2f}".format(rd1_correct_num*100/total_query_mag_num)))
recovery_str_rd2   = '%s/%s(%s)' % (rd2_correct_num,  total_query_mag_num, float("{0:.2f}".format(rd2_correct_num*100/total_query_mag_num)))
recovery_str_both  = '%s/%s(%s)' % (correct_num_both, total_query_mag_num, float("{0:.2f}".format(correct_num_both*100/total_query_mag_num)))
accuracy_str_rd1   = '%s/%s(%s)' % (rd1_correct_num,  (len(linked_mag_dict_rd1) - len(rd1_unknown_gnm)), float("{0:.2f}".format(rd1_correct_num*100/(len(linked_mag_dict_rd1) - len(rd1_unknown_gnm)))))
accuracy_str_rd2 = '0/0(0)'
if (len(linked_mag_dict_rd2) - len(rd2_unknown_gnm)) > 0:
    accuracy_str_rd2   = '%s/%s(%s)' % (rd2_correct_num,  (len(linked_mag_dict_rd2) - len(rd2_unknown_gnm)), float("{0:.2f}".format(rd2_correct_num*100/(len(linked_mag_dict_rd2) - len(rd2_unknown_gnm)))))
accuracy_str_both  = '%s/%s(%s)' % (correct_num_both, (total_linked_mag_num - unknown_gnm_both_num), float("{0:.2f}".format(correct_num_both*100/(total_linked_mag_num - unknown_gnm_both_num))))

correct_link_both = correct_link_rd1 + correct_link_rd2
wrong_link_both   = wrong_link_rd1 + wrong_link_rd2
unknown_link_both = unknown_link_rd1 + unknown_link_rd2

accuracy_str_rd1_link_level   = '%s/%s(%s)' % (correct_link_rd1,  (total_linkage_num_rd1 - unknown_link_rd1), float("{0:.2f}".format(correct_link_rd1*100/(total_linkage_num_rd1 - unknown_link_rd1))))
accuracy_str_rd2_link_level = '0/0(0)'
if (total_linkage_num_rd2 - unknown_link_rd2) > 0:
    accuracy_str_rd2_link_level   = '%s/%s(%s)' % (correct_link_rd2,  (total_linkage_num_rd2 - unknown_link_rd2), float("{0:.2f}".format(correct_link_rd2*100/(total_linkage_num_rd2 - unknown_link_rd2))))
accuracy_str_both_link_level  = '%s/%s(%s)' % (correct_link_both, (total_linkage_num_both - unknown_link_both), float("{0:.2f}".format(correct_link_both*100/(total_linkage_num_both - unknown_link_both))))

prefix = os.path.basename(MarkerMAG_linkages).split('_identified_linkages_genome_level')[0]
print('%s\tRound\t|\tLink\tYes\tNA\tNo\tAccuracy\t|\tMAG\tYes\tNA\tNo\tY/N\tRecovery\tAccuracy' % prefix)
print('%s\tRd_1\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd1,  correct_link_rd1,  unknown_link_rd1,  wrong_link_rd1,  accuracy_str_rd1_link_level,  len(linked_mag_dict_rd1), rd1_correct_num,  len(rd1_unknown_gnm), len(rd1_wrong_gnm), len(rd1_ambiguous_gnm), recovery_str_rd1, accuracy_str_rd1))
print('%s\tRd_2\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd2,  correct_link_rd2,  unknown_link_rd2,  wrong_link_rd2,  accuracy_str_rd2_link_level,  len(linked_mag_dict_rd2), rd2_correct_num,  len(rd2_unknown_gnm), len(rd2_wrong_gnm), len(rd2_ambiguous_gnm), recovery_str_rd2, accuracy_str_rd2))
print('%s\tBoth\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_both, correct_link_both, unknown_link_both, wrong_link_both, accuracy_str_both_link_level, total_linked_mag_num,      correct_num_both, unknown_gnm_both_num, wrong_gnm_both_num, ambiguous_gnm_both_num, recovery_str_both, accuracy_str_both))
