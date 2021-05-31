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

wd = '/Users/songweizhi/Desktop/assess_linkages'

########## reference to cluster ##########

drep_ani_cutoff             = 97
drep_cdb_file               = '%s/file_in/Cdb_%s.csv'                                       % (wd, drep_ani_cutoff)
ref_to_strain_file          = '%s/file_in/ref_to_strain.txt'                                % wd

########## bin to reference ##########

parse_blastn_bin_vs_ref     = False  # True or False
blastn_bin_vs_ref           = '%s/file_in/bin_vs_ref.tab'                                   % wd
iden_cutoff                 = 99.5
aln_len_cutoff              = 1500
cov_q_cutoff                = 90
min_match_length            = 524288  # 100 Kbp  102400
mag_metadata                = '%s/file_in/MAG_metadata.txt'                                 % wd

'''
GI_0414_specific	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0414_specific	Rd_2	|	22	22	0	0	22/22(100.0)	|	16	16	0	0	0	16/97(16.49)	16/16(100.0)
GI_0414_specific	Both	|	72	71	1	0	71/71(100.0)	|	48	48	0	0	0	48/97(49.48)	48/48(100.0)


'''
########## 16S to reference ##########

perform_blastn_16s_vs_refs  = False  # True or False
combined_GI_ref_16S         = '%s/file_in/combined_GI_ref_16S.ffn'           % wd
matam_16s_seqs              = '%s/file_in/GI_128_16S_0.999.fasta'            % wd
matam_16s_blastn            = '%s/file_in/GI_128_16S_0.999_vs_ref.tab'       % wd
iden_cutoff_16s             = 99.5  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 97

########## assessment results ##########

MarkerMAG_linkages          = '%s/GI_0528_128_60_60_stats_combined_filtered.txt'             % wd
mlen = '45'
#MarkerMAG_linkages          = '/Users/songweizhi/Desktop/test_filter/GI_0524_128_%s_%s_stats_combined_filtered.txt' % (mlen, mlen)

MarkerMAG_linkages          = '/Users/songweizhi/Desktop/tunning_rd1/stats_combined_filtered.txt'
linkages_from_rd1           = True


'''

3_GI_subsample_25_484	Refined_95	43	S1	Correct
3_GI_subsample_50_869	Refined_95	42	S1	Correct
3_GI_subsample_50_811	Refined_95	42	S1	Correct
3_GI_subsample_75_1381	Refined_95	42	S1	Correct


3_GI_subsample_100_1336	Refined_95	15	S1	Wrong
3_GI_subsample_50_650	Refined_95	15	S1	Wrong
3_GI_subsample_75_1159	Refined_95	15	S1	Wrong
3_GI_subsample_75_1644	Refined_95	15	S1	Wrong





# new
GI_0524_128_45_45_stats_combined_filtered.txt	Rd_1	|	128	101	14	13	101/114(88.6)	|	40	30	3	3	4	30/97(30.93)	30/37(81.08)
GI_0524_128_55_55_stats_combined_filtered.txt	Rd_1	|	135	108	14	13	108/121(89.26)	|	35	31	1	2	1	31/97(31.96)	31/34(91.18)
GI_0524_128_60_60_stats_combined_filtered.txt	Rd_1	|	102	95	7	0	95/95(100.0)	|	30	29	1	0	0	29/97(29.9)	    29/29(100.0)
GI_0524_128_65_65_stats_combined_filtered.txt	Rd_1	|	69	65	4	0	65/65(100.0)	|	28	27	1	0	0	27/97(27.84)	27/27(100.0)


# old
GI_0524_128_45_45_stats_combined_filtered.txt	Rd_1	|	124	96	14	14	96/110(87.27)	|	36	27	3	4	2	27/97(27.84)	27/33(81.82)
GI_0524_128_55_55_stats_combined_filtered.txt	Rd_1	|	140	96	20	24	96/120(80.0)	|	31	28	1	1	1	28/97(28.87)	28/30(93.33)
GI_0524_128_60_60_stats_combined_filtered.txt	Rd_1	|	106	96	8	2	96/98(97.96)	|	29	27	1	0	1	27/97(27.84)	27/28(96.43)
GI_0524_128_65_65_stats_combined_filtered.txt	Rd_1	|	68	62	6	0	62/62(100.0)	|	27	26	1	0	0	26/97(26.8)	    26/26(100.0)


3_GI_subsample_75_1	Refined_28	10	S1	Correct

3_GI_subsample_75_3263	Refined_28	42	S1	Wrong
3_GI_subsample_100_4409	Refined_28	42	S1	Unknown
3_GI_subsample_25_928	Refined_28	35	S1	Unknown

3_GI_subsample_100_3600	Refined_36	203	S1	Unknown
MarkerGene__3_GI_subsample_100_1395,GenomicSeq__Refined_36,221









GI_0524_128_60_60_stats_combined_filtered.txt	            Rd_1	|	106	96	8	2	96/98(97.96)	|	29	27	1	0	1	27/97(27.84)	27/28(96.43)
GI_0524_128_55_55_stats_combined_filtered.txt	            Rd_1	|	140	96	20	24	96/120(80.0)	|	31	28	1	1	1	28/97(28.87)	28/30(93.33)
GI_0524_128_45_45_no_clp_check2_stats_combined_filtered.txt	Rd_1	|	107	91	12	4	91/95(95.79)	|	35	29	3	2	1	29/97(29.9)	    29/32(90.62)
GI_0524_128_45_45_no_clp_check2_stats_combined_filtered.txt	Rd_1	|	135	110	9	16	110/126(87.3)	|	33	26	2	4	1	26/97(26.8)	    26/31(83.87)









3_GI_subsample_5_85	Refined_59	156	S1	Correct

3_GI_subsample_75_1344	Refined_59	57	S1	Wrong
3_GI_subsample_75_1054	Refined_59	57	S1	Wrong
3_GI_subsample_50_708	Refined_59	57	S1	Wrong
3_GI_subsample_75_1194	Refined_59	50	S1	Wrong
3_GI_subsample_50_649	Refined_59	50	S1	Wrong
3_GI_subsample_100_1406	Refined_59	50	S1	Wrong








GI_0508_mis2_45_45	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0508_mis2_45_45	Rd_1	|	98	91	6	1	91/92(98.91)	|	31	30	0	0	1	30/97(30.93)	30/31(96.77)
GI_0508_mis2_45_45	Rd_2	|	18	17	1	0	17/17(100.0)	|	9	9	0	0	0	9/97(9.28)	    9/9(100.0)
GI_0508_mis2_45_45	Both	|	116	108	7	1	108/109(99.08)	|	40	39	0	0	1	39/97(40.21)	39/40(97.5)



3_GI_subsample_5_85 Refined_59  161
3_GI_subsample_5_85 Refined_52  160

3_GI_subsample_5_85___Refined_52___NODE_19882_length_2622_cov_0.501804
S11_9543565.2   2/39  should be ignored!

can S11_9543565.2 mapped to 3_GI_subsample_5_85 with very sensitive mode?


cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0508_mis2_45_45_MarkerMAG_wd/GI_0508_mis2_45_45_step_1_wd
grep -E '@|S11_9543565.' GI_0508_mis2_45_45_input_reads_to_16S_best_match.sam > input_reads_to_16S_best_match.sam
grep -E '@|S11_9543565.' GI_0508_mis2_45_45_unmapped_mates_to_16s_reformat.sam > unmapped_mates_to_16S.sam

grep -E '@|S4_602526.' /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd/GI_0507_mis2_45_45_input_reads_to_16S_best_match.sam > sub.sam
grep -E '@|S4_602526.' GI_0507_mis2_45_45_clipping_parts_best_match.sam > sub.sam

GI_0508_mis2_45_45	Rd_1	|	98	91	6	1	91/92(98.91)	|	31	30	0	0	1	30/97(30.93)	30/31(96.77)
GI_0508_mis2_45_45	Rd_2	|	18	17	1	0	17/17(100.0)	|	9	8	1	0	0	8/97(8.25)	    8/8(100.0)
GI_0508_mis2_45_45	Both	|	116	108	7	1	108/109(99.08)	|	40	38	1	0	1	38/97(39.18)	38/39(97.44)

GI_0508_mis2_45_45	Rd_1	|	98	91	6	1	91/92(98.91)	|	31	30	0	0	1	30/97(30.93)	30/31(96.77)
GI_0508_mis2_45_45	Rd_2	|	18	17	1	0	17/17(100.0)	|	9	9	0	0	0	9/97(9.28)	    9/9(100.0)
GI_0508_mis2_45_45	Both	|	116	108	7	1	108/109(99.08)	|	40	39	0	0	1	39/97(40.21)	39/40(97.5)

GI_0507_mis2_45_45	Rd_1	|	103	90	6	7	90/97(92.78)	|	32	30	0	1	1	30/97(30.93)	30/32(93.75)
GI_0507_mis2_45_45	Rd_2	|	16	16	0	0	16/16(100.0)	|	9	9	0	0	0	9/97(9.28)	    9/9(100.0)
GI_0507_mis2_45_45	Both	|	119	106	6	7	106/113(93.81)	|	41	39	0	1	1	39/97(40.21)	39/41(95.12)

GI_0420_good	    Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0420_good	    Rd_2	|	21	19	0	2	19/21(90.48)	|	17	15	0	2	0	15/97(15.46)	15/17(88.24)
GI_0420_good	    Both	|	74	71	1	2	71/73(97.26)	|	51	49	0	2	0	49/97(50.52)	49/51(96.08)


GI_0505_mis2_75_45_drop_short_M	Rd_1	|	83	74	4	5	74/79(93.67)	|	28	26	0	1	1	26/97(26.8)	26/28(92.86)
GI_0505_mis2_75_45_drop_short_M	Rd_2	|	4	4	0	0	4/4(100.0)	    |	3	3	0	0	0	3/97(3.09)	3/3(100.0)
GI_0505_mis2_75_45_drop_short_M	Both	|	87	78	4	5	78/83(93.98)	|	31	29	0	1	1	29/97(29.9)	29/31(93.55)

GI_0505_mis2_75_45_keep_short_M	Rd_1	|	75	70	4	1	70/71(98.59)	|	28	27	0	0	1	27/97(27.84)	27/28(96.43)
GI_0505_mis2_75_45_keep_short_M	Rd_2	|	36	13	1	22	13/35(37.14)	|	9	7	0	2	0	7/97(7.22)	    7/9(77.78)
GI_0505_mis2_75_45_keep_short_M	Both	|	111	83	5	23	83/106(78.3)	|	37	34	0	2	1	34/97(35.05)	34/37(91.89)

GI_0505_mis2_75_45	Rd_1	|	83	78	4	1	78/79(98.73)	|	28	27	0	0	1	27/97(27.84)	27/28(96.43)
GI_0505_mis2_75_45	Rd_2	|	11	10	1	0	10/10(100.0)	|	7	6	1	0	0	6/97(6.19)	    6/6(100.0)
GI_0505_mis2_75_45	Both	|	94	88	5	1	88/89(98.88)	|	35	33	1	0	1	33/97(34.02)	33/34(97.06)

GI_0505_mis2_75_45	Rd_1	|	83	78	4	1	78/79(98.73)	|	28	27	0	0	1	27/97(27.84)	27/28(96.43)
GI_0505_mis2_75_45	Rd_2	|	11	10	1	0	10/10(100.0)	|	7	6	1	0	0	6/97(6.19)	    6/6(100.0)
GI_0505_mis2_75_45	Both	|	94	88	5	1	88/89(98.88)	|	35	33	1	0	1	33/97(34.02)	33/34(97.06)

GI_0504_mis2_75_45	    Rd_1	|	85	74	4	7	74/81(91.36)	|	28	25	0	2	1	25/97(25.77)	25/28(89.29)
GI_0504_mis2_75_45	    Rd_2	|	11	10	1	0	10/10(100.0)	|	7	6	1	0	0	6/97(6.19)	    6/6(100.0)
GI_0504_mis2_75_45	    Both	|	96	84	5	7	84/91(92.31)	|	35	31	1	2	1	31/97(31.96)	31/34(91.18)

GI_0504_mis1.5_75_45	Rd_1	|	96	75	16	5	75/80(93.75)	|	28	25	0	2	1	25/97(25.77)	25/28(89.29)
GI_0504_mis1.5_75_45	Rd_2	|	7	6	1	0	6/6(100.0)	    |	5	4	1	0	0	4/97(4.12)	    4/4(100.0)
GI_0504_mis1.5_75_45	Both	|	103	81	17	5	81/86(94.19)	|	33	29	1	2	1	29/97(29.9)	    29/32(90.62)

GI_0504_mis1_75_45	    Rd_1	|	95	67	16	12	67/79(84.81)	|	24	19	0	4	1	19/97(19.59)	19/24(79.17)
GI_0504_mis1_75_45	    Rd_2	|	11	11	0	0	11/11(100.0)	|	8	8	0	0	0	8/97(8.25)	    8/8(100.0)
GI_0504_mis1_75_45	    Both	|	106	78	16	12	78/90(86.67)	|	32	27	0	4	1	27/97(27.84)	27/32(84.38)






GI_0503_mis1_75_75	Rd_1	|	56	50	4	2	50/52(96.15)	|	10	9	0	0	1	9/97(9.28)	    9/10(90.0)
GI_0503_mis1_75_75	Rd_2	|	8	8	0	0	8/8(100.0)	    |	5	5	0	0	0	5/97(5.15)	    5/5(100.0)
GI_0503_mis1_75_75	Both	|	64	58	4	2	58/60(96.67)	|	15	14	0	0	1	14/97(14.43)	14/15(93.33)

GI_0503_mis1_75_45	Rd_1	|	127	64	7	56	64/120(53.33)	|	25	19	0	5	1	19/97(19.59)	19/25(76.0)
GI_0503_mis1_75_45	Rd_2	|	3	3	0	0	3/3(100.0)	    |	3	3	0	0	0	3/97(3.09)	    3/3(100.0)
GI_0503_mis1_75_45	Both	|	130	67	7	56	67/123(54.47)	|	28	22	0	5	1	22/97(22.68)	22/28(78.57)

GI_0503_mis1_45_45	Rd_1	|	225	85	11	129	85/214(39.72)	|	31	22	0	8	1	22/97(22.68)	22/31(70.97)
GI_0503_mis1_45_45	Rd_2	|	7	7	0	0	7/7(100.0)	    |	6	6	0	0	0	6/97(6.19)	    6/6(100.0)
GI_0503_mis1_45_45	Both	|	232	92	11	129	92/221(41.63)	|	37	28	0	8	1	28/97(28.87)	28/37(75.68)


GI_0420_good	    Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0420_good	    Rd_2	|	21	19	0	2	19/21(90.48)	|	17	15	0	2	0	15/97(15.46)	15/17(88.24)
GI_0420_good	    Both	|	74	71	1	2	71/73(97.26)	|	51	49	0	2	0	49/97(50.52)	49/51(96.08)




GI_0503_mis2_75_75	Rd_1	|	84	80	4	0	80/80(100.0)	|	23	23	0	0	0	23/97(23.71)	23/23(100.0)
GI_0503_mis2_75_75	Rd_2	|	10	10	0	0	10/10(100.0)	|	6	6	0	0	0	6/97(6.19)	    6/6(100.0)
GI_0503_mis2_75_75	Both	|	94	90	4	0	90/90(100.0)	|	29	29	0	0	0	29/97(29.9)	    29/29(100.0)

GI_0503_mis2_75_45	Rd_1	|	141	63	7	71	63/134(47.01)	|	27	20	0	6	1	20/97(20.62)	20/27(74.07)
GI_0503_mis2_75_45	Rd_2	|	13	11	0	2	11/13(84.62)	|	7	6	0	1	0	6/97(6.19)	    6/7(85.71)
GI_0503_mis2_75_45	Both	|	154	74	7	73	74/147(50.34)	|	34	26	0	7	1	26/97(26.8)	    26/34(76.47)

GI_0503_mis2_45_45	Rd_1	|	244	87	18	139	87/226(38.5)	|	34	24	0	9	1	24/97(24.74)	24/34(70.59)
GI_0503_mis2_45_45	Rd_2	|	16	15	1	0	15/15(100.0)	|	11	10	1	0	0	10/97(10.31)	10/10(100.0)
GI_0503_mis2_45_45	Both	|	260	102	19	139	102/241(42.32)	|	45	34	1	9	1	34/97(35.05)	34/44(77.27)

GI_0503_mis2	Rd_1	|	117	76	17	24	76/100(76.0)	|	22	20	0	2	0	20/97(20.62)	20/22(90.91)
GI_0503_mis2	Rd_2	|	13	7	0	6	7/13(53.85)	|	6	5	0	1	0	5/97(5.15)	5/6(83.33)
GI_0503_mis2	Both	|	130	83	17	30	83/113(73.45)	|	28	25	0	3	0	25/97(25.77)	25/28(89.29)

3_GI_subsample_100_1336	Refined_4	338
3_GI_subsample_75_988	Refined_4	338
3_GI_subsample_25_763	Refined_4	338
3_GI_subsample_75_1644	Refined_4	338
3_GI_subsample_75_69	Refined_4	338
3_GI_subsample_100_2066	Refined_4	338
3_GI_subsample_100_3364	Refined_4	338
3_GI_subsample_25_345	Refined_4	338
3_GI_subsample_75_1008	Refined_4	338
3_GI_subsample_10_173	Refined_4	338
3_GI_subsample_25_565	Refined_4	338
3_GI_subsample_75_2502	Refined_4	338
3_GI_subsample_75_1159	Refined_4	306

GI_0501_very_specific_clp	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	Accuracy
GI_0501_very_specific_clp	Rd_1	|	59	50	4	5	50/55(90.91)	|	11	9	0	1	1	9/97(9.28)	9/11(81.82)
GI_0501_very_specific_clp	Rd_2	|	39	18	2	19	18/37(48.65)	|	13	10	1	2	0	10/97(10.31)	10/12(83.33)
GI_0501_very_specific_clp	Both	|	98	68	6	24	68/92(73.91)	|	24	19	1	3	1	19/97(19.59)	19/23(82.61)

GI_0420_bias20	Rd_1	|	54	52	2	0	52/52(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0420_bias20	Rd_2	|	20	19	0	1	19/20(95.0)	    |	16	15	0	1	0	15/97(15.46)	15/16(93.75)
GI_0420_bias20	Both	|	74	71	2	1	71/72(98.61)	|	51	49	1	1	0	49/97(50.52)	49/50(98.0)

GI_0420_bias30	Rd_1	|	54	52	2	0	52/52(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0420_bias30	Rd_2	|	18	17	0	1	17/18(94.44)	|	14	13	0	1	0	13/97(13.4)	    13/14(92.86)
GI_0420_bias30	Both	|	72	69	2	1	69/70(98.57)	|	49	47	1	1	0	47/97(48.45)	47/48(97.92)

GI_0420_bias40	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0420_bias40	Rd_2	|	17	17	0	0	17/17(100.0)	|	13	13	0	0	0	13/97(13.4)	    13/13(100.0)
GI_0420_bias40	Both	|	70	69	1	0	69/69(100.0)	|	47	47	0	0	0	47/97(48.45)	47/47(100.0)

GI_0420_bias50	Rd_1	|	54	52	2	0	52/52(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0420_bias50	Rd_2	|	14	14	0	0	14/14(100.0)	|	11	11	0	0	0	11/97(11.34)	11/11(100.0)
GI_0420_bias50	Both	|	68	66	2	0	66/66(100.0)	|	46	45	1	0	0	45/97(46.39)	45/45(100.0)


GI_0420_godddd	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0420_godddd	Rd_2	|	21	19	0	2	19/21(90.48)	|	17	15	0	2	0	15/97(15.46)	15/17(88.24)
GI_0420_godddd	Both	|	74	71	1	2	71/73(97.26)	|	51	49	0	2	0	49/97(50.52)	49/51(96.08)

GI_0420_meta	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0420_meta	Rd_2	|	21	19	0	2	19/21(90.48)	|	17	15	0	2	0	15/97(15.46)	15/17(88.24)
GI_0420_meta	Both	|	74	71	1	2	71/73(97.26)	|	51	49	0	2	0	49/97(50.52)	49/51(96.08)

GI_0420_careful	Rd_1	|	54	52	2	0	52/52(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0420_careful	Rd_2	|	17	16	0	1	16/17(94.12)	|	12	11	0	1	0	11/97(11.34)	11/12(91.67)
GI_0420_careful	Both	|	71	68	2	1	68/69(98.55)	|	47	45	1	1	0	45/97(46.39)	45/46(97.83)

GI_0419_mira	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0419_mira	Rd_2	|	20	19	0	1	19/20(95.0)	    |	17	16	0	1	0	16/97(16.49)	16/17(94.12)
GI_0419_mira	Both	|	73	71	1	1	71/72(98.61)	|	51	50	0	1	0	50/97(51.55)	50/51(98.04)

GI_0419_mira_assemble_clp	Rd_1	|	54	52	2	0	52/52(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0419_mira_assemble_clp	Rd_2	|	21	18	1	2	18/20(90.0)	    |	18	15	1	2	0	15/97(15.46)	15/17(88.24)
GI_0419_mira_assemble_clp	Both	|	75	70	3	2	70/72(97.22)	|	53	49	2	2	0	49/97(50.52)	49/51(96.08)

GI_0418_mis2_1_no_clp	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0418_mis2_1_no_clp	Rd_2	|	16	15	0	1	15/16(93.75)	|	9	8	0	1	0	8/97(8.25)	    8/9(88.89)
GI_0418_mis2_1_no_clp	Both	|	69	67	1	1	67/68(98.53)	|	43	42	0	1	0	42/97(43.3)	    42/43(97.67)

GI_0418_mis2_clp	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0418_mis2_clp	Rd_2	|	21	16	0	5	16/21(76.19)	|	14	12	0	2	0	12/97(12.37)	12/14(85.71)
GI_0418_mis2_clp	Both	|	74	68	1	5	68/73(93.15)	|	48	46	0	2	0	46/97(47.42)	46/48(95.83)

GI_0418_mis2_no_clp	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0418_mis2_no_clp	Rd_2	|	19	18	0	1	18/19(94.74)	|	15	14	0	1	0	14/97(14.43)	14/15(93.33)
GI_0418_mis2_no_clp	Both	|	72	70	1	1	70/71(98.59)	|	49	48	0	1	0	48/97(49.48)	48/49(97.96)

GI_0418_mis1	Rd_1	|	40	39	1	0	39/39(100.0)	|	26	26	0	0	0	26/97(26.8)	    26/26(100.0)
GI_0418_mis1	Rd_2	|	14	13	0	1	13/14(92.86)	|	10	9	0	1	0	9/97(9.28)	    9/10(90.0)
GI_0418_mis1	Both	|	54	52	1	1	52/53(98.11)	|	36	35	0	1	0	35/97(36.08)	35/36(97.22)

GI_0418_mis2	Rd_1	|	54	52	2	0	52/52(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0418_mis2	Rd_2	|	19	18	0	1	18/19(94.74)	|	15	14	0	1	0	14/97(14.43)	14/15(93.33)
GI_0418_mis2	Both	|	73	70	2	1	70/71(98.59)	|	50	48	1	1	0	48/97(49.48)	48/49(97.96)

GI_0418_mis3	Rd_1	|	52	51	1	0	51/51(100.0)	|	35	35	0	0	0	35/97(36.08)	35/35(100.0)
GI_0418_mis3	Rd_2	|	24	18	1	5	18/23(78.26)	|	15	12	1	2	0	12/97(12.37)	12/14(85.71)
GI_0418_mis3	Both	|	76	69	2	5	69/74(93.24)	|	50	47	1	2	0	47/97(48.45)	47/49(95.92)




GI_0418_very_specific	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0418_very_specific	Rd_1	|	40	39	1	0	39/39(100.0)	|	25	25	0	0	0	25/97(25.77)	25/25(100.0)
GI_0418_very_specific	Rd_2	|	10	10	0	0	10/10(100.0)	|	7	7	0	0	0	7/97(7.22)	    7/7(100.0)
GI_0418_very_specific	Both	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)

GI_0418_specific	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0418_specific	Rd_1	|	51	50	1	0	50/50(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0418_specific	Rd_2	|	22	20	0	2	20/22(90.91)	|	18	16	0	2	0	16/97(16.49)	16/18(88.89)
GI_0418_specific	Both	|	73	70	1	2	70/72(97.22)	|	50	48	0	2	0	48/97(49.48)	48/50(96.0)


GI_0416_specific_Mrd2_60	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0416_specific_Mrd2_60	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0416_specific_Mrd2_60	Rd_2	|	21	20	0	1	20/21(95.24)	|	17	16	0	1	0	16/97(16.49)	16/17(94.12)
GI_0416_specific_Mrd2_60	Both	|	71	69	1	1	69/70(98.57)	|	49	48	0	1	0	48/97(49.48)	48/49(97.96)

GI_0416_specific_Mrd2_75	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0416_specific_Mrd2_75	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0416_specific_Mrd2_75	Rd_2	|	21	20	0	1	20/21(95.24)	|	17	16	0	1	0	16/97(16.49)	16/17(94.12)
GI_0416_specific_Mrd2_75	Both	|	71	69	1	1	69/70(98.57)	|	49	48	0	1	0	48/97(49.48)	48/49(97.96)

GI_0415_specific_min_M_pct45_mink75	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0415_specific_min_M_pct45_mink75	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0415_specific_min_M_pct45_mink75	Rd_2	|	11	10	0	1	10/11(90.91)	|	8	7	0	1	0	7/97(7.22)	    7/8(87.5)
GI_0415_specific_min_M_pct45_mink75	Both	|	61	59	1	1	59/60(98.33)	|	40	39	0	1	0	39/97(40.21)	39/40(97.5)


GI_0415_specific_min_M_pct45_mink49_careful	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0415_specific_min_M_pct45_mink49_careful	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0415_specific_min_M_pct45_mink49_careful	Rd_2	|	20	19	0	1	19/20(95.0)	    |	16	15	0	1	0	15/97(15.46)	15/16(93.75)
GI_0415_specific_min_M_pct45_mink49_careful	Both	|	70	68	1	1	68/69(98.55)	|	48	47	0	1	0	47/97(48.45)	47/48(97.92)







GI_0415_default2	    Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0415_default2	    Rd_1	|	53	51	2	0	51/51(100.0)	|	35	34	1	0	0	34/97(35.05)	34/34(100.0)
GI_0415_default2	    Rd_2	|	19	18	0	1	18/19(94.74)	|	15	14	0	1	0	14/97(14.43)	14/15(93.33)
GI_0415_default2	    Both	|	72	69	2	1	69/70(98.57)	|	50	48	1	1	0	48/97(49.48)	48/49(97.96)

GI_0415_specific2	    Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0415_specific2	    Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0415_specific2	    Rd_2	|	21	20	0	1	20/21(95.24)	|	17	16	0	1	0	16/97(16.49)	16/17(94.12)
GI_0415_specific2	    Both	|	71	69	1	1	69/70(98.57)	|	49	48	0	1	0	48/97(49.48)	48/49(97.96)

# 3_GI_subsample_75_963	Refined_28	260	S2	Wrong

GI_0415_very_specific2	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0415_very_specific2	Rd_1	|	38	37	1	0	37/37(100.0)	|	25	25	0	0	0	25/97(25.77)	25/25(100.0)
GI_0415_very_specific2	Rd_2	|	9	9	0	0	9/9(100.0)	    |	7	7	0	0	0	7/97(7.22)	    7/7(100.0)
GI_0415_very_specific2	Both	|	47	46	1	0	46/46(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)






            	        Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy

GI_0414_default	        Rd_1	|	58	56	1	1	56/57(98.25)	|	39	38	0	1	0	38/97(39.18)	38/39(97.44)
GI_0414_default	        Rd_2	|	15	13	1	1	13/14(92.86)	|	13	11	1	1	0	11/97(11.34)	11/12(91.67)
GI_0414_default	        Both	|	73	69	2	2	69/71(97.18)	|	52	49	1	2	0	49/97(50.52)	49/51(96.08)

GI_0414_specific	    Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0414_specific	    Rd_2	|	22	22	0	0	22/22(100.0)	|	16	16	0	0	0	16/97(16.49)	16/16(100.0)
GI_0414_specific	    Both	|	72	71	1	0	71/71(100.0)	|	48	48	0	0	0	48/97(49.48)	48/48(100.0)

GI_0414_very_specific	Rd_1	|	50	49	1	0	49/49(100.0)	|	32	32	0	0	0	32/97(32.99)	32/32(100.0)
GI_0414_very_specific	Rd_2	|	21	21	0	0	21/21(100.0)	|	15	15	0	0	0	15/97(15.46)	15/15(100.0)
GI_0414_very_specific	Both	|	71	70	1	0	70/70(100.0)	|	47	47	0	0	0	47/97(48.45)	47/47(100.0)


3_GI_subsample_100_1437	Refined_31	141	S1	Correct
3_GI_subsample_75_1079	Refined_31	103	S1	Correct
3_GI_subsample_75_1075	Refined_31	96	S1	Correct
3_GI_subsample_75_1076	Refined_31	78	S1	Correct
3_GI_subsample_75_1078	Refined_31	77	S1	Correct
3_GI_subsample_50_663	Refined_31	75	S1	Correct
3_GI_subsample_75_1074	Refined_31	68	S1	Correct
3_GI_subsample_25_402	Refined_31	40	S1	Correct
3_GI_subsample_50_712	Refined_31	38	S1	Wrong
3_GI_subsample_75_1055	Refined_31	33	S1	Wrong
3_GI_subsample_75_1162	Refined_31	27	S1	Wrong
3_GI_subsample_75_1057	Refined_31	25	S1	Wrong
3_GI_subsample_25_431	Refined_31	22	S1	Wrong
3_GI_subsample_100_1418	Refined_31	19	S1	Wrong
3_GI_subsample_25_598	Refined_31	14	S1	Wrong
3_GI_subsample_100_2329	Refined_31	12	S1	Wrong
3_GI_subsample_50_1728	Refined_31	8	S1	Unknown
3_GI_subsample_100_3714	Refined_31	8	S1	Unknown
3_GI_subsample_50_723	Refined_31	7	S1	Wrong







3_GI_subsample_10_82	Refined_19	2370	S2	Correct
3_GI_subsample_100_563	Refined_19	2224	S2	Correct
3_GI_subsample_50_207	Refined_19	2219	S2	Correct
3_GI_subsample_100_1450	Refined_53	497	S2	Correct
3_GI_subsample_100_1692	Refined_37	428	S2	Correct
3_GI_subsample_75_487	Refined_1	289	S2	Correct
3_GI_subsample_75_963	Refined_28	260	S2	Wrong       (got this in the new round) 
3_GI_subsample_100_577	Refined_8	119	S2	Correct
3_GI_subsample_100_579	Refined_8	107	S2	Correct
3_GI_subsample_100_578	Refined_8	103	S2	Correct
3_GI_subsample_100_1740	Refined_45	77	S2	Correct
3_GI_subsample_50_234	Refined_24	60	S2	Correct
3_GI_subsample_100_1404	Refined_38	56	S2	Correct
3_GI_subsample_50_647	Refined_38	50	S2	Correct
3_GI_subsample_100_1342	Refined_46	48	S2	Correct
3_GI_subsample_100_331	Refined_15	48	S2	Correct
3_GI_subsample_100_878	Refined_15	42	S2	Correct
3_GI_subsample_50_811	Refined_32	30	S2	Correct
3_GI_subsample_100_2288	Refined_89	28	S2	Correct
3_GI_subsample_75_2796	Refined_73	27	S2	Correct
3_GI_subsample_75_746	Refined_49	21	S2	Correct
3_GI_subsample_100_457	Refined_54	13	S2	Correct
3_GI_subsample_75_1549	Refined_26	13	S2	Correct


=========================================================================================

3_GI_subsample_75_963	Refined_28	256	S2  Wrong       (got this in the new round) 
NODE_4309_length_14678_cov_0.076146


0414
MarkerGene__3_GI_subsample_10_82,GenomicSeq__Refined_19,2370    =
MarkerGene__3_GI_subsample_100_563,GenomicSeq__Refined_19,2224  =
MarkerGene__3_GI_subsample_50_207,GenomicSeq__Refined_19,2219   =
MarkerGene__3_GI_subsample_100_1450,GenomicSeq__Refined_53,497  =
MarkerGene__3_GI_subsample_100_1692,GenomicSeq__Refined_37,428  =
MarkerGene__3_GI_subsample_75_487,GenomicSeq__Refined_1,289     = good
MarkerGene__3_GI_subsample_100_577,GenomicSeq__Refined_8,119    = good
MarkerGene__3_GI_subsample_100_579,GenomicSeq__Refined_8,107    =
MarkerGene__3_GI_subsample_100_578,GenomicSeq__Refined_8,103    !!! not in below
MarkerGene__3_GI_subsample_100_1740,GenomicSeq__Refined_45,77   =



0415
MarkerGene__3_GI_subsample_10_82,GenomicSeq__Refined_19,2384    =
MarkerGene__3_GI_subsample_50_207,GenomicSeq__Refined_19,2230   =
MarkerGene__3_GI_subsample_100_563,GenomicSeq__Refined_19,2227  =
MarkerGene__3_GI_subsample_100_1450,GenomicSeq__Refined_53,502  =
MarkerGene__3_GI_subsample_100_1692,GenomicSeq__Refined_37,435  =
MarkerGene__3_GI_subsample_75_487,GenomicSeq__Refined_1,291     = good
MarkerGene__3_GI_subsample_75_963,GenomicSeq__Refined_28,256    !!!
MarkerGene__3_GI_subsample_100_577,GenomicSeq__Refined_8,131    = good
MarkerGene__3_GI_subsample_100_579,GenomicSeq__Refined_8,113    =
MarkerGene__3_GI_subsample_100_1740,GenomicSeq__Refined_45,80   =







'''


########## script ##########

pwd_plot_sankey_R = '%s/file_in/get_sankey_plot.R' % wd

############################################### define file/folder name ################################################

# bin to reference
bin_vs_ref_txt              = '%s/bin_vs_ref_imag%s.txt'                    % (wd, iden_cutoff)
stats_bin_to_ref_txt        = '%s/stats_bin_to_ref_imag%s.txt'              % (wd, iden_cutoff)
stats_ref_to_bin_txt        = '%s/stats_ref_to_bin_imag%s.txt'              % (wd, iden_cutoff)
stats_bin_to_cluster_txt    = '%s/stats_bin_to_cluster_ani%s_imag%s.txt'    % (wd, drep_ani_cutoff, iden_cutoff)
stats_cluster_to_bin_txt    = '%s/stats_cluster_to_bin_ani%s_imag%s.txt'    % (wd, drep_ani_cutoff, iden_cutoff)

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

print('cluster_to_ref_dict: %s' % len(cluster_to_ref_dict) )

################################################### bin to reference ###################################################

bin_ref_connector           = '__|__'

# get ref_to_strain_dict
ref_to_strain_dict = {}
for ref in open(ref_to_strain_file):
    ref_split = ref.strip().split('\t')
    ref_to_strain_dict[ref_split[0]] = ref_split[1]

if parse_blastn_bin_vs_ref is True:
    bin_vs_ref_dict = {}
    for match in open(blastn_bin_vs_ref):
        match_split = match.strip().split('\t')
        query = match_split[0]
        query_genome = '_'.join(query.split('_')[:2])
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
        ref_to_bin_txt_handle.write('%s\t%s\t%s\n' % (each_ref, ','.join(ref_to_bin_dict[each_ref]), ref_to_strain_dict[each_ref]))
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
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict[each_ref_1]))

                for each_mag_cluster in matched_cluster_mag:
                    current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                    for each_ref_2 in current_mag_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict[each_ref_2]))
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
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict[each_ref_1]))
            for each_mag_cluster in matched_cluster_mag:
                current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                for each_ref_2 in current_mag_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict[each_ref_2]))
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


########################################################################################################################

'''

########################################################################################################################

# use 80_20 in round 2
GI_0412_specific_spades_M2_85_diff_80_20	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0412_specific_spades_M2_85_diff_80_20	Rd_2	|	23	    21	1	1	21/22(95.45)|	18	16	1	1	0	16/97(16.49)	16/17(94.12)
GI_0412_specific_spades_M2_85_diff_75_20	Rd_2	|	26	    23	1	2	23/25(92.0)	|	17	15	1	1	0	15/97(15.46)	15/16(93.75)
GI_0412_specific_spades_M2_85_diff_50_20	Rd_2	|	42	    37	1	4	37/41(90.24)|	17	15	1	1	0	15/97(15.46)	15/16(93.75)

########################################################################################################################

'''
'''


3_GI_subsample_75_1644	Refined_4	958	S1	Wrong
3_GI_subsample_75_1159	Refined_4	878	S1	Wrong

3_GI_subsample_100_1336	Refined_4	958	S1	Correct
3_GI_subsample_25_345	Refined_4	958	S1	Unknown
3_GI_subsample_75_2502	Refined_4	958	S1	Correct
3_GI_subsample_25_763	Refined_4	958	S1	Correct
3_GI_subsample_100_2066	Refined_4	958	S1	Correct
3_GI_subsample_75_69	Refined_4	958	S1	Correct
3_GI_subsample_10_173	Refined_4	958	S1	Unknown
3_GI_subsample_100_3364	Refined_4	958	S1	Correct
3_GI_subsample_75_1008	Refined_4	958	S1	Correct
3_GI_subsample_25_565	Refined_4	958	S1	Correct
3_GI_subsample_75_988	Refined_4	958	S1	Unknown
3_GI_subsample_100_2322	Refined_4	929	S1	Correct
3_GI_subsample_50_2029	Refined_4	927	S1	Unknown
3_GI_subsample_25_913	Refined_4	926	S1	Unknown
3_GI_subsample_75_2954	Refined_4	919	S1	Unknown
3_GI_subsample_75_2742	Refined_4	919	S1	Correct
3_GI_subsample_50_1564	Refined_4	898	S1	Correct
3_GI_subsample_75_2408	Refined_4	888	S1	Correct
3_GI_subsample_100_3244	Refined_4	888	S1	Correct
3_GI_subsample_25_762	Refined_4	888	S1	Correct
3_GI_subsample_100_3779	Refined_4	849	S1	Correct
3_GI_subsample_75_2776	Refined_4	836	S1	Correct

Refined_4               C46_1

3_GI_subsample_75_1644  C45_0
3_GI_subsample_75_1159  C44_1

3_GI_subsample_75_1644	3_GI_subsample_100_1336	98.442	1476
3_GI_subsample_75_1159	3_GI_subsample_75_2502	97.736	1369
3_GI_subsample_75_1159	3_GI_subsample_25_763	99.139  697
3_GI_subsample_75_1159	3_GI_subsample_100_2066	97.450	1412


3_GI_subsample_75_1159
3_GI_subsample_75_1159	3_GI_subsample_100_1336	98.060	1495
3_GI_subsample_75_1159	3_GI_subsample_75_2502	97.736	1369
3_GI_subsample_75_1159	3_GI_subsample_75_1644	99.119


3_GI_subsample_75_1644	Refined_4	958	S1	Wrong
3_GI_subsample_75_1644	3_GI_subsample_100_1336	98.442	1476
3_GI_subsample_75_1644	3_GI_subsample_75_2502	98.028	1369	24
3_GI_subsample_75_1644	3_GI_subsample_25_763	99.283	697	
3_GI_subsample_75_1644	3_GI_subsample_100_2066	97.678	1421
3_GI_subsample_75_1644	3_GI_subsample_75_69	98.252	1373

3_GI_subsample_75_988	3_GI_subsample_100_1336	99.094	1214
3_GI_subsample_75_988	3_GI_subsample_75_2502	98.848	1215
3_GI_subsample_75_988	3_GI_subsample_25_763	99.283	697
3_GI_subsample_75_988	3_GI_subsample_100_2066	98.847
3_GI_subsample_75_988	3_GI_subsample_75_69	99.012
3_GI_subsample_75_988	3_GI_subsample_100_3364	98.848






'''
print(cluster_to_bin_dict['C45_0'])  # 'Refined_47', 'Refined_50'
print(cluster_to_bin_dict['C44_1'])  # 'Refined_56', 'Refined_34'
