import os
import os
import os
import glob
import random
import argparse
import numpy as np
from Bio import SeqIO
import multiprocessing as mp
from itertools import groupby
from operator import itemgetter
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


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

wd = '/Users/songweizhi/Desktop/assess_linkages_Oral'

########## reference to cluster ##########

drep_ani_cutoff             = 97
drep_cdb_file               = '%s/file_in/Cdb_%s.csv'                                       % (wd, drep_ani_cutoff)
# ref_to_strain_file          = '%s/file_in/ref_to_strain.txt'                                % wd

########## bin to reference ##########

parse_blastn_bin_vs_ref     = False  # True or False
blastn_bin_vs_ref           = '%s/file_in/bin_vs_ref.tab'                                   % wd
iden_cutoff                 = 99.5
aln_len_cutoff              = 1500
cov_q_cutoff                = 90
min_match_length            = 524288  # 100 Kbp  102400
# mag_metadata                = '%s/file_in/MAG_metadata.txt'                                 % wd

########## 16S to reference ##########

perform_blastn_16s_vs_refs  = False  # True or False
combined_GI_ref_16S         = '%s/file_in/combined_Oral_ref_16S.ffn'                     % wd
matam_16s_seqs              = '%s/file_in/CAMI_Oral_138_16S_0.999.polished.fa'           % wd
matam_16s_blastn            = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_vs_ref.tab'   % wd
iden_cutoff_16s             = 99.0  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 87
#matam_16s_seqs              = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_min1200.QC.fa'           % wd

########## assessment results ##########

MarkerMAG_linkages          = '%s/Oral_0909_iden995_linkages_by_genome.txt'             % wd
linkages_from_rd1           = False


# MarkerMAG_linkages          = '/Users/songweizhi/Desktop/tunning_rd1_Oral/stats_combined_filtered.txt'
# MarkerMAG_linkages          = '/Users/songweizhi/Desktop/tunning_rd2_Oral/stats_GapFilling_gnm_filtered.txt'
# linkages_from_rd1           = True

'''

Oral_0909_iden999_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Accuracy
Oral_0909_iden999_linkages_by_genome.txt	Rd_1	|	28	27	1	0	27/27(100.0)	|	10	10	0	0	0	10/10(100.0)
Oral_0909_iden999_linkages_by_genome.txt	Rd_2	|	8	8	0	0	8/8(100.0)	    |	5	5	0	0	0	5/5(100.0)
Oral_0909_iden999_linkages_by_genome.txt	Both	|	36	35	1	0	35/35(100.0)	|	15	15	0	0	0	15/15(100.0)

Oral_0909_iden995_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Accuracy
Oral_0909_iden995_linkages_by_genome.txt	Rd_1	|	10	9	1	0	9/9(100.0)	    |	9	9	0	0	0	9/9(100.0)
Oral_0909_iden995_linkages_by_genome.txt	Rd_2	|	7	7	0	0	7/7(100.0)	    |	6	6	0	0	0	6/6(100.0)
Oral_0909_iden995_linkages_by_genome.txt	Both	|	17	16	1	0	16/16(100.0)	|	15	15	0	0	0	15/15(100.0)

Oral_0909_iden990_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Accuracy
Oral_0909_iden990_linkages_by_genome.txt	Rd_1	|	9	9	0	0	9/9(100.0)	    |	9	9	0	0	0	9/9(100.0)
Oral_0909_iden990_linkages_by_genome.txt	Rd_2	|	5	4	1	0	4/4(100.0)	    |	5	4	1	0	0	4/4(100.0)
Oral_0909_iden990_linkages_by_genome.txt	Both	|	14	13	1	0	13/13(100.0)	|	14	13	1	0	0	13/13(100.0)


Oral_0729_45_45_min1200_mismatch2_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Rd_1	|	17	16	1	0	16/16(100.0)	|	7	7	0	0	0	7/7(100.0)
Oral_0729_45_45_min1200_mismatch2_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Rd_2	|	19	19	0	0	19/19(100.0)	|	8	8	0	0	0	8/8(100.0)
Oral_0729_45_45_min1200_mismatch2_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Both	|	36	35	1	0	35/35(100.0)	|	15	15	0	0	0	15/15(100.0)

Oral_0729_45_45_min1200_mismatch1_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Rd_1	|	14	13	1	0	13/13(100.0)	|	5	5	0	0	0	5/5(100.0)
Oral_0729_45_45_min1200_mismatch1_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Rd_2	|	11	10	0	1	10/11(90.91)	|	3	2	0	1	0	3/3(100.0)
Oral_0729_45_45_min1200_mismatch1_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Both	|	25	23	1	1	23/24(95.83)	|	8	7	0	1	0	8/8(100.0)

Oral_0729_45_45_min1200_mismatch0_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Rd_1	|	14	13	1	0	13/13(100.0)	|	5	5	0	0	0	5/5(100.0)
Oral_0729_45_45_min1200_mismatch0_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Rd_2	|	11	10	0	1	10/11(90.91)	|	2	1	0	1	0	2/2(100.0)
Oral_0729_45_45_min1200_mismatch0_diff80_iden99_GoodMAGs_linkages_by_genome.txt	Both	|	25	23	1	1	23/24(95.83)	|	7	6	0	1	0	7/7(100.0)



Oral_0729_45_45_min1200_mismatch2_diff80_iden99_linkages_by_genome.txt	Rd_1	|	35	34	1	0	34/34(100.0)	|	10	10	0	0	0	10/10(100.0)
Oral_0729_45_45_min1200_mismatch2_diff80_iden99_linkages_by_genome.txt	Rd_2	|	25	25	0	0	25/25(100.0)	|	9	9	0	0	0	9/9(100.0)
Oral_0729_45_45_min1200_mismatch2_diff80_iden99_linkages_by_genome.txt	Both	|	60	59	1	0	59/59(100.0)	|	19	19	0	0	0	19/19(100.0)

Oral_0729_45_45_min1200_mismatch1_diff80_iden99_linkages_by_genome.txt	Rd_1	|	22	21	1	0	21/21(100.0)	|	7	7	0	0	0	7/7(100.0)
Oral_0729_45_45_min1200_mismatch1_diff80_iden99_linkages_by_genome.txt	Rd_2	|	11	10	0	1	10/11(90.91)	|	3	3	0	0	0	3/3(100.0)
Oral_0729_45_45_min1200_mismatch1_diff80_iden99_linkages_by_genome.txt	Both	|	33	31	1	1	31/32(96.88)	|	10	10	0	0	0	10/10(100.0)

Oral_0729_45_45_min1200_mismatch0_diff80_iden99_linkages_by_genome.txt	Rd_1	|	21	20	1	0	20/20(100.0)	|	6	6	0	0	0	6/6(100.0)
Oral_0729_45_45_min1200_mismatch0_diff80_iden99_linkages_by_genome.txt	Rd_2	|	11	10	0	1	10/11(90.91)	|	2	1	0	1	0	1/2(50.0)
Oral_0729_45_45_min1200_mismatch0_diff80_iden99_linkages_by_genome.txt	Both	|	32	30	1	1	30/31(96.77)	|	8	7	0	1	0	7/8(87.5)






Oral_0719_45_45_min1200_linkages_by_genome.txt	Rd_1	|	59	58	1	0	58/58(100.0)	|	14	14	0	0	0	14/14(100.0)
Oral_0719_45_45_min1200_linkages_by_genome.txt	Rd_2	|	11	10	0	1	10/11(90.91)	|	7	6	0	1	0	6/7(85.71)
Oral_0719_45_45_min1200_linkages_by_genome.txt	Both	|	70	68	1	1	68/69(98.55)	|	21	20	0	1	0	20/21(95.24)

Oral_0719_60_60_min1200_linkages_by_genome.txt	Rd_1	|	43	42	1	0	42/42(100.0)	|	13	13	0	0	0	13/13(100.0)
Oral_0719_60_60_min1200_linkages_by_genome.txt	Rd_2	|	7	6	0	1	6/7(85.71)	|	4	3	0	1	0	3/4(75.0)
Oral_0719_60_60_min1200_linkages_by_genome.txt	Both	|	50	48	1	1	48/49(97.96)	|	17	16	0	1	0	16/17(94.12)

Oral_0719_45_45_min1500_linkages_by_genome.txt	Rd_1	|	46	46	0	0	46/46(100.0)	|	13	13	0	0	0	13/13(100.0)
Oral_0719_45_45_min1500_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	        |	0	0	0	0	0	0/0(0)
Oral_0719_45_45_min1500_linkages_by_genome.txt	Both	|	46	46	0	0	46/46(100.0)	|	13	13	0	0	0	13/13(100.0)

Oral_0719_45_45_min1500_diff10_linkages_by_genome.txt	Rd_1	|	52	52	0	0	52/52(100.0)	|	13	13	0	0	0	13/13(100.0)
Oral_0719_45_45_min1500_diff10_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	        |	0	0	0	0	0	0/0(0)
Oral_0719_45_45_min1500_diff10_linkages_by_genome.txt	Both	|	52	52	0	0	52/52(100.0)	|	13	13	0	0	0	13/13(100.0)

Oral_0719_60_60_min1500_linkages_by_genome.txt	Rd_1	|	36	36	0	0	36/36(100.0)	|	12	12	0	0	0	12/12(100.0)
Oral_0719_60_60_min1500_linkages_by_genome.txt	Rd_2	|	8	8	0	0	8/8(100.0)	    |	5	5	0	0	0	5/5(100.0)
Oral_0719_60_60_min1500_linkages_by_genome.txt	Both	|	44	44	0	0	44/44(100.0)	|	17	17	0	0	0	17/17(100.0)

Oral_0718_45_45_min1200_linkages_by_genome.txt	Rd_1	|	60	59	1	0	59/59(100.0)	|	15	15	0	0	0	15/15(100.0)
Oral_0718_45_45_min1200_linkages_by_genome.txt	Rd_2	|	11	10	0	1	10/11(90.91)	|	7	6	0	1	0	6/7(85.71)
Oral_0718_45_45_min1200_linkages_by_genome.txt	Both	|	71	69	1	1	69/70(98.57)	|	22	21	0	1	0	21/22(95.45)

Oral_0718_60_60_min1500_linkages_by_genome.txt	Rd_1	|	37	37	0	0	37/37(100.0)	|	13	13	0	0	0	13/13(100.0)
Oral_0718_60_60_min1500_linkages_by_genome.txt	Rd_2	|	8	5	0	3	5/8(62.5)	|	5	4	0	1	0	4/5(80.0)
Oral_0718_60_60_min1500_linkages_by_genome.txt	Both	|	45	42	0	3	42/45(93.33)	|	18	17	0	1	0	17/18(94.44)

Oral_0717_45_45_min1200_linkages_by_genome.txt	Rd_1	|	60	59	1	0	59/59(100.0)	|	15	15	0	0	0	15/15(100.0)
Oral_0717_45_45_min1200_linkages_by_genome.txt	Rd_2	|	12	10	0	2	10/12(83.33)	|	8	6	0	2	0	6/8(75.0)
Oral_0717_45_45_min1200_linkages_by_genome.txt	Both	|	72	69	1	2	69/71(97.18)	|	23	21	0	2	0	21/23(91.3)

Oral_0715_45_45_min1200_linkages_by_genome.txt	Rd_1	|	60	59	1	0	59/59(100.0)	|	15	15	0	0	0	15/15(100.0)
Oral_0715_45_45_min1200_linkages_by_genome.txt	Rd_2	|	12	10	0	2	10/12(83.33)	|	8	6	0	2	0	7/8(87.5)
Oral_0715_45_45_min1200_linkages_by_genome.txt	Both	|	72	69	1	2	69/71(97.18)	|	23	21	0	2	0	22/23(95.7)

Oral_0715_60_60_min1200_linkages_by_genome.txt	Rd_1	|	44	43	1	0	43/43(100.0)	|	14	14	0	0	0	14/14(100.0)
Oral_0715_60_60_min1200_linkages_by_genome.txt	Rd_2	|	8	7	0	1	7/8(87.5)	    |	5	4	0	1	0	5/5(100.0)
Oral_0715_60_60_min1200_linkages_by_genome.txt	Both	|	52	50	1	1	50/51(98.04)	|	19	18	0	1	0	19/19(100.0)

Oral_0715_45_45_min1200_linkages_by_genome.txt	Rd_1	|	60	59	1	0	59/59(100.0)	|	15	15	0	0	0	15/15(100.0)
Oral_0715_45_45_min1200_linkages_by_genome.txt	Rd_2	|	12	10	0	2	10/12(83.33)	|	8	6	0	2	0	6/8(75.0)
Oral_0715_45_45_min1200_linkages_by_genome.txt	Both	|	72	69	1	2	69/71(97.18)	|	23	21	0	2	0	21/23(91.3)

Oral_0715_60_60_min1200_linkages_by_genome.txt	Rd_1	|	44	43	1	0	43/43(100.0)	|	14	14	0	0	0	14/14(100.0)
Oral_0715_60_60_min1200_linkages_by_genome.txt	Rd_2	|	8	7	0	1	7/8(87.5)	    |	5	4	0	1	0	4/5(80.0)
Oral_0715_60_60_min1200_linkages_by_genome.txt	Both	|	52	50	1	1	50/51(98.04)	|	19	18	0	1	0	18/19(94.74)

Oral_0715_60_60_min1500_linkages_by_genome.txt	Rd_1	|	37	37	0	0	37/37(100.0)	|	13	13	0	0	0	13/13(100.0)
Oral_0715_60_60_min1500_linkages_by_genome.txt	Rd_2	|	8	5	0	3	5/8(62.5)	    |	5	4	0	1	0	4/5(80.0)
Oral_0715_60_60_min1500_linkages_by_genome.txt	Both	|	45	42	0	3	42/45(93.33)	|	18	17	0	1	0	17/18(94.44)

good!!!
Oral_0713_60_60_min1500_linkages_by_genome.txt	Rd_1	|	45	40	0	5	40/45(88.89)	|	15	14	0	1	0	14/15(93.33)
Oral_0713_60_60_min1500_linkages_by_genome.txt	Rd_2	|	10	10	0	0	10/10(100.0)	|	3	3	0	0	0	3/3(100.0)
Oral_0713_60_60_min1500_linkages_by_genome.txt	Both	|	55	50	0	5	50/55(90.91)	|	18	17	0	1	0	17/18(94.44)
good!!!
Oral_0713_60_60_min1200_linkages_by_genome.txt	Rd_1	|	49	43	1	5	43/48(89.58)	|	15	14	0	1	0	14/15(93.33)
Oral_0713_60_60_min1200_linkages_by_genome.txt	Rd_2	|	4	3	0	1	3/4(75.0)	    |	3	2	0	1	0	2/3(66.67)
Oral_0713_60_60_min1200_linkages_by_genome.txt	Both	|	53	46	1	6	46/52(88.46)	|	18	16	0	2	0	16/18(88.89)

Oral_0713_45_45_min1200_linkages_by_genome.txt	Rd_1	|	93	81	1	11	81/92(88.04)	|	22	18	0	4	0	18/22(81.82)
Oral_0713_45_45_min1200_linkages_by_genome.txt	Rd_2	|	15	5	0	10	5/15(33.33)	    |	7	4	0	3	0	4/7(57.14)
Oral_0713_45_45_min1200_linkages_by_genome.txt	Both	|	108	86	1	21	86/107(80.37)	|	29	22	0	7	0	22/29(75.86)

Oral_0708_60_60_min1200_mask_linkages_by_genome.txt	Rd_1	|	50	43	1	6	43/49(87.76)	|	16	14	0	2	0	14/16(87.5)
Oral_0708_60_60_min1200_mask_linkages_by_genome.txt	Rd_2	|	16	12	0	4	12/16(75.0)	    |	4	2	0	2	0	2/4(50.0)
Oral_0708_60_60_min1200_mask_linkages_by_genome.txt	Both	|	66	55	1	10	55/65(84.62)	|	20	16	0	4	0	16/20(80.0)
Oral_0708_60_60_min1500_mask_linkages_by_genome.txt	Rd_1	|	46	40	0	6	40/46(86.96)	|	16	14	0	2	0	14/16(87.5)
Oral_0708_60_60_min1500_mask_linkages_by_genome.txt	Rd_2	|	9	1	0	8	1/9(11.11)	    |	2	1	0	1	0	1/2(50.0)
Oral_0708_60_60_min1500_mask_linkages_by_genome.txt	Both	|	55	41	0	14	41/55(74.55)	|	18	15	0	3	0	15/18(83.33)
Oral_0705_60_60_min1500_linkages_by_genome.txt	Rd_1	|	45	42	0	3	42/45(93.33)	|	16	15	0	1	0	15/16(93.75)
Oral_0705_60_60_min1500_linkages_by_genome.txt	Rd_2	|	5	5	0	0	5/5(100.0)	    |	2	2	0	0	0	2/2(100.0)
Oral_0705_60_60_min1500_linkages_by_genome.txt	Both	|	50	47	0	3	47/50(94.0)	    |	18	17	0	1	0	17/18(94.44)
Oral_0705_60_60_min1400_linkages_by_genome.txt	Rd_1	|	38	37	0	1	37/38(97.37)	|	15	14	0	1	0	14/15(93.33)
Oral_0705_60_60_min1400_linkages_by_genome.txt	Rd_2	|	17	1	0	16	1/17(5.88)	    |	3	0	0	3	0	0/3(0.0)
Oral_0705_60_60_min1400_linkages_by_genome.txt	Both	|	55	38	0	17	38/55(69.09)	|	18	14	0	4	0	14/18(77.78)
Oral_0622_60_60_polish_min1200_99_127_linkages_by_genome.txt	Rd_1	|	44	39	1	4	39/43(90.7)	    |	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1200_99_127_linkages_by_genome.txt	Rd_2	|	20	6	0	14	6/20(30.0)	    |	3	1	0	2	0	1/3(33.33)
Oral_0622_60_60_polish_min1200_99_127_linkages_by_genome.txt	Both	|	64	45	1	18	45/63(71.43)	|	19	16	0	3	0	16/19(84.21)
Oral_0622_60_60_polish_new_linkages_by_genome.txt	Rd_1	|	46	44	1	1	44/45(97.78)	|	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_new_linkages_by_genome.txt	Rd_2	|	23	9	0	14	9/23(39.13)	    |	4	2	0	2	0	2/4(50.0)
Oral_0622_60_60_polish_new_linkages_by_genome.txt	Both	|	69	53	1	15	53/68(77.94)	|	20	17	0	3	0	17/20(85.0)
Oral_0622_60_60_polish_min1200_75_127_linkages_by_genome.txt	Rd_1	|	45	43	1	1	43/44(97.73)	|	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1200_75_127_linkages_by_genome.txt	Rd_2	|	23	7	0	16	7/23(30.43)	    |	4	1	0	3	0	1/4(25.0)
Oral_0622_60_60_polish_min1200_75_127_linkages_by_genome.txt	Both	|	68	50	1	17	50/67(74.63)	|	20	16	0	4	0	16/20(80.0)
Oral_0622_60_60_polish_min1200_linkages_by_genome.txt	Rd_1	|	50	44	1	5	44/49(89.8)	    |	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1200_linkages_by_genome.txt	Rd_2	|	14	11	0	3	11/14(78.57)	|	3	1	0	2	0	1/3(33.33)
Oral_0622_60_60_polish_min1200_linkages_by_genome.txt	Both	|	64	55	1	8	55/63(87.3)	    |	19	16	0	3	0	16/19(84.21)
Oral_0622_60_60_polish_min1000_linkages_by_genome.txt	Rd_1	|	48	42	1	5	42/47(89.36)	|	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1000_linkages_by_genome.txt	Rd_2	|	25	8	0	17	8/25(32.0)	    |	5	2	0	3	0	2/5(40.0)
Oral_0622_60_60_polish_min1000_linkages_by_genome.txt	Both	|	73	50	1	22	50/72(69.44)	|	21	17	0	4	0	17/21(80.95)
Oral_0622_60_60_polish_min900_linkages_by_genome.txt	Rd_1	|	48	45	1	2	45/47(95.74)	|	16	14	0	1	1	14/16(87.5)
Oral_0622_60_60_polish_min900_linkages_by_genome.txt	Rd_2	|	33	14	0	19	14/33(42.42)	|	7	3	0	4	0	3/7(42.86)
Oral_0622_60_60_polish_min900_linkages_by_genome.txt	Both	|	81	59	1	21	59/80(73.75)	|	23	17	0	5	1	17/23(73.91)
Oral_0622_60_60_polish_linkages_by_genome.txt	Rd_1	|	63	52	8	3	52/55(94.55)	|	20	16	1	2	1	16/19(84.21)
Oral_0622_60_60_polish_linkages_by_genome.txt	Rd_2	|	14	10	0	4	10/14(71.43)	|	6	3	0	3	0	3/6(50.0)
Oral_0622_60_60_polish_linkages_by_genome.txt	Both	|	77	62	8	7	62/69(89.86)	|	26	19	1	5	1	19/25(76.0)

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
    if ('MarkerGene\tGenomicSeq\tLinkage\tRound' in each_linkage) or ('MarkerGene,GenomicSeq,Number' in each_linkage):
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
                    MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tCorrect\n' % (id_16s, id_mag, link_num))

                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Correct']
                else:
                    linkage_assessment_dict[id_mag].append('Correct')

            else:
                linkage_num_wrong += 1

                if linkages_from_rd1 is False:
                    MarkerMAG_linkages_assessed_handle.write('%s\tWrong\n' % each_linkage.strip())
                else:
                    MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tWrong\n' % (id_16s, id_mag, link_num))

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
                MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tUnknown\n' % (id_16s, id_mag, link_num))

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
        if linked_rd == 'Rd1':
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
        if linked_rd == 'Rd2':
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
print('%s\tRound\t|\tLink\tYes\tNA\tNo\tAccuracy\t|\tMAG\tYes\tNA\tNo\tY/N\tAccuracy' % prefix)
print('%s\tRd_1\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd1,  correct_link_rd1,  unknown_link_rd1,  wrong_link_rd1,  accuracy_str_rd1_link_level,  len(linked_mag_dict_rd1), rd1_correct_num,  len(rd1_unknown_gnm), len(rd1_wrong_gnm), len(rd1_ambiguous_gnm), accuracy_str_rd1))
print('%s\tRd_2\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd2,  correct_link_rd2,  unknown_link_rd2,  wrong_link_rd2,  accuracy_str_rd2_link_level,  len(linked_mag_dict_rd2), rd2_correct_num,  len(rd2_unknown_gnm), len(rd2_wrong_gnm), len(rd2_ambiguous_gnm), accuracy_str_rd2))
print('%s\tBoth\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_both, correct_link_both, unknown_link_both, wrong_link_both, accuracy_str_both_link_level, total_linked_mag_num,      correct_num_both, unknown_gnm_both_num, wrong_gnm_both_num, ambiguous_gnm_both_num, accuracy_str_both))


########################################################################################################################

'''

########################################################################################################################

CAMI_Oral_subsample_50_1551	Oral_57	2123	Rd1	Wrong

CAMI_Oral_subsample_50_1551     C125_1,C125_3,C126_1,C126_2,C126_3
Oral_57                         C126_5

# blast between 16S from Oral_57 ref gnm and Matam assemblies
OTU_97.23667.1_16S	CAMI_Oral_subsample_50_571	100.000	1542	0	0	1	1542	5	1546	0.0	2848

# blast between CAMI_Oral_subsample_50_1551 and 16S from ref gnms

BioSAK iTOL -ColorRange -lg leaf_group_raw.txt -lt Identity -out leaf_group.txt 

########################################################################################################################

'''
# cluster_id = 'C161_2'
# print('cluster to ref  : %s\t%s\t%s' % (cluster_id, len(cluster_to_ref_dict.get(cluster_id, 'NA')), cluster_to_ref_dict.get(cluster_id, 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % (cluster_id, len(cluster_to_bin_dict.get(cluster_id, 'NA')), cluster_to_bin_dict.get(cluster_id, 'NA')))
# print('cluster to Matam: %s\t%s\t%s' % (cluster_id, len(cluster_to_matam_16s_dict.get(cluster_id, 'NA')), cluster_to_matam_16s_dict.get(cluster_id, 'NA')))
# #print('cluster to MAG: %s\t%s' % ('C126_5', cluster_to_bin_dict.get('C126_5', 'NA')))

#
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_1', len(cluster_to_bin_dict.get('C126_1', 'NA')), cluster_to_bin_dict.get('C126_1', 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_2', len(cluster_to_bin_dict.get('C126_2', 'NA')), cluster_to_bin_dict.get('C126_2', 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_3', len(cluster_to_bin_dict.get('C126_3', 'NA')), cluster_to_bin_dict.get('C126_3', 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_4', len(cluster_to_bin_dict.get('C126_4', 'NA')), cluster_to_bin_dict.get('C126_4', 'NA')))
#


# print('%s\t%s' % ('CAMI_Oral_subsample_50_735', matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_735']))
# print('%s\t%s' % ('CAMI_Oral_subsample_50_818', matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_818']))
# print('%s\t%s' % ('CAMI_Oral_subsample_50_734', matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_734']))
# print('%s\t%s' % ('CAMI_Oral_subsample_10_308', matam_16s_to_cluster_dict['CAMI_Oral_subsample_10_308']))
# print('%s\t%s' % ('CAMI_Oral_subsample_50_279', matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_279']))

# print('%s\t%s' % ('CAMI_Oral_subsample_50_894', matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_894']))
# print('%s\t%s' % ('CAMI_Oral_subsample_50_884', matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_884']))
# print('%s\t%s' % ('CAMI_Oral_subsample_25_465', matam_16s_to_cluster_dict['CAMI_Oral_subsample_25_465']))
# print(bin_to_cluster_dict['Oral_28'])
# print(bin_to_cluster_dict['Oral_42'])
# print(bin_to_cluster_dict['Oral_43'])
# print(bin_to_cluster_dict['Oral_44'])
# print(bin_to_cluster_dict['Oral_46'])
# print(bin_to_cluster_dict['Oral_47'])
# print(bin_to_cluster_dict['Oral_53'])
# print(bin_to_cluster_dict['Oral_57'])
# print(bin_to_cluster_dict['Oral_65'])
# print(bin_to_cluster_dict['Oral_70'])
# print(bin_to_cluster_dict['Oral_72'])
# print(bin_to_cluster_dict['Oral_75'])


#
# print(cluster_to_matam_16s_dict)
# print(len(cluster_to_matam_16s_dict))



from Bio import SeqIO

matam_16s_seqs_qc       = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_min1200.QC.fa'            % wd
good_quality_MAG_txt    = '/Users/songweizhi/Desktop/get_depth/Oral_refined_MAGs_complete50_contain5.txt'
MAG_with_matam_16s      = '/Users/songweizhi/Desktop/get_depth/Oral_MAG_with_matam_16s.txt'

passed_qc_16s_set = set()
for each_16s in SeqIO.parse(matam_16s_seqs_qc, 'fasta'):
    passed_qc_16s_set.add(each_16s.id)


MAG_with_matam_16s_handle = open(MAG_with_matam_16s, 'w')
mag_num = 0
mag_with_matam_16s_num = 0
for each_mag in open(good_quality_MAG_txt):
    mag_id = each_mag.strip()
    mag_clusters = bin_to_cluster_dict.get(mag_id, [])

    mag_with_matam_16s = False
    for each_cluster in mag_clusters:
        current_cluster_16s = cluster_to_matam_16s_dict.get(each_cluster, [])
        for each_16s in current_cluster_16s:
            if each_16s in passed_qc_16s_set:
                mag_with_matam_16s = True
    mag_num += 1
    if mag_with_matam_16s is True:
        mag_with_matam_16s_num += 1
        MAG_with_matam_16s_handle.write(each_mag)
MAG_with_matam_16s_handle.close()

print('Good-quality MAGs: %s'                % mag_num)
print('Good-quality MAGs with Matam 16S: %s' % mag_with_matam_16s_num)


########################################################################################################################

def get_scatter_plot(num_list_1, num_list_2, png_file):

    num_arrary_1 = np.array(num_list_1)
    num_arrary_2 = np.array(num_list_2)

    fig = plt.figure(figsize=(6, 6))
    plt.margins(0)

    plt.scatter(num_arrary_1, num_arrary_2)
    plt.xlabel("Estimated copy number", fontsize=12)
    plt.ylabel("User provided copy number", fontsize=12)

    # set axis range
    plt.xlim(0, round(max(num_arrary_1) + 1))
    plt.ylim(0, round(max(num_arrary_2) + 1))

    # add fit line
    coeffs = np.polyfit(num_arrary_1, num_arrary_2, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    plt.plot(num_arrary_1, slope * num_arrary_1 + intercept)

    # get R-squared value
    p = np.poly1d(coeffs)
    yhat = p(num_arrary_1)
    ybar = np.sum(num_arrary_2)/len(num_arrary_2)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((num_arrary_2 - ybar)**2)
    r_squared = ssreg/sstot

    title_str = 'y = %sx + %s, R-squared = %s' % (float("{0:.4f}".format(slope)), float("{0:.4f}".format(intercept)), float("{0:.4f}".format(r_squared)))
    plt.title(title_str)

    # save plot
    plt.tight_layout()
    plt.savefig(png_file)
    plt.close()


# input
ref_gnm_cp_num      = '/Users/songweizhi/Desktop/cp_num_Oral/oral_799_ref_genomes_16S_stats.txt'
estimated_cp_num    = '/Users/songweizhi/Desktop/cp_num_Oral/MAG_mis2_aln75___16S_mis2_aln50___ignore_end150_lowest45_highest5_all_estimation.txt'
assessment_png      = '/Users/songweizhi/Desktop/cp_num_Oral/Oral_assessment.png'

# output
provided_mag_16s_cp_num = '/Users/songweizhi/Desktop/cp_num_Oral/provided_mag_16s_cp_num.txt'


# read in reference cp num
ref_to_cp_num_dict = {}
for each_line in open(ref_gnm_cp_num):
    if not each_line.startswith('Genome	Copies(16S)'):
        each_line_split = each_line.strip().split('\t')
        ref_to_cp_num_dict[each_line_split[0]] = int(each_line_split[1])

# read in assessments
mag_to_cp_num_dict = {}
for each_line in open(estimated_cp_num):
    if not each_line.startswith('MAG	Copy_number'):
        each_line_split = each_line.strip().split('\t')
        mag_to_cp_num_dict[each_line_split[0]] = float(each_line_split[1])

estimated_cp_num_list = []
ref_gnm_cp_num_list = []
provided_mag_16s_cp_num_handle = open(provided_mag_16s_cp_num, 'w')
provided_mag_16s_cp_num_handle.write('MAG\tCopy_num_mean\tCopy_num_sep\n')
for each_mag in mag_to_cp_num_dict:

    mag_estimated_cp_num = mag_to_cp_num_dict.get(each_mag)
    mag_clusters = bin_to_cluster_dict.get(each_mag)

    cluster_ref_list = []
    for each_cluster in mag_clusters:
        current_cluster_refs = cluster_to_ref_dict.get(each_cluster)
        for each_ref in current_cluster_refs:
            cluster_ref_list.append(each_ref)

    cluster_ref_cp_num_list = []
    cluster_ref_cp_num_set = set()
    for each_ref in cluster_ref_list:
        cluster_ref_cp_num_list.append(ref_to_cp_num_dict.get(each_ref, 0))
        cluster_ref_cp_num_set.add(ref_to_cp_num_dict.get(each_ref, 0))
    if len(cluster_ref_cp_num_set) == 1:
        ref_mean_cp_num = cluster_ref_cp_num_list[0]
    else:
        ref_mean_cp_num = sum(cluster_ref_cp_num_list)/len(cluster_ref_cp_num_list)
        ref_mean_cp_num = float("{0:.2f}".format(ref_mean_cp_num))

    #print('%s\t%s\t%s\t%s\t%s\t%s' % (each_mag, mag_estimated_cp_num, ref_mean_cp_num, '+'.join([str(i) for i in cluster_ref_cp_num_list]), ','.join(mag_clusters), ','.join(cluster_ref_list)))
    print('%s\t%s\t%s' % (each_mag, mag_estimated_cp_num, ref_mean_cp_num))

    estimated_cp_num_list.append(mag_estimated_cp_num)
    ref_gnm_cp_num_list.append(ref_mean_cp_num)

    provided_mag_16s_cp_num_handle.write('%s\t%s\t(%s)/%s\t%s\n' % (each_mag, ref_mean_cp_num, '+'.join([str(i) for i in cluster_ref_cp_num_list]), len(cluster_ref_cp_num_list), ','.join(cluster_ref_list)))

provided_mag_16s_cp_num_handle.close()


get_scatter_plot(estimated_cp_num_list, ref_gnm_cp_num_list, assessment_png)

