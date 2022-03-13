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

wd = '/Users/songweizhi/Desktop/assess_linkages'

########## reference to cluster ##########

drep_ani_cutoff             = 95
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

########## 16S to reference ##########

perform_blastn_16s_vs_refs  = False  # True or False
combined_GI_ref_16S         = '%s/file_in/combined_GI_ref_16S.ffn'           % wd
matam_16s_seqs              = '%s/file_in/GI_128_16S_0.999.fasta'            % wd
matam_16s_blastn            = '%s/file_in/GI_128_16S_0.999_vs_ref.tab'       % wd
iden_cutoff_16s             = 99.5  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 97
# matam_16s_seqs              = '%s/file_in/GI_128_16S_0.999.QC.fasta'            % wd
# matam_16s_blastn            = '%s/file_in/GI_128_16S_0.999_QC_vs_ref.tab'       % wd

########## assessment results ##########

MarkerMAG_linkages          = '%s/GI_0909_iden990_linkages_by_genome.txt'             % wd


# linkages_from_rd1         = True
linkages_from_rd1           = False

########################################################################################################################

'''

GI_0909_iden999_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0909_iden999_linkages_by_genome.txt	Rd_1	|	61	61	0	0	61/61(100.0)	|	20	20	0	0	0	20/97(20.62)	20/20(100.0)
GI_0909_iden999_linkages_by_genome.txt	Rd_2	|	4	4	0	0	4/4(100.0)	    |	4	4	0	0	0	4/97(4.12)	    4/4(100.0)
GI_0909_iden999_linkages_by_genome.txt	Both	|	65	65	0	0	65/65(100.0)	|	24	24	0	0	0	24/97(24.74)	24/24(100.0)

GI_0909_iden995_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0909_iden995_linkages_by_genome.txt	Rd_1	|	21	20	0	1	20/21(95.24)	|	19	18	0	0	1	18/97(18.56)	18/19(94.74)
GI_0909_iden995_linkages_by_genome.txt	Rd_2	|	3	3	0	0	3/3(100.0)	    |	3	3	0	0	0	3/97(3.09)	    3/3(100.0)
GI_0909_iden995_linkages_by_genome.txt	Both	|	24	23	0	1	23/24(95.83)	|	22	21	0	0	1	21/97(21.65)	21/22(95.45)

GI_0909_iden990_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0909_iden990_linkages_by_genome.txt	Rd_1	|	20	20	0	0	20/20(100.0)	|	20	20	0	0	0	20/97(20.62)	20/20(100.0)
GI_0909_iden990_linkages_by_genome.txt	Rd_2	|	3	3	0	0	3/3(100.0)	    |	3	3	0	0	0	3/97(3.09)	    3/3(100.0)
GI_0909_iden990_linkages_by_genome.txt	Both	|	23	23	0	0	23/23(100.0)	|	23	23	0	0	0	23/97(23.71)	23/23(100.0)



GI_0819_linkages_by_genome.txt	        Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0819_linkages_by_genome.txt	        Rd_1	|	66	66	0	0	66/66(100.0)	|	20	20	0	0	0	20/97(20.62)	20/20(100.0)
GI_0819_linkages_by_genome.txt	        Rd_2	|	4	3	0	1	3/4(75.0)	    |	4	3	0	1	0	3/97(3.09)	    3/4(75.0)
GI_0819_linkages_by_genome.txt	        Both	|	70	69	0	1	69/70(98.57)	|	24	23	0	1	0	23/97(23.71)	23/24(95.83)

GI_0830_iden99_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0830_iden99_linkages_by_genome.txt	Rd_1	|	20	20	0	0	20/20(100.0)	|	20	20	0	0	0	20/97(20.62)	20/20(100.0)
GI_0830_iden99_linkages_by_genome.txt	Rd_2	|	3	3	0	0	3/3(100.0)	    |	3	3	0	0	0	3/97(3.09)	    3/3(100.0)
GI_0830_iden99_linkages_by_genome.txt	Both	|	23	23	0	0	23/23(100.0)	|	23	23	0	0	0	23/97(23.71)	23/23(100.0)

GI_0726_128_45_45_1200_diff80_iden99_mismatch0_GoodMAGs_linkages_by_genome.txt	Rd_1	|	38	38	0	0	38/38(100.0)	|	9	9	0	0	0	9/97(9.28)	9/9(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch0_GoodMAGs_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	|	0	0	0	0	0	0/97(0.0)	0/0(0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch0_GoodMAGs_linkages_by_genome.txt	Both	|	38	38	0	0	38/38(100.0)	|	9	9	0	0	0	9/97(9.28)	9/9(100.0)

GI_0726_128_45_45_1200_diff80_iden99_mismatch1_GoodMAGs_linkages_by_genome.txt	Rd_1	|	39	39	0	0	39/39(100.0)	|	11	11	0	0	0	11/97(11.34)	11/11(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch1_GoodMAGs_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	|	0	0	0	0	0	0/97(0.0)	0/0(0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch1_GoodMAGs_linkages_by_genome.txt	Both	|	39	39	0	0	39/39(100.0)	|	11	11	0	0	0	11/97(11.34)	11/11(100.0)

GI_0726_128_45_45_1200_diff80_iden99_mismatch2_GoodMAGs_linkages_by_genome.txt	Rd_1	|	59	59	0	0	59/59(100.0)	|	15	15	0	0	0	15/97(15.46)	15/15(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch2_GoodMAGs_linkages_by_genome.txt	Rd_2	|	6	6	0	0	6/6(100.0)	|	6	6	0	0	0	6/97(6.19)	6/6(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch2_GoodMAGs_linkages_by_genome.txt	Both	|	65	65	0	0	65/65(100.0)	|	21	21	0	0	0	21/97(21.65)	21/21(100.0)


GI_0726_128_45_45_1200_diff80_iden99_mismatch2_linkages_by_genome.txt	Rd_1	|	78	78	0	0	78/78(100.0)	|	20	20	0	0	0	20/97(20.62)	20/20(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch2_linkages_by_genome.txt	Rd_2	|	7	7	0	0	7/7(100.0)	|	6	6	0	0	0	6/97(6.19)	6/6(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch2_linkages_by_genome.txt	Both	|	85	85	0	0	85/85(100.0)	|	26	26	0	0	0	26/97(26.8)	26/26(100.0)

GI_0726_128_45_45_1200_diff80_iden99_mismatch1_linkages_by_genome.txt	Rd_1	|	55	55	0	0	55/55(100.0)	|	15	15	0	0	0	15/97(15.46)	15/15(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch1_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	|	0	0	0	0	0	0/97(0.0)	0/0(0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch1_linkages_by_genome.txt	Both	|	55	55	0	0	55/55(100.0)	|	15	15	0	0	0	15/97(15.46)	15/15(100.0)


# diff80, good, but decreased recovery with mismatch0
GI_0726_128_45_45_1200_diff80_iden99_mismatch0_linkages_by_genome.txt	Rd_1	|	52	52	0	0	52/52(100.0)	|	13	13	0	0	0	13/97(13.4)	13/13(100.0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch0_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	        |	0	0	0	0	0	0/97(0.0)	0/0(0)
GI_0726_128_45_45_1200_diff80_iden99_mismatch0_linkages_by_genome.txt	Both	|	52	52	0	0	52/52(100.0)	|	13	13	0	0	0	13/97(13.4)	13/13(100.0)

# diff10, not good
GI_0726_128_45_45_1200_diff10_iden99_mismatch0_linkages_by_genome.txt	Rd_1	|	86	85	0	1	85/86(98.84)	|	13	12	0	0	1	12/97(12.37)	12/13(92.31)
GI_0726_128_45_45_1200_diff10_iden99_mismatch0_linkages_by_genome.txt	Rd_2	|	0	0	0	0	0/0(0)	        |	0	0	0	0	0	0/97(0.0)	0/0(0)
GI_0726_128_45_45_1200_diff10_iden99_mismatch0_linkages_by_genome.txt	Both	|	86	85	0	1	85/86(98.84)	|	13	12	0	0	1	12/97(12.37)	12/13(92.31)


GI_0719_128_45_45_1200_linkages_by_genome.txt	Rd_1	|	59	59	0	0	59/59(100.0)	|	19	19	0	0	0	19/97(19.59)	19/19(100.0)
GI_0719_128_45_45_1200_linkages_by_genome.txt	Rd_2	|	9	9	0	0	9/9(100.0)	    |	7	7	0	0	0	7/97(7.22)	    7/7(100.0)
GI_0719_128_45_45_1200_linkages_by_genome.txt	Both	|	68	68	0	0	68/68(100.0)	|	26	26	0	0	0	26/97(26.8)	    26/26(100.0)

GI_0719_128_45_45_1200_diff10_linkages_by_genome.txt	Rd_1	|	108	96	0	12	96/108(88.89)	|	20	18	0	0	2	18/97(18.56)	18/20(90.0)
GI_0719_128_45_45_1200_diff10_linkages_by_genome.txt	Rd_2	|	14	8	3	3	8/11(72.73)	    |	6	5	0	1	0	5/97(5.15)	    5/6(83.33)
GI_0719_128_45_45_1200_diff10_linkages_by_genome.txt	Both	|	122	104	3	15	104/119(87.39)	|	26	23	0	1	2	23/97(23.71)	23/26(88.46)


GI_0719_128_60_60_1200_linkages_by_genome.txt	Rd_1	|	56	56	0	0	56/56(100.0)	|	18	18	0	0	0	18/97(18.56)	18/18(100.0)
GI_0719_128_60_60_1200_linkages_by_genome.txt	Rd_2	|	6	6	0	0	6/6(100.0)	    |	5	5	0	0	0	5/97(5.15)	    5/5(100.0)
GI_0719_128_60_60_1200_linkages_by_genome.txt	Both	|	62	62	0	0	62/62(100.0)	|	23	23	0	0	0	23/97(23.71)	23/23(100.0)

GI_0719_128_60_60_1200_diff10_linkages_by_genome.txt	Rd_1	|	88	87	0	1	87/88(98.86)	|	19	18	0	0	1	18/97(18.56)	18/19(94.74)
GI_0719_128_60_60_1200_diff10_linkages_by_genome.txt	Rd_2	|	12	8	2	2	8/10(80.0)	    |	5	5	0	0	0	5/97(5.15)	    5/5(100.0)
GI_0719_128_60_60_1200_diff10_linkages_by_genome.txt	Both	|	100	95	2	3	95/98(96.94)	|	24	23	0	0	1	23/97(23.71)	23/24(95.83)


GI_0717_128_45_45_1500_linkages_by_genome.txt	Rd_1	|	11	11	0	0	11/11(100.0)	|	9	9	0	0	0	9/97(9.28)	9/9(100.0)
GI_0717_128_45_45_1500_linkages_by_genome.txt	Rd_2	|	2	0	0	2	0/2(0.0)	|	1	0	0	1	0	0/97(0.0)	0/1(0.0)
GI_0717_128_45_45_1500_linkages_by_genome.txt	Both	|	13	11	0	2	11/13(84.62)	|	10	9	0	1	0	9/97(9.28)	9/10(90.0)

GI_0717_128_60_60_1500_linkages_by_genome.txt	Rd_1	|	10	10	0	0	10/10(100.0)	|	8	8	0	0	0	8/97(8.25)	8/8(100.0)
GI_0717_128_60_60_1500_linkages_by_genome.txt	Rd_2	|	3	0	0	3	0/3(0.0)	|	2	0	0	2	0	0/97(0.0)	0/2(0.0)
GI_0717_128_60_60_1500_linkages_by_genome.txt	Both	|	13	10	0	3	10/13(76.92)	|	10	8	0	2	0	8/97(8.25)	8/10(80.0)


GI_0715_128_45_45_linkages_by_genome.txt	Rd_1	|	94	87	4	3	87/90(96.67)	|	25	23	0	1	1	23/97(23.71)	23/25(92.0)
GI_0715_128_45_45_linkages_by_genome.txt	Rd_2	|	7	7	0	0	7/7(100.0)	    |	5	5	0	0	0	5/97(5.15)	    5/5(100.0)
GI_0715_128_45_45_linkages_by_genome.txt	Both	|	101	94	4	3	94/97(96.91)	|	30	28	0	1	1	28/97(28.87)	28/30(93.33)

GI_0715_128_60_60_linkages_by_genome.txt	Rd_1	|	61	61	0	0	61/61(100.0)	|	21	21	0	0	0	21/97(21.65)	21/21(100.0)
GI_0715_128_60_60_linkages_by_genome.txt	Rd_2	|	6	6	0	0	6/6(100.0)	    |	5	5	0	0	0	5/97(5.15)	    5/5(100.0)
GI_0715_128_60_60_linkages_by_genome.txt	Both	|	67	67	0	0	67/67(100.0)	|	26	26	0	0	0	26/97(26.8)	    26/26(100.0)




GI_0630_128_60_60_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	Accuracy
GI_0630_128_60_60_linkages_by_genome.txt	Rd_1	|	82	79	3	0	79/79(100.0)	|	26	26	0	0	0	26/97(26.8)	26/26(100.0)
GI_0630_128_60_60_linkages_by_genome.txt	Rd_2	|	3	3	0	0	3/3(100.0)	    |	3	3	0	0	0	3/97(3.09)	3/3(100.0)
GI_0630_128_60_60_linkages_by_genome.txt	Both	|	85	82	3	0	82/82(100.0)	|	29	29	0	0	0	29/97(29.9)	29/29(100.0)

GI_0602_128_60_60_good_linkages_by_genome.txt	Rd_1	|	80	74	6	0	74/74(100.0)	|	23	23	0	0	0	23/97(23.71)	23/23(100.0)
GI_0602_128_60_60_good_linkages_by_genome.txt	Rd_2	|	14	13	1	0	13/13(100.0)	|	6	5	1	0	0	5/97(5.15)	5/5(100.0)
GI_0602_128_60_60_good_linkages_by_genome.txt	Both	|	94	87	7	0	87/87(100.0)	|	29	28	1	0	0	28/97(28.87)	28/28(100.0)


GI_0602_128_60_60_linkages_by_genome.txt	Round	|	Link	Yes	NA	No	Accuracy	|	MAG	Yes	NA	No	Y/N	Recovery	    Accuracy
GI_0602_128_60_60_linkages_by_genome.txt	Rd_1	|	101	95	6	0	95/95(100.0)	|	30	30	0	0	0	30/97(30.93)	30/30(100.0)
GI_0602_128_60_60_linkages_by_genome.txt	Rd_2	|	15	14	1	0	14/14(100.0)	|	7	6	1	0	0	6/97(6.19)	    6/6(100.0)
GI_0602_128_60_60_linkages_by_genome.txt	Both	|	116	109	7	0	109/109(100.0)	|	37	36	1	0	0	36/97(37.11)	36/36(100.0)

GI_0420_bias40	Rd_1	|	53	52	1	0	52/52(100.0)	|	34	34	0	0	0	34/97(35.05)	34/34(100.0)
GI_0420_bias40	Rd_2	|	17	17	0	0	17/17(100.0)	|	13	13	0	0	0	13/97(13.4)	    13/13(100.0)
GI_0420_bias40	Both	|	70	69	1	0	69/69(100.0)	|	47	47	0	0	0	47/97(48.45)	47/47(100.0)

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
                MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tUnknown\n' % (id_16s, id_mag, link_num))

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
print('%s\tRound\t|\tLink\tYes\tNA\tNo\tAccuracy\t|\tMAG\tYes\tNA\tNo\tY/N\tRecovery\tAccuracy' % prefix)
print('%s\tRd_1\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd1,  correct_link_rd1,  unknown_link_rd1,  wrong_link_rd1,  accuracy_str_rd1_link_level,  len(linked_mag_dict_rd1), rd1_correct_num,  len(rd1_unknown_gnm), len(rd1_wrong_gnm), len(rd1_ambiguous_gnm), recovery_str_rd1, accuracy_str_rd1))
print('%s\tRd_2\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd2,  correct_link_rd2,  unknown_link_rd2,  wrong_link_rd2,  accuracy_str_rd2_link_level,  len(linked_mag_dict_rd2), rd2_correct_num,  len(rd2_unknown_gnm), len(rd2_wrong_gnm), len(rd2_ambiguous_gnm), recovery_str_rd2, accuracy_str_rd2))
print('%s\tBoth\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_both, correct_link_both, unknown_link_both, wrong_link_both, accuracy_str_both_link_level, total_linked_mag_num,      correct_num_both, unknown_gnm_both_num, wrong_gnm_both_num, ambiguous_gnm_both_num, recovery_str_both, accuracy_str_both))


########################################################################################################################

'''
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4256179/
A ANI cutoff score of >95% indicates that they belong to the same species. 
'''

from Bio import SeqIO

matam_16s_seqs_qc       = '%s/file_in/GI_128_16S_0.999.QC.fasta'            % wd
good_quality_MAG_txt    = '/Users/songweizhi/Desktop/get_depth/GI_refined_bins_complete50_contain5.txt'
MAG_with_matam_16s      = '/Users/songweizhi/Desktop/get_depth/GI_MAG_with_matam_16s.txt'

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


good_quality_MAG_txt    = '%s/GI_refined_bins_complete50_contain5.txt' % '/Users/songweizhi/Desktop/get_depth'
good_quality_MAG_set = set()
for each_mag in open(good_quality_MAG_txt):
    good_quality_MAG_set.add(each_mag.strip())

bin_to_cluster_dict_50_5 = {}
for each_match in open(stats_bin_to_cluster_txt):
    each_match_split = each_match.strip().split('\t')
    if each_match_split[0] in good_quality_MAG_set:
        bin_to_cluster_dict_50_5[each_match_split[0]] = {i for i in each_match_split[1].split(',')}

cluster_to_bin_dict_50_5 = {}
for each_bin in bin_to_cluster_dict_50_5:
    matched_clusters =  bin_to_cluster_dict_50_5[each_bin]
    for matched_cluster in matched_clusters:
        if matched_cluster not in cluster_to_bin_dict_50_5:
            cluster_to_bin_dict_50_5[matched_cluster] = {each_bin}
        else:
            cluster_to_bin_dict_50_5[matched_cluster].add(each_bin)


########################################################################################################################

ref_gnm_cp_num          = '/Users/songweizhi/Desktop/cp_num_GI/all_ref_16S_stats.txt'
provided_mag_16s_cp_num = '/Users/songweizhi/Desktop/cp_num_GI/provided_mag_16s_cp_num.txt'

# read in reference cp num
ref_to_cp_num_dict = {}
for each_line in open(ref_gnm_cp_num):
    if not each_line.startswith('Genome	Copies(16S)'):
        each_line_split = each_line.strip().split('\t')
        ref_to_cp_num_dict[each_line_split[0]] = int(each_line_split[1])

estimated_cp_num_list = []
ref_gnm_cp_num_list = []
provided_mag_16s_cp_num_handle = open(provided_mag_16s_cp_num, 'w')
provided_mag_16s_cp_num_handle.write('MAG\tCopy_num_mean\tCopy_num_sep\n')
for each_mag in bin_to_cluster_dict_50_5:
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

    ref_gnm_cp_num_list.append(ref_mean_cp_num)
    provided_mag_16s_cp_num_handle.write('%s\t%s\t(%s)/%s\n' % (each_mag, ref_mean_cp_num, '+'.join([str(i) for i in cluster_ref_cp_num_list]), len(cluster_ref_cp_num_list)))

provided_mag_16s_cp_num_handle.close()


########################################################################################################################

print('\n--------------------\n')

c99_16s = '/Users/songweizhi/Desktop/GI_128_16S_0.999_polished_min1200bp_c99.0.fasta'
mag_id  = 'Refined_4'


c99_16s_id_set = set()
for each_seq in SeqIO.parse(c99_16s, 'fasta'):
    c99_16s_id_set.add(each_seq.id)


print('MAG id: %s' % (mag_id))
mag_clusters = bin_to_cluster_dict[mag_id]
print('Cluster of %s: %s' % (mag_id, mag_clusters))

for each_cluster in mag_clusters:

    cluster_ref_gnm_set = cluster_to_ref_dict[each_cluster]
    print('Reference genome in cluster %s: %s' % (each_cluster, cluster_ref_gnm_set))


    cluster_16s_set = cluster_to_matam_16s_dict[each_cluster]
    print('Matam 16S in cluster %s: %s' % (each_cluster, cluster_16s_set))

    for each_16s in cluster_16s_set:
        if each_16s in c99_16s_id_set:
            print(each_16s)



'''
MAG id: Refined_4
Cluster of Refined_4: {'C46_1'}
Reference genome in cluster C46_1: ['OTU_97.1182.0', 'OTU_97.15599.0', 'OTU_97.16157.0', 'OTU_97.326.0', 'OTU_97.365.0']
OTU_97.1182.0	12
OTU_97.15599.0	9
OTU_97.16157.0	9
OTU_97.326.0	11
OTU_97.365.0	8

3_GI_subsample_100_1336	3_GI_subsample_100_3361	99.142	1398	10	2	41	1436	1	1398	0.0	2461	1493	1398



'''
