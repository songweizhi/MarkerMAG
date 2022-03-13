import os

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


###################################################### file in/out #####################################################

wd                  = '/Users/songweizhi/Desktop/MarkerMAG_wd/4_Kelp'
bin_id_file         = '%s/bin_id.txt'                                                                   % wd
bin_id_file_50_5    = '%s/bin_id_50_5.txt'                                                              % wd
barrnap_output      = '%s/BH_ER_050417_16S.txt'                                                         % wd
checkm_output       = '%s/BH_ER_050417_refined_MAG_qualities.txt'                                       % wd
mean_depth_file     = '%s/BH_ER_050417_refined_bins_mean_depth.txt'                                     % wd

# r95
taxonomy_mag        = '%s/BH_ER_050417.bac120.summary.tsv'                                              % wd
taxonomy_16s        = '%s/Kelp_SILVA138_id99_assembled_16S_uclust_0.999_vs_GTDB_classifications.txt'    % wd
taxonomy_16s_blca   = '%s/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.fasta.blca.out'                 % wd
gtdb_gnm_to_ssu_txt = '%s/GTDB_ssu_all_r95_gnm_to_ssu.txt'                                              % wd

# r202
taxonomy_mag        = '%s/BH_ER_050417.bac120.summary.r202.tsv'                                         % wd
taxonomy_16s        = '%s/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.QC_vs_GTDB_r202_classifications.txt'    % wd
taxonomy_16s_blca   = '%s/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.QC_blca_vs_gtdb_r202.txt'                 % wd
gtdb_gnm_to_ssu_txt = '%s/GTDB_ssu_all_r202_gnm_to_ssu.txt'                                             % wd

linkage_file = '/Users/songweizhi/Desktop/assess_linkages_Kelp/Kelp_0907_linkages_by_genome.txt'


'''
=============================
Linkage accuracy: Kelp_0907_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	13	0	13	0	13/13(100.0)
class	13	0	12	0	12/13(92.31)
order	13	0	12	0	12/13(92.31)
family	13	0	10	0	10/13(76.92)
genus	13	0	7	2	9/13(69.23)
=============================
linked MAGs: 14

=============================
Linkage accuracy: Kelp_0726_45_45_1200_iden99_mismatch0_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	14	0	14	0	14/14(100.0)
class	14	0	12	0	12/14(85.71)
order	14	0	12	0	12/14(85.71)
family	14	0	10	0	10/14(71.43)
genus	14	0	8	2	10/14(71.43)
=============================
linked MAGs: 11


=============================
Linkage accuracy: Kelp_0726_45_45_1200_iden99_mismatch1_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	19	0	19	0	19/19(100.0)
class	19	0	17	0	17/19(89.47)
order	19	0	17	0	17/19(89.47)
family	19	0	15	0	15/19(78.95)
genus	19	0	10	2	12/19(63.16)
=============================
linked MAGs: 14


=============================
Linkage accuracy: Kelp_0726_45_45_1200_iden99_mismatch2_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	19	0	19	0	19/19(100.0)
class	19	0	17	0	17/19(89.47)
order	19	0	17	0	17/19(89.47)
family	19	0	15	0	15/19(78.95)
genus	19	0	10	2	12/19(63.16)
=============================
linked MAGs: 14








=============================
Linkage accuracy: Kelp_0721_60_60_900_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	38	0	38	0	38/38(100.0)
class	38	0	36	0	36/38(94.74)
order	38	0	34	0	34/38(89.47)
family	38	0	32	0	32/38(84.21)
genus	38	0	19	8	27/38(71.05)
=============================
linked MAGs: 24


=============================
Linkage accuracy: Kelp_0721_45_45_900_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	38	0	38	0	38/38(100.0)
class	38	0	36	0	36/38(94.74)
order	38	0	34	0	34/38(89.47)
family	38	0	32	0	32/38(84.21)
genus	38	0	20	7	27/38(71.05)
=============================
linked MAGs: 22



Linkage accuracy: Kelp_0602_60_60_linkages_by_genome
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	59	1	58	0	58/58(100.0)
class	59	1	56	0	56/58(96.55)
order	59	1	54	0	54/58(93.1)
family	59	1	49	0	49/58(84.48)
genus	59	1	32	12	44/58(75.86)

Linkage accuracy: Kelp_0413_very_sensitive_identified_linkages_genome_level
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	44	3	41	0	41/41(100.0)
class	44	3	40	0	40/41(97.56)
order	44	3	37	0	37/41(90.24)
family	44	3	34	0	34/41(82.93)
genus	44	3	23	6	29/41(70.73)
linked MAGs: 31

Linkage accuracy: Kelp_0413_sensitive_identified_linkages_genome_level
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	45	3	42	0	42/42(100.0)
class	45	3	41	0	41/42(97.62)
order	45	3	38	0	38/42(90.48)
family	45	3	35	0	35/42(83.33)
genus	45	3	22	8	30/42(71.43)
linked MAGs: 31

Linkage accuracy: Kelp_0413_default_identified_linkages_genome_level
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	39	2	37	0	37/37(100.0)
class	39	2	36	0	36/37(97.3)
order	39	2	33	0	33/37(89.19)
family	39	2	30	0	30/37(81.08)
genus	39	2	19	7	26/37(70.27)
linked MAGs: 29

Linkage accuracy: Kelp_0413_specific_identified_linkages_genome_level
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	35	4	31	0	31/31(100.0)
class	35	4	30	0	30/31(96.77)
order	35	4	27	0	27/31(87.1)
family	35	4	24	0	24/31(77.42)
genus	35	4	15	6	21/31(67.74)
linked MAGs: 25

Linkage accuracy: Kelp_0413_very_specific_identified_linkages_genome_level
Rank	All	NA	Correct	Unknown	Total(Accuracy)
phylum	33	2	31	0	31/31(100.0)
class	33	2	30	0	30/31(96.77)
order	33	2	27	0	27/31(87.1)
family	33	2	24	0	24/31(77.42)
genus	33	2	15	7	22/31(70.97)
linked MAGs: 24

'''


################################################### define file name ###################################################

linkage_file_path, linkage_file_basename, linkage_file_extension = sep_path_basename_ext(linkage_file)
linkage_file_with_assessment        = '%s/%s_with_assessment%s'         % (linkage_file_path, linkage_file_basename, linkage_file_extension)
linkage_file_with_assessment_by_mag = '%s/%s_with_assessment_by_MAG%s'  % (linkage_file_path, linkage_file_basename, linkage_file_extension)


################################################### read in metadata ###################################################

gtdb_gnm_to_ssu_dict = {}
for each_gnm in open(gtdb_gnm_to_ssu_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    ssu_id_list = each_gnm_split[1].split(',')
    gtdb_gnm_to_ssu_dict[gnm_id] = ssu_id_list

bin_id_list = []
for each_bin_id in open(bin_id_file):
    bin_id_list.append(each_bin_id.strip())


bin_id_list_50_5 = []
for each_bin_id in open(bin_id_file_50_5):
    bin_id_list_50_5.append(each_bin_id.strip())

bin_closest_ref_dict = {}
bin_taxon_dict = {}
for each_bin_taxon in open(taxonomy_mag):
    if not each_bin_taxon.startswith('user_genome	classification'):
        each_bin_taxon_split = each_bin_taxon.strip().split('\t')
        closest_ref = each_bin_taxon_split[7]
        bin_taxon_dict[each_bin_taxon_split[0]] = each_bin_taxon_split[1]
        bin_closest_ref_dict[each_bin_taxon_split[0]] = closest_ref

mag_completeness_dict = {}
for each_mag_quality in open(checkm_output):
    each_mag_quality_split = each_mag_quality.strip().split(',')
    mag_completeness_dict[each_mag_quality_split[0]] = each_mag_quality_split[1]

bin_depth_dict = {}
for each_bin_depth in open(mean_depth_file):
    each_bin_depth_split = each_bin_depth.strip().split('\t')
    bin_depth_dict[each_bin_depth_split[0]] = each_bin_depth_split[2]

barrnap_op_dict = {}
for each_barrnap_op in open(barrnap_output):
    if not each_barrnap_op.startswith('Genome	16S	Length'):
        each_barrnap_op_split = each_barrnap_op.strip().split('\t')
        barrnap_op_dict[each_barrnap_op_split[0]] = each_barrnap_op_split[2]

s16_taxon_dict = {}
s16_iden_dict = {}
s16_aln_dict = {}
for each_16s_taxon in open(taxonomy_16s):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    id_16s = each_16s_taxon_split[0]
    taxon_best_hit = each_16s_taxon_split[4].split(' [')[0]
    iden_with_best_match = each_16s_taxon_split[2]
    aln_with_best_match = each_16s_taxon_split[3]
    s16_taxon_dict[id_16s] = taxon_best_hit
    s16_iden_dict[id_16s]  = iden_with_best_match
    s16_aln_dict[id_16s]   = aln_with_best_match

s16_taxon_blca_dict = {}
for each_16s_taxon in open(taxonomy_16s_blca):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    s16_taxon_blca_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1]

s16_taxon_blca_dict_formatted = {}
for each_16s in s16_taxon_blca_dict:
    taxon_blca_raw = s16_taxon_blca_dict[each_16s]
    formatted_taxon_str = 'Unclassified'
    if taxon_blca_raw != 'Unclassified':
        taxon_blca_raw_split_1 = taxon_blca_raw.strip().split(':')[1:]
        formatted_taxon_list = []
        for each_str in taxon_blca_raw_split_1:
            each_str_split = each_str.split(';')
            taxon_with_confidence = '%s(%s)' % (each_str_split[0], each_str_split[1][:5])
            formatted_taxon_list.append(taxon_with_confidence)
        formatted_taxon_str = ';'.join(formatted_taxon_list)
    s16_taxon_blca_dict_formatted[each_16s] = formatted_taxon_str


###################################################### by linkage ######################################################

linkage_file_with_assessment_handle = open(linkage_file_with_assessment, 'w')
total_linkage = 0
correct_linkage_p = 0
linkage_with_unknown_p = 0
correct_linkage_c = 0
linkage_with_unknown_c = 0
correct_linkage_o = 0
linkage_with_unknown_o = 0
correct_linkage_f = 0
linkage_with_unknown_f = 0
correct_linkage_g = 0
linkage_with_unknown_g = 0
correct_linkage_s = 0
linkage_with_unknown_s = 0
unknown_mag_or_16s_taxon = 0
for each_linkage in open(linkage_file):
    if each_linkage.startswith('MarkerGene	GenomicSeq'):
        linkage_file_with_assessment_handle.write('%s\tPhylum\tClass\tOrder\tFamily\tGenus\tIdentity\tAln_len\tclosest_ref\tclosest_ref_16s_num\tMAG_taxon\t16S_taxon\t16S_taxon_BLCA\n' % each_linkage.strip())
    else:
        each_linkage_split = each_linkage.strip().split('\t')
        line_to_write = each_linkage.strip()
        s16_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        taxon_mag = bin_taxon_dict.get(mag_id, 'NA')
        closest_ref = bin_closest_ref_dict.get(mag_id, 'NA')
        closest_ref_16s_num = len(gtdb_gnm_to_ssu_dict.get(closest_ref, []))
        taxon_16s = s16_taxon_dict.get(s16_id, 'NA')
        iden_16s = s16_iden_dict.get(s16_id, 'NA')
        aln_len_16s = s16_aln_dict.get(s16_id, 'NA')
        taxon_16s_blca = s16_taxon_blca_dict_formatted.get(s16_id, 'NA')
        taxon_mag_split = taxon_mag.split(';')
        taxon_16s_split = taxon_16s.split(';')
        total_linkage += 1

        if (taxon_mag == 'NA') or (taxon_16s == 'NA'):
            unknown_mag_or_16s_taxon += 1
            line_to_write += '\tna\tna\tna\tna\tna'
        else:
            taxon_mag_p = taxon_mag_split[1]
            taxon_mag_c = taxon_mag_split[2]
            taxon_mag_o = taxon_mag_split[3]
            taxon_mag_f = taxon_mag_split[4]
            taxon_mag_g = taxon_mag_split[5]
            taxon_mag_s = taxon_mag_split[6]
            taxon_16s_p = taxon_16s_split[1]
            taxon_16s_c = taxon_16s_split[2]
            taxon_16s_o = taxon_16s_split[3]
            taxon_16s_f = taxon_16s_split[4]
            taxon_16s_g = taxon_16s_split[5]
            taxon_16s_s = taxon_16s_split[6]
            already_wrong = False

            if taxon_mag_p == taxon_16s_p:
                correct_linkage_p += 1
                line_to_write += '\t1'
            else:
                if taxon_mag_p == 'p__':
                    linkage_with_unknown_p += 1
                    line_to_write += '\tna'
                else:
                    already_wrong = True
                    line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_c == taxon_16s_c:
                    correct_linkage_c += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_c == 'c__':
                        linkage_with_unknown_c += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_o == taxon_16s_o:
                    correct_linkage_o += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_o == 'o__':
                        linkage_with_unknown_o += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_f == taxon_16s_f:
                    correct_linkage_f += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_f == 'f__':
                        linkage_with_unknown_f += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_g == taxon_16s_g:
                    correct_linkage_g += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_g == 'g__':
                        linkage_with_unknown_g += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

        line_to_write += '\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (iden_16s, aln_len_16s, closest_ref, closest_ref_16s_num, taxon_mag, taxon_16s, taxon_16s_blca)
        linkage_file_with_assessment_handle.write('%s\n' % line_to_write)


if total_linkage > 0:

    print(total_linkage)
    print()
    print('=============================')
    accuracy_p = float("{0:.2f}".format((correct_linkage_p + linkage_with_unknown_p) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_c = float("{0:.2f}".format((correct_linkage_c + linkage_with_unknown_c) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_o = float("{0:.2f}".format((correct_linkage_o + linkage_with_unknown_o) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_f = float("{0:.2f}".format((correct_linkage_f + linkage_with_unknown_f) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_g = float("{0:.2f}".format((correct_linkage_g + linkage_with_unknown_g) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))

    print('Linkage accuracy: %s' % linkage_file_basename)
    print('Rank\tAll\tNA\tCorrect\tUnknown\tTotal(Accuracy)')
    print('phylum\t%s\t%s\t%s\t%s\t%s/%s(%s)' % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_p, linkage_with_unknown_p, (correct_linkage_p + linkage_with_unknown_p), (total_linkage - unknown_mag_or_16s_taxon), accuracy_p))
    print('class\t%s\t%s\t%s\t%s\t%s/%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_c, linkage_with_unknown_c, (correct_linkage_c + linkage_with_unknown_c), (total_linkage - unknown_mag_or_16s_taxon), accuracy_c))
    print('order\t%s\t%s\t%s\t%s\t%s/%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_o, linkage_with_unknown_o, (correct_linkage_o + linkage_with_unknown_o), (total_linkage - unknown_mag_or_16s_taxon), accuracy_o))
    print('family\t%s\t%s\t%s\t%s\t%s/%s(%s)' % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_f, linkage_with_unknown_f, (correct_linkage_f + linkage_with_unknown_f), (total_linkage - unknown_mag_or_16s_taxon), accuracy_f))
    print('genus\t%s\t%s\t%s\t%s\t%s/%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_g, linkage_with_unknown_g, (correct_linkage_g + linkage_with_unknown_g), (total_linkage - unknown_mag_or_16s_taxon), accuracy_g))
    print('=============================')

######################################################## by MAG ########################################################

linked_mag_list = set()
for each_linkage in open(linkage_file):
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
        each_linkage_split = each_linkage.strip().split('\t')
        mag_id = each_linkage_split[1]
        linked_mag_list.add(mag_id)

linkage_file_with_assessment_by_mag_handle = open(linkage_file_with_assessment_by_mag, 'w')
linkage_file_with_assessment_by_mag_handle.write('MAG\tdepth\t16S\tcompleteness\tlinked\ttaxon\n')
for each_mag in sorted(bin_id_list_50_5):

    mag_depth = '0'
    if each_mag in bin_depth_dict:
        mag_depth = bin_depth_dict[each_mag]

    mag_taxon = 'unknown'
    if each_mag in bin_taxon_dict:
        mag_taxon = bin_taxon_dict[each_mag]

    with_16s_in_mag = 'no'
    if each_mag in barrnap_op_dict:
        with_16s_in_mag = barrnap_op_dict[each_mag]

    linked_to_16s = 'no'
    if each_mag in linked_mag_list:
        linked_to_16s = 'yes'

    mag_cpl = 0
    if each_mag in mag_completeness_dict:
        mag_cpl = mag_completeness_dict[each_mag]

    linkage_file_with_assessment_by_mag_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (each_mag, mag_depth, with_16s_in_mag, mag_cpl, linked_to_16s, mag_taxon))
linkage_file_with_assessment_by_mag_handle.close()

print('linked MAGs: %s' % len(linked_mag_list))

'''

Kelp_SILVA138_id99_subsample_75_306		Refined_20	42	S1	1	0	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__HTCC2089;g__SYFY01;s__	d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Erythrobacter_A;s__Erythrobacter_A_sp002895025	Bacteria(100.0);Proteobacteria(100.0);Alphaproteobacteria(88.5);Sphingomonadales(88.5);Sphingomonadaceae(88.5);Erythrobacter_A(88.5);Erythrobacter_A sp002895025(88.5)
Kelp_SILVA138_id99_subsample_100_365	Refined_20	23	S1	1	0	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__HTCC2089;g__SYFY01;s__	d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Sphingomonadales;f__Sphingomonadaceae;g__Erythrobacter_A;s__Erythrobacter_A_sp002895025	Bacteria(100.0);Proteobacteria(100.0);Alphaproteobacteria(85.0);Sphingomonadales(85.0);Sphingomonadaceae(85.0);Erythrobacter_A(85.0);Erythrobacter_A sp002895025(85.0)
Kelp_SILVA138_id99_subsample_25_212		Refined_20	11	S1	1	1	1	1	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__HTCC2089;g__SYFY01;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__HTCC2089;g__SZUA-521;s__SZUA-521_sp003248125	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(55.50);Pseudomonadales(55.50);HTCC2089(55.50);SZUA-521(55.50);SZUA-521 sp003248125(55.50)

Kelp_SILVA138_id99_subsample_75_668		Refined_26	91	S1	1	1	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__GCA-1730015;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Sedimenticolaceae;g__QGON01;s__QGON01_sp003660235	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Chromatiales(81.33);Sedimenticolaceae(81.33);QGON01(81.33);QGON01 sp003660235(80.99)
Kelp_SILVA138_id99_subsample_50_553		Refined_26	20	S1	1	1	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__GCA-1730015;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Sedimenticolaceae;g__QGON01;s__QGON01_sp003660235	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Granulosicoccales(78.33);Granulosicoccaceae(78.33);Granulosicoccus(78.33);Granulosicoccus antarcticus(78.33)
Kelp_SILVA138_id99_subsample_10_229		Refined_26	42	S1	1	1	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__GCA-1730015;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Sedimenticolaceae;g__QGON01;s__QGON01_sp003660235	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Chromatiales(95.75);Sedimenticolaceae(95.75);QGON01(95.75);QGON01 sp003660235(91.50)
Kelp_SILVA138_id99_subsample_75_314		Refined_26	17	S1	1	1	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__GCA-1730015;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Sedimenticolaceae;g__QGON01;s__QGON01_sp003660235	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Chromatiales(82.0);Sedimenticolaceae(82.0);QGON01(82.0);QGON01 sp003660235(45.0)
Kelp_SILVA138_id99_subsample_25_393		Refined_26	9	S1	1	1	0	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__GCA-1730015;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Sedimenticolaceae;g__QGON01;s__QGON01_sp003660235	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Granulosicoccales(51.0);Granulosicoccaceae(51.0);Granulosicoccus(51.0);Granulosicoccus antarcticus(51.0)
Kelp_SILVA138_id99_subsample_100_789	Refined_26	13	S1	1	1	1	1	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__GCA-1730015;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Granulosicoccales;f__Granulosicoccaceae;g__Granulosicoccus;s__Granulosicoccus_antarcticus	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Granulosicoccales(97.0);Granulosicoccaceae(97.0);Granulosicoccus(97.0);Granulosicoccus antarcticus(97.0)

Kelp_SILVA138_id99_subsample_10_136		Refined_57	154	S1	1	1	1	0	0	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Porticoccaceae;g__;s__	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Cellvibrionaceae;g__Teredinibacter;s__Teredinibacter_sp000966245	Bacteria(100.0);Proteobacteria(100.0);Gammaproteobacteria(100.0);Pseudomonadales(100.0);Cellvibrionaceae(78.0);Teredinibacter(60.0);Teredinibacter sp003634075(60.0)

Kelp_SILVA138_id99_subsample_50_115		Refined_19	44	S1	1	1	1	0	0	d__Bacteria;p__Planctomycetota;c__UBA1135;o__UBA1135;f__GCA-002686595;g__;s__	d__Bacteria;p__Planctomycetota;c__UBA1135;o__UBA1135;f__UBA1135;g__GCA-2705055;s__GCA-2705055_sp002705055	Unclassified

Kelp_SILVA138_id99_subsample_50_418		Refined_31	28	S2	1	1	1	0	0	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__SHLQ01;g__SHLQ01;s__	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__MedAcidi-G1;g__S20-B6;s__S20-B6_sp002346745	Bacteria(100.0);Actinobacteriota(100.0);Acidimicrobiia(100.0);Acidimicrobiales(100.0);MedAcidi-G1(100.0);MedAcidi-G1(60.5);MedAcidi-G1 sp002713545(60.5)
Kelp_SILVA138_id99_subsample_75_575		Refined_31	28	S2	1	1	1	0	0	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__SHLQ01;g__SHLQ01;s__	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__MedAcidi-G1;g__S20-B6;s__S20-B6_sp002699725	Bacteria(100.0);Actinobacteriota(100.0);Acidimicrobiia(100.0);Acidimicrobiales(100.0);MedAcidi-G1(100.0);UBA3125(78.33);UBA3125 sp002687745(78.33)
Kelp_SILVA138_id99_subsample_100_694	Refined_31	28	S2	1	1	1	0	0	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__SHLQ01;g__SHLQ01;s__	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__MedAcidi-G1;g__S20-B6;s__S20-B6_sp002346745	Bacteria(100.0);Actinobacteriota(100.0);Acidimicrobiia(100.0);Acidimicrobiales(100.0);MedAcidi-G1(87.0);MedAcidi-G1(81.0);MedAcidi-G1 sp003214465(81.0)
Kelp_SILVA138_id99_subsample_75_493		Refined_31	3	S2	1	1	1	0	0	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__SHLQ01;g__SHLQ01;s__	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__MedAcidi-G1;g__MedAcidi-G1;s__MedAcidi-G1_sp002713545	Bacteria(100.0);Actinobacteriota(100.0);Acidimicrobiia(100.0);Acidimicrobiales(100.0);MedAcidi-G1(100.0);MedAcidi-G1(72.33);MedAcidi-G1 sp002713545(72.33)
Kelp_SILVA138_id99_subsample_50_487		Refined_31	8	S2	1	1	1	0	0	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__SHLQ01;g__SHLQ01;s__	d__Bacteria;p__Actinobacteriota;c__Acidimicrobiia;o__Acidimicrobiales;f__MedAcidi-G1;g__S20-B6;s__S20-B6_sp002346745	Bacteria(100.0);Actinobacteriota(100.0);Acidimicrobiia(100.0);Acidimicrobiales(100.0);MedAcidi-G1(57.5);MedAcidi-G3(42.5);MedAcidi-G3 sp000817105(42.5)


'''