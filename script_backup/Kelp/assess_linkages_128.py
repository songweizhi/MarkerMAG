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
taxonomy_mag        = '%s/BH_ER_050417.bac120.summary.tsv'                              % wd
taxonomy_16s        = '%s/Kelp_SILVA128_id95_assembled_16S_uclust_0.999_classified.txt' % wd
taxonomy_16s_blca   = '%s/Kelp_SILVA128_id95_assembled_16S_uclust_0.999.fasta.blca.out' % wd
bin_id_file         = '%s/bin_id.txt'                                                   % wd
barrnap_output      = '%s/BH_ER_050417_16S.txt'                                         % wd
checkm_output       = '%s/BH_ER_050417_refined_MAG_qualities.txt'                       % wd
mean_depth_file     = '%s/BH_ER_050417_refined_bins_mean_depth.txt'                     % wd

linkage_file        = '%s/file_in/Kelp_70_3_128_combined_linkages.txt'                  % wd


'''
Linkage accuracy: Kelp_NewCigar_70_mis1_combined_linkages
Rank	Links	NA	Correct	Unknown	Total(Accuracy)
phylum	37	0	37	0	37(100.0)
class	37	0	34	0	34(91.89)
order	37	0	30	0	30(81.08)
family	37	0	28	0	28(75.68)
genus	37	0	18	7	25(67.57)

Linkage accuracy: Kelp_NewCigar_70_mis2_combined_linkages
Rank	Links	NA	Correct	Unknown	Total(Accuracy)
phylum	46	0	46	0	46(100.0)
class	46	0	43	0	43(93.48)
order	46	0	39	0	39(84.78)
family	46	0	37	0	37(80.43)
genus	46	0	27	7	34(73.91)

Linkage accuracy: Kelp_NewCigar_70_mis3_combined_linkages
Rank	Links	NA	Correct	Unknown	Total(Accuracy)
phylum	50	0	50	0	50(100.0)
class	50	0	47	0	47(94.0)
order	50	0	43	0	43(86.0)
family	50	0	41	0	41(82.0)
genus	50	0	31	7	38(76.0)

        70_mis1         70_mis2         70_mis3         70_mis4         70_mis5
phylum	37/37(100.0)	46/46(100.0)	50/50(100.0)	52/52(100.0)    54/54(100.0)
class	34/37(91.89)	43/46(93.48)	47/50(94.0)	    49/52(94.23)    51/54(94.44)
order	30/37(81.08)	39/46(84.78)	43/50(86.0) 	45/52(86.54)    47/54(87.04)
family	28/37(75.68)	37/46(80.43)	41/50(82.0) 	43/52(82.69)    45/54(83.33)
genus	18/37(48.65)	27/46(58.7)	    31/50(62.0) 	32/52(61.54)    34/54(62.96)

80_mis1
phylum	33/33(100.0)
class	31/33(93.94)
order	27/33(81.82)
family	25/33(75.76)
genus	17/33(51.52)

80_mis2
phylum	42/42(100.0)
class	40/42(95.24)
order	36/42(85.71)
family	34/42(80.95)
genus	24/42(57.14)


80_mis3
phylum	44/44(100.0)
class	42/44(95.45)
order	38/44(86.36)
family	36/44(81.82)
genus	26/44(59.09)

60_mis1
phylum	38/38(100.0)
class	35/38(92.11)
order	31/38(81.58)
family	29/38(76.32)
genus	18/38(47.37)

60_mis2
phylum	48/48(100.0)
class	45/48(93.75)
order	41/48(85.42)
family	39/48(81.25)
genus	29/48(60.42)

60_mis3
phylum	46/46(100.0)
class	43/46(93.48)
order	39/46(84.78)
family	37/46(80.43)
genus	27/46(58.7)

'''

################################################### define file name ###################################################

linkage_file_path, linkage_file_basename, linkage_file_extension = sep_path_basename_ext(linkage_file)
linkage_file_with_assessment        = '%s/%s_with_assessment%s'         % (linkage_file_path, linkage_file_basename, linkage_file_extension)
linkage_file_with_assessment_by_mag = '%s/%s_with_assessment_by_MAG%s'  % (linkage_file_path, linkage_file_basename, linkage_file_extension)


################################################### read in metadata ###################################################


bin_id_list = []
for each_bin_id in open(bin_id_file):
    bin_id_list.append(each_bin_id.strip())

bin_taxon_dict = {}
for each_bin_taxon in open(taxonomy_mag):
    if not each_bin_taxon.startswith('user_genome	classification'):
        each_bin_taxon_split = each_bin_taxon.strip().split('\t')
        bin_taxon_dict[each_bin_taxon_split[0]] = each_bin_taxon_split[1]

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
for each_16s_taxon in open(taxonomy_16s):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    s16_taxon_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1].split(' [')[0]

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
        linkage_file_with_assessment_handle.write('%s\tPhylum\tClass\tOrder\tFamily\tGenus\tMAG_taxon\t16S_taxon\t16S_taxon_BLCA\n' % each_linkage.strip())
    else:
        each_linkage_split = each_linkage.strip().split('\t')
        line_to_write = each_linkage.strip()
        s16_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        taxon_mag = bin_taxon_dict.get(mag_id, 'NA')
        taxon_16s = s16_taxon_dict.get(s16_id, 'NA')
        taxon_16s_blca = s16_taxon_blca_dict_formatted.get(s16_id, 'NA')
        taxon_mag_split = taxon_mag.split(';')
        taxon_16s_split = taxon_16s.split(';')
        total_linkage += 1

        if (taxon_mag == 'NA') or (taxon_16s == 'NA'):
            unknown_mag_or_16s_taxon += 1
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

        line_to_write += '\t%s\t%s\t%s' % (taxon_mag, taxon_16s, taxon_16s_blca)
        linkage_file_with_assessment_handle.write('%s\n' % line_to_write)


if total_linkage > 0:

    print('=============================')
    accuracy_p = float("{0:.2f}".format((correct_linkage_p + linkage_with_unknown_p) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_c = float("{0:.2f}".format((correct_linkage_c + linkage_with_unknown_c) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_o = float("{0:.2f}".format((correct_linkage_o + linkage_with_unknown_o) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_f = float("{0:.2f}".format((correct_linkage_f + linkage_with_unknown_f) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_g = float("{0:.2f}".format((correct_linkage_g + linkage_with_unknown_g) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))

    print('Linkage accuracy: %s' % linkage_file_basename)
    print('Rank\tLinks\tNA\tCorrect\tUnknown\tTotal(Accuracy)')
    print('phylum\t%s\t%s\t%s\t%s\t%s(%s)' % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_p, linkage_with_unknown_p, (correct_linkage_p + linkage_with_unknown_p), accuracy_p))
    print('class\t%s\t%s\t%s\t%s\t%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_c, linkage_with_unknown_c, (correct_linkage_c + linkage_with_unknown_c), accuracy_c))
    print('order\t%s\t%s\t%s\t%s\t%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_o, linkage_with_unknown_o, (correct_linkage_o + linkage_with_unknown_o), accuracy_o))
    print('family\t%s\t%s\t%s\t%s\t%s(%s)' % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_f, linkage_with_unknown_f, (correct_linkage_f + linkage_with_unknown_f), accuracy_f))
    print('genus\t%s\t%s\t%s\t%s\t%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_g, linkage_with_unknown_g, (correct_linkage_g + linkage_with_unknown_g), accuracy_g))
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
for each_mag in sorted(bin_id_list):

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
