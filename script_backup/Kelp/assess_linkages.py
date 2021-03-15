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
taxonomy_mag        = '%s/BH_ER_050417.bac120.summary.tsv'                          % wd
taxonomy_16s        = '%s/BH_ER_050417_assembled_16S_uclust_0.999_classified.txt'   % wd
taxonomy_16s_blca   = '%s/BH_ER_050417_assembled_16S_uclust_0.999.fasta.blca.out'   % wd
bin_id_file         = '%s/bin_id.txt'                                               % wd
barrnap_output      = '%s/BH_ER_050417_16S.txt'                                     % wd
checkm_output       = '%s/BH_ER_050417_refined_MAG_qualities.txt'                   % wd
mean_depth_file     = '%s/BH_ER_050417_refined_bins_mean_depth.txt'                 % wd

linkage_file        = '%s/Kelp_0.999_combined_linkages.txt'                         % wd
linkage_file        = '%s/Kelp_0.999_aa_combined_linkages.txt'                      % wd
linkage_file        = '%s/Kelp_0.999_bbmap_combined_linkages.txt'                   % wd
linkage_file        = '%s/Kelp_0.999_bbmap_all_steps_combined_linkages.txt'         % wd
linkage_file        = '%s/Kelp_0.999_bbmap_60_40_combined_linkages.txt'             % wd
linkage_file        = '%s/Kelp_NewCigar_combined_linkages.txt'                      % wd
linkage_file        = '%s/Kelp_NewCigar_combined_linkages.txt'                      % wd
linkage_file        = '%s/Kelp_NewCigar_combined_linkages.txt'                      % wd
linkage_file        = '%s/file_in/Kelp_NewCigar_70_mis5_combined_linkages.txt'      % wd


'''
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
correct_linkage_c = 0
correct_linkage_o = 0
correct_linkage_f = 0
correct_linkage_g = 0
unknown_mag_16s_taxon = 0
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

        if (taxon_mag != 'NA') and (taxon_16s != 'NA'):

            if taxon_mag_split[1] == taxon_16s_split[1]:
                correct_linkage_p += 1
                line_to_write += '\t1'
            else:
                line_to_write += '\t0'
            if taxon_mag_split[2] == taxon_16s_split[2]:
                correct_linkage_c += 1
                line_to_write += '\t1'
            else:
                line_to_write += '\t0'
            if taxon_mag_split[3] == taxon_16s_split[3]:
                correct_linkage_o += 1
                line_to_write += '\t1'
            else:
                line_to_write += '\t0'
            if taxon_mag_split[4] == taxon_16s_split[4]:
                correct_linkage_f += 1
                line_to_write += '\t1'
            else:
                line_to_write += '\t0'

            if taxon_mag_split[5] == taxon_16s_split[5]:
                correct_linkage_g += 1
                line_to_write += '\t1'
            else:
                line_to_write += '\t0'

            total_linkage += 1

        else:
            unknown_mag_16s_taxon += 1
            line_to_write += '\t0\t0\t0\t0\t0'

        line_to_write += '\t%s\t%s\t%s' % (taxon_mag, taxon_16s, taxon_16s_blca)
        linkage_file_with_assessment_handle.write('%s\n' % line_to_write)

linkage_file_with_assessment_handle.close()

print('=============================')

print('Linkage accuracy:')
print('phylum\t%s/%s(%s)' % (correct_linkage_p, total_linkage, float("{0:.2f}".format(correct_linkage_p*100/total_linkage))))
print('class\t%s/%s(%s)'  % (correct_linkage_c, total_linkage, float("{0:.2f}".format(correct_linkage_c*100/total_linkage))))
print('order\t%s/%s(%s)'  % (correct_linkage_o, total_linkage, float("{0:.2f}".format(correct_linkage_o*100/total_linkage))))
print('family\t%s/%s(%s)' % (correct_linkage_f, total_linkage, float("{0:.2f}".format(correct_linkage_f*100/total_linkage))))
print('genus\t%s/%s(%s)'  % (correct_linkage_g, total_linkage, float("{0:.2f}".format(correct_linkage_g*100/total_linkage))))
print('Unclassified MAG/16S: %s' % unknown_mag_16s_taxon)

print('=============================')


'''
Linkage accuracy (Bowtie aa):
phylum	39/39(100.0)
class	39/39(100.0)
order	35/39(89.74)
family	31/39(79.49)
genus	23/39(58.97)
Unclassified MAG/16S: 0

Linkage accuracy (bbmap):
phylum	46/46(100.0)
class	44/46(95.65)
order	40/46(86.96)
family	37/46(80.43)
genus	28/46(60.87)
Unclassified MAG/16S: 1

Linkage accuracy (bbmap all step):
phylum	49/49(100.0)
class	49/49(100.0)
order	45/49(91.84)
family	43/49(87.76)
genus	32/49(65.31)
Unclassified MAG/16S: 2


Linkage accuracy (60_40):
phylum	53/53(100.0)
class	53/53(100.0)
order	49/53(92.45)
family	47/53(88.68)
genus	35/53(66.04)
Unclassified MAG/16S: 3




'''

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

'''

Refined_33	11.06	no	83.3	no	d__Bacteria;p__Myxococcota;c__Polyangia;o__Haliangiales;f__Haliangiaceae;g__;s__

BH_ER_050417_subsample_100_98	d__Bacteria;p__Myxococcota;c__Polyangia;o__Haliangiales;f__Haliangiaceae;g__Haliangium;s__Haliangium_ochraceum
BH_ER_050417_subsample_75_86	d__Bacteria;p__Myxococcota;c__Polyangia;o__Haliangiales;f__Haliangiaceae;g__Haliangium;s__Haliangium_ochraceum

'''