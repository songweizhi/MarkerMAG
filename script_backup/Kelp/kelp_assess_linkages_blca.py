
wd                  = '/Users/songweizhi/Desktop/MarkerMAG_wd/4_Kelp'
linkage_file        = '%s/Kelp_0.999_combined_linkages.txt'                        % wd
linkage_file        = '%s/Kelp_0.999_aa_combined_linkages.txt'                     % wd
taxonomy_mag        = '%s/BH_ER_050417.bac120.summary.tsv'                         % wd
taxonomy_16s_blca   = '%s/BH_ER_050417_assembled_16S_uclust_0.999.fasta.blca.out'  % wd


bin_taxon_dict = {}
for each_bin_taxon in open(taxonomy_mag):
    if not each_bin_taxon.startswith('user_genome	classification'):
        each_bin_taxon_split = each_bin_taxon.strip().split('\t')
        bin_taxon_dict[each_bin_taxon_split[0]] = each_bin_taxon_split[1]

s16_taxon_dict = {}
for each_16s_taxon in open(taxonomy_16s_blca):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    s16_taxon_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1].split(' [')[0]

print(s16_taxon_dict)


# total_linkage = 0
# correct_linkage_p = 0
# correct_linkage_c = 0
# correct_linkage_o = 0
# correct_linkage_f = 0
# correct_linkage_g = 0
for each_linkage in open(linkage_file):
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
        each_linkage_split = each_linkage.strip().split('\t')
        s16_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        taxon_mag = bin_taxon_dict.get(mag_id, 'NA')
        taxon_16s = s16_taxon_dict.get(s16_id, 'NA')
        print(each_linkage.strip())
        print(taxon_mag)
        print(taxon_16s)
        print()




        # taxon_mag_split = taxon_mag.split(';')
        # taxon_16s_split = taxon_16s.split(';')
        #
        # if taxon_mag_split[1] == taxon_16s_split[1]:
        #     correct_linkage_p += 1
        # if taxon_mag_split[2] == taxon_16s_split[2]:
        #     correct_linkage_c += 1
        # if taxon_mag_split[3] == taxon_16s_split[3]:
        #     correct_linkage_o += 1
        # if taxon_mag_split[4] == taxon_16s_split[4]:
        #     correct_linkage_f += 1
        # if taxon_mag_split[5] == taxon_16s_split[5]:
        #     correct_linkage_g += 1
        # total_linkage += 1


# print('Number of correct linkages:')
# print('phylum\t%s/%s' % (correct_linkage_p, total_linkage))
# print('class\t%s/%s' % (correct_linkage_c, total_linkage))
# print('order\t%s/%s' % (correct_linkage_o, total_linkage))
# print('family\t%s/%s' % (correct_linkage_f, total_linkage))
# print('genus\t%s/%s' % (correct_linkage_g, total_linkage))
#


