
wd                  = '/Users/songweizhi/Desktop/MarkerMAG_wd/4_Kelp'
linkage_file        = '%s/Kelp_0.999_combined_linkages.txt'                        % wd
#linkage_file        = '%s/Kelp_0.999_aa_combined_linkages.txt'                     % wd
linkage_file        = '%s/Kelp_0.999_bbmap_combined_linkages.txt'                     % wd
taxonomy_mag        = '%s/BH_ER_050417.bac120.summary.tsv'                         % wd
taxonomy_16s        = '%s/BH_ER_050417_assembled_16S_uclust_0.999_classified.txt'  % wd
#taxonomy_16s_blca   = '%s/BH_ER_050417_assembled_16S_uclust_0.999.fasta.blca.out'  % wd
with_blca           = True

bin_taxon_dict = {}
for each_bin_taxon in open(taxonomy_mag):
    if not each_bin_taxon.startswith('user_genome	classification'):
        each_bin_taxon_split = each_bin_taxon.strip().split('\t')
        bin_taxon_dict[each_bin_taxon_split[0]] = each_bin_taxon_split[1]

s16_taxon_dict = {}
for each_16s_taxon in open(taxonomy_16s):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    s16_taxon_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1].split(' [')[0]

# s16_taxon_blca_dict = {}
# for each_16s_taxon in open(taxonomy_16s_blca):
#     each_16s_taxon_split = each_16s_taxon.strip().split('\t')
#     s16_taxon_blca_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1]


total_linkage = 0
correct_linkage_p = 0
correct_linkage_c = 0
correct_linkage_o = 0
correct_linkage_f = 0
correct_linkage_g = 0
for each_linkage in open(linkage_file):
    if each_linkage.startswith('MarkerGene	GenomicSeq'):
        print('%s\tPhylum\tClass\tOrder\tFamily\tGenus\tMAG_taxon\t16S_taxon' % each_linkage.strip())
    else:
        each_linkage_split = each_linkage.strip().split('\t')
        line_to_write = each_linkage.strip()
        s16_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        taxon_mag = bin_taxon_dict.get(mag_id, 'NA')
        taxon_16s = s16_taxon_dict.get(s16_id, 'NA')
        #taxon_16s_blca = s16_taxon_blca_dict.get(s16_id, 'NA')
        taxon_mag_split = taxon_mag.split(';')
        taxon_16s_split = taxon_16s.split(';')
        #taxon_16s_blca_split = taxon_16s_blca.split(';')
        #print(taxon_16s)
        print(each_linkage.strip())
        print(taxon_mag_split)
        print(taxon_16s_split)
        #print('%s\t%s' % (len(taxon_16s_split), len(taxon_16s_blca_split)))

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
                pass
                # print(each_linkage.strip())
                # print('MAG\t%s' % taxon_mag)
                # print('16S\t%s' % taxon_16s)
                # print('blca\t%s' % taxon_16s_blca)
                # print()

            if taxon_mag_split[5] == taxon_16s_split[5]:
                correct_linkage_g += 1
                line_to_write += '\t1'
            else:
                line_to_write += '\t0'

        total_linkage += 1

        line_to_write += '\t%s\t%s' % (taxon_mag, taxon_16s)

        print()


print('Linkage accuracy:')
print('phylum\t%s/%s(%s)' % (correct_linkage_p, total_linkage, float("{0:.2f}".format(correct_linkage_p*100/total_linkage))))
print('class\t%s/%s(%s)' % (correct_linkage_c, total_linkage, float("{0:.2f}".format(correct_linkage_c*100/total_linkage))))
print('order\t%s/%s(%s)' % (correct_linkage_o, total_linkage, float("{0:.2f}".format(correct_linkage_o*100/total_linkage))))
print('family\t%s/%s(%s)' % (correct_linkage_f, total_linkage, float("{0:.2f}".format(correct_linkage_f*100/total_linkage))))
print('genus\t%s/%s(%s)' % (correct_linkage_g, total_linkage, float("{0:.2f}".format(correct_linkage_g*100/total_linkage))))

'''
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_0.999_aa_MarkerMAG_wd/Kelp_0.999_aa_step_1_wd/BH_ER_050417_assembled_16S_uclust_0.999_index

# currently used
bowtie2 -x BH_ER_050417_assembled_16S_uclust_0.999 -1 ../../../Kelp_R1.fasta -2 ../../../Kelp_R2.fasta -S test_1.sam -f -a --local --no-unal --quiet --threads 12

#
bowtie2 -x BH_ER_050417_assembled_16S_uclust_0.999 -1 ../../../Kelp_R1.fasta -2 ../../../Kelp_R2.fasta -S test_1_very_sensitive.sam -f -a --local --no-unal --quiet --threads 12 --very-sensitive-local
bowtie2 -x BH_ER_050417_assembled_16S_uclust_0.999 -1 ../../../Kelp_R1.fasta -2 ../../../Kelp_R2.fasta -S test_1_very_fast.sam -f -a --local --no-unal --quiet --threads 12 --very-fast-local


phylum	13/13(100.0)
class	13/13(100.0)
order	10/13(76.92)
family	7/13(53.85)
genus	5/13(38.46)

Linkage accuracy:
phylum	39/39(100.0)
class	39/39(100.0)
order	35/39(89.74)
family	31/39(79.49)
genus	23/39(58.97)

'''