import os
import glob


wd = '/Users/songweizhi/Desktop/subset_hc_redas'
reads_mapping_tsv = '%s/mapping_subset.tsv' % wd
refs_with_16s_txt = '%s/ref_genomes_hc_with_16S.txt'    % wd


refs_with_16s_set = set()
for each_ref in open(refs_with_16s_txt):
    refs_with_16s_set.add(each_ref.strip()[4:])

n = 0
for each_read in open(reads_mapping_tsv):
    each_read_split = each_read.strip().split('\t')
    print(each_read_split)
    tax_id = each_read_split[2]

    if tax_id in refs_with_16s_set:
        print(each_read_split)
        n += 1
print(n)



