
wd      = '/Users/songweizhi/Desktop/subset_hc_reads'
wd      = '.'
sample  = 'S001'

reads_mapping_tsv   = '%s/H_%s__insert_180_gs_read_mapping.tsv'     % (wd, sample)
refs_with_16s_txt   = '%s/ref_genomes_hc_with_16S.txt'              % wd
reads_to_keep_file  = '%s/reads_to_keep_%s.txt'                     % (wd, sample)


refs_with_16s_set = set()
for each_ref in open(refs_with_16s_txt):
    refs_with_16s_set.add(each_ref.strip()[4:])


reads_to_keep_file_handle = open(reads_to_keep_file, 'w')
for each_read in open(reads_mapping_tsv):
    each_read_split = each_read.strip().split('\t')
    read_id = each_read_split[0]
    tax_id = each_read_split[2]
    if tax_id in refs_with_16s_set:
        reads_to_keep_file_handle.write('%s/1\n' % read_id)
        reads_to_keep_file_handle.write('%s/2\n' % read_id)
reads_to_keep_file_handle.close()

