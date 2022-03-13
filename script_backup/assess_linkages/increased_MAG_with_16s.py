
def get_increased_mag_with_16s(barrnap_op_txt, len_16s_cutoff, markermag_op):

    mag_with_16S_set = set()
    for each in open(barrnap_op_txt):
        each_split = each.strip().split('\t')
        if not each.startswith('Genome\t16S'):
            if len(each_split) > 1:
                mag_id = each_split[0]
                len_16s = int(each_split[2])
                if len_16s >= len_16s_cutoff:
                    mag_with_16S_set.add(mag_id)

    all_linked_gnm_set= set()
    increased_mag_with_16s_set = set()
    for each_linkage in open(markermag_op):
        each_linkage_split = each_linkage.strip().split('\t')
        if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Round'):
            gnm_id= each_linkage_split[1]
            all_linked_gnm_set.add(gnm_id)
            if gnm_id not in mag_with_16S_set:
                increased_mag_with_16s_set.add(gnm_id)


    return mag_with_16S_set, all_linked_gnm_set, increased_mag_with_16s_set


wd                  = '/Users/songweizhi/Desktop/increased_MAGs_with_16S'
len_16s_cutoff      = 1200

# kelp
barrnap_op_txt      = '%s/BH_ER_050417_refined_bins_complete50.0_contain5.0_16S.txt'    % wd
markermag_op        = '%s/Kelp_0726_45_45_1200_iden99_mismatch2_linkages_by_genome.txt' % wd

# oil field
barrnap_op_txt      = '%s/OF_16S.txt'                                                   % wd
markermag_op        = '%s/OF_0726_45_45_1200_iden99_mismatch2_linkages_by_genome.txt'   % wd

# pig gut
barrnap_op_txt      = '%s/pig_322_MAGs_16S.txt'                                         % wd
markermag_op        = '%s/Pig_0726_45_45_1200_iden99_mismatch2_linkages_by_genome.txt'  % wd

mags_with_16s, all_linked_gnm, increased_mag_with_16s = get_increased_mag_with_16s(barrnap_op_txt, len_16s_cutoff, markermag_op)
print('Number of MAGs with 16S longer than %s bp:\t%s'  % (len_16s_cutoff, len(mags_with_16s)))
print('Number of MAGs linked by MarkerMAG:\t%s'         % len(all_linked_gnm))
print('Increased number of MAGs with 16S:\t%s'          % len(increased_mag_with_16s))
