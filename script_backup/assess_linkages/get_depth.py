
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


wd                      = '/Users/songweizhi/Desktop/get_depth'
len_16s_cutoff          = 1200

# MBARC26
barrnap_op_txt          = '%s/MBARC26_refined_bins_50_5_16S.txt'                                            % wd
markermag_op            = '%s/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_linkages_by_genome.txt'    % wd
good_quality_MAG_txt    = '%s/MBARC26_refined_bins_50_5.txt'                                                % wd
mag_depth_txt           = '%s/MBARC26_refined_bins_50_5_depth.txt'                                          % wd
MAG_with_matam_16s      = '%s/MBARC26_MAG_with_matam_16s.txt'                                               % wd
MBARC26 = True


# # GI
# barrnap_op_txt          = '%s/GI_refined_bins_complete50_contain5_16S.txt'    % wd
# markermag_op            = '%s/GI_0726_128_45_45_1200_diff80_iden99_mismatch2_GoodMAGs_linkages_by_genome.txt' % wd
# good_quality_MAG_txt    = '%s/GI_refined_bins_complete50_contain5.txt'                      % wd
# mag_depth_txt           = '%s/GI_bin_mean_depth.txt'                                        % wd
# MAG_with_matam_16s      = '%s/GI_MAG_with_matam_16s.txt' % wd

# # Oral
# barrnap_op_txt          = '%s/Oral_refined_MAGs_complete50_contain5_16S.txt'    % wd
# markermag_op            = '%s/Oral_0729_45_45_min1200_mismatch2_diff80_iden99_GoodMAGs_linkages_by_genome.txt' % wd
# good_quality_MAG_txt    = '%s/Oral_refined_MAGs_complete50_contain5.txt'                      % wd
# mag_depth_txt           = '%s/Oral_MAGs_depth.txt'                                        % wd
# MAG_with_matam_16s      = '%s/Oral_MAG_with_matam_16s.txt' % wd

# # Kelp
# barrnap_op_txt          = '%s/Kelp_BH_ER_050417_refined_bins_complete50.0_contain5.0_16S.txt'    % wd
# markermag_op            = '%s/Kelp_0726_45_45_1200_iden99_mismatch2_linkages_by_genome.txt' % wd
# good_quality_MAG_txt    = '%s/Kelp_BH_ER_050417_refined_bins_complete50_contain5.txt'                      % wd
# mag_depth_txt           = '%s/Kelp_BH_ER_050417_refined_bins_mean_depth.txt'                                        % wd
#
# # Oil
# barrnap_op_txt          = '%s/OF_16S.txt'    % wd
# markermag_op            = '%s/OF_0726_45_45_1200_iden99_mismatch2_linkages_by_genome.txt' % wd
# good_quality_MAG_txt    = '%s/OF_42_MAGs.txt'                      % wd
# mag_depth_txt           = '%s/OF_42_MAGs_mean_depth.txt'                                        % wd
#
# # Pig
# barrnap_op_txt          = '%s/pig_322_MAGs_16S.txt'    % wd
# markermag_op            = '%s/Pig_0726_45_45_1200_iden99_mismatch2_linkages_by_genome.txt' % wd
# good_quality_MAG_txt    = '%s/pig_322_MAGs.txt'                      % wd
# mag_depth_txt           = '%s/pig_322_MAGs_mean_depth.txt'                                        % wd


MAG_with_matam_16s_set = set()
for each_mag in open(MAG_with_matam_16s):
    MAG_with_matam_16s_set.add(each_mag.strip())

good_quality_MAG_set = set()
for each_mag in open(good_quality_MAG_txt):
    good_quality_MAG_set.add(each_mag.strip())

mag_depth_dict = {}
for each_mag_depth in open(mag_depth_txt):
    each_mag_depth_split = each_mag_depth.strip().split('\t')
    if each_mag_depth_split[0] in good_quality_MAG_set:
        mag_depth_dict[each_mag_depth_split[0]] = float(each_mag_depth_split[2])

mags_with_16s, all_linked_gnm, increased_mag_with_16s = get_increased_mag_with_16s(barrnap_op_txt, len_16s_cutoff, markermag_op)

print('MAG\tDepth\tMatam_16S\tLinked')
for each_mag in good_quality_MAG_set:
    mag_depth = mag_depth_dict.get(each_mag, 'NA')

    linked_to_16s = 'no'
    if each_mag in all_linked_gnm:
        linked_to_16s = 'yes'

    with_matam_16s = 'no'
    if each_mag in MAG_with_matam_16s_set:
        with_matam_16s = 'yes'

    print('%s\t%s\t%s\t%s' % (each_mag, mag_depth, with_matam_16s, linked_to_16s))

