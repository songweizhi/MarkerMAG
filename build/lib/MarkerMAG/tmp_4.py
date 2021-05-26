#
# linkage_file       = '/Users/songweizhi/Desktop/assess_linkages_hc/hc_0512_tmp_stats_combined_filtered_with_assessment_ani97_imag99.5_i16S99.5.txt'
# clp_pct_ratio_file = '/Users/songweizhi/Desktop/assess_linkages_hc/rd1_clp_pct_diff.txt'
#
# filtered_marker_to_gnm_set = set()
# for each_link in open(linkage_file):
#     each_link_split = each_link.strip().split('\t')
#     marker_id = each_link_split[0]
#     genome_id = each_link_split[1]
#     marker_to_gnm_key = '%s___%s' % (marker_id, genome_id)
#     filtered_marker_to_gnm_set.add(marker_to_gnm_key)
#
#
# for each_link in open(clp_pct_ratio_file):
#     each_link_split = each_link.strip().split('\t')
#     marker_id = each_link_split[0]
#     genome_id = each_link_split[1]
#     marker_to_gnm_key = '%s___%s' % (marker_id, genome_id)
#     if marker_to_gnm_key in filtered_marker_to_gnm_set:
#         if genome_id in ['cami_hc_3', 'cami_hc_13', 'cami_hc_14', 'cami_hc_101', 'cami_hc_38', 'cami_hc_127', 'cami_hc_11']:
#             print(each_link.strip())
#             pass
#         else:
#             pass
#


reads_num = 315957101
subsample_rate = 10
reads_to_subsample = round(reads_num/subsample_rate)
print(reads_num/subsample_rate)
print(reads_to_subsample)