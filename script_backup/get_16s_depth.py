#
# # calculate depth for input 16S
# mean_depth_dict_16s = {}
# if min_16s_gnm_multiple > 0:
#
#     if depth_file_16s is None:
#
#         print('Round 1: calculating depth for 16S')
#         # get mean depth for 16S sequences
#         mean_depth_dict_16s, s16_len_dict = get_ctg_mean_depth_by_samtools_coverage_global(True, input_16s_qc,
#                                                                                            reads_file_r1_subset,
#                                                                                            reads_file_r2_subset, '',
#                                                                                            subsample_rate_for_depth_estimation,
#                                                                                            num_threads)
#
#         # write out 16s depth
#         depth_file_16s_handle = open(depth_file_16s_calculated, 'w')
#         for s16 in mean_depth_dict_16s:
#             depth_file_16s_handle.write('%s\t%s\n' % (s16, mean_depth_dict_16s[s16]))
#         depth_file_16s_handle.close()
#     else:
#         print('Round 1: read in provided 16S depth')
#         # read in depth and store in mean_depth_dict_gnm here
#         for each_16s_depth in open(depth_file_16s):
#             each_16s_depth_split = each_16s_depth.strip().split('\t')
#             s16_id = each_16s_depth_split[0]
#             s16_depth = float(each_16s_depth_split[1])
#             mean_depth_dict_16s[s16_id] = s16_depth