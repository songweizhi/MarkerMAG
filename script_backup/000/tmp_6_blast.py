

def get_mag_gc_bias_worker(argument_list):

    mag                                 = argument_list[0]
    current_mag_ctg_set                 = argument_list[1]
    current_mag_pos_depth_dict          = argument_list[2]
    ctg_seq_dict                        = argument_list[3]
    ctg_len_dict                        = argument_list[4]
    ctg_len_cutoff                      = argument_list[5]
    current_mag_region_to_ignore_dict   = argument_list[6]
    ignore_ends_len_mag                 = argument_list[7]
    window_len                          = argument_list[8]
    mag_gc_content_dict                 = argument_list[9]
    gnm_to_linked_16s_gc_content_dict   = argument_list[10]
    current_mag_gc_bias_txt             = argument_list[11]
    current_mag_gc_bias_png             = argument_list[12]
    current_mag_depth_GC_content_file   = argument_list[13]

    print('Processing %s' % mag)
    current_mag_total_mapped_read_len = 0
    current_mag_total_counted_ctg_len = 0
    current_mag_gc_content_to_depth_list_dict = {}
    for each_ctg in current_mag_ctg_set:
        current_ctg_pos_depth_dict = current_mag_pos_depth_dict.get(each_ctg, {})
        current_ctg_len = ctg_len_dict.get(each_ctg, 0)
        if (len(current_ctg_pos_depth_dict) > 0) and (current_ctg_len >= ctg_len_cutoff):
            current_ctg_seq = ctg_seq_dict[each_ctg]
            current_ctg_regions_to_ignore = current_mag_region_to_ignore_dict.get(each_ctg, set())
            print('current_ctg_regions_to_ignore')
            print(current_ctg_regions_to_ignore)
            window_start_pos = ignore_ends_len_mag + 1
            while window_start_pos <= (current_ctg_len - (ignore_ends_len_mag + window_len) + 1):
                window_end_pos = window_start_pos + window_len - 1

                if (str(window_start_pos) not in current_ctg_regions_to_ignore) and (
                        str(window_end_pos) not in current_ctg_regions_to_ignore):
                    if (str(window_start_pos) in current_ctg_pos_depth_dict) and (
                            str(window_end_pos) in current_ctg_pos_depth_dict):

                        # get sequence of current window
                        if window_start_pos < (len(current_ctg_seq) - window_len + 1):
                            window_seq = current_ctg_seq[(window_start_pos - 1):window_end_pos]
                        else:
                            window_seq = current_ctg_seq[(window_start_pos - 1):]
                        window_seq_upper = window_seq.upper()

                        masked_base_num = window_seq_upper.count('N')
                        if masked_base_num <= 4:
                            window_pos_list = list(range(window_start_pos, (window_end_pos + 1)))
                            window_pos_depth_list = [current_ctg_pos_depth_dict[str(pos)] for pos in window_pos_list]
                            window_mean_depth = sum(window_pos_depth_list) / len(window_pos_depth_list)
                            window_gc_content = (window_seq_upper.count('G') + window_seq_upper.count('C')) * 100 / len(
                                window_seq)
                            # print('%s-%s\t%sx\t%s%s' % (window_start_pos, window_end_pos, window_mean_depth, window_gc_content, '%'))

                            current_mag_total_mapped_read_len += window_mean_depth
                            current_mag_total_counted_ctg_len += 1

                            # when the length of sliding window is 100 bp
                            window_gc_content = int(window_gc_content)

                            if window_gc_content not in current_mag_gc_content_to_depth_list_dict:
                                current_mag_gc_content_to_depth_list_dict[window_gc_content] = [window_mean_depth]
                            else:
                                current_mag_gc_content_to_depth_list_dict[window_gc_content].append(window_mean_depth)
                window_start_pos += 1

    # get current_mag_average_global_coverage
    current_mag_average_global_coverage = current_mag_total_mapped_read_len / current_mag_total_counted_ctg_len
    current_mag_GG_content = mag_gc_content_dict.get(mag, 'NA')

    # write out current_mag_depth_GC_content_file
    current_mag_depth_GC_content_file_handle = open(current_mag_depth_GC_content_file, 'w')
    # current_mag_depth_GC_content_file_handle.write('MAG\tDepth(x)\tGC(%)\n')
    current_mag_depth_GC_content_file_handle.write('%s\t%s\t%s\n' % (mag, float("{0:.2f}".format(current_mag_average_global_coverage)), current_mag_GG_content))
    current_mag_depth_GC_content_file_handle.close()

    # write out current_mag_gc_bias_txt
    current_mag_gc_bias_txt_handle = open(current_mag_gc_bias_txt, 'w')
    current_mag_gc_bias_txt_handle.write('GC\tCoverage\tNormalized_coverage\tNumber\tProportion\n')
    for each_gc_group in sorted([i for i in current_mag_gc_content_to_depth_list_dict]):
        current_gc_group_depth_list = current_mag_gc_content_to_depth_list_dict[each_gc_group]
        current_gc_group_proportion = len(current_gc_group_depth_list) * 100 / current_mag_total_counted_ctg_len
        current_gc_group_mean_depth = sum(current_gc_group_depth_list) / len(current_gc_group_depth_list)
        current_gc_group_mean_depth_normalized = current_gc_group_mean_depth / current_mag_average_global_coverage
        current_mag_gc_bias_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_gc_group,
                                                                       float("{0:.2f}".format(
                                                                           current_gc_group_mean_depth)),
                                                                       float("{0:.2f}".format(
                                                                           current_gc_group_mean_depth_normalized)),
                                                                       len(current_gc_group_depth_list),
                                                                       float("{0:.4f}".format(
                                                                           current_gc_group_proportion))))
    current_mag_gc_bias_txt_handle.close()

    # plot GC bias
    linked_16s_gc_content_list = gnm_to_linked_16s_gc_content_dict[mag]
    plot_gc_bias(current_mag_gc_bias_txt, current_mag_GG_content, linked_16s_gc_content_list,
                 current_mag_average_global_coverage, current_mag_gc_bias_png)

