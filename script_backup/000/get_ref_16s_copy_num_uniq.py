
ref_gnm_renamed_16S_id_txt  = '/Users/songweizhi/Desktop/666/ref_gnm_renamed_16S_id.txt'
depth_file                  = '/Users/songweizhi/Desktop/666/000/uniq/reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted_depth.txt'
top_num                     = 13


ref_gnm_to_renamed_16S_id_dict = {}
for each_line in open(ref_gnm_renamed_16S_id_txt):
    gnm_id = each_line.strip().split('_')[0]
    if gnm_id not in ref_gnm_to_renamed_16S_id_dict:
        ref_gnm_to_renamed_16S_id_dict[gnm_id] = {each_line.strip()}
    else:
        ref_gnm_to_renamed_16S_id_dict[gnm_id].add(each_line.strip())


ref_to_pos_depth_dict = {}
for eac_pos in open(depth_file):
    eac_pos_split = eac_pos.strip().split('\t')
    ref_id = eac_pos_split[0]
    bp_depth = int(eac_pos_split[2])
    if ref_id not in ref_to_pos_depth_dict:
        ref_to_pos_depth_dict[ref_id] = [bp_depth]
    else:
        ref_to_pos_depth_dict[ref_id].append(bp_depth)


ref_gnm_depth_list_dict = {}
for each_ref in ref_to_pos_depth_dict:
    ref_gnm = each_ref.split('_')[0]
    pos_depth_list = ref_to_pos_depth_dict[each_ref]
    pos_depth_list_sorted = sorted(pos_depth_list)[::-1]
    top_depth_list = pos_depth_list_sorted[:top_num]
    mean_top_depth = sum(top_depth_list)/len(top_depth_list)
    mean_top_depth = float("{0:.2f}".format(mean_top_depth))
    if ref_gnm not in ref_gnm_depth_list_dict:
        ref_gnm_depth_list_dict[ref_gnm] = [mean_top_depth]
    else:
        ref_gnm_depth_list_dict[ref_gnm].append(mean_top_depth)


ref_gnm_mean_depth_dict = {}
for each_ref_gnm in ref_gnm_depth_list_dict:
    ref_gnm_depth_list = ref_gnm_depth_list_dict[each_ref_gnm]
    ref_gnm_mean_depth = sum(ref_gnm_depth_list)/len(ref_gnm_depth_list)
    ref_gnm_mean_depth = float("{0:.2f}".format(ref_gnm_mean_depth))
    ref_gnm_mean_depth_dict[each_ref_gnm] = ref_gnm_mean_depth

print(ref_gnm_mean_depth_dict)
