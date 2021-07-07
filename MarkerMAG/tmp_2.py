
stats_GapFilling_ctg = '/Users/songweizhi/Desktop/stats_GapFilling_ctg.txt'
marker_to_ctg_gnm_Key_connector = '------'
free_living_16s_to_ctg_linkage_dict_to_use = {}
for each_ctg_level_link in open(stats_GapFilling_ctg):
    if ',Number\n' not in each_ctg_level_link:
        each_ctg_level_link_split = each_ctg_level_link.split(',')
        marker_id = each_ctg_level_link_split[0]
        ctg_id = each_ctg_level_link_split[1]
        link_num = int(each_ctg_level_link_split[2])
        current_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, ctg_id)
        free_living_16s_to_ctg_linkage_dict_to_use[current_key] = link_num
