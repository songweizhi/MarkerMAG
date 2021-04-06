
final_op = ''
final_op_ctg_level = ''


# summarize linkages at contig level
mock_final_op_ctg_level_handle = open(final_op_ctg_level, 'w')
mock_final_op_ctg_level_handle.write('Marker___Genome(total)\tContig\tPaired\tClipping\tStep\n')
for each_linkage in open(final_op):
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
        each_linkage_split = each_linkage.strip().split('\t')
        marker_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        total_link_num = int(each_linkage_split[2])
        link_step = each_linkage_split[3]

        # first go through link num dict by paired reads
        counted_16s_to_ctg_key = set()
        for each_paired_link in marker_to_ctg_link_num_dict_pair:
            paired_link_16s_id = each_paired_link.split(marker_to_ctg_Key_connector_str)[0]
            paired_link_ctg_id = each_paired_link.split(marker_to_ctg_Key_connector_str)[1]
            paired_link_ctg_id_no_gnm = paired_link_ctg_id.split(gnm_ctg_connector)[1]
            paired_link_gnm_id = paired_link_ctg_id.split(gnm_ctg_connector)[0]
            if (paired_link_16s_id == marker_id) and (paired_link_gnm_id == mag_id):
                current_pair_link_num = marker_to_ctg_link_num_dict_pair[each_paired_link]
                current_clip_link_num = marker_to_ctg_link_num_dict_clip.get(each_paired_link, 0)
                mock_final_op_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t%s\t%s\n' % (paired_link_16s_id, paired_link_gnm_id, total_link_num, paired_link_ctg_id_no_gnm, current_pair_link_num, current_clip_link_num, link_step))
                counted_16s_to_ctg_key.add(each_paired_link)

        # then go through link num dict by clipping reads
        for each_clip_link in marker_to_ctg_link_num_dict_clip:
            if each_clip_link not in counted_16s_to_ctg_key:
                clip_link_16s_id = each_clip_link.split(marker_to_ctg_Key_connector_str)[0]
                clip_link_ctg_id = each_clip_link.split(marker_to_ctg_Key_connector_str)[1]
                clip_link_ctg_id_no_gnm = clip_link_ctg_id.split(gnm_ctg_connector)[1]
                clip_link_gnm_id = clip_link_ctg_id.split(gnm_ctg_connector)[0]
                if (clip_link_16s_id == marker_id) and (clip_link_gnm_id == mag_id):
                    current_pair_link_num = marker_to_ctg_link_num_dict_pair.get(each_clip_link, 0)
                    current_clip_link_num = marker_to_ctg_link_num_dict_clip[each_clip_link]
                    mock_final_op_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t%s\t%s\n' % (clip_link_16s_id, clip_link_gnm_id, total_link_num, clip_link_ctg_id_no_gnm, current_pair_link_num, current_clip_link_num, link_step))
                    counted_16s_to_ctg_key.add(each_clip_link)
mock_final_op_ctg_level_handle.close()


