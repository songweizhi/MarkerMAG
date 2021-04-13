
def get_GapFilling_stats_by_assembly(free_living_16s_ref_file,
                                     free_living_ctg_ref_file,
                                     mini_assembly_to_16s_reads,
                                     mini_assembly_to_ctg_reads,
                                     ctg_level_min_link,
                                     mini_assembly_to_16s_ctg_connector,
                                     gnm_to_ctg_connector,
                                     marker_to_ctg_gnm_Key_connector,
                                     stats_GapFilling_ctg,
                                     stats_GapFilling_gnm):

    max_within_cate_diff_pct = 80
    max_between_cate_diff_pct = 20

    round2_free_living_16s_ref_dict = {}
    for free_living_read_16s in open(free_living_16s_ref_file):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_16s_id = free_living_read_16s_split[0]
            read_16s_refs = free_living_read_16s_split[1].split(',')
            round2_free_living_16s_ref_dict[read_16s_id] = read_16s_refs

    round2_free_living_ctg_ref_dict = {}
    for free_living_read_ctg in open(free_living_ctg_ref_file):
        free_living_read_ctg_split = free_living_read_ctg.strip().split('\t')
        read_ctg_id = free_living_read_ctg_split[0]
        read_ctg_refs = free_living_read_ctg_split[1].split(',')
        read_ctg_refs_no_suffix = []
        for each_read_ctg_ref in read_ctg_refs:
            if each_read_ctg_ref[-2:] in ['_l', '_r']:
                each_read_ctg_ref_no_suffix = each_read_ctg_ref[:-2]
                read_ctg_refs_no_suffix.append(each_read_ctg_ref_no_suffix)
        round2_free_living_ctg_ref_dict[read_ctg_id] = read_ctg_refs_no_suffix

    mini_assembly_to_16s_dict = {}
    for each_mini_assembly in open(mini_assembly_to_16s_reads):
        mini_assembly_split = each_mini_assembly.strip().split('\t')
        mini_assembly_id = mini_assembly_split[0]
        mini_assembly_mapped_reads = mini_assembly_split[1].split(',')
        for each_mapped_read in mini_assembly_mapped_reads:
            mapped_read_16s_refs = round2_free_living_16s_ref_dict.get(each_mapped_read, [])
            for each_mapped_read_16s_ref in mapped_read_16s_refs:
                mini_assembly_to_16s_key = '%s%s%s' % (mini_assembly_id, mini_assembly_to_16s_ctg_connector, each_mapped_read_16s_ref)
                if mini_assembly_to_16s_key not in mini_assembly_to_16s_dict:
                    mini_assembly_to_16s_dict[mini_assembly_to_16s_key] = 1
                else:
                    mini_assembly_to_16s_dict[mini_assembly_to_16s_key] += 1

    mini_assembly_to_ctg_dict = {}
    for each_mini_assembly in open(mini_assembly_to_ctg_reads):
        mini_assembly_split = each_mini_assembly.strip().split('\t')
        mini_assembly_id = mini_assembly_split[0]
        mini_assembly_mapped_reads = mini_assembly_split[1].split(',')
        for each_mapped_read in mini_assembly_mapped_reads:
            mapped_read_ctg_refs = round2_free_living_ctg_ref_dict.get(each_mapped_read, [])
            for each_mapped_read_ctg_ref in mapped_read_ctg_refs:
                mini_assembly_to_ctg_key = '%s%s%s' % (mini_assembly_id, mini_assembly_to_16s_ctg_connector, each_mapped_read_ctg_ref)
                if mini_assembly_to_ctg_key not in mini_assembly_to_ctg_dict:
                    mini_assembly_to_ctg_dict[mini_assembly_to_ctg_key] = 1
                else:
                    mini_assembly_to_ctg_dict[mini_assembly_to_ctg_key] += 1

    mini_assembly_to_16s_dict_reformatted = {}
    for each in mini_assembly_to_16s_dict:
        mini_assembly_id = each.split(mini_assembly_to_16s_ctg_connector)[0]
        seq_16s_id = each.split(mini_assembly_to_16s_ctg_connector)[1]
        linkage_num = mini_assembly_to_16s_dict[each]
        seq_16s_with_num = '%s__num__%s' % (seq_16s_id, linkage_num)
        if linkage_num >= ctg_level_min_link:
            if mini_assembly_id not in mini_assembly_to_16s_dict_reformatted:
                mini_assembly_to_16s_dict_reformatted[mini_assembly_id] = {seq_16s_with_num}
            else:
                mini_assembly_to_16s_dict_reformatted[mini_assembly_id].add(seq_16s_with_num)

    mini_assembly_to_ctg_dict_reformatted = {}
    for each in mini_assembly_to_ctg_dict:
        mini_assembly_id = each.split(mini_assembly_to_16s_ctg_connector)[0]
        ctg_id = each.split(mini_assembly_to_16s_ctg_connector)[1]
        linkage_num = mini_assembly_to_ctg_dict[each]
        ctg_with_num = '%s__num__%s' % (ctg_id, linkage_num)
        if linkage_num >= ctg_level_min_link:
            if mini_assembly_id not in mini_assembly_to_ctg_dict_reformatted:
                mini_assembly_to_ctg_dict_reformatted[mini_assembly_id] = {ctg_with_num}
            else:
                mini_assembly_to_ctg_dict_reformatted[mini_assembly_id].add(ctg_with_num)

    mini_assembly_linked_both = set(mini_assembly_to_16s_dict_reformatted).intersection(mini_assembly_to_ctg_dict_reformatted)

    stats_GapFilling_ctg_handle = open(stats_GapFilling_ctg, 'w')
    stats_GapFilling_gnm_dict = {}
    for each_mini_assembly in mini_assembly_linked_both:
        linked_16s = mini_assembly_to_16s_dict_reformatted[each_mini_assembly]
        linked_ctg = mini_assembly_to_ctg_dict_reformatted[each_mini_assembly]
        linked_16s_num_list = [int(i.split('__num__')[1]) for i in linked_16s]
        linked_ctg_num_list = [int(i.split('__num__')[1]) for i in linked_ctg]
        linked_16s_num_max = max(linked_16s_num_list)
        linked_ctg_num_max = max(linked_ctg_num_list)

        for each_linked_16s in linked_16s:
            linked_16s_id = each_linked_16s.split('__num__')[0]
            linked_16s_num = int(each_linked_16s.split('__num__')[1])
            for each_linked_ctg in linked_ctg:
                linked_ctg_id = each_linked_ctg.split('__num__')[0]
                linked_gnm_id = linked_ctg_id.split(gnm_to_ctg_connector)[0]
                linked_ctg_num = int(each_linked_ctg.split('__num__')[1])
                if linked_gnm_id in ['TR', 'PS']:
                    print('%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_gnm_id, linked_ctg_num))

        if (min(linked_16s_num_max, linked_ctg_num_max) * 100 / max(linked_16s_num_max, linked_ctg_num_max)) >= max_between_cate_diff_pct:

            linked_16s_filtered = [i for i in linked_16s if int(i.split('__num__')[1])*100/linked_16s_num_max >= max_within_cate_diff_pct]
            linked_ctg_filtered = [i for i in linked_ctg if int(i.split('__num__')[1])*100/linked_ctg_num_max >= max_within_cate_diff_pct]

            for each_linked_16s in linked_16s_filtered:
                linked_16s_id = each_linked_16s.split('__num__')[0]
                linked_16s_num = int(each_linked_16s.split('__num__')[1])
                for each_linked_ctg in linked_ctg_filtered:
                    linked_ctg_id = each_linked_ctg.split('__num__')[0]
                    linked_gnm_id = linked_ctg_id.split(gnm_to_ctg_connector)[0]
                    linked_ctg_num = int(each_linked_ctg.split('__num__')[1])
                    stats_GapFilling_ctg_handle.write('%s\t%s\t%s\n' % (linked_16s_id, linked_ctg_id, (linked_16s_num + linked_ctg_num)))
                    marker_to_gnm_key = '%s%s%s' % (linked_16s_id, marker_to_ctg_gnm_Key_connector, linked_gnm_id)
                    if marker_to_gnm_key not in stats_GapFilling_gnm_dict:
                        stats_GapFilling_gnm_dict[marker_to_gnm_key] = (linked_16s_num + linked_ctg_num)
                    else:
                        stats_GapFilling_gnm_dict[marker_to_gnm_key] += (linked_16s_num + linked_ctg_num)


    stats_GapFilling_ctg_handle.close()

    stats_GapFilling_gnm_handle = open(stats_GapFilling_gnm, 'w')
    stats_GapFilling_gnm_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_16s_to_gnm in stats_GapFilling_gnm_dict:
        each_16s_to_gnm_split = each_16s_to_gnm.split(marker_to_ctg_gnm_Key_connector)
        id_16s = each_16s_to_gnm_split[0]
        id_gnm = each_16s_to_gnm_split[1]
        linkage_num = stats_GapFilling_gnm_dict[each_16s_to_gnm]
        stats_GapFilling_gnm_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (id_16s, id_gnm, linkage_num))
    stats_GapFilling_gnm_handle.close()


step_2_wd                   = '/Users/songweizhi/Desktop/step_2_wd'
free_living_16s_ref_file    = '%s/round2_free_living_16s_refs.txt'  % step_2_wd
free_living_ctg_ref_file    = '%s/round2_free_living_ctg_refs.txt'  % step_2_wd
mini_assembly_to_16s_reads  = '%s/mini_assembly_to_16s_reads.txt'   % step_2_wd
mini_assembly_to_ctg_reads  = '%s/mini_assembly_to_ctg_reads.txt'   % step_2_wd
stats_GapFilling_ctg        = '%s/stats_GapFilling_ctg.txt'         % step_2_wd
stats_GapFilling_gnm        = '%s/stats_GapFilling_gnm.txt'         % step_2_wd
ctg_level_min_link = 3
mini_assembly_to_16s_ctg_connector = '___Mini___'
gnm_to_ctg_connector = '___C___'
marker_to_ctg_gnm_Key_connector = '___M___'


get_GapFilling_stats_by_assembly(free_living_16s_ref_file,
                                 free_living_ctg_ref_file,
                                 mini_assembly_to_16s_reads,
                                 mini_assembly_to_ctg_reads,
                                 ctg_level_min_link,
                                 mini_assembly_to_16s_ctg_connector,
                                 gnm_to_ctg_connector,
                                 marker_to_ctg_gnm_Key_connector,
                                 stats_GapFilling_ctg,
                                 stats_GapFilling_gnm)


'''
3_GI_subsample_10_82	Refined_19	2369	S2	Correct
3_GI_subsample_100_563	Refined_19	2224	S2	Correct
3_GI_subsample_50_207	Refined_19	2219	S2	Correct
3_GI_subsample_50_890	Refined_6	525	S2	Correct
3_GI_subsample_100_1450	Refined_53	496	S2	Correct
3_GI_subsample_100_1692	Refined_37	428	S2	Correct
3_GI_subsample_75_487	Refined_1	289	S2	Correct
3_GI_subsample_100_577	Refined_8	129	S2	Correct
3_GI_subsample_50_241	Refined_61	113	S2	Unknown
3_GI_subsample_100_579	Refined_8	113	S2	Correct
3_GI_subsample_100_1740	Refined_45	75	S2	Correct
3_GI_subsample_50_234	Refined_24	60	S2	Correct
3_GI_subsample_100_1404	Refined_38	56	S2	Correct
3_GI_subsample_50_647	Refined_38	50	S2	Correct
3_GI_subsample_100_1342	Refined_46	48	S2	Correct
3_GI_subsample_100_331	Refined_15	48	S2	Correct
3_GI_subsample_100_878	Refined_15	42	S2	Correct
3_GI_subsample_75_1022	Refined_28	41	S2	Wrong
3_GI_subsample_100_2288	Refined_89	28	S2	Correct
3_GI_subsample_75_2796	Refined_73	27	S2	Correct
3_GI_subsample_75_746	Refined_49	21	S2	Correct
3_GI_subsample_75_1549	Refined_26	13	S2	Correct
3_GI_subsample_100_457	Refined_54	13	S2	Correct
'''
