
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
    max_link_nun_dict_16s = {}
    for each in mini_assembly_to_16s_dict:
        mini_assembly_id = each.split(mini_assembly_to_16s_ctg_connector)[0]
        seq_16s_id = each.split(mini_assembly_to_16s_ctg_connector)[1]
        linkage_num = mini_assembly_to_16s_dict[each]
        seq_16s_with_num = '%s__num__%s' % (seq_16s_id, linkage_num)
        if linkage_num >= ctg_level_min_link:

            # add to mini_assembly_to_16s_dict_reformatted
            if mini_assembly_id not in mini_assembly_to_16s_dict_reformatted:
                mini_assembly_to_16s_dict_reformatted[mini_assembly_id] = {seq_16s_with_num}
            else:
                mini_assembly_to_16s_dict_reformatted[mini_assembly_id].add(seq_16s_with_num)

            # add to max_link_nun_dict_16s
            if seq_16s_id not in max_link_nun_dict_16s:
                max_link_nun_dict_16s[seq_16s_id] = linkage_num
            else:
                if linkage_num > max_link_nun_dict_16s[seq_16s_id]:
                    max_link_nun_dict_16s[seq_16s_id] = linkage_num

    mini_assembly_to_ctg_dict_reformatted = {}
    max_link_nun_dict_ctg = {}
    for each in mini_assembly_to_ctg_dict:
        mini_assembly_id = each.split(mini_assembly_to_16s_ctg_connector)[0]
        ctg_id = each.split(mini_assembly_to_16s_ctg_connector)[1]
        linkage_num = mini_assembly_to_ctg_dict[each]
        ctg_with_num = '%s__num__%s' % (ctg_id, linkage_num)
        if linkage_num >= ctg_level_min_link:

            # add to mini_assembly_to_ctg_dict_reformatted
            if mini_assembly_id not in mini_assembly_to_ctg_dict_reformatted:
                mini_assembly_to_ctg_dict_reformatted[mini_assembly_id] = {ctg_with_num}
            else:
                mini_assembly_to_ctg_dict_reformatted[mini_assembly_id].add(ctg_with_num)

            # add to max_link_nun_dict_ctg
            if ctg_id not in max_link_nun_dict_ctg:
                max_link_nun_dict_ctg[ctg_id] = linkage_num
            else:
                if linkage_num > max_link_nun_dict_ctg[ctg_id]:
                    max_link_nun_dict_ctg[ctg_id] = linkage_num

    mini_assembly_linked_both = set(mini_assembly_to_16s_dict_reformatted).intersection(mini_assembly_to_ctg_dict_reformatted)

    print(max_link_nun_dict_16s['3_GI_subsample_75_1022'])
    print(max_link_nun_dict_ctg['Refined_28___C___NODE_4309_length_14678_cov_0.076146'])


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
                if linked_gnm_id in ['Refined_28']:
                    #print('%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_gnm_id, linked_ctg_num))
                    pass
                if linked_16s_id in ['3_GI_subsample_75_1022']:
                    #print('%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_ctg_id, linked_ctg_num))
                    pass

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

                    linked_16s_max_num = max_link_nun_dict_16s[linked_16s_id]
                    linked_ctg_max_num = max_link_nun_dict_ctg[linked_ctg_id]
                    # print('linked_16s_max_num: %s' % linked_16s_max_num)
                    # print('linked_ctg_max_num: %s' % linked_ctg_max_num)
                    # if linked_gnm_id in ['Refined_28']:
                    #     print('%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_gnm_id, linked_ctg_num))
                    #
                    # if linked_16s_id in ['3_GI_subsample_75_1022']:
                    #     print('%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_ctg_id, linked_ctg_num))
                    linked_16s_num_pct_by_max = linked_16s_num*100/linked_16s_max_num
                    linked_ctg_num_pct_by_max = linked_ctg_num*100/linked_ctg_max_num

                    if (linked_16s_num_pct_by_max >= 50) and (linked_ctg_num_pct_by_max >= 50):
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



step_2_wd                   = '/Users/songweizhi/Desktop/step_2_wd_GI'
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

NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1022(64)	Refined_28(10)
NODE_375_length_663_cov_0.011194	3_GI_subsample_75_1022(25)	Refined_65(3)
NODE_375_length_663_cov_0.011194	3_GI_subsample_75_1022(25)	Refined_28(16)



NODE_375_length_663_cov_0.011194	3_GI_subsample_100_991(16)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_100_1043(16)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_50_761(7)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_50_741(7)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_75_1(17)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_50_743(5)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_75_1267(5)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_75_1022(25)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_50_708(3)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_100_1392(19)	Refined_28(16)
NODE_375_length_663_cov_0.011194	3_GI_subsample_75_2092(7)	Refined_28(16)
NODE_328_length_709_cov_0.273196	3_GI_subsample_100_578(4)	Refined_28(6)
NODE_328_length_709_cov_0.273196	3_GI_subsample_100_577(4)	Refined_28(6)
NODE_328_length_709_cov_0.273196	3_GI_subsample_100_1755(4)	Refined_28(6)
NODE_328_length_709_cov_0.273196	3_GI_subsample_100_1450(322)	Refined_28(6)
NODE_328_length_709_cov_0.273196	3_GI_subsample_50_672(122)	Refined_28(6)
NODE_328_length_709_cov_0.273196	3_GI_subsample_100_579(4)	Refined_28(6)
NODE_667_length_449_cov_0.080745	3_GI_subsample_100_1740(31)	Refined_28(5)
NODE_667_length_449_cov_0.080745	3_GI_subsample_100_2291(21)	Refined_28(5)
NODE_667_length_449_cov_0.080745	3_GI_subsample_50_739(3)	Refined_28(5)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_743(152)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1538(12)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_741(299)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1517(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_725(34)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_10_219(23)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_711(56)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_2065(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1129(16)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_742(55)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1554(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_436(19)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1587(21)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_735(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1658(151)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_692(6)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_649(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_424(43)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1519(6)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_757(50)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_727(249)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1170(22)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_761(283)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_100_1491(38)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_695(42)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_692(151)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1267(247)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1181(36)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_693(4)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_694(38)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_734(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_724(121)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_434(54)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_447(20)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_788(18)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_738(25)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_449(13)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_10_211(4)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_739(264)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1009(9)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_438(6)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_740(73)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_709(3)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_729(37)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_435(8)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_728(64)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_710(38)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_444(13)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1143(21)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1111(35)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_25_443(12)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_726(45)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1344(17)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1203(11)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_708(22)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_699(5)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_50_730(58)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_5_117(35)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1295(25)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_1022(64)	Refined_28(10)
NODE_3887_length_301_cov_0.597701	3_GI_subsample_75_2092(8)	Refined_28(10)

'''
