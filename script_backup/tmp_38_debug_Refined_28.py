
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


    ######################### debug report #########################

    check_list_gnm  = []
    check_list_16s  = []

    check_list_gnm = ['Refined_51', 'Refined_91']
    check_list_16s = ['3_GI_subsample_50_882', '3_GI_subsample_75_963']

    report_all      = False

    '''
    # 3_GI_subsample_75_963
    NODE_4066_length_298_cov_0.076023	3_GI_subsample_75_963(20)	Refined_28___C___NODE_4309_length_14678_cov_0.076146(6)
    NODE_660_length_449_cov_0.232919	3_GI_subsample_75_963(8)	Refined_84___C___NODE_14371_length_3638_cov_0.021646(3)
    NODE_660_length_449_cov_0.232919	3_GI_subsample_75_963(8)	Refined_91___C___NODE_3602_length_17682_cov_0.132270(6)
    NODE_12_length_1555_cov_0.052521	3_GI_subsample_75_963(160)	Refined_51___C___NODE_7061_length_8073_cov_0.110244(51)
    NODE_12_length_1555_cov_0.052521	3_GI_subsample_75_963(160)	Refined_28___C___NODE_4309_length_14678_cov_0.076146(96)
    
    # Refined_28
 
    '''

    ################################################################

    # for report
    error_report_linked_gnm_list = []
    error_report_linked_16s_list = []
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
                error_report_str_linked_gnm = '%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_gnm_id, linked_ctg_num)
                error_report_str_linked_16s = '%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_ctg_id, linked_ctg_num)

                if report_all is True:
                    if linked_gnm_id in check_list_gnm:
                        error_report_linked_gnm_list.append(error_report_str_linked_gnm)
                    if linked_16s_id in check_list_16s:
                        error_report_linked_16s_list.append(error_report_str_linked_16s)


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
                    error_report_str_linked_gnm = '%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_gnm_id, linked_ctg_num)
                    error_report_str_linked_16s = '%s\t%s(%s)\t%s(%s)' % (each_mini_assembly, linked_16s_id, linked_16s_num, linked_ctg_id, linked_ctg_num)

                    if report_all is False:
                        if linked_gnm_id in check_list_gnm:
                            error_report_linked_gnm_list.append(error_report_str_linked_gnm)
                        if linked_16s_id in check_list_16s:
                            error_report_linked_16s_list.append(error_report_str_linked_16s)

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

    ######################### debug report #########################


    if check_list_gnm != []:
        print('Mini-assembly linked to:\t%s' % ', '.join(check_list_gnm))
        for each_link_gnm in error_report_linked_gnm_list:
            print(each_link_gnm)
        print()

    if check_list_16s != []:
        print('Mini-assembly linked to:\t%s' % ', '.join(check_list_16s))
        for each_link_16s in error_report_linked_16s_list:
            print(each_link_16s)
        print()

    ################################################################

    stats_GapFilling_gnm_handle = open(stats_GapFilling_gnm, 'w')
    stats_GapFilling_gnm_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_16s_to_gnm in stats_GapFilling_gnm_dict:
        each_16s_to_gnm_split = each_16s_to_gnm.split(marker_to_ctg_gnm_Key_connector)
        id_16s = each_16s_to_gnm_split[0]
        id_gnm = each_16s_to_gnm_split[1]
        linkage_num = stats_GapFilling_gnm_dict[each_16s_to_gnm]
        stats_GapFilling_gnm_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (id_16s, id_gnm, linkage_num))
    stats_GapFilling_gnm_handle.close()

'''
14
Mini-assembly linked to:	Refined_28
NODE_419_length_586_cov_0.065359	3_GI_subsample_25_431(14)	Refined_28(45)
NODE_4801_length_285_cov_0.082278	3_GI_subsample_50_581(26)	Refined_28(6)
NODE_4801_length_285_cov_0.082278	3_GI_subsample_75_961(29)	Refined_28(6)
NODE_2373_length_336_cov_0.062201	3_GI_subsample_50_782(11)	Refined_28(5)
NODE_2373_length_336_cov_0.062201	3_GI_subsample_100_1646(13)	Refined_28(5)

15
Mini-assembly linked to:	Refined_28
NODE_4066_length_298_cov_0.076023	3_GI_subsample_75_961(29)	Refined_28(6)
NODE_4066_length_298_cov_0.076023	3_GI_subsample_50_581(26)	Refined_28(6)
NODE_12_length_1555_cov_0.052521	3_GI_subsample_75_963(160)	Refined_28(96)
NODE_2978_length_321_cov_0.000000	3_GI_subsample_100_1646(13)	Refined_28(5)
NODE_2978_length_321_cov_0.000000	3_GI_subsample_50_782(11)	Refined_28(5)
NODE_417_length_586_cov_0.089325	3_GI_subsample_25_431(14)	Refined_28(45)

'''


step_2_wd                   = '/Users/songweizhi/Desktop/51_91'
free_living_16s_ref_file    = '%s/round2_free_living_16s_refs.txt'      % step_2_wd
free_living_ctg_ref_file    = '%s/round2_free_living_ctg_refs.txt'      % step_2_wd
mini_assembly_to_16s_reads  = '%s/mini_assembly_to_16s_reads.txt'       % step_2_wd
mini_assembly_to_ctg_reads  = '%s/mini_assembly_to_ctg_reads.txt'       % step_2_wd
stats_GapFilling_ctg        = '%s/stats_GapFilling_ctg.txt'             % step_2_wd
stats_GapFilling_gnm        = '%s/stats_GapFilling_gnm.txt'             % step_2_wd
ctg_level_min_link          = 3

mini_assembly_to_16s_ctg_connector  = '___Mini___'
gnm_to_ctg_connector                = '___C___'
marker_to_ctg_gnm_Key_connector     = '___M___'


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

3_GI_subsample_75_963(169)  --- NODE_2265_length_348_cov_0.257919   ---     Refined_51___C___NODE_7061_length_8073_cov_0.110244(44)
3_GI_subsample_75_963(125)  --- NODE_126_length_934_cov_0.154895    ---     Refined_51___C___NODE_7061_length_8073_cov_0.110244(48)

3_GI_subsample_50_882(9)	--- NODE_3017_length_328_cov_0.144279   ---     Refined_91___C___NODE_3602_length_17682_cov_0.132270(33)



'''









'''
contig id:  NODE_11_length_1555_cov_0.052521
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0418_mis2_no_clp_MarkerMAG_wd/GI_0418_mis2_no_clp_step_2_wd
BioSAK select_seq -seq mini_assembly_SPAdes_wd/scaffolds.fasta -id NODE_11_length_1555_cov_0.052521.txt -out NODE_11_length_1555_cov_0.052521.fa -option 1
BioSAK select_seq -seq /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/3_GI_Matam16S_wd/3_GI_assembled_16S_uclust_0.999.fasta -id 3_GI_subsample_75_963.txt -out 3_GI_subsample_75_963.fa -option 1
BioSAK select_seq -seq /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_refined_bins/Refined_28.fasta -id NODE_4309_length_14678_cov_0.076146.txt -out NODE_4309_length_14678_cov_0.076146.fa -option 1


cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0418_mis2_no_clp_MarkerMAG_wd/GI_0418_mis2_no_clp_step_2_wd
blastn -query 3_GI_subsample_75_963.fa -subject /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/ref_genomes_renamed.fa -out 3_GI_subsample_75_963_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn
blastn -query NODE_11_length_1555_cov_0.052521.fa -subject /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/ref_genomes_renamed.fa -out NODE_11_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn
blastn -query NODE_4309_length_14678_cov_0.076146.fa -subject /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/ref_genomes_renamed.fa -out NODE_4309_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn


# ctg matched to (not mis-binned ctg)
CP002770.1  OTU_97.3402.0
Refined_28	OTU_97.3402.0
NODE_4309_length_14678_cov_0.076146	    CP002770.1	    99.986	14678	2	0	1	14678	701936	687259	0.0	26462	14678	3601386

# 16S matched to 
3_GI_subsample_75_963	FR745875.1	100.000	1494	0	0	1	1494	3691677	3690184	0.0	2695	1494	3781509
3_GI_subsample_75_963	CP001056.1	100.000	1494	0	0	1	1494	3691664	3690171	0.0	2695	1494	3800327
3_GI_subsample_75_963	CP001078.1	99.933	1494	1	0	1	1494	3543149	3541656	0.0	2691	1494	3659644
3_GI_subsample_75_963	CP010521.1	99.933	1494	1	0	1	1494	3495403	3493910	0.0	2691	1494	3611898
3_GI_subsample_75_963	CP010520.1	99.933	1494	1	0	1	1494	3495402	3493909	0.0	2691	1494	3611897
3_GI_subsample_75_963	CP006903.1	99.866	1494	2	0	1	1494	2392237	2390744	0.0	2686	1494	3874462

# mini-assembly matched to
NODE_11_length_1555_cov_0.052521	FR745875.1	99.486	1555	8	0	1	1555	3704931	3703377	0.0	2769	1555	3781509
NODE_11_length_1555_cov_0.052521	CP001056.1	99.550	1555	7	0	1	1555	9725	11279	0.0	2773	1555	3800327
NODE_11_length_1555_cov_0.052521	CP001078.1	99.550	1555	7	0	1	1555	3424658	3423104	0.0	2773	1555	3659644
NODE_11_length_1555_cov_0.052521	CP010521.1	99.550	1555	7	0	1	1555	3377284	3375730	0.0	2773	1555	3611898
NODE_11_length_1555_cov_0.052521	CP010520.1	99.550	1555	7	0	1	1555	3377283	3375729	0.0	2773	1555	3611897
NODE_11_length_1555_cov_0.052521	CP006903.1	99.357	1555	9	1	1	1555	2392131	2390578	0.0	2756	1555	3874462










'''


