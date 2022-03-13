import os

ref_genome_metadata         = '/Users/songweizhi/Desktop/666/reference_genome_metadata.txt'
gnm_id_txt                  = '/Users/songweizhi/Desktop/gnm_id.txt'
barrnap_op_16s_location_txt = '/Users/songweizhi/Desktop/reference_genomes_renamed_16S.txt'
depth_file_folder           = '/Users/songweizhi/Desktop/depth_files'
ref_gnm_folder              = '/Users/songweizhi/Desktop/reference_genomes_renamed'
markermag_op                = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_linkages_by_genome.txt'
flk_len                     = 10000


linked_16s_set = set()
gnm_to_linked_16s_dict = {}
for each_linkage in open(markermag_op):
    each_linkage_split = each_linkage.strip().split('\t')
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Round'):
        id_16s = each_linkage_split[0]
        id_gnm = each_linkage_split[1]
        linked_16s_set.add(id_16s)
        if id_gnm not in gnm_to_linked_16s_dict:
            gnm_to_linked_16s_dict[id_gnm] = {id_16s}
        else:
            gnm_to_linked_16s_dict[id_gnm].add(id_16s)

reported_ref_16s_copy_num_dict = {}
calculated_ref_16s_copy_num_dict = {}
ref_16s_total_depth_dict = {}
ref_gnm_depth_dict = {}
for each_ref in open(ref_genome_metadata):
    if not each_ref.startswith('Genome\t'):
        each_ref_split = each_ref.strip().split('\t')
        ref_id = each_ref_split[0]
        reported_copy_num = int(each_ref_split[1])

        if each_ref_split[2] == 'NA':
            calculated_copy_num = each_ref_split[2]
        else:
            calculated_copy_num = float(each_ref_split[2])

        ref_16s_total_depth = float(each_ref_split[3])
        ref_gnm_depth = float(each_ref_split[4])
        reported_ref_16s_copy_num_dict[ref_id] = reported_copy_num
        calculated_ref_16s_copy_num_dict[ref_id] = calculated_copy_num
        ref_16s_total_depth_dict[ref_id] = ref_16s_total_depth
        ref_gnm_depth_dict[ref_id] = ref_gnm_depth


gnm_to_16s_len_dict = {}
s16_ctg_to_gnm_dict = {}
gnm_to_16s_ctg_dict = {}
ctg_to_16s_region_dict = {}
ctg_to_16s_region_dict_sep = {}
for each_line in open(barrnap_op_16s_location_txt):
    if 'Length(bp)' not in each_line:
        each_line_split = each_line.strip().split('\t')
        gnm_id = each_line_split[0]
        len_16s = int(each_line_split[2])
        ctg_id = each_line_split[3]
        pos_start = int(each_line_split[4])
        pos_end = int(each_line_split[5])
        bp_list = list(range(pos_start, pos_end))

        if gnm_id not in gnm_to_16s_len_dict:
            gnm_to_16s_len_dict[gnm_id] = [len_16s]
        else:
            gnm_to_16s_len_dict[gnm_id].append(len_16s)

        if ctg_id not in s16_ctg_to_gnm_dict:
            s16_ctg_to_gnm_dict[ctg_id] = gnm_id

        if gnm_id not in gnm_to_16s_ctg_dict:
            gnm_to_16s_ctg_dict[gnm_id] = {ctg_id}
        else:
            gnm_to_16s_ctg_dict[gnm_id].add(ctg_id)

        if ctg_id not in ctg_to_16s_region_dict:
            ctg_to_16s_region_dict[ctg_id] = set()
        for each_bp in bp_list:
            ctg_to_16s_region_dict[ctg_id].add(each_bp)

        if ctg_id not in ctg_to_16s_region_dict_sep:
            ctg_to_16s_region_dict_sep[ctg_id] = []
        ctg_to_16s_region_dict_sep[ctg_id].append([pos_start, pos_end])


print('Genome\tMean_depth_16S\tDepth_gnm\tCalculated\tReported')
for each_gnm in open(gnm_id_txt):
    gnm_id = each_gnm.strip()

    if gnm_id in gnm_to_linked_16s_dict:

        current_gnm_16s_ctg = gnm_to_16s_ctg_dict[gnm_id]
        current_gnm_depth_file = '%s/%s_report_best_mis0_sorted_depth.txt' % (depth_file_folder, gnm_id)
        current_gnm_total_16s_depth = 0
        current_gnm_total_16s_len = 0
        for each_pos in open(current_gnm_depth_file):
            each_pos_split = each_pos.strip().split('\t')
            ctg_id = each_pos_split[0]
            ctg_pos = int(each_pos_split[1])
            pos_depth = int(each_pos_split[2])
            if (ctg_id in current_gnm_16s_ctg):
                current_ctg_16s_pos_set = ctg_to_16s_region_dict[ctg_id]
                if ctg_pos in current_ctg_16s_pos_set:
                    current_gnm_total_16s_depth += pos_depth
                    current_gnm_total_16s_len += 1

        if current_gnm_total_16s_len == 0:
            current_gnm_mean_16s_depth = 'NA'
        else:
            current_gnm_mean_16s_depth = current_gnm_total_16s_depth/current_gnm_total_16s_len
            current_gnm_mean_16s_depth = float("{0:.2f}".format(current_gnm_mean_16s_depth))

        ref_gnm_depth = ref_gnm_depth_dict[gnm_id]


        current_gnm_mean_16s_len = sum(gnm_to_16s_len_dict[gnm_id])/len(gnm_to_16s_len_dict[gnm_id])
        calculated_ref_16s_copy_num = (current_gnm_total_16s_depth/current_gnm_mean_16s_len)/ref_gnm_depth
        calculated_ref_16s_copy_num = float("{0:.2f}".format(calculated_ref_16s_copy_num))

        reported_ref_16s_copy_num = reported_ref_16s_copy_num_dict[gnm_id]

        print('%s\t%s\t%s\t%s\t%s' % (gnm_id, current_gnm_mean_16s_depth, ref_gnm_depth, calculated_ref_16s_copy_num, reported_ref_16s_copy_num))


# for each_gnm in gnm_to_linked_16s_dict:
#     gnm_id = each_gnm.strip()
#     gnm_16s_ctg_set = gnm_to_16s_ctg_dict[gnm_id]
#     for each_16s_ctg in gnm_16s_ctg_set:
#         ctg_16s_region_list = ctg_to_16s_region_dict_sep[each_16s_ctg]
#         for each_16s_region in ctg_16s_region_list:
#             left_end = each_16s_region[0]
#             right_end = each_16s_region[1]
#             plot_file_name = '%s__%s__%s-%s' % (gnm_id, each_16s_ctg, left_end, right_end)
#             pwd_ref_gnm = '%s/%s' % (ref_gnm_folder, gnm_id)
#             pwd_depth_file = '%s/%s_report_best_mis0_sorted_depth.txt' % (depth_file_folder, gnm_id)
#             plot_depth_cmd = 'BioSAK plot_sam_depth -r %s.fna -d %s -i %s -s %s -e %s -k 300 -l %s,%s' % (pwd_ref_gnm, pwd_depth_file, each_16s_ctg, (left_end-flk_len), (right_end+flk_len), left_end, right_end)
#             os.system(plot_depth_cmd)
