
metabat_16S_depth   = '/Users/songweizhi/Desktop/666/000/reference_genomes_renamed_16S_depth.txt'
ref_depth_txt       = '/Users/songweizhi/Desktop/666/ref_depth.txt'
ref_16S_id_txt      = '/Users/songweizhi/Desktop/666/ref_16S_id.txt'

# output
ref_genome_metadata = '/Users/songweizhi/Desktop/666/reference_genome_metadata.txt'


ref_to_16s_dict = {}
for each_line in open(ref_16S_id_txt):
    gnm_id = each_line.strip().split('_')[0]
    if gnm_id not in ref_to_16s_dict:
        ref_to_16s_dict[gnm_id] = {each_line.strip()}
    else:
        ref_to_16s_dict[gnm_id].add(each_line.strip())


ref_depth_dict = {}
for each_ref_depth in open(ref_depth_txt):
    each_ref_depth_split = each_ref_depth.strip().split('\t')
    ref_id = each_ref_depth_split[0].split('_')[0]
    ref_depth = float(each_ref_depth_split[1])
    ref_depth_dict[ref_id] = ref_depth


gnm_total_weighted_16s_depth_dict = {}
gnm_16s_len_dict = {}
for each_line in open(metabat_16S_depth):
    if not each_line.startswith('contigName'):
        each_line_split = each_line.strip().split('\t')
        id_16s = each_line_split[0]
        id_gnm = id_16s.split('_16S_')[0]
        len_16s = int(each_line_split[1])
        depth_16s = float(each_line_split[2])
        depth_16s_weighted = len_16s*depth_16s
        if id_gnm not in gnm_total_weighted_16s_depth_dict:
            gnm_total_weighted_16s_depth_dict[id_gnm] = depth_16s_weighted
        else:
            gnm_total_weighted_16s_depth_dict[id_gnm] += depth_16s_weighted
        if id_gnm not in gnm_16s_len_dict:
            gnm_16s_len_dict[id_gnm] = [len_16s]
        else:
            gnm_16s_len_dict[id_gnm].append(len_16s)


ref_genome_metadata_handle = open(ref_genome_metadata, 'w')
ref_genome_metadata_handle.write('Genome\tReported_ref_16S_copy_num\tCalculated_ref_16S_copy_num\tCalculated_ref_16S_total_depth\tRef_depth\n')
for each_gnm in gnm_total_weighted_16s_depth_dict:
    gnm_total_weighted_16s_depth = gnm_total_weighted_16s_depth_dict[each_gnm]
    gnm_16s_len_list = gnm_16s_len_dict[each_gnm]
    gnm_16s_mean_len = sum(gnm_16s_len_list)/len(gnm_16s_len_list)
    gnm_16s_total_depth = gnm_total_weighted_16s_depth/gnm_16s_mean_len
    gnm_16s_total_depth = float("{0:.2f}".format(gnm_16s_total_depth))
    ref_gnm_depth = ref_depth_dict[each_gnm]
    reported_ref_16s_cp_num = len(ref_to_16s_dict[each_gnm])
    if ref_gnm_depth == 0:
        ref_16s_cp_num = 'NA'
    else:
        ref_16s_cp_num = gnm_16s_total_depth/ref_gnm_depth
        ref_16s_cp_num = float("{0:.2f}".format(ref_16s_cp_num))
    ref_genome_metadata_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_gnm, reported_ref_16s_cp_num, ref_16s_cp_num, gnm_16s_total_depth, ref_gnm_depth))
ref_genome_metadata_handle.close()
