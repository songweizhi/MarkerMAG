
depth_file_mag = '/Users/songweizhi/Desktop/666/hc_0513_with_depth_global_mean_depth_gnm.txt'
depth_file_16s = '/Users/songweizhi/Desktop/666/hc_0513_with_depth_global_mean_depth_16s.txt'
linkage_file   = '/Users/songweizhi/Desktop/666/hc_0514_stats_combined_filtered_with_assessment_ani97_imag99.5_i16S99.5.txt'


mean_depth_dict_gnm = {}
for each_mag_depth in open(depth_file_mag):
    each_mag_depth_split = each_mag_depth.strip().split('\t')
    mag_id = each_mag_depth_split[0]
    mag_depth = float(each_mag_depth_split[1])
    mean_depth_dict_gnm[mag_id] = mag_depth

mean_depth_dict_16s = {}
for each_16s_depth in open(depth_file_16s):
    each_16s_depth_split = each_16s_depth.strip().split('\t')
    s16_id = each_16s_depth_split[0]
    s16_depth = float(each_16s_depth_split[1])
    mean_depth_dict_16s[s16_id] = s16_depth

for each_link in open(linkage_file):
    each_link_split = each_link.strip().split('\t')
    id_16s = each_link_split[0]
    id_mag = each_link_split[1]

    depth_16s = mean_depth_dict_16s[id_16s]
    depth_mag = mean_depth_dict_gnm[id_mag]
    ratio = depth_16s/depth_mag
    depth_16s = float("{0:.2f}".format(depth_16s))
    depth_mag = float("{0:.2f}".format(depth_mag))
    ratio = float("{0:.2f}".format(ratio))

    print('%s\t%s\t%s\t%s' % (each_link.strip(), depth_16s, depth_mag, ratio))



