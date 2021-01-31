
def uniq_element_with_num(list_in):
    list_uniq = []
    for i in list_in:
        if i not in list_uniq:
            list_uniq.append(i)

    list_uniq_with_num = []
    for j in list_uniq:
        list_uniq_with_num.append('%s(%s)' % (j, list_in.count(j)))

    return list_uniq_with_num


bin_len_depth_file  = '/Users/songweizhi/Desktop/333/GI_bin_mean_depth.txt'
bin_quality_file    = '/Users/songweizhi/Desktop/333/GI_refined_bins_quality_refmt.txt'
bin_to_ref_file     = '/Users/songweizhi/Desktop/333/bin_vs_ref_imag99.9.txt'
ref_to_16s_num_file = '/Users/songweizhi/Desktop/333/GI_16S_16S_stats.txt'
bin_to_16s_num_file = '/Users/songweizhi/Desktop/333/GI_refined_bins_16S_stats.txt'
assess_results      = '/Users/songweizhi/Desktop/333/MarkerMAG_Mita_test_combined_linkages_with_assessment_ani95_imag99.9_i16S99.5.txt'


bin_id_list = set()
bin_to_len_dict = {}
bin_to_depth_dict = {}
for each_mag in open(bin_len_depth_file):
    each_mag_split = each_mag.strip().split('\t')
    mag_id = '.'.join(each_mag_split[0].split('.')[:-1])
    mag_len = each_mag_split[1]
    mag_depth = each_mag_split[2]
    bin_to_len_dict[mag_id] = mag_len
    bin_to_depth_dict[mag_id] = mag_depth
    bin_id_list.add(mag_id)

bin_to_cpl_dict = {}
for each_mag in open(bin_quality_file):
    each_mag_split = each_mag.strip().split(',')
    mag_id = each_mag_split[0]
    mag_cpl = each_mag_split[1]
    bin_to_cpl_dict[mag_id] = mag_cpl

bin_to_linkage_dict = {}
for each_linkage in open(assess_results):
    each_linkage_split = each_linkage.strip().split('\t')
    mag_id = each_linkage_split[1]
    link_assess = each_linkage_split[4]
    if mag_id not in bin_to_linkage_dict:
        bin_to_linkage_dict[mag_id] = [link_assess]
    else:
        bin_to_linkage_dict[mag_id].append(link_assess)

bin_to_ref_dict = {}
for each_mag in open(bin_to_ref_file):
    each_mag_split = each_mag.strip().split(',')
    mag_id = each_mag_split[0]
    mag_ref = each_mag_split[1]
    if mag_id not in bin_to_ref_dict:
        bin_to_ref_dict[mag_id] = {mag_ref}
    else:
        bin_to_ref_dict[mag_id].add(mag_ref)

ref_to_16s_num_dict = {}
for each_ref in open(ref_to_16s_num_file):
    each_ref_split = each_ref.strip().split('\t')
    ref_id = each_ref_split[0]
    num_16s = int(each_ref_split[1])
    ref_to_16s_num_dict[ref_id] = num_16s

bin_to_16s_num_dict = {}
for each_bin in open(bin_to_16s_num_file):
    each_bin_split = each_bin.strip().split('\t')
    bin_id = each_bin_split[0]
    num_16s_bin = int(each_bin_split[1])
    bin_to_16s_num_dict[bin_id] = num_16s_bin


print('MAG\tDepth\tCompleteness\tMatched_ref\tMatched_16S\tLinked')
for mag in bin_id_list:

    mag_linkage = ['no']
    if mag in bin_to_linkage_dict:
        mag_linkage = bin_to_linkage_dict[mag]

    matched_ref_num = 0
    matched_16s_num = 0
    if mag in bin_to_ref_dict:
        matched_ref_num = len(bin_to_ref_dict[mag])
        matched_16s_num = sum([ref_to_16s_num_dict[ref] for ref in bin_to_ref_dict[mag]])

    if len(mag_linkage) == 1:
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (mag, bin_to_depth_dict[mag], bin_to_cpl_dict[mag], bin_to_16s_num_dict[mag], matched_ref_num, matched_16s_num, mag_linkage[0]))
    else:
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (mag, bin_to_depth_dict[mag], bin_to_cpl_dict[mag], bin_to_16s_num_dict[mag], matched_ref_num, matched_16s_num, ','.join(sorted(uniq_element_with_num(mag_linkage)))))

