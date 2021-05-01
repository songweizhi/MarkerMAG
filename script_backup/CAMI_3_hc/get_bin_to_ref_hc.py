
bin_to_ref_txt   = '/Users/songweizhi/Desktop/assess_linkages_hc/stats_bin_to_ref_imag99.5.txt'
ref_with_16s_txt = '/Users/songweizhi/Desktop/assess_linkages_hc/file_in/ref_genomes_hc_with_16S.txt'


refs_with_16s_set = set()
for each_ref in open(ref_with_16s_txt):
    refs_with_16s_set.add(each_ref.strip()[4:])


n = 0
for eac_line in open(bin_to_ref_txt):
    eac_line_split = eac_line.strip().split('\t')
    bin_id = eac_line_split[0]
    ref_list = eac_line_split[1].split(',')
    ref_list_onyl_id = {i[4:] for i in ref_list}
    # print(eac_line_split)
    # print(bin_id)
    # print(ref_list_onyl_id)
    # print()

    ref_with_16s = False
    for each_ref in ref_list_onyl_id:
        if each_ref in refs_with_16s_set:
            ref_with_16s = True

    if ref_with_16s is False:
        print(bin_id)
        n += 1

print(n)


'''
cami_hc_110
cami_hc_115
cami_hc_116
cami_hc_119
cami_hc_122
cami_hc_133
cami_hc_150
cami_hc_154
cami_hc_53
cami_hc_64
cami_hc_80
cami_hc_87
cami_hc_99
'''
