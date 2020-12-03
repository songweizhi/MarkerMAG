
def overlap_between_list(list_1, list_2):

    overlap = False
    for list_1_element in list_1:
        if list_1_element in list_2:
            overlap = True

    return overlap

iden = 99.9
wd                      = '/Users/songweizhi/Desktop/ttt_mag'
ref_to_strain_file      = '%s/ref_to_strain.txt'                % wd
bin_to_ref_txt          = '%s/stats_bin_to_ref_%s.txt'          % (wd, iden)
s16_to_ref_txt          = '%s/stats_16s_to_ref.txt'             % wd
MarkerMAG_linkages      = '%s/CAMI2_HMP_combined_linkages.txt'  % wd
#MarkerMAG_linkages     = '%s/CAMI2_HMP_s2m30_combined_linkages.txt'  % wd
total_query_mag_num     = 97

wrong_linkages_txt      = '%s/wrong_linkages_%s.txt'           % (wd, iden)
unknown_linkages_txt    = '%s/unknown_linkages_%s.txt'         % (wd, iden)


# get ref_to_strain_dict
ref_to_strain_dict = {}
for ref in open(ref_to_strain_file):
    ref_split = ref.strip().split('\t')
    ref_to_strain_dict[ref_split[0]] = ref_split[1]


# get bin_to_ref_dict
bin_to_ref_dict = {}
for each_bin in open(bin_to_ref_txt):
    each_bin_split = each_bin.strip().split('\t')
    bin_id = each_bin_split[0]
    matched_refs = each_bin_split[1].split(',')
    bin_to_ref_dict[bin_id] = matched_refs


# get s16_to_ref_dict
s16_to_ref_dict = {}
for each_16s in open(s16_to_ref_txt):
    each_16s_split = each_16s.strip().split('\t')
    s16_id = each_16s_split[0]
    matched_refs = each_16s_split[1].split(',')
    s16_to_ref_dict[s16_id] = matched_refs


wrong_linkages_txt_handle = open(wrong_linkages_txt, 'w')
unknown_linkages_txt_handle = open(unknown_linkages_txt, 'w')
linkage_num_right = 0
linkage_num_wrong = 0
linkage_num_unknown = 0
for each_linkage in open(MarkerMAG_linkages):
    if not each_linkage.startswith('MarkerGene\tGenomicSeq\tLinkage\tStep'):
        each_linkage_split = each_linkage.strip().split('\t')
        id_16s = each_linkage_split[0]
        id_mag = each_linkage_split[1]
        if (id_16s in s16_to_ref_dict) and (id_mag in bin_to_ref_dict):
            matched_refs_16s = s16_to_ref_dict[id_16s]
            matched_refs_mag = bin_to_ref_dict[id_mag]
            if overlap_between_list(matched_refs_16s, matched_refs_mag) is True:
                linkage_num_right += 1
            else:
                linkage_num_wrong += 1
                wrong_linkages_txt_handle.write(each_linkage)
                for i in matched_refs_16s:
                    wrong_linkages_txt_handle.write('%s\t%s\t%s\n' % ('16S', i, ref_to_strain_dict[i]))
                for j in matched_refs_mag:
                    wrong_linkages_txt_handle.write('%s\t%s\t%s\n' % ('MAG', j, ref_to_strain_dict[j]))
                wrong_linkages_txt_handle.write('\n')
        else:
            linkage_num_unknown += 1
            unknown_linkages_txt_handle.write(each_linkage)
            if id_16s in s16_to_ref_dict:
                for m in s16_to_ref_dict[id_16s]:
                    unknown_linkages_txt_handle.write('16S\t%s\t%s\n' % (m, ref_to_strain_dict[m]))
            else:
                unknown_linkages_txt_handle.write('16S\n')
            if id_mag in bin_to_ref_dict:
                for n in bin_to_ref_dict[id_mag]:
                    unknown_linkages_txt_handle.write('MAG\t%s\t%s\n' % (n, ref_to_strain_dict[n]))
            else:
                unknown_linkages_txt_handle.write('MAG\n')
            unknown_linkages_txt_handle.write('\n')
wrong_linkages_txt_handle.close()
unknown_linkages_txt_handle.close()


link_recovery = linkage_num_right*100/total_query_mag_num
link_accuracy = linkage_num_right*100/(linkage_num_right + linkage_num_wrong)
link_recovery = float("{0:.2f}".format(link_recovery))
link_accuracy = float("{0:.2f}".format(link_accuracy))
recovery_str  = '%s/%s(%s)' % (linkage_num_right, total_query_mag_num, link_recovery)
accuracy_str  = '%s/%s(%s)' % (linkage_num_right, (linkage_num_right + linkage_num_wrong), link_accuracy)

print('%s\tRecovery\tAccuracy\tUnknown' % 'Linkage')
print('%s\t%s\t%s\t%s' % ('Linkage', recovery_str, accuracy_str, linkage_num_unknown))


'''



'''
