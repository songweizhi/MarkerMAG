import os


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def overlap_between_list(list_1, list_2):

    overlap = False
    for list_1_element in list_1:
        if list_1_element in list_2:
            overlap = True

    return overlap


########################################################################################################################

wd = '/Users/songweizhi/Desktop/MarkerMAG_wd/3_CAMI_GI/assess_linkages'

########## reference to cluster ##########

drep_ani_cutoff             = 95
drep_cdb_file               = '%s/file_in/Cdb_%s.csv'                                       % (wd, drep_ani_cutoff)
ref_to_strain_file          = '%s/file_in/ref_to_strain.txt'                                % wd

########## bin to reference ##########

parse_blastn_bin_vs_ref     = False  # True or False
blastn_bin_vs_ref           = '%s/file_in/bin_vs_ref.tab'                                   % wd
iden_cutoff                 = 99
aln_len_cutoff              = 1000
cov_q_cutoff                = 80
min_match_length            = 10240  # 10 Kbp
mag_metadata                = '%s/file_in/MAG_metadata.txt'                                 % wd

########## 16S to reference ##########

matam_16s_blastn            = '%s/file_in/3_GI_assembled_16S_uclust_0.999_vs_ref.tab'       % wd
iden_cutoff_16s             = 99.5  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 97

########## script ##########

pwd_plot_sankey_R = '%s/file_in/get_sankey_plot.R' % wd

########## assessment results ##########

MarkerMAG_linkages          = '%s/GI_3_70_10_5_combined_linkages.txt'       % wd
# Linkage	156/97(160.82)	156/174(89.66)	21

MarkerMAG_linkages          = '%s/GI_3_80_10_5_combined_linkages.txt'       % wd
# Linkage	151/97(155.67)	151/166(90.96)	19

MarkerMAG_linkages          = '%s/GI_0_70_10_5_combined_linkages.txt'       % wd

'''

GI_0_70_10_5    Linkage	89/97(91.75)	89/95(93.68)	9	Correctly_linked_MAG(29)
GI_1_70_10_5    Linkage	105/97(108.25)	105/108(97.22)	9	Correctly_linked_MAG(34)
GI_2_70_10_5    Linkage	136/97(140.21)	136/150(90.67)	17	Correctly_linked_MAG(41)
GI_3_70_10_5    Linkage	156/97(160.82)	156/174(89.66)	21	Correctly_linked_MAG(43)


GI_1_80_10_5    Linkage	102/97(105.15)	102/106(96.23)	8
GI_3_80_10_5    Linkage	151/97(155.67)	151/166(90.96)	19

'''




############################################### define file/folder name ################################################

# bin to reference
cluster_to_ref_file         = '%s/cluster_to_ref.txt'                       % wd
bin_vs_ref_txt              = '%s/bin_vs_ref_imag%s.txt'                    % (wd, iden_cutoff)
stats_bin_to_ref_txt        = '%s/stats_bin_to_ref_imag%s.txt'              % (wd, iden_cutoff)
stats_ref_to_bin_txt        = '%s/stats_ref_to_bin_imag%s.txt'              % (wd, iden_cutoff)
stats_bin_to_cluster_txt    = '%s/stats_bin_to_cluster_ani%s_imag%s.txt'    % (wd, drep_ani_cutoff, iden_cutoff)
stats_cluster_to_bin_txt    = '%s/stats_cluster_to_bin_ani%s_imag%s.txt'    % (wd, drep_ani_cutoff, iden_cutoff)

# assessment results
linkage_file_path, linkage_file_basename, linkage_file_extension = sep_path_basename_ext(MarkerMAG_linkages)
MarkerMAG_linkages_assessed = '%s/%s_with_assessment_ani%s_imag%s_i16S%s.txt'   % (wd, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)
wrong_linkages_txt          = '%s/%s_wrong_ani%s_imag%s_i16S%s.txt'             % (wd, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)
unknown_linkages_txt        = '%s/%s_unknown_ani%s_imag%s_i16S%s.txt'           % (wd, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)


################################################# reference to cluster #################################################

cluster_to_ref_dict = {}
ref_to_cluster_dict = {}
for each_ref in open(drep_cdb_file):
    if not each_ref.startswith('genome,secondary_cluster'):
        each_ref_split = each_ref.strip().split(',')
        ref_file_name = each_ref_split[0]
        ref_file_name_no_ext = '.'.join(ref_file_name.split('.')[:-1])
        ref_cluster = 'C' + each_ref_split[1]
        ref_to_cluster_dict[ref_file_name_no_ext] = ref_cluster
        if ref_cluster not in cluster_to_ref_dict:
            cluster_to_ref_dict[ref_cluster] = [ref_file_name_no_ext]
        else:
            cluster_to_ref_dict[ref_cluster].append(ref_file_name_no_ext)

cluster_to_ref_file_handle = open(cluster_to_ref_file, 'w')
for each_cluster in cluster_to_ref_dict:
    cluster_to_ref_file_handle.write('%s\t%s\n' % (each_cluster, ','.join(cluster_to_ref_dict[each_cluster])))
cluster_to_ref_file_handle.close()


################################################### bin to reference ###################################################

bin_ref_connector           = '__|__'

# get ref_to_strain_dict
ref_to_strain_dict = {}
for ref in open(ref_to_strain_file):
    ref_split = ref.strip().split('\t')
    ref_to_strain_dict[ref_split[0]] = ref_split[1]

if parse_blastn_bin_vs_ref is True:
    bin_vs_ref_dict = {}
    for match in open(blastn_bin_vs_ref):
        match_split = match.strip().split('\t')
        query = match_split[0]
        query_genome = '_'.join(query.split('_')[:2])
        subject = match_split[1]
        subject_genome = '_'.join(subject.split('_')[:2])
        iden = float(match_split[2])
        aln_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        coverage_q = float(aln_len) * 100 / float(query_len)
        if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff) and (coverage_q >= cov_q_cutoff):
            key_bin_ref = '%s%s%s' % (query_genome, bin_ref_connector, subject_genome)
            if key_bin_ref not in bin_vs_ref_dict:
                bin_vs_ref_dict[key_bin_ref] = aln_len
            else:
                bin_vs_ref_dict[key_bin_ref] += aln_len

    bin_vs_ref_dict_filtered = {}
    for each_key in bin_vs_ref_dict:
        if bin_vs_ref_dict[each_key] >= min_match_length:
            bin_vs_ref_dict_filtered[each_key] = bin_vs_ref_dict[each_key]

    # write out linkages
    linkage_txt_handle = open(bin_vs_ref_txt, 'w')
    linkage_txt_handle.write('Bin,Ref,Length\n')
    bin_to_ref_dict = {}
    ref_to_bin_dict = {}
    bin_to_cluster_dict = {}
    cluster_to_bin_dict = {}
    for each_linkage in bin_vs_ref_dict_filtered:
        linkage_split = each_linkage.split(bin_ref_connector)
        bin_id = linkage_split[0]
        ref_id = linkage_split[1]
        ref_cluster = ref_to_cluster_dict[ref_id]

        if bin_id not in bin_to_ref_dict:
            bin_to_ref_dict[bin_id] = {ref_id}
        else:
            bin_to_ref_dict[bin_id].add(ref_id)

        if bin_id not in bin_to_cluster_dict:
            bin_to_cluster_dict[bin_id] = {ref_cluster}
        else:
            bin_to_cluster_dict[bin_id].add(ref_cluster)

        if ref_id not in ref_to_bin_dict:
            ref_to_bin_dict[ref_id] = {bin_id}
        else:
            ref_to_bin_dict[ref_id].add(bin_id)

        if ref_cluster not in cluster_to_bin_dict:
            cluster_to_bin_dict[ref_cluster] = {bin_id}
        else:
            cluster_to_bin_dict[ref_cluster].add(bin_id)

        linkage_txt_handle.write('%s,%s,%s\n' % (bin_id, ref_id, bin_vs_ref_dict_filtered[each_linkage]))
    linkage_txt_handle.close()

    # visualize
    cmd_sankey_bin_to_ref = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, bin_vs_ref_txt, 600, 1800)
    os.system(cmd_sankey_bin_to_ref)

    bin_to_ref_txt_handle = open(stats_bin_to_ref_txt, 'w')
    for each_bin in bin_to_ref_dict:
        bin_to_ref_txt_handle.write('%s\t%s\n' % (each_bin, ','.join(bin_to_ref_dict[each_bin])))
    bin_to_ref_txt_handle.close()

    ref_to_bin_txt_handle = open(stats_ref_to_bin_txt, 'w')
    for each_ref in ref_to_bin_dict:
        ref_to_bin_txt_handle.write('%s\t%s\t%s\n' % (each_ref, ','.join(ref_to_bin_dict[each_ref]), ref_to_strain_dict[each_ref]))
    ref_to_bin_txt_handle.close()

    stats_bin_to_cluster_txt_handle = open(stats_bin_to_cluster_txt, 'w')
    for each_bin in bin_to_cluster_dict:
        stats_bin_to_cluster_txt_handle.write('%s\t%s\n' % (each_bin, ','.join(bin_to_cluster_dict[each_bin])))
    stats_bin_to_cluster_txt_handle.close()

    stats_cluster_to_bin_txt_handle = open(stats_cluster_to_bin_txt, 'w')
    for each_ref in cluster_to_bin_dict:
        stats_cluster_to_bin_txt_handle.write('%s\t%s\n' % (each_ref, ','.join(cluster_to_bin_dict[each_ref])))
    stats_cluster_to_bin_txt_handle.close()

bin_to_cluster_dict = {}
for each_match in open(stats_bin_to_cluster_txt):
    each_match_split = each_match.strip().split('\t')
    bin_to_cluster_dict[each_match_split[0]] = {i for i in each_match_split[1].split(',')}


cluster_to_bin_dict = {}
for each_bin in bin_to_cluster_dict:
    matched_clusters =  bin_to_cluster_dict[each_bin]
    for matched_cluster in matched_clusters:
        if matched_cluster not in cluster_to_bin_dict:
            cluster_to_bin_dict[matched_cluster] = {each_bin}
        else:
            cluster_to_bin_dict[matched_cluster].add(each_bin)

################################################### 16S to reference ###################################################

# get matam_16s_to_cluster_dict
matam_16s_to_cluster_dict = {}
for match in open(matam_16s_blastn):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_gnm = subject.split('_16S_')[0]
    subject_cluster = ref_to_cluster_dict[subject_gnm]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    if (iden >= iden_cutoff_16s) and (aln_len > aln_len_cutoff_16s) and (coverage_q >= cov_q_cutoff_16s):
        if query not in matam_16s_to_cluster_dict:
            matam_16s_to_cluster_dict[query] = {subject_cluster}
        else:
            matam_16s_to_cluster_dict[query].add(subject_cluster)

cluster_to_matam_16s_dict = {}
for matam_16s in matam_16s_to_cluster_dict:
    matched_clusters =  matam_16s_to_cluster_dict[matam_16s]
    for matched_cluster in matched_clusters:
        if matched_cluster not in cluster_to_matam_16s_dict:
            cluster_to_matam_16s_dict[matched_cluster] = {matam_16s}
        else:
            cluster_to_matam_16s_dict[matched_cluster].add(matam_16s)


###################################################### assessment ######################################################

MarkerMAG_linkages_assessed_handle = open(MarkerMAG_linkages_assessed, 'w')
wrong_linkages_txt_handle = open(wrong_linkages_txt, 'w')
unknown_linkages_txt_handle = open(unknown_linkages_txt, 'w')

linkage_mag_right = set()
linkage_num_right = 0
linkage_num_wrong = 0
linkage_num_unknown = 0
linkage_assessment_dict = {}
for each_linkage in open(MarkerMAG_linkages):
    if each_linkage.startswith('MarkerGene\tGenomicSeq\tLinkage\tStep'):
        MarkerMAG_linkages_assessed_handle.write('MarkerGene\tGenomicSeq\tLinkage\tStep\tAssessment\n')
    else:
        each_linkage_split = each_linkage.strip().split('\t')
        id_16s = each_linkage_split[0]
        id_mag = each_linkage_split[1]
        key_16s_mag = '%s___%s' % (id_16s, id_mag)

        matched_cluster_16s = {}
        if id_16s in matam_16s_to_cluster_dict:
            matched_cluster_16s = matam_16s_to_cluster_dict[id_16s]

        matched_cluster_mag = {}
        if id_mag in bin_to_cluster_dict:
            matched_cluster_mag = bin_to_cluster_dict[id_mag]

        if (matched_cluster_16s != {}) and (matched_cluster_mag != {}):

            if overlap_between_list(matched_cluster_mag, matched_cluster_16s) is True:
                linkage_num_right += 1
                linkage_mag_right.add(id_mag)
                MarkerMAG_linkages_assessed_handle.write('%s\tCorrect\n' % each_linkage.strip())
                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Correct']
                else:
                    linkage_assessment_dict[id_mag].append('Correct')

            else:
                print(matched_cluster_mag)
                print(matched_cluster_16s)
                print()
                linkage_num_wrong += 1
                MarkerMAG_linkages_assessed_handle.write('%s\tWrong\n' % each_linkage.strip())
                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Wrong']
                else:
                    linkage_assessment_dict[id_mag].append('Wrong')

                for each_16s_cluster in matched_cluster_16s:
                    current_16s_cluster_ref = cluster_to_ref_dict[each_16s_cluster]
                    for each_ref_1 in current_16s_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict[each_ref_1]))

                for each_mag_cluster in matched_cluster_mag:
                    current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                    for each_ref_2 in current_mag_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict[each_ref_2]))
                wrong_linkages_txt_handle.write('\n')
        else:
            linkage_num_unknown += 1
            MarkerMAG_linkages_assessed_handle.write('%s\tUnknown\n' % each_linkage.strip())
            if id_mag not in linkage_assessment_dict:
                linkage_assessment_dict[id_mag] = ['Unknown']
            else:
                linkage_assessment_dict[id_mag].append('Unknown')

            for each_16s_cluster in matched_cluster_16s:
                current_16s_cluster_ref = cluster_to_ref_dict[each_16s_cluster]
                for each_ref_1 in current_16s_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict[each_ref_1]))
            for each_mag_cluster in matched_cluster_mag:
                current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                for each_ref_2 in current_mag_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict[each_ref_2]))
            unknown_linkages_txt_handle.write('\n')
wrong_linkages_txt_handle.close()
unknown_linkages_txt_handle.close()
MarkerMAG_linkages_assessed_handle.close()


# report
link_recovery = linkage_num_right*100/total_query_mag_num
link_accuracy = linkage_num_right*100/(linkage_num_right + linkage_num_wrong)
link_recovery = float("{0:.2f}".format(link_recovery))
link_accuracy = float("{0:.2f}".format(link_accuracy))
recovery_str  = '%s/%s(%s)' % (linkage_num_right, total_query_mag_num, link_recovery)
accuracy_str  = '%s/%s(%s)' % (linkage_num_right, (linkage_num_right + linkage_num_wrong), link_accuracy)

print('%s\tRecovery\tAccuracy\tUnknown' % 'Linkage')
print('%s\t%s\t%s\t%s\tCorrectly_linked_MAG(%s)' % ('Linkage', recovery_str, accuracy_str, linkage_num_unknown, len(linkage_mag_right)))



########################################################################################################################

# linkage_assessment_stats_dict = {}
# for each_assess_key in linkage_assessment_dict:
#     assess_result = linkage_assessment_dict[each_assess_key]
#     assess_result_sorted = sorted(assess_result)
#     assess_result_uniq = []
#     for each_a in assess_result_sorted:
#         if each_a not in assess_result_uniq:
#             assess_result_uniq.append(each_a)
#     assess_stats = []
#     for uniq_a in assess_result_uniq:
#         assess_stats.append('%s(%s)' % (uniq_a, assess_result_sorted.count(uniq_a)))
#     linkage_assessment_stats_dict[each_assess_key] = assess_stats
#
#
# for each_mag in open(mag_metadata):
#     if not each_mag.startswith('MAG	Depth	Completeness	16S_in_MAG	Matched_ref'):
#         each_mag_split = each_mag.strip().split('\t')
#         mag_id = each_mag_split[0]
#         mag_cluster = {'unknown'}
#         if mag_id in bin_to_cluster_dict:
#             mag_cluster = bin_to_cluster_dict[mag_id]
#
#         link_result = {'no'}
#         if mag_id in linkage_assessment_stats_dict:
#             link_result = linkage_assessment_stats_dict[mag_id]
#
#         Matam16S_matched = set()
#         for each_cluster in mag_cluster:
#             if each_cluster in cluster_to_matam_16s_dict:
#                 for each_16s in cluster_to_matam_16s_dict[each_cluster]:
#                     Matam16S_matched.add(each_16s)
#
#         print('%s\t%s\t%s\t%s' % (each_mag.strip(), ','.join(mag_cluster), ','.join(link_result), len(Matam16S_matched)))
#


# # print('MAG\tDepth\tCompleteness\t16S_in_MAG\tMatched_ref\tCluster\tLinkage\tMatam16S')
# # print(linkage_assessment_stats_dict)
# # print(cluster_to_matam_16s_dict)
#
# for bin in bin_to_cluster_dict:
#     print('%s\t%s' % (bin, bin_to_cluster_dict[bin]))

