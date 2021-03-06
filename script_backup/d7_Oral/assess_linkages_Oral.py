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


def gnm_level_asessment(linked_mag_dict_rd1):

    rd1_correct_num = 0
    rd1_wrong_gnm = set()
    rd1_unknown_gnm = set()
    rd1_ambiguous_gnm = set()
    for each_rd_1_mag in linked_mag_dict_rd1:
        mag_assess = linked_mag_dict_rd1[each_rd_1_mag]
        if len(mag_assess) == 1:
            if mag_assess == {'Correct'}:
                rd1_correct_num += 1
            elif mag_assess == {'Wrong'}:
                rd1_wrong_gnm.add(each_rd_1_mag)
            elif mag_assess == {'Unknown'}:
                rd1_unknown_gnm.add(each_rd_1_mag)
        elif len(mag_assess) == 2:
            if ('Correct' in mag_assess) and ('Unknown' in mag_assess):
                rd1_correct_num += 1
            if ('Correct' in mag_assess) and ('Wrong' in mag_assess):
                rd1_ambiguous_gnm.add(each_rd_1_mag)
            if ('Wrong' in mag_assess) and ('Unknown' in mag_assess):
                rd1_wrong_gnm.add(each_rd_1_mag)
        elif len(mag_assess) == 3:
            rd1_ambiguous_gnm.add(each_rd_1_mag)

    return rd1_correct_num, rd1_unknown_gnm, rd1_wrong_gnm, rd1_ambiguous_gnm


########################################################################################################################

wd = '/Users/songweizhi/Desktop/assess_linkages_Oral'

########## reference to cluster ##########

drep_ani_cutoff             = 97
drep_cdb_file               = '%s/file_in/Cdb_%s.csv'                                       % (wd, drep_ani_cutoff)
# ref_to_strain_file          = '%s/file_in/ref_to_strain.txt'                                % wd

########## bin to reference ##########

parse_blastn_bin_vs_ref     = False  # True or False
blastn_bin_vs_ref           = '%s/file_in/bin_vs_ref.tab'                                   % wd
iden_cutoff                 = 99.5
aln_len_cutoff              = 1500
cov_q_cutoff                = 90
min_match_length            = 524288  # 100 Kbp  102400
# mag_metadata                = '%s/file_in/MAG_metadata.txt'                                 % wd

########## 16S to reference ##########

perform_blastn_16s_vs_refs  = False  # True or False
combined_GI_ref_16S         = '%s/file_in/combined_Oral_ref_16S.ffn'                     % wd
matam_16s_seqs              = '%s/file_in/CAMI_Oral_138_16S_0.999.polished.fa'           % wd
matam_16s_blastn            = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_vs_ref.tab'   % wd
iden_cutoff_16s             = 99.5  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 87

########## assessment results ##########

MarkerMAG_linkages          = '%s/Oral_0708_60_60_min1200_mask_linkages_by_genome.txt'             % wd
linkages_from_rd1           = False


# MarkerMAG_linkages          = '/Users/songweizhi/Desktop/tunning_rd2/stats_GapFilling_gnm_filtered.txt'
# linkages_from_rd1           = True

'''

Oral_0708_60_60_min1200_mask_linkages_by_genome.txt	Rd_1	|	50	43	1	6	43/49(87.76)	|	16	14	0	2	0	14/16(87.5)
Oral_0708_60_60_min1200_mask_linkages_by_genome.txt	Rd_2	|	16	12	0	4	12/16(75.0)	    |	4	2	0	2	0	2/4(50.0)
Oral_0708_60_60_min1200_mask_linkages_by_genome.txt	Both	|	66	55	1	10	55/65(84.62)	|	20	16	0	4	0	16/20(80.0)

Oral_0708_60_60_min1500_mask_linkages_by_genome.txt	Rd_1	|	46	40	0	6	40/46(86.96)	|	16	14	0	2	0	14/16(87.5)
Oral_0708_60_60_min1500_mask_linkages_by_genome.txt	Rd_2	|	9	1	0	8	1/9(11.11)	    |	2	1	0	1	0	1/2(50.0)
Oral_0708_60_60_min1500_mask_linkages_by_genome.txt	Both	|	55	41	0	14	41/55(74.55)	|	18	15	0	3	0	15/18(83.33)


Oral_0705_60_60_min1500_linkages_by_genome.txt	Rd_1	|	45	42	0	3	42/45(93.33)	|	16	15	0	1	0	15/16(93.75)
Oral_0705_60_60_min1500_linkages_by_genome.txt	Rd_2	|	5	5	0	0	5/5(100.0)	    |	2	2	0	0	0	2/2(100.0)
Oral_0705_60_60_min1500_linkages_by_genome.txt	Both	|	50	47	0	3	47/50(94.0)	    |	18	17	0	1	0	17/18(94.44)

Oral_0705_60_60_min1400_linkages_by_genome.txt	Rd_1	|	38	37	0	1	37/38(97.37)	|	15	14	0	1	0	14/15(93.33)
Oral_0705_60_60_min1400_linkages_by_genome.txt	Rd_2	|	17	1	0	16	1/17(5.88)	    |	3	0	0	3	0	0/3(0.0)
Oral_0705_60_60_min1400_linkages_by_genome.txt	Both	|	55	38	0	17	38/55(69.09)	|	18	14	0	4	0	14/18(77.78)

Oral_0622_60_60_polish_min1200_99_127_linkages_by_genome.txt	Rd_1	|	44	39	1	4	39/43(90.7)	    |	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1200_99_127_linkages_by_genome.txt	Rd_2	|	20	6	0	14	6/20(30.0)	    |	3	1	0	2	0	1/3(33.33)
Oral_0622_60_60_polish_min1200_99_127_linkages_by_genome.txt	Both	|	64	45	1	18	45/63(71.43)	|	19	16	0	3	0	16/19(84.21)

Oral_0622_60_60_polish_new_linkages_by_genome.txt	Rd_1	|	46	44	1	1	44/45(97.78)	|	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_new_linkages_by_genome.txt	Rd_2	|	23	9	0	14	9/23(39.13)	    |	4	2	0	2	0	2/4(50.0)
Oral_0622_60_60_polish_new_linkages_by_genome.txt	Both	|	69	53	1	15	53/68(77.94)	|	20	17	0	3	0	17/20(85.0)

Oral_0622_60_60_polish_min1200_75_127_linkages_by_genome.txt	Rd_1	|	45	43	1	1	43/44(97.73)	|	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1200_75_127_linkages_by_genome.txt	Rd_2	|	23	7	0	16	7/23(30.43)	    |	4	1	0	3	0	1/4(25.0)
Oral_0622_60_60_polish_min1200_75_127_linkages_by_genome.txt	Both	|	68	50	1	17	50/67(74.63)	|	20	16	0	4	0	16/20(80.0)

Oral_0622_60_60_polish_min1200_linkages_by_genome.txt	Rd_1	|	50	44	1	5	44/49(89.8)	    |	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1200_linkages_by_genome.txt	Rd_2	|	14	11	0	3	11/14(78.57)	|	3	1	0	2	0	1/3(33.33)
Oral_0622_60_60_polish_min1200_linkages_by_genome.txt	Both	|	64	55	1	8	55/63(87.3)	    |	19	16	0	3	0	16/19(84.21)

Oral_0622_60_60_polish_min1000_linkages_by_genome.txt	Rd_1	|	48	42	1	5	42/47(89.36)	|	16	15	0	1	0	15/16(93.75)
Oral_0622_60_60_polish_min1000_linkages_by_genome.txt	Rd_2	|	25	8	0	17	8/25(32.0)	    |	5	2	0	3	0	2/5(40.0)
Oral_0622_60_60_polish_min1000_linkages_by_genome.txt	Both	|	73	50	1	22	50/72(69.44)	|	21	17	0	4	0	17/21(80.95)

Oral_0622_60_60_polish_min900_linkages_by_genome.txt	Rd_1	|	48	45	1	2	45/47(95.74)	|	16	14	0	1	1	14/16(87.5)
Oral_0622_60_60_polish_min900_linkages_by_genome.txt	Rd_2	|	33	14	0	19	14/33(42.42)	|	7	3	0	4	0	3/7(42.86)
Oral_0622_60_60_polish_min900_linkages_by_genome.txt	Both	|	81	59	1	21	59/80(73.75)	|	23	17	0	5	1	17/23(73.91)

Oral_0622_60_60_polish_linkages_by_genome.txt	Rd_1	|	63	52	8	3	52/55(94.55)	|	20	16	1	2	1	16/19(84.21)
Oral_0622_60_60_polish_linkages_by_genome.txt	Rd_2	|	14	10	0	4	10/14(71.43)	|	6	3	0	3	0	3/6(50.0)
Oral_0622_60_60_polish_linkages_by_genome.txt	Both	|	77	62	8	7	62/69(89.86)	|	26	19	1	5	1	19/25(76.0)

'''

########## script ##########

pwd_plot_sankey_R = '%s/file_in/get_sankey_plot.R' % wd

############################################### define file/folder name ################################################

# bin to reference
bin_vs_ref_txt              = '%s/bin_vs_ref_imag%s.txt'                    % (wd, iden_cutoff)
stats_bin_to_ref_txt        = '%s/stats_bin_to_ref_imag%s.txt'              % (wd, iden_cutoff)
stats_ref_to_bin_txt        = '%s/stats_ref_to_bin_imag%s.txt'              % (wd, iden_cutoff)
stats_bin_to_cluster_txt    = '%s/stats_bin_to_cluster_ani%s_imag%s.txt'    % (wd, drep_ani_cutoff, iden_cutoff)
stats_cluster_to_bin_txt    = '%s/stats_cluster_to_bin_ani%s_imag%s.txt'    % (wd, drep_ani_cutoff, iden_cutoff)

# assessment results
linkage_file_path, linkage_file_basename, linkage_file_extension = sep_path_basename_ext(MarkerMAG_linkages)
MarkerMAG_linkages_assessed = '%s/%s_with_assessment_ani%s_imag%s_i16S%s.txt'   % (linkage_file_path, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)
wrong_linkages_txt          = '%s/%s_wrong_ani%s_imag%s_i16S%s.txt'             % (linkage_file_path, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)
unknown_linkages_txt        = '%s/%s_unknown_ani%s_imag%s_i16S%s.txt'           % (linkage_file_path, linkage_file_basename, drep_ani_cutoff, iden_cutoff, iden_cutoff_16s)


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


################################################### bin to reference ###################################################

bin_ref_connector           = '__|__'

# get ref_to_strain_dict
ref_to_strain_dict = {}
# for ref in open(ref_to_strain_file):
#     ref_split = ref.strip().split('\t')
#     ref_to_strain_dict[ref_split[0]] = ref_split[1]

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
        ref_to_bin_txt_handle.write('%s\t%s\t%s\n' % (each_ref, ','.join(ref_to_bin_dict[each_ref]), ref_to_strain_dict.get(each_ref, 'NA')))
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

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 1'
blast_cmd = 'blastn -query %s -subject %s -out %s %s' % (matam_16s_seqs, combined_GI_ref_16S, matam_16s_blastn, blast_parameters)
if perform_blastn_16s_vs_refs is True:
    os.system(blast_cmd)

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
linkage_num_right = 0
linkage_num_wrong = 0
linkage_num_unknown = 0
linkage_assessment_dict = {}
for each_linkage in open(MarkerMAG_linkages):
    if ('MarkerGene\tGenomicSeq\tLinkage\tRound' in each_linkage) or ('MarkerGene,GenomicSeq,Number' in each_linkage):
        MarkerMAG_linkages_assessed_handle.write('MarkerGene\tGenomicSeq\tLinkage\tStep\tAssessment\n')
    else:
        id_16s = ''
        id_mag = ''
        link_num = ''
        if linkages_from_rd1 is False:
            each_linkage_split = each_linkage.strip().split('\t')
            id_16s = each_linkage_split[0]
            id_mag = each_linkage_split[1]
            link_num = each_linkage_split[2]
        else:
            each_linkage_split = each_linkage.strip().split(',')
            id_16s = each_linkage_split[0][12:]
            id_mag = each_linkage_split[1][12:]
            link_num = each_linkage_split[2]

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

                if linkages_from_rd1 is False:
                    MarkerMAG_linkages_assessed_handle.write('%s\tCorrect\n' % each_linkage.strip())
                else:
                    MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tCorrect\n' % (id_16s, id_mag, link_num))

                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Correct']
                else:
                    linkage_assessment_dict[id_mag].append('Correct')

            else:
                linkage_num_wrong += 1

                if linkages_from_rd1 is False:
                    MarkerMAG_linkages_assessed_handle.write('%s\tWrong\n' % each_linkage.strip())
                else:
                    MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tWrong\n' % (id_16s, id_mag, link_num))

                if id_mag not in linkage_assessment_dict:
                    linkage_assessment_dict[id_mag] = ['Wrong']
                else:
                    linkage_assessment_dict[id_mag].append('Wrong')

                for each_16s_cluster in matched_cluster_16s:
                    current_16s_cluster_ref = cluster_to_ref_dict[each_16s_cluster]
                    for each_ref_1 in current_16s_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict.get(each_ref_1, 'NA')))

                for each_mag_cluster in matched_cluster_mag:
                    current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                    for each_ref_2 in current_mag_cluster_ref:
                        wrong_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_wrong, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict.get(each_ref_2, 'NA')))
                wrong_linkages_txt_handle.write('\n')
        else:
            linkage_num_unknown += 1

            if linkages_from_rd1 is False:
                MarkerMAG_linkages_assessed_handle.write('%s\tUnknown\n' % each_linkage.strip())
            else:
                MarkerMAG_linkages_assessed_handle.write('%s\t%s\t%s\tRd1\tUnknown\n' % (id_16s, id_mag, link_num))

            if id_mag not in linkage_assessment_dict:
                linkage_assessment_dict[id_mag] = ['Unknown']
            else:
                linkage_assessment_dict[id_mag].append('Unknown')

            for each_16s_cluster in matched_cluster_16s:
                current_16s_cluster_ref = cluster_to_ref_dict[each_16s_cluster]
                for each_ref_1 in current_16s_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, '16S', each_16s_cluster, each_ref_1, ref_to_strain_dict.get(each_ref_1, 'NA')))
            for each_mag_cluster in matched_cluster_mag:
                current_mag_cluster_ref = cluster_to_ref_dict[each_mag_cluster]
                for each_ref_2 in current_mag_cluster_ref:
                    unknown_linkages_txt_handle.write('%s\t%s___%s\t%s\t%s\t%s\t%s\n' % (linkage_num_unknown, id_16s, id_mag, 'MAG', each_mag_cluster, each_ref_2, ref_to_strain_dict.get(each_ref_2, 'NA')))
            unknown_linkages_txt_handle.write('\n')
wrong_linkages_txt_handle.close()
unknown_linkages_txt_handle.close()
MarkerMAG_linkages_assessed_handle.close()


# get assess stats at mag level
linked_mag_dict_rd1 = {}
linked_mag_dict_rd2 = {}
total_linkage_num_rd1 = 0
total_linkage_num_rd2 = 0
correct_link_rd1 = 0
correct_link_rd2 = 0
wrong_link_rd1 = 0
wrong_link_rd2 = 0
unknown_link_rd1 = 0
unknown_link_rd2 = 0
for each_link in open(MarkerMAG_linkages_assessed):
    if not each_link.startswith('MarkerGene	GenomicSeq	Linkage	Step	Assessment'):
        each_link_split = each_link.strip().split('\t')
        gnm_id = each_link_split[1]
        linked_rd = each_link_split[3]
        assessment = each_link_split[4]
        if linked_rd == 'Rd1':
            total_linkage_num_rd1 += 1

            if assessment == 'Correct':
                correct_link_rd1 += 1
            if assessment == 'Wrong':
                wrong_link_rd1 += 1
            if assessment == 'Unknown':
                unknown_link_rd1 += 1

            if gnm_id not in linked_mag_dict_rd1:
                linked_mag_dict_rd1[gnm_id] = {assessment}
            else:
                linked_mag_dict_rd1[gnm_id].add(assessment)
        if linked_rd == 'Rd2':
            total_linkage_num_rd2 += 1

            if assessment == 'Correct':
                correct_link_rd2 += 1
            if assessment == 'Wrong':
                wrong_link_rd2 += 1
            if assessment == 'Unknown':
                unknown_link_rd2 += 1

            if gnm_id not in linked_mag_dict_rd1:
                linked_mag_dict_rd2[gnm_id] = {assessment}
            else:
                linked_mag_dict_rd2[gnm_id].add(assessment)


total_linkage_num_both = total_linkage_num_rd1 + total_linkage_num_rd2
total_linked_mag_num = len(linked_mag_dict_rd1) + len(linked_mag_dict_rd2)

rd1_correct_num, rd1_unknown_gnm, rd1_wrong_gnm, rd1_ambiguous_gnm = gnm_level_asessment(linked_mag_dict_rd1)
rd2_correct_num, rd2_unknown_gnm, rd2_wrong_gnm, rd2_ambiguous_gnm = gnm_level_asessment(linked_mag_dict_rd2)

correct_num_both        = rd1_correct_num + rd2_correct_num
unknown_gnm_both_num    = len(rd1_unknown_gnm) + len(rd2_unknown_gnm)
wrong_gnm_both_num      = len(rd1_wrong_gnm) + len(rd2_wrong_gnm)
ambiguous_gnm_both_num  = len(rd1_ambiguous_gnm) + len(rd2_ambiguous_gnm)

recovery_str_rd1   = '%s/%s(%s)' % (rd1_correct_num,  total_query_mag_num, float("{0:.2f}".format(rd1_correct_num*100/total_query_mag_num)))
recovery_str_rd2   = '%s/%s(%s)' % (rd2_correct_num,  total_query_mag_num, float("{0:.2f}".format(rd2_correct_num*100/total_query_mag_num)))
recovery_str_both  = '%s/%s(%s)' % (correct_num_both, total_query_mag_num, float("{0:.2f}".format(correct_num_both*100/total_query_mag_num)))
accuracy_str_rd1   = '%s/%s(%s)' % (rd1_correct_num,  (len(linked_mag_dict_rd1) - len(rd1_unknown_gnm)), float("{0:.2f}".format(rd1_correct_num*100/(len(linked_mag_dict_rd1) - len(rd1_unknown_gnm)))))
accuracy_str_rd2 = '0/0(0)'
if (len(linked_mag_dict_rd2) - len(rd2_unknown_gnm)) > 0:
    accuracy_str_rd2   = '%s/%s(%s)' % (rd2_correct_num,  (len(linked_mag_dict_rd2) - len(rd2_unknown_gnm)), float("{0:.2f}".format(rd2_correct_num*100/(len(linked_mag_dict_rd2) - len(rd2_unknown_gnm)))))
accuracy_str_both  = '%s/%s(%s)' % (correct_num_both, (total_linked_mag_num - unknown_gnm_both_num), float("{0:.2f}".format(correct_num_both*100/(total_linked_mag_num - unknown_gnm_both_num))))

correct_link_both = correct_link_rd1 + correct_link_rd2
wrong_link_both   = wrong_link_rd1 + wrong_link_rd2
unknown_link_both = unknown_link_rd1 + unknown_link_rd2

accuracy_str_rd1_link_level   = '%s/%s(%s)' % (correct_link_rd1,  (total_linkage_num_rd1 - unknown_link_rd1), float("{0:.2f}".format(correct_link_rd1*100/(total_linkage_num_rd1 - unknown_link_rd1))))
accuracy_str_rd2_link_level = '0/0(0)'
if (total_linkage_num_rd2 - unknown_link_rd2) > 0:
    accuracy_str_rd2_link_level   = '%s/%s(%s)' % (correct_link_rd2,  (total_linkage_num_rd2 - unknown_link_rd2), float("{0:.2f}".format(correct_link_rd2*100/(total_linkage_num_rd2 - unknown_link_rd2))))
accuracy_str_both_link_level  = '%s/%s(%s)' % (correct_link_both, (total_linkage_num_both - unknown_link_both), float("{0:.2f}".format(correct_link_both*100/(total_linkage_num_both - unknown_link_both))))

prefix = os.path.basename(MarkerMAG_linkages).split('_identified_linkages_genome_level')[0]
print('%s\tRound\t|\tLink\tYes\tNA\tNo\tAccuracy\t|\tMAG\tYes\tNA\tNo\tY/N\tAccuracy' % prefix)
print('%s\tRd_1\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd1,  correct_link_rd1,  unknown_link_rd1,  wrong_link_rd1,  accuracy_str_rd1_link_level,  len(linked_mag_dict_rd1), rd1_correct_num,  len(rd1_unknown_gnm), len(rd1_wrong_gnm), len(rd1_ambiguous_gnm), accuracy_str_rd1))
print('%s\tRd_2\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_rd2,  correct_link_rd2,  unknown_link_rd2,  wrong_link_rd2,  accuracy_str_rd2_link_level,  len(linked_mag_dict_rd2), rd2_correct_num,  len(rd2_unknown_gnm), len(rd2_wrong_gnm), len(rd2_ambiguous_gnm), accuracy_str_rd2))
print('%s\tBoth\t|\t%s\t%s\t%s\t%s\t%s\t|\t%s\t%s\t%s\t%s\t%s\t%s' % (prefix, total_linkage_num_both, correct_link_both, unknown_link_both, wrong_link_both, accuracy_str_both_link_level, total_linked_mag_num,      correct_num_both, unknown_gnm_both_num, wrong_gnm_both_num, ambiguous_gnm_both_num, accuracy_str_both))


########################################################################################################################

'''

########################################################################################################################

CAMI_Oral_subsample_50_1551	Oral_57	2123	Rd1	Wrong

CAMI_Oral_subsample_50_1551     C125_1,C125_3,C126_1,C126_2,C126_3
Oral_57                         C126_5

# blast between 16S from Oral_57 ref gnm and Matam assemblies
OTU_97.23667.1_16S	CAMI_Oral_subsample_50_571	100.000	1542	0	0	1	1542	5	1546	0.0	2848

# blast between CAMI_Oral_subsample_50_1551 and 16S from ref gnms

BioSAK iTOL -ColorRange -lg leaf_group_raw.txt -lt Identity -out leaf_group.txt 

########################################################################################################################

'''
# cluster_id = 'C159_1'
# print('cluster to ref  : %s\t%s\t%s' % (cluster_id, len(cluster_to_ref_dict.get(cluster_id, 'NA')), cluster_to_ref_dict.get(cluster_id, 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % (cluster_id, len(cluster_to_bin_dict.get(cluster_id, 'NA')), cluster_to_bin_dict.get(cluster_id, 'NA')))
# print('cluster to Matam: %s\t%s\t%s' % (cluster_id, len(cluster_to_matam_16s_dict.get(cluster_id, 'NA')), cluster_to_matam_16s_dict.get(cluster_id, 'NA')))
# #print('cluster to MAG: %s\t%s' % ('C126_5', cluster_to_bin_dict.get('C126_5', 'NA')))

# print(bin_to_cluster_dict['Oral_70'])
#
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_1', len(cluster_to_bin_dict.get('C126_1', 'NA')), cluster_to_bin_dict.get('C126_1', 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_2', len(cluster_to_bin_dict.get('C126_2', 'NA')), cluster_to_bin_dict.get('C126_2', 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_3', len(cluster_to_bin_dict.get('C126_3', 'NA')), cluster_to_bin_dict.get('C126_3', 'NA')))
# print('cluster to MAG  : %s\t%s\t%s' % ('C126_4', len(cluster_to_bin_dict.get('C126_4', 'NA')), cluster_to_bin_dict.get('C126_4', 'NA')))
#
# print(matam_16s_to_cluster_dict['CAMI_Oral_subsample_50_767'])
# print(cluster_to_matam_16s_dict['C118_4'])


'''
OTU_97.16455.0
OTU_97.18445.1
OTU_97.30379.0
OTU_97.3393.0
OTU_97.41293.0
OTU_97.41776.0
OTU_97.42505.0
OTU_97.43450.0
OTU_97.44958.0
OTU_97.45053.1
OTU_97.45058.0
OTU_97.4910.0




'''
