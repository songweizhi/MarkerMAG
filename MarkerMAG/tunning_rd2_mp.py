import os
import pandas as pd


def blast_results_to_pairwise_16s_iden_dict(blastn_output, align_len_cutoff, cov_cutoff):
    pairwise_iden_dict = {}
    for match in open(blastn_output):
        match_split = match.strip().split('\t')
        query = match_split[0]
        subject = match_split[1]
        iden = float(match_split[2])
        align_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        coverage_q = float(align_len) * 100 / float(query_len)
        coverage_s = float(align_len) * 100 / float(subject_len)

        if (align_len >= align_len_cutoff) and (query != subject) and (coverage_q >= cov_cutoff) and (
                coverage_s >= cov_cutoff):
            query_to_subject_key = '__|__'.join(sorted([query, subject]))
            if query_to_subject_key not in pairwise_iden_dict:
                pairwise_iden_dict[query_to_subject_key] = iden
            else:
                if iden > pairwise_iden_dict[query_to_subject_key]:
                    pairwise_iden_dict[query_to_subject_key] = iden

    return pairwise_iden_dict


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = ''

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def sort_csv_by_col(file_in, file_out, col_header):
    df_in = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, within_gnm_linkage_num_diff, file_out):

    # get MarkerGene_to_GenomicSeq_dict
    MarkerGene_to_GenomicSeq_dict = {}
    for each_linkage in open(file_in):
        if not ((each_linkage.startswith('MarkerGene,GenomicSeq,Number')) or (each_linkage.startswith('MarkerGene,MiniAssembly,Number'))):
            each_linkage_split = each_linkage.strip().split(',')
            MarkerGene_id = each_linkage_split[0][12:]
            GenomicSeq_id = each_linkage_split[1][12:]
            linkage_num = int(each_linkage_split[2])
            if linkage_num > 1:
                if MarkerGene_id not in MarkerGene_to_GenomicSeq_dict:
                    MarkerGene_to_GenomicSeq_dict[MarkerGene_id] = {GenomicSeq_id}
                else:
                    MarkerGene_to_GenomicSeq_dict[MarkerGene_id].add(GenomicSeq_id)

    file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
    file_in_sorted = '%s/%s_sorted%s' % (file_in_path, file_in_basename, file_in_extension)

    # sort file in
    sort_csv_by_col(file_in, file_in_sorted, sort_by_col_header)

    # fileter linkage
    gnm_max_link_num_dict = {}
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    GenomicSeq_best_marker_dict = {}
    for each_match in open(file_in_sorted):
        if (each_match.startswith('MarkerGene,GenomicSeq,Number')) or (each_match.startswith('MarkerGene,MiniAssembly,Number')):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            MarkerGene = match_split[0][12:]
            GenomicSeq = match_split[1][12:]
            linkage_num = int(match_split[2])

            current_min_linkage = min_linkages_for_uniq_linked_16s
            if MarkerGene in MarkerGene_to_GenomicSeq_dict:
                if len(MarkerGene_to_GenomicSeq_dict[MarkerGene]) > 1:
                    current_min_linkage = min_linkages

            if linkage_num >= current_min_linkage:
                if MarkerGene not in MarkerGene_with_assignment:

                    # consider depth
                    if min_16s_gnm_multiple > 0:

                        # get marker and genome depth
                        MarkerGene_depth = marker_gene_depth_dict.get(MarkerGene, 'na')
                        GenomicSeq_depth = genomic_seq_depth_dict.get(GenomicSeq, 'na')
                        marker_genome_depth_ratio = 'na'
                        if (MarkerGene_depth != 'na') and (GenomicSeq_depth != 'na'):
                            if GenomicSeq_depth > 0:
                                marker_genome_depth_ratio = MarkerGene_depth / GenomicSeq_depth

                        if marker_genome_depth_ratio >= min_16s_gnm_multiple:
                            if GenomicSeq not in GenomicSeq_best_marker_dict:
                                GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                                gnm_max_link_num_dict[GenomicSeq] = linkage_num
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)
                            else:
                                # get identity with best marker
                                current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                                key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))
                                iden_with_best_marker = pairwise_16s_iden_dict.get(key_str, 0)
                                if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                                    gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
                                    if (linkage_num * 100 / gnm_max_link_num) >= within_gnm_linkage_num_diff:
                                        file_out_handle.write(each_match)
                                        MarkerGene_with_assignment.add(MarkerGene)
                                    else:
                                        MarkerGene_with_assignment.add(MarkerGene)
                    # ignore depth
                    else:
                        if GenomicSeq not in GenomicSeq_best_marker_dict:
                            GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                            gnm_max_link_num_dict[GenomicSeq] = linkage_num
                            file_out_handle.write(each_match)
                            MarkerGene_with_assignment.add(MarkerGene)
                        else:
                            # get identity with best marker
                            current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                            key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))
                            iden_with_best_marker = pairwise_16s_iden_dict.get(key_str, 0)
                            if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                                gnm_max_link_num = gnm_max_link_num_dict[GenomicSeq]
                                if (linkage_num*100/gnm_max_link_num) >= within_gnm_linkage_num_diff:
                                    file_out_handle.write(each_match)
                                    MarkerGene_with_assignment.add(MarkerGene)
                                else:
                                    MarkerGene_with_assignment.add(MarkerGene)
    file_out_handle.close()


def filter_linkages_iteratively_mini_assembly_to_ctg(file_in_sorted, min_linkages, file_out):

    # do mini-assemblies assigned to the same mag need to have roughly the same number of linkages? think about this later
    mag_ctg_max_link_num_dict = {}
    mini_assembly_to_mag_dict = {}
    file_out_handle = open(file_out, 'w')
    mini_assembly_with_assignment = set()
    for each_match in open(file_in_sorted):
        if each_match.startswith('MiniAssembly,GenomicSeq,Number'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            mini_assembly = match_split[0]
            mag_ctg_id = match_split[1]
            mag_id = mag_ctg_id.split('___C___')[0]

            linkage_num = int(match_split[2])
            if linkage_num >= min_linkages:
                if mini_assembly not in mini_assembly_with_assignment:
                    if mag_ctg_id not in mag_ctg_max_link_num_dict:
                        mag_ctg_max_link_num_dict[mag_ctg_id] = linkage_num
                        file_out_handle.write(each_match)
                        mini_assembly_to_mag_dict[mini_assembly] = mag_id
                        mini_assembly_with_assignment.add(mini_assembly)
                    else:
                        ratio_with_best_assignment = linkage_num/(mag_ctg_max_link_num_dict[mag_ctg_id])
                        if ratio_with_best_assignment >= 0.8:
                            file_out_handle.write(each_match)
                            mini_assembly_to_mag_dict[mini_assembly] = mag_id
                            mini_assembly_with_assignment.add(mini_assembly)
                        else:
                            mini_assembly_with_assignment.add(mini_assembly)
    file_out_handle.close()

    return mini_assembly_to_mag_dict


def get_GapFilling_stats_by_assembly_separately(free_living_16s_ref_file, free_living_ctg_ref_file,
                                                mini_assembly_to_16s_reads, mini_assembly_to_ctg_reads, ctg_level_min_link,
                                                mini_assembly_to_16s_ctg_connector, gnm_to_ctg_connector, marker_to_ctg_gnm_Key_connector,
                                                stats_mini_assembly_to_ctg, stats_mini_assembly_to_ctg_sorted, stats_mini_assembly_to_ctg_filtered,
                                                stats_GapFilling_ctg, stats_GapFilling_gnm):

    ########## link mini-assembly to MAGs ##########

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

    stats_mini_assembly_to_ctg_handle = open(stats_mini_assembly_to_ctg, 'w')
    stats_mini_assembly_to_ctg_handle.write('MiniAssembly,GenomicSeq,Number\n')
    for each_mini_assembly_to_ctg in mini_assembly_to_ctg_dict:
        link_num = mini_assembly_to_ctg_dict[each_mini_assembly_to_ctg]
        id_mini_assembly = each_mini_assembly_to_ctg.split(mini_assembly_to_16s_ctg_connector)[0]
        id_ctg = each_mini_assembly_to_ctg.split(mini_assembly_to_16s_ctg_connector)[1]
        stats_mini_assembly_to_ctg_handle.write('%s,%s,%s\n' % (id_mini_assembly, id_ctg, link_num))
    stats_mini_assembly_to_ctg_handle.close()

    # sort  and filter
    sort_csv_by_col(stats_mini_assembly_to_ctg, stats_mini_assembly_to_ctg_sorted, 'Number')
    os.remove(stats_mini_assembly_to_ctg)
    mini_assembly_to_mag_dict = filter_linkages_iteratively_mini_assembly_to_ctg(stats_mini_assembly_to_ctg_sorted, 3, stats_mini_assembly_to_ctg_filtered)

    # visualize linkages


    ########## link 16S to mini-assemblies with MAG assignment ##########

    round2_free_living_16s_ref_dict = {}
    for free_living_read_16s in open(free_living_16s_ref_file):
        free_living_read_16s_split = free_living_read_16s.strip().split('\t')
        if len(free_living_read_16s_split) > 1:
            read_16s_id = free_living_read_16s_split[0]
            read_16s_refs = free_living_read_16s_split[1].split(',')
            round2_free_living_16s_ref_dict[read_16s_id] = read_16s_refs

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

    stats_mini_assembly_to_16s_handle = open(stats_GapFilling_ctg, 'w')
    stats_mini_assembly_to_16s_handle.write('MarkerGene,MiniAssembly,Number\n')
    marker_to_gnm_link_num_dict_rd2 = {}
    for each_mini_assembly_to_16s in mini_assembly_to_16s_dict:
        link_num = mini_assembly_to_16s_dict[each_mini_assembly_to_16s]
        id_mini_assembly = each_mini_assembly_to_16s.split(mini_assembly_to_16s_ctg_connector)[0]
        id_16s = each_mini_assembly_to_16s.split(mini_assembly_to_16s_ctg_connector)[1]
        mini_assembly_mag = mini_assembly_to_mag_dict.get(id_mini_assembly, None)
        if link_num >= ctg_level_min_link:
            if mini_assembly_mag != None:
                stats_mini_assembly_to_16s_handle.write('%s,%s%s%s,%s\n' % (id_16s, mini_assembly_mag, gnm_to_ctg_connector, id_mini_assembly, link_num))
                marker_to_gnm_key_rd2 = '%s%s%s' % (id_16s, marker_to_ctg_gnm_Key_connector, mini_assembly_mag)
                if marker_to_gnm_key_rd2 not in marker_to_gnm_link_num_dict_rd2:
                    marker_to_gnm_link_num_dict_rd2[marker_to_gnm_key_rd2] = link_num
                else:
                    marker_to_gnm_link_num_dict_rd2[marker_to_gnm_key_rd2] += link_num
    stats_mini_assembly_to_16s_handle.close()

    # write out linkages at genome level
    stats_GapFilling_gnm_handle = open(stats_GapFilling_gnm, 'w')
    stats_GapFilling_gnm_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_linkage in marker_to_gnm_link_num_dict_rd2:
        stats_GapFilling_gnm_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
            each_linkage.split(marker_to_ctg_gnm_Key_connector)[0],
            each_linkage.split(marker_to_ctg_gnm_Key_connector)[1],
            marker_to_gnm_link_num_dict_rd2[each_linkage]))
    stats_GapFilling_gnm_handle.close()


####################################################### file in ########################################################

step_2_wd                                       = '/Users/songweizhi/Desktop/tunning_rd2'
free_living_16s_ref_file                        = '%s/file_in/round2_free_living_16s_refs.txt'                          % step_2_wd
free_living_ctg_ref_file                        = '%s/file_in/round2_free_living_ctg_refs.txt'                          % step_2_wd
mini_assembly_to_16s_reads                      = '%s/file_in/mini_assembly_to_16s_reads.txt'                           % step_2_wd
mini_assembly_to_ctg_reads                      = '%s/file_in/mini_assembly_to_ctg_reads.txt'                           % step_2_wd
blast_results_all_vs_all_16s                    = '%s/file_in/Oral_0622_60_60_polish_new_16S_all_vs_all_blastn.tab'     % step_2_wd
ctg_level_min_link                              = 3
within_gnm_linkage_num_diff                     = 80
max_mini_assembly_link_num_diff_between_ctg_16s = 10
min_aln_16s                                     = 500
min_cov_16s                                     = 30
min_iden_16s                                    = 98
min_link_num                                    = 8
marker_to_ctg_gnm_Key_connector                 = '___M___'
gnm_to_ctg_connector                            = '___C___'
mini_assembly_to_16s_ctg_connector              = '___Mini___'
read_to_marker_connector                        = '___r___'
mean_depth_dict_gnm                             = {}
mean_depth_dict_16s                             = {}
min_16s_gnm_multiple                            = 0

####################################################### file out #######################################################

stats_mini_assembly_to_ctg                      = '%s/stats_mini_assembly_to_ctg.txt'                                   % step_2_wd
stats_mini_assembly_to_ctg_sorted               = '%s/stats_mini_assembly_to_ctg_sorted.txt'                            % step_2_wd
stats_mini_assembly_to_ctg_filtered             = '%s/stats_mini_assembly_to_ctg_filtered.txt'                          % step_2_wd
stats_GapFilling_ctg                            = '%s/stats_GapFilling_ctg.txt'                                         % step_2_wd
stats_GapFilling_file                           = '%s/stats_GapFilling_gnm.txt'                                         % step_2_wd
stats_GapFilling_file_filtered                  = '%s/stats_GapFilling_gnm_filtered.txt'                                % step_2_wd

########################################################################################################################

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

get_GapFilling_stats_by_assembly_separately(free_living_16s_ref_file,
                                 free_living_ctg_ref_file,
                                 mini_assembly_to_16s_reads,
                                 mini_assembly_to_ctg_reads,
                                 ctg_level_min_link,
                                 mini_assembly_to_16s_ctg_connector,
                                 gnm_to_ctg_connector,
                                 marker_to_ctg_gnm_Key_connector,
                                 stats_mini_assembly_to_ctg,
                                 stats_mini_assembly_to_ctg_sorted,
                                 stats_mini_assembly_to_ctg_filtered,
                                 stats_GapFilling_ctg,
                                 stats_GapFilling_file)

filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm,
                            mean_depth_dict_16s, min_16s_gnm_multiple, min_iden_16s, min_link_num, min_link_num,
                            within_gnm_linkage_num_diff, stats_GapFilling_file_filtered)


os.system('python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/d7_Oral/assess_linkages_Oral.py')


