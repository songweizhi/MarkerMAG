import os
import pandas as pd
from Bio import SeqIO


def sort_csv_by_col(file_in, file_out, col_header):
    df_in        = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


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

        if (align_len >= align_len_cutoff) and (query != subject) and (coverage_q >= cov_cutoff) and (coverage_s >= cov_cutoff):
            query_to_subject_key = '__|__'.join(sorted([query, subject]))
            if query_to_subject_key not in pairwise_iden_dict:
                pairwise_iden_dict[query_to_subject_key] = iden
            else:
                if iden > pairwise_iden_dict[query_to_subject_key]:
                    pairwise_iden_dict[query_to_subject_key] = iden

    return pairwise_iden_dict


def get_max_clp_and_index(r1_cigar_list, r2_cigar_list):
    r1_cigar_list_split = [cigar_splitter(i) for i in r1_cigar_list]
    r2_cigar_list_split = [cigar_splitter(i) for i in r2_cigar_list]

    r1_cigar_list_split_only_clp = []
    for each_r1_cigar_split in r1_cigar_list_split:
        clp_len_l = 0
        if each_r1_cigar_split[0][-1] in ['S', 's']:
            clp_len_l = int(each_r1_cigar_split[0][:-1])
        clp_len_r = 0
        if each_r1_cigar_split[-1][-1] in ['S', 's']:
            clp_len_r = int(each_r1_cigar_split[-1][:-1])
        r1_cigar_list_split_only_clp.append([clp_len_l, clp_len_r])

    r2_cigar_list_split_only_clp = []
    for each_r2_cigar_split in r2_cigar_list_split:
        clp_len_l = 0
        if each_r2_cigar_split[0][-1] in ['S', 's']:
            clp_len_l = int(each_r2_cigar_split[0][:-1])
        clp_len_r = 0
        if each_r2_cigar_split[-1][-1] in ['S', 's']:
            clp_len_r = int(each_r2_cigar_split[-1][:-1])
        r2_cigar_list_split_only_clp.append([clp_len_l, clp_len_r])

    cigar_list_split_only_clp_r1_r2 = [r1_cigar_list_split_only_clp, r2_cigar_list_split_only_clp]

    max_value = 0
    max_value_index = ''
    for num_list_1 in cigar_list_split_only_clp_r1_r2[0]:
        if num_list_1[0] > max_value:
            max_value = num_list_1[0]
            max_value_index = 'r1_l'
        if num_list_1[1] > max_value:
            max_value = num_list_1[1]
            max_value_index = 'r1_r'
    for num_list_2 in cigar_list_split_only_clp_r1_r2[1]:
        if num_list_2[0] > max_value:
            max_value = num_list_2[0]
            max_value_index = 'r2_l'
        if num_list_2[1] > max_value:
            max_value = num_list_2[1]
            max_value_index = 'r2_r'

    return max_value, max_value_index


class MappingRecord:

    def __init__(self):

        self.r1_seq = ''
        self.r2_seq = ''
        self.r1_refs = dict()
        self.r2_refs = dict()

        self.qualified_reads           = False
        self.consider_r1_unmapped_mate = False
        self.consider_r1_clipping_part = False
        self.consider_r2_unmapped_mate = False
        self.consider_r2_clipping_part = False

        self.r1_filtered_refs = set()
        self.r2_filtered_refs = set()

        self.r1_clipping_seq = ''
        self.r2_clipping_seq = ''

        self.unmapped_r1_refs = set()
        self.unmapped_r2_refs = set()

        self.clipping_r1_refs = set()
        self.clipping_r2_refs = set()


def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if (each_element.isalpha() is True) or (each_element == '='):
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


def get_cigar_stats(cigar_splitted):

    # aligned_len: M I X =
    # clipping_len: S
    # mismatch_len: X I D
    # mismatch_pct = mismatch_len / aligned_len
    # aligned_pct  = aligned_len  / (aligned_len + clipping_len)
    # clipping_pct = clipping_len / (aligned_len + clipping_len)

    aligned_len = 0
    clipping_len = 0
    mismatch_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]

        # get aligned_len
        if each_part_cate in {'M', 'm', 'I', 'i', 'X', 'x', '='}:
            aligned_len += each_part_len

        # get clipping_len
        if each_part_cate in ['S', 's']:
            clipping_len += each_part_len

        # get mismatch_len
        if each_part_cate in {'I', 'i', 'X', 'x', 'D', 'd'}:
            mismatch_len += each_part_len

    aligned_pct  = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
    clipping_pct = float("{0:.2f}".format(clipping_len * 100 / (aligned_len + clipping_len)))
    mismatch_pct = float("{0:.2f}".format(mismatch_len * 100 / (aligned_len)))

    return aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, file_out):

    # get MarkerGene_to_GenomicSeq_dict
    MarkerGene_to_GenomicSeq_dict = {}
    for each_linkage in open(file_in):
        if not each_linkage.startswith('MarkerGene,GenomicSeq,Number'):
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
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    GenomicSeq_best_marker_dict = {}
    for each_match in open(file_in_sorted):
        if each_match.startswith('MarkerGene,GenomicSeq,Number'):
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

                # consider depth
                if min_16s_gnm_multiple > 0:
                    MarkerGene_depth = marker_gene_depth_dict[MarkerGene]
                    GenomicSeq_depth = genomic_seq_depth_dict[GenomicSeq]
                    if (MarkerGene_depth/GenomicSeq_depth) >= min_16s_gnm_multiple:
                        if MarkerGene not in MarkerGene_with_assignment:

                            if GenomicSeq not in GenomicSeq_best_marker_dict:
                                GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)
                            else:
                                current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                                key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))

                                iden_with_best_marker = 0
                                if key_str in pairwise_16s_iden_dict:
                                    iden_with_best_marker = pairwise_16s_iden_dict[key_str]

                                if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                                    file_out_handle.write(each_match)
                                    MarkerGene_with_assignment.add(MarkerGene)
                # ignore depth
                else:
                    if MarkerGene not in MarkerGene_with_assignment:
                        if GenomicSeq not in GenomicSeq_best_marker_dict:
                            GenomicSeq_best_marker_dict[GenomicSeq] = MarkerGene
                            file_out_handle.write(each_match)
                            MarkerGene_with_assignment.add(MarkerGene)
                        else:
                            current_GenomicSeq_best_marker = GenomicSeq_best_marker_dict[GenomicSeq]
                            key_str = '__|__'.join(sorted([MarkerGene, current_GenomicSeq_best_marker]))

                            iden_with_best_marker = 0
                            if key_str in pairwise_16s_iden_dict:
                                iden_with_best_marker = pairwise_16s_iden_dict[key_str]

                            if iden_with_best_marker >= within_genome_16s_divergence_cutoff:
                                file_out_handle.write(each_match)
                                MarkerGene_with_assignment.add(MarkerGene)

    file_out_handle.close()

    # remove tmp file
    # os.remove(file_in_sorted)


########################################################################################################################

# file in
pwd_samfile                         = '/Users/songweizhi/Desktop/new_algorithm_Kelp/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.sam'
min_M_pct                           = 20
min_M_len                           = 30
max_mis_pct                         = 3
min_clp_len                         = 35

# file out
pwd_samfile_reads                   = '/Users/songweizhi/Desktop/new_algorithm_Kelp/Kelp_SILVA138_id99_assembled_16S_uclust_0.999_reads.fa'
clipping_part_seq                   = '/Users/songweizhi/Desktop/new_algorithm_Kelp/clipping_parts.fa'
unmapped_mates                      = '/Users/songweizhi/Desktop/new_algorithm_Kelp/unmapped_mates.fa'
clipping_part_seq_sam               = '/Users/songweizhi/Desktop/new_algorithm_Kelp/clipping_parts.sam'
unmapped_mates_sam                  = '/Users/songweizhi/Desktop/new_algorithm_Kelp/unmapped_mates.sam'
stats_combined                      = '/Users/songweizhi/Desktop/new_algorithm_Kelp/stats_combined.txt'
stats_combined_filtered             = '/Users/songweizhi/Desktop/new_algorithm_Kelp/stats_combined_filtered.txt'
blast_results_all_vs_all_16s        = '/Users/songweizhi/Desktop/new_algorithm_Kelp/16S_all_vs_all_blastn.tab'
mock_final_op                       = '/Users/songweizhi/Desktop/new_algorithm_Kelp/mock_final_op.txt'

marker_to_ctg_Key_connector_str = '___M___'
marker_to_gnm_Key_connector_str = '___M___'
gnm_ctg_connector = '___'


########################################################################################################################

pwd_samfile_reads_handle = open(pwd_samfile_reads, 'w')
samfile_reads_wrote = set()
MappingRecord_dict = {}
for each_read in open(pwd_samfile):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_seq = each_read_split[9]
        cigar = each_read_split[5]
        if read_id not in samfile_reads_wrote:
            pwd_samfile_reads_handle.write('>%s\n' % read_id)
            pwd_samfile_reads_handle.write('%s\n' % read_seq)
            samfile_reads_wrote.add(read_id)
        if cigar != '*':
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = each_read_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

            if (aligned_len >= min_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= (max_mis_pct)):
                if read_id_base not in MappingRecord_dict:
                    MappingRecord_dict[read_id_base] = MappingRecord()
                if read_strand == '1':
                    MappingRecord_dict[read_id_base].r1_refs[ref_id_with_pos] = cigar
                    MappingRecord_dict[read_id_base].r1_seq = read_seq
                if read_strand == '2':
                    MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar
                    MappingRecord_dict[read_id_base].r2_seq = read_seq
pwd_samfile_reads_handle.close()


for each_mp in MappingRecord_dict:
    current_mp_record = MappingRecord_dict[each_mp]
    current_mp_r1_seq  = current_mp_record.r1_seq
    current_mp_r2_seq  = current_mp_record.r2_seq
    current_mp_r1_refs = current_mp_record.r1_refs
    current_mp_r2_refs = current_mp_record.r2_refs

    # only r1 mapped
    if (current_mp_r1_refs != {}) and (current_mp_r2_refs == {}):

        for matched_ref in current_mp_r1_refs:

            matched_ref_no_pos = matched_ref.split('_pos_')[0]
            current_match_cigar = current_mp_r1_refs[matched_ref]
            current_match_cigar_split = cigar_splitter(current_match_cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(current_match_cigar_split)

            # consider as paired
            current_mp_record.qualified_reads = True
            current_mp_record.consider_r1_unmapped_mate = True
            current_mp_record.r1_filtered_refs.add(matched_ref_no_pos)

            # consider as clipping mapped
            if clipping_len >= min_clp_len:

                if (current_match_cigar_split[0][-1] in ['S', 's']) and (current_match_cigar_split[-1][-1] not in ['S', 's']):
                    current_mp_record.qualified_reads = True
                    current_mp_record.consider_r1_clipping_part = True
                    current_clipping_part_len = int(current_match_cigar_split[0][:-1])
                    if current_clipping_part_len > len(current_mp_record.r1_clipping_seq):
                        current_clipping_seq = current_mp_r1_seq[:int(current_match_cigar_split[0][:-1])]
                        current_mp_record.r1_clipping_seq = current_clipping_seq

                elif (current_match_cigar_split[0][-1] not in ['S', 's']) and (current_match_cigar_split[-1][-1] in ['S', 's']):
                    current_mp_record.qualified_reads = True
                    current_mp_record.consider_r1_clipping_part = True
                    current_clipping_part_len = int(current_match_cigar_split[-1][:-1])
                    if current_clipping_part_len > len(current_mp_record.r1_clipping_seq):
                        current_clipping_seq = current_mp_r1_seq[-int(current_match_cigar_split[-1][:-1]):]
                        current_mp_record.r1_clipping_seq = current_clipping_seq

                elif (current_match_cigar_split[0][-1] in ['S', 's']) and (current_match_cigar_split[-1][-1] in ['S', 's']):

                    clipping_len_l = int(current_match_cigar_split[0][:-1])
                    clipping_len_r = int(current_match_cigar_split[-1][:-1])
                    clipping_seq_l = current_mp_r1_seq[:int(current_match_cigar_split[0][:-1])]
                    clipping_seq_r = current_mp_r1_seq[-int(current_match_cigar_split[-1][:-1]):]

                    # keep left ['75S', '125=', '1X', '7=', '1X', '7=', '3S']
                    if (clipping_len_l >= min_clp_len) and (clipping_len_r <= 3):
                        current_mp_record.qualified_reads = True
                        current_mp_record.consider_r1_clipping_part = True
                        if clipping_len_l > len(current_mp_record.r1_clipping_seq):
                            current_mp_record.r1_clipping_seq = clipping_seq_l

                    # keep right ['2S', '125=', '1X', '7=', '1X', '7=', '56S']
                    if (clipping_len_l <= 3) and (clipping_len_r >= min_clp_len):
                        current_mp_record.qualified_reads = True
                        current_mp_record.consider_r1_clipping_part = True
                        if clipping_len_r > len(current_mp_record.r1_clipping_seq):
                            current_mp_record.r1_clipping_seq = clipping_seq_r

    # only r2 mapped
    elif (current_mp_r1_refs == {}) and (current_mp_r2_refs != {}):

        for matched_ref in current_mp_r2_refs:

            matched_ref_no_pos = matched_ref.split('_pos_')[0]
            current_match_cigar = current_mp_r2_refs[matched_ref]
            current_match_cigar_split = cigar_splitter(current_match_cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(current_match_cigar_split)

            # consider as paired
            current_mp_record.qualified_reads = True
            current_mp_record.consider_r2_unmapped_mate = True
            current_mp_record.r2_filtered_refs.add(matched_ref_no_pos)

            # consider as clipping mapped
            if clipping_len >= min_clp_len:

                if (current_match_cigar_split[0][-1] in ['S', 's']) and (current_match_cigar_split[-1][-1] not in ['S', 's']):
                    current_mp_record.qualified_reads = True
                    current_mp_record.consider_r2_clipping_part = True
                    current_clipping_part_len = int(current_match_cigar_split[0][:-1])
                    if current_clipping_part_len > len(current_mp_record.r2_clipping_seq):
                        current_clipping_seq = current_mp_r2_seq[:int(current_match_cigar_split[0][:-1])]
                        current_mp_record.r2_clipping_seq = current_clipping_seq

                elif (current_match_cigar_split[0][-1] not in ['S', 's']) and (current_match_cigar_split[-1][-1] in ['S', 's']):
                    current_mp_record.qualified_reads = True
                    current_mp_record.consider_r2_clipping_part = True
                    current_clipping_part_len = int(current_match_cigar_split[-1][:-1])
                    if current_clipping_part_len > len(current_mp_record.r2_clipping_seq):
                        current_clipping_seq = current_mp_r2_seq[-int(current_match_cigar_split[-1][:-1]):]
                        current_mp_record.r2_clipping_seq = current_clipping_seq

                elif (current_match_cigar_split[0][-1] in ['S', 's']) and (current_match_cigar_split[-1][-1] in ['S', 's']):

                    clipping_len_l = int(current_match_cigar_split[0][:-1])
                    clipping_len_r = int(current_match_cigar_split[-1][:-1])
                    clipping_seq_l = current_mp_r2_seq[:int(current_match_cigar_split[0][:-1])]
                    clipping_seq_r = current_mp_r2_seq[-int(current_match_cigar_split[-1][:-1]):]

                    # keep left ['75S', '125=', '1X', '7=', '1X', '7=', '3S']
                    if (clipping_len_l >= min_clp_len) and (clipping_len_r <= 3):
                        current_mp_record.qualified_reads = True
                        current_mp_record.consider_r2_clipping_part = True
                        if clipping_len_l > len(current_mp_record.r2_clipping_seq):
                            current_mp_record.r2_clipping_seq = clipping_seq_l

                    # keep right ['2S', '125=', '1X', '7=', '1X', '7=', '56S']
                    if (clipping_len_l <= 3) and (clipping_len_r >= min_clp_len):
                        current_mp_record.qualified_reads = True
                        current_mp_record.consider_r2_clipping_part = True
                        if clipping_len_r > len(current_mp_record.r2_clipping_seq):
                            current_mp_record.r2_clipping_seq = clipping_seq_r

    # both of r1 and r2 mapped
    else:
        clipping_mapped_r1_len = 0
        for matched_r1_ref in current_mp_r1_refs:
            if ('S' in current_mp_r1_refs[matched_r1_ref]) or ('s' in current_mp_r1_refs[matched_r1_ref]):
                for cigmar_element in cigar_splitter(current_mp_r1_refs[matched_r1_ref]):
                    if ('S' in cigmar_element) or ('s' in cigmar_element):
                        current_clipping_len = int(cigmar_element[:-1])
                        if current_clipping_len > clipping_mapped_r1_len:
                            clipping_mapped_r1_len = current_clipping_len

        clipping_mapped_r2_len = 0
        for matched_r2_ref in current_mp_r2_refs:
            if ('S' in current_mp_r2_refs[matched_r2_ref]) or ('s' in current_mp_r2_refs[matched_r2_ref]):
                for cigmar_element in cigar_splitter(current_mp_r2_refs[matched_r2_ref]):
                    if ('S' in cigmar_element) or ('s' in cigmar_element):
                        current_clipping_len = int(cigmar_element[:-1])
                        if current_clipping_len > clipping_mapped_r2_len:
                            clipping_mapped_r2_len = current_clipping_len

        if (clipping_mapped_r1_len > min_clp_len) or (clipping_mapped_r2_len > min_clp_len):

            if (len(current_mp_r1_refs) == 1) and (len(current_mp_r2_refs) == 1):
                r1_ref_with_pos = list(current_mp_r1_refs.keys())[0]
                r2_ref_with_pos = list(current_mp_r1_refs.keys())[0]
                r1_ref_no_pos = r1_ref_with_pos.split('_pos_')[0]
                r2_ref_no_pos = r2_ref_with_pos.split('_pos_')[0]
                r1_cigar = list(current_mp_r1_refs.values())[0]
                r2_cigar = list(current_mp_r2_refs.values())[0]
                r1_cigar_split = cigar_splitter(r1_cigar)
                r2_cigar_split = cigar_splitter(r2_cigar)

                if r1_ref_no_pos == r2_ref_no_pos:

                    current_mp_record.qualified_reads = True
                    longest_clipping_len = max([clipping_mapped_r1_len, clipping_mapped_r2_len])

                    if (r1_cigar_split[0][-1] == 'S') and (int(r1_cigar_split[0][:-1]) == longest_clipping_len):
                        longest_clipping_seq = current_mp_r1_seq[:int(r1_cigar_split[0][:-1])]
                        current_mp_record.consider_r1_clipping_part = True
                        current_mp_record.r1_clipping_seq = longest_clipping_seq
                        current_mp_record.r1_filtered_refs.add(r1_ref_no_pos)

                    elif (r1_cigar_split[-1][-1] == 'S') and (int(r1_cigar_split[-1][:-1]) == longest_clipping_len):
                        longest_clipping_seq = current_mp_r1_seq[-int(r1_cigar_split[-1][:-1]):]
                        current_mp_record.consider_r1_clipping_part = True
                        current_mp_record.r1_clipping_seq = longest_clipping_seq
                        current_mp_record.r1_filtered_refs.add(r1_ref_no_pos)

                    elif (r2_cigar_split[0][-1] == 'S') and (int(r2_cigar_split[0][:-1]) == longest_clipping_len):
                        longest_clipping_seq = current_mp_r2_seq[:int(r2_cigar_split[0][:-1])]
                        current_mp_record.consider_r2_clipping_part = True
                        current_mp_record.r2_clipping_seq = longest_clipping_seq
                        current_mp_record.r2_filtered_refs.add(r2_ref_no_pos)

                    elif (r2_cigar_split[-1][-1] == 'S') and (int(r2_cigar_split[-1][:-1]) == longest_clipping_len):
                        longest_clipping_seq = current_mp_r2_seq[-int(r2_cigar_split[-1][:-1]):]
                        current_mp_record.consider_r2_clipping_part = True
                        current_mp_record.r2_clipping_seq = longest_clipping_seq
                        current_mp_record.r2_filtered_refs.add(r2_ref_no_pos)
                else:
                    pass  # ignore

            # at least one mate mapped to multiple refs
            else:
                r1_ref_list_with_pos = list(current_mp_r1_refs.keys())
                r2_ref_list_with_pos = list(current_mp_r2_refs.keys())
                r1_ref_list_no_pos = {ref1.split('_pos_')[0] for ref1 in r1_ref_list_with_pos}
                r2_ref_list_no_pos = {ref2.split('_pos_')[0] for ref2 in r2_ref_list_with_pos}
                r1_cigar_list = list(current_mp_r1_refs.values())
                r2_cigar_list = list(current_mp_r2_refs.values())

                if r1_ref_list_no_pos == r2_ref_list_no_pos:
                    current_mp_record.qualified_reads = True
                    max_clp, max_clp_location = get_max_clp_and_index(r1_cigar_list, r2_cigar_list)
                    if max_clp_location == 'r1_l':
                        current_mp_record.consider_r1_clipping_part = True
                        current_mp_record.r1_clipping_seq = current_mp_record.r1_seq[:max_clp]
                        current_mp_record.r1_filtered_refs = r1_ref_list_no_pos
                    if max_clp_location == 'r1_r':
                        current_mp_record.consider_r1_clipping_part = True
                        current_mp_record.r1_clipping_seq = current_mp_record.r1_seq[-max_clp:]
                        current_mp_record.r1_filtered_refs = r1_ref_list_no_pos
                    if max_clp_location == 'r2_l':
                        current_mp_record.consider_r2_clipping_part = True
                        current_mp_record.r2_clipping_seq = current_mp_record.r2_seq[:max_clp]
                        current_mp_record.r2_filtered_refs = r2_ref_list_no_pos
                    if max_clp_location == 'r2_r':
                        current_mp_record.consider_r2_clipping_part = True
                        current_mp_record.r2_clipping_seq = current_mp_record.r2_seq[-max_clp:]
                        current_mp_record.r2_filtered_refs = r2_ref_list_no_pos

                # important !!!
                else:
                    current_mp_r1_refs_no_pos = [i.split('_pos_')[0] for i in current_mp_r1_refs]
                    current_mp_r2_refs_no_pos = [i.split('_pos_')[0] for i in current_mp_r2_refs]
                    r1_r2_shared_refs = set(current_mp_r1_refs_no_pos).intersection(current_mp_r2_refs_no_pos)

                    if len(r1_r2_shared_refs) > 0:

                        # get shared refs mean mismatch pct
                        shared_refs_mismatch_list = []
                        r1_ref_to_mismatch_pct_dict = {}
                        for r1_ref in current_mp_r1_refs:
                            r1_ref_no_pos = r1_ref.split('_pos_')[0]
                            r1_ref_cigar = current_mp_r1_refs[r1_ref]
                            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r1_ref_cigar))
                            r1_ref_to_mismatch_pct_dict[r1_ref] = mismatch_pct
                            if r1_ref_no_pos in r1_r2_shared_refs:
                                shared_refs_mismatch_list.append(mismatch_pct)

                        r2_ref_to_mismatch_pct_dict = {}
                        for r2_ref in current_mp_r2_refs:
                            r2_ref_no_pos = r2_ref.split('_pos_')[0]
                            r2_ref_cigar = current_mp_r2_refs[r2_ref]
                            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r2_ref_cigar))
                            r2_ref_to_mismatch_pct_dict[r2_ref] = mismatch_pct
                            if r2_ref_no_pos in r1_r2_shared_refs:
                                shared_refs_mismatch_list.append(mismatch_pct)

                        shared_refs_mean_mismatch_pct = sum(shared_refs_mismatch_list)/len(shared_refs_mismatch_list)

                        current_mp_filtered_refs_r1_r2_no_pos = set()
                        current_mp_r1_refs_filtered = {}
                        for r1_ref in current_mp_r1_refs:
                            r1_ref_no_pos = r1_ref.split('_pos_')[0]
                            if r1_ref_no_pos in r1_r2_shared_refs:
                                current_mp_r1_refs_filtered[r1_ref] = current_mp_r1_refs[r1_ref]
                                current_mp_filtered_refs_r1_r2_no_pos.add(r1_ref_no_pos)
                            else:
                                r1_ref_mismatch = r1_ref_to_mismatch_pct_dict[r1_ref]
                                if r1_ref_mismatch <= shared_refs_mean_mismatch_pct:
                                    current_mp_r1_refs_filtered[r1_ref] = current_mp_r1_refs[r1_ref]
                                    current_mp_filtered_refs_r1_r2_no_pos.add(r1_ref_no_pos)

                        current_mp_r2_refs_filtered = {}
                        for r2_ref in current_mp_r2_refs:
                            r2_ref_no_pos = r2_ref.split('_pos_')[0]
                            if r2_ref_no_pos in r1_r2_shared_refs:
                                current_mp_r2_refs_filtered[r2_ref] = current_mp_r2_refs[r2_ref]
                                current_mp_filtered_refs_r1_r2_no_pos.add(r2_ref_no_pos)
                            else:
                                r2_ref_mismatch = r2_ref_to_mismatch_pct_dict[r2_ref]
                                if r2_ref_mismatch <= shared_refs_mean_mismatch_pct:
                                    current_mp_r2_refs_filtered[r2_ref] = current_mp_r2_refs[r2_ref]
                                    current_mp_filtered_refs_r1_r2_no_pos.add(r2_ref_no_pos)

                        r1_cigar_list_filtered = list(current_mp_r1_refs_filtered.values())
                        r2_cigar_list_filtered = list(current_mp_r2_refs_filtered.values())
                        max_clp, max_clp_location = get_max_clp_and_index(r1_cigar_list_filtered, r2_cigar_list_filtered)
                        if max_clp >= min_clp_len:
                            current_mp_record.qualified_reads = True

                            if max_clp_location == 'r1_l':
                                current_mp_record.consider_r1_clipping_part = True
                                current_mp_record.r1_clipping_seq = current_mp_record.r1_seq[:max_clp]
                                current_mp_record.r1_filtered_refs = current_mp_filtered_refs_r1_r2_no_pos

                            if max_clp_location == 'r1_r':
                                current_mp_record.consider_r1_clipping_part = True
                                current_mp_record.r1_clipping_seq = current_mp_record.r1_seq[-max_clp:]
                                current_mp_record.r1_filtered_refs = current_mp_filtered_refs_r1_r2_no_pos

                            if max_clp_location == 'r2_l':
                                current_mp_record.consider_r2_clipping_part = True
                                current_mp_record.r2_clipping_seq = current_mp_record.r2_seq[:max_clp]
                                current_mp_record.r2_filtered_refs = current_mp_filtered_refs_r1_r2_no_pos

                            if max_clp_location == 'r2_r':
                                current_mp_record.consider_r2_clipping_part = True
                                current_mp_record.r2_clipping_seq = current_mp_record.r2_seq[-max_clp:]
                                current_mp_record.r2_filtered_refs = current_mp_filtered_refs_r1_r2_no_pos
                    else:
                        pass  # worth to look at

        else:
            pass  # to check, most likely to ignore
#     if (MappingRecord_dict[each_mp].qualified_reads is True):
#
#         print('Summary:')
#         print('qualified_reads: %s'                             % current_mp_record.qualified_reads)
#         print('consider_r1_unmapped_mate: %s'                   % current_mp_record.consider_r1_unmapped_mate)
#         print('consider_r1_clipping_part: %s'                   % current_mp_record.consider_r1_clipping_part)
#         print('consider_r2_unmapped_mate: %s'                   % current_mp_record.consider_r2_unmapped_mate)
#         print('consider_r2_clipping_part: %s'                   % current_mp_record.consider_r2_clipping_part)
#         print('r1_clipping_seq: %s'                             % current_mp_record.r1_clipping_seq)
#         print('r2_clipping_seq: %s'                             % current_mp_record.r2_clipping_seq)
#         print('r1_filtered_refs: %s'                            % current_mp_record.r1_filtered_refs)
#         print('r2_filtered_refs: %s'                            % current_mp_record.r2_filtered_refs)
#     print()


# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    if MappingRecord_dict[each_mp].qualified_reads is False:
        MappingRecord_dict.pop(each_mp)



clipping_part_seq_handle = open(clipping_part_seq, 'w')
unmapped_reads_to_extract = set()
for qualified_read in MappingRecord_dict:

    r1_name = '%s.1' % qualified_read
    r2_name = '%s.2' % qualified_read
    read_mr = MappingRecord_dict[qualified_read]

    # print('%s\tconsider_r1_unmapped_mate: %s' % (qualified_read, read_mr.consider_r1_unmapped_mate))
    # print('%s\tconsider_r1_clipping_part: %s' % (qualified_read, read_mr.consider_r1_clipping_part))
    # print('%s\tconsider_r2_unmapped_mate: %s' % (qualified_read, read_mr.consider_r2_unmapped_mate))
    # print('%s\tconsider_r2_clipping_part: %s' % (qualified_read, read_mr.consider_r2_clipping_part))
    # print('%s\tr1_clipping_seq: %s' % (qualified_read, read_mr.r1_clipping_seq))
    # print('%s\tr2_clipping_seq: %s' % (qualified_read, read_mr.r2_clipping_seq))
    # print('%s\tr1_filtered_refs: %s' % (qualified_read, read_mr.r1_filtered_refs))
    # print('%s\tr2_filtered_refs: %s' % (qualified_read, read_mr.r2_filtered_refs))
    # print()

    if read_mr.consider_r1_unmapped_mate is True:
        unmapped_reads_to_extract.add(r2_name)
    if read_mr.consider_r2_unmapped_mate is True:
        unmapped_reads_to_extract.add(r1_name)
    if read_mr.consider_r1_clipping_part is True:
        clipping_part_seq_handle.write('>%s\n' % r1_name)
        clipping_part_seq_handle.write('%s\n'  % read_mr.r1_clipping_seq)
    if read_mr.consider_r2_clipping_part is True:
        clipping_part_seq_handle.write('>%s\n' % r2_name)
        clipping_part_seq_handle.write('%s\n'  % read_mr.r2_clipping_seq)
clipping_part_seq_handle.close()


# extract unmapped_reads
unmapped_mates_handle = open(unmapped_mates, 'w')
for each_read in SeqIO.parse(pwd_samfile_reads, 'fasta'):
    if each_read.id in unmapped_reads_to_extract:
        unmapped_mates_handle.write('>%s\n' % each_read.id)
        unmapped_mates_handle.write('%s\n'  % str(each_read.seq))
unmapped_mates_handle.close()

########################################################################################################################

# mapping with bbmap
pwd_bbmap_exe                       = 'bbmap.sh'
blast_db                            = '/srv/scratch/z5039045/MarkerMAG_wd/BBAY/BBAY_0205_70_3_99.9_MarkerMAG_wd/BBAY_0205_70_3_99.9_step_1_wd/CF_refined_bins_db/CF_refined_bins_combined.fa'
unmapped_mates_file                 = 'unmapped_mates.fa'
clipping_parts_file                 = 'clipping_parts.fa'
pwd_samfile_to_mag_unmapped         = 'unmapped_mates.sam'
pwd_samfile_to_mag_clipping         = 'clipping_parts.sam'
bbmap_parameter                     = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % ('12', '10')
pwd_samfile_to_mag_stderr_unmapped  = 'unmapped_mates_bbmap_stderr.txt'
pwd_samfile_to_mag_stderr_clipping  = 'clipping_parts_bbmap_stderr.txt'
bbmap_index_and_mapping_cmd_to_mag_paired   = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, blast_db, unmapped_mates_file, pwd_samfile_to_mag_unmapped, bbmap_parameter, pwd_samfile_to_mag_stderr_unmapped)
bbmap_index_and_mapping_cmd_to_mag_clipping = '%s ref=%s in=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, blast_db, clipping_parts_file, pwd_samfile_to_mag_clipping, bbmap_parameter, pwd_samfile_to_mag_stderr_clipping)

#print(bbmap_index_and_mapping_cmd_to_mag_paired)
#print(bbmap_index_and_mapping_cmd_to_mag_clipping)

'''
cd /srv/scratch/z5039045/MarkerMAG_wd/BBAY/BBAY_0205_70_3_99.9_MarkerMAG_wd
bbmap.sh ref=/srv/scratch/z5039045/MarkerMAG_wd/BBAY/BBAY_0205_70_3_99.9_MarkerMAG_wd/BBAY_0205_70_3_99.9_step_1_wd/CF_refined_bins_db/CF_refined_bins_combined.fa in=unmapped_mates.fa outm=unmapped_mates.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g 2> unmapped_mates_bbmap_stderr.txt
bbmap.sh ref=/srv/scratch/z5039045/MarkerMAG_wd/BBAY/BBAY_0205_70_3_99.9_MarkerMAG_wd/BBAY_0205_70_3_99.9_step_1_wd/CF_refined_bins_db/CF_refined_bins_combined.fa in=clipping_parts.fa outm=clipping_parts.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g 2> clipping_parts_bbmap_stderr.txt

cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_70_3_MarkerMAG_wd
bbmap.sh ref=/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_70_3_MarkerMAG_wd/Kelp_70_3_step_1_wd/BH_ER_050417_refined_bins_db/BH_ER_050417_refined_bins_combined.fa in=unmapped_mates.fa outm=unmapped_mates.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g 2> unmapped_mates_bbmap_stderr.txt
bbmap.sh ref=/srv/scratch/z5039045/MarkerMAG_wd/Kelp/Kelp_70_3_MarkerMAG_wd/Kelp_70_3_step_1_wd/BH_ER_050417_refined_bins_db/BH_ER_050417_refined_bins_combined.fa in=clipping_parts.fa outm=clipping_parts.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g 2> clipping_parts_bbmap_stderr.txt
'''

########################################################################################################################

for each_read in open(unmapped_mates_sam):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_seq = each_read_split[9]
        cigar = each_read_split[5]
        if cigar != '*':
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = each_read_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

            if (aligned_len >= min_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= max_mis_pct):
                if read_strand == '1':
                    MappingRecord_dict[read_id_base].unmapped_r1_refs.add(ref_id)
                if read_strand == '2':
                    MappingRecord_dict[read_id_base].unmapped_r2_refs.add(ref_id)

                # print(each_read.strip())
                #print(cigar_splitted)
                # print('%s\tconsider_r1_unmapped_mate: %s'   % (read_id, MappingRecord_dict[read_id_base].consider_r1_unmapped_mate))
                # print('%s\tconsider_r2_unmapped_mate: %s'   % (read_id, MappingRecord_dict[read_id_base].consider_r2_unmapped_mate))
                # print('%s\tr1_filtered_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].r1_filtered_refs))
                # print('%s\tr2_filtered_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].r2_filtered_refs))
                # print('%s\tunmapped_r1_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].unmapped_r1_refs))
                # print('%s\tunmapped_r2_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].unmapped_r2_refs))
                # print()


for each_read in open(clipping_part_seq_sam):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_seq = each_read_split[9]
        cigar = each_read_split[5]
        if cigar != '*':
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = each_read_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

            if (aligned_len >= min_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= max_mis_pct):
                if read_strand == '1':
                    MappingRecord_dict[read_id_base].clipping_r1_refs.add(ref_id)
                if read_strand == '2':
                    MappingRecord_dict[read_id_base].clipping_r2_refs.add(ref_id)

                # print(each_read.strip())
                # print(cigar_splitted)
                # print('%s\tconsider_r1_clipping_part: %s'   % (read_id, MappingRecord_dict[read_id_base].consider_r1_clipping_part))
                # print('%s\tconsider_r2_clipping_part: %s'   % (read_id, MappingRecord_dict[read_id_base].consider_r2_clipping_part))
                # print('%s\tr1_filtered_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].r1_filtered_refs))
                # print('%s\tr2_filtered_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].r2_filtered_refs))
                # print('%s\tclipping_r1_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].clipping_r1_refs))
                # print('%s\tclipping_r2_refs: %s'            % (read_id, MappingRecord_dict[read_id_base].clipping_r2_refs))
                # print()


########################################################################################################################

marker_to_ctg_link_num_dict_pair = {}
marker_to_ctg_link_num_dict_clip = {}
for qualified_read in MappingRecord_dict:

    read_mr = MappingRecord_dict[qualified_read]

    if (len(read_mr.unmapped_r2_refs) > 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) == 0):

        for r1_ref in read_mr.r1_filtered_refs:
            for unmapped_r2_ref in read_mr.unmapped_r2_refs:
                marker_to_ctg_key = '%s%s%s' % (r1_ref, marker_to_ctg_Key_connector_str, unmapped_r2_ref)
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1

    elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) > 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) == 0):

        for r2_ref in read_mr.r2_filtered_refs:
            for unmapped_r1_ref in read_mr.unmapped_r1_refs:
                marker_to_ctg_key = '%s%s%s' % (r2_ref, marker_to_ctg_Key_connector_str, unmapped_r1_ref)
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1

    elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) > 0) and (len(read_mr.clipping_r2_refs) == 0):

        for r1_ref in read_mr.r1_filtered_refs:
            for r1_clipping_ref in read_mr.clipping_r1_refs:
                marker_to_ctg_key = '%s%s%s' % (r1_ref, marker_to_ctg_Key_connector_str, r1_clipping_ref)
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1

    elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) > 0):

        for r2_ref in read_mr.r2_filtered_refs:
            for r2_clipping_ref in read_mr.clipping_r2_refs:
                marker_to_ctg_key = '%s%s%s' % (r2_ref, marker_to_ctg_Key_connector_str, r2_clipping_ref)
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1

    elif (len(read_mr.unmapped_r2_refs) > 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) > 0) and (len(read_mr.clipping_r2_refs) == 0):

        shared_refs = set(read_mr.unmapped_r2_refs).intersection(read_mr.clipping_r1_refs)

        for r1_ref in read_mr.r1_filtered_refs:
            for shared_ref in shared_refs:
                marker_to_ctg_key = '%s%s%s' % (r1_ref, marker_to_ctg_Key_connector_str, shared_ref)

                # add to marker_to_ctg_link_num_dict_pair
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1

                # add to marker_to_ctg_link_num_dict_clip
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1

    elif (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) > 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) > 0):

        shared_refs = set(read_mr.unmapped_r1_refs).intersection(read_mr.clipping_r2_refs)

        for r2_ref in read_mr.r2_filtered_refs:
            for shared_ref in shared_refs:
                marker_to_ctg_key = '%s%s%s' % (r2_ref, marker_to_ctg_Key_connector_str, shared_ref)

                # add to marker_to_ctg_link_num_dict_pair
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_pair:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_pair[marker_to_ctg_key] += 1

                # add to marker_to_ctg_link_num_dict_clip
                if marker_to_ctg_key not in marker_to_ctg_link_num_dict_clip:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] = 1
                else:
                    marker_to_ctg_link_num_dict_clip[marker_to_ctg_key] += 1

    # elif (len(read_mr.unmapped_r2_refs) > 0) or (len(read_mr.unmapped_r1_refs) > 0) or (len(read_mr.clipping_r1_refs) > 0) or (len(read_mr.clipping_r2_refs) > 0):
    #
    #     print('\n%s\tr1_filtered_refs:\t%s'         % (qualified_read, read_mr.r1_filtered_refs))
    #     print('%s\tr2_filtered_refs:\t%s\n'         % (qualified_read, read_mr.r2_filtered_refs))
    #     print('%s\tconsider_r1_unmapped_mate:\t%s'  % (qualified_read, read_mr.consider_r1_unmapped_mate))
    #     print('%s\tunmapped_r2_refs:\t%s\n'         % (qualified_read, read_mr.unmapped_r2_refs))
    #     print('%s\tconsider_r2_unmapped_mate:\t%s'  % (qualified_read, read_mr.consider_r2_unmapped_mate))
    #     print('%s\tunmapped_r1_refs:\t%s\n'         % (qualified_read, read_mr.unmapped_r1_refs))
    #     print('%s\tconsider_r1_clipping_part:\t%s'  % (qualified_read, read_mr.consider_r1_clipping_part))
    #     print('%s\tclipping_r1_refs:\t%s\n'         % (qualified_read, read_mr.clipping_r1_refs))
    #     print('%s\tconsider_r2_clipping_part:\t%s'  % (qualified_read, read_mr.consider_r2_clipping_part))
    #     print('%s\tclipping_r2_refs:\t%s\n'         % (qualified_read, read_mr.clipping_r2_refs))
    #     print('---------------\n')


# combine clip and pair dict
marker_to_ctg_link_num_combined = {}
for each_marker_to_ctg_key in marker_to_ctg_link_num_dict_pair:
    marker_to_ctg_link_num_combined[each_marker_to_ctg_key] = marker_to_ctg_link_num_dict_pair[each_marker_to_ctg_key]
for each_marker_to_ctg_key in marker_to_ctg_link_num_dict_clip:
    if each_marker_to_ctg_key not in marker_to_ctg_link_num_combined:
        marker_to_ctg_link_num_combined[each_marker_to_ctg_key] = marker_to_ctg_link_num_dict_clip[each_marker_to_ctg_key]
    else:
        marker_to_ctg_link_num_combined[each_marker_to_ctg_key] += marker_to_ctg_link_num_dict_clip[each_marker_to_ctg_key]


marker_to_gnm_link_num_combined = {}
for each_marker_to_ctg_key in marker_to_ctg_link_num_combined:

    marker_id = each_marker_to_ctg_key.split(marker_to_ctg_Key_connector_str)[0]
    ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_Key_connector_str)[1]
    gnm_id = ctg_id.split(gnm_ctg_connector)[0]
    marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_gnm_Key_connector_str, gnm_id)

    if marker_to_gnm_key not in marker_to_gnm_link_num_combined:
        marker_to_gnm_link_num_combined[marker_to_gnm_key] = marker_to_ctg_link_num_combined[each_marker_to_ctg_key]
    else:
        marker_to_gnm_link_num_combined[marker_to_gnm_key] += marker_to_ctg_link_num_combined[each_marker_to_ctg_key]


sankey_file_in_clipping_handle = open(stats_combined, 'w')
sankey_file_in_clipping_handle.write('MarkerGene,GenomicSeq,Number\n')
for each_clipping in marker_to_gnm_link_num_combined:
    sankey_file_in_clipping_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (each_clipping.split(marker_to_gnm_Key_connector_str)[0], each_clipping.split(marker_to_gnm_Key_connector_str)[1], marker_to_gnm_link_num_combined[each_clipping]))
sankey_file_in_clipping_handle.close()


pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, 500, 30)

filter_linkages_iteratively(stats_combined, 'Number', pairwise_16s_iden_dict, {}, {}, 0, 98, 10, 5, stats_combined_filtered)

mock_final_op_handle = open(mock_final_op, 'w')
mock_final_op_handle.write('MarkerGene	GenomicSeq	Linkage	Step\n')
for each_line in open(stats_combined_filtered):
    if not each_line.startswith('MarkerGene,GenomicSeq,Number'):
        each_line_split = each_line.strip().split(',')
        mock_final_op_handle.write('%s\t%s\t%s\tS1\n' % (each_line_split[0][12:], each_line_split[1][12:], each_line_split[2]))
mock_final_op_handle.close()

