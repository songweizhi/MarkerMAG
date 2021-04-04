import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc


def sam_flag_to_rc(flag_value):

    read_rced = 'na'
    if flag_value != '':
        binary_flag = "{0:b}".format(int(flag_value))
        binary_flag_len = len(str(binary_flag))
        binary_flag_polished = '0' * (12 - binary_flag_len) + str(binary_flag)

        if binary_flag_polished[7] == '0':
            read_rced = False
        if binary_flag_polished[7] == '1':
            read_rced = True

    return read_rced


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

    # get the best cigar
    best_cigar = ''
    if max_value_index == 'r1_l':
        if best_cigar == '':
            for each_cigar in r1_cigar_list:
                if (each_cigar.startswith('%sS' % max_value)) or (each_cigar.startswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r1_r':
        if best_cigar == '':
            for each_cigar in r1_cigar_list:
                if (each_cigar.endswith('%sS' % max_value)) or (each_cigar.endswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r2_l':
        if best_cigar == '':
            for each_cigar in r2_cigar_list:
                if (each_cigar.startswith('%sS' % max_value)) or (each_cigar.startswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r2_r':
        if best_cigar == '':
            for each_cigar in r2_cigar_list:
                if (each_cigar.endswith('%sS' % max_value)) or (each_cigar.endswith('%ss' % max_value)):
                    best_cigar = each_cigar

    return best_cigar, max_value, max_value_index


class MappingRecord:

    #  sequences store in r1_seq and r2_seq should NOT been reverse complemented

    def __init__(self):

        self.r1_seq = ''
        self.r2_seq = ''
        self.r1_seq_qual = '*'
        self.r2_seq_qual = '*'

        self.r1_refs = dict()
        self.r2_refs = dict()

        self.r1_cigar_to_flag = dict()
        self.r2_cigar_to_flag = dict()

        self.r1_longest_clp_cigar = ''
        self.r1_longest_clp_falg  = ''
        self.r2_longest_clp_cigar = ''
        self.r2_longest_clp_falg  = ''

        self.qualified_reads           = False
        self.consider_round_2          = False
        self.consider_r1_unmapped_mate = False
        self.consider_r1_clipping_part = False
        self.consider_r2_unmapped_mate = False
        self.consider_r2_clipping_part = False

        self.r1_filtered_refs = set()
        self.r2_filtered_refs = set()

        self.r1_clipping_seq = ''
        self.r2_clipping_seq = ''

        self.r1_clipping_seq_qual = '*'
        self.r2_clipping_seq_qual = '*'

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
min_clp_len                         = 30
min_clp_M_len                       = 25
max_mis_pct                         = 3
min_clp_len_round2                  = 30
end_seq_len                         = 1000

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
mock_final_op_ctg_level             = '/Users/songweizhi/Desktop/new_algorithm_Kelp/mock_final_op_ctg_level.txt'

marker_to_ctg_Key_connector_str = '___M___'
marker_to_gnm_Key_connector_str = '___M___'
gnm_ctg_connector = '___'


########################################################################################################################
print('read in sam file')

print('rparse mapping results')
MappingRecord_dict = {}
for each_read in open(pwd_samfile):
    if not each_read.startswith('@'):
        store_read_seq = False
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        read_flag = int(each_read_split[1])
        read_seq = each_read_split[9]
        read_seq_qual = each_read_split[10]
        cigar = each_read_split[5]

        if cigar != '*':
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
                    MappingRecord_dict[read_id_base].r1_cigar_to_flag[cigar] = read_flag
                    store_read_seq = True
                if read_strand == '2':
                    MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar
                    MappingRecord_dict[read_id_base].r2_cigar_to_flag[cigar] = read_flag
                    store_read_seq = True
            else:
                store_read_seq = True
        else:
            store_read_seq = True

        # store_read_seq into dict
        if store_read_seq is True:

            # turn back if read reverse complemented
            read_rc = sam_flag_to_rc(read_flag)
            read_seq_to_store = read_seq
            read_seq_qual_to_store = read_seq_qual
            if read_rc is True:
                read_seq_to_store = get_rc(read_seq)
                read_seq_qual_to_store = read_seq_qual[::-1]

            if read_id_base not in MappingRecord_dict:
                MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if MappingRecord_dict[read_id_base].r1_seq == '':
                    MappingRecord_dict[read_id_base].r1_seq = read_seq_to_store
                    MappingRecord_dict[read_id_base].r1_seq_qual = read_seq_qual_to_store
            if read_strand == '2':
                if MappingRecord_dict[read_id_base].r2_seq == '':
                    MappingRecord_dict[read_id_base].r2_seq = read_seq_to_store
                    MappingRecord_dict[read_id_base].r2_seq_qual = read_seq_qual_to_store


for each_mp in MappingRecord_dict:

    current_mp_record           = MappingRecord_dict[each_mp]
    current_mp_r1_seq           = current_mp_record.r1_seq
    current_mp_r2_seq           = current_mp_record.r2_seq
    current_mp_r1_seq_qual      = current_mp_record.r1_seq_qual
    current_mp_r2_seq_qual      = current_mp_record.r2_seq_qual
    current_mp_r1_refs          = current_mp_record.r1_refs
    current_mp_r2_refs          = current_mp_record.r2_refs
    current_mp_r1_refs_no_pos   = {i.split('_pos_')[0] for i in current_mp_r1_refs}
    current_mp_r2_refs_no_pos   = {i.split('_pos_')[0] for i in current_mp_r2_refs}
    current_mp_all_refs_no_pos  = current_mp_r1_refs_no_pos.union(current_mp_r2_refs_no_pos)
    r1_cigar_list               = list(current_mp_r1_refs.values())
    r2_cigar_list               = list(current_mp_r2_refs.values())
    best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, r2_cigar_list)

    best_cigar_flag = ''
    if max_value_index in ['r1_l', 'r1_r']:
        best_cigar_flag         = current_mp_record.r1_cigar_to_flag.get(best_cigar, '')
    if max_value_index in ['r2_l', 'r2_r']:
        best_cigar_flag = current_mp_record.r2_cigar_to_flag.get(best_cigar, '')

    best_cigar_rc           = sam_flag_to_rc(best_cigar_flag)

    # only r1 mapped
    if (current_mp_r1_refs != {}) and (current_mp_r2_refs == {}):

        # consider as paired
        current_mp_record.qualified_reads = True
        current_mp_record.consider_r1_unmapped_mate = True
        current_mp_record.r1_filtered_refs = current_mp_r1_refs_no_pos

        # consider as clipping mapped
        if max_value >= min_clp_len:

            current_clipping_seq = ''
            current_clipping_seq_qual = '*'
            if best_cigar_rc is False:
                if max_value_index == 'r1_l':
                    current_clipping_seq = current_mp_record.r1_seq[:max_value]
                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r1_seq_qual[:max_value]

                if max_value_index == 'r1_r':
                    current_clipping_seq = current_mp_record.r1_seq[-max_value:]
                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r1_seq_qual[-max_value:]

            if best_cigar_rc is True:
                r1_seq_rc = get_rc(current_mp_record.r1_seq)
                r1_seq_rc_qual = current_mp_record.r1_seq_qual[::-1]

                if max_value_index == 'r1_l':
                    current_clipping_seq_rc = r1_seq_rc[:max_value]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_rc_qual = r1_seq_rc_qual[:max_value]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

                if max_value_index == 'r1_r':
                    current_clipping_seq_rc = r1_seq_rc[-max_value:]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_rc_qual = r1_seq_rc_qual[-max_value:]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

            current_mp_record.consider_r1_clipping_part = True
            current_mp_record.r1_clipping_seq           = current_clipping_seq
            current_mp_record.r1_clipping_seq_qual      = current_clipping_seq_qual

    # only r2 mapped
    elif (current_mp_r1_refs == {}) and (current_mp_r2_refs != {}):

        # consider as paired
        current_mp_record.qualified_reads = True
        current_mp_record.consider_r2_unmapped_mate = True
        current_mp_record.r2_filtered_refs = current_mp_r2_refs_no_pos

        # consider as clipping mapped
        if max_value >= min_clp_len:

            current_clipping_seq = ''
            current_clipping_seq_qual = '*'
            if best_cigar_rc is False:
                if max_value_index == 'r2_l':
                    current_clipping_seq = current_mp_record.r2_seq[:max_value]
                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r2_seq_qual[:max_value]

                if max_value_index == 'r2_r':
                    current_clipping_seq = current_mp_record.r2_seq[-max_value:]
                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r2_seq_qual[-max_value:]

            if best_cigar_rc is True:
                r2_seq_rc = get_rc(current_mp_record.r2_seq)
                r2_seq_rc_qual = current_mp_record.r2_seq_qual[::-1]

                if max_value_index == 'r2_l':
                    current_clipping_seq_rc = r2_seq_rc[:max_value]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_rc_qual = r2_seq_rc_qual[:max_value]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

                if max_value_index == 'r2_r':
                    current_clipping_seq_rc = r2_seq_rc[-max_value:]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_rc_qual = r2_seq_rc_qual[-max_value:]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

            current_mp_record.consider_r2_clipping_part = True
            current_mp_record.r2_clipping_seq           = current_clipping_seq
            current_mp_record.r2_clipping_seq_qual      = current_clipping_seq_qual



    # both of r1 and r2 mapped
    else:
        if max_value >= min_clp_len:

            current_mp_record.qualified_reads = True

            # print('%s r1_seq: %s' % (each_mp, current_mp_record.r1_seq))
            # print('%s r1_seq: %s' % (each_mp, current_mp_record.r1_seq_qual))
            # print('%s r2_seq: %s' % (each_mp, current_mp_record.r2_seq))
            # print('%s r2_seq: %s' % (each_mp, current_mp_record.r2_seq_qual))
            # print('%s r1_refs: %s' % (each_mp, current_mp_record.r1_refs))
            # print('%s r2_refs: %s' % (each_mp, current_mp_record.r2_refs))
            # #print('%s all_refs: %s' % (each_mp, current_mp_all_refs_no_pos))
            # #print('%s r1_cigar_to_flag: %s' % (each_mp, current_mp_record.r1_cigar_to_flag))
            # #print('%s r2_cigar_to_flag: %s' % (each_mp, current_mp_record.r2_cigar_to_flag))
            # #print('r1_cigar_list: %s' % r1_cigar_list)
            # #print('r2_cigar_list: %s' % r2_cigar_list)
            # print('best_cigar: %s' % best_cigar)
            # #print('best_cigar_flag: %s' % best_cigar_flag)
            # print('best_cigar_rc: %s' % best_cigar_rc)
            # print('best_cigar_len: %s' % max_value)
            # print('best_cigar_loc: %s' % max_value_index)

            # for clipping part, the refs from its mate are also considered !!!
            current_clipping_seq = ''
            current_clipping_seq_qual = '*'
            if best_cigar_rc is False:

                if max_value_index == 'r1_l':
                    current_clipping_seq = current_mp_record.r1_seq[:max_value]
                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r1_seq_qual[:max_value]

                if max_value_index == 'r1_r':
                    current_clipping_seq = current_mp_record.r1_seq[-max_value:]
                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r1_seq_qual[-max_value:]

                if max_value_index == 'r2_l':
                    current_clipping_seq = current_mp_record.r2_seq[:max_value]
                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r2_seq_qual[:max_value]

                if max_value_index == 'r2_r':
                    current_clipping_seq = current_mp_record.r2_seq[-max_value:]
                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_qual = current_mp_record.r2_seq_qual[-max_value:]

            if best_cigar_rc is True:

                r1_seq_rc       = get_rc(current_mp_record.r1_seq)
                r2_seq_rc       = get_rc(current_mp_record.r2_seq)
                r1_seq_rc_qual  = current_mp_record.r1_seq_qual[::-1]
                r2_seq_rc_qual  = current_mp_record.r2_seq_qual[::-1]

                if max_value_index == 'r1_l':
                    current_clipping_seq_rc = r1_seq_rc[:max_value]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_rc_qual = r1_seq_rc_qual[:max_value]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

                if max_value_index == 'r1_r':
                    current_clipping_seq_rc = r1_seq_rc[-max_value:]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r1_seq_qual != '*':
                        current_clipping_seq_rc_qual = r1_seq_rc_qual[-max_value:]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

                if max_value_index == 'r2_l':
                    current_clipping_seq_rc = r2_seq_rc[:max_value]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_rc_qual = r2_seq_rc_qual[:max_value]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

                if max_value_index == 'r2_r':
                    current_clipping_seq_rc = r2_seq_rc[-max_value:]
                    current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if current_mp_record.r2_seq_qual != '*':
                        current_clipping_seq_rc_qual = r2_seq_rc_qual[-max_value:]
                        current_clipping_seq_qual = current_clipping_seq_rc_qual[::-1]

            if max_value_index in ['r1_l', 'r1_r']:
                current_mp_record.consider_r1_clipping_part = True
                current_mp_record.r1_clipping_seq = current_clipping_seq
                current_mp_record.r1_clipping_seq_qual = current_clipping_seq_qual
                current_mp_record.r1_filtered_refs = current_mp_all_refs_no_pos
            if max_value_index in ['r2_l', 'r2_r']:
                current_mp_record.consider_r2_clipping_part = True
                current_mp_record.r2_clipping_seq = current_clipping_seq
                current_mp_record.r2_clipping_seq_qual = current_clipping_seq_qual
                current_mp_record.r2_filtered_refs = current_mp_all_refs_no_pos

            # print('r1_clipping_seq: %s' % current_mp_record.r1_clipping_seq)
            # print('r1_clipping_seq: %s' % current_mp_record.r1_clipping_seq_qual)
            # print('r2_clipping_seq: %s' % current_mp_record.r2_clipping_seq)
            # print('r2_clipping_seq: %s' % current_mp_record.r2_clipping_seq_qual)
            # print('consider_r1_unmapped_mate: %s' % current_mp_record.consider_r1_unmapped_mate)
            # print('consider_r1_clipping_part: %s' % current_mp_record.consider_r1_clipping_part)
            # print('consider_r2_unmapped_mate: %s' % current_mp_record.consider_r2_unmapped_mate)
            # print('consider_r2_clipping_part: %s' % current_mp_record.consider_r2_clipping_part)
            # print('r1_filtered_refs: %s' % current_mp_record.r1_filtered_refs)
            # print('r2_filtered_refs: %s' % current_mp_record.r2_filtered_refs)
            # print()


# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    if MappingRecord_dict[each_mp].qualified_reads is False:
        MappingRecord_dict.pop(each_mp)


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

            if (aligned_len >= min_clp_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= max_mis_pct):
                if read_strand == '1':
                    MappingRecord_dict[read_id_base].clipping_r1_refs.add(ref_id)
                if read_strand == '2':
                    MappingRecord_dict[read_id_base].clipping_r2_refs.add(ref_id)


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

# get number of linkages at genome level
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

filter_linkages_iteratively(stats_combined, 'Number', pairwise_16s_iden_dict, {}, {}, 0, 98, 5, 3, stats_combined_filtered)


########################################################################################################################

mock_final_op_handle = open(mock_final_op, 'w')
mock_final_op_handle.write('MarkerGene	GenomicSeq	Linkage	Step\n')
for each_line in open(stats_combined_filtered):
    if not each_line.startswith('MarkerGene,GenomicSeq,Number'):
        each_line_split = each_line.strip().split(',')
        mock_final_op_handle.write('%s\t%s\t%s\tS1\n' % (each_line_split[0][12:], each_line_split[1][12:], each_line_split[2]))
mock_final_op_handle.close()


# summarize linkages at contig level
mock_final_op_ctg_level_handle = open(mock_final_op_ctg_level, 'w')
mock_final_op_ctg_level_handle.write('Marker___Genome(total)\tContig\tPaired\tClipping\tStep\n')
for each_linkage in open(mock_final_op):
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
        each_linkage_split = each_linkage.strip().split('\t')
        marker_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        total_link_num = int(each_linkage_split[2])
        link_step = each_linkage_split[3]

        # first go through link num dict by paired reads
        counted_16s_to_ctg_key = set()
        for each_paired_link in marker_to_ctg_link_num_dict_pair:
            paired_link_16s_id = each_paired_link.split(marker_to_ctg_Key_connector_str)[0]
            paired_link_ctg_id = each_paired_link.split(marker_to_ctg_Key_connector_str)[1]
            paired_link_ctg_id_no_gnm = paired_link_ctg_id.split(gnm_ctg_connector)[1]
            paired_link_gnm_id = paired_link_ctg_id.split(gnm_ctg_connector)[0]
            if (paired_link_16s_id == marker_id) and (paired_link_gnm_id == mag_id):
                current_pair_link_num = marker_to_ctg_link_num_dict_pair[each_paired_link]
                current_clip_link_num = marker_to_ctg_link_num_dict_clip.get(each_paired_link, 0)
                mock_final_op_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t%s\t%s\n' % (paired_link_16s_id, paired_link_gnm_id, total_link_num, paired_link_ctg_id_no_gnm, current_pair_link_num, current_clip_link_num, link_step))
                counted_16s_to_ctg_key.add(each_paired_link)

        # then go through link num dict by clipping reads
        for each_clip_link in marker_to_ctg_link_num_dict_clip:
            if each_clip_link not in counted_16s_to_ctg_key:
                clip_link_16s_id = each_clip_link.split(marker_to_ctg_Key_connector_str)[0]
                clip_link_ctg_id = each_clip_link.split(marker_to_ctg_Key_connector_str)[1]
                clip_link_ctg_id_no_gnm = clip_link_ctg_id.split(gnm_ctg_connector)[1]
                clip_link_gnm_id = clip_link_ctg_id.split(gnm_ctg_connector)[0]
                if (clip_link_16s_id == marker_id) and (clip_link_gnm_id == mag_id):
                    current_pair_link_num = marker_to_ctg_link_num_dict_pair.get(each_clip_link, 0)
                    current_clip_link_num = marker_to_ctg_link_num_dict_clip[each_clip_link]
                    mock_final_op_ctg_level_handle.write('%s___%s(%s)\t%s\t%s\t%s\t%s\n' % (clip_link_16s_id, clip_link_gnm_id, total_link_num, clip_link_ctg_id_no_gnm, current_pair_link_num, current_clip_link_num, link_step))
                    counted_16s_to_ctg_key.add(each_clip_link)
mock_final_op_ctg_level_handle.close()


########################################################################################################################
####################################################### round 2 ########################################################
########################################################################################################################

#################### extract sequences flanking 16S ends ####################

link_stats_combined_filtered_s1 = stats_combined_filtered
round2_fq   = False

free_living_ext = 'fa'
if round2_fq is True:
    free_living_ext = 'fq'


free_living_16s_R1 = '/Users/songweizhi/Desktop/new_algorithm_Kelp/round2_free_living_16s_R1.%s' % free_living_ext
free_living_16s_R2 = '/Users/songweizhi/Desktop/new_algorithm_Kelp/round2_free_living_16s_R2.%s' % free_living_ext
free_living_16s_UP = '/Users/songweizhi/Desktop/new_algorithm_Kelp/round2_free_living_16s_UP.%s' % free_living_ext

# get linked marker genes and genomic sequences in step 1
linked_marker_gene_set = set()
linked_genomic_seq_set = set()
for each_link in open(link_stats_combined_filtered_s1):
    if not each_link.startswith('MarkerGene,GenomicSeq,Number'):
        each_link_split = each_link.strip().split(',')
        linked_marker_gene_set.add(each_link_split[0][12:])
        linked_genomic_seq_set.add(each_link_split[1][12:])


free_living_16s_R1_handle = open(free_living_16s_R1, 'w')
free_living_16s_R2_handle = open(free_living_16s_R2, 'w')
free_living_16s_UP_handle = open(free_living_16s_UP, 'w')
for qualified_read in MappingRecord_dict:
    read_mr = MappingRecord_dict[qualified_read]

    if (len(read_mr.unmapped_r2_refs) == 0) and (len(read_mr.unmapped_r1_refs) == 0) and (len(read_mr.clipping_r1_refs) == 0) and (len(read_mr.clipping_r2_refs) == 0):

        ref_been_linked = False
        for r1_ref in read_mr.r1_filtered_refs:
            if r1_ref in linked_marker_gene_set:
                ref_been_linked = True
        for r2_ref in read_mr.r2_filtered_refs:
            if r2_ref in linked_marker_gene_set:
                ref_been_linked = True

        if ref_been_linked is False:
            read_mr.consider_round_2 = True

            if (read_mr.consider_r1_unmapped_mate is True) and (read_mr.consider_r1_clipping_part is True):

                if round2_fq is False:
                    # write out R1 fa
                    free_living_16s_R1_handle.write('>%s.1\n' % qualified_read)
                    free_living_16s_R1_handle.write('%s\n' % read_mr.r1_clipping_seq)
                    # write out R2 fa
                    free_living_16s_R2_handle.write('>%s.2\n' % qualified_read)
                    free_living_16s_R2_handle.write('%s\n' % get_rc(read_mr.r2_seq))
                else:
                    # write out R1 fq
                    free_living_16s_R1_handle.write('@%s.1\n' % qualified_read)
                    free_living_16s_R1_handle.write('%s\n' % read_mr.r1_clipping_seq)
                    free_living_16s_R1_handle.write('+\n')
                    free_living_16s_R1_handle.write('%s\n' % read_mr.r1_clipping_seq_qual)
                    # write out R2 fq
                    free_living_16s_R2_handle.write('@%s.2\n' % qualified_read)
                    free_living_16s_R2_handle.write('%s\n' % get_rc(read_mr.r2_seq))
                    free_living_16s_R2_handle.write('+\n')
                    free_living_16s_R2_handle.write('%s\n' % read_mr.r2_seq_qual[::-1])

            elif (read_mr.consider_r2_unmapped_mate is True) and (read_mr.consider_r2_clipping_part is True):

                if round2_fq is False:
                    # write out R1 fa
                    free_living_16s_R1_handle.write('>%s.1\n' % qualified_read)
                    free_living_16s_R1_handle.write('%s\n' % read_mr.r1_seq)
                    # write out R2 fa
                    free_living_16s_R2_handle.write('>%s.2\n' % qualified_read)
                    free_living_16s_R2_handle.write('%s\n' % get_rc(read_mr.r2_clipping_seq))
                else:
                    # write out R1 fq
                    free_living_16s_R1_handle.write('@%s.1\n' % qualified_read)
                    free_living_16s_R1_handle.write('%s\n' % read_mr.r1_seq)
                    free_living_16s_R1_handle.write('+\n')
                    free_living_16s_R1_handle.write('%s\n' % read_mr.r1_seq_qual)
                    # write out R2 fq
                    free_living_16s_R2_handle.write('@%s.2\n' % qualified_read)
                    free_living_16s_R2_handle.write('%s\n' % get_rc(read_mr.r2_clipping_seq))
                    free_living_16s_R2_handle.write('+\n')
                    free_living_16s_R2_handle.write('%s\n' % read_mr.r2_clipping_seq_qual[::-1])

            else:
                if read_mr.consider_r1_unmapped_mate is True:

                    if round2_fq is False:
                        free_living_16s_UP_handle.write('>%s.2\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % get_rc(read_mr.r2_seq))
                    else:
                        free_living_16s_UP_handle.write('@%s.2\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % get_rc(read_mr.r2_seq))
                        free_living_16s_UP_handle.write('+\n')
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r2_seq_qual[::-1])

                elif read_mr.consider_r2_unmapped_mate is True:

                    if round2_fq is False:
                        free_living_16s_UP_handle.write('>%s.1\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_seq)
                    else:
                        free_living_16s_UP_handle.write('@%s.1\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_seq)
                        free_living_16s_UP_handle.write('+\n')
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_seq_qual)

                elif read_mr.consider_r1_clipping_part is True:

                    if round2_fq is False:
                        free_living_16s_UP_handle.write('>%s.1\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_clipping_seq)
                    else:
                        free_living_16s_UP_handle.write('@%s.1\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_clipping_seq)
                        free_living_16s_UP_handle.write('+\n')
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r1_clipping_seq_qual)

                elif read_mr.consider_r2_clipping_part is True:

                    if round2_fq is False:
                        free_living_16s_UP_handle.write('>%s.2\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % get_rc(read_mr.r2_clipping_seq))
                    else:
                        free_living_16s_UP_handle.write('@%s.2\n' % qualified_read)
                        free_living_16s_UP_handle.write('%s\n' % get_rc(read_mr.r2_clipping_seq))
                        free_living_16s_UP_handle.write('+\n')
                        free_living_16s_UP_handle.write('%s\n' % read_mr.r2_clipping_seq_qual[::-1])

free_living_16s_R1_handle.close()
free_living_16s_R2_handle.close()
free_living_16s_UP_handle.close()

#################### extract sequences flanking ctg ends ####################

combined_1st_round_unlinked_mags_sam = '/Users/songweizhi/Desktop/round2/round_1_unlinked_gnm.sam'

print('read in combined_1st_round_unlinked_mags_sam')
# parse sam file
round_2_MappingRecord_dict = {}
for each_line in open(combined_1st_round_unlinked_mags_sam):
    each_line_split = each_line.strip().split('\t')
    if not each_line.startswith('@'):
        store_read_seq = False
        read_id = each_line_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        read_flag = int(each_line_split[1])
        cigar = each_line_split[5]
        read_seq = each_line_split[9]
        read_seq_qual = each_line_split[10]
        if cigar != '*':
            ref_id = each_line_split[2]
            ref_pos = each_line_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
            if mismatch_pct <= max_mis_pct:

                if read_id_base not in round_2_MappingRecord_dict:
                    round_2_MappingRecord_dict[read_id_base] = MappingRecord()

                if read_strand == '1':
                    round_2_MappingRecord_dict[read_id_base].r1_refs[ref_id_with_pos] = cigar
                    round_2_MappingRecord_dict[read_id_base].r1_cigar_to_flag[cigar] = read_flag

                if read_strand == '2':
                    round_2_MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar
                    round_2_MappingRecord_dict[read_id_base].r2_cigar_to_flag[cigar] = read_flag

                if clipping_len >= min_clp_len_round2:
                    store_read_seq = True
            else:
                store_read_seq = True
        else:
            store_read_seq = True

        # store_read_seq into dict
        if store_read_seq is True:

            # turn back if read reverse complemented
            read_rc = sam_flag_to_rc(read_flag)
            read_seq_to_store = read_seq
            read_seq_qual_to_store = read_seq_qual
            if read_rc is True:
                read_seq_to_store = get_rc(read_seq)
                read_seq_qual_to_store = read_seq_qual[::-1]

            if read_id_base not in round_2_MappingRecord_dict:
                round_2_MappingRecord_dict[read_id_base] = MappingRecord()
            if read_strand == '1':
                if round_2_MappingRecord_dict[read_id_base].r1_seq == '':
                    round_2_MappingRecord_dict[read_id_base].r1_seq = read_seq_to_store
                    round_2_MappingRecord_dict[read_id_base].r1_seq_qual = read_seq_qual_to_store
            if read_strand == '2':
                if round_2_MappingRecord_dict[read_id_base].r2_seq == '':
                    round_2_MappingRecord_dict[read_id_base].r2_seq = read_seq_to_store
                    round_2_MappingRecord_dict[read_id_base].r2_seq_qual = read_seq_qual_to_store

# parse round_2_MappingRecord_dict
for read_basename in round_2_MappingRecord_dict.copy():
    read_mr = round_2_MappingRecord_dict[read_basename]
    r1_ref_cigar_list = list(read_mr.r1_refs.values())
    r2_ref_cigar_list = list(read_mr.r2_refs.values())
    best_cigar, max_clp, max_clp_location = get_max_clp_and_index(r1_ref_cigar_list, r2_ref_cigar_list)

    if (read_mr.r1_refs == {}) and (read_mr.r2_refs != {}):

        if len(read_mr.r2_refs) == 1:

            r2_ref_cigar = r2_ref_cigar_list[0]
            r2_ref_flag = read_mr.r2_cigar_to_flag[r2_ref_cigar]
            r2_ref_cigar_rc = sam_flag_to_rc(r2_ref_flag)

            # consider the unmapped mate only
            if max_clp < min_clp_len_round2:
                read_mr.consider_round_2 = True
                read_mr.consider_r2_unmapped_mate = True

            # consider both of unmapped mate and clipping part
            else:
                r2_ref_no_pos = list(read_mr.r2_refs.keys())[0].split('_pos_')[0]
                r2_ref_pos = int(list(read_mr.r2_refs.keys())[0].split('_pos_')[1])

                if (r2_ref_no_pos[-1]) == (max_clp_location[-1]) == 'l':
                    read_mr.consider_round_2 = True
                    read_mr.consider_r2_unmapped_mate = True
                    read_mr.consider_r2_clipping_part = True

                    if r2_ref_cigar_rc is False:
                        read_mr.r2_clipping_seq = read_mr.r2_seq[:max_clp]
                        if read_mr.r2_seq_qual != '*':
                            read_mr.r2_clipping_seq_qual = read_mr.r2_seq_qual[:max_clp]

                    if r2_ref_cigar_rc is True:
                        r2_seq_rc = get_rc(read_mr.r2_seq)
                        r2_clipping_seq_rc = r2_seq_rc[:max_clp]
                        read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)
                        if read_mr.r2_seq_qual != '*':
                            r2_seq_rc_qual = read_mr.r2_seq_qual[::-1]
                            r2_clipping_seq_rc_qual = r2_seq_rc_qual[:max_clp]
                            read_mr.r2_clipping_seq_qual = r2_clipping_seq_rc_qual[::-1]

                elif (r2_ref_no_pos[-1]) == (max_clp_location[-1]) == 'r':
                    read_mr.consider_round_2 = True
                    read_mr.consider_r2_unmapped_mate = True
                    read_mr.consider_r2_clipping_part = True

                    if r2_ref_cigar_rc is False:
                        read_mr.r2_clipping_seq = read_mr.r2_seq[-max_clp:]
                        if read_mr.r2_seq_qual != '*':
                            read_mr.r2_clipping_seq_qual = read_mr.r2_seq_qual[-max_clp:]

                    if r2_ref_cigar_rc is True:
                        r2_seq_rc = get_rc(read_mr.r2_seq)
                        r2_clipping_seq_rc = r2_seq_rc[-max_clp:]
                        read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)
                        if read_mr.r2_seq_qual != '*':
                            r2_seq_rc_qual = read_mr.r2_seq_qual[::-1]
                            r2_clipping_seq_rc_qual = r2_seq_rc_qual[-max_clp:]
                            read_mr.r2_clipping_seq_qual = r2_clipping_seq_rc_qual[::-1]

                else:  # mapped to unwanted end, ignore
                    round_2_MappingRecord_dict.pop(read_basename)

        else:  # r2 mapped to multiple refs, ignore
            round_2_MappingRecord_dict.pop(read_basename)

    elif (read_mr.r1_refs != {}) and (read_mr.r2_refs == {}):

        if len(read_mr.r1_refs) == 1:

            r1_ref_cigar = r1_ref_cigar_list[0]
            r1_ref_flag = read_mr.r1_cigar_to_flag[r1_ref_cigar]
            r1_ref_cigar_rc = sam_flag_to_rc(r1_ref_flag)

            # consider the unmapped mate only
            if max_clp < min_clp_len_round2:
                read_mr.consider_round_2 = True
                read_mr.consider_r1_unmapped_mate = True

            # consider both of unmapped mate and clipping part
            else:
                r1_ref_no_pos = list(read_mr.r1_refs.keys())[0].split('_pos_')[0]
                r1_ref_pos = int(list(read_mr.r1_refs.keys())[0].split('_pos_')[1])

                if (r1_ref_no_pos[-1]) == (max_clp_location[-1]) == 'l':
                    read_mr.consider_round_2 = True
                    read_mr.consider_r1_unmapped_mate = True
                    read_mr.consider_r1_clipping_part = True

                    if r1_ref_cigar_rc is False:
                        read_mr.r1_clipping_seq = read_mr.r1_seq[:max_clp]
                        if read_mr.r1_seq_qual != '*':
                            read_mr.r1_clipping_seq_qual = read_mr.r1_seq_qual[:max_clp]

                    if r1_ref_cigar_rc is True:
                        r1_seq_rc = get_rc(read_mr.r1_seq)
                        r1_clipping_seq_rc = r1_seq_rc[:max_clp]
                        read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)
                        if read_mr.r1_seq_qual != '*':
                            r1_seq_rc_qual = read_mr.r1_seq_qual[::-1]
                            r1_clipping_seq_rc_qual = r1_seq_rc_qual[:max_clp]
                            read_mr.r1_clipping_seq_qual = r1_clipping_seq_rc_qual[::-1]

                elif (r1_ref_no_pos[-1]) == (max_clp_location[-1]) == 'r':
                    read_mr.consider_round_2 = True
                    read_mr.consider_r1_unmapped_mate = True
                    read_mr.consider_r1_clipping_part = True

                    if r1_ref_cigar_rc is False:
                        read_mr.r1_clipping_seq = read_mr.r1_seq[-max_clp:]
                        if read_mr.r1_seq_qual != '*':
                            read_mr.r1_clipping_seq_qual = read_mr.r1_seq_qual[-max_clp:]

                    if r1_ref_cigar_rc is True:
                        r1_seq_rc = get_rc(read_mr.r1_seq)
                        r1_clipping_seq_rc = r1_seq_rc[-max_clp:]
                        read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)
                        if read_mr.r1_seq_qual != '*':
                            r1_seq_rc_qual = read_mr.r1_seq_qual[::-1]
                            r1_clipping_seq_rc_qual = r1_seq_rc_qual[-max_clp:]
                            read_mr.r1_clipping_seq_qual = r1_clipping_seq_rc_qual[::-1]

                else:  # mapped to unwanted end, ignore
                    round_2_MappingRecord_dict.pop(read_basename)

        else:  # r1 mapped to multiple refs, ignore
            round_2_MappingRecord_dict.pop(read_basename)

    elif (read_mr.r1_refs != {}) and (read_mr.r2_refs != {}):
        if max_clp >= min_M_len:
            if (len(read_mr.r1_refs) == 1) and (len(read_mr.r2_refs) == 1):
                r1_ref_no_pos = list(read_mr.r1_refs.keys())[0].split('_pos_')[0]
                r2_ref_no_pos = list(read_mr.r2_refs.keys())[0].split('_pos_')[0]
                r1_ref_pos = int(list(read_mr.r1_refs.keys())[0].split('_pos_')[1])
                r2_ref_pos = int(list(read_mr.r2_refs.keys())[0].split('_pos_')[1])

                if r1_ref_no_pos == r2_ref_no_pos:
                    if (r1_ref_no_pos[-1]) == (r2_ref_no_pos[-1]) == (max_clp_location[-1]):

                        if (max_clp_location == 'r1_l') and (r1_ref_pos <= 5):
                            read_mr.consider_round_2 = True
                            read_mr.consider_r1_clipping_part = True
                            best_cigar_flag = read_mr.r1_cigar_to_flag[best_cigar]
                            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                            if best_cigar_rc is False:
                                read_mr.r1_clipping_seq = read_mr.r1_seq[:max_clp]
                                if read_mr.r1_seq_qual != '*':
                                    read_mr.r1_clipping_seq_qual = read_mr.r1_seq_qual[:max_clp]
                            if best_cigar_rc is True:
                                r1_seq_rc = get_rc(read_mr.r1_seq)
                                r1_clipping_seq_rc = r1_seq_rc[:max_clp]
                                read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)
                                if read_mr.r1_seq_qual != '*':
                                    r1_seq_rc_qual = read_mr.r1_seq_qual[::-1]
                                    r1_clipping_seq_rc_qual = r1_seq_rc_qual[:max_clp]
                                    read_mr.r1_clipping_seq_qual = r1_clipping_seq_rc_qual[::-1]

                        elif (max_clp_location == 'r2_l') and (r2_ref_pos <= 5):
                            read_mr.consider_round_2 = True
                            read_mr.consider_r2_clipping_part = True
                            best_cigar_flag = read_mr.r2_cigar_to_flag[best_cigar]
                            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                            if best_cigar_rc is False:
                                read_mr.r2_clipping_seq = read_mr.r2_seq[:max_clp]
                                if read_mr.r2_seq_qual != '*':
                                    read_mr.r2_clipping_seq_qual = read_mr.r2_seq_qual[:max_clp]
                            if best_cigar_rc is True:
                                r2_seq_rc = get_rc(read_mr.r2_seq)
                                r2_clipping_seq_rc = r2_seq_rc[:max_clp]
                                read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)
                                if read_mr.r2_seq_qual != '*':
                                    r2_seq_rc_qual = read_mr.r2_seq_qual[::-1]
                                    r2_clipping_seq_rc_qual = r2_seq_rc_qual[:max_clp]
                                    read_mr.r2_clipping_seq_qual = r2_clipping_seq_rc_qual[::-1]

                        elif (max_clp_location == 'r1_r') and (r1_ref_pos >= (end_seq_len / 2)):
                            read_mr.consider_round_2 = True
                            read_mr.consider_r1_clipping_part = True
                            best_cigar_flag = read_mr.r1_cigar_to_flag[best_cigar]
                            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                            if best_cigar_rc is False:
                                read_mr.r1_clipping_seq = read_mr.r1_seq[-max_clp:]
                                if read_mr.r1_seq_qual != '*':
                                    read_mr.r1_clipping_seq_qual = read_mr.r1_seq_qual[-max_clp:]
                            if best_cigar_rc is True:
                                r1_seq_rc = get_rc(read_mr.r1_seq)
                                r1_clipping_seq_rc = r1_seq_rc[-max_clp:]
                                read_mr.r1_clipping_seq = get_rc(r1_clipping_seq_rc)
                                if read_mr.r1_seq_qual != '*':
                                    r1_seq_rc_qual = read_mr.r1_seq_qual[::-1]
                                    r1_clipping_seq_rc_qual = r1_seq_rc_qual[-max_clp:]
                                    read_mr.r1_clipping_seq_qual = r1_clipping_seq_rc_qual[::-1]

                        elif (max_clp_location == 'r2_r') and (r2_ref_pos >= (end_seq_len / 2)):
                            read_mr.consider_round_2 = True
                            read_mr.consider_r2_clipping_part = True
                            best_cigar_flag = read_mr.r2_cigar_to_flag[best_cigar]
                            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                            if best_cigar_rc is False:
                                read_mr.r2_clipping_seq = read_mr.r2_seq[-max_clp:]
                                if read_mr.r2_seq_qual != '*':
                                    read_mr.r2_clipping_seq_qual = read_mr.r2_seq_qual[-max_clp:]
                            if best_cigar_rc is True:
                                r2_seq_rc = get_rc(read_mr.r2_seq)
                                r2_clipping_seq_rc = r2_seq_rc[-max_clp:]
                                read_mr.r2_clipping_seq = get_rc(r2_clipping_seq_rc)
                                if read_mr.r2_seq_qual != '*':
                                    r2_seq_rc_qual = read_mr.r2_seq_qual[::-1]
                                    r2_clipping_seq_rc_qual = r2_seq_rc_qual[-max_clp:]
                                    read_mr.r2_clipping_seq_qual = r2_clipping_seq_rc_qual[::-1]

                        else:  # not too many of them, ignore now, maybe worth check later,
                            round_2_MappingRecord_dict.pop(read_basename)
                            # print('%s\tr1_refs:\t%s\t%s' % (read_basename, read_mr.r1_refs, read_mr.r1_seq))
                            # print('%s\tr2_refs:\t%s\t%s' % (read_basename, read_mr.r2_refs, read_mr.r2_seq))

                    else:  # mapped to unwanted end, ignore
                        round_2_MappingRecord_dict.pop(read_basename)

                else:  # r1 and r2 mapped to different refs
                    # not too many of them, ignore now
                    round_2_MappingRecord_dict.pop(read_basename)

            else:  # r1 or r2 mapped to multiple refs
                # not too many of them, ignore now
                round_2_MappingRecord_dict.pop(read_basename)

        else:  # ignore and remove element from dict
            round_2_MappingRecord_dict.pop(read_basename)

    else:  # ignore and remove element from dict
        round_2_MappingRecord_dict.pop(read_basename)


########################################################################################################################

print('read in sam_file_mini_assembly_combined')

sam_file_mini_assembly_combined     = '/Users/songweizhi/Desktop/round2/scaffolds_combined.sam'

# remove reads mapped to multiple miniassembly? check later
gap_seq_to_reads_dict = {}
for each_read in open(sam_file_mini_assembly_combined):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            ref_id = each_read_split[2]
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
            if (aligned_len >= min_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= (max_mis_pct)):
                if ref_id not in gap_seq_to_reads_dict:
                    gap_seq_to_reads_dict[ref_id] = [read_id]
                else:
                    gap_seq_to_reads_dict[ref_id].append(read_id)


for gap_seq in gap_seq_to_reads_dict:
    gap_seq_mapped_reads = gap_seq_to_reads_dict[gap_seq]
    gap_seq_mapped_reads_linked_to_16s = {}
    gap_seq_mapped_reads_linked_to_ctg = {}

    print('%s\t%s' % (gap_seq, gap_seq_mapped_reads))

    for mapped_read in gap_seq_mapped_reads:
        mapped_read_basename = mapped_read[:-2]
        mapped_read_strand = mapped_read[-1]
        print(mapped_read)
        print(mapped_read_basename)
        print(mapped_read_strand)

    print()



