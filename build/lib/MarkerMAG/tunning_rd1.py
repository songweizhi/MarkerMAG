import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime
from Bio.SeqRecord import SeqRecord


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


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


class MappingRecord:

    def __init__(self):

        #################### overall ####################

        self.qualified_reads = False

        #################### round 1 16s ####################

        self.consider_r1_unmapped_mate = False
        self.consider_r2_unmapped_mate = False

        self.r1_16s_ref_dict = dict()
        self.r2_16s_ref_dict = dict()

        self.r1_16s_refs_lowest_mismatch = None
        self.r2_16s_refs_lowest_mismatch = None

        self.r1_16s_refs_no_ignored = dict()
        self.r2_16s_refs_no_ignored = dict()
        self.shared_16s_refs_no_ignored = dict()

        self.both_mapped_to_16s = False

        #################### round 1 ctg ####################

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()

        self.r1_ctg_refs_lowest_mismatch = None
        self.r2_ctg_refs_lowest_mismatch = None

        self.r1_ctg_refs_no_ignored = dict()
        self.r2_ctg_refs_no_ignored = dict()
        self.shared_ctg_refs_no_ignored = dict()

        self.matched_to_ctg = False

        #################### round 2 ####################

        self.qualified_reads_rd2 = False

        self.r1_ctg_ref_dict_rd2 = dict()
        self.r2_ctg_ref_dict_rd2 = dict()

        self.r1_ctg_refs_lowest_mismatch_rd2 = None
        self.r2_ctg_refs_lowest_mismatch_rd2 = None

        #################### round 2 mini_assembly ####################

        self.r1_mini_ref_dict = dict()
        self.r2_mini_ref_dict = dict()

        self.r1_mini_refs_lowest_mismatch = None
        self.r2_mini_refs_lowest_mismatch = None

        self.r1_mini_refs_no_ignored = dict()
        self.r2_mini_refs_no_ignored = dict()
        self.shared_mini_refs_no_ignored = dict()


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


def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc


def get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len):

    mismatch_set_all_cigar = set()
    mismatch_set_long_M_cigars = set()
    for each_cigar in r1_ref_cigar_set:
        aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
        mismatch_set_all_cigar.add(mismatch_pct)
        if aligned_len >= min_M_len:
            mismatch_set_long_M_cigars.add(mismatch_pct)

    min_mismatch = 'NA'
    if len(mismatch_set_all_cigar) > 0:
        min_mismatch = min(mismatch_set_all_cigar)
        if len(mismatch_set_long_M_cigars) > 0:
            min_mismatch = min(mismatch_set_long_M_cigars)

    return min_mismatch


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


def check_cigar_quality(cigar_str, mismatch_cutoff, min_M_len, ref_pos, ref_len):
    r2_ctg_ref_cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        r2_aligned_len_ctg, r2_aligned_pct_ctg, r2_clipping_len_ctg, r2_clipping_pct_ctg, r2_mismatch_pct_ctg = get_cigar_stats(
            r2_ctg_ref_cigar_splitted)
        if r2_mismatch_pct_ctg <= mismatch_cutoff:
            # check aligned length
            if r2_aligned_len_ctg >= min_M_len:
                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar_str) or ('s' in cigar_str):
                    clip_in_middle = True
                    if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + r2_aligned_len_ctg - 1) == ref_len:
                            clip_in_middle = False

                if clip_in_middle is False:
                    qualified_cigar = True

    return qualified_cigar


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def sort_csv_by_col(file_in, file_out, col_header):
    df_in        = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


def filter_linkages_iteratively(file_in, sort_by_col_header, pairwise_16s_iden_dict, genomic_seq_depth_dict, marker_gene_depth_dict, min_16s_gnm_multiple, within_genome_16s_divergence_cutoff, min_linkages, min_linkages_for_uniq_linked_16s, within_gnm_linkage_num_diff, file_out):

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
    gnm_max_link_num_dict = {}
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


##################################################### file in ####################################################

step_1_wd = '/Users/songweizhi/Desktop/tunning_rd1'
input_16s_polished                          = '%s/file_in/MBARC26_SILVA138_polished.fasta'                          % step_1_wd
combined_input_gnms                         = '%s/file_in/Refined_refined_bins_renamed_combined.fa'                 % step_1_wd
input_reads_to_16s_sam_sorted               = '%s/file_in/MBARC26_0521_input_reads_to_16S_reformatted_sorted.sam'   % step_1_wd
input_reads_to_16s_sam_sorted               = '%s/file_in/MBARC26_0521_input_reads_to_16S_reformatted_sorted_subset.sam'   % step_1_wd
blast_results_all_vs_all_16s                = '%s/file_in/MBARC26_0521_16S_all_vs_all_blastn.tab'                   % step_1_wd
rd1_extracted_to_gnm_sam_reformatted_sorted = '%s/file_in/rd1_extracted_to_gnm_reformatted_sorted.sam'              % step_1_wd
min_M_len_16s                               = 45
min_M_len_ctg                               = 45
min_M_pct                                   = 35
mismatch_cutoff                             = 2
read_to_marker_connector                    = '___r___'
marker_to_ctg_gnm_Key_connector             = '___M___'
gnm_to_ctg_connector                        = '___C___'
min_aln_16s = 500
min_cov_16s = 30
mean_depth_dict_gnm = {}
mean_depth_dict_16s = {}
min_16s_gnm_multiple = 0
min_iden_16s = 98
min_link_num = 8
within_gnm_linkage_num_diff = 80
time_format = '[%Y-%m-%d %H:%M:%S]'

rd1_clp_pct_diff_txt                        = '%s/rd1_clp_pct_diff.txt'             % step_1_wd
rd1_clp_pct_diff_txt_to_ignore              = '%s/rd1_clp_pct_diff_to_ignore.txt'   % step_1_wd
link_stats_combined                         = '%s/stats_combined.txt'               % step_1_wd
link_stats_combined_filtered_s1             = '%s/stats_combined_filtered.txt'      % step_1_wd
linking_reads_rd1                           = '%s/linking_reads_rd1.txt'            % step_1_wd


##################################################### read in sam file ####################################################

print('%s %s' % ((datetime.now().strftime(time_format)), 'read in sam file'))

# get marker len dict
marker_len_dict = {}
for each_marker_record in SeqIO.parse(input_16s_polished, 'fasta'):
    marker_len_dict[each_marker_record.id] = len(each_marker_record.seq)

processed_num = 0
MappingRecord_dict = {}
to_extract_read_base = set()
current_read_base = ''
current_read_base_r1_16s_ref_dict = dict()
current_read_base_r2_16s_ref_dict = dict()
with open(input_reads_to_16s_sam_sorted) as input_reads_to_16s_sam_opened:
    for each_read in input_reads_to_16s_sam_opened:
        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            cigar = each_read_split[5]
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])

            if current_read_base == '':
                current_read_base = read_id_base

                if cigar != '*':
                    if read_strand == '1':
                        current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                    if read_strand == '2':
                        current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}

            elif read_id_base == current_read_base:

                if cigar != '*':
                    if read_strand == '1':
                        if ref_id not in current_read_base_r1_16s_ref_dict:
                            current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                        else:
                            current_read_base_r1_16s_ref_dict[ref_id][ref_pos] = cigar
                    if read_strand == '2':
                        if ref_id not in current_read_base_r2_16s_ref_dict:
                            current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}
                        else:
                            current_read_base_r2_16s_ref_dict[ref_id][ref_pos] = cigar
            else:
                ################################### analysis previous read refs ####################################

                ########## get lowest mismatch for r1/r2 16s refs ##########

                # get r1_ref_cigar_set
                r1_ref_cigar_set = set()
                for each_pos_dict in current_read_base_r1_16s_ref_dict.values():
                    each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                    r1_ref_cigar_set.update(each_pos_dict_values)

                # get r2_ref_cigar_set
                r2_ref_cigar_set = set()
                for each_pos_dict in current_read_base_r2_16s_ref_dict.values():
                    each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                    r2_ref_cigar_set.update(each_pos_dict_values)

                r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_16s)
                r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_16s)

                refs_to_ignore = set()

                ########## filter r1 16s refs ##########

                r1_16s_refs_passed_qc = {}
                for r1_16s_ref in current_read_base_r1_16s_ref_dict:
                    r1_matched_pos_dict = current_read_base_r1_16s_ref_dict[r1_16s_ref]

                    # one read need to mapped to one 16S only for one time
                    if len(r1_matched_pos_dict) > 1:
                        refs_to_ignore.add(r1_16s_ref)
                    else:
                        r1_16s_ref_pos = list(r1_matched_pos_dict.keys())[0]
                        r1_16s_ref_cigar = r1_matched_pos_dict[r1_16s_ref_pos]
                        r1_16s_ref_cigar_splitted = cigar_splitter(r1_16s_ref_cigar)

                        # check both end clip
                        both_end_clp = check_both_ends_clipping(r1_16s_ref_cigar_splitted)
                        if both_end_clp is True:
                            refs_to_ignore.add(r1_16s_ref)
                        else:
                            # check mismatch
                            r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                r1_16s_ref_cigar_splitted)
                            if r1_ref_min_mismatch == 'NA':
                                refs_to_ignore.add(r1_16s_ref)
                            elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                                refs_to_ignore.add(r1_16s_ref)
                            else:
                                # check aligned length
                                if r1_aligned_len < min_M_len_16s:
                                    refs_to_ignore.add(r1_16s_ref)
                                else:
                                    # check if clp in the middle
                                    clip_in_middle = False
                                    if ('S' in r1_16s_ref_cigar) or ('s' in r1_16s_ref_cigar):
                                        clip_in_middle = True
                                        if (r1_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                r1_16s_ref_pos == 1):
                                            clip_in_middle = False
                                        if (r1_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                            if (r1_16s_ref_pos + r1_aligned_len - 1) == marker_len_dict[r1_16s_ref]:
                                                clip_in_middle = False

                                    # exclude the ref if clp in the middle is True
                                    if clip_in_middle is True:
                                        refs_to_ignore.add(r1_16s_ref)
                                    else:
                                        r1_16s_refs_passed_qc[r1_16s_ref] = [r1_16s_ref_cigar]

                ########## filter r2 16s refs ##########

                r2_16s_refs_passed_qc = {}
                for r2_16s_ref in current_read_base_r2_16s_ref_dict:
                    r2_matched_pos_dict = current_read_base_r2_16s_ref_dict[r2_16s_ref]

                    # one read need to mapped to one 16S only once
                    if len(r2_matched_pos_dict) > 1:
                        refs_to_ignore.add(r2_16s_ref)
                    else:
                        r2_16s_ref_pos = list(r2_matched_pos_dict.keys())[0]
                        r2_16s_ref_cigar = r2_matched_pos_dict[r2_16s_ref_pos]
                        r2_16s_ref_cigar_splitted = cigar_splitter(r2_16s_ref_cigar)

                        # check both end clip
                        both_end_clp = check_both_ends_clipping(r2_16s_ref_cigar_splitted)
                        if both_end_clp is True:
                            refs_to_ignore.add(r2_16s_ref)
                        else:
                            # check mismatch
                            r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                r2_16s_ref_cigar_splitted)
                            if r2_ref_min_mismatch == 'NA':
                                refs_to_ignore.add(r2_16s_ref)
                            elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                                refs_to_ignore.add(r2_16s_ref)
                            else:
                                # check aligned length
                                if r2_aligned_len < min_M_len_16s:
                                    refs_to_ignore.add(r2_16s_ref)
                                else:
                                    # check if clp in the middle
                                    clip_in_middle = False
                                    if ('S' in r2_16s_ref_cigar) or ('s' in r2_16s_ref_cigar):
                                        clip_in_middle = True
                                        if (r2_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                r2_16s_ref_pos == 1):
                                            clip_in_middle = False
                                        if (r2_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                            if (r2_16s_ref_pos + r2_aligned_len - 1) == marker_len_dict[r2_16s_ref]:
                                                clip_in_middle = False

                                    # exclude the ref if clp in the middle is True
                                    if clip_in_middle is True:
                                        refs_to_ignore.add(r2_16s_ref)
                                    else:
                                        r2_16s_refs_passed_qc[r2_16s_ref] = [r2_16s_ref_cigar]

                #################################

                r1_16s_refs_no_ignored = {key: value for key, value in r1_16s_refs_passed_qc.items() if
                                          key not in refs_to_ignore}
                r2_16s_refs_no_ignored = {key: value for key, value in r2_16s_refs_passed_qc.items() if
                                          key not in refs_to_ignore}

                # no mate has no_ignored alignments
                if (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) == 0):
                    pass

                # only r1 has no_ignored alignments
                elif (len(r1_16s_refs_no_ignored) > 0) and (len(r2_16s_refs_no_ignored) == 0):

                    if current_read_base not in MappingRecord_dict:
                        MappingRecord_dict[current_read_base] = MappingRecord()

                    MappingRecord_dict[current_read_base].qualified_reads = True
                    MappingRecord_dict[current_read_base].consider_r1_unmapped_mate = True
                    MappingRecord_dict[current_read_base].r1_16s_refs_no_ignored = r1_16s_refs_no_ignored
                    MappingRecord_dict[current_read_base].r1_16s_ref_dict = current_read_base_r1_16s_ref_dict
                    to_extract_read_base.add(current_read_base)

                # only r2 has no_ignored alignments
                elif (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) > 0):

                    if current_read_base not in MappingRecord_dict:
                        MappingRecord_dict[current_read_base] = MappingRecord()

                    MappingRecord_dict[current_read_base].qualified_reads = True
                    MappingRecord_dict[current_read_base].consider_r2_unmapped_mate = True
                    MappingRecord_dict[current_read_base].r2_16s_refs_no_ignored = r2_16s_refs_no_ignored
                    MappingRecord_dict[current_read_base].r2_16s_ref_dict = current_read_base_r2_16s_ref_dict

                    to_extract_read_base.add(current_read_base)

                # both r1 and r2 have no_ignored alignments
                else:
                    shared_16s_ref_dict = {key: [r1_16s_refs_no_ignored[key][0], r2_16s_refs_no_ignored[key][0]] for key
                                           in set(r1_16s_refs_no_ignored).intersection(set(r2_16s_refs_no_ignored))}
                    if len(shared_16s_ref_dict) > 0:

                        if current_read_base not in MappingRecord_dict:
                            MappingRecord_dict[current_read_base] = MappingRecord()

                        MappingRecord_dict[current_read_base].qualified_reads = True
                        MappingRecord_dict[current_read_base].both_mapped_to_16s = True
                        MappingRecord_dict[current_read_base].shared_16s_refs_no_ignored = shared_16s_ref_dict
                        MappingRecord_dict[current_read_base].r1_16s_ref_dict = current_read_base_r1_16s_ref_dict
                        MappingRecord_dict[current_read_base].r2_16s_ref_dict = current_read_base_r2_16s_ref_dict
                        to_extract_read_base.add(current_read_base)

                    '''
                    qualified_reads
                    consider_r1_unmapped_mate
                    consider_r2_unmapped_mate
                    both_mapped_to_16s
                    r1_16s_ref_dict
                    r2_16s_ref_dict
                    r1_16s_refs_no_ignored
                    r2_16s_refs_no_ignored
                    shared_16s_refs_no_ignored
                    '''
                # if current_read_base in MappingRecord_dict:
                #     print('%s\tr1_16s_ref_dict\t%s'             % (current_read_base, MappingRecord_dict[current_read_base].r1_16s_ref_dict))
                #     print('%s\tr2_16s_ref_dict\t%s'             % (current_read_base, MappingRecord_dict[current_read_base].r2_16s_ref_dict))
                #     print('%s\tr1_16s_refs_no_ignored\t%s'      % (current_read_base, MappingRecord_dict[current_read_base].r1_16s_refs_no_ignored))
                #     print('%s\tr2_16s_refs_no_ignored\t%s'      % (current_read_base, MappingRecord_dict[current_read_base].r2_16s_refs_no_ignored))
                #     print('%s\tshared_16s_refs_no_ignored\t%s'  % (current_read_base, MappingRecord_dict[current_read_base].shared_16s_refs_no_ignored))
                #     print()
                #
                #     print()


                ########################################### reset values ###########################################

                current_read_base = read_id_base
                current_read_base_r1_16s_ref_dict = dict()
                current_read_base_r2_16s_ref_dict = dict()

                if cigar != '*':
                    if read_strand == '1':
                        current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                    if read_strand == '2':
                        current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}

current_read_base = 'MBARC26_251146'
print('%s\tr1_16s_ref_dict\t%s'             % (current_read_base, MappingRecord_dict[current_read_base].r1_16s_ref_dict))
print('%s\tr2_16s_ref_dict\t%s'             % (current_read_base, MappingRecord_dict[current_read_base].r2_16s_ref_dict))
print('%s\tr1_16s_refs_no_ignored\t%s'      % (current_read_base, MappingRecord_dict[current_read_base].r1_16s_refs_no_ignored))
print('%s\tr2_16s_refs_no_ignored\t%s'      % (current_read_base, MappingRecord_dict[current_read_base].r2_16s_refs_no_ignored))
print('%s\tshared_16s_refs_no_ignored\t%s'  % (current_read_base, MappingRecord_dict[current_read_base].shared_16s_refs_no_ignored))
print()


######################################### read mapping results of rd1 extracted mates into mp dict  #########################################

print('%s %s' % ((datetime.now().strftime(time_format)), 'read in sam file'))

ctg_len_dict = {}
for each_ctg_record in SeqIO.parse(combined_input_gnms, 'fasta'):
    ctg_len_dict[each_ctg_record.id] = len(each_ctg_record.seq)


current_read_base = ''
current_read_base_r1_ctg_ref_dict = dict()
current_read_base_r2_ctg_ref_dict = dict()
with open(rd1_extracted_to_gnm_sam_reformatted_sorted) as rd1_extracted_to_gnm_sam_reformatted_opened:
    for each_read in rd1_extracted_to_gnm_sam_reformatted_opened:
        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            cigar = each_read_split[5]
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])

            if current_read_base == '':
                current_read_base = read_id_base
                if cigar != '*':
                    if read_strand == '1':
                        current_read_base_r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                    if read_strand == '2':
                        current_read_base_r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}

            elif read_id_base == current_read_base:
                if cigar != '*':
                    if read_strand == '1':
                        if ref_id not in current_read_base_r1_ctg_ref_dict:
                            current_read_base_r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                        else:
                            current_read_base_r1_ctg_ref_dict[ref_id][ref_pos] = cigar

                    if read_strand == '2':
                        if ref_id not in current_read_base_r2_ctg_ref_dict:
                            current_read_base_r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                        else:
                            current_read_base_r2_ctg_ref_dict[ref_id][ref_pos] = cigar
            else:
                ################################### analysis previous read refs ####################################

                refs_to_ignore_ctg = set()

                ########## filter r1 ctg refs ##########

                r1_ctg_refs_passed_qc = {}
                for r1_ctg_ref in current_read_base_r1_ctg_ref_dict:
                    r1_ctg_ref_matched_pos_dict = current_read_base_r1_ctg_ref_dict[r1_ctg_ref]

                    # one read need to mapped to one ctg only for one time
                    if len(r1_ctg_ref_matched_pos_dict) > 1:
                        refs_to_ignore_ctg.add(r1_ctg_ref)
                    else:
                        r1_ctg_ref_pos = list(r1_ctg_ref_matched_pos_dict.keys())[0]
                        r1_ctg_ref_cigar = r1_ctg_ref_matched_pos_dict[r1_ctg_ref_pos]
                        r1_ctg_ref_len = ctg_len_dict[r1_ctg_ref]
                        qualified_cigar = check_cigar_quality(r1_ctg_ref_cigar, mismatch_cutoff, min_M_len_ctg,
                                                              r1_ctg_ref_pos, r1_ctg_ref_len)
                        if qualified_cigar is True:
                            r1_ctg_refs_passed_qc[r1_ctg_ref] = [r1_ctg_ref_cigar]
                        else:
                            refs_to_ignore_ctg.add(r1_ctg_ref)

                ########## filter r2 ctg refs ##########

                r2_ctg_refs_passed_qc = {}
                for r2_ctg_ref in current_read_base_r2_ctg_ref_dict:
                    r2_ctg_ref_matched_pos_dict = current_read_base_r2_ctg_ref_dict[r2_ctg_ref]

                    # one read need to mapped to one ctg only for one time
                    if len(r2_ctg_ref_matched_pos_dict) > 1:
                        refs_to_ignore_ctg.add(r2_ctg_ref)
                    else:
                        r2_ctg_ref_pos = list(r2_ctg_ref_matched_pos_dict.keys())[0]
                        r2_ctg_ref_cigar = r2_ctg_ref_matched_pos_dict[r2_ctg_ref_pos]
                        r2_ctg_ref_len = ctg_len_dict.get(r2_ctg_ref, 0)
                        qualified_cigar = check_cigar_quality(r2_ctg_ref_cigar, mismatch_cutoff, min_M_len_ctg,
                                                              r2_ctg_ref_pos, r2_ctg_ref_len)
                        if qualified_cigar is True:
                            r2_ctg_refs_passed_qc[r2_ctg_ref] = [r2_ctg_ref_cigar]
                        else:
                            refs_to_ignore_ctg.add(r2_ctg_ref)

                ####################################################################################################

                r1_ctg_refs_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}
                r2_ctg_refs_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}

                # no mate matched to ctg
                if (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) == 0):
                    pass

                # only r1 matched to ctg
                elif (len(r1_ctg_refs_no_ignored) > 0) and (len(r2_ctg_refs_no_ignored) == 0):

                    if current_read_base not in MappingRecord_dict:
                        MappingRecord_dict[current_read_base] = MappingRecord()

                    MappingRecord_dict[current_read_base].matched_to_ctg = True
                    MappingRecord_dict[current_read_base].r1_ctg_refs_no_ignored = r1_ctg_refs_no_ignored
                    MappingRecord_dict[current_read_base].r1_ctg_ref_dict = current_read_base_r1_ctg_ref_dict
                # only r2 matched to ctg
                elif (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) > 0):

                    if current_read_base not in MappingRecord_dict:
                        MappingRecord_dict[current_read_base] = MappingRecord()

                    MappingRecord_dict[current_read_base].matched_to_ctg = True
                    MappingRecord_dict[current_read_base].r2_ctg_refs_no_ignored = r2_ctg_refs_no_ignored
                    MappingRecord_dict[current_read_base].r2_ctg_ref_dict = current_read_base_r2_ctg_ref_dict

                # both r1 and r2 matched to ctg
                else:
                    shared_ctg_ref_set = {key: [r1_ctg_refs_no_ignored[key][0], r2_ctg_refs_no_ignored[key][0]] for key in set(r1_ctg_refs_no_ignored).intersection(set(r2_ctg_refs_no_ignored))}
                    if len(shared_ctg_ref_set) > 0:

                        if current_read_base not in MappingRecord_dict:
                            MappingRecord_dict[current_read_base] = MappingRecord()

                        MappingRecord_dict[current_read_base].matched_to_ctg = True
                        MappingRecord_dict[current_read_base].shared_ctg_refs_no_ignored = shared_ctg_ref_set
                        MappingRecord_dict[current_read_base].r1_ctg_ref_dict = current_read_base_r1_ctg_ref_dict
                        MappingRecord_dict[current_read_base].r2_ctg_ref_dict = current_read_base_r2_ctg_ref_dict

                ########################################### reset values ###########################################

                current_read_base = read_id_base
                current_read_base_r1_ctg_ref_dict = dict()
                current_read_base_r2_ctg_ref_dict = dict()

                if cigar != '*':
                    if read_strand == '1':
                        current_read_base_r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                    if read_strand == '2':
                        current_read_base_r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}

##################################################### get linkages from MappingRecord_dict #####################################################

print('%s %s' % ((datetime.now().strftime(time_format)), 'get linkages from MappingRecord_dict'))

marker_to_ctg_linkage_cigar_dict_16s_side = {}
marker_to_ctg_linkage_cigar_dict_ctg_side = {}
marker_to_ctg_linkage_num_dict = {}
marker_to_ctg_linking_reads_dict = {}
for qualified_read in MappingRecord_dict:

    r1_16s_refs     = MappingRecord_dict[qualified_read].r1_16s_refs_no_ignored
    r2_16s_refs     = MappingRecord_dict[qualified_read].r2_16s_refs_no_ignored
    shared_16s_refs = MappingRecord_dict[qualified_read].shared_16s_refs_no_ignored
    r1_ctg_refs     = MappingRecord_dict[qualified_read].r1_ctg_refs_no_ignored
    r2_ctg_refs     = MappingRecord_dict[qualified_read].r2_ctg_refs_no_ignored
    shared_ctg_refs = MappingRecord_dict[qualified_read].shared_ctg_refs_no_ignored

    # combine dict
    combined_16s_refs_dict = {**r1_16s_refs, **r2_16s_refs, **shared_16s_refs}
    combined_ctg_refs_dict = {**r1_ctg_refs, **r2_ctg_refs, **shared_ctg_refs}

    if (len(combined_16s_refs_dict) > 0) and (len(combined_ctg_refs_dict) > 0):
        for each_16s_ref in combined_16s_refs_dict:
            current_cigar_list_16s_side = combined_16s_refs_dict[each_16s_ref]
            for each_ctg_ref in combined_ctg_refs_dict:
                current_cigar_list_ctg_side = combined_ctg_refs_dict[each_ctg_ref]

                marker_to_ctg_key = '%s%s%s' % (each_16s_ref, marker_to_ctg_gnm_Key_connector, each_ctg_ref)

                # initialize key values if not in the dict
                if marker_to_ctg_key not in marker_to_ctg_linkage_num_dict:
                    marker_to_ctg_linkage_num_dict[marker_to_ctg_key] = 0
                    marker_to_ctg_linking_reads_dict[marker_to_ctg_key] = set()
                    marker_to_ctg_linkage_cigar_dict_16s_side[marker_to_ctg_key] = []
                    marker_to_ctg_linkage_cigar_dict_ctg_side[marker_to_ctg_key] = []

                # add values to dict
                marker_to_ctg_linkage_num_dict[marker_to_ctg_key] += 1
                marker_to_ctg_linking_reads_dict[marker_to_ctg_key].add(qualified_read)
                for each in current_cigar_list_16s_side:
                    marker_to_ctg_linkage_cigar_dict_16s_side[marker_to_ctg_key].append(each)
                for each in current_cigar_list_ctg_side:
                    marker_to_ctg_linkage_cigar_dict_ctg_side[marker_to_ctg_key].append(each)





print('\nFinally got:')
print(marker_to_ctg_linking_reads_dict['DA_m2___M___DA___C___NODE_537_length_39828_cov_427.956453'])
print(marker_to_ctg_linkage_num_dict['DA_m2___M___DA___C___NODE_537_length_39828_cov_427.956453'])
print(marker_to_ctg_linkage_cigar_dict_16s_side['DA_m2___M___DA___C___NODE_537_length_39828_cov_427.956453'])
print(marker_to_ctg_linkage_cigar_dict_ctg_side['DA_m2___M___DA___C___NODE_537_length_39828_cov_427.956453'])
print()





# calculate clp_pct_ratio
rd1_clp_pct_diff_txt_handle = open(rd1_clp_pct_diff_txt, 'w')
rd1_clp_pct_diff_txt_to_ignore_handle = open(rd1_clp_pct_diff_txt_to_ignore, 'w')
rd1_clp_pct_diff_txt_handle.write('Marker\tGenome\tContig\tclp_pct_ctg\tclp_pct_16s\tRatio\n')
rd1_clp_pct_diff_txt_to_ignore_handle.write('Marker\tGenome\tContig\tclp_pct_ctg\tclp_pct_16s\tRatio\n')
linkages_to_ignore = set()
for each_link in marker_to_ctg_linkage_num_dict:
    if marker_to_ctg_linkage_num_dict[each_link] >= 3:
        each_link_split = each_link.split(marker_to_ctg_gnm_Key_connector)
        marker_id = each_link_split[0]
        gnm_id = each_link_split[1].split(gnm_to_ctg_connector)[0]
        ctg_id = each_link_split[1].split(gnm_to_ctg_connector)[1]

        linkage_cigar_16s_side_all = marker_to_ctg_linkage_cigar_dict_16s_side[each_link]
        linkage_cigar_16s_side_clp = [i for i in linkage_cigar_16s_side_all if (('S' in i) or ('s' in i))]
        linkage_cigar_16s_side_clp_num = len(linkage_cigar_16s_side_clp)

        linkage_cigar_ctg_side_all = marker_to_ctg_linkage_cigar_dict_ctg_side[each_link]
        linkage_cigar_ctg_side_clp = [i for i in linkage_cigar_ctg_side_all if (('S' in i) or ('s' in i))]
        linkage_cigar_ctg_side_clp_num = len(linkage_cigar_ctg_side_clp)

        # if each_link == 'DA_m2___M___DA___C___NODE_537_length_39828_cov_427.956453':
        #     print('%s\tlinkage_cigar_16s_side_all: %s' % (each_link, linkage_cigar_16s_side_all))
        #     print('%s\tlinkage_cigar_16s_side_clp: %s' % (each_link, linkage_cigar_16s_side_clp))
        #     print('%s\tlinkage_cigar_ctg_side_all: %s' % (each_link, linkage_cigar_ctg_side_all))
        #     print('%s\tlinkage_cigar_ctg_side_clp: %s' % (each_link, linkage_cigar_ctg_side_clp))
        #     print()

        mean_clp_num = (linkage_cigar_16s_side_clp_num + linkage_cigar_ctg_side_clp_num)/2

        to_ignore = False
        if (linkage_cigar_16s_side_clp_num == 0) and (linkage_cigar_ctg_side_clp_num >= 10):
            to_ignore = True
        elif (linkage_cigar_16s_side_clp_num >= 10) and (linkage_cigar_ctg_side_clp_num == 0):
            to_ignore = True
        elif (linkage_cigar_16s_side_clp_num > 0) and (linkage_cigar_ctg_side_clp_num > 0):
            if mean_clp_num <= 20:
                if not (0.33 <= (linkage_cigar_16s_side_clp_num/linkage_cigar_ctg_side_clp_num) <= 3):
                    to_ignore = True
            # elif 20 < mean_clp_num <= 50:
            #     if not (0.5<= (linkage_cigar_16s_side_clp_num/linkage_cigar_ctg_side_clp_num) <= 2):
            #         to_ignore = True
            elif mean_clp_num > 20:
                if not (0.5 <= (linkage_cigar_16s_side_clp_num/linkage_cigar_ctg_side_clp_num) <= 2):
                    to_ignore = True

        # clp_pct_16s_side = 'na'
        # if len(linkage_cigar_16s_side_all) > 0:
        #     clp_pct_16s_side = len(linkage_cigar_16s_side_clp) * 100 / len(linkage_cigar_16s_side_all)
        #     clp_pct_16s_side = float("{0:.2f}".format(clp_pct_16s_side))
        #
        # clp_pct_ctg_side = 'na'
        # if len(linkage_cigar_ctg_side_all) > 0:
        #     clp_pct_ctg_side = len(linkage_cigar_ctg_side_clp) * 100 / len(linkage_cigar_ctg_side_all)
        #     clp_pct_ctg_side = float("{0:.2f}".format(clp_pct_ctg_side))
        #
        # clp_pct_ratio = 'na'
        # if (clp_pct_16s_side != 'na') and (clp_pct_ctg_side != 'na'):
        #     if clp_pct_16s_side > 0:
        #         clp_pct_ratio = float("{0:.2f}".format(clp_pct_ctg_side / clp_pct_16s_side))
        #
        # to_ignore = False
        # if clp_pct_ratio != 'na':
        #     if (clp_pct_ctg_side >= clp_pct_ctg_side_max_num) and (clp_pct_ratio >= clp_pct_ratio_cutoff):
        #         to_ignore = True
        # else:
        #     if (clp_pct_ctg_side != 'na') and (clp_pct_16s_side != 'na'):
        #         if clp_pct_ctg_side >= clp_pct_ctg_side_max_num:
        #             to_ignore = True

        if to_ignore is True:
            linkages_to_ignore.add(each_link)
            rd1_clp_pct_diff_txt_to_ignore_handle.write('%s\t%s\t%s\t%s/%s\t%s/%s\n'  % (marker_id, gnm_id, ctg_id, len(linkage_cigar_ctg_side_clp), len(linkage_cigar_ctg_side_all), len(linkage_cigar_16s_side_clp),len(linkage_cigar_16s_side_all)))
        else:
            rd1_clp_pct_diff_txt_handle.write('%s\t%s\t%s\t%s/%s\t%s/%s\n'            % (marker_id, gnm_id, ctg_id, len(linkage_cigar_ctg_side_clp), len(linkage_cigar_ctg_side_all), len(linkage_cigar_16s_side_clp), len(linkage_cigar_16s_side_all)))
rd1_clp_pct_diff_txt_handle.close()
rd1_clp_pct_diff_txt_to_ignore_handle.close()

# ignore linkages if did not pass the clp pct ratio test
marker_to_ctg_linkage_num_dict_min3_passed_ratio_check = {}
for each_link in marker_to_ctg_linkage_num_dict:
    if marker_to_ctg_linkage_num_dict[each_link] >= 3:
        if each_link not in linkages_to_ignore:
            marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_link] = marker_to_ctg_linkage_num_dict[
                each_link]

# get number of linkages at genome level
marker_to_gnm_link_num = {}
for each_marker_to_ctg_key in marker_to_ctg_linkage_num_dict_min3_passed_ratio_check:
    marker_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[0]
    ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[1]
    gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
    marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
    if marker_to_gnm_key not in marker_to_gnm_link_num:
        marker_to_gnm_link_num[marker_to_gnm_key] = marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[
            each_marker_to_ctg_key]
    else:
        marker_to_gnm_link_num[marker_to_gnm_key] += marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[
            each_marker_to_ctg_key]

# write out linkages at genome level
sankey_file_in_handle = open(link_stats_combined, 'w')
sankey_file_in_handle.write('MarkerGene,GenomicSeq,Number\n')
for each_linkage in marker_to_gnm_link_num:
    sankey_file_in_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
        each_linkage.split(marker_to_ctg_gnm_Key_connector)[0],
        each_linkage.split(marker_to_ctg_gnm_Key_connector)[1],
        marker_to_gnm_link_num[each_linkage]))
sankey_file_in_handle.close()

# write out linking reads
linking_reads_rd1_handle = open(linking_reads_rd1, 'w')
for each_link in marker_to_ctg_linkage_num_dict:
    current_linking_reads = marker_to_ctg_linking_reads_dict[each_link]
    linking_reads_rd1_handle.write('%s\t%s\n' % (each_link, ','.join(current_linking_reads)))
linking_reads_rd1_handle.close()


####################################### filter_linkages_iteratively ########################################

print('%s %s' % ((datetime.now().strftime(time_format)), 'filter_linkages_iteratively'))

pairwise_16s_iden_dict = blast_results_to_pairwise_16s_iden_dict(blast_results_all_vs_all_16s, min_aln_16s, min_cov_16s)

filter_linkages_iteratively(link_stats_combined, 'Number', pairwise_16s_iden_dict,
                            mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple,
                            min_iden_16s, min_link_num,
                            min_link_num, within_gnm_linkage_num_diff, link_stats_combined_filtered_s1)

