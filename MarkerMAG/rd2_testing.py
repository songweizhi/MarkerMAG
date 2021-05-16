import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
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

    aligned_pct = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
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

        self.r1_16s_ref_dict = dict()
        self.r2_16s_ref_dict = dict()

        self.r1_16s_refs_lowest_mismatch = None
        self.r2_16s_refs_lowest_mismatch = None

        self.r1_16s_refs_no_ignored = dict()
        self.r2_16s_refs_no_ignored = dict()
        self.shared_16s_refs_no_ignored = dict()

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()
        self.r1_ctg_refs_no_ignored = dict()
        self.r2_ctg_refs_no_ignored = dict()
        self.shared_ctg_refs_no_ignored = dict()

        self.r1_longest_clp_cigar = ''
        self.r1_longest_clp_falg = ''
        self.r2_longest_clp_cigar = ''
        self.r2_longest_clp_falg = ''

        self.qualified_reads = False
        self.matched_to_ctg = False
        self.both_mapped_to_16s = False
        self.consider_r1_unmapped_mate = False
        self.consider_r1_clipping_part = False
        self.consider_r2_unmapped_mate = False
        self.consider_r2_clipping_part = False

        self.consider_round_2 = False

        self.r1_filtered_refs = set()
        self.r2_filtered_refs = set()
        self.unmapped_r1_refs = set()
        self.unmapped_r2_refs = set()
        self.clipping_r1_refs = set()
        self.clipping_r2_refs = set()
        self.unmapped_r1_refs_with_pos = set()
        self.unmapped_r2_refs_with_pos = set()
        self.clipping_r1_refs_with_pos = set()
        self.clipping_r2_refs_with_pos = set()


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


def check_cigar_all_clp(cigar_list):
    cigar_all_clp = True
    for each_cigar in cigar_list:
        if ('S' not in each_cigar) and ('s' not in each_cigar):
            cigar_all_clp = False

    return cigar_all_clp


def get_min_max_cigar_S(cigar_list):

    # get cigar_S_set
    cigar_S_set = set()
    for each_cigar in cigar_list:
        if ('S' not in each_cigar) and ('s' not in each_cigar):
            cigar_S_set.add(0)
        else:
            cigar_splitted = cigar_splitter(each_cigar)

            if cigar_splitted[0][-1] in ['S', 's']:
                cigar_S_set.add(int(cigar_splitted[0][:-1]))

            if cigar_splitted[-1][-1] in ['S', 's']:
                cigar_S_set.add(int(cigar_splitted[1][:-1]))

    # get min_cigar_S and max_cigar_S
    min_cigar_S = 0
    max_cigar_S = 0
    if len(cigar_S_set) > 0:
        min_cigar_S = min(cigar_S_set)
        max_cigar_S = max(cigar_S_set)

    return min_cigar_S, max_cigar_S


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


########################################################################################################################

# do not remove both ends clip alignments, important !!!!!!

wd                                              = '/Users/songweizhi/Desktop/rd2_testing'
input_reads_to_16s_sam_best_match               = '%s/MBARC26_0515_input_reads_to_16S_reformatted.sam'  % wd
rd1_extracted_to_gnm_sam_reformatted_best_match = '%s/rd1_extracted_to_gnm_reformatted.sam'             % wd

min_M_len_16s                                   = 45
min_M_len_ctg                                   = 45
min_M_pct                                       = 35
mismatch_cutoff                                 = 2
read_to_marker_connector                        = '___r___'
marker_to_ctg_gnm_Key_connector                 = '___M___'
gnm_to_ctg_connector                            = '___C___'
report_interval                                 = 25000
clp_pct_ctg_side_max_num                        = 65
clp_pct_ratio_cutoff                            = 3.5
link_stats_combined                             = '%s/hc_0509_stats_combined.txt'               % wd
rd1_clp_pct_diff_txt                            = '%s/rd1_clp_pct_diff.txt'                     % wd
rd1_clp_pct_diff_txt_to_ignore                  = '%s/rd1_clp_pct_diff_to_ignore.txt'           % wd
linking_reads_rd1                               = '%s/linking_reads_rd1.txt'                    % wd
link_stats_combined_filtered_s1                 = '%s/MBARC26_0515_stats_combined_filtered.txt' % wd

# on katana
# input_reads_to_16s_sam_best_match               = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0511_very_sensitive_MarkerMAG_wd/hc_0511_very_sensitive_step_1_wd/hc_0511_very_sensitive_input_reads_to_16S_reformatted.sam'
# rd1_extracted_to_gnm_sam_reformatted_best_match = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_tmp_MarkerMAG_wd_up/hc_0512_tmp_step_1_wd/rd1_extracted_to_gnm_reformatted.sam'
# link_stats_combined                             = 'linkage_stats_combined.txt'
# rd1_clp_pct_diff_txt                            = 'rd1_clp_pct_diff.txt'

reads_file_r1_fasta = ''
reads_file_r2_fasta = ''
rd2_to_extract_flking_16s_r1_id = '%s/rd2_to_extract_flking_16s_r1_id.txt' % wd
rd2_to_extract_flking_16s_r2_id = '%s/rd2_to_extract_flking_16s_r2_id.txt' % wd
rd2_extracted_flking_16s_r1_seq_tmp = '%s/rd2_extracted_flking_16s_r1_tmp.fa' % wd
rd2_extracted_flking_16s_r2_seq_tmp = '%s/rd2_extracted_flking_16s_r2_tmp.fa' % wd
rd2_extracted_flking_16s_r1_seq = '%s/rd2_extracted_flking_16s_r1.fa' % wd
rd2_extracted_flking_16s_r2_seq = '%s/rd2_extracted_flking_16s_r2.fa' % wd

########################################################################################################################

marker_len_dict = {}
MappingRecord_dict = {}
for each_read in open(input_reads_to_16s_sam_best_match):
    each_read_split = each_read.strip().split('\t')

    if each_read.startswith('@'):
        marker_id = ''
        marker_len = 0
        for each_element in each_read_split:
            if each_element.startswith('SN:'):
                marker_id = each_element[3:]
            if each_element.startswith('LN:'):
                marker_len = int(each_element[3:])
        marker_len_dict[marker_id] = marker_len
    else:
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])

            if read_id_base not in MappingRecord_dict:
                MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if ref_id not in MappingRecord_dict[read_id_base].r1_16s_ref_dict:
                    MappingRecord_dict[read_id_base].r1_16s_ref_dict[ref_id] = {ref_pos:cigar}
                else:
                    MappingRecord_dict[read_id_base].r1_16s_ref_dict[ref_id][ref_pos] = cigar

            if read_strand == '2':
                if ref_id not in MappingRecord_dict[read_id_base].r2_16s_ref_dict:
                    MappingRecord_dict[read_id_base].r2_16s_ref_dict[ref_id] = {ref_pos:cigar}
                else:
                    MappingRecord_dict[read_id_base].r2_16s_ref_dict[ref_id][ref_pos] = cigar


##################################################### parse MappingRecord_dict ####################################################

processed_num = 0
to_extract_read_base = set()
for each_mp in MappingRecord_dict:

    current_mp_record = MappingRecord_dict[each_mp]
    refs_to_ignore = set()

    ########## get lowest mismatch for r1/r2 16s refs ##########

    # get r1_ref_cigar_set
    r1_ref_cigar_set = set()
    for each_pos_dict in current_mp_record.r1_16s_ref_dict.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r1_ref_cigar_set.update(each_pos_dict_values)

    # get r2_ref_cigar_set
    r2_ref_cigar_set = set()
    for each_pos_dict in current_mp_record.r2_16s_ref_dict.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r2_ref_cigar_set.update(each_pos_dict_values)

    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_16s)
    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_16s)
    MappingRecord_dict[each_mp].r1_16s_refs_lowest_mismatch = r1_ref_min_mismatch
    MappingRecord_dict[each_mp].r2_16s_refs_lowest_mismatch = r2_ref_min_mismatch

    ########## filter r1 16s refs ##########

    r1_16s_refs_passed_qc = {}
    for r1_16s_ref in current_mp_record.r1_16s_ref_dict:
        r1_matched_pos_dict = current_mp_record.r1_16s_ref_dict[r1_16s_ref]

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
                            if (r1_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_16s_ref_pos == 1):
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
    for r2_16s_ref in current_mp_record.r2_16s_ref_dict:
        r2_matched_pos_dict = current_mp_record.r2_16s_ref_dict[r2_16s_ref]

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
                            if (r2_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_16s_ref_pos == 1):
                                clip_in_middle = False
                            if (r2_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r2_16s_ref_pos + r2_aligned_len - 1) == marker_len_dict[r2_16s_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            refs_to_ignore.add(r2_16s_ref)
                        else:
                            r2_16s_refs_passed_qc[r2_16s_ref] = [r2_16s_ref_cigar]

    ####################################################################################################

    r1_16s_refs_no_ignored = {key: value for key, value in r1_16s_refs_passed_qc.items() if key not in refs_to_ignore}
    r2_16s_refs_no_ignored = {key: value for key, value in r2_16s_refs_passed_qc.items() if key not in refs_to_ignore}

    # no mate has no_ignored alignments
    if (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) == 0):
        pass

    # only r1 has no_ignored alignments
    elif (len(r1_16s_refs_no_ignored) > 0) and (len(r2_16s_refs_no_ignored) == 0):

        MappingRecord_dict[each_mp].qualified_reads = True
        MappingRecord_dict[each_mp].consider_r1_unmapped_mate = True
        MappingRecord_dict[each_mp].r1_16s_refs_no_ignored = r1_16s_refs_no_ignored
        to_extract_read_base.add(each_mp)

        # if the shortest S longer than min_M_ctg_len, also consider as clipping mapped
        r1_16s_refs_no_ignored_cigar_list = list(r1_16s_refs_no_ignored.values())
        r1_min_cigar_S, r1_max_cigar_S = get_min_max_cigar_S(r1_16s_refs_no_ignored_cigar_list)
        if r1_min_cigar_S >= min_M_len_ctg:
            MappingRecord_dict[each_mp].consider_r1_clipping_part = True
            to_extract_read_base.add(each_mp)

    # only r2 has no_ignored alignments
    elif (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) > 0):

        MappingRecord_dict[each_mp].qualified_reads = True
        MappingRecord_dict[each_mp].consider_r2_unmapped_mate = True
        MappingRecord_dict[each_mp].r2_16s_refs_no_ignored = r2_16s_refs_no_ignored
        to_extract_read_base.add(each_mp)

        # if the shortest S longer than min_M_ctg_len, also consider as clipping mapped
        r2_16s_refs_no_ignored_cigar_list = list(r2_16s_refs_no_ignored.values())
        r2_min_cigar_S, r2_max_cigar_S = get_min_max_cigar_S(r2_16s_refs_no_ignored_cigar_list)
        if r2_min_cigar_S >= min_M_len_ctg:
            MappingRecord_dict[each_mp].consider_r2_clipping_part = True
            to_extract_read_base.add(each_mp)

    # both r1 and r2 have no_ignored alignments
    else:
        shared_16s_ref_dict = {key: [r1_16s_refs_no_ignored[key][0], r2_16s_refs_no_ignored[key][0]] for key in set(r1_16s_refs_no_ignored).intersection(set(r2_16s_refs_no_ignored))}
        if len(shared_16s_ref_dict) > 0:
            MappingRecord_dict[each_mp].qualified_reads = True
            MappingRecord_dict[each_mp].both_mapped_to_16s = True
            MappingRecord_dict[each_mp].shared_16s_refs_no_ignored = shared_16s_ref_dict
            to_extract_read_base.add(each_mp)

    processed_num += 1
    if (processed_num % report_interval == 0):
        processed_pct = processed_num*100/len(MappingRecord_dict)
        processed_pct = float("{0:.2f}".format(processed_pct))
        print('Processed %sk (%s%s)' % (int(processed_num/1000), processed_pct, '%'))


# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    if MappingRecord_dict[each_mp].qualified_reads is False:
        MappingRecord_dict.pop(each_mp)


######################################### read mapping results of rd1 extracted mates into mp dict  #########################################

ctg_len_dict = {}
for each_read in open(rd1_extracted_to_gnm_sam_reformatted_best_match):
    each_read_split = each_read.strip().split('\t')

    if each_read.startswith('@'):
        ctg_id = ''
        ctg_len = 0
        for each_element in each_read_split:
            if each_element.startswith('SN:'):
                ctg_id = each_element[3:]
            if each_element.startswith('LN:'):
                ctg_len = int(each_element[3:])
        ctg_len_dict[ctg_id] = ctg_len
    else:
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])

            if read_id_base not in MappingRecord_dict:
                MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if ref_id not in MappingRecord_dict[read_id_base].r1_ctg_ref_dict:
                    MappingRecord_dict[read_id_base].r1_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                else:
                    MappingRecord_dict[read_id_base].r1_ctg_ref_dict[ref_id][ref_pos] = cigar

            if read_strand == '2':
                if ref_id not in MappingRecord_dict[read_id_base].r2_ctg_ref_dict:
                    MappingRecord_dict[read_id_base].r2_ctg_ref_dict[ref_id] = {ref_pos: cigar}
                else:
                    MappingRecord_dict[read_id_base].r2_ctg_ref_dict[ref_id][ref_pos] = cigar


##################################################### parse MappingRecord_dict ####################################################

processed_num = 0
for each_mp in MappingRecord_dict:

    # print('%s\tr1_ctg_ref_dict\t%s'       % (each_mp, current_mp_record.r1_ctg_ref_dict))
    # print('%s\tr2_ctg_ref_dict\t%s'       % (each_mp, current_mp_record.r2_ctg_ref_dict))

    current_mp_record = MappingRecord_dict[each_mp]
    refs_to_ignore_ctg = set()

    ########## filter r1 ctg refs ##########

    r1_ctg_refs_passed_qc = {}
    for r1_ctg_ref in current_mp_record.r1_ctg_ref_dict:
        r1_ctg_ref_matched_pos_dict = current_mp_record.r1_ctg_ref_dict[r1_ctg_ref]

        # one read need to mapped to one ctg only for one time
        if len(r1_ctg_ref_matched_pos_dict) > 1:
            refs_to_ignore_ctg.add(r1_ctg_ref)
        else:
            r1_ctg_ref_pos = list(r1_ctg_ref_matched_pos_dict.keys())[0]
            r1_ctg_ref_cigar = r1_ctg_ref_matched_pos_dict[r1_ctg_ref_pos]
            r1_ctg_ref_cigar_splitted = cigar_splitter(r1_ctg_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r1_ctg_ref_cigar_splitted)
            if both_end_clp is True:
                refs_to_ignore_ctg.add(r1_ctg_ref)
            else:
                # check mismatch
                r1_aligned_len_ctg, r1_aligned_pct_ctg, r1_clipping_len_ctg, r1_clipping_pct_ctg, r1_mismatch_pct_ctg = get_cigar_stats(r1_ctg_ref_cigar_splitted)
                if r1_mismatch_pct_ctg > mismatch_cutoff:
                    refs_to_ignore_ctg.add(r1_ctg_ref)
                else:
                    # check aligned length
                    if r1_aligned_len_ctg < min_M_len_ctg:
                        refs_to_ignore_ctg.add(r1_ctg_ref)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r1_ctg_ref_cigar) or ('s' in r1_ctg_ref_cigar):
                            clip_in_middle = True
                            if (r1_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_ctg_ref_pos == 1):
                                clip_in_middle = False
                            if (r1_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r1_ctg_ref_pos + r1_aligned_len_ctg - 1) == ctg_len_dict[r1_ctg_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            refs_to_ignore_ctg.add(r1_ctg_ref)
                        else:
                            r1_ctg_refs_passed_qc[r1_ctg_ref] = [r1_ctg_ref_cigar]

    ########## filter r2 ctg refs ##########

    r2_ctg_refs_passed_qc = {}
    for r2_ctg_ref in current_mp_record.r2_ctg_ref_dict:
        r2_ctg_ref_matched_pos_dict = current_mp_record.r2_ctg_ref_dict[r2_ctg_ref]

        # one read need to mapped to one ctg only for one time
        if len(r2_ctg_ref_matched_pos_dict) > 1:
            refs_to_ignore_ctg.add(r2_ctg_ref)
        else:
            r2_ctg_ref_pos = list(r2_ctg_ref_matched_pos_dict.keys())[0]
            r2_ctg_ref_cigar = r2_ctg_ref_matched_pos_dict[r2_ctg_ref_pos]
            r2_ctg_ref_cigar_splitted = cigar_splitter(r2_ctg_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
            if both_end_clp is True:
                refs_to_ignore_ctg.add(r2_ctg_ref)
            else:
                # check mismatch
                r2_aligned_len_ctg, r2_aligned_pct_ctg, r2_clipping_len_ctg, r2_clipping_pct_ctg, r2_mismatch_pct_ctg = get_cigar_stats(r2_ctg_ref_cigar_splitted)
                if r2_mismatch_pct_ctg > mismatch_cutoff:
                    refs_to_ignore_ctg.add(r2_ctg_ref)
                else:
                    # check aligned length
                    if r2_aligned_len_ctg < min_M_len_ctg:
                        refs_to_ignore_ctg.add(r2_ctg_ref)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r2_ctg_ref_cigar) or ('s' in r2_ctg_ref_cigar):
                            clip_in_middle = True
                            if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_ctg_ref_pos == 1):
                                clip_in_middle = False
                            if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r2_ctg_ref_pos + r2_aligned_len_ctg - 1) == ctg_len_dict[r2_ctg_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            refs_to_ignore_ctg.add(r2_ctg_ref)
                        else:
                            r2_ctg_refs_passed_qc[r2_ctg_ref] = [r2_ctg_ref_cigar]

    ####################################################################################################

    r1_ctg_refs_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}
    r2_ctg_refs_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in refs_to_ignore_ctg}

    # no mate matched to ctg
    if (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) == 0):
        pass

    # only r1 matched to ctg
    elif (len(r1_ctg_refs_no_ignored) > 0) and (len(r2_ctg_refs_no_ignored) == 0):
        MappingRecord_dict[each_mp].matched_to_ctg = True
        MappingRecord_dict[each_mp].r1_ctg_refs_no_ignored = r1_ctg_refs_no_ignored

    # only r2 matched to ctg
    elif (len(r1_ctg_refs_no_ignored) == 0) and (len(r2_ctg_refs_no_ignored) > 0):
        MappingRecord_dict[each_mp].matched_to_ctg = True
        MappingRecord_dict[each_mp].r2_ctg_refs_no_ignored = r2_ctg_refs_no_ignored

    # both r1 and r2 matched to ctg
    else:
        shared_ctg_ref_set = {key: [r1_ctg_refs_no_ignored[key][0], r2_ctg_refs_no_ignored[key][0]] for key in set(r1_ctg_refs_no_ignored).intersection(set(r2_ctg_refs_no_ignored))}
        if len(shared_ctg_ref_set) > 0:
            MappingRecord_dict[each_mp].matched_to_ctg = True
            MappingRecord_dict[each_mp].shared_ctg_refs_no_ignored = shared_ctg_ref_set

    processed_num += 1
    if (processed_num % report_interval == 0):
        processed_pct = processed_num*100/len(MappingRecord_dict)
        processed_pct = float("{0:.2f}".format(processed_pct))
        print('Processed %sk (%s%s)' % (int(processed_num/1000), processed_pct, '%'))


##################################################### get linkages from MappingRecord_dict #####################################################

marker_to_ctg_linkage_cigar_dict_16s_side = {}
marker_to_ctg_linkage_cigar_dict_ctg_side = {}
marker_to_ctg_linkage_num_dict = {}
marker_to_ctg_linking_reads_dict = {}
for qualified_read in MappingRecord_dict:

    read_mr = MappingRecord_dict[qualified_read]
    r1_16s_refs = read_mr.r1_16s_refs_no_ignored
    r2_16s_refs = read_mr.r2_16s_refs_no_ignored
    shared_16s_refs = read_mr.shared_16s_refs_no_ignored
    r1_ctg_refs = read_mr.r1_ctg_refs_no_ignored
    r2_ctg_refs = read_mr.r2_ctg_refs_no_ignored
    shared_ctg_refs = read_mr.shared_ctg_refs_no_ignored

    # combine dict
    combined_16s_refs_dict = {**r1_16s_refs, **r2_16s_refs, **shared_16s_refs}
    combined_ctg_refs_dict = {**r1_ctg_refs, **r2_ctg_refs, **shared_ctg_refs}

    if (len(combined_16s_refs_dict) > 0) and (len(combined_ctg_refs_dict) > 0):
        for each_16s_ref in combined_16s_refs_dict:
            current_cigar_list_16s_side = combined_16s_refs_dict[each_16s_ref]
            for each_ctg_ref in combined_ctg_refs_dict:
                current_cigar_list_ctg_side = combined_ctg_refs_dict[each_ctg_ref]
                marker_to_ctg_key = '%s%s%s' % (each_16s_ref, marker_to_ctg_gnm_Key_connector, each_ctg_ref)

                if marker_to_ctg_key not in marker_to_ctg_linkage_num_dict:
                    marker_to_ctg_linkage_num_dict[marker_to_ctg_key] = 1
                    marker_to_ctg_linking_reads_dict[marker_to_ctg_key] = {qualified_read}
                    marker_to_ctg_linkage_cigar_dict_16s_side[marker_to_ctg_key] = current_cigar_list_16s_side
                    marker_to_ctg_linkage_cigar_dict_ctg_side[marker_to_ctg_key] = current_cigar_list_ctg_side
                else:
                    marker_to_ctg_linkage_num_dict[marker_to_ctg_key] += 1
                    marker_to_ctg_linking_reads_dict[marker_to_ctg_key].add(qualified_read)
                    marker_to_ctg_linkage_cigar_dict_16s_side[marker_to_ctg_key] += current_cigar_list_16s_side
                    marker_to_ctg_linkage_cigar_dict_ctg_side[marker_to_ctg_key] += current_cigar_list_ctg_side

# calculate clp_pct_ratio
linkages_to_ignore = set()
for each_link in marker_to_ctg_linkage_num_dict:
    if marker_to_ctg_linkage_num_dict[each_link] >= 3:
        each_link_split = each_link.split(marker_to_ctg_gnm_Key_connector)
        marker_id = each_link_split[0]
        gnm_id = each_link_split[1].split(gnm_to_ctg_connector)[0]
        ctg_id = each_link_split[1].split(gnm_to_ctg_connector)[1]

        linkage_cigar_16s_side_all = marker_to_ctg_linkage_cigar_dict_16s_side[each_link]
        linkage_cigar_ctg_side_all = marker_to_ctg_linkage_cigar_dict_ctg_side[each_link]
        linkage_cigar_16s_side_clp = [i for i in linkage_cigar_16s_side_all if (('S' in i) or ('s' in i))]
        linkage_cigar_ctg_side_clp = [i for i in linkage_cigar_ctg_side_all if (('S' in i) or ('s' in i))]

        clp_pct_16s_side = 'na'
        if len(linkage_cigar_16s_side_all) > 0:
            clp_pct_16s_side = len(linkage_cigar_16s_side_clp) * 100 / len(linkage_cigar_16s_side_all)
            clp_pct_16s_side = float("{0:.2f}".format(clp_pct_16s_side))

        clp_pct_ctg_side = 'na'
        if len(linkage_cigar_ctg_side_all) > 0:
            clp_pct_ctg_side = len(linkage_cigar_ctg_side_clp) * 100 / len(linkage_cigar_ctg_side_all)
            clp_pct_ctg_side = float("{0:.2f}".format(clp_pct_ctg_side))

        clp_pct_ratio = 'na'
        if (clp_pct_16s_side != 'na') and (clp_pct_ctg_side != 'na'):
            if clp_pct_16s_side > 0:
                clp_pct_ratio = float("{0:.2f}".format(clp_pct_ctg_side / clp_pct_16s_side))

        to_ignore = False
        if clp_pct_ratio != 'na':
            if (clp_pct_ctg_side >= clp_pct_ctg_side_max_num) and (clp_pct_ratio >= clp_pct_ratio_cutoff):
                to_ignore = True
        else:
            if (clp_pct_ctg_side != 'na') and (clp_pct_16s_side != 'na'):
                if clp_pct_ctg_side >= clp_pct_ctg_side_max_num:
                    to_ignore = True

        if to_ignore is True:
            linkages_to_ignore.add(each_link)

# ignore linkages if did not pass the clp pct ratio test
marker_to_ctg_linkage_num_dict_min3_passed_ratio_check = {}
for each_link in marker_to_ctg_linkage_num_dict:
    if marker_to_ctg_linkage_num_dict[each_link] >= 3:
        if each_link not in linkages_to_ignore:
            marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_link] = marker_to_ctg_linkage_num_dict[each_link]


##########################################################################################################

# get linked marker genes and genomic sequences in step 1
linked_marker_gene_set = set()
linked_genomic_seq_set = set()
linked_16s_to_gnm_set = set()
for each_link in open(link_stats_combined_filtered_s1):
    if not each_link.startswith('MarkerGene,GenomicSeq,Number'):
        each_link_split = each_link.strip().split(',')
        id_16s = each_link_split[0][12:]
        id_gnm = each_link_split[1][12:]
        linked_marker_gene_set.add(id_16s)
        linked_genomic_seq_set.add(id_gnm)
        linked_16s_to_gnm_set.add('%s%s%s' % (id_16s, marker_to_ctg_gnm_Key_connector, id_gnm))

# get all rd1 linking reads
all_linking_reads_for_rd1_linkages = set()
for each_link in marker_to_ctg_linkage_num_dict_min3_passed_ratio_check:
    current_linking_reads = marker_to_ctg_linking_reads_dict[each_link]
    marker_to_gnm_key = each_link.split(gnm_to_ctg_connector)[0]
    if marker_to_gnm_key in linked_16s_to_gnm_set:
        all_linking_reads_for_rd1_linkages.update(current_linking_reads)

# get the id to extract
rd2_read_to_extract_flanking_16s = set()
for each_mp in MappingRecord_dict.copy():
    if each_mp in all_linking_reads_for_rd1_linkages:
        MappingRecord_dict.pop(each_mp)
    else:
        if len(MappingRecord_dict[each_mp].shared_16s_refs_no_ignored) > 0:
            MappingRecord_dict.pop(each_mp)
        elif (len(MappingRecord_dict[each_mp].r1_16s_refs_no_ignored) > 0) or (len(MappingRecord_dict[each_mp].r2_16s_refs_no_ignored) > 0):
            rd2_read_to_extract_flanking_16s.add(each_mp)

# write out id of linking reads for extraction
with open(rd2_to_extract_flking_16s_r1_id, 'w') as rd2_to_extract_flking_16s_r1_id_handle:
    rd2_to_extract_flking_16s_r1_id_handle.write('%s\n' % '\n'.join(sorted([('%s.1' % i) for i in rd2_read_to_extract_flanking_16s])))
with open(rd2_to_extract_flking_16s_r2_id, 'w') as rd2_to_extract_flking_16s_r2_id_handle:
    rd2_to_extract_flking_16s_r2_id_handle.write('%s\n' % '\n'.join(sorted([('%s.2' % i) for i in rd2_read_to_extract_flanking_16s])))

# extract reads with seqtk
seqtk_extract_cmd_rd2_flk_16s_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_to_extract_flking_16s_r1_id, rd2_extracted_flking_16s_r1_seq)
seqtk_extract_cmd_rd2_flk_16s_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_to_extract_flking_16s_r2_id, rd2_extracted_flking_16s_r2_seq)
#report_and_log((seqtk_extract_cmd_rd2_flk_16s_r1), pwd_log_file, True)
#report_and_log((seqtk_extract_cmd_rd2_flk_16s_r2), pwd_log_file, True)
#os.system(seqtk_extract_cmd_rd2_flk_16s_r1)
#os.system(seqtk_extract_cmd_rd2_flk_16s_r2)

# read extracted read sequences into dict
extract_rd2_flking_16s_read_seq_dict = {}
for extracted_r1 in SeqIO.parse(rd2_extracted_flking_16s_r1_seq_tmp, 'fasta'):
    extract_rd2_flking_16s_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
for extracted_r2 in SeqIO.parse(rd2_extracted_flking_16s_r2_seq_tmp, 'fasta'):
    extract_rd2_flking_16s_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)

# write out paired in the same order
rd2_extracted_flking_16s_r1_handle = open(rd2_extracted_flking_16s_r1_seq, 'w')
rd2_extracted_flking_16s_r2_handle = open(rd2_extracted_flking_16s_r2_seq, 'w')
for each_read_base in rd2_read_to_extract_flanking_16s:
    current_r1 = '%s.1' % each_read_base
    current_r2 = '%s.2' % each_read_base
    current_r1_seq = extract_rd2_flking_16s_read_seq_dict.get(current_r1, '')
    current_r2_seq = extract_rd2_flking_16s_read_seq_dict.get(current_r2, '')
    rd2_extracted_flking_16s_r1_handle.write('>%s\n' % current_r1)
    rd2_extracted_flking_16s_r1_handle.write('%s\n' % current_r1_seq)
    rd2_extracted_flking_16s_r2_handle.write('>%s\n' % current_r2)
    rd2_extracted_flking_16s_r2_handle.write('%s\n' % current_r2_seq)
rd2_extracted_flking_16s_r1_handle.close()
rd2_extracted_flking_16s_r2_handle.close()







