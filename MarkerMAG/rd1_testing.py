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
'''
------------------------------------------------------------------------------------------------------------

# cami_hc_SILVA138_id99_50_subsample_100_3569
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_tmp_MarkerMAG_wd_up/hc_0512_tmp_step_1_wd
grep -E '@|cami_hc_SILVA138_id99_50_subsample_100_3569' /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0511_very_sensitive_MarkerMAG_wd/hc_0511_very_sensitive_step_1_wd/hc_0511_very_sensitive_input_reads_to_16S_reformatted.sam > input_reads_to_16S_subset.sam

# cami_hc_14___NODE_13073_length_6967_cov_0.639487
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_tmp_MarkerMAG_wd_up/hc_0512_tmp_step_1_wd
grep -E '@|cami_hc_14___C___NODE_13073_length_6967_cov_0.639487' rd1_extracted_to_gnm_reformatted.sam > rd1_extracted_to_gnm_subset.sam

------------------------------------------------------------------------------------------------------------

'''
# do not remove both ends clip alignments, important !!!!!!

wd                                              = '/Users/songweizhi/Desktop/consider_all'
input_reads_to_16s_sam_best_match               = '%s/input_reads_to_16S_subset.sam'            % wd
rd1_extracted_to_gnm_sam_reformatted_best_match = '%s/rd1_extracted_to_gnm_subset.sam'          % wd

min_M_len_16s                                   = 45
min_M_len_ctg                                   = 45
min_M_pct                                       = 35
mismatch_cutoff                                 = 2
read_to_marker_connector                        = '___r___'
marker_to_ctg_gnm_Key_connector                 = '___M___'
gnm_to_ctg_connector                            = '___C___'
report_interval                                 = 10000

link_stats_combined                             = '%s/hc_0509_stats_combined.txt'               % wd
rd1_clp_pct_diff_txt                            = '%s/rd1_clp_pct_diff.txt'                     % wd
rd1_clp_pct_diff_txt_to_ignore                  = '%s/rd1_clp_pct_diff_to_ignore.txt'           % wd
linking_reads_rd1                               = '%s/linking_reads_rd1.txt'                 % wd

# on katana
# input_reads_to_16s_sam_best_match               = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0511_very_sensitive_MarkerMAG_wd/hc_0511_very_sensitive_step_1_wd/hc_0511_very_sensitive_input_reads_to_16S_reformatted.sam'
# rd1_extracted_to_gnm_sam_reformatted_best_match = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0512_tmp_MarkerMAG_wd_up/hc_0512_tmp_step_1_wd/rd1_extracted_to_gnm_reformatted.sam'
# link_stats_combined                             = 'linkage_stats_combined.txt'
# rd1_clp_pct_diff_txt                            = 'rd1_clp_pct_diff.txt'


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

    read_mr         = MappingRecord_dict[qualified_read]
    r1_16s_refs     = read_mr.r1_16s_refs_no_ignored
    r2_16s_refs     = read_mr.r2_16s_refs_no_ignored
    shared_16s_refs = read_mr.shared_16s_refs_no_ignored
    r1_ctg_refs     = read_mr.r1_ctg_refs_no_ignored
    r2_ctg_refs     = read_mr.r2_ctg_refs_no_ignored
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
            if (clp_pct_ctg_side >= 65) and (clp_pct_ratio >= 5):
                to_ignore = True

        if to_ignore is True:
            linkages_to_ignore.add(each_link)
            rd1_clp_pct_diff_txt_to_ignore_handle.write('%s\t%s\t%s\t%s(%s/%s)\t%s(%s/%s)\t%s\n' % (marker_id, gnm_id, ctg_id, clp_pct_ctg_side, len(linkage_cigar_ctg_side_clp), len(linkage_cigar_ctg_side_all), clp_pct_16s_side, len(linkage_cigar_16s_side_clp), len(linkage_cigar_16s_side_all), clp_pct_ratio))
        else:
            rd1_clp_pct_diff_txt_handle.write('%s\t%s\t%s\t%s(%s/%s)\t%s(%s/%s)\t%s\n' % (marker_id, gnm_id, ctg_id, clp_pct_ctg_side, len(linkage_cigar_ctg_side_clp), len(linkage_cigar_ctg_side_all), clp_pct_16s_side, len(linkage_cigar_16s_side_clp), len(linkage_cigar_16s_side_all), clp_pct_ratio))
rd1_clp_pct_diff_txt_handle.close()
rd1_clp_pct_diff_txt_to_ignore_handle.close()

# ignore linkages if did not pass the clp pct ratio test
marker_to_ctg_linkage_num_dict_min3_passed_ratio_check = {}
for each_link in marker_to_ctg_linkage_num_dict:
    if marker_to_ctg_linkage_num_dict[each_link] >= 3:
        if each_link not in linkages_to_ignore:
            marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_link] = marker_to_ctg_linkage_num_dict[each_link]

# get number of linkages at genome level
marker_to_gnm_link_num = {}
for each_marker_to_ctg_key in marker_to_ctg_linkage_num_dict_min3_passed_ratio_check:
    marker_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[0]
    ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[1]
    gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
    marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
    if marker_to_gnm_key not in marker_to_gnm_link_num:
        marker_to_gnm_link_num[marker_to_gnm_key] = marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_marker_to_ctg_key]
    else:
        marker_to_gnm_link_num[marker_to_gnm_key] += marker_to_ctg_linkage_num_dict_min3_passed_ratio_check[each_marker_to_ctg_key]

# write out linkages at genome level
sankey_file_in_handle = open(link_stats_combined, 'w')
sankey_file_in_handle.write('MarkerGene,GenomicSeq,Number\n')
for each_linkage in marker_to_gnm_link_num:
    sankey_file_in_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
    each_linkage.split(marker_to_ctg_gnm_Key_connector)[0], each_linkage.split(marker_to_ctg_gnm_Key_connector)[1],
    marker_to_gnm_link_num[each_linkage]))
sankey_file_in_handle.close()

# write out linking reads
linking_reads_rd1_handle = open(linking_reads_rd1, 'w')
for each_link in marker_to_ctg_linkage_num_dict:
    marker_id = each_link.split(marker_to_ctg_gnm_Key_connector)[0]
    ctg_id = each_link.split(marker_to_ctg_gnm_Key_connector)[1]
    gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
    current_linking_reads = marker_to_ctg_linking_reads_dict[each_link]
    linking_reads_rd1_handle.write('%s\t%s\n' % (each_link, ','.join(current_linking_reads)))
linking_reads_rd1_handle.close()


