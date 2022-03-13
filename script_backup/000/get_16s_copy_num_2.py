#!/usr/bin/env python3

import os
import random


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


def check_cigar_quality(cigar_str, consider_ref_min_mismatch, ref_min_mismatch, mismatch_cutoff, min_M_len, ref_pos, marker_len):

    # check the following:
    # 1. both end clip  2. mismatch     3. aligned length   4. clp in the middle    5. ref_min_mismatch

    cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        aln_len, aln_pct, clp_len, clp_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

        # check mismatch
        passed_mismatch_check = False
        if consider_ref_min_mismatch is False:
            if mismatch_pct <= mismatch_cutoff:
                passed_mismatch_check = True
        else:
            if ref_min_mismatch != 'NA':
                if (mismatch_pct <= ref_min_mismatch) and (mismatch_pct <= mismatch_cutoff):
                    passed_mismatch_check = True

        if passed_mismatch_check is True:
            # check aligned length
            if aln_len >= min_M_len:
                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar_str) or ('s' in cigar_str):
                    clip_in_middle = True
                    if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + aln_len - 1) == marker_len:
                            clip_in_middle = False

                if clip_in_middle is False:
                    qualified_cigar = True

    return qualified_cigar


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


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


def r12_16s_ref_dict_to_str(r2_16s_ref_dict):

    top_dict_str = 'na'

    if r2_16s_ref_dict != {}:
        top_dict_list = []
        for sub_dict_key in r2_16s_ref_dict:
            sub_dict_value = r2_16s_ref_dict[sub_dict_key]
            bottom_dict_list = []
            for bottom_dict_key in sub_dict_value:
                bottom_dict_value = sub_dict_value[bottom_dict_key]
                bottom_dict_list.append('%s:%s' % (bottom_dict_key, bottom_dict_value))
            bottom_dict_str = '%s' % (';'.join(bottom_dict_list))
            top_dict_list.append('%s:::%s' % (sub_dict_key, bottom_dict_str))

        top_dict_str = ';;;'.join(top_dict_list)

    return top_dict_str


def get_r12_16s_ref_dict_from_str(r2_16s_ref_dict_str):

    ref_dict_from_str = {}
    if r2_16s_ref_dict_str != 'na':
        r2_16s_ref_dict_str_split = r2_16s_ref_dict_str.split(';;;')

        ref_dict_from_str = {}
        for each_sub_dict_str in r2_16s_ref_dict_str_split:
            each_sub_dict_str_split = each_sub_dict_str.split(':::')
            sub_dict_key = each_sub_dict_str_split[0]
            sub_dict_str = each_sub_dict_str_split[1]
            sub_dict_str_split = sub_dict_str.split(';')
            bottom_dict = dict()
            for each_bottom_dict_str in sub_dict_str_split:
                bottom_dict_str_split = each_bottom_dict_str.split(':')
                bottom_dict_key = int(bottom_dict_str_split[0])
                bottom_dict_value = bottom_dict_str_split[1]
                bottom_dict[bottom_dict_key] = bottom_dict_value
            ref_dict_from_str[sub_dict_key] = bottom_dict

    return ref_dict_from_str


def no_ignored_dict_to_str(r2_16s_refs_no_ignored_dict):
    top_dict_str = 'na'

    if r2_16s_refs_no_ignored_dict != {}:
        top_dict_list = []
        for each_ref in r2_16s_refs_no_ignored_dict:
            each_ref_cigar_list = r2_16s_refs_no_ignored_dict[each_ref]
            each_ref_cigar_str = ','.join(each_ref_cigar_list)
            top_dict_list.append('%s:%s' % (each_ref, each_ref_cigar_str))
        top_dict_str = ';'.join(top_dict_list)

    return top_dict_str


def get_no_ignored_dict_from_str(r2_16s_refs_no_ignored_dict_str):

    no_ignored_dict_from_str = {}
    if r2_16s_refs_no_ignored_dict_str != 'na':
        r2_16s_refs_no_ignored_dict_str_split = r2_16s_refs_no_ignored_dict_str.split(';')
        for each_sub_dict_str in r2_16s_refs_no_ignored_dict_str_split:
            each_sub_dict_str_split = each_sub_dict_str.split(':')
            ref_id = each_sub_dict_str_split[0]
            cigar_list_str = each_sub_dict_str_split[1]
            cigar_list = cigar_list_str.split(',')
            no_ignored_dict_from_str[ref_id] = cigar_list

    return no_ignored_dict_from_str


def parse_sam16s_worker(argument_list):

    sorted_sam          = argument_list[0]
    MappingRecord_file  = argument_list[1]
    min_M_len_16s       = argument_list[2]
    mismatch_cutoff     = argument_list[3]
    marker_len_dict     = argument_list[4]

    MappingRecord_file_handle = open(MappingRecord_file, 'w')
    current_read_base = ''
    current_read_base_r1_16s_ref_dict = dict()
    current_read_base_r2_16s_ref_dict = dict()
    with open(sorted_sam) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
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

                    current_read_base___qualified_reads             = 0
                    current_read_base___consider_r1_unmapped_mate   = 0
                    current_read_base___consider_r2_unmapped_mate   = 0
                    current_read_base___both_mapped_to_16s          = 0
                    current_read_base___r1_16s_ref_dict             = dict()
                    current_read_base___r2_16s_ref_dict             = dict()
                    current_read_base___r1_16s_refs_no_ignored      = dict()
                    current_read_base___r2_16s_refs_no_ignored      = dict()
                    current_read_base___shared_16s_refs_no_ignored  = dict()

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
                    r1_16s_refs_passed_qc_with_pos = {}
                    for r1_16s_ref in current_read_base_r1_16s_ref_dict:
                        r1_matched_pos_dict = current_read_base_r1_16s_ref_dict[r1_16s_ref]

                        # one read need to mapped to one 16S only for one time
                        if len(r1_matched_pos_dict) > 1:
                            refs_to_ignore.add(r1_16s_ref)
                        else:
                            r1_16s_ref_pos = list(r1_matched_pos_dict.keys())[0]
                            r1_16s_ref_cigar = r1_matched_pos_dict[r1_16s_ref_pos]

                            r1_16s_ref_qualified_cigar = check_cigar_quality(r1_16s_ref_cigar, True, r1_ref_min_mismatch,
                                                                             mismatch_cutoff, min_M_len_16s, r1_16s_ref_pos,
                                                                             marker_len_dict[r1_16s_ref])

                            # exclude the ref if cigar not qualified
                            if r1_16s_ref_qualified_cigar is False:
                                refs_to_ignore.add(r1_16s_ref)
                            else:
                                r1_16s_refs_passed_qc[r1_16s_ref] = [r1_16s_ref_cigar]
                                r1_16s_refs_passed_qc_with_pos[r1_16s_ref] = {r1_16s_ref_pos: r1_16s_ref_cigar}

                    ########## filter r2 16s refs ##########

                    r2_16s_refs_passed_qc = {}
                    r2_16s_refs_passed_qc_with_pos = {}
                    for r2_16s_ref in current_read_base_r2_16s_ref_dict:
                        r2_matched_pos_dict = current_read_base_r2_16s_ref_dict[r2_16s_ref]

                        # one read need to mapped to one 16S only once
                        if len(r2_matched_pos_dict) > 1:
                            refs_to_ignore.add(r2_16s_ref)
                        else:
                            r2_16s_ref_pos = list(r2_matched_pos_dict.keys())[0]
                            r2_16s_ref_cigar = r2_matched_pos_dict[r2_16s_ref_pos]

                            r2_16s_ref_qualified_cigar = check_cigar_quality(r2_16s_ref_cigar, True, r2_ref_min_mismatch,
                                                                             mismatch_cutoff, min_M_len_16s, r2_16s_ref_pos,
                                                                             marker_len_dict[r2_16s_ref])

                            # exclude the ref if cigar not qualified
                            if r2_16s_ref_qualified_cigar is False:
                                refs_to_ignore.add(r2_16s_ref)
                            else:
                                r2_16s_refs_passed_qc[r2_16s_ref] = [r2_16s_ref_cigar]
                                r2_16s_refs_passed_qc_with_pos[r2_16s_ref] = {r2_16s_ref_pos: r2_16s_ref_cigar}

                    #################################

                    r1_16s_refs_no_ignored          = {key: value for key, value in r1_16s_refs_passed_qc.items() if key not in refs_to_ignore}
                    r2_16s_refs_no_ignored          = {key: value for key, value in r2_16s_refs_passed_qc.items() if key not in refs_to_ignore}
                    r1_16s_refs_no_ignored_with_pos = {key: value for key, value in r1_16s_refs_passed_qc_with_pos.items() if key not in refs_to_ignore}
                    r2_16s_refs_no_ignored_with_pos = {key: value for key, value in r2_16s_refs_passed_qc_with_pos.items() if key not in refs_to_ignore}

                    # no mate has no_ignored alignments
                    if (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) == 0):
                        pass

                    # only r1 has no_ignored alignments
                    elif (len(r1_16s_refs_no_ignored) > 0) and (len(r2_16s_refs_no_ignored) == 0):

                        current_read_base___qualified_reads = 1
                        current_read_base___consider_r1_unmapped_mate = 1
                        current_read_base___r1_16s_refs_no_ignored = r1_16s_refs_no_ignored
                        current_read_base___r1_16s_ref_dict = current_read_base_r1_16s_ref_dict

                    # only r2 has no_ignored alignments
                    elif (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) > 0):

                        current_read_base___qualified_reads = 1
                        current_read_base___consider_r2_unmapped_mate = 1
                        current_read_base___r2_16s_refs_no_ignored = r2_16s_refs_no_ignored
                        current_read_base___r2_16s_ref_dict = current_read_base_r2_16s_ref_dict

                    # both r1 and r2 have no_ignored alignments
                    else:
                        shared_16s_ref_dict = {key: [r1_16s_refs_no_ignored[key][0], r2_16s_refs_no_ignored[key][0]] for key
                                               in set(r1_16s_refs_no_ignored).intersection(set(r2_16s_refs_no_ignored))}
                        if len(shared_16s_ref_dict) > 0:

                            current_read_base___qualified_reads = 1
                            current_read_base___both_mapped_to_16s = 1
                            current_read_base___shared_16s_refs_no_ignored = shared_16s_ref_dict
                            current_read_base___r1_16s_ref_dict = current_read_base_r1_16s_ref_dict
                            current_read_base___r2_16s_ref_dict = current_read_base_r2_16s_ref_dict

                    if current_read_base___qualified_reads == 1:
                        current_read_base___r1_16s_ref_dict_str            = r12_16s_ref_dict_to_str(current_read_base___r1_16s_ref_dict)
                        current_read_base___r2_16s_ref_dict_str            = r12_16s_ref_dict_to_str(current_read_base___r2_16s_ref_dict)
                        current_read_base___r1_16s_refs_no_ignored_str     = no_ignored_dict_to_str(current_read_base___r1_16s_refs_no_ignored)
                        current_read_base___r2_16s_refs_no_ignored_str     = no_ignored_dict_to_str(current_read_base___r2_16s_refs_no_ignored)
                        current_read_base___shared_16s_refs_no_ignored_str = no_ignored_dict_to_str(current_read_base___shared_16s_refs_no_ignored)
                        MappingRecord_file_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (current_read_base,
                                                                                                      current_read_base___qualified_reads,
                                                                                                      current_read_base___consider_r1_unmapped_mate,
                                                                                                      current_read_base___consider_r2_unmapped_mate,
                                                                                                      current_read_base___both_mapped_to_16s,
                                                                                                      current_read_base___r1_16s_ref_dict_str,
                                                                                                      current_read_base___r2_16s_ref_dict_str,
                                                                                                      current_read_base___r1_16s_refs_no_ignored_str,
                                                                                                      current_read_base___r2_16s_refs_no_ignored_str,
                                                                                                      current_read_base___shared_16s_refs_no_ignored_str))


                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_16s_ref_dict = dict()
                    current_read_base_r2_16s_ref_dict = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}

    MappingRecord_file_handle.close()


###################################### quality filter 16S sequences ######################################

ref_genome_metadata     = '/Users/songweizhi/Desktop/666/reference_genome_metadata.txt'
markermag_op            = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_linkages_by_genome.txt'
mag_depth_txt           = '/Users/songweizhi/Desktop/666/MBARC26_refined_bins_50_5_depth.txt'
both_pair_mapped        = True
min_M_len_16s           = 70
mismatch_cutoff         = 2

sam_basename            = 'reads_16s_to_16s_sam_all'

pwd_sam_file                        = '/Users/songweizhi/Desktop/666/000/%s.sam'                                     % sam_basename


pwd_sam_mp_file                     = '/Users/songweizhi/Desktop/666/000/%s_all_mis%s_minM%s_mr.txt'                 % (sam_basename, mismatch_cutoff, min_M_len_16s)
pwd_sam_file_filtered               = '/Users/songweizhi/Desktop/666/000/%s_all_mis%s_minM%s_filtered.sam'           % (sam_basename, mismatch_cutoff, min_M_len_16s)
pwd_sam_file_filtered_random        = '/Users/songweizhi/Desktop/666/000/%s_all_mis%s_minM%s_filtered_random.sam'    % (sam_basename, mismatch_cutoff, min_M_len_16s)
if both_pair_mapped is True:
    pwd_sam_mp_file                 = '/Users/songweizhi/Desktop/666/000/%s_both_mis%s_minM%s_mr.txt'                % (sam_basename, mismatch_cutoff, min_M_len_16s)
    pwd_sam_file_filtered           = '/Users/songweizhi/Desktop/666/000/%s_both_mis%s_minM%s_filtered.sam'          % (sam_basename, mismatch_cutoff, min_M_len_16s)
    pwd_sam_file_filtered_random    = '/Users/songweizhi/Desktop/666/000/%s_both_mis%s_minM%s_filtered_random.sam'   % (sam_basename, mismatch_cutoff, min_M_len_16s)


############################################# read in ref_gnm/MAG metadata #############################################

reported_ref_16s_copy_num_dict = {}
calculated_ref_16s_copy_num_dict = {}
ref_16s_total_depth_dict = {}
ref_gnm_depth_dict = {}
for each_ref in open(ref_genome_metadata):
    if not each_ref.startswith('Genome\t'):
        each_ref_split = each_ref.strip().split('\t')
        ref_id = each_ref_split[0]
        reported_copy_num = int(each_ref_split[1])

        if each_ref_split[2] == 'NA':
            calculated_copy_num = each_ref_split[2]
        else:
            calculated_copy_num = float(each_ref_split[2])

        ref_16s_total_depth = float(each_ref_split[3])
        ref_gnm_depth = float(each_ref_split[4])
        reported_ref_16s_copy_num_dict[ref_id] = reported_copy_num
        calculated_ref_16s_copy_num_dict[ref_id] = calculated_copy_num
        ref_16s_total_depth_dict[ref_id] = ref_16s_total_depth
        ref_gnm_depth_dict[ref_id] = ref_gnm_depth


mag_depth_dict = {}
for each_mag_depth in open(mag_depth_txt):
    if not each_mag_depth.startswith('MAG	Length(bp)	Depth'):
        each_mag_depth_split = each_mag_depth.strip().split('\t')
        mag_depth_dict[each_mag_depth_split[0]] = float(each_mag_depth_split[2])


########################################################################################################################

marker_len_dict = {}
read_len_dict = {}
for each_line in open(pwd_sam_file):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('@'):
        mini_assembly_id = ''
        mini_assembly_len = 0
        for each_element in each_line_split:
            if each_element.startswith('SN:'):
                mini_assembly_id = each_element[3:]
            if each_element.startswith('LN:'):
                mini_assembly_len = int(each_element[3:])
        marker_len_dict[mini_assembly_id] = mini_assembly_len


if os.path.isfile(pwd_sam_mp_file) is False:
    parse_sam16s_worker([pwd_sam_file, pwd_sam_mp_file, min_M_len_16s, mismatch_cutoff, marker_len_dict])


# reads filter results into dict
MappingRecord_dict = {}
with open(pwd_sam_mp_file) as each_mp_file_opened:
    for each_read_base in each_mp_file_opened:
        each_read_base_split = each_read_base.strip().split('\t')
        current_read_base___id = each_read_base_split[0]
        current_read_base___qualified_reads = each_read_base_split[1]
        current_read_base___consider_r1_unmapped_mate = each_read_base_split[2]
        current_read_base___consider_r2_unmapped_mate = each_read_base_split[3]
        current_read_base___both_mapped_to_16s = each_read_base_split[4]
        current_read_base___r1_16s_ref_dict = get_r12_16s_ref_dict_from_str(each_read_base_split[5])
        current_read_base___r2_16s_ref_dict = get_r12_16s_ref_dict_from_str(each_read_base_split[6])
        current_read_base___r1_16s_refs_no_ignored = get_no_ignored_dict_from_str(each_read_base_split[7])
        current_read_base___r2_16s_refs_no_ignored = get_no_ignored_dict_from_str(each_read_base_split[8])
        current_read_base___shared_16s_refs_no_ignored = get_no_ignored_dict_from_str(each_read_base_split[9])

        MappingRecord_dict[current_read_base___id] = MappingRecord()
        MappingRecord_dict[current_read_base___id].qualified_reads = True

        if current_read_base___consider_r1_unmapped_mate == '1':
            MappingRecord_dict[current_read_base___id].consider_r1_unmapped_mate = True
        if current_read_base___consider_r2_unmapped_mate == '1':
            MappingRecord_dict[current_read_base___id].consider_r2_unmapped_mate = True
        if current_read_base___both_mapped_to_16s == '1':
            MappingRecord_dict[current_read_base___id].both_mapped_to_16s = True

        MappingRecord_dict[current_read_base___id].r1_16s_ref_dict = current_read_base___r1_16s_ref_dict
        MappingRecord_dict[current_read_base___id].r2_16s_ref_dict = current_read_base___r2_16s_ref_dict
        MappingRecord_dict[current_read_base___id].r1_16s_refs_no_ignored = current_read_base___r1_16s_refs_no_ignored
        MappingRecord_dict[current_read_base___id].r2_16s_refs_no_ignored = current_read_base___r2_16s_refs_no_ignored
        MappingRecord_dict[current_read_base___id].shared_16s_refs_no_ignored = current_read_base___shared_16s_refs_no_ignored


# filter sam file
pwd_sam_file_filtered_handle = open(pwd_sam_file_filtered, 'w')
for each_line in open(pwd_sam_file):
    if each_line.startswith('@'):
        pwd_sam_file_filtered_handle.write(each_line)
    else:
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        read_base = read_id[:-2]
        read_strand = read_id[-1]
        ref_id = each_line_split[2]
        if read_base in MappingRecord_dict:

            to_write = False

            if both_pair_mapped is True:
                if (ref_id in MappingRecord_dict[read_base].shared_16s_refs_no_ignored):
                    to_write = True
            else:
                if read_strand == '1':
                    if (ref_id in MappingRecord_dict[read_base].r1_16s_refs_no_ignored):
                        to_write = True
                if read_strand == '2':
                    if (ref_id in MappingRecord_dict[read_base].r2_16s_refs_no_ignored):
                        to_write = True
            if to_write is True:
                pwd_sam_file_filtered_handle.write(each_line)
pwd_sam_file_filtered_handle.close()


linked_16s_set = set()
gnm_to_linked_16s_dict = {}
linked_16s_to_gnm_dict = {}
for each_linkage in open(markermag_op):
    each_linkage_split = each_linkage.strip().split('\t')
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Round'):
        id_16s = each_linkage_split[0]
        id_gnm = each_linkage_split[1]
        linked_16s_set.add(id_16s)
        if id_gnm not in gnm_to_linked_16s_dict:
            gnm_to_linked_16s_dict[id_gnm] = {id_16s}
        else:
            gnm_to_linked_16s_dict[id_gnm].add(id_16s)
        linked_16s_to_gnm_dict[id_16s] = id_gnm


read_to_ref_dict = {}
ref_len_dict = {}
read_len_dict = {}
for each_line in open(pwd_sam_file_filtered):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('@'):
        mini_assembly_id = ''
        mini_assembly_len = 0
        for each_element in each_line_split:
            if each_element.startswith('SN:'):
                mini_assembly_id = each_element[3:]
            if each_element.startswith('LN:'):
                mini_assembly_len = int(each_element[3:])
        ref_len_dict[mini_assembly_id] = mini_assembly_len
    else:
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        read_len = len(each_line_split[9])

        #if ref_id in linked_16s_set:
        if read_id not in read_to_ref_dict:
            read_to_ref_dict[read_id] = {ref_id}
        else:
            read_to_ref_dict[read_id].add(ref_id)
        read_len_dict[read_id] = read_len


ref_to_read_dict_overall = {}
ref_to_read_dict_random_one = {}
read_to_ref_dict_random_one = {}
for each_read in read_to_ref_dict:
    matched_ref_list = [i for i in read_to_ref_dict[each_read]]
    random_ref = random.sample(matched_ref_list, 1)[0]
    read_to_ref_dict_random_one[each_read] = random_ref
    if random_ref not in ref_to_read_dict_random_one:
        ref_to_read_dict_random_one[random_ref] = {each_read}
    else:
        ref_to_read_dict_random_one[random_ref].add(each_read)

    for each_matched_ref in matched_ref_list:
        if each_matched_ref not in ref_to_read_dict_overall:
            ref_to_read_dict_overall[each_matched_ref] = {each_read}
        else:
            ref_to_read_dict_overall[each_matched_ref].add(each_read)


# filter sam file
pwd_sam_file_filtered_random_handle = open(pwd_sam_file_filtered_random, 'w')
for each_line in open(pwd_sam_file_filtered):
    if each_line.startswith('@'):
        pwd_sam_file_filtered_random_handle.write(each_line)
    else:
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        random_ref = read_to_ref_dict_random_one.get(read_id, 'NA')
        if ref_id == random_ref:
            pwd_sam_file_filtered_random_handle.write(each_line)
pwd_sam_file_filtered_random_handle.close()


ref_depth_dict = {}
ref_cp_num_dict = {}
for each_ref in ref_to_read_dict_random_one:
    if each_ref in linked_16s_set:
        matched_reads_list = ref_to_read_dict_random_one[each_ref]
        ref_len = ref_len_dict[each_ref]
        matched_reads_total_len = 0
        for each_read in matched_reads_list:
            matched_reads_total_len += read_len_dict.get(each_read)

        ref_gnm_depth = matched_reads_total_len / ref_len
        ref_gnm_depth = float("{0:.2f}".format(ref_gnm_depth))
        ref_depth_dict[each_ref] = ref_gnm_depth
        linked_mag = linked_16s_to_gnm_dict[each_ref]
        linked_mag_depth = mag_depth_dict[linked_mag]

        ref_cp_num = ref_gnm_depth / linked_mag_depth
        ref_cp_num = float("{0:.2f}".format(ref_cp_num))
        ref_cp_num_dict[each_ref] = ref_cp_num

gnm_total_cp_num_dict_sep = {}
gnm_total_cp_num_dict = {}
for each_ref in ref_cp_num_dict:
    ref_gnm = each_ref.split('_')[0]
    ref_cp_num = ref_cp_num_dict[each_ref]

    if ref_gnm not in gnm_total_cp_num_dict:
        gnm_total_cp_num_dict[ref_gnm] = ref_cp_num
        gnm_total_cp_num_dict_sep[ref_gnm] = [ref_cp_num]
    else:
        gnm_total_cp_num_dict[ref_gnm] += ref_cp_num
        gnm_total_cp_num_dict_sep[ref_gnm].append(ref_cp_num)


print('Genome\tMAG\tREF(calculated)\tREF(reported)\tSep')
for each_mag in gnm_to_linked_16s_dict:
    mag_16s_copy_num = gnm_total_cp_num_dict[each_mag]
    mag_16s_copy_num = float("{0:.2f}".format(mag_16s_copy_num))
    mag_16s_copy_num_sep = gnm_total_cp_num_dict_sep[each_mag]
    mag_16s_copy_num_sep_str = '+'.join([str(i) for i in mag_16s_copy_num_sep])
    calculated_ref_16s_copy_num = calculated_ref_16s_copy_num_dict[each_mag]
    reported_ref_16s_copy_num = reported_ref_16s_copy_num_dict[each_mag]

    print('%s\t%s\t%s\t%s\t%s' % (each_mag, mag_16s_copy_num, calculated_ref_16s_copy_num, reported_ref_16s_copy_num, mag_16s_copy_num_sep_str))

print('\n')


############
print('Genome\tCopy_num(MAG)\tEstimated_copy_num(ref)\tReported_copy_num(ref)\tmag_16s_depth\tref_16s_depth\tmag_depth\tref_gnm_depth')
for each_mag in gnm_to_linked_16s_dict:
    linked_16s_set = gnm_to_linked_16s_dict[each_mag]
    linked_16s_len_list = [ref_len_dict[s16] for s16 in linked_16s_set]
    linked_16s_mean_len = sum(linked_16s_len_list)/len(linked_16s_len_list)
    linked_16s_total_len = sum(linked_16s_len_list)
    current_mag_linked_16s_read_set = set()
    for each_linked_16s in linked_16s_set:
        current_16s_mapped_reads = ref_to_read_dict_overall.get(each_linked_16s, [])
        for i in current_16s_mapped_reads:
            current_mag_linked_16s_read_set.add(i)

    current_mag_linked_16s_total_len = 0
    for read_len_id in current_mag_linked_16s_read_set:
        current_mag_linked_16s_total_len += read_len_dict[read_len_id]

    ref_gnm_depth = ref_gnm_depth_dict[each_mag]
    calculated_ref_16s_copy_num = calculated_ref_16s_copy_num_dict[each_mag]
    reported_ref_16s_copy_num = reported_ref_16s_copy_num_dict[each_mag]
    ref_gnm_16s_depth = ref_16s_total_depth_dict[each_mag]

    mag_16s_depth = current_mag_linked_16s_total_len / linked_16s_mean_len
    mag_16s_depth = float("{0:.2f}".format(mag_16s_depth))
    mag_depth = mag_depth_dict[each_mag]
    mag_16s_copy_num_per_gnm = mag_16s_depth / mag_depth
    mag_16s_copy_num_per_gnm = float("{0:.2f}".format(mag_16s_copy_num_per_gnm))

    #print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_mag, copy_num, calculated_ref_16s_copy_num, reported_ref_16s_copy_num, mag_16s_depth, ref_16s_depth, mag_depth, ref_depth))
    print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_mag, mag_16s_copy_num_per_gnm, calculated_ref_16s_copy_num, reported_ref_16s_copy_num, mag_16s_depth, ref_gnm_16s_depth, mag_depth, ref_gnm_depth))

