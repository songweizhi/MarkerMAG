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

        #################### overall ####################

        self.qualified_reads = False

        #################### round 1 16s ####################

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

        #################################################

        self.r1_longest_clp_cigar = ''
        self.r1_longest_clp_falg = ''
        self.r2_longest_clp_cigar = ''
        self.r2_longest_clp_falg = ''

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

wd                                              = '/Users/songweizhi/Desktop/mini_testing'
sam_file_mini_assembly_reformatted              = '%s/scaffolds_bowtie_reformatted.sam'  % wd
rd2_read_to_extract_flanking_16s                = {'MBARC26_20140984'}
to_extract_read_base_rd2_ctg                    = set()

min_M_len_16s                                   = 45
min_M_len_ctg                                   = 45
min_M_len_mini                                  = 45
min_M_pct                                       = 35
mismatch_cutoff                                 = 2
read_to_marker_connector                        = '___r___'
marker_to_ctg_gnm_Key_connector                 = '___M___'
gnm_to_ctg_connector                            = '___C___'
report_interval                                 = 25000
clp_pct_ctg_side_max_num                        = 65
clp_pct_ratio_cutoff                            = 3.5


mini_assembly_to_16s_reads  = '%s/mini_assembly_to_16s_reads.txt' % wd
mini_assembly_to_ctg_reads  = '%s/mini_assembly_to_ctg_reads.txt' % wd
#stats_GapFilling_ctg        = '%s/stats_GapFilling_ctg.txt' % wd


########################################################################################################################

mini_assembly_len_dict = {}
mini_assembly_mp_dict = {}
for each_read in open(sam_file_mini_assembly_reformatted):
    each_read_split = each_read.strip().split('\t')
    if each_read.startswith('@'):
        mini_assembly_id = ''
        mini_assembly_len = 0
        for each_element in each_read_split:
            if each_element.startswith('SN:'):
                mini_assembly_id = each_element[3:]
            if each_element.startswith('LN:'):
                mini_assembly_len = int(each_element[3:])
        mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
    else:
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])

            if read_id_base not in mini_assembly_mp_dict:
                mini_assembly_mp_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if ref_id not in mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict:
                    mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                else:
                    mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict[ref_id][ref_pos] = cigar

            if read_strand == '2':
                if ref_id not in mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict:
                    mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict[ref_id] = {ref_pos: cigar}
                else:
                    mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict[ref_id][ref_pos] = cigar

##################################################### parse MappingRecord_dict ####################################################

mini_assembly_to_16s_reads_dict = dict()
mini_assembly_to_ctg_reads_dict = dict()
processed_num = 0
to_extract_read_base = set()
for each_mp in mini_assembly_mp_dict:

    current_mp_record = mini_assembly_mp_dict[each_mp]
    mini_refs_to_ignore = set()

    ########## get lowest mismatch for r1/r2 16s refs ##########

    # get r1_ref_cigar_set
    r1_ref_cigar_set = set()
    for each_pos_dict in current_mp_record.r1_mini_ref_dict.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r1_ref_cigar_set.update(each_pos_dict_values)

    # get r2_ref_cigar_set
    r2_ref_cigar_set = set()
    for each_pos_dict in current_mp_record.r2_mini_ref_dict.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r2_ref_cigar_set.update(each_pos_dict_values)

    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_mini)
    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_mini)
    mini_assembly_mp_dict[each_mp].r1_mini_refs_lowest_mismatch = r1_ref_min_mismatch
    mini_assembly_mp_dict[each_mp].r2_mini_refs_lowest_mismatch = r2_ref_min_mismatch

    ########## filter r1 mini refs ##########

    r1_mini_refs_passed_qc = {}
    for r1_mini_ref in current_mp_record.r1_mini_ref_dict:
        r1_matched_pos_dict = current_mp_record.r1_mini_ref_dict[r1_mini_ref]

        # one read need to mapped to one 16S only for one time
        if len(r1_matched_pos_dict) > 1:
            mini_refs_to_ignore.add(r1_mini_ref)
        else:
            r1_mini_ref_pos = list(r1_matched_pos_dict.keys())[0]
            r1_mini_ref_cigar = r1_matched_pos_dict[r1_mini_ref_pos]
            r1_mini_ref_cigar_splitted = cigar_splitter(r1_mini_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r1_mini_ref_cigar_splitted)
            if both_end_clp is True:
                mini_refs_to_ignore.add(r1_mini_ref)
            else:
                # check mismatch
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(r1_mini_ref_cigar_splitted)
                if r1_ref_min_mismatch == 'NA':
                    mini_refs_to_ignore.add(r1_mini_ref)
                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                    mini_refs_to_ignore.add(r1_mini_ref)
                else:
                    # check aligned length
                    if r1_aligned_len < min_M_len_mini:
                        mini_refs_to_ignore.add(r1_mini_ref)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r1_mini_ref_cigar) or ('s' in r1_mini_ref_cigar):
                            clip_in_middle = True
                            if (r1_mini_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_mini_ref_pos == 1):
                                clip_in_middle = False
                            if (r1_mini_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r1_mini_ref_pos + r1_aligned_len - 1) == mini_assembly_len_dict[r1_mini_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            mini_refs_to_ignore.add(r1_mini_ref)
                        else:
                            r1_mini_refs_passed_qc[r1_mini_ref] = [r1_mini_ref_cigar]

    ########## filter r2 mini refs ##########

    r2_mini_refs_passed_qc = {}
    for r2_mini_ref in current_mp_record.r2_mini_ref_dict:
        r2_matched_pos_dict = current_mp_record.r2_mini_ref_dict[r2_mini_ref]

        # one read need to mapped to one 16S only once
        if len(r2_matched_pos_dict) > 1:
            mini_refs_to_ignore.add(r2_mini_ref)
        else:
            r2_mini_ref_pos = list(r2_matched_pos_dict.keys())[0]
            r2_mini_ref_cigar = r2_matched_pos_dict[r2_mini_ref_pos]
            r2_mini_ref_cigar_splitted = cigar_splitter(r2_mini_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r2_mini_ref_cigar_splitted)
            if both_end_clp is True:
                mini_refs_to_ignore.add(r2_mini_ref)
            else:
                # check mismatch
                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(r2_mini_ref_cigar_splitted)
                if r2_ref_min_mismatch == 'NA':
                    mini_refs_to_ignore.add(r2_mini_ref)
                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                    mini_refs_to_ignore.add(r2_mini_ref)
                else:
                    # check aligned length
                    if r2_aligned_len < min_M_len_mini:
                        mini_refs_to_ignore.add(r2_mini_ref)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r2_mini_ref_cigar) or ('s' in r2_mini_ref_cigar):
                            clip_in_middle = True
                            if (r2_mini_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_mini_ref_pos == 1):
                                clip_in_middle = False
                            if (r2_mini_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r2_mini_ref_pos + r2_aligned_len - 1) == mini_assembly_len_dict[r2_mini_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            mini_refs_to_ignore.add(r2_mini_ref)
                        else:
                            r2_mini_refs_passed_qc[r2_mini_ref] = [r2_mini_ref_cigar]

    ####################################################################################################

    r1_mini_refs_no_ignored = {key: value for key, value in r1_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
    r2_mini_refs_no_ignored = {key: value for key, value in r2_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
    shared_mini_ref_dict = {key: [r1_mini_refs_no_ignored[key][0], r2_mini_refs_no_ignored[key][0]] for key in set(r1_mini_refs_no_ignored).intersection(set(r2_mini_refs_no_ignored))}

    # both r1 and r2 mapped to the same mini assembly
    if len(shared_mini_ref_dict) > 0:
        #print('\n---------------------------------\n')
        #print('%s.1\tr1_mini_refs_no_ignored\t%s'   % ((each_mp, r1_mini_refs_no_ignored)))
        #print('%s.2\tr2_mini_refs_no_ignored\t%s'   % ((each_mp, r2_mini_refs_no_ignored)))
        print('%s\tshared_16s_ref_dict\t%s' % ((each_mp, shared_mini_ref_dict)))

        if (each_mp in rd2_read_to_extract_flanking_16s) and (each_mp not in to_extract_read_base_rd2_ctg):
            for each_mini_assembly in shared_mini_ref_dict:
                if each_mini_assembly not in mini_assembly_to_16s_reads_dict:
                    mini_assembly_to_16s_reads_dict[each_mini_assembly] = {each_mp}
                else:
                    mini_assembly_to_16s_reads_dict[each_mini_assembly].add(each_mp)

        elif (each_mp not in rd2_read_to_extract_flanking_16s) and (each_mp in to_extract_read_base_rd2_ctg):
            for each_mini_assembly in shared_mini_ref_dict:
                if each_mini_assembly not in mini_assembly_to_ctg_reads_dict:
                    mini_assembly_to_ctg_reads_dict[each_mini_assembly] = {each_mp}
                else:
                    mini_assembly_to_ctg_reads_dict[each_mini_assembly].add(each_mp)


# write out mini_assembly_to_16s_reads
mini_assembly_to_16s_reads_handle = open(mini_assembly_to_16s_reads, 'w')
for each_mini_assembly in mini_assembly_to_16s_reads_dict:
    mini_assembly_to_16s_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join([i for i in mini_assembly_to_16s_reads_dict[each_mini_assembly]])))
mini_assembly_to_16s_reads_handle.close()

# write out mini_assembly_to_ctg_reads
mini_assembly_to_ctg_reads_handle = open(mini_assembly_to_ctg_reads, 'w')
for each_mini_assembly in mini_assembly_to_ctg_reads_dict:
    mini_assembly_to_ctg_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join([i for i in mini_assembly_to_ctg_reads_dict[each_mini_assembly]])))
mini_assembly_to_ctg_reads_handle.close()



