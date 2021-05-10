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
        self.r1_16s_refs_no_ignored = set()
        self.r2_16s_refs_no_ignored = set()
        self.shared_16s_refs_no_ignored = set()

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()
        self.r1_ctg_refs_no_ignored = set()
        self.r2_ctg_refs_no_ignored = set()
        self.shared_ctg_refs_no_ignored = set()


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


def keep_best_matches_in_sam_keep_short_M(sam_in, min_M_len, sam_out):

    # get read_to_cigar_dict
    read_to_cigar_dict = {}
    for each_line in open(sam_in):
        each_line_split = each_line.strip().split('\t')
        if not each_line.startswith('@'):
            read_id = each_line_split[0]
            cigar = each_line_split[5]
            if cigar != '*':
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    if read_id not in read_to_cigar_dict:
                        read_to_cigar_dict[read_id] = {cigar}
                    else:
                        read_to_cigar_dict[read_id].add(cigar)

    # get min_mismatch for each read
    read_min_mismatch_dict = {}
    for each_read in read_to_cigar_dict:
        read_mismatch_set_all_M = set()
        read_mismatch_set_long_M = set()
        for each_cigar in read_to_cigar_dict[each_read]:
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
            read_mismatch_set_all_M.add(mismatch_pct)
            if aligned_len >= min_M_len:
                read_mismatch_set_long_M.add(mismatch_pct)
        read_min_mismatch = min(read_mismatch_set_all_M)
        if len(read_mismatch_set_long_M) > 0:
            read_min_mismatch = min(read_mismatch_set_long_M)
        read_min_mismatch_dict[each_read] = read_min_mismatch

    sam_file_best_match_handle = open(sam_out, 'w')
    for each_line in open(sam_in):
        if each_line.startswith('@'):
            sam_file_best_match_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            cigar = each_line_split[5]
            if cigar == '*':
                sam_file_best_match_handle.write(each_line)
            else:
                read_id = each_line_split[0]
                cigar_split = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_split)
                    if mismatch_pct <= (read_min_mismatch_dict[read_id] * 1.5):
                        sam_file_best_match_handle.write(each_line)
    sam_file_best_match_handle.close()


########################################################################################################################

# do not remove both ends clip alignments, important !!!!!!

wd                                  = '/Users/songweizhi/Desktop/consider_all'
input_reads_to_16s_sam_best_match   = '%s/hc_0509_input_reads_to_16S_best_match_subset.sam' % wd
min_M_len_16s                       = 45
min_M_len_ctg                       = 45
min_M_pct                           = 35
mismatch_cutoff                     = 2
read_to_marker_connector            = '___r___'
marker_to_ctg_gnm_Key_connector     = '___M___'
gnm_to_ctg_connector = '___C___'
num_threads = 12
combined_input_gnms_no_ext = 'cami_hc_refined_bins_combined'
input_r1_fasta      = 'cami_hc_combined_R1.fasta'
input_r2_fasta      = 'cami_hc_combined_R2.fasta'
marker_gene_seqs    = '%s/cami_hc_SILVA138_uclust_0.999_polished.fa'    % wd
combined_input_gnms = '%s/cami_hc_refined_bins_combined.fa'             % wd

rd1_r1_to_extract                               = '%s/rd1_r1_to_extract.txt'                            % wd
rd1_r2_to_extract                               = '%s/rd1_r2_to_extract.txt'                            % wd
rd1_extracted_all_r1                            = '%s/rd1_extracted_all_r1.fasta'                       % wd
rd1_extracted_all_r2                            = '%s/rd1_extracted_all_r2.fasta'                       % wd
rd1_extracted_p_r1                              = '%s/rd1_extracted_R1.fasta'                           % wd
rd1_extracted_p_r2                              = '%s/rd1_extracted_R2.fasta'                           % wd
rd1_extracted_up                                = '%s/rd1_extracted_UP.fasta'                           % wd
rd1_extracted_to_gnm_sam                        = '%s/rd1_extracted_to_gnm.sam'                         % wd
rd1_extracted_to_gnm_sam_log                    = '%s/rd1_extracted_to_gnm.sam.log'                     % wd
rd1_extracted_to_gnm_sam_reformatted            = '%s/rd1_extracted_to_gnm_reformatted.sam'             % wd
rd1_extracted_to_gnm_sam_reformat_log           = '%s/rd1_extracted_to_gnm_reformat.log'                % wd
rd1_extracted_to_gnm_sam_reformatted_best_match = '%s/rd1_extracted_to_gnm_reformatted_best_match.sam'  % wd

link_stats_combined             = '%s/hc_0509_stats_combined.txt'       % wd
linking_reads_tab               = '%s/linking_reads.txt'                % wd
linking_reads_r1_txt            = '%s/rd1_linking_reads_R1.txt'         % wd
linking_reads_r2_txt            = '%s/rd1_linking_reads_R2.txt'         % wd
linking_reads_r1_fasta          = '%s/rd1_linking_reads_R1.fasta'       % wd
linking_reads_r2_fasta          = '%s/rd1_linking_reads_R2.fasta'       % wd
linked_ends_rd1_txt             = '%s/linked_ends_rd1.txt'              % wd
linked_contigs_txt              = '%s/linked_contigs_rd1.txt'               % wd
linked_contigs_fasta            = '%s/linked_contigs_rd1.fasta'               % wd

mafft_seq_folder                = '%s/mafft_seq_folder' % wd
end_ctg_len_for_mafft           = 1000
gap_N_num                       = 50


'''
seqtk subseq cami_hc_combined_R1.fasta rd1_r1_to_extract.txt > rd1_extracted_all_r1.fasta
seqtk subseq cami_hc_combined_R2.fasta rd1_r2_to_extract.txt > rd1_extracted_all_r2.fasta

seqtk subseq cami_hc_combined_R1.fasta rd1_r1_to_extract.txt > rd1_extracted_r1.fasta
seqtk subseq cami_hc_combined_R2.fasta rd1_r2_to_extract.txt > rd1_extracted_r2.fasta
reformat.sh in=rd1_extracted_to_gnm.sam out=rd1_extracted_to_gnm_reformatted.sam sam=1.4 2> rd1_extracted_to_gnm_reformat.log

seqtk subseq cami_hc_combined_R1.fasta rd1_linking_reads_R1.txt > rd1_linking_reads_R1.fasta
seqtk subseq cami_hc_combined_R2.fasta rd1_linking_reads_R2.txt > rd1_linking_reads_R2.fasta

'''

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

n = 0
r1_to_extract = set()
r2_to_extract = set()
for each_mp in MappingRecord_dict:

    current_mp_record = MappingRecord_dict[each_mp]
    refs_to_ignore = set()

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
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(r1_16s_ref_cigar_splitted)
                if r1_mismatch_pct > mismatch_cutoff:
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
                            r1_16s_refs_passed_qc[r1_16s_ref] = r1_16s_ref_cigar

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
                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(r2_16s_ref_cigar_splitted)
                if r2_mismatch_pct > mismatch_cutoff:
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
                            r2_16s_refs_passed_qc[r2_16s_ref] = r2_16s_ref_cigar


    ####################################################################################################

    #print('%s\tr1_16s_ref_dict\t%s'       % (each_mp, current_mp_record.r1_16s_ref_dict))
    #print('%s\tr1_16s_refs_passed_qc\t%s' % (each_mp, r1_16s_refs_passed_qc))
    #print('%s\tr1_16s_refs_to_ignore\t%s' % (each_mp, r1_16s_refs_to_ignore))
    #print(len(current_mp_record.r1_16s_ref_dict) == len(r1_16s_refs_passed_qc) + len(r1_16s_refs_to_ignore))

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
        r2_to_extract.add('%s.2' % each_mp)

        # if the shortest S longer than min_M_ctg_len, also consider as clipping mapped
        r1_16s_refs_no_ignored_cigar_list = list(r1_16s_refs_no_ignored.values())
        r1_min_cigar_S, r1_max_cigar_S = get_min_max_cigar_S(r1_16s_refs_no_ignored_cigar_list)
        if r1_min_cigar_S >= min_M_len_ctg:
            MappingRecord_dict[each_mp].consider_r1_clipping_part = True
            r1_to_extract.add('%s.1' % each_mp)

    # only r2 has no_ignored alignments
    elif (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) > 0):

        MappingRecord_dict[each_mp].qualified_reads = True
        MappingRecord_dict[each_mp].consider_r2_unmapped_mate = True
        MappingRecord_dict[each_mp].r2_16s_refs_no_ignored = r2_16s_refs_no_ignored
        r1_to_extract.add('%s.1' % each_mp)

        # if the shortest S longer than min_M_ctg_len, also consider as clipping mapped
        r2_16s_refs_no_ignored_cigar_list = list(r2_16s_refs_no_ignored.values())
        r2_min_cigar_S, r2_max_cigar_S = get_min_max_cigar_S(r2_16s_refs_no_ignored_cigar_list)
        if r2_min_cigar_S >= min_M_len_ctg:
            MappingRecord_dict[each_mp].consider_r2_clipping_part = True
            r2_to_extract.add('%s.2' % each_mp)

    # both r1 and r2 have no_ignored alignments
    else:
        shared_16s_ref_set = set(r1_16s_refs_no_ignored).intersection(set(r2_16s_refs_no_ignored))
        if len(shared_16s_ref_set) > 0:
            MappingRecord_dict[each_mp].qualified_reads = True
            MappingRecord_dict[each_mp].both_mapped_to_16s = True
            MappingRecord_dict[each_mp].shared_16s_refs_no_ignored = shared_16s_ref_set
            r1_to_extract.add('%s.1' % each_mp)
            r2_to_extract.add('%s.2' % each_mp)

        # else:
        #     r1_16s_refs_no_ignored_shared = {key: value for key, value in r1_16s_refs_no_ignored.items() if key in shared_16s_ref_set}
        #     r2_16s_refs_no_ignored_shared = {key: value for key, value in r2_16s_refs_no_ignored.items() if key in shared_16s_ref_set}
        #     #print('%s\tr1_16s_ref_dict\t%s' % (each_mp, current_mp_record.r1_16s_ref_dict))
        #     #print('%s\tr2_16s_ref_dict\t%s' % (each_mp, current_mp_record.r2_16s_ref_dict))
        #     #print()
        #     print('%s\tr1_16s_refs_no_ignored\t%s' % (each_mp, r1_16s_refs_no_ignored))
        #     print('%s\tr2_16s_refs_no_ignored\t%s' % (each_mp, r2_16s_refs_no_ignored))
        #     print()
        #     print(list(r1_16s_refs_no_ignored.values()))
        #     print(list(r2_16s_refs_no_ignored.values()))
        #     print()
        #     #print('%s\tr1_16s_refs_no_ignored_shared\t%s' % (each_mp, r1_16s_refs_no_ignored_shared))
        #     #print('%s\tr2_16s_refs_no_ignored_shared\t%s' % (each_mp, r2_16s_refs_no_ignored_shared))
        #     print('%s\tshared_ref_set\t%s' % (each_mp, shared_16s_ref_set))
        #     print()
        #     #print('%s\trefs_to_ignore\t%s' % (each_mp, refs_to_ignore))
        #     print()
        #     print('----------------------')
        #     print()
        #     n += 1
        #     # what if one fully matched, one partially mapped?

print('n: %s' % n)



# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    if MappingRecord_dict[each_mp].qualified_reads is False:
        MappingRecord_dict.pop(each_mp)

# write out id of reads to extract
with open(rd1_r1_to_extract, 'w') as rd1_r1_to_extract_handle:
    rd1_r1_to_extract_handle.write('%s\n' % '\n'.join(sorted([i for i in r1_to_extract])))
with open(rd1_r2_to_extract, 'w') as rd1_r2_to_extract_handle:
    rd1_r2_to_extract_handle.write('%s\n' % '\n'.join(sorted([i for i in r2_to_extract])))

# extract reads with seqtk
seqtk_extract_cmd_rd1_r1 = 'seqtk subseq %s %s > %s' % (input_r1_fasta, rd1_r1_to_extract, rd1_extracted_all_r1)
seqtk_extract_cmd_rd1_r2 = 'seqtk subseq %s %s > %s' % (input_r2_fasta, rd1_r2_to_extract, rd1_extracted_all_r2)
print(seqtk_extract_cmd_rd1_r1)
print(seqtk_extract_cmd_rd1_r2)
#os.system(seqtk_extract_cmd_rd1_r1)
#os.system(seqtk_extract_cmd_rd1_r2)

# read extracted read sequences into dict
extract_read_seq_dict = {}
for extracted_r1 in SeqIO.parse(rd1_extracted_all_r1, 'fasta'):
    extract_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
for extracted_r2 in SeqIO.parse(rd1_extracted_all_r2, 'fasta'):
    extract_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)

# write out paired and unpaired reads separately
reads_to_extract_base_set = {'.'.join(i.split('.')[:-1]) for i in r1_to_extract.union(r2_to_extract)}
rd1_extracted_p_r1_handle = open(rd1_extracted_p_r1, 'w')
rd1_extracted_p_r2_handle = open(rd1_extracted_p_r2, 'w')
rd1_extracted_up_handle = open(rd1_extracted_up, 'w')
rd1_extracted_paired = set()
rd1_extracted_unpaired = set()
for each_read_base in sorted(list(reads_to_extract_base_set)):
    current_r1 = '%s.1' % each_read_base
    current_r2 = '%s.2' % each_read_base
    current_r1_seq = extract_read_seq_dict.get(current_r1, '')
    current_r2_seq = extract_read_seq_dict.get(current_r2, '')

    if (current_r1 in r1_to_extract) and (current_r2 in r2_to_extract):
        rd1_extracted_p_r1_handle.write('>%s\n' % current_r1)
        rd1_extracted_p_r1_handle.write('%s\n'  % current_r1_seq)
        rd1_extracted_p_r2_handle.write('>%s\n' % current_r2)
        rd1_extracted_p_r2_handle.write('%s\n'  % current_r2_seq)
        rd1_extracted_paired.add(each_read_base)

    if (current_r1 in r1_to_extract) and (current_r2 not in r2_to_extract):
        rd1_extracted_up_handle.write('>%s\n' % current_r1)
        rd1_extracted_up_handle.write('%s\n'  % current_r1_seq)
        rd1_extracted_unpaired.add(current_r1)

    if (current_r1 not in r1_to_extract) and (current_r2 in r2_to_extract):
        rd1_extracted_up_handle.write('>%s\n' % current_r2)
        rd1_extracted_up_handle.write('%s\n'  % current_r2_seq)
        rd1_extracted_unpaired.add(current_r2)

rd1_extracted_p_r1_handle.close()
rd1_extracted_p_r2_handle.close()
rd1_extracted_up_handle.close()


############################# map extracted reads to combined input genomes #############################

bowtie_cmd_rd1_extracted_to_mag = ''
if (len(rd1_extracted_paired) > 0) and (len(rd1_extracted_unpaired) > 0):
    bowtie_cmd_rd1_extracted_to_mag = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p %s -f --local --all --no-unal 2> %s' % (combined_input_gnms_no_ext, rd1_extracted_p_r1, rd1_extracted_p_r2, rd1_extracted_up, rd1_extracted_to_gnm_sam, num_threads, rd1_extracted_to_gnm_sam_log)
if (len(rd1_extracted_paired) > 0) and (len(rd1_extracted_unpaired) == 0):
    bowtie_cmd_rd1_extracted_to_mag = 'bowtie2 -x %s -1 %s -2 %s -S %s -p %s -f --local --all --no-unal 2> %s'       % (combined_input_gnms_no_ext, rd1_extracted_p_r1, rd1_extracted_p_r2, rd1_extracted_to_gnm_sam, num_threads, rd1_extracted_to_gnm_sam_log)
if (len(rd1_extracted_paired) == 0) and (len(rd1_extracted_unpaired) > 0):
    bowtie_cmd_rd1_extracted_to_mag = 'bowtie2 -x %s -U %s -S %s -p %s -f --local --all --no-unal 2> %s'             % (combined_input_gnms_no_ext, rd1_extracted_up, rd1_extracted_to_gnm_sam, num_threads, rd1_extracted_to_gnm_sam_log)

print(bowtie_cmd_rd1_extracted_to_mag)
#os.system(bowtie_cmd_rd1_extracted_to_mag)

bbmap_reformat_cmd_rd1_extracted_to_mag = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (rd1_extracted_to_gnm_sam, rd1_extracted_to_gnm_sam_reformatted, rd1_extracted_to_gnm_sam_reformat_log)
print(bbmap_reformat_cmd_rd1_extracted_to_mag)
#os.system(bbmap_reformat_cmd_rd1_extracted_to_mag)

keep_best_matches_in_sam_keep_short_M(rd1_extracted_to_gnm_sam_reformatted, min_M_len_ctg, rd1_extracted_to_gnm_sam_reformatted_best_match)


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

m = 0
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
                            r1_ctg_refs_passed_qc[r1_ctg_ref] = r1_ctg_ref_cigar

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
                            r2_ctg_refs_passed_qc[r2_ctg_ref] = r2_ctg_ref_cigar

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
        shared_ctg_ref_set = set(r1_ctg_refs_no_ignored).intersection(set(r2_ctg_refs_no_ignored))
        if len(shared_ctg_ref_set) > 0:
            MappingRecord_dict[each_mp].matched_to_ctg = True
            MappingRecord_dict[each_mp].shared_ctg_refs_no_ignored = shared_ctg_ref_set


##################################################### get linkages from MappingRecord_dict #####################################################

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

    combined_16s_refs = set.union(set(r1_16s_refs), set(r2_16s_refs), set(shared_16s_refs))
    combined_ctg_refs = set.union(set(r1_ctg_refs), set(r2_ctg_refs), set(shared_ctg_refs))

    for each_16s_ref in combined_16s_refs:
        for each_ctg_ref in combined_ctg_refs:
            marker_to_ctg_key = '%s%s%s' % (each_16s_ref, marker_to_ctg_gnm_Key_connector, each_ctg_ref)

            if marker_to_ctg_key not in marker_to_ctg_linkage_num_dict:
                marker_to_ctg_linkage_num_dict[marker_to_ctg_key] = 1
                marker_to_ctg_linking_reads_dict[marker_to_ctg_key] = {qualified_read}
            else:
                marker_to_ctg_linkage_num_dict[marker_to_ctg_key] += 1
                marker_to_ctg_linking_reads_dict[marker_to_ctg_key].add(qualified_read)

# remove linkages less than 3
marker_to_ctg_linkage_num_min3_dict = {}
for each_key in marker_to_ctg_linkage_num_dict:
    if marker_to_ctg_linkage_num_dict[each_key] >= 3:
        marker_to_ctg_linkage_num_min3_dict[each_key] = marker_to_ctg_linkage_num_dict[each_key]

# get number of linkages at genome level
marker_to_gnm_link_num = {}
for each_marker_to_ctg_key in marker_to_ctg_linkage_num_min3_dict:
    marker_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[0]
    ctg_id = each_marker_to_ctg_key.split(marker_to_ctg_gnm_Key_connector)[1]
    gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
    marker_to_gnm_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, gnm_id)
    if marker_to_gnm_key not in marker_to_gnm_link_num:
        marker_to_gnm_link_num[marker_to_gnm_key] = marker_to_ctg_linkage_num_min3_dict[each_marker_to_ctg_key]
    else:
        marker_to_gnm_link_num[marker_to_gnm_key] += marker_to_ctg_linkage_num_min3_dict[each_marker_to_ctg_key]

sankey_file_in_handle = open(link_stats_combined, 'w')
sankey_file_in_handle.write('MarkerGene,GenomicSeq,Number\n')
for each_linkage in marker_to_gnm_link_num:
    sankey_file_in_handle.write('MarkerGene__%s,GenomicSeq__%s,%s\n' % (
    each_linkage.split(marker_to_ctg_gnm_Key_connector)[0], each_linkage.split(marker_to_ctg_gnm_Key_connector)[1],
    marker_to_gnm_link_num[each_linkage]))
sankey_file_in_handle.close()


########################################################################################################################
######################################### get linking reads for visualization ##########################################
########################################################################################################################

if os.path.isdir(mafft_seq_folder) is True:
    os.system('rm -r %s' % mafft_seq_folder)
os.mkdir(mafft_seq_folder)

ctgs_to_extract = set()
all_linking_reads_base_set = set()
linking_reads_txt_handle = open(linking_reads_tab, 'w')
for marker_to_ctg in marker_to_ctg_linkage_num_min3_dict:
    marker_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[0]
    ctg_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[1]
    gnm_id = ctg_id.split(gnm_to_ctg_connector)[0]
    linking_reads = marker_to_ctg_linking_reads_dict[marker_to_ctg]
    all_linking_reads_base_set.update(linking_reads)
    ctgs_to_extract.add(ctg_id)
    linking_reads_txt_handle.write('%s\t%s\t%s\n' % (marker_id, ctg_id, ','.join(linking_reads)))
linking_reads_txt_handle.close()

# write out id of linking reads for extraction
with open(linking_reads_r1_txt, 'w') as linking_reads_r1_txt_handle:
    linking_reads_r1_txt_handle.write('%s\n' % '\n'.join(sorted([('%s.1' % i) for i in all_linking_reads_base_set])))
with open(linking_reads_r2_txt, 'w') as linking_reads_r2_txt_handle:
    linking_reads_r2_txt_handle.write('%s\n' % '\n'.join(sorted([('%s.2' % i) for i in all_linking_reads_base_set])))

# extract linking reads with seqtk
seqtk_extract_cmd_rd1_linking_r1 = 'seqtk subseq %s %s > %s' % (input_r1_fasta, linking_reads_r1_txt, linking_reads_r1_fasta)
seqtk_extract_cmd_rd1_linking_r2 = 'seqtk subseq %s %s > %s' % (input_r2_fasta, linking_reads_r2_txt, linking_reads_r2_fasta)
print(seqtk_extract_cmd_rd1_linking_r1)
print(seqtk_extract_cmd_rd1_linking_r2)
#os.system(seqtk_extract_cmd_rd1_linking_r1)
#os.system(seqtk_extract_cmd_rd1_linking_r2)

# subset combined genome file
linked_contigs_txt_handle = open(linked_contigs_txt, 'w')
linked_contigs_txt_handle.write('\n'.join(ctgs_to_extract) + '\n')
linked_contigs_txt_handle.close()
subset_linked_ctgs_cmd = 'seqtk subseq %s %s > %s' % (combined_input_gnms, linked_contigs_txt, linked_contigs_fasta)
os.system(subset_linked_ctgs_cmd)

# read sequence of linked contigs into dict
linked_ctg_seq_dict = {}
for linked_ctg in SeqIO.parse(linked_contigs_fasta, 'fasta'):
    linked_ctg_seq_dict[linked_ctg.id] = str(linked_ctg.seq)

# read sequence of linking reads into dict
linking_read_seq_dict = {}
for linking_r1 in SeqIO.parse(linking_reads_r1_fasta, 'fasta'):
    linking_read_seq_dict[linking_r1.id] = str(linking_r1.seq)
for linking_r2 in SeqIO.parse(linking_reads_r2_fasta, 'fasta'):
    linking_read_seq_dict[linking_r2.id] = str(linking_r2.seq)

# read marker sequences into dict
marker_seq_dict = {}
for each_16s in SeqIO.parse(marker_gene_seqs, 'fasta'):
    marker_seq_dict[each_16s.id] = str(each_16s.seq)

# write out sequences for each linkage
concatenate_dict = {}
concatenated_ref_id_dict = {}
concatenate_pos_dict = {}
for marker_to_ctg in marker_to_ctg_linkage_num_min3_dict:
    linking_reads = marker_to_ctg_linking_reads_dict[marker_to_ctg]
    marker_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[0]
    marker_seq = marker_seq_dict[marker_id]
    ctg_id = marker_to_ctg.split(marker_to_ctg_gnm_Key_connector)[1]
    contig_seq = linked_ctg_seq_dict[ctg_id]
    reads_file_base_tmp = marker_to_ctg.replace(marker_to_ctg_gnm_Key_connector, '___')
    reads_file_base = reads_file_base_tmp.replace(gnm_to_ctg_connector, '___')
    os.mkdir('%s/%s' % (mafft_seq_folder, reads_file_base))

    linked_16s_seq_file = '%s/%s/%s.fa'    % (mafft_seq_folder, reads_file_base, marker_id)
    reads_file_r1       = '%s/%s/%s_R1.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)
    reads_file_r2       = '%s/%s/%s_R2.fa' % (mafft_seq_folder, reads_file_base, reads_file_base)

    # write out sequences of linking reads
    marker_pos_list = []
    contig_pos_list = []
    reads_file_r1_handle = open(reads_file_r1, 'w')
    reads_file_r2_handle = open(reads_file_r2, 'w')
    for each_linking_read in linking_reads:
        linking_r1_id  = '%s.1' % each_linking_read
        linking_r1_seq = linking_read_seq_dict[linking_r1_id]
        linking_r2_id  = '%s.2' % each_linking_read
        linking_r2_seq = linking_read_seq_dict[linking_r2_id]
        reads_file_r1_handle.write('>%s\n' % linking_r1_id)
        reads_file_r1_handle.write('%s\n' % linking_r1_seq)
        reads_file_r2_handle.write('>%s\n' % linking_r2_id)
        reads_file_r2_handle.write('%s\n' % linking_r2_seq)

        # get matched position on makrer and contig
        marker_pos_r1 = list(MappingRecord_dict[each_linking_read].r1_16s_ref_dict.get(marker_id, dict()).keys())
        marker_pos_r2 = list(MappingRecord_dict[each_linking_read].r2_16s_ref_dict.get(marker_id, dict()).keys())
        contig_pos_r1 = list(MappingRecord_dict[each_linking_read].r1_ctg_ref_dict.get(ctg_id, dict()).keys())
        contig_pos_r2 = list(MappingRecord_dict[each_linking_read].r2_ctg_ref_dict.get(ctg_id, dict()).keys())
        if len(marker_pos_r1) == 1:
            marker_pos_list.append(marker_pos_r1[0])
        if len(marker_pos_r2) == 1:
            marker_pos_list.append(marker_pos_r2[0])
        if len(contig_pos_r1) == 1:
            contig_pos_list.append(contig_pos_r1[0])
        if len(contig_pos_r2) == 1:
            contig_pos_list.append(contig_pos_r2[0])

    reads_file_r1_handle.close()
    reads_file_r2_handle.close()

    # get marker linked end
    marker_pos_median = np.median(marker_pos_list)
    linked_end_marker = 'middle'
    if marker_pos_median <= (marker_len_dict[marker_id] / 3):
        linked_end_marker = 'left'
    elif marker_pos_median >= (marker_len_dict[marker_id] * 2 / 3):
        linked_end_marker = 'right'

    # get contig linked end
    contig_pos_median = np.median(contig_pos_list)
    contig_pos_middle = int(round(float(contig_pos_median)))
    linked_end_contig = 'middle(%s)' % contig_pos_middle
    if contig_pos_median <= 500:
        linked_end_contig = 'left'
    elif (ctg_len_dict[ctg_id] - contig_pos_median) <= 500:
        linked_end_contig = 'right'

    # get contig end sequence
    if linked_end_contig == 'left':
        contig_seq_for_mafft = contig_seq[:end_ctg_len_for_mafft]
    elif linked_end_contig == 'right':
        contig_seq_for_mafft = contig_seq[-end_ctg_len_for_mafft:]
    else:
        left_end_pos = 0
        if (contig_pos_middle - end_ctg_len_for_mafft) > 0:
            left_end_pos = contig_pos_middle - end_ctg_len_for_mafft
        right_end_pos = contig_pos_middle + end_ctg_len_for_mafft
        if right_end_pos > ctg_len_dict[ctg_id]:
            right_end_pos = ctg_len_dict[ctg_id]
        contig_seq_for_mafft = contig_seq[left_end_pos:(right_end_pos - 1)]

    # concatenate 16s and contig sequences
    to_concatenate = False
    concatenated_seq_id = ''
    concatenated_seq = ''
    concatenate_pos = 0
    if linked_end_contig in ['left', 'right']:
        if (linked_end_marker == 'right') and (linked_end_contig == 'left'):
            concatenated_seq_id = 'Marker___NNN___Contig'
            concatenated_seq = '%s%s%s' % (marker_seq, 'N'*gap_N_num, contig_seq_for_mafft)
            to_concatenate = True
            concatenate_pos = len(marker_seq) + round(gap_N_num/2)

        if (linked_end_marker == 'left') and (linked_end_contig == 'right'):
            concatenated_seq_id = 'Contig___NNN___Marker'
            concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N'*gap_N_num, marker_seq)
            to_concatenate = True
            concatenate_pos = len(contig_seq_for_mafft) + round(gap_N_num/2)

        if (linked_end_marker == 'left') and (linked_end_contig == 'left'):
            marker_seq_rc = get_rc(marker_seq)
            concatenated_seq_id = 'Marker_RC___NNN___Contig'
            concatenated_seq = '%s%s%s' % (marker_seq_rc, 'N'*gap_N_num, contig_seq_for_mafft)
            to_concatenate = True
            concatenate_pos = len(marker_seq_rc) + round(gap_N_num/2)

        if (linked_end_marker == 'right') and (linked_end_contig == 'right'):
            marker_seq_rc = get_rc(marker_seq)
            concatenated_seq_id = 'Contig___NNN___Marker_RC'
            concatenated_seq = '%s%s%s' % (contig_seq_for_mafft, 'N'*gap_N_num, marker_seq_rc)
            to_concatenate = True
            concatenate_pos = len(contig_seq_for_mafft) + round(gap_N_num/2)

    concatenate_dict[reads_file_base] = to_concatenate
    concatenated_ref_id_dict[reads_file_base] = concatenated_seq_id
    concatenate_pos_dict[reads_file_base] = concatenate_pos

    # write out sequences
    pwd_seq_file_cbd = '%s/%s/%s_cbd.fasta' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_16s = '%s/%s/%s_16s.fasta' % (mafft_seq_folder, reads_file_base, reads_file_base)
    pwd_seq_file_ctg = '%s/%s/%s_ctg.fasta' % (mafft_seq_folder, reads_file_base, reads_file_base)
    if to_concatenate is True:
        pwd_seq_file_cbd_handle = open(pwd_seq_file_cbd, 'w')
        pwd_seq_file_cbd_handle.write('>%s\n' % concatenated_seq_id)
        pwd_seq_file_cbd_handle.write('%s\n' % concatenated_seq)
        pwd_seq_file_cbd_handle.close()
    else:
        # write out 16s sequence
        pwd_seq_file_16s_handle = open(pwd_seq_file_16s, 'w')
        pwd_seq_file_16s_handle.write('>Marker\n')
        pwd_seq_file_16s_handle.write('%s\n' % marker_seq)
        pwd_seq_file_16s_handle.close()
        # write out ctg sequence
        pwd_seq_file_ctg_handle = open(pwd_seq_file_ctg, 'w')
        pwd_seq_file_ctg_handle.write('>Contig\n')
        pwd_seq_file_ctg_handle.write('%s\n' % contig_seq_for_mafft)
        pwd_seq_file_ctg_handle.close()


########## align and visualize ##########

# prepare arguments for mapping_worker
list_for_mapping_worker = []
for each_marker_to_ctg in concatenate_dict:
    list_for_mapping_worker.append([mafft_seq_folder, each_marker_to_ctg, concatenate_dict[each_marker_to_ctg], concatenated_ref_id_dict[each_marker_to_ctg], concatenate_pos_dict.get(each_marker_to_ctg, 0), mafft_seq_folder_local])

# run mapping_worker with multiprocessing
pool = mp.Pool(processes=num_threads)
pool.map(mapping_worker, list_for_mapping_worker)
pool.close()
pool.join()
