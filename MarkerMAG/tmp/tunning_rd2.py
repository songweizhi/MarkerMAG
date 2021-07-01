import os
from Bio import SeqIO


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


####################################################################################################################

rd1_unlinked_mags_sam_bowtie_reformat = '/Users/songweizhi/Desktop/tunning_rd2/round_1_unlinked_gnm_bowtie_reformatted.sam'
min_M_len_ctg = 30
mismatch_cutoff = 2
report_interval = 100000
'''
# should keep both, but only kept one !!!
MBARC26_13313152	{'TC___C___NODE_809_length_22129_cov_127.040953_l': {1: '66S55=1X8='}, 'TC___C___NODE_189_length_107521_cov_122.830635_l': {1: '76S45=1X8='}}
MBARC26_13313152.1	TC___C___NODE_809_length_22129_cov_127.040953_l	{1: '66S55=1X8='}
'''

reads_file_r1_fasta = ''
reads_file_r2_fasta = ''
rd2_to_extract_flking_ctg_r1_id     = 'rd2_to_extract_flking_ctg_r1_id.txt'
rd2_to_extract_flking_ctg_r2_id     = 'rd2_to_extract_flking_ctg_r2_id.txt'
rd2_extracted_flking_ctg_r1_seq_tmp = 'rd2_extracted_flking_ctg_r1_tmp.fa'
rd2_extracted_flking_ctg_r2_seq_tmp = 'rd2_extracted_flking_ctg_r2_tmp.fa'
rd2_extracted_flking_ctg_r1_seq     = 'rd2_extracted_flking_ctg_r1.fa'
rd2_extracted_flking_ctg_r2_seq     = 'rd2_extracted_flking_ctg_r2.fa'


####################################################################################################################

round_2_ctg_end_seq_len_dict = {}
round_2_MappingRecord_dict = {}
for each_line in open(rd1_unlinked_mags_sam_bowtie_reformat):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('@'):
        rd2_ref_id = ''
        rd2_ref_len = 0
        for each_element in each_line_split:
            if each_element.startswith('SN:'):
                rd2_ref_id = each_element[3:]
            if each_element.startswith('LN:'):
                rd2_ref_len = int(each_element[3:])
        round_2_ctg_end_seq_len_dict[rd2_ref_id] = rd2_ref_len
    else:
        cigar = each_line_split[5]
        if cigar != '*':
            read_id = each_line_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_line_split[2]
            ref_pos = int(each_line_split[3])

            if read_id_base not in round_2_MappingRecord_dict:
                round_2_MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                if ref_id not in round_2_MappingRecord_dict[read_id_base].r1_ctg_ref_dict_rd2:
                    round_2_MappingRecord_dict[read_id_base].r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                else:
                    round_2_MappingRecord_dict[read_id_base].r1_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar

            if read_strand == '2':
                if ref_id not in round_2_MappingRecord_dict[read_id_base].r2_ctg_ref_dict_rd2:
                    round_2_MappingRecord_dict[read_id_base].r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                else:
                    round_2_MappingRecord_dict[read_id_base].r2_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar


####################################################################################################

n = 0
processed_num = 0
to_extract_read_base_rd2_ctg = set()
for each_mp in round_2_MappingRecord_dict.copy():

    ctg_mp_record = round_2_MappingRecord_dict[each_mp]
    ctg_refs_to_ignore_rd2 = set()

    ########## get lowest mismatch for r1/r2 ctg refs ##########

    # get r1_ref_cigar_set
    r1_ref_cigar_set = set()
    for each_pos_dict in ctg_mp_record.r1_ctg_ref_dict_rd2.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r1_ref_cigar_set.update(each_pos_dict_values)

    # get r2_ref_cigar_set
    r2_ref_cigar_set = set()
    for each_pos_dict in ctg_mp_record.r2_ctg_ref_dict_rd2.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r2_ref_cigar_set.update(each_pos_dict_values)

    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_ctg)
    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_ctg)
    round_2_MappingRecord_dict[each_mp].r1_ctg_refs_lowest_mismatch_rd2 = r1_ref_min_mismatch
    round_2_MappingRecord_dict[each_mp].r2_ctg_refs_lowest_mismatch_rd2 = r2_ref_min_mismatch

    ########## filter r1 16s refs ##########

    r1_ctg_refs_passed_qc = {}
    for r1_ctg_ref_rd2 in ctg_mp_record.r1_ctg_ref_dict_rd2:
        r1_matched_pos_dict = ctg_mp_record.r1_ctg_ref_dict_rd2[r1_ctg_ref_rd2]
        if len(r1_matched_pos_dict) > 1:
            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
        else:
            r1_ctg_ref_pos = list(r1_matched_pos_dict.keys())[0]
            r1_ctg_ref_cigar = r1_matched_pos_dict[r1_ctg_ref_pos]
            r1_ctg_ref_cigar_splitted = cigar_splitter(r1_ctg_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r1_ctg_ref_cigar_splitted)
            if both_end_clp is True:
                ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
            else:
                # check mismatch
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(r1_ctg_ref_cigar_splitted)
                if r1_ref_min_mismatch == 'NA':
                    ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                    ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                else:
                    # check aligned length
                    if r1_aligned_len < min_M_len_ctg:
                        ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r1_ctg_ref_cigar) or ('s' in r1_ctg_ref_cigar):
                            clip_in_middle = True
                            if (r1_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_ctg_ref_pos == 1):
                                clip_in_middle = False
                            if (r1_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r1_ctg_ref_pos + r1_aligned_len - 1) == round_2_ctg_end_seq_len_dict[r1_ctg_ref_rd2]:
                                    clip_in_middle = False

                        if clip_in_middle is True:
                            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                        else:
                            r1_ctg_refs_passed_qc[r1_ctg_ref_rd2] = [r1_ctg_ref_cigar]

    ########## filter r2 16s refs ##########

    r2_ctg_refs_passed_qc = {}
    for r2_ctg_ref_rd2 in ctg_mp_record.r2_ctg_ref_dict_rd2:
        r2_matched_pos_dict = ctg_mp_record.r2_ctg_ref_dict_rd2[r2_ctg_ref_rd2]
        if len(r2_matched_pos_dict) > 1:
            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
        else:
            r2_ctg_ref_pos = list(r2_matched_pos_dict.keys())[0]
            r2_ctg_ref_cigar = r2_matched_pos_dict[r2_ctg_ref_pos]
            r2_ctg_ref_cigar_splitted = cigar_splitter(r2_ctg_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
            if both_end_clp is True:
                ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
            else:
                # check mismatch
                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(r2_ctg_ref_cigar_splitted)
                if r2_ref_min_mismatch == 'NA':
                    ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                    ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                else:
                    # check aligned length
                    if r2_aligned_len < min_M_len_ctg:
                        ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r2_ctg_ref_cigar) or ('s' in r2_ctg_ref_cigar):
                            clip_in_middle = True
                            if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_ctg_ref_pos == 1):
                                clip_in_middle = False
                            if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r2_ctg_ref_pos + r2_aligned_len - 1) == round_2_ctg_end_seq_len_dict[r2_ctg_ref_rd2]:
                                    clip_in_middle = False

                        if clip_in_middle is True:
                            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                        else:
                            r2_ctg_refs_passed_qc[r2_ctg_ref_rd2] = [r2_ctg_ref_cigar]

    ####################################################################################################

    r1_ctg_refs_rd2_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in ctg_refs_to_ignore_rd2}
    r2_ctg_refs_rd2_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in ctg_refs_to_ignore_rd2}

    # only r1 has no_ignored alignments
    if (len(r1_ctg_refs_rd2_no_ignored) > 0) and (len(r2_ctg_refs_rd2_no_ignored) == 0):
        round_2_MappingRecord_dict[each_mp].qualified_reads_rd2 = True
        to_extract_read_base_rd2_ctg.add(each_mp)

    # only r2 has no_ignored alignments
    if (len(r1_ctg_refs_rd2_no_ignored) == 0) and (len(r2_ctg_refs_rd2_no_ignored) > 0):
        round_2_MappingRecord_dict[each_mp].qualified_reads_rd2 = True
        to_extract_read_base_rd2_ctg.add(each_mp)



# extract reads with seqtk
seqtk_extract_cmd_rd1_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_to_extract_flking_ctg_r1_id, rd2_extracted_flking_ctg_r1_seq)
seqtk_extract_cmd_rd1_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_to_extract_flking_ctg_r2_id, rd2_extracted_flking_ctg_r2_seq)
#report_and_log((seqtk_extract_cmd_rd1_r1), pwd_log_file, True)
#report_and_log((seqtk_extract_cmd_rd1_r2), pwd_log_file, True)
#os.system(seqtk_extract_cmd_rd1_r1)
#os.system(seqtk_extract_cmd_rd1_r2)

# read extracted read sequences into dict
extract_rd2_flking_ctg_read_seq_dict = {}
for extracted_r1 in SeqIO.parse(rd2_extracted_flking_ctg_r1_seq_tmp, 'fasta'):
    extract_rd2_flking_ctg_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
for extracted_r2 in SeqIO.parse(rd2_extracted_flking_ctg_r2_seq_tmp, 'fasta'):
    extract_rd2_flking_ctg_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)

# write out paired in the same order
rd2_extracted_flking_ctg_r1_handle = open(rd2_extracted_flking_ctg_r1_seq, 'w')
rd2_extracted_flking_ctg_r2_handle = open(rd2_extracted_flking_ctg_r2_seq, 'w')
for each_read_base in to_extract_read_base_rd2_ctg:
    current_r1 = '%s.1' % each_read_base
    current_r2 = '%s.2' % each_read_base
    current_r1_seq = extract_rd2_flking_ctg_read_seq_dict.get(current_r1, '')
    current_r2_seq = extract_rd2_flking_ctg_read_seq_dict.get(current_r2, '')
    rd2_extracted_flking_ctg_r1_handle.write('>%s\n' % current_r1)
    rd2_extracted_flking_ctg_r1_handle.write('%s\n' % current_r1_seq)
    rd2_extracted_flking_ctg_r2_handle.write('>%s\n' % current_r2)
    rd2_extracted_flking_ctg_r2_handle.write('%s\n' % current_r2_seq)
rd2_extracted_flking_ctg_r1_handle.close()
rd2_extracted_flking_ctg_r2_handle.close()


