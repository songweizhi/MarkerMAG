import os
import glob
from Bio import SeqIO
import multiprocessing as mp


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


def parse_sam_gnm_worker(arguments_list):
    rd1_unlinked_mags_sam_bowtie_reformat_sorted = arguments_list[0]
    free_living_ctg_ref_file = arguments_list[1]
    min_M_len_ctg = arguments_list[2]
    mismatch_cutoff = arguments_list[3]
    round_2_ctg_end_seq_len_dict = arguments_list[4]

    free_living_ctg_ref_file_handle = open(free_living_ctg_ref_file, 'w')
    to_extract_read_base_rd2_ctg = set()
    current_read_base = ''
    current_read_base_r1_ctg_ref_dict_rd2 = dict()
    current_read_base_r2_ctg_ref_dict_rd2 = dict()
    with open(rd1_unlinked_mags_sam_bowtie_reformat_sorted) as rd1_unlinked_mags_sam_bowtie_reformat_sorted_opened:
        for each_line in rd1_unlinked_mags_sam_bowtie_reformat_sorted_opened:
            if not each_line.startswith('@'):
                each_line_split = each_line.strip().split('\t')
                cigar = each_line_split[5]
                read_id = each_line_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_line_split[2]
                ref_pos = int(each_line_split[3])

                if current_read_base == '':
                    current_read_base = read_id_base

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}

                elif read_id_base == current_read_base:

                    if cigar != '*':
                        if read_strand == '1':
                            if ref_id not in current_read_base_r1_ctg_ref_dict_rd2:
                                current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r1_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar
                        if read_strand == '2':
                            if ref_id not in current_read_base_r2_ctg_ref_dict_rd2:
                                current_read_base_r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r2_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar
                else:
                    ################################### analysis previous read refs ####################################

                    ctg_refs_to_ignore_rd2 = set()

                    ########## get lowest mismatch for r1/r2 ctg refs ##########

                    # get r1_ref_cigar_set
                    r1_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r1_ctg_ref_dict_rd2.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r1_ref_cigar_set.update(each_pos_dict_values)

                    # get r2_ref_cigar_set
                    r2_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r2_ctg_ref_dict_rd2.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r2_ref_cigar_set.update(each_pos_dict_values)

                    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_ctg)
                    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_ctg)

                    ########## filter r1 ctg refs ##########
                    r1_ctg_refs_passed_qc = {}
                    for r1_ctg_ref_rd2 in current_read_base_r1_ctg_ref_dict_rd2:
                        r1_matched_pos_dict = current_read_base_r1_ctg_ref_dict_rd2[r1_ctg_ref_rd2]
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
                                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                    r1_ctg_ref_cigar_splitted)
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
                                            if (r1_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r1_ctg_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r1_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r1_ctg_ref_pos + r1_aligned_len - 1) == \
                                                        round_2_ctg_end_seq_len_dict[
                                                            r1_ctg_ref_rd2]:
                                                    clip_in_middle = False

                                        if clip_in_middle is True:
                                            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                        else:
                                            r1_ctg_refs_passed_qc[r1_ctg_ref_rd2] = [r1_ctg_ref_cigar]

                    ########## filter r2 ctg refs ##########
                    r2_ctg_refs_passed_qc = {}
                    for r2_ctg_ref_rd2 in current_read_base_r2_ctg_ref_dict_rd2:
                        r2_matched_pos_dict = current_read_base_r2_ctg_ref_dict_rd2[r2_ctg_ref_rd2]
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
                                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                    r2_ctg_ref_cigar_splitted)
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
                                            if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r2_ctg_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r2_ctg_ref_pos + r2_aligned_len - 1) == \
                                                        round_2_ctg_end_seq_len_dict[r2_ctg_ref_rd2]:
                                                    clip_in_middle = False

                                        if clip_in_middle is True:
                                            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                        else:
                                            r2_ctg_refs_passed_qc[r2_ctg_ref_rd2] = [r2_ctg_ref_cigar]

                    ####################################################################################################

                    r1_ctg_refs_rd2_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if
                                                  key not in ctg_refs_to_ignore_rd2}
                    r2_ctg_refs_rd2_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if
                                                  key not in ctg_refs_to_ignore_rd2}

                    # only r1 has no_ignored alignments
                    if (len(r1_ctg_refs_rd2_no_ignored) > 0) and (len(r2_ctg_refs_rd2_no_ignored) == 0):
                        current_read_base___qualified_reads_rd2 = True
                        current_read_base___r1_ctg_ref_dict_rd2 = current_read_base_r1_ctg_ref_dict_rd2
                        to_extract_read_base_rd2_ctg.add(current_read_base)
                        free_living_ctg_ref_file_handle.write(
                            '%s\t%s\n' % (current_read_base, ','.join(r1_ctg_refs_rd2_no_ignored)))

                    # only r2 has no_ignored alignments
                    if (len(r1_ctg_refs_rd2_no_ignored) == 0) and (len(r2_ctg_refs_rd2_no_ignored) > 0):
                        current_read_base___qualified_reads_rd2 = True
                        current_read_base___r2_ctg_ref_dict_rd2 = current_read_base_r2_ctg_ref_dict_rd2
                        to_extract_read_base_rd2_ctg.add(current_read_base)
                        free_living_ctg_ref_file_handle.write(
                            '%s\t%s\n' % (current_read_base, ','.join(r2_ctg_refs_rd2_no_ignored)))

                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_ctg_ref_dict_rd2 = dict()
                    current_read_base_r2_ctg_ref_dict_rd2 = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}

    free_living_ctg_ref_file_handle.close()


####################################################################################################################

if __name__ == '__main__':


    wd = '/Users/songweizhi/Desktop/tunning_rd2'
    combined_1st_round_unlinked_mag_end_seq      = '%s/round_1_unlinked_gnm_end_500bp.fa' % wd
    rd1_unlinked_mags_sam_bowtie_reformat_sorted = '%s/subset.sam' % wd
    rd1_unlinked_mags_sam_split_folder           = '%s/round_1_unlinked_gnm_bowtie_reformatted_split' % wd
    free_living_ctg_ref_file                     = '%s/round2_free_living_ctg_refs.txt' % wd
    free_living_ctg_ref_folder_split              = '%s/round2_free_living_ctg_refs_split' % wd
    line_num_per_file = 10000
    min_M_len_ctg = 30
    mismatch_cutoff = 2
    num_threads = 4
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
    for each_ctg_end_record in SeqIO.parse(combined_1st_round_unlinked_mag_end_seq, 'fasta'):
        round_2_ctg_end_seq_len_dict[each_ctg_end_record.id] = len(each_ctg_end_record.seq)

    if os.path.isdir(rd1_unlinked_mags_sam_split_folder) is True:
        os.system('rm -r %s' % rd1_unlinked_mags_sam_split_folder)
    os.mkdir(rd1_unlinked_mags_sam_split_folder)

    split_gnm_sam_cmd = 'split -l %s %s %s/splitted_sam_ ' % (line_num_per_file, rd1_unlinked_mags_sam_bowtie_reformat_sorted, rd1_unlinked_mags_sam_split_folder)
    os.system(split_gnm_sam_cmd)

    # get splitted sam file list
    splitted_gnm_sam_file_list = [os.path.basename(file_name) for file_name in glob.glob('%s/splitted_sam_*' % rd1_unlinked_mags_sam_split_folder)]


    ####################################################################################################################

    if os.path.isdir(free_living_ctg_ref_folder_split) is True:
        os.system('rm -r %s' % free_living_ctg_ref_folder_split)
    os.mkdir(free_living_ctg_ref_folder_split)

    # prepare lol for mp worker
    list_for_parse_sam_gnm_worker = []
    splitted_sam_mp_file_set = set()
    for splitted_gnm_sam_file in splitted_gnm_sam_file_list:
        pwd_splitted_gnm_sam_file                       = '%s/%s'                           % (rd1_unlinked_mags_sam_split_folder, splitted_gnm_sam_file)
        pwd_splitted_gnm_sam_free_living_ctg_ref_file   = '%s/%s_free_living_ctg_refs.txt'  % (free_living_ctg_ref_folder_split, splitted_gnm_sam_file)
        splitted_sam_mp_file_set.add(pwd_splitted_gnm_sam_free_living_ctg_ref_file)
        list_for_parse_sam_gnm_worker.append([pwd_splitted_gnm_sam_file,
                                             pwd_splitted_gnm_sam_free_living_ctg_ref_file,
                                             min_M_len_ctg,
                                             mismatch_cutoff,
                                             round_2_ctg_end_seq_len_dict])

    pool_parse_sam_gnm = mp.Pool(processes=num_threads)
    pool_parse_sam_gnm.map(parse_sam_gnm_worker, list_for_parse_sam_gnm_worker)
    pool_parse_sam_gnm.close()
    pool_parse_sam_gnm.join()

    # combine free_living_ctg_ref_files
    os.system('cat %s > %s' % (' '.join(splitted_sam_mp_file_set), free_living_ctg_ref_file))

    ####################################################################################################


    # # extract reads with seqtk
    # seqtk_extract_cmd_rd1_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_to_extract_flking_ctg_r1_id, rd2_extracted_flking_ctg_r1_seq)
    # seqtk_extract_cmd_rd1_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_to_extract_flking_ctg_r2_id, rd2_extracted_flking_ctg_r2_seq)
    # #report_and_log((seqtk_extract_cmd_rd1_r1), pwd_log_file, True)
    # #report_and_log((seqtk_extract_cmd_rd1_r2), pwd_log_file, True)
    # #os.system(seqtk_extract_cmd_rd1_r1)
    # #os.system(seqtk_extract_cmd_rd1_r2)
    #
    # # read extracted read sequences into dict
    # extract_rd2_flking_ctg_read_seq_dict = {}
    # for extracted_r1 in SeqIO.parse(rd2_extracted_flking_ctg_r1_seq_tmp, 'fasta'):
    #     extract_rd2_flking_ctg_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
    # for extracted_r2 in SeqIO.parse(rd2_extracted_flking_ctg_r2_seq_tmp, 'fasta'):
    #     extract_rd2_flking_ctg_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)
    #
    # # write out paired in the same order
    # to_extract_read_base_rd2_ctg = {}
    # rd2_extracted_flking_ctg_r1_handle = open(rd2_extracted_flking_ctg_r1_seq, 'w')
    # rd2_extracted_flking_ctg_r2_handle = open(rd2_extracted_flking_ctg_r2_seq, 'w')
    # for each_read_base in to_extract_read_base_rd2_ctg:
    #     current_r1 = '%s.1' % each_read_base
    #     current_r2 = '%s.2' % each_read_base
    #     current_r1_seq = extract_rd2_flking_ctg_read_seq_dict.get(current_r1, '')
    #     current_r2_seq = extract_rd2_flking_ctg_read_seq_dict.get(current_r2, '')
    #     rd2_extracted_flking_ctg_r1_handle.write('>%s\n' % current_r1)
    #     rd2_extracted_flking_ctg_r1_handle.write('%s\n' % current_r1_seq)
    #     rd2_extracted_flking_ctg_r2_handle.write('>%s\n' % current_r2)
    #     rd2_extracted_flking_ctg_r2_handle.write('%s\n' % current_r2_seq)
    # rd2_extracted_flking_ctg_r1_handle.close()
    # rd2_extracted_flking_ctg_r2_handle.close()
    #

