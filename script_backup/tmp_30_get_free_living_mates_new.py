from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


wd          = '/Users/songweizhi/Desktop/round2'
sam_file    = '%s/step_1_unlinked_combined_gnms_ends_1000bp.sam'    % wd
sam_file    = '%s/subset.sam'                                       % wd
round2_fq   = False


free_living_ext = 'fa'
if round2_fq is True:
    free_living_ext = 'fq'

free_living_ctg_R1 = '/Users/songweizhi/Desktop/new_algorithm_Kelp/round2_free_living_ctg_R1.%s' % free_living_ext
free_living_ctg_R2 = '/Users/songweizhi/Desktop/new_algorithm_Kelp/round2_free_living_ctg_R2.%s' % free_living_ext
free_living_ctg_UP = '/Users/songweizhi/Desktop/new_algorithm_Kelp/round2_free_living_ctg_UP.%s' % free_living_ext
min_M_len                               = 30
min_clp_len_round2                      = 30
max_mis_pct                             = 3
end_seq_len                             = 1000  # bp


round_2_MappingRecord_dict = {}
for each_line in open(sam_file):
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
free_living_ctg_R1_handle = open(free_living_ctg_R1, 'w')
free_living_ctg_R2_handle = open(free_living_ctg_R2, 'w')
free_living_ctg_UP_handle = open(free_living_ctg_UP, 'w')
for read_basename in round_2_MappingRecord_dict.copy():
    read_mr = round_2_MappingRecord_dict[read_basename]
    r1_ref_cigar_list = list(read_mr.r1_refs.values())
    r2_ref_cigar_list = list(read_mr.r2_refs.values())
    best_cigar, max_clp, max_clp_location = get_max_clp_and_index(r1_ref_cigar_list, r2_ref_cigar_list)

    if (read_mr.r1_refs == {}) and (read_mr.r2_refs != {}):

        if len(read_mr.r2_refs) == 1:

            r2_ref_cigar    = r2_ref_cigar_list[0]
            r2_ref_flag     = read_mr.r2_cigar_to_flag[r2_ref_cigar]
            r2_ref_cigar_rc = sam_flag_to_rc(r2_ref_flag)

            # consider the unmapped mate only
            if max_clp < min_clp_len_round2:
                read_mr.consider_round_2 = True
                read_mr.consider_r2_unmapped_mate = True

                # write out sequence
                if round2_fq is False:
                    free_living_ctg_UP_handle.write('>%s.1\n' % read_basename)
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_seq)
                else:
                    free_living_ctg_UP_handle.write('@%s.1\n' % read_basename)
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_seq)
                    free_living_ctg_UP_handle.write('+\n')
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_seq_qual)

            # consider both of unmapped mate and clipping part
            else:
                r2_ref_no_pos = list(read_mr.r2_refs.keys())[0].split('_pos_')[0]
                r2_ref_pos    = int(list(read_mr.r2_refs.keys())[0].split('_pos_')[1])

                if (r2_ref_no_pos[-1]) == (max_clp_location[-1]) == 'l':
                    read_mr.consider_round_2            = True
                    read_mr.consider_r2_unmapped_mate   = True
                    read_mr.consider_r2_clipping_part   = True

                    if r2_ref_cigar_rc is False:
                        read_mr.r2_clipping_seq          = read_mr.r2_seq[:max_clp]
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
                    read_mr.consider_round_2            = True
                    read_mr.consider_r2_unmapped_mate   = True
                    read_mr.consider_r2_clipping_part   = True

                    if r2_ref_cigar_rc is False:
                        read_mr.r2_clipping_seq          = read_mr.r2_seq[-max_clp:]
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

                # write out sequence
                if read_mr.consider_round_2 is True:

                    if round2_fq is False:
                        # write out R1 fa
                        free_living_ctg_R1_handle.write('>%s.1\n' % read_basename)
                        free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_seq)
                        # write out R2 fa
                        free_living_ctg_R2_handle.write('>%s.2\n' % read_basename)
                        free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_clipping_seq)
                    else:
                        # write out R1 fq
                        free_living_ctg_R1_handle.write('@%s.1\n' % read_basename)
                        free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_seq)
                        free_living_ctg_R1_handle.write('+\n')
                        free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_seq_qual)
                        # write out R2 fq
                        free_living_ctg_R2_handle.write('@%s.2\n' % read_basename)
                        free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_clipping_seq)
                        free_living_ctg_R2_handle.write('+\n')
                        free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_clipping_seq_qual)

                # print('r1_clipping_seq: %s' % read_mr.r1_clipping_seq)
                # print('r1_clipping_seq: %s' % read_mr.r1_clipping_seq_qual)
                # print('r2_clipping_seq: %s' % read_mr.r2_clipping_seq)
                # print('r2_clipping_seq: %s' % read_mr.r2_clipping_seq_qual)
                # print()

        else:  # r2 mapped to multiple refs, ignore
            round_2_MappingRecord_dict.pop(read_basename)

    elif (read_mr.r1_refs != {}) and (read_mr.r2_refs == {}):

        if len(read_mr.r1_refs) == 1:

            r1_ref_cigar    = r1_ref_cigar_list[0]
            r1_ref_flag     = read_mr.r1_cigar_to_flag[r1_ref_cigar]
            r1_ref_cigar_rc = sam_flag_to_rc(r1_ref_flag)

            # consider the unmapped mate only
            if max_clp < min_clp_len_round2:
                read_mr.consider_round_2 = True
                read_mr.consider_r1_unmapped_mate = True

                # write out sequence
                if round2_fq is False:
                    free_living_ctg_UP_handle.write('>%s.2\n' % read_basename)
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_seq)
                else:
                    free_living_ctg_UP_handle.write('@%s.2\n' % read_basename)
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_seq)
                    free_living_ctg_UP_handle.write('+\n')
                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_seq_qual)

            # consider both of unmapped mate and clipping part
            else:
                r1_ref_no_pos = list(read_mr.r1_refs.keys())[0].split('_pos_')[0]
                r1_ref_pos    = int(list(read_mr.r1_refs.keys())[0].split('_pos_')[1])

                if (r1_ref_no_pos[-1]) == (max_clp_location[-1]) == 'l':
                    read_mr.consider_round_2            = True
                    read_mr.consider_r1_unmapped_mate   = True
                    read_mr.consider_r1_clipping_part   = True

                    if r1_ref_cigar_rc is False:
                        read_mr.r1_clipping_seq          = read_mr.r1_seq[:max_clp]
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
                    read_mr.consider_round_2            = True
                    read_mr.consider_r1_unmapped_mate   = True
                    read_mr.consider_r1_clipping_part   = True

                    if r1_ref_cigar_rc is False:
                        read_mr.r1_clipping_seq          = read_mr.r1_seq[-max_clp:]
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

                # write out sequence
                if read_mr.consider_round_2 is True:

                    if round2_fq is False:
                        # write out R1 fa
                        free_living_ctg_R1_handle.write('>%s.1\n' % read_basename)
                        free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_clipping_seq)
                        # write out R2 fa
                        free_living_ctg_R2_handle.write('>%s.2\n' % read_basename)
                        free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_seq)
                    else:
                        # write out R1 fq
                        free_living_ctg_R1_handle.write('@%s.1\n' % read_basename)
                        free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_clipping_seq)
                        free_living_ctg_R1_handle.write('+\n')
                        free_living_ctg_R1_handle.write('%s\n' % read_mr.r1_clipping_seq_qual)
                        # write out R2 fq
                        free_living_ctg_R2_handle.write('@%s.2\n' % read_basename)
                        free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_seq)
                        free_living_ctg_R2_handle.write('+\n')
                        free_living_ctg_R2_handle.write('%s\n' % read_mr.r2_seq_qual)

        else:  # r1 mapped to multiple refs, ignore
            round_2_MappingRecord_dict.pop(read_basename)

    elif (read_mr.r1_refs != {}) and (read_mr.r2_refs != {}):
        if max_clp >= min_M_len:
            if (len(read_mr.r1_refs) == 1) and (len(read_mr.r2_refs) == 1):
                r1_ref_no_pos = list(read_mr.r1_refs.keys())[0].split('_pos_')[0]
                r2_ref_no_pos = list(read_mr.r2_refs.keys())[0].split('_pos_')[0]
                r1_ref_pos    = int(list(read_mr.r1_refs.keys())[0].split('_pos_')[1])
                r2_ref_pos    = int(list(read_mr.r2_refs.keys())[0].split('_pos_')[1])
                r1_ref_cigar  = list(read_mr.r1_refs.values())[0]
                r2_ref_cigar  = list(read_mr.r2_refs.values())[0]

                if r1_ref_no_pos == r2_ref_no_pos:
                    if (r1_ref_no_pos[-1]) == (r2_ref_no_pos[-1]) == (max_clp_location[-1]):

                        # print('%s r1_seq: %s'  % (read_basename, read_mr.r1_seq))
                        # print('%s r2_seq: %s'  % (read_basename, read_mr.r2_seq))
                        # print('%s r1_refs: %s' % (read_basename, read_mr.r1_refs))
                        # print('%s r2_refs: %s' % (read_basename, read_mr.r2_refs))
                        # print('%s r1_cigar_list: %s' % (read_basename, r1_ref_cigar_list))
                        # print('%s r1_cigar_to_flag: %s' % (read_basename, read_mr.r1_cigar_to_flag))
                        # print('%s r2_cigar_list: %s' % (read_basename, r2_ref_cigar_list))
                        # print('%s r2_cigar_to_flag: %s' % (read_basename, read_mr.r2_cigar_to_flag))
                        # print()
                        #print(r2_ref_flag)
                        #print(r2_ref_cigar_rc)

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

                        elif (max_clp_location == 'r1_r') and (r1_ref_pos >= (end_seq_len/2)):
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

                        elif (max_clp_location == 'r2_r') and (r2_ref_pos >= (end_seq_len/2)):
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

                        # write out sequence
                        if read_mr.consider_round_2 is True:

                            if read_mr.consider_r1_clipping_part is True:
                                if round2_fq is False:
                                    free_living_ctg_UP_handle.write('>%s.1\n' % read_basename)
                                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_clipping_seq)
                                else:
                                    free_living_ctg_UP_handle.write('@%s.1\n' % read_basename)
                                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_clipping_seq)
                                    free_living_ctg_UP_handle.write('+\n')
                                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r1_clipping_seq_qual)

                            if read_mr.consider_r2_clipping_part is True:
                                if round2_fq is False:
                                    free_living_ctg_UP_handle.write('>%s.2\n' % read_basename)
                                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_clipping_seq)
                                else:
                                    free_living_ctg_UP_handle.write('@%s.2\n' % read_basename)
                                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_clipping_seq)
                                    free_living_ctg_UP_handle.write('+\n')
                                    free_living_ctg_UP_handle.write('%s\n' % read_mr.r2_clipping_seq_qual)

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

free_living_ctg_R1_handle.close()
free_living_ctg_R2_handle.close()
free_living_ctg_UP_handle.close()


