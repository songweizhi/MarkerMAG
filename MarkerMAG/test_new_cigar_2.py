
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


def get_cigar_matched_and_mismatch_pct(cigar):

    cigar_splitted = cigar_splitter(cigar)

    total_len_with_s = 0
    total_len_without_s = 0
    matched_seq_len = 0
    mismatched_seq_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]

        # get total len with/without s
        if each_part_cate in {'M', 'I', '=', 'X'}:
            total_len_with_s += each_part_len
            total_len_without_s += each_part_len

        # get total len with s
        if each_part_cate == 'S':
            total_len_with_s += each_part_len

        # get matched part len
        if each_part_cate == '=':
            matched_seq_len += each_part_len

        # get mismatched part len
        if each_part_cate in {'I', 'X'}:
            mismatched_seq_len += each_part_len

    matched_pct    = float("{0:.2f}".format(matched_seq_len*100/total_len_with_s))
    mismatch_pct = float("{0:.2f}".format(mismatched_seq_len*100/total_len_without_s))

    return matched_pct, mismatch_pct


def paired_blast_results_to_dict_by_mapping(unmapped_paired_reads_mapping_results):

    cigar_M_pct_min_value = 50
    mismatch_pct_max_value = 1

    unmapped_paired_reads_to_ctg_dict_by_mapping = {}

    ref_len_dict = {}
    for unmapped_read in open(unmapped_paired_reads_mapping_results):
        unmapped_read_split = unmapped_read.strip().split('\t')

        # get ref len dict
        if unmapped_read.startswith('@'):
            ref_id = ''
            ref_len = 0
            for each_element in unmapped_read_split:
                if each_element.startswith('SN:'):
                    ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_len = int(each_element[3:])
            ref_len_dict[ref_id] = ref_len

        else:
            cigar = unmapped_read_split[5]
            if cigar != '*':
                unmapped_read_split = unmapped_read.strip().split('\t')
                read_id = unmapped_read_split[0]
                ref_id = unmapped_read_split[2]
                ref_len = ref_len_dict[ref_id]
                ref_id_with_prefix = 'GenomicSeq__%s' % unmapped_read_split[2]
                ref_pos = int(unmapped_read_split[3])
                read_seq = unmapped_read_split[9]
                read_len = len(read_seq)
                cigar_splitted = cigar_splitter(cigar)
                qualified_unmapped_read = False

                # e.g. 189=
                if ('=' in cigar) and (len(cigar_splitted) == 1):
                    qualified_unmapped_read = True

                elif len(cigar_splitted) == 2:

                    # e.g. ['44S', '176=']
                    if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == '='):
                        cigar_M_pct = int(cigar_splitted[1][:-1]) * 100 / read_len
                        if (ref_pos == 1) and (cigar_M_pct >= cigar_M_pct_min_value):
                            qualified_unmapped_read = True

                    # e.g. ['154=', '66S']
                    if (cigar_splitted[0][-1] == '=') and (cigar_splitted[1][-1] == 'S'):
                        cigar_M_pct = int(cigar_splitted[0][:-1]) * 100 / read_len
                        matched_to_bp = ref_pos + int(cigar_splitted[0][:-1]) - 1
                        if (matched_to_bp == ref_len) and (cigar_M_pct >= cigar_M_pct_min_value):
                            qualified_unmapped_read = True


                elif len(cigar_splitted) == 3:

                    # e.g. ['134=', '1X', '22='], ['215=', '1X', '4=']
                    if (cigar_splitted[0][-1] == '=') and (cigar_splitted[2][-1] == '='):
                        mismatch_pct = int(cigar_splitted[1][:-1])*100/read_len
                        if mismatch_pct <= mismatch_pct_max_value:
                            qualified_unmapped_read = True

                elif len(cigar_splitted) == 4:

                    # e.g. ['118=', '1X', '99=', '1S'], ['8=', '1D', '127=', '84S']
                    if (cigar_splitted[0][-1] == '=') and (cigar_splitted[2][-1] == '=') and (cigar_splitted[3][-1] == 'S'):
                        cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1]))*100/read_len
                        if ((ref_pos + int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1])) == ref_len) and (int(cigar_splitted[1][:-1]) <= 2):
                            if cigar_M_pct >= cigar_M_pct_min_value:
                                qualified_unmapped_read = True

                    # e.g. ['51S', '158=', '1X', '10='] and ref_pos = 1
                    elif (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == '=') and (cigar_splitted[3][-1] == '='):
                        cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1]))*100/read_len
                        if (ref_pos == 1) and (int(cigar_splitted[2][:-1]) <= 1) and (cigar_M_pct >= cigar_M_pct_min_value):
                            qualified_unmapped_read = True

                # add to dict
                if qualified_unmapped_read is True:
                    if read_id not in unmapped_paired_reads_to_ctg_dict_by_mapping:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id] = [ref_id_with_prefix]
                    else:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id].append(ref_id_with_prefix)

    return unmapped_paired_reads_to_ctg_dict_by_mapping


def paired_blast_results_to_dict_by_mapping_new(unmapped_paired_reads_mapping_results):

    cigar_M_pct_min_value = 50
    mismatch_pct_max_value = 1

    unmapped_paired_reads_to_ctg_dict_by_mapping = {}

    for unmapped_read in open(unmapped_paired_reads_mapping_results):

        # get ref len dict
        if not unmapped_read.startswith('@'):
            unmapped_read_split = unmapped_read.strip().split('\t')
            cigar = unmapped_read_split[5]
            if cigar != '*':
                read_id = unmapped_read_split[0]
                ref_id_with_prefix = 'GenomicSeq__%s' % unmapped_read_split[2]
                cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)

                if (cigar_match_pct >= cigar_M_pct_min_value) and (cigar_mismatch_pct <= mismatch_pct_max_value):

                    if read_id not in unmapped_paired_reads_to_ctg_dict_by_mapping:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id] = [ref_id_with_prefix]
                    else:
                        unmapped_paired_reads_to_ctg_dict_by_mapping[read_id].append(ref_id_with_prefix)

    return unmapped_paired_reads_to_ctg_dict_by_mapping



pwd_samfile_to_mag = '/Users/songweizhi/Desktop/BH_ER_050417_refined_bins_combined.sam'
unmapped_paired_reads_to_ctg_dict = paired_blast_results_to_dict_by_mapping(pwd_samfile_to_mag)
print(len(unmapped_paired_reads_to_ctg_dict))

unmapped_paired_reads_to_ctg_dict = paired_blast_results_to_dict_by_mapping_new(pwd_samfile_to_mag)
print(len(unmapped_paired_reads_to_ctg_dict))



