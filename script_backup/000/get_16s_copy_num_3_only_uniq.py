
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


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


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


def check_cigar_quality(cigar_str, mismatch_cutoff, min_M_len, ref_pos, marker_len):

    # check the following:
    # 1. both end clip  2. mismatch     3. aligned length   4. clp in the middle

    cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        aln_len, aln_pct, clp_len, clp_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

        # check mismatch
        passed_mismatch_check = False
        if mismatch_pct <= mismatch_cutoff:
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


min_M_len_16s                       = 45
mismatch_cutoff                     = 2
sam_basename                        = 'reads_16s_to_16s_sam_all'
pwd_sam_file                        = '/Users/songweizhi/Desktop/666/000/%s.sam'                                     % sam_basename
pwd_sam_file_filtered               = '/Users/songweizhi/Desktop/666/000/%s_all_mis%s_minM%s_filtered.sam'           % (sam_basename, mismatch_cutoff, min_M_len_16s)


read_to_ref_dict = {}
marker_len_dict = {}
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
    else:
        cigar_str = each_line_split[5]
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        ref_pos = int(each_line_split[3])
        r1_16s_ref_qualified_cigar = check_cigar_quality(cigar_str, mismatch_cutoff, min_M_len_16s, ref_pos, marker_len_dict[ref_id])

        if r1_16s_ref_qualified_cigar is True:

            if read_id not in read_to_ref_dict:
                read_to_ref_dict[read_id] = [ref_id]
            else:
                read_to_ref_dict[read_id].append(ref_id)

read_to_ref_dict_uniq = {}
for each_read in read_to_ref_dict:
    matched_ref_list = read_to_ref_dict[each_read]
    if len(matched_ref_list) == 1:
        read_to_ref_dict_uniq[each_read] = matched_ref_list[0]


# filter sam file
pwd_sam_file_filtered_handle = open(pwd_sam_file_filtered, 'w')
for each_line in open(pwd_sam_file):
    if each_line.startswith('@'):
        pwd_sam_file_filtered_handle.write(each_line)
    else:
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        cigar_str = each_line_split[5]
        if read_id in read_to_ref_dict_uniq:
            if ref_id == read_to_ref_dict_uniq[read_id]:
                if ('S' not in cigar_str) and ('s' not in cigar_str):
                    pwd_sam_file_filtered_handle.write(each_line)
pwd_sam_file_filtered_handle.close()


print(len(read_to_ref_dict))
print(len(read_to_ref_dict_uniq))
