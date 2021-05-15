import os


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


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


def remove_both_ends_clp(sam_in, sam_out):
    sam_out_handle = open(sam_out, 'w')
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)
        else:
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    sam_out_handle.write(each_read)
    sam_out_handle.close()


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


def remove_clp_in_middle(sam_in, sam_out):

    sam_out_handle = open(sam_out, 'w')

    marker_len_dict = {}
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)

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

            # check if clp in the middle
            if ('S' not in cigar) and ('s' not in cigar):
                sam_out_handle.write(each_read)
            else:
                ref_id = each_read_split[2]
                ref_pos = int(each_read_split[3])
                cigar_splitted = cigar_splitter(cigar)
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)

                clip_in_middle = True
                if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                    clip_in_middle = False
                if (cigar_splitted[-1][-1] in ['S', 's']):
                    if (ref_pos + r1_aligned_len - 1) == marker_len_dict[ref_id]:
                        clip_in_middle = False

                if clip_in_middle is False:
                    sam_out_handle.write(each_read)

    sam_out_handle.close()


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


def remove_high_mismatch(sam_in, mismatch_cutoff, sam_out):
    sam_out_handle = open(sam_out, 'w')
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)
        else:
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                cigar_splitted = cigar_splitter(cigar)
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)
                if r1_mismatch_pct <= mismatch_cutoff:
                    sam_out_handle.write(each_read)
    sam_out_handle.close()



input_sam       = '5_261.sam'
#input_sam       = 'linked_to_hc_99.sam'
mismatch_cutoff = 2

sam_file_path, sam_file_basename, sam_file_extension = sep_path_basename_ext(input_sam)
sam_file_one_end_clp                                        = '%s/%s_one_end_clp.sam'                                               % (sam_file_path, sam_file_basename)
sam_file_one_end_clp_no_middle                              = '%s/%s_one_end_clp_no_middle.sam'                                     % (sam_file_path, sam_file_basename)
sam_file_one_end_clp_reformatted                            = '%s/%s_one_end_clp_no_middle_reformatted.sam'                         % (sam_file_path, sam_file_basename)
sam_file_one_end_clp_reformatted_log                        = '%s/%s_one_end_clp_no_middle_reformatted.log'                         % (sam_file_path, sam_file_basename)
sam_file_one_end_clp_reformatted_best_match                 = '%s/%s_one_end_clp_no_middle_reformatted_best_match.sam'              % (sam_file_path, sam_file_basename)
sam_file_one_end_clp_reformatted_best_match_low_mismatch    = '%s/%s_one_end_clp_no_middle_reformatted_best_match_low_mismatch.sam' % (sam_file_path, sam_file_basename)
sam_file_sorted                                             = '%s/%s_sorted.sam'                                                    % (sam_file_path, sam_file_basename)
coverage_file                                               = '%s/%s_cov.txt'                                                       % (sam_file_path, sam_file_basename)


# filter mapping
remove_both_ends_clp(input_sam, sam_file_one_end_clp)
remove_clp_in_middle(sam_file_one_end_clp, sam_file_one_end_clp_no_middle)
bbmap_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_one_end_clp_no_middle, sam_file_one_end_clp_reformatted, sam_file_one_end_clp_reformatted_log)
os.system(bbmap_reformat_cmd)

keep_best_matches_in_sam_keep_short_M(sam_file_one_end_clp_reformatted, 35, sam_file_one_end_clp_reformatted_best_match)

remove_high_mismatch(sam_file_one_end_clp_reformatted_best_match, mismatch_cutoff, sam_file_one_end_clp_reformatted_best_match_low_mismatch)

# sort mapping
cmd_samtools_sort = 'samtools sort %s -o %s' % (sam_file_one_end_clp_reformatted_best_match_low_mismatch, sam_file_sorted)
os.system(cmd_samtools_sort)

# get mean depth
cmd_samtools_coverage = 'samtools --ff 4 coverage %s -o %s' % (sam_file_sorted, coverage_file)
os.system(cmd_samtools_coverage)



'''
#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
cami_hc_SILVA138_id99_50_subsample_5_261	1	1560	4587	1560	100	286.179	40	255
cami_hc_SILVA138_id99_50_subsample_5_261	1	1560	4587	1560	100	286.179	40	255

cami_hc_SILVA138_id99_50_subsample_5_261	1	1560	4117	1560	100	255.337	40	0.958

cami_hc_SILVA138_id99_50_subsample_100_894	1	1573	4130	1573	100	245.599	40	0.926
cami_hc_SILVA138_id99_50_subsample_100_793	1	1573	3743	1573	100	228.688	40	0.944
cami_hc_SILVA138_id99_50_subsample_50_608	1	1573	4061	1573	100	248.905	40	0.944
cami_hc_SILVA138_id99_50_subsample_100_891	1	1573	3115	1572	99.9364	189.449	40	1.01
cami_hc_SILVA138_id99_50_subsample_100_1029	1	1569	4388	1569	100	270.631	40	0.948
cami_hc_SILVA138_id99_50_subsample_100_1027	1	1569	4168	1569	100	256.609	40	0.946
cami_hc_SILVA138_id99_50_subsample_5_263	1	1560	4336	1560	100	269.605	40	0.961
cami_hc_SILVA138_id99_50_subsample_50_660	1	1478	2903	1478	100	197.816	40	1.02




cami_hc_SILVA138_id99_50_subsample_5_261	1	1560	609	1553	99.5513	37.4974	40	0.969

cami_hc_SILVA138_id99_50_subsample_100_894	1	1573	532	1573	100	30.5683	40	0.914
cami_hc_SILVA138_id99_50_subsample_100_793	1	1573	494	1571	99.8729	29.417	40	0.947
cami_hc_SILVA138_id99_50_subsample_50_608	1	1573	592	1573	100	36.1087	40	0.938
cami_hc_SILVA138_id99_50_subsample_100_891	1	1573	408	1569	99.7457	25.9282	40	1.53
cami_hc_SILVA138_id99_50_subsample_100_1029	1	1569	665	1569	100	41.144	40	0.959
cami_hc_SILVA138_id99_50_subsample_100_1027	1	1569	582	1569	100	36.334	40	0.966
cami_hc_SILVA138_id99_50_subsample_5_263	1	1560	620	1560	100	39.7391	40	1
cami_hc_SILVA138_id99_50_subsample_50_660	1	1478	349	1478	100	24.6698	40	1.11


'''