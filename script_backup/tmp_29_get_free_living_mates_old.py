import os
from Bio import SeqIO
from datetime import datetime


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


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
        if each_part_cate in {'I', 'X', 'D'}:
            mismatched_seq_len += each_part_len

    matched_pct    = float("{0:.2f}".format(matched_seq_len*100/total_len_with_s))
    mismatch_pct = float("{0:.2f}".format(mismatched_seq_len*100/total_len_without_s))

    return matched_pct, mismatch_pct


def get_free_living_mate(ref_in, reads_r1, reads_r2, end_seq_len, max_gap_to_end, num_threads, pwd_bbmap_exe, bbmap_memory, perfect_match_min_cigar_M_pct, global_max_mismatch_pct, pwd_log_file, keep_quiet):

    ref_in_path, ref_in_basename, ref_in_ext = sep_path_basename_ext(ref_in)

    ref_subset      = '%s/%s_ends_%sbp%s'                % (ref_in_path, ref_in_basename, end_seq_len, ref_in_ext)
    sam_file        = '%s/%s_ends_%sbp.sam'              % (ref_in_path, ref_in_basename, end_seq_len)
    bbmap_stderr    = '%s/%s_ends_%sbp_bbmap_stderr.txt' % (ref_in_path, ref_in_basename, end_seq_len)

    # get ref seqs subset
    report_and_log(('Round 2: map reads to %s, get ends sequences' % ref_in_basename), pwd_log_file, keep_quiet)
    ref_subset_handle = open(ref_subset, 'w')
    for ref_seq in SeqIO.parse(ref_in, 'fasta'):

        ref_seq_id = ref_seq.id
        ref_seq_len = len(ref_seq.seq)

        if ref_seq_len < end_seq_len * 2:
            ref_subset_handle.write('>%s\n' % ref_seq_id)
            ref_subset_handle.write('%s\n' % ref_seq.seq)
        else:
            ref_seq_left_end_id = '%s_l' % ref_seq_id
            ref_seq_right_end_id = '%s_r' % ref_seq_id
            ref_seq_left_end = ref_seq.seq[:end_seq_len]
            ref_seq_right_end = ref_seq.seq[-end_seq_len:]

            # write out left end
            ref_subset_handle.write('>%s\n' % ref_seq_left_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_left_end)

            # write out right end
            ref_subset_handle.write('>%s\n' % ref_seq_right_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_right_end)
    ref_subset_handle.close()

    # mapping with bbmap
    report_and_log(('Round 2: map reads to %s, mapping' % ref_in_basename), pwd_log_file, keep_quiet)
    bbmap_parameter_round2 = 'local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=%s -Xmx%sg' % (num_threads, bbmap_memory)
    bbmap_cmd_round2 = '%s ref=%s in=%s in2=%s outm=%s %s 2> %s' % (pwd_bbmap_exe, ref_subset, reads_r1, reads_r2, sam_file, bbmap_parameter_round2, bbmap_stderr)
    os.system(bbmap_cmd_round2)

    # parse sam file
    report_and_log(('Round 2: map reads to %s, filter sam' % ref_in_basename), pwd_log_file, keep_quiet)
    all_mapped_reads = set()
    qualified_reads_dict = {}
    qualified_reads_to_ref_dict = {}
    ref_in_sam_len_dict = {}
    for each_line in open(sam_file):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('@'):
            ref_id = ''
            ref_len = 0
            for each_element in each_line_split:
                if each_element.startswith('SN:'):
                    ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_len = int(each_element[3:])
            ref_in_sam_len_dict[ref_id] = ref_len
        else:
            cigar = each_line_split[5]
            if cigar != '*':
                read_id = each_line_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_line_split[2]
                ref_len = ref_in_sam_len_dict[ref_id]
                ref_pos = int(each_line_split[3])
                cigar_splitted = cigar_splitter(cigar)
                cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)
                all_mapped_reads.add(read_id)

                # get aligned length
                aligned_len = 0
                for each_section in cigar_splitted:
                    each_section_len = int(each_section[:-1])
                    each_section_cate = each_section[-1]
                    if each_section_cate in {'D', '=', 'X'}:
                        aligned_len += each_section_len

                # filter mapped reads
                qualified_mapping = False
                if (cigar_match_pct >= perfect_match_min_cigar_M_pct) and (cigar_mismatch_pct <= global_max_mismatch_pct):

                    # check left end for contig's left end sequence
                    if ref_id[-2:] == '_l':
                        if ref_pos <= max_gap_to_end:
                            qualified_mapping = True

                    # check right end for contig's right end sequence
                    elif ref_id[-2:] == '_r':
                        if cigar_splitted[0][-1] != 'S':
                            ref_pos_end = ref_pos + aligned_len
                            if (ref_len - ref_pos_end) <= max_gap_to_end:
                                qualified_mapping = True

                    # check both ends for contigs without subset
                    else:
                        if (ref_pos <= max_gap_to_end) or ((ref_len - ref_pos - aligned_len) <= max_gap_to_end):
                            qualified_mapping = True

                # proceed with qualified mappings
                if qualified_mapping is True:
                    qualified_reads_to_ref_dict[read_id] = ref_id

                    if read_id_base not in qualified_reads_dict:
                        qualified_reads_dict[read_id_base] = [read_strand]
                    else:
                        qualified_reads_dict[read_id_base].append(read_strand)

    reads_to_extract_to_ref_dict = {}
    for qualified_read in qualified_reads_dict:
        read_strand = qualified_reads_dict[qualified_read]
        if len(read_strand) == 1:
            mapped_mate = ''
            mate_to_extract = ''
            if read_strand == ['1']:
                mapped_mate = '%s.1' % (qualified_read)
                mate_to_extract = '%s.2' % (qualified_read)
            if read_strand == ['2']:
                mapped_mate = '%s.2' % (qualified_read)
                mate_to_extract = '%s.1' % (qualified_read)

            if mate_to_extract not in all_mapped_reads:
                reads_to_extract_to_ref_dict[mate_to_extract] = qualified_reads_to_ref_dict[mapped_mate]

    return reads_to_extract_to_ref_dict

