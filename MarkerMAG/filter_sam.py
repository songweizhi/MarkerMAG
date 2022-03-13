import argparse


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


def remove_high_mismatch(sam_in, aln_len_cutoff, mismatch_cutoff, sam_out):

    sam_out_handle = open(sam_out, 'w')
    ref_len_dict = {}
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
            ref_len_dict[marker_id] = marker_len

        else:
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                cigar_splitted = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitted)
                if both_ends_clp is False:
                    r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)
                    if r1_mismatch_pct <= mismatch_cutoff:

                        if r1_aligned_len >= aln_len_cutoff:

                            # check if clp in middle
                            if ('S' not in cigar) and ('s' not in cigar):
                                sam_out_handle.write(each_read)
                            else:
                                clip_in_middle = True
                                if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                                    clip_in_middle = False
                                if (cigar_splitted[-1][-1] in ['S', 's']):
                                    if (ref_pos + r1_aligned_len - 1) == ref_len_dict[ref_id]:
                                        clip_in_middle = False

                                if clip_in_middle is False:
                                    sam_out_handle.write(each_read)
    sam_out_handle.close()


# sam_file                       = 'FP_report_best.sam'
# sam_file_filtered              = 'FP_report_best_mis0.sam'
# aln_len_cutoff                 = 70
# mismatch_cutoff                = 0
#

parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help='input file')
parser.add_argument('-out', required=True, help='output file')
parser.add_argument('-mm', required=True, type=int, help='mismatch cutoff')
parser.add_argument('-aln', required=True, type=int, help='alignment length cutoff')
args = vars(parser.parse_args())

sam_file                       = args['in']
sam_file_filtered              = args['out']
mismatch_cutoff                = args['mm']
aln_len_cutoff                 = args['aln']

remove_high_mismatch(sam_file, aln_len_cutoff, mismatch_cutoff, sam_file_filtered)
