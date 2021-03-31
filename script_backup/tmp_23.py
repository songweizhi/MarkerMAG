
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


pwd_samfile                         = '/Users/songweizhi/Desktop/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.sam'
clipping_reads_not_matched_part_seq = '/Users/songweizhi/Desktop/clipping_reads_not_matched_part.fa'
perfect_match_min_cigar_M_pct       = 70
global_max_mismatch_pct             = 3


# export clipping mapped reads and perfectly mapped reads
all_mapped_reads_set = set()
clipping_mapped_reads_list = set()
clipping_reads_mapped_part_dict = {}
perfectly_mapped_reads_dict = {}
clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
for each_read in open(pwd_samfile):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
            ref_pos = int(each_read_split[3])
            read_seq = each_read_split[9]
            cigar_splitted = cigar_splitter(cigar)
            cigar_match_pct, cigar_mismatch_pct = get_cigar_matched_and_mismatch_pct(cigar)
            read_id_with_ref_pos = '%s__x__%s__x__%s' % (read_id, ref_id, ref_pos)
            all_mapped_reads_set.add(read_id)

            # treat_as_full_match and store into dict
            if (cigar_match_pct >= perfect_match_min_cigar_M_pct) and (
                    cigar_mismatch_pct <= global_max_mismatch_pct):

                if read_id_base not in perfectly_mapped_reads_dict:
                    perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
                else:
                    if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                        perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                    else:
                        perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

            # treat as clipping matched
            if (cigar_match_pct < perfect_match_min_cigar_M_pct) and (
                    cigar_mismatch_pct <= global_max_mismatch_pct):

                # only one end of cigar is S
                if ((cigar_splitted[0][-1] == 'S') and (cigar_splitted[-1][-1] != 'S')) or (
                        (cigar_splitted[0][-1] != 'S') and (cigar_splitted[-1][-1] == 'S')):

                    # if clipped at left
                    if cigar_splitted[0][-1] == 'S':
                        read_seq_clipped = read_seq[:int(cigar_splitted[0][:-1])]

                        # write out the clipped part
                        clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id_with_ref_pos)
                        clipping_reads_not_matched_part_seq_handle.write(read_seq_clipped + '\n')

                        # store the matching info of aligned part
                        if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                            clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                        else:
                            clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                    # if clipped at right
                    if cigar_splitted[-1][-1] == 'S':
                        read_seq_clipped = read_seq[-int(cigar_splitted[-1][:-1]):]

                        # write out the clipped part
                        clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id_with_ref_pos)
                        clipping_reads_not_matched_part_seq_handle.write(read_seq_clipped + '\n')

                        # store the matching info of aligned part
                        if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                            clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                        else:
                            clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                    clipping_mapped_reads_list.add(read_id)

clipping_reads_not_matched_part_seq_handle.close()


perfectly_mapped_read_singleton_dict = {}
for perfectly_mapped_read in perfectly_mapped_reads_dict:
    current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
    if len(current_value) == 1:
        perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value

print(len(perfectly_mapped_read_singleton_dict))
