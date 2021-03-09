
def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if each_element.isalpha() is True:
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


sam_file_mini_assembly  = '/Users/songweizhi/Desktop/scaffolds.sam'

gap_seq_to_reads_dict = {}
for each_line in open(sam_file_mini_assembly):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        cigar = each_line_split[5]
        cigar_splitted = cigar_splitter(cigar)

        qualified_mapping = False
        if (len(cigar_splitted) == 1) and (cigar[-1] == 'M'):
            qualified_mapping = True
            print(each_line_split)
        # allow one mismatch
        elif (len(cigar_splitted) == 3) and (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M'):
            left_m_len = int(cigar_splitted[0][:-1])
            mismatch_len = int(cigar_splitted[1][:-1])
            right_m_len = int(cigar_splitted[2][:-1])
            read_len = left_m_len + right_m_len + mismatch_len
            left_m_pct = left_m_len*100/read_len
            right_m_pct = right_m_len*100/read_len
            if (left_m_len >= 30) and (right_m_len >= 30) and (left_m_pct >= 25) and (right_m_pct >= 25) and (mismatch_len <= 1):
                qualified_mapping = True

            else:
                pass
                #print(each_line_split)

        # elif len(cigar_splitted) == 4:



        if qualified_mapping is True:
            if ref_id not in gap_seq_to_reads_dict:
                gap_seq_to_reads_dict[ref_id] = [read_id]
            else:
                gap_seq_to_reads_dict[ref_id].append(read_id)




# print(gap_seq_to_reads_dict)
# print(len(gap_seq_to_reads_dict))
