
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
        if (len(cigar_splitted) == 1) and (cigar[-1] == 'M'):
            if ref_id not in gap_seq_to_reads_dict:
                gap_seq_to_reads_dict[ref_id] = [read_id]
            else:
                gap_seq_to_reads_dict[ref_id].append(read_id)

        elif (len(cigar_splitted) == 3):


            if () and ():
                print(cigar_splitted)



        # allow one mismatch
        # elif

print(len(gap_seq_to_reads_dict))
