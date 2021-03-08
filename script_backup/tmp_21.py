

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


pwd_samfile = '/Users/songweizhi/Desktop/BH_ER_050417_assembled_16S_uclust_0.999.sam'

all_mapped_reads_set = set()
for each_read in open(pwd_samfile):

    if not each_read.startswith('@'):
        print(each_read.strip())
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        ref_id = each_read_split[2]
        ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
        ref_pos = int(each_read_split[3])
        cigar = each_read_split[5]
        read_seq = each_read_split[9]
        cigar_splitted = cigar_splitter(cigar)
        read_id_with_ref_pos = '%s__x__%s__x__%s' % (read_id, ref_id, ref_pos)
        all_mapped_reads_set.add(read_id)
        treat_as_full_match = False

        # for perfectly mapped reads, e.g. 150M
        if ('M' in cigar) and (len(cigar_splitted) == 1):
            treat_as_full_match = True

        # for perfectly mapped reads, e.g. 120M-1D-80M
        if len(cigar_splitted) == 3:
            if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M'):
                mismatch_ratio = int(cigar_splitted[1][:-1]) * 100 / len(read_seq)
                M_pct_l = int(cigar_splitted[0][:-1]) * 100 / len(read_seq)
                M_pct_r = int(cigar_splitted[2][:-1]) * 100 / len(read_seq)
                if (M_pct_l >= 20) and (M_pct_r >= 20) and (mismatch_ratio <= 0.5):
                    treat_as_full_match = True

        if ('S' in cigar) and (len(cigar_splitted) == 2):
            cigar_M_len = 0
            cigar_S_len = 0
            if cigar_splitted[0][-1] == 'M':
                cigar_M_len = int(cigar_splitted[0][:-1])
                cigar_S_len = int(cigar_splitted[1][:-1])
            if cigar_splitted[1][-1] == 'M':
                cigar_M_len = int(cigar_splitted[1][:-1])
                cigar_S_len = int(cigar_splitted[0][:-1])

            cigar_M_pct = cigar_M_len * 100 / (cigar_M_len + cigar_S_len)
            cigar_S_pct = cigar_S_len * 100 / (cigar_M_len + cigar_S_len)

