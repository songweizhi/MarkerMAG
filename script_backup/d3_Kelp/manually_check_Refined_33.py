
'''

Refined_33	11.06	no	83.3	no	d__Bacteria;p__Myxococcota;c__Polyangia;o__Haliangiales;f__Haliangiaceae;g__;s__

BH_ER_050417_subsample_100_98	d__Bacteria;p__Myxococcota;c__Polyangia;o__Haliangiales;f__Haliangiaceae;g__Haliangium;s__Haliangium_ochraceum
BH_ER_050417_subsample_75_86	d__Bacteria;p__Myxococcota;c__Polyangia;o__Haliangiales;f__Haliangiaceae;g__Haliangium;s__Haliangium_ochraceum

'''

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
perfect_match_min_cigar_M_pct = 50
perfect_match_max_cigar_S_pct = 50
mismatch_ratio_max_value = 1


# export clipping mapped reads and perfectly mapped reads
all_mapped_reads_set = set()
clipping_mapped_reads_list = set()
clipping_reads_mapped_part_dict = {}
perfectly_mapped_reads_dict = {}
for each_read in open(pwd_samfile):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        cigar = each_read_split[5]
        if cigar != '*':

            read_id = each_read_split[0]
            ref_id = each_read_split[2]

            if ref_id in ['BH_ER_050417_subsample_100_98', 'BH_ER_050417_subsample_75_86']:

                # if read_id in ['Kelp_31369131.1', 'Kelp_31369131.2']:


                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
                ref_pos = int(each_read_split[3])
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
                    #print(each_read_split)

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
                    # print(cigar_M_pct)
                    # print(cigar_S_pct)

                    # for clipping reads with unmapped part >= min_cigar_S
                    if (cigar_M_pct >= 30) and (cigar_S_pct >= 30):
                        read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                        read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                        if cigar_splitted[0][-1] == 'M':

                            # store the match info of mapped part
                            if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                                clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                            else:
                                clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                        if cigar_splitted[1][-1] == 'M':

                            # store the match info of mapped part
                            if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                                clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                            else:
                                clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                        clipping_mapped_reads_list.add(read_id)

                    # for clipping reads with unmapped part < min_cigar_S
                    # treat these reads as perfectly mapped reads
                    if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
                        if read_id_base not in perfectly_mapped_reads_dict:
                            perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
                        else:
                            if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                                perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                            else:
                                perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

                elif ('S' in cigar) and (len(cigar_splitted) > 2):

                    if ((cigar_splitted[0][-1] == 'M') and (cigar_splitted[-1][-1] == 'S')) or (
                            (cigar_splitted[0][-1] == 'S') and (cigar_splitted[-1][-1] == 'M')):

                        if len(cigar_splitted) == 4:

                            # e.g. 30M-1D-50M-40S
                            if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (
                                    cigar_splitted[3][-1] == 'S'):
                                mismatch_ratio = int(cigar_splitted[1][:-1]) * 100 / len(read_seq)
                                if (mismatch_ratio <= mismatch_ratio_max_value) and (int(cigar_splitted[2][:-1]) >= 30):
                                    cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1])) * 100 / len(
                                        read_seq)
                                    cigar_S_pct = (int(cigar_splitted[3][:-1])) * 100 / len(read_seq)
                                    if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                                            cigar_S_pct < perfect_match_max_cigar_S_pct):
                                        treat_as_full_match = True

                            # e.g. 40S-30M-1D-50M
                            if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (
                                    cigar_splitted[3][-1] == 'M'):
                                mismatch_ratio = int(cigar_splitted[2][:-1]) * 100 / len(read_seq)
                                if (mismatch_ratio <= mismatch_ratio_max_value) and (int(cigar_splitted[1][:-1]) >= 30):
                                    cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1])) * 100 / len(
                                        read_seq)
                                    cigar_S_pct = (int(cigar_splitted[0][:-1])) * 100 / len(read_seq)
                                    if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                                            cigar_S_pct < perfect_match_max_cigar_S_pct):
                                        treat_as_full_match = True

                        # elif len(cigar_splitted) == 6:
                        #     if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (cigar_splitted[4][-1] == 'M') and (cigar_splitted[5][-1] == 'S'):
                        #         mismatch_ratio = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1])) * 100 / len(read_seq)
                        #         if mismatch_ratio <= mismatch_ratio_max_value:
                        #             cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1]) + int(cigar_splitted[4][:-1])) * 100 / len(read_seq)
                        #             cigar_S_pct = (int(cigar_splitted[5][:-1])) * 100 / len(read_seq)
                        #             if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                        #                     cigar_S_pct < perfect_match_max_cigar_S_pct):
                        #                 treat_as_full_match = True
                        #
                        #     if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (cigar_splitted[3][-1] == 'M') and (cigar_splitted[5][-1] == 'M'):
                        #         mismatch_ratio = (int(cigar_splitted[2][:-1]) + int(cigar_splitted[4][:-1])) * 100 / len(read_seq)
                        #         if mismatch_ratio <= mismatch_ratio_max_value:
                        #             cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1]) + int(cigar_splitted[5][:-1])) * 100 / len(read_seq)
                        #             cigar_S_pct = (int(cigar_splitted[0][:-1])) * 100 / len(read_seq)
                        #             if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
                        #                 treat_as_full_match = True

                if treat_as_full_match is True:
                    if read_id_base not in perfectly_mapped_reads_dict:
                        perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
                    else:
                        if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                            perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                        else:
                            perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)


########################################## extract reads with multiprocessing ##########################################

perfectly_mapped_read_singleton_dict = {}
for perfectly_mapped_read in perfectly_mapped_reads_dict:
    current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
    if len(current_value) == 1:
        perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value

for each in perfectly_mapped_read_singleton_dict:
    print('%s\t%s' % (each, perfectly_mapped_read_singleton_dict[each]))

print(len(perfectly_mapped_read_singleton_dict))

# Kelp_31369131.2
print(perfectly_mapped_reads_dict.get('Kelp_31369131'))

'''

Kelp_43912804.2
'''