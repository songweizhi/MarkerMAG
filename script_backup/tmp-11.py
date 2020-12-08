
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


unmapped_paired_reads_blastn    = '/Users/songweizhi/Desktop/unmapped_paired_reads_blast.txt'
reads_iden_cutoff               = 100  # s1_ri
reads_cov_cutoff                = 90   # s1_rc


pwd_samfile = '/Users/songweizhi/Desktop/BH_ER_050417_all_depth_assemblies.dereplicated.sam'
perfect_match_min_cigar_M_pct       = 90
perfect_match_max_cigar_S_pct       = 10
clipping_reads_not_matched_part_seq = '/Users/songweizhi/Desktop/test.fa'
error_rate_cutoff = 0.005


# export clipping mapped reads and perfectly mapped reads
perfect_num = 0
all_mapped_reads_set = set()
clipping_mapped_reads_list = set()
clipping_reads_mapped_part_dict = {}
perfectly_mapped_reads_dict = {}
clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
for each_read in open(pwd_samfile):

    qualified_perfect_match = False

    if not each_read.startswith('@'):
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

        # for perfectly mapped reads
        if ('M' in cigar) and (len(cigar_splitted) == 1):
            if read_id_base not in perfectly_mapped_reads_dict:
                perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
            else:
                if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                    perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                else:
                    perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

            perfect_num += 1

        # for clipping mapped reads
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

            # for clipping reads with unmapped part >= min_cigar_S
            if (cigar_M_pct >= 0.3) and (cigar_S_pct >= 0.3):
                read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                if cigar_splitted[0][-1] == 'M':

                    # write out the sequence of unmapped part
                    clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id_with_ref_pos)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_right + '\n')

                    # store the match info of mapped part
                    if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                if cigar_splitted[1][-1] == 'M':

                    # write out the sequence of unmapped part
                    clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id_with_ref_pos)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_left + '\n')

                    # store the match info of mapped part
                    if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                clipping_mapped_reads_list.add(read_id)

            # for clipping reads with cigar_S_pct < perfect_match_max_cigar_S_pct
            # treat these reads as perfectly mapped reads
            if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
                qualified_perfect_match = True


        # if ('S' in cigar) and (len(cigar_splitted) == 3):
        #     if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (cigar_splitted[2][-1] == 'S'):
        #         cigar_M_len = int(cigar_splitted[1][:-1])
        #         cigar_S_len_l = int(cigar_splitted[0][:-1])
        #         cigar_S_len_r = int(cigar_splitted[2][:-1])
        #         cigar_M_pct = cigar_M_len * 100 / (cigar_M_len + cigar_S_len_l + cigar_S_len_r)
        #         cigar_S_pct = (cigar_S_len_l + cigar_S_len_r) * 100 / (cigar_M_len + cigar_S_len_l + cigar_S_len_r)
        #         if (cigar_S_len_l <= 1) or (cigar_S_len_r <= 1):
        #             if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
        #                 qualified_perfect_match = True


        if ('S' in cigar) and (len(cigar_splitted) == 4):

            if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (cigar_splitted[3][-1] == 'M'):
                cigar_m1_len = int(cigar_splitted[1][:-1])
                cigar_m2_len = int(cigar_splitted[3][:-1])
                mismatch_len = int(cigar_splitted[2][:-1])
                total_M = cigar_m1_len + cigar_m2_len
                error_rate = mismatch_len / (total_M)
                cigar_S_len = int(cigar_splitted[0][:-1])
                cigar_M_pct = total_M * 100 / (total_M + cigar_S_len)
                cigar_S_pct = cigar_S_len * 100 / (total_M + cigar_S_len)
                if (cigar_m1_len >= 30) and (cigar_m2_len >= 30) and (error_rate <= error_rate_cutoff):
                    if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
                        qualified_perfect_match = True

            if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (cigar_splitted[3][-1] == 'S'):
                cigar_m1_len = int(cigar_splitted[0][:-1])
                cigar_m2_len = int(cigar_splitted[2][:-1])
                mismatch_len = int(cigar_splitted[1][:-1])
                total_M = cigar_m1_len + cigar_m2_len
                error_rate = mismatch_len / (total_M)
                cigar_S_len = int(cigar_splitted[3][:-1])
                cigar_M_pct = total_M * 100 / (total_M + cigar_S_len)
                cigar_S_pct = cigar_S_len * 100 / (total_M + cigar_S_len)
                if (cigar_m1_len >= 30) and (cigar_m2_len >= 30) and (error_rate <= error_rate_cutoff):
                    if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
                        qualified_perfect_match = True

        # add to dict
        if qualified_perfect_match is True:
            if read_id_base not in perfectly_mapped_reads_dict:
                perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
            else:
                if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                    perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                else:
                    perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

clipping_reads_not_matched_part_seq_handle.close()

bothed = 0
perfectly_mapped_read_singleton_dict = {}
for perfectly_mapped_read in perfectly_mapped_reads_dict:
    current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
    if len(current_value) == 1:
        perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value
    if len(current_value) == 2:
        bothed += 1


# get the id of paired reads to extract
solely_perfectly_mapped_reads_r1 = set()
solely_perfectly_mapped_reads_r2 = set()
for perfectly_mapped_read in perfectly_mapped_read_singleton_dict:
    current_value = perfectly_mapped_read_singleton_dict[perfectly_mapped_read]
    strand = list(current_value.keys())[0]
    if strand == '1':
        r2_to_extract = '%s.2' % perfectly_mapped_read
        if r2_to_extract not in all_mapped_reads_set:
        #if r2_to_extract not in clipping_mapped_reads_list:
            solely_perfectly_mapped_reads_r2.add(r2_to_extract)

    elif strand == '2':
        r1_to_extract = '%s.1' % perfectly_mapped_read
        if r1_to_extract not in all_mapped_reads_set:
        #if r1_to_extract not in clipping_mapped_reads_list:
            solely_perfectly_mapped_reads_r1.add(r1_to_extract)


print(len(solely_perfectly_mapped_reads_r1) + len(solely_perfectly_mapped_reads_r2))
# default: 1318
# new with 0.005 error rate: 1504


