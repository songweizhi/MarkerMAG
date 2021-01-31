def paired_blast_results_to_dict_by_mapping(unmapped_paired_reads_mapping_results):

    unmapped_paired_reads_to_ctg_dict_by_mapping = {}

    ref_len_dict = {}
    for unmapped_read in open(unmapped_paired_reads_mapping_results):

        # get ref len dict
        if unmapped_read.startswith('@'):
            unmapped_read_split = unmapped_read.strip().split('\t')
            ref_id = ''
            ref_len = 0
            for each_element in unmapped_read_split:
                if each_element.startswith('SN:'):
                    ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_len = int(each_element[3:])
            ref_len_dict[ref_id] = ref_len

        else:
            qualified_unmapped_read = False
            unmapped_read_split = unmapped_read.strip().split('\t')
            read_id = unmapped_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = unmapped_read_split[2]
            ref_len = ref_len_dict[ref_id]
            ref_id_with_prefix = 'GenomicSeq__%s' % unmapped_read_split[2]
            ref_pos = int(unmapped_read_split[3])
            cigar = unmapped_read_split[5]
            read_seq = unmapped_read_split[9]
            read_len = len(read_seq)
            cigar_splitted = cigar_splitter(cigar)

            # e.g. 189M
            if ('M' in cigar) and (len(cigar_splitted) == 1):
                qualified_unmapped_read = True

            elif len(cigar_splitted) == 2:

                # e.g. 139S61M
                if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M'):
                    cigar_M_pct = int(cigar_splitted[1][:-1]) * 100 / read_len
                    if (ref_pos == 1) and (cigar_M_pct >= 25):
                        qualified_unmapped_read = True

                # e.g. 147M53S
                if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[1][-1] == 'S'):
                    cigar_M_pct = int(cigar_splitted[0][:-1]) * 100 / read_len
                    matched_to_bp = ref_pos + int(cigar_splitted[0][:-1]) - 1
                    if (matched_to_bp == ref_len) and (cigar_M_pct >= 25):
                        qualified_unmapped_read = True

            elif len(cigar_splitted) == 3:

                # e.g. 121M1D66M, 121M1I66M
                if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M'):
                    mismatch_pct = int(cigar_splitted[1][:-1])*100/read_len
                    if mismatch_pct <= 1:
                        qualified_unmapped_read = True

                # e.g. 11S81M3S, not sure, not considered yet
                elif cigar_splitted[1][-1] == 'M':
                    min_mismatch = min(int(cigar_splitted[0][:-1]), int(cigar_splitted[2][:-1]))
                    if min_mismatch <= 3:
                        pass
            elif len(cigar_splitted) == 4:

                # e.g. ['181M', '1D', '8M', '11S'], ['181M', '2D', '8M', '11S']
                if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (cigar_splitted[3][-1] == 'S'):
                    cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1]))*100/read_len
                    if ((ref_pos + int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1])) == ref_len) and (int(cigar_splitted[1][:-1]) <= 2):
                        if cigar_M_pct >= 25:
                            qualified_unmapped_read = True

                # e.g. ['24S', '95M', '1D', '27M'] and ref_pos = 1
                elif (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (cigar_splitted[3][-1] == 'M'):
                    cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1]))*100/read_len
                    if (ref_pos == 1) and (int(cigar_splitted[2][:-1]) <= 1) and (cigar_M_pct >= 25):
                        qualified_unmapped_read = True

            # add to dict
            if qualified_unmapped_read is True:
                if read_id not in unmapped_paired_reads_to_ctg_dict_by_mapping:
                    unmapped_paired_reads_to_ctg_dict_by_mapping[read_id] = [ref_id_with_prefix]
                else:
                    unmapped_paired_reads_to_ctg_dict_by_mapping[read_id].append(ref_id_with_prefix)


    return unmapped_paired_reads_to_ctg_dict_by_mapping
