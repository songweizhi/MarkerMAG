
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

sam_file_mini_assembly_combined     = '/Users/songweizhi/Desktop/round2/scaffolds_combined.sam'
min_M_pct                           = 20
min_M_len                           = 30
min_clp_len                         = 30
min_clp_M_len                       = 25
max_mis_pct                         = 3
MappingRecord_dict = {}
round_2_MappingRecord_dict = {}


# remove reads mapped to multiple miniassembly? check later
gap_seq_to_reads_dict = {}
for each_read in open(sam_file_mini_assembly_combined):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            ref_id  = each_read_split[2]
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
            if (aligned_len >= min_M_len) and (aligned_pct >= min_M_pct) and (mismatch_pct <= (max_mis_pct)):
                if ref_id not in gap_seq_to_reads_dict:
                    gap_seq_to_reads_dict[ref_id] = [read_id]
                else:
                    gap_seq_to_reads_dict[ref_id].append(read_id)

stats_GapFilling_file_16s = ''
stats_GapFilling_file_ctg = ''
stats_GapFilling_file_16s_handle = open(stats_GapFilling_file_16s, 'w')
stats_GapFilling_file_ctg_handle = open(stats_GapFilling_file_ctg, 'w')
stats_GapFilling_file_16s_handle.write('Gap_seq,s16,Number\n')
stats_GapFilling_file_ctg_handle.write('Gap_seq,ctg,Number\n')
for gap_seq in gap_seq_to_reads_dict:
    gap_seq_mapped_reads = gap_seq_to_reads_dict[gap_seq]
    gap_seq_to_16s_linkage_dict = {}
    gap_seq_to_ctg_linkage_dict = {}
    gap_seq_to_gnm_linkage_dict = {}
    for mapped_read in gap_seq_mapped_reads:
        mapped_read_basename = mapped_read[:-2]
        mapped_read_strand   = mapped_read[-1]

        mapped_read_mate_or_clp_refs_16s = set()
        if mapped_read_basename in MappingRecord_dict:
            read_mr_round_1 = MappingRecord_dict[mapped_read_basename]
            if read_mr_round_1.consider_round_2 is True:
                if mapped_read_strand == '1':
                    if read_mr_round_1.consider_r1_unmapped_mate is True:
                        mapped_read_mate_or_clp_refs_16s = mapped_read_mate_or_clp_refs_16s.union(read_mr_round_1.r2_filtered_refs)
                    if read_mr_round_1.consider_r1_clipping_part is True:
                        mapped_read_mate_or_clp_refs_16s = mapped_read_mate_or_clp_refs_16s.union(read_mr_round_1.r1_filtered_refs)
                if mapped_read_strand == '2':
                    if read_mr_round_1.consider_r2_unmapped_mate is True:
                        mapped_read_mate_or_clp_refs_16s = mapped_read_mate_or_clp_refs_16s.union(read_mr_round_1.r1_filtered_refs)
                    if read_mr_round_1.consider_r2_clipping_part is True:
                        mapped_read_mate_or_clp_refs_16s = mapped_read_mate_or_clp_refs_16s.union(read_mr_round_1.r2_filtered_refs)

        for each_mate_or_clp_ref_16s in mapped_read_mate_or_clp_refs_16s:
            if each_mate_or_clp_ref_16s not in gap_seq_to_16s_linkage_dict:
                gap_seq_to_16s_linkage_dict[each_mate_or_clp_ref_16s] = 1
            else:
                gap_seq_to_16s_linkage_dict[each_mate_or_clp_ref_16s] += 1

        mapped_read_mate_or_clp_refs_ctg = set()
        if mapped_read_basename in round_2_MappingRecord_dict:
            read_mr_round_2 = round_2_MappingRecord_dict[mapped_read_basename]
            if mapped_read_strand == '1':
                if read_mr_round_2.consider_r1_unmapped_mate is True:
                    mapped_read_mate_or_clp_refs_ctg = mapped_read_mate_or_clp_refs_ctg.union(read_mr_round_2.r2_filtered_refs)
                if read_mr_round_2.consider_r1_clipping_part is True:
                    mapped_read_mate_or_clp_refs_ctg = mapped_read_mate_or_clp_refs_ctg.union(read_mr_round_2.r1_filtered_refs)
            if mapped_read_strand == '2':
                if read_mr_round_2.consider_r2_unmapped_mate is True:
                    mapped_read_mate_or_clp_refs_ctg = mapped_read_mate_or_clp_refs_ctg.union(read_mr_round_2.r1_filtered_refs)
                if read_mr_round_2.consider_r2_clipping_part is True:
                    mapped_read_mate_or_clp_refs_ctg = mapped_read_mate_or_clp_refs_ctg.union(read_mr_round_2.r2_filtered_refs)

        for each_mate_or_clp_ref_ctg in mapped_read_mate_or_clp_refs_ctg:
            if each_mate_or_clp_ref_ctg not in gap_seq_to_ctg_linkage_dict:
                gap_seq_to_ctg_linkage_dict[each_mate_or_clp_ref_ctg] = 1
            else:
                gap_seq_to_ctg_linkage_dict[each_mate_or_clp_ref_ctg] += 1

        for each_16s in gap_seq_to_16s_linkage_dict:
            stats_GapFilling_file_16s_handle.write('%s,%s,%s\n' % (gap_seq, each_16s, gap_seq_to_16s_linkage_dict[each_16s]))
        for each_ctg in gap_seq_to_ctg_linkage_dict:
            stats_GapFilling_file_ctg_handle.write('%s,%s,%s\n' % (gap_seq, each_ctg, gap_seq_to_ctg_linkage_dict[each_ctg]))

stats_GapFilling_file_16s_handle.close()
stats_GapFilling_file_ctg_handle.close()




