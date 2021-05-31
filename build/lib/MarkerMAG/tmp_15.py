
def parse_sam_gnm_worker(arguments_list):

    rd1_unlinked_mags_sam_bowtie_reformat_sorted    = arguments_list[0]
    free_living_ctg_ref_file                        = arguments_list[1]
    min_M_len_ctg                                   = arguments_list[2]
    mismatch_cutoff                                 = arguments_list[3]
    round_2_ctg_end_seq_len_dict                    = arguments_list[4]

    free_living_ctg_ref_file_handle = open(free_living_ctg_ref_file, 'w')
    to_extract_read_base_rd2_ctg = set()
    current_read_base = ''
    current_read_base_r1_ctg_ref_dict_rd2 = dict()
    current_read_base_r2_ctg_ref_dict_rd2 = dict()
    with open(rd1_unlinked_mags_sam_bowtie_reformat_sorted) as rd1_unlinked_mags_sam_bowtie_reformat_sorted_opened:
        for each_line in rd1_unlinked_mags_sam_bowtie_reformat_sorted_opened:
            if not each_line.startswith('@'):
                each_line_split = each_line.strip().split('\t')
                cigar = each_line_split[5]
                read_id = each_line_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_line_split[2]
                ref_pos = int(each_line_split[3])

                if current_read_base == '':
                    current_read_base = read_id_base

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}

                elif read_id_base == current_read_base:

                    if cigar != '*':
                        if read_strand == '1':
                            if ref_id not in current_read_base_r1_ctg_ref_dict_rd2:
                                current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r1_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar
                        if read_strand == '2':
                            if ref_id not in current_read_base_r2_ctg_ref_dict_rd2:
                                current_read_base_r2_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r2_ctg_ref_dict_rd2[ref_id][ref_pos] = cigar
                else:
                    ################################### analysis previous read refs ####################################

                    ctg_refs_to_ignore_rd2 = set()

                    ########## get lowest mismatch for r1/r2 ctg refs ##########

                    # get r1_ref_cigar_set
                    r1_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r1_ctg_ref_dict_rd2.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r1_ref_cigar_set.update(each_pos_dict_values)

                    # get r2_ref_cigar_set
                    r2_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r2_ctg_ref_dict_rd2.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r2_ref_cigar_set.update(each_pos_dict_values)

                    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_ctg)
                    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_ctg)

                    ########## filter r1 ctg refs ##########
                    r1_ctg_refs_passed_qc = {}
                    for r1_ctg_ref_rd2 in current_read_base_r1_ctg_ref_dict_rd2:
                        r1_matched_pos_dict = current_read_base_r1_ctg_ref_dict_rd2[r1_ctg_ref_rd2]
                        if len(r1_matched_pos_dict) > 1:
                            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                        else:
                            r1_ctg_ref_pos = list(r1_matched_pos_dict.keys())[0]
                            r1_ctg_ref_cigar = r1_matched_pos_dict[r1_ctg_ref_pos]
                            r1_ctg_ref_cigar_splitted = cigar_splitter(r1_ctg_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r1_ctg_ref_cigar_splitted)
                            if both_end_clp is True:
                                ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                            else:
                                # check mismatch
                                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                    r1_ctg_ref_cigar_splitted)
                                if r1_ref_min_mismatch == 'NA':
                                    ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                                    ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                else:
                                    # check aligned length
                                    if r1_aligned_len < min_M_len_ctg:
                                        ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r1_ctg_ref_cigar) or ('s' in r1_ctg_ref_cigar):
                                            clip_in_middle = True
                                            if (r1_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_ctg_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r1_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r1_ctg_ref_pos + r1_aligned_len - 1) == round_2_ctg_end_seq_len_dict[
                                                    r1_ctg_ref_rd2]:
                                                    clip_in_middle = False

                                        if clip_in_middle is True:
                                            ctg_refs_to_ignore_rd2.add(r1_ctg_ref_rd2)
                                        else:
                                            r1_ctg_refs_passed_qc[r1_ctg_ref_rd2] = [r1_ctg_ref_cigar]

                    ########## filter r2 ctg refs ##########
                    r2_ctg_refs_passed_qc = {}
                    for r2_ctg_ref_rd2 in current_read_base_r2_ctg_ref_dict_rd2:
                        r2_matched_pos_dict = current_read_base_r2_ctg_ref_dict_rd2[r2_ctg_ref_rd2]
                        if len(r2_matched_pos_dict) > 1:
                            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                        else:
                            r2_ctg_ref_pos = list(r2_matched_pos_dict.keys())[0]
                            r2_ctg_ref_cigar = r2_matched_pos_dict[r2_ctg_ref_pos]
                            r2_ctg_ref_cigar_splitted = cigar_splitter(r2_ctg_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
                            if both_end_clp is True:
                                ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                            else:
                                # check mismatch
                                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                    r2_ctg_ref_cigar_splitted)
                                if r2_ref_min_mismatch == 'NA':
                                    ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                                    ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                else:
                                    # check aligned length
                                    if r2_aligned_len < min_M_len_ctg:
                                        ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r2_ctg_ref_cigar) or ('s' in r2_ctg_ref_cigar):
                                            clip_in_middle = True
                                            if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_ctg_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r2_ctg_ref_pos + r2_aligned_len - 1) == round_2_ctg_end_seq_len_dict[r2_ctg_ref_rd2]:
                                                    clip_in_middle = False

                                        if clip_in_middle is True:
                                            ctg_refs_to_ignore_rd2.add(r2_ctg_ref_rd2)
                                        else:
                                            r2_ctg_refs_passed_qc[r2_ctg_ref_rd2] = [r2_ctg_ref_cigar]

                    ####################################################################################################

                    r1_ctg_refs_rd2_no_ignored = {key: value for key, value in r1_ctg_refs_passed_qc.items() if key not in ctg_refs_to_ignore_rd2}
                    r2_ctg_refs_rd2_no_ignored = {key: value for key, value in r2_ctg_refs_passed_qc.items() if key not in ctg_refs_to_ignore_rd2}

                    # only r1 has no_ignored alignments
                    if (len(r1_ctg_refs_rd2_no_ignored) > 0) and (len(r2_ctg_refs_rd2_no_ignored) == 0):
                        current_read_base___qualified_reads_rd2 = True
                        current_read_base___r1_ctg_ref_dict_rd2 = current_read_base_r1_ctg_ref_dict_rd2
                        to_extract_read_base_rd2_ctg.add(current_read_base)
                        free_living_ctg_ref_file_handle.write('%s\t%s\n' % (current_read_base, ','.join(r1_ctg_refs_rd2_no_ignored)))

                    # only r2 has no_ignored alignments
                    if (len(r1_ctg_refs_rd2_no_ignored) == 0) and (len(r2_ctg_refs_rd2_no_ignored) > 0):
                        current_read_base___qualified_reads_rd2 = True
                        current_read_base___r2_ctg_ref_dict_rd2 = current_read_base_r2_ctg_ref_dict_rd2
                        to_extract_read_base_rd2_ctg.add(current_read_base)
                        free_living_ctg_ref_file_handle.write('%s\t%s\n' % (current_read_base, ','.join(r2_ctg_refs_rd2_no_ignored)))

                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_ctg_ref_dict_rd2 = dict()
                    current_read_base_r2_ctg_ref_dict_rd2 = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r1_ctg_ref_dict_rd2[ref_id] = {ref_pos: cigar}

    free_living_ctg_ref_file_handle.close()

