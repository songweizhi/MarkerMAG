import os
import multiprocessing as mp


def parse_sam16s_worker(argument_list):

    sorted_sam                  = argument_list[0]
    MappingRecord_dict          = argument_list[1]
    min_M_len_16s = argument_list[2]
    mismatch_cutoff = argument_list[3]
    marker_len_dict = argument_list[4]


    to_extract_read_base = set()
    current_read_base = ''
    current_read_base_r1_16s_ref_dict = dict()
    current_read_base_r2_16s_ref_dict = dict()
    with open(sorted_sam) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if not each_read.startswith('@'):
                each_read_split = each_read.strip().split('\t')
                cigar = each_read_split[5]
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_read_split[2]
                ref_pos = int(each_read_split[3])

                if current_read_base == '':
                    current_read_base = read_id_base

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}

                elif read_id_base == current_read_base:

                    if cigar != '*':
                        if read_strand == '1':
                            if ref_id not in current_read_base_r1_16s_ref_dict:
                                current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r1_16s_ref_dict[ref_id][ref_pos] = cigar
                        if read_strand == '2':
                            if ref_id not in current_read_base_r2_16s_ref_dict:
                                current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}
                            else:
                                current_read_base_r2_16s_ref_dict[ref_id][ref_pos] = cigar
                else:
                    ################################### analysis previous read refs ####################################

                    ########## get lowest mismatch for r1/r2 16s refs ##########

                    # get r1_ref_cigar_set
                    r1_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r1_16s_ref_dict.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r1_ref_cigar_set.update(each_pos_dict_values)

                    # get r2_ref_cigar_set
                    r2_ref_cigar_set = set()
                    for each_pos_dict in current_read_base_r2_16s_ref_dict.values():
                        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                        r2_ref_cigar_set.update(each_pos_dict_values)

                    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_16s)
                    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_16s)

                    refs_to_ignore = set()

                    ########## filter r1 16s refs ##########

                    r1_16s_refs_passed_qc = {}
                    for r1_16s_ref in current_read_base_r1_16s_ref_dict:
                        r1_matched_pos_dict = current_read_base_r1_16s_ref_dict[r1_16s_ref]

                        # one read need to mapped to one 16S only for one time
                        if len(r1_matched_pos_dict) > 1:
                            refs_to_ignore.add(r1_16s_ref)
                        else:
                            r1_16s_ref_pos = list(r1_matched_pos_dict.keys())[0]
                            r1_16s_ref_cigar = r1_matched_pos_dict[r1_16s_ref_pos]
                            r1_16s_ref_cigar_splitted = cigar_splitter(r1_16s_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r1_16s_ref_cigar_splitted)
                            if both_end_clp is True:
                                refs_to_ignore.add(r1_16s_ref)
                            else:
                                # check mismatch
                                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                                    r1_16s_ref_cigar_splitted)
                                if r1_ref_min_mismatch == 'NA':
                                    refs_to_ignore.add(r1_16s_ref)
                                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                                    refs_to_ignore.add(r1_16s_ref)
                                else:
                                    # check aligned length
                                    if r1_aligned_len < min_M_len_16s:
                                        refs_to_ignore.add(r1_16s_ref)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r1_16s_ref_cigar) or ('s' in r1_16s_ref_cigar):
                                            clip_in_middle = True
                                            if (r1_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r1_16s_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r1_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r1_16s_ref_pos + r1_aligned_len - 1) == marker_len_dict[r1_16s_ref]:
                                                    clip_in_middle = False

                                        # exclude the ref if clp in the middle is True
                                        if clip_in_middle is True:
                                            refs_to_ignore.add(r1_16s_ref)
                                        else:
                                            r1_16s_refs_passed_qc[r1_16s_ref] = [r1_16s_ref_cigar]

                    ########## filter r2 16s refs ##########

                    r2_16s_refs_passed_qc = {}
                    for r2_16s_ref in current_read_base_r2_16s_ref_dict:
                        r2_matched_pos_dict = current_read_base_r2_16s_ref_dict[r2_16s_ref]

                        # one read need to mapped to one 16S only once
                        if len(r2_matched_pos_dict) > 1:
                            refs_to_ignore.add(r2_16s_ref)
                        else:
                            r2_16s_ref_pos = list(r2_matched_pos_dict.keys())[0]
                            r2_16s_ref_cigar = r2_matched_pos_dict[r2_16s_ref_pos]
                            r2_16s_ref_cigar_splitted = cigar_splitter(r2_16s_ref_cigar)

                            # check both end clip
                            both_end_clp = check_both_ends_clipping(r2_16s_ref_cigar_splitted)
                            if both_end_clp is True:
                                refs_to_ignore.add(r2_16s_ref)
                            else:
                                # check mismatch
                                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                                    r2_16s_ref_cigar_splitted)
                                if r2_ref_min_mismatch == 'NA':
                                    refs_to_ignore.add(r2_16s_ref)
                                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                                    refs_to_ignore.add(r2_16s_ref)
                                else:
                                    # check aligned length
                                    if r2_aligned_len < min_M_len_16s:
                                        refs_to_ignore.add(r2_16s_ref)
                                    else:
                                        # check if clp in the middle
                                        clip_in_middle = False
                                        if ('S' in r2_16s_ref_cigar) or ('s' in r2_16s_ref_cigar):
                                            clip_in_middle = True
                                            if (r2_16s_ref_cigar_splitted[0][-1] in ['S', 's']) and (
                                                    r2_16s_ref_pos == 1):
                                                clip_in_middle = False
                                            if (r2_16s_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                                if (r2_16s_ref_pos + r2_aligned_len - 1) == marker_len_dict[r2_16s_ref]:
                                                    clip_in_middle = False

                                        # exclude the ref if clp in the middle is True
                                        if clip_in_middle is True:
                                            refs_to_ignore.add(r2_16s_ref)
                                        else:
                                            r2_16s_refs_passed_qc[r2_16s_ref] = [r2_16s_ref_cigar]

                    #################################

                    r1_16s_refs_no_ignored = {key: value for key, value in r1_16s_refs_passed_qc.items() if
                                              key not in refs_to_ignore}
                    r2_16s_refs_no_ignored = {key: value for key, value in r2_16s_refs_passed_qc.items() if
                                              key not in refs_to_ignore}

                    # no mate has no_ignored alignments
                    if (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) == 0):
                        pass

                    # only r1 has no_ignored alignments
                    elif (len(r1_16s_refs_no_ignored) > 0) and (len(r2_16s_refs_no_ignored) == 0):

                        MappingRecord_dict[current_read_base].qualified_reads = True
                        MappingRecord_dict[current_read_base].consider_r1_unmapped_mate = True
                        MappingRecord_dict[current_read_base].r1_16s_refs_no_ignored = r1_16s_refs_no_ignored
                        MappingRecord_dict[current_read_base].r1_16s_ref_dict = current_read_base_r1_16s_ref_dict
                        to_extract_read_base.add(current_read_base)

                    # only r2 has no_ignored alignments
                    elif (len(r1_16s_refs_no_ignored) == 0) and (len(r2_16s_refs_no_ignored) > 0):

                        MappingRecord_dict[current_read_base].qualified_reads = True
                        MappingRecord_dict[current_read_base].consider_r2_unmapped_mate = True
                        MappingRecord_dict[current_read_base].r2_16s_refs_no_ignored = r2_16s_refs_no_ignored
                        MappingRecord_dict[current_read_base].r2_16s_ref_dict = current_read_base_r2_16s_ref_dict
                        to_extract_read_base.add(current_read_base)

                    # both r1 and r2 have no_ignored alignments
                    else:
                        shared_16s_ref_dict = {key: [r1_16s_refs_no_ignored[key][0], r2_16s_refs_no_ignored[key][0]] for key
                                               in set(r1_16s_refs_no_ignored).intersection(set(r2_16s_refs_no_ignored))}
                        if len(shared_16s_ref_dict) > 0:

                            MappingRecord_dict[current_read_base].qualified_reads = True
                            MappingRecord_dict[current_read_base].both_mapped_to_16s = True
                            MappingRecord_dict[current_read_base].shared_16s_refs_no_ignored = shared_16s_ref_dict
                            MappingRecord_dict[current_read_base].r1_16s_ref_dict = current_read_base_r1_16s_ref_dict
                            MappingRecord_dict[current_read_base].r2_16s_ref_dict = current_read_base_r2_16s_ref_dict
                            to_extract_read_base.add(current_read_base)

                    ########################################### reset values ###########################################

                    current_read_base = read_id_base
                    current_read_base_r1_16s_ref_dict = dict()
                    current_read_base_r2_16s_ref_dict = dict()

                    if cigar != '*':
                        if read_strand == '1':
                            current_read_base_r1_16s_ref_dict[ref_id] = {ref_pos: cigar}
                        if read_strand == '2':
                            current_read_base_r2_16s_ref_dict[ref_id] = {ref_pos: cigar}


