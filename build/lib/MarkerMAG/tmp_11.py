
def cigar_splitter():
    pass

def check_both_ends_clipping(aaa):
    pass

def cigar_splitter(cifar):
    return []

def get_cigar_stats():
    pass

aa = 0




def check_cigar_quality(cigar_str, mismatch_cutoff, min_M_len, ref_pos, ref_len):

    r2_ctg_ref_cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(r2_ctg_ref_cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        r2_aligned_len_ctg, r2_aligned_pct_ctg, r2_clipping_len_ctg, r2_clipping_pct_ctg, r2_mismatch_pct_ctg = get_cigar_stats(r2_ctg_ref_cigar_splitted)
        if r2_mismatch_pct_ctg <= mismatch_cutoff:
            # check aligned length
            if r2_aligned_len_ctg >= min_M_len:
                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar_str) or ('s' in cigar_str):
                    clip_in_middle = True
                    if (r2_ctg_ref_cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (r2_ctg_ref_cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + r2_aligned_len_ctg - 1) == ref_len:
                            clip_in_middle = False

                if clip_in_middle is False:
                    qualified_cigar = True

    return qualified_cigar


