import pandas as pd
from Bio import SeqIO


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


def get_cigar_aln_len(cigar_splitted):
    # aligned_len: M I X =
    aligned_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]
        # get aligned_len
        if each_part_cate in {'M', 'm', 'I', 'i', 'X', 'x', '='}:
            aligned_len += each_part_len
    return aligned_len


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


class MappingRecord:

    def __init__(self):

        #################### overall ####################

        self.qualified_reads = False

        #################### round 1 16s ####################

        self.consider_r1_unmapped_mate = False
        self.consider_r2_unmapped_mate = False

        self.r1_16s_ref_dict = dict()
        self.r2_16s_ref_dict = dict()

        self.r1_16s_refs_lowest_mismatch = None
        self.r2_16s_refs_lowest_mismatch = None

        self.r1_16s_refs_no_ignored = dict()
        self.r2_16s_refs_no_ignored = dict()
        self.shared_16s_refs_no_ignored = dict()

        self.both_mapped_to_16s = False

        #################### round 1 ctg ####################

        self.r1_ctg_ref_dict = dict()
        self.r2_ctg_ref_dict = dict()

        self.r1_ctg_refs_lowest_mismatch = None
        self.r2_ctg_refs_lowest_mismatch = None

        self.r1_ctg_refs_no_ignored = dict()
        self.r2_ctg_refs_no_ignored = dict()
        self.shared_ctg_refs_no_ignored = dict()

        self.matched_to_ctg = False

        #################### round 2 ####################

        self.qualified_reads_rd2 = False

        self.r1_ctg_ref_dict_rd2 = dict()
        self.r2_ctg_ref_dict_rd2 = dict()

        self.r1_ctg_refs_lowest_mismatch_rd2 = None
        self.r2_ctg_refs_lowest_mismatch_rd2 = None

        #################### round 2 mini_assembly ####################

        self.r1_mini_ref_dict = dict()
        self.r2_mini_ref_dict = dict()

        self.r1_mini_refs_lowest_mismatch = None
        self.r2_mini_refs_lowest_mismatch = None

        self.r1_mini_refs_no_ignored = dict()
        self.r2_mini_refs_no_ignored = dict()
        self.shared_mini_refs_no_ignored = dict()


class LinkingRecord:

    def __init__(self):

        self.linked_seq_l = ''
        self.linked_seq_r = ''

        self.linked_seq_len_l = 0
        self.linked_seq_len_r = 0

        self.linking_reads_base = []

        self.linking_reads_l = []
        self.linking_reads_r = []

        self.linking_cigar_l = []
        self.linking_cigar_r = []

        self.linking_pos_l = []
        self.linking_pos_r = []

        self.min_dist_to_end_l = []
        self.min_dist_to_end_r = []


def get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len):

    mismatch_set_all_cigar = set()
    mismatch_set_long_M_cigars = set()
    for each_cigar in r1_ref_cigar_set:
        aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
        mismatch_set_all_cigar.add(mismatch_pct)
        if aligned_len >= min_M_len:
            mismatch_set_long_M_cigars.add(mismatch_pct)

    min_mismatch = 'NA'
    if len(mismatch_set_all_cigar) > 0:
        min_mismatch = min(mismatch_set_all_cigar)
        if len(mismatch_set_long_M_cigars) > 0:
            min_mismatch = min(mismatch_set_long_M_cigars)

    return min_mismatch


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


def check_cigar_quality(cigar_str, consider_ref_min_mismatch, ref_min_mismatch, mismatch_cutoff, min_M_len, ref_pos, marker_len):

    # check the following:
    # 1. both end clip  2. mismatch     3. aligned length   4. clp in the middle    5. ref_min_mismatch

    cigar_splitted = cigar_splitter(cigar_str)
    qualified_cigar = False

    # check both end clip
    both_end_clp = check_both_ends_clipping(cigar_splitted)
    if both_end_clp is False:
        # check mismatch
        aln_len, aln_pct, clp_len, clp_pct, mismatch_pct = get_cigar_stats(cigar_splitted)

        # check mismatch
        passed_mismatch_check = False
        if consider_ref_min_mismatch is False:
            if mismatch_pct <= mismatch_cutoff:
                passed_mismatch_check = True
        else:
            if ref_min_mismatch != 'NA':
                if (mismatch_pct <= ref_min_mismatch) and (mismatch_pct <= mismatch_cutoff):
                    passed_mismatch_check = True

        if passed_mismatch_check is True:
            # check aligned length
            if aln_len >= min_M_len:
                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar_str) or ('s' in cigar_str):
                    clip_in_middle = True
                    if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + aln_len - 1) == marker_len:
                            clip_in_middle = False

                if clip_in_middle is False:
                    qualified_cigar = True

    return qualified_cigar


def get_qualified_ref_dict(ref_dict, ref_len_dict, ref_min_mismatch, min_M_len, mismatch_cutoff, refs_to_ignore):

    refs_passed_qc = {}
    refs_passed_qc_with_pos = {}

    for each_ref in ref_dict:

        # one read can only mapped to on ref once
        matched_pos_dict = ref_dict[each_ref]
        if len(matched_pos_dict) > 1:
            refs_to_ignore.add(each_ref)
        else:
            ref_pos = list(matched_pos_dict.keys())[0]
            ref_cigar = matched_pos_dict[ref_pos]
            qualified_cigar = check_cigar_quality(ref_cigar, True, ref_min_mismatch, mismatch_cutoff, min_M_len, ref_pos, ref_len_dict[each_ref])

            if qualified_cigar is False:
                refs_to_ignore.add(each_ref)
            else:
                refs_passed_qc[each_ref] = [ref_cigar]
                refs_passed_qc_with_pos[each_ref] = {ref_pos: ref_cigar}

    return refs_passed_qc, refs_passed_qc_with_pos


def get_short_cigar_pct(cigar_list, short_M_len):
    short_cigar_num = 0
    for each_cigar in cigar_list:
        aln_len, aln_pct, clp_len, clp_pct, mis_pct = get_cigar_stats(cigar_splitter(each_cigar))
        if aln_len <= short_M_len:
            short_cigar_num += 1
    short_cigar_pct = short_cigar_num * 100 / len(cigar_list)
    short_cigar_pct = float("{0:.2f}".format(short_cigar_pct))
    return short_cigar_pct


def get_matched_pos_num_pct(linking_cigar_list, linking_pos_list, ref_len):

    matched_pos = set()
    for (each_linking_cigar, each_linking_pos) in zip(linking_cigar_list, linking_pos_list):
        current_cigar_aligned_len = get_cigar_aln_len(cigar_splitter(each_linking_cigar))
        matched_seq_end = each_linking_pos + current_cigar_aligned_len
        pos_list = list(range(each_linking_pos, matched_seq_end))
        matched_pos.update(pos_list)

    matched_pos_num = len(matched_pos)
    matched_pos_pct = matched_pos_num * 100 / ref_len
    matched_pos_pct = float("{0:.2f}".format(matched_pos_pct))

    return matched_pos_num, matched_pos_pct


def sort_csv_by_col(file_in, file_out, col_header):
    df_in        = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


def filter_linkages_iteratively_mini_assembly_to_ctg(file_in_sorted, min_linkages, file_out):

    # do mini-assemblies assigned to the same mag need to have roughly the same number of linkages? think about this later
    mag_ctg_max_link_num_dict = {}
    mini_assembly_to_mag_dict = {}
    mini_assembly_to_ctg_dict = {}
    file_out_handle = open(file_out, 'w')
    mini_assembly_with_assignment = set()
    mag_ctg_set_with_linked_mini_assembly = set()
    for each_match in open(file_in_sorted):
        if each_match.startswith('MiniAssembly,GenomicSeq,Number'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            mini_assembly = match_split[0]
            mag_ctg_id = match_split[1]
            mag_id = mag_ctg_id.split('___C___')[0]

            linkage_num = int(match_split[2])
            if linkage_num >= min_linkages:
                if mini_assembly not in mini_assembly_with_assignment:
                    if mag_ctg_id not in mag_ctg_max_link_num_dict:
                        mag_ctg_max_link_num_dict[mag_ctg_id] = linkage_num
                        file_out_handle.write(each_match)
                        mini_assembly_to_ctg_dict[mini_assembly] = mag_ctg_id
                        mini_assembly_to_mag_dict[mini_assembly] = mag_id
                        mini_assembly_with_assignment.add(mini_assembly)
                        mag_ctg_set_with_linked_mini_assembly.add(mag_ctg_id)
                    else:
                        ratio_with_best_assignment = linkage_num / (mag_ctg_max_link_num_dict[mag_ctg_id])
                        if ratio_with_best_assignment >= 0.8:
                            file_out_handle.write(each_match)
                            mini_assembly_to_ctg_dict[mini_assembly] = mag_ctg_id
                            mini_assembly_to_mag_dict[mini_assembly] = mag_id
                            mini_assembly_with_assignment.add(mini_assembly)
                            mag_ctg_set_with_linked_mini_assembly.add(mag_ctg_id)
                        else:
                            mini_assembly_with_assignment.add(mini_assembly)
    file_out_handle.close()

    return mini_assembly_to_ctg_dict, mini_assembly_to_mag_dict, mag_ctg_set_with_linked_mini_assembly


def check_clp_cigar_pct_diff(cigar_list_l, cigar_list_r, diff_cutoff):
    uneven_clp_cigar_pct = 'passed'

    clp_cigar_num_l = 0
    for cigar_l in cigar_list_l:
        if 'S' in cigar_l:
            clp_cigar_num_l += 1

    clp_cigar_num_r = 0
    for cigar_r in cigar_list_r:
        if 'S' in cigar_r:
            clp_cigar_num_r += 1

    clp_cigar_pct_l = clp_cigar_num_l * 100 / len(cigar_list_l)
    clp_cigar_pct_r = clp_cigar_num_r * 100 / len(cigar_list_r)
    clp_cigar_pct_l = float("{0:.2f}".format(clp_cigar_pct_l))
    clp_cigar_pct_r = float("{0:.2f}".format(clp_cigar_pct_r))

    if abs(clp_cigar_pct_l - clp_cigar_pct_r) > diff_cutoff:
        uneven_clp_cigar_pct = 'failed'

    return uneven_clp_cigar_pct, clp_cigar_pct_l, clp_cigar_pct_r


########################################################################################################################

min_M_len_mini                                      = 75
mismatch_cutoff                                     = 2
short_M_len                                         = 75
linked_to_ctg_end_cigar_num_cutoff                  = 1
linked_to_ctg_end_cigar_pct_cutoff                  = 20
consider_as_low_linking_reads_num                   = 100
max_short_cigar_pct_cutoff_linking_reads_num_high   = 75
max_short_cigar_pct_cutoff_linking_reads_num_low    = 85
end_seq_len                                         = 1000
uneven_clp_cigar_pct_cutoff                         = 25
mini_assembly_to_16s_ctg_connector                  = '___Mini___'

# input file
step_2_wd                                   = '/Users/songweizhi/Desktop/tunning_mini_to_mag'
mini_assemblies                             = '%s/file_in/scaffolds.fasta'                                  % step_2_wd
sam_file_mini_assembly_reformatted_sorted   = '%s/file_in/scaffolds_bowtie_sorted.sam'                      % step_2_wd
combined_1st_round_unlinked_mag_end_seq     = '%s/file_in/round_1_unlinked_gnm_end_%sbp.fa'                 % (step_2_wd, end_seq_len)
free_living_ctg_ref_file_with_pos_cigar     = '%s/file_in/round2_free_living_ctg_refs_with_pos_cigar.txt'   % step_2_wd
free_living_ctg_ref_file_no_linked          = '%s/file_in/round2_free_living_ctg_refs_no_linked.txt'        % step_2_wd

# output file
stats_mini_assembly_to_ctg_qc_report        = '%s/stats_mini_assembly_to_ctg_qc_report.txt'                 % step_2_wd
stats_mini_assembly_to_ctg                  = '%s/stats_mini_assembly_to_ctg.txt'                           % step_2_wd
stats_mini_assembly_to_ctg_sorted           = '%s/stats_mini_assembly_to_ctg_sorted.txt'                    % step_2_wd
stats_mini_assembly_to_ctg_filtered         = '%s/stats_mini_assembly_to_ctg_filtered.txt'                  % step_2_wd

# NODE_21_length_536_cov_0.956416,Oral_70___C___NODE_5490_length_8437_cov_0.862575_l,67

########################################################################################################################

mini_assembly_len_dict = {}
MappingRecord_dict_mini = {}
current_read_base = ''
current_read_base_r1_mini_ref_dict = dict()
current_read_base_r2_mini_ref_dict = dict()
with open(sam_file_mini_assembly_reformatted_sorted) as sam_file_mini_assembly_reformatted_sorted_opened:
    for each_line in sam_file_mini_assembly_reformatted_sorted_opened:
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('@'):
            mini_assembly_id = ''
            mini_assembly_len = 0
            for each_element in each_line_split:
                if each_element.startswith('SN:'):
                    mini_assembly_id = each_element[3:]
                if each_element.startswith('LN:'):
                    mini_assembly_len = int(each_element[3:])
            mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
        else:
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
                        current_read_base_r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                    if read_strand == '2':
                        current_read_base_r2_mini_ref_dict[ref_id] = {ref_pos: cigar}

            elif read_id_base == current_read_base:

                if cigar != '*':
                    if read_strand == '1':
                        if ref_id not in current_read_base_r1_mini_ref_dict:
                            current_read_base_r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                        else:
                            current_read_base_r1_mini_ref_dict[ref_id][ref_pos] = cigar
                    if read_strand == '2':
                        if ref_id not in current_read_base_r2_mini_ref_dict:
                            current_read_base_r2_mini_ref_dict[ref_id] = {ref_pos: cigar}
                        else:
                            current_read_base_r2_mini_ref_dict[ref_id][ref_pos] = cigar
            else:
                ################################### analysis previous read refs ####################################

                mini_refs_to_ignore = set()

                ########## get lowest mismatch for r1/r2 mini refs ##########

                # get r1_ref_cigar_set
                r1_ref_cigar_set = set()
                for each_pos_dict in current_read_base_r1_mini_ref_dict.values():
                    each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                    r1_ref_cigar_set.update(each_pos_dict_values)

                # get r2_ref_cigar_set
                r2_ref_cigar_set = set()
                for each_pos_dict in current_read_base_r2_mini_ref_dict.values():
                    each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
                    r2_ref_cigar_set.update(each_pos_dict_values)

                r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_mini)
                r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_mini)

                ########## filter mini refs for r1 and r2 ##########

                r1_mini_refs_passed_qc, r1_mini_refs_passed_qc_with_pos = get_qualified_ref_dict(
                    current_read_base_r1_mini_ref_dict, mini_assembly_len_dict, r1_ref_min_mismatch,
                    min_M_len_mini, mismatch_cutoff, mini_refs_to_ignore)

                r2_mini_refs_passed_qc, r2_mini_refs_passed_qc_with_pos = get_qualified_ref_dict(
                    current_read_base_r2_mini_ref_dict, mini_assembly_len_dict, r2_ref_min_mismatch,
                    min_M_len_mini, mismatch_cutoff, mini_refs_to_ignore)

                ####################################################################################################

                r1_mini_refs_no_ignored = {key: value for key, value in r1_mini_refs_passed_qc.items() if
                                           key not in mini_refs_to_ignore}
                r2_mini_refs_no_ignored = {key: value for key, value in r2_mini_refs_passed_qc.items() if
                                           key not in mini_refs_to_ignore}

                # no mate has no_ignored alignments
                if (len(r1_mini_refs_no_ignored) == 0) and (len(r2_mini_refs_no_ignored) == 0):
                    pass

                # only r1 has no_ignored alignments
                elif (len(r1_mini_refs_no_ignored) > 0) and (len(r2_mini_refs_no_ignored) == 0):
                    if current_read_base not in MappingRecord_dict_mini:
                        MappingRecord_dict_mini[current_read_base] = MappingRecord()
                    MappingRecord_dict_mini[current_read_base].r1_mini_ref_dict = current_read_base_r1_mini_ref_dict
                    MappingRecord_dict_mini[current_read_base].r1_mini_refs_no_ignored = r1_mini_refs_no_ignored

                # only r2 has no_ignored alignments
                elif (len(r1_mini_refs_no_ignored) == 0) and (len(r2_mini_refs_no_ignored) > 0):
                    if current_read_base not in MappingRecord_dict_mini:
                        MappingRecord_dict_mini[current_read_base] = MappingRecord()
                    MappingRecord_dict_mini[current_read_base].r2_mini_ref_dict = current_read_base_r2_mini_ref_dict
                    MappingRecord_dict_mini[current_read_base].r2_mini_refs_no_ignored = r2_mini_refs_no_ignored

                # both r1 and r2 have no_ignored alignments
                else:
                    shared_mini_refs_no_ignored = {
                        key: [r1_mini_refs_no_ignored[key][0], r2_mini_refs_no_ignored[key][0]] for key in
                        set(r1_mini_refs_no_ignored).intersection(set(r2_mini_refs_no_ignored))}
                    if len(shared_mini_refs_no_ignored) > 0:
                        if current_read_base not in MappingRecord_dict_mini:
                            MappingRecord_dict_mini[current_read_base] = MappingRecord()
                        MappingRecord_dict_mini[current_read_base].r1_mini_ref_dict = current_read_base_r1_mini_ref_dict
                        MappingRecord_dict_mini[current_read_base].r2_mini_ref_dict = current_read_base_r2_mini_ref_dict
                        MappingRecord_dict_mini[
                            current_read_base].shared_mini_refs_no_ignored = shared_mini_refs_no_ignored

                ########################################### reset values ###########################################

                current_read_base = read_id_base
                current_read_base_r1_mini_ref_dict = dict()
                current_read_base_r2_mini_ref_dict = dict()

                if cigar != '*':
                    if read_strand == '1':
                        current_read_base_r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                    if read_strand == '2':
                        current_read_base_r2_mini_ref_dict[ref_id] = {ref_pos: cigar}

#################### read sequences into dict ####################

# read sequence of mini_assembly into dict
mini_assembly_seq_dict = {}
for linked_ctg in SeqIO.parse(mini_assemblies, 'fasta'):
    mini_assembly_seq_dict[linked_ctg.id] = str(linked_ctg.seq)

# read sequence of unlinked mag end seq into dict
unlinked_mag_end_seq_dict = {}
for linked_ctg in SeqIO.parse(combined_1st_round_unlinked_mag_end_seq, 'fasta'):
    unlinked_mag_end_seq_dict[linked_ctg.id] = str(linked_ctg.seq)

############################################ link mini-assembly to MAGs ############################################

# read ctg side aln_pos and cigar into dict
mini_ctg_side_pos_dict = dict()
mini_ctg_side_cigar_dict = dict()
for each_read_to_ctg_ref in open(free_living_ctg_ref_file_with_pos_cigar):
    each_read_to_ctg_ref_split = each_read_to_ctg_ref.strip().split('\t')
    read_id = each_read_to_ctg_ref_split[0]
    ctg_refs = each_read_to_ctg_ref_split[1].split(',')
    for each_ctg_ref in ctg_refs:
        each_ctg_ref_split = each_ctg_ref.split('__pc__')
        ctg_ref_id = each_ctg_ref_split[0]
        ctg_ref_pos = int(each_ctg_ref_split[1])
        ctg_ref_cigar = each_ctg_ref_split[2]
        read_to_ctg_ref_key = '%s__ctg__%s' % (read_id, ctg_ref_id)
        mini_ctg_side_pos_dict[read_to_ctg_ref_key] = ctg_ref_pos
        mini_ctg_side_cigar_dict[read_to_ctg_ref_key] = ctg_ref_cigar

mini_to_mag_LinkingRecord_dict = {}
for free_living_read_ctg in open(free_living_ctg_ref_file_no_linked):
    free_living_read_ctg_split = free_living_read_ctg.strip().split('\t')
    read_id = free_living_read_ctg_split[0]
    read_id_base = '.'.join(read_id.split('.')[:-1])
    read_strand = read_id.split('.')[-1]
    read_ctg_refs = free_living_read_ctg_split[1].split(',')
    read_mini_mp = MappingRecord_dict_mini.get(read_id_base, None)
    if read_mini_mp is not None:

        read_mini_refs = set()
        for mini_ref in read_mini_mp.shared_mini_refs_no_ignored:
            read_mini_refs.add(mini_ref)
        if read_strand == '1':
            for mini_ref in read_mini_mp.r1_mini_refs_no_ignored:
                read_mini_refs.add(mini_ref)
        if read_strand == '2':
            for mini_ref in read_mini_mp.r2_mini_refs_no_ignored:
                read_mini_refs.add(mini_ref)

        for each_mini_ref in read_mini_refs:
            for each_ctg_ref in read_ctg_refs:
                mini_to_ctg_key_with_l_r = '%s%s%s' % (each_mini_ref, mini_assembly_to_16s_ctg_connector, each_ctg_ref)
                mini_ref_len = mini_assembly_len_dict[each_mini_ref]
                ctg_ref_len = len(unlinked_mag_end_seq_dict[each_ctg_ref])
                if mini_to_ctg_key_with_l_r not in mini_to_mag_LinkingRecord_dict:
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r] = LinkingRecord()
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_l = each_mini_ref
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_r = each_ctg_ref
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_len_l = mini_ref_len
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linked_seq_len_r = ctg_ref_len

                # get mini side cigar here
                current_mini_mp = MappingRecord_dict_mini[read_id_base]
                current_linking_cigar = ''
                current_linking_pos = ''
                if read_strand == '1':
                    current_linking_cigar = list(current_mini_mp.r1_mini_ref_dict[each_mini_ref].values())[0]
                    current_linking_pos = list(current_mini_mp.r1_mini_ref_dict[each_mini_ref].keys())[0]
                if read_strand == '2':
                    current_linking_cigar = list(current_mini_mp.r2_mini_ref_dict[each_mini_ref].values())[0]
                    current_linking_pos = list(current_mini_mp.r2_mini_ref_dict[each_mini_ref].keys())[0]

                # get min dist to mini end
                current_linking_cigar_aln_len = get_cigar_aln_len(cigar_splitter(current_linking_cigar))
                current_linking_cigar_dist_to_left = current_linking_pos - 1
                current_linking_cigar_dist_to_right = mini_ref_len - current_linking_pos - current_linking_cigar_aln_len + 1
                min_dist_to_mini_end = min(current_linking_cigar_dist_to_left, current_linking_cigar_dist_to_right)
                mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_reads_base.append(read_id_base)
                mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_cigar_l.append(current_linking_cigar)
                mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_pos_l.append(current_linking_pos)
                mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].min_dist_to_end_l.append(min_dist_to_mini_end)

                # get ctg side cigar here
                read_to_ctg_ref_key = '%s__ctg__%s' % (read_id, each_ctg_ref)
                current_linking_cigar_ctg_side = mini_ctg_side_cigar_dict.get(read_to_ctg_ref_key, None)
                current_linking_pos_ctg_side = mini_ctg_side_pos_dict.get(read_to_ctg_ref_key, None)
                if current_linking_cigar_ctg_side is not None:
                    # get min dist to ctg end
                    current_linking_cigar_ctg_side_aln_len = get_cigar_aln_len(
                        cigar_splitter(current_linking_cigar_ctg_side))
                    current_linking_cigar_ctg_side_dist_to_left = current_linking_pos_ctg_side - 1
                    current_linking_cigar_ctg_side_dist_to_right = ctg_ref_len - current_linking_pos_ctg_side - current_linking_cigar_ctg_side_aln_len + 1
                    min_dist_to_mini_end_ctg_side = min(current_linking_cigar_ctg_side_dist_to_left,
                                                        current_linking_cigar_ctg_side_dist_to_right)
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_cigar_r.append(
                        current_linking_cigar_ctg_side)
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].linking_pos_r.append(
                        current_linking_pos_ctg_side)
                    mini_to_mag_LinkingRecord_dict[mini_to_ctg_key_with_l_r].min_dist_to_end_r.append(
                        min_dist_to_mini_end_ctg_side)

link_to_check = ['NODE_299_length_313_cov_1.121053___Mini___Oral_48___C___NODE_18526_length_3167_cov_0.810855_r',
                 'NODE_552_length_278_cov_0.225806___Mini___Oral_35___C___NODE_7603_length_6638_cov_0.038857_l',
                 'NODE_380_length_299_cov_0.613636___Mini___Oral_50___C___NODE_8973_length_5826_cov_0.485699_r',
                 'NODE_157_length_409_cov_0.146853___Mini___Oral_54___C___NODE_5721_length_8207_cov_0.036015_r']


all_linking_reads_base_set_rd2_mini_to_ctg = set()
stats_mini_assembly_to_ctg_qc_report_handle = open(stats_mini_assembly_to_ctg_qc_report, 'w')
stats_mini_assembly_to_ctg_qc_report_handle.write('Mini\tContig\tLinkages\tLinked_to_mini_end(num pct)\tLinked_to_ctg_end(num pct)\tshort_cigar_pct_mini\tshort_cigar_pct_ctg\tclp_cigar_pct\n')
stats_mini_assembly_to_ctg_handle = open(stats_mini_assembly_to_ctg, 'w')
stats_mini_assembly_to_ctg_handle.write('MiniAssembly,GenomicSeq,Number\n')
mini_assembly_to_ctg_dict_with_l_r_read_base_min3 = {}
for each_link in mini_to_mag_LinkingRecord_dict.copy():
    current_linkage_linking_read_base = mini_to_mag_LinkingRecord_dict[each_link].linking_reads_base
    if len(current_linkage_linking_read_base) < 3:
        mini_to_mag_LinkingRecord_dict.pop(each_link)
    else:
        id_mini_assembly = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_l
        id_ctg = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_r
        current_linkage_cigar_mini_side = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_l
        current_linkage_cigar_ctg_side = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_r

        # get pct of short aligned cigar
        short_cigar_pct_mini = get_short_cigar_pct(current_linkage_cigar_mini_side, short_M_len)
        short_cigar_pct_ctg = get_short_cigar_pct(current_linkage_cigar_ctg_side, short_M_len)

        # check num/pct of reads linked to mini end
        min_dist_list_to_mini_end = mini_to_mag_LinkingRecord_dict[each_link].min_dist_to_end_l
        linked_to_mini_end_cigar_num = 0
        for each_min_dist in min_dist_list_to_mini_end:
            if each_min_dist <= 10:
                linked_to_mini_end_cigar_num += 1
        linked_to_mini_end_cigar_pct = linked_to_mini_end_cigar_num * 100 / len(min_dist_list_to_mini_end)
        linked_to_mini_end_cigar_pct = float("{0:.2f}".format(linked_to_mini_end_cigar_pct))

        # check num/pct of reads linked to ctg end
        min_dist_list_to_ctg_end = mini_to_mag_LinkingRecord_dict[each_link].min_dist_to_end_r
        linked_to_ctg_end_cigar_num = 0
        for each_min_dist in min_dist_list_to_ctg_end:
            if each_min_dist <= 10:
                linked_to_ctg_end_cigar_num += 1
        linked_to_ctg_end_cigar_pct = linked_to_ctg_end_cigar_num * 100 / len(min_dist_list_to_ctg_end)
        linked_to_ctg_end_cigar_pct = float("{0:.2f}".format(linked_to_ctg_end_cigar_pct))

        uneven_clp_cigar_pct, clp_cigar_pct_mini_side, clp_cigar_pct_ctg_side = check_clp_cigar_pct_diff(current_linkage_cigar_mini_side, current_linkage_cigar_ctg_side, uneven_clp_cigar_pct_cutoff)
        uneven_clp_cigar_pct_report_str = '%s(%s_vs_%s)' % (uneven_clp_cigar_pct, clp_cigar_pct_mini_side, clp_cigar_pct_ctg_side)

        if each_link in link_to_check:
            print(each_link)
            print('current_linkage_cigar_mini_side (%s): %s' % (len(current_linkage_cigar_mini_side), current_linkage_cigar_mini_side))
            print('current_linkage_cigar_ctg_side (%s): %s' % (len(current_linkage_cigar_ctg_side), current_linkage_cigar_ctg_side))
            print('uneven_clp_cigar_pct: %s' % uneven_clp_cigar_pct_report_str)

        linked_to_mini_end = False
        if (linked_to_mini_end_cigar_num > 0) and (linked_to_mini_end_cigar_pct >= 5):
            linked_to_mini_end = True

        linked_to_ctg_end = False
        if (linked_to_ctg_end_cigar_num >= linked_to_ctg_end_cigar_num_cutoff) and (linked_to_ctg_end_cigar_pct >= linked_to_ctg_end_cigar_pct_cutoff):
            linked_to_ctg_end = True

        max_short_cigar_pct_cutoff_to_use_mini = max_short_cigar_pct_cutoff_linking_reads_num_high
        if len(current_linkage_cigar_mini_side) < consider_as_low_linking_reads_num:
            max_short_cigar_pct_cutoff_to_use_mini = max_short_cigar_pct_cutoff_linking_reads_num_low

        max_short_cigar_pct_cutoff_to_use_ctg = max_short_cigar_pct_cutoff_linking_reads_num_high
        if len(current_linkage_cigar_ctg_side) < consider_as_low_linking_reads_num:
            max_short_cigar_pct_cutoff_to_use_ctg = max_short_cigar_pct_cutoff_linking_reads_num_low

        short_cigar_pct_mini_bool = False
        if short_cigar_pct_mini < max_short_cigar_pct_cutoff_to_use_mini:
            short_cigar_pct_mini_bool = True

        short_cigar_pct_ctg_bool = False
        if short_cigar_pct_ctg < max_short_cigar_pct_cutoff_to_use_ctg:
            short_cigar_pct_ctg_bool = True

        # get bp and pct of matched positions
        linking_cigar_list_mini = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_l
        linking_pos_list_mini = mini_to_mag_LinkingRecord_dict[each_link].linking_pos_l
        ref_len_mini = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_len_l
        matched_pos_mini_num, matched_pos_mini_pct = get_matched_pos_num_pct(linking_cigar_list_mini, linking_pos_list_mini, ref_len_mini)

        linking_cigar_list_ctg = mini_to_mag_LinkingRecord_dict[each_link].linking_cigar_r
        linking_pos_list_ctg = mini_to_mag_LinkingRecord_dict[each_link].linking_pos_r
        ref_len_ctg = mini_to_mag_LinkingRecord_dict[each_link].linked_seq_len_r
        matched_pos_ctg_num, matched_pos_ctg_pct = get_matched_pos_num_pct(linking_cigar_list_ctg, linking_pos_list_ctg, ref_len_ctg)

        matched_region_passed_qc_mini = False
        if (matched_pos_mini_num >= 300) or (matched_pos_mini_pct >= 50):
            matched_region_passed_qc_mini = True

        matched_region_passed_qc_ctg = False
        if matched_pos_ctg_num >= 300:
            matched_region_passed_qc_ctg = True

        # write out qc
        stats_mini_assembly_to_ctg_qc_report_handle.write('%s\t%s\t%s\t%s(%s %s)__v__%s(%sbp %s)\t%s(%s %s)__v__%s(%sbp)\t%s(%s)\t%s(%s)\t%s\n' % (
                id_mini_assembly, id_ctg, len(current_linkage_linking_read_base),
                linked_to_mini_end, linked_to_mini_end_cigar_num, linked_to_mini_end_cigar_pct,
                matched_region_passed_qc_mini, matched_pos_mini_num, matched_pos_mini_pct,
                linked_to_ctg_end, linked_to_ctg_end_cigar_num, linked_to_ctg_end_cigar_pct,
                matched_region_passed_qc_ctg, matched_pos_ctg_num,
                short_cigar_pct_mini_bool, short_cigar_pct_mini,
                short_cigar_pct_ctg_bool, short_cigar_pct_ctg,
                uneven_clp_cigar_pct_report_str))
        if each_link in link_to_check:
            print('%s\t%s\t%s\t%s(%s %s)__v__%s(%sbp %s)\t%s(%s %s)__v__%s(%sbp)\t%s(%s)\t%s(%s)\t%s\n' % (
                    id_mini_assembly, id_ctg, len(current_linkage_linking_read_base),
                    linked_to_mini_end, linked_to_mini_end_cigar_num, linked_to_mini_end_cigar_pct,
                    matched_region_passed_qc_mini, matched_pos_mini_num, matched_pos_mini_pct,
                    linked_to_ctg_end, linked_to_ctg_end_cigar_num, linked_to_ctg_end_cigar_pct,
                    matched_region_passed_qc_ctg, matched_pos_ctg_num,
                    short_cigar_pct_mini_bool, short_cigar_pct_mini,
                    short_cigar_pct_ctg_bool, short_cigar_pct_ctg,
                    uneven_clp_cigar_pct_report_str))

        if uneven_clp_cigar_pct == 'passed':
            if (short_cigar_pct_mini_bool is True) and (short_cigar_pct_ctg_bool is True):
                if ((linked_to_mini_end is True) or (matched_region_passed_qc_mini is True)) and (linked_to_ctg_end is True):
                    mini_assembly_to_ctg_dict_with_l_r_read_base_min3[each_link] = current_linkage_linking_read_base
                    all_linking_reads_base_set_rd2_mini_to_ctg.update(current_linkage_linking_read_base)
                    stats_mini_assembly_to_ctg_handle.write('%s,%s,%s\n' % (id_mini_assembly, id_ctg, len(current_linkage_linking_read_base)))
stats_mini_assembly_to_ctg_handle.close()
stats_mini_assembly_to_ctg_qc_report_handle.close()

# sort and filter
sort_csv_by_col(stats_mini_assembly_to_ctg, stats_mini_assembly_to_ctg_sorted, 'Number')
mini_assembly_to_ctg_dict, mini_assembly_to_mag_dict, mag_ctg_set_with_linked_mini_assembly = filter_linkages_iteratively_mini_assembly_to_ctg(
    stats_mini_assembly_to_ctg_sorted, 3,
    stats_mini_assembly_to_ctg_filtered)
