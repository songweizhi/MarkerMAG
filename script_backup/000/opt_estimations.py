import numpy as np
from Bio import SeqIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def vxtractor_output_to_dict(vxtractor_op, seq_len_dict, end_len_to_ignore, ignore_region_end_pct, min_region_len):

    seq_to_check = ''

    # initialize var con region id list
    var_region_id_list     = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9']
    con_region_id_list     = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10']
    var_con_region_id_list = ['C1', 'V1', 'C2', 'V2', 'C3', 'V3', 'C4', 'V4', 'C5', 'V5', 'C6', 'V6', 'C7', 'V7', 'C8', 'V8', 'C9', 'V9', 'C10']

    # read in variable regions
    seq_var_con_region_dict = {}
    for each_line in open(vxtractor_op):
        if not each_line.startswith(('Command:', 'Options:', 'Sequence,')):
            each_line_no_quote_symbol = ''
            for each_char in each_line.strip():
                if each_char not in ['"', '\'', '']:
                    each_line_no_quote_symbol += each_char

            each_line_split = each_line_no_quote_symbol.strip().split(',')
            seq_id = each_line_split[0]
            seq_var_con_region_dict[seq_id] = {}
            seq_var_regions = each_line_split[37:46]
            # print('%s\t%s' % (seq_id, '\t'.join(seq_var_regions)))

            for (var_id, var_range) in zip(var_region_id_list, seq_var_regions):
                var_range_to_store = 'na'
                if var_range not in ['', 'notfound']:
                    if 'wrongorder' not in var_range:
                        if 'HMM=' not in var_range:
                            var_range_to_store = sorted([int(i) for i in var_range.split('-')])
                seq_var_con_region_dict[seq_id][var_id] = var_range_to_store
                # print('%s\t%s' % (var_id, var_range))

    # get conserved regions
    for each_seq in seq_var_con_region_dict:
        current_var_con_dict = seq_var_con_region_dict.get(each_seq, dict())
        for each_con_index in range(1, 11):

            ##### get C1 #####
            if each_con_index == 1:
                if current_var_con_dict['V1'] == 'na':
                    seq_var_con_region_dict[each_seq]['C1'] = 'na'
                else:
                    c1_left = 1
                    c1_right = current_var_con_dict['V1'][0] - 1
                    if c1_right >= 1:
                        seq_var_con_region_dict[each_seq]['C1'] = [c1_left, c1_right]
                    else:
                        seq_var_con_region_dict[each_seq]['C1'] = 'na'

            ##### get C2 - C9 #####
            elif 1 < each_con_index < 10:
                con_id = 'C%s' % each_con_index
                var_id_left = 'V%s' % (each_con_index - 1)
                var_id_right = 'V%s' % each_con_index
                if (current_var_con_dict[var_id_left] == 'na') or (current_var_con_dict[var_id_right] == 'na'):
                    seq_var_con_region_dict[each_seq][con_id] = 'na'
                else:
                    con_left = current_var_con_dict[var_id_left][1] + 1
                    con_right = current_var_con_dict[var_id_right][0] - 1
                    if con_right - con_left >= 1:
                        seq_var_con_region_dict[each_seq][con_id] = [con_left, con_right]
                    else:
                        seq_var_con_region_dict[each_seq][con_id] = 'na'

            ##### get C10 #####
            else:
                if current_var_con_dict['V9'] == 'na':
                    seq_var_con_region_dict[each_seq]['C10'] = 'na'
                else:
                    c10_left = current_var_con_dict['V9'][1] + 1
                    c10_right = seq_len_dict.get(each_seq, 0)
                    if c10_right - c10_left >= 1:
                        seq_var_con_region_dict[each_seq]['C10'] = [c10_left, c10_right]
                    else:
                        seq_var_con_region_dict[each_seq]['C10'] = 'na'

        ##############################
        # print('\n--------------------\n')
        # for each_var_con in var_con_region_id_list:
        #     var_con_range = 'na'
        #     if each_var_con[0] == 'V':
        #         var_con_range = current_var_con_dict.get(each_var_con, 'na')
        #     if each_var_con[0] == 'C':
        #         var_con_range = seq_var_con_region_dict[each_seq].get(each_var_con, 'na')
        #     print('%s\t%s\t%s' % (each_seq, each_var_con, var_con_range))

    # filter seq_var_con_region_dict
    seq_var_con_region_dict_filtered = {}
    for each_seq in seq_var_con_region_dict:
        current_seq_len = seq_len_dict.get(each_seq, 0)
        seq_var_con_region_dict_filtered[each_seq] = {}
        current_var_con_region_dict = seq_var_con_region_dict[each_seq]
        # print('%s\t%s' % (each_seq, current_var_con_region_dict))
        for each_con_var in var_con_region_id_list:
            con_var_range = current_var_con_region_dict[each_con_var]

            if con_var_range == 'na':
                seq_var_con_region_dict_filtered[each_seq][each_con_var] = 'na'
            else:
                pos_l = con_var_range[0]
                pos_r = con_var_range[1]
                region_len = pos_r - pos_l + 1
                region_end_len_to_ignore = round(region_len * ignore_region_end_pct / 100)
                pos_l_no_ignored = pos_l + region_end_len_to_ignore
                pos_r_no_ignored = pos_r - region_end_len_to_ignore

                region_to_consider = False
                pos_l_no_ignored_no_seq_end = pos_l_no_ignored
                pos_r_no_ignored_no_seq_end = pos_r_no_ignored
                if pos_l_no_ignored < end_len_to_ignore:
                    pos_l_no_ignored_no_seq_end = end_len_to_ignore + 1
                    if (pos_r_no_ignored - pos_l_no_ignored_no_seq_end) >= min_region_len:
                        region_to_consider = True
                elif pos_r_no_ignored > (current_seq_len - end_len_to_ignore):
                    pos_r_no_ignored_no_seq_end = current_seq_len - end_len_to_ignore
                    if (pos_r_no_ignored_no_seq_end - pos_l_no_ignored) >= min_region_len:
                        region_to_consider = True
                else:
                    if (pos_r_no_ignored - pos_l_no_ignored + 1) >= min_region_len:
                        region_to_consider = True

                if region_to_consider is False:
                    if each_seq == seq_to_check:
                        print('%s\t%s\t%s-%s\t%s\t%s-%s\t%s\t%s' % (each_seq, each_con_var, end_len_to_ignore, (current_seq_len - end_len_to_ignore), con_var_range, pos_l_no_ignored_no_seq_end, pos_r_no_ignored_no_seq_end, (pos_r_no_ignored_no_seq_end - pos_l_no_ignored_no_seq_end + 1), 'bad'))
                    seq_var_con_region_dict_filtered[each_seq][each_con_var] = 'na'
                else:
                    if each_seq == seq_to_check:
                        print('%s\t%s\t%s-%s\t%s\t%s-%s\t%s\t%s' % (each_seq, each_con_var, end_len_to_ignore, (current_seq_len - end_len_to_ignore), con_var_range, pos_l_no_ignored_no_seq_end, pos_r_no_ignored_no_seq_end, (pos_r_no_ignored_no_seq_end - pos_l_no_ignored_no_seq_end + 1), 'good'))
                    seq_var_con_region_dict_filtered[each_seq][each_con_var] = [pos_l_no_ignored_no_seq_end, pos_r_no_ignored_no_seq_end]

    # ignore region if its two flankings are na
    seq_var_con_region_dict_filtered2 = {}
    for each_seq in seq_var_con_region_dict_filtered:
        current_var_con_region_dict = seq_var_con_region_dict_filtered[each_seq]
        current_var_con_region_dict2 = {}

        # check vars
        for each_var in var_region_id_list:
            var_index = int(each_var[1:])
            flk_con_l_id = 'C%s' % var_index
            flk_con_r_id = 'C%s' % (var_index + 1)
            flk_con_l = current_var_con_region_dict.get(flk_con_l_id, 'na')
            flk_con_r = current_var_con_region_dict.get(flk_con_r_id, 'na')
            if (flk_con_l == 'na') and (flk_con_r == 'na'):
                current_var_con_region_dict2[each_var] = 'na'
            else:
                current_var_con_region_dict2[each_var] = current_var_con_region_dict[each_var]

        # check cons
        for each_con in con_region_id_list:
            con_index = int(each_con[1:])
            flk_var_l_id = 'V%s' % (con_index - 1)
            flk_var_r_id = 'V%s' % (con_index)
            flk_var_l = current_var_con_region_dict.get(flk_var_l_id, 'na')
            flk_var_r = current_var_con_region_dict.get(flk_var_r_id, 'na')
            if (flk_var_l == 'na') and (flk_var_r == 'na'):
                current_var_con_region_dict2[each_con] = 'na'
            else:
                current_var_con_region_dict2[each_con] = current_var_con_region_dict[each_con]

        seq_var_con_region_dict_filtered2[each_seq] = current_var_con_region_dict2

    seq_var_con_combined_pos_dict = {}
    for each_seq in seq_var_con_region_dict_filtered2:
        current_var_con_region_dict = seq_var_con_region_dict_filtered2[each_seq]
        seq_var_con_combined_pos_dict[each_seq] = {}
        seq_var_con_combined_pos_dict[each_seq]['Con'] = set()
        seq_var_con_combined_pos_dict[each_seq]['Var'] = set()
        for each_con_var in current_var_con_region_dict:
            current_var_con_range = current_var_con_region_dict[each_con_var]
            if current_var_con_range != 'na':

                # get var pos list
                if each_con_var in var_region_id_list:
                    current_var_pos_list = list(range((current_var_con_range[0] - 1), current_var_con_range[1]))
                    for each_var_pos in current_var_pos_list:
                        seq_var_con_combined_pos_dict[each_seq]['Var'].add(each_var_pos)

                # get con pos list
                if each_con_var in con_region_id_list:
                    current_con_pos_list = list(range((current_var_con_range[0] - 1), current_var_con_range[1]))
                    for each_con_pos in current_con_pos_list:
                        seq_var_con_combined_pos_dict[each_seq]['Con'].add(each_con_pos)

    return seq_var_con_combined_pos_dict


def opt_estimation(depth_file_all, depth_file_linked, pos_cov_file_all, pos_cov_file_linked, vxtractor_op, end_len_to_ignore, ignore_region_end_pct, min_region_len_no_ignored, seq_len_dict, depth_file_opt):

    # read in mean depth (all 16S)
    marker_to_mag_dict = {}
    mean_depth_mag_dict = {}
    mean_depth_16s_all = {}
    mean_depth_16s_all_norm = {}
    for each_16s in open(depth_file_all):
        if not each_16s.startswith('MAG	16S	MAG_depth	16S_depth'):
            each_16s_split = each_16s.strip().split('\t')
            id_mag = each_16s_split[0]
            id_16s = each_16s_split[1]
            mag_depth = float(each_16s_split[2])
            mean_depth = float(each_16s_split[3])
            mean_depth_norm = float(each_16s_split[4])
            mean_depth_mag_dict[id_mag] = mag_depth
            mean_depth_16s_all[id_16s] = mean_depth
            mean_depth_16s_all_norm[id_16s] = mean_depth_norm
            marker_to_mag_dict[id_16s] = id_mag

    # read in mean depth (linked 16S)
    mean_depth_16s_linked = {}
    mean_depth_16s_linked_norm = {}
    for each_16s in open(depth_file_linked):
        if not each_16s.startswith('MAG	16S	MAG_depth	16S_depth'):
            each_16s_split = each_16s.strip().split('\t')
            id_16s = each_16s_split[1]
            mean_depth = float(each_16s_split[3])
            mean_depth_norm = float(each_16s_split[4])
            mean_depth_16s_linked[id_16s] = mean_depth
            mean_depth_16s_linked_norm[id_16s] = mean_depth_norm

    # read in depth by bp (all 16S)
    depth_by_bp_dict_all = {}
    for each_bp_all in open(pos_cov_file_all):
        each_bp_all_split = each_bp_all.strip().split('\t')
        ref_if = each_bp_all_split[0]
        pos = each_bp_all_split[1]
        depth = int(each_bp_all_split[2])
        if ref_if not in depth_by_bp_dict_all:
            depth_by_bp_dict_all[ref_if] = {}
        depth_by_bp_dict_all[ref_if][pos] = depth

    # read in depth by bp (linked 16S)
    depth_by_bp_dict_linked = {}
    for each_bp_all in open(pos_cov_file_linked):
        each_bp_all_split = each_bp_all.strip().split('\t')
        ref_if = each_bp_all_split[0]
        pos = each_bp_all_split[1]
        depth = int(each_bp_all_split[2])
        if ref_if not in depth_by_bp_dict_linked:
            depth_by_bp_dict_linked[ref_if] = {}
        depth_by_bp_dict_linked[ref_if][pos] = depth

    seq_var_con_combined_pos_dict = vxtractor_output_to_dict(vxtractor_op, seq_len_dict, end_len_to_ignore, ignore_region_end_pct, min_region_len_no_ignored)

    estimation_with_depth_opt_handle = open(depth_file_opt, 'w')
    estimation_with_depth_opt_handle.write('MAG\t16S\tCoverage(var/con)(linked)\tCoverage(var/con)(all)\tCoverage_opt\tCopies_opt\tCopies_opt_norm\n')
    for each_seq in seq_var_con_combined_pos_dict:
        seq_var_pos_list = sorted([i for i in seq_var_con_combined_pos_dict[each_seq]['Var']])
        seq_con_pos_list = sorted([i for i in seq_var_con_combined_pos_dict[each_seq]['Con']])

        seq_var_pos_total_depth_all = 0
        for each_var_pos in seq_var_pos_list:
            var_pos_depth_all = depth_by_bp_dict_all.get(each_seq, dict()).get(str(each_var_pos), 0)
            seq_var_pos_total_depth_all += var_pos_depth_all

        seq_con_pos_total_depth_all = 0
        for each_con_pos in seq_con_pos_list:
            con_pos_depth_all = depth_by_bp_dict_all.get(each_seq, dict()).get(str(each_con_pos), 0)
            seq_con_pos_total_depth_all += con_pos_depth_all

        seq_var_pos_total_depth_linked = 0
        for each_var_pos in seq_var_pos_list:
            var_pos_depth_linked = depth_by_bp_dict_linked.get(each_seq, dict()).get(str(each_var_pos), 0)
            seq_var_pos_total_depth_linked += var_pos_depth_linked

        seq_con_pos_total_depth_linked = 0
        for each_con_pos in seq_con_pos_list:
            con_pos_depth_linked = depth_by_bp_dict_linked.get(each_seq, dict()).get(str(each_con_pos), 0)
            seq_con_pos_total_depth_linked += con_pos_depth_linked

        var_mean_depth_all = 0
        var_mean_depth_linked = 0
        if len(seq_var_pos_list) != 0:
            var_mean_depth_all = seq_var_pos_total_depth_all/len(seq_var_pos_list)
            var_mean_depth_linked = seq_var_pos_total_depth_linked/len(seq_var_pos_list)

        con_mean_depth_all = 0
        con_mean_depth_linked = 0
        if len(seq_con_pos_list) != 0:
            con_mean_depth_all = seq_con_pos_total_depth_all/len(seq_con_pos_list)
            con_mean_depth_linked = seq_con_pos_total_depth_linked/len(seq_con_pos_list)

        var_to_con_cov_ratio_all = 0
        if con_mean_depth_all > 0:
            var_to_con_cov_ratio_all = var_mean_depth_all/con_mean_depth_all
            var_to_con_cov_ratio_all = float("{0:.2f}".format(var_to_con_cov_ratio_all))

        var_to_con_cov_ratio_linked = 0
        if con_mean_depth_linked > 0:
            var_to_con_cov_ratio_linked = var_mean_depth_linked/con_mean_depth_linked
            var_to_con_cov_ratio_linked = float("{0:.2f}".format(var_to_con_cov_ratio_linked))

        if (var_to_con_cov_ratio_linked != 0) and (var_to_con_cov_ratio_all != 0):

            seq_cov_by_linked      = mean_depth_16s_linked.get(each_seq, 0)
            seq_cov_by_linked_norm = mean_depth_16s_linked_norm.get(each_seq, 0)
            seq_cov_by_all         = mean_depth_16s_all.get(each_seq, 0)
            seq_cov_by_all_norm    = mean_depth_16s_all_norm.get(each_seq, 0)
            linked_mag             = marker_to_mag_dict[each_seq]
            mag_cov                = mean_depth_mag_dict[linked_mag]

            # chose the higher one
            if (0.85 <= var_to_con_cov_ratio_linked <= 1.15) and (0.85 <= var_to_con_cov_ratio_all <= 1.15):
                seq_cov_opt = max([seq_cov_by_linked, seq_cov_by_all])
                seq_cov_opt_norm = max([seq_cov_by_linked_norm, seq_cov_by_all_norm])

            # reads from other genome mapped to conserved regions
            elif (var_to_con_cov_ratio_linked < 0.85) and (var_to_con_cov_ratio_all < 0.85):
                seq_cov_opt = seq_cov_by_all*var_to_con_cov_ratio_all
                seq_cov_opt_norm = seq_cov_by_all_norm*var_to_con_cov_ratio_all

            # reads from conserved regions mapped to other genomes
            elif (var_to_con_cov_ratio_linked > 1.15) and (var_to_con_cov_ratio_all > 1.15):
                seq_cov_opt = seq_cov_by_all/var_to_con_cov_ratio_all
                seq_cov_opt_norm = seq_cov_by_all_norm/var_to_con_cov_ratio_all
            else:
                distance_to_1_linked = abs(var_to_con_cov_ratio_linked - 1)
                distance_to_1_all = abs(var_to_con_cov_ratio_all - 1)
                seq_cov_opt = seq_cov_by_linked/var_to_con_cov_ratio_linked
                seq_cov_opt_norm = seq_cov_by_linked_norm/var_to_con_cov_ratio_linked
                if distance_to_1_linked > distance_to_1_all:
                    seq_cov_opt = seq_cov_by_all/var_to_con_cov_ratio_all
                    seq_cov_opt_norm = seq_cov_by_all_norm/var_to_con_cov_ratio_all
                seq_cov_opt = float("{0:.2f}".format(seq_cov_opt))
                seq_cov_opt_norm = float("{0:.2f}".format(seq_cov_opt_norm))

            seq_cp_num_opt = float("{0:.2f}".format(seq_cov_opt/mag_cov))
            seq_cp_num_opt_norm = float("{0:.2f}".format(seq_cov_opt_norm/mag_cov))
            estimation_with_depth_opt_handle.write('%s\t%s\t%s(%s)\t%s(%s)\t%s\t%s\t%s\n' % (linked_mag, each_seq, seq_cov_by_linked, var_to_con_cov_ratio_linked, seq_cov_by_all, var_to_con_cov_ratio_all, seq_cov_opt, seq_cp_num_opt, seq_cp_num_opt_norm))
    estimation_with_depth_opt_handle.close()


def get_scatter_plot(num_list_1, num_list_2, png_file):

    num_arrary_1 = np.array(num_list_1)
    num_arrary_2 = np.array(num_list_2)

    fig = plt.figure(figsize=(6, 6))
    plt.margins(0)

    plt.scatter(num_arrary_1, num_arrary_2)
    plt.xlabel("Estimated copy number", fontsize=12)
    plt.ylabel("User provided copy number", fontsize=12)

    # set axis range
    max_value = max([max(num_list_1), max(num_list_2)])
    plt.xlim(0, round(max_value + 1))
    plt.ylim(0, round(max_value + 1))

    # add fit line
    coeffs = np.polyfit(num_arrary_1, num_arrary_2, 1)
    slope = coeffs[0]
    intercept = coeffs[1]
    plt.plot(num_arrary_1, slope * num_arrary_1 + intercept)

    # get R-squared value
    p = np.poly1d(coeffs)
    yhat = p(num_arrary_1)
    ybar = np.sum(num_arrary_2)/len(num_arrary_2)
    ssreg = np.sum((yhat-ybar)**2)
    sstot = np.sum((num_arrary_2 - ybar)**2)
    r_squared = ssreg/sstot

    title_str = 'y = %sx + %s, R-squared = %s' % (float("{0:.4f}".format(slope)), float("{0:.4f}".format(intercept)), float("{0:.4f}".format(r_squared)))
    plt.title(title_str)

    # save plot
    plt.tight_layout()
    plt.savefig(png_file)
    plt.close()


# # MBARC26
# depth_by_bp_file_all            = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_16S_2_50_PairedOnly_depth_by_bp_all.txt'
# depth_by_bp_file_linked         = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_16S_2_50_PairedOnly_depth_by_bp_linked.txt'
# seq_16s_file                    = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_SILVA138_polished_polished_min1200bp_c99.0.fasta'
# vxtractor_op                    = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_SILVA138_polished_polished_min1200bp_c99_variable_region.csv'
# estimation_with_depth_all       = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_with_depth_all.txt'
# estimation_with_depth_linked    = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_with_depth_linked.txt'
# provided_cp_num_txt             = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_ref_16S_cp_num.txt'
# estimation_with_depth_opt       = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_optimized.txt'
# estimation_by_gnm               = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.txt'
# scatter_plot                    = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.png'
# scatter_plot_norm               = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG_norm.png'
#

# # GI
# depth_by_bp_file_all            = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_16S_2_50_PairedOnly_depth_by_bp_all.txt'
# depth_by_bp_file_linked         = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_16S_2_50_PairedOnly_depth_by_bp_linked.txt'
# seq_16s_file                    = '/Users/songweizhi/Desktop/get_var_region_wd/GI_128_16S_0.999_polished_min1200bp_c99.0.fasta'
# vxtractor_op                    = '/Users/songweizhi/Desktop/get_var_region_wd/GI_128_16S_0.999_polished_min1200bp_c99_variable_region_as_bacteria.csv'
# estimation_with_depth_all       = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_with_depth_all.txt'
# estimation_with_depth_linked    = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_with_depth_linked.txt'
# provided_cp_num_txt             = '/Users/songweizhi/Desktop/get_var_region_wd/GI_provided_mag_16s_cp_num.txt'
# estimation_with_depth_opt       = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_optimized.txt'
# estimation_by_gnm               = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.txt'
# scatter_plot                    = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.png'
# scatter_plot_norm               = '/Users/songweizhi/Desktop/get_var_region_wd/GI_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG_norm.png'


# # Oral
# depth_by_bp_file_all            = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_16S_2_50_PairedOnly_depth_by_bp_all.txt'
# depth_by_bp_file_linked         = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_16S_2_50_PairedOnly_depth_by_bp_linked.txt'
# seq_16s_file                    = '/Users/songweizhi/Desktop/get_var_region_wd/CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.0.fa'
# vxtractor_op                    = '/Users/songweizhi/Desktop/get_var_region_wd/CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99_variable_region_as_bacteria.csv'
# estimation_with_depth_all       = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_with_depth_all.txt'
# estimation_with_depth_linked    = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_with_depth_linked.txt'
# provided_cp_num_txt             = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_provided_mag_16s_cp_num.txt'
# estimation_with_depth_opt       = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_optimized.txt'
# estimation_by_gnm               = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.txt'
# scatter_plot                    = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.png'
# scatter_plot_norm               = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG_norm.png'



# Oral
depth_by_bp_file_all            = '/Users/songweizhi/Desktop/opt/test_sub20_16S_pos_cov_all.txt'
depth_by_bp_file_linked         = '/Users/songweizhi/Desktop/opt/test_sub20_16S_pos_cov_linked.txt'
seq_16s_file                    = '/Users/songweizhi/Desktop/opt/CAMI_Oral_138_16S_0.999.polished_min1200_polished_min1200bp_c99.0.fa'
vxtractor_op                    = '/Users/songweizhi/Desktop/opt/test_3_25_25_16S_vxtractor.csv'
estimation_with_depth_all       = '/Users/songweizhi/Desktop/opt/test_3_25_25_MAG_16S_cov_all.txt'
estimation_with_depth_linked    = '/Users/songweizhi/Desktop/opt/test_3_25_25_MAG_16S_cov_linked.txt'
provided_cp_num_txt             = '/Users/songweizhi/Desktop/get_var_region_wd/Oral_provided_mag_16s_cp_num.txt'

estimation_with_depth_opt       = '/Users/songweizhi/Desktop/opt/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_optimized.txt'

estimation_by_gnm               = '/Users/songweizhi/Desktop/opt/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.txt'
scatter_plot                    = '/Users/songweizhi/Desktop/opt/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG.png'
scatter_plot_norm               = '/Users/songweizhi/Desktop/opt/Oral_0906_MAGs_2_50__16S_2_50__ign_ed150_ld25.0_hd25.0__PairedOnly_estimation_by_MAG_norm.png'



end_len_to_ignore           = 200
ignore_region_end_pct       = 25
min_region_len_no_ignored   = 25

marker_len_dict = {}
for each_marker in SeqIO.parse(seq_16s_file, 'fasta'):
    marker_id = each_marker.id
    marker_seq = str(each_marker.seq).upper()
    marker_len_dict[marker_id] = len(marker_seq)


opt_estimation(estimation_with_depth_all, estimation_with_depth_linked, depth_by_bp_file_all, depth_by_bp_file_linked, vxtractor_op, end_len_to_ignore, ignore_region_end_pct, min_region_len_no_ignored, marker_len_dict, estimation_with_depth_opt)


mag_cp_num_dict = {}
mag_cp_num_dict_norm = {}
for each_16s in open(estimation_with_depth_opt):
    if not each_16s.startswith('MAG	16S	Coverage(var/con)(linked)'):
        each_16s_split = each_16s.strip().split('\t')
        mag_id = each_16s_split[0]
        cp_num = float(each_16s_split[5])
        cp_num_norm = float(each_16s_split[6])
        if mag_id not in mag_cp_num_dict:
            mag_cp_num_dict[mag_id] = cp_num
            mag_cp_num_dict_norm[mag_id] = cp_num_norm
        else:
            mag_cp_num_dict[mag_id] += cp_num
            mag_cp_num_dict_norm[mag_id] += cp_num_norm


reported_ref_16s_copy_num_dict = {}
if provided_cp_num_txt is not None:
    for each_ref in open(provided_cp_num_txt):
        if not each_ref.startswith('Genome\t'):
            each_ref_split = each_ref.strip().split('\t')
            ref_id = each_ref_split[0]
            reported_copy_num = float(each_ref_split[1])
            reported_ref_16s_copy_num_dict[ref_id] = reported_copy_num


estimated_16s_cp_num_list = []
estimated_16s_cp_num_list_norm = []
reported_16s_cp_num_list = []
estimation_by_gnm_handle = open(estimation_by_gnm, 'w')
if provided_cp_num_txt is None:
    estimation_by_gnm_handle.write('MAG\tCopies\tCopies_norm\n')
else:
    estimation_by_gnm_handle.write('MAG\tCopies\tCopies_norm\tProvided\n')
for each_mag in mag_cp_num_dict:
    mag_cp_num = mag_cp_num_dict.get(each_mag, 'na')
    mag_cp_num_norm = mag_cp_num_dict_norm.get(each_mag, 'na')
    mag_cp_num_provided = reported_ref_16s_copy_num_dict.get(each_mag, 'na')

    # write out
    if provided_cp_num_txt is None:
        estimation_by_gnm_handle.write('%s\t%s\t%s\n' % (each_mag, mag_cp_num, mag_cp_num_norm))
    else:
        estimation_by_gnm_handle.write(
            '%s\t%s\t%s\t%s\n' % (each_mag, mag_cp_num, mag_cp_num_norm, mag_cp_num_provided))

    # get cp num list
    estimated_16s_cp_num_list.append(mag_cp_num)
    estimated_16s_cp_num_list_norm.append(mag_cp_num_norm)
    reported_16s_cp_num_list.append(mag_cp_num_provided)
estimation_by_gnm_handle.close()

get_scatter_plot(estimated_16s_cp_num_list, reported_16s_cp_num_list, scatter_plot)
get_scatter_plot(estimated_16s_cp_num_list_norm, reported_16s_cp_num_list, scatter_plot_norm)
