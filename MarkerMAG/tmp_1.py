from Bio import SeqIO


Get_COG_boxplot = '''
cd /Users/songweizhi/Desktop
python3 ~/PycharmProjects/BioSAK/BioSAK/boxplot_matrix_COG.py -in COG_func_stats -out COG_func_stats_df.txt -skip_1st_row -with_functional_description -in_percent
# important !!!
# please make sure the calculated proportions of genes in each COG category among detected HGTs are in the last row in COG_func_stats_df.txt, as values in the last row will be used to position the coloured triangles and squares in the plot.
Rscript ~/PycharmProjects/BioSAK/BioSAK/COG_boxplot_last1row.R -i COG_func_stats_df.txt -o COG_func_stats_df.png
# step_1
python3 boxplot_matrix_COG.py -in COG_func_stats -out COG_func_stats_df.txt -skip_1st_row -with_functional_description -in_percent
# step_2, important !!!
# please make sure sumarries for detected HGTs are in the last row in COG_func_stats_df.txt, as values in the last row will be used to position the coloured triangles and squares in the plot.
# step_3
Rscript COG_boxplot_last1row.R -i COG_func_stats_df.txt -o COG_func_stats_df.png
'''


def get_regions_to_ignore_from_barrnap_output(combined_barrnap_gff, ctg_len_dict):

    ctg_ignore_region_dict = {}

    for each_line in open(combined_barrnap_gff):
        if (not each_line.startswith('#')) and ('16S_rRNA' in each_line):
            each_line_split = each_line.strip().split('\t')
            ctg_id = each_line_split[0]
            ctg_len = ctg_len_dict[ctg_id]
            start_pos = int(each_line_split[3])
            end_pos = int(each_line_split[4])
            len_16s = end_pos - start_pos + 1
            left_gap = start_pos - 1
            right_gap = ctg_len - end_pos - 1

            if left_gap <= 50:
                if ctg_id not in ctg_ignore_region_dict:
                    ctg_ignore_region_dict[ctg_id] = {'left_end'}
                else:
                    ctg_ignore_region_dict[ctg_id].add('left_end')

            if right_gap <= 50:
                if ctg_id not in ctg_ignore_region_dict:
                    ctg_ignore_region_dict[ctg_id] = {'right_end'}
                else:
                    ctg_ignore_region_dict[ctg_id].add('right_end')

    return ctg_ignore_region_dict


def get_unlinked_mag_end_seq(ref_in, ref_in_end_seq, end_seq_len, ctg_ignore_region_dict_rd1):
    ctg_ignore_region_dict_rd2 = dict()

    # get ref seqs subset
    ref_subset_handle = open(ref_in_end_seq, 'w')
    for ref_seq in SeqIO.parse(ref_in, 'fasta'):

        ref_seq_id = ref_seq.id
        ref_seq_len = len(ref_seq.seq)

        if ref_seq_len < end_seq_len * 2:
            ref_subset_handle.write('>%s\n' % ref_seq_id)
            ref_subset_handle.write('%s\n' % ref_seq.seq)

            # add to ctg_ignore_region_dict_rd2
            if ref_seq_id in ctg_ignore_region_dict_rd1:
                ctg_ignore_region_dict_rd2[ref_seq_id] = ctg_ignore_region_dict_rd1[ref_seq_id]
        else:
            ref_seq_left_end_id = '%s_l' % ref_seq_id
            ref_seq_right_end_id = '%s_r' % ref_seq_id
            ref_seq_left_end = ref_seq.seq[:end_seq_len]
            ref_seq_right_end = ref_seq.seq[-end_seq_len:]

            # write out left end
            ref_subset_handle.write('>%s\n' % ref_seq_left_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_left_end)

            # write out right end
            ref_subset_handle.write('>%s\n' % ref_seq_right_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_right_end)

            # add to ctg_ignore_region_dict_rd2
            if ref_seq_id in ctg_ignore_region_dict_rd1:
                current_seq_to_ignore_ends = ctg_ignore_region_dict_rd1[ref_seq_id]
                for end_to_ignore in current_seq_to_ignore_ends:
                    if end_to_ignore == 'left_end':
                        ctg_ignore_region_dict_rd2[ref_seq_left_end_id] = {'left_end'}
                    if end_to_ignore == 'right_end':
                        ctg_ignore_region_dict_rd2[ref_seq_right_end_id] = {'right_end'}
    ref_subset_handle.close()

    return ctg_ignore_region_dict_rd2


step_2_wd                                   = '/Users/songweizhi/Desktop/tunning_rd2'
combined_input_gnms                         = '%s/file_in/3_Oral_refined_MAGs_combined.fa'                  % step_2_wd
combined_input_gnms_no_ending_16s           = '%s/file_in/3_Oral_refined_MAGs_combined_no_ending_16s.fa'    % step_2_wd
end_seq_len = 500
combined_1st_round_unlinked_mags            = '%s/file_in/round_1_unlinked_gnm.fa'                          % step_2_wd
combined_1st_round_unlinked_mag_end_seq     = '%s/file_in/round_1_unlinked_gnm_end_%sbp.fa'                 % (step_2_wd, end_seq_len)
combined_barrnap_gff                        = '%s/file_in/combined_barrnap.gff'                             % step_2_wd


ctg_len_dict = {}
for each_ctg_record in SeqIO.parse(combined_input_gnms, 'fasta'):
    ctg_len_dict[each_ctg_record.id] = len(each_ctg_record.seq)

ctg_ended_16s_to_keep_part_dict = dict()
for each_line in open(combined_barrnap_gff):
    if (not each_line.startswith('#')) and ('16S_rRNA' in each_line):
        each_line_split = each_line.strip().split('\t')
        ctg_id = each_line_split[0]
        ctg_len = ctg_len_dict[ctg_id]
        start_pos = int(each_line_split[3])
        end_pos = int(each_line_split[4])
        len_16s = end_pos - start_pos + 1
        left_gap = start_pos - 1
        right_gap = ctg_len - end_pos - 1
        if left_gap <= 100:
            ctg_ended_16s_to_keep_part_dict[ctg_id] = '%s-%s' % ((end_pos + 1), ctg_len)
        elif right_gap <= 100:
            ctg_ended_16s_to_keep_part_dict[ctg_id] = '%s-%s' % (1, start_pos)

combined_input_gnms_no_ending_16s_handle = open(combined_input_gnms_no_ending_16s, 'w')
for each_ctg in SeqIO.parse(combined_input_gnms, 'fasta'):
    ctg_id = each_ctg.id
    ctg_seq = str(each_ctg.seq)

    if ctg_id not in ctg_ended_16s_to_keep_part_dict:
        combined_input_gnms_no_ending_16s_handle.write('>%s\n' % ctg_id)
        combined_input_gnms_no_ending_16s_handle.write('%s\n' % ctg_seq)
    else:
        keep_region_str = ctg_ended_16s_to_keep_part_dict[ctg_id].split('-')
        keep_region_l = int(keep_region_str[0])
        keep_region_r = int(keep_region_str[1])
        combined_input_gnms_no_ending_16s_handle.write('>%s\n' % ctg_id)
        combined_input_gnms_no_ending_16s_handle.write('%s\n' % ctg_seq[(keep_region_l - 1): keep_region_r])

combined_input_gnms_no_ending_16s_handle.close()

# ctg_ignore_region_dict = get_regions_to_ignore_from_barrnap_output(combined_barrnap_gff, ctg_len_dict)
# print(ctg_ignore_region_dict)
#
# ctg_ignore_region_dict_rd2 = get_unlinked_mag_end_seq(combined_1st_round_unlinked_mags,
#                                                       combined_1st_round_unlinked_mag_end_seq, end_seq_len,
#                                                       ctg_ignore_region_dict)
