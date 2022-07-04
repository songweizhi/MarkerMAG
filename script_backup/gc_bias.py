from Bio import SeqIO


'''
GC Bias Report
https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/GCBiasReport_fDG.htm
'''


gnm_id = 'FP'
gnm_seq_file    = '/Users/songweizhi/Desktop/gc_bias/%s.fna'                                        % gnm_id
depth_file      = '/Users/songweizhi/Desktop/gc_bias/%s_global_report_one_mis0_sorted_depth.txt'    % gnm_id
gc_bias_txt     = '/Users/songweizhi/Desktop/gc_bias/%s_gc_bias.txt'                                % gnm_id
window_len      = 100
ignore_ends_len = 300

sequence_id = ''
if gnm_id == 'FP':
    sequence_id = 'CP003260.1'
if gnm_id == 'DA':
    sequence_id = 'CP003639.1'
if gnm_id == 'NO':
    sequence_id = 'CP003929.1'

# get sequence
sequence_str = str(SeqIO.read(gnm_seq_file, 'fasta').seq)
print(len(sequence_str))

# read depth into dict
seq_pos_depth_dict = {}
seq_overall_depth_dict = {}
for each_bp in open(depth_file):
    each_bp_split = each_bp.strip().split('\t')
    seq_id = each_bp_split[0]
    seq_pos = each_bp_split[1]
    pos_depth = int(each_bp_split[2])

    if seq_id not in seq_pos_depth_dict:
        seq_pos_depth_dict[seq_id] = dict()
    seq_pos_depth_dict[seq_id][seq_pos] = pos_depth

    if seq_id not in seq_overall_depth_dict:
        seq_overall_depth_dict[seq_id] = 0
    seq_overall_depth_dict[seq_id] += pos_depth


average_global_coverage = seq_overall_depth_dict[sequence_id]/len(sequence_str)


window_start_pos = ignore_ends_len + 1
gc_content_to_depth_list_dict = {}
while window_start_pos <= (len(sequence_str) - (ignore_ends_len + window_len) + 1):
    window_end_pos = window_start_pos + window_len - 1

    # get sequence of current window
    if window_start_pos < (len(sequence_str) - window_len + 1):
        window_seq = sequence_str[(window_start_pos - 1):window_end_pos]
    else:
        window_seq = sequence_str[(window_start_pos - 1):]
    window_seq_upper = window_seq.upper()

    masked_base_num = window_seq_upper.count('N')
    if masked_base_num <= 4:
        window_pos_list = list(range(window_start_pos, (window_end_pos + 1)))
        window_pos_depth_list = [seq_pos_depth_dict[sequence_id][str(pos)] for pos in window_pos_list]
        window_mean_depth = sum(window_pos_depth_list)/len(window_pos_depth_list)
        window_gc_content = (window_seq_upper.count('G') + window_seq_upper.count('C'))*100/len(window_seq)
        #print('%s-%s\t%s\t%s' % (window_start_pos, window_end_pos, window_mean_depth, window_gc_content))

        # when the length of sliding window is 100 bp
        window_gc_content = int(window_gc_content)

        if window_gc_content not in gc_content_to_depth_list_dict:
            gc_content_to_depth_list_dict[window_gc_content] = [window_mean_depth]
        else:
            gc_content_to_depth_list_dict[window_gc_content].append(window_mean_depth)

    window_start_pos += 1


gc_bias_txt_handle = open(gc_bias_txt, 'w')
gc_bias_txt_handle.write('GC\tCoverage\tNormalized_coverage\tNumber\tProportion\n')
for each_gc_group in sorted([i for i in gc_content_to_depth_list_dict]):
    current_gc_group_depth_list = gc_content_to_depth_list_dict[each_gc_group]
    current_gc_group_proportion = len(current_gc_group_depth_list)*100/(len(sequence_str) - window_len)
    current_gc_group_mean_depth = sum(current_gc_group_depth_list)/len(current_gc_group_depth_list)
    current_gc_group_mean_depth_normalized = current_gc_group_mean_depth/average_global_coverage
    gc_bias_txt_handle.write('%s\t%s\t%s\t%s\t%s\n' % (each_gc_group,
                                                       float("{0:.2f}".format(current_gc_group_mean_depth)),
                                                       float("{0:.2f}".format(current_gc_group_mean_depth_normalized)),
                                                       len(current_gc_group_depth_list),
                                                       float("{0:.4f}".format(current_gc_group_proportion))))
gc_bias_txt_handle.close()


print('average_global_coverage: %s' % average_global_coverage)
print(len(sequence_str))

