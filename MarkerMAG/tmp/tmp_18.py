from Bio import SeqIO


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

            # print(each_line_split)
            print('%s(%sbp)\t%s\t%s\t%s' % (ctg_id, ctg_len, len_16s, left_gap, right_gap))
            # print()

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


combined_gnm_file       = '/Users/songweizhi/Desktop/GI_refined_bins_combined.fa'
combined_barrnap_gff    = '/Users/songweizhi/Desktop/combined.gff'


ctg_len_dict = {}
for each_ctg in SeqIO.parse(combined_gnm_file, 'fasta'):
    ctg_len_dict[each_ctg.id] = len(each_ctg.seq)


ctg_ignore_region_dict = get_regions_to_ignore_from_barrnap_output(combined_barrnap_gff, ctg_len_dict)
print(ctg_ignore_region_dict)



