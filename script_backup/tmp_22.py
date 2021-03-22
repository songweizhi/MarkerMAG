from Bio import SeqIO


ref_in      = '/Users/songweizhi/Desktop/000/BH_ER_050417_refined_bins_combined.fa'
ref_subset  = '/Users/songweizhi/Desktop/000/BH_ER_050417_refined_bins_combined_subset.fa'
end_seq_len = 1800


# get ref seqs subset
ref_subset_len_dict = {}
ref_subset_handle = open(ref_subset, 'w')
for ref_seq in SeqIO.parse(ref_in, 'fasta'):
    ref_seq_id = ref_seq.id
    ref_seq_len = len(ref_seq.seq)
    if ref_seq_len < end_seq_len * 3:
        ref_subset_handle.write('>%s\n' % ref_seq_id)
        ref_subset_handle.write('%s\n' % ref_seq.seq)
        ref_subset_len_dict[ref_seq_id] = ref_seq_len
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
ref_subset_handle.close()


