import os
from Bio import SeqIO, AlignIO


seq_in_16s      = '/Users/songweizhi/Desktop/00/GI_128_16S_0.999.fasta'
min_cov_id = 1
sliding_window = 25


seq_in_16s_aln          = '/Users/songweizhi/Desktop/00/GI_128_16S_0.999.fasta.aln'
seq_in_16s_filtered     = '/Users/songweizhi/Desktop/00/GI_128_16S_0.999_filtered.fasta'
alignment_file_out      = '/Users/songweizhi/Desktop/00/GI_128_16S_0.999_filtered.fasta.aln'


seq_in_16s_mafft_cmd = 'mafft --quiet --retree 1 %s > %s' % (seq_in_16s, seq_in_16s_aln)
os.system(seq_in_16s_mafft_cmd)

alignment_in = AlignIO.read(seq_in_16s_aln, "fasta")
total_col_num = alignment_in.get_alignment_length()
seq_id_list = [seq.id for seq in alignment_in]
seqs_to_remove = set()
col_index = 0
while col_index < (total_col_num - sliding_window + 1):
    current_column = [str(seg.seq) for seg in alignment_in[:, col_index:(col_index + sliding_window)]]
    non_dash_number = len(current_column) - current_column.count(sliding_window*'-')
    non_gap_percent = float("{0:.3f}".format(non_dash_number*100/len(seq_id_list)))
    if non_gap_percent < min_cov_id:
        row_index = 0
        for each_base in current_column:
            if each_base != (sliding_window*'-'):
                seqs_to_remove.add(seq_id_list[row_index])
            row_index += 1
    col_index += 1


seq_in_16s_filtered_handle = open(seq_in_16s_filtered, 'w')
for read_seq in SeqIO.parse(seq_in_16s, 'fasta'):
    if read_seq.id not in seqs_to_remove:
        seq_in_16s_filtered_handle.write('>%s\n' % read_seq.id)
        seq_in_16s_filtered_handle.write('%s\n' % str(read_seq.seq))
seq_in_16s_filtered_handle.close()


filtered_16s_mafft_cmd = 'mafft --quiet --retree 1 %s > %s' % (seq_in_16s_filtered, alignment_file_out)
os.system(filtered_16s_mafft_cmd)

print(seqs_to_remove)
print(len(seqs_to_remove))
