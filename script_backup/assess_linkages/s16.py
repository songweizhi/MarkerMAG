from Bio import SeqIO

seq_16s_MBARC26 = '/Users/songweizhi/Desktop/get_depth/MBARC26_SILVA138_polished.QC.fasta'
seq_16s_GI      = '/Users/songweizhi/Desktop/get_depth/GI_128_16S_0.999.QC.fasta'
seq_16s_Oral    = '/Users/songweizhi/Desktop/get_depth/CAMI_Oral_138_16S_0.999.polished_min1200.QC.fa'


seq_len_list_MBARC26 = []
for each_seq in SeqIO.parse(seq_16s_MBARC26, 'fasta'):
    seq_len_list_MBARC26.append(str(len(each_seq.seq)))
print(','.join(seq_len_list_MBARC26))

seq_len_list_GI = []
for each_seq in SeqIO.parse(seq_16s_GI, 'fasta'):
    seq_len_list_GI.append(str(len(each_seq.seq)))
print(','.join(seq_len_list_GI))


seq_len_list_Oral = []
for each_seq in SeqIO.parse(seq_16s_Oral, 'fasta'):
    seq_len_list_Oral.append(str(len(each_seq.seq)))
print(','.join(seq_len_list_Oral))


