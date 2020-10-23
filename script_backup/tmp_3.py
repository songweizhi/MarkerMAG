from Bio import SeqIO

file_1         = '/Users/songweizhi/Desktop/combined_all_depth_assemblies_iden99.5_uniq_default_centroids..fasta'
file_2         = '/Users/songweizhi/Desktop/combined_all_depth_assemblies_iden99.5_uniq.fasta'


file_1_id_list = []
for seq in SeqIO.parse(file_1, 'fasta'):
    file_1_id_list.append(seq.id)

file_2_id_list = []
for seq in SeqIO.parse(file_2, 'fasta'):
    file_2_id_list.append(seq.id)


c = set(file_1_id_list).intersection(file_2_id_list)



print(len(file_1_id_list))
print(len(file_2_id_list))
print(len(c))
