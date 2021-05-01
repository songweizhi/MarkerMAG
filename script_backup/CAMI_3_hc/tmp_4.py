import os
from Bio import SeqIO


def usearch_uc(uc_file):

    cluster_to_seq_dict = {}
    seq_to_cluster_dict = {}
    for each_line in open(uc_file):
        each_line_split = each_line.strip().split('\t')
        cluster_id = 'cluster_%s' % each_line_split[1]
        seq_id = each_line_split[8].split(' ')[0]

        if cluster_id not in cluster_to_seq_dict:
            cluster_to_seq_dict[cluster_id] = {seq_id}
        else:
            cluster_to_seq_dict[cluster_id].add(seq_id)

        seq_to_cluster_dict[seq_id] = cluster_id

    return cluster_to_seq_dict, seq_to_cluster_dict


uc_file  = '/Users/songweizhi/Desktop/usearch.uc'
seq_file = '/Users/songweizhi/Desktop/cami_hc_SILVA138_id99_assembled_16S_uclust_0.999.fasta'
cluster_folder = '/Users/songweizhi/Desktop/cluster_folder'


if os.path.isdir(cluster_folder) is True:
    os.system('rm -r %s' % cluster_folder)
os.system('mkdir %s' % cluster_folder)


cluster_to_seq_dict, seq_to_cluster_dict = usearch_uc(uc_file)


for each_seq in SeqIO.parse(seq_file, 'fasta'):
    seq_id = each_seq.id
    seq_cluster = seq_to_cluster_dict[seq_id]
    cluster_seq_file = '%s/%s.fa' % (cluster_folder, seq_cluster)
    with open(cluster_seq_file, 'a') as cluster_seq_file_handle:
        cluster_seq_file_handle.write('>%s\n' % seq_id)
        cluster_seq_file_handle.write('%s\n' % str(each_seq.seq))


print(cluster_to_seq_dict)
print(seq_to_cluster_dict)
print(seq_to_cluster_dict['cami_hc_SILVA138_id99_75_subsample_75_1629'])

for each_cluster in cluster_to_seq_dict:
    cluster_seqs = cluster_to_seq_dict[each_cluster]
    pwd_cluster_seq_file = '%s/%s.fa' % (cluster_folder, each_cluster)
    pwd_cluster_aln_file = '%s/%s.fa.aln' % (cluster_folder, each_cluster)
    pwd_cluster_aln_file_1line = '%s/%s.fa.oneline.aln' % (cluster_folder, each_cluster)

    if len(cluster_seqs) > 1:
        mafft_cmd = 'mafft --quiet --retree 1 %s > %s' % (pwd_cluster_seq_file, pwd_cluster_aln_file)
        os.system(mafft_cmd)
        one_line_aln_cmd = 'BioSAK OneLineAln -in %s -out %s -upper' % (pwd_cluster_aln_file, pwd_cluster_aln_file_1line)
        os.system(one_line_aln_cmd)


def polish_16s():
    pass


