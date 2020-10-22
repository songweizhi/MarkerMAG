#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO


uclust_16S_usage = '''
========================== uclust_16S example commands ==========================

module load usearch/10.0.240
python3 uclust_16S.py -in combined_16S.fasta -i 1 -out dereplicated_16S.fasta

=================================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def parse_uclust_output(seq_file_in, uclust_output_table, seq_file_out, cluster_to_member_file):

    seq_len_dict = {}
    for seq_record in SeqIO.parse(seq_file_in, 'fasta'):
        seq_len_dict[seq_record.id] = len(seq_record.seq)

    cluster_id_set = set()
    cluster_to_rep_seq_dict = {}
    cluster_to_seq_member_dict = {}
    for each_line in open(uclust_output_table):
        each_line_split = each_line.strip().split('\t')
        cluster_id = each_line_split[1]
        seq_id = each_line_split[8].split(' ')[0]
        cluster_id_set.add(int(cluster_id))

        if cluster_id not in cluster_to_rep_seq_dict:
            cluster_to_rep_seq_dict[cluster_id] = seq_id
            cluster_to_seq_member_dict[cluster_id] = {seq_id}
        else:
            if seq_len_dict[seq_id] > seq_len_dict[cluster_to_rep_seq_dict[cluster_id]]:
                cluster_to_rep_seq_dict[cluster_id] = seq_id
            cluster_to_seq_member_dict[cluster_id].add(seq_id)

    # write out cluster sequence members
    cluster_to_member_file_handle = open(cluster_to_member_file, 'w')
    for each_cluster in sorted([i for i in cluster_id_set]):
        cluster_to_member_file_handle.write('Cluster_%s\t%s\n' % (each_cluster, ','.join(sorted([i for i in cluster_to_seq_member_dict[str(each_cluster)]]))))
    cluster_to_member_file_handle.close()

    # write out the longest sequence in each cluster
    rep_seq_id_set = set()
    for i in cluster_to_rep_seq_dict:
        rep_seq_id_set.add(cluster_to_rep_seq_dict[i])

    rep_seq_file_handle = open(seq_file_out, 'w')
    for each_seq in SeqIO.parse(seq_file_in, 'fasta'):
        if each_seq.id in rep_seq_id_set:
            rep_seq_file_handle.write('>%s\n' % each_seq.id)
            rep_seq_file_handle.write('%s\n' % str(each_seq.seq))
    rep_seq_file_handle.close()


parser = argparse.ArgumentParser()

parser.add_argument('-in',  required=True,                         help='fasta file')
parser.add_argument('-i',   required=False, type=float, default=1, help='Identity cutoff (0-1), default: 1')
parser.add_argument('-out', required=True,                         help='file out')

args = vars(parser.parse_args())
fasta_in    = args['in']
iden_cutoff = args['i']
fasta_out   = args['out']

# define tmp file name
fasta_out_path, fasta_out_basename, fasta_out_extension = sep_path_basename_ext(fasta_out)
uclust_output_table             = '%s/%s.uc'                     % (fasta_out_path, fasta_out_basename)
default_centroids_to_be_ignored = '%s/%s_default_centroids.%s'   % (fasta_out_path, fasta_out_basename, fasta_out_extension)
cluster_to_member_file          = '%s/%s.uc.reorganised.txt'     % (fasta_out_path, fasta_out_basename)

# run cluster_fast
uclust_cmd = 'usearch -cluster_fast %s -id %s -centroids %s -uc %s -quiet' % (fasta_in, iden_cutoff, default_centroids_to_be_ignored, uclust_output_table)
os.system(uclust_cmd)

# extract the longest sequence from each cluster
parse_uclust_output(fasta_in, uclust_output_table, fasta_out, cluster_to_member_file)

# report
print('Dereplicated sequences exported to: %s' % fasta_out)
