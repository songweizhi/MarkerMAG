import os
from Bio import SeqIO


# file in
wd                      = '/Users/songweizhi/Desktop/get_fake_bins_wd'
combined_ref            = '%s/combined_refs.fna'                         % wd
ctg_depth_file          = '%s/ctg_depth.txt'                             % wd
file_in_16s             = '%s/16S_gene_list.txt'                         % wd


ref_id_set = set()
genome_to_seq_dict = {}
seq_len_dict = {}
total_ref_seq = 0
for seq_record in SeqIO.parse(combined_ref, 'fasta'):
    seq_genome = seq_record.id.split('_')[0]
    seq_len_dict[seq_record.id] = len(seq_record.seq)
    if seq_genome not in genome_to_seq_dict:
        genome_to_seq_dict[seq_genome] = [seq_record.id]
    else:
        genome_to_seq_dict[seq_genome].append(seq_record.id)


seq_depth_dict = {}
for each_depth in open(ctg_depth_file):
    if not each_depth.startswith('contigName'):
        each_depth_split = each_depth.strip().split('\t')
        seq_depth_dict[each_depth_split[0]] = float(each_depth_split[1])


genome_to_mean_depth_dict = {}
for genome in genome_to_seq_dict:
    genome_total_length = 0
    genome_total_depth = 0
    for each_seq in genome_to_seq_dict[genome]:
        genome_total_length += seq_len_dict[each_seq]
        genome_total_depth += ((seq_len_dict[each_seq])*(seq_depth_dict[each_seq]))
    average_depth = genome_total_depth/genome_total_length
    average_depth = float("{0:.2f}".format(average_depth))
    genome_to_mean_depth_dict[genome] = average_depth


for each_16s in open(file_in_16s):
    each_16s_genome = each_16s.strip().split('_')[0]
    print('%s\t%s' % (each_16s.strip(), genome_to_mean_depth_dict[each_16s_genome]))

# cd /Users/songweizhi/Desktop/get_fake_bins_wd
# BioSAK iTOL -SimpleBar -lv 16S_gene_to_depth.txt -scale 0-500-1000-1500 -lt Depth -out SimpleBar_depth.txt

print()
for gnm in genome_to_mean_depth_dict:
    print('%s\t%s' % (gnm, genome_to_mean_depth_dict[gnm]))

