import os
import glob
from Bio import SeqIO


ref_16s_seq = '/Users/songweizhi/Desktop/MarkerMAG_wd/2_MBARC26/combined_16S.ffn'
linkage_file_dir = '/Users/songweizhi/Desktop/linkage_files'

linkage_file_re   = '%s/*_identified_linkages_genome_level.txt' % linkage_file_dir
linkage_file_list = [os.path.basename(file_name) for file_name in glob.glob(linkage_file_re)]


gnm_16s_num_dict = {}
for each_seq in SeqIO.parse(ref_16s_seq, 'fasta'):
    seq_id = each_seq.id
    gnm_id = seq_id.split('_')[0]
    if gnm_id not in gnm_16s_num_dict:
        gnm_16s_num_dict[gnm_id] = 1
    else:
        gnm_16s_num_dict[gnm_id] += 1


ref_gnm_list = sorted([i for i in gnm_16s_num_dict])
ref_16s_num_list = []
for each_ref in ref_gnm_list:
    each_ref_linked_16s_num = gnm_16s_num_dict.get(each_ref, 0)
    ref_16s_num_list.append(each_ref_linked_16s_num)
ref_16s_num_list_str = [str(i) for i in ref_16s_num_list]

print('Reference_id\t%s' % '\t'.join(ref_gnm_list))
print('Reference\t%s' % '\t'.join(ref_16s_num_list_str))


for each_linkage_file in sorted(linkage_file_list):

    pwd_linkage_file = '%s/%s' % (linkage_file_dir, each_linkage_file)
    linkage_file_prefix = each_linkage_file.split('_identified_linkages_genome_level')[0]

    linked_gnm_16s_num_dict = {}
    for each_match in open(pwd_linkage_file):
        if not each_match.startswith('MarkerGene,GenomicSeq,Number'):
            match_split = each_match.strip().split('\t')
            id_16s = match_split[0]
            id_16s_gnm = id_16s.split('_')[0]
            id_gnm = match_split[1]
            if id_16s_gnm == id_gnm:
                if id_gnm not in linked_gnm_16s_num_dict:
                    linked_gnm_16s_num_dict[id_gnm] = 1
                else:
                    linked_gnm_16s_num_dict[id_gnm] += 1

    linked_16s_num_list = []
    for each_ref in ref_gnm_list:
        each_ref_linked_16s_num = linked_gnm_16s_num_dict.get(each_ref, 0)
        linked_16s_num_list.append(each_ref_linked_16s_num)
    linked_16s_num_list_str = [str(i) for i in linked_16s_num_list]
    print('%s\t%s' % (linkage_file_prefix, '\t'.join(linked_16s_num_list_str)))
