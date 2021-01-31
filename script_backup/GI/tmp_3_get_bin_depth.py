import os
import glob
from Bio import SeqIO


ctg_len_depth_file  = '/Users/songweizhi/Desktop/GI_depth.txt'
genome_folder       = '/Users/songweizhi/Desktop/GI_refined_bins'
genome_ext          = 'fasta'
gnm_mean_depth_file = '/Users/songweizhi/Desktop/GI_bin_mean_depth.txt'


# ctg_len_depth_file  = '/Users/songweizhi/Desktop/BH_ER_050417_depth.txt'
# genome_folder       = '/Users/songweizhi/Desktop/BH_ER_050417_refined_bins'
# genome_ext          = 'fasta'
# gnm_mean_depth_file = '/Users/songweizhi/Desktop/BH_ER_050417_refined_bins_mean_depth.txt'


# get ctg_len_dict and ctg_depth_dict
ctg_len_dict = {}
ctg_depth_dict = {}
for each_line in open(ctg_len_depth_file):
    if not each_line.startswith('contigName	contigLen	totalAvgDepth'):
        each_line_split = each_line.strip().split('\t')
        ctg_id = each_line_split[0]
        ctg_len = int(each_line_split[1])
        ctg_depth = float(each_line_split[2])
        ctg_len_dict[ctg_id] = ctg_len
        ctg_depth_dict[ctg_id] = ctg_depth


file_re = '%s/*.%s' % (genome_folder, genome_ext)
file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]

gnm_mean_depth_dict = {}
gnm_mean_depth_file_handle = open(gnm_mean_depth_file, 'w')
for genome in file_list:
    pwd_genome = '%s/%s' % (genome_folder, genome)
    gnm_total_len = 0
    gnm_total_depth = 0
    for each_ctg in SeqIO.parse(pwd_genome, 'fasta'):
        gnm_total_len += len(each_ctg.seq)
        gnm_total_depth += (len(each_ctg.seq)*ctg_depth_dict[each_ctg.id])
    gnm_mean_depth = gnm_total_depth/gnm_total_len
    gnm_mean_depth = float("{0:.2f}".format(gnm_mean_depth))
    gnm_mean_depth_file_handle.write('%s\t%s\t%s\n' % (genome, gnm_total_len, gnm_mean_depth))
    gnm_mean_depth_dict[genome] = gnm_mean_depth
gnm_mean_depth_file_handle.close()


#
# linkaged_bin_list = []
# for each_link in open('/Users/songweizhi/Desktop/CAMI2_HMP_combined_linkages.txt'):
#     linkaged_bin_list.append(each_link.strip().split('\t')[1])
# print(linkaged_bin_list)
#
# print(gnm_mean_depth_dict)
#
# for each_gnm in gnm_mean_depth_dict:
#
#     gnm_no_ext= '.'.join(each_gnm.split('.')[:-1])
#
#     gnm_linked = 'no'
#     if gnm_no_ext in linkaged_bin_list:
#         gnm_linked = 'yes'
#
#     print('%s\t%s\t%s' % (gnm_no_ext, gnm_mean_depth_dict[each_gnm], gnm_linked))
#
