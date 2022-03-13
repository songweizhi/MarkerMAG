import os
import glob
from Bio import SeqIO


ref_gnm_folder  = '/Users/songweizhi/Desktop/reference_genomes_renamed'
ref_gnm_ext     = 'fna'
ref_depth_txt   = '/Users/songweizhi/Desktop/666/ref_depth.txt'
read_len        = 130

ref_depth_dict = {}
for each_ref_16s in open(ref_depth_txt):
    each_ref_16s_split = each_ref_16s.strip().split('\t')
    gnm_id = each_ref_16s_split[0].split('_')[0]
    gnm_depth = float(each_ref_16s_split[1])
    ref_depth_dict[gnm_id] = gnm_depth
print(ref_depth_dict)


gnm_file_re = '%s/*.%s' % (ref_gnm_folder, ref_gnm_ext)
gnm_file_list = [os.path.basename(file_name) for file_name in glob.glob(gnm_file_re)]

print(gnm_file_list)

gnm_total_len_dict = {}
for each_gnm in gnm_file_list:

    gnm_id = '.'.join(each_gnm.split('.')[:-1])
    pwd_gnm_file = '%s/%s' % (ref_gnm_folder, each_gnm)

    current_gnm_total_len = 0
    for each_seq in SeqIO.parse(pwd_gnm_file, 'fasta'):
        current_gnm_total_len += len(each_seq.seq)

    gnm_depth = ref_depth_dict[gnm_id]

    total_bp_to_simulate = current_gnm_total_len * gnm_depth

    total_num_to_simulate = round(total_bp_to_simulate/(read_len * 2))

    #print('%s\t%s\t%s\t%s' % (gnm_id, current_gnm_total_len, gnm_depth, total_num_to_simulate))



    simulate_cmd = 'BioSAK Reads_simulator -p %s -r %s/%s -l 130 -i 200 -split -n %s' % (gnm_id, '/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/reference_genomes_renamed', each_gnm, total_num_to_simulate)
    print(simulate_cmd)





















