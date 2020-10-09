import os
import glob
import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument('-g', required=True, help='genome folder')
parser.add_argument('-x', required=True, help='genome extension')
parser.add_argument('-a', required=True, help='abundance file')
parser.add_argument('-d', required=True, type=int, help='sequencing depth')
parser.add_argument('-l', required=True, type=int, help='read length')
parser.add_argument('-i', required=True, type=int, help='insert size')
parser.add_argument('-GemSIM', required=False, action="store_true", help='print GemSIM commands')

args = vars(parser.parse_args())

genome_folder       = args['g']
genome_ext          = args['x']
genome_abundance    = args['a']
sequencing_depth    = args['d']
read_len            = args['l']
insert_size         = args['i']
run_GemSIM          = args['GemSIM']

# genome_folder       = '/Users/songweizhi/Desktop/selected_genomes_renamed_no_plasmid'
# genome_abundance    = '/Users/songweizhi/Desktop/abundance_3.txt'
# genome_ext          = 'fna'
# sequencing_depth    = 25
# read_len            = 150
# insert_size         = 200

'''
python3 /srv/scratch/z5039045/MarkerMAG_wd/z_simulate_reads.py -g /srv/scratch/z5039045/MarkerMAG_wd/genome_selection/selected_genomes_renamed_no_plasmid -x fna -a abundance_even.txt -d 7 -l 300 -i 200



'''


genome_re = '%s/*.%s' % (genome_folder, genome_ext)
genome_list = [os.path.basename(file_name) for file_name in glob.glob(genome_re)]


genome_size_dict = {}
for genome_in in genome_list:
    pwd_genome_file = '%s/%s' % (genome_folder, genome_in)
    current_genome_len = 0
    for seq_record in SeqIO.parse(pwd_genome_file, 'fasta'):
        current_genome_len += len(seq_record.seq)
    genome_size_dict[genome_in] = current_genome_len


genome_abundance_dict = {}
for genome_abund in open(genome_abundance):
    genome_abund_split = genome_abund.strip().split('\t')
    genome_abundance_dict[genome_abund_split[0]] = float(genome_abund_split[1])


genome_size_total = 0
for genome in genome_size_dict:
    genome_size_total += genome_size_dict[genome]


total_reads_len_to_simulate = genome_size_total*sequencing_depth
total_reads_pair_to_simulate = total_reads_len_to_simulate//(2*read_len)


paires_to_sim_dict_by_size_by_abund = {}
paired_reads_to_sim_by_size_by_abund_total = 0
for genome in genome_abundance_dict:
    current_genome_size = genome_size_dict[genome]
    current_genome_abun = genome_abundance_dict[genome]
    paired_reads_to_sim_by_size_by_abund = total_reads_pair_to_simulate * (current_genome_size/genome_size_total) * current_genome_abun
    paires_to_sim_dict_by_size_by_abund[genome] = paired_reads_to_sim_by_size_by_abund
    paired_reads_to_sim_by_size_by_abund_total += paired_reads_to_sim_by_size_by_abund


paires_to_sim_dict_by_size_by_abund_norm = {}
for paires_to_sim in paires_to_sim_dict_by_size_by_abund:
    paires_to_sim_norm = paires_to_sim_dict_by_size_by_abund[paires_to_sim] * (total_reads_pair_to_simulate/paired_reads_to_sim_by_size_by_abund_total)
    paires_to_sim_dict_by_size_by_abund_norm[paires_to_sim] = round(paires_to_sim_norm)


for genome_to_sim in paires_to_sim_dict_by_size_by_abund_norm:

    genome_to_sim_basename = '.'.join(genome_to_sim.split('.')[:-1])
    if run_GemSIM is True:
        Reads_simulator_cmd = 'python /srv/scratch/z5039045/Softwares/GemSIM_v1.6/GemReads.py -r %s -n %s -m ../../phix/phix_2020_i5_p.gzip -o %s_%sx -q 33 -u d -l 151 -p -c\nmv rd.log %s_rd.log' % (genome_to_sim, paires_to_sim_dict_by_size_by_abund_norm[genome_to_sim], genome_to_sim_basename, sequencing_depth, genome_to_sim_basename)

    else:
        Reads_simulator_cmd = 'BioSAK Reads_simulator -r %s/%s -n %s -l %s -i %s -split' % (genome_folder, genome_to_sim, paires_to_sim_dict_by_size_by_abund_norm[genome_to_sim], read_len, insert_size)

    print(Reads_simulator_cmd)


print('total_reads_pair_to_simulate %sx: %s' % (sequencing_depth, total_reads_pair_to_simulate))
