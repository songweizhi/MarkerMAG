import os
import glob
from Bio import SeqIO


genome_folder       = '/Users/songweizhi/Desktop/selected_genomes_renamed_no_plasmid'
genome_abundance    = '/Users/songweizhi/Desktop/abundance_3.txt'
genome_ext          = 'fna'
sequencing_depth    = 25
read_len            = 150
insert_size         = 200


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
    Reads_simulator_cmd = 'BioSAK Reads_simulator -r %s -n %s -l %s -i %s -split' % (genome_to_sim, paires_to_sim_dict_by_size_by_abund_norm[genome_to_sim], read_len, insert_size)
    print(Reads_simulator_cmd)



abun_1 = '''
BioSAK Reads_simulator -r s1.fna -n 275510 -l 150 -i 200 -split
BioSAK Reads_simulator -r s2.fna -n 276512 -l 150 -i 200 -split
BioSAK Reads_simulator -r s3.fna -n 789988 -l 150 -i 200 -split
BioSAK Reads_simulator -r s4.fna -n 667154 -l 150 -i 200 -split
BioSAK Reads_simulator -r s5.fna -n 779904 -l 150 -i 200 -split
BioSAK Reads_simulator -r g1.fna -n 376005 -l 150 -i 200 -split
BioSAK Reads_simulator -r g2.fna -n 142261 -l 150 -i 200 -split
BioSAK Reads_simulator -r g3.fna -n 827388 -l 150 -i 200 -split
BioSAK Reads_simulator -r g4.fna -n 75433 -l 150 -i 200 -split
BioSAK Reads_simulator -r g5.fna -n 586525 -l 150 -i 200 -split
BioSAK Reads_simulator -r f1.fna -n 172421 -l 150 -i 200 -split
BioSAK Reads_simulator -r f2.fna -n 490120 -l 150 -i 200 -split
BioSAK Reads_simulator -r f3.fna -n 263348 -l 150 -i 200 -split
BioSAK Reads_simulator -r f4.fna -n 291691 -l 150 -i 200 -split
BioSAK Reads_simulator -r f5.fna -n 210338 -l 150 -i 200 -split
BioSAK Reads_simulator -r o1.fna -n 97120 -l 150 -i 200 -split
BioSAK Reads_simulator -r o2.fna -n 337343 -l 150 -i 200 -split
BioSAK Reads_simulator -r o3.fna -n 312352 -l 150 -i 200 -split
BioSAK Reads_simulator -r o4.fna -n 440042 -l 150 -i 200 -split
BioSAK Reads_simulator -r o5.fna -n 340160 -l 150 -i 200 -split
BioSAK Reads_simulator -r c1.fna -n 49574 -l 150 -i 200 -split
BioSAK Reads_simulator -r c2.fna -n 348565 -l 150 -i 200 -split
BioSAK Reads_simulator -r c3.fna -n 312528 -l 150 -i 200 -split
BioSAK Reads_simulator -r p1.fna -n 4480 -l 150 -i 200 -split
BioSAK Reads_simulator -r p2.fna -n 133356 -l 150 -i 200 -split
BioSAK Reads_simulator -r p3.fna -n 195077 -l 150 -i 200 -split
BioSAK Reads_simulator -r p4.fna -n 379490 -l 150 -i 200 -split
BioSAK Reads_simulator -r p5.fna -n 780155 -l 150 -i 200 -split

'''


abun_2 ='''
BioSAK Reads_simulator -r s1.fna -n 897956 -l 150 -i 200 -split
BioSAK Reads_simulator -r s2.fna -n 738443 -l 150 -i 200 -split
BioSAK Reads_simulator -r s3.fna -n 857948 -l 150 -i 200 -split
BioSAK Reads_simulator -r s4.fna -n 629486 -l 150 -i 200 -split
BioSAK Reads_simulator -r s5.fna -n 607319 -l 150 -i 200 -split
BioSAK Reads_simulator -r g1.fna -n 985937 -l 150 -i 200 -split
BioSAK Reads_simulator -r g2.fna -n 28969 -l 150 -i 200 -split
BioSAK Reads_simulator -r g3.fna -n 424690 -l 150 -i 200 -split
BioSAK Reads_simulator -r g4.fna -n 27429 -l 150 -i 200 -split
BioSAK Reads_simulator -r g5.fna -n 444313 -l 150 -i 200 -split
BioSAK Reads_simulator -r f1.fna -n 129326 -l 150 -i 200 -split
BioSAK Reads_simulator -r f2.fna -n 388281 -l 150 -i 200 -split
BioSAK Reads_simulator -r f3.fna -n 579996 -l 150 -i 200 -split
BioSAK Reads_simulator -r f4.fna -n 9712 -l 150 -i 200 -split
BioSAK Reads_simulator -r f5.fna -n 348945 -l 150 -i 200 -split
BioSAK Reads_simulator -r o1.fna -n 219486 -l 150 -i 200 -split
BioSAK Reads_simulator -r o2.fna -n 510318 -l 150 -i 200 -split
BioSAK Reads_simulator -r o3.fna -n 241223 -l 150 -i 200 -split
BioSAK Reads_simulator -r o4.fna -n 238579 -l 150 -i 200 -split
BioSAK Reads_simulator -r o5.fna -n 59394 -l 150 -i 200 -split
BioSAK Reads_simulator -r c1.fna -n 118448 -l 150 -i 200 -split
BioSAK Reads_simulator -r c2.fna -n 95731 -l 150 -i 200 -split
BioSAK Reads_simulator -r c3.fna -n 7835 -l 150 -i 200 -split
BioSAK Reads_simulator -r p1.fna -n 25624 -l 150 -i 200 -split
BioSAK Reads_simulator -r p2.fna -n 98856 -l 150 -i 200 -split
BioSAK Reads_simulator -r p3.fna -n 64123 -l 150 -i 200 -split
BioSAK Reads_simulator -r p4.fna -n 1092338 -l 150 -i 200 -split
BioSAK Reads_simulator -r p5.fna -n 84140 -l 150 -i 200 -split

'''


abun_3 ='''
BioSAK Reads_simulator -r s1.fna -n 801741 -l 150 -i 200 -split
BioSAK Reads_simulator -r s2.fna -n 94795 -l 150 -i 200 -split
BioSAK Reads_simulator -r s3.fna -n 559157 -l 150 -i 200 -split
BioSAK Reads_simulator -r s4.fna -n 584592 -l 150 -i 200 -split
BioSAK Reads_simulator -r s5.fna -n 40707 -l 150 -i 200 -split
BioSAK Reads_simulator -r g1.fna -n 119535 -l 150 -i 200 -split
BioSAK Reads_simulator -r g2.fna -n 375285 -l 150 -i 200 -split
BioSAK Reads_simulator -r g3.fna -n 96564 -l 150 -i 200 -split
BioSAK Reads_simulator -r g4.fna -n 349913 -l 150 -i 200 -split
BioSAK Reads_simulator -r g5.fna -n 883424 -l 150 -i 200 -split
BioSAK Reads_simulator -r f1.fna -n 144196 -l 150 -i 200 -split
BioSAK Reads_simulator -r f2.fna -n 503550 -l 150 -i 200 -split
BioSAK Reads_simulator -r f3.fna -n 123655 -l 150 -i 200 -split
BioSAK Reads_simulator -r f4.fna -n 150978 -l 150 -i 200 -split
BioSAK Reads_simulator -r f5.fna -n 20798 -l 150 -i 200 -split
BioSAK Reads_simulator -r o1.fna -n 46860 -l 150 -i 200 -split
BioSAK Reads_simulator -r o2.fna -n 560325 -l 150 -i 200 -split
BioSAK Reads_simulator -r o3.fna -n 528863 -l 150 -i 200 -split
BioSAK Reads_simulator -r o4.fna -n 546313 -l 150 -i 200 -split
BioSAK Reads_simulator -r o5.fna -n 270647 -l 150 -i 200 -split
BioSAK Reads_simulator -r c1.fna -n 362995 -l 150 -i 200 -split
BioSAK Reads_simulator -r c2.fna -n 420119 -l 150 -i 200 -split
BioSAK Reads_simulator -r c3.fna -n 540118 -l 150 -i 200 -split
BioSAK Reads_simulator -r p1.fna -n 167556 -l 150 -i 200 -split
BioSAK Reads_simulator -r p2.fna -n 178746 -l 150 -i 200 -split
BioSAK Reads_simulator -r p3.fna -n 35521 -l 150 -i 200 -split
BioSAK Reads_simulator -r p4.fna -n 808233 -l 150 -i 200 -split
BioSAK Reads_simulator -r p5.fna -n 639659 -l 150 -i 200 -split

'''



