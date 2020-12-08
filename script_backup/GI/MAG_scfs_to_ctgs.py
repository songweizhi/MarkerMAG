import os
import glob
from Bio import SeqIO


genome_folder       = '/Users/songweizhi/Desktop/BH_ER_050417_refined_bins'
genome_ext          = 'fasta'
genome_folder_out   = '/Users/songweizhi/Desktop/BH_ER_050417_refined_bins_noNs'


file_re = '%s/*.%s' % (genome_folder, genome_ext)
file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]


for genome in file_list:
    pwd_genome = '%s/%s' % (genome_folder, genome)
    pwd_genome_out = '%s/%s' % (genome_folder_out, genome)
    pwd_genome_out_tmp = '%s/%s.tmp' % (genome_folder_out, genome)
    pwd_genome_out_tmp_handle = open(pwd_genome_out_tmp, 'w')
    for seq_record in SeqIO.parse(pwd_genome, 'fasta'):
        pwd_genome_out_tmp_handle.write('>%s\n' % seq_record.id)
        pwd_genome_out_tmp_handle.write('%s\n' % seq_record.seq)
    pwd_genome_out_tmp_handle.close()


    split_cmd = 'perl /Users/songweizhi/Desktop/split.scaffolds.to.contigs.pl -i %s -o %s -m 2000' % (pwd_genome_out_tmp, pwd_genome_out)
    os.system(split_cmd)

    os.system('rm %s' % pwd_genome_out_tmp)

