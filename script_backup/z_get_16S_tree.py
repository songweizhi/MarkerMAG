import os
from Bio import SeqIO


combined_ffn_file   = '/Users/songweizhi/Desktop/combined.ffn'
seq_file_out        = '/Users/songweizhi/Desktop/combined_16S.ffn'
seq_file_aln        = '/Users/songweizhi/Desktop/combined_16S.aln'
seq_file_tree       = '/Users/songweizhi/Desktop/combined_16S.tree'
all_vs_all_16s      = '/Users/songweizhi/Desktop/combined_16S_all_vs_all.txt'


##################################################### get 16S tree #####################################################

seq_file_out_handle = open(seq_file_out, 'w')
for ffn_record in SeqIO.parse(combined_ffn_file, 'fasta'):
    if '16S ribosomal RNA' in ffn_record.description:
        seq_file_out_handle.write('>%s\n' % ffn_record.id)
        seq_file_out_handle.write('%s\n' % str(ffn_record.seq))
seq_file_out_handle.close()


mafft_cmd = 'mafft --quiet --retree 1 %s > %s' % (seq_file_out, seq_file_aln)
fasttree_cmd = '/Users/songweizhi/Softwares/FastTree/bin/FastTree -quiet %s > %s 2>/dev/null' % (seq_file_aln, seq_file_tree)
os.system(mafft_cmd)
os.system(fasttree_cmd)


################################################## get 16S divergence ##################################################

blastn_cmd = 'blastn -query %s -subject %s -out %s -outfmt 6' % (seq_file_out, seq_file_out, all_vs_all_16s)
os.system(blastn_cmd)


max_inter_genome = 0
min_intra_genome = 100
for line in open(all_vs_all_16s):

    line_split = line.strip().split('\t')
    query_genome    = line_split[0].split('_')[0]
    subject_genome  = line_split[1].split('_')[0]
    iden            = float(line_split[2])
    aln_len         = int(line_split[3])

    if aln_len >= 1300:
        if query_genome != subject_genome:
            if iden > max_inter_genome:
                max_inter_genome = iden
        if query_genome == subject_genome:
            if iden < min_intra_genome:
                min_intra_genome = iden

print('max_inter_genome: %s' % max_inter_genome)
print('min_intra_genome: %s' % min_intra_genome)

# max_inter_genome: 99.282
# min_intra_genome: 99.016

# MBARC26
# max_inter_genome: 98.765
# min_intra_genome: 98.04
