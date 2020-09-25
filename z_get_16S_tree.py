import os
from Bio import SeqIO

combined_ffn_file   = '/Users/songweizhi/Desktop/combined.ffn'
seq_file_out        = '/Users/songweizhi/Desktop/combined_16S.ffn'
seq_file_aln        = '/Users/songweizhi/Desktop/combined_16S.aln'
seq_file_tree       = '/Users/songweizhi/Desktop/combined_16S.tree'

seq_file_out_handle = open(seq_file_out, 'w')
for ffn_record in SeqIO.parse(combined_ffn_file, 'fasta'):
    if '16S ribosomal RNA' in ffn_record.description:
        seq_file_out_handle.write('>%s\n' % ffn_record.id)
        print(ffn_record.id)
        seq_file_out_handle.write('%s\n' % str(ffn_record.seq))
seq_file_out_handle.close()

mafft_cmd = 'mafft --quiet --retree 1 %s > %s' % (seq_file_out, seq_file_aln)
fasttree_cmd = '/Users/songweizhi/Softwares/FastTree/bin/FastTree -quiet %s > %s 2>/dev/null' % (seq_file_aln, seq_file_tree)
os.system(mafft_cmd)
os.system(fasttree_cmd)

