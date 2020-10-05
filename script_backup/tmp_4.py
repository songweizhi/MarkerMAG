from Bio import SeqIO

file_in  = '/Users/songweizhi/Desktop/combined.ffn'
file_out = '/Users/songweizhi/Desktop/combined_16S.ffn'

file_out_handle = open(file_out, 'w')
for gene_record in SeqIO.parse(file_in, 'fasta'):

    if '16S ribosomal RNA' in gene_record.description:
        file_out_handle.write('>%s\n'% gene_record.id)
        file_out_handle.write('%s\n' % gene_record.seq)
file_out_handle.close()
