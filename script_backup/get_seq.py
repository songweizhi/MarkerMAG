from Bio import SeqIO

for seq_record in SeqIO.parse('combined_CF_R2.fasta', 'fasta'):

    if seq_record.id == 'CF3_54230053.2':
        print(seq_record.id)
        print(str(seq_record.seq))



