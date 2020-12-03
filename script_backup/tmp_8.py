from Bio import SeqIO


fasta_in  = '/Users/songweizhi/Desktop/CB_ER_070716.fasta'
fasta_out = '/Users/songweizhi/Desktop/CB_ER_070716_renamed.fasta'

fasta_in  = 'combined_16s_reads.fa'
fasta_out = 'combined_16s_reads_renamed.fa'


fasta_out_handle = open(fasta_out, 'w')
for read_record in SeqIO.parse(fasta_in, 'fasta'):
    read_record_description = read_record.description
    read_record_description_split = read_record_description.split(' ')
    read_record_new_id = ''
    if read_record_description_split[1][0] == '1':
        read_record_new_id = '%s.1' % read_record.id
    if read_record_description_split[1][0] == '2':
        read_record_new_id = '%s.2' % read_record.id
    fasta_out_handle.write('>%s\n' % read_record_new_id)
    fasta_out_handle.write('%s\n' % read_record.seq)
fasta_out_handle.close()
