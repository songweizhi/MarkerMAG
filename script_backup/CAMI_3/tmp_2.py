from Bio import SeqIO

op_handle = open('combined_hc_refs_uniq.fa', 'w')
wrote_seq_set = set()
for each_seq in SeqIO.parse('combined_hc_refs.fa', 'fasta'):
    seq_id = each_seq.id
    if seq_id not in wrote_seq_set:
        SeqIO.write(each_seq, op_handle, 'fasta')
        wrote_seq_set.add(seq_id)
op_handle.close()

