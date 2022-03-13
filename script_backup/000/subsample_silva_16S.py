from Bio import SeqIO

fasta_in  = '/Users/songweizhi/Desktop/subset.fa'
fasta_out = '/Users/songweizhi/Desktop/subsampled.fa'

fasta_in  = 'SILVA_138.1_SSURef_NR99_tax_silva_NR99.complete.fasta'
fasta_out = 'subsample_o.fasta'


fasta_out_handle = open(fasta_out, 'w')
wrote_taxon_set = set()
for each_seq in SeqIO.parse(fasta_in, 'fasta'):
    seq_des = each_seq.description
    taxon_str = ' '.join(seq_des.split(' ')[1:])

    if (taxon_str.startswith('Bacteria')) or (taxon_str.startswith('Archaea')):
        taxon_list = taxon_str.split(';')
        if len(taxon_list) >= 4:
            seq_d = taxon_list[0]
            seq_p = taxon_list[1]
            seq_c = taxon_list[2]
            seq_o = taxon_list[3]
            if seq_o not in wrote_taxon_set:
                fasta_out_handle.write('>%s\n' % seq_des)
                fasta_out_handle.write('%s\n' % str(each_seq.seq))
                wrote_taxon_set.add(seq_o)
                wrote_taxon_set.add(seq_c)
        else:
            if len(taxon_list) == 3:
                if taxon_list[2] not in wrote_taxon_set:
                    fasta_out_handle.write('>%s\n' % seq_des)
                    fasta_out_handle.write('%s\n' % str(each_seq.seq))
                    wrote_taxon_set.add(taxon_list[2])
            else:
             print(taxon_str)

fasta_out_handle.close()
print(len(wrote_taxon_set))
