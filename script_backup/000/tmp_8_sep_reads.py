import os
from Bio import SeqIO


def sep_paired_and_singleton_reads(seq_in, seq_out_r1, seq_out_r2, seq_out_singleton):

    input_seq_ext = os.path.splitext(seq_in)[1]
    input_seq_fmt = 'fasta'
    if ('q' in input_seq_ext) or ('Q' in input_seq_ext):
        input_seq_fmt = 'fastq'

    fmt_to_provide_for_write = 'fasta-2line'
    if input_seq_fmt == 'fastq':
        fmt_to_provide_for_write = 'fastq'

    reads_pair_dict = {}
    for read_record in SeqIO.parse(seq_in, input_seq_fmt):
        read_id_base = '.'.join(read_record.id.split('.')[:-1])
        read_strand = read_record.id.split('.')[-1]
        if read_id_base not in reads_pair_dict:
            reads_pair_dict[read_id_base] = {read_strand}
        else:
            reads_pair_dict[read_id_base].add(read_strand)

    read_list_paired = set()
    read_list_singleton = set()
    for read_base in reads_pair_dict:
        if len(reads_pair_dict[read_base]) == 1:
            read_list_singleton.add(read_base)
        if len(reads_pair_dict[read_base]) == 2:
            read_list_paired.add(read_base)

    fasta_out_r1_handle = open(seq_out_r1, 'w')
    fasta_out_r2_handle = open(seq_out_r2, 'w')
    fasta_out_singleton_handle = open(seq_out_singleton, 'w')
    for read_record in SeqIO.parse(seq_in, input_seq_fmt):
        read_id_base = '.'.join(read_record.id.split('.')[:-1])
        read_strand = read_record.id.split('.')[-1]
        if read_id_base in read_list_singleton:
            SeqIO.write(read_record, fasta_out_singleton_handle, fmt_to_provide_for_write)
        if read_id_base in read_list_paired:
            if read_strand == '1':
                SeqIO.write(read_record, fasta_out_r1_handle, fmt_to_provide_for_write)
            if read_strand == '2':
                SeqIO.write(read_record, fasta_out_r2_handle, fmt_to_provide_for_write)
    fasta_out_r1_handle.close()
    fasta_out_r2_handle.close()
    fasta_out_singleton_handle.close()


seq_in              = '/Users/songweizhi/Desktop/111/CAMI_Oral_16S_reads.fasta'
seq_out_r1          = '/Users/songweizhi/Desktop/111/CAMI_Oral_16S_reads_R1.fasta'
seq_out_r2          = '/Users/songweizhi/Desktop/111/CAMI_Oral_16S_reads_R2.fasta'
seq_out_singleton   = '/Users/songweizhi/Desktop/111/CAMI_Oral_16S_reads_UP.fasta'

sep_paired_and_singleton_reads(seq_in, seq_out_r1, seq_out_r2, seq_out_singleton)



