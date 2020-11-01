#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def export_dna_record(gene_seq, gene_id, gene_description, output_handle):
    seq_object = Seq(gene_seq, IUPAC.unambiguous_dna)
    seq_record = SeqRecord(seq_object)
    seq_record.id = gene_id
    seq_record.description = gene_description
    SeqIO.write(seq_record, output_handle, 'fasta')


# file in
rm_16s_parser = argparse.ArgumentParser(description='', usage='')
rm_16s_parser.add_argument('-gbk', required=True, type=str, help='gbk file')
args = vars(rm_16s_parser.parse_args())
gbk_file = args['gbk']


# define file name
gbk_file_path, gbk_file_basename, gbk_file_extension = sep_path_basename_ext(gbk_file)
fasta_out = '%s/%s_no_16S.fasta' % (gbk_file_path, gbk_file_basename)


pos_list_16s = [1]
genome_seq_str = ''
for record in SeqIO.parse(gbk_file, 'genbank'):
    genome_seq_str = str(record.seq)
    for gene in record.features:
        if 'locus_tag' in gene.qualifiers:
            if gene.qualifiers['product'] == ['16S ribosomal RNA']:
                start_pos = gene.location.start
                end_pos = gene.location.end
                pos_list_16s.append(int(start_pos))
                pos_list_16s.append(int(end_pos))

pos_list_16s.append(len(genome_seq_str))

pos_list_16s = sorted(pos_list_16s)

fasta_out_handle = open(fasta_out, 'w')
n = 1
while n <= len(pos_list_16s):

    between_region_index = n//2 + 1

    if pos_list_16s[n-1] == 1:
        between_region_start = 1
    else:
        between_region_start = pos_list_16s[n-1] + 1

    if pos_list_16s[n] == len(genome_seq_str):
        between_region_end = len(genome_seq_str)
    else:
        between_region_end = pos_list_16s[n] - 1

    between_region_seq     = genome_seq_str[between_region_start - 1: between_region_end]
    between_region_seq_id  = '%s_%s'   % (gbk_file_basename, between_region_index)
    between_region_seq_des = '%s-%sbp' % (between_region_start, between_region_end)

    export_dna_record(between_region_seq, between_region_seq_id, between_region_seq_des, fasta_out_handle)

    n += 2

fasta_out_handle.close()


'''

cd /Users/songweizhi/Desktop
python3 ~/PycharmProjects/MarkerMAG/script_backup/rm_16S.py -gbk c3.gbk

module load python/3.7.3
cd /srv/scratch/z5039045/MarkerMAG_wd/genome_selection_3/gbk_files
python3 /srv/scratch/z5039045/MarkerMAG_wd/rm_16S.py c1.gbk
python3 /srv/scratch/z5039045/MarkerMAG_wd/rm_16S.py c2.gbk
python3 /srv/scratch/z5039045/MarkerMAG_wd/rm_16S.py c3.gbk

'''