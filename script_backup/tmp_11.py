import os
import glob
import shutil
import argparse
import pandas as pd
from Bio import SeqIO
from time import sleep
import multiprocessing as mp
from datetime import datetime
from distutils.spawn import find_executable


def split_list(list_in, subset_num):

    list_in_formatted = [i for i in list_in]

    # get the number of element per subset
    file_num_per_folder = round(len(list_in_formatted) / subset_num)

    n = 1
    lol_out = []
    while n <= subset_num:

        if n < subset_num:
            current_subset_elements = {i for i in list_in_formatted[(file_num_per_folder * (n - 1)):(file_num_per_folder * n)]}
            lol_out.append(current_subset_elements)
        else:
            current_subset_elements = {i for i in list_in_formatted[(file_num_per_folder * (n - 1)):]}
            lol_out.append(current_subset_elements)

        n += 1

    return lol_out


def extract_reads_worker(argument_list):

    reads_file_in    = argument_list[0]
    reads_to_extract = argument_list[1]
    reads_file_out   = argument_list[2]

    reads_file_out_handle = open(reads_file_out, 'w')
    for read_record in SeqIO.parse(reads_file_in, 'fasta'):
        if read_record.id in reads_to_extract:
            reads_file_out_handle.write('>%s\n' % read_record.id)
            reads_file_out_handle.write('%s\n' % read_record.seq)
    reads_file_out_handle.close()


def extracted_reads_with_multiprocessing(reads_r1, reads_r2, r1_to_extract, r2_to_extract, output_folder, num_threads):
    solely_perfectly_mapped_reads_r1_splitted = split_list(r1_to_extract, num_threads // 2)
    solely_perfectly_mapped_reads_r2_splitted = split_list(r2_to_extract, num_threads // 2)

    argument_list_for_extract_reads_worker = []
    extract_reads_file_index_r1 = 1
    for reads_subset_r1 in solely_perfectly_mapped_reads_r1_splitted:
        current_output_file = '%s/extract_r1_subset_%s.fasta' % (output_folder, extract_reads_file_index_r1)
        argument_list_for_extract_reads_worker.append([reads_r1, reads_subset_r1, current_output_file])
        extract_reads_file_index_r1 += 1

    extract_reads_file_index_r2 = 1
    for reads_subset_r2 in solely_perfectly_mapped_reads_r2_splitted:
        current_output_file = '%s/extract_r2_subset_%s.fasta' % (output_folder, extract_reads_file_index_r2)
        argument_list_for_extract_reads_worker.append([reads_r2, reads_subset_r2, current_output_file])
        extract_reads_file_index_r2 += 1

    # extract reads with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(extract_reads_worker, argument_list_for_extract_reads_worker)
    pool.close()
    pool.join()


fq_1 = '/Users/songweizhi/Desktop/rename_fq/o1_R1_Q30_P.fastq'
fq_2 = '/Users/songweizhi/Desktop/rename_fq/o1_R2_Q30_P.fastq'

for read_record in SeqIO.parse(fq_1, 'fastq'):
    print(read_record.id)

    read_record.id = '000'

    print(read_record.id)


def fq2fa(fq_file, fa_file):
    pass
