#!/usr/bin/env python3

# Copyright (C) 2020, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com or t.thomas@unsw.edu.au

# MarkerMAG is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MarkerMAG is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import os
import argparse
from Bio import SeqIO
import multiprocessing as mp


rename_reads_usage = '''
=================== rename_reads example commands ===================

MarkerMAG rename_reads -r1 R1.fasta -r2 R2.fasta -p Soil -t 2
MarkerMAG rename_reads -r1 R1.fastq -r2 R2.fastq -p Soil -fq -t 2

Note: The order of paired reads in the two files must be the same.

=====================================================================
'''


def sep_path_basename_ext(file_in):

    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def rename_reads_worker(arg_list):

    file_in         = arg_list[0]
    file_out        = arg_list[1]
    output_prefix   = arg_list[2]
    file_format     = arg_list[3]
    read_strand     = arg_list[4]

    file_out_handle = open(file_out, 'w')
    read_index = 1
    for read_record in SeqIO.parse(file_in, file_format):
        if file_format == 'fasta':
            file_out_handle.write('>%s_%s.%s\n' % (output_prefix, read_index, read_strand))
            file_out_handle.write('%s\n' % read_record.seq)
        if file_format == 'fastq':
            read_record.id = '%s_%s.%s' % (output_prefix, read_index, read_strand)
            read_record.description = ''
            SeqIO.write(read_record, file_out_handle, file_format)
        read_index += 1
    file_out_handle.close()


def rename_reads(args):

    reads_r1        = args['r1']
    reads_r2        = args['r2']
    output_prefix   = args['p']
    fastq_format    = args['fq']
    num_threads     = args['t']

    # define name of output files
    format_to_parse = 'fasta'
    output_file_r1 = '%s_R1.fasta' % output_prefix
    output_file_r2 = '%s_R2.fasta' % output_prefix
    if fastq_format is True:
        format_to_parse = 'fastq'
        output_file_r1 = '%s_R1.fastq' % output_prefix
        output_file_r2 = '%s_R2.fastq' % output_prefix

    if os.path.isfile(output_file_r1) is True:
        print('%s detected, please remove your existing file or specify a different prefix' % output_file_r1)
        exit()
    if os.path.isfile(output_file_r2) is True:
        print('%s detected, please remove your existing file or specify a different prefix' % output_file_r2)
        exit()

    pool = mp.Pool(processes=num_threads)
    pool.map(rename_reads_worker, [[reads_r1, output_file_r1, output_prefix, format_to_parse, 1], [reads_r2, output_file_r2, output_prefix, format_to_parse, 2]])
    pool.close()
    pool.join()


if __name__ == '__main__':

    rename_reads_parser = argparse.ArgumentParser(description='', usage=rename_reads_usage)
    rename_reads_parser.add_argument('-r1',  required=True, type=str,               help='forward reads, fasta format')
    rename_reads_parser.add_argument('-r2',  required=True, type=str,               help='reverse reads, fasta format')
    rename_reads_parser.add_argument('-p',   required=True, type=str,               help='prefix of output file and read id')
    rename_reads_parser.add_argument('-fq',  required=False, action="store_true",   help='rename reads with quality score')
    rename_reads_parser.add_argument('-t',   required=False, type=int, default=1,   help='number of threads, default: 1')
    args = vars(rename_reads_parser.parse_args())
    rename_reads(args)
