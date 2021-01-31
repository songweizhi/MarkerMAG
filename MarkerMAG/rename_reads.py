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


rename_reads_usage = '''
============================== rename_reads example commands ==============================

MarkerMAG rename_reads -r1 R1.fasta -r2 R2.fasta -p Soil
MarkerMAG rename_reads -r1 R1.fastq -r2 R2.fastq -p Soil -fq

Note: The order of paired reads in the two files must be the same.

python3 ~/PycharmProjects/MarkerMAG/MarkerMAG/rename_reads.py -r1 o1_R1_Q30_P.fastq  -r2 o1_R2_Q30_P.fastq  -p Test -fq



===========================================================================================
'''


def sep_path_basename_ext(file_in):

    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def rename_reads(args):

    reads_r1        = args['r1']
    reads_r2        = args['r2']
    output_prefix   = args['p']
    fastq_format    = args['fq']


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

    output_file_r1_handle = open(output_file_r1, 'w')
    r1_index = 1
    for r1_record in SeqIO.parse(reads_r1, format_to_parse):
        if format_to_parse == 'fasta':
            output_file_r1_handle.write('>%s_%s.1\n' % (output_prefix, r1_index))
            output_file_r1_handle.write('%s\n' % r1_record.seq)
        if format_to_parse == 'fastq':
            r1_record.id = '%s_%s.1' % (output_prefix, r1_index)
            r1_record.description = ''
            SeqIO.write(r1_record, output_file_r1_handle, format_to_parse)
        r1_index += 1
    output_file_r1_handle.close()

    output_file_r2_handle = open(output_file_r2, 'w')
    r2_index = 1
    for r2_record in SeqIO.parse(reads_r2, format_to_parse):
        if format_to_parse == 'fasta':
            output_file_r2_handle.write('>%s_%s.2\n' % (output_prefix, r2_index))
            output_file_r2_handle.write('%s\n' % r2_record.seq)
        if format_to_parse == 'fastq':
            r2_record.id = '%s_%s.2' % (output_prefix, r2_index)
            r2_record.description = ''
            SeqIO.write(r2_record, output_file_r2_handle, format_to_parse)
        r2_index += 1
    output_file_r2_handle.close()


if __name__ == '__main__':

    rename_reads_parser = argparse.ArgumentParser(description='', usage=rename_reads_usage)
    rename_reads_parser.add_argument('-r1',  required=True, type=str,               help='forward reads, fasta format')
    rename_reads_parser.add_argument('-r2',  required=True, type=str,               help='reverse reads, fasta format')
    rename_reads_parser.add_argument('-p',   required=True, type=str,               help='prefix of output file and read id')
    rename_reads_parser.add_argument('-fq',  required=False, action="store_true",   help='rename reads with quality score')
    args = vars(rename_reads_parser.parse_args())

    rename_reads(args)

