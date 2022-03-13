#!/usr/bin/env python3

import os
import argparse


subsample_reads_usage = '''
================== subsample_reads example commands ==================

MarkerMAG subsample_reads -r1 R1.fasta -r2 R2.fasta -ratio 0.1
MarkerMAG subsample_reads -r1 R1.fasta -r2 R2.fasta -ratio 0.1,0.3,0.5
MarkerMAG subsample_reads -r1 R1.fastq -r2 R2.fastq -ratio 0.1,0.3,0.5

======================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def subsample_reads(args):


    r1_file         = args['r1']
    r2_file         = args['r2']
    subsample_ratio = args['ratio']
    usearch_exe     = args['usearch']

    r1_file_path, r1_file_basename, r1_file_extension = sep_path_basename_ext(r1_file)
    r2_file_path, r2_file_basename, r2_file_extension = sep_path_basename_ext(r2_file)

    subsample_step_list = [float(i) for i in subsample_ratio.split(',')]

    for subsample_step in subsample_step_list:

        r1_file_out = '%s_%s%s' % (r1_file_basename, subsample_step, r1_file_extension)
        r2_file_out = '%s_%s%s' % (r2_file_basename, subsample_step, r2_file_extension)

        sunsample_cmd = ''
        if r1_file_extension != r2_file_extension:
            print('Different reads formats detected, program exited!')
            exit()
        elif ('q' in r1_file_extension) and ('q' in r2_file_extension):
            sunsample_cmd = '%s -fastx_subsample %s -reverse %s -fastqout %s -output2 %s -sample_pct %s' % (usearch_exe, r1_file, r2_file, r1_file_out, r2_file_out, subsample_step*100)
            os.system(sunsample_cmd)
        else:
            sunsample_cmd_r1 = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, r1_file, r1_file_out, subsample_step * 100)
            sunsample_cmd_r2 = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, r2_file, r2_file_out, subsample_step * 100)
            os.system(sunsample_cmd_r1)
            os.system(sunsample_cmd_r2)


if __name__ == '__main__':

    subsample_reads_parser = argparse.ArgumentParser(description='', usage=subsample_reads_usage)

    subsample_reads_parser.add_argument('-r1',      required=True, type=str,            help='forward reads')
    subsample_reads_parser.add_argument('-r2',      required=True, type=str,            help='reverse reads')
    subsample_reads_parser.add_argument('-ratio',   required=True, type=str,            help='subsample ratio, 0-1')
    subsample_reads_parser.add_argument('-usearch', required=False, default='usearch',  help='path to usearch executable file, default: usearch')

    args = vars(subsample_reads_parser.parse_args())
    subsample_reads(args)
