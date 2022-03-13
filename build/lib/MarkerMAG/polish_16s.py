#!/usr/bin/env python3

import os
import glob
import shutil
import argparse
from Bio import SeqIO
from datetime import datetime
from distutils.spawn import find_executable


polish_16s_usage = '''
===================== polish_16s example commands =====================

MarkerMAG polish_16s -in Matam_output.fa -out Matam_output_polished.fa

=======================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def polish_16s(args):

    file_in = args['in']
    file_out_ffn = args['out']

    file_out_path, file_out_base, file_out_ext = sep_path_basename_ext(file_out_ffn)

    barrnap_stdout   = '%s/%s.log'    % (file_out_path, file_out_base)
    file_out_gff     = '%s/%s.gff'    % (file_out_path, file_out_base)
    file_out_ffn_tmp = '%s/%s_tmp%s' % (file_out_path, file_out_base, file_out_ext)

    barrnap_cmd = 'barrnap --quiet -o %s %s 2> %s > %s' % (file_out_ffn_tmp, file_in, barrnap_stdout, file_out_gff)
    os.system(barrnap_cmd)

    wrote_id = []
    file_out_ffn_handle = open(file_out_ffn, 'w')
    for each_16s in SeqIO.parse(file_out_ffn_tmp, 'fasta'):
        seq_id = each_16s.id
        if seq_id.startswith('16S_rRNA::'):
            seq_id_polished = seq_id[10:].split(':')[0]

            if seq_id_polished not in wrote_id:
                file_out_ffn_handle.write('>%s\n' % seq_id_polished)
                file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
                wrote_id.append(seq_id_polished)
            else:
                file_out_ffn_handle.write('>%s_%s\n' % (seq_id_polished, (wrote_id.count(seq_id_polished) + 1)))
                file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
                wrote_id.append(seq_id_polished)

    file_out_ffn_handle.close()

    os.system('rm %s' % file_out_ffn_tmp)
    os.system('rm %s.fai' % file_in)


if __name__ == '__main__':

    matam_16s_parser = argparse.ArgumentParser()

    matam_16s_parser.add_argument('-in',  required=True, help='input 16S sequences')
    matam_16s_parser.add_argument('-out', required=True, help='output sequences')

    args = vars(matam_16s_parser.parse_args())
    polish_16s(args)
