import os
import argparse
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def polish_16s(file_in, file_out_ffn):

    file_out_path, file_out_base, file_out_ext = sep_path_basename_ext(file_out_ffn)

    file_out_gff     = '%s/%s.gff'    % (file_out_path, file_out_base)
    file_out_ffn_tmp = '%s/%s_tmp.%s' % (file_out_path, file_out_base, file_out_ext)

    barrnap_cmd = 'barrnap --quiet -o %s %s > %s' % (file_out_ffn_tmp, file_in, file_out_gff)
    os.system(barrnap_cmd)

    file_out_ffn_handle = open(file_out_ffn, 'w')
    for each_16s in SeqIO.parse(file_out_ffn_tmp, 'fasta'):
        seq_id = each_16s.id
        if seq_id.startswith('16S_rRNA::'):
            seq_id_polished = seq_id[10:].split(':')[0]
            file_out_ffn_handle.write('>%s\n' % seq_id_polished)
            file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
    file_out_ffn_handle.close()

    os.system('rm %s' % file_out_ffn_tmp)
    os.system('rm %s.fai' % file_in)



file_in = 'trim_16s.fasta'
file_out_ffn = 'trimed_16s.fasta'
polish_16s(file_in, file_out_ffn)
