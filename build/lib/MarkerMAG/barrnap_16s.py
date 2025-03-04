import os
import glob
import shutil
import argparse
from Bio import SeqIO
import multiprocessing as mp

barrnap_16s_usage = '''
========= barrnap_16s example commands =========

MarkerMAG barrnap_16s -g MAG_1.fa
MarkerMAG barrnap_16s -g MAG_files -x fa -t 6

================================================
'''

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


def barrnap_16s(args):

    genome_folder   = args['g']
    genome_ext      = args['x']
    num_threads     = args['t']
    force_overwrite = args['force']

    genome_folder_no_path = genome_folder
    if '/' in genome_folder:
        genome_folder_no_path = genome_folder.split('/')[-1]

    output_prefix = ''
    if os.path.isdir(genome_folder) is True:
        genome_folder_no_path = genome_folder
        if '/' in genome_folder:
            genome_folder_no_path = genome_folder.split('/')[-1]
        output_prefix = genome_folder_no_path
    elif os.path.isfile(genome_folder) is True:
        gnm_file_path, gnm_file_basename, gnm_file_extension = sep_path_basename_ext(genome_folder)
        output_prefix = gnm_file_basename
    else:
        print('No folder/file detected, program exited!')
        exit()

    # define file name
    barrnap_16s_wd          = '%s_barrnap_16s_wd'               % output_prefix
    barrnap_cmd_txt         = '%s/barrnap_commadns.txt'         % barrnap_16s_wd
    barrnap_op_foler        = '%s/%s_barrnap_outputs'           % (barrnap_16s_wd, output_prefix)
    output_seq_16s_folder   = '%s/%s_16S_seq'                   % (barrnap_16s_wd, output_prefix)
    combined_16s_seqs       = '%s/%s_16S.fasta'                 % (barrnap_16s_wd, output_prefix)
    output_table            = '%s/%s_16S.txt'                   % (barrnap_16s_wd, output_prefix)
    output_table_stats      = '%s/%s_16S_stats.txt'             % (barrnap_16s_wd, output_prefix)

    # create folder
    if (os.path.isdir(barrnap_16s_wd) is True) and (force_overwrite is False):
        print('Output folder detected, program exited: %s' % barrnap_16s_wd)
        exit()
    else:
        force_create_folder(barrnap_16s_wd)
        os.mkdir(barrnap_op_foler)

    genome_file_list = []
    if os.path.isdir(genome_folder) is True:
        genome_file_re = '%s/*.%s' % (genome_folder, genome_ext)
        genome_file_list = [os.path.basename(file_name) for file_name in glob.glob(genome_file_re)]
        if len(genome_file_list) == 0:
            print('No query genome detected, program exited!')
            exit()
    elif os.path.isfile(genome_folder) is True:
        genome_file_list = [genome_folder_no_path]
    else:
        print('No query genome detected, program exited!')
        exit()

    # prepare commands for running barrnap
    argument_list_for_barrnap = []
    for genome in genome_file_list:
        gnm_id = '.'.join(genome.split('.')[:-1])
        pwd_genome = '%s/%s' % (genome_folder, genome)
        if os.path.isfile(genome_folder) is True:
            pwd_genome = genome_folder
            gnm_path, gnm_id, gnm_ext = sep_path_basename_ext(genome_folder)

        pwd_op_ffn = '%s/%s.ffn' % (barrnap_op_foler, gnm_id)
        pwd_op_gff = '%s/%s.gff' % (barrnap_op_foler, gnm_id)
        pwd_op_log = '%s/%s.log' % (barrnap_op_foler, gnm_id)
        barrnap_cmd = 'barrnap --quiet --outseq %s %s > %s 2> %s' % (pwd_op_ffn, pwd_genome, pwd_op_gff, pwd_op_log)
        argument_list_for_barrnap.append(barrnap_cmd)

    # write out barrnap commands
    with open(barrnap_cmd_txt, 'w') as barrnap_cmd_txt_handle:
        barrnap_cmd_txt_handle.write('\n'.join(argument_list_for_barrnap))

    # run barrnap with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(os.system, argument_list_for_barrnap)
    pool.close()
    pool.join()

    # create folder
    os.mkdir(output_seq_16s_folder)

    # parse gff files
    genome_failed_barrnap = set()
    output_table_handle = open(output_table, 'w')
    output_table_handle.write('Genome\t16S\tLength(bp)\tLocation\tStart\tEnd\n')
    output_table_stats_handle = open(output_table_stats, 'w')
    output_table_stats_handle.write('Genome\tcopy_number\n')
    for genome in genome_file_list:

        current_genome_id           = '.'.join(genome.split('.')[:-1])
        current_genome_gff          = '%s/%s.gff'       % (barrnap_op_foler, current_genome_id)
        current_genome_ffn          = '%s/%s.ffn'       % (barrnap_op_foler, current_genome_id)
        current_genome_16s_seq_file = '%s/%s_16S.ffn'   % (output_seq_16s_folder, current_genome_id)
        current_genome_fai_file     = '%s/%s.%s.fai'    % (genome_folder, current_genome_id, genome_ext)

        current_genome_16s_id_dict = {}
        seq_index_num = 1
        if os.path.isfile(current_genome_gff):
            for each_line in open(current_genome_gff):
                if (not each_line.startswith('#')) and ('16S_rRNA' in each_line):
                    each_line_split = each_line.strip().split('\t')
                    ctg_id = each_line_split[0]
                    start_pos = int(each_line_split[3])
                    end_pos = int(each_line_split[4])
                    length = end_pos - start_pos + 1
                    renamed_id_16s = '%s_16S_%s' % (genome, seq_index_num)
                    output_table_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (genome, renamed_id_16s, length, ctg_id, start_pos, end_pos))
                    current_genome_16s_id_dict[str(end_pos)] = renamed_id_16s
                    seq_index_num += 1
        else:
            genome_failed_barrnap.add(genome)

        # extract only 16S sequences
        current_genome_16s_seq_file_handle = open(current_genome_16s_seq_file, 'w')
        for seq_record in SeqIO.parse(current_genome_ffn, 'fasta'):
            seq_id = seq_record.id
            if seq_id.startswith('16S_rRNA'):
                seq_end_pos = seq_id.split(':')[3].split('(')[0].split('-')[-1]
                seq_id_new = current_genome_16s_id_dict[seq_end_pos]
                current_genome_16s_seq_file_handle.write('>%s\n' % seq_id_new)
                current_genome_16s_seq_file_handle.write('%s\n' % seq_record.seq)
        current_genome_16s_seq_file_handle.close()

        # remove index files
        if os.path.isfile(current_genome_fai_file):
            os.system('rm %s' % current_genome_fai_file)

        # write out copy number
        if len(current_genome_16s_id_dict) > 0:
            output_table_stats_handle.write('%s\t%s\n' % (genome, len(current_genome_16s_id_dict)))
        else:
            output_table_stats_handle.write('%s\t%s\n' % (genome, 0))

    output_table_handle.close()
    output_table_stats_handle.close()

    # combine 16S sequences
    os.system('cat %s/*.ffn > %s' % (output_seq_16s_folder, combined_16s_seqs))
    os.system('rm -r %s' % output_seq_16s_folder)

    # report failed genomes
    if len(genome_failed_barrnap) > 0:
        print('Failed to run barrnap on the following genomes: %s' % ','.join(genome_failed_barrnap))

    print('Output table exported to %s' % barrnap_16s_wd)
    print('Done!')


if __name__ == '__main__':

    barrnap_16s_parser = argparse.ArgumentParser(description='', usage=barrnap_16s_usage)
    barrnap_16s_parser.add_argument('-g',       required=True, type=str,               help='genome folder/file')
    barrnap_16s_parser.add_argument('-x',       required=False, default='fasta',       help='genome file extension, default: fasta')
    barrnap_16s_parser.add_argument('-t',       required=False, type=int, default=1,   help='number of threads, default: 1')
    barrnap_16s_parser.add_argument('-force',   required=False, action="store_true",   help='force overwrite existing results')
    args = vars(barrnap_16s_parser.parse_args())
    barrnap_16s(args)
