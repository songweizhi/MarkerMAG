#!/usr/bin/env python3

import os
import glob
import shutil
import argparse
from Bio import SeqIO
from datetime import datetime
from distutils.spawn import find_executable


matam_16s_usage = '''
=================================== matam_16s example commands ===================================

MarkerMAG matam_16s -p Test -r1 R1.fastq -r2 R2.fastq -t 12 -d path/to/SILVA_138_SSURef_NR99
MarkerMAG matam_16s -p Test -r16s 16S_reads.fastq -t 12 -d path/to/SILVA_138_SSURef_NR99

# please refer to here for the preparation of input files
https://github.com/songweizhi/MarkerMAG/blob/master/doc/README_matam_16s.md

==================================================================================================
'''


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def str_to_num_list(nums_str):

    subsample_pct_list = []
    for pct_value in [str(float(i)) for i in nums_str.split(',')]:
        if pct_value[-2:] == '.0':
            subsample_pct_list.append(int(float(pct_value)))
        else:
            subsample_pct_list.append(float(pct_value))

    return sorted(subsample_pct_list)


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


def single_line_seq(fasta_in, fasta_out):
    fasta_out_handle = open(fasta_out, 'w')
    for seq_record in SeqIO.parse(fasta_in, 'fasta'):
        fasta_out_handle.write('>%s\n' % seq_record.id)
        fasta_out_handle.write('%s\n' % seq_record.seq)
    fasta_out_handle.close()


def subsample_paired_and_singleton_reads(paired_r1, paired_r2, singleton_in, subsample_pct, paired_r1_out, paired_r2_out, singleton_out, usearch_exe):

    # check input file format
    input_seq_ext = os.path.splitext(paired_r1)[1]

    # define tmp files
    paired_r1_out_tmp = '%s.tmp' % paired_r1_out
    paired_r2_out_tmp = '%s.tmp' % paired_r2_out
    singleton_out_tmp = '%s.tmp' % singleton_out

    # prepare commands
    subsample_cmd_paired_r1 = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, paired_r1, paired_r1_out_tmp, subsample_pct)
    subsample_cmd_paired_r2 = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, paired_r2, paired_r2_out_tmp, subsample_pct)
    subsample_cmd_singleton = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, singleton_in, singleton_out_tmp, subsample_pct)
    if ('q' in input_seq_ext) or ('Q' in input_seq_ext):
        subsample_cmd_paired_r1 = '%s -fastx_subsample %s -fastqout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, paired_r1, paired_r1_out, subsample_pct)
        subsample_cmd_paired_r2 = '%s -fastx_subsample %s -fastqout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, paired_r2, paired_r2_out, subsample_pct)
        subsample_cmd_singleton = '%s -fastx_subsample %s -fastqout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, singleton_in, singleton_out, subsample_pct)

    # execute commands
    os.system(subsample_cmd_paired_r1)
    os.system(subsample_cmd_paired_r2)
    os.system(subsample_cmd_singleton)

    if ('q' not in input_seq_ext) and ('Q' not in input_seq_ext):

        # single line sequences
        single_line_seq(paired_r1_out_tmp, paired_r1_out)
        single_line_seq(paired_r2_out_tmp, paired_r2_out)
        single_line_seq(singleton_out_tmp, singleton_out)

        # remove tmp files
        os.remove(paired_r1_out_tmp)
        os.remove(paired_r2_out_tmp)
        os.remove(singleton_out_tmp)


def subsample_sortmerna_output(sortmerna_op, subsample_pct, sortmerna_op_subsampled, usearch_exe, seqtk_exe):

    # define file name
    fasta_in_path,  fasta_in_basename,  fasta_in_ext  = sep_path_basename_ext(sortmerna_op)
    fasta_out_path, fasta_out_basename, fasta_out_ext = sep_path_basename_ext(sortmerna_op_subsampled)

    fasta_in_paired_r1               = '%s/%s_%s_paired_r1%s'            % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)
    fasta_in_paired_r2               = '%s/%s_%s_paired_r2%s'            % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)
    fasta_in_singleton               = '%s/%s_%s_singleton%s'            % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)
    fasta_in_paired_r1_subsampled    = '%s/%s_%s_paired_r1_subsampled%s' % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)
    fasta_in_paired_r2_subsampled    = '%s/%s_%s_paired_r2_subsampled%s' % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)
    fasta_in_singleton_subsampled    = '%s/%s_%s_singleton_subsampled%s' % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)
    fasta_in_paired_subsampled       = '%s/%s_%s_paired_subsampled%s'    % (fasta_out_path, fasta_out_basename, fasta_in_basename, fasta_out_ext)

    # separate paired and singleton reads
    sep_paired_and_singleton_reads(sortmerna_op, fasta_in_paired_r1, fasta_in_paired_r2, fasta_in_singleton)

    # subsample paired and singleton reads
    subsample_paired_and_singleton_reads(fasta_in_paired_r1, fasta_in_paired_r2, fasta_in_singleton,
                                         subsample_pct,
                                         fasta_in_paired_r1_subsampled, fasta_in_paired_r2_subsampled, fasta_in_singleton_subsampled,
                                         usearch_exe)

    # combine subsampled paired reads
    combine_paired_fasta_cmd = '%s mergepe %s %s > %s' % (seqtk_exe, fasta_in_paired_r1_subsampled, fasta_in_paired_r2_subsampled, fasta_in_paired_subsampled)
    os.system(combine_paired_fasta_cmd)

    # combine subsampled paired and singleton reads
    combine_all_subsample_reads_cmd = 'cat %s %s > %s' % (fasta_in_paired_subsampled, fasta_in_singleton_subsampled, sortmerna_op_subsampled)
    os.system(combine_all_subsample_reads_cmd)

    # remove tmp files
    os.remove(fasta_in_paired_r1)
    os.remove(fasta_in_paired_r2)
    os.remove(fasta_in_singleton)
    os.remove(fasta_in_paired_r1_subsampled)
    os.remove(fasta_in_paired_r2_subsampled)
    os.remove(fasta_in_singleton_subsampled)
    os.remove(fasta_in_paired_subsampled)


def prefix_seq(seq_in, prefix, seq_out):

    seq_out_handle = open(seq_out, 'w')
    for seq_record in SeqIO.parse(seq_in, 'fasta'):
        seq_id_new = '%s_%s' % (prefix, seq_record.id)
        seq_out_handle.write('>%s\n' % seq_id_new)
        seq_out_handle.write('%s\n' % str(seq_record.seq))
    seq_out_handle.close()


def parse_uclust_output(uclust_output_table, cluster_to_member_file):

    cluster_id_set = set()
    cluster_to_seq_member_dict = {}
    for each_line in open(uclust_output_table):
        each_line_split = each_line.strip().split('\t')
        cluster_id = each_line_split[1]
        seq_id = each_line_split[8].split(' ')[0]
        cluster_id_set.add(int(cluster_id))

        if cluster_id not in cluster_to_seq_member_dict:
            cluster_to_seq_member_dict[cluster_id] = {seq_id}
        else:
            cluster_to_seq_member_dict[cluster_id].add(seq_id)

    # write out cluster sequence members
    cluster_to_member_file_handle = open(cluster_to_member_file, 'w')
    for each_cluster in sorted([i for i in cluster_id_set]):
        cluster_to_member_file_handle.write('Cluster_%s\t%s\n' % (each_cluster, ','.join(sorted([i for i in cluster_to_seq_member_dict[str(each_cluster)]]))))
    cluster_to_member_file_handle.close()


def matam_16s(args):

    ###################################################### file in #####################################################

    # file in
    output_prefix                   = args['p']
    reads_file_r1                   = args['r1']
    reads_file_r2                   = args['r2']
    input_16S_reads                 = args['r16s']
    subsample_pcts                  = args['pct']
    matam_ref                       = args['d']
    uclust_iden_cutoff              = args['i']
    num_threads                     = args['t']
    mem_mb                          = args['mem']
    force_overwrite                 = args['force']
    matam_assembly_script           = args['matam']
    usearch_exe                     = args['usearch']
    seqtk_exe                       = args['seqtk']
    keep_quiet                      = args['quiet']


    ################################################ check dependencies ################################################

    program_list = [usearch_exe, seqtk_exe, matam_assembly_script]

    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ', '.join(not_detected_programs))
        exit()


    ####################################################################################################################

    extract_16s_reads = True
    if (reads_file_r1 is None) and (reads_file_r2 is None) and (input_16S_reads is not None):
        if os.path.isfile(input_16S_reads) is True:
            extract_16s_reads = False
        else:
            print('%s not found, program exited!' % os.path.basename(input_16S_reads))
            exit()
    elif (reads_file_r1 is not None) and (reads_file_r2 is not None) and (input_16S_reads is None):
        if os.path.isfile(reads_file_r1) is False:
            print('%s not found, program exited!' % os.path.basename(reads_file_r1))
            exit()
        if os.path.isfile(reads_file_r2) is False:
            print('%s not found, program exited!' % os.path.basename(reads_file_r2))
            exit()
    else:
        print('please either provide paired-end reads (-r1 and -r2) or 16S reads (-r16s) extracted from the paired-end reads as input, program exited!')
        exit()

    input_reads_ext = ''
    if extract_16s_reads is True:
        input_reads_ext = os.path.splitext(reads_file_r1)[1]
    else:
        input_reads_ext = os.path.splitext(input_16S_reads)[1]


    ####################################################################################################################

    matam16s_wd                             = '%s_Matam16S_wd'                          % output_prefix
    log_file                                = '%s/%s_matam_16s.log'                     % (matam16s_wd, output_prefix)
    combined_r1_r2                          = '%s/%s_combined_R1_R2%s'                  % (matam16s_wd, output_prefix, input_reads_ext)
    extracted_16S_reads                     = '%s/%s_16S_reads%s'                       % (matam16s_wd, output_prefix, input_reads_ext)
    combined_all_depth_matam_assemblies     = '%s/%s_assembled_16S_unclustered.fasta'   % (matam16s_wd, output_prefix)
    uclust_output_uc                        = '%s/%s_assembled_16S_uclust_%s.uc'        % (matam16s_wd, output_prefix, uclust_iden_cutoff)
    uclust_output_txt                       = '%s/%s_assembled_16S_uclust_%s.txt'       % (matam16s_wd, output_prefix, uclust_iden_cutoff)
    uclust_output_fasta                     = '%s/%s_assembled_16S_uclust_%s.fasta'     % (matam16s_wd, output_prefix, uclust_iden_cutoff)

    # create folder
    if os.path.isdir(matam16s_wd) is False:
        os.mkdir(matam16s_wd)
    else:
        if force_overwrite is False:
            print('Output folder detected, program exited: %s' % matam16s_wd)
            exit()
        else:
            os.system('rm -r %s' % matam16s_wd)
            os.mkdir(matam16s_wd)

    if extract_16s_reads is True:

        # combine R1 and R2
        report_and_log(('Combining R1 and R2 with seqtk'), log_file, keep_quiet)
        combined_r1_r2_cmd = '%s mergepe %s %s > %s' % (seqtk_exe, reads_file_r1, reads_file_r2, combined_r1_r2)
        report_and_log(combined_r1_r2_cmd, log_file, keep_quiet)
        os.system(combined_r1_r2_cmd)

        # Extract 16S reads
        report_and_log(('Extracting 16S reads with Matam'), log_file, keep_quiet)
        matam_filter_cmd = '%s -i %s -o %s/%s_get_16S_reads_wd --cpu %s --max_memory %s -v --filter_only -d %s' % (matam_assembly_script, combined_r1_r2, matam16s_wd, output_prefix, num_threads, mem_mb, matam_ref)
        report_and_log(matam_filter_cmd, log_file, keep_quiet)
        os.system(matam_filter_cmd)

        # get file name of extracted 16S reads

        matam_16S_reads_re = '%s/%s_get_16S_reads_wd/workdir/*%s' % (matam16s_wd, output_prefix, input_reads_ext)
        matam_16S_reads_list = glob.glob(matam_16S_reads_re)
        if len(matam_16S_reads_list) == 0:
            print('16S reads not found, please check if 16S reads extraction step finished successfully.')
            exit()

        matam_16S_reads    = matam_16S_reads_list[0]
        os.system('mv %s %s' % (matam_16S_reads, extracted_16S_reads))
    else:
        extracted_16S_reads = input_16S_reads


    # subsample 16S reads and assemble
    renamed_matam_assembly_list = []
    for subsample_pct in str_to_num_list(subsample_pcts):

        report_and_log(('Subsample RNA reads at %s percent' % subsample_pct), log_file, keep_quiet)

        subsample_reads_file = '%s/%s_16S_reads_subset_%s%s'    % (matam16s_wd, output_prefix, subsample_pct, input_reads_ext)
        matam_output_folder  = '%s/%s_16S_reads_subset_%s_Matam_wd' % (matam16s_wd, output_prefix, subsample_pct)

        # subsample
        subsample_sortmerna_output(extracted_16S_reads, subsample_pct, subsample_reads_file, usearch_exe, seqtk_exe)

        # assemble with Matam
        report_and_log(('Assembling subsampled reads with Matam'), log_file, keep_quiet)

        matam_cmd = '%s -d %s -i %s --cpu %s --max_memory %s -v -o %s' % (matam_assembly_script, matam_ref, subsample_reads_file, num_threads, mem_mb, matam_output_folder)
        report_and_log(matam_cmd, log_file, keep_quiet)
        os.system(matam_cmd)

        matam_assemblies          = '%s/workdir/scaffolds.NR.min_500bp.fa'          % (matam_output_folder)
        matam_assemblies_prefixed = '%s/workdir/scaffolds.NR.min_500bp.prefixed.fa' % (matam_output_folder)
        seq_prefix                = '%s_subsample_%s'                               % (output_prefix, subsample_pct)

        if os.path.isfile(matam_assemblies) is True:
            report_and_log(('Adding prefix to Matam assemblies'), log_file, keep_quiet)
            prefix_seq(matam_assemblies, seq_prefix, matam_assemblies_prefixed)
            renamed_matam_assembly_list.append(matam_assemblies_prefixed)
        else:
            report_and_log(('No 16S rRNA gene sequence reconstructed at current depth!'), log_file, keep_quiet)

    # combine Matam outputs
    report_and_log(('Combine Matam assemblies at all depth'), log_file, keep_quiet)
    combine_cmd = 'cat %s > %s' % (' '.join(renamed_matam_assembly_list), combined_all_depth_matam_assemblies)
    report_and_log(combine_cmd, log_file, keep_quiet)
    os.system(combine_cmd)

    # dereplicate combined assemblies with Usearch
    report_and_log(('Cluster combined Matam assemblies with UCLUST at %s identity cutoff' % uclust_iden_cutoff), log_file, keep_quiet)
    uclust_cmd = '%s -cluster_fast %s -id %s -centroids %s -uc %s -sort length -quiet' % (usearch_exe, combined_all_depth_matam_assemblies, uclust_iden_cutoff, uclust_output_fasta, uclust_output_uc)
    report_and_log(uclust_cmd, log_file, keep_quiet)
    os.system(uclust_cmd)

    # get readable cluster results
    parse_uclust_output(uclust_output_uc, uclust_output_txt)

    # remove tmp file
    os.system('rm %s' % combined_r1_r2)

    # report
    report_and_log(('SortMeRNA identified 16S reads exported to %s' % os.path.basename(extracted_16S_reads)), log_file, keep_quiet)
    report_and_log(('Dereplicated Matam assemblies exported to %s'  % os.path.basename(uclust_output_fasta)), log_file, keep_quiet)
    report_and_log(('Done!'), log_file, keep_quiet)


if __name__ == '__main__':

    matam_16s_parser = argparse.ArgumentParser()
    matam_16s_parser.add_argument('-p',              required=True,  metavar='',                                            help='output prefix')
    matam_16s_parser.add_argument('-r1',             required=False, metavar='', default=None,                              help='forward reads')
    matam_16s_parser.add_argument('-r2',             required=False, metavar='', default=None,                              help='reverse reads')
    matam_16s_parser.add_argument('-r16s',           required=False, metavar='', default=None,                              help='extracted 16S reads, specify if you already have 16S reads identified from metagenomic dataset')
    matam_16s_parser.add_argument('-pct',            required=True,  metavar='', type=str, default='1,5,10,25,50,75,100',   help='subsample percentage, must be integer, between 1-100, deafault: 1,5,10,25,50,75,100')
    matam_16s_parser.add_argument('-d',              required=False, metavar='', type=str,                                  help='Matam ref db, controls "-d" in Matam')
    matam_16s_parser.add_argument('-i',              required=False, metavar='', type=float, default=0.999,                 help='cluster identity cutoff, do not provide value lower than 0.99, default: 0.999')
    matam_16s_parser.add_argument('-t',              required=False, metavar='', type=int, default=1,                       help='number of threads, default: 1')
    matam_16s_parser.add_argument('-mem',            required=False, metavar='', type=int, default=10240,                   help='maximum memory to be allocated to Matam (in MB), controls "--max_memory" in Matam, default: 10240')
    matam_16s_parser.add_argument('-force',          required=False, action="store_true",                                   help='force overwrite existing results')
    matam_16s_parser.add_argument('-quiet',          required=False, action="store_true",                                   help='not report progress')
    matam_16s_parser.add_argument('-matam',          required=False, metavar='', type=str, default='matam_assembly.py',     help='path to matam_assembly.py, default: matam_assembly.py')
    matam_16s_parser.add_argument('-seqtk',          required=False, metavar='', type=str, default='seqtk',                 help='path to seqtk executable file, default: seqtk')
    matam_16s_parser.add_argument('-usearch',        required=False, metavar='', type=str, default='usearch',               help='path to usearch executable file, default: usearch')
    args = vars(matam_16s_parser.parse_args())
    matam_16s(args)
