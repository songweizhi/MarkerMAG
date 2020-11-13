#!/usr/bin/env python3

import os
import shutil
import argparse
from Bio import SeqIO
from datetime import datetime
from distutils.spawn import find_executable


matam_16s_usage = '''
=================================== matam_16s example commands ===================================

MarkerMAG matam_16s -p Test -in combined_R1_R2.fasta -pct 0.1,0.5,1,5,10,25,50,75 -i 0.995 -t 12 -force -ref path/to/SILVA_128_SSURef_NR95 -matam_assembly path/to/matam_assembly.py -sortmerna /path/to/sortmerna

==================================================================================================
'''


def report_and_log(message_for_report, log_file, keep_quiet):

    time_format = '[%Y-%m-%d %H:%M:%S]'
    with open(log_file, 'a') as log_handle:
        log_handle.write('%s %s\n' % ((datetime.now().strftime(time_format)), message_for_report))

    if keep_quiet is False:
        print('%s %s' % ((datetime.now().strftime(time_format)), message_for_report))


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)
    os.mkdir(folder_to_create)


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

    return subsample_pct_list


def sep_paired_and_singleton_reads(fasta_in, fasta_out_r1, fasta_out_r2, fasta_out_singleton):
    reads_pair_dict = {}
    for read_record in SeqIO.parse(fasta_in, 'fasta'):
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

    fasta_out_r1_handle = open(fasta_out_r1, 'w')
    fasta_out_r2_handle = open(fasta_out_r2, 'w')
    fasta_out_singleton_handle = open(fasta_out_singleton, 'w')

    for read_record in SeqIO.parse(fasta_in, 'fasta'):

        read_id_base = '.'.join(read_record.id.split('.')[:-1])
        read_strand = read_record.id.split('.')[-1]

        if read_id_base in read_list_singleton:
            fasta_out_singleton_handle.write('>%s\n' % read_record.id)
            fasta_out_singleton_handle.write('%s\n' % str(read_record.seq))

        if read_id_base in read_list_paired:

            if read_strand == '1':
                fasta_out_r1_handle.write('>%s\n' % read_record.id)
                fasta_out_r1_handle.write('%s\n' % str(read_record.seq))

            if read_strand == '2':
                fasta_out_r2_handle.write('>%s\n' % read_record.id)
                fasta_out_r2_handle.write('%s\n' % str(read_record.seq))

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

    # define tmp files
    paired_r1_out_tmp = '%s.tmp' % paired_r1_out
    paired_r2_out_tmp = '%s.tmp' % paired_r2_out
    singleton_out_tmp = '%s.tmp' % singleton_out

    # prepare commands
    subsample_cmd_paired_r1 = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, paired_r1, paired_r1_out_tmp, subsample_pct)
    subsample_cmd_paired_r2 = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, paired_r2, paired_r2_out_tmp, subsample_pct)
    subsample_cmd_singleton = '%s -fastx_subsample %s -fastaout %s -sample_pct %s -randseed 1 -quiet' % (usearch_exe, singleton_in, singleton_out_tmp, subsample_pct)

    # execute commands
    os.system(subsample_cmd_paired_r1)
    os.system(subsample_cmd_paired_r2)
    os.system(subsample_cmd_singleton)

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
    fasta_in_path, fasta_in_basename, fasta_in_ext = sep_path_basename_ext(sortmerna_op)
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


def parse_uclust_output(seq_file_in, uclust_output_table, seq_file_out, cluster_to_member_file):

    seq_len_dict = {}
    for seq_record in SeqIO.parse(seq_file_in, 'fasta'):
        seq_len_dict[seq_record.id] = len(seq_record.seq)

    cluster_id_set = set()
    cluster_to_rep_seq_dict = {}
    cluster_to_seq_member_dict = {}
    for each_line in open(uclust_output_table):
        each_line_split = each_line.strip().split('\t')
        cluster_id = each_line_split[1]
        seq_id = each_line_split[8].split(' ')[0]
        cluster_id_set.add(int(cluster_id))

        if cluster_id not in cluster_to_rep_seq_dict:
            cluster_to_rep_seq_dict[cluster_id] = seq_id
            cluster_to_seq_member_dict[cluster_id] = {seq_id}
        else:
            if seq_len_dict[seq_id] > seq_len_dict[cluster_to_rep_seq_dict[cluster_id]]:
                cluster_to_rep_seq_dict[cluster_id] = seq_id
            cluster_to_seq_member_dict[cluster_id].add(seq_id)

    # write out cluster sequence members
    cluster_to_member_file_handle = open(cluster_to_member_file, 'w')
    for each_cluster in sorted([i for i in cluster_id_set]):
        cluster_to_member_file_handle.write('Cluster_%s\t%s\n' % (each_cluster, ','.join(sorted([i for i in cluster_to_seq_member_dict[str(each_cluster)]]))))
    cluster_to_member_file_handle.close()

    # write out the longest sequence in each cluster
    rep_seq_id_set = set()
    for i in cluster_to_rep_seq_dict:
        rep_seq_id_set.add(cluster_to_rep_seq_dict[i])

    rep_seq_file_handle = open(seq_file_out, 'w')
    for each_seq in SeqIO.parse(seq_file_in, 'fasta'):
        if each_seq.id in rep_seq_id_set:
            rep_seq_file_handle.write('>%s\n' % each_seq.id)
            rep_seq_file_handle.write('%s\n' % str(each_seq.seq))
    rep_seq_file_handle.close()


def matam_16s(args):

    ###################################################### file in #####################################################

    # file in
    output_prefix                   = args['p']
    reads_file_r1                   = args['r1']
    reads_file_r2                   = args['r2']
    subsample_pcts                  = args['pct']
    matam_ref                       = args['ref']
    uclust_iden_cutoff              = args['i']
    num_threads                     = args['t']
    force_overwrite                 = args['force']
    matam_assembly_script           = args['matam_assembly']
    sortmerna_exe                   = args['sortmerna']
    usearch_exe                     = args['usearch']
    seqtk_exe                       = args['seqtk']
    keep_quiet                      = args['quiet']


    ################################################ check dependencies ################################################

    program_list = [usearch_exe, seqtk_exe, sortmerna_exe, matam_assembly_script]

    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ', '.join(not_detected_programs))
        exit()


    ####################################################################################################################

    subsample_pct_list  = str_to_num_list(subsample_pcts)
    sortmerna_ref       = '%s.clustered.fasta,%s.clustered' % (matam_ref, matam_ref)
    matam16s_wd         = '%s_Matam16S_wd'                  % (output_prefix)
    log_file            = '%s/%s_matam_16s.log'             % (matam16s_wd, output_prefix)
    combined_r1_r2      = '%s/%s_combined_R1_R2.fasta'      % (matam16s_wd, output_prefix)
    sortmerna_op_fasta  = '%s/%s.fasta'                     % (matam16s_wd, output_prefix)
    sortmerna_stdout    = '%s/%s.SortMeRNA_stdout.txt'      % (matam16s_wd, output_prefix)


    # create folder
    if (os.path.isdir(matam16s_wd) is True) and (force_overwrite is False):
        print('Output folder detected, program exited: %s' % matam16s_wd)
        exit()
    else:
        force_create_folder(matam16s_wd)

    # combine R1 and R2
    report_and_log(('Combining the forward and reverse reads'), log_file, keep_quiet)
    combined_r1_r2_cmd = '%s mergepe %s %s > %s' % (seqtk_exe, reads_file_r1, reads_file_r2, combined_r1_r2)
    os.system(combined_r1_r2_cmd)

    # run SortMeRNA
    report_and_log(('Extracting rRNA reads with SortMeRNA from %s' % combined_r1_r2), log_file, keep_quiet)

    sortmerna_cmd = '%s --ref %s --reads %s --aligned %s/%s --fastx --sam --blast "1" --log --best 10 --min_lis 10 -e 1.00e-05 -a %s -v > %s' % (sortmerna_exe, sortmerna_ref, combined_r1_r2, matam16s_wd, output_prefix, num_threads, sortmerna_stdout)
    report_and_log(sortmerna_cmd, log_file, True)
    os.system(sortmerna_cmd)

    # subsample SortMeRNA output and assemble
    renamed_matam_assembly_list = []
    for subsample_pct in subsample_pct_list:

        report_and_log(('Subsample RNA reads at %s percent' % subsample_pct), log_file, keep_quiet)

        subsample_reads_file = '%s/%s_subsample_%s.fasta'    % (matam16s_wd, output_prefix, subsample_pct)
        matam_output_folder  = '%s/%s_subsample_%s_Matam_wd' % (matam16s_wd, output_prefix, subsample_pct)

        # subsample
        subsample_sortmerna_output(sortmerna_op_fasta, subsample_pct, subsample_reads_file, usearch_exe, seqtk_exe)

        # assemble with Matam
        report_and_log(('Assembling subsampled reads with Matam'), log_file, keep_quiet)

        matam_cmd = 'python3 %s -d %s -i %s --cpu %s --max_memory 100000 -v -o %s' % (matam_assembly_script, matam_ref, subsample_reads_file, num_threads, matam_output_folder)
        report_and_log(matam_cmd, log_file, True)
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


    # assemble with Matam without subsample
    report_and_log(('Assembling nonsubsampled RNA reads'), log_file, keep_quiet)
    matam_output_folder_no_subsample   = '%s/%s_nonsubsampled_Matam_wd'                                 % (matam16s_wd, output_prefix)
    matam_cmd_no_subsample             = 'python3 %s -d %s -i %s --cpu %s --max_memory 100000 -v -o %s' % (matam_assembly_script, matam_ref, sortmerna_op_fasta, num_threads, matam_output_folder_no_subsample)
    report_and_log(matam_cmd_no_subsample, log_file, True)
    os.system(matam_cmd_no_subsample)


    # prefix no subsampled Matam assemblies
    matam_assemblies_no_subsample          = '%s/workdir/scaffolds.NR.min_500bp.fa'          % (matam_output_folder_no_subsample)
    matam_assemblies_no_subsample_prefixed = '%s/workdir/scaffolds.NR.min_500bp.prefixed.fa' % (matam_output_folder_no_subsample)
    seq_prefix                             = '%s_no_subsample'                               % (output_prefix)
    prefix_seq(matam_assemblies_no_subsample, seq_prefix, matam_assemblies_no_subsample_prefixed)
    renamed_matam_assembly_list.append(matam_assemblies_no_subsample_prefixed)


    # combine Matam outputs
    report_and_log(('Combine Matam assemblies at all depth'), log_file, keep_quiet)
    combined_all_depth_matam_assemblies    = '%s/%s_all_depth_assemblies.fasta' % (matam16s_wd, output_prefix)
    combine_cmd = 'cat %s > %s' % (' '.join(renamed_matam_assembly_list), combined_all_depth_matam_assemblies)
    report_and_log(combine_cmd, log_file, True)
    os.system(combine_cmd)


    # dereplicate combined assemblies with Usearch
    report_and_log(('Dereplicate combined Matam assemblies at %s identity cutoff' % uclust_iden_cutoff), log_file, keep_quiet)
    uclust_output_fasta             = '%s/%s_all_depth_assemblies.dereplicated.fasta'                       % (matam16s_wd, output_prefix)
    uclust_output_table             = '%s/%s_all_depth_assemblies.uc'                                       % (matam16s_wd, output_prefix)
    cluster_to_member_file          = '%s/%s_all_depth_assemblies.uc.reorganised.txt'                       % (matam16s_wd, output_prefix)
    default_centroids_to_be_ignored = '%s/%s_all_depth_assemblies_default_centroids.fasta'                  % (matam16s_wd, output_prefix)
    uclust_cmd                      = '%s -cluster_fast %s -id %s -centroids %s -uc %s -sort length -quiet' % (usearch_exe, combined_all_depth_matam_assemblies, uclust_iden_cutoff, default_centroids_to_be_ignored, uclust_output_table)
    report_and_log(uclust_cmd, log_file, True)
    os.system(uclust_cmd)


    # parse_uclust_output
    parse_uclust_output(combined_all_depth_matam_assemblies, uclust_output_table, uclust_output_fasta, cluster_to_member_file)


    # report
    report_and_log(('Dereplicated Matam assemblies exported to %s' % uclust_output_fasta), log_file, keep_quiet)


if __name__ == '__main__':

    matam_16s_parser = argparse.ArgumentParser()

    matam_16s_parser.add_argument('-p',                 required=True,                                          help='output prefix')
    matam_16s_parser.add_argument('-r1',                required=True,                                          help='paired reads r1')
    matam_16s_parser.add_argument('-r2',                required=True,                                          help='paired reads r2')
    matam_16s_parser.add_argument('-pct',               required=True,  type=str, default='1,5,10,25,50,75',    help='subsample percentage, deafault: 1,5,10,25,50,75')
    matam_16s_parser.add_argument('-ref',               required=False, type=str,                               help='Path to Matam reference database')
    matam_16s_parser.add_argument('-i',                 required=False, type=float, default=0.995,              help='cluster identity cutoff (0-1), default: 0.995')
    matam_16s_parser.add_argument('-t',                 required=False, type=int, default=1,                    help='number of threads, default: 1')
    matam_16s_parser.add_argument('-force',             required=False, action="store_true",                    help='force overwrite existing results')
    matam_16s_parser.add_argument('-quiet',             required=False, action="store_true",                    help='not report progress')
    matam_16s_parser.add_argument('-matam_assembly',    required=False, type=str, default='matam_assembly.py',  help='path to matam_assembly.py, default: matam_assembly.py')
    matam_16s_parser.add_argument('-sortmerna',         required=False, type=str, default='sortmerna',          help='path to sortmerna executable file, default: sortmerna')
    matam_16s_parser.add_argument('-seqtk',             required=False, type=str, default='seqtk',              help='path to seqtk executable file, default: seqtk')
    matam_16s_parser.add_argument('-usearch',           required=False, type=str, default='usearch',            help='path to usearch executable file, default: usearch')

    args = vars(matam_16s_parser.parse_args())

    matam_16s(args)

'''

To do:
1. support to customize parameters for Matam
2. silent Matam

'''
