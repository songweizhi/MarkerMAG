import os
import glob
import shutil
import argparse
from Bio import SeqIO
from time import sleep
from datetime import datetime
from BioSAK.BioSAK_config import config_dict


link_16S_MAG_usage = '''
=========================================== link_16S_MAG example commands ===========================================

1. rename paired reads
link_16S_MAG assumes the id of paired reads in a format of XXXX.1 and XXXX.2. The only difference of their names is 
the last character. "BioSAK rename_reads_for_Reago" do this job.
 
BioSAK link_16S_MAG -p Test -r1 ST13_R1.fasta -r2 ST13_R2.fasta -m ST13_16S_seqs.fa -mag ST13_ctgs.fasta -t 4
BioSAK link_16S_MAG -p Test -r1 ST13_R1.fasta -r2 ST13_R2.fasta -m ST13_16S_seqs.fa -mag ST13_MAGs -x fa -t 4

=====================================================================================================================
'''

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


def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if each_element.isalpha() is True:
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


def remove_unmapped_reads(sam_file, sam_file_only_mapped):

    sam_file_only_mapped_handle = open(sam_file_only_mapped, 'w')

    for each_read in open(sam_file):
        if not each_read.startswith('@'):
            ref_id = each_read.strip().split('\t')[2]
            cigar = each_read.strip().split('\t')[5]

            if (ref_id != '*') and (cigar != '*'):
                sam_file_only_mapped_handle.write(each_read)

    sam_file_only_mapped_handle.close()


def blast_results_to_dict(blastn_results, iden_cutoff, query_cov_cutoff):
    query_to_subject_list_dict = {}

    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        subject = blast_hit_split[1]
        subject_with_prefix = 'GenomicSeq__%s' % subject
        iden = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_len = int(blast_hit_split[12])
        coverage_q = float(align_len) * 100 / float(query_len)
        if (iden >= iden_cutoff) and (coverage_q >= query_cov_cutoff):
            if query not in query_to_subject_list_dict:
                query_to_subject_list_dict[query] = [subject_with_prefix]
            else:
                query_to_subject_list_dict[query].append(subject_with_prefix)

    return query_to_subject_list_dict


def stats_dict_to_sankey_file_in(clipping_stats_dict, paired_stats_dict, value_type, signal_cutoff, sankey_file_in_clipping, sankey_file_in_paired, sankey_file_in_combined):

    file_header = ''
    if value_type == 'reads_number':
        file_header = 'GenomicSeq,MarkerGene,Number\n'
    if value_type == 'reads_length':
        file_header = 'GenomicSeq,MarkerGene,Length_bp\n'

    contig_to_plot_paired = set()
    contig_to_plot_clipping = set()
    contig_to_plot_combined = set()

    # prepare input file for plot of clipping mapped reads
    sankey_file_in_combined_handle = open(sankey_file_in_combined, 'w')
    sankey_file_in_clipping_handle = open(sankey_file_in_clipping, 'w')
    sankey_file_in_combined_handle.write(file_header)
    sankey_file_in_clipping_handle.write(file_header)
    for each_clipping in clipping_stats_dict:
        if clipping_stats_dict[each_clipping] >= signal_cutoff:
            sankey_file_in_clipping_handle.write('%s,%s\n' % (','.join(each_clipping.split('_|_')), clipping_stats_dict[each_clipping]))
            sankey_file_in_combined_handle.write('%s,%s\n' % (','.join(each_clipping.split('_|_')), clipping_stats_dict[each_clipping]))
            contig_to_plot_clipping.add(each_clipping.split('_|_')[1])
            contig_to_plot_combined.add(each_clipping.split('_|_')[1])
    sankey_file_in_clipping_handle.close()

    # prepare input file for plot of paired reads
    sankey_file_in_paired_handle = open(sankey_file_in_paired, 'w')
    sankey_file_in_paired_handle.write(file_header)
    for each_paired in paired_stats_dict:
        if paired_stats_dict[each_paired] >= signal_cutoff:
            sankey_file_in_paired_handle.write('%s,%s\n' % (','.join(each_paired.split('_|_')), paired_stats_dict[each_paired]))
            sankey_file_in_combined_handle.write('%s,%s\n' % (','.join(each_paired.split('_|_')), paired_stats_dict[each_paired]))
            contig_to_plot_paired.add(each_paired.split('_|_')[1])
            contig_to_plot_combined.add(each_paired.split('_|_')[1])
    sankey_file_in_paired_handle.close()
    sankey_file_in_combined_handle.close()

    return len(contig_to_plot_paired), len(contig_to_plot_clipping), len(contig_to_plot_combined)


def rename_seq(ctg_file_in, ctg_file_out, prefix):

    ctg_file_out_handle = open(ctg_file_out, 'w')
    for Seq_record in SeqIO.parse(ctg_file_in, 'fasta'):
        Seq_record.id = '%s___%s' % (prefix, Seq_record.id)
        SeqIO.write(Seq_record, ctg_file_out_handle, 'fasta')
    ctg_file_out_handle.close()


def link_16S_MAG(args, config_dict):

    ###################################################### file in/out #####################################################

    # file in
    output_prefix                       = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    genomic_assemblies                  = args['g']
    mag_folder                          = args['mag']
    mag_file_extension                  = args['x']
    marker_gene_seqs                    = args['m']
    min_cigar_M                         = args['cigarM']
    min_cigar_S                         = args['cigarS']
    iden_cutoff                         = args['iden']
    cov_cutoff                          = args['cov']
    min_num_to_plot                     = args['num']
    min_len_to_plot                     = args['len']
    skip_mapping                         = args['skip_mapping']
    skip_blastn                          = args['skip_blastn']
    force_overwrite                     = args['force']
    keep_quiet                          = args['quiet']
    num_threads                         = args['t']
    keep_temp                           = args['tmp']

    pwd_plot_sankey_R                   = config_dict['get_sankey_plot_R']
    pwd_bowtie2_exe                     = config_dict['bowtie2']
    pwd_bowtie2_build_exe               = config_dict['bowtie2_build']
    pwd_makeblastdb_exe                 = config_dict['makeblastdb']
    pwd_blastn_exe                      = config_dict['blastn']


    ################################################ check dependencies ################################################


    ############################################# create working directory #############################################

    # create working directory
    working_directory = '%s_link_16S_and_MAG_wd' % output_prefix
    if (os.path.isdir(working_directory) is True) and (force_overwrite is False):
        print('Working directory detected, program exited!')
        exit()
    else:
        force_create_folder(working_directory)


    ######################## check genomic sequence type and prepare files for making blast db #########################

    blast_db         = ''
    genomic_seq_type = ''  # ctg or mag

    # check the type of input genomic sequences
    if (genomic_assemblies is not None) and (mag_folder is None):
        genomic_seq_type = 'ctg'
        metagenomic_assemblies_file_path, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension = sep_path_basename_ext(genomic_assemblies)
        blast_db_dir = '%s/%s_%s_blast_db' % (working_directory, output_prefix, metagenomic_assemblies_file_basename)
        blast_db     = '%s/%s%s'           % (blast_db_dir, metagenomic_assemblies_file_basename, metagenomic_assemblies_file_extension)
        os.mkdir(blast_db_dir)
        os.system('cp %s %s/' % (genomic_assemblies, blast_db_dir))

    elif (genomic_assemblies is None) and (mag_folder is not None):
        genomic_seq_type = 'mag'
        mag_folder_name     = mag_folder.split('/')[-1]
        blast_db_dir        = '%s/%s_blast_db'            % (working_directory, mag_folder_name)
        renamed_mag_folder  = '%s/%s_blast_db/%s_renamed' % (working_directory, mag_folder_name, mag_folder_name)
        os.mkdir(blast_db_dir)
        os.mkdir(renamed_mag_folder)

        # get input mag file list
        mag_file_re = '%s/*%s' % (mag_folder, mag_file_extension)
        mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
        if len(mag_file_list) == 0:
            print('No MAG detected, program exited!')
            exit()

        # add mag id to its sequences
        for mag_in in mag_file_list:
            pwd_mag_in      = '%s/%s' % (mag_folder, mag_in)
            pwd_mag_renamed = '%s/%s' % (renamed_mag_folder, mag_in)
            mag_basename    = '.'.join(mag_in.split('.')[:-1])
            rename_seq(pwd_mag_in, pwd_mag_renamed, mag_basename)

        # combine renamed MAGs
        blast_db = '%s/%s_combined.fa' % (blast_db_dir, mag_folder_name)
        os.system('cat %s/*%s > %s' % (renamed_mag_folder, mag_file_extension, blast_db))

    else:
        print('Please provide genomic sequences either as raw assemblies (-g) or as MAGs (-mag)')
        exit()


    ########################################### define folder and file name ############################################

    marker_gene_seqs_file_path, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension = sep_path_basename_ext(marker_gene_seqs)

    pwd_log_file                                = '%s/%s_%s.log.txt'                                                        % (working_directory, output_prefix, datetime.now().strftime('%Y-%m-%d_%Hh-%Mm-%Ss_%f'))
    bowtie_index_dir                            = '%s/%s_%s_index'                                                          % (working_directory, output_prefix, marker_gene_seqs_file_basename)
    pwd_samfile                                 = '%s/%s.sam'                                                               % (working_directory, marker_gene_seqs_file_basename)
    clipping_reads_matched_part                 = '%s/%s_clipping_matched_part.txt'                                         % (working_directory, output_prefix)
    clipping_reads_not_matched_part_seq         = '%s/%s_clipping_not_matched_part_seq.fasta'                               % (working_directory, output_prefix)
    clipping_reads_not_matched_part_seq_blastn  = '%s/%s_clipping_not_matched_part_seq_blast.txt'                           % (working_directory, output_prefix)
    clipping_reads_match_profile                = '%s/%s_match_profile_clipping.txt'                                        % (working_directory, output_prefix)
    extracted_reads_file                        = '%s/%s_paired_unmapped.fasta'                                             % (working_directory, output_prefix)
    extracted_reads_blastn                      = '%s/%s_paired_unmapped_blast.txt'                                         % (working_directory, output_prefix)
    paired_reads_match_profile                  = '%s/%s_match_profile_paired.txt'                                          % (working_directory, output_prefix)
    clipping_stats_num                          = '%s/%s_plot_cigarM%s_cigarS%s_iden%s_cov%s_MinNum%s_clipping.txt'         % (working_directory, output_prefix, min_cigar_M, min_cigar_S, iden_cutoff, cov_cutoff, min_num_to_plot)
    paired_stats_num                            = '%s/%s_plot_cigarM%s_cigarS%s_iden%s_cov%s_MinNum%s_paired.txt'           % (working_directory, output_prefix, min_cigar_M, min_cigar_S, iden_cutoff, cov_cutoff, min_num_to_plot)
    combined_stats_num                          = '%s/%s_plot_cigarM%s_cigarS%s_iden%s_cov%s_MinNum%s_combined.txt'         % (working_directory, output_prefix, min_cigar_M, min_cigar_S, iden_cutoff, cov_cutoff, min_num_to_plot)
    clipping_stats_len                          = '%s/%s_plot_cigarM%s_cigarS%s_iden%s_cov%s_MinLen%sbp_clipping.txt'       % (working_directory, output_prefix, min_cigar_M, min_cigar_S, iden_cutoff, cov_cutoff, min_len_to_plot)
    paired_stats_len                            = '%s/%s_plot_cigarM%s_cigarS%s_iden%s_cov%s_MinLen%sbp_paired.txt'         % (working_directory, output_prefix, min_cigar_M, min_cigar_S, iden_cutoff, cov_cutoff, min_len_to_plot)
    combined_stats_len                          = '%s/%s_plot_cigarM%s_cigarS%s_iden%s_cov%s_MinLen%sbp_combined.txt'       % (working_directory, output_prefix, min_cigar_M, min_cigar_S, iden_cutoff, cov_cutoff, min_len_to_plot)


    ######################################## map reads to marker gene sequences ########################################

    if skip_mapping is False:

        # copy marker gene sequence file to index folder
        report_and_log(('Indexing marker gene sequences for mapping'), pwd_log_file, keep_quiet)
        os.mkdir(bowtie_index_dir)
        os.system('cp %s %s/' % (marker_gene_seqs, bowtie_index_dir))

        # run mapping
        report_and_log(('Mapping reads to marker gene sequences'), pwd_log_file, keep_quiet)

        sleep(1)
        report_and_log(('If you see any warning messages starting with "Use of uninitialized value" during bowtie mapping, please ignore them :)'), pwd_log_file, keep_quiet)
        sleep(1)

        bowtie2_index_ref_cmd = '%s -f %s/%s%s %s/%s --quiet --threads %s' % (pwd_bowtie2_build_exe, bowtie_index_dir, marker_gene_seqs_file_basename, marker_gene_seqs_file_extension, bowtie_index_dir, marker_gene_seqs_file_basename, num_threads)
        bowtie2_mapping_cmd   = '%s -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s' % (pwd_bowtie2_exe,  bowtie_index_dir, marker_gene_seqs_file_basename, reads_file_r1, reads_file_r2, pwd_samfile, num_threads)
        os.system(bowtie2_index_ref_cmd)
        os.system(bowtie2_mapping_cmd)

        sleep(1)
        report_and_log(('If you see any warning messages starting with "Use of uninitialized value" during bowtie mapping, please ignore them :)'), pwd_log_file, keep_quiet)
        sleep(1)


    ##################################################### extract reads ####################################################

    report_and_log(('Extract unmapped paired reads and unmapped part of clipping mapped reads'), pwd_log_file, keep_quiet)

    # export clipping mapped reads and perfectly mapped reads
    clipping_mapped_reads_list = set()
    clipping_reads_mapped_part_dict = {}
    perfectly_mapped_reads_dict = {}
    clipping_reads_matched_part_handle = open(clipping_reads_matched_part, 'w')
    clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
    for each_read in open(pwd_samfile):

        if not each_read.startswith('@'):
            each_read_split = each_read.strip().split('\t')
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
            ref_pos = int(each_read_split[3])
            cigar = each_read_split[5]
            read_seq = each_read_split[9]
            cigar_splitted  = cigar_splitter(cigar)

            # for perfectly mapped reads
            if ('M' in cigar) and (len(cigar_splitted) == 1):
                if read_id_base not in perfectly_mapped_reads_dict:
                    perfectly_mapped_reads_dict[read_id_base] = {read_strand:[ref_id_with_prefix]}
                else:
                    if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                        perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                    else:
                        perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

            # for clipping mapped reads
            if ('S' in cigar) and (len(cigar_splitted) == 2):
                cigar_M_len = 0
                cigar_S_len = 0
                split_pos = 0
                if cigar_splitted[0][-1] == 'M':
                    cigar_M_len = int(cigar_splitted[0][:-1])
                    cigar_S_len = int(cigar_splitted[1][:-1])
                    split_pos = ref_pos + cigar_M_len
                if cigar_splitted[1][-1] == 'M':
                    cigar_M_len = int(cigar_splitted[1][:-1])
                    cigar_S_len = int(cigar_splitted[0][:-1])
                    split_pos = ref_pos

                if (cigar_M_len >= min_cigar_M) and (cigar_S_len >= min_cigar_S):
                    read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                    read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                    if cigar_splitted[0][-1] == 'M':
                        clipping_reads_matched_part_handle.write('%s_l\t%s\t%s\n' % (read_id, ref_id, split_pos))
                        clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id)
                        clipping_reads_not_matched_part_seq_handle.write(read_seq_right + '\n')
                        if ('%s_l' % read_id) not in clipping_reads_mapped_part_dict:
                            clipping_reads_mapped_part_dict[('%s_l' % read_id)] = [ref_id_with_prefix]
                        else:
                            clipping_reads_mapped_part_dict[('%s_l' % read_id)].append(ref_id_with_prefix)

                    if cigar_splitted[1][-1] == 'M':
                        clipping_reads_matched_part_handle.write('%s_r\t%s\t%s\n' % (read_id, ref_id, split_pos))
                        clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id)
                        clipping_reads_not_matched_part_seq_handle.write(read_seq_left + '\n')
                        if ('%s_r' % read_id) not in clipping_reads_mapped_part_dict:
                            clipping_reads_mapped_part_dict[('%s_r' % read_id)] = [ref_id_with_prefix]
                        else:
                            clipping_reads_mapped_part_dict[('%s_r' % read_id)].append(ref_id_with_prefix)

                    clipping_mapped_reads_list.add(read_id)
    clipping_reads_not_matched_part_seq_handle.close()

    perfectly_mapped_read_singleton_dict = {}
    for perfectly_mapped_read in perfectly_mapped_reads_dict:
        current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
        if len(current_value) == 1:
            perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value

    # get the id of paired reads to extract
    solely_perfectly_mapped_reads_r1_test = []
    solely_perfectly_mapped_reads_r2_test = []
    for perfectly_mapped_read in perfectly_mapped_read_singleton_dict:
        current_value = perfectly_mapped_read_singleton_dict[perfectly_mapped_read]
        strand = list(current_value.keys())[0]
        if strand == '1':
            r2_to_extract = '%s.2' % perfectly_mapped_read
            if r2_to_extract not in clipping_mapped_reads_list:
                solely_perfectly_mapped_reads_r2_test.append(r2_to_extract)
        if strand == '2':
            r1_to_extract = '%s.1' % perfectly_mapped_read
            if r1_to_extract not in clipping_mapped_reads_list:
                solely_perfectly_mapped_reads_r1_test.append(r1_to_extract)

    # extract reads
    if skip_mapping is False:
        extracted_reads_file_handle = open(extracted_reads_file, 'w')
        for read_r1 in SeqIO.parse(reads_file_r1, 'fasta'):
            if read_r1.id in solely_perfectly_mapped_reads_r1_test:
                extracted_reads_file_handle.write('>%s\n' % read_r1.id)
                extracted_reads_file_handle.write('%s\n' % read_r1.seq)
        for read_r2 in SeqIO.parse(reads_file_r2, 'fasta'):
            if read_r2.id in solely_perfectly_mapped_reads_r2_test:
                extracted_reads_file_handle.write('>%s\n' % read_r2.id)
                extracted_reads_file_handle.write('%s\n' % read_r2.seq)
        extracted_reads_file_handle.close()

    # store reads length into dict
    reads_len_dict = {}
    for read_r1 in SeqIO.parse(reads_file_r1, 'fasta'):
        reads_len_dict[read_r1.id] = len(read_r1.seq)
    for read_r2 in SeqIO.parse(reads_file_r2, 'fasta'):
        reads_len_dict[read_r2.id] = len(read_r2.seq)


    ############################# run blast between extracted reads and metagenomic assemblies #############################

    # run blast
    if skip_blastn is False:

        # copy genomic assemblies to blast db folder

        # run blastn
        blast_parameters    = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % num_threads
        makeblastdb_cmd     = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blast_db)
        blastn_cmd_paired   = '%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, extracted_reads_file, blast_db, extracted_reads_blastn, blast_parameters)
        blastn_cmd_clipping = '%s -query %s -db %s -out %s %s' % (pwd_blastn_exe, clipping_reads_not_matched_part_seq, blast_db, clipping_reads_not_matched_part_seq_blastn, blast_parameters)

        report_and_log(('Making blastn database'), pwd_log_file, keep_quiet)
        os.system(makeblastdb_cmd)

        report_and_log(('Running blastn for unmapped paired reads'), pwd_log_file, keep_quiet)
        os.system(blastn_cmd_paired)

        report_and_log(('Running blastn for unmapped parts of clipping mapped reads'), pwd_log_file, keep_quiet)
        os.system(blastn_cmd_clipping)


    ######################################### parse blast results for paired reads #########################################

    report_and_log(('Parsing blast results for paired reads'), pwd_log_file, keep_quiet)

    # filter blast results for paired and clipping mapped reads
    unmapped_paired_reads_to_ctg_dict = blast_results_to_dict(extracted_reads_blastn, iden_cutoff, cov_cutoff)
    clipping_parts_to_ctg_dict        = blast_results_to_dict(clipping_reads_not_matched_part_seq_blastn, iden_cutoff, cov_cutoff)
    paired_stats_dict_num = {}
    paired_stats_dict_len = {}
    paired_reads_match_profile_handle = open(paired_reads_match_profile, 'w')
    paired_reads_match_profile_handle.write('ID\tR1\tR2\tPairedTotalLength(bp)\n')
    for unmapped_read in unmapped_paired_reads_to_ctg_dict:

        unmapped_read_base = '.'.join(unmapped_read.split('.')[:-1])
        unmapped_read_strand = unmapped_read.split('.')[-1]

        unmapped_read_base_r1_matched = []
        unmapped_read_base_r2_matched = []
        if unmapped_read_strand == '1':
            unmapped_read_base_r1_matched = unmapped_paired_reads_to_ctg_dict[unmapped_read]
            if unmapped_read_base in perfectly_mapped_read_singleton_dict:
                unmapped_read_base_r2_matched = perfectly_mapped_read_singleton_dict[unmapped_read_base]['2']
        if unmapped_read_strand == '2':
            unmapped_read_base_r2_matched = unmapped_paired_reads_to_ctg_dict[unmapped_read]
            if unmapped_read_base in perfectly_mapped_read_singleton_dict:
                unmapped_read_base_r1_matched = perfectly_mapped_read_singleton_dict[unmapped_read_base]['1']

        if (unmapped_read_base_r1_matched != []) and (unmapped_read_base_r2_matched != []):
            for r1 in unmapped_read_base_r1_matched:
                for r2 in unmapped_read_base_r2_matched:

                    # write out to file
                    paired_total_length = reads_len_dict[unmapped_read_base + '.1'] + reads_len_dict[unmapped_read_base + '.2']
                    paired_reads_match_profile_handle.write('%s\t%s\t%s\t%s\n' % (unmapped_read_base, r1, r2, paired_total_length))

                    # store in dict
                    paired_key = '_|_'.join(sorted([r1, r2])[::-1])
                    if genomic_seq_type == 'mag':
                        paired_key = '___'.join(paired_key.split('___')[:-1])
                    if paired_key not in paired_stats_dict_len:
                        paired_stats_dict_len[paired_key] = paired_total_length
                        paired_stats_dict_num[paired_key] = 1
                    else:
                        paired_stats_dict_len[paired_key] += paired_total_length
                        paired_stats_dict_num[paired_key] += 1

    paired_reads_match_profile_handle.close()


    #################################### parse blast results for clipping mapped reads #####################################

    report_and_log(('Parsing blast results for clipping mapped reads'), pwd_log_file, keep_quiet)

    clipping_reads_match_profile_handle = open(clipping_reads_match_profile, 'w')
    clipping_reads_match_profile_handle.write('ID\tLeft\tRight\tReadLength(bp)\n')
    clipping_stats_dict_num = {}
    clipping_stats_dict_len = {}
    for clipping_mapped_read in clipping_reads_mapped_part_dict:

        clipping_mapped_read_id = clipping_mapped_read[:-2]
        mapped_part = clipping_mapped_read[-1]

        clipping_mapped_read_matches_l = []
        clipping_mapped_read_matches_r = []
        if mapped_part == 'l':
            clipping_mapped_read_matches_l = clipping_reads_mapped_part_dict[clipping_mapped_read]
            if ('%s_r' % clipping_mapped_read_id) in clipping_parts_to_ctg_dict:
                clipping_mapped_read_matches_r = clipping_parts_to_ctg_dict[('%s_r' % clipping_mapped_read_id)]
        if mapped_part == 'r':
            clipping_mapped_read_matches_r = clipping_reads_mapped_part_dict[clipping_mapped_read]
            if ('%s_l' % clipping_mapped_read_id) in clipping_parts_to_ctg_dict:
                clipping_mapped_read_matches_l = clipping_parts_to_ctg_dict[('%s_l' % clipping_mapped_read_id)]

        if (clipping_mapped_read_matches_l != []) and (clipping_mapped_read_matches_r != []):

            for l in clipping_mapped_read_matches_l:
                for r in clipping_mapped_read_matches_r:

                    # write out to file
                    clipping_reads_match_profile_handle.write('%s\t%s\t%s\t%s\n' % (clipping_mapped_read_id, l, r, reads_len_dict[clipping_mapped_read_id]))

                    # store in dict
                    clipping_key = '_|_'.join(sorted([l, r])[::-1])
                    if genomic_seq_type == 'mag':
                        clipping_key = '___'.join(clipping_key.split('___')[:-1])

                    if clipping_key not in clipping_stats_dict_len:
                        clipping_stats_dict_len[clipping_key] = reads_len_dict[clipping_mapped_read_id]
                        clipping_stats_dict_num[clipping_key] = 1
                    else:
                        clipping_stats_dict_len[clipping_key] += reads_len_dict[clipping_mapped_read_id]
                        clipping_stats_dict_num[clipping_key] += 1

    clipping_reads_match_profile_handle.close()


    ######################################################### plot #########################################################

    report_and_log(('Plotting'), pwd_log_file, keep_quiet)

    # prepare input file for sankey plot
    ctgs_to_plot_clipping_num, ctgs_to_plot_paired_num, ctgs_to_plot_combined_num = stats_dict_to_sankey_file_in(clipping_stats_dict_num, paired_stats_dict_num, 'reads_number', min_num_to_plot, clipping_stats_num, paired_stats_num, combined_stats_num)
    ctgs_to_plot_clipping_len, ctgs_to_plot_paired_len, ctgs_to_plot_combined_len = stats_dict_to_sankey_file_in(clipping_stats_dict_len, paired_stats_dict_len, 'reads_length', min_len_to_plot, clipping_stats_len, paired_stats_len, combined_stats_len)

    # get plot height
    plot_height_clipping_num = 450 if ctgs_to_plot_clipping_num <= 30 else ctgs_to_plot_clipping_num * 15
    plot_height_paired_num   = 450 if ctgs_to_plot_paired_num   <= 30 else ctgs_to_plot_paired_num   * 15
    plot_height_combined_num = 450 if ctgs_to_plot_combined_num <= 30 else ctgs_to_plot_combined_num * 15
    plot_height_clipping_len = 450 if ctgs_to_plot_clipping_len <= 30 else ctgs_to_plot_clipping_len * 15
    plot_height_paired_len   = 450 if ctgs_to_plot_paired_len   <= 30 else ctgs_to_plot_paired_len   * 15
    plot_height_combined_len = 450 if ctgs_to_plot_combined_len <= 30 else ctgs_to_plot_combined_len * 15

    # prepare commands
    cmd_sankey_clipping_num = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, clipping_stats_num, 600, plot_height_clipping_num)
    cmd_sankey_clipping_len = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, clipping_stats_len, 600, plot_height_clipping_len)
    cmd_sankey_paired_num   = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, paired_stats_num,   600, plot_height_paired_num)
    cmd_sankey_paired_len   = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, paired_stats_len,   600, plot_height_paired_len)
    cmd_sankey_combined_num = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, combined_stats_num, 600, plot_height_combined_num)
    cmd_sankey_combined_len = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, combined_stats_len, 600, plot_height_combined_len)

    # plot
    #os.system(cmd_sankey_clipping_num)
    #os.system(cmd_sankey_clipping_len)
    #os.system(cmd_sankey_paired_num)
    #os.system(cmd_sankey_paired_len)
    os.system(cmd_sankey_combined_num)
    os.system(cmd_sankey_combined_len)


    ################################################### remove tmp files ###################################################

    if keep_temp is False:

        report_and_log(('Removing temporary files'), pwd_log_file, keep_quiet)
        os.remove(pwd_samfile)
        os.remove(clipping_reads_matched_part)
        os.remove(clipping_reads_not_matched_part_seq)
        os.remove(clipping_reads_not_matched_part_seq_blastn)
        os.remove(extracted_reads_file)
        os.remove(extracted_reads_blastn)
    os.remove(clipping_stats_num)
    os.remove(paired_stats_num)
    os.remove(clipping_stats_len)
    os.remove(paired_stats_len)


    report_and_log(('Done!'), pwd_log_file, keep_quiet)


######################################################### main #########################################################

if __name__ == '__main__':

    link_16S_and_MAG_parser = argparse.ArgumentParser()

    # arguments for link_16S_and_MAG
    link_16S_and_MAG_parser.add_argument('-p',              required=True,                           help='Output prefix')
    link_16S_and_MAG_parser.add_argument('-r1',             required=True,                           help='Paired reads r1')
    link_16S_and_MAG_parser.add_argument('-r2',             required=True,                           help='Paired reads r2')
    link_16S_and_MAG_parser.add_argument('-m',              required=True,                           help='Marker gene sequences')
    link_16S_and_MAG_parser.add_argument('-g',              required=False, default=None,            help='Genomic sequences')
    link_16S_and_MAG_parser.add_argument('-mag',            required=False, default=None,            help='Metagenome-assembled-genome (MAG) folder')
    link_16S_and_MAG_parser.add_argument('-x',              required=False, default='fasta',         help='MAG file extension, default: fasta')
    link_16S_and_MAG_parser.add_argument('-iden',           required=False, type=float, default=100, help='identity cutoff, default: 100')
    link_16S_and_MAG_parser.add_argument('-cov',            required=False, type=float, default=100, help='coverage cutoff, default: 100')
    link_16S_and_MAG_parser.add_argument('-cigarM',         required=False, type=int, default=50,    help='cigarM cutoff, default: 50')
    link_16S_and_MAG_parser.add_argument('-cigarS',         required=False, type=int, default=50,    help='cigarS cutoff, default: 50')
    link_16S_and_MAG_parser.add_argument('-num',            required=False, type=int, default=5,     help='minimum number of reads for a link to plot, default: 5')
    link_16S_and_MAG_parser.add_argument('-len',            required=False, type=int, default=500,   help='minimum length of reads for a link to plot, default: 500')
    link_16S_and_MAG_parser.add_argument('-t',              required=False, type=int, default=1,     help='number of threads, default: 1')
    link_16S_and_MAG_parser.add_argument('-quiet',          required=False, action="store_true",     help='not report progress')
    link_16S_and_MAG_parser.add_argument('-force',          required=False, action="store_true",     help='force overwrite existing results')
    link_16S_and_MAG_parser.add_argument('-tmp',            required=False, action="store_true",     help='keep temporary files')
    link_16S_and_MAG_parser.add_argument('-skip_mapping',   required=False, action="store_true",     help='skip mapping')
    link_16S_and_MAG_parser.add_argument('-skip_blastn',    required=False, action="store_true",     help='skip blastn')

    args = vars(link_16S_and_MAG_parser.parse_args())
    link_16S_MAG(args, config_dict)


To_do = '''

1. enable multiple placement for a read?
2. where does the paired read of the clipping mapped read mapped to? (should take into consideration!!!)
3. check whether the final output is blank
4. how to incorporate the taxonomy of MAGs and 16S sequences
5. rename reads as xxx.1 and xxx.2
6. contig specific signal cutoff
7. the effect of sequencing depth, insert size and read length
8. metric for describing accuracy and completeness


Torsten's comments:
1. get the best match in a iterative way
2. divergence degree of 16S sequences in a genome
3. set a cutoff of 20 is not a good approach
4. separate mate reads and clipping mapped reads
5. put 16S sequences into clusters


export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH
cd /Users/songweizhi/Desktop/link_16S_MAG
python3 ~/PycharmProjects/BioSAK/BioSAK/link_16S_MAG.py -p Shan -r1 file_in/ST13_R1_P_renamed.fasta -r2 file_in/ST13_R2_P_renamed.fasta -m file_in/scaffolds.NR.min_500bp.abd.fa -g file_in/scaffold.fa -force -t 4

export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH
cd /Users/songweizhi/Desktop/link_16S_MAG
python3 ~/PycharmProjects/BioSAK/BioSAK/link_16S_MAG.py -p Shan -r1 file_in/ST13_R1_P_renamed.fasta -r2 file_in/ST13_R2_P_renamed.fasta -m file_in/scaffolds.NR.min_500bp.abd.fa -mag file_in/ST13_MAGs -x fasta -force -t 4
BioSAK link_16S_MAG -p Shan -r1 file_in/ST13_R1_P_renamed.fasta -r2 file_in/ST13_R2_P_renamed.fasta -m file_in/scaffolds.NR.min_500bp.abd.fa -mag file_in/ST13_MAGs -x fasta -force -t 4

export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH
cd /Users/songweizhi/Desktop/link_16S_MAG
python3 ~/PycharmProjects/BioSAK/BioSAK/link_16S_MAG.py -p Madeup -r1 madeup_file_in/reads_R1.fasta -r2 madeup_file_in/reads_R2.fasta -m madeup_file_in/16S_seq.fasta -g madeup_file_in/scaffold.fa -force -t 4 -len 1 -num 1

'''