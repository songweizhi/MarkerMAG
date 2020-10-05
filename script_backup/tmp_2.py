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
#from MarkerMAG.MarkerMAG_config import config_dict


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


pwd_samfile                         = '/Users/songweizhi/Desktop/test_bowtie/f5_16S_a.sam'
clipping_reads_matched_part         = '/Users/songweizhi/Desktop/test_bowtie/clipping_reads_matched_part.txt'
clipping_reads_not_matched_part_seq = '/Users/songweizhi/Desktop/test_bowtie/clipping_reads_not_matched_part.fasta'
min_cigar_M = 30
min_cigar_S = 30

reads_file_r1 = '/Users/songweizhi/Desktop/test_bowtie/f5_R1.fasta'
reads_file_r2 = '/Users/songweizhi/Desktop/test_bowtie/f5_R2.fasta'
unmapped_paired_reads_folder = '/Users/songweizhi/Desktop/test_bowtie/unmapped_paired_reads_folder'
unmapped_paired_reads_file = '/Users/songweizhi/Desktop/test_bowtie/unmapped_paired_reads.fasta'
blast_db = '/Users/songweizhi/Desktop/test_bowtie/f5.fna'

unmapped_paired_reads_blastn = '/Users/songweizhi/Desktop/test_bowtie/unmapped_paired_reads_blastn.tab'
clipping_reads_not_matched_part_seq_blastn = '/Users/songweizhi/Desktop/test_bowtie/clipping_reads_not_matched_part_blastn.tab'

paired_reads_match_profile = '/Users/songweizhi/Desktop/test_bowtie/paired_reads_match_profile.txt'


'''
cd /Users/songweizhi/Desktop/test_bowtie
export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH

BioSAK Reads_simulator -r f5.fna -n 200000 -l 150 -i 200 -split

bowtie2-build -f f5_16S.ffn f5_16S 
bowtie2 -x f5_16S -1 f5_R1.fasta -2 f5_R2.fasta -S f5_16S_a.sam --local --no-unal -f -a 

'''


##################################################### extract reads ####################################################

clipping_mapped_reads_list = set()
perfectly_mapped_reads_dict = {}
clipping_reads_mapped_part_dict = {}
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
        cigar_splitted = cigar_splitter(cigar)

        read_id_with_ref_pos = '%s__c__%s__c__%s' % (read_id, ref_id, ref_pos)

        # for perfectly mapped reads (allow mismatch?)
        if ('M' in cigar) and (len(cigar_splitted) == 1):
            if read_id_base not in perfectly_mapped_reads_dict:
                perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
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

                    # write out the sequence of unmapped part
                    clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id_with_ref_pos)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_right + '\n')

                    # store the match info of mapped part
                    if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                if cigar_splitted[1][-1] == 'M':

                    # write out the sequence of unmapped part
                    clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id_with_ref_pos)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_left + '\n')

                    # store the match info of mapped part
                    if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                clipping_mapped_reads_list.add(read_id)

clipping_reads_not_matched_part_seq_handle.close()


########################################## extract reads with multiprocessing ##########################################

perfectly_mapped_read_singleton_dict = {}
for perfectly_mapped_read in perfectly_mapped_reads_dict:
    current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
    if len(current_value) == 1:
        perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value

# get the id of unmapped paired reads to extract
solely_perfectly_mapped_reads_r1 = set()
solely_perfectly_mapped_reads_r2 = set()
for perfectly_mapped_read in perfectly_mapped_read_singleton_dict:
    current_value = perfectly_mapped_read_singleton_dict[perfectly_mapped_read]
    strand = list(current_value.keys())[0]

    if strand == '1':
        r2_to_extract = '%s.2' % perfectly_mapped_read
        if r2_to_extract not in clipping_mapped_reads_list:
            solely_perfectly_mapped_reads_r2.add(r2_to_extract)

    if strand == '2':
        r1_to_extract = '%s.1' % perfectly_mapped_read
        if r1_to_extract not in clipping_mapped_reads_list:
            solely_perfectly_mapped_reads_r1.add(r1_to_extract)


# extract reads
#os.mkdir(unmapped_paired_reads_folder)

# extract reads with multiprocessing
#extracted_reads_with_multiprocessing(reads_file_r1, reads_file_r2, solely_perfectly_mapped_reads_r1, solely_perfectly_mapped_reads_r2, unmapped_paired_reads_folder, 4)

# combine extracted reads
#os.system('cat %s/*.fasta > %s' % (unmapped_paired_reads_folder, unmapped_paired_reads_file))


############################# run blast between extracted reads and metagenomic assemblies #############################

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 4'

# run blastn
makeblastdb_cmd     = 'makeblastdb -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (blast_db)
blastn_cmd_paired   = 'blastn -query %s -db %s -out %s %s'                          % (unmapped_paired_reads_file, blast_db, unmapped_paired_reads_blastn, blast_parameters)
blastn_cmd_clipping = 'blastn -query %s -db %s -out %s %s'                          % (clipping_reads_not_matched_part_seq, blast_db, clipping_reads_not_matched_part_seq_blastn, blast_parameters)
#
# os.system(makeblastdb_cmd)
# os.system(blastn_cmd_paired)
# os.system(blastn_cmd_clipping)
#

######################################### parse blast results for paired reads #########################################

reads_iden_cutoff = 100
reads_cov_cutoff = 100
unmapped_paired_reads_to_ctg_dict = blast_results_to_dict(unmapped_paired_reads_blastn, reads_iden_cutoff, reads_cov_cutoff)
clipping_parts_to_ctg_dict = blast_results_to_dict(clipping_reads_not_matched_part_seq_blastn, reads_iden_cutoff, reads_cov_cutoff)


# combine match info of r1 and r2
paired_stats_dict_num = {}
paired_reads_match_profile_handle = open(paired_reads_match_profile, 'w')
paired_reads_match_profile_handle.write('ID\tR1\tR2\n')
for unmapped_read in unmapped_paired_reads_to_ctg_dict:
    unmapped_read_base = '.'.join(unmapped_read.split('.')[:-1])
    unmapped_read_strand = unmapped_read.split('.')[-1]

    print(unmapped_read)
    print(unmapped_read_base)
    print(unmapped_read_strand)
    print()

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
                paired_reads_match_profile_handle.write('%s\t%s\t%s\n' % (unmapped_read_base, r1, r2))

                # store in dict
                paired_key = '_|_'.join(sorted([r1, r2])[::-1])
                if genomic_seq_type == 'mag':
                    paired_key = '___'.join(paired_key.split('___')[:-1])
                if paired_key not in paired_stats_dict_num:
                    paired_stats_dict_num[paired_key] = 1
                else:
                    paired_stats_dict_num[paired_key] += 1

paired_reads_match_profile_handle.close()


print(unmapped_paired_reads_to_ctg_dict)
print(perfectly_mapped_read_singleton_dict)


# rm rmp file
#os.system('rm -r %s' % unmapped_paired_reads_folder)



