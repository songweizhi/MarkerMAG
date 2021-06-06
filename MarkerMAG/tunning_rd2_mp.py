import os
import glob
import shutil
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
from datetime import datetime
import plotly.graph_objects as go
from Bio.SeqRecord import SeqRecord
from distutils.spawn import find_executable


####################################################################################################################
############################################### second round linking ###############################################
####################################################################################################################

step_2_wd = '/Users/songweizhi/Desktop/tunning_rd2'
min_M_len_ctg = 60
end_seq_len = 500
marker_to_ctg_gnm_Key_connector = '___M___'
gnm_to_ctg_connector = '___C___'
mini_assembly_to_16s_ctg_connector = '___Mini___'
read_to_marker_connector = '___r___'
ctg_level_min_link = 3
end_ctg_len_for_mafft = 1000
keep_short_M = True
gap_N_num = 50
report_interval = 25000
clp_pct_ctg_side_max_num = 65
clp_pct_ratio_cutoff = 3.5
# mismatch cutoff for filtering matches between unmapped mates and clipping reads against contig end
mismatch_ctg_ends = 1
subsample_rate_for_depth_estimation = 0.1  # between 0 and 1
min_M_len_mini = min_M_len_ctg
min_M_len_16s = 45
min_M_len_ctg = 45
min_M_pct = 35
mismatch_cutoff = 2
min_aln_16s = 500
min_cov_16s = 30
mean_depth_dict_gnm = {}
mean_depth_dict_16s = {}
min_16s_gnm_multiple = 0
min_iden_16s = 98
min_link_num = 8
within_gnm_linkage_num_diff = 80
time_format = '[%Y-%m-%d %H:%M:%S]'
num_threads = 4
pwd_makeblastdb_exe = 'makeblastdb'
pwd_blastn_exe = 'blastn'
pwd_bowtie2_build_exe = 'bowtie2-build'
pwd_bowtie2_exe = 'bowtie2'
pwd_samtools_exe = 'samtools'
pwd_bbmap_exe = 'bbmap.sh'
pwd_spades_exe = 'spades.py'
seqtk_exe = 'seqtk'

combined_1st_round_unlinked_mags                = '%s/round_1_unlinked_gnm.fa'                      % step_2_wd
combined_1st_round_unlinked_mag_end_seq         = '%s/round_1_unlinked_gnm_end_%sbp.fa'             % (step_2_wd, end_seq_len)
rd1_unlinked_mag_end_seq_no_ext                 = '%s/round_1_unlinked_gnm_end_%sbp'                % (step_2_wd, end_seq_len)
rd1_unlinked_mags_sam_bowtie_log                = '%s/round_1_unlinked_gnm_bowtie.log'              % step_2_wd
rd1_unlinked_mags_sam_bowtie                    = '%s/round_1_unlinked_gnm_bowtie.sam'              % step_2_wd
rd1_unlinked_mags_sam_bowtie_reformat           = '%s/round_1_unlinked_gnm_bowtie_reformatted.sam'  % step_2_wd
rd1_unlinked_mags_sam_bowtie_reformat_log       = '%s/round_1_unlinked_gnm_bowtie_reformat.log'     % step_2_wd
rd1_unlinked_mags_sam_bowtie_reformat_sorted    = '%s/round_1_unlinked_gnm_bowtie_reformatted_sorted.sam'  % step_2_wd
rd1_unlinked_mags_sam_line_num                  = '%s/round_1_unlinked_gnm_line_num.txt'            % step_2_wd
rd1_unlinked_mags_sam_split_folder              = '%s/round_1_unlinked_gnm_split'                   % step_2_wd
rd1_unlinked_mags_sam_MappingRecord_folder      = '%s/round_1_unlinked_gnm_MappingRecord'           % step_2_wd
stats_GapFilling_file                           = '%s/stats_GapFilling_gnm.txt'                     % step_2_wd
stats_GapFilling_file_filtered                  = '%s/stats_GapFilling_gnm_filtered.txt'            % step_2_wd
free_living_16s_ref_file                        = '%s/round2_free_living_16s_refs.txt'              % step_2_wd
free_living_ctg_ref_file                        = '%s/round2_free_living_ctg_refs.txt'              % step_2_wd
free_living_all                                 = '%s/round2_free_living_all.fa'                    % step_2_wd
free_living_all_fq_r1                           = '%s/round2_free_living_all_R1.fastq'              % step_2_wd
free_living_all_fq_r2                           = '%s/round2_free_living_all_R2.fastq'              % step_2_wd
free_living_all_fq                              = '%s/round2_free_living_all.fastq'                 % step_2_wd
spades_wd                                       = '%s/mini_assembly_SPAdes_wd'                      % step_2_wd
spades_log                                      = '%s/SPAdes_stdout.txt'                            % step_2_wd
mira_manifest                                   = '%s/mira_manifest.txt'                            % step_2_wd
mira_stdout                                     = '%s/mira_stdout.txt'                              % step_2_wd
mini_assembly_to_16s_reads                      = '%s/mini_assembly_to_16s_reads.txt'               % step_2_wd
mini_assembly_to_ctg_reads                      = '%s/mini_assembly_to_ctg_reads.txt'               % step_2_wd
stats_GapFilling_ctg                            = '%s/stats_GapFilling_ctg.txt'                     % step_2_wd
sam_file_mini_assembly                          = '%s/scaffolds_bowtie.sam'                         % step_2_wd
sam_file_mini_assembly_log                      = '%s/scaffolds_bowtie.log'                         % step_2_wd
sam_file_mini_assembly_reformatted              = '%s/scaffolds_bowtie_reformatted.sam'             % step_2_wd
sam_file_mini_assembly_reformatted_log          = '%s/scaffolds_bowtie_reformatted.log'             % step_2_wd
rd2_to_extract_flking_16s_r1_id                 = '%s/rd2_to_extract_flking_16s_r1_id.txt'          % step_2_wd
rd2_to_extract_flking_16s_r2_id                 = '%s/rd2_to_extract_flking_16s_r2_id.txt'          % step_2_wd
rd2_extracted_flking_16s_r1_seq_tmp             = '%s/rd2_extracted_flking_16s_r1_tmp.fa'           % step_2_wd
rd2_extracted_flking_16s_r2_seq_tmp             = '%s/rd2_extracted_flking_16s_r2_tmp.fa'           % step_2_wd
rd2_to_extract_flking_ctg_r1_id                 = '%s/rd2_to_extract_flking_ctg_r1_id.txt'          % step_2_wd
rd2_to_extract_flking_ctg_r2_id                 = '%s/rd2_to_extract_flking_ctg_r2_id.txt'          % step_2_wd
rd2_extracted_flking_ctg_r1_seq_tmp             = '%s/rd2_extracted_flking_ctg_r1_tmp.fa'           % step_2_wd
rd2_extracted_flking_ctg_r2_seq_tmp             = '%s/rd2_extracted_flking_ctg_r2_tmp.fa'           % step_2_wd
rd2_extracted_r1_combined                       = '%s/rd2_extracted_r1_combined.fa'                 % step_2_wd
rd2_extracted_r2_combined                       = '%s/rd2_extracted_r2_combined.fa'                 % step_2_wd


##################################################### read in sam file ####################################################

round_2_ctg_end_seq_len_dict = {}
for each_ctg_end_record in SeqIO.parse(combined_1st_round_unlinked_mag_end_seq, 'fasta'):
    round_2_ctg_end_seq_len_dict[each_ctg_end_record.id] = len(each_ctg_end_record.seq)


# get splitted sam file list
splitted_gnm_sam_file_list = [os.path.basename(file_name) for file_name in glob.glob('%s/splitted_sam_*' % rd1_unlinked_mags_sam_split_folder)]

# prepare lol for mp worker
list_for_parse_sam_gnm_worker = []
splitted_sam_mp_file_set = set()
for splitted_gnm_sam_file in splitted_gnm_sam_file_list:
    pwd_splitted_gnm_sam_file                       = '%s/%s'                           % (rd1_unlinked_mags_sam_split_folder, splitted_gnm_sam_file)
    pwd_splitted_gnm_sam_free_living_ctg_ref_file   = '%s/%s_free_living_ctg_refs.txt'  % (rd1_unlinked_mags_sam_MappingRecord_folder, splitted_gnm_sam_file)
    splitted_sam_mp_file_set.add(pwd_splitted_gnm_sam_free_living_ctg_ref_file)
    list_for_parse_sam_gnm_worker.append([pwd_splitted_gnm_sam_file,
                                         pwd_splitted_gnm_sam_free_living_ctg_ref_file,
                                         min_M_len_ctg,
                                         mismatch_cutoff,
                                         round_2_ctg_end_seq_len_dict])

pool_parse_sam_gnm = mp.Pool(processes=num_threads)
pool_parse_sam_gnm.map(parse_sam_gnm_worker, list_for_parse_sam_gnm_worker)
pool_parse_sam_gnm.close()
pool_parse_sam_gnm.join()

os.system('rm -r %s' % rd1_unlinked_mags_sam_split_folder)

# combine free_living_ctg_ref_files
os.system('cat %s > %s' % (' '.join(splitted_sam_mp_file_set), free_living_ctg_ref_file))

##################################################### extract reads flanking contig ends ####################################################

to_extract_read_base_rd2_ctg = set()
for each_read_base in open(free_living_ctg_ref_file):
    to_extract_read_base_rd2_ctg.add(each_read_base.strip().split('\t')[0])

# write out id of linking reads for extraction
with open(rd2_to_extract_flking_ctg_r1_id, 'w') as rd2_to_extract_flking_ctg_r1_id_handle:
    rd2_to_extract_flking_ctg_r1_id_handle.write('%s\n' % '\n'.join(sorted([('%s.1' % i) for i in to_extract_read_base_rd2_ctg])))
with open(rd2_to_extract_flking_ctg_r2_id, 'w') as rd2_to_extract_flking_ctg_r2_id_handle:
    rd2_to_extract_flking_ctg_r2_id_handle.write('%s\n' % '\n'.join(sorted([('%s.2' % i) for i in to_extract_read_base_rd2_ctg])))

# extract reads with seqtk
seqtk_extract_cmd_rd2_flk_ctg_r1 = 'seqtk subseq %s %s > %s' % (reads_file_r1_fasta, rd2_to_extract_flking_ctg_r1_id, rd2_extracted_flking_ctg_r1_seq_tmp)
seqtk_extract_cmd_rd2_flk_ctg_r2 = 'seqtk subseq %s %s > %s' % (reads_file_r2_fasta, rd2_to_extract_flking_ctg_r2_id, rd2_extracted_flking_ctg_r2_seq_tmp)
os.system(seqtk_extract_cmd_rd2_flk_ctg_r1)
os.system(seqtk_extract_cmd_rd2_flk_ctg_r2)


##################################################### combine reads for assembly ####################################################

# read extracted read sequences into dict
extract_rd2_read_seq_dict = {}
for extracted_r1 in SeqIO.parse(rd2_extracted_flking_16s_r1_seq_tmp, 'fasta'):
    extract_rd2_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
for extracted_r2 in SeqIO.parse(rd2_extracted_flking_16s_r2_seq_tmp, 'fasta'):
    extract_rd2_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)
for extracted_r1 in SeqIO.parse(rd2_extracted_flking_ctg_r1_seq_tmp, 'fasta'):
    extract_rd2_read_seq_dict[extracted_r1.id] = str(extracted_r1.seq)
for extracted_r2 in SeqIO.parse(rd2_extracted_flking_ctg_r2_seq_tmp, 'fasta'):
    extract_rd2_read_seq_dict[extracted_r2.id] = str(extracted_r2.seq)

# write out paired in the same order
rd2_extracted_r1_combined_handle = open(rd2_extracted_r1_combined, 'w')
rd2_extracted_r2_combined_handle = open(rd2_extracted_r2_combined, 'w')
for each_read_base in set.union(rd2_read_to_extract_flanking_16s, to_extract_read_base_rd2_ctg):
    current_r1 = '%s.1' % each_read_base
    current_r2 = '%s.2' % each_read_base
    current_r1_seq = extract_rd2_read_seq_dict.get(current_r1, '')
    current_r2_seq = extract_rd2_read_seq_dict.get(current_r2, '')
    rd2_extracted_r1_combined_handle.write('>%s\n' % current_r1)
    rd2_extracted_r1_combined_handle.write('%s\n' % current_r1_seq)
    rd2_extracted_r2_combined_handle.write('>%s\n' % current_r2)
    rd2_extracted_r2_combined_handle.write('%s\n' % current_r2_seq)
rd2_extracted_r1_combined_handle.close()
rd2_extracted_r2_combined_handle.close()


####################################################################################################################
######################################### second round linking by assembly #########################################
####################################################################################################################

# assemble
spades_cmd = '%s --only-assembler --meta -1 %s -2 %s -o %s -t %s -k 55,75,99,127 > %s' % (pwd_spades_exe, rd2_extracted_r1_combined, rd2_extracted_r2_combined, spades_wd, num_threads, spades_log)
spades_cmd = '%s --only-assembler --meta -1 %s -2 %s -o %s -t %s -k 55,75,99 > %s' % (pwd_spades_exe, rd2_extracted_r1_combined, rd2_extracted_r2_combined, spades_wd, num_threads, spades_log)
os.system(spades_cmd)
mini_assemblies = '%s/scaffolds.fasta' % spades_wd

##################################################### mapping reads to mini_assemblies ####################################################

# index miniassembly
mini_assemblies_no_ext = '.'.join(mini_assemblies.split('.')[:-1])
bowtie_build_mini_assemblies_cmd = 'bowtie2-build --quiet --threads %s -f %s %s' % (num_threads, mini_assemblies, mini_assemblies_no_ext)
os.system(bowtie_build_mini_assemblies_cmd)

# mapping with bowtie
bowtie_cmd_miniassembly = 'bowtie2 -x %s -U %s,%s -S %s -p %s -f %s 2> %s'    % (mini_assemblies_no_ext, rd2_extracted_r1_combined, rd2_extracted_r2_combined, sam_file_mini_assembly, num_threads, bowtie_parameter_mini_assembly, sam_file_mini_assembly_log)
os.system(bowtie_cmd_miniassembly)

# reformat cigar string
mini_assembly_reformat_cmd = 'reformat.sh in=%s out=%s sam=1.4 2> %s' % (sam_file_mini_assembly, sam_file_mini_assembly_reformatted, sam_file_mini_assembly_reformatted_log)
os.system(mini_assembly_reformat_cmd)


##################################################### read in sam file ####################################################

mini_assembly_len_dict = {}
mini_assembly_mp_dict = {}
with open(sam_file_mini_assembly_reformatted) as sam_file_mini_assembly_reformatted_opened:
    for each_read in sam_file_mini_assembly_reformatted_opened:
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            mini_assembly_id = ''
            mini_assembly_len = 0
            for each_element in each_read_split:
                if each_element.startswith('SN:'):
                    mini_assembly_id = each_element[3:]
                if each_element.startswith('LN:'):
                    mini_assembly_len = int(each_element[3:])
            mini_assembly_len_dict[mini_assembly_id] = mini_assembly_len
        else:
            cigar = each_read_split[5]
            if cigar != '*':
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]
                ref_id = each_read_split[2]
                ref_pos = int(each_read_split[3])

                if read_id_base not in mini_assembly_mp_dict:
                    mini_assembly_mp_dict[read_id_base] = MappingRecord()

                if read_strand == '1':
                    if ref_id not in mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict:
                        mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict[ref_id] = {ref_pos: cigar}
                    else:
                        mini_assembly_mp_dict[read_id_base].r1_mini_ref_dict[ref_id][ref_pos] = cigar

                if read_strand == '2':
                    if ref_id not in mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict:
                        mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict[ref_id] = {ref_pos: cigar}
                    else:
                        mini_assembly_mp_dict[read_id_base].r2_mini_ref_dict[ref_id][ref_pos] = cigar


#################################################### parse sam file ####################################################

# get mini assembly to 16s reads dict
mini_assembly_to_16s_reads_dict = {}
mini_assembly_to_ctg_reads_dict = {}
for each_mp in mini_assembly_mp_dict:

    current_mp_record = mini_assembly_mp_dict[each_mp]
    mini_refs_to_ignore = set()

    ########## get lowest mismatch for r1/r2 16s refs ##########

    # get r1_ref_cigar_set
    r1_ref_cigar_set = set()
    for each_pos_dict in current_mp_record.r1_mini_ref_dict.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r1_ref_cigar_set.update(each_pos_dict_values)

    # get r2_ref_cigar_set
    r2_ref_cigar_set = set()
    for each_pos_dict in current_mp_record.r2_mini_ref_dict.values():
        each_pos_dict_values = {each_pos_dict[i] for i in each_pos_dict}
        r2_ref_cigar_set.update(each_pos_dict_values)

    r1_ref_min_mismatch = get_min_mismatch_from_cigar_list(r1_ref_cigar_set, min_M_len_mini)
    r2_ref_min_mismatch = get_min_mismatch_from_cigar_list(r2_ref_cigar_set, min_M_len_mini)
    mini_assembly_mp_dict[each_mp].r1_mini_refs_lowest_mismatch = r1_ref_min_mismatch
    mini_assembly_mp_dict[each_mp].r2_mini_refs_lowest_mismatch = r2_ref_min_mismatch

    ########## filter r1 mini refs ##########

    r1_mini_refs_passed_qc = {}
    for r1_mini_ref in current_mp_record.r1_mini_ref_dict:
        r1_matched_pos_dict = current_mp_record.r1_mini_ref_dict[r1_mini_ref]

        # one read need to mapped to one 16S only for one time
        if len(r1_matched_pos_dict) > 1:
            mini_refs_to_ignore.add(r1_mini_ref)
        else:
            r1_mini_ref_pos = list(r1_matched_pos_dict.keys())[0]
            r1_mini_ref_cigar = r1_matched_pos_dict[r1_mini_ref_pos]
            r1_mini_ref_cigar_splitted = cigar_splitter(r1_mini_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r1_mini_ref_cigar_splitted)
            if both_end_clp is True:
                mini_refs_to_ignore.add(r1_mini_ref)
            else:
                # check mismatch
                r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(
                    r1_mini_ref_cigar_splitted)
                if r1_ref_min_mismatch == 'NA':
                    mini_refs_to_ignore.add(r1_mini_ref)
                elif (r1_mismatch_pct > r1_ref_min_mismatch) or (r1_mismatch_pct > mismatch_cutoff):
                    mini_refs_to_ignore.add(r1_mini_ref)
                else:
                    # check aligned length
                    if r1_aligned_len < min_M_len_mini:
                        mini_refs_to_ignore.add(r1_mini_ref)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r1_mini_ref_cigar) or ('s' in r1_mini_ref_cigar):
                            clip_in_middle = True
                            if (r1_mini_ref_cigar_splitted[0][-1] in ['S', 's']) and (r1_mini_ref_pos == 1):
                                clip_in_middle = False
                            if (r1_mini_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r1_mini_ref_pos + r1_aligned_len - 1) == mini_assembly_len_dict[r1_mini_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            mini_refs_to_ignore.add(r1_mini_ref)
                        else:
                            r1_mini_refs_passed_qc[r1_mini_ref] = [r1_mini_ref_cigar]

    ########## filter r2 mini refs ##########

    r2_mini_refs_passed_qc = {}
    for r2_mini_ref in current_mp_record.r2_mini_ref_dict:
        r2_matched_pos_dict = current_mp_record.r2_mini_ref_dict[r2_mini_ref]

        # one read need to mapped to one 16S only once
        if len(r2_matched_pos_dict) > 1:
            mini_refs_to_ignore.add(r2_mini_ref)
        else:
            r2_mini_ref_pos = list(r2_matched_pos_dict.keys())[0]
            r2_mini_ref_cigar = r2_matched_pos_dict[r2_mini_ref_pos]
            r2_mini_ref_cigar_splitted = cigar_splitter(r2_mini_ref_cigar)

            # check both end clip
            both_end_clp = check_both_ends_clipping(r2_mini_ref_cigar_splitted)
            if both_end_clp is True:
                mini_refs_to_ignore.add(r2_mini_ref)
            else:
                # check mismatch
                r2_aligned_len, r2_aligned_pct, r2_clipping_len, r2_clipping_pct, r2_mismatch_pct = get_cigar_stats(
                    r2_mini_ref_cigar_splitted)
                if r2_ref_min_mismatch == 'NA':
                    mini_refs_to_ignore.add(r2_mini_ref)
                elif (r2_mismatch_pct > r2_ref_min_mismatch) or (r2_mismatch_pct > mismatch_cutoff):
                    mini_refs_to_ignore.add(r2_mini_ref)
                else:
                    # check aligned length
                    if r2_aligned_len < min_M_len_mini:
                        mini_refs_to_ignore.add(r2_mini_ref)
                    else:
                        # check if clp in the middle
                        clip_in_middle = False
                        if ('S' in r2_mini_ref_cigar) or ('s' in r2_mini_ref_cigar):
                            clip_in_middle = True
                            if (r2_mini_ref_cigar_splitted[0][-1] in ['S', 's']) and (r2_mini_ref_pos == 1):
                                clip_in_middle = False
                            if (r2_mini_ref_cigar_splitted[-1][-1] in ['S', 's']):
                                if (r2_mini_ref_pos + r2_aligned_len - 1) == mini_assembly_len_dict[r2_mini_ref]:
                                    clip_in_middle = False

                        # exclude the ref if clp in the middle is True
                        if clip_in_middle is True:
                            mini_refs_to_ignore.add(r2_mini_ref)
                        else:
                            r2_mini_refs_passed_qc[r2_mini_ref] = [r2_mini_ref_cigar]

    ####################################################################################################

    r1_mini_refs_no_ignored = {key: value for key, value in r1_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
    r2_mini_refs_no_ignored = {key: value for key, value in r2_mini_refs_passed_qc.items() if key not in mini_refs_to_ignore}
    shared_mini_ref_dict = {key: [r1_mini_refs_no_ignored[key][0], r2_mini_refs_no_ignored[key][0]] for key in set(r1_mini_refs_no_ignored).intersection(set(r2_mini_refs_no_ignored))}

    # both r1 and r2 mapped to the same mini assembly
    if len(shared_mini_ref_dict) > 0:

        if (each_mp in rd2_read_to_extract_flanking_16s) and (each_mp not in to_extract_read_base_rd2_ctg):
            for each_mini_assembly in shared_mini_ref_dict:
                if each_mini_assembly not in mini_assembly_to_16s_reads_dict:
                    mini_assembly_to_16s_reads_dict[each_mini_assembly] = {each_mp}
                else:
                    mini_assembly_to_16s_reads_dict[each_mini_assembly].add(each_mp)

        elif (each_mp not in rd2_read_to_extract_flanking_16s) and (each_mp in to_extract_read_base_rd2_ctg):
            for each_mini_assembly in shared_mini_ref_dict:
                if each_mini_assembly not in mini_assembly_to_ctg_reads_dict:
                    mini_assembly_to_ctg_reads_dict[each_mini_assembly] = {each_mp}
                else:
                    mini_assembly_to_ctg_reads_dict[each_mini_assembly].add(each_mp)

# write out mini_assembly_to_16s_reads
mini_assembly_to_16s_reads_handle = open(mini_assembly_to_16s_reads, 'w')
for each_mini_assembly in mini_assembly_to_16s_reads_dict:
    mini_assembly_to_16s_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join(mini_assembly_to_16s_reads_dict[each_mini_assembly])))
mini_assembly_to_16s_reads_handle.close()

# write out mini_assembly_to_ctg_reads
mini_assembly_to_ctg_reads_handle = open(mini_assembly_to_ctg_reads, 'w')
for each_mini_assembly in mini_assembly_to_ctg_reads_dict:
    mini_assembly_to_ctg_reads_handle.write('%s\t%s\n' % (each_mini_assembly, ','.join(mini_assembly_to_ctg_reads_dict[each_mini_assembly])))
mini_assembly_to_ctg_reads_handle.close()


############################################# get_GapFilling_stats #############################################

get_GapFilling_stats_by_assembly(free_living_16s_ref_file,
                                 free_living_ctg_ref_file,
                                 mini_assembly_to_16s_reads,
                                 mini_assembly_to_ctg_reads,
                                 ctg_level_min_link,
                                 mini_assembly_to_16s_ctg_connector,
                                 gnm_to_ctg_connector,
                                 marker_to_ctg_gnm_Key_connector,
                                 within_gnm_linkage_num_diff,
                                 max_mini_assembly_link_num_diff_between_ctg_16s,
                                 stats_GapFilling_ctg,
                                 stats_GapFilling_file)

filter_linkages_iteratively(stats_GapFilling_file, 'Number', pairwise_16s_iden_dict, mean_depth_dict_gnm, mean_depth_dict_16s, min_16s_gnm_multiple, min_iden_16s, min_link_num, min_link_num, within_gnm_linkage_num_diff, stats_GapFilling_file_filtered)

free_living_16s_to_ctg_linkage_dict_to_use = {}
for each_ctg_level_link in open(stats_GapFilling_ctg):
    each_ctg_level_link_split = each_ctg_level_link.split('\t')
    marker_id = each_ctg_level_link_split[0]
    ctg_id = each_ctg_level_link_split[1]
    link_num = int(each_ctg_level_link_split[2])
    current_key = '%s%s%s' % (marker_id, marker_to_ctg_gnm_Key_connector, ctg_id)
    free_living_16s_to_ctg_linkage_dict_to_use[current_key] = link_num

