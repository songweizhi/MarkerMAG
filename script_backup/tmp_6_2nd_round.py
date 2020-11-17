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


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_free_living_mate(ref_in, reads_r1, reads_r2, end_seq_len, minCigarM, max_gap_to_end, bowtie_build_exe, bowtie2_exe, num_threads):

    ref_in_path, ref_in_basename, ref_in_ext = sep_path_basename_ext(ref_in)

    ref_subset = '%s/%s_ends_%sbp%s'   % (ref_in_path, ref_in_basename, end_seq_len, ref_in_ext)
    sam_file   = '%s/%s_ends_%sbp.sam' % (ref_in_path, ref_in_basename, end_seq_len)

    # get ref seqs subset
    ref_subset_len_dict = {}
    ref_subset_handle = open(ref_subset, 'w')
    for ref_seq in SeqIO.parse(ref_in, 'fasta'):
        ref_seq_id = ref_seq.id
        ref_seq_len = len(ref_seq.seq)
        if ref_seq_len < end_seq_len*3:
            ref_subset_handle.write('>%s\n' % ref_seq_id)
            ref_subset_handle.write('%s\n' % ref_seq.seq)
            ref_subset_len_dict[ref_seq_id] = ref_seq_len
        else:
            ref_seq_left_end_id     = '%s_l' % ref_seq_id
            ref_seq_right_end_id    = '%s_r' % ref_seq_id
            ref_seq_left_end        = ref_seq.seq[:end_seq_len]
            ref_seq_right_end       = ref_seq.seq[-end_seq_len:]

            # write out left end
            ref_subset_handle.write('>%s\n' % ref_seq_left_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_left_end)

            # write out right end
            ref_subset_handle.write('>%s\n' % ref_seq_right_end_id)
            ref_subset_handle.write('%s\n' % ref_seq_right_end)
    ref_subset_handle.close()

    # index ref seq subset
    ref_subset_no_ext = '.'.join(ref_subset.split('.')[:-1])
    bowtie2_index_ref_cmd = '%s -f %s %s --quiet --threads %s' % (bowtie_build_exe, ref_subset, ref_subset_no_ext, num_threads)
    os.system(bowtie2_index_ref_cmd)

    # mapping
    bowtie2_mapping_cmd = '%s -x %s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s' % (bowtie2_exe, ref_subset_no_ext, reads_r1, reads_r2, sam_file, num_threads)
    os.system(bowtie2_mapping_cmd)

    # parse sam file
    all_mapped_reads = set()
    qualified_reads_dict = {}
    qualified_reads_to_ref_dict = {}
    for each_line in open(sam_file):
        if not each_line.startswith('@'):
            each_line_split = each_line.strip().split('\t')
            read_id         = each_line_split[0]
            read_id_base    = '.'.join(read_id.split('.')[:-1])
            read_strand     = read_id.split('.')[-1]
            ref_id          = each_line_split[2]
            ref_pos         = int(each_line_split[3])
            cigar           = each_line_split[5]
            cigar_splitted  = cigar_splitter(cigar)

            all_mapped_reads.add(read_id)

            qualified_mapping = False
            if ref_id[-2:] == '_l':
                if ref_pos <= max_gap_to_end:
                    if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and(int(cigar[:-1]) >= minCigarM):
                        qualified_mapping = True

                    elif (len(cigar_splitted) == 2) and (cigar_splitted[-1][-1] == 'M') and (cigar_splitted[0][-1] == 'S') and (int(cigar_splitted[-1][:-1]) >= minCigarM) and (ref_pos == 1):
                        qualified_mapping = True

            elif ref_id[-2:] == '_r':

                if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                    ref_pos_end = ref_pos + int(cigar[:-1])

                    if (end_seq_len - ref_pos_end) <= max_gap_to_end:
                        qualified_mapping = True

                elif (len(cigar_splitted) == 2) and (cigar_splitted[0][-1] == 'M') and (cigar_splitted[1][-1] == 'S') and (int(cigar_splitted[0][:-1]) >= minCigarM):
                    if (ref_pos + int(cigar_splitted[0][:-1]) - 1) == end_seq_len:
                        qualified_mapping = True

            else:
                ref_len = ref_subset_len_dict[ref_id]
                if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                    ref_pos_end = ref_pos + int(cigar[:-1])

                    # left side
                    if ref_pos <= max_gap_to_end:
                        qualified_mapping = True

                    # right side
                    elif (ref_len - ref_pos_end) <= max_gap_to_end:
                        qualified_mapping = True

                if len(cigar_splitted) == 2:

                    # left side
                    if (cigar_splitted[-1][-1] == 'M') and (cigar_splitted[0][-1] == 'S') and (int(cigar_splitted[-1][:-1]) >= minCigarM) and (ref_pos == 1):
                        qualified_mapping = True

                    # right side
                    elif (cigar_splitted[-0][-1] == 'M') and (cigar_splitted[1][-1] == 'S') and (int(cigar_splitted[0][:-1]) >= minCigarM):
                        if (ref_pos + int(cigar_splitted[0][:-1]) - 1) == ref_len:
                            qualified_mapping = True

            if qualified_mapping is True:
                qualified_reads_to_ref_dict[read_id] = ref_id

                if read_id_base not in qualified_reads_dict:
                    qualified_reads_dict[read_id_base] = [read_strand]
                else:
                    qualified_reads_dict[read_id_base].append(read_strand)


    reads_to_extract_to_ref_dict = {}
    for qualified_read in qualified_reads_dict:
        read_strand = qualified_reads_dict[qualified_read]
        if len(read_strand) == 1:

            mapped_mate = ''
            mate_to_extract = ''
            if read_strand == ['1']:
                mapped_mate     = '%s.1' % (qualified_read)
                mate_to_extract = '%s.2' % (qualified_read)
            if read_strand == ['2']:
                mapped_mate     = '%s.2' % (qualified_read)
                mate_to_extract = '%s.1' % (qualified_read)

            if mate_to_extract not in all_mapped_reads:
                reads_to_extract_to_ref_dict[mate_to_extract] = qualified_reads_to_ref_dict[mapped_mate]

    return reads_to_extract_to_ref_dict


########################################################################################################################

wd = '/Users/songweizhi/Desktop/777'

marker_gene_seqs_1st_round_unlinked    = '%s/rest_16S.fasta' % wd
combined_1st_round_unlinked_mags              = '%s/rest_mags.fna' % wd
reads_file_r1       = '%s/MBARC26_R1_0.05.fasta'    % wd
reads_file_r2       = '%s/MBARC26_R2_0.05.fasta'    % wd
end_seq_len         = 3000
minCigarM           = 50
max_gap_to_end      = 100
num_threads         = 4

pwd_bowtie2_build_exe    = '/Users/songweizhi/Softwares/bowtie2/bowtie2-build'
pwd_bowtie2_exe         = '/Users/songweizhi/Softwares/bowtie2/bowtie2'
pwd_spades_exe          = '/Users/songweizhi/Softwares/SPAdes-3.14.1/bin/spades.py'

pwd_bowtie2_build_exe    = 'bowtie2-build'
pwd_bowtie2_exe         = 'bowtie2'
pwd_spades_exe          = 'spades.py'


# define output
free_living_mate_gnm        = '%s/free_living_mate_ctg.txt'             % wd
free_living_mate_16s        = '%s/free_living_mate_16s.txt'             % wd
extracted_reads_folder      = '%s/free_living_reads'                    % wd
extracted_reads_cbd         = '%s/free_living_read_combined.fasta'      % wd
spades_wd                   = '%s/combined_free_living_reads_SPAdes_wd' % wd
spades_assemblies           = '%s/scaffolds.fasta'                      % spades_wd
sam_file                    = '%s/scaffolds.sam'                        % wd


########################################################################################################################

reads_to_extract_to_ref_dict_gnm = get_free_living_mate(combined_1st_round_unlinked_mags, reads_file_r1, reads_file_r2, end_seq_len, minCigarM, max_gap_to_end, pwd_bowtie2_build_exe, pwd_bowtie2_exe, num_threads)
reads_to_extract_to_ref_dict_16s = get_free_living_mate(marker_gene_seqs_1st_round_unlinked, reads_file_r1, reads_file_r2, end_seq_len, minCigarM, max_gap_to_end, pwd_bowtie2_build_exe, pwd_bowtie2_exe, num_threads)

# write out free_living_mate_gnm
free_living_mate_gnm_handle = open(free_living_mate_gnm, 'w')
for fl_read_gnm in reads_to_extract_to_ref_dict_gnm:
    free_living_mate_gnm_handle.write('%s\t%s\n' % (fl_read_gnm, reads_to_extract_to_ref_dict_gnm[fl_read_gnm]))
free_living_mate_gnm_handle.close()

# write out free_living_mate_16s
free_living_mate_16s_handle = open(free_living_mate_16s, 'w')
for fl_read_16s in reads_to_extract_to_ref_dict_16s:
    free_living_mate_16s_handle.write('%s\t%s\n' % (fl_read_16s, reads_to_extract_to_ref_dict_16s[fl_read_16s]))
free_living_mate_16s_handle.close()

# reads_to_extract_to_ref_dict_gnm = {}
# for i in open(free_living_mate_gnm):
#     i_split = i.strip().split('\t')
#     reads_to_extract_to_ref_dict_gnm[i_split[0]] = i_split[1]
#
# reads_to_extract_to_ref_dict_16s = {}
# for j in open(free_living_mate_16s):
#     j_split = j.strip().split('\t')
#     reads_to_extract_to_ref_dict_16s[j_split[0]] = j_split[1]


####################################################  extract reads ####################################################

os.mkdir(extracted_reads_folder)

extract_list_gnm = set()
extract_list_16s = set()
extract_list_combined_r1 = set()
extract_list_combined_r2 = set()
for fl_mate_gnm in reads_to_extract_to_ref_dict_gnm:
    extract_list_gnm.add(fl_mate_gnm)
    if fl_mate_gnm[-2:] == '.1':
        extract_list_combined_r1.add(fl_mate_gnm)
    if fl_mate_gnm[-2:] == '.2':
        extract_list_combined_r2.add(fl_mate_gnm)
for fl_mate_16s in reads_to_extract_to_ref_dict_16s:
    extract_list_16s.add(fl_mate_16s)
    if fl_mate_16s[-2:] == '.1':
        extract_list_combined_r1.add(fl_mate_16s)
    if fl_mate_16s[-2:] == '.2':
        extract_list_combined_r2.add(fl_mate_16s)


# remove overlap
extract_list_combined_r1_no_overlap = set()
for r1 in extract_list_combined_r1:
    if (r1 in extract_list_gnm) and (r1 in extract_list_16s):
        pass
    else:
        extract_list_combined_r1_no_overlap.add(r1)

extract_list_combined_r2_no_overlap = set()
for r2 in extract_list_combined_r2:
    if (r2 in extract_list_gnm) and (r2 in extract_list_16s):
        pass
    else:
        extract_list_combined_r2_no_overlap.add(r2)

# extract reads with multiprocessing
extracted_reads_with_multiprocessing(reads_file_r1, reads_file_r2, extract_list_combined_r1_no_overlap, extract_list_combined_r2_no_overlap, extracted_reads_folder, num_threads)

# combine extracted reads
os.system('cat %s/*.fasta > %s' % (extracted_reads_folder, extracted_reads_cbd))
os.system('rm -r %s' % extracted_reads_folder)


####################################################### assemble #######################################################

spades_cmd = '%s -s %s -o %s -t %s -k 21,33,55,75,99,127 --only-assembler' % (pwd_spades_exe, extracted_reads_cbd, spades_wd, num_threads)
os.system(spades_cmd)


######################################################## mapping #######################################################

spades_assemblies_no_ext = '.'.join(spades_assemblies.split('.')[:-1])
index_ref_cmd = '%s -f %s %s' % (pwd_bowtie2_build_exe, spades_assemblies, spades_assemblies_no_ext)
mapping_cmd = '%s -x %s -U %s -S %s -p %s -f' % (pwd_bowtie2_exe, spades_assemblies_no_ext, extracted_reads_cbd, sam_file, num_threads)
os.system(index_ref_cmd)
os.system(mapping_cmd)


#################################################### parse sam file ####################################################

short_ctg_to_reads_dict = {}
for each_line in open(sam_file):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id         = each_line_split[0]
        read_id_base    = '.'.join(read_id.split('.')[:-1])
        read_strand     = read_id.split('.')[-1]
        ref_id          = each_line_split[2]
        ref_pos         = int(each_line_split[3])
        cigar           = each_line_split[5]
        cigar_splitted  = cigar_splitter(cigar)
        if (len(cigar_splitted) == 1) and (cigar[-1] == 'M'):
            if ref_id not in short_ctg_to_reads_dict:
                short_ctg_to_reads_dict[ref_id] = [read_id]
            else:
                short_ctg_to_reads_dict[ref_id].append(read_id)


for short_ctg in short_ctg_to_reads_dict:
    short_ctg_mapped_reads = short_ctg_to_reads_dict[short_ctg]
    for mapped_read in short_ctg_mapped_reads:

        mapped_read_mate_ref = ''
        if mapped_read in reads_to_extract_to_ref_dict_gnm:
            mapped_read_mate_ref = reads_to_extract_to_ref_dict_gnm[mapped_read]
        if mapped_read in reads_to_extract_to_ref_dict_16s:
            mapped_read_mate_ref = reads_to_extract_to_ref_dict_16s[mapped_read]

        print('%s\t%s' % (short_ctg, mapped_read_mate_ref))

    print()

'''
export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH
export PATH=/Users/songweizhi/Softwares/samtools-1.11/bin:$PATH
cd /Users/songweizhi/Desktop/777
cat free_living_read_16s.fasta free_living_read_gnm.fasta > free_living_read_cbd.fasta
bowtie2-build -f combined_refs.fa combined_refs
bowtie2 -x combined_refs -U free_living_read_16s.fasta -S free_living_read_16s.sam -p 4 -f
bowtie2 -x combined_refs -U free_living_read_gnm.fasta -S free_living_read_gnm.sam -p 4 -f
bowtie2 -x combined_refs -U free_living_read_cbd.fasta -S free_living_read_cbd.sam -p 4 -f
samtools view -bS free_living_read_16s.sam -o free_living_read_16s.bam
samtools view -bS free_living_read_gnm.sam -o free_living_read_gnm.bam
samtools view -bS free_living_read_cbd.sam -o free_living_read_cbd.bam
samtools sort free_living_read_16s.bam -o free_living_read_16s_sorted.bam
samtools sort free_living_read_gnm.bam -o free_living_read_gnm_sorted.bam
samtools sort free_living_read_cbd.bam -o free_living_read_cbd_sorted.bam
samtools index free_living_read_16s_sorted.bam
samtools index free_living_read_gnm_sorted.bam
samtools index free_living_read_cbd_sorted.bam

module load spades/3.14.0
cd /srv/scratch/z5039045/MarkerMAG_wd/new_algorithm
spades.py -s free_living_read_cbd.fasta -o free_living_read_cbd_spades -t 4 --only-assembler
cd /srv/scratch/z5039045/MarkerMAG_wd/new_algorithm/free_living_read_cbd_spades
cd /Users/songweizhi/Desktop/777
bowtie2-build -f scaffolds.fasta scaffolds
bowtie2 -x scaffolds -U free_living_read_cbd.fasta -S scaffolds.sam -p 4 -f
samtools view -bS scaffolds.sam -o scaffolds.bam
samtools sort scaffolds.bam -o scaffolds_sorted.bam
samtools index scaffolds_sorted.bam

cd /Users/songweizhi/Desktop/777
/Users/songweizhi/Softwares/SPAdes-3.14.1/bin/spades.py -s combined_free_living_reads.fasta -o combined_free_living_reads_spades_wd -t 4 --only-assembler

# on Katana
module load python/3.7.3
source ~/mypython3env/bin/activate
module load spades/3.14.0
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/new_algorithm
python3 tmp_6_2nd_round.py

'''
