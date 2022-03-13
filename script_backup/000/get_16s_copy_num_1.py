
def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if (each_element.isalpha() is True) or (each_element == '='):
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


def get_cigar_stats(cigar_splitted):

    # aligned_len: M I X =
    # clipping_len: S
    # mismatch_len: X I D
    # mismatch_pct = mismatch_len / aligned_len
    # aligned_pct  = aligned_len  / (aligned_len + clipping_len)
    # clipping_pct = clipping_len / (aligned_len + clipping_len)

    aligned_len = 0
    clipping_len = 0
    mismatch_len = 0
    for each_part in cigar_splitted:
        each_part_len = int(each_part[:-1])
        each_part_cate = each_part[-1]

        # get aligned_len
        if each_part_cate in {'M', 'm', 'I', 'i', 'X', 'x', '='}:
            aligned_len += each_part_len

        # get clipping_len
        if each_part_cate in ['S', 's']:
            clipping_len += each_part_len

        # get mismatch_len
        if each_part_cate in {'I', 'i', 'X', 'x', 'D', 'd'}:
            mismatch_len += each_part_len

    aligned_pct  = float("{0:.2f}".format(aligned_len * 100 / (aligned_len + clipping_len)))
    clipping_pct = float("{0:.2f}".format(clipping_len * 100 / (aligned_len + clipping_len)))
    mismatch_pct = float("{0:.2f}".format(mismatch_len * 100 / (aligned_len)))

    return aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct


def check_both_ends_clipping(cigar_splitted):

    both_ends_clipping = False
    if len(cigar_splitted) >= 3:
        if (cigar_splitted[0][-1] in ['S', 's']) and (cigar_splitted[-1][-1] in ['S', 's']):
            both_ends_clipping = True

    return both_ends_clipping


def remove_high_mismatch(sam_in, aln_len_cutoff, mismatch_cutoff, sam_out):

    sam_out_handle = open(sam_out, 'w')
    ref_len_dict = {}
    for each_read in open(sam_in):
        each_read_split = each_read.strip().split('\t')
        if each_read.startswith('@'):
            sam_out_handle.write(each_read)

            marker_id = ''
            marker_len = 0
            for each_element in each_read_split:
                if each_element.startswith('SN:'):
                    marker_id = each_element[3:]
                if each_element.startswith('LN:'):
                    marker_len = int(each_element[3:])
            ref_len_dict[marker_id] = marker_len

        else:
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])
            cigar = each_read_split[5]
            if cigar == '*':
                sam_out_handle.write(each_read)
            else:
                cigar_splitted = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitted)
                if both_ends_clp is False:
                    r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)
                    if r1_mismatch_pct <= mismatch_cutoff:

                        if r1_aligned_len >= aln_len_cutoff:

                            # check if clp in middle
                            if ('S' not in cigar) and ('s' not in cigar):
                                sam_out_handle.write(each_read)
                            else:
                                clip_in_middle = True
                                if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                                    clip_in_middle = False
                                if (cigar_splitted[-1][-1] in ['S', 's']):
                                    if (ref_pos + r1_aligned_len - 1) == ref_len_dict[ref_id]:
                                        clip_in_middle = False

                                if clip_in_middle is False:
                                    sam_out_handle.write(each_read)
    sam_out_handle.close()


def keep_best_matches_in_sam(sam_in, sam_out):
    # get read_to_cigar_dict
    read_to_cigar_dict = {}
    for each_line in open(sam_in):
        each_line_split = each_line.strip().split('\t')
        if not each_line.startswith('@'):
            read_id = each_line_split[0]
            cigar = each_line_split[5]
            if cigar != '*':
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    if read_id not in read_to_cigar_dict:
                        read_to_cigar_dict[read_id] = {cigar}
                    else:
                        read_to_cigar_dict[read_id].add(cigar)

    # get min_mismatch for each read
    read_min_mismatch_dict = {}
    for each_read in read_to_cigar_dict:
        read_mismatch_set = set()
        for each_cigar in read_to_cigar_dict[each_read]:
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(each_cigar))
            read_mismatch_set.add(mismatch_pct)
        read_min_mismatch = min(read_mismatch_set)
        read_min_mismatch_dict[each_read] = read_min_mismatch

    sam_file_best_match_handle = open(sam_out, 'w')
    for each_line in open(sam_in):
        if each_line.startswith('@'):
            sam_file_best_match_handle.write(each_line)
        else:
            each_line_split = each_line.strip().split('\t')
            read_id = each_line_split[0]
            cigar = each_line_split[5]
            if cigar == '*':
                sam_file_best_match_handle.write(each_line)
            else:
                cigar_split = cigar_splitter(cigar)
                both_ends_clp = check_both_ends_clipping(cigar_splitter(cigar))
                if both_ends_clp is False:
                    aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_split)
                    if mismatch_pct <= (read_min_mismatch_dict[read_id]):

                        sam_file_best_match_handle.write(each_line)
    sam_file_best_match_handle.close()


bowtie_parameter                    = '--xeq --local --all --no-unal -N 1 -L 30'
pwd_bowtie2_exe                     = 'bowtie2'
pwd_samtools_exe                    = 'samtools'
reads_16s                           = '/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta'
input_16s_qc_no_ext                 = '/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_MarkerMAG_wd/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_rd1_wd/input_16S/MBARC26_SILVA138_polished.QC'
reads_16s_to_16s_sam                = 'reads_16s_to_16s_sam.sam'
reads_16s_to_16s_log                = 'reads_16s_to_16s_sam.log'
reads_16s_to_16s_sam_filtered       = 'reads_16s_to_16s_sam_filtered.sam'
mag_depth_txt                       = '/Users/songweizhi/Desktop/666/MBARC26_refined_bins_50_5_depth.txt'
ref_depth_txt                       = '/Users/songweizhi/Desktop/666/ref_depth.txt'
markermag_op                        = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_linkages_by_genome.txt'
ref_16S_id_txt                      = '/Users/songweizhi/Desktop/666/ref_16S_id.txt'
num_threads                         = 12
aln_len_cutoff                      = 70
ref_genome_metadata                 = '/Users/songweizhi/Desktop/666/reference_genome_metadata.txt'


bowtie_read_to_16s_cmd  = '%s -x %s -U %s -S %s -p %s -f %s 2> %s'          % (pwd_bowtie2_exe, input_16s_qc_no_ext, reads_16s, reads_16s_to_16s_sam, num_threads, bowtie_parameter, reads_16s_to_16s_log)
print(bowtie_read_to_16s_cmd)

reads_16s_to_16s_sam                = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam.sam'
reads_16s_to_16s_sam_filtered_mis0  = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_filtered_mis0.sam'
reads_16s_to_16s_sam_filtered_mis1  = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_filtered_mis1.sam'
reads_16s_to_16s_sam_filtered_mis2  = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_filtered_mis2.sam'

# reads_16s_to_16s_sam                = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_best.sam'
# reads_16s_to_16s_sam_filtered_mis0  = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_best_filtered_mis0.sam'
# reads_16s_to_16s_sam_filtered_mis1  = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_best_filtered_mis1.sam'
# reads_16s_to_16s_sam_filtered_mis2  = '/Users/songweizhi/Desktop/666/reads_16s_to_16s_sam_best_filtered_mis2.sam'

# remove_high_mismatch(reads_16s_to_16s_sam, aln_len_cutoff, 0, reads_16s_to_16s_sam_filtered_mis0)
# remove_high_mismatch(reads_16s_to_16s_sam, aln_len_cutoff, 1, reads_16s_to_16s_sam_filtered_mis1)
# remove_high_mismatch(reads_16s_to_16s_sam, aln_len_cutoff, 2, reads_16s_to_16s_sam_filtered_mis2)


reads_to_16s_sam_filtered                   = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted.sam'
reads_to_16s_sam_filtered_mis0              = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted_mis0.sam'
reads_to_16s_sam_filtered_mis0_best_match   = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted_mis0_best_match.sam'
reads_to_16s_sam_filtered_mis1              = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted_mis1.sam'
reads_to_16s_sam_filtered_mis1_best_match   = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted_mis1_best_match.sam'
reads_to_16s_sam_filtered_mis2              = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted_mis2.sam'
reads_to_16s_sam_filtered_mis2_best_match   = '/Users/songweizhi/Desktop/666/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_input_reads_to_16S_sorted_mis2_best_match.sam'

# remove_high_mismatch(reads_to_16s_sam_filtered, aln_len_cutoff, 2, reads_to_16s_sam_filtered_mis2)
# keep_best_matches_in_sam(reads_to_16s_sam_filtered_mis2, reads_to_16s_sam_filtered_mis2_best_match)

#remove_high_mismatch(reads_to_16s_sam_filtered, aln_len_cutoff, 1, reads_to_16s_sam_filtered_mis1)
#keep_best_matches_in_sam(reads_to_16s_sam_filtered_mis1, reads_to_16s_sam_filtered_mis1_best_match)

# remove_high_mismatch(reads_to_16s_sam_filtered, aln_len_cutoff, 0, reads_to_16s_sam_filtered_mis0)
# keep_best_matches_in_sam(reads_to_16s_sam_filtered_mis0, reads_to_16s_sam_filtered_mis0_best_match)

pwd_sam_file_filtered           = '/Users/songweizhi/Desktop/666/000/reads_16s_to_16s_sam_all_filtered_mis%s_minM%s.sam'        % (0, 70)

pwd_sam_file_filtered_random    = '/Users/songweizhi/Desktop/666/000/reads_16s_to_16s_sam_all_filtered_mis%s_minM%s_random.sam' % (0, 70)


sam_basename = 'reads_16s_to_16s_sam_all_sorted'
#sam_basename = 'MBARC26_SILVA138_polished'
#sam_basename = 'MBARC26_SILVA138_polished_best'
pwd_sam_file                    = '/Users/songweizhi/Desktop/666/000/%s.sam'                                % sam_basename


min_M_len_16s                   = 70
mismatch_cutoff                 = 0
pwd_sam_mp_file                 = '/Users/songweizhi/Desktop/666/000/%s_mr_mis%s_minM%s.txt'                % (sam_basename, mismatch_cutoff, min_M_len_16s)
pwd_sam_file_filtered           = '/Users/songweizhi/Desktop/666/000/%s_filtered_mis%s_minM%s.sam'          % (sam_basename, mismatch_cutoff, min_M_len_16s)
pwd_sam_file_filtered_random    = '/Users/songweizhi/Desktop/666/000/%s_filtered_mis%s_minM%s_random.sam'   % (sam_basename, mismatch_cutoff, min_M_len_16s)


ref_to_read_dict = {}
ref_len_dict = {}
read_len_dict = {}
for each_line in open(pwd_sam_file_filtered_random):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('@'):
        mini_assembly_id = ''
        mini_assembly_len = 0
        for each_element in each_line_split:
            if each_element.startswith('SN:'):
                mini_assembly_id = each_element[3:]
            if each_element.startswith('LN:'):
                mini_assembly_len = int(each_element[3:])
        ref_len_dict[mini_assembly_id] = mini_assembly_len
    else:
        read_id = each_line_split[0]
        ref_id = each_line_split[2]
        read_len = len(each_line_split[9])
        if ref_id not in ref_to_read_dict:
            ref_to_read_dict[ref_id] = {read_id}
        else:
            ref_to_read_dict[ref_id].add(read_id)
        read_len_dict[read_id] = read_len


mag_depth_dict = {}
for each_mag_depth in open(mag_depth_txt):
    if not each_mag_depth.startswith('MAG	Length(bp)	Depth'):
        each_mag_depth_split = each_mag_depth.strip().split('\t')
        mag_depth_dict[each_mag_depth_split[0]] = float(each_mag_depth_split[2])


linked_16s_set = set()
gnm_to_linked_16s_dict = {}
for each_linkage in open(markermag_op):
    each_linkage_split = each_linkage.strip().split('\t')
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Round'):
        id_16s = each_linkage_split[0]
        id_gnm = each_linkage_split[1]
        linked_16s_set.add(id_16s)
        if id_gnm not in gnm_to_linked_16s_dict:
            gnm_to_linked_16s_dict[id_gnm] = {id_16s}
        else:
            gnm_to_linked_16s_dict[id_gnm].add(id_16s)


reported_ref_16s_copy_num_dict = {}
calculated_ref_16s_copy_num_dict = {}
ref_16s_total_depth_dict = {}
ref_depth_dict = {}
for each_ref in open(ref_genome_metadata):
    if not each_ref.startswith('Genome\t'):
        each_ref_split = each_ref.strip().split('\t')
        ref_id = each_ref_split[0]
        reported_copy_num = int(each_ref_split[1])

        if each_ref_split[2] == 'NA':
            calculated_copy_num = each_ref_split[2]
        else:
            calculated_copy_num = float(each_ref_split[2])

        ref_16s_total_depth = float(each_ref_split[3])
        ref_depth = float(each_ref_split[4])
        reported_ref_16s_copy_num_dict[ref_id] = reported_copy_num
        calculated_ref_16s_copy_num_dict[ref_id] = calculated_copy_num
        ref_16s_total_depth_dict[ref_id] = ref_16s_total_depth
        ref_depth_dict[ref_id] = ref_depth


print('Genome\tMAG\tREF\tmag_16s_depth\tref_16s_depth\tmag_depth\tref_depth')
for each_mag in gnm_to_linked_16s_dict:
    linked_16s_set = gnm_to_linked_16s_dict[each_mag]
    linked_16s_len_list = [ref_len_dict[s16] for s16 in linked_16s_set]
    linked_16s_mean_len = sum(linked_16s_len_list)/len(linked_16s_len_list)
    linked_16s_total_len = sum(linked_16s_len_list)
    current_mag_linked_16s_read_set = set()
    for each_linked_16s in linked_16s_set:
        current_16s_mapped_reads = ref_to_read_dict.get(each_linked_16s, [])
        for i in current_16s_mapped_reads:
            current_mag_linked_16s_read_set.add(i)

    current_mag_linked_16s_total_len = 0
    for read_len_id in current_mag_linked_16s_read_set:
        current_mag_linked_16s_total_len += read_len_dict[read_len_id]

    ref_depth = ref_depth_dict[each_mag]
    reported_ref_16s_copy_num = reported_ref_16s_copy_num_dict[each_mag]
    calculated_ref_16s_copy_num = calculated_ref_16s_copy_num_dict[each_mag]
    ref_16s_depth = ref_16s_total_depth_dict[each_mag]

    mag_16s_depth = current_mag_linked_16s_total_len / linked_16s_mean_len
    mag_16s_depth = float("{0:.2f}".format(mag_16s_depth))
    mag_depth = mag_depth_dict[each_mag]
    copy_num = mag_16s_depth / mag_depth
    copy_num = float("{0:.2f}".format(copy_num))

    print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (each_mag, copy_num, calculated_ref_16s_copy_num, reported_ref_16s_copy_num, mag_16s_depth, ref_16s_depth, mag_depth, ref_depth))




'''
module unload python
module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_MarkerMAG_wd/reads_16s_to_16s_sam
bowtie2 -x /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_MarkerMAG_wd/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_rd1_wd/input_16S/MBARC26_SILVA138_polished.QC -U /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta -S reads_16s_to_16s_sam.sam -p 12 -f --xeq --local --all --no-unal -N 1 -L 30 2> reads_16s_to_16s_sam.log
bowtie2 -x /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_MarkerMAG_wd/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_rd1_wd/input_16S/MBARC26_SILVA138_polished.QC -U /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta -S reads_16s_to_16s_sam_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30 2> reads_16s_to_16s_sam.log

cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26
bowtie2-build MBARC26_SILVA138_polished.fasta MBARC26_SILVA138_polished
bowtie2 -x MBARC26_SILVA138_polished -U /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta -S MBARC26_SILVA138_polished.sam -p 12 -f --xeq --local --all --no-unal -N 1 -L 30
bowtie2 -x MBARC26_SILVA138_polished -U /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta -S MBARC26_SILVA138_polished_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30

/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_polished.fasta

module unload python
module load python/3.7.3
source ~/mypython3env/bin/activate
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/reference_genomes_combined
bowtie2-build -f reference_genomes_combined.fna reference_genomes_combined
bowtie2 -x reference_genomes_combined -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S reference_genomes_combined_best_match.sam -p 16 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x reference_genomes_combined -U ../MBARC26_simulated_R1.fa,../MBARC26_simulated_R1.fa -S reference_genomes_combined_simulated_report_best.sam -p 16 -f --xeq --local --no-unal -N 1 -L 30


module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/FP
bowtie2-build -f FP.fna FP
bowtie2 -x FP -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S FP_report_best.sam -p 16 -f --xeq --local --no-unal -N 1 -L 30


module load python/3.7.3
source ~/mypython3env/bin/activate
module load samtools/1.10
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/reference_genomes_combined
BioSAK sam2bam -sam reference_genomes_combined_best_match_mis0.sam
BioSAK sam2bam -sam reference_genomes_combined_best_match_mis1.sam


module load bowtie/2.3.5.1
module load samtools/1.10
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/get_16S_sequencning_depth
bowtie2 -x reference_genomes_renamed_16S -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S reference_genomes_renamed_16S_report_all.sam -p 16 -f --xeq --local --all --no-unal -N 1 -L 30
samtools sort -n -O sam --threads 16 -o reference_genomes_renamed_16S_report_all_sorted.sam reference_genomes_renamed_16S_report_all.sam 

samtools sort -n -O sam --threads 16 -o reads_16s_to_16s_sam_all_sorted.sam reads_16s_to_16s_sam_all.sam 


samtools depth depth_old_sorted.bam > depth_old_sorted.depth

cd /Users/songweizhi/Desktop/666/000/uniq

samtools sort -O sam --threads 4 -o reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted.sam reads_16s_to_16s_sam_all_all_mis0_minM45_filtered.sam 
samtools sort -O sam --threads 4 -o reads_16s_to_16s_sam_all_all_mis1_minM45_filtered_sorted.sam reads_16s_to_16s_sam_all_all_mis1_minM45_filtered.sam 
samtools sort -O sam --threads 4 -o reads_16s_to_16s_sam_all_all_mis2_minM45_filtered_sorted.sam reads_16s_to_16s_sam_all_all_mis2_minM45_filtered.sam 

samtools depth reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted.sam > reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted_depth.txt
samtools depth reads_16s_to_16s_sam_all_all_mis1_minM45_filtered_sorted.sam > reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted_depth.txt
samtools depth reads_16s_to_16s_sam_all_all_mis2_minM45_filtered_sorted.sam > reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted_depth.txt

cd /Users/songweizhi/Desktop/666/000/uniq
samtools sort -O sam --threads 4 -o reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted.sam reads_16s_to_16s_sam_all_all_mis0_minM45_filtered.sam 
samtools depth reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted.sam > reads_16s_to_16s_sam_all_all_mis0_minM45_filtered_sorted_depth.txt



samtools depth reference_genomes_combined_best_match_mis0_sorted.bam > reference_genomes_combined_best_match_mis0_sorted_depth.txt

'''


'''

module load bowtie/2.3.5.1
module load samtools/1.10
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_0727_45_45_min1200_mismatch2_iden99_diff80_MarkerMAG_wd
bowtie2-build MBARC26_SILVA138_polished.QC.fasta MBARC26_SILVA138_polished.QC
bowtie2 -x MBARC26_SILVA138_polished.QC -U /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta -S MBARC26_SILVA138_polished.QC_report_all.sam -p 12 -f --xeq --local --all --no-unal -N 1 -L 30

samtools sort -n -O sam --threads 16 -o MBARC26_SILVA138_polished.QC_report_all_sorted_by_read.sam MBARC26_SILVA138_polished.QC_report_all.sam


module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/map_to_ref_gnm_sep
bowtie2-build -f AS.fna	AS
bowtie2-build -f CA.fna	CA
bowtie2-build -f CG.fna	CG
bowtie2-build -f CP.fna	CP
bowtie2-build -f DA.fna	DA
bowtie2-build -f DG.fna	DG
bowtie2-build -f DM.fna	DM
bowtie2-build -f EC.fna	EC
bowtie2-build -f EV.fna	EV
bowtie2-build -f FA.fna	FA
bowtie2-build -f FP.fna	FP
bowtie2-build -f FP_MAG.fna	FP_MAG
bowtie2-build -f HB.fna	HB
bowtie2-build -f HR.fna	HR
bowtie2-build -f HT.fna	HT
bowtie2-build -f MS.fna	MS
bowtie2-build -f ND.fna	ND
bowtie2-build -f NG.fna	NG
bowtie2-build -f NO.fna	NO
bowtie2-build -f OU.fna	OU
bowtie2-build -f PS.fna	PS
bowtie2-build -f SB.fna	SB
bowtie2-build -f SP.fna	SP
bowtie2-build -f SR.fna	SR
bowtie2-build -f SS.fna	SS
bowtie2-build -f TC.fna	TC
bowtie2-build -f TR.fna	TR


module load bowtie/2.3.5.1
module load samtools/1.10
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/get_16S_cp_num
bowtie2-build CAMI_Oral_138_16S_0.999.polished_min1200.qualified.fa CAMI_Oral_138_16S_0.999.polished_min1200.qualified
bowtie2 -x CAMI_Oral_138_16S_0.999.polished_min1200.qualified -U /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_16S_reads.fasta -S CAMI_Oral_138_16S_0.999.polished_min1200.qualified.sam -p 16 -f --xeq --no-unal -N 1 -L 30
python3 filter_sam.py -in CAMI_Oral_138_16S_0.999.polished_min1200.qualified.sam -out CAMI_Oral_138_16S_0.999.polished_min1200.qualified_mis2_aln50.sam -mm 2 -aln 50

bowtie2 -x linked_16s -U /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_16S_reads.fasta -S linked_16s_report_best.sam -p 16 -f --xeq --no-unal -N 1 -L 30
python3 filter_sam.py -in linked_16s_report_best.sam -out linked_16s_report_best_mis2_aln50.sam -mm 2 -aln 50


bowtie2-build DA.fna DA
bowtie2-build NO.fna NO
bowtie2-build FP.fna FP
bowtie2-build FP_16S.fa FP_16S
bowtie2-build FP_16S_with_pacbio.fa FP_16S_with_pacbio
bowtie2-build FP_16S_pacbio.fa FP_16S_pacbio
bowtie2-build FP_chr_1631000_1647000.fa FP_chr_1631000_1647000
bowtie2-build reference_genomes_with_prefix_cbd.fna reference_genomes_with_prefix_cbd

bowtie2 -x FP_16S -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_report_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x FP_16S -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_report_all.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30 --all
bowtie2 -x FP_16S_with_pacbio -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_with_pacbio_report_all.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30 --all
bowtie2 -x FP_16S_pacbio -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_pacbio_report_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x FP_16S_pacbio -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_pacbio_report_all.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30 --all
bowtie2 -x FP_chr_1631000_1647000 -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S FP_chr_1631000_1647000_report_one.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x FP_chr_1631000_1647000 -U ../MBARC26_SILVA138_id99_16S_reads.fasta -S FP_chr_1631000_1647000_report_one_2.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x FP -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S FP_report_one.sam -p 16 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x FP -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S FP_global_report_one.sam -p 16 -f --xeq --no-unal -N 1 -L 30
bowtie2 -x DA -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S DA_global_report_one.sam -p 16 -f --xeq --no-unal -N 1 -L 30
bowtie2 -x NO -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S NO_global_report_one.sam -p 16 -f --xeq --no-unal -N 1 -L 30
bowtie2 -x reference_genomes_with_prefix_cbd -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S reference_genomes_global_report_one.sam -p 16 -f --xeq --no-unal -N 1 -L 30
bowtie2 -x FP_MAG -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S FP_MAG_global_report_one.sam -p 16 -f --xeq --no-unal -N 1 -L 30

python3 filter_sam.py -in FP_16S_report_best.sam -out FP_16S_report_best_mis5.sam -mm 5 -aln 50
python3 filter_sam.py -in FP_16S_report_all.sam -out FP_16S_report_all_mis5.sam -mm 5 -aln 50
python3 filter_sam.py -in FP_16S_with_pacbio_report_all.sam -out FP_16S_with_pacbio_report_all_mis5.sam -mm 5 -aln 50
python3 filter_sam.py -in FP_16S_pacbio_report_best.sam -out FP_16S_pacbio_report_best_mis5.sam -mm 5 -aln 50
python3 filter_sam.py -in FP_16S_pacbio_report_all.sam -out FP_16S_pacbio_report_all_mis5.sam -mm 5 -aln 50
python3 filter_sam.py -in FP_chr_1631000_1647000_report_one.sam -out FP_chr_1631000_1647000_report_one_mis0.sam -mm 0 -aln 50
python3 filter_sam.py -in FP_global_report_one.sam -out FP_global_report_one_mis0.sam -mm 0 -aln 50
python3 filter_sam.py -in DA_global_report_one.sam -out DA_global_report_one_mis0.sam -mm 0 -aln 50
python3 filter_sam.py -in NO_global_report_one.sam -out NO_global_report_one_mis0.sam -mm 0 -aln 50
python3 filter_sam.py -in reference_genomes_global_report_one.sam -out reference_genomes_global_report_one_mis0.sam -mm 0 -aln 50

samtools sort -O sam --threads 12 -o FP_global_report_one_mis0_sorted.sam FP_global_report_one_mis0.sam 
samtools sort -O sam --threads 12 -o DA_global_report_one_mis0_sorted.sam DA_global_report_one_mis0.sam 
samtools sort -O sam --threads 12 -o NO_global_report_one_mis0_sorted.sam NO_global_report_one_mis0.sam 
samtools sort -O sam --threads 12 -o reference_genomes_global_report_one_mis0_sorted.sam reference_genomes_global_report_one_mis0.sam 

samtools depth -a FP_global_report_one_mis0_sorted.sam > FP_global_report_one_mis0_sorted_depth.txt
samtools depth -a DA_global_report_one_mis0_sorted.sam > DA_global_report_one_mis0_sorted_depth.txt
samtools depth -a NO_global_report_one_mis0_sorted.sam > NO_global_report_one_mis0_sorted_depth.txt
samtools depth -a reference_genomes_global_report_one_mis0_sorted.sam > reference_genomes_global_report_one_mis0_sorted_depth.txt


samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis0_50bp.sam 
samtools sort -O sam --threads 12 -o combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted.sam combined_prefixed_MAGs_global_report_one_mis0_100bp.sam 

samtools depth -a combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis0_50bp_sorted_depth.txt
samtools depth -a combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted.sam > combined_prefixed_MAGs_global_report_one_mis0_100bp_sorted_depth.txt


python3 filter_sam.py -in FP_MAG_global_report_one.sam -out FP_MAG_global_report_one_mis0.sam -mm 0 -aln 50
samtools sort -O sam --threads 12 -o FP_MAG_global_report_one_mis0_sorted.sam FP_MAG_global_report_one_mis0.sam 
samtools depth -a FP_MAG_global_report_one_mis0_sorted.sam > FP_MAG_global_report_one_mis0_sorted_depth.txt


bbmap.sh ref=FP_16S.fa in=MBARC26_SILVA138_id99_16S_reads.fasta outm=FP_16S_report_all_bbmap.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g
bbmap.sh ref=FP_16S.fa in=MBARC26_SILVA138_id99_16S_reads.fasta outm=FP_16S_report_best_bbmap.sam local=t nodisk=t ambiguous=best keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g

bowtie2-build -f FP_16S_ref.fa FP_16S_ref
bowtie2-build -f FP_16S_matam.fa FP_16S_matam
bowtie2 -x FP_16S_ref -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_ref_report_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30
bowtie2 -x FP_16S_matam -U MBARC26_R1.fasta,MBARC26_R2.fasta -S FP_16S_matam_report_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30
python3 filter_sam.py -in FP_16S_ref_report_best.sam -out FP_16S_ref_report_best_mis5.sam -mm 5 -aln 70
python3 filter_sam.py -in FP_16S_matam_report_best.sam -out FP_16S_matam_report_best_mis5.sam -mm 5 -aln 70

/srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_SILVA138_id99_Matam16S_wd/MBARC26_SILVA138_id99_16S_reads.fasta

cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/try_blast
blastn -query MBARC26_SILVA138_id99_16S_reads.fasta -db seq_16s_overall.fa -out a.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 16
blastn -query MBARC26_SILVA138_id99_16S_reads.fasta -db seq_16s_overall.fa -out a.tab -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -num_threads 16

module load java/8u201-jdk
module load bbmap/38.51
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26
bbmap.sh ref=Kelp_16S.fasta in=Kelp_R1_subset.fastq in2=Kelp_R2_subset.fastq outm=reads_to_16S_fq.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g 2> reads_to_16S_fq_bbmap_stderr.txt
bbmap.sh ref=Kelp_16S.fasta in=Kelp_R1_subset.fasta in2=Kelp_R2_subset.fasta outm=reads_to_16S_fa.sam local=t nodisk=t ambiguous=all keepnames=t saa=f trd=t silent=true threads=12 -Xmx10g 2> reads_to_16S_fa_bbmap_stderr.txt


'''