from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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


class MappingRecord:

    #  sequences store in r1_seq and r2_seq should NOT been reverse complemented

    def __init__(self):

        self.r1_seq = ''
        self.r2_seq = ''

        self.r1_refs = dict()
        self.r2_refs = dict()

        self.r1_cigar_to_flag = dict()
        self.r2_cigar_to_flag = dict()

        self.r1_longest_clp_cigar = ''
        self.r1_longest_clp_falg  = ''
        self.r2_longest_clp_cigar = ''
        self.r2_longest_clp_falg  = ''

        self.qualified_reads           = False
        self.consider_round_2          = False
        self.consider_r1_unmapped_mate = False
        self.consider_r1_clipping_part = False
        self.consider_r2_unmapped_mate = False
        self.consider_r2_clipping_part = False

        self.r1_filtered_refs = set()
        self.r2_filtered_refs = set()

        self.r1_clipping_seq = ''
        self.r2_clipping_seq = ''

        self.unmapped_r1_refs = set()
        self.unmapped_r2_refs = set()

        self.clipping_r1_refs = set()
        self.clipping_r2_refs = set()

        self.unmapped_r1_refs_with_pos = set()
        self.unmapped_r2_refs_with_pos = set()

        self.clipping_r1_refs_with_pos = set()
        self.clipping_r2_refs_with_pos = set()


def sam_flag_to_rc(flag_value):

    read_rced = 'na'
    if flag_value != '':
        binary_flag = "{0:b}".format(int(flag_value))
        binary_flag_len = len(str(binary_flag))
        binary_flag_polished = '0' * (12 - binary_flag_len) + str(binary_flag)

        if binary_flag_polished[7] == '0':
            read_rced = False
        if binary_flag_polished[7] == '1':
            read_rced = True

    return read_rced


def get_rc(seq_in):
    seq_in_rc = str(SeqRecord(Seq(seq_in)).reverse_complement().seq)
    return seq_in_rc


def get_max_clp_and_index(r1_cigar_list, r2_cigar_list):

    r1_cigar_list_split = [cigar_splitter(i) for i in r1_cigar_list]
    r2_cigar_list_split = [cigar_splitter(i) for i in r2_cigar_list]

    r1_cigar_list_split_only_clp = []
    for each_r1_cigar_split in r1_cigar_list_split:
        clp_len_l = 0
        if each_r1_cigar_split[0][-1] in ['S', 's']:
            clp_len_l = int(each_r1_cigar_split[0][:-1])
        clp_len_r = 0
        if each_r1_cigar_split[-1][-1] in ['S', 's']:
            clp_len_r = int(each_r1_cigar_split[-1][:-1])
        r1_cigar_list_split_only_clp.append([clp_len_l, clp_len_r])

    r2_cigar_list_split_only_clp = []
    for each_r2_cigar_split in r2_cigar_list_split:
        clp_len_l = 0
        if each_r2_cigar_split[0][-1] in ['S', 's']:
            clp_len_l = int(each_r2_cigar_split[0][:-1])
        clp_len_r = 0
        if each_r2_cigar_split[-1][-1] in ['S', 's']:
            clp_len_r = int(each_r2_cigar_split[-1][:-1])
        r2_cigar_list_split_only_clp.append([clp_len_l, clp_len_r])

    cigar_list_split_only_clp_r1_r2 = [r1_cigar_list_split_only_clp, r2_cigar_list_split_only_clp]

    max_value = 0
    max_value_index = ''
    for num_list_1 in cigar_list_split_only_clp_r1_r2[0]:
        if num_list_1[0] > max_value:
            max_value = num_list_1[0]
            max_value_index = 'r1_l'
        if num_list_1[1] > max_value:
            max_value = num_list_1[1]
            max_value_index = 'r1_r'
    for num_list_2 in cigar_list_split_only_clp_r1_r2[1]:
        if num_list_2[0] > max_value:
            max_value = num_list_2[0]
            max_value_index = 'r2_l'
        if num_list_2[1] > max_value:
            max_value = num_list_2[1]
            max_value_index = 'r2_r'

    # get the best cigar
    best_cigar = ''
    if max_value_index == 'r1_l':
        if best_cigar == '':
            for each_cigar in r1_cigar_list:
                if (each_cigar.startswith('%sS' % max_value)) or (each_cigar.startswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r1_r':
        if best_cigar == '':
            for each_cigar in r1_cigar_list:
                if (each_cigar.endswith('%sS' % max_value)) or (each_cigar.endswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r2_l':
        if best_cigar == '':
            for each_cigar in r2_cigar_list:
                if (each_cigar.startswith('%sS' % max_value)) or (each_cigar.startswith('%ss' % max_value)):
                    best_cigar = each_cigar
    elif max_value_index == 'r2_r':
        if best_cigar == '':
            for each_cigar in r2_cigar_list:
                if (each_cigar.endswith('%sS' % max_value)) or (each_cigar.endswith('%ss' % max_value)):
                    best_cigar = each_cigar

    return best_cigar, max_value, max_value_index


wd = '/Users/songweizhi/Desktop/tunning'
#input_reads_to_16s_sam_best_match       = '%s/GI_0503_mis1_75_45_input_reads_to_16S_best_match.sam' % wd
input_reads_to_16s_sam_best_match       = '%s/input_reads_to_16S_best_match.sam'                                                             % wd
sam_best_match_unmapped_mates_seq_file  = '%s/GI_0507_mis2_45_45_input_reads_to_16S_best_match_unmapped_mates.fa'   % wd
unmapped_mates_seq_file                 = '%s/unmapped_mates.fa'                                                    % wd
clipping_parts_seq_file                 = '%s/clipping_parts.fa'                                                    % wd
min_M_len_16s                           = 45
min_M_len_ctg                           = 45
min_M_pct                               = 35
mismatch_cutoff                         = 2
unmapped_mates_to_16s_sam               = '%s/unmapped_mates_to_16S.sam'                            % wd
read_to_marker_connector                = '___r___'


'''
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U GI_0507_mis2_45_45_unmapped_mates.fa -S unmapped_mates_to_16S_fast.sam -p 12 -f --local --all --no-unal --very-fast-local
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U GI_0507_mis2_45_45_unmapped_mates.fa -S unmapped_mates_to_16S_slow.sam -p 12 -f --local --all --no-unal
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U GI_0507_mis2_45_45_unmapped_mates.fa -S unmapped_mates_to_16S_very_slow.sam -p 12 -f --local --all --no-unal --very-sensitive-local


S11_9543565.2.fa


# can S11_9543565.2 mapped to 3_GI_subsample_5_85 with very sensitive mode?
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0508_mis2_45_45_MarkerMAG_wd/GI_0508_mis2_45_45_step_1_wd
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U S11_9543565.2.fa -S S11_9543565.2_to_16S_fast.sam -p 12 -f --local --all --no-unal --very-fast-local
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U S11_9543565.2.fa -S S11_9543565.2_to_16S_slow.sam -p 12 -f --local --all --no-unal
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U S11_9543565.2.fa -S S11_9543565.2_to_16S_very_slow.sam -p 12 -f --local --all --no-unal --very-sensitive-local

bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U S11_9543565.2.fa -S S11_9543565.2_to_16S_very_slow_N1.sam -p 12 -f --local --all --no-unal --very-sensitive-local -N 1
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U S11_9543565.2.fa -S S11_9543565.2_to_16S_slow_N1.sam -p 12 -f --local --all --no-unal -N 1


grep S11_4013158.1 unmapped_mates_to_16S_fast.sam > subset_fast.sam
grep S11_4013158.1 unmapped_mates_to_16S_slow.sam > subset_slow.sam
grep S11_4013158.1 unmapped_mates_to_16S_very_slow.sam > subset_very_slow.sam


S11_4013158.1   3/61
S11_4013158.1   3/61
S11_4013158.1   3/61
S11_4013158.1   3/61



                mismatch
S11_4013158.1   3/61
S2_12033538.1   3/68
S5_7471726.1    2/55
S10_11295486.1  1/36

Marker              read             mismatch
3_GI_subsample_5_85 S2_12033538.1    3/68        should be ignored! 
#3_GI_subsample_5_85 S11_4013158.1   3/61        should be ignored! 
#3_GI_subsample_5_85 S5_7471726.1    2/55        should be ignored! 
#3_GI_subsample_5_85 S10_11295486.1  1/36        should be ignored! 

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd
grep -E '@|S11_4013158.' GI_0507_mis2_45_45_input_reads_to_16S_best_match.sam > input_reads_to_16S_best_match.sam
grep -E '@|S11_4013158.' GI_0507_mis2_45_45_unmapped_mates_to_16s_reformat.sam > unmapped_mates_to_16S.sam

map S2_12033538.1 to 3_GI_subsample_5_85  why didn't get match?

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/0_the_wrong_one
bowtie2-build --quiet -f 3_GI_subsample_5_85.fasta 3_GI_subsample_5_85
bowtie2 -x 3_GI_subsample_5_85 -U S2_12033538.1.fa -S S2_12033538.1.sam -p 6 -f --local --all --no-unal

bowtie2 -x ref -U S2_12033538.1.fa -S S2_12033538.1_concate.sam -p 6 -f --local --all --no-unal


index_ref_cmd = '
 %s' % (ref_seq, ref_index)
bowtie2_cmd = ''
if (len(paired_base_set) > 0) and (len(unpaired_read_set) > 0):
    bowtie2_cmd   = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p 6 -f --local --all --no-unal' % (ref_index, op_r1, op_r2, op_up, op_sam)
if (len(paired_base_set) > 0) and (len(unpaired_read_set) == 0):
    bowtie2_cmd   = 'bowtie2 -x %s -1 %s -2 %s -S %s -p 6 -f --local --all --no-unal' % (ref_index, op_r1, op_r2, op_sam)
if (len(paired_base_set) == 0) and (len(unpaired_read_set) > 0):
    bowtie2_cmd   = '' % (ref_index, op_up, op_sam)
    


S4_602526.2
both mates mapped to the same 16s with very good quality, should be ignored    3_GI_subsample_75_2814


cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd
grep -E '@|S4_602526.' GI_0507_mis2_45_45_input_reads_to_16S_best_match.sam > input_reads_to_16S_best_match.sam
grep -E '@|S4_602526.' GI_0507_mis2_45_45_unmapped_mates_to_16s_reformat.sam > unmapped_mates_to_16S.sam


grep -E '@|S4_602526.' /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd/GI_0507_mis2_45_45_input_reads_to_16S_best_match.sam > sub.sam

grep -E '@|S4_602526.' GI_0507_mis2_45_45_clipping_parts_best_match.sam > sub.sam

grep 'S4_602526.2	419	3_GI_subsample_100_1396' /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd/GI_0507_mis2_45_45_input_reads_to_16S.sam

grep 'S4_602526.2	3_GI_subsample_75_2814' /srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/GI_0507_mis2_45_45_MarkerMAG_wd.1/GI_0507_mis2_45_45_step_1_wd/GI_0507_mis2_45_45_unmapped_mates_to_16s_reformat.sam


S4_8191380.2
S4_8191380.1    high mismatch, ignored, should not consider this one as unmapped mate

seqtk subseq GI_R1.fasta to_extract.txt > to_extract_R1.fa
seqtk subseq GI_R2.fasta to_extract.txt > to_extract_R2.fa
bowtie2 -x 3_GI_assembled_16S_uclust_0.999_index/3_GI_assembled_16S_uclust_0.999 -U GI_0506_mis2_45_45_unmapped_mates.fa -S unmapped_mates_to_16s.sam -p 12 -f --local --all --no-unal --very-fast-local
reformat.sh in=unmapped_mates_to_16s.sam out=s_unmapped_mates_to_16s_reformat.sam sam=1.4
grep 3_GI_subsample_25_818 s_unmapped_mates_to_16s_reformat.sam > unmapped_mates_to_16s_subset.sam


grep 3_GI_subsample_75_2802 s_unmapped_mates_to_16s_reformat.sam > unmapped_mates_to_16s_subset2.sam
grep 3_GI_subsample_75_2814 s_unmapped_mates_to_16s_reformat.sam > unmapped_mates_to_16s_subset3.sam


cat unmapped_mates_to_16s_subset.sam unmapped_mates_to_16s_subset2.sam unmapped_mates_to_16s_subset3.sam > unmapped_mates_to_16s_subset_combined.sam




grep -E '@|3_GI_subsample_25_818|3_GI_subsample_75_2802|3_GI_subsample_75_2814' GI_0506_mis2_45_45_input_reads_to_16S_best_match.sam > test3.sam

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0505_mis2_75_45_keep_short_M_MarkerMAG_wd/hc_0505_mis2_75_45_keep_short_M_step_1_wd
grep cami_hc_SILVA138_id99_75_subsample_75_2741 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_1.sam
grep cami_hc_SILVA138_id99_50_subsample_10_876 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_2.sam
grep cami_hc_SILVA138_id99_50_subsample_50_1792 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S_best_match.sam > sub_3.sam
cat sub_1.sam sub_2.sam sub_3.sam > sub.sam

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/hc_0505_mis2_75_45_keep_short_M_MarkerMAG_wd/hc_0505_mis2_75_45_keep_short_M_step_1_wd
grep cami_hc_SILVA138_id99_75_subsample_75_2741 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S.sam > sub_1_all.sam
grep cami_hc_SILVA138_id99_50_subsample_10_876 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S.sam > sub_2_all.sam
grep cami_hc_SILVA138_id99_50_subsample_50_1792 hc_0505_mis2_75_45_keep_short_M_input_reads_to_16S.sam > sub_3_all.sam
cat sub_1_all.sam sub_2_all.sam sub_3_all.sam > sub_all.sam


'''

MappingRecord_dict = {}
marker_len_dict = {}
for each_read in open(input_reads_to_16s_sam_best_match):
    each_read_split = each_read.strip().split('\t')

    if each_read.startswith('@'):
        marker_id = ''
        marker_len = 0
        for each_element in each_read_split:
            if each_element.startswith('SN:'):
                marker_id = each_element[3:]
            if each_element.startswith('LN:'):
                marker_len = int(each_element[3:])
        marker_len_dict[marker_id] = marker_len
    else:
        read_id = each_read_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        read_flag = int(each_read_split[1])
        read_seq = each_read_split[9]
        cigar = each_read_split[5]

        if cigar != '*':
            ref_id = each_read_split[2]
            ref_pos = int(each_read_split[3])
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            cigar_splitted = cigar_splitter(cigar)
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
            both_ends_clp = check_both_ends_clipping(cigar_splitted)

            if (aligned_pct >= min_M_pct) and (mismatch_pct <= mismatch_cutoff) and (both_ends_clp is False):

                # check if clp in the middle
                clip_in_middle = False
                if ('S' in cigar) or ('s' in cigar):
                    clip_in_middle = True
                    if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
                        clip_in_middle = False
                    if (cigar_splitted[-1][-1] in ['S', 's']):
                        if (ref_pos + aligned_len - 1) == marker_len_dict[ref_id]:
                            clip_in_middle = False

                # if not clp in the middle
                if clip_in_middle is False:
                    if read_id_base not in MappingRecord_dict:
                        MappingRecord_dict[read_id_base] = MappingRecord()
                    if read_strand == '1':
                        MappingRecord_dict[read_id_base].r1_refs[ref_id_with_pos] = cigar
                        MappingRecord_dict[read_id_base].r1_cigar_to_flag[cigar] = read_flag
                    if read_strand == '2':
                        MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar
                        MappingRecord_dict[read_id_base].r2_cigar_to_flag[cigar] = read_flag

        # store_read_seq into dict
        if read_id_base not in MappingRecord_dict:
            MappingRecord_dict[read_id_base] = MappingRecord()

        if read_strand == '1':
            if MappingRecord_dict[read_id_base].r1_seq == '':

                # turn back if read reverse complemented
                read_rc = sam_flag_to_rc(read_flag)
                r1_seq_to_store = read_seq
                if read_rc is True:
                    r1_seq_to_store = get_rc(read_seq)

                MappingRecord_dict[read_id_base].r1_seq = r1_seq_to_store

        if read_strand == '2':
            if MappingRecord_dict[read_id_base].r2_seq == '':

                # turn back if read reverse complemented
                read_rc = sam_flag_to_rc(read_flag)
                r2_seq_to_store = read_seq
                if read_rc is True:
                    r2_seq_to_store = get_rc(read_seq)

                MappingRecord_dict[read_id_base].r2_seq = r2_seq_to_store


# add sequences of unmapped mates to mp dict
for each_read in SeqIO.parse(sam_best_match_unmapped_mates_seq_file, 'fasta'):
    read_id = str(each_read.id)
    read_basename = '.'.join(read_id.split('.')[:-1])
    read_strand = read_id.split('.')[-1]

    if read_basename in MappingRecord_dict:
        if read_strand == '1':
            MappingRecord_dict[read_basename].r1_seq = str(each_read.seq)
        if read_strand == '2':
            MappingRecord_dict[read_basename].r2_seq = str(each_read.seq)


##################################################### parse MappingRecord_dict ####################################################

for each_mp in MappingRecord_dict:
    current_mp_record = MappingRecord_dict[each_mp]
    print('%s\tcurrent_mp_r1_refs\t%s' % (each_mp, current_mp_record.r1_refs))
    print('%s\tcurrent_mp_r2_refs\t%s' % (each_mp, current_mp_record.r2_refs))

    # only r1 mapped
    if (current_mp_record.r1_refs != {}) and (current_mp_record.r2_refs == {}):

        #print(each_mp)
        #print(current_mp_record.r1_refs)
        #print('%s\tcurrent_mp_r1_refs\t%s' % (each_mp, current_mp_record.r1_refs))
        #print('%s\tcurrent_mp_r2_refs\t%s' % (each_mp, current_mp_r2_refs))

        r1_refs_filtered_by_M_len = {}
        for r1_ref in current_mp_record.r1_refs:
            r1_ref_cigar = current_mp_record.r1_refs[r1_ref]
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r1_ref_cigar))
            if aligned_len >= min_M_len_16s:
                r1_refs_filtered_by_M_len[r1_ref] = current_mp_record.r1_refs[r1_ref]
        #print('%s\tr1_refs_filtered\t%s' % (each_mp, r1_refs_filtered_by_M_len))
        #print()
        #print('r1_refs_filtered_by_M_len: %s' % r1_refs_filtered_by_M_len)

        if len(r1_refs_filtered_by_M_len) == 0:
            current_mp_record.r1_refs = {}
        else:
            # consider unmapped mate
            current_mp_record.qualified_reads = True
            current_mp_record.consider_r1_unmapped_mate = True
            current_mp_record.r1_refs = r1_refs_filtered_by_M_len
            current_mp_record.r1_filtered_refs = {i.split('_pos_')[0] for i in r1_refs_filtered_by_M_len}

            # if full length match not found, consider as clipping mapped
            r1_cigar_list = list(r1_refs_filtered_by_M_len.values())
            full_length_match_found_r1 = False
            for each_cigar_r1 in r1_cigar_list:
                if ('S' not in each_cigar_r1) and ('s' not in each_cigar_r1):
                    full_length_match_found_r1 = True

            if full_length_match_found_r1 is False:
                best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, [])
                best_cigar_flag = current_mp_record.r1_cigar_to_flag.get(best_cigar, '')
                best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                if max_value >= min_M_len_ctg:
                    current_mp_record.consider_r1_clipping_part = True

                    current_clipping_seq = ''
                    if best_cigar_rc is False:
                        if max_value_index == 'r1_l':
                            current_clipping_seq = current_mp_record.r1_seq[:max_value]
                        if max_value_index == 'r1_r':
                            current_clipping_seq = current_mp_record.r1_seq[-max_value:]
                    else:
                        r1_seq_rc = get_rc(current_mp_record.r1_seq)
                        if max_value_index == 'r1_l':
                            current_clipping_seq_rc = r1_seq_rc[:max_value]
                            current_clipping_seq = get_rc(current_clipping_seq_rc)
                        if max_value_index == 'r1_r':
                            current_clipping_seq_rc = r1_seq_rc[-max_value:]
                            current_clipping_seq = get_rc(current_clipping_seq_rc)
                    current_mp_record.r1_clipping_seq = current_clipping_seq

    # only r2 mapped
    elif (current_mp_record.r1_refs == {}) and (current_mp_record.r2_refs != {}):

        r2_refs_filtered_by_M_len = {}
        for r2_ref in current_mp_record.r2_refs:
            r2_ref_cigar = current_mp_record.r2_refs[r2_ref]
            aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitter(r2_ref_cigar))
            if aligned_len >= min_M_len_16s:
                r2_refs_filtered_by_M_len[r2_ref] = current_mp_record.r2_refs[r2_ref]

        if len(r2_refs_filtered_by_M_len) == 0:
            current_mp_record.r2_refs = {}
        else:
            # consider as paired
            current_mp_record.qualified_reads = True
            current_mp_record.consider_r2_unmapped_mate = True
            current_mp_record.r2_refs = r2_refs_filtered_by_M_len
            current_mp_record.r2_filtered_refs =  {i.split('_pos_')[0] for i in r2_refs_filtered_by_M_len}

            # if full length match not found, consider as clipping mapped
            r2_cigar_list = list(r2_refs_filtered_by_M_len.values())
            full_length_match_found_r2 = False
            for each_cigar_r2 in r2_cigar_list:
                if ('S' not in each_cigar_r2) and ('s' not in each_cigar_r2):
                    full_length_match_found_r2 = True

            if full_length_match_found_r2 is False:

                best_cigar, max_value, max_value_index = get_max_clp_and_index([], r2_cigar_list)
                best_cigar_flag = current_mp_record.r2_cigar_to_flag.get(best_cigar, '')
                best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

                if max_value >= min_M_len_ctg:
                    current_mp_record.consider_r2_clipping_part = True

                    current_clipping_seq = ''
                    if best_cigar_rc is False:
                        if max_value_index == 'r2_l':
                            current_clipping_seq = current_mp_record.r2_seq[:max_value]
                        if max_value_index == 'r2_r':
                            current_clipping_seq = current_mp_record.r2_seq[-max_value:]
                    else:
                        r2_seq_rc = get_rc(current_mp_record.r2_seq)
                        if max_value_index == 'r2_l':
                            current_clipping_seq_rc = r2_seq_rc[:max_value]
                            current_clipping_seq = get_rc(current_clipping_seq_rc)
                        if max_value_index == 'r2_r':
                            current_clipping_seq_rc = r2_seq_rc[-max_value:]
                            current_clipping_seq = get_rc(current_clipping_seq_rc)

                    current_mp_record.r2_clipping_seq = current_clipping_seq

    # both of r1 and r2 mapped
    elif (current_mp_record.r1_refs != {}) and (current_mp_record.r2_refs != {}):

        r1_cigar_list = list(current_mp_record.r1_refs.values())
        r2_cigar_list = list(current_mp_record.r2_refs.values())

        r1_cigar_all_clp = True
        for each_cigar_r1 in r1_cigar_list:
            if ('S' not in each_cigar_r1) and ('s' not in each_cigar_r1):
                r1_cigar_all_clp = False

        r2_cigar_all_clp = True
        for each_cigar_r2 in r2_cigar_list:
            if ('S' not in each_cigar_r2) and ('s' not in each_cigar_r2):
                r2_cigar_all_clp = False

        best_cigar = ''
        max_value = 0
        max_value_index = ''
        if (r1_cigar_all_clp is True) and (r2_cigar_all_clp is False):
            best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, [])
        if (r1_cigar_all_clp is False) and (r2_cigar_all_clp is True):
            best_cigar, max_value, max_value_index = get_max_clp_and_index([], r2_cigar_list)
        if (r1_cigar_all_clp is True) and (r2_cigar_all_clp is True):
            best_cigar, max_value, max_value_index = get_max_clp_and_index(r1_cigar_list, r2_cigar_list)

        if max_value >= min_M_len_ctg:

            best_cigar_flag = ''
            if max_value_index in ['r1_l', 'r1_r']:
                best_cigar_flag = current_mp_record.r1_cigar_to_flag.get(best_cigar, '')
            if max_value_index in ['r2_l', 'r2_r']:
                best_cigar_flag = current_mp_record.r2_cigar_to_flag.get(best_cigar, '')
            best_cigar_rc = sam_flag_to_rc(best_cigar_flag)

            current_mp_r1_refs_no_pos = {i.split('_pos_')[0] for i in current_mp_record.r1_refs}
            current_mp_r2_refs_no_pos = {i.split('_pos_')[0] for i in current_mp_record.r2_refs}
            current_mp_shared_refs_no_pos = current_mp_r1_refs_no_pos.intersection(current_mp_r2_refs_no_pos)

            if len(current_mp_shared_refs_no_pos) > 0:
                current_mp_record.qualified_reads = True

                # for clipping part, only consider shared refs
                current_clipping_seq = ''
                if best_cigar_rc is False:

                    if max_value_index == 'r1_l':
                        current_clipping_seq = current_mp_record.r1_seq[:max_value]

                    if max_value_index == 'r1_r':
                        current_clipping_seq = current_mp_record.r1_seq[-max_value:]

                    if max_value_index == 'r2_l':
                        current_clipping_seq = current_mp_record.r2_seq[:max_value]

                    if max_value_index == 'r2_r':
                        current_clipping_seq = current_mp_record.r2_seq[-max_value:]
                else:
                    r1_seq_rc = get_rc(current_mp_record.r1_seq)
                    r2_seq_rc = get_rc(current_mp_record.r2_seq)

                    if max_value_index == 'r1_l':
                        current_clipping_seq_rc = r1_seq_rc[:max_value]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if max_value_index == 'r1_r':
                        current_clipping_seq_rc = r1_seq_rc[-max_value:]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if max_value_index == 'r2_l':
                        current_clipping_seq_rc = r2_seq_rc[:max_value]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                    if max_value_index == 'r2_r':
                        current_clipping_seq_rc = r2_seq_rc[-max_value:]
                        current_clipping_seq = get_rc(current_clipping_seq_rc)

                if max_value_index in ['r1_l', 'r1_r']:
                    current_mp_record.consider_r1_clipping_part = True
                    current_mp_record.r1_clipping_seq = current_clipping_seq
                    current_mp_record.r1_filtered_refs = current_mp_shared_refs_no_pos
                if max_value_index in ['r2_l', 'r2_r']:
                    current_mp_record.consider_r2_clipping_part = True
                    current_mp_record.r2_clipping_seq = current_clipping_seq
                    current_mp_record.r2_filtered_refs = current_mp_shared_refs_no_pos



# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    if MappingRecord_dict[each_mp].qualified_reads is False:
        MappingRecord_dict.pop(each_mp)


# write out sequences
unmapped_mates_handle = open(unmapped_mates_seq_file, 'w')
clipping_part_seq_handle = open(clipping_parts_seq_file, 'w')
for qualified_read in MappingRecord_dict:

    r1_name = '%s.1' % qualified_read
    r2_name = '%s.2' % qualified_read
    read_mr = MappingRecord_dict[qualified_read]

    if read_mr.consider_r1_unmapped_mate is True:
        unmapped_mates_handle.write('>%s\n' % r2_name)
        unmapped_mates_handle.write('%s\n' % read_mr.r2_seq)

    if read_mr.consider_r2_unmapped_mate is True:
        unmapped_mates_handle.write('>%s\n' % r1_name)
        unmapped_mates_handle.write('%s\n' % read_mr.r1_seq)

    if read_mr.consider_r1_clipping_part is True:
        clipping_part_seq_handle.write('>%s\n' % r1_name)
        clipping_part_seq_handle.write('%s\n' % read_mr.r1_clipping_seq)

    if read_mr.consider_r2_clipping_part is True:
        clipping_part_seq_handle.write('>%s\n' % r2_name)
        clipping_part_seq_handle.write('%s\n' % read_mr.r2_clipping_seq)

    # if qualified_read == 'S4_8191380':
    #     print('to extract:')
    #     if read_mr.qualified_reads is True:
    #         #if len(read_mr.r1_filtered_refs) > 0:
    #         print('%s.1\t%s' % (qualified_read, read_mr.r1_filtered_refs))
    #         #if len(read_mr.r2_filtered_refs) > 0:
    #         print('%s.2\t%s' % (qualified_read, read_mr.r2_filtered_refs))
unmapped_mates_handle.close()
clipping_part_seq_handle.close()


# read in unmapped_mates_to_16s mapping results
unmapped_mates_to_16s_dict = {}
for each_match in open(unmapped_mates_to_16s_sam):
    if not each_match.startswith('@'):
        each_match_split = each_match.strip().split('\t')
        read_id = each_match_split[0]
        ref_id = each_match_split[2]
        cigar = each_match_split[5]
        read_to_ref_key = '%s%s%s' % (read_id, read_to_marker_connector, ref_id)
        if read_to_ref_key not in unmapped_mates_to_16s_dict:
            unmapped_mates_to_16s_dict[read_to_ref_key] = {cigar}
        else:
            unmapped_mates_to_16s_dict[read_to_ref_key].add(cigar)

# for each_read_to_ref_key in unmapped_mates_to_16s_dict:
#     if len(unmapped_mates_to_16s_dict[each_read_to_ref_key]) > 1:
#         print(unmapped_mates_to_16s_dict[each_read_to_ref_key])


unmapped_mates_to_remove = set()
for each_mp in MappingRecord_dict.copy():
    mp_record = MappingRecord_dict[each_mp]

    if mp_record.consider_r1_unmapped_mate is True:

        if each_mp == 'S4_602526':
            print('%s\tr1_refs\t%s'             % (each_mp, mp_record.r1_refs))
            print('%s\tr1_filtered_refs\t%s'    % (each_mp, mp_record.r1_filtered_refs))
            print('%s\tr2_refs\t%s'             % (each_mp, mp_record.r2_refs))
            print('%s\tr2_filtered_refs\t%s'    % (each_mp, mp_record.r2_filtered_refs))

        for r1_ref_with_pos in mp_record.r1_refs:
            r1_ref_no_pos = r1_ref_with_pos.split('_pos_')[0]
            r2_to_ref_key     = '%s.2%s%s' % (each_mp, read_to_marker_connector, r1_ref_no_pos)
            r2_to_ref_matches = unmapped_mates_to_16s_dict.get(r2_to_ref_key, [])

            # if each_mp == 'S9_7256286':
            #     print(r1_to_ref_matches)

            if len(r2_to_ref_matches) == 1:
                cigar_splitted = cigar_splitter([i for i in r2_to_ref_matches][0])
                both_ends_clp = check_both_ends_clipping(cigar_splitted)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                #if (mismatch_pct > mismatch_cutoff) or (both_ends_clp is True) or ((mismatch_pct <= mismatch_cutoff) and (clipping_len <= min_M_len_ctg)):
                if (mismatch_pct > mismatch_cutoff) or (both_ends_clp is True):
                    # MappingRecord_dict[each_mp].r1_refs.pop(r1_ref_with_pos)
                    MappingRecord_dict[each_mp].r1_filtered_refs.remove(r1_ref_no_pos)

        if len(MappingRecord_dict[each_mp].r1_filtered_refs) == 0:
            MappingRecord_dict[each_mp].consider_r1_unmapped_mate = False
            unmapped_mates_to_remove.add('%s.2' % each_mp)


    if mp_record.consider_r2_unmapped_mate is True:

        # if each_mp == 'S9_7256286':
        #     print('%s\tr1_refs\t%s'             % (each_mp, mp_record.r1_refs))
        #     print('%s\tr1_filtered_refs\t%s'    % (each_mp, mp_record.r1_filtered_refs))
        #     print('%s\tr2_refs\t%s'             % (each_mp, mp_record.r2_refs))
        #     print('%s\tr2_filtered_refs\t%s'    % (each_mp, mp_record.r2_filtered_refs))

        for r2_ref_with_pos in mp_record.r2_refs:
            r2_ref_no_pos = r2_ref_with_pos.split('_pos_')[0]
            r1_to_ref_key     = '%s.1%s%s' % (each_mp, read_to_marker_connector, r2_ref_no_pos)
            r1_to_ref_matches = unmapped_mates_to_16s_dict.get(r1_to_ref_key, [])

            # if each_mp == 'S9_7256286':
            #     print(r1_to_ref_matches)

            if len(r1_to_ref_matches) == 1:
                cigar_splitted = cigar_splitter([i for i in r1_to_ref_matches][0])
                both_ends_clp = check_both_ends_clipping(cigar_splitted)
                aligned_len, aligned_pct, clipping_len, clipping_pct, mismatch_pct = get_cigar_stats(cigar_splitted)
                #if (mismatch_pct > mismatch_cutoff) or (both_ends_clp is True) or ((mismatch_pct <= mismatch_cutoff) and (clipping_len <= min_M_len_ctg)):
                if (mismatch_pct > mismatch_cutoff) or (both_ends_clp is True):
                    MappingRecord_dict[each_mp].r2_filtered_refs.remove(r2_ref_no_pos)

        if len(MappingRecord_dict[each_mp].r2_filtered_refs) == 0:
            MappingRecord_dict[each_mp].consider_r2_unmapped_mate = False
            unmapped_mates_to_remove.add('%s.1' % each_mp)

print(len(MappingRecord_dict))

# remove unqualified mapping record from dict
for each_mp in MappingRecord_dict.copy():
    consider_r1_unmapped_mate = MappingRecord_dict[each_mp].consider_r1_unmapped_mate
    consider_r1_clipping_part = MappingRecord_dict[each_mp].consider_r1_clipping_part
    consider_r2_unmapped_mate = MappingRecord_dict[each_mp].consider_r2_unmapped_mate
    consider_r2_clipping_part = MappingRecord_dict[each_mp].consider_r2_clipping_part
    if (consider_r1_unmapped_mate is False) and (consider_r1_clipping_part is False) and (consider_r2_unmapped_mate is False) and (consider_r2_clipping_part is False):
        MappingRecord_dict.pop(each_mp)

print(len(MappingRecord_dict))


'''
S9_7256286.1

S4_8191380.2
S4_8191380.1    high mismatch, ignored, should not consider this one as unmapped mate
'''

print('\nunmapped_mates_to_remove (%s)' % len(unmapped_mates_to_remove))
print(unmapped_mates_to_remove)

# print('%s\tr1_refs\t%s' % ('S4_602526', MappingRecord_dict['S4_602526'].r1_refs))
# print('%s\tr1_filtered_refs\t%s' % ('S4_602526', MappingRecord_dict['S4_602526'].r1_filtered_refs))
# print('%s\tr2_refs\t%s' % ('S4_602526', MappingRecord_dict['S4_602526'].r2_refs))
# print('%s\tr2_filtered_refs\t%s' % ('S4_602526', MappingRecord_dict['S4_602526'].r2_filtered_refs))
#
