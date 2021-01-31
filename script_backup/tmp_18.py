
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


def paired_blast_results_to_dict(blastn_results, iden_cutoff, query_cov_cutoff):

    query_to_subject_list_dict = {}
    for blast_hit in open(blastn_results):
        blast_hit_split = blast_hit.strip().split('\t')
        query = blast_hit_split[0]
        subject = blast_hit_split[1]
        subject_with_prefix = 'GenomicSeq__%s' % subject
        iden = float(blast_hit_split[2])
        align_len = int(blast_hit_split[3])
        query_len = int(blast_hit_split[12])
        subject_len = int(blast_hit_split[13])
        coverage_q = float(align_len) * 100 / float(query_len)

        # for perfect hits
        if (iden >= iden_cutoff) and (coverage_q == 100):
            if query not in query_to_subject_list_dict:
                query_to_subject_list_dict[query] = [subject_with_prefix]
            else:
                query_to_subject_list_dict[query].append(subject_with_prefix)

        # for nearly perfect hits
        elif (iden >= iden_cutoff) and (query_cov_cutoff <= coverage_q < 100):
            s_l = sorted([int(blast_hit_split[8]), int(blast_hit_split[9])])[0]
            s_r = sorted([int(blast_hit_split[8]), int(blast_hit_split[9])])[1]
            subject_min_gap = min([s_l, (subject_len - s_r)])
            if subject_min_gap <= 5:
                if query not in query_to_subject_list_dict:
                    query_to_subject_list_dict[query] = [subject_with_prefix]
                else:
                    query_to_subject_list_dict[query].append(subject_with_prefix)

    return query_to_subject_list_dict


########################################################################################################################

pwd_samfile                     = '/Users/songweizhi/Desktop/tuning_wd/BH_ER_050417_assembled_16S_uclust_0.995.sam'
unmapped_paired_reads_blastn    = '/Users/songweizhi/Desktop/tuning_wd/unmapped_paired_reads_blast.txt'
paired_reads_match_profile      = '/Users/songweizhi/Desktop/tuning_wd/match_profile_paired.txt'

reads_iden_cutoff               = 100
reads_cov_cutoff                = 90
perfect_match_min_cigar_M_pct   = 90
perfect_match_max_cigar_S_pct   = 10
genomic_seq_type                = 'mag'


pwd_samfile_to_mag                      = '/Users/songweizhi/Desktop/tuning_wd/BH_ER_050417_refined_bins_combined.sam'
paired_reads_match_profile_by_mapping   = '/Users/songweizhi/Desktop/tuning_wd/match_profile_paired_by_mapping.txt'


##################################################### extract reads ####################################################


# export clipping mapped reads and perfectly mapped reads
all_mapped_reads_set = set()
clipping_mapped_reads_list = set()
clipping_reads_mapped_part_dict = {}
perfectly_mapped_reads_dict = {}
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
        read_id_with_ref_pos = '%s__x__%s__x__%s' % (read_id, ref_id, ref_pos)
        all_mapped_reads_set.add(read_id)

        # for perfectly mapped reads
        if ('M' in cigar) and (len(cigar_splitted) == 1):
            if read_id_base not in perfectly_mapped_reads_dict:
                perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
            else:
                if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                    perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                else:
                    perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

        if ('S' in cigar) and (len(cigar_splitted) == 2):
            cigar_M_len = 0
            cigar_S_len = 0
            if cigar_splitted[0][-1] == 'M':
                cigar_M_len = int(cigar_splitted[0][:-1])
                cigar_S_len = int(cigar_splitted[1][:-1])
            if cigar_splitted[1][-1] == 'M':
                cigar_M_len = int(cigar_splitted[1][:-1])
                cigar_S_len = int(cigar_splitted[0][:-1])

            cigar_M_pct = cigar_M_len * 100 / (cigar_M_len + cigar_S_len)
            cigar_S_pct = cigar_S_len * 100 / (cigar_M_len + cigar_S_len)

            # for clipping reads with unmapped part >= min_cigar_S
            if (cigar_M_pct >= 30) and (cigar_S_pct >= 30):
                read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                if cigar_splitted[0][-1] == 'M':

                    # write out the sequence of unmapped part

                    # store the match info of mapped part
                    if ('%s_l' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                if cigar_splitted[1][-1] == 'M':

                    # write out the sequence of unmapped part

                    # store the match info of mapped part
                    if ('%s_r' % read_id_with_ref_pos) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id_with_ref_pos)].append(ref_id_with_prefix)

                clipping_mapped_reads_list.add(read_id)

            # for clipping reads with unmapped part < min_cigar_S
            # treat these reads as perfectly mapped reads
            elif (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (cigar_S_pct < perfect_match_max_cigar_S_pct):
                if read_id_base not in perfectly_mapped_reads_dict:
                    perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
                else:
                    if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                        perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                    else:
                        perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

        elif ('S' in cigar) and (len(cigar_splitted) > 2):

            treat_as_full_match = False
            if ((cigar_splitted[0][-1] == 'M') and (cigar_splitted[-1][-1] == 'S')) or (
                    (cigar_splitted[0][-1] == 'S') and (cigar_splitted[-1][-1] == 'M')):

                if len(cigar_splitted) == 4:

                    if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (
                            cigar_splitted[3][-1] == 'S'):
                        mismatch_ratio = int(cigar_splitted[1][:-1]) * 100 / len(read_seq)
                        if mismatch_ratio <= 1:
                            cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1])) * 100 / len(
                                read_seq)
                            cigar_S_pct = (int(cigar_splitted[3][:-1])) * 100 / len(read_seq)
                            if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                                    cigar_S_pct < perfect_match_max_cigar_S_pct):
                                treat_as_full_match = True

                    if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (
                            cigar_splitted[3][-1] == 'M'):
                        mismatch_ratio = int(cigar_splitted[2][:-1]) * 100 / len(read_seq)
                        if mismatch_ratio <= 1:
                            cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1])) * 100 / len(
                                read_seq)
                            cigar_S_pct = (int(cigar_splitted[0][:-1])) * 100 / len(read_seq)
                            if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                                    cigar_S_pct < perfect_match_max_cigar_S_pct):
                                treat_as_full_match = True

                elif len(cigar_splitted) == 6:

                    if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (
                            cigar_splitted[4][-1] == 'M') and (cigar_splitted[5][-1] == 'S'):
                        mismatch_ratio = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1])) * 100 / len(
                            read_seq)
                        if mismatch_ratio <= 1:
                            cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1]) + int(
                                cigar_splitted[4][:-1])) * 100 / len(read_seq)
                            cigar_S_pct = (int(cigar_splitted[5][:-1])) * 100 / len(read_seq)
                            if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                                    cigar_S_pct < perfect_match_max_cigar_S_pct):
                                treat_as_full_match = True

                    if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (
                            cigar_splitted[3][-1] == 'M') and (cigar_splitted[5][-1] == 'M'):
                        mismatch_ratio = (int(cigar_splitted[2][:-1]) + int(cigar_splitted[4][:-1])) * 100 / len(
                            read_seq)
                        if mismatch_ratio <= 1:
                            cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1]) + int(
                                cigar_splitted[5][:-1])) * 100 / len(read_seq)
                            cigar_S_pct = (int(cigar_splitted[0][:-1])) * 100 / len(read_seq)
                            if (cigar_M_pct >= perfect_match_min_cigar_M_pct) and (
                                    cigar_S_pct < perfect_match_max_cigar_S_pct):
                                treat_as_full_match = True

            if treat_as_full_match is True:
                if read_id_base not in perfectly_mapped_reads_dict:
                    perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
                else:
                    if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                        perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                    else:
                        perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

perfectly_mapped_read_singleton_dict = {}
for perfectly_mapped_read in perfectly_mapped_reads_dict:
    current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
    if len(current_value) == 1:
        perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value

######################################### parse blast results for paired reads #########################################

# filter blast results for paired reads
unmapped_paired_reads_to_ctg_dict = paired_blast_results_to_dict(unmapped_paired_reads_blastn, reads_iden_cutoff, reads_cov_cutoff)

def paired_blast_results_to_dict_by_mapping(unmapped_paired_reads_mapping_results):

    unmapped_paired_reads_to_ctg_dict_by_mapping = {}

    ref_len_dict = {}
    for unmapped_read in open(unmapped_paired_reads_mapping_results):

        # get ref len dict
        if unmapped_read.startswith('@'):
            unmapped_read_split = unmapped_read.strip().split('\t')
            ref_id = ''
            ref_len = 0
            for each_element in unmapped_read_split:
                if each_element.startswith('SN:'):
                    ref_id = each_element[3:]
                if each_element.startswith('LN:'):
                    ref_len = int(each_element[3:])
            ref_len_dict[ref_id] = ref_len

        else:
            qualified_unmapped_read = False
            unmapped_read_split = unmapped_read.strip().split('\t')
            read_id = unmapped_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = unmapped_read_split[2]
            ref_len = ref_len_dict[ref_id]
            ref_id_with_prefix = 'GenomicSeq__%s' % unmapped_read_split[2]
            ref_pos = int(unmapped_read_split[3])
            cigar = unmapped_read_split[5]
            read_seq = unmapped_read_split[9]
            read_len = len(read_seq)
            cigar_splitted = cigar_splitter(cigar)

            # e.g. 189M
            if ('M' in cigar) and (len(cigar_splitted) == 1):
                qualified_unmapped_read = True

            elif len(cigar_splitted) == 2:

                # e.g. 139S61M
                if (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M'):
                    cigar_M_pct = int(cigar_splitted[1][:-1]) * 100 / read_len
                    if (ref_pos == 1) and (cigar_M_pct >= 25):
                        qualified_unmapped_read = True

                # e.g. 147M53S
                if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[1][-1] == 'S'):
                    cigar_M_pct = int(cigar_splitted[0][:-1]) * 100 / read_len
                    matched_to_bp = ref_pos + int(cigar_splitted[0][:-1]) - 1
                    if (matched_to_bp == ref_len) and (cigar_M_pct >= 25):
                        qualified_unmapped_read = True

            elif len(cigar_splitted) == 3:

                # e.g. 121M1D66M, 121M1I66M
                if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M'):
                    mismatch_pct = int(cigar_splitted[1][:-1])*100/read_len
                    if mismatch_pct <= 1:
                        qualified_unmapped_read = True

                # e.g. 11S81M3S, not sure, not considered yet
                elif cigar_splitted[1][-1] == 'M':
                    min_mismatch = min(int(cigar_splitted[0][:-1]), int(cigar_splitted[2][:-1]))
                    if min_mismatch <= 3:
                        pass
            elif len(cigar_splitted) == 4:

                # e.g. ['181M', '1D', '8M', '11S'], ['181M', '2D', '8M', '11S']
                if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M') and (cigar_splitted[3][-1] == 'S'):
                    cigar_M_pct = (int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1]))*100/read_len
                    if ((ref_pos + int(cigar_splitted[0][:-1]) + int(cigar_splitted[2][:-1])) == ref_len) and (int(cigar_splitted[1][:-1]) <= 2):
                        if cigar_M_pct >= 25:
                            qualified_unmapped_read = True

                # e.g. ['24S', '95M', '1D', '27M'] and ref_pos = 1
                elif (cigar_splitted[0][-1] == 'S') and (cigar_splitted[1][-1] == 'M') and (cigar_splitted[3][-1] == 'M'):
                    cigar_M_pct = (int(cigar_splitted[1][:-1]) + int(cigar_splitted[3][:-1]))*100/read_len
                    if (ref_pos == 1) and (int(cigar_splitted[2][:-1]) <= 1) and (cigar_M_pct >= 25):
                        qualified_unmapped_read = True

            # add to dict
            if qualified_unmapped_read is True:
                if read_id not in unmapped_paired_reads_to_ctg_dict_by_mapping:
                    unmapped_paired_reads_to_ctg_dict_by_mapping[read_id] = [ref_id_with_prefix]
                else:
                    unmapped_paired_reads_to_ctg_dict_by_mapping[read_id].append(ref_id_with_prefix)


    return unmapped_paired_reads_to_ctg_dict_by_mapping


unmapped_paired_reads_to_ctg_dict_by_mapping = paired_blast_results_to_dict_by_mapping(pwd_samfile_to_mag)

print('unmapped_paired_reads_to_ctg_dict')
print(unmapped_paired_reads_to_ctg_dict)
print(len(unmapped_paired_reads_to_ctg_dict))
print()
print('unmapped_paired_reads_to_ctg_dict_by_mapping')
print(unmapped_paired_reads_to_ctg_dict_by_mapping)
print(len(unmapped_paired_reads_to_ctg_dict_by_mapping))




paired_stats_dict_num = {}
paired_reads_match_profile_handle = open(paired_reads_match_profile, 'w')
paired_reads_match_profile_handle.write('ID\tR1\tR2\n')
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


'''
cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp/BH_ER_050417_Mira_MarkerMAG_wd/BH_ER_050417_Mira_step_1_wd/BH_ER_050417_refined_bins_db
bowtie2-build -f BH_ER_050417_refined_bins_combined.fa BH_ER_050417_refined_bins_combined --quiet --threads 12
bowtie2 -x BH_ER_050417_refined_bins_combined -U ../unmapped_paired_reads.fasta -S BH_ER_050417_refined_bins_combined.sam -f --local --no-unal --quiet --threads 12


178
205
209
210

'''