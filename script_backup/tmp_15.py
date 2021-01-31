import os

genomic_seq_type = 'mag'

############################# run blast between extracted reads and metagenomic assemblies #############################

# run blastn
makeblastdb_cmd     = '%s -in %s -dbtype nucl -parse_seqids -logfile /dev/null' % (pwd_makeblastdb_exe, blast_db)
blastn_cmd_paired   = '%s -query %s -db %s -out %s %s'                          % (pwd_blastn_exe, unmapped_paired_reads_file, blast_db, unmapped_paired_reads_blastn, blast_parameters)
blastn_cmd_clipping = '%s -query %s -db %s -out %s %s'                          % (pwd_blastn_exe, clipping_reads_not_matched_part_seq, blast_db, clipping_reads_not_matched_part_seq_blastn, blast_parameters)

os.system(makeblastdb_cmd)
os.system(blastn_cmd_paired)


######################################### parse blast results for paired reads #########################################


# filter blast results for paired reads
unmapped_paired_reads_to_ctg_dict = paired_blast_results_to_dict(unmapped_paired_reads_blastn, reads_iden_cutoff, reads_cov_cutoff)

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

