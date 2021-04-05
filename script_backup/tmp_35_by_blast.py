import os


wd = '/Users/songweizhi/Desktop/round_2_by_blast'
round_2_free_living_ctg             = '%s/round2_free_living_ctg.fa'    % wd
round_2_free_living_16s             = '%s/round2_free_living_16s.fa'    % wd
round_2_free_living_blast_result    = '%s/free_living_reads_blastn.tab' % wd
round_2_align_len_cutoff = 50
round_2_iden_cutoff = 99.5
round_2_free_living_16s_to_ctg_connector = '__f__'
round2_free_living_16s_ref_file = '%s/round2_free_living_16s_refs.txt' % wd
round2_free_living_ctg_ref_file = '%s/round2_free_living_ctg_refs.txt' % wd
gnm_ctg_connector = '___'

stats_GapFilling_file = '%s/stats_GapFilling.txt' % wd


blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads %s' % 4
makeblastdb_cmd  = 'makeblastdb -in %s -dbtype nucl -parse_seqids' % round_2_free_living_ctg
blastn_cmd       = 'blastn -query %s -db %s -out %s %s' % (round_2_free_living_16s, round_2_free_living_ctg, round_2_free_living_blast_result, blast_parameters)
#os.system(makeblastdb_cmd)
#os.system(blastn_cmd)

'''
cat round2_free_living_16s_R1.fa round2_free_living_16s_R2.fa round2_free_living_16s_UP.fa > round2_free_living_16s.fa
cat round2_free_living_ctg_R1.fa round2_free_living_ctg_R2.fa round2_free_living_ctg_UP.fa > round2_free_living_ctg.fa
'''


round2_free_living_16s_ref_dict = {}
for free_living_read_16s in open(round2_free_living_16s_ref_file):
    free_living_read_16s_split = free_living_read_16s.strip().split('\t')
    read_16s_id = free_living_read_16s_split[0]
    read_16s_refs = free_living_read_16s_split[1].split(',')
    round2_free_living_16s_ref_dict[read_16s_id] = read_16s_refs

round2_free_living_ctg_ref_dict = {}
for free_living_read_ctg in open(round2_free_living_ctg_ref_file):
    free_living_read_ctg_split = free_living_read_ctg.strip().split('\t')
    read_ctg_id = free_living_read_ctg_split[0]
    read_ctg_refs = free_living_read_ctg_split[1].split(',')
    round2_free_living_ctg_ref_dict[read_ctg_id] = read_ctg_refs

round_2_min_aln_len = 50
round_2_min_iden = 99.9
round_2_min_cov = 35
m = 0
free_living_16s_to_ctg_linkage_dict = {}
free_living_16s_to_gnm_linkage_dict = {}
for each_hit in open(round_2_free_living_blast_result):
    match_split = each_hit.strip().split('\t')
    query       = match_split[0]
    subject     = match_split[1]
    identity    = float(match_split[2])
    align_len   = int(match_split[3])
    query_len   = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q  = float(align_len) * 100 / float(query_len)
    coverage_s  = float(align_len) * 100 / float(subject_len)
    qstart      = int(match_split[6])
    qend        = int(match_split[7])
    sstart      = int(match_split[8])
    send        = int(match_split[9])

    if (align_len >= round_2_min_aln_len) and (identity >= round_2_min_iden) and (coverage_q >= round_2_min_cov) and (coverage_s >= round_2_min_cov):

        if ((1 in [qstart, qend]) or (query_len in [qstart, qend])) and ((1 in [sstart, send]) or (subject_len in [sstart, send])):

            print(each_hit.strip())
            print(qstart)
            print(qend)
            print(sstart)
            print(send)
            print()
            m += 1


        query_16s_refs = round2_free_living_16s_ref_dict.get(query, [])
        subject_ctg_refs = round2_free_living_ctg_ref_dict.get(subject, [])
        for each_query_ref in query_16s_refs:
            for each_subject_ref in subject_ctg_refs:

                if each_subject_ref[-2:] in ['_l', '_r']:
                    each_subject_ref = each_subject_ref[:-2]

                subject_ref_gnm = each_subject_ref.split(gnm_ctg_connector)[0]

                q_ref_to_s_ref_key = '%s%s%s' % (each_query_ref, round_2_free_living_16s_to_ctg_connector, each_subject_ref)
                q_ref_to_s_ref_gnm_key = '%s%s%s' % (each_query_ref, round_2_free_living_16s_to_ctg_connector, subject_ref_gnm)

                if q_ref_to_s_ref_key not in free_living_16s_to_ctg_linkage_dict:
                    free_living_16s_to_ctg_linkage_dict[q_ref_to_s_ref_key] = 1
                else:
                    free_living_16s_to_ctg_linkage_dict[q_ref_to_s_ref_key] += 1

                if q_ref_to_s_ref_gnm_key not in free_living_16s_to_gnm_linkage_dict:
                    free_living_16s_to_gnm_linkage_dict[q_ref_to_s_ref_gnm_key] = 1
                else:
                    free_living_16s_to_gnm_linkage_dict[q_ref_to_s_ref_gnm_key] += 1

print(m)