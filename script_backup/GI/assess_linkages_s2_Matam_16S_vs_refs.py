import os
import glob


wd                  = '/Users/songweizhi/Desktop/ttt_mag'
blast_results_16s   = '%s/16S_vs_ref.tab'                   % wd
iden_cutoff         = 99.5
aln_len_cutoff      = 300
s16_to_ref_txt      = '%s/stats_16s_to_ref.txt'             % wd


s16_to_ref_dict = {}
for match in open(blast_results_16s):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_genome = subject.split('_16S_')[0]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff):
        if query not in s16_to_ref_dict:
            s16_to_ref_dict[query] = {subject_genome}
        else:
            s16_to_ref_dict[query].add(subject_genome)


s16_to_ref_txt_handle = open(s16_to_ref_txt, 'w')
for each_16s in s16_to_ref_dict:
    s16_to_ref_txt_handle.write('%s\t%s\n' % (each_16s, ','.join(s16_to_ref_dict[each_16s])))
s16_to_ref_txt_handle.close()
