import os
from Bio import SeqIO


###################################################### file in/out #####################################################

wd = '/Users/songweizhi/Desktop/ttt_dRep97'

# file in
drep_cdb_file               = '%s/Cdb.csv'                                      % wd
combined_GI_ref_16S         = '%s/combined_GI_ref_16S.ffn'                      % wd
matam_16s_seqs              = '%s/3_GI_assembled_16S_uclust_0.995.fasta'        % wd
iden_cutoff                 = 99.5
iden_cutoff                 = 100
aln_len_cutoff              = 500
cov_q_cutoff                = 90

# file out
r16s_to_cluster_file        = '%s/16S_to_cluster.txt'                           % wd
rename_16s_file             = '%s/16S_rename.txt'                               % wd
combined_GI_ref_16S_renamed = '%s/combined_GI_ref_16S_renamed.ffn'              % wd

# tmp files
matam_16s_blastn            = '%s/3_GI_assembled_16S_uclust_0.995_vs_refs.tab'  % wd


########################################################################################################################

# get new cluster id
cluster_id_default_to_new_dict = {}
cluster_id_index = 1
for each_mag in open(drep_cdb_file):
    if not each_mag.startswith('genome,secondary_cluster'):
        mag_cluster = each_mag.strip().split(',')[1]
        if mag_cluster not in cluster_id_default_to_new_dict:
            new_id = 'C%s' % cluster_id_index
            cluster_id_index += 1
            cluster_id_default_to_new_dict[mag_cluster] = new_id


cluster_to_mag_dict = {}
mag_to_cluster_dict = {}
for each_mag in open(drep_cdb_file):
    if not each_mag.startswith('genome,secondary_cluster'):
        each_mag_split = each_mag.strip().split(',')
        mag_file_name = each_mag_split[0]
        mag_file_name_no_ext = '.'.join(mag_file_name.split('.')[:-1])
        mag_cluster_default = each_mag_split[1]
        mag_cluster_new = cluster_id_default_to_new_dict[mag_cluster_default]
        mag_to_cluster_dict[mag_file_name_no_ext] = mag_cluster_new
        if mag_cluster_new not in cluster_to_mag_dict:
            cluster_to_mag_dict[mag_cluster_new] = [mag_file_name_no_ext]
        else:
            cluster_to_mag_dict[mag_cluster_new].append(mag_file_name_no_ext)


r16s_to_cluster_file_handle = open(r16s_to_cluster_file, 'w')
s16_to_cluster_dict = {}
cluster_to_s16_dict = {}
for record_16s in SeqIO.parse(combined_GI_ref_16S, 'fasta'):
    s16_genome = record_16s.id.split('_16S_')[0]
    s16_genome_cluster = mag_to_cluster_dict[s16_genome]
    r16s_to_cluster_file_handle.write('%s\t%s\n' % (record_16s.id, s16_genome_cluster))
    s16_to_cluster_dict[record_16s.id] = s16_genome_cluster
    if s16_genome_cluster not in cluster_to_s16_dict:
        cluster_to_s16_dict[s16_genome_cluster] = [record_16s.id]
    else:
        cluster_to_s16_dict[s16_genome_cluster].append(record_16s.id)
r16s_to_cluster_file_handle.close()


print('cluster_id_default_to_new_dict')
print(cluster_id_default_to_new_dict)
print()
print('cluster_to_mag_dict')
print(cluster_to_mag_dict)
print()
print('mag_to_cluster_dict')
print(mag_to_cluster_dict)
print()
print('cluster_to_s16_dict')
print(cluster_to_s16_dict)
print()
print('s16_to_cluster_dict')
print(s16_to_cluster_dict)
print()


blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 1'
blast_cmd = 'blastn -query %s -subject %s -out %s %s' % (matam_16s_seqs, combined_GI_ref_16S, matam_16s_blastn, blast_parameters)
#os.system(blast_cmd)

# parse blast results
matam_16s_to_cluster_dict = {}
for match in open(matam_16s_blastn):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff) and (coverage_q >= cov_q_cutoff):
        subject_cluster = s16_to_cluster_dict[subject]
        if query not in matam_16s_to_cluster_dict:
            matam_16s_to_cluster_dict[query] = {subject_cluster}
        else:
            matam_16s_to_cluster_dict[query].add(subject_cluster)

m = 0
n = 0
for matam_16s in matam_16s_to_cluster_dict:
    matched_clusters = matam_16s_to_cluster_dict[matam_16s]
    if len(matched_clusters) == 1:
        m += 1
    if len(matched_clusters) > 1:
        n += 1

        print('%s\t%s' % (matam_16s, matched_clusters))

print(m)
print(n)

# C63   d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus; s__Bacillus subtilis
'''
C58: ['OTU_97.1969.0']  195103
C59: ['OTU_97.38225.0'] 1502
C60: ['OTU_97.11107.0'] 1502
C61: ['OTU_97.9303.0']  195102
C62: ['OTU_97.29691.0'] 1502
C63: ['OTU_97.207.0']   289380
C64: ['OTU_97.31703.1'] 1502
'''

