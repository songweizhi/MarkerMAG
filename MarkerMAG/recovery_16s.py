import os
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def overlap_between_list(list_1, list_2):

    overlap = False
    for list_1_element in list_1:
        if list_1_element in list_2:
            overlap = True

    return overlap


def gnm_level_asessment(linked_mag_dict_rd1):

    rd1_correct_num = 0
    rd1_wrong_gnm = set()
    rd1_unknown_gnm = set()
    rd1_ambiguous_gnm = set()
    for each_rd_1_mag in linked_mag_dict_rd1:
        mag_assess = linked_mag_dict_rd1[each_rd_1_mag]
        if len(mag_assess) == 1:
            if mag_assess == {'Correct'}:
                rd1_correct_num += 1
            elif mag_assess == {'Wrong'}:
                rd1_wrong_gnm.add(each_rd_1_mag)
            elif mag_assess == {'Unknown'}:
                rd1_unknown_gnm.add(each_rd_1_mag)
        elif len(mag_assess) == 2:
            if ('Correct' in mag_assess) and ('Unknown' in mag_assess):
                rd1_correct_num += 1
            if ('Correct' in mag_assess) and ('Wrong' in mag_assess):
                rd1_ambiguous_gnm.add(each_rd_1_mag)
            if ('Wrong' in mag_assess) and ('Unknown' in mag_assess):
                rd1_wrong_gnm.add(each_rd_1_mag)
        elif len(mag_assess) == 3:
            rd1_ambiguous_gnm.add(each_rd_1_mag)

    return rd1_correct_num, rd1_unknown_gnm, rd1_wrong_gnm, rd1_ambiguous_gnm


########################################################################################################################

wd = '/Users/songweizhi/Desktop/assess_linkages_Oral'

########## reference to cluster ##########

drep_ani_cutoff             = 97
drep_cdb_file               = '%s/file_in/Cdb_%s.csv'                                       % (wd, drep_ani_cutoff)


########## 16S to reference ##########

perform_blastn_16s_vs_refs  = False  # True or False
combined_GI_ref_16S         = '%s/file_in/combined_Oral_ref_16S.ffn'                     % wd
matam_16s_seqs              = '%s/file_in/CAMI_Oral_138_16S_0.999.polished.fa'           % wd
matam_16s_blastn            = '%s/file_in/CAMI_Oral_138_16S_0.999.polished_vs_ref.tab'   % wd
iden_cutoff_16s             = 99.5  # 99.3 (best), 99.5
aln_len_cutoff_16s          = 500
cov_q_cutoff_16s            = 70
total_query_mag_num         = 87

matam_len_cutoff_list       = [500, 900, 1000, 1200, 1500]

'''
module load blast+/2.10.1
cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/oral_799_ref_genomes_barrnap_16s_wd
makeblastdb -in combined_Oral_ref_16S.ffn -dbtype nucl -parse_seqids -logfile /dev/null

cd /srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral
blastn -query CAMI_Oral_138_16S_0.999.polished.fa -db oral_799_ref_genomes_barrnap_16s_wd/combined_Oral_ref_16S.ffn -out CAMI_Oral_138_16S_0.999.polished_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12
blastn -query CAMI_Oral_138_16S_0.999.polished_min900.fa -db oral_799_ref_genomes_barrnap_16s_wd/combined_Oral_ref_16S.ffn -out CAMI_Oral_138_16S_0.999.polished_min900_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12
blastn -query CAMI_Oral_138_16S_0.999.polished_min1000.fa -db oral_799_ref_genomes_barrnap_16s_wd/combined_Oral_ref_16S.ffn -out CAMI_Oral_138_16S_0.999.polished_min1000_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12
blastn -query CAMI_Oral_138_16S_0.999.polished_min1100.fa -db oral_799_ref_genomes_barrnap_16s_wd/combined_Oral_ref_16S.ffn -out CAMI_Oral_138_16S_0.999.polished_min1100_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12
blastn -query CAMI_Oral_138_16S_0.999.polished_min1200.fa -db oral_799_ref_genomes_barrnap_16s_wd/combined_Oral_ref_16S.ffn -out CAMI_Oral_138_16S_0.999.polished_min1200_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12
blastn -query CAMI_Oral_138_16S_0.999.polished_min1300.fa -db oral_799_ref_genomes_barrnap_16s_wd/combined_Oral_ref_16S.ffn -out CAMI_Oral_138_16S_0.999.polished_min1300_vs_ref.tab -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 12
'''


################################################# reference to cluster #################################################

cluster_to_ref_dict = {}
ref_to_cluster_dict = {}
for each_ref in open(drep_cdb_file):
    if not each_ref.startswith('genome,secondary_cluster'):
        each_ref_split = each_ref.strip().split(',')
        ref_file_name = each_ref_split[0]
        ref_file_name_no_ext = '.'.join(ref_file_name.split('.')[:-1])
        ref_cluster = 'C' + each_ref_split[1]
        ref_to_cluster_dict[ref_file_name_no_ext] = ref_cluster
        if ref_cluster not in cluster_to_ref_dict:
            cluster_to_ref_dict[ref_cluster] = [ref_file_name_no_ext]
        else:
            cluster_to_ref_dict[ref_cluster].append(ref_file_name_no_ext)


################################################### 16S to reference ###################################################

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 1'
blast_cmd = 'blastn -query %s -subject %s -out %s %s' % (matam_16s_seqs, combined_GI_ref_16S, matam_16s_blastn, blast_parameters)
if perform_blastn_16s_vs_refs is True:
    os.system(blast_cmd)


for matam_len_cutoff in matam_len_cutoff_list:

    # get matam_16s_to_cluster_dict
    matam_16s_to_cluster_dict = {}
    for match in open(matam_16s_blastn):
        match_split = match.strip().split('\t')
        query = match_split[0]
        subject = match_split[1]
        subject_gnm = subject.split('_16S_')[0]
        subject_cluster = ref_to_cluster_dict[subject_gnm]
        iden = float(match_split[2])
        aln_len = int(match_split[3])
        query_len = int(match_split[12])
        subject_len = int(match_split[13])
        coverage_q = float(aln_len) * 100 / float(query_len)

        if query_len >= matam_len_cutoff:
            if (iden >= iden_cutoff_16s) and (aln_len > aln_len_cutoff_16s) and (coverage_q >= cov_q_cutoff_16s):
                if query not in matam_16s_to_cluster_dict:
                    matam_16s_to_cluster_dict[query] = {subject_cluster}
                else:
                    matam_16s_to_cluster_dict[query].add(subject_cluster)


    cluster_to_matam_16s_dict = {}
    for matam_16s in matam_16s_to_cluster_dict:
        matched_clusters =  matam_16s_to_cluster_dict[matam_16s]
        for matched_cluster in matched_clusters:
            if matched_cluster not in cluster_to_matam_16s_dict:
                cluster_to_matam_16s_dict[matched_cluster] = {matam_16s}
            else:
                cluster_to_matam_16s_dict[matched_cluster].add(matam_16s)


    print('%s\t%s/%s' % (matam_len_cutoff, len(cluster_to_matam_16s_dict), len(cluster_to_ref_dict)))

