import os
from Bio import SeqIO


# file in
wd                  = '/Users/songweizhi/Desktop/get_fake_bins_wd'
spades_ctg          = '%s/scaffolds_min2000_spades.fa'               % wd
combined_ref        = '%s/combined_refs.fna'                         % wd
fake_bin_folder     = '%s/fake_bins'                                 % wd
iden_cutoff         = 99.5
coverage_q_cutoff   = 90


# tmp files
blast_op         = '%s/spades_scaffolds_min2000_vs_refs.tab'         % wd
blast_op_BestHit = '%s/spades_scaffolds_min2000_vs_refs_BestHit.tab' % wd


# makeblastdb, run blastn and keep best hit
makeblastdb_cmd   = 'makeblastdb -in %s -dbtype nucl -parse_seqids' % combined_ref
blastn_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 4'
blastn_cmd        = 'blastn -query %s -db %s -out %s %s' % (spades_ctg, combined_ref, blast_op, blastn_parameters)
get_best_hit_cmd  = 'BioSAK BestHit -i %s -o %s' % (blast_op, blast_op_BestHit)
# os.system(makeblastdb_cmd)
# os.system(blastn_cmd)
# os.system(get_best_hit_cmd)


ref_id_set = set()
ref_len_dict = {}
total_ref_seq = 0
for seq_record in SeqIO.parse(combined_ref, 'fasta'):
    seq_genome = seq_record.id.split('_')[0]
    ref_id_set.add(seq_genome)
    total_ref_seq += len(seq_record.seq)
    if seq_genome not in ref_len_dict:
        ref_len_dict[seq_genome] = len(seq_record.seq)
    else:
        ref_len_dict[seq_genome] += len(seq_record.seq)


total_matched_ctg = 0
mag_size_dict = {}
ctg_to_genome_dict = {}
for each_hit in open(blast_op_BestHit):
    match_split = each_hit.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_genome = subject.split('_')[0]
    iden = float(match_split[2])
    align_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    query_bin_name = '_'.join(query.split('_')[:-1])
    subject_bin_name = '_'.join(subject.split('_')[:-1])
    coverage_q = float(align_len) * 100 / float(query_len)
    coverage_s = float(align_len) * 100 / float(subject_len)
    if (iden >= iden_cutoff) and (coverage_q >= coverage_q_cutoff):
        total_matched_ctg += query_len
        ctg_to_genome_dict[query] = subject_genome
        if subject_genome not in mag_size_dict:
            mag_size_dict[subject_genome] = query_len
        else:
            mag_size_dict[subject_genome] += query_len


print('MAG\tLength(Kbp)\tCompleteness(%)')
for each_ref in sorted([i for i in ref_id_set]):
    ref_len = ref_len_dict[each_ref]
    mag_len = 0
    if each_ref in mag_size_dict:
        mag_len = mag_size_dict[each_ref]
    mag_len_kbp = float("{0:.2f}".format(mag_len/(1024)))
    recovery = mag_len*100/ref_len
    recovery = float("{0:.2f}".format(recovery))
    print('%s\t%s\t%s' % (each_ref, mag_len_kbp, recovery))


if os.path.isdir(fake_bin_folder) is True:
    os.system('rm -r %s' % fake_bin_folder)
os.mkdir(fake_bin_folder)


total_ctg_len = 0
for ctg_record in SeqIO.parse(spades_ctg, 'fasta'):
    total_ctg_len += len(ctg_record.seq)
    if ctg_record.id in ctg_to_genome_dict:
        assignment = ctg_to_genome_dict[ctg_record.id]
        fake_bin_file = '%s/%s.fa' % (fake_bin_folder, assignment)
        fake_bin_file_handle = open(fake_bin_file, 'a')
        fake_bin_file_handle.write('>%s\n' % ctg_record.id)
        fake_bin_file_handle.write('%s\n'  % ctg_record.seq)
        fake_bin_file_handle.close()


total_ref_seq_Mbp       = total_ref_seq/(1024*1024)
total_ref_seq_Mbp       = float("{0:.2f}".format(total_ref_seq_Mbp))
total_ctg_len_Mbp       = total_ctg_len/(1024*1024)
total_ctg_len_Mbp       = float("{0:.2f}".format(total_ctg_len_Mbp))
total_matched_ctg_Mbp   = total_matched_ctg/(1024*1024)
total_matched_ctg_Mbp   = float("{0:.2f}".format(total_matched_ctg_Mbp))
print('Total length of reference genomes: %s Mbp' % total_ref_seq_Mbp)
print('Total length of contigs >= 2000bp: %s Mbp' % total_ctg_len_Mbp)
print('Total length of matched contigs: %s Mbp'   % total_matched_ctg_Mbp)
print('Fake bins exported to: %s'                 % fake_bin_folder)

'''

MAG	Length(Kbp)	Completeness(%)
AS	3884.6	86.46
CA	1937.62	52.9
CG	2525.75	78.15
CP	1603.3	50.41
DA	3579.45	73.44
DG	4286.35	90.4
DM	4132.4	86.83
EC	3325.88	73.37
EV	3633.37	66.34
FA	2114.1	60.08
FP	2052.16	97.0
HB	2.01	0.06
HR	2651.31	84.21
HT	3291.99	87.71
MS	2923.64	80.44
ND	0.0	0.0
NG	3081.81	83.3
NO	3483.43	82.68
OU	289.21	14.43
PS	3302.12	73.5
SB	1512.16	34.72
SP	1697.96	93.86
SR	2348.46	76.16
SS	3697.88	81.36
TC	3189.34	74.98
TR	1310.49	25.67
Total length of reference genomes: 100.13 Mbp
Total length of contigs >= 2000bp: 88.51 Mbp
Total length of matched contigs: 64.31 Mbp
Fake bins exported to: /Users/songweizhi/Desktop/get_fake_bins_wd/fake_bins

'''