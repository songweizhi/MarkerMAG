import os

# file in
strain = 'CP'
wd                  = '/Users/songweizhi/Desktop/distance/%s' % strain
fake_bin            = '%s/bin_%s.fna'                       % (wd, strain)
ref_genome          = '%s/ref_%s.fna'                       % (wd, strain)
ref_genome_gff      = '%s/%s.gff'                           % (wd, strain)
ctg_cov_cutoff      = 90
iden_cutoff         = 99.5

# file out
output_table        = '%s/%s_distance.txt'                  % (wd, strain)


# tmp file
blast_op            = '%s/%s_bin_vs_ref.tab'               % (wd, strain)
blast_op_best_hit   = '%s/%s_bin_vs_ref_best_hit.tab'      % (wd, strain)


blastn_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn'
blastn_cmd        = 'blastn -query %s -subject %s -out %s %s' % (fake_bin, ref_genome, blast_op, blastn_parameters)
os.system(blastn_cmd)

get_best_hit_cmd  = 'BioSAK BestHit -i %s -o %s' % (blast_op, blast_op_best_hit)
os.system(get_best_hit_cmd)


total_hit_num = 0
qualified_hit_num = 0
ref_len_dict = {}
recovered_ref_dict = {}
for each_hit in open(blast_op_best_hit):
    match_split = each_hit.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_genome = subject.split('_')[0]
    iden = float(match_split[2])
    align_len = int(match_split[3])
    sstart = int(match_split[8])
    send = int(match_split[9])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    query_bin_name = '_'.join(query.split('_')[:-1])
    subject_bin_name = '_'.join(subject.split('_')[:-1])
    coverage_q = float(align_len) * 100 / float(query_len)
    coverage_s = float(align_len) * 100 / float(subject_len)
    total_hit_num += 1
    if (coverage_q >= ctg_cov_cutoff) and (iden >= iden_cutoff):
        qualified_hit_num += 1
        current_recovered_ref_region = '-'.join([str(i) for i in sorted([sstart, send])])
        if subject not in recovered_ref_dict:
            recovered_ref_dict[subject] = [current_recovered_ref_region]
        else:
            recovered_ref_dict[subject].append(current_recovered_ref_region)
        if subject not in ref_len_dict:
            ref_len_dict[subject] = subject_len

rna_16s_region_dict = {}
for each_line in open(ref_genome_gff):
    if not each_line.startswith('#'):
        if ';' in each_line:
            line_split = each_line.strip().split('\t')
            ctg_id = line_split[0]
            gene_product = line_split[-1].split('=')[-1]
            if gene_product == '16S ribosomal RNA':
                current_16s_region = '-'.join([str(i) for i in sorted([int(line_split[3]), int(line_split[4])])])
                if ctg_id not in rna_16s_region_dict:
                    rna_16s_region_dict[ctg_id] = [current_16s_region]
                else:
                    rna_16s_region_dict[ctg_id].append(current_16s_region)


output_table_handle = open(output_table, 'w')
for ref_seq in recovered_ref_dict:
    ref_seq_recovered = recovered_ref_dict[ref_seq]
    for ctg_region in ref_seq_recovered:
        output_table_handle.write('%s\tcontig\t%s\t%s\n' % (ref_seq, ctg_region.split('-')[0], ctg_region.split('-')[1]))
for ref_seq in rna_16s_region_dict:
    ref_16s_recovered = rna_16s_region_dict[ref_seq]
    for region_16s in ref_16s_recovered:
        output_table_handle.write('%s\trRNA16S\t%s\t%s\n' % (ref_seq, region_16s.split('-')[0], region_16s.split('-')[1]))
output_table_handle.close()

