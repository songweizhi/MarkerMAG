max_inter_genome = 0
min_intra_genome = 100
for line in open('/Users/songweizhi/Desktop/combined_16S_all_vs_all.txt'):

    line_split = line.strip().split('\t')
    query_genome    = line_split[0].split('_')[0]
    subject_genome  = line_split[1].split('_')[0]
    iden            = float(line_split[2])
    aln_len         = int(line_split[3])

    if aln_len >= 1300:
        if query_genome != subject_genome:
            if iden > max_inter_genome:
                max_inter_genome = iden
        if query_genome == subject_genome:
            if iden < min_intra_genome:
                min_intra_genome = iden

print('max_inter_genome: %s' % max_inter_genome)
print('min_intra_genome: %s' % min_intra_genome)
# max_inter_genome: 99.282
# min_intra_genome: 99.016
