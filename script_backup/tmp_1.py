
all_vs_all_16s = '/Users/songweizhi/Desktop/combined_16S_all_vs_all.txt.txt'

min_intra_genome_iden_dict = {}
for line in open(all_vs_all_16s):
    line_split     = line.strip().split('\t')
    query_genome   = line_split[0].split('_')[0]
    subject_genome = line_split[1].split('_')[0]
    iden           = float(line_split[2])
    aln_len        = int(line_split[3])
    # if (query_genome == 'p5') and (subject_genome == 'p5') and (line_split[0] != line_split[1]):
    #     print('%s\t%s\t%s' % (line_split[0], line_split[1], iden))


    if query_genome == subject_genome:
        print(line_split)
        print(query_genome)
        print(subject_genome)
        if query_genome not in min_intra_genome_iden_dict:
            min_intra_genome_iden_dict[query_genome] = iden
        else:
            if iden < min_intra_genome_iden_dict[query_genome]:
                min_intra_genome_iden_dict[query_genome] = iden

print(min_intra_genome_iden_dict)

for each in min_intra_genome_iden_dict:

    print('%s\t%s' % (each, min_intra_genome_iden_dict[each]))