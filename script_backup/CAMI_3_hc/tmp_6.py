
r1_cigar_list = ['38S37=', '39S36=', '48S27=', '72=3S']

r1_cigar_all_clp = True
for each_cigar_r1 in r1_cigar_list:
    if ('S' not in each_cigar_r1) and ('s' not in each_cigar_r1):
        r1_cigar_all_clp = False

print(r1_cigar_all_clp)