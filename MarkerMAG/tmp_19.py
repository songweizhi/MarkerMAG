
rd2_read_to_extract_flanking_16s_r1_up = {1, 2, 3}
rd2_read_to_extract_flanking_16s_r2_up = {3, 4, 5}
rd2_read_to_extract_flanking_ctg_r1_up = {5, 6, 7}
rd2_read_to_extract_flanking_ctg_r2_up = {7, 8, 9}

rd2_read_to_extract_flanking_both_r1_up = set.union(rd2_read_to_extract_flanking_16s_r1_up, rd2_read_to_extract_flanking_ctg_r1_up)
rd2_read_to_extract_flanking_both_r2_up = set.union(rd2_read_to_extract_flanking_16s_r2_up, rd2_read_to_extract_flanking_ctg_r2_up)




print(rd2_read_to_extract_flanking_both_r1_up)
print(rd2_read_to_extract_flanking_both_r2_up)
print(rd2_read_to_extract_flanking_16s_r1_up)
print(rd2_read_to_extract_flanking_16s_r2_up)
print(rd2_read_to_extract_flanking_ctg_r1_up)
print(rd2_read_to_extract_flanking_ctg_r2_up)


print(round(9/2))