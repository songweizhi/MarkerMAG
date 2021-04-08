
sub_num_dict = {}
for each_line in open('/Users/songweizhi/Desktop/STY_Merged_OTU01-Transporter.tbl'):
    sub_name = each_line.strip().split('\t')[2]
    if sub_name not in sub_num_dict:
        sub_num_dict[sub_name] = 1
    else:
        sub_num_dict[sub_name] += 1

for each_sub in sub_num_dict:
    print('%s\t%s' % (each_sub, sub_num_dict[each_sub]))
