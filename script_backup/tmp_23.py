
ko_id_to_num_dict = {}
for each in open('/Users/songweizhi/Desktop/DP_kegg.txt'):
    each_split = each.strip().split('\t')
    if len(each_split) > 1:
        print(each_split)

        if each_split[1] not in ko_id_to_num_dict:
            ko_id_to_num_dict[each_split[1]] = 1
        else:
            ko_id_to_num_dict[each_split[1]] += 1

print(ko_id_to_num_dict)

total_num = 0
for each_ko in ko_id_to_num_dict:
    total_num += ko_id_to_num_dict[each_ko]
print(total_num)
