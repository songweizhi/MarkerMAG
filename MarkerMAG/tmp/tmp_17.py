


name_list = ['a', 'b', 'c', 'd']
value_list = [4, 7, 2, 0]

name_list = ['a']
value_list = [4]

sorted_best_16s_list = [[seq_id, mean_iden] for mean_iden, seq_id in sorted(zip(value_list, name_list), reverse=True)]

best_matched_marker = sorted_best_16s_list[0][0]
best_matched_marker_iden = sorted_best_16s_list[0][1]


print(best_matched_marker)
print(best_matched_marker_iden)
