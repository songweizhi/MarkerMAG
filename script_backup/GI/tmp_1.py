
markermag_linkage = '/Users/songweizhi/Desktop/CAMI2_HMP_combined_linkages.txt'
GI_bin_depth_file = '/Users/songweizhi/Desktop/GI_bin_mean_depth.txt'

linked_mag_set = set()
for each_linkage in open(markermag_linkage):
    each_linkage_split = each_linkage.strip().split('\t')
    linked_mag = each_linkage_split[1]
    linked_mag_set.add(linked_mag)
print(linked_mag_set)

for mag_depth in open(GI_bin_depth_file):
    mag_depth_split = mag_depth.strip().split('\t')
    mag_id = '.'.join(mag_depth_split[0].split('.')[:-1])

    if mag_id in linked_mag_set:
        print('%s\tyes' % (mag_depth.strip()))
    else:
        print('%s\tno' % (mag_depth.strip()))

print(len(linked_mag_set))