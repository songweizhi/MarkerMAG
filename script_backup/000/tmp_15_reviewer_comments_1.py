

'''

cd /Users/songweizhi/Desktop/MarkerMAG_Matam
/Users/songweizhi/Software/usearch/usearch -cluster_fast ref_16S.fasta -id 0.99 -centroids ref_16S_uclust_0.99.fasta -uc ref_16S_uclust_0.99.uc -sort length -quiet

BioSAK BestHit -i unclustered_vs_refs.tab -o unclustered_vs_refs_BestHit.txt
BioSAK BestHit -i uclust_vs_refs.tab -o uclust_vs_refs_BestHit.txt

'''


# read in ref to cluster info
cluster_set = set()
ref_seq_to_cluster_dict = {}
for each in open('/Users/songweizhi/Desktop/MarkerMAG_Matam/ref_16S_uclust_0.99.txt'):
    each_split = each.strip().split('\t')
    cluster_id = each_split[0]
    ref_seq_list = each_split[1].split(',')
    cluster_set.add(cluster_id)
    for each_ref_seq in ref_seq_list:
        ref_seq_to_cluster_dict[each_ref_seq] = cluster_id
print(ref_seq_to_cluster_dict)


subset_set = set()
dict_of_dict = {}
for each in open('/Users/songweizhi/Desktop/MarkerMAG_Matam/unclustered_vs_refs_BestHit_aln1200.txt'):
    each_split = each.strip().split('\t')
    matam_16s_id = each_split[0]
    subset_rate = matam_16s_id.split('_')[-2]
    subset_rate = int(subset_rate)
    subset_set.add(subset_rate)
    ref_seq_id = each_split[1]
    ref_seq_cluster = ref_seq_to_cluster_dict[ref_seq_id]
    if subset_rate not in dict_of_dict:
        dict_of_dict[subset_rate] = {}
    if ref_seq_cluster not in dict_of_dict[subset_rate]:
        dict_of_dict[subset_rate][ref_seq_cluster] = 1
    else:
        dict_of_dict[subset_rate][ref_seq_cluster] += 1
print(dict_of_dict)


for each in open('/Users/songweizhi/Desktop/MarkerMAG_Matam/uclust_vs_refs_BestHit_aln1200.txt'):
    each_split = each.strip().split('\t')
    matam_16s_id = each_split[0]
    subset_rate = 'combined'
    ref_seq_id = each_split[1]
    ref_seq_cluster = ref_seq_to_cluster_dict[ref_seq_id]
    if subset_rate not in dict_of_dict:
        dict_of_dict[subset_rate] = {}
    if ref_seq_cluster not in dict_of_dict[subset_rate]:
        dict_of_dict[subset_rate][ref_seq_cluster] = 1
    else:
        dict_of_dict[subset_rate][ref_seq_cluster] += 1
print(dict_of_dict)


cluster_list_sorted = sorted([i for i in cluster_set])
subset_list_sorted  = sorted([i for i in subset_set])
subset_list_sorted.append('combined')

print('%s\t%s' % ('Cluster', '\t'.join([str(i) for i in cluster_list_sorted])))

for each_subset in subset_list_sorted:

    current_count_list = []
    for each_cluster in cluster_list_sorted:
        current_count_list.append(dict_of_dict[each_subset].get(each_cluster, '0'))

    print('%s\t%s' % (each_subset, '\t'.join([str(i) for i in current_count_list])))

