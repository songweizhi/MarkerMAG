
blast_best_hits = '/Users/songweizhi/Desktop/BH_ER_050417_vs_GTDB_r95_BestHit.tab'
ref_taxon_file  = '/Users/songweizhi/Desktop/GTDB_ssu_all_r95.txt'
query_taxon     = '/Users/songweizhi/Desktop/query_taxon.txt'

ref_taxon_dict = {}
for each_ref_seq in open(ref_taxon_file):
    each_ref_seq_split = each_ref_seq.strip().split(' ')
    ref_seq_id = each_ref_seq_split[0]
    ref_seq_taxon = ' '.join(each_ref_seq_split[1:])
    ref_taxon_dict[ref_seq_id] = ref_seq_taxon


query_to_hit_dict = {}
for each_hit in open(blast_best_hits):
    each_hit_split = each_hit.strip().split('\t')
    query_to_hit_dict[each_hit_split[0]] = each_hit_split[1]


query_taxon_handle = open(query_taxon, 'w')
for each_query in query_to_hit_dict:
    query_taxon = ref_taxon_dict[query_to_hit_dict[each_query]]
    query_taxon_handle.write('%s\t%s\n' % (each_query, query_taxon))
query_taxon_handle.close()




