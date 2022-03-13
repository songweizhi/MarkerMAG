import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


BH_ER_050417_vs_GTDB_txt       = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_zOTUs_vs_GTDB_r202_top1_classifications.txt'
BH_ER_050417_vs_SILVA_txt      = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_zOTUs_vs_SILVA_top1_classifications.txt'
BH_ER_050417_vs_GTDB_BLCA_txt  = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_zOTUs_vs_GTDB_BLCA_reformatted_1.txt'
BH_ER_050417_vs_SILVA_BLCA_txt = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_zOTUs_vs_SILVA_BLCA_reformatted_1.txt'
taxon_level                    = 'species'
confidence_score_cutoff        = 80

plot_file                      = '/Users/songweizhi/Desktop/taxon_16S/GTDB_vs_SILVA_confidence%s_%s.png' % (confidence_score_cutoff, taxon_level)


taxon_rank_index = {'domain': 0, 'phylum': 1, 'class': 2, 'order': 3, 'family': 4, 'genus': 5, 'species': 6}

best_hit_iden_dict_gtdb = {}
for each_otu in open(BH_ER_050417_vs_GTDB_txt):
    if not each_otu.startswith('Query	Reference	Identity'):
        each_otu_split = each_otu.strip().split('\t')
        otu_id = each_otu_split[0]
        best_hit_iden = float(each_otu_split[2])
        best_hit_iden_dict_gtdb[otu_id] = best_hit_iden

best_hit_iden_dict_silva = {}
for each_otu in open(BH_ER_050417_vs_SILVA_txt):
    if not each_otu.startswith('Query	Reference	Identity'):
        each_otu_split = each_otu.strip().split('\t')
        otu_id = each_otu_split[0]
        best_hit_iden = float(each_otu_split[2])
        best_hit_iden_dict_silva[otu_id] = best_hit_iden

similarity_list_gtdb = []
similarity_list_gtdb_lt_90 = 0
for each_otu in open(BH_ER_050417_vs_GTDB_BLCA_txt):
    each_otu_split = each_otu.strip().split('\t')
    otu_id = each_otu_split[0]
    taxon_list = each_otu_split[1].split(';')
    if len(taxon_list) == 7:
        taxon_to_process = taxon_list[taxon_rank_index[taxon_level]]
        confidence_score = float(taxon_to_process.split('(')[1][:-1])
        if confidence_score >= confidence_score_cutoff:
            similarity_with_best_hit =best_hit_iden_dict_gtdb.get(otu_id, 'na')
            if similarity_with_best_hit != 'na':
                similarity_list_gtdb.append(similarity_with_best_hit)
                if similarity_with_best_hit >= 90:
                    similarity_list_gtdb_lt_90 += 1

similarity_list_silva = []
similarity_list_silva_lt_90 = 0
for each_otu in open(BH_ER_050417_vs_SILVA_BLCA_txt):
    each_otu_split = each_otu.strip().split('\t')
    otu_id = each_otu_split[0]
    taxon_list = each_otu_split[1].split(';')
    if len(taxon_list) == 7:
        taxon_to_process = taxon_list[taxon_rank_index[taxon_level]]
        confidence_score = float(taxon_to_process.split('(')[-1][:-1])
        if confidence_score >= confidence_score_cutoff:
            similarity_with_best_hit = best_hit_iden_dict_gtdb.get(otu_id, 'na')
            if similarity_with_best_hit != 'na':
                similarity_list_silva.append(similarity_with_best_hit)
                if similarity_with_best_hit >= 90:
                    similarity_list_silva_lt_90 += 1

min_value = round(min([min(similarity_list_silva), min(similarity_list_gtdb)]) - 1)
tick_value_list = list(range(min_value, 101, 1))

bin_num_1 = list(range(round(min(similarity_list_silva)), 101, 1))
bin_num_2 = list(range(round(min(similarity_list_gtdb)), 101, 1))

plt.figure(figsize=(8,6))
plt.hist(similarity_list_silva, bins=len(bin_num_1)-1, alpha=0.5, label="SILVA (%s)" % len(similarity_list_silva))
plt.hist(similarity_list_gtdb, bins=len(bin_num_2)-1, alpha=0.5, label="GTDB (%s)" % len(similarity_list_gtdb))
plt.xticks(tick_value_list, size=15)
plt.yticks(size=15)
plt.title(("Confidence score = %s at %s level" % (confidence_score_cutoff, taxon_level)), fontsize=16)
plt.xlabel("Similarity with best blast hit (%)", size=16)
plt.ylabel("Number of query sequences", size=16)
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig(plot_file)
plt.close()

print('sequence similarity >=90 (GTDB): %s' % (similarity_list_gtdb_lt_90/len(similarity_list_gtdb)))
print('sequence similarity >=90 (SILVA): %s' % (similarity_list_silva_lt_90/len(similarity_list_silva)))
