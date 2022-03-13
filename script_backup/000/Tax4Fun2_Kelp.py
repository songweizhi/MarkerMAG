import os
import numpy as np
import pandas as pd
from scipy import stats
import numpy as np
from statistics import median
from scipy.stats import mannwhitneyu
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def fun_table_to_dict(fun_table, col_list):

    ko_set = set()
    metagenome_fun_dict = {}
    for each_line in open(fun_table):
        if each_line.startswith('KO'):
            for each_sample in col_list:
                metagenome_fun_dict[each_sample] = {}
        else:
            each_line_split = each_line.strip().split('\t')
            ko_id = each_line_split[0]
            abun_list = each_line_split[1:]
            ko_set.add(ko_id)
            for (sample, ko_abun) in zip(col_list, abun_list):
                metagenome_fun_dict[sample][ko_id] = float(ko_abun)

    return metagenome_fun_dict, ko_set


def get_unused_otu_seq_list(logfile2):
    sample_id_counted = []
    unused_otu = []
    unused_seq = []
    for each in open(logfile2):
        each_split = each.strip().split(': ')
        if len(each_split) == 2:
            sample_id = each_split[0]
            sample_value = float(each_split[1])
            if sample_id not in sample_id_counted:
                unused_otu.append(sample_value)
                sample_id_counted.append(sample_id)
            else:
                unused_seq.append(sample_value)

    return unused_otu, unused_seq


################################################# subsample OTU table ##################################################

# 57884-57886, 57903-57905 and 57940-57942
# # BH_ER_050417: ['57884', '57885', '57886']
# # BH_ER_110816: ['57903', '57904', '57905']
# # CB_ER_080217: ['57940', '57941', '57942']
#
# # file in
# otu_seq             = '/Users/songweizhi/Desktop/Tax4Fun2/zOTUs_nc.fasta'
# otu_table           = '/Users/songweizhi/Desktop/Tax4Fun2/otu_table2.txt'
# sample_id           = ['57884', '57885', '57886', '57903', '57904', '57905', '57940', '57941', '57942']
#
# # file out
# otu_table_subset    = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_otu_table.txt'
# otu_subset_id       = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_otu_id.txt'
# otu_seq_subset      = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_zOTUs_nc.fasta'
#
# otu_subset_id_handle = open(otu_subset_id, 'w')
# otu_table_subset_handle = open(otu_table_subset, 'w')
# otu_table_subset_handle.write('ID\t%s\n' % '\t'.join(sample_id))
# sample_index_dict = {}
# for each_line in open(otu_table):
#     each_line_split = each_line.strip().split('\t')
#     if each_line.startswith('34830'):
#         for each_sample in sample_id:
#             sample_index = each_line_split.index(each_sample) + 1
#             sample_index_dict[each_sample] = sample_index
#     else:
#         extract_element_list = []
#         for each_sample in sample_id:
#             sample_index = sample_index_dict.get(each_sample)
#             sample_value = each_line_split[sample_index]
#             extract_element_list.append(int(sample_value))
#         if extract_element_list != ([0]*len(sample_id)):
#             otu_table_subset_handle.write('%s\t%s\n' % (each_line_split[0], '\t'.join([str(i) for i in extract_element_list])))
#             otu_subset_id_handle.write('%s\n' % each_line_split[0])
# otu_table_subset_handle.close()
# otu_subset_id_handle.close()
#
# select_seq_cmd = 'BioSAK select_seq -seq %s -id %s -out %s -option 1' % (otu_seq, otu_subset_id, otu_seq_subset)
# os.system(select_seq_cmd)


############################################ get metagenome functional table ###########################################

# # file in
# metagenome_vs_KEGG_best_hits_folder = '/Users/songweizhi/Desktop/Tax4Fun2/metagenome_vs_KEGG_best_hits'
# sample_list = ['BH_ER_050417_1', 'BH_ER_050417_2', 'BH_ER_050417_3', 'BH_ER_110816_1', 'BH_ER_110816_2', 'BH_ER_110816_3', 'CB_ER_080217_1', 'CB_ER_080217_2', 'CB_ER_080217_3']
#
# # file out
# metagenome_fun_tab = '/Users/songweizhi/Desktop/Tax4Fun2/9_samples_metagenome_fun_profile.txt'
#
# # read in blast results
# overall_identified_ko_set = set()
# sample_ko_count_dict = {}
# sample_total_ko_dict = {}
# for sample_with_rep in sample_list:
#     diamond_op = '%s/%s_vs_Tax4Fun2_KEGG_best_hit.tab' % (metagenome_vs_KEGG_best_hits_folder, sample_with_rep)
#     current_sample_ko_count_dict = {}
#     current_sample_total_ko = 0
#     if os.path.isfile(diamond_op) is True:
#         current_sample_total_ko = 0
#         for each_line in open(diamond_op):
#             each_line_split = each_line.strip().split('\t')
#             subject = each_line_split[1]
#             subject_list = [subject]
#             if len(subject) > 6:
#                 subject_list = subject.split(':')
#             for each_ko in subject_list:
#                 overall_identified_ko_set.add(each_ko)
#                 current_sample_total_ko += 1
#                 if each_ko not in current_sample_ko_count_dict:
#                     current_sample_ko_count_dict[each_ko] = 1
#                 else:
#                     current_sample_ko_count_dict[each_ko] += 1
#     sample_ko_count_dict[sample_with_rep] = current_sample_ko_count_dict
#     sample_total_ko_dict[sample_with_rep] = current_sample_total_ko
#
# metagenome_fun_tab_handle = open(metagenome_fun_tab, 'w')
# metagenome_fun_tab_handle.write('KO\t%s\n' % '\t'.join(sample_list))
# for each_ko in sorted([i for i in overall_identified_ko_set]):
#     current_ko_freq_list = []
#     for sample_with_rep in sample_list:
#         ko_count = sample_ko_count_dict[sample_with_rep].get(each_ko, 0)
#         samlple_total_ko = sample_total_ko_dict[sample_with_rep]
#         ko_freq = 0
#         if samlple_total_ko > 0:
#             ko_freq = 0
#             if ko_count > 0:
#                 ko_freq = ko_count/samlple_total_ko
#         current_ko_freq_list.append(str(ko_freq))
#     metagenome_fun_tab_handle.write('%s\t%s\n' % (each_ko, '\t'.join(current_ko_freq_list)))
# metagenome_fun_tab_handle.close()

########################################################################################################################

fun_table_metagenome                            = '/Users/songweizhi/Desktop/Tax4Fun2/9_samples_metagenome_fun_profile.txt'
fun_table_Tax4Fun2_default                      = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_Tax4Fun2_output_default_db/functional_prediction_CN_T.txt'
fun_table_Tax4Fun2_with_prior_16S               = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_Tax4Fun2_output_MAG_with_prior_16S/functional_prediction.txt'
fun_table_Tax4Fun2_with_prior_and_linked_16S    = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_Tax4Fun2_output_MAG_with_prior_and_linked_16S/functional_prediction.txt'
sample_list                                     = ['BH_ER_050417_1', 'BH_ER_050417_2', 'BH_ER_050417_3', 'BH_ER_110816_1', 'BH_ER_110816_2', 'BH_ER_110816_3', 'CB_ER_080217_1', 'CB_ER_080217_2', 'CB_ER_080217_3']

metagenome_fun_dict, metagenome_identified_ko_set                           = fun_table_to_dict(fun_table_metagenome, sample_list)
Tax4Fun2_default_fun_dict, Tax4Fun2_default_identified_ko_set               = fun_table_to_dict(fun_table_Tax4Fun2_default, sample_list)
Tax4Fun2_with_prior_16S_fun_dict, Tax4Fun2_with_prior_16S_identified_ko_set = fun_table_to_dict(fun_table_Tax4Fun2_with_prior_16S, sample_list)
Tax4Fun2_with_prior_and_linked_16S_fun_dict, Tax4Fun2_with_prior_and_linke_16S_identified_ko_set = fun_table_to_dict(fun_table_Tax4Fun2_with_prior_and_linked_16S, sample_list)

# combine ko sets
all_identified_ko_set = set.union(metagenome_identified_ko_set, Tax4Fun2_default_identified_ko_set, Tax4Fun2_with_prior_16S_identified_ko_set, Tax4Fun2_with_prior_and_linke_16S_identified_ko_set)
all_identified_ko_list = sorted([i for i in all_identified_ko_set])

print('Sample\tDefault\tWith_prior_16S\tWith_linked_16S')
rho_list_default = []
rho_list_with_prior_16S = []
rho_list_with_prior_and_linked_16S = []
for each_sample in sample_list:

    fun_dict_metagenome                 = metagenome_fun_dict[each_sample]
    fun_dict_Tax4Fun2_default           = Tax4Fun2_default_fun_dict[each_sample]
    fun_dict_Tax4Fun2_with_prior_16S    = Tax4Fun2_with_prior_16S_fun_dict[each_sample]
    fun_dict_Tax4Fun2_with_prior_and_linked_16S    = Tax4Fun2_with_prior_and_linked_16S_fun_dict[each_sample]

    abun_list_metagenome = []
    abun_list_Tax4Fun2_default = []
    abun_list_Tax4Fun2_with_prior_16S = []
    abun_list_Tax4Fun2_with_prior_and_linked_16S = []
    for each_ko in all_identified_ko_list:
        abun_list_metagenome.append(fun_dict_metagenome.get(each_ko, 0))
        abun_list_Tax4Fun2_default.append(fun_dict_Tax4Fun2_default.get(each_ko, 0))
        abun_list_Tax4Fun2_with_prior_16S.append(fun_dict_Tax4Fun2_with_prior_16S.get(each_ko, 0))
        abun_list_Tax4Fun2_with_prior_and_linked_16S.append(fun_dict_Tax4Fun2_with_prior_and_linked_16S.get(each_ko, 0))

    default_rho, default_pval                                       = stats.spearmanr(abun_list_metagenome, abun_list_Tax4Fun2_default)
    with_prior_16S_rho, with_prior_16S_pval                         = stats.spearmanr(abun_list_metagenome, abun_list_Tax4Fun2_with_prior_16S)
    with_prior_and_linked_16S_rho, with_prior_and_linked_16S_pval   = stats.spearmanr(abun_list_metagenome, abun_list_Tax4Fun2_with_prior_and_linked_16S)

    print('%s\t%s\t%s\t%s'         % (each_sample, float("{0:.4f}".format(default_rho)), float("{0:.4f}".format(with_prior_16S_rho)), float("{0:.4f}".format(with_prior_and_linked_16S_rho))))
    rho_list_default.append(default_rho)
    rho_list_with_prior_16S.append(with_prior_16S_rho)
    rho_list_with_prior_and_linked_16S.append(with_prior_and_linked_16S_rho)


########################################################################################################################

logfile2_default                    = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_Tax4Fun2_output_default_db/logfile2_CN_T.txt'
logfile2_with_prior_16S             = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_Tax4Fun2_output_MAG_with_prior_16S/logfile2.txt'
logfile2_with_prior_and_linked_16S  = '/Users/songweizhi/Desktop/Tax4Fun2/9_sample_Tax4Fun2_output_MAG_with_prior_and_linked_16S/logfile2.txt'
output_plot                         = '/Users/songweizhi/Desktop/Tax4Fun2_stats.svg'

# get unused_otu and unused_seq lists
unused_otu_default,                   unused_seq_default                   = get_unused_otu_seq_list(logfile2_default)
unused_otu_with_prior_16S,            unused_seq_with_prior_16S            = get_unused_otu_seq_list(logfile2_with_prior_16S)
unused_otu_with_prior_and_linked_16S, unused_seq_with_prior_and_linked_16S = get_unused_otu_seq_list(logfile2_with_prior_and_linked_16S)

num_list_of_list = []
# add rho
num_list_of_list.append(rho_list_default)
num_list_of_list.append(rho_list_with_prior_16S)
num_list_of_list.append(rho_list_with_prior_and_linked_16S)
# add unused_otu
num_list_of_list.append(unused_otu_default)
num_list_of_list.append(unused_otu_with_prior_16S)
num_list_of_list.append(unused_otu_with_prior_and_linked_16S)
# add unused_seq
num_list_of_list.append(unused_seq_default)
num_list_of_list.append(unused_seq_with_prior_16S)
num_list_of_list.append(unused_seq_with_prior_and_linked_16S)

color_list = ['wheat', 'lightblue', 'lightgreen', 'wheat', 'lightblue', 'lightgreen', 'wheat', 'lightblue', 'lightgreen']

box_plot = plt.boxplot(num_list_of_list, showfliers=False, patch_artist=True, whiskerprops=dict(color='grey', linewidth=1), capprops=dict(color='grey'))
plt.xticks([])
plt.ylabel('Spearman corfficient / fraction unused')
plt.title('Correlation to metagenome        Unused OTUs                  Unused sequences    ', fontsize=9)
plt.axvline(x=3.5, c='black', linestyle='dashed', dashes=(5, 10), linewidth=0.5)
plt.axvline(x=6.5, c='black', linestyle='dashed', dashes=(5, 10), linewidth=0.5)

# set the color pf box
for box, color in zip(box_plot['boxes'], color_list):
    box.set(linewidth=0)
    box.set_facecolor(color)
    box.set_alpha(0.85)
    box.set_zorder(1)

# add dots
pos_index = 1
for num_list in num_list_of_list:
    plt.plot(np.random.normal(pos_index, 0.075, len(num_list)), num_list, '.', alpha=0.8, marker='o', color='dimgrey', markersize=3, markeredgewidth=0, zorder=2)
    pos_index += 1

plt.savefig(output_plot, format='svg')
plt.close()


rho_U1_default_vs_prior, rho_p_default_vs_prior = mannwhitneyu(rho_list_default, rho_list_with_prior_16S)
rho_U1_prior_vs_linked,  rho_p_prior_vs_linked  = mannwhitneyu(rho_list_with_prior_16S, rho_list_with_prior_and_linked_16S)
otu_U1_default_vs_prior, otu_p_default_vs_prior = mannwhitneyu(unused_otu_default, unused_otu_with_prior_16S)
otu_U1_prior_vs_linked,  otu_p_prior_vs_linked  = mannwhitneyu(unused_otu_with_prior_16S, unused_otu_with_prior_and_linked_16S)
seq_U1_default_vs_prior, seq_p_default_vs_prior = mannwhitneyu(unused_seq_default, unused_seq_with_prior_16S)
seq_U1_prior_vs_linked,  seq_p_prior_vs_linked  = mannwhitneyu(unused_seq_with_prior_16S, unused_seq_with_prior_and_linked_16S)

print('Coefficient(Median):\t%s\t%s\t%s'    % (median(rho_list_default),   median(rho_list_with_prior_16S),   median(rho_list_with_prior_and_linked_16S)))
print('Unused OTUs(Median):\t%s\t%s\t%s'    % (median(unused_otu_default), median(unused_otu_with_prior_16S), median(unused_otu_with_prior_and_linked_16S)))
print('Unused seqs(Median):\t%s\t%s\t%s'    % (median(unused_seq_default), median(unused_seq_with_prior_16S), median(unused_seq_with_prior_and_linked_16S)))
print('rho default_vs_prior:\t%s\t%s'       % (rho_U1_default_vs_prior, rho_p_default_vs_prior))
print('rho prior_vs_linked:\t%s\t%s'        % (rho_U1_prior_vs_linked,  rho_p_prior_vs_linked))
print('otu default_vs_prior:\t%s\t%s'       % (otu_U1_default_vs_prior, otu_p_default_vs_prior))
print('otu prior_vs_linked:\t%s\t%s'        % (otu_U1_prior_vs_linked,  otu_p_prior_vs_linked))
print('seq default_vs_prior:\t%s\t%s'       % (seq_U1_default_vs_prior, seq_p_default_vs_prior))
print('seq prior_vs_linked:\t%s\t%s'        % (seq_U1_prior_vs_linked,  seq_p_prior_vs_linked))

