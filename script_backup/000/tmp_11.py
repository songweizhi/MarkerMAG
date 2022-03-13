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


################################################# subsample OTU table ##################################################

# 57884-57886, 57903-57905 and 57940-57942
# BH_ER_050417: ['57884', '57885', '57886']
# BH_ER_110816: ['57903', '57904', '57905']
# CB_ER_080217: ['57940', '57941', '57942']

# file in
otu_seq             = '/Users/songweizhi/Desktop/taxon_16S/zOTUs_nc.fasta'
otu_table           = '/Users/songweizhi/Desktop/taxon_16S/otu_table2.txt'
sample_id           = ['57893', '57894', '57895', '57900', '57901', '57902', '57906', '57907', '57908', '57890', '57891', '57892', '57903', '57904', '57905', '57896', '57897', '57909', '57910', '57898', '57899', '57884', '57885', '57886', '57887', '57888', '57889', '57911', '57912', '57924', '57925', '57926', '57929', '57930', '57931', '57918', '57919', '57920', '57927', '57928', '57921', '57922', '57932', '57933', '57923', '57913', '57914', '57915', '57916', '57917', '57943', '57944', '57945', '57938', '57939', '57952', '57949', '57950', '57951', '57940', '57941', '57942', '57934', '57935', '57936', '57946', '57947', '57948', '57953', '57954', '57955', '57971', '57961', '57962', '57963', '57969', '57964', '57965', '57966', '57972', '57973', '57974', '57967', '57968', '57956', '57958', '57959', '57960']
sample_id           = ['57884', '57885', '57886']

# file out
otu_table_subset    = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_otu_table.txt'
otu_subset_id       = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_otu_id.txt'
otu_seq_subset      = '/Users/songweizhi/Desktop/taxon_16S/BH_ER_050417_zOTUs.fasta'

otu_subset_id_handle = open(otu_subset_id, 'w')
otu_table_subset_handle = open(otu_table_subset, 'w')
otu_table_subset_handle.write('ID\t%s\n' % '\t'.join(sample_id))
sample_index_dict = {}
for each_line in open(otu_table):
    each_line_split = each_line.strip().split('\t')
    if each_line.startswith('34830'):
        for each_sample in sample_id:
            sample_index = each_line_split.index(each_sample) + 1
            sample_index_dict[each_sample] = sample_index
    else:
        extract_element_list = []
        for each_sample in sample_id:
            sample_index = sample_index_dict.get(each_sample)
            sample_value = each_line_split[sample_index]
            extract_element_list.append(int(sample_value))
        if extract_element_list != ([0]*len(sample_id)):
            otu_table_subset_handle.write('%s\t%s\n' % (each_line_split[0], '\t'.join([str(i) for i in extract_element_list])))
            otu_subset_id_handle.write('%s\n' % each_line_split[0])
otu_table_subset_handle.close()
otu_subset_id_handle.close()

select_seq_cmd = 'BioSAK select_seq -seq %s -id %s -out %s -option 1' % (otu_seq, otu_subset_id, otu_seq_subset)
os.system(select_seq_cmd)
