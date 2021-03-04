
wd = '/Users/songweizhi/Desktop'
bin_id_file     = '%s/bin_id.txt'                               % wd
bin_taxon_file  = '%s/BH_ER_050417.bac120.summary.tsv'          % wd
mean_depth_file = '%s/BH_ER_050417_refined_bins_mean_depth.txt' % wd

link_output     = '%s/Kelp_0202_combined_linkages.txt'          % wd
barrnap_output  = '%s/BH_ER_050417_16S.txt'                     % wd
checkm_output   = '%s/BH_ER_050417_refined_MAG_qualities.2.txt' % wd


bin_id_list = []
for each_bin_id in open(bin_id_file):
    bin_id_list.append(each_bin_id.strip())


bin_taxon_dict = {}
for each_bin_taxon in open(bin_taxon_file):
    if not each_bin_taxon.startswith('user_genome	classification'):
        each_bin_taxon_split = each_bin_taxon.strip().split('\t')
        bin_taxon_dict[each_bin_taxon_split[0]] = each_bin_taxon_split[1]

bin_depth_dict = {}
for each_bin_depth in open(mean_depth_file):
    each_bin_depth_split = each_bin_depth.strip().split('\t')
    bin_depth_dict[each_bin_depth_split[0]] = each_bin_depth_split[2]



linked_mag_list = set()
for each_linkage in open(link_output):
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
        each_linkage_split = each_linkage.strip().split('\t')
        mag_id = each_linkage_split[1]
        linked_mag_list.add(mag_id)


barrnap_op_dict = {}
for each_barrnap_op in open(barrnap_output):
    if not each_barrnap_op.startswith('Genome	16S	Length'):
        each_barrnap_op_split = each_barrnap_op.strip().split('\t')
        barrnap_op_dict[each_barrnap_op_split[0]] = each_barrnap_op_split[2]


mag_completeness_dict = {}
for each_mag_quality in open(checkm_output):
    each_mag_quality_split = each_mag_quality.strip().split(',')
    mag_completeness_dict[each_mag_quality_split[0]] = each_mag_quality_split[1]


print('MAG\t16S\tcompleteness\tlinked\ttaxon')
for each_mag in sorted(bin_id_list):

    mag_depth = '0'
    if each_mag in bin_depth_dict:
        mag_depth = bin_depth_dict[each_mag]

    mag_taxon = 'unknown'
    if each_mag in bin_taxon_dict:
        mag_taxon = bin_taxon_dict[each_mag]

    with_16s_in_mag = 'no'
    if each_mag in barrnap_op_dict:
        with_16s_in_mag = barrnap_op_dict[each_mag]

    linked_to_16s = 'no'
    if each_mag in linked_mag_list:
        linked_to_16s = 'yes'

    mag_cpl = 0
    if each_mag in mag_completeness_dict:
        mag_cpl = mag_completeness_dict[each_mag]

    print('%s\t%s\t%s\t%s\t%s\t%s' % (each_mag, mag_depth, with_16s_in_mag, mag_cpl, linked_to_16s, mag_taxon))







