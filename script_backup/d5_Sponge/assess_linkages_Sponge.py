import os
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


###################################################### file in/out #####################################################

wd                  = '/Users/songweizhi/Desktop/Sponge'
taxonomy_mag        = '%s/Sponge_MAGs.summary.tsv'                                                       % wd
taxonomy_16s        = '%s/SRR9841444_assembled_16S_uclust_0.999_vs_GTDB_classifications.txt'             % wd
taxonomy_16s_blca   = '%s/SRR9841444_assembled_16S_uclust_0.999_blca_vs_gtdb.txt'                        % wd
bin_id_file         = '%s/bin_id.txt'                                                                   % wd
barrnap_output      = '%s/Sponge_MAGs_16S.txt'                                                                   % wd
checkm_output       = '%s/Sponge_MAGs_quality_formatted.txt'                                                                             % wd
mean_depth_file     = '%s/Sponge_MAGs_mean_depth.txt'                                                                             % wd
gtdb_gnm_to_ssu_txt = '%s/GTDB_ssu_all_r95_gnm_to_ssu.txt'                                              % wd
seq_file_16s        = '%s/SRR9841444_assembled_16S_uclust_0.999.fasta'                                   % wd
linkage_file        = '%s/Sponge_linkages_by_genome.txt'                                                    % wd

'''

'''


################################################### define file name ###################################################

linkage_file_path, linkage_file_basename, linkage_file_extension = sep_path_basename_ext(linkage_file)
linkage_file_with_assessment        = '%s/%s_with_assessment%s'         % (linkage_file_path, linkage_file_basename, linkage_file_extension)
linkage_file_with_assessment_by_mag = '%s/%s_with_assessment_by_MAG%s'  % (linkage_file_path, linkage_file_basename, linkage_file_extension)

################################################### read in metadata ###################################################

len_dict_16s = {}
for each_16s in SeqIO.parse(seq_file_16s, 'fasta'):
    len_dict_16s[each_16s.id] = len(each_16s.seq)

gtdb_gnm_to_ssu_dict = {}
for each_gnm in open(gtdb_gnm_to_ssu_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    ssu_id_list = each_gnm_split[1].split(',')
    gtdb_gnm_to_ssu_dict[gnm_id] = ssu_id_list

bin_id_list = []
for each_bin_id in open(bin_id_file):
    bin_id_list.append(each_bin_id.strip())

bin_closest_ref_dict = {}
bin_taxon_dict = {}
for each_bin_taxon in open(taxonomy_mag):
    if not each_bin_taxon.startswith('user_genome	classification'):
        each_bin_taxon_split = each_bin_taxon.strip().split('\t')
        closest_ref = each_bin_taxon_split[7]
        bin_taxon_dict[each_bin_taxon_split[0]] = each_bin_taxon_split[1]
        bin_closest_ref_dict[each_bin_taxon_split[0]] = closest_ref

mag_completeness_dict = {}
for each_mag_quality in open(checkm_output):
    each_mag_quality_split = each_mag_quality.strip().split(',')
    mag_completeness_dict[each_mag_quality_split[0]] = each_mag_quality_split[1]

bin_depth_dict = {}
for each_bin_depth in open(mean_depth_file):
    each_bin_depth_split = each_bin_depth.strip().split('\t')
    bin_depth_dict[each_bin_depth_split[0]] = each_bin_depth_split[2]

barrnap_op_dict = {}
for each_barrnap_op in open(barrnap_output):
    if not each_barrnap_op.startswith('Genome\t16S\tLength(bp)'):
        each_barrnap_op_split = each_barrnap_op.strip().split('\t')
        barrnap_op_dict[each_barrnap_op_split[0]] = each_barrnap_op_split[2]

s16_taxon_dict = {}
s16_iden_dict = {}
s16_aln_dict = {}
for each_16s_taxon in open(taxonomy_16s):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    id_16s = each_16s_taxon_split[0]
    taxon_best_hit = each_16s_taxon_split[4].split(' [')[0]
    iden_with_best_match = each_16s_taxon_split[2]
    aln_with_best_match = each_16s_taxon_split[3]
    s16_taxon_dict[id_16s] = taxon_best_hit
    s16_iden_dict[id_16s]  = iden_with_best_match
    s16_aln_dict[id_16s]   = aln_with_best_match

s16_taxon_blca_dict = {}
for each_16s_taxon in open(taxonomy_16s_blca):
    each_16s_taxon_split = each_16s_taxon.strip().split('\t')
    s16_taxon_blca_dict[each_16s_taxon_split[0]] = each_16s_taxon_split[1]

s16_taxon_blca_dict_formatted = {}
for each_16s in s16_taxon_blca_dict:
    taxon_blca_raw = s16_taxon_blca_dict[each_16s]
    formatted_taxon_str = 'Unclassified'
    if taxon_blca_raw != 'Unclassified':
        taxon_blca_raw_split_1 = taxon_blca_raw.strip().split(':')[1:]
        formatted_taxon_list = []
        for each_str in taxon_blca_raw_split_1:
            each_str_split = each_str.split(';')
            taxon_with_confidence = '%s(%s)' % (each_str_split[0], each_str_split[1][:5])
            formatted_taxon_list.append(taxon_with_confidence)
        formatted_taxon_str = ';'.join(formatted_taxon_list)
    s16_taxon_blca_dict_formatted[each_16s] = formatted_taxon_str


###################################################### by linkage ######################################################

linkage_file_with_assessment_handle = open(linkage_file_with_assessment, 'w')
total_linkage = 0
correct_linkage_p = 0
linkage_with_unknown_p = 0
correct_linkage_c = 0
linkage_with_unknown_c = 0
correct_linkage_o = 0
linkage_with_unknown_o = 0
correct_linkage_f = 0
linkage_with_unknown_f = 0
correct_linkage_g = 0
linkage_with_unknown_g = 0
correct_linkage_s = 0
linkage_with_unknown_s = 0
unknown_mag_or_16s_taxon = 0
for each_linkage in open(linkage_file):
    if each_linkage.startswith('MarkerGene	GenomicSeq'):
        linkage_file_with_assessment_handle.write('%s\tPhylum\tClass\tOrder\tFamily\tGenus\tIdentity\tAln_len\tclosest_ref\tclosest_ref_16s_num\tMAG_taxon\t16S_taxon\t16S_taxon_BLCA\n' % each_linkage.strip())
    else:
        each_linkage_split = each_linkage.strip().split('\t')
        line_to_write = each_linkage.strip()
        s16_id = each_linkage_split[0]
        mag_id = each_linkage_split[1]
        taxon_mag = bin_taxon_dict.get(mag_id, 'NA')
        closest_ref = bin_closest_ref_dict.get(mag_id, 'NA')
        closest_ref_16s_num = len(gtdb_gnm_to_ssu_dict.get(closest_ref, []))
        taxon_16s = s16_taxon_dict.get(s16_id, 'NA')
        iden_16s = s16_iden_dict.get(s16_id, 'NA')
        aln_len_16s = s16_aln_dict.get(s16_id, 'NA')
        taxon_16s_blca = s16_taxon_blca_dict_formatted.get(s16_id, 'NA')
        taxon_mag_split = taxon_mag.split(';')
        taxon_16s_split = taxon_16s.split(';')
        total_linkage += 1

        if (taxon_mag == 'NA') or (taxon_16s == 'NA'):
            unknown_mag_or_16s_taxon += 1
            line_to_write += '\tna\tna\tna\tna\tna'
        else:
            taxon_mag_p = taxon_mag_split[1]
            taxon_mag_c = taxon_mag_split[2]
            taxon_mag_o = taxon_mag_split[3]
            taxon_mag_f = taxon_mag_split[4]
            taxon_mag_g = taxon_mag_split[5]
            taxon_mag_s = taxon_mag_split[6]
            taxon_16s_p = taxon_16s_split[1]
            taxon_16s_c = taxon_16s_split[2]
            taxon_16s_o = taxon_16s_split[3]
            taxon_16s_f = taxon_16s_split[4]
            taxon_16s_g = taxon_16s_split[5]
            taxon_16s_s = taxon_16s_split[6]
            already_wrong = False

            if taxon_mag_p == taxon_16s_p:
                correct_linkage_p += 1
                line_to_write += '\t1'
            else:
                if taxon_mag_p == 'p__':
                    linkage_with_unknown_p += 1
                    line_to_write += '\tna'
                else:
                    already_wrong = True
                    line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_c == taxon_16s_c:
                    correct_linkage_c += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_c == 'c__':
                        linkage_with_unknown_c += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_o == taxon_16s_o:
                    correct_linkage_o += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_o == 'o__':
                        linkage_with_unknown_o += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_f == taxon_16s_f:
                    correct_linkage_f += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_f == 'f__':
                        linkage_with_unknown_f += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

            if already_wrong is False:
                if taxon_mag_g == taxon_16s_g:
                    correct_linkage_g += 1
                    line_to_write += '\t1'
                else:
                    if taxon_mag_g == 'g__':
                        linkage_with_unknown_g += 1
                        line_to_write += '\tna'
                    else:
                        already_wrong = True
                        line_to_write += '\t0'
            else:
                line_to_write += '\t0'

        line_to_write += '\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (iden_16s, aln_len_16s, closest_ref, closest_ref_16s_num, taxon_mag, taxon_16s, taxon_16s_blca)
        linkage_file_with_assessment_handle.write('%s\n' % line_to_write)


if total_linkage > 0:

    print()
    print('==========================================================')
    accuracy_p = float("{0:.2f}".format((correct_linkage_p + linkage_with_unknown_p) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_c = float("{0:.2f}".format((correct_linkage_c + linkage_with_unknown_c) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_o = float("{0:.2f}".format((correct_linkage_o + linkage_with_unknown_o) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_f = float("{0:.2f}".format((correct_linkage_f + linkage_with_unknown_f) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))
    accuracy_g = float("{0:.2f}".format((correct_linkage_g + linkage_with_unknown_g) * 100 / (total_linkage - unknown_mag_or_16s_taxon)))

    print('Linkage accuracy: %s' % linkage_file_basename)
    print('Rank\tAll\tNA\tCorrect\tUnknown\tTotal(Accuracy)')
    print('phylum\t%s\t%s\t%s\t%s\t%s/%s(%s)' % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_p, linkage_with_unknown_p, (correct_linkage_p + linkage_with_unknown_p), (total_linkage - unknown_mag_or_16s_taxon), accuracy_p))
    print('class\t%s\t%s\t%s\t%s\t%s/%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_c, linkage_with_unknown_c, (correct_linkage_c + linkage_with_unknown_c), (total_linkage - unknown_mag_or_16s_taxon), accuracy_c))
    print('order\t%s\t%s\t%s\t%s\t%s/%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_o, linkage_with_unknown_o, (correct_linkage_o + linkage_with_unknown_o), (total_linkage - unknown_mag_or_16s_taxon), accuracy_o))
    print('family\t%s\t%s\t%s\t%s\t%s/%s(%s)' % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_f, linkage_with_unknown_f, (correct_linkage_f + linkage_with_unknown_f), (total_linkage - unknown_mag_or_16s_taxon), accuracy_f))
    print('genus\t%s\t%s\t%s\t%s\t%s/%s(%s)'  % (total_linkage, unknown_mag_or_16s_taxon, correct_linkage_g, linkage_with_unknown_g, (correct_linkage_g + linkage_with_unknown_g), (total_linkage - unknown_mag_or_16s_taxon), accuracy_g))

######################################################## by MAG ########################################################

mag_to_linked_16s_dict = {}
for each_linkage in open(linkage_file):
    if not each_linkage.startswith('MarkerGene	GenomicSeq	Linkage	Step'):
        each_linkage_split = each_linkage.strip().split('\t')
        mag_id = each_linkage_split[1]
        if mag_id not in mag_to_linked_16s_dict:
            mag_to_linked_16s_dict[mag_id] = {each_linkage_split[0]}
        else:
            mag_to_linked_16s_dict[mag_id].add(each_linkage_split[0])

linkage_file_with_assessment_by_mag_handle = open(linkage_file_with_assessment_by_mag, 'w')
linkage_file_with_assessment_by_mag_handle.write('MAG\tdepth\t16S\tcompleteness\tlinked\ttaxon\n')
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
    if each_mag in mag_to_linked_16s_dict:
        linked_to_16s = 'yes'

    linked_16s = mag_to_linked_16s_dict.get(each_mag, [''])
    linked_16s_len = [len_dict_16s.get(i, '0') for i in linked_16s]
    linked_16s_len_str = ';'.join([str(j) for j in linked_16s_len])

    mag_cpl = 0
    if each_mag in mag_completeness_dict:
        mag_cpl = mag_completeness_dict[each_mag]

    linkage_file_with_assessment_by_mag_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (each_mag, mag_depth, with_16s_in_mag, mag_cpl, linked_to_16s, linked_16s_len_str, mag_taxon))
linkage_file_with_assessment_by_mag_handle.close()

if total_linkage > 0:

    print('Total linkages:\t%s'  % total_linkage)
    print('Linked MAGs:\t%s'     % len(mag_to_linked_16s_dict))
    print('==========================================================')


