import os
import pandas as pd
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def sort_csv_by_col(file_in, file_out, col_header):
    df_in        = pd.read_csv(file_in)
    df_in_sorted = df_in.sort_values(by=[col_header], ascending=False)
    df_in_sorted.to_csv(file_out, index=False)


def filter_linkages_iteratively(file_in, sort_by_col_header, file_out):

    file_in_path, file_in_basename, file_in_extension = sep_path_basename_ext(file_in)
    file_in_sorted = '%s/%s_sorted%s' % (file_in_path, file_in_basename, file_in_extension)

    # sort file in
    sort_csv_by_col(file_in, file_in_sorted, sort_by_col_header)

    # get GenomicSeq_max_num_dict
    GenomicSeq_max_num_dict = {}
    for each_linkage in open(file_in):
        if not each_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            linkage_split = each_linkage.strip().split(',')
            if linkage_split[1][12:] not in GenomicSeq_max_num_dict:
                GenomicSeq_max_num_dict[linkage_split[1][12:]] = int(linkage_split[2])
            else:
                if int(linkage_split[2]) > GenomicSeq_max_num_dict[linkage_split[1][12:]]:
                    GenomicSeq_max_num_dict[linkage_split[1][12:]] = int(linkage_split[2])

    # fileter linkage
    file_out_handle = open(file_out, 'w')
    MarkerGene_with_assignment = set()
    for each_match in open(file_in_sorted):
        if each_match.startswith('MarkerGene,GenomicSeq,Number'):
            file_out_handle.write(each_match)
        else:
            match_split = each_match.strip().split(',')
            MarkerGene = match_split[0][12:]
            if MarkerGene not in MarkerGene_with_assignment:
                MarkerGene_with_assignment.add(MarkerGene)
                file_out_handle.write(each_match)
    file_out_handle.close()

    # remove tmp file
    os.remove(file_in_sorted)


def combine_paired_and_clipping_linkages(paired_linkages, clipping_linkages, file_out_summary, file_out_intersect_linkages):

    # file in:   file_in_paired    and  file_in_clipping
    # file out:  file_out_summary  and  file_out_intersection

    combined_paired_and_clipping_keys = set()

    # read in paired linkages
    paired_linkages_dict = {}
    for paired_linkage in open(paired_linkages):
        if not paired_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            paired_linkage_split = paired_linkage.strip().split(',')
            paired_key = '%s__|__%s' % (paired_linkage_split[0], paired_linkage_split[1])
            paired_value = int(paired_linkage_split[2])
            paired_linkages_dict[paired_key] = paired_value
            combined_paired_and_clipping_keys.add(paired_key)

    # read in clipping linkages
    clipping_linkages_dict = {}
    for clipping_linkage in open(clipping_linkages):
        if not clipping_linkage.startswith('MarkerGene,GenomicSeq,Number'):
            clipping_linkage_split = clipping_linkage.strip().split(',')
            clipping_key = '%s__|__%s' % (clipping_linkage_split[0], clipping_linkage_split[1])
            clipping_value = int(clipping_linkage_split[2])
            clipping_linkages_dict[clipping_key] = clipping_value
            combined_paired_and_clipping_keys.add(clipping_key)

    # combine paired and clipping linkages
    file_out_summary_handle = open(file_out_summary, 'w')
    file_out_intersect_linkages_handle = open(file_out_intersect_linkages, 'w')
    file_out_summary_handle.write('MarkerGene\tGenomicSeq\tPaired\tClipping\tIntersection\n')
    file_out_intersect_linkages_handle.write('MarkerGene,GenomicSeq,Number\n')
    for each_key in combined_paired_and_clipping_keys:

        current_key_paired_value = 0
        if each_key in paired_linkages_dict:
            current_key_paired_value = paired_linkages_dict[each_key]

        current_key_clipping_value = 0
        if each_key in clipping_linkages_dict:
            current_key_clipping_value = clipping_linkages_dict[each_key]

        current_key_combined = 0
        intersection = 'no'
        if (current_key_paired_value > 0) and (current_key_clipping_value > 0):
            current_key_combined = current_key_paired_value + current_key_clipping_value
            intersection = 'yes'

        # write out
        file_out_summary_handle.write('%s\t%s\t%s\t%s\n' % ('\t'.join(each_key.split('__|__')), current_key_paired_value, current_key_clipping_value, intersection))
        if current_key_combined > 0:
            file_out_intersect_linkages_handle.write('%s,%s\n' % (','.join(each_key.split('__|__')) , current_key_combined))

    file_out_summary_handle.close()
    file_out_intersect_linkages_handle.close()


def get_accuracy(file_in, marker_num):

    linkage_num_total = 0
    linkage_num_correct = 0
    recovered_markers = set()
    for each_match in open(file_in):
        if not each_match.startswith('MarkerGene,GenomicSeq,Number'):
            match_split = each_match.strip().split(',')
            linkage_num = int(match_split[2])
            MarkerGene_genome = match_split[0][12:][:2]
            GenomicSeq_genome = match_split[1][12:]

            linkage_num_total += linkage_num
            if MarkerGene_genome == GenomicSeq_genome:
                linkage_num_correct += linkage_num
                recovered_markers.add(match_split[0][12:])

    marker_recovery = float("{0:.4f}".format(len(recovered_markers)/marker_num))
    link_accuracy = float("{0:.4f}".format(linkage_num_correct/linkage_num_total))

    return marker_recovery, link_accuracy


# file in
wd                      = '/Users/songweizhi/Desktop/test_in'
marker_seq_file         = '%s/combined_16S.ffn'                                                  % wd
file_in_paired          = '%s/EvenDepth25x_stats_M30_S30_iden100_cov100_n1_paired.txt'           % wd
file_in_clipping        = '%s/EvenDepth25x_stats_M30_S30_iden100_cov100_n1_clipping.txt'         % wd

# file out
file_out_combined_table = '%s/EvenDepth25x_stats_M30_S30_iden100_cov100_n1_combined_table.txt'   % wd
file_out_intersection   = '%s/EvenDepth25x_stats_M30_S30_iden100_cov100_n1_intersection.txt'     % wd


# get the number of marker genes
marker_num = 0
for marker_seq_record in SeqIO.parse(marker_seq_file, 'fasta'):
    marker_num += 1

# define output file name
file_in_paired_path, file_in_paired_basename, file_in_paired_extension = sep_path_basename_ext(file_in_paired)
file_in_clipping_path, file_in_clipping_basename, file_in_clipping_extension = sep_path_basename_ext(file_in_clipping)
file_out_paired_filtered   = '%s/%s_filtered.txt' % (file_in_paired_path, file_in_paired_basename)
file_out_clipping_filtered = '%s/%s_filtered.txt' % (file_in_clipping_path, file_in_clipping_basename)

# filter paired and clipping linkages
filter_linkages_iteratively(file_in_paired, 'Number', file_out_paired_filtered)
filter_linkages_iteratively(file_in_clipping, 'Number', file_out_clipping_filtered)

# combine_paired_and_clipping_linkages and get summary table
combine_paired_and_clipping_linkages(file_out_paired_filtered, file_out_clipping_filtered, file_out_combined_table, file_out_intersection)

# get recovery and accuracy
marker_recovery_paired, link_accuracy_paired             = get_accuracy(file_out_paired_filtered, marker_num)
marker_recovery_clipping, link_accuracy_clipping         = get_accuracy(file_out_clipping_filtered, marker_num)
marker_recovery_intersection, link_accuracy_intersection = get_accuracy(file_out_intersection, marker_num)


print('sample\trec_p\tacc_p\trec_c\tacc_c\trec_i\tacc_i')
print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('Sample', marker_recovery_paired, link_accuracy_paired, marker_recovery_clipping, link_accuracy_clipping, marker_recovery_intersection, link_accuracy_intersection))


#get_sanky_plot_cmd = 'Rscript ~/PycharmProjects/LinkMAG16S/get_sankey_plot.R -f %s -x 500 -y 800' % file_out
#os.system(get_sanky_plot_cmd)


'''
the number of linkage either via mate reads or clipping mapped reads were
between each pair of genomic sequences and 16S rRNA gene sequences were sorted according to the number of linking reads (either paired reads or clipping mapped reads) between them in descending order.
16S rRNA gene sequences already linked to other genomic sequences were ignored in later interations.


Recovery: the percentage of correctly recovered introduced 16S rRNA gene sequences.
Accuracy: the percentage of correct linkages in all identified linkage.


# with 16S iden cutoff 99.5
Depth	rec_p	acc_p	rec_c	acc_c	rec_i	acc_i
1x	    0.6897	0.9947	0.5402	0.9686	0.4138	1.0
2x		0.908	1.0		0.7241	0.9932	0.6667	1.0
3x	    0.8966	0.9716	0.8391	1.0		0.7586	1.0
5x	    0.954	1.0		0.9425	0.9963	0.931	1.0
25x	    0.954	1.0		0.9425	0.9948	0.931	1.0
50x	    0.954	1.0		0.9655	1.0		0.954	1.0
100x	0.9425	1.0		0.9425	0.9951	0.8966	1.0


# with 16S iden cutoff 99
Depth	rec_p	acc_p	rec_c	acc_c	rec_i	acc_i
1x	    0.7241	1.0	    0.5862	0.9709	0.4598	1.0
2x		0.9655	1.0	    0.7586	0.9933	0.7356	1.0
3x	    0.931	0.9724	0.8851	1.0 	0.8161	1.0
5x	    1.0	    1.0	    0.977	0.9963	0.977	1.0
25x	    1.0	    1.0	    0.9885	0.9969	0.9885	1.0
50x	    1.0	    1.0	    0.9885	0.9983	0.9885	1.0
100x	1.0	    1.0	    0.9885	0.9973	0.9885	1.0


'''




################################################# rank level assessment ################################################

# seq_id_file_16S    = '%s/16S_seq_id.txt'                                              % wd

# rank_to_16s_num_dict = {'c': 0, 'f': 0, 'g': 0, 'o': 0, 'p': 0, 's': 0}
# for seq_id in open(seq_id_file_16S):
#     if seq_id.startswith('s'):
#         rank_to_16s_num_dict['s'] += 1
#     if seq_id.startswith('g'):
#         rank_to_16s_num_dict['g'] += 1
#     if seq_id.startswith('f'):
#         rank_to_16s_num_dict['f'] += 1
#     if seq_id.startswith('o'):
#         rank_to_16s_num_dict['o'] += 1
#     if seq_id.startswith('c'):
#         rank_to_16s_num_dict['c'] += 1
#     if seq_id.startswith('p'):
#         rank_to_16s_num_dict['p'] += 1
#
# linkage_num_total   =  {'c': 0, 'f': 0, 'g': 0, 'o': 0, 'p': 0, 's': 0}
# linkage_num_correct =  {'c': 0, 'f': 0, 'g': 0, 'o': 0, 'p': 0, 's': 0}
# recovered_markers = {'c': set(), 'f': set(), 'g': set(), 'o': set(), 'p': set(), 's': set()}
# for each_match in open(file_out):
#     if not each_match.startswith('MarkerGene,GenomicSeq,Number'):
#         match_split = each_match.strip().split(',')
#         linkage_num = int(match_split[2])
#         MarkerGene_genome = match_split[0][12:][:2]
#         GenomicSeq_genome = match_split[1][12:]
#         linkage_num_total[MarkerGene_genome[0]] += 1
#         if MarkerGene_genome == GenomicSeq_genome:
#             linkage_num_correct[MarkerGene_genome[0]] += 1
#             recovered_markers[MarkerGene_genome[0]].add(match_split[0][12:])
#
# print('Rank\tRecovery\tAccuracy')
# for each_rank in ['s', 'g', 'f', 'o', 'c', 'p'][::-1]:
#     current_rank_linkage_num_total     = linkage_num_total[each_rank]
#     current_rank_linkage_num_correct   = linkage_num_correct[each_rank]
#     current_rank_introduced_marker_num = rank_to_16s_num_dict[each_rank]
#     current_rank_recovered_markers     = recovered_markers[each_rank]
#     current_rank_accuracy              = float("{0:.4f}".format(current_rank_linkage_num_correct/current_rank_linkage_num_total))
#     current_rank__recovery             = float("{0:.4f}".format(len(current_rank_recovered_markers)/current_rank_introduced_marker_num))
#     print('%s\t%s\t%s\t%s' % (file_in, each_rank, current_rank__recovery, current_rank_accuracy))
#

########################################################################################################################

