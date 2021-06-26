import os
import glob
from Bio import SeqIO


def get_accuracy(file_in, marker_num):

    linkage_num_total = 0
    linkage_num_correct = 0
    recovered_markers = set()
    for each_match in open(file_in):
        if not each_match.startswith('MarkerGene\tGenomicSeq\tLinkage\tStep'):
            match_split = each_match.strip().split('\t')
            linkage_num = int(match_split[2])
            MarkerGene_genome = match_split[0].split('_')[0]
            GenomicSeq_genome = match_split[1]
            linkage_num_total += linkage_num
            if MarkerGene_genome == GenomicSeq_genome:
                linkage_num_correct += linkage_num
                recovered_markers.add(match_split[0])

    marker_recovery = float("{0:.2f}".format(len(recovered_markers)*100/marker_num))

    link_accuracy = 0
    if linkage_num_total > 0:
        link_accuracy = float("{0:.2f}".format(linkage_num_correct*100/linkage_num_total))
        print('%s/%s' % (linkage_num_correct, linkage_num_total))

    marker_recovery = '%s/%s(%s)' % (len(recovered_markers), marker_num, marker_recovery)

    print('marker_recovery')
    print(marker_recovery)

    return marker_recovery, link_accuracy, recovered_markers


def get_accuracy_by_genome(file_in, mag_folder, mag_file_extension):

    # get MAG file list
    mag_file_re             = '%s/*%s' % (mag_folder, mag_file_extension)
    mag_file_list           = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
    mag_file_list_no_ext    = {'.'.join(i.split('.')[:-1]) for i in mag_file_list}

    genome_with_right_16s_assignment_tmp = set()
    genome_with_wrong_16s_assignment = set()
    for each_match in open(file_in):
        if not each_match.startswith('MarkerGene\tGenomicSeq\tLinkage\tStep'):
            match_split = each_match.strip().split('\t')
            MarkerGene_genome = match_split[0].split('_')[0]
            GenomicSeq_genome = match_split[1]

            if GenomicSeq_genome == MarkerGene_genome:
                genome_with_right_16s_assignment_tmp.add(GenomicSeq_genome)
            else:
                genome_with_wrong_16s_assignment.add(GenomicSeq_genome)


    genome_with_right_16s_assignment_always = []
    genome_without_right_16s_assignment = []
    for input_genome in mag_file_list_no_ext:
        if (input_genome in genome_with_right_16s_assignment_tmp) and (input_genome not in genome_with_wrong_16s_assignment):
            genome_with_right_16s_assignment_always.append(input_genome)
        else:
            genome_without_right_16s_assignment.append(input_genome)


    marker_gene_assignment_rate = float("{0:.2f}".format(len(genome_with_right_16s_assignment_always)*100/len(mag_file_list_no_ext)))

    marker_gene_assignment_accuracy = 0
    if (len(genome_with_right_16s_assignment_always) + len(genome_with_wrong_16s_assignment)) > 0:
        marker_gene_assignment_accuracy = float("{0:.2f}".format(len(genome_with_right_16s_assignment_always)*100/(len(genome_with_right_16s_assignment_always) + len(genome_with_wrong_16s_assignment))))
    marker_gene_assignment_rate = '%s/%s(%s)' % (len(genome_with_right_16s_assignment_always), len(mag_file_list_no_ext), marker_gene_assignment_rate)

    print('marker_gene_assignment_accuracy')
    print('%s/%s' % ())
    return marker_gene_assignment_rate, marker_gene_assignment_accuracy, genome_with_right_16s_assignment_always, genome_without_right_16s_assignment


def get_unrecovered_markers(marker_all, marker_recovered):
    unrecovered_markers = []
    for each_marker in marker_all:
        if each_marker not in marker_recovered:
            unrecovered_markers.append(each_marker)
    return sorted(unrecovered_markers)


########################################################################################################################

wd                          = '/Users/songweizhi/Desktop/MarkerMAG_wd/2_MBARC26'
marker_gene_seqs            = '%s/MBARC26_all_depth_assemblies_uclust_0.999_99.9_500bp.fasta'  % wd
mag_folder                  = '%s/Refined_refined_bins_renamed'                     % wd
mag_file_extension          = 'fna'
combined_linkage_file       = '%s/MBARC26_mis_4_combined_linkages.txt' % wd


########################################################################################################################

marker_id_set = set()
for marker_seq_record in SeqIO.parse(marker_gene_seqs, 'fasta'):
    marker_id_set.add(marker_seq_record.id)

# get recovery and accuracy
recovery_combined, accuracy_combined, recovered_combined = get_accuracy(combined_linkage_file, len(marker_id_set))

# get unrecovered markers
unrecovered_markers_paired = get_unrecovered_markers(marker_id_set, recovered_combined)
unrecovered_markers_paired_str = 'Unrecovered(%s):%s' % (len(unrecovered_markers_paired), ','.join(sorted([i for i in unrecovered_markers_paired])))

# assessment by genome
assign_rate, assign_accuracy, right_assign, wrong_assign = get_accuracy_by_genome(combined_linkage_file, mag_folder, mag_file_extension)
unrecovered_paired_report_str = 'Unrecovered(%s):%s' % (len(wrong_assign), ','.join(sorted([i for i in wrong_assign])))

# report
print('By\tRecovery\tAccuracy\tUnrecovered')
print('Marker\t%s\t%s\t%s' % (recovery_combined, accuracy_combined, unrecovered_markers_paired_str))
print('Genome\t%s\t%s\t%s' % (assign_rate, assign_accuracy, unrecovered_paired_report_str))
