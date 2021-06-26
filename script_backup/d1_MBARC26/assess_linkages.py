import os
import glob
from Bio import SeqIO

'''


MBARC26_mis_1_10_5_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	    Unlinked:AS,CG,CP,HR,PS,TR
MBARC26_mis_2_10_5_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	    Unlinked:AS,CG,CP,HR,PS,TR
MBARC26_mis_3_10_5_combined_linkages.txt	Genome	16/23(69.57)	16/16(100.0)	Wrong:	    Unlinked:AS,CG,CP,EC,HR,PS,TR
MBARC26_mis_4_10_5_combined_linkages.txt	Genome	18/23(78.26)	18/18(100.0)	Wrong:	    Unlinked:CG,CP,HR,PS,TR
MBARC26_mis_5_10_5_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	    Unlinked:AS,CG,CP,HR,PS,TR
MBARC26_mis_6_10_5_combined_linkages.txt	Genome	16/23(69.57)	16/17(94.12)	Wrong:EC	Unlinked:AS,CG,CP,HR,PS,TR
MBARC26_mis_7_10_5_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	    Unlinked:CG,CP,EC,HR,PS,TR
MBARC26_mis_8_10_5_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	    Unlinked:CG,CP,EC,HR,PS,TR
MBARC26_mis_9_10_5_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	    Unlinked:CG,CP,EC,HR,PS,TR

MBARC26_mis_3_10_5_spades_combined_linkages.txt	Genome	17/23(73.91)	17/17(100.0)	Wrong:	Unlinked:AS,CG,CP,EC,PS,TR
MBARC26_mis_4_10_5_spades_combined_linkages.txt	Genome	18/23(78.26)	18/18(100.0)	Wrong:	Unlinked:CG,CP,EC,PS,TR
MBARC26_mis_5_10_5_spades_combined_linkages.txt	Genome	18/23(78.26)	18/18(100.0)	Wrong:	Unlinked:CG,CP,EC,PS,TR

'''

########################################################################################################################

wd                          = '/Users/songweizhi/Desktop/MarkerMAG_wd/2_MBARC26'
marker_gene_seqs            = '%s/MBARC26_all_depth_assemblies_uclust_0.999_99.9_500bp.fasta'   % wd
mag_folder                  = '%s/Refined_refined_bins_renamed'                                 % wd
mag_file_extension          = 'fna'
combined_linkage_file       = '%s/MBARC26_mis_6_10_5_combined_linkages.txt'                          % wd

'''

# SB and AS formed one bin,

[2021-04-07 22:01:23] MBARC26_0407_very_sensitive	Genome	19/23(82.61)	90.48	Unrecovered(4):AS,CP,CA,PS
[2021-04-07 21:49:24] MBARC26_0407_sensitive	    Genome	19/23(82.61)	90.48	Unrecovered(4):AS,CP,CA,PS
[2021-04-07 21:58:59] MBARC26_0407_default	        Genome	21/23(91.3)	    95.45	Unrecovered(2):AS,CP
[2021-04-07 21:53:56] MBARC26_0407_specific	        Genome	21/23(91.3)	    95.45	Unrecovered(2):AS,CP
[2021-04-07 21:50:10] MBARC26_0407_very_specific	Genome	20/23(86.96)	95.24	Unrecovered(3):AS,CP,TR
[2021-04-07 21:54:39] MBARC26_0407_super_specific	Genome	19/23(82.61)	95.0	Unrecovered(4):AS,CP,HR,TR

'''

########################################################################################################################

marker_id_set = set()
for marker_seq_record in SeqIO.parse(marker_gene_seqs, 'fasta'):
    marker_id_set.add(marker_seq_record.id)

linkage_file_basename = os.path.basename(combined_linkage_file)

mag_file_re             = '%s/*%s' % (mag_folder, mag_file_extension)
mag_file_list           = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
mag_file_list_no_ext    = {'.'.join(i.split('.')[:-1]) for i in mag_file_list}

markers_with_wrong_linkage = set()
markers_with_right_linkage = set()
genomes_with_right_linkage = set()
genomes_with_wrong_linkage = set()
all_linked_genomes = set()
total_linkages = 0
right_linkages = 0
for each_match in open(combined_linkage_file):
    if not each_match.startswith('MarkerGene\tGenomicSeq\tLinkage\tStep'):
        total_linkages += 1
        match_split = each_match.strip().split('\t')
        linkage_num = int(match_split[2])
        MarkerGene_genome = match_split[0].split('_')[0]
        GenomicSeq_genome = match_split[1]
        all_linked_genomes.add(GenomicSeq_genome)
        if MarkerGene_genome == GenomicSeq_genome:
            right_linkages += 1
            genomes_with_right_linkage.add(GenomicSeq_genome)
            markers_with_right_linkage.add(match_split[0])
        else:
            genomes_with_wrong_linkage.add(GenomicSeq_genome)
            markers_with_wrong_linkage.add(match_split[0])

recovery_by_marker_str      = '%s/%s(%s)'               % (len(markers_with_right_linkage), len(marker_id_set), float("{0:.2f}".format(len(markers_with_right_linkage)*100/len(marker_id_set))))
accuracy_by_marker_str      = '%s/%s(%s)'               % (right_linkages, total_linkages, float("{0:.2f}".format(right_linkages*100/total_linkages)))
unrecovered_by_marker_str   = 'Wrong:%s\tUnlinked:'     % (','.join(markers_with_wrong_linkage))
recovery_by_genome_str      = '%s/%s(%s)'               % (len(genomes_with_right_linkage), len(mag_file_list_no_ext), float("{0:.2f}".format(len(genomes_with_right_linkage)*100/len(mag_file_list_no_ext))))
accuracy_by_genome_str      = '%s/%s(%s)'               % (len(genomes_with_right_linkage), len(all_linked_genomes), float("{0:.2f}".format(len(genomes_with_right_linkage)*100/len(all_linked_genomes))))
unrecovered_by_genome_str   = 'Wrong:%s\tUnlinked:%s'   % (','.join(genomes_with_wrong_linkage), ','.join(sorted([i for i in (mag_file_list_no_ext - genomes_with_right_linkage - genomes_with_wrong_linkage)])))

print('Linkages\tBy\tRecovery\tAccuracy\tUnrecovered')
print('%s\tMarker\t%s\t%s\t%s' % (linkage_file_basename, recovery_by_marker_str, accuracy_by_marker_str, unrecovered_by_marker_str))
print('%s\tGenome\t%s\t%s\t%s' % (linkage_file_basename, recovery_by_genome_str, accuracy_by_genome_str, unrecovered_by_genome_str))
