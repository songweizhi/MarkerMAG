import os
import glob
from Bio import SeqIO


genomic_assemblies          = ''
renamed_mag_folder          = '/Users/songweizhi/Desktop/777/Refined_refined_bins_renamed_renamed'
mag_file_extension          = 'fna'
marker_gene_seqs            = '/Users/songweizhi/Desktop/777/MBARC26_all_depth_assemblies.dereplicated_99.5_500bp.fasta'
link_stats_combined_table   = '/Users/songweizhi/Desktop/777/MBARC26_link_Matam_RealBins_500bp_depth_d0_mpl10_stats_combined_table.txt'
genomic_seq_type            = 'mag'

marker_gene_seqs_1st_round_unlinked = '/Users/songweizhi/Desktop/777/MBARC26_all_depth_assemblies.dereplicated_99.5_500bp_1st_round_unlinked.fasta'
combined_1st_round_unlinked_mags    = '/Users/songweizhi/Desktop/777/combined_1st_round_unlinked_mags.fasta'
combined_1st_round_unlinked_ctgs    = '/Users/songweizhi/Desktop/777/combined_1st_round_unlinked_ctgs.fasta'


# get linked marker genes and genomic sequences in 1st round
linked_marker_gene_set = set()
linked_genomic_seq_set = set()
for each_link in open(link_stats_combined_table):
    if not each_link.startswith('MarkerGene	GenomicSeq	Paired	Clipping'):
        each_link_split = each_link.strip().split('\t')
        linked_marker_gene_set.add(each_link_split[0])
        linked_genomic_seq_set.add(each_link_split[1])


# get the sequence of unlinked marker genes
marker_gene_seqs_1st_round_unlinked_handle = open(marker_gene_seqs_1st_round_unlinked, 'w')
for marker_gene_record in SeqIO.parse(marker_gene_seqs, 'fasta'):
    if marker_gene_record.id not in linked_marker_gene_set:
        marker_gene_seqs_1st_round_unlinked_handle.write('>%s\n' % marker_gene_record.id)
        marker_gene_seqs_1st_round_unlinked_handle.write('%s\n' % marker_gene_record.seq)
marker_gene_seqs_1st_round_unlinked_handle.close()


# get the sequence of unlinked genomic seqs
if genomic_seq_type == 'mag':

    # put all renamed mag into list
    renamed_gnm_re   = '%s/*.%s' % (renamed_mag_folder, mag_file_extension)
    renamed_gnm_list = [os.path.basename(file_name) for file_name in glob.glob(renamed_gnm_re)]
    renamed_gnm_list_no_ext = ['.'.join(i.split('.')[:-1]) for i in renamed_gnm_list]

    # keep only unlinked mags
    unlinked_mag_list_with_pwd = []
    for renamed_mag in renamed_gnm_list_no_ext:
        if renamed_mag not in linked_genomic_seq_set:
            pwd_renamed_mag = '%s/%s.%s' % (renamed_mag_folder, renamed_mag, mag_file_extension)
            unlinked_mag_list_with_pwd.append(pwd_renamed_mag)

    # combine unlinked mags
    cat_cmd = 'cat %s > %s' % (' '.join(unlinked_mag_list_with_pwd), combined_1st_round_unlinked_mags)
    os.system(cat_cmd)


# get the sequence of unlinked metagenomic assemblies
if genomic_seq_type == 'ctg':
    combined_1st_round_unlinked_ctgs_handle = open(combined_1st_round_unlinked_ctgs, 'w')
    for ctg_record in SeqIO.parse(genomic_assemblies, 'fasta'):
        if ctg_record.id not in linked_genomic_seq_set:
            SeqIO.write(ctg_record, combined_1st_round_unlinked_ctgs_handle, 'fasta')
    combined_1st_round_unlinked_ctgs_handle.close()

