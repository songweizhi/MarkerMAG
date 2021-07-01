from Bio import SeqIO


def prefix_seq(seq_in, prefix, seq_out):

    seq_out_handle = open(seq_out, 'w')
    for seq_record in SeqIO.parse(seq_in, 'fasta'):
        seq_id_new = '%s_%s' % (prefix, seq_record.id)
        seq_out_handle.write('>%s\n' % seq_id_new)
        seq_out_handle.write('%s\n' % str(seq_record.seq))
    seq_out_handle.close()


matam_assemblies            = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_Matam16S_wd/CAMI_Oral_16S_reads_subset_50_Matam_wd/workdir/scaffolds.NR.min_500bp.fa'
matam_assemblies_prefixed   = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/CAMI_Oral_Matam16S_wd/CAMI_Oral_16S_reads_subset_50_Matam_wd/workdir/scaffolds.NR.min_500bp.prefixed.fa'
seq_prefix                  = '%s_subsample_%s' % ('CAMI_Oral', 50)

prefix_seq(matam_assemblies, seq_prefix, matam_assemblies_prefixed)

