from Bio import SeqIO

fa_file = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/CAMI2_GI_mplu5_longkmer_55-127_2nd_MarkerMAG_wd/CAMI2_GI_mplu5_longkmer_55-127_2nd_step_2_wd/free_living_read_combined.fasta'
fq_r1 = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/1_filtered_reads/GI_R1.fastq'
fq_r2 = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/1_filtered_reads/GI_R2.fastq'
fq_out = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI2_HMP/CAMI2_GI_mplu5_longkmer_55-127_2nd_MarkerMAG_wd/CAMI2_GI_mplu5_longkmer_55-127_2nd_step_2_wd/free_living_read_combined.fastq'


reads_to_extract = set()
for fa_read_id in SeqIO.parse(fa_file, 'fasta'):
    reads_to_extract.add(fa_read_id.id)

fq_out_handle = open(fq_out, 'w')
for r1_record in SeqIO.parse(fq_r1, 'fastq'):
    if r1_record.id in reads_to_extract:
        SeqIO.write(r1_record, fq_out_handle, 'fastq')
for r2_record in SeqIO.parse(fq_r2, 'fastq'):
    if r2_record.id in reads_to_extract:
        SeqIO.write(r2_record, fq_out_handle, 'fastq')
fq_out_handle.close()


