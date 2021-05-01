from Bio import SeqIO

# ref_gnm_id_txt = 'ref_gnm_id.txt'
#
# for each_ref in open(ref_gnm_id_txt):
#     gnm_id = each_ref.strip()
#     pwd_gnm_old = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/ref_genomes_hc_04_20/%s' % gnm_id
#     pwd_gnm_new = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/ref_genomes_hc_04_20_uniq/%s' % gnm_id
#
#     wrote_ctgs = set()
#     pwd_gnm_new_handle = open(pwd_gnm_new, 'w')
#     for each_ctg in SeqIO.parse(pwd_gnm_old, 'fasta'):
#         if each_ctg.id not in wrote_ctgs:
#             SeqIO.write(each_ctg, pwd_gnm_new_handle, 'fasta')
#             wrote_ctgs.add(each_ctg.id)
#     pwd_gnm_new_handle.close()

fa_file = '/Users/songweizhi/Desktop/NODE_345_reads.fasta'

fa_file_r1 = '/Users/songweizhi/Desktop/NODE_345_reads_R1.fasta'
fa_file_r2 = '/Users/songweizhi/Desktop/NODE_345_reads_R2.fasta'
fa_file_up = '/Users/songweizhi/Desktop/NODE_345_reads_UP.fasta'


seq_id_set = set()
for each_read in SeqIO.parse(fa_file, 'fasta'):
    seq_id_set.add(each_read.id)

print(seq_id_set)

seq_id_set_no_strand = {'.'.join(i.split('.')[:-1]) for i in seq_id_set}

print(seq_id_set_no_strand)

paired_reads = set()
for each_base_name in seq_id_set_no_strand:
    r1 = '%s.1' % each_base_name
    r2 = '%s.2' % each_base_name
    if (r1 in seq_id_set) and (r2 in seq_id_set):
        paired_reads.add(r1)
        paired_reads.add(r2)

fa_file_r1_handle = open(fa_file_r1, 'w')
fa_file_r2_handle = open(fa_file_r2, 'w')
fa_file_up_handle = open(fa_file_up, 'w')
for each_read in SeqIO.parse(fa_file, 'fasta'):
    if each_read.id in paired_reads:
        if each_read.id[-1] == '1':
            fa_file_r1_handle.write('>%s\n' % each_read.id)
            fa_file_r1_handle.write('%s\n' % each_read.seq)
        if each_read.id[-1] == '2':
            fa_file_r2_handle.write('>%s\n' % each_read.id)
            fa_file_r2_handle.write('%s\n' % each_read.seq)
    else:
        fa_file_up_handle.write('>%s\n' % each_read.id)
        fa_file_up_handle.write('%s\n' % each_read.seq)
fa_file_r1_handle.close()
fa_file_r2_handle.close()
fa_file_up_handle.close()