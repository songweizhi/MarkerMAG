import os
from Bio import SeqIO

# wd = '/Users/songweizhi/Desktop/lc'
# lc_30_profile               = '/Users/songweizhi/Desktop/lc/30_profile.tsv'
# complete_fna_gi_acc_tax_tsv = '/Users/songweizhi/Desktop/lc/complete.fna.gi_acc_tax.tsv'
# bacteria_fna                = 'bacteria.fna'
# ref_gnm_folder = '%s/ref_genomes' % wd

complete_fna_gi_acc_tax_tsv         = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/refseq/complete.fna.gi_acc_tax.tsv'
bacteria_fna_gi_acc_tax_tsv         = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/genomes/bacteria/bacteria.fna.gi_acc_tax.tsv'
bacteria_draft_fna_gi_acc_tax_tsv   = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/genomes/bacteria_draft/bacteria_draft.fna.gi_acc_tax.tsv'
complete_fna                        = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/refseq/complete.fna'
bacteria_fna                        = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/genomes/bacteria/bacteria.fna'
bacteria_draft_fna                  = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/genomes/bacteria_draft/bacteria_draft.fna'


lc_30_profile                       = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/2_Medium_Complexity/225_profile_d1_180bp.tsv'
lc_30_profile                       = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/1_Low_Complexity/30_profile.tsv'
lc_30_profile                       = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/450_profile_d1.tsv'
ref_gnm_folder                      = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/ref_genomes_hc_04_20'

'''

grep 'NZ_KB912814.1' /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/genomes/bacteria/bacteria.fna.id.lines.txt
grep 'NZ_KB912814.1' /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/genomes/bacteria_draft/bacteria_draft.fna.id.lines.txt
grep 'NZ_KB912814.1' /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/PROCESSED_NCBI/refseq/complete.fna.id.txt

'''


missing_tax_lc = ['1131272.1', '1230476.1']

os.mkdir(ref_gnm_folder)
if os.path.isdir(ref_gnm_folder) is True:
    os.system('rm -r %s' % ref_gnm_folder)
os.mkdir(ref_gnm_folder)

select_ref_tax_id_set = set()
for each_line in open(lc_30_profile):
    each_line_split = each_line.strip().split('\t')
    if len(each_line_split) > 1:
        if each_line_split[1] == 'strain':
                select_ref_tax_id_set.add(each_line_split[0])

select_ref_tax_id_set_no_version = {i.split('.')[0] for i in select_ref_tax_id_set}


#################### get gi_to_tax_dict ####################

gi_to_tax_dict = {}

for each_line in open(complete_fna_gi_acc_tax_tsv):
    each_line_split = each_line.strip().split('\t')
    id_gi  = each_line_split[0]
    id_acc = each_line_split[1]
    id_tax = each_line_split[2]
    if id_tax in select_ref_tax_id_set_no_version:
        gi_to_tax_dict[id_gi] = id_tax

for each_line in open(bacteria_fna_gi_acc_tax_tsv):
    each_line_split = each_line.strip().split('\t')
    id_gi  = each_line_split[0]
    id_acc = each_line_split[1]
    id_tax = each_line_split[2]
    if id_tax in select_ref_tax_id_set_no_version:
        gi_to_tax_dict[id_gi] = id_tax

for each_line in open(bacteria_draft_fna_gi_acc_tax_tsv):
    each_line_split = each_line.strip().split('\t')
    id_gi  = each_line_split[0]
    id_acc = each_line_split[1]
    id_tax = each_line_split[2]
    if id_tax in select_ref_tax_id_set_no_version:
        gi_to_tax_dict[id_gi] = id_tax


############################################################

for each_seq in SeqIO.parse(bacteria_fna, 'fasta'):
    gi_id = each_seq.id.split('|')[1]
    if gi_id in gi_to_tax_dict:
        tax_id = gi_to_tax_dict[gi_id]
        pwd_gnm_file = '%s/tax_%s.fa' % (ref_gnm_folder, tax_id)
        with open(pwd_gnm_file, 'a') as gnm_file_handle:
            gnm_file_handle.write('>%s\n' % each_seq.id)
            gnm_file_handle.write('%s\n' % str(each_seq.seq))

for each_seq in SeqIO.parse(bacteria_draft_fna, 'fasta'):
    gi_id = each_seq.id.split('|')[1]
    if gi_id in gi_to_tax_dict:
        tax_id = gi_to_tax_dict[gi_id]
        pwd_gnm_file = '%s/tax_%s.fa' % (ref_gnm_folder, tax_id)
        with open(pwd_gnm_file, 'a') as gnm_file_handle:
            gnm_file_handle.write('>%s\n' % each_seq.id)
            gnm_file_handle.write('%s\n' % str(each_seq.seq))

for each_seq in SeqIO.parse(complete_fna, 'fasta'):
    gi_id = each_seq.id.split('|')[1]
    if gi_id in gi_to_tax_dict:
        tax_id = gi_to_tax_dict[gi_id]
        pwd_gnm_file = '%s/tax_%s.fa' % (ref_gnm_folder, tax_id)
        with open(pwd_gnm_file, 'a') as gnm_file_handle:
            gnm_file_handle.write('>%s\n' % each_seq.id)
            gnm_file_handle.write('%s\n' % str(each_seq.seq))

