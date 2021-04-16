import os
from Bio import SeqIO


wd                  = '/Users/songweizhi/Desktop/CAMI_1'
combined_gnm_file   = '%s/S_S001__genomes_30__insert_180_gsa_anonymous.fasta'   % wd
ctg_to_gnm_file     = '%s/S_S001__genomes_30__insert_180_gsa_mapping.tsv'       % wd

ref_gnm_folder = '%s/ref_genomes' % wd


os.mkdir(ref_gnm_folder)
if os.path.isdir(ref_gnm_folder) is True:
    os.system('rm -r %s' % ref_gnm_folder)
os.mkdir(ref_gnm_folder)


ctg_to_gnm_dict = {}
ctg_to_seq_id = {}
for each_ctg in open(ctg_to_gnm_file):
    each_ctg_split = each_ctg.strip().split('\t')
    ctg_id = each_ctg_split[0]
    gnm_id = each_ctg_split[1]
    ctg_id_raw = each_ctg_split[3]
    ctg_to_seq_id[ctg_id] = ctg_id_raw
    ctg_to_gnm_dict[ctg_id] = gnm_id


for each_seq in SeqIO.parse(combined_gnm_file, 'fasta'):

    seq_id = each_seq.id
    seq_id_raw = ctg_to_seq_id[seq_id]
    seq_gnm = ctg_to_gnm_dict[seq_id]

    pwd_gnm_file = '%s/%s.fa' % (ref_gnm_folder, seq_gnm)
    print(pwd_gnm_file)
    with open(pwd_gnm_file, 'a') as gnm_file_handle:
        gnm_file_handle.write('>%s\n' % seq_id_raw)
        gnm_file_handle.write('%s\n' % str(each_seq.seq))



