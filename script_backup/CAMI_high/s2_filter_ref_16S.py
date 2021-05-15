import os
import glob
from Bio import SeqIO


wd                  = '/Users/songweizhi/Desktop/CAMI_high'
gnm_folder          = '%s/ref_gnm_596'                                                  % wd
ffn_folder          = '%s/ref_gnm_596_barrnap_outputs'                                  % wd
ref_596_genomes     = '%s/goldstandard_high_pool.filtered.profile.stain.596_ref.txt'    % wd
len_cutoff_16s      = 1000
len_cutoff_ctg      = 5000


qualified_16s_seq    = '%s/ref_16S_%sbp_%sbp.ffn'         % (wd, len_cutoff_16s, len_cutoff_ctg)
qualified_gnm_folder = '%s/ref_with_16S_%sbp_%sbp'        % (wd, len_cutoff_16s, len_cutoff_ctg)


genomes_to_tax_dict = {}
for each_ref in open(ref_596_genomes):
    each_ref_split = each_ref.strip().split('\t')
    tax_id = each_ref_split[0]
    cami_genome_id = each_ref_split[5]
    genomes_to_tax_dict[cami_genome_id] = tax_id


if os.path.isdir(qualified_gnm_folder) is True:
    os.system('rm -r %s' % qualified_gnm_folder)
os.system('mkdir %s' % qualified_gnm_folder)

ffn_file_re = '%s/*.ffn' % ffn_folder
ffn_file_list = [os.path.basename(file_name) for file_name in glob.glob(ffn_file_re)]

qualified_16s_seq_handle = open(qualified_16s_seq, 'w')
qualified_gnm_set = set()
for each_ffn in ffn_file_list:

    gnm_id = '.'.join(each_ffn.split('.')[:-1])

    pwd_gnm = '%s/%s.fna'   % (gnm_folder, gnm_id)
    pwd_ffn = '%s/%s'       % (ffn_folder, each_ffn)

    seq_len_dict = {}
    for each_ctg in SeqIO.parse(pwd_gnm, 'fasta'):
        seq_len_dict[each_ctg.id] = len(str(each_ctg.seq))

    index_16s = 1
    for each_rna in SeqIO.parse(pwd_ffn, 'fasta'):
        rna_id = each_rna.id
        if '16S_rRNA' in rna_id:
            rna_id_split = rna_id.split(':')
            located_ctg = rna_id_split[2]
            rna_len = len(str(each_rna.seq))
            ctg_len = seq_len_dict[located_ctg]
            if (rna_len >= len_cutoff_16s) and (ctg_len >= len_cutoff_ctg):
                qualified_gnm_set.add(gnm_id)

                # move genome to qualified folder
                os.system('cp %s %s/' % (pwd_gnm, qualified_gnm_folder))
                qualified_16s_seq_handle.write('>%s_16S_%s\n' % (gnm_id, index_16s))
                qualified_16s_seq_handle.write('%s\n' % str(each_rna.seq))
                index_16s += 1

qualified_16s_seq_handle.close()

# write out
qualified_gnm_id_file = '%s/ref_with_16S_%sbp_%sbp_%s.txt' % (wd, len_cutoff_16s, len_cutoff_ctg, len(qualified_gnm_set))
qualified_gnm_id_file_handle = open(qualified_gnm_id_file, 'w')
#squalified_gnm_id_file_handle.write('cami_genome_id\ttax_id\n')
for qualified_gnm in sorted([i for i in qualified_gnm_set]):
    #qualified_gnm_id_file_handle.write('%s\t%s\n' % (qualified_gnm, genomes_to_tax_dict[qualified_gnm]))
    qualified_gnm_id_file_handle.write('%s\n' % qualified_gnm)
qualified_gnm_id_file_handle.close()