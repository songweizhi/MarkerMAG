import os


abundance_df_file       = '/Users/songweizhi/Desktop/gt_samples.txt'
genome_to_id_file       = '/Users/songweizhi/Desktop/genome_to_id.tsv'

abundance_df_file       = 'gt_samples.txt'
genome_to_id_file       = 'raw_files/genome_to_id.tsv'

genome_folder           = 'ref_genomes'
genome_folder_renamed   = 'ref_genomes_renamed'

genome_folder_subset    = 'ref_genomes_GI'

#################################################### rename genomes ####################################################

# # read in genome id into dict
# id_to_name_dict = {}
# for each_gnm in open(genome_to_id_file):
#     each_gnm_split = each_gnm.strip().split()
#     id_to_name_dict[each_gnm_split[0]] = each_gnm_split[1]
#
# # rename genome
# for each_line in open(abundance_df_file):
#     each_line_split = each_line.strip().split('\t')
#     genome_id = each_line_split[0]
#     genome_name = id_to_name_dict[genome_id]
#
#     cp_cmd = 'cp %s/%s %s/%s.fa' % (genome_folder, genome_name, genome_folder_renamed, genome_id)
#     os.system(cp_cmd)


#################################################### subset genomes ####################################################

n = 0
for each_line in open(abundance_df_file):
    each_line_split = each_line.strip().split('\t')
    genome_id = each_line_split[0]
    if each_line_split[1:] != ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0']:
        cp_cmd = 'cp %s/%s.fa %s/' % (genome_folder_renamed, genome_id, genome_folder_subset)
        print(cp_cmd)
        os.system(cp_cmd)
        n += 1
print(n)

