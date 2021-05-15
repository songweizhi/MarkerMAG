import os

wd = '/Users/songweizhi/Desktop/CAMI_high'
wd = '.'
ref_596_genomes   = '%s/goldstandard_high_pool.filtered.profile.stain.596_ref.txt' % wd
source_genome_txt = '%s/source_genomes.txt' % wd


ref_596_genomes_set = set()
for each_ref in open(ref_596_genomes):
    each_ref_split = each_ref.strip().split('\t')
    cami_genome_id = each_ref_split[5]
    ref_596_genomes_set.add(cami_genome_id)

for each_gnm_file in open(source_genome_txt):
    file_name = each_gnm_file.strip()
    gnm_id = ''
    if '_run' in file_name:
        gnm_id = file_name.split('_run')[0]
        #print(gnm_id)
    elif '.gt1kb.' in file_name:
        gnm_id = file_name.split('.gt1kb.')[0]
    else:
        gnm_id = '.'.join(file_name.split('.')[:-1])
    if gnm_id in ref_596_genomes_set:
        rename_cmd = 'cp /srv/scratch/z5039045/MarkerMAG_wd/CAMI_high/source_genomes/%s /srv/scratch/z5039045/MarkerMAG_wd/CAMI_high/ref_gnm_596/%s.fna' % (file_name, gnm_id)
        os.system(rename_cmd)
