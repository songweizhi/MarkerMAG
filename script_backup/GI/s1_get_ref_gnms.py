import os
import glob
from random import randrange


wd = '/Users/songweizhi/Desktop/CAMI_oral'
#wd = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/0_raw_files'
genome_to_id_tsv = '%s/genome_to_id.tsv' % wd
file_re = '%s/abundance/*.tsv' % wd
file_list = glob.glob(file_re)
fake_quality = '%s/oral_799_ref_genomes_fake_quality.txt' % wd


genome_id_to_name_dict = {}
for each_ref in open(genome_to_id_tsv):
    each_ref_split = each_ref.strip().split('\t')
    genome_id_to_name_dict[each_ref_split[0]] = each_ref_split[1]



fake_quality_handle = open(fake_quality, 'w')
fake_quality_handle.write('genome,completeness,contamination\n')
ref_id_set = set()
for each_file in file_list:
    for each_line in open(each_file):
        gnm_id = each_line.split('\t')[0]
        gnm_name = genome_id_to_name_dict.get(gnm_id, 'na')
        abun = float(each_line.split('\t')[1])
        if abun > 0:
            cp_cmd = 'cp 0_raw_files/genomes/%s oral_799_ref_genomes/%s.fa' % (gnm_name, gnm_id)
            #print(cp_cmd)
            #os.system(cp_cmd)
            print(abun)
            fake_quality_handle.write('%s.fa,%s,%s\n' % (gnm_id, randrange(95, 100), randrange(0, 3)))

            ref_id_set.add(gnm_id)
print(len(ref_id_set))

fake_quality_handle.close()

