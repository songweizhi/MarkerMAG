import os
import glob


sample_list = ['Airways', 'Gastrointestinal_tract', 'Oral', 'Skin', 'Urogenital_tract']


for each_sample in sample_list:
    total_abun_file = '/Users/songweizhi/Desktop/ref_total_abun_%s.txt' % each_sample
    abundance_folder = '/Users/songweizhi/Desktop/abundance_%s' % each_sample
    abundance_file_re = '%s/*.tsv' % abundance_folder
    abundance_file_list = glob.glob(abundance_file_re)
    current_sample_ref_set = set()
    ref_total_abun_dict = {}
    for abundance_file in abundance_file_list:
        for each_ref in open(abundance_file):
            each_ref_split = each_ref.strip().split('\t')
            ref_id = each_ref_split[0]
            ref_abun = float(each_ref_split[1])
            if ref_abun > 0:
                current_sample_ref_set.add(ref_id)

                if ref_id not in ref_total_abun_dict:
                    ref_total_abun_dict[ref_id] = ref_abun
                else:
                    ref_total_abun_dict[ref_id] += ref_abun

    total_abun_file_handle = open(total_abun_file, 'w')
    for each_red_abun in ref_total_abun_dict:
        total_abun_file_handle.write('%s\t%s\n' % (each_red_abun, ref_total_abun_dict[each_red_abun]))
    total_abun_file_handle.close()

    print('%s\t%s' % (len(current_sample_ref_set), each_sample))



# wd = '/Users/songweizhi/Desktop/CAMI_oral'
# #wd = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI_Oral/0_raw_files'
# genome_to_id_tsv = '%s/genome_to_id.tsv' % wd
# file_re = '%s/abundance/*.tsv' % wd
# file_list = glob.glob(file_re)
# fake_quality = '%s/oral_799_ref_genomes_fake_quality.txt' % wd
#
#
# genome_id_to_name_dict = {}
# for each_ref in open(genome_to_id_tsv):
#     each_ref_split = each_ref.strip().split('\t')
#     genome_id_to_name_dict[each_ref_split[0]] = each_ref_split[1]
#
#
#
# fake_quality_handle = open(fake_quality, 'w')
# fake_quality_handle.write('genome,completeness,contamination\n')
# ref_id_set = set()
# for each_file in file_list:
#     for each_line in open(each_file):
#         gnm_id = each_line.split('\t')[0]
#         gnm_name = genome_id_to_name_dict.get(gnm_id, 'na')
#         abun = float(each_line.split('\t')[1])
#         if abun > 0:
#             cp_cmd = 'cp 0_raw_files/genomes/%s oral_799_ref_genomes/%s.fa' % (gnm_name, gnm_id)
#             #print(cp_cmd)
#             #os.system(cp_cmd)
#             print(abun)
#             fake_quality_handle.write('%s.fa,%s,%s\n' % (gnm_id, randrange(95, 100), randrange(0, 3)))
#
#             ref_id_set.add(gnm_id)
# print(len(ref_id_set))
#
# fake_quality_handle.close()

