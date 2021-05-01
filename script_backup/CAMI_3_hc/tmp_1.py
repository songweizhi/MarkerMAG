import os

bin_id_file = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/cami_hc_refined_bins_id.txt'
bin_id_file = '/srv/scratch/z5039045/MarkerMAG_wd/CAMI1/ref_id.txt'
#bin_id_file = '/Users/songweizhi/Desktop/cami_hc_refined_bins_id.txt'

for each in open(bin_id_file):
    file_name = each.strip()
    bin_id = file_name.split('.')[0]
    #add_prefix_cmd = 'BioSAK rename_seq -in /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/3_High_Complexity/cami_hc_refined_bins/%s -prefix %s' % (file_name, bin_id)
    add_prefix_cmd = 'BioSAK rename_seq -in /srv/scratch/z5039045/MarkerMAG_wd/CAMI1/ref_genomes_hc_04_20/%s -prefix %s' % (file_name, bin_id)
    print(add_prefix_cmd)
    os.system(add_prefix_cmd)


