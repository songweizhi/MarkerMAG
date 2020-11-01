import os
import glob


mag_folder    = '/Users/songweizhi/Desktop/Refined_refined_bins_renamed'
mag_ext       = 'fna'
combined_16s  = '/Users/songweizhi/Desktop/combined_16S.ffn'
combined_mags = '/Users/songweizhi/Desktop/combined_MAGs.fna'


# rename ctgs in mags
mag_file_re = '%s/*.%s' % (mag_folder, mag_ext)
mag_file_list = [os.path.basename(file_name) for file_name in glob.glob(mag_file_re)]
for mag_file in mag_file_list:
    mag_file_no_ext = '.'.join(mag_file.split('.')[:-1])
    rename_cmd = 'BioSAK rename_seq -in %s/%s -prefix %s' % (mag_folder, mag_file, mag_file_no_ext)
    os.system(rename_cmd)


# combine renamed mags
combine_cmd = 'cat %s/*_renamed.fna > %s' % (mag_folder, combined_mags)
remove_renamed_mags_cmd = 'rm %s/*_renamed.fna' % mag_folder
os.system(combine_cmd)
os.system(remove_renamed_mags_cmd)




