import os
from Bio import SeqIO


###################################################### file in/out #####################################################

wd = '/Users/songweizhi/Desktop/ttt_mag'

# file in
drep_cdb_file               = '%s/Cdb.csv'                         % wd
fastANI_output              = '%s/GI_fastANI.txt'                  % wd
ani_cutoff                  = 99
mag_folder                  = '%s/GI_refined_bins'                  % wd
mag_ext                     = 'fasta'

# file out
rename_mag_file             = '%s/mag_rename.txt'                  % wd
mag_folder_renamed          = '%s/GI_refined_bins_renamed'         % wd


########################################################################################################################

cluster_to_ref_dict = {}
ref_to_cluster_dict = {}
for each_ref in open(drep_cdb_file):
    if not each_ref.startswith('genome,secondary_cluster'):
        each_ref_split = each_ref.strip().split(',')
        ref_file_name = each_ref_split[0]
        ref_file_name_no_ext = '.'.join(ref_file_name.split('.')[:-1])
        ref_cluster = 'C' + each_ref_split[1]
        ref_to_cluster_dict[ref_file_name_no_ext] = ref_cluster
        if ref_cluster not in cluster_to_ref_dict:
            cluster_to_ref_dict[ref_cluster] = [ref_file_name_no_ext]
        else:
            cluster_to_ref_dict[ref_cluster].append(ref_file_name_no_ext)

mag_id_set = set()
ref_to_mag_dict = {}
mag_to_ref_dict = {}
cluster_to_mag_dict = {}
mag_to_cluster_dict = {}
for line in open(fastANI_output):
    line_split = line.strip().split('\t')
    query_name = line_split[0]
    if '/' in query_name:
        query_name = query_name.split('/')[-1]
    ref_name = line_split[1]
    if '/' in ref_name:
        ref_name = ref_name.split('/')[-1]

    query_name = '.'.join(query_name.split('.')[:-1])
    ref_name = '.'.join(ref_name.split('.')[:-1])
    ref_cluster = ref_to_cluster_dict[ref_name]

    mag_id_set.add(query_name)

    ani_value = float(line_split[2])
    if ani_value >= ani_cutoff:
        if ref_name not in ref_to_mag_dict:
            ref_to_mag_dict[ref_name] = {query_name}
        else:
            ref_to_mag_dict[ref_name].add(query_name)
        if query_name not in mag_to_ref_dict:
            mag_to_ref_dict[query_name] = {ref_name}
        else:
            mag_to_ref_dict[query_name].add(ref_name)

        if ref_cluster not in cluster_to_mag_dict:
            cluster_to_mag_dict[ref_cluster] = {query_name}
        else:
            cluster_to_mag_dict[ref_cluster].add(query_name)
        if query_name not in mag_to_cluster_dict:
            mag_to_cluster_dict[query_name] = {ref_cluster}
        else:
            mag_to_cluster_dict[query_name].add(ref_cluster)


mags_assigned_to_multi_ref_genome = set()
for mag in mag_to_ref_dict:
    if len(mag_to_ref_dict[mag]) > 1:
        mags_assigned_to_multi_ref_genome.add(mag)

mags_assigned_to_multi_clusters = set()
for mag in mag_to_cluster_dict:
    if len(mag_to_cluster_dict[mag]) > 1:
        mags_assigned_to_multi_clusters.add(mag)


if len(mags_assigned_to_multi_ref_genome) == 0:
    print('Good! no MAG assigned to multiple reference genomes')
else:
    print('Caution! found MAG assigned to multiple reference genomes: %s' % ','.join(mags_assigned_to_multi_ref_genome))

if len(mags_assigned_to_multi_clusters) == 0:
    print('Good! no MAG assigned to multiple clusters')
else:
    print('Caution! found MAG assigned to multiple clusters: %s' % ','.join(mags_assigned_to_multi_clusters))


if len(mags_assigned_to_multi_clusters) == 0:

    # get rename dict
    rename_mag_file_handle = open(rename_mag_file, 'w')
    rename_dict_mag = {}
    for each_cluster in cluster_to_mag_dict:
        current_cluster_mags = [i for i in cluster_to_mag_dict[each_cluster]]
        if len(current_cluster_mags) == 1:
            current_mag_new_id = '%s_bin' % each_cluster
            rename_dict_mag[current_cluster_mags[0]] = current_mag_new_id
            rename_mag_file_handle.write('%s\t%s\n' % (current_cluster_mags[0], current_mag_new_id))
        else:
            index_mag = 1
            for current_mag in current_cluster_mags:
                current_mag_new_id = '%s_bin%s' % (each_cluster, index_mag)
                rename_dict_mag[current_mag] = current_mag_new_id
                rename_mag_file_handle.write('%s\t%s\n' % (current_mag, current_mag_new_id))
                index_mag += 1
    rename_mag_file_handle.close()

    # rename mags
    for mag in rename_dict_mag:
        rename_cmd = 'cp %s/%s.%s %s/%s.fasta' % (mag_folder, mag, mag_ext, mag_folder_renamed, rename_dict_mag[mag])
        os.system(rename_cmd)

    mag_without_assignment = set()
    for mag in mag_id_set:
        if mag not in mag_to_cluster_dict:
            mag_without_assignment.add(mag)
    if len(mag_without_assignment) > 0:
        print('The following %s MAGs not assigned to any cluster: %s' % (len(mag_without_assignment), ','.join(mag_without_assignment)))
