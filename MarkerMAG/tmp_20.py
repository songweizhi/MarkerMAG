from Bio import SeqIO

gtdb_ssu_fna = '/Users/songweizhi/Desktop/GTDB_ssu_all_r95.fna'
gtdb_gnm_to_ssu = '/Users/songweizhi/Desktop/GTDB_ssu_all_r95_gnm_to_ssu.txt'

gnm_to_ssu_dict= {}
for each_ssu in SeqIO.parse(gtdb_ssu_fna, 'fasta'):
    ssu_id = each_ssu.id
    assembly_id = ssu_id.split('~')[0][3:]
    if assembly_id not in gnm_to_ssu_dict:
        gnm_to_ssu_dict[assembly_id] = {ssu_id}
    else:
        gnm_to_ssu_dict[assembly_id].add(ssu_id)


print(gnm_to_ssu_dict)
gtdb_gnm_to_ssu_handle = open(gtdb_gnm_to_ssu, 'w')
for each_gnm in gnm_to_ssu_dict:
    gtdb_gnm_to_ssu_handle.write('%s\t%s\n' % (each_gnm, ','.join(gnm_to_ssu_dict[each_gnm])))



gtdb_gnm_to_ssu_handle.close()