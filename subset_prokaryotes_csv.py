
def subset_prokaryotes_csv(file_in, id_file_in, prokaryotes_csv, file_out):

    num_to_id_dict = {}
    GTDB_genome_id_list = set()

    if id_file_in is True:
        for each_id in open(file_in):

            NCBI_ID_without_version = each_id.strip()
            if '.' in each_id.strip():
                NCBI_ID_without_version = '.'.join(each_id.strip().split('.')[:-1])

            NCBI_ID_no_version_only_number = NCBI_ID_without_version[4:]

            GTDB_genome_id_list.add(NCBI_ID_no_version_only_number)
            num_to_id_dict[NCBI_ID_no_version_only_number] = each_id.strip()

    else:
        for each in open(file_in):
            if not each.startswith('"ID",'):
                each_split = each.strip().split(',')

                NCBI_ID_with_version = each_split[0][1:-1]
                NCBI_ID_without_version = '.'.join(NCBI_ID_with_version.split('.')[:-1])
                NCBI_ID_no_version_only_number = NCBI_ID_without_version[4:]

                GTDB_genome_id_list.add(NCBI_ID_no_version_only_number)
                num_to_id_dict[NCBI_ID_no_version_only_number] = NCBI_ID_with_version


    GenBank_genomes_subset_handle = open(file_out, 'w')
    found_in_GenBank_genome_id_list = set()
    for GenBank_genome in open(prokaryotes_csv):
        if not GenBank_genome.startswith('#Organism Name'):
            GenBank_genome_split = GenBank_genome.strip().split(',')
            GenBank_genome_id_with_version = GenBank_genome_split[5][1:-1]
            GenBank_genome_id_without_version = '.'.join(GenBank_genome_id_with_version.split('.')[:-1])
            GenBank_genome_id_no_version_only_number = GenBank_genome_id_without_version[4:]
            if GenBank_genome_id_no_version_only_number in GTDB_genome_id_list:
                found_in_GenBank_genome_id_list.add(GenBank_genome_id_no_version_only_number)
                GenBank_genomes_subset_handle.write(GenBank_genome)
    GenBank_genomes_subset_handle.close()

    for each in GTDB_genome_id_list:
        if each not in found_in_GenBank_genome_id_list:
            print('not found: %s' % num_to_id_dict[each])


# file_in         = '/Users/songweizhi/Desktop/select_strain/GTDB-Pseudomonas_E.csv'
file_in         = '/Users/songweizhi/Desktop/accession_id.txt'
id_file_in      = True
prokaryotes_csv = '/Users/songweizhi/Desktop/select_strain/prokaryotes.csv'
file_out        = '/Users/songweizhi/Desktop/GenBank_g_test.csv'

subset_prokaryotes_csv(file_in, id_file_in, prokaryotes_csv, file_out)

