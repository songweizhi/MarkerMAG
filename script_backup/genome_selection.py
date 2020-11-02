

def subset_prokaryotes_csv(GTDB_genomes, GenBank_genomes, GenBank_genomes_subset):

    num_to_full_id_dict = {}
    GTDB_genome_id_list = set()
    for each in open(GTDB_genomes):
        if not each.startswith('"ID",'):
            each_split = each.strip().split(',')
            NCBI_ID_with_version = each_split[0][1:-1]
            NCBI_ID_without_version = '.'.join(NCBI_ID_with_version.split('.')[:-1])
            NCBI_ID_no_version_only_number = NCBI_ID_without_version[3:]
            GTDB_genome_id_list.add(NCBI_ID_no_version_only_number)
            num_to_full_id_dict[NCBI_ID_no_version_only_number] = NCBI_ID_with_version

    GenBank_genomes_subset_handle = open(GenBank_genomes_subset, 'w')
    found_in_GenBank_genome_id_list = set()
    for GenBank_genome in open(GenBank_genomes):
        if not GenBank_genome.startswith('#Organism Name'):
            GenBank_genome_split = GenBank_genome.strip().split(',')
            GenBank_genome_id_with_version = GenBank_genome_split[5][1:-1]
            GenBank_genome_id_without_version = '.'.join(GenBank_genome_id_with_version.split('.')[:-1])
            GenBank_genome_id_no_version_only_number = GenBank_genome_id_without_version[3:]
            if GenBank_genome_id_no_version_only_number in GTDB_genome_id_list:
                found_in_GenBank_genome_id_list.add(GenBank_genome_id_no_version_only_number)
                GenBank_genomes_subset_handle.write(GenBank_genome)
    GenBank_genomes_subset_handle.close()

    for each in GTDB_genome_id_list:
        if each not in found_in_GenBank_genome_id_list:
            print('not found: %s' % num_to_full_id_dict[each])


def subset_prokaryotes_csv_by_accession(GenBank_genomes, accession_file, GenBank_genomes_subset):

    num_to_full_id_dict = {}
    genomes_in_accession = set()
    for each in open(accession_file):
        if not each.startswith('"ID",'):
            each_split = each.strip().split(',')
            NCBI_ID_with_version = each_split[0][1:-1]
            NCBI_ID_without_version = '.'.join(NCBI_ID_with_version.split('.')[:-1])
            NCBI_ID_no_version_only_number = NCBI_ID_without_version[3:]
            genomes_in_accession.add(NCBI_ID_no_version_only_number)
            num_to_full_id_dict[NCBI_ID_no_version_only_number] = NCBI_ID_with_version

    GenBank_genomes_subset_handle = open(GenBank_genomes_subset, 'w')
    found_in_GenBank_genome_id_list = set()
    for GenBank_genome in open(GenBank_genomes):
        if not GenBank_genome.startswith('#Organism Name'):
            GenBank_genome_split = GenBank_genome.strip().split(',')
            GenBank_genome_id_with_version = GenBank_genome_split[5][1:-1]
            GenBank_genome_id_without_version = '.'.join(GenBank_genome_id_with_version.split('.')[:-1])
            GenBank_genome_id_no_version_only_number = GenBank_genome_id_without_version[3:]
            if GenBank_genome_id_no_version_only_number in genomes_in_accession:
                found_in_GenBank_genome_id_list.add(GenBank_genome_id_no_version_only_number)
                GenBank_genomes_subset_handle.write(GenBank_genome)
    GenBank_genomes_subset_handle.close()

    for each in genomes_in_accession:
        if each not in found_in_GenBank_genome_id_list:
            print('not found: %s' % num_to_full_id_dict[each])

GenBank_genomes                                 = '/Users/songweizhi/Desktop/select_strain/prokaryotes.csv'


GTDB_Pseudomonas_csv                            = '/Users/songweizhi/Desktop/select_strain/GTDB-Pseudomonas.csv'
GTDB_o__Pseudomonadales_csv                     = '/Users/songweizhi/Desktop/select_strain/GTDB-o__Pseudomonadales.csv'
GTDB_o__Pseudomonadales_no_Pseudomonadaceae_csv = '/Users/songweizhi/Desktop/select_strain/GTDB-o__Pseudomonadales_no_Pseudomonadaceae.csv'
GTDB_f__Pseudomonadaceae_no_Pseudomonas_E_csv   = '/Users/songweizhi/Desktop/select_strain/GTDB-f__Pseudomonadaceae_no_Pseudomonas_E.csv'
GTDB_f__Pseudomonadaceae_csv                    = '/Users/songweizhi/Desktop/select_strain/GTDB-f__Pseudomonadaceae.csv'
GTDB_g__Pseudomonas_E_csv                       = '/Users/songweizhi/Desktop/select_strain/GTDB-Pseudomonas_E.csv'
GenBank_genomes_subset_g__Pseudomonas_E         = '/Users/songweizhi/Desktop/select_strain/GenBank_g__Pseudomonas_E_subset.csv'
GenBank_genomes_subset_f__Pseudomonadaceae_no_Pseudomonas_E   = '/Users/songweizhi/Desktop/select_strain/GenBank_f__Pseudomonadaceae_subset_no_Pseudomonas_E.csv'
GenBank_genomes_subset_o__Pseudomonadales_no_Pseudomonadaceae = '/Users/songweizhi/Desktop/select_strain/GenBank_o__Pseudomonadales_no_Pseudomonadaceae.csv'


GTDB_Pseudomonas_E_csv_handle = open(GTDB_g__Pseudomonas_E_csv, 'w')
for genome in open(GTDB_Pseudomonas_csv):
    genome_split = genome.strip().split(',')
    GTDB_Taxonomy = genome_split[3]
    if 'g__Pseudomonas_E' in GTDB_Taxonomy:
        GTDB_Pseudomonas_E_csv_handle.write(genome)
GTDB_Pseudomonas_E_csv_handle.close()


subset_prokaryotes_csv(GTDB_g__Pseudomonas_E_csv, GenBank_genomes, GenBank_genomes_subset_g__Pseudomonas_E)
subset_prokaryotes_csv(GTDB_f__Pseudomonadaceae_no_Pseudomonas_E_csv, GenBank_genomes, GenBank_genomes_subset_f__Pseudomonadaceae_no_Pseudomonas_E)
subset_prokaryotes_csv(GTDB_o__Pseudomonadales_no_Pseudomonadaceae_csv, GenBank_genomes, GenBank_genomes_subset_o__Pseudomonadales_no_Pseudomonadaceae)


# for line in open(GenBank_genomes_subset_g__Pseudomonas_E):
#     line_split = line.strip().split(',')
#     if line_split[6] in ['"Complete"', '" Chromosome"']:
#         pass
#         #print(line.strip())
#         #print()


# for line_f in open(GTDB_o__Pseudomonadales_csv):
#     line_f_split = line_f.strip().split(',')
#     GTDB_Taxonomy = line_f_split[3]
#
#     if 'f__Pseudomonadaceae' not in GTDB_Taxonomy:
#         print(line_f.strip())
#



# for line in open(GenBank_genomes_subset_o__Pseudomonadales_no_Pseudomonadaceae):
#     line_split = line.strip().split(',')
#     if line_split[6] in ['"Complete"', '" Chromosome"']:
#
#         print(line.strip())
#         #print()




'''

cd /Users/songweizhi/Desktop/select_strain
BioSAK dwnld_GenBank_genome -csv GenBank_Pseudomonas_E_subset.csv -fna -name -t 4
BioSAK dwnld_GenBank_genome -csv GenBank_Pseudomonas_E_subset_completed.csv -fna -name -t 4

BioSAK dwnld_GenBank_genome -csv GenBank_f__Pseudomonadaceae_subset_no_Pseudomonas_E_completed.csv -fna -name -t 4


'''

x= 0
for line3 in open('/Users/songweizhi/Desktop/select_strain/bac120_metadata_r95.tsv'):
    line3_split = line3.strip().split('\t')
    accession = line3_split[0]
    gtdb_taxonomy = line3_split[16]
    ncbi_assembly_level = line3_split[45]
    ncbi_type_material_designation = line3_split[84]
    ssu_count = line3_split[91]

    if ('d__Bacteria' in gtdb_taxonomy) and (ncbi_assembly_level not in ['Contig', 'Scaffold']) and ('p__Proteobacteria' not in gtdb_taxonomy) and (ncbi_type_material_designation == 'assembly from type material') and (line3_split[17] == 'type strain of species'):
            #and ('g__Pseudomonas_E' not in gtdb_taxonomy) and (ncbi_assembly_level not in ['Contig', 'Scaffold']) and ('g__Pseudomonas_D' in gtdb_taxonomy):
            #
            #


        print('%s\t%s' % (accession, gtdb_taxonomy))
        #print(line3_split[17])
        #print(ncbi_type_material_designation)
        x += 1
print(x)