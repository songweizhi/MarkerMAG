
gtdb_metadata      = '/Users/songweizhi/Desktop/select_strain/bac120_metadata_r95.tsv'
get_genome_from    = ['d__Bacteria']
taxon_to_exculde   = ['p__Proteobacteria']
complete_genome    = True
ncbi_type_material = True
gtdb_type_strain   = True


x = 0
qualified_genome_num = 0
for genome_record in open(gtdb_metadata):
    genome_record_split             = genome_record.strip().split('\t')
    accession                       = genome_record_split[0]
    gtdb_taxonomy                   = genome_record_split[16]
    ncbi_assembly_level             = genome_record_split[45]
    ncbi_type_material_designation  = genome_record_split[84]
    ssu_count                       = genome_record_split[91]

    if ('d__Bacteria' in gtdb_taxonomy) and (ncbi_assembly_level not in ['Contig', 'Scaffold']) and ('p__Proteobacteria' not in gtdb_taxonomy) and (ncbi_type_material_designation == 'assembly from type material') and (genome_record_split[17] == 'type strain of species'):
        print('%s\t%s' % (accession, gtdb_taxonomy))
        qualified_genome_num += 1

print('Number of qualified genome is %s' % qualified_genome_num)





