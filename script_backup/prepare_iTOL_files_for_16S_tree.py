import os

wd               = '/Users/songweizhi/Desktop/MarkerMAG_wd/Figures/16S_tree'
gene_list_file   = '%s/gene_list.txt'                                           % wd
rank_list        = ['d', 'p', 'c', 'o', 'f', 'g', 's']
ranks_with_color = ['d', 'p', 'c', 'o', 'f', 'g']

gene_list = []
for gene in open(gene_list_file):
    gene_list.append(gene.strip())


for taxon_rank in rank_list:
    print(taxon_rank)
    genome_to_taxon     = '%s/16S_grouping_rank_%s.txt'                 % (wd, taxon_rank)
    gene_to_taxon       = '%s/16S_grouping_rank_%s_gene_level.txt'      % (wd, taxon_rank)
    gene_to_taxon_itol  = '%s/16S_grouping_rank_%s_gene_level_iTOL.txt' % (wd, taxon_rank)
    group_to_color      = '%s/16S_group_color_%s.txt'                   % (wd, taxon_rank)

    genome_to_taxon_dict = {}
    for genome_taxon in open(genome_to_taxon):
        genome_taxon_split = genome_taxon.strip().split('\t')
        genome_to_taxon_dict[genome_taxon_split[0]] = '_'.join(genome_taxon_split[1].split('_'))

    gene_to_taxon_handle = open(gene_to_taxon, 'w')
    for gene in sorted(gene_list):
        gene_to_taxon_handle.write('%s\t%s\n' % (gene, genome_to_taxon_dict[gene[:2]]))
    gene_to_taxon_handle.close()

    BioSAK_cmd = 'BioSAK iTOL -ColorStrip -lg %s -lt %s -out %s' % (gene_to_taxon, taxon_rank, gene_to_taxon_itol)
    if taxon_rank in ranks_with_color:
        BioSAK_cmd = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt %s -out %s' % (gene_to_taxon, group_to_color, taxon_rank, gene_to_taxon_itol)

    if taxon_rank == 's':
        #BioSAK_cmd = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt %s -out %s' % (gene_to_taxon, group_to_color, taxon_rank, gene_to_taxon_itol)
        BioSAK_cmd = 'BioSAK iTOL -ColorRange -lg %s -gc %s -lt %s -out %s' % (gene_to_taxon, group_to_color, taxon_rank, gene_to_taxon_itol)

    os.system(BioSAK_cmd)

    os.remove(gene_to_taxon)
