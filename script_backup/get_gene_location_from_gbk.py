#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO

'''
Example command
python3 get_gene_location_from_gbk.py -h
python3 get_gene_location_from_gbk.py -id gene_id.txt -gbk Test_pcofg_prodigal_output -o gene_location.txt
'''

parser = argparse.ArgumentParser()
parser.add_argument('-id',  required=True, help='gene id file, one id per line')
parser.add_argument('-gbk', required=True, help='folder holds gbk files, e.g. prodigal_output folder from MetaCHIP')
parser.add_argument('-o',   required=True, help='output file')
args = vars(parser.parse_args())

gene_id_file        = args['id']
gbk_folder          = args['gbk']
gene_location_file  = args['o']


genome_to_gene_dict = dict()
for each_gene in open(gene_id_file):
    gene_id = each_gene.strip()
    genome_id = '_'.join(gene_id.split('_')[:-1])
    if genome_id not in genome_to_gene_dict:
        genome_to_gene_dict[genome_id] = {gene_id}
    else:
        genome_to_gene_dict[genome_id].add(gene_id)


gene_location_file_handle = open(gene_location_file, 'w')
gene_location_file_handle.write('Gene\tGenome\tContig\n' )
for each_gnm in genome_to_gene_dict:

    gnm_gbk = '%s/%s.gbk' % (gbk_folder, each_gnm)
    current_gnm_genes = genome_to_gene_dict[each_gnm]

    if os.path.isfile(gnm_gbk) is False:
        for each_current_gene in current_gnm_genes:
            gene_location_file_handle.write('%s\t%s\tgbk_not_found\n' % (each_current_gene, each_gnm))
    else:
        gene_to_ctg_dict = dict()
        for ctg_record in SeqIO.parse(gnm_gbk, "genbank"):
            ctg_record_id = ctg_record.id
            for gene_record in ctg_record.features:
                if 'locus_tag' in gene_record.qualifiers:
                    gene_record_id = gene_record.qualifiers['locus_tag'][0]
                    if gene_record_id in current_gnm_genes:
                        gene_to_ctg_dict[gene_record_id] = ctg_record_id
        # write out
        for each_gene_id in current_gnm_genes:
            gene_location_file_handle.write('%s\t%s\t%s\n' % (each_gene_id, each_gnm, gene_to_ctg_dict.get(each_gene_id, 'gene_not_found_in_gbk')))

gene_location_file_handle.close()
