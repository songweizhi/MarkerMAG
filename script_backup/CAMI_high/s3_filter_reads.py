import argparse
from Bio import SeqIO

# BINID in gs_read_mapping_1.binning refers to cami_genome_id

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, type=int, help='sample_index')
args = vars(parser.parse_args())
sample_index = args['i']

wd = '.'
qualified_gnm_id_file   = '%s/ref_with_16S_1000bp_5000bp_416.txt'   % wd
gs_read_mapping_file    = '%s/gs_read_mapping_%s.binning'           % (wd, sample_index)
r1_in                   = '%s/RH_S00%s_R1.fa'                       % (wd, sample_index)
r2_in                   = '%s/RH_S00%s_R2.fa'                       % (wd, sample_index)
r1_out                  = '%s/RH_S00%s_R1_filtered.fa'              % (wd, sample_index)
r2_out                  = '%s/RH_S00%s_R2_filtered.fa'              % (wd, sample_index)

qualified_gnm_id_set = set()
for each_ref in open(qualified_gnm_id_file):
    qualified_gnm_id_set.add(each_ref.strip())

read_to_source_dict = {}
for each_read in open(gs_read_mapping_file):
    if (not each_read.startswith('@')) and (len(each_read) > 1):
        each_read_split = each_read.strip().split('\t')
        read_base = each_read_split[0]
        ref_id = each_read_split[1]
        read_to_source_dict[read_base] = ref_id

r1_out_handle = open(r1_out, 'w')
for each_read in SeqIO.parse(r1_in, 'fasta'):
    read_id = each_read.id
    read_id_base = read_id.split('/')[0]
    read_source = read_to_source_dict.get(read_id_base, 'na')
    if read_source in qualified_gnm_id_set:
        r1_out_handle.write('>%s\n' % read_id)
        r1_out_handle.write('%s\n' % str(each_read.seq))
r1_out_handle.close()

r2_out_handle = open(r2_out, 'w')
for each_read in SeqIO.parse(r2_in, 'fasta'):
    read_id = each_read.id
    read_id_base = read_id.split('/')[0]
    read_source = read_to_source_dict.get(read_id_base, 'na')
    if read_source in qualified_gnm_id_set:
        r2_out_handle.write('>%s\n' % read_id)
        r2_out_handle.write('%s\n' % str(each_read.seq))
r2_out_handle.close()
