import os
import argparse
from Bio import SeqIO


def sep_path_basename_ext(file_in):
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'
    file_basename, file_extension = os.path.splitext(file_name)
    return file_path, file_basename, file_extension


parser = argparse.ArgumentParser()
parser.add_argument('-ref',   required=True, help='reference')
parser.add_argument('-reads', required=True, help='reads')
args = vars(parser.parse_args())
ref_seq = args['ref']
combined_reads_in = args['reads']


ref_path, ref_basename, ref_extension = sep_path_basename_ext(ref_seq)
reads_in_path, reads_in_basename, reads_in_extension = sep_path_basename_ext(combined_reads_in)
op_r1       = '%s/%s_R1.fa' % (reads_in_path, reads_in_basename)
op_r2       = '%s/%s_R2.fa' % (reads_in_path, reads_in_basename)
op_up       = '%s/%s_UP.fa' % (reads_in_path, reads_in_basename)
ref_index   = '%s/%s'       % (ref_path, ref_basename)
op_sam      = '%s/%s.sam'   % (ref_path, ref_basename)


read_id_set = set()
read_base_set = set()
read_seq_dict = {}
for each_read in SeqIO.parse(combined_reads_in, 'fasta'):
    read_id = each_read.id
    read_base = '.'.join(read_id.split('.')[:-1])
    read_id_set.add(read_id)
    read_base_set.add(read_base)
    read_seq_dict[read_id] = str(each_read.seq)

paired_base_set = set()
unpaired_read_set = set()
for each_base in read_base_set:
    base_r1 = '%s.1' % each_base
    base_r2 = '%s.2' % each_base
    if (base_r1 in read_id_set) and (base_r2 in read_id_set):
        paired_base_set.add(each_base)
    elif (base_r1 in read_id_set) and (base_r2 not in read_id_set):
        unpaired_read_set.add(base_r1)
    elif (base_r1 not in read_id_set) and (base_r2 in read_id_set):
        unpaired_read_set.add(base_r2)

if len(paired_base_set) > 0:
    op_r1_handle = open(op_r1, 'w')
    op_r2_handle = open(op_r2, 'w')
    for each_paired_base in paired_base_set:
        paired_base_r1 = '%s.1' % each_paired_base
        paired_base_r2 = '%s.2' % each_paired_base
        paired_base_r1_seq = read_seq_dict[paired_base_r1]
        paired_base_r2_seq = read_seq_dict[paired_base_r2]
        op_r1_handle.write('>%s\n' % paired_base_r1)
        op_r1_handle.write('%s\n'  % paired_base_r1_seq)
        op_r2_handle.write('>%s\n' % paired_base_r2)
        op_r2_handle.write('%s\n'  % paired_base_r2_seq)
    op_r1_handle.close()
    op_r2_handle.close()

if len(unpaired_read_set) > 0:
    op_up_handle = open(op_up, 'w')
    for each_unpaired_read in unpaired_read_set:
        op_up_handle.write('>%s\n' % each_unpaired_read)
        op_up_handle.write('%s\n'  % read_seq_dict[each_unpaired_read])
    op_up_handle.close()

index_ref_cmd = 'bowtie2-build --quiet -f %s %s' % (ref_seq, ref_index)
bowtie2_cmd = ''
if (len(paired_base_set) > 0) and (len(unpaired_read_set) > 0):
    bowtie2_cmd   = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p 6 -f --local --all --no-unal' % (ref_index, op_r1, op_r2, op_up, op_sam)
if (len(paired_base_set) > 0) and (len(unpaired_read_set) == 0):
    bowtie2_cmd   = 'bowtie2 -x %s -1 %s -2 %s -S %s -p 6 -f --local --all --no-unal' % (ref_index, op_r1, op_r2, op_sam)
if (len(paired_base_set) == 0) and (len(unpaired_read_set) > 0):
    bowtie2_cmd   = 'bowtie2 -x %s -U %s -S %s -p 6 -f --local --all --no-unal' % (ref_index, op_up, op_sam)
if bowtie2_cmd != '':
    os.system('module load bowtie')
    os.system(index_ref_cmd)
    os.system(bowtie2_cmd)
else:
    print('No reads found!')
