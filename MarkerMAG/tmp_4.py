import os
import argparse
from Bio import SeqIO



def mapping_worker(argument_list):

    vis_folder          = argument_list[0]
    each_marker_to_ctg  = argument_list[1]
    concatenated        = argument_list[2]

    pwd_seq_file_cbd        = '%s/%s/%s_cbd.fasta'      % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_16s        = '%s/%s/%s_16s.fasta'      % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_ctg        = '%s/%s/%s_ctg.fasta'      % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_reads      = '%s/%s/%s_reads.fasta'    % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_reads_r1   = '%s/%s/%s_reads_R1.fasta' % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_reads_r2   = '%s/%s/%s_reads_R2.fasta' % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_reads_up   = '%s/%s/%s_reads_UP.fasta' % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_cbd_index  = '%s/%s/%s_cbd'            % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_16s_index  = '%s/%s/%s_16s'            % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_ctg_index  = '%s/%s/%s_ctg'            % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_cbd_sam    = '%s/%s/%s_cbd.sam'        % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_16s_sam    = '%s/%s/%s_16s.sam'        % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)
    pwd_seq_file_ctg_sam    = '%s/%s/%s_ctg.sam'        % (vis_folder, each_marker_to_ctg, each_marker_to_ctg)

    # sep R1, R2 and unpaired reads
    read_id_set = set()
    read_base_set = set()
    read_seq_dict = {}
    for each_read in SeqIO.parse(pwd_seq_file_reads, 'fasta'):
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
        op_r1_handle = open(pwd_seq_file_reads_r1, 'w')
        op_r2_handle = open(pwd_seq_file_reads_r2, 'w')
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
        op_up_handle = open(pwd_seq_file_reads_up, 'w')
        for each_unpaired_read in unpaired_read_set:
            op_up_handle.write('>%s\n' % each_unpaired_read)
            op_up_handle.write('%s\n'  % read_seq_dict[each_unpaired_read])
        op_up_handle.close()


    if concatenated is True:
        index_ref_cmd = 'bowtie2-build --quiet -f %s %s' % (pwd_seq_file_cbd, pwd_seq_file_cbd_index)
        bowtie2_cmd = ''
        if (len(paired_base_set) > 0) and (len(unpaired_read_set) > 0):
            bowtie2_cmd   = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p 1 -f --local --all --no-unal' % (pwd_seq_file_cbd_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_reads_up, pwd_seq_file_cbd_sam)
        if (len(paired_base_set) > 0) and (len(unpaired_read_set) == 0):
            bowtie2_cmd   = 'bowtie2 -x %s -1 %s -2 %s -S %s -p 1 -f --local --all --no-unal'       % (pwd_seq_file_cbd_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_cbd_sam)
        if (len(paired_base_set) == 0) and (len(unpaired_read_set) > 0):
            bowtie2_cmd   = 'bowtie2 -x %s -U %s -S %s -p 1 -f --local --all --no-unal'             % (pwd_seq_file_cbd_index, pwd_seq_file_reads_up, pwd_seq_file_cbd_sam)
        if bowtie2_cmd != '':
            os.system(index_ref_cmd)
            os.system(bowtie2_cmd)
    else:
        index_ref_cmd_16s = 'bowtie2-build --quiet -f %s %s' % (pwd_seq_file_16s, pwd_seq_file_16s_index)
        index_ref_cmd_ctg = 'bowtie2-build --quiet -f %s %s' % (pwd_seq_file_ctg, pwd_seq_file_ctg_index)
        os.system(index_ref_cmd_16s)
        os.system(index_ref_cmd_ctg)

        bowtie2_cmd_16s = ''
        bowtie2_cmd_ctg = ''
        if (len(paired_base_set) > 0) and (len(unpaired_read_set) > 0):
            bowtie2_cmd_16s = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p 6 -f --local --all --no-unal' % (pwd_seq_file_16s_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_reads_up, pwd_seq_file_16s_sam)
            bowtie2_cmd_ctg = 'bowtie2 -x %s -1 %s -2 %s -U %s -S %s -p 6 -f --local --all --no-unal' % (pwd_seq_file_ctg_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_reads_up, pwd_seq_file_ctg_sam)

        if (len(paired_base_set) > 0) and (len(unpaired_read_set) == 0):
            bowtie2_cmd_16s = 'bowtie2 -x %s -1 %s -2 %s -S %s -p 6 -f --local --all --no-unal'       % (pwd_seq_file_16s_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_16s_sam)
            bowtie2_cmd_ctg = 'bowtie2 -x %s -1 %s -2 %s -S %s -p 6 -f --local --all --no-unal'       % (pwd_seq_file_ctg_index, pwd_seq_file_reads_r1, pwd_seq_file_reads_r2, pwd_seq_file_ctg_sam)

        if (len(paired_base_set) == 0) and (len(unpaired_read_set) > 0):
            bowtie2_cmd_16s = 'bowtie2 -x %s -U %s -S %s -p 6 -f --local --all --no-unal'             % (pwd_seq_file_16s_index, pwd_seq_file_reads_up, pwd_seq_file_16s_sam)
            bowtie2_cmd_ctg = 'bowtie2 -x %s -U %s -S %s -p 6 -f --local --all --no-unal'             % (pwd_seq_file_ctg_index, pwd_seq_file_reads_up, pwd_seq_file_ctg_sam)


        if bowtie2_cmd_16s != '':
            os.system(bowtie2_cmd_16s)
        if bowtie2_cmd_ctg != '':
            os.system(bowtie2_cmd_ctg)

    # remove tmp files

