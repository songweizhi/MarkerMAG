import os
import numpy as np


def get_read_num_and_length(reads_file, tmp_file_location):

    reads_file_line_num = '%s/R1_line_num.txt'  % (tmp_file_location)
    reads_file_sub1000  = '%s/R1_sub1000.fasta' % (tmp_file_location)

    # get the number of paired reads
    os.system('wc -l %s > %s' % (reads_file, reads_file_line_num))
    paired_reads_num = int(int(open(reads_file_line_num).readline().strip().split(' ')[0])/2)

    # subsample 1000 reads
    os.system('seqtk sample -s100 %s 1000 > %s' % (reads_file, reads_file_sub1000))

    read_len_list = []
    for each_seq in open(reads_file_sub1000):
        if each_seq[0] not in ['>', '@', '+']:
            read_len_list.append(len(each_seq.strip()))
    read_len_median = np.median(read_len_list)

    os.system('rm %s' % reads_file_line_num)
    os.system('rm %s' % reads_file_sub1000)

    return paired_reads_num, read_len_median
