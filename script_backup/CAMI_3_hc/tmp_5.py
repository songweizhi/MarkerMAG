import os

def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension

def get_unmapped_mates_seq(sam_file, input_r1_fasta, input_r2_fasta, extracted_seq_file):

    sam_path, sam_basename, sam_ext = sep_path_basename_ext(sam_file)
    reads_to_extract_r1_txt = '%s/%s_unmapped_mates_R1.txt' % (sam_path, sam_basename)
    reads_to_extract_r2_txt = '%s/%s_unmapped_mates_R2.txt' % (sam_path, sam_basename)
    reads_to_extract_r1_fa  = '%s/%s_unmapped_mates_R1.fa'  % (sam_path, sam_basename)
    reads_to_extract_r2_fa  = '%s/%s_unmapped_mates_R2.fa'  % (sam_path, sam_basename)

    prescreening_qualified_reads_set = set()
    prescreening_qualified_reads_base_set = set()
    for each_read in open(sam_file):
        each_read_split = each_read.strip().split('\t')
        if not each_read.startswith('@'):
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            prescreening_qualified_reads_set.add(read_id)
            prescreening_qualified_reads_base_set.add(read_id_base)

    # get id files
    reads_to_extract_r1 = set()
    reads_to_extract_r2 = set()
    for each_read_base in prescreening_qualified_reads_base_set:
        read_r1 = '%s.1' % each_read_base
        read_r2 = '%s.2' % each_read_base
        if (read_r1 not in prescreening_qualified_reads_set) and (read_r2 in prescreening_qualified_reads_set):
            reads_to_extract_r1.add(read_r1)
        if (read_r1 in prescreening_qualified_reads_set) and (read_r2 not in prescreening_qualified_reads_set):
            reads_to_extract_r2.add(read_r2)

    reads_to_extract_r1_txt_handle = open(reads_to_extract_r1_txt, 'w')
    reads_to_extract_r1_txt_handle.write('%s\n' % '\n'.join(reads_to_extract_r1))
    reads_to_extract_r1_txt_handle.close()

    reads_to_extract_r2_txt_handle = open(reads_to_extract_r2_txt, 'w')
    reads_to_extract_r2_txt_handle.write('%s\n' % '\n'.join(reads_to_extract_r2))
    reads_to_extract_r2_txt_handle.close()

    seqtk_extract_r1_read_cmd = 'seqtk subseq %s %s > %s' % (input_r1_fasta, reads_to_extract_r1_txt, reads_to_extract_r1_fa)
    seqtk_extract_r2_read_cmd = 'seqtk subseq %s %s > %s' % (input_r2_fasta, reads_to_extract_r2_txt, reads_to_extract_r2_fa)
    os.system(seqtk_extract_r1_read_cmd)
    os.system(seqtk_extract_r2_read_cmd)

    os.system('cat %s %s > %s' % (reads_to_extract_r1_fa, reads_to_extract_r2_fa, extracted_seq_file))

    # rm tmp files
    os.system('rm %s' % reads_to_extract_r1_txt)
    os.system('rm %s' % reads_to_extract_r2_txt)
    os.system('rm %s' % reads_to_extract_r1_fa)
    os.system('rm %s' % reads_to_extract_r2_fa)


sam_file = '/Users/songweizhi/Desktop/subset_best_match.sam'

reads_to_extract_r1_txt = '/Users/songweizhi/Desktop/subset_best_match_R1.txt'
reads_to_extract_r2_txt = '/Users/songweizhi/Desktop/subset_best_match_R2.txt'
reads_to_extract_r1_fa = '/Users/songweizhi/Desktop/subset_best_match_R1.fa'
reads_to_extract_r2_fa = '/Users/songweizhi/Desktop/subset_best_match_R2.fa'
