import os
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def cigar_splitter(cigar):

    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if each_element.isalpha() is True:
            letter_pos_list.append(n)
        n += 1

    # split cigar
    index = 0
    cigar_splitted = []
    while index <= len(letter_pos_list) - 1:
        if index == 0:
            cigar_splitted.append(cigar[:(letter_pos_list[index] + 1)])
        else:
            cigar_splitted.append(cigar[(letter_pos_list[index - 1] + 1):(letter_pos_list[index] + 1)])
        index += 1

    return cigar_splitted


sam_file        = '/Users/songweizhi/Desktop/step_1_unlinked_marker_genes_ends_3000bp.sam'
ref_in          = '/Users/songweizhi/Desktop/step_1_unlinked_marker_genes_ends_3000bp.fasta'
max_gap_to_end  = 300
minCigarM       = 50
end_seq_len     = 3000


# get ref seqs subset
ref_subset_len_dict = {}
for ref_seq in SeqIO.parse(ref_in, 'fasta'):
    ref_seq_id = ref_seq.id
    ref_seq_len = len(ref_seq.seq)
    if ref_seq_len < end_seq_len * 3:
        ref_subset_len_dict[ref_seq_id] = ref_seq_len
    else:
        ref_seq_left_end_id = '%s_l' % ref_seq_id
        ref_seq_right_end_id = '%s_r' % ref_seq_id
        ref_seq_left_end = ref_seq.seq[:end_seq_len]
        ref_seq_right_end = ref_seq.seq[-end_seq_len:]


# parse sam file
all_mapped_reads = set()
qualified_reads_dict = {}
qualified_reads_to_ref_dict = {}
for each_line in open(sam_file):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        cigar = each_line_split[5]

        if cigar != '*':
            read_id = each_line_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_line_split[2]
            ref_pos = int(each_line_split[3])
            cigar_splitted = cigar_splitter(cigar)
            all_mapped_reads.add(read_id)
            qualified_mapping = False

            if ref_id[-2:] == '_l':
                if ref_pos <= max_gap_to_end:
                    if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                        qualified_mapping = True


                    elif (len(cigar_splitted) == 2) and (cigar_splitted[-1][-1] == 'M') and (
                            cigar_splitted[0][-1] == 'S') and (int(cigar_splitted[-1][:-1]) >= minCigarM) and (
                            ref_pos == 1):
                        qualified_mapping = True


            elif ref_id[-2:] == '_r':

                if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                    ref_pos_end = ref_pos + int(cigar[:-1])
                    if (end_seq_len - ref_pos_end) <= max_gap_to_end:
                        qualified_mapping = True

                elif (len(cigar_splitted) == 2) and (cigar_splitted[0][-1] == 'M') and (
                        cigar_splitted[1][-1] == 'S') and (int(cigar_splitted[0][:-1]) >= minCigarM):
                    if (ref_pos + int(cigar_splitted[0][:-1]) - 1) == end_seq_len:
                        qualified_mapping = True

            else:
                ref_len = ref_subset_len_dict[ref_id]
                if (len(cigar_splitted) == 1) and (cigar[-1] == 'M') and (int(cigar[:-1]) >= minCigarM):
                    ref_pos_end = ref_pos + int(cigar[:-1])

                    # left side: 220M
                    if ref_pos <= max_gap_to_end:
                        qualified_mapping = True
                        # good
                        print(each_line_split)


                    # right side: 116M
                    elif (ref_len - ref_pos_end) <= max_gap_to_end:
                        qualified_mapping = True
                        # good

                if len(cigar_splitted) == 2:

                    # left side: ['6S', '214M']
                    if (cigar_splitted[-1][-1] == 'M') and (cigar_splitted[0][-1] == 'S') and (int(cigar_splitted[-1][:-1]) >= minCigarM) and (ref_pos == 1):
                        qualified_mapping = True
                        # good

                    # right side: ['212M', '8S']
                    elif (cigar_splitted[-0][-1] == 'M') and (cigar_splitted[1][-1] == 'S') and (int(cigar_splitted[0][:-1]) >= minCigarM):
                        if (ref_pos + int(cigar_splitted[0][:-1]) - 1) == ref_len:
                            qualified_mapping = True
                            # good
                # allow one mismatch
                if len(cigar_splitted) == 3:
                    if (cigar_splitted[0][-1] == 'M') and (cigar_splitted[2][-1] == 'M'):
                        left_m_len = int(cigar_splitted[0][:-1])
                        mismatch_len = int(cigar_splitted[1][:-1])
                        right_m_len = int(cigar_splitted[2][:-1])
                        read_len = left_m_len + right_m_len + mismatch_len
                        left_m_pct = left_m_len * 100 / read_len
                        right_m_pct = right_m_len * 100 / read_len
                        if (left_m_len >= 30) and (right_m_len >= 30) and (left_m_pct >= 25) and (right_m_pct >= 25) and (mismatch_len == 1):
                            qualified_mapping = True

            if qualified_mapping is True:
                qualified_reads_to_ref_dict[read_id] = ref_id

                if read_id_base not in qualified_reads_dict:
                    qualified_reads_dict[read_id_base] = [read_strand]
                else:
                    qualified_reads_dict[read_id_base].append(read_strand)

reads_to_extract_to_ref_dict = {}
for qualified_read in qualified_reads_dict:
    read_strand = qualified_reads_dict[qualified_read]
    if len(read_strand) == 1:

        mapped_mate = ''
        mate_to_extract = ''
        if read_strand == ['1']:
            mapped_mate = '%s.1' % (qualified_read)
            mate_to_extract = '%s.2' % (qualified_read)
        if read_strand == ['2']:
            mapped_mate = '%s.2' % (qualified_read)
            mate_to_extract = '%s.1' % (qualified_read)

        if mate_to_extract not in all_mapped_reads:
            reads_to_extract_to_ref_dict[mate_to_extract] = qualified_reads_to_ref_dict[mapped_mate]


print(reads_to_extract_to_ref_dict)
print(len(reads_to_extract_to_ref_dict))

'''
module load java/8u201-jdk
module load bbmap/38.51
cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp/test_cigar_fmt
bbmap.sh ref=test.fasta in=../Kelp_R1.fasta in2=../Kelp_R2.fasta outm=test_1.3.sam local=t nodisk=t ambiguous=all sam=1.3 keepnames=t saa=f silent=true threads=6 -Xmx10g
bbmap.sh ref=test.fasta in=../Kelp_R1.fasta in2=../Kelp_R2.fasta outm=test_1.4.sam local=t nodisk=t ambiguous=all keepnames=t saa=f silent=true threads=6 -Xmx10g



'''