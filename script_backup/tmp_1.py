

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


pwd_samfile                         = '/Users/songweizhi/Desktop/test_bowtie/f5_16S.sam'
clipping_reads_matched_part         = '/Users/songweizhi/Desktop/test_bowtie/clipping_reads_matched_part.txt'
clipping_reads_not_matched_part_seq = '/Users/songweizhi/Desktop/test_bowtie/clipping_reads_not_matched_part.fasta'
min_cigar_M = 30
min_cigar_S = 30

'''
cd /Users/songweizhi/Desktop/test_bowtie
export PATH=/Users/songweizhi/Softwares/bowtie2:$PATH

BioSAK Reads_simulator -r f5.fna -n 200000 -l 150 -i 200 -split

bowtie2-build -f f5_16S.ffn f5_16S 
bowtie2 -x f5_16S -1 f5_R1.fasta -2 f5_R2.fasta -S f5_16S.sam --local --no-unal -f 
bowtie2 -x f5_16S -1 f5_R1.fasta -2 f5_R2.fasta -S f5_16S_a.sam --local --no-unal -f -a 

'''

# export clipping mapped reads and perfectly mapped reads
clipping_mapped_reads_list = set()
clipping_reads_mapped_part_dict = {}
perfectly_mapped_reads_dict = {}
clipping_reads_matched_part_handle = open(clipping_reads_matched_part, 'w')
clipping_reads_not_matched_part_seq_handle = open(clipping_reads_not_matched_part_seq, 'w')
for each_read in open(pwd_samfile):

    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_id_base = '.'.join(read_id.split('.')[:-1])
        read_strand = read_id.split('.')[-1]
        ref_id = each_read_split[2]
        ref_id_with_prefix = 'MarkerGene__%s' % each_read_split[2]
        ref_pos = int(each_read_split[3])
        cigar = each_read_split[5]
        read_seq = each_read_split[9]
        cigar_splitted = cigar_splitter(cigar)

        # for perfectly mapped reads
        if ('M' in cigar) and (len(cigar_splitted) == 1):
            if read_id_base not in perfectly_mapped_reads_dict:
                perfectly_mapped_reads_dict[read_id_base] = {read_strand: [ref_id_with_prefix]}
            else:
                if read_strand not in perfectly_mapped_reads_dict[read_id_base]:
                    perfectly_mapped_reads_dict[read_id_base][read_strand] = [ref_id_with_prefix]
                else:
                    perfectly_mapped_reads_dict[read_id_base][read_strand].append(ref_id_with_prefix)

        # for clipping mapped reads
        if ('S' in cigar) and (len(cigar_splitted) == 2):
            cigar_M_len = 0
            cigar_S_len = 0
            split_pos = 0
            if cigar_splitted[0][-1] == 'M':
                cigar_M_len = int(cigar_splitted[0][:-1])
                cigar_S_len = int(cigar_splitted[1][:-1])
                split_pos = ref_pos + cigar_M_len
            if cigar_splitted[1][-1] == 'M':
                cigar_M_len = int(cigar_splitted[1][:-1])
                cigar_S_len = int(cigar_splitted[0][:-1])
                split_pos = ref_pos

            if (cigar_M_len >= min_cigar_M) and (cigar_S_len >= min_cigar_S):
                read_seq_left = read_seq[: int(cigar_splitted[0][:-1])]
                read_seq_right = read_seq[-int(cigar_splitted[1][:-1]):]

                if cigar_splitted[0][-1] == 'M':
                    clipping_reads_matched_part_handle.write('%s_l\t%s\t%s\n' % (read_id, ref_id, split_pos))
                    clipping_reads_not_matched_part_seq_handle.write('>%s_r\n' % read_id)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_right + '\n')
                    if ('%s_l' % read_id) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_l' % read_id)].append(ref_id_with_prefix)

                if cigar_splitted[1][-1] == 'M':
                    clipping_reads_matched_part_handle.write('%s_r\t%s\t%s\n' % (read_id, ref_id, split_pos))
                    clipping_reads_not_matched_part_seq_handle.write('>%s_l\n' % read_id)
                    clipping_reads_not_matched_part_seq_handle.write(read_seq_left + '\n')
                    if ('%s_r' % read_id) not in clipping_reads_mapped_part_dict:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id)] = [ref_id_with_prefix]
                    else:
                        clipping_reads_mapped_part_dict[('%s_r' % read_id)].append(ref_id_with_prefix)

                clipping_mapped_reads_list.add(read_id)
clipping_reads_not_matched_part_seq_handle.close()

print(perfectly_mapped_reads_dict)
print(clipping_reads_mapped_part_dict)


perfectly_mapped_read_singleton_dict = {}
for perfectly_mapped_read in perfectly_mapped_reads_dict:
    current_value = perfectly_mapped_reads_dict[perfectly_mapped_read]
    if len(current_value) == 1:
        perfectly_mapped_read_singleton_dict[perfectly_mapped_read] = current_value


print(perfectly_mapped_read_singleton_dict)