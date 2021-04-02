
pwd_samfile = '/Users/songweizhi/Desktop/input_reads_to_16S.sam'
pwd_samfile = 'round_1_unlinked_gnm.sam'


id_to_seq_dict = {}
for each_read in open(pwd_samfile):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        read_id = each_read_split[0]
        read_seq = each_read_split[9]
        if read_id not in id_to_seq_dict:
            id_to_seq_dict[read_id] = {read_seq}
        else:
            id_to_seq_dict[read_id].add(read_seq)

for each_read in id_to_seq_dict:
    if len(id_to_seq_dict[each_read]) > 1:
        for each_seq in id_to_seq_dict[each_read]:
            print('%s\t%s' % (each_read, each_seq))
        print()
