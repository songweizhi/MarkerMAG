import argparse


def keep_only_both_mapped_reads(sam_in, sam_out):
    mapping_dic = {}
    with open(sam_in) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if not each_read.startswith('@'):
                each_read_split = each_read.strip().split('\t')
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                read_strand = read_id.split('.')[-1]

                if read_id_base not in mapping_dic:
                    mapping_dic[read_id_base] = {read_strand}
                else:
                    mapping_dic[read_id_base].add(read_strand)

    mapping_dic_both_mapped = {}
    for read_base in mapping_dic:
        if len(mapping_dic[read_base]) == 2:
            mapping_dic_both_mapped[read_base] = mapping_dic[read_base]

    sam_out_handle = open(sam_out, 'w')
    with open(sam_in) as sorted_sam_opened:
        for each_read in sorted_sam_opened:
            if each_read.startswith('@'):
                sam_out_handle.write(each_read)
            else:
                each_read_split = each_read.strip().split('\t')
                read_id = each_read_split[0]
                read_id_base = '.'.join(read_id.split('.')[:-1])
                if read_id_base in mapping_dic_both_mapped:
                    sam_out_handle.write(each_read)
    sam_out_handle.close()


parser = argparse.ArgumentParser()
parser.add_argument('-in', required=True, help='input file')
parser.add_argument('-out', required=True, help='output file')
args = vars(parser.parse_args())


sam_in               = args['in']
sam_out              = args['out']


keep_only_both_mapped_reads(sam_in, sam_out)