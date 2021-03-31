
class MappingRecord:

    def __init__(self):
        self.r1_refs = dict()
        self.r2_refs = dict()


#
# demo_mr_1 = MappingRecord()
#
# demo_mr_1.r1_refs['ref_1'] = '220='
# demo_mr_1.r1_refs['ref_2'] = '111=1X108='
#
#
# print(demo_mr_1.r1_refs)


MappingRecord_dict = {}
demo_sam = '/Users/songweizhi/Desktop/Kelp_SILVA138_id99_assembled_16S_uclust_0.999.sam'
for each_read in open(demo_sam):
    if not each_read.startswith('@'):
        each_read_split = each_read.strip().split('\t')
        cigar = each_read_split[5]
        if cigar != '*':
            read_id = each_read_split[0]
            read_id_base = '.'.join(read_id.split('.')[:-1])
            read_strand = read_id.split('.')[-1]
            ref_id = each_read_split[2]
            ref_pos = each_read_split[3]
            ref_id_with_pos = '%s_pos_%s' % (ref_id, ref_pos)
            read_seq = each_read_split[9]
            #print(each_read.strip())
            #print(read_strand)
            if read_id_base not in MappingRecord_dict:
                MappingRecord_dict[read_id_base] = MappingRecord()

            if read_strand == '1':
                MappingRecord_dict[read_id_base].r1_refs[ref_id_with_pos] = cigar
            if read_strand == '2':
                MappingRecord_dict[read_id_base].r2_refs[ref_id_with_pos] = cigar



#print(MappingRecord_dict)

#print('\n################\n')

for each_mp in MappingRecord_dict:

    print('%s\tr1\t%s' % (each_mp, MappingRecord_dict[each_mp].r1_refs))
    print('%s\tr2\t%s' % (each_mp, MappingRecord_dict[each_mp].r2_refs))
    print()









'''
MappingRecord_example = MappingRecord
MappingRecord_example.r1 = 'aaa'

print(MappingRecord_example.r1)
print(MappingRecord_example.r2)

'''








