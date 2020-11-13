import os
from Bio import SeqIO


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


free_living_mate_gnm = '/Users/songweizhi/Desktop/777/free_living_mate_ctg.txt'
free_living_mate_16s = '/Users/songweizhi/Desktop/777/free_living_mate_16s.txt'
sam_file             = '/Users/songweizhi/Desktop/777/scaffolds.sam'


reads_to_extract_to_ref_dict_gnm = {}
for i in open(free_living_mate_gnm):
    i_split = i.strip().split('\t')
    reads_to_extract_to_ref_dict_gnm[i_split[0]] = i_split[1]


reads_to_extract_to_ref_dict_16s = {}
for j in open(free_living_mate_16s):
    j_split = j.strip().split('\t')
    reads_to_extract_to_ref_dict_16s[j_split[0]] = j_split[1]


short_ctg_to_reads_dict = {}
for each_line in open(sam_file):
    if not each_line.startswith('@'):
        each_line_split = each_line.strip().split('\t')
        read_id         = each_line_split[0]
        read_id_base    = '.'.join(read_id.split('.')[:-1])
        read_strand     = read_id.split('.')[-1]
        ref_id          = each_line_split[2]
        ref_pos         = int(each_line_split[3])
        cigar           = each_line_split[5]
        cigar_splitted  = cigar_splitter(cigar)
        if (len(cigar_splitted) == 1) and (cigar[-1] == 'M'):
            if ref_id not in short_ctg_to_reads_dict:
                short_ctg_to_reads_dict[ref_id] = [read_id]
            else:
                short_ctg_to_reads_dict[ref_id].append(read_id)


for short_ctg in short_ctg_to_reads_dict:
    short_ctg_mapped_reads = short_ctg_to_reads_dict[short_ctg]

    short_ctg_mapped_reads_mate_linkage_dict = {}
    marker_mate = False
    contig_mate = False
    for mapped_read in short_ctg_mapped_reads:

        mapped_read_mate_ref = ''
        if mapped_read in reads_to_extract_to_ref_dict_gnm:
            mapped_read_mate_ref = reads_to_extract_to_ref_dict_gnm[mapped_read]
        if mapped_read in reads_to_extract_to_ref_dict_16s:
            mapped_read_mate_ref = reads_to_extract_to_ref_dict_16s[mapped_read]

        # check mate mapped to marker gene or genomic sequence
        if 'NODE' in mapped_read_mate_ref:
            contig_mate = True
        if '_m' in mapped_read_mate_ref:
            marker_mate = True


        if mapped_read_mate_ref not in short_ctg_mapped_reads_mate_linkage_dict:
            short_ctg_mapped_reads_mate_linkage_dict[mapped_read_mate_ref] = 1
        else:
            short_ctg_mapped_reads_mate_linkage_dict[mapped_read_mate_ref] += 1

    #print(short_ctg)
    if (len(short_ctg_mapped_reads_mate_linkage_dict) > 1) and (marker_mate is True) and (contig_mate is True):
        print('%s\t%s' % (short_ctg, short_ctg_mapped_reads_mate_linkage_dict))

'''
SP  NODE_8_length_727_cov_18.005952	{'SP___NODE_320_length_67508_cov_39.369843_l': 12, 'SP_m1': 93, 'SP___NODE_65_length_228742_cov_39.546017_l': 20, 'SP___NODE_424_length_51123_cov_38.390616_r': 17, 'SP___NODE_424_length_51123_cov_38.390616_l': 1, 'SP___NODE_107_length_162778_cov_40.127653_l': 25}
FA  NODE_4_length_788_cov_61.869031	{'FA___NODE_21_length_420594_cov_73.571978_l': 160, 'FA___NODE_614_length_33921_cov_65.224827_r': 93, 'FA_m1': 301, 'FA___NODE_124_length_146339_cov_79.710495_l': 87}
TR  NODE_3_length_848_cov_14.722573	{'TR_m1': 98, 'TR___NODE_4_length_1143814_cov_37.243114_l': 38, 'TR___NODE_869_length_20205_cov_36.143970_l': 41, 'HT_m1': 1}
HR  NODE_14_length_606_cov_9.152450	{'HR___NODE_261_length_83213_cov_25.070576_l': 16, 'HR_m1': 61}
PS  NODE_15_length_589_cov_28.988764	{'PS_m1': 134, 'PS___NODE_14_length_483430_cov_23.135121_l': 31, 'PS___NODE_412_length_52305_cov_23.412325_r': 17, 'PS___NODE_42_length_297492_cov_24.514163_l': 22, 'PS___NODE_53_length_264444_cov_22.977193_l': 17, 'AS_m4': 1, 'EC_m3': 3, 'AS_m1': 1, 'AS_m2': 2}
NO  NODE_6_length_759_cov_23.151989	{'NO_m1': 127, 'NO___NODE_658_length_30562_cov_23.290884_l': 70, 'NO___NODE_1253_length_11103_cov_32.877896_r': 31, 'NO___NODE_829_length_21405_cov_27.166698_r': 16}
CA  NODE_12_length_625_cov_36.142105	{'CA___NODE_7_length_713150_cov_99.400944_r': 72, 'CA___NODE_379_length_56028_cov_93.717578_l': 79, 'CA_m1': 138}
HT  NODE_34_length_518_cov_7.732181	{'HT___NODE_556_length_38425_cov_18.514751_r': 7, 'HT_m1': 24, 'HT___NODE_611_length_34055_cov_17.717882_r': 11, 'HT___NODE_1741_length_5855_cov_20.533103': 7, 'HT___NODE_141_length_134940_cov_20.167787_l': 2}
CG  NODE_58_length_459_cov_9.569307	{'CG_m1': 6, 'CG___NODE_363_length_58443_cov_9.958519_l': 4, 'CG_m2': 7, 'CG___NODE_409_length_52773_cov_10.642968_r': 6, 'CG___NODE_1017_length_15878_cov_10.670100_l': 7, 'CG___NODE_100_length_171566_cov_11.802584_r': 1}


NODE_11_length_643_cov_10.226190	{'EC___NODE_1930_length_4753_cov_4.611324': 4, 'AS_m1': 12, 'AS_m4': 6, 'EC___NODE_1271_length_10862_cov_8.428056_l': 2, 'SB_m1': 5, 'AS_m2': 12, 'AS___NODE_723_length_26105_cov_20.427946_r': 10, 'EC_m2': 3, 'PS_m1': 4, 'AS___NODE_1508_length_7793_cov_16.639183': 6, 'AS_m3': 1, 'EC_m3': 2, 'AS___NODE_1896_length_4980_cov_15.570355': 1}
NODE_165_length_329_cov_3.003650	{'PS___NODE_272_length_79871_cov_21.607623_r': 12, 'DG_m2': 1, 'PS___NODE_97_length_175325_cov_21.731432_l': 1}
NODE_421_length_212_cov_1.248408	{'NO_m1': 1, 'AS___NODE_731_length_25859_cov_13.804488_l': 1, 'SP___NODE_449_length_48541_cov_41.452976_l': 1}
NODE_371_length_227_cov_2.180233	{'EC___NODE_1091_length_14060_cov_7.360443_l': 4, 'NO_m1': 1}
NODE_169_length_327_cov_4.584559	{'AS___NODE_1508_length_7793_cov_16.639183': 3, 'SB_m1': 1, 'AS_m2': 4, 'AS_m1': 1, 'AS_m4': 1}
NODE_400_length_219_cov_2.121951	{'AS___NODE_2479_length_2837_cov_17.326384': 2, 'DG_m2': 1}
NODE_436_length_207_cov_3.046053	{'SP_m1': 3, 'HT_m1': 3, 'AS___NODE_2553_length_2616_cov_14.632175': 1}

NODE_427_length_210_cov_2.890323	{'AS___NODE_1151_length_12921_cov_15.667418_r': 7, 'AS___NODE_519_length_41173_cov_13.496303_r': 1, 'FA_m1': 1}
NODE_83_length_425_cov_2.010811	{'AS___NODE_1592_length_7069_cov_9.316082': 10, 'SB_m1': 1}


'''