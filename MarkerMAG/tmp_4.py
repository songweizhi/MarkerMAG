import os
import argparse
from Bio import SeqIO


def cigar_splitter(cigar):
    # get the position of letters
    letter_pos_list = []
    n = 0
    for each_element in cigar:
        if (each_element.isalpha() is True) or (each_element == '='):
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


# def sep_path_basename_ext(file_in):
#
#     # separate path and file name
#     file_path, file_name = os.path.split(file_in)
#     if file_path == '':
#         file_path = '.'
#
#     # separate file basename and extension
#     file_basename, file_extension = os.path.splitext(file_name)
#
#     return file_path, file_basename, file_extension
#
#
#
# def polish_16s(file_in, file_out_ffn):
#
#     file_out_path, file_out_base, file_out_ext = sep_path_basename_ext(file_out_ffn)
#
#     barrnap_stdout   = '%s/%s.log'    % (file_out_path, file_out_base)
#     file_out_gff     = '%s/%s.gff'    % (file_out_path, file_out_base)
#     file_out_ffn_tmp = '%s/%s_tmp%s' % (file_out_path, file_out_base, file_out_ext)
#
#     barrnap_cmd = 'barrnap --quiet -o %s %s 2> %s > %s' % (file_out_ffn_tmp, file_in, barrnap_stdout, file_out_gff)
#     os.system(barrnap_cmd)
#
#     wrote_id = []
#     file_out_ffn_handle = open(file_out_ffn, 'w')
#     for each_16s in SeqIO.parse(file_out_ffn_tmp, 'fasta'):
#         seq_id = each_16s.id
#         if seq_id.startswith('16S_rRNA::'):
#             seq_id_polished = seq_id[10:].split(':')[0]
#
#             if seq_id_polished not in wrote_id:
#                 file_out_ffn_handle.write('>%s\n' % seq_id_polished)
#                 file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
#                 wrote_id.append(seq_id_polished)
#             else:
#                 file_out_ffn_handle.write('>%s_%s\n' % (seq_id_polished, (wrote_id.count(seq_id_polished) + 1)))
#                 file_out_ffn_handle.write('%s\n' % str(each_16s.seq))
#                 wrote_id.append(seq_id_polished)
#
#     file_out_ffn_handle.close()
#
#     #os.system('rm %s' % file_out_ffn_tmp)
#     #os.system('rm %s.fai' % file_in)
#
#
# file_in      = 'cami_hc_SILVA138_uclust_0.999.fa'
# file_out_ffn = 'cami_hc_SILVA138_uclust_0.999_polished.fa'
# polish_16s(file_in, file_out_ffn)

cigar_List = {'cami_hc_SILVA138_id99_50_subsample_100_534': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_1461': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_1872': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1654': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1336': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1478': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1489': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_1497': '46=54S', 'cami_hc_SILVA138_id99_75_subsample_75_1057': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_1679': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_1107': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1556': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1187': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1368': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1440': '47=53S', 'cami_hc_SILVA138_id99_75_subsample_75_1152': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1430': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_25_740': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1413': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_461': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_580': '47=53S', 'cami_hc_SILVA138_id99_75_subsample_75_638': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1439': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1543': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1372': '47=53S', 'cami_hc_SILVA138_id99_75_subsample_75_1255': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1479': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1298': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1557': '47=53S', 'cami_hc_SILVA138_id99_75_subsample_75_1129': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1438': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1330': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_700': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1831': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1373': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1452': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1251': '47=53S', 'cami_hc_SILVA138_id99_75_subsample_75_1097': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_970': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_1268': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_972': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_10_292': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_1456': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_50_1075': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_2969': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_25_249': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_319': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_321': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_10_177': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_50_1197': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_50_253': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_10_172': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_25_245': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_100_881': '47=53S', 'cami_hc_SILVA138_id99_75_subsample_75_1313': '47=53S', 'cami_hc_SILVA138_id99_50_subsample_25_479': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_732': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_728': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_3623': '46=54S', 'cami_hc_SILVA138_id99_50_subsample_100_2284': '46=54S'}

def check_cigar_all_clp(cigar_list):
    cigar_all_clp = True
    for each_cigar in cigar_list:
        if ('S' not in each_cigar) and ('s' not in each_cigar):
            cigar_all_clp = False

    return cigar_all_clp

print(check_cigar_all_clp(cigar_List.values()))

cigar_List = []
cigar_List = ['1S99=', '1S99=', '1S99=', '6S94=']

def get_min_max_cigar_S(cigar_list):

    # get cigar_S_set
    cigar_S_set = set()
    for each_cigar in cigar_list:
        if ('S' not in each_cigar) and ('s' not in each_cigar):
            cigar_S_set.add(0)
        else:
            cigar_splitted = cigar_splitter(each_cigar)

            if cigar_splitted[0][-1] in ['S', 's']:
                cigar_S_set.add(int(cigar_splitted[0][:-1]))

            if cigar_splitted[-1][-1] in ['S', 's']:
                cigar_S_set.add(int(cigar_splitted[1][:-1]))

    # get min_cigar_S and max_cigar_S
    min_cigar_S = 0
    max_cigar_S = 0
    if len(cigar_S_set) > 0:
        min_cigar_S = min(cigar_S_set)
        max_cigar_S = max(cigar_S_set)

    return min_cigar_S, max_cigar_S


set1 = {1, 2}
set2 = {1, 3}
set3 = {1, 5}

combined = set.union(set1, set2, set3)
print(combined)