#!/usr/bin/env python3

# Copyright (C) 2020, Weizhi Song, Torsten Thomas.
# songwz03@gmail.com or t.thomas@unsw.edu.au

# MarkerMAG is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MarkerMAG is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import argparse
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '../script_backup'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


parser = argparse.ArgumentParser()

parser.add_argument('-matam',   required=True,                              help='Matam assemblies')
parser.add_argument('-barrnap', required=True,                              help='Barrnap predicted 16S rRNA gene sequences')
parser.add_argument('-iden',    required=False, type=float, default=99.5,   help='Identity cutoff, default: 99.5')
parser.add_argument('-aln',     required=False, type=int,   default=1500,   help='Alignment length cutoff, default: 1500')
parser.add_argument('-noblast', required=False, action="store_true",        help='Skip blastn')

args = vars(parser.parse_args())
matam_16s_seqs   = args['matam']
barrnap_16s_seqs = args['barrnap']
iden_cutoff      = args['iden']
aln_len_cutoff   = args['aln']
skip_blastn      = args['noblast']


matam_16s_seqs_file_path, matam_16s_seqs_file_basename, matam_16s_seqs_file_extension = sep_path_basename_ext(matam_16s_seqs)
matam_16s_blastn = '%s/%s_%s_%sbp_blastn.tab' % (matam_16s_seqs_file_path, matam_16s_seqs_file_basename, iden_cutoff, aln_len_cutoff)
output_file      = '%s/%s_%s_%sbp.fasta' % (matam_16s_seqs_file_path, matam_16s_seqs_file_basename, iden_cutoff, aln_len_cutoff)

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 1'
blast_cmd = 'blastn -query %s -subject %s -out %s %s' % (matam_16s_seqs, barrnap_16s_seqs, matam_16s_blastn, blast_parameters)
if skip_blastn is False:
    os.system(blast_cmd)


query_to_subject_dict = {}
subject_to_query_dict = {}
query_set = set()
subject_set = set()
for match in open(matam_16s_blastn):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_genome = subject.split('_')[0]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    coverage_s = float(aln_len) * 100 / float(subject_len)

    if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff):

        if query not in query_to_subject_dict:
            query_to_subject_dict[query] = {subject_genome}
        else:
            query_to_subject_dict[query].add(subject_genome)

        if subject_genome not in subject_to_query_dict:
            subject_to_query_dict[subject_genome] = {query}
        else:
            subject_to_query_dict[subject_genome].add(query)

        query_set.add(query)
        subject_set.add(subject)


multiple_assign = False
for i in query_to_subject_dict:
    if len(query_to_subject_dict[i]) > 1:
        multiple_assign = True


# get rename dict
rename_dict = {}
for genome in subject_to_query_dict:
    n = 1
    for m16s in subject_to_query_dict[genome]:
        if m16s not in rename_dict:
            rename_dict[m16s] = {'%s_m%s' % (genome, n)}
        else:
            rename_dict[m16s].add('%s_m%s' % (genome, n))
        n += 1


# rename and write out
output_file_handle = open(output_file, 'w')
for m16s_record in SeqIO.parse(matam_16s_seqs, 'fasta'):
    if m16s_record.id in rename_dict:
        new_name_list = rename_dict[m16s_record.id]
        for new_name in new_name_list:
            output_file_handle.write('>%s matam_%s\n' % (new_name, m16s_record.id) )
            output_file_handle.write('%s\n' % m16s_record.seq)
output_file_handle.close()


os.remove(matam_16s_blastn)


# report
print('Identity cutoff: %s' % iden_cutoff)
print('Alignment length cutoff: %s' % aln_len_cutoff)
print('Matched Matam assemblies: %s' % len(query_set))
print('Matched 16S sequences: %s' % len(subject_set))
if multiple_assign is True:
    print('Some/one Matam assembly matched to multiple genomes!')
else:
    print('Good, no Matam assembly matched to multiple genomes!')


#               query   subject
# 1500, 100:    19      42
# 1500, 99.9:   22      52
# 1500, 99.5:   24      66

# 1200, 100:    28      55
# 1200, 99.9:   34      77
# 1200, 99.5:   36      83

# 1000, 100:    30      58
# 1000, 99.9:   36      80
# 1000, 99.5:   38      86

'''
cd /Users/songweizhi/Desktop/ttt
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -matam final_assembly.fa -barrnap combined_16S.ffn

cd /Users/songweizhi/Desktop/111
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -matam 5x_scaffolds.NR.min_500bp.fa -barrnap combined_16S.ffn


cd /Users/songweizhi/Desktop/111
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -matam 10x_scaffolds.NR.min_500bp.fa  -barrnap combined_16S.ffn

python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -matam 10x_scaffolds.NR.min_500bp.fa  -barrnap combined_16S.ffn -iden 99.5 -aln 1000

'''