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
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


parser = argparse.ArgumentParser()

parser.add_argument('-m',       required=True,                              help='Matam assemblies')
parser.add_argument('-c',       required=True,                              help='Control, e.g. Barrnap predicted 16S rRNA gene sequences')
parser.add_argument('-i',       required=False, type=float, default=99.5,   help='Identity cutoff, default: 99.5')
parser.add_argument('-l',       required=False, type=int,   default=1500,   help='Alignment length cutoff, default: 1500')
parser.add_argument('-noblast', required=False, action="store_true",        help='Skip blastn')

args = vars(parser.parse_args())
matam_16s_seqs   = args['m']
barrnap_16s_seqs = args['c']
iden_cutoff      = args['i']
aln_len_cutoff   = args['l']
skip_blastn      = args['noblast']


matam_16s_seqs_file_path, matam_16s_seqs_file_basename, matam_16s_seqs_file_extension = sep_path_basename_ext(matam_16s_seqs)
matam_16s_blastn = '%s/%s_%s_%sbp_blastn.tab' % (matam_16s_seqs_file_path, matam_16s_seqs_file_basename, iden_cutoff, aln_len_cutoff)
output_file      = '%s/%s_%s_%sbp.fasta' % (matam_16s_seqs_file_path, matam_16s_seqs_file_basename, iden_cutoff, aln_len_cutoff)
output_file_txt  = '%s/%s_%s_%sbp.txt'   % (matam_16s_seqs_file_path, matam_16s_seqs_file_basename, iden_cutoff, aln_len_cutoff)

blast_parameters = '-evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -task blastn -num_threads 1'
blast_cmd = 'blastn -query %s -subject %s -out %s %s' % (matam_16s_seqs, barrnap_16s_seqs, matam_16s_blastn, blast_parameters)
if skip_blastn is False:
    os.system(blast_cmd)


# get the number of control sequences
genome_id_set = set()
control_seq_id_list = set()
for seq_record in SeqIO.parse(barrnap_16s_seqs, 'fasta'):
    control_seq_id_list.add(seq_record.id)
    genome_id_set.add(seq_record.id.split('_16S_')[0])


query_to_subject_dict = {}
query_to_subject_genome_dict = {}
subject_genome_to_query_dict = {}
matched_matam_set = set()
matched_control_set = set()
matched_control_set_genome_level = set()
for match in open(matam_16s_blastn):
    match_split = match.strip().split('\t')
    query = match_split[0]
    subject = match_split[1]
    subject_genome = subject.split('_16S_')[0]
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    coverage_s = float(aln_len) * 100 / float(subject_len)

    if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff):

        # get query_to_subject_dict
        if query not in query_to_subject_dict:
            query_to_subject_dict[query] = {subject}
        else:
            query_to_subject_dict[query].add(subject)

        # get query_to_subject_genome_dict
        if query not in query_to_subject_genome_dict:
            query_to_subject_genome_dict[query] = {subject_genome}
        else:
            query_to_subject_genome_dict[query].add(subject_genome)

        # get subject_genome_to_query_dict
        if subject_genome not in subject_genome_to_query_dict:
            subject_genome_to_query_dict[subject_genome] = {query}
        else:
            subject_genome_to_query_dict[subject_genome].add(query)

        matched_matam_set.add(query)
        matched_control_set.add(subject)
        matched_control_set_genome_level.add(subject_genome)


multiple_assign = False
for i in query_to_subject_genome_dict:
    if len(query_to_subject_genome_dict[i]) > 1:
        multiple_assign = True


# get rename dict
rename_dict = {}
for genome in subject_genome_to_query_dict:
    n = 1
    for m16s in subject_genome_to_query_dict[genome]:
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


output_file_txt_handle = open(output_file_txt, 'w')
for matched_assembly in rename_dict:
    new_names = rename_dict[matched_assembly]
    matched_to_list = query_to_subject_dict[matched_assembly]
    for new_name in new_names:
        output_file_txt_handle.write('%s\t%s\n' % (new_name, ','.join(matched_to_list)))
output_file_txt_handle.close()


######################################################## report ########################################################

# get unrecovered genome
unrecovered_genome_list = []
for genome_id in genome_id_set:
    if genome_id not in matched_control_set_genome_level:
        unrecovered_genome_list.append(genome_id)

# get unrecovered control set
unrecovered_control_set = set()
for i in control_seq_id_list:
    if i not in matched_control_set:
        unrecovered_control_set.add(i)

unrecovered_str = 'Unrecovered(%s):%s' % (len(unrecovered_control_set), ','.join(sorted([i for i in unrecovered_control_set])))
unrecovered_str_genome_level = 'Unrecovered(%s):%s' % (len(unrecovered_genome_list), ','.join(sorted(unrecovered_genome_list)))

print('Iden\tLen\tControl\tMatam\tUnrecovered')
print('Marker\t%s%s\t%sbp\t%s/%s\t%s\t%s' % (iden_cutoff, '%', aln_len_cutoff, len(matched_control_set), len(control_seq_id_list), len(matched_matam_set), unrecovered_str))
print('Genome\t%s%s\t%sbp\t%s/%s\t%s\t%s' % (iden_cutoff, '%', aln_len_cutoff, len(matched_control_set_genome_level), len(genome_id_set), len(matched_matam_set), unrecovered_str_genome_level))

if multiple_assign is True:
    print('Some/one Matam assembly matched to multiple genomes!')
else:
    print('Good, no Matam assembly matched to multiple genomes!')


########################################################################################################################


'''
cd /Users/songweizhi/Desktop/111
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -c /srv/scratch/z5039045/MarkerMAG_wd/genome_selection_3/combined_16S.ffn -m scaffolds.NR.min_500bp.abd.fa -i 99.9 -a 1300

python3 /srv/scratch/z5039045/MarkerMAG_wd/rename_matam_assemblies.py -c /srv/scratch/z5039045/MarkerMAG_wd/genome_selection_3/combined_16S.ffn -m scaffolds.NR.min_500bp.abd.fa -i 99.9 -a 1300

cd /Users/songweizhi/Desktop/999
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -c combined_16S.ffn -m combined_all_depth_assemblies_iden99.5_uniq.fasta -i 99.5 -l 1300

'''