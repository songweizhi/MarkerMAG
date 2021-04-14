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
    genome_id_set.add('_'.join(seq_record.id.split('_')[:-1]))


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
    subject_genome = subject.split('_')[0]
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
multiple_assignmeng_dict = {}
for i in query_to_subject_genome_dict:
    if len(query_to_subject_genome_dict[i]) > 1:
        multiple_assignmeng_dict[i] = query_to_subject_genome_dict[i]
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
    for each in multiple_assignmeng_dict:
        print('%s\t%s' % (each, ','.join(multiple_assignmeng_dict[each])))
else:
    print('Good, no Matam assembly matched to multiple genomes!')


########################################################################################################################


'''
cd /Users/songweizhi/Desktop/111
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -c /srv/scratch/z5039045/MarkerMAG_wd/genome_selection_3/combined_16S.ffn -m scaffolds.NR.min_500bp.abd.fa -i 99.9 -a 1300

python3 /srv/scratch/z5039045/MarkerMAG_wd/rename_matam_assemblies.py -c /srv/scratch/z5039045/MarkerMAG_wd/genome_selection_3/combined_16S.ffn -m scaffolds.NR.min_500bp.abd.fa -i 99.9 -a 1300

cd /Users/songweizhi/Desktop/999
python3 /Users/songweizhi/PycharmProjects/MarkerMAG/script_backup/rename_matam_assemblies.py -c combined_16S.ffn -m combined_all_depth_assemblies_iden99.5_uniq.fasta -i 99.5 -l 1300


DG_m4	DA	634	S1

MarkerGene,GenomicSeq,Number
MarkerGene__FP_m2,GenomicSeq__FP,1014
MarkerGene__DM_m4,GenomicSeq__DM,644
MarkerGene__DG_m4,GenomicSeq__DA,634
MarkerGene__DM_m2,GenomicSeq__DM,633
MarkerGene__DM_m8,GenomicSeq__DM,581
MarkerGene__DA_m5,GenomicSeq__DA,550
MarkerGene__MS_m1,GenomicSeq__MS,499
MarkerGene__DM_m7,GenomicSeq__DM,459
MarkerGene__DG_m7,GenomicSeq__DG,319


DG_m4	DG_04004	99.586	1449	5	1	1	1449	81	1528	0.0	2641
DG_m4	DG_01205	99.241	1449	11	0	1	1449	81	1529	0.0	2615
DG_m4	DG_01680	99.580	1430	6	0	20	1449	152	1581	0.0	2608
DG_m4	DG_00008	99.034	1449	14	0	1	1449	81	1529	0.0	2599
DG_m4	DG_03978	99.231	1430	11	0	20	1449	152	1581	0.0	2580
DG_m4	DG_01668	98.281	1454	15	8	1	1449	81	1529	0.0	2538
DG_m4	DG_01491	98.395	1433	17	5	20	1449	152	1581	0.0	2514
DG_m4	DG_03995	99.102	1337	12	0	20	1356	152	1488	0.0	2403

DG_m4	DA_03472	86.305	1475	149	28	20	1449	100	1566	0.0	1555
DG_m4	DA_02778	86.305	1475	149	28	20	1449	191	1657	0.0	1555
DG_m4	DA_00795	86.305	1475	149	28	20	1449	100	1566	0.0	1555
DG_m4	DA_00160	86.305	1475	149	28	20	1449	100	1566	0.0	1555
DG_m4	DA_00117	86.305	1475	149	28	20	1449	100	1566	0.0	1555
DG_m4	DA_03192	86.237	1475	150	29	20	1449	100	1566	0.0	1550
DG_m4	DA_02224	86.237	1475	150	28	20	1449	100	1566	0.0	1550
DG_m4	DA_00092	86.237	1475	150	29	20	1449	100	1566	0.0	1550
DG_m4	DA_03459	86.111	1476	150	28	20	1449	100	1566	0.0	1539

NODE_798_length_22778_cov_480.292743


to remove:
MBARC26_63374472.1	81	DG_m4	1367	16	14S17=1X2=1X64=30S	DA_m5	1193	0	CCCGCAAGGGAGCTAGCCGTCGAAGGTGGGGCCGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCTAAGGAGAACGGTTTAGAGCTTAGGCTTTAA	*	NM:i:2	AM:i:16	NH:i:1



wrong:
MBARC26_63374472.1	81	DG_m4	1367	16	14S17=1X2=1X64=30S	DA_m5	1193	0	CCCGCAAGGGAGCTAGCCGTCGAAGGTGGGGCCGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCTAAGGAGAACGGTTTAGAGCTTAGGCTTTAA	*	NM:i:2	AM:i:16	NH:i:1
MBARC26_63374472.2	161	DA_m5	1193	3	123=	DG_m4	1367	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	XT:A:R	NM:i:0	AM:i:3	NH:i:3
MBARC26_63374472.2	353	DA_m4	1095	3	123=	DG_m4	1367	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	NM:i:0	AM:i:3	NH:i:3
MBARC26_63374472.2	353	DA_m3	1040	3	123=	DG_m4	1367	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	NM:i:0	AM:i:3	NH:i:3

MBARC26_48021260.2	145	DG_m4	1373	13	11=1X2=1X64=50S	DA_m5	1111	0	CGAAGGTGGGGCCGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCTAAGGAGAACGGTTTAGAGCTTAGGCTTTAAACGAACATCCTATGGGTCGA	*	NM:i:2	AM:i:13	NH:i:1
MBARC26_48021260.1	97	DA_m5	1111	3	130=	DG_m4	1373	0	ATAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATC	*	XT:A:R	NM:i:0	AM:i:3	NH:i:3
MBARC26_48021260.1	353	DA_m4	1013	3	130=	DG_m4	1373	0	ATAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATC	*	NM:i:0	AM:i:3	NH:i:3
MBARC26_48021260.1	353	DA_m3	958	3	130=	DG_m4	1373	0	ATAGACTGCCGGTGACAAACCGGAGGAAGGTGGGGATGACGTCAAATCATCATGCCCCTTATGTCCTGGGCTACACACGTGCTACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATC	*	NM:i:0	AM:i:3	NH:i:3

MBARC26_47678810.1	83	DG_m4	1367	16	12S17=1X2=1X64=33S	=	1235	-217	CGCAAGGGAGCTAGCCGTCGAAGGTGGGGCCGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCTAAGGAGAACGGTTTAGAGCTTAGGCTTTAAACG	*	NM:i:2	AM:i:16	NH:i:1
MBARC26_47678810.2	163	DG_m4	1235	44	92=	=	1367	217	CATGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAG	*	NM:i:0	AM:i:44	NH:i:5
MBARC26_47678810.2	353	DG_m8	1351	44	92=	DG_m4	1367	0	CATGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAG	*	NM:i:0	AM:i:44	NH:i:5
MBARC26_47678810.2	353	DM_m1	740	44	92=	DG_m4	1367	0	CATGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAG	*	NM:i:0	AM:i:44	NH:i:5
MBARC26_47678810.2	353	DM_m4	1316	44	92=	DG_m4	1367	0	CATGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAG	*	NM:i:0	AM:i:44	NH:i:5
MBARC26_47678810.2	353	DG_m9	1299	44	92=	DG_m4	1367	0	CATGAAGTCGGAATCGCTAGTAATCGCAGGTCAGCATACTGCGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACGAAAG	*	NM:i:0	AM:i:44	NH:i:5

MBARC26_9544440.1	83	DG_m4	1380	14	4=1X2=1X64=41S	=	1334	-118	GGGGCCGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGGTGCGGCTGGATCACCTCCTTTCTAAGGAGAACGGTTTAGAGCTTAGGCTTTAAACGAACATCCT	*	NM:i:2	AM:i:14	NH:i:1
MBARC26_9544440.2	163	DG_m4	1334	25	17=1X5=3X1=2I3=3X17=1X2=1X41=	=	1380	118	CACCCGAAGCCGGTGAGGTAACCCGCAAGGGAGCTAGCCGTCGAAGGTGGGGCCGATAATTGGGGTGAAGTCGTAACAAGGTAGCCGTATCGGAAGG	*	NM:i:11	AM:i:25	NH:i:1




previous
MBARC26_63374472.1	69	DA_m1	1233	0	*	=	1233	0	TTAAAGCCTAAGCTCTAAACCGTTCTCCTTAGAAAGGAGGTGATCCAGCCGCACCTTCCGATACGGCTACCTTGTTACGACTTCACCCCAATTATCGGCCCCACCTTCGACGGCTAGCTCCCTTGCGGG	*
MBARC26_63374472.2	137	DA_m1	1233	3	123=	=	1233	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	XT:A:R	NM:i:0	AM:i:0	NH:i:5
MBARC26_63374472.2	329	DA_m3	1173	3	123=	=	1173	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	NM:i:0	AM:i:0	NH:i:5
MBARC26_63374472.2	329	DA_m2	1332	3	123=	=	1332	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	NM:i:0	AM:i:0	NH:i:5
MBARC26_63374472.2	329	DA_m7	1241	3	123=	=	1241	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	NM:i:0	AM:i:0	NH:i:5
MBARC26_63374472.2	329	DA_m6	1233	3	123=	=	1233	0	TACAATGGCCGGTACAGACGGAAGCGAAGCCGCGAGGTGAAGCCAATCCGAGAAAGCCGGTCTCAGTTCGGATTGCAGGCTGCAACTCGCCTGCATGAAGTCGGAATCGCTAGTAATCGCAGG	*	NM:i:0	AM:i:0	NH:i:5

'''