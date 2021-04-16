import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def vis_msa_with_mview(mview_exe, aln_in, page_title, html_out):
    # coloring: any,identity,mismatch,consensus,group
    mview_parameter_str = '-in fasta -moltype dna -colormap CLUSTAL_NUC -coloring any -css on -html head -ruler off -label0 -label4 -label5'
    mview_cmd = '%s %s -title %s %s > %s' % (mview_exe, mview_parameter_str, page_title, aln_in, html_out)
    os.system(mview_cmd)


msa_file                = '/Users/songweizhi/Desktop/vis_msa/combined.aln'
msa_file_updated        = '/Users/songweizhi/Desktop/vis_msa/combined_updated.aln'
msa_file_updated_html   = '/Users/songweizhi/Desktop/vis_msa/combined_updated.html'
ref_id                  = 'ref'
mview_exe               = '/Users/songweizhi/bin/mview'
html_page_title         = 'Reads_linking_16S_and_contig'


mapped_reads_dict = {}
ref_seq = ''
for each_seq in AlignIO.read(msa_file, "fasta"):
    seq_id = each_seq.id
    seq_seq = str(each_seq.seq).upper()
    if seq_id == ref_id:
        ref_seq = seq_seq.upper()
    else:
        seq_id_base = '.'.join(seq_id.split('.')[:-1])
        seq_id_strand = seq_id.split('.')[-1]
        if seq_id_base not in mapped_reads_dict:
            if seq_id_strand == '1':
                mapped_reads_dict[seq_id_base] = [seq_seq, '']
            if seq_id_strand == '2':
                mapped_reads_dict[seq_id_base] = ['', seq_seq]
        else:
            if seq_id_strand == '1':
                mapped_reads_dict[seq_id_base][0] = seq_seq
            if seq_id_strand == '2':
                mapped_reads_dict[seq_id_base][1] = seq_seq


seq_record_list = []

# add ref_seq to seq_record_list
align_record_updated = MultipleSeqAlignment([])
ref_seq_split_by_n = ref_seq.split('N')
ref_seq_split_by_n_updated = []
for segment in ref_seq_split_by_n:
    if ('-' in segment) and (segment == '-' * len(segment)):
        segment = 'N' * len(segment)
    ref_seq_split_by_n_updated.append(segment)
ref_seq_updated = 'N'.join(ref_seq_split_by_n_updated)
seq_record_list.append(SeqRecord(Seq(ref_seq_updated), id=ref_id, description='Reference'))

# add break line
seq_record_list.append(SeqRecord(Seq('='*len(ref_seq_updated)), id='#', description=''))

# add paired reads to seq_record_list
singleton_dict = {}
overlapping_reads_dict = {}
for each_read_base in mapped_reads_dict:

    r1_seq = mapped_reads_dict[each_read_base][0]
    r2_seq = mapped_reads_dict[each_read_base][1]

    if (r1_seq == '') or (r2_seq == ''):
        if not ((r1_seq == '') and (r2_seq == '')):
            singleton_dict[each_read_base] = mapped_reads_dict[each_read_base]
    else:
        overlapping_bps = 0
        for bp1, bp2 in zip(r1_seq, r2_seq):
            if (bp1 != '-') and (bp2 != '-'):
                overlapping_bps += 1

        if overlapping_bps > 0:
            overlapping_reads_dict[each_read_base] = mapped_reads_dict[each_read_base]
        else:
            merge_r1_r2 = ''
            for bp1, bp2 in zip(r1_seq, r2_seq):
                if (bp1 == '-') and (bp2 == '-'):
                    merge_r1_r2 += '-'
                if (bp1 != '-') and (bp2 == '-'):
                    merge_r1_r2 += bp1
                if (bp1 == '-') and (bp2 != '-'):
                    merge_r1_r2 += bp2
            seq_record_list.append(SeqRecord(Seq(merge_r1_r2), id=each_read_base, description='In Pair'))

# add break line
seq_record_list.append(SeqRecord(Seq('='*len(ref_seq_updated)), id='#', description=''))

# add singleton to seq_record_list
for each_singleton in singleton_dict:
    r1_id  = '%s.1' % each_singleton
    r2_id  = '%s.2' % each_singleton
    r1_seq = singleton_dict[each_singleton][0]
    r2_seq = singleton_dict[each_singleton][1]
    if (r1_seq != '') and (r2_seq == ''):
        seq_record_list.append(SeqRecord(Seq(r1_seq), id=r1_id, description='Unpaired'))
    if (r1_seq == '') and (r2_seq != ''):
        seq_record_list.append(SeqRecord(Seq(r2_seq), id=r2_id, description='Unpaired'))

# add break line
seq_record_list.append(SeqRecord(Seq('='*len(ref_seq_updated)), id='#', description=''))

# add overlapped reads to seq_record_list
# to be added

# add break line
seq_record_list.append(SeqRecord(Seq('='*len(ref_seq_updated)), id='#', description=''))

# get updated msa record
align_record_updated = MultipleSeqAlignment(seq_record_list)

# write out updated msa
msa_file_updated_handle = open(msa_file_updated, 'w')
AlignIO.write(align_record_updated, msa_file_updated_handle, 'fasta')
msa_file_updated_handle.close()

# visualize update msa
vis_msa_with_mview(mview_exe, msa_file_updated, html_page_title, msa_file_updated_html)

