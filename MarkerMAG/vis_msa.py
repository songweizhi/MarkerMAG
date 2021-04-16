import os
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def mview_linkage(msa_file, plot_title, mview_exe):

    msa_path, msa_basename, msa_ext = sep_path_basename_ext(msa_file)
    msa_file_mview         = '%s/%s_MView.aln' % (msa_path, msa_basename)
    msa_file_mviewd_html   = '%s/%s_MView.html' % (msa_path, msa_basename)

    mapped_reads_dict = {}
    ref_id = ''
    ref_seq = ''
    current_line = 0
    for each_seq in AlignIO.read(msa_file, "fasta"):
        seq_id = each_seq.id
        seq_seq = str(each_seq.seq).upper()
        if current_line == 0:
            ref_id = seq_id
            ref_seq = seq_seq
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

        current_line += 1

    seq_record_list = []

    # add ref_seq to seq_record_list
    align_record_mview = MultipleSeqAlignment([])
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
    align_record_mview = MultipleSeqAlignment(seq_record_list)

    # write out updated msa
    msa_file_mview_handle = open(msa_file_mview, 'w')
    AlignIO.write(align_record_mview, msa_file_mview_handle, 'fasta')
    msa_file_mview_handle.close()

    # visualize update msa
    # coloring: any,identity,mismatch,consensus,group
    mview_parameter_str = '-in fasta -moltype dna -colormap CLUSTAL_NUC -coloring any -css on -html head -ruler off -label0 -label4 -label5'
    mview_cmd = '%s %s -title %s %s > %s' % (mview_exe, mview_parameter_str, plot_title, msa_file_mview, msa_file_mviewd_html)
    os.system(mview_cmd)


msa_file                = '/Users/songweizhi/Desktop/vis_msa/combined.aln'
mview_exe               = '/Users/songweizhi/bin/mview'
html_page_title         = 'Reads_linking_16S_and_contig'

mview_linkage(msa_file, html_page_title, mview_exe)
