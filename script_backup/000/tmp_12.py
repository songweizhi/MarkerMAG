import os
from Bio import SeqIO


def keep_best_blast_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def run_vxtractor_on_input_16s(input_16s, silva_ref_seq, get_16s_domain_wd, num_threads, vxtractor_perl, hmm_folder_bac,
                               hmm_folder_arc, vxtractor_op_combined):
    pwd_silva_16s = '%s/SILVA_16S.fa' % get_16s_domain_wd
    makeblastdb_log = '%s/SILVA_16S_makeblastdb.log' % get_16s_domain_wd
    blast_op = '%s/Input_16S_vs_SILVA.tab' % get_16s_domain_wd
    blast_op_best_hit = '%s/Input_16S_vs_SILVA_best_hit.tab' % get_16s_domain_wd
    input_16s_bac = '%s/Input_16S_bac.fa' % get_16s_domain_wd
    input_16s_arc = '%s/Input_16S_arc.fa' % get_16s_domain_wd
    vxtractor_log_bac = '%s/vxtractor_bac.log' % get_16s_domain_wd
    vxtractor_log_arc = '%s/vxtractor_arc.log' % get_16s_domain_wd
    vxtractor_op_csv_bac = '%s/vxtractor_op_bac.csv' % get_16s_domain_wd
    vxtractor_op_csv_arc = '%s/vxtractor_op_arc.csv' % get_16s_domain_wd
    vxtractor_op_fa_bac = '%s/vxtractor_op_bac.fa' % get_16s_domain_wd
    vxtractor_op_fa_arc = '%s/vxtractor_op_arc.fa' % get_16s_domain_wd

    os.system('cp %s %s' % (silva_ref_seq, pwd_silva_16s))

    # read in ref taxon
    silva_seq_domain_dict = {}
    for each_seq in SeqIO.parse(pwd_silva_16s, 'fasta'):
        seq_description = each_seq.description
        seq_des_split = seq_description.split(' ')
        seq_id = seq_des_split[0]
        taxon_list = ' '.join(seq_description.split(' ')[1:]).split(';')
        silva_seq_domain_dict[seq_id] = taxon_list[0]

    make_blastdb_cmd = 'makeblastdb -in %s -dbtype nucl -parse_seqids -logfile %s' % (pwd_silva_16s, makeblastdb_log)
    os.system(make_blastdb_cmd)
    blastn_cmd = 'blastn -query %s -db %s -out %s -outfmt 6 -num_threads %s' % (input_16s, pwd_silva_16s, blast_op, num_threads)
    os.system(blastn_cmd)

    keep_best_blast_hit(blast_op, blast_op_best_hit)

    query_domain_dict = {}
    for each_hit in open(blast_op_best_hit):
        each_hit_split = each_hit.strip().split('\t')
        query_id = each_hit_split[0]
        ref_id = each_hit_split[1]
        ref_domain = silva_seq_domain_dict[ref_id]
        query_domain_dict[query_id] = ref_domain

    bac_16s_num = 0
    arc_16s_num = 0
    input_16s_bac_handle = open(input_16s_bac, 'w')
    input_16s_arc_handle = open(input_16s_arc, 'w')
    for each_16s in SeqIO.parse(input_16s, 'fasta'):
        seq_domain = query_domain_dict.get(each_16s.id, 'Bacteria')
        if seq_domain == 'Bacteria':
            input_16s_bac_handle.write('>%s\n' % each_16s.id)
            input_16s_bac_handle.write('%s\n' % str(each_16s.seq))
            bac_16s_num += 1
        if seq_domain == 'Archaea':
            input_16s_arc_handle.write('>%s\n' % each_16s.id)
            input_16s_arc_handle.write('%s\n' % str(each_16s.seq))
            arc_16s_num += 1
    input_16s_bac_handle.close()
    input_16s_arc_handle.close()

    vxtractor_cmd_bac = 'perl %s -a -h %s -c %s -o %s %s 2> %s' % (
    vxtractor_perl, hmm_folder_bac, vxtractor_op_csv_bac, vxtractor_op_fa_bac, input_16s_bac, vxtractor_log_bac)
    vxtractor_cmd_arc = 'perl %s -a -h %s -c %s -o %s %s 2> %s' % (
    vxtractor_perl, hmm_folder_arc, vxtractor_op_csv_arc, vxtractor_op_fa_arc, input_16s_arc, vxtractor_log_arc)
    if bac_16s_num > 0:
        os.system(vxtractor_cmd_bac)
    if arc_16s_num > 0:
        os.system(vxtractor_cmd_arc)

    cat_cmd = ''
    if (bac_16s_num > 0) and (arc_16s_num > 0):
        cat_cmd = 'cat %s %s > %s' % (vxtractor_op_csv_bac, vxtractor_op_csv_arc, vxtractor_op_combined)
    elif bac_16s_num > 0:
        cat_cmd = 'cat %s > %s' % (vxtractor_op_csv_bac, vxtractor_op_combined)
    elif arc_16s_num > 0:
        cat_cmd = 'cat %s > %s' % (vxtractor_op_csv_arc, vxtractor_op_combined)
    os.system(cat_cmd)


input_16s               = '/Users/songweizhi/Desktop/get_var_region_wd/MBARC26_SILVA138_polished_polished_min1200bp_c99.0.fasta'
silva_ref_seq           = '/Users/songweizhi/Desktop/get_var_region_wd/SILVA_138.1_one_seq_per_order.fasta'
get_16s_domain_wd       = '/Users/songweizhi/Desktop/get_var_region_wd/get_16s_domain_wd'
num_threads             = 4
vxtractor_perl          = '/Users/songweizhi/Softwares/vxtractor/vxtractor.pl'
hmm_folder_bac          = '/Users/songweizhi/Softwares/vxtractor/HMMs/SSU/bacteria'
hmm_folder_arc          = '/Users/songweizhi/Softwares/vxtractor/HMMs/SSU/archaea'
vxtractor_op_combined   = '/Users/songweizhi/Desktop/get_var_region_wd/Input_16S_vxtractor.csv'

run_vxtractor_on_input_16s(input_16s, silva_ref_seq, get_16s_domain_wd, num_threads, vxtractor_perl, hmm_folder_bac, hmm_folder_arc, vxtractor_op_combined)

