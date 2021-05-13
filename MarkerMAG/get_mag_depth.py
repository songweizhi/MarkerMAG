import os
import argparse


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def get_ctg_mean_depth_by_samtools_coverage(index_ref, ref_seq, reads_r1, reads_r2, reads_unpaired, num_threads):

    ref_seq_file_path, ref_seq_file_basename, ref_seq_file_extension = sep_path_basename_ext(ref_seq)

    sam_file        = '%s/%s.sam'          % (ref_seq_file_path, ref_seq_file_basename)
    sam_file_sorted = '%s/%s_sorted.sam'   % (ref_seq_file_path, ref_seq_file_basename)
    coverage_file   = '%s/%s_cov.txt'      % (ref_seq_file_path, ref_seq_file_basename)

    # build reference index
    cmd_bowtie2_build   = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, ref_seq, ref_seq_file_path, ref_seq_file_basename)
    if index_ref is True:
        os.system(cmd_bowtie2_build)

    # mapping
    if reads_unpaired == '':
        cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -1 %s -2 %s -S %s -p %s --all -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
    else:
        cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -1 %s -2 %s -U %s -S %s -p %s --all -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, reads_unpaired, sam_file, num_threads)
    os.system(cmd_bowtie2_mapping)

    # sort mapping
    cmd_samtools_sort = 'samtools sort %s -o %s' % (sam_file, sam_file_sorted)
    os.system(cmd_samtools_sort)

    # get mean depth
    cmd_samtools_coverage = 'samtools coverage %s -o %s' % (sam_file_sorted, coverage_file)
    os.system(cmd_samtools_coverage)

    # remove sam files
    os.system('rm %s' % sam_file)
    os.system('rm %s' % sam_file_sorted)

    # store mean depth into dict
    mean_depth_dict_ctg = {}
    ctg_len_dict = {}
    for each_ctg_depth in open(coverage_file):
        if not each_ctg_depth.startswith('#'):
            ctg_depth_split = each_ctg_depth.strip().split('\t')
            ctg_id = ctg_depth_split[0]
            ctg_len = int(ctg_depth_split[2])
            ctg_depth = float(ctg_depth_split[6])
            mean_depth_dict_ctg[ctg_id] = ctg_depth
            ctg_len_dict[ctg_id] = ctg_len

    return mean_depth_dict_ctg, ctg_len_dict


def get_mag_depth(args):

    output_prefix                       = args['p']
    reads_file_r1                       = args['r1']
    reads_file_r2                       = args['r2']
    mag_folder                          = args['mag']
    num_threads                         = args['t']
    force_overwrite                     = args['force']


    ############################################# create working directory #############################################

    # create working directory
    working_directory = '%s_get_MAG_depth_wd' % output_prefix

    if (os.path.isdir(working_directory) is True) and (force_overwrite is False):
        print('Working directory detected, program exited!')
        exit()
    else:
        force_create_folder(working_directory)

    step_1_wd = '%s/%s_step_1_wd' % (working_directory, output_prefix)
    step_2_wd = '%s/%s_step_2_wd' % (working_directory, output_prefix)
    os.mkdir(step_1_wd)



    print('Round 1: calculating depth for genomes in %s' % mag_folder)

    mean_depth_dict_gnm = {}

    gnm_to_ctg_connector                            = '___C___'



    # get mean depth for contig
    mean_depth_dict_ctg, ctg_len_dict = get_ctg_mean_depth_by_samtools_coverage(True, combined_input_gnms, reads_file_r1, reads_file_r2, '', num_threads)

    # write out ctg depth
    depth_file_ctg_handle = open(depth_file_ctg, 'w')
    for ctg in mean_depth_dict_ctg:
        depth_file_ctg_handle.write('%s\t%s\n' % (ctg, mean_depth_dict_ctg[ctg]))
    depth_file_ctg_handle.close()

    # get mean_depth_dict_gnm
    # if genomic_seq_type == 'mag':
    gnm_len_total_depth_dict = {}
    for ctg in mean_depth_dict_ctg:
        ctg_genome = ctg.split(gnm_to_ctg_connector)[0]
        ctg_len = ctg_len_dict[ctg]
        ctg_depth = mean_depth_dict_ctg[ctg]
        ctg_total_depth = ctg_depth * ctg_len
        if ctg_genome not in gnm_len_total_depth_dict:
            gnm_len_total_depth_dict[ctg_genome] = [ctg_len, ctg_total_depth]
        else:
            gnm_len_total_depth_dict[ctg_genome][0] += ctg_len
            gnm_len_total_depth_dict[ctg_genome][1] += ctg_total_depth

    for each_gnm in gnm_len_total_depth_dict:
        gnm_len = gnm_len_total_depth_dict[each_gnm][0]
        gnm_total_depth = gnm_len_total_depth_dict[each_gnm][1]
        gnm_mean_depth = float("{0:.6f}".format(gnm_total_depth / gnm_len))
        mean_depth_dict_gnm[each_gnm] = gnm_mean_depth

    # write out gnm depth
    depth_file_gnm_handle = open(depth_file_gnm, 'w')
    for gnm in mean_depth_dict_gnm:
        depth_file_gnm_handle.write('%s\t%s\n' % (gnm, mean_depth_dict_gnm[gnm]))
    depth_file_gnm_handle.close()


######################################################### main #########################################################

if __name__ == '__main__':

    link_16s_parser = argparse.ArgumentParser(description='Linking MAGs with marker genes', usage=link_Marker_MAG_usage)

    # specify argument group
    link_16s_parser_input_files = link_16s_parser.add_argument_group("input files")
    link_16s_parser_16s         = link_16s_parser.add_argument_group("16S rRNA gene related parameters")
    link_16s_parser_both_rds    = link_16s_parser.add_argument_group("parameters for linking (round 1 and 2)")
    link_16s_parser_rd1         = link_16s_parser.add_argument_group("parameters for linking (round 1)")
    link_16s_parser_rd2         = link_16s_parser.add_argument_group("parameters for linking (round 2)")
    link_16s_parser_preset      = link_16s_parser.add_argument_group("preset parameters, decide automatically if not specified")
    link_16s_parser_dependency  = link_16s_parser.add_argument_group("provide if dependencies are not in your system path")
    link_16s_parser_others      = link_16s_parser.add_argument_group("program settings")
    link_16s_parser_debug       = link_16s_parser.add_argument_group("for debugging, do NOT specify")

    # input files
    link_16s_parser_input_files.add_argument('-p',          required=False, metavar='',             default=default_prefix, help='output prefix, (default: MyRun_SystemTime)')
    link_16s_parser_input_files.add_argument('-r1',         required=True,  metavar='',                                     help='paired reads r1 (fasta format)')
    link_16s_parser_input_files.add_argument('-r2',         required=True,  metavar='',                                     help='paired reads r2 (fasta format)')
    link_16s_parser_input_files.add_argument('-marker',     required=True,  metavar='',                                     help='marker gene sequences')
    #link_16s_parser_input_files.add_argument('-g',          required=False, metavar='',             default=None,           help='genomic sequences')
    link_16s_parser_input_files.add_argument('-mag',        required=False, metavar='',             default=None,           help='metagenome-assembled-genome (MAG) folder')
    link_16s_parser_input_files.add_argument('-x',          required=False, metavar='',             default='fasta',        help='MAG file extension, (default: %(default)s)')
    link_16s_parser_input_files.add_argument('-depth',      required=False, metavar='', type=float, default=0,              help='minimum depth multiple between 16S and  genomic sequences, a value of no higher than 0.2 is recommended, (default: %(default)s)')
    link_16s_parser_input_files.add_argument('-r16s',       required=False, metavar='',                                     help='16S reads for depth calculation')
    link_16s_parser_input_files.add_argument('-no_polish',  required=False, action="store_true",                            help='skip polishing 16S before linking')

    # 16S rRNA gene related parameters
    link_16s_parser_16s.add_argument('-min_iden_16s',       required=False, metavar='', type=float, default=98,             help='minimum similarity for 16S sequences to be assigned to the same genome, (default: %(default)s)')
    link_16s_parser_16s.add_argument('-min_cov_16s',        required=False, metavar='', type=float, default=30,             help='coverage cutoff for calculating pairwise 16S similarity, (default: %(default)s)')
    link_16s_parser_16s.add_argument('-min_aln_16s',        required=False, metavar='', type=int,   default=500,            help='alignment length cutoff for calculating pairwise 16S similarity, (default: %(default)s)')

    # parameters for both rounds linking
    link_16s_parser_both_rds.add_argument('-mismatch',      required=False, metavar='', type=float, default=1,              help='maximum mismatch percentage, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_len_16s', required=False, metavar='', type=int,   default=75,             help='minimum length aligned to 16S, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_len_ctg', required=False, metavar='', type=int,   default=45,             help='minimum length aligned to ctg, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_M_pct',     required=False, metavar='', type=float, default=35,             help='minimum aligned percentage, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-min_link',      required=False, metavar='', type=int,   default=8,              help='minimum number of linkages to report, (default: %(default)s)')
    link_16s_parser_both_rds.add_argument('-link_num_diff', required=False, metavar='', type=float, default=80,             help='within_gnm_linkage_num_diff, (default: %(default)s)')

    # parameters for 1st round linking
    #link_16s_parser_rd1.add_argument('-min_clp_len',        required=False, metavar='', type=int,   default=45,             help='minimum clipping sequence length (bp), (default: %(default)s)')
    #link_16s_parser_rd1.add_argument('-min_clp_M_len',      required=False, metavar='', type=int,   default=35,             help='minimum aligned clipping sequence length (bp), (default: %(default)s)')

    # parameters for 2nd round linking
    #link_16s_parser_rd2.add_argument('-min_overlap_iden',   required=False, metavar='', type=float, default=100,            help='min_overlap_iden, (default: %(default)s)')
    #link_16s_parser_rd2.add_argument('-min_overlap_cov',    required=False, metavar='', type=float, default=50,             help='min_overlap_cov, (default: %(default)s)')
    #link_16s_parser_rd2.add_argument('-min_overlap_len',    required=False, metavar='', type=int,   default=50,             help='min_overlap_len, (default: %(default)s)')
    #link_16s_parser_rd2.add_argument('-min_overlap_num',    required=False, metavar='', type=int,   default=10,             help='minimum number of overlapping reads for a linkages to be reported, (default: %(default)s)')
    link_16s_parser_rd2.add_argument('-assemble_clp',       required=False, action="store_true",                            help='use clipping mapped reads for mini-assembly')
    link_16s_parser_rd2.add_argument('-mira_tmp',           required=False, default=None,                                   help='tmp dir for mira')
    link_16s_parser_rd2.add_argument('-spades',             required=False, action="store_true",                            help='run spades, instead of Mira')
    link_16s_parser_rd2.add_argument('-link_bias_rd2',      required=False, metavar='', type=float, default=40,             help='max_mini_assembly_link_num_diff_between_ctg_16s, (default: %(default)s)')

    # # preset parameters
    # link_16s_parser_preset.add_argument('-very_sensitive',  required=False, action="store_true",                            help='for greater sensitivity, shortcut for  "min_overlap_iden 99.5 min_overlap_cov 25 min_overlap_len 50 min_overlap_num 3"')
    # link_16s_parser_preset.add_argument('-sensitive',       required=False, action="store_true",                            help='for better sensitivity, shortcut for   "min_overlap_iden 99.5 min_overlap_cov 35 min_overlap_len 50 min_overlap_num 5"')
    # link_16s_parser_preset.add_argument('-specific',        required=False, action="store_true",                            help='for better specificity, shortcut for   "min_overlap_iden 100  min_overlap_cov 55 min_overlap_len 50 min_overlap_num 8"')
    # link_16s_parser_preset.add_argument('-very_specific',   required=False, action="store_true",                            help='for greater specificity, shortcut for  "min_overlap_iden 100  min_overlap_cov 75 min_overlap_len 50 min_overlap_num 10"')

    # program settings
    link_16s_parser_others.add_argument('-bbmap_mem',       required=False, metavar='', type=int,   default=10,             help='bbmap memory allocation (in gigabyte), (default: %(default)s)')
    link_16s_parser_others.add_argument('-t',               required=False, metavar='', type=int,   default=1,              help='number of threads, (default: %(default)s)')
    link_16s_parser_others.add_argument('-tmp',             required=False, action="store_true",                            help='keep temporary files')
    link_16s_parser_others.add_argument('-quiet',           required=False, action="store_true",                            help='not report progress')
    link_16s_parser_others.add_argument('-force',           required=False, action="store_true",                            help='force overwrite existing results')

    # dependency related

    # for debugging
    link_16s_parser_debug.add_argument('-test_mode',        required=False, action="store_true",                            help='only for debugging, do not provide')
    link_16s_parser_debug.add_argument('-rd2_only',         required=False, action="store_true",                            help='run round 2 only')
    link_16s_parser_debug.add_argument('-filtered_sam',     required=False, default=None,                                   help='filtered_sam')

    args = vars(link_16s_parser.parse_args())
    get_mag_depth(args)
