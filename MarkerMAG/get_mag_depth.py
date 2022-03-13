#
# def sep_paired_and_singleton_reads(fasta_in, fasta_out_r1, fasta_out_r2, fasta_out_singleton):
#     reads_pair_dict = {}
#     for read_record in SeqIO.parse(fasta_in, 'fasta'):
#         read_id_base = '.'.join(read_record.id.split('.')[:-1])
#         read_strand = read_record.id.split('.')[-1]
#         if read_id_base not in reads_pair_dict:
#             reads_pair_dict[read_id_base] = {read_strand}
#         else:
#             reads_pair_dict[read_id_base].add(read_strand)
#
#     read_list_paired = set()
#     read_list_singleton = set()
#     for read_base in reads_pair_dict:
#         if len(reads_pair_dict[read_base]) == 1:
#             read_list_singleton.add(read_base)
#         if len(reads_pair_dict[read_base]) == 2:
#             read_list_paired.add(read_base)
#
#     fasta_out_r1_handle = open(fasta_out_r1, 'w')
#     fasta_out_r2_handle = open(fasta_out_r2, 'w')
#     fasta_out_singleton_handle = open(fasta_out_singleton, 'w')
#
#     for read_record in SeqIO.parse(fasta_in, 'fasta'):
#
#         read_id_base = '.'.join(read_record.id.split('.')[:-1])
#         read_strand = read_record.id.split('.')[-1]
#
#         if read_id_base in read_list_singleton:
#             fasta_out_singleton_handle.write('>%s\n' % read_record.id)
#             fasta_out_singleton_handle.write('%s\n' % str(read_record.seq))
#
#         if read_id_base in read_list_paired:
#
#             if read_strand == '1':
#                 fasta_out_r1_handle.write('>%s\n' % read_record.id)
#                 fasta_out_r1_handle.write('%s\n' % str(read_record.seq))
#
#             if read_strand == '2':
#                 fasta_out_r2_handle.write('>%s\n' % read_record.id)
#                 fasta_out_r2_handle.write('%s\n' % str(read_record.seq))
#
#     fasta_out_r1_handle.close()
#     fasta_out_r2_handle.close()
#     fasta_out_singleton_handle.close()
#
#
# def get_read_num_and_length(reads_file, tmp_file_location, seqtk_exe):
#
#     reads_file_line_num = '%s/R1_line_num.txt'  % (tmp_file_location)
#     reads_file_sub1000  = '%s/R1_sub1000.fasta' % (tmp_file_location)
#
#     # get the number of paired reads
#     os.system('wc -l %s > %s' % (reads_file, reads_file_line_num))
#     paired_reads_num = int(int(open(reads_file_line_num).readline().strip().split(' ')[0]) / 2)
#     if reads_file[-1] in ['Q', 'q']:
#         paired_reads_num = int(int(open(reads_file_line_num).readline().strip().split(' ')[0])/4)
#
#     # subsample 1000 reads
#     os.system('%s sample -s100 %s 1000 > %s' % (seqtk_exe, reads_file, reads_file_sub1000))
#
#     read_len_list = []
#     for each_seq in open(reads_file_sub1000):
#         if each_seq[0] not in ['>', '@', '+']:
#             read_len_list.append(len(each_seq.strip()))
#
#     read_len_median = np.median(read_len_list)
#     read_len_max    = np.max(read_len_list)
#
#     os.system('rm %s' % reads_file_line_num)
#     os.system('rm %s' % reads_file_sub1000)
#
#     return paired_reads_num, read_len_median, read_len_max
#
#
# def remove_clp_in_middle(sam_in, sam_out):
#
#     sam_out_handle = open(sam_out, 'w')
#
#     marker_len_dict = {}
#     for each_read in open(sam_in):
#         each_read_split = each_read.strip().split('\t')
#         if each_read.startswith('@'):
#             sam_out_handle.write(each_read)
#
#             marker_id = ''
#             marker_len = 0
#             for each_element in each_read_split:
#                 if each_element.startswith('SN:'):
#                     marker_id = each_element[3:]
#                 if each_element.startswith('LN:'):
#                     marker_len = int(each_element[3:])
#             marker_len_dict[marker_id] = marker_len
#         else:
#             cigar = each_read_split[5]
#
#             # check if clp in the middle
#             if ('S' not in cigar) and ('s' not in cigar):
#                 sam_out_handle.write(each_read)
#             else:
#                 ref_id = each_read_split[2]
#                 ref_pos = int(each_read_split[3])
#                 cigar_splitted = cigar_splitter(cigar)
#                 r1_aligned_len, r1_aligned_pct, r1_clipping_len, r1_clipping_pct, r1_mismatch_pct = get_cigar_stats(cigar_splitted)
#
#                 clip_in_middle = True
#                 if (cigar_splitted[0][-1] in ['S', 's']) and (ref_pos == 1):
#                     clip_in_middle = False
#                 if (cigar_splitted[-1][-1] in ['S', 's']):
#                     if (ref_pos + r1_aligned_len - 1) == marker_len_dict[ref_id]:
#                         clip_in_middle = False
#
#                 if clip_in_middle is False:
#                     sam_out_handle.write(each_read)
#
#     sam_out_handle.close()
#
#
# def get_ctg_mean_depth_by_samtools_coverage(index_ref, ref_seq, reads_r1, reads_r2, reads_unpaired, subsample_rate, num_threads):
#
#     ref_seq_file_path, ref_seq_file_basename, ref_seq_file_extension = sep_path_basename_ext(ref_seq)
#
#     sam_file                                                = '%s/%s.sam'                                               % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_one_end_clp                                    = '%s/%s_one_end_clp.sam'                                   % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_one_end_clp_no_middle                          = '%s/%s_one_end_clp_no_middle.sam'                         % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_one_end_clp_no_middle_best_match               = '%s/%s_one_end_clp_no_middle_best_match.sam'              % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_one_end_clp_no_middle_best_match_low_mismatch  = '%s/%s_one_end_clp_no_middle_best_match_low_mismatch.sam' % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_sorted                                         = '%s/%s_sorted.sam'                                        % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_sorted_log                                         = '%s/%s_sorted.log'                                        % (ref_seq_file_path, ref_seq_file_basename)
#     coverage_file                                           = '%s/%s_cov.txt'                                           % (ref_seq_file_path, ref_seq_file_basename)
#
#     # build reference index
#     cmd_bowtie2_build   = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, ref_seq, ref_seq_file_path, ref_seq_file_basename)
#     if index_ref is True:
#         os.system(cmd_bowtie2_build)
#
#     # mapping
#     if reads_unpaired == '':
#         cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -U %s,%s -S %s -p %s --xeq --local --all --no-unal -N 1 -L 30 -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
#     else:
#         cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -U %s,%s,%s -S %s -p %s --xeq --local --all --no-unal -N 1 -L 30 -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, reads_unpaired, sam_file, num_threads)
#     os.system(cmd_bowtie2_mapping)
#
#     # filter mapping
#     remove_both_ends_clp(sam_file, sam_file_one_end_clp)
#     remove_clp_in_middle(sam_file_one_end_clp, sam_file_one_end_clp_no_middle)
#     keep_best_matches_in_sam_keep_short_M(sam_file_one_end_clp_no_middle, 35, sam_file_one_end_clp_no_middle_best_match)
#     remove_high_mismatch(sam_file_one_end_clp_no_middle_best_match, 2, sam_file_one_end_clp_no_middle_best_match_low_mismatch)
#
#     # sort mapping
#     cmd_samtools_sort = 'samtools sort %s -o %s 2> %s' % (sam_file_one_end_clp_no_middle_best_match_low_mismatch, sam_file_sorted, sam_file_sorted_log)
#     os.system(cmd_samtools_sort)
#
#     # get mean depth
#     cmd_samtools_coverage = 'samtools coverage --ff 4 %s -o %s' % (sam_file_sorted, coverage_file)
#     os.system(cmd_samtools_coverage)
#
#     # remove sam files
#     os.system('rm %s' % sam_file)
#     # os.system('rm %s' % sam_file_sorted)
#
#     # store mean depth into dict
#     mean_depth_dict_ctg = {}
#     ctg_len_dict = {}
#     for each_ctg_depth in open(coverage_file):
#         if not each_ctg_depth.startswith('#'):
#             ctg_depth_split = each_ctg_depth.strip().split('\t')
#             ctg_id = ctg_depth_split[0]
#             ctg_len = int(ctg_depth_split[2])
#             ctg_depth = float(ctg_depth_split[6]) * (1 / subsample_rate)
#             mean_depth_dict_ctg[ctg_id] = ctg_depth
#             ctg_len_dict[ctg_id] = ctg_len
#
#     return mean_depth_dict_ctg, ctg_len_dict
#
#
# def get_ctg_mean_depth_by_samtools_coverage_global(index_ref, ref_seq, reads_r1, reads_r2, reads_unpaired, subsample_rate, num_threads):
#
#     ref_seq_file_path, ref_seq_file_basename, ref_seq_file_extension = sep_path_basename_ext(ref_seq)
#
#     sam_file                          = '%s/%s.sam'                         % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_best_match               = '%s/%s_best_match.sam'              % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_best_match_low_mismatch  = '%s/%s_best_match_low_mismatch.sam' % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_sorted                   = '%s/%s_sorted.sam'                  % (ref_seq_file_path, ref_seq_file_basename)
#     sam_file_sorted_log                   = '%s/%s_sorted.log'                  % (ref_seq_file_path, ref_seq_file_basename)
#     coverage_file                     = '%s/%s_cov.txt'                     % (ref_seq_file_path, ref_seq_file_basename)
#
#     # build reference index
#     cmd_bowtie2_build   = 'bowtie2-build --quiet --threads %s -f %s %s/%s' % (num_threads, ref_seq, ref_seq_file_path, ref_seq_file_basename)
#     if index_ref is True:
#         os.system(cmd_bowtie2_build)
#
#     # mapping
#     cmd_bowtie2_mapping = 'bowtie2 -x %s/%s -U %s,%s -S %s -p %s --all --no-unal -N 1 -L 30 -f --quiet' % (ref_seq_file_path, ref_seq_file_basename, reads_r1, reads_r2, sam_file, num_threads)
#     os.system(cmd_bowtie2_mapping)
#
#     # filter mapping
#     keep_best_matches_in_sam_keep_short_M(sam_file, 35, sam_file_best_match)
#     remove_high_mismatch(sam_file_best_match, 2, sam_file_best_match_low_mismatch)
#
#     # sort mapping
#     cmd_samtools_sort = 'samtools sort %s -o %s 2> %s' % (sam_file_best_match_low_mismatch, sam_file_sorted, sam_file_sorted_log)
#     os.system(cmd_samtools_sort)
#
#     # get mean depth
#     cmd_samtools_coverage = 'samtools coverage --ff 4 %s -o %s' % (sam_file_sorted, coverage_file)
#     os.system(cmd_samtools_coverage)
#
#     # remove sam files
#     os.system('rm %s' % sam_file)
#     # os.system('rm %s' % sam_file_sorted)
#
#     # store mean depth into dict
#     mean_depth_dict_ctg = {}
#     ctg_len_dict = {}
#     for each_ctg_depth in open(coverage_file):
#         if not each_ctg_depth.startswith('#'):
#             ctg_depth_split = each_ctg_depth.strip().split('\t')
#             ctg_id = ctg_depth_split[0]
#             ctg_len = int(ctg_depth_split[2])
#             ctg_depth = float(ctg_depth_split[6]) * (1 / subsample_rate)
#             mean_depth_dict_ctg[ctg_id] = ctg_depth
#             ctg_len_dict[ctg_id] = ctg_len
#
#     return mean_depth_dict_ctg, ctg_len_dict
#
#
# reads_file_r1_subset = '%s/input_R1_subset.fa' % mag_folder_in_wd
# reads_file_r2_subset = '%s/input_R2_subset.fa' % mag_folder_in_wd
#
# mean_depth_dict_gnm = {}
# if min_16s_gnm_multiple > 0:
#
#     if depth_file_mag is None:
#
#         report_and_log(
#             ('Round 1: depth info will be considered (but not provided!) for linking MAGs and 16S rRNA genes'),
#             pwd_log_file, keep_quiet)
#         report_and_log((
#                            'Round 1: depth estimation is time consuming, if you already have them, specify with -depth_16s and -depth_mag'),
#                        pwd_log_file, keep_quiet)
#         report_and_log(('Round 1: calculating MAG depth now, be patient!'), pwd_log_file, keep_quiet)
#
#         # get the number of paired reads
#         input_r1_line_num_file = '%s/R1_line_num.txt' % (mag_folder_in_wd)
#         input_r2_line_num_file = '%s/R2_line_num.txt' % (mag_folder_in_wd)
#         os.system('wc -l %s > %s' % (reads_file_r1, input_r1_line_num_file))
#         os.system('wc -l %s > %s' % (reads_file_r2, input_r2_line_num_file))
#         paired_r1_num = int(int(open(input_r1_line_num_file).readline().strip().split(' ')[0]) / 2)
#         paired_r2_num = int(int(open(input_r2_line_num_file).readline().strip().split(' ')[0]) / 2)
#         if reads_file_r1[-1] in ['Q', 'q']:
#             paired_r1_num = int(int(open(input_r1_line_num_file).readline().strip().split(' ')[0]) / 2)
#             paired_r2_num = int(int(open(input_r2_line_num_file).readline().strip().split(' ')[0]) / 2)
#
#         if paired_r1_num != paired_r2_num:
#             print('Inconsistent number of reads found in r1 and r2, program exited!')
#             exit()
#
#         # get the number of reads paired to subset
#         to_extract_reads_num = round(paired_r1_num * subsample_rate_for_depth_estimation)
#
#         # remember to use the same random seed to keep pairing
#         subsample_r1_cmd = '%s sample -s100 %s %s > %s' % (
#         seqtk_exe, reads_file_r1, to_extract_reads_num, reads_file_r1_subset)
#         subsample_r2_cmd = '%s sample -s100 %s %s > %s' % (
#         seqtk_exe, reads_file_r2, to_extract_reads_num, reads_file_r2_subset)
#
#         # subsample with multiprocessing
#         pool = mp.Pool(processes=2)
#         pool.map(os.system, [subsample_r1_cmd, subsample_r2_cmd])
#         pool.close()
#         pool.join()
#
#         # get mean depth for contig
#         mean_depth_dict_ctg, ctg_len_dict = get_ctg_mean_depth_by_samtools_coverage_global(True, combined_input_gnms,
#                                                                                            reads_file_r1_subset,
#                                                                                            reads_file_r2_subset, '',
#                                                                                            subsample_rate_for_depth_estimation,
#                                                                                            num_threads)
#
#         # write out ctg depth
#         depth_file_ctg_handle = open(depth_file_ctg, 'w')
#         for ctg in mean_depth_dict_ctg:
#             depth_file_ctg_handle.write('%s\t%s\n' % (ctg, mean_depth_dict_ctg[ctg]))
#         depth_file_ctg_handle.close()
#
#         # get mean_depth_dict_gnm
#         gnm_len_total_depth_dict = {}
#         for ctg in mean_depth_dict_ctg:
#             ctg_genome = ctg.split(gnm_to_ctg_connector)[0]
#             ctg_len = ctg_len_dict[ctg]
#             ctg_depth = mean_depth_dict_ctg[ctg]
#             ctg_total_depth = ctg_depth * ctg_len
#             if ctg_genome not in gnm_len_total_depth_dict:
#                 gnm_len_total_depth_dict[ctg_genome] = [ctg_len, ctg_total_depth]
#             else:
#                 gnm_len_total_depth_dict[ctg_genome][0] += ctg_len
#                 gnm_len_total_depth_dict[ctg_genome][1] += ctg_total_depth
#
#         for each_gnm in gnm_len_total_depth_dict:
#             gnm_len = gnm_len_total_depth_dict[each_gnm][0]
#             gnm_total_depth = gnm_len_total_depth_dict[each_gnm][1]
#             gnm_mean_depth = float("{0:.6f}".format(gnm_total_depth / gnm_len))
#             mean_depth_dict_gnm[each_gnm] = gnm_mean_depth
#
#         # write out gnm depth
#         depth_file_gnm_handle = open(depth_file_gnm, 'w')
#         for gnm in mean_depth_dict_gnm:
#             depth_file_gnm_handle.write('%s\t%s\n' % (gnm, mean_depth_dict_gnm[gnm]))
#         depth_file_gnm_handle.close()
#     else:
#         report_and_log(('Round 1: read in provided MAG depth'), pwd_log_file, keep_quiet)
#
#         # read in depth and store in mean_depth_dict_gnm here
#         for each_mag_depth in open(depth_file_mag):
#             each_mag_depth_split = each_mag_depth.strip().split('\t')
#             mag_id = each_mag_depth_split[0]
#             mag_depth = float(each_mag_depth_split[1])
#             mean_depth_dict_gnm[mag_id] = mag_depth
#
