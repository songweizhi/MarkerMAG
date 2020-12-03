import os
import glob


genome_folder       = '/Users/songweizhi/Desktop/ttt_mag/GI_refined_bins'
genome_ext          = 'fasta'
folder_out          = '/Users/songweizhi/Desktop/ttt_mag/GI_refined_bins_test'

wd                  = '/Users/songweizhi/Desktop'
wd                  = '.'

ref_to_strain_file  = '%s/ref_to_strain.txt'      % wd
blast_results       = '%s/bin_vs_ref.tab'         % wd
linkage_txt         = '%s/stats_linkage.txt'      % wd
bin_to_ref_txt      = '%s/stats_bin_to_ref.txt'   % wd
ref_to_bin_txt      = '%s/stats_ref_to_bin.txt'   % wd
iden_cutoff         = 100
aln_len_cutoff      = 1500
cov_q_cutoff        = 90
min_match_length    = 102400  # 100 Kbp
bin_ref_connector   = '__|__'
pwd_plot_sankey_R   = '/Users/songweizhi/PycharmProjects/MarkerMAG/MarkerMAG/get_sankey_plot.R'
pwd_plot_sankey_R   = '%s/get_sankey_plot.R' % wd


# get ref_to_strain_dict
ref_to_strain_dict = {}
for ref in open(ref_to_strain_file):
    ref_split = ref.strip().split('\t')
    ref_to_strain_dict[ref_split[0]] = ref_split[1]


file_re = '%s/*.%s' % (genome_folder, genome_ext)
file_list = [os.path.basename(file_name) for file_name in glob.glob(file_re)]
for genome in file_list:
    genome_no_ext          = '.'.join(genome.split('.')[:-1])
    pwd_genome             = '%s/%s' % (genome_folder, genome)
    pwd_genome_renamed_tmp = '%s/%s_renamed.%s' % (genome_folder, genome_no_ext, genome_ext)
    pwd_genome_renamed     = '%s/%s' % (folder_out, genome)
    #os.system('BioSAK rename_seq -in %s -prefix %s' % (pwd_genome, genome_no_ext))
    #os.system('mv %s %s' % (pwd_genome_renamed_tmp, pwd_genome_renamed))

linkage_dict = {}
for match in open(blast_results):
    match_split = match.strip().split('\t')
    query = match_split[0]
    query_genome = '_'.join(query.split('_')[:2])
    subject = match_split[1]
    subject_genome = '_'.join(subject.split('_')[:2])
    iden = float(match_split[2])
    aln_len = int(match_split[3])
    query_len = int(match_split[12])
    subject_len = int(match_split[13])
    coverage_q = float(aln_len) * 100 / float(query_len)
    if (iden >= iden_cutoff) and (aln_len > aln_len_cutoff) and (coverage_q >= cov_q_cutoff):
        key_bin_ref = '%s%s%s' % (query_genome, bin_ref_connector, subject_genome)
        if key_bin_ref not in linkage_dict:
            linkage_dict[key_bin_ref] = aln_len
        else:
            linkage_dict[key_bin_ref] += aln_len


linkage_dict_filtered = {}
for each_key in linkage_dict:
    if linkage_dict[each_key] >= min_match_length:
        linkage_dict_filtered[each_key] = linkage_dict[each_key]


# write out linkages
linkage_txt_handle = open(linkage_txt, 'w')
linkage_txt_handle.write('Bin,Ref,Length\n')
bin_to_ref_dict = {}
ref_to_bin_dict = {}
for each_linkage in linkage_dict_filtered:
    linkage_split = each_linkage.split(bin_ref_connector)
    bin_id = linkage_split[0]
    ref_id = linkage_split[1]

    if bin_id not in bin_to_ref_dict:
        bin_to_ref_dict[bin_id] = {ref_id}
    else:
        bin_to_ref_dict[bin_id].add(ref_id)

    if ref_id not in ref_to_bin_dict:
        ref_to_bin_dict[ref_id] = {bin_id}
    else:
        ref_to_bin_dict[ref_id].add(bin_id)

    linkage_txt_handle.write('%s,%s,%s\n' % (bin_id, ref_id, linkage_dict_filtered[each_linkage]))
linkage_txt_handle.close()

# visualize
cmd_sankey_paired = 'Rscript %s -f %s -x %s -y %s' % (pwd_plot_sankey_R, linkage_txt, 600, 1800)
os.system(cmd_sankey_paired)


bin_to_ref_txt_handle = open(bin_to_ref_txt, 'w')
for each_bin in bin_to_ref_dict:
    bin_to_ref_txt_handle.write('%s\t%s\n' % (each_bin, ','.join(bin_to_ref_dict[each_bin])))
bin_to_ref_txt_handle.close()

ref_to_bin_txt_handle = open(ref_to_bin_txt, 'w')
for each_ref in ref_to_bin_dict:
    ref_to_bin_txt_handle.write('%s\t%s\t%s\n' % (each_ref, ','.join(ref_to_bin_dict[each_ref]), ref_to_strain_dict[each_ref]))
ref_to_bin_txt_handle.close()



'''
Refined_19  40  Bordetella pertussis genomes 
Refined_25  3   Ruminiclostridium thermocellum
Refined_4   4   Clostridium butyricum
Refined_52  3   Clostridium botulinum
Refined_59  2   Clostridium botulinum B
Refined_7   4   Clostridium pasteurianum
Refined_29  2   [Clostridium] stercorarium subsp. stercorarium DSM 8532
Refined_60  2   Erysipelothrix rhusiopathiae
Refined_21	2   Desulfovibrio vulgaris

OTU_97.5710.0	Refined_82	CP002660.1 Clostridium acetobutylicum DSM 1731
OTU_97.5710.0	Refined_30	CP002660.1 Clostridium acetobutylicum DSM 1731
OTU_97.16157.1	Refined_30	AE001437.1 Clostridium acetobutylicum ATCC 824

OTU_97.2555.1	Refined_95	CP006763.1 Clostridium autoethanogenum DSM 10061
OTU_97.2555.1	Refined_91	CP006763.1 Clostridium autoethanogenum DSM 10061
OTU_97.2555.1	Refined_32	CP006763.1 Clostridium autoethanogenum DSM 10061
OTU_97.5374.1	Refined_32	CP012395.1 Clostridium autoethanogenum DSM 10061

OTU_97.2496.0	Refined_96,Refined_1	CP011319.1 Janthinobacterium sp. 1_2014MBL_MicDiv
OTU_97.5222.0	Refined_20,Refined_40	CP014028.1 Achromobacter xylosoxidans strain FDAARGOS_150
OTU_97.25906.0	Refined_27,Refined_61,Refined_88	CP016172.1 Bordetella flabilis strain AU10664
OTU_97.326.1	Refined_51,Refined_31	CP006721.1 Clostridium saccharobutylicum DSM 13864
OTU_97.577.0	Refined_39,Refined_43	CP009933.1 Clostridium scatologenes strain ATCC 25775
OTU_97.40665.0	Refined_83,Refined_5	AE015928.1 Bacteroides thetaiotaomicron VPI-5482
OTU_97.740.0	Refined_68,Refined_69	AP012331.1 Bifidobacterium scardovii JCM 12489 = DSM 13734 DNA
OTU_97.156.0	Refined_6,Refined_93	AM902716.1 Bordetella petrii strain DSM 12804
OTU_97.158.0	Refined_86,Refined_9	CP012334.1 Bordetella sp. H567

OTU_97.2472.0	Refined_11	CP002160.1 Clostridium cellulovorans 743B
OTU_97.680.0	Refined_12	CP009687.1 Clostridium aceticum strain DSM 1496
OTU_97.1070.0	Refined_13	CP003221.1 Desulfovibrio africanus str. Walvis Bay
OTU_97.7000.0	Refined_14	CP016440.1 Bordetella pseudohinzii strain HI4681
OTU_97.3427.0	Refined_16	CP012938.1 Bacteroides ovatus strain ATCC 8483
OTU_97.21930.0	Refined_17	LN877293.1 Bacteroides fragilis genome assembly BFBE1.1
OTU_97.4796.0	Refined_18	CP001998.1 Coraliomargarita akajimensis DSM 45221
OTU_97.19208.1	Refined_22	CP001726.1 Eggerthella lenta DSM 2243
OTU_97.20224.0	Refined_24	CP002530.1 Bacteroides salanitronis DSM 18170
OTU_97.3402.0	Refined_28	CP002770.1 Desulfotomaculum kuznetsovii DSM 6115
OTU_97.216.0	Refined_2	CP000468.1 Escherichia coli APEC O1
OTU_97.68.0		Refined_33	CP003040.1 Roseburia hominis A2-183
OTU_97.26691.0	Refined_36	AP012212.1 Clostridium sp. SY8519 DNA
OTU_97.427.0	Refined_37	CP015409.1 Akkermansia muciniphila strain YL44
OTU_97.3188.0	Refined_38	CP001684.1 Slackia heliotrinireducens DSM 20476
OTU_97.22911.0	Refined_3	CP013232.1 Collimonas fungivorans strain Ter6
OTU_97.9303.0	Refined_42	BA000016.3 Clostridium perfringens str. 13 DNA
OTU_97.463.0	Refined_53	CP009302.1 Coriobacteriaceae bacterium 68-1-3
OTU_97.29331.0	Refined_58	CP001850.2 Mageeibacillus indolicus UPII9-5
OTU_97.3311.0	Refined_62	CP017151.1 Lactobacillus fermentum strain NCC2970
OTU_97.2657.0	Refined_64	LN515532.1 Porphyromonadaceae bacterium ING2-E5B
OTU_97.740.1	Refined_71	CP002220.1 Bifidobacterium bifidum S17
OTU_97.40660.0	Refined_76	CP002122.1 Prevotella melaninogenica ATCC 25845
OTU_97.30976.0	Refined_78	CP017037.1 Dialister pneumosintes strain F0677
OTU_97.41827.1	Refined_8	CP002109.1 [Clostridium] saccharolyticum WM1
OTU_97.1154.0	Refined_97	HE965803.1 Bordetella parapertussis Bpp5

'''

