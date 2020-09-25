import os

# extract path to the config file
pwd_config_file = os.path.realpath(__file__)
config_file_path = '/'.join(pwd_config_file.split('/')[:-1])

# specify full path to corresponding executables at the right side of colon
config_dict = {'config_file_path'       : config_file_path,
               'prodigal'               : 'prodigal',
               'hmmsearch'              : 'hmmsearch',
               'hmmfetch'               : 'hmmfetch',
               'hmmalign'               : 'hmmalign',
               'hmmstat'                : 'hmmstat',
               'mafft'                  : 'mafft',
               'bowtie2'                : 'bowtie2',
               'bowtie2_build'          : 'bowtie2-build',
               'blastp'                 : 'blastp',
               'blastn'                 : 'blastn',
               'makeblastdb'            : 'makeblastdb',
               'fasttree'               : 'FastTree',
               'ranger_mac'             : '%s/Ranger-DTL-Dated.mac'           % config_file_path,
               'ranger_linux'           : '%s/Ranger-DTL-Dated.linux'         % config_file_path,
               'path_to_hmm'            : '%s/MetaCHIP_phylo.hmm'             % config_file_path,
               'circos_HGT_R'           : '%s/MetaCHIP_circos_HGT.R'          % config_file_path,
               'label_tree_R'           : '%s/label_tree.R'                   % config_file_path,
               'cdd2cog_perl'           : '%s/cdd2cog.pl'                     % config_file_path,
               'get_sankey_plot_R'      : '%s/get_sankey_plot.R'              % config_file_path,
               'ko00001_keg'            : '%s/ko00001.keg'                    % config_file_path,
               'MetaCyc_rxns_with_ec'   : '%s/MetaCyc_reactions_with_ec.txt'  % config_file_path
               }

