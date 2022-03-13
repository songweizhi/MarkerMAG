import os

pwd_config_file = os.path.realpath(__file__)
config_file_path = '/'.join(pwd_config_file.split('/')[:-1])

config_dict = {'config_file_path'   : config_file_path,
               'silva_order_refs'   : '%s/SILVA_16S_order.fasta'        % config_file_path,
               'vxtractor'          : '%s/vxtractor/vxtractor.pl'       % config_file_path,
               'hmm_bac'            : '%s/vxtractor/HMMs/SSU_bacteria'  % config_file_path,
               'hmm_arc'            : '%s/vxtractor/HMMs/SSU_archaea'   % config_file_path}
