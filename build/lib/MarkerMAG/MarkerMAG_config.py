import os

pwd_config_file = os.path.realpath(__file__)
config_file_path = '/'.join(pwd_config_file.split('/')[:-1])
config_dict = {'config_file_path': config_file_path}
