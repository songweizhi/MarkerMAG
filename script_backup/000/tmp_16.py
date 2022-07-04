import os

cms_1 = 'spades.py --only-assembler -s Demo_MarkerMAG_wd/Demo_rd2_wd/rd2_read_to_extract_flanking_both_R12_up.fa -o Demo_MarkerMAG_wd/Demo_rd2_wd/mini_assembly_SPAdes_wd -t 6 -k 75,99,123'
cms_2 = 'spades.py --only-assembler -s Demo_MarkerMAG_wd/Demo_rd2_wd/rd2_read_to_extract_flanking_both_R12_up.fa -o Demo_MarkerMAG_wd/Demo_rd2_wd/mini_assembly_SPAdes_wd -t 6 -k 75,99,123 -m 1024'

try:
    os.system(cms_1)
except:
    os.system(cms_2)


