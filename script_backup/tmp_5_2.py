import os

js_header = '''#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=90gb
#PBS -l walltime=71:59:00
#PBS -j oe
#PBS -M weizhi.song@unsw.edu.au
#PBS -m ae

module load python/3.6.5
source ~/mypython365env/bin/activate
module load java/7u51
module load gcc/7.3.0
module load sparsehash/2.0.3
module load matam/1.5.3
module load samtools/1.9
module load usearch/11.0.667
module load seqtk/20190219
export PATH=/home/z5039045/anaconda3/bin:$PATH

cd /srv/scratch/z5039045/Kelp_16S
'''

Kelp_sample_id = '/Users/songweizhi/Desktop/Kelp_sample_id.txt'
js_folder = '/Users/songweizhi/Desktop/js_matam_16s'


co_assembly_dict = {}
for co_assembly in open(Kelp_sample_id):
    co_assembly_split = co_assembly.strip().split('\t')
    if co_assembly_split[0] not in co_assembly_dict:
        co_assembly_dict[co_assembly_split[0]] = [co_assembly_split[1]]
    else:
        co_assembly_dict[co_assembly_split[0]].append(co_assembly_split[1])


for each_co_assembly in co_assembly_dict:
    #print(each_co_assembly)

    js_file = '%s/js_matam_%s.sh' % (js_folder, each_co_assembly)
    js_file_handle = open(js_file, 'w')
    js_file_handle.write(js_header)
    js_file_handle.write('MarkerMAG matam_16s -p %s -r1 %s_R1.fasta -r2 %s_R2.fasta -pct 0.1,0.5,1,5,10,25,50,75 -ref /srv/scratch/z5039045/DB/Matam/SILVA_128_SSURef_NR95 -i 0.995 -t 16 -force -matam_assembly /home/z5039045/anaconda3/pkgs/matam-v1.5.3-0/bin/matam_assembly.py -sortmerna /home/z5039045/anaconda3/pkgs/matam-v1.5.3-0/opt/matam-v1.5.3/sortmerna/sortmerna\n' % (each_co_assembly, each_co_assembly, each_co_assembly))
    js_file_handle.close()

