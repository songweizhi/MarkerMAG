
js_header = '''#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=30gb
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -M weizhi.song@unsw.edu.au
#PBS -m ae

module load samtools/1.10
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/Tara_SOC/2_mapping
'''

js_folder = '/Users/songweizhi/Desktop/Tara_SOC_js_mapping'


for each_sample in ['ERR598945', 'ERR599059', 'ERR3587138', 'ERR3587191', 'ERR599090', 'ERR599176', 'ERR3587139', 'ERR3587192', 'ERR599104', 'ERR599121', 'ERR3587140', 'ERR3587193']:

    pwd_js = '%s/js_%s.sh' % (js_folder, each_sample)

    r1_qc_p  = '%s_R1.fasta'   % each_sample
    r2_qc_p  = '%s_R2.fasta'   % each_sample
    cmd_mapping = 'bowtie2 -x TARA_SOC_RAW_min2500 -1 ../1_raw_reads/%s -2 ../1_raw_reads/%s -S %s.sam -f -p 12\n' % (r1_qc_p, r2_qc_p, each_sample)
    cmd_sam2bam = 'samtools view -bS %s.sam -o %s.bam\n' % (each_sample, each_sample)
    cmd_sort    = 'samtools sort %s.bam -o %s_sorted.bam\n' % (each_sample, each_sample)
    cmd_index   = 'samtools index %s_sorted.bam\n' % each_sample
    cmd_rm_sam  = 'rm %s.sam\n' % each_sample
    cmd_rm_bam  = 'rm %s.bam\n' % each_sample

    pwd_js_handle = open(pwd_js, 'w')
    pwd_js_handle.write(js_header)
    pwd_js_handle.write(cmd_mapping)
    pwd_js_handle.write(cmd_sam2bam)
    pwd_js_handle.write(cmd_sort)
    pwd_js_handle.write(cmd_index)
    pwd_js_handle.write(cmd_rm_sam)
    pwd_js_handle.write(cmd_rm_bam)
    pwd_js_handle.close()
