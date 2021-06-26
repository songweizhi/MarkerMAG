
js_header = '''#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -l mem=30gb
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -M weizhi.song@unsw.edu.au
#PBS -m ae

module load python/3.7.3
source ~/mypython3env/bin/activate
module load java/8u121
module load fastqc/0.11.8
module load idba/1.1.3 
cd /srv/scratch/z5039045/MarkerMAG_wd/Tara_SOC/1_raw_reads
'''

js_folder = '/Users/songweizhi/Desktop/Tara_SOC_js_qc'


for each_sample in ['ERR598945', 'ERR599059', 'ERR3587138', 'ERR3587191', 'ERR599090', 'ERR599176', 'ERR3587139', 'ERR3587192', 'ERR599104', 'ERR599121', 'ERR3587140', 'ERR3587193']:

    pwd_js = '%s/js_%s.sh' % (js_folder, each_sample)

    cmd_gunzip_1 = 'gunzip %s_1.fastq.gz\n' % each_sample
    cmd_gunzip_2 = 'gunzip %s_2.fastq.gz\n' % each_sample

    r1_raw   = '%s_1.fastq'     % each_sample
    r2_raw   = '%s_2.fastq'     % each_sample
    r1_qc_p  = '%s_1_P.fastq'   % each_sample
    r1_qc_up = '%s_1_UP.fastq'  % each_sample
    r2_qc_p  = '%s_2_P.fastq'   % each_sample
    r2_qc_up = '%s_2_UP.fastq'  % each_sample

    cmd_qc_raw      = 'fastqc %s %s\n' % (r1_raw, r2_raw)
    cmd_trimmomatic = 'java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar PE %s %s %s %s %s %s ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/TruSeq3-PE-2.fa:2:30:10 CROP:99 HEADCROP:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:6:25 MINLEN:50\n' % (r1_raw, r2_raw, r1_qc_p, r1_qc_up, r2_qc_p, r2_qc_up)
    cmd_qc_p        = 'fastqc %s %s\n' % (r1_qc_p, r2_qc_p)
    cmd_rename_fq   = 'MarkerMAG rename_reads -r1 %s -r2 %s -p %s -fq -t 2\n' % (r1_qc_p, r2_qc_p, each_sample)
    cmd_fq2fa_1     = 'fq2fa %s_R1.fastq %s_R1.fasta\n' % (each_sample, each_sample)
    cmd_fq2fa_2     = 'fq2fa %s_R2.fastq %s_R2.fasta\n' % (each_sample, each_sample)
    cmd_gzip_1      = 'gzip %s_1.fastq\n' % each_sample
    cmd_gzip_2      = 'gzip %s_2.fastq\n' % each_sample

    pwd_js_handle = open(pwd_js, 'w')
    pwd_js_handle.write(js_header)
    pwd_js_handle.write(cmd_gunzip_1)
    pwd_js_handle.write(cmd_gunzip_2)
    pwd_js_handle.write(cmd_qc_raw)
    pwd_js_handle.write(cmd_trimmomatic)
    pwd_js_handle.write(cmd_qc_p)
    pwd_js_handle.write(cmd_rename_fq)
    pwd_js_handle.write(cmd_fq2fa_1)
    pwd_js_handle.write(cmd_fq2fa_2)
    pwd_js_handle.write(cmd_gzip_1)
    pwd_js_handle.write(cmd_gzip_2)
    pwd_js_handle.close()
