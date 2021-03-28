
js_header = '''#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=30gb
#PBS -l walltime=11:59:00
#PBS -j oe
#PBS -M weizhi.song@unsw.edu.au
#PBS -m ae

cd /srv/scratch/z5039045/MarkerMAG_wd/Tara_SOC/1_raw_reads
'''

js_folder = '/Users/songweizhi/Desktop/Tara_SOC_js_wget'


for each_sample in ['wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR598/ERR598945/ERR598945_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR598/ERR598945/ERR598945_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599059/ERR599059_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599059/ERR599059_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/008/ERR3587138/ERR3587138_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/008/ERR3587138/ERR3587138_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/001/ERR3587191/ERR3587191_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/001/ERR3587191/ERR3587191_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599090/ERR599090_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599090/ERR599090_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599176/ERR599176_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599176/ERR599176_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/009/ERR3587139/ERR3587139_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/009/ERR3587139/ERR3587139_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/002/ERR3587192/ERR3587192_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/002/ERR3587192/ERR3587192_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599104/ERR599104_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599104/ERR599104_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599121/ERR599121_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR599/ERR599121/ERR599121_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/000/ERR3587140/ERR3587140_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/000/ERR3587140/ERR3587140_2.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/003/ERR3587193/ERR3587193_1.fastq.gz', 'wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/003/ERR3587193/ERR3587193_2.fastq.gz']:

    sample_id = each_sample.split('/')[-1].split('.')[0]
    pwd_js = '%s/js_%s.sh' % (js_folder, sample_id)
    pwd_js_handle = open(pwd_js, 'w')
    pwd_js_handle.write(js_header)
    pwd_js_handle.write(each_sample)
    pwd_js_handle.close()
