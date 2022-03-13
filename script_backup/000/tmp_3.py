

for each_gnm in open('/Users/songweizhi/Desktop/gnm_id.txt'):

    gnm_id = each_gnm.strip()
    mapping_cmd = 'bowtie2 -x %s -U ../../MBARC26_R1.fasta,../../MBARC26_R2.fasta -S %s_report_best_default.sam -p 12 -f --xeq --local --no-unal' % (gnm_id, gnm_id)
    #mapping_cmd = 'bowtie2 -x %s -U ../MBARC26_R1.fasta,../MBARC26_R2.fasta -S %s_report_best.sam -p 12 -f --xeq --local --no-unal -N 1 -L 30' % (gnm_id, gnm_id)
    filter_cmd = 'python3 filter_sam.py -in %s_report_best_default.sam -out %s_report_best_default_mis5.sam -mm 5 -aln 50' % (gnm_id, gnm_id)
    print(mapping_cmd)
    #print(filter_cmd)

    #print(mapping_cmd)
    #print(filter_cmd)
    #concate_cmd = '%s; %s' % (mapping_cmd, filter_cmd)
    #print(concate_cmd)
    sort_cmd = 'samtools sort -O sam --threads 12 -o %s_report_best_mis2_sorted.sam %s_report_best_mis2.sam' % (gnm_id, gnm_id)
    #print(sort_cmd)
    get_depth_cmd = 'samtools depth %s_report_best_mis2_sorted.sam > %s_report_best_mis2_sorted_depth.txt' % (gnm_id, gnm_id)
    #print(get_depth_cmd)

''':key


python3 filter_sam.py -in DA_report_best.sam -out DA_report_best_mis3.sam -mm 3 -aln 50
python3 filter_sam.py -in DA_report_best.sam -out DA_report_best_mis4.sam -mm 4 -aln 50
python3 filter_sam.py -in DA_report_best.sam -out DA_report_best_mis5.sam -mm 5 -aln 50


python3 filter_sam.py -in FP_report_best.sam -out FP_report_best_mis5.sam -mm 5 -aln 50




bowtie2 -x FP -U ../../MBARC26_R1.fasta,../../MBARC26_R2.fasta -S FP_report_all_default.sam -p 12 -f --all --xeq --local --no-unal
python3 filter_sam.py -in FP_report_all_default.sam -out FP_report_all_default_mis5.sam -mm 5 -aln 50

'''