

for each_sample in ['BH_ER_050417', 'BH_ER_110816', 'CB_ER_080217']:
    for each_rep in [1, 2, 3]:

        samp_with_rep = '%s_%s' % (each_sample, each_rep)
        trim_cmd = 'java -jar /apps/trimmomatic/0.38/trimmomatic-0.38.jar PE %s_pairedForward.fastq %s_pairedReverse.fastq %s_R1_P.fastq %s_R1_UP.fastq %s_R2_P.fastq %s_R2_UP.fastq ILLUMINACLIP:/apps/trimmomatic/0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:25 TRAILING:25 CROP:225 HEADCROP:5 SLIDINGWINDOW:6:20 MINLEN:50' % (samp_with_rep, samp_with_rep, samp_with_rep, samp_with_rep, samp_with_rep, samp_with_rep)
        #fastqc_cmd = 'fastqc %s_R1_P.fastq %s_R2_P.fastq' % (samp_with_rep, samp_with_rep)
        #print(trim_cmd)
        #print(fastqc_cmd)

        #fq2fa_cmd = 'fq2fa %s_R1_P.fastq %s_R1_P.fasta\nfq2fa %s_R2_P.fastq %s_R2_P.fasta\n' % (samp_with_rep, samp_with_rep, samp_with_rep, samp_with_rep)
        #print(fq2fa_cmd)

        # spades_cmd = 'spades.py --meta -1 %s_%s_pairedForward.fastq -2 %s_%s_pairedReverse.fastq -o %s_%s_spades_wd -t 12 -k 55,75,99,127' % (each_sample, each_rep, each_sample, each_rep, each_sample, each_rep)
        # print(spades_cmd)

        #spades_cmd = 'spades.py --only-assembler --meta -1 %s_%s_R1_P.fasta -2 %s_%s_R2_P.fasta -o %s_%s_spades_wd -t 12 -k 55,75,99,127' % (each_sample, each_rep, each_sample, each_rep, each_sample, each_rep)
        #print(spades_cmd)

        #prokka_cmd = 'prokka --force --compliant --cpus 12 --kingdom Bacteria --prefix %s_%s --locustag %s_%s --outdir %s_%s_prokka_wd %s_%s_scaffolds.fasta' % (each_sample, each_rep, each_sample, each_rep, each_sample, each_rep, each_sample, each_rep)
        #print(prokka_cmd)

        diamond_cmd = 'diamond blastp -q %s_%s.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out %s_%s_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10' % (each_sample, each_rep, each_sample, each_rep)
        #print(diamond_cmd)

        best_hit_cmd = 'BioSAK BestHit -i %s_%s_vs_Tax4Fun2_KEGG.tab -o %s_%s_vs_Tax4Fun2_KEGG_best_hit.tab' % (each_sample, each_rep, each_sample, each_rep)
        print(best_hit_cmd)
        #
        # print()

# spades.py --meta -1 Kelp_R1.fastq -2 Kelp_R2.fastq -o BH_ER_050417_SPAdes_wd -t 12 -k 21,33,55,75,99,127
# diamond blastp -q BH_ER_050417.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_050417_vs_Tax4Fun2_KEGG.tab --outfmt 6 --block-size 1 --threads 16 -e 1e-10
# BioSAK BestHit -i BH_ER_050417_vs_Tax4Fun2_KEGG.tab -o BH_ER_050417_vs_Tax4Fun2_KEGG_best_hit.tab
