
# annotate metagenome (Tax4Fun2 cutoff: '-e 1e-10')
module load python/3.7.3
source ~/mypython3env/bin/activate
module load diamond/0.9.24
cd /srv/scratch/z5039045/MarkerMAG_wd/Tax4Fun2_wd/BH_ER_050417_prokka_wd
diamond blastp -q BH_ER_050417.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_050417_vs_Tax4Fun2_KEGG.tab --outfmt 6 --block-size 1 --threads 16 -e 1e-10
BioSAK BestHit -i BH_ER_050417_vs_Tax4Fun2_KEGG.tab -o BH_ER_050417_vs_Tax4Fun2_KEGG_best_hit.tab


module load python/3.7.3
source ~/mypython3env/bin/activate
module load diamond/0.9.24
cd /srv/scratch/z5039045/MarkerMAG_wd/Kelp_Tax4Fun2/spades_wd2

diamond blastp -q BH_ER_050417_1_prokka_wd/BH_ER_050417_1.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_050417_1_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i BH_ER_050417_1_vs_Tax4Fun2_KEGG.tab -o BH_ER_050417_1_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q BH_ER_050417_2_prokka_wd/BH_ER_050417_2.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_050417_2_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i BH_ER_050417_2_vs_Tax4Fun2_KEGG.tab -o BH_ER_050417_2_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q BH_ER_050417_3_prokka_wd/BH_ER_050417_3.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_050417_3_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i BH_ER_050417_3_vs_Tax4Fun2_KEGG.tab -o BH_ER_050417_3_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q BH_ER_110816_1_prokka_wd/BH_ER_110816_1.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_110816_1_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i BH_ER_110816_1_vs_Tax4Fun2_KEGG.tab -o BH_ER_110816_1_vs_Tax4Fun2_KEGG_best_hit.tab

#diamond blastp -q BH_ER_110816_2_prokka_wd/BH_ER_110816_2.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_110816_2_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
#BioSAK BestHit -i BH_ER_110816_2_vs_Tax4Fun2_KEGG.tab -o BH_ER_110816_2_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q BH_ER_110816_3_prokka_wd/BH_ER_110816_3.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out BH_ER_110816_3_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i BH_ER_110816_3_vs_Tax4Fun2_KEGG.tab -o BH_ER_110816_3_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q CB_ER_080217_1_prokka_wd/CB_ER_080217_1.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out CB_ER_080217_1_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i CB_ER_080217_1_vs_Tax4Fun2_KEGG.tab -o CB_ER_080217_1_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q CB_ER_080217_2_prokka_wd/CB_ER_080217_2.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out CB_ER_080217_2_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i CB_ER_080217_2_vs_Tax4Fun2_KEGG.tab -o CB_ER_080217_2_vs_Tax4Fun2_KEGG_best_hit.tab

diamond blastp -q CB_ER_080217_3_prokka_wd/CB_ER_080217_3.faa --db /srv/scratch/z5039045/DB/KEGG_Tax4Fun2/prokaryotes.dmnd --out CB_ER_080217_3_vs_Tax4Fun2_KEGG.tab --outfmt 6 --threads 16 -e 1e-10
BioSAK BestHit -i CB_ER_080217_3_vs_Tax4Fun2_KEGG.tab -o CB_ER_080217_3_vs_Tax4Fun2_KEGG_best_hit.tab



################ get copy number of 16S ################

# MAG with 16S but not linked
Refined_23	1.33  1
Refined_19	0.75  1
Refined_50	1.13  1
Refined_5	  1.11  1

# MAG with 16S and got linked
Refined_57	1     1
Refined_34	1.02  1
Refined_20	1.18  1
Refined_37	1.88  2

# MAG without 16S but got linked
Refined_4	  4.26  4
Refined_32	2.77  3
Refined_17	3.5   4
Refined_47	1.05  1
Refined_56	1.66  2
Refined_27	5.12  5
Refined_26	0.81  1
Refined_48	0.21  1
Refined_3	  2.85  3

Refined_4	  Kelp_SILVA138_id99_subsample_100_784	117.86(1.07)	115.09(1.08)	117.86	4.28	4.26
Refined_32	Kelp_SILVA138_id99_subsample_75_625	124.87(1.04)	124.87(1.04)	124.87	2.82	2.77
Refined_17	Kelp_SILVA138_id99_subsample_75_664	134.17(1.11)	133.8(1.12)	134.17	3.48	3.5
Refined_47	Kelp_SILVA138_id99_subsample_50_390	160.28(1.32)	159.17(1.33)	119.6766917293233	1.05	1.05
Refined_56	Kelp_SILVA138_id99_subsample_25_300	272.1(0.83)	182.26(0.74)	134.8724	1.67	1.66
Refined_27	Kelp_SILVA138_id99_subsample_75_277	698.33(0.93)	699.79(0.93)	699.79	5.12	5.12
Refined_26	Kelp_SILVA138_id99_subsample_100_789	260.62(1.01)	248.5(1.02)	260.62	0.82	0.81
Refined_48	Kelp_SILVA138_id99_subsample_100_531	12.06(0.44)	8.28(0.47)	3.8915999999999995	0.21	0.21
Refined_3	  Kelp_SILVA138_id99_subsample_100_644	124.34(0.99)	124.18(0.99)	124.34	2.86	2.85



########################################################












