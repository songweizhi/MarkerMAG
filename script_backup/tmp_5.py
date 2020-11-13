import os


'''
# subsample reads 
module load usearch/11.0.667
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26
MarkerMAG subsample_reads -r1 MBARC26_R1.fasta -r2 MBARC26_R2.fasta -ratio 0.05 -usearch /srv/scratch/z5039045/Softwares/usearch/usearch11.0.667_i86linux64

# map reads to reference
module load samtools/1.10
module load bowtie/2.3.5.1
cd /srv/scratch/z5039045/MarkerMAG_wd/new_algorithm
bowtie2 -x MBARC26_fake_bins_combined -1 MBARC26_R1_0.05.fasta -2 MBARC26_R2_0.05.fasta -S MBARC26.sam -p 12 -f --quiet
samtools view -bS MBARC26.sam -o MBARC26.bam
samtools sort MBARC26.bam -o MBARC26_sorted.bam
samtools index MBARC26_sorted.bam




module load samtools/1.10
cd /srv/scratch/z5039045/MarkerMAG_wd/MBARC26/MBARC26_link_Matam_FakeBins_500bp_depth_d0.2_MarkerMAG_wd/MBARC26_fake_bins_db
bowtie2 -x %s/%s -1 %s -2 %s -S %s -f --local --no-unal --quiet --threads %s
samtools view -f 8 -F 4 -b foo.bam > foo.filtered.bam

-f 8 -F 4   :   yield reads that are mapped with unmapped mates.
-F 8 -f 4   :   yield unmapped reads with mapped mates, but flagstat doesn't count mapped mates unless the read being looked at is also mapped (thus, the singleton metric with -f 8 -F 4).

-f INT   only include reads with all  of the FLAGs in INT present [0]
-F INT   only include reads with none of the FLAGS in INT present [0]
4   segment unmapped
8   next segment in the template unmapped

-f 4
-f 8
-F 4
-F 8

samtools view -f 8 -F 4 -b MBARC26_fake_bins_combined_sorted.bam > MBARC26_fake_bins_combined_sorted.f8.F4.bam
samtools view -F 8 -f 4 -b MBARC26_fake_bins_combined_sorted.bam > MBARC26_fake_bins_combined_sorted.F8.f4.bam



cd /srv/scratch/z5039045/MarkerMAG_wd/new_algorithm
samtools view -f 8 -F 4 -b MBARC26_sorted.bam > MBARC26_sorted.f8.F4.bam
samtools view -F 8 -f 4 -b MBARC26_sorted.bam > MBARC26_sorted.F8.f4.bam

samtools view -h -o MBARC26_sorted.f8.F4.sam MBARC26_sorted.f8.F4.bam
samtools view -h -o MBARC26_sorted.F8.f4.sam MBARC26_sorted.F8.f4.bam

'''





