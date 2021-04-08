
Manual for the `matam_16s` module
---

The reconstruction of 16S rRNA genes by Matam is highly affected by sequencing depth ([ref](to/be/added)), we thus recommend to 
   run Matam on reads subsets subsampled at different percentage and combine assemblies at all depth, followed by dereplication.


The following command extracts 16S rRNA reads from `combined_paired_reads.fasta` and subsample at percentage of `1, 5, 10, 25, 50, 75 and 100`.
16S rRNA genes reconstructed from all subsets will be combined and clustered at identity cut-off of `99.9%`.
The longest sequence from each cluster will be kept.

Please refer to [here](demo_files/README_Matam.md) for running Matam using the latest SILVA SSU database (v138.1).
    
   # convert fastq files fasta files (e.g. with idba's fq2fa)
   fq2fa R1.fastq R1.fasta
   fq2fa R2.fastq R2.fasta
   MarkerMAG matam_16s -p Soil -r1 R1.fasta -r2 R2.fasta -pct 1,5,10,25,50,75,100 -i 0.999 -ref /srv/scratch/z5039045/DB/SILVA/SILVA_138_1_SSURef_NR99_id99/SILVA_138.1_SSURef_NR99_tax_silva_NR99 -t 12 -force
















Prepare Matam database with the latest SILVA SSU sequences (v138.1)
---

1. Download SILVA SSU sequences (v138.1)

       # specify a location where you want to store the db files
       matam_db_folder='/srv/scratch/z5039045/DB/Matam'

       # download the SILVA SSU sequence file to the specified folder and decompress it
       cd $matam_db_folder
       wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/README.txt
       wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
       gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz

1. Format SILVA SSU sequences with Matam

       matam_db_preprocessing.py --clustering_id_threshold 0.99 --max_memory 30000 --cpu 12 -v -i SILVA_138.1_SSURef_NR99_tax_silva.fasta -d SILVA_138_1_SSURef_NR99_id99

1. Run Matam with generated SILVA SSU database

       matam_assembly.py -i filtered_reads_R1_R2.fasta -o Matam_outputs -d $matam_db_folder/SILVA_138_1_SSURef_NR99_id99/SILVA_138.1_SSURef_NR99_tax_silva_NR99 -v --cpu 12 --max_memory 30000 

For UNSW Katana users
---

    # Module needed 
    module load python/3.6.5
    module load java/7u51
    module load gcc/8.4.0
    module load sparsehash/2.0.3
    module load matam/1.5.3
    module load samtools/1.9
    matam_assembly.py -h
