
Prepare SILVA SSU database for Matam
---

1. Download SILVA SSU sequences (version 138.1)

       wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/README.txt
       wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
       gunzip SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz  

1. Format SILVA SSU sequences with Matam

       matam_db_preprocessing.py --clustering_id_threshold 0.99 --max_memory 30000 --cpu 6 -v -i SILVA_138.1_SSURef_NR99_tax_silva.fasta -d SILVA_138_1_SSURef_NR99_id99

1. Run Matam with generated SILVA SSU database

       matam_assembly.py -i filtered_reads.fasta -o Matam_outputs -d SILVA_138_1_SSURef_NR99_id99/SILVA_138.1_SSURef_NR99_tax_silva_NR99 --filter_only -v --cpu 6 --max_memory 30000 

For UNSW Katana users
---

    # Module needed 
    module load python/3.6.5
    module load java/7u51
    module load gcc/8.4.0
    module load sparsehash/2.0.3
    module load matam/1.5.3
    module load samtools/1.9

