

### Notes

1. All reads in the R1.fastq and R2.fastq must be in pair and their orders in the two files must be the same :warning:.

1. The sizes of the output files might be slightly smaller than the input files, 
which might due to the shorter length of read names in the renamed files.

### Example commands

      # in fasta format
      MarkerMAG rename_reads -r1 R1.fasta -r2 R2.fasta -p soil -t 2
      
      # in fastq format
      MarkerMAG rename_reads -r1 R1.fastq -r2 R2.fastq -p soil -fq -t 2


### Format of renamed reads

    [prefix]_R1.fasta
    
      >[prefix]_1.1
      ATGCATGCATGCATGCATGC
      >[prefix]_2.1
      ATGCATGCATGCATGCATGC
      ...
    
    [prefix]_R2.fasta
    
      >[prefix]_1.2
      ATGCATGCATGCATGCATGC
      >[prefix]_2.2
      ATGCATGCATGCATGCATGC
      ...

