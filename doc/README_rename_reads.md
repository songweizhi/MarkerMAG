
Manual for the `rename_reads` module
---

### Notes

+ :warning: All reads in the forward and reverse reads files **must be in pair** and their orders in the two files **must be the same**.

+ The sizes of the output files might be slightly smaller than the input files, 
   which might due to the shorter length of read names in the output files.

   You can use `wc -l input_R1.fasta` and `wc -l output_R1.fasta` to get the number of lines in the input and output files and see if the numbers are the same.   
   You can also use `head R1.fasta` or `tail R1.fasta` to print on the screen the first/last few lines of a file.


### Example commands

+ In fasta format

      MarkerMAG rename_reads -r1 R1.fasta -r2 R2.fasta -p soil -t 2
      
+ In fastq format

      MarkerMAG rename_reads -r1 R1.fastq -r2 R2.fastq -p soil -fq -t 2


### Format of renamed reads

+ [prefix]_R1.fasta
    
      >[prefix]_1.1
      ATGCATGCATGCATGCATGC
      >[prefix]_2.1
      ATGCATGCATGCATGCATGC
      ...
    
+ [prefix]_R2.fasta
    
      >[prefix]_1.2
      ATGCATGCATGCATGCATGC
      >[prefix]_2.2
      ATGCATGCATGCATGCATGC
      ...
