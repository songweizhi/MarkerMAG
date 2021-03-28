
## MarkerMAG (linking MAGs with 16S rRNA marker genes)

![logo](images/MarkerMAG_logo.jpg) 

[![pypi licence](https://img.shields.io/pypi/l/MarkerMAG.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version](https://img.shields.io/pypi/v/MarkerMAG.svg)](https://pypi.python.org/pypi/MarkerMAG) 


Publication
---
+ In preparation
+ **Weizhi Song** (songwz03[at]gmail.com)
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia


Dependencies
---
 
 Dependencies are module-specific, please see details below:
 
+ `link`: 
  [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/),
  [Mira v5rc2](https://github.com/bachev/mira) (default) or [SPAdes](https://github.com/ablab/spades), 
  [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download),
  [optparse](https://cran.r-project.org/web/packages/optparse/index.html) (R) and 
  [googleVis](https://cran.r-project.org/web/packages/googleVis/index.html) (R)

+ `matam_16s`: 
  [SortMeRNA](https://github.com/biocore/sortmerna), 
  [MATAM](https://github.com/bonsai-team/matam),
  [Usearch](https://www.drive5.com/usearch/) and 
  [seqtk](https://github.com/lh3/seqtk)

+ `uclust_16S`: 
  [Usearch](https://www.drive5.com/usearch/)

+ `subsample_reads`: 
  [Usearch](https://www.drive5.com/usearch/)


How to install:
---

MarkerMAG is implemented in python3, you can install it with pip3:

    # install with 
    pip3 install MarkerMAG
        
    # upgrade with 
    pip3 install --upgrade MarkerMAG


MarkerMAG modules:
---

+ Main module

      link             ->  link MAGs with marker genes
    
+ Supplementary modules

      rename_reads     ->  rename paired reads 
      matam_16s        ->  Assemble 16S rRNA genes with Matam, including subsample and dereplication
      uclust_16s       ->  cluster marker genes with Usearch
      barrnap_16s      ->  identify 16S gene sequences with Barrnap
      subsample_reads  ->  subsample reads with Usearch


Important Notes :warning:
---

1. MarkerMAG assumes the id of paired reads in the format of `XXXX.1` and `XXXX.2`. The only difference is the last character.
   You can rename your reads with MarkerMAG's `rename_reads` module. 
   
    :warning: All reads in the R1.fastq and R2.fastq must be in pair and their orders in the two files must be the same.

       MarkerMAG rename_reads -r1 R1.fastq -r2 R2.fastq -p soil -fq -t 2
        
       # output files and format of renamed reads id:
       # Soil_R1.fasta: soil_1.1, soil_2.1, soil_3.1 ...
       # Soil_R2.fasta: soil_1.2, soil_2.2, soil_3.2 ...

1. The reconstruction of 16S rRNA genes by Matam is highly affected by sequencing depth ([ref](to/be/added)), we thus recommend to 
   run Matam on reads subsets subsampled at different percentage and combine assemblies at all depth, followed by dereplication.

   The following command extracts 16S rRNA reads from `combined_paired_reads.fasta` and subsample at percentage of `1, 5, 10, 25, 50, 75 and 100`.
   16S rRNA genes reconstructed from all subsets will be combined and clustered at identity cut-off of `99.9%`.
   The longest sequence from each cluster will be kept.
   
   Please refer to [here](demo_files/README_Matam.md) for running Matam using an updated SILVA SSU database (version 138.1).
        
       # convert fastq files fasta files (e.g. with idba's fq2fa)
       fq2fa R1.fastq R1.fasta
       fq2fa R2.fastq R2.fasta
       MarkerMAG matam_16s -p Soil -r1 R1.fasta -r2 R2.fasta -pct 1,5,10,25,50,75,100 -i 0.999 -ref /srv/scratch/z5039045/DB/SILVA/SILVA_138_1_SSURef_NR99_id99/SILVA_138.1_SSURef_NR99_tax_silva_NR99 -t 12 -force

2. :warning: All MAGs derived from a set should be included in MarkerMAG run. (more details need to be added)


How to run:
---

+ Link 16S rRNA gene sequences with MAGs: 

      MarkerMAG link -p Soil -r1 R1.fastq -r2 R2.fastq -marker Soil_16S_uclust_0.999.fasta -mag refined_MAG -x fasta -t 12 -tmp -force


Output files:
---

1. Linkage table

    | MarkerGene | Genome | Linkages | Step |
    |:---:|:---:|:---:|:---:|
    | Soil_subsample_10_56   | Refined_MAG_26| 65| S1 |
    | Soil_subsample_25_66   | Refined_MAG_19| 23| S1 |
    | Soil_subsample_100_322 | Refined_MAG_47| 10| S1 |
    | Soil_subsample_75_284  | Refined_MAG_5 | 5 | S1 |
    | Soil_subsample_100_563 | Refined_MAG_26| 3 | S1 |
    | Soil_subsample_100_133 | Refined_MAG_42| 30| S2 |
    | Soil_subsample_100_133 | Refined_MAG_42| 11| S2 |
    | Soil_subsample_100_262 | Refined_MAG_31| 7 | S2 |

1. Visualization of linkages
![linkages](images/linkage.png) 
