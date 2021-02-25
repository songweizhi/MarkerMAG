
## MarkerMAG (link MAGs with marker genes)

![logo](images/MarkerMAG_logo.jpg) 

[![pypi licence](https://img.shields.io/pypi/l/MarkerMAG.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version](https://img.shields.io/pypi/v/MarkerMAG.svg)](https://pypi.python.org/pypi/MarkerMAG) 


Contact
---

+ **Weizhi Song**, Postdoctoral Researcher
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia
+ E-mail: songwz03@gmail.com


Dependencies
---
 
+ `link`: 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
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

MarkerMAG can be installed via `pip3`:

    # First-time installation
    pip3 install MarkerMAG
        
    # upgrade
    pip3 install --upgrade MarkerMAG


Important Notes :warning:
---

1. MarkerMAG assumes the id of paired reads in the format of `XXXX.1` and `XXXX.2`. The only difference is the last character.
   You can rename your reads with MarkerMAG's `rename_reads` module. 
   
   :warning: Please note that all reads in R1.fastq and R2.fastq must be in pair and their orders in the two files must be the same.

       MarkerMAG rename_reads -r1 R1.fastq -r2 R2.fastq -p soil -fq -t 2
        
       # output files and format of renamed reads id:
       # Soil_R1.fasta: soil_1.1, soil_2.1, soil_3.1 ...
       # Soil_R2.fasta: soil_1.2, soil_2.2, soil_3.2 ...

1. The reconstruction of 16S rRNA genes by Matam is highly affected by sequencing depth ([ref](to/be/added)), we thus recommend to 
   run Matam on reads subsets subsampled at different percentage and combine assemblies at all depth, followed by dereplication.

   The following command extracts 16S rRNA reads from `combined_paired_reads.fasta` and subsample at percentage of `1, 5, 10, 25, 50 and 75`.
   16S rRNA genes reconstructed from all subsets will be combined and clustered at identity cut-off of `99.5%`.
   The longest sequence from each cluster will be kept.  
        
       # convert fastq files fasta files with idba
       fq2fa R1.fastq R1.fasta
       fq2fa R2.fastq R2.fasta
       MarkerMAG matam_16s -p Soil -r1 R1.fasta -r2 R2.fasta -pct 1,5,10,25,50,75,100 -i 0.995 -t 12 -force -ref /srv/scratch/z5039045/DB/Matam/SILVA_128_SSURef_NR95

How to run:
---

+ Link 16S rRNA gene sequences with MAGs: 

      MarkerMAG link -p Soil -r1 R1.fastq -r2 R2.fastq -marker Soil_16S_uclust_0.995.fasta -mag refined_MAG -x fasta -t 12 -tmp -force


Output files:
---

1. Linkage table

    | Marker | Genome | Paired | Clipping |
    |:---:|:---:|:---:|:---:|
    | g4_00414 | bin_g4 | 196 | 139 |
    | o3_02626 | bin_o3 | 100 | 81 |
    | s4_04216 | bin_s4 | 97 | 39 |
    | s4_00580 | bin_s4 | 84 | 41 |
    | o2_01394 | bin_o2 | 58 | 0 |


1. Visualization of linkages
![linkages](images/linkage.png) 
