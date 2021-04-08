
## MarkerMAG (linking MAGs with 16S rRNA marker genes)

![logo](doc/images/MarkerMAG_logo.jpg) 

[![pypi licence](https://img.shields.io/pypi/l/MarkerMAG.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version](https://img.shields.io/pypi/v/MarkerMAG.svg)](https://pypi.python.org/pypi/MarkerMAG) 


Publication
---
+ In preparation
+ **Weizhi Song** (songwz03[at]gmail.com)
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia


MarkerMAG modules:
---

1. Main module

    + `link`: linking MAGs with 16S rRNA marker genes
    
1. Supplementary modules

    + `rename_reads`: rename paired reads ([manual](doc/README_rename_reads.md))
    + `matam_16s`: assemble 16S rRNA genes with Matam, including subsample and dereplication ([manual](doc/README_matam_16s.md))
    + `uclust_16s`: cluster marker genes with Usearch ([manual](doc/README_uclust_16s.md))
    + `barrnap_16s`: identify 16S gene sequences with Barrnap ([manual](doc/README_barrnap_16s.md))
    + `subsample_reads`: subsample reads with Usearch ([manual](doc/README_subsample_reads.md))


Dependencies
---
 
 Dependencies of MarkerMAG are module-specific.
 
+ `link`: 
  [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/) and 
  [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

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



Important Notes :warning:
---

1. MarkerMAG assumes the id of paired reads in the format of `XXXX.1` and `XXXX.2`. The only difference is the last character.
   You can rename your reads with MarkerMAG's `rename_reads` module ([manual](doc/README_rename_reads.md)). 
   
1. Although you can use your preferred tool to reconstruct 16S rRNA gene sequences from the metagenomic dataset, 
   MarkerMAG does have a supplementary module (`matam_16s`) to reconstruct 16S using Matam. 
   Please refer to the manual [here](doc/README_matam_16s.md) if you want to use it.

1. :warning: All MAGs derived from a metagenomic dataset should be included in MarkerMAG run. (more details need to be added)


How to run:
---

+ Link 16S rRNA gene sequences with MAGs: 

      MarkerMAG link -p Soil -r1 R1.fastq -r2 R2.fastq -marker Soil_16S_uclust_0.999.fasta -mag refined_MAG -x fasta -t 12 -tmp -force

+ Preset parameters: 

      very_sensitive:  -min_clp_len 30 -min_clp_M_len 20 -s1_mpl 5  -s1_mplu 3  -min_M_len 30 -min_M_pct 20 -mismatch 3 -min_overlap_iden 99.9 -min_overlap_cov 25 -min_overlap_len 50 -min_overlap_num 3
      sensitive:       -min_clp_len 30 -min_clp_M_len 20 -s1_mpl 5  -s1_mplu 3  -min_M_len 30 -min_M_pct 20 -mismatch 3 -min_overlap_iden 99.9 -min_overlap_cov 30 -min_overlap_len 50 -min_overlap_num 5        
      default:         -min_clp_len 30 -min_clp_M_len 20 -s1_mpl 10 -s1_mplu 5  -min_M_len 30 -min_M_pct 25 -mismatch 3 -min_overlap_iden 99.9 -min_overlap_cov 35 -min_overlap_len 50 -min_overlap_num 5
      specific:        -min_clp_len 30 -min_clp_M_len 20 -s1_mpl 10 -s1_mplu 5  -min_M_len 30 -min_M_pct 30 -mismatch 2 -min_overlap_iden 100  -min_overlap_cov 55 -min_overlap_len 50 -min_overlap_num 8
      very_specific:   -min_clp_len 30 -min_clp_M_len 20 -s1_mpl 10 -s1_mplu 10 -min_M_len 30 -min_M_pct 35 -mismatch 1 -min_overlap_iden 100  -min_overlap_cov 75 -min_overlap_len 50 -min_overlap_num 10
      super_specific:  -min_clp_len 30 -min_clp_M_len 20 -s1_mpl 10 -s1_mplu 10 -min_M_len 30 -min_M_pct 35 -mismatch 1 -min_overlap_iden 100  -min_overlap_cov 85 -min_overlap_len 50 -min_overlap_num 10



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
