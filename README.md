
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


 
+ Dependencies for the `link` module:
  [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), 
  [seqtk](https://github.com/lh3/seqtk), 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
  [Samtools](http://www.htslib.org) and 
  [metaSPAdes](https://cab.spbu.ru/software/meta-spades/)

+ Dependencies for other supplementary modules can be found from their own manual page.
 
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

+ BioSAK has been tested on Linux/Mac, but NOT on Windows.
+ MarkerMAG is implemented in python3, you can install it with pip3:

      # install with 
      pip3 install MarkerMAG
        
      # upgrade with 
      pip3 install --upgrade MarkerMAG

+ If you clone the repository directly off GitHub you might end up with a version that is still under development.


Notes 
---

1. :warning: MarkerMAG assumes the id of paired reads in the format of `XXXX.1` and `XXXX.2`. The only difference is the last character.
   You can rename your reads with MarkerMAG's `rename_reads` module ([manual](doc/README_rename_reads.md)). 
   
1. Although you can use your preferred tool to reconstruct 16S rRNA gene sequences from the metagenomic dataset, 
   MarkerMAG does have a supplementary module (`matam_16s`) to reconstruct 16S using Matam. 
   Please refer to the manual [here](doc/README_matam_16s.md) if you want to use it.
   
1. :warning: All MAGs derived from a metagenomic dataset should be included in MarkerMAG run. (more details need to be added)


How to run:
---

+ Make sure you have bbmap, BLAST+ and Spades in your path.

+ Link 16S rRNA gene sequences with MAGs: 

      MarkerMAG link -p Soil -r1 R1.fasta -r2 R2.fasta -marker Soil_16S_uclust_0.999.fasta -mag refined_MAG -x fasta -t 12 -tmp -force


Output files:
---

+ #### Summary of identified linkages at genome level

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

+ #### Summary of identified linkages at contig level

    |Marker___Genome(total_number_of_linkages)	|Contig	|Paired	|Clipping	|Overlapped	|Step|
    |:---:|:---:|:---:|:---:|:---:|:---:|
    |Matam_100_310___cami_MAG_1(284)|	Contig_C12361|	72	|1|	0|	S1|
    |Matam_100_310___cami_MAG_1(284)|	Contig_C15695|	72	|1|	0|	S1|
    |Matam_100_310___cami_MAG_1(284)|	Contig_C52142|	61	|0|	0|	S1|
    |Matam_100_284___cami_MAG_16(109)|	Contig_C28768|	81	|28|	0|	S1|
    |Matam_75_135___cami_MAG_10(57)|	Contig_C4223|	44	|8	|0|	S1|
    |Matam_75_135___cami_MAG_10(68)|	Contig_C44260|	32|	3	|0|	S1|
    |Matam_100_162___cami_MAG_10(60)|	Contig_C4223|	42	|4|	0|	S1|
    |Matam_100_162___cami_MAG_10(56)|	Contig_C44260|	51	|3|	0|	S1|
    |Matam_50_12___cami_MAG_10(51)|	Contig_C4223|	42|	4	|0	|S1|
    |Matam_50_12___cami_MAG_10(42)|	Contig_C44260|	51|	3	|0	|S1|
    |Matam_50_142___cami_MAG_15(59)|	Contig_C685|	54|	5	|0|	S1|
    |Matam_100_235___cami_MAG_15(49)|	Contig_C685|	45|	4	|0|	S1|
    |Matam_75_12___cami_MAG_35(80)|	Contig_C65|	0|	0|	80	|S2|
    |Matam_100_25___cami_MAG_35(69)|	Contig_C65|	0|	0|	69	|S2|
    |Matam_100_23___cami_MAG_35(39)|	Contig_C65|	0|	0|	39	|S2|

    ![linkages](doc/images/linkages_plot.png)

+ #### Visualization of linkages
    ![linkages](doc/images/linking_reads.png)
   