
## MarkerMAG: linking MAGs with 16S rRNA marker genes using paired-end short reads

[![pypi licence](https://img.shields.io/pypi/l/MarkerMAG.svg)](https://opensource.org/licenses/gpl-3.0.html)
[![pypi version](https://img.shields.io/pypi/v/MarkerMAG.svg)](https://pypi.python.org/pypi/MarkerMAG) 


Publication
---

+ **Song WZ**, Zhang S, Thomas T* (2021) MarkerMAG: linking metagenome-assembled genomes (MAGs) with 16S rRNA marker genes using paired-end short reads (under review)
+ Contact: Dr. Weizhi Song (songwz03@gmail.com), Prof. Torsten Thomas (t.thomas@unsw.edu.au)
+ Center for Marine Science & Innovation, University of New South Wales, Sydney, Australia


Updates
---

+ 2022-03-12 - A [demo dataset](https://doi.org/10.5281/zenodo.6466784) (together with command) has been prepared! You can use it to check if MarkerMAG is installed successfully on your system.


How it works
---

+ Workflow of MarkerMAG
![linkages](doc/images/MarkerMAG_workflow.png)

+ GC content bias
  
  Read coverage of MAGs and their linked 16S rRNA genes might be biased by guanine-cytosine (GC) content [[Reference](https://doi.org/10.1093/nar/gks001)].
    Read coverage are weighted by GC content bias before estimating the copy number of 16S rRNA genes in MAGs. 
    GC content bias is calculated as described [here](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/GCBiasReport_fDG.htm). 
    An example of GC content bias from the MBARC-26 dataset that we used for benchmarking MarkerMAG is [here](doc/README_GC_bias.md).


MarkerMAG modules
---

1. Main module

    + `link`: linking MAGs with 16S rRNA marker genes
    
1. Supplementary modules

    + `rename_reads`: rename paired reads ([manual](doc/README_rename_reads.md))
    + `matam_16s`: assemble 16S rRNA genes with Matam ([manual](doc/README_matam_16s.md))
    + `barrnap_16s`: identify 16S rRNA genes from genomes/MAGs with Barrnap ([manual](doc/README_barrnap_16s.md))


How to install
---

+ MarkerMAG has been tested on Linux and MacOS, but NOT on Windows.


+ MarkerMAG is implemented in [python3](https://www.python.org), It can be installed with pip. 
  Software dependencies need to be in your system path. 
  Dependencies for the `link` module include 
  [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), 
  [Barrnap](https://github.com/tseemann/barrnap), 
  [seqtk](https://github.com/lh3/seqtk), 
  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), 
  [Samtools](http://www.htslib.org), 
  [HMMER](http://hmmer.org), 
  [metaSPAdes](https://cab.spbu.ru/software/meta-spades/) and 
  [Usearch](https://www.drive5.com/usearch/).
  Dependencies for the supplementary modules are provided in their corresponding manual page.
  
      # install with 
      pip3 install MarkerMAG
      
      # install a specific version of MarkerMAG (e.g. 1.1.24)
      pip3 install MarkerMAG==1.1.24
        
      # upgrade with 
      pip3 install --upgrade MarkerMAG

+ A Conda package for MarkerMAG is also available, which will install the third-party dependencies for you.
  However, you'll need to install usearch on your own as it's not available in conda due to license issue.

      conda create -n MarkerMAG -c bioconda -c songweizhi MarkerMAG
      conda activate MarkerMAG
      MarkerMAG -h

+ :warning: If you clone the repository directly off GitHub you might end up with a version that is still under development.


+ [Here](doc/README_example_cmds.md) are some example commands for UNSW Katana users.


How to run
---

+ MarkerMAG’s input consists of 
   1. A set of user-provided MAGs
   2. A set of 16S rRNA gene sequences (either user-provided or generated with the `matam_16s` module) 
   3. The **quality-filtered** paired-end reads used to generate the data above. 
      If the input reads are provided in fastq format, MarkerMAG will first convert them into fasta format.
   
+ :warning: MarkerMAG is designed to work with paired short-read data (i.e. Illumina). It assumes the id of reads in pair in the format of `XXXX.1` and `XXXX.2`. The only difference is the last character.
   You can rename your reads with MarkerMAG's `rename_reads` module ([manual](doc/README_rename_reads.md)). 

+ Although you can use your preferred tool to reconstruct 16S rRNA gene sequences from the metagenomic dataset, 
   MarkerMAG does have a supplementary module (`matam_16s`) to reconstruct 16S rRNA genes. 
   Please refer to the manual [here](doc/README_matam_16s.md) if you want to give it a go.

+ Link 16S rRNA gene sequences with MAGs ([demo dataset](https://doi.org/10.5281/zenodo.6466784)): 

      MarkerMAG link -p Demo -r1 demo_R1.fasta -r2 demo_R2.fasta -marker demo_16S.fasta -mag demo_MAGs -x fa -t 12


Output files
---

1. Summary of identified linkages at genome level:

    | Marker | MAG | Linkage | Round |
    |:---:|:---:|:---:|:---:|
    | matam_16S_7   | MAG_6 | 181| Rd1 |
    | matam_16S_12  | MAG_9 | 102| Rd1 |
    | matam_16S_6   | MAG_59| 55 | Rd2 |

2. Summary of identified linkages at contig level (with figure):

    |Marker___MAG (linkages)	|Contig	        |Round_1	|Round_2	|
    |:---:|:---:|:---:|:---:|
    |matam_16S_7___MAG_6(181)	            |Contig_1799	|176	    |0          |
    |matam_16S_7___MAG_6(181)	            |Contig_1044	|5	        |0          |
    |matam_16S_12___MAG_9(102)	            |Contig_840	    |102	    |0          |
    |matam_16S_6___MAG_59(39)	            |Contig_171	    |0	        |55         |

   ![linkages](doc/images/linkages_plot_2.png)

3. Copy number of linked 16S rRNA genes.


4. Visualization of individual linkage.
  
   MarkerMAG supports the visualization of identified linkages (needs [Tablet](https://ics.hutton.ac.uk/tablet/)). 
   Output files for visualization ([example](doc/vis_folder)) can be found in the [Prefix]_linkage_visualization_rd1/2 folders. 
   You can visualize how the linking reads are aligned to MAG contig and 16S rRNA gene by double-clicking the corresponding ".tablet" file. 
   Fifty Ns are added between the linked MAG contig and 16S rRNA gene.
 
   ![linkages](doc/images/linking_reads.png)
 
   *If you saw error message from Tablet that says input files format can not be understood, 
   please refer to [here](https://github.com/cropgeeks/tablet/issues/15) for a potential solution.


