
+ Install MarkerMAG on Katana with Python virtual environment

      module load python/3.7.3
      mkdir ~/mypython3env_MarkerMAG
      python3 -m venv --system-site-packages ~/mypython3env_MarkerMAG
      source ~/mypython3env_MarkerMAG/bin/activate
       
      # install with 
      pip3 install MarkerMAG
        
      # upgrade with 
      pip3 install --upgrade MarkerMAG


+ To run MarkerMAG, just run the following commands to activate python's virtual environment and load needed modules.

      module load python/3.7.3
      source ~/mypython3env_MarkerMAG/bin/activate
      module load R/4.0.2
      module load blast+/2.11.0
      module load bowtie/2.3.5.1
      module load samtools/1.10
      module load spades/3.14.0
      module load gcc/8.4.0
      module load boost/1.73.0-gcc8   
      module load java/8u201-jdk
      module load bbmap/38.51
      module load seqtk/20190219
      module load mafft/7.407
      module load perl/5.28.0
      module load hmmer/3.2.1
      module load bedtools/2.27.1
      module load barrnap/0.9
      module load usearch/11.0.667
      cd /srv/scratch/zID
      MarkerMAG link -p Soil -marker Soil_16S.fa -mag Soil_MAGs -x fa -r1 R1.fa -r2 R2.fa -t 12
      