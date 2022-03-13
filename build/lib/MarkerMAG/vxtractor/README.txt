V-Xtractor: An open-source, high-throughput software tool to identify and extract hypervariable regions of small and large subunit ribosomal RNA gene sequences of bacterial, archaeal, and fungal origin. 

Source code available at:
http://www.microbiome.ch/web/Tools.html

Version:	2.1
Purpose:	Extraction of variable subregions from SSU or LSU rRNA sequences
Copyright:	Hartmann et al. 2010
Contact:	Martin Hartmann (contact(at)microbiome.ch)
Affiliation:	Swiss Federal Research Institute WSL, Forest Soils and Biogeochemistry, 8903 Birmensdorf, Switzerland
Programmer:	Charles Howes (vxtractor(at)ch.pkts.ca)

The full installation instructions can be found in the Users Guide.
Quick installation instructions can be found below.


1) Perl and HMMER

Perl is available at http://www.perl.org/

Installing HMMER version 3
note: V-Xtractor 2.1 requires HMMER version 3. An older version of V-Xtractor running on HMMER version 2 can be provided on request.

 - download			http://hmmer.janelia.org/#download  
 - uncompress/unpack		tar zxf hmmer-3.0.tar.gz
 - move to new directory	cd hmmer-3.0
 - configure			./configure
 - build			make
 - automated tests		make check
 - automated install		make install

Precompiled binaries are available for some operating systems. They can be found in the "binaries/" directory after the unpacking the tarball.
Installing the package takes nothing more than moving these binaries wherever you want them (e.g. /usr/local/bin).


2) V-Xtractor

 - download	at http://www.microbiome.ch/web/Tools.html
 - unpack	"unzip vxtractor.pl"


3) Input file and running V-Xtractor

Use sequences in FASTA formatted input file. Copy FASTA file to directory containing vxtractor.pl. Move into this directory with "cd path/". Typing "perl vxtractor.pl" lists all options of the program.

Example:
perl vxtractor.pl -a -r .V1-V3. -h HMMs/SSU/bacteria/ -o out.fasta  in.fasta
-- this will extract V1 through V3, from bacterial small subunit rRNA gene sequences found in the file in.fasta and save the results to out.fasta, checking the correct order of V1, V2, and V3.
