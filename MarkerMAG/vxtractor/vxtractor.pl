#!/usr/bin/perl -w
# Program name: V-Xtractor
# Version: 2.1
# Purpose: Extraction of variable subregions from SSU and LSU rRNA sequences
# Copyright (c) Hartmann et al. 2010
# Contact: Martin Hartmann (contact(at)microbiome.ch)
# Swiss Federal Research Institute WSL, Forest Soils and Biogeochemistry, 
# 8903 Birmensdorf, Switzerland
# Programmer: Charles Howes (vxtractor(at)ch.pkts.ca)
# 
# Citation: Hartmann M, Howes CG, Abarenkov K, Mohn WW, Nilsson RH (2010)
# V-Xtractor: An open-source, high-throughput software tool to identify
# and extract hypervariable regions of small subunit (16S/18S) ribosomal
# RNA gene sequences. Journal of Microbiological Methods 83(2): 250-253.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have receive a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>
#
# Description:
# Hidden Markov Models (HMM) can be used to detect conserved regions
# on SSU or LSU rRNA sequences (for example) and can then be used to extract
# interjacent segments.  We have created HMMs so that each hypervariable
# region (SSU V1-V9 or LSU V1-V16) is flanked by two upstream (or left) and two 
# downstream (or right) HMMs.  At each position where necessary, a long and a
# short HMM was designed, i.e. V1leftlong.HMM and V1leftshort.HMM. For each
# gene (SSU or LSU), the HMMs are stored in three separate directories for the
# three taxonomic groups: bacteria, archaea, and fungi.

print STDERR "\nV-Xtractor v. 2.1. Copyright (c) Hartmann et al. 2010.\n";
print STDERR "\n";

use strict;
use Getopt::Long;
use POSIX qw(strftime);
use Carp qw(cluck);

# Record the full command line for logging purposes:
my $COMMAND=join(" ",$0,@ARGV);

# Parse command line:
my $OUTFILE;
my $EVALUE=0.01;
my $SCORE=0;
my $HMMDIR="HMMs/bacteria";
my @REGIONS;
my $OUTCSV;
my $INCLUDEHMM;
my $ALPHABET;
my $USEBITSCORE;
my $DEBUG;

my $result=GetOptions(
  "alphabet"=>\$ALPHABET,
  "bitscore"=>\$USEBITSCORE,
  "csv:s"=>\$OUTCSV,
  "debug:i"=>\$DEBUG,
  "evalue:f"=>\$EVALUE,
  "hmmdir:s"=>\$HMMDIR,
  "includehmm:s"=>\$INCLUDEHMM,
  "output=s"=>\$OUTFILE,
  "region:s"=>\@REGIONS,
  "score:f"=>\$SCORE,
);
my @INFILES=@ARGV;

# Checking for errors and usage message:
if (!@INFILES or !defined $OUTFILE or !$result) {
  die("Usage: $0 [-a] [-b] [-d] [-e evalue] [-s score] [-r region] [-i (long|short)]
   [-h hmmdirectory] [-c csvoutput] [-o outputfile] inputfile [inputfiles...]

  This program will analyze each sequence in each input file, looking
  for the HMMs in the hmm directory.

  Options:
    -o outputfile: Write the HMM region information to a FASTA file
    -c csvoutput: Write the HMM region information to a CSV file

    -h hmmdirectory: The directory containing HMM files named
       V[1-x]leftlong.HMM   V[1-x]leftshort.HMM
       V[1-x]rightlong.HMM  V[1-x]rightshort.HMM

    -r region: The regions to extract, in the following format:
      -r V1       -- the V1 region only
      -r .V1-V2.  -- the region from the left of V1 to the right of V2
      -r V1.-.V2  -- the region from the right of V1 to the left of V2
      -r .V3-.V7  -- the region from the left of V3 to the left of V7

    -i (long|short): Include HMM regions in the fasta output (default: exclude)
      The long or short HMM region will be chosen where applicable.

    -b: Use bitscore instead of evalue threshold (only use one or the other)
    -e evalue: Set the global evalue threshold (default: $EVALUE)
    -s score: Set the global score threshold (default: $SCORE)

    -a: Check that HMMs occur in alphabetical order in each sequence

    Example:
    $0 -a -r .V1-V3. -h HMMs/SSU/bacteria/ -o out.fasta  in.fasta
    --this will extract V1 through V3, for SSU bacteria, from the file in.fasta
    and save the results to out.fasta, checking correct order of V1, V2, and V3.
");
}

# Look for the specified input files:
if (my @m=grep {!-f} @INFILES) {
  die("These files were not found: @m\n");
}

# Check to see if an output file was specified:
if (!defined $OUTCSV and !defined $OUTFILE) {
  die("Either -c or -o must be specified, or else there's no point in
running this program.\n");
}

# Check that the -i option is either 'short' or 'long':
if (defined $INCLUDEHMM) {
  $INCLUDEHMM=lc($INCLUDEHMM);
  if ($INCLUDEHMM !~ m/^(long|short)$/i) {die("The -i option can be either 'long' or 'short', not '$INCLUDEHMM'\n");}
}

# Testing for the hmmscan program:
#my     $HMMSCAN=`which hmmscan`; chomp($HMMSCAN);
# if(!-f $HMMSCAN || !-x _) { die("$HMMSCAN isn't an executable file!"); }
my $HMMSCAN="hmmscan";  # Hope it's in the path.

# Getting the list of HMM files as a hash:
my $HMMcount=0;
opendir(DIR,$HMMDIR) or die("$HMMDIR: $!");
my %HMMS=map {$HMMcount++;$_=>1} grep {-f && m/[.]HMM$/i} map {"$HMMDIR/$_"} readdir(DIR);
closedir(DIR);
my $badness=0;
foreach my $n (sort keys %HMMS) {
  if ($n!~ m/^(\S+)(left|right)(short|long)[.]HMM$/i) {
    warn("The filename '$n' doesn't match the pattern '???(left|right)(short|long)[.]HMM'\n");
    $badness++;
  }
}
if ($badness) {exit 2;}
if ($HMMcount==0) {die("No HMMs were found in $HMMDIR");}

# Get the e-values and scores from each HMM file.
# The hmmscan program can be given an e-value and a score cutoff on
# the command line with -e and -T; we're storing this information
# inside the HMM files in a DESC line.  Format:
# DESC  evalue=0.005  score=0
my %HMMevalue=map {$_=>getevalue($_)} keys %HMMS;
#print "HMMevalue:\n",map {"$_=>$HMMevalue{$_}\n"} sort keys %HMMevalue;
my %HMMscore=map {$_=>getscore($_)} keys %HMMS;
#print "HMMscore:\n",map {"$_=>$HMMscore{$_}\n"} sort keys %HMMscore;

# Turn the HMM filenames into more useful names by stripping the
# directory and suffix (in HMMbasenames) and also the left/right
# long/short words (in HMMnames).  Each hash contains the full
# filename(s) of the unstripped files.
my %HMMbasenames=map {basename($_,".HMM")=>$_} keys %HMMS;
my %HMMnames;
map {$HMMnames{basename($_,"(left|right)(long|short)[.]HMM")}.=" ".$_} keys %HMMS;
#print "HMMnames:\n",join("\n",sort keys %HMMnames),"\n";

# Error checking to make sure that each region has a left and right
# long and short hmm, or else it will be hard to define the regions:
foreach my $x (keys %HMMnames) {
  my %set=map {basename($_,"[.]HMM")=>1} grep {/./} split(/ /,$HMMnames{$x});
  if (!defined $set{"${x}leftshort"}) {die("${x}leftshort.HMM was not found\n");}
  if (!defined $set{"${x}leftlong"}) {die("${x}leftlong.HMM was not found\n");}
  if (!defined $set{"${x}rightshort"}) {die("${x}rightshort.HMM was not found\n");}
  if (!defined $set{"${x}rightlong"}) {die("${x}rightlong.HMM was not found\n");}
}

# Optimization: sort the keys only once, and store them for later use:
my @HMMS=sort keys %HMMS;
my @REGLIST=sort keys %HMMnames;

# If no regions were specified on the command line, return *all* of them.
# If the regions do not appear in the sequence in the same order as
# they appear in @REGLIST (alphabetical order), there will be an error
# about things being out of order if -a is given. 
if (grep {/^$/} @REGIONS) {
  die("Error: -r needs to be followed by an argument.\n");
}
if (@REGIONS==0) {@REGIONS=@REGLIST;}
my %REGTODO=map {$_=>1} map {/([^.-]+)/g} @REGIONS; 
my @REGTODO=sort keys %REGTODO;

my @HEADERS=grep { /^(.*)(left|right)(short|long)$/;defined $REGTODO{$1}} sort keys %HMMbasenames;

# Find size of all input files, for the progress bar later:
my $INSIZE=0;
map {(-e $_) && ($INSIZE+= -s $_)} @INFILES; # I'm using map more and more

# Progress bar variables:
my $last=time();
my $start=time();
my $offset=0;

# Prepare the output files:
if (defined $OUTCSV) {
  open(OUTCSV,">$OUTCSV") or die("$OUTCSV: $!");
  print OUTCSV join(",","Command:",$COMMAND),"\n";
  print OUTCSV "Options:,";
  print OUTCSV $INCLUDEHMM?"Include $INCLUDEHMM HMM in fasta":"Don't include HMM in fasta";
  print OUTCSV "\n";
  print OUTCSV join(",","Sequence",@HEADERS,@REGTODO),"\n";
}
if (defined $OUTFILE) {open(OUT,">$OUTFILE") or die("$OUTFILE: $!");}

my $SUMMARY;

# Process the input files:
map {my $aa=$_;processfile($aa);$offset+=-s $aa} @INFILES;

# Close everything that's open:
if (defined $OUTFILE) {close(OUT);}
if (defined $OUTCSV) {close(OUTCSV);}

# Blank line at the end:
print STDERR "\n";

#----------------------------------------------------------------------
# Process a file:
sub processfile {
  my $file=$_[0];
  if (!-e $file) {die("$file doesn't exist.\n");}
  if (!-f $file) {die("$file isn't a file.\n");}
  open(IN,"<$file") or die("$file: $!\n");
  do {
    local $/="\n>";  # Record separator is '\n>'
    while (<IN>) {
      my $v=$_;
      $v=~s/>\n?$//s;
      process_sequence($file,$v);
    };
  };
  close(IN);
}

# Process a sequence:
sub process_sequence {
  my ($file,$seq)=@_;

  # Update the progress meter:
  $last=progress(tell(IN)+$offset,$INSIZE,$start,$last);

  # Remove carriage returns from windows text files:
  $seq=~s/\r//g;

  # Tidy-up and separate input data:
  $seq=~s/^>?([^\n]+)\n//; # Remove name/description from sequence
  my $name=$1;  # Store it here
  $seq=~tr/\n\r//d; # Delete newlines of any kind from sequence

  # Search for HMMs.  This only looks for HMMs in @REGTODO:
  my %result=map {$_=>matchhmm($_,$name,$seq)} grep {/([^\/]+)(left|right)(short|long)[.]HMM$/;defined $REGTODO{$1}} @HMMS;

  $SUMMARY=join("",map {"$_=>$result{$_}\n"} sort keys %result);

  # Calculate the regions, stick into %result
  map {$result{$_}=getstartend($seq,\%result,"${_}left","${_}right",$file,$name);} @REGTODO;

  # Check that the regions found are in the same order as @REGIONS:
  my $error="";
  if (defined $ALPHABET) {
    my $last=0;
    my $lastregion="";
    foreach my $s (@REGLIST) {
      if (!defined $result{$s}) {next;} # Was not scanned
      if ($result{$s} eq "notfound") {next;}
      if ($result{$s} !~ m/^(\d+)-(\d+)/) {die("Got '$result{$s}' for region $s in sequence $name\n");}
      if ($1<$last) {
        $error.=" $s<$lastregion";
		warn("$file:$name: The regions $lastregion-$s are in the wrong order: $last>$1\n");
		#die("$file:$name: The regions $lastregion-$s are in the wrong order: $last>$1\n$SUMMARY");
      }
      $last=$1; # It was $2, but that was too stringent
      $lastregion=$s;
    }
    if ($error ne "") {$error="Error:$error";}
  }

  # Write HMM matches using an over-complicated expression:
  if (defined $OUTCSV) {
    print OUTCSV join(",",map({"\"'$_\""} $name,@result{(grep {/([^\/]+)(left|right)(short|long)[.]HMM$/;defined $REGTODO{$1}} @HMMS),@REGTODO},$error)),"\n";
  }

  # Extract sequence data for fasta file:
  foreach my $x (@REGIONS) {
    my $region="";
    if ($x=~m/^([^.-]+)$/)   # -r V1
      {$region=getstartend($seq,\%result,"$1left","$1right",$file,$name,1);}

    if ($x=~m/^[.]([^.-]+)-[.]([^.-]+)$/)  # -r .V1-.V2
      {$region=getstartend($seq,\%result,"$1left","$2left",$file,$name,1);}
    if ($x=~m/^[.]([^.-]+)-([^.-]+)[.]$/)  # -r .V1-V2.
      {$region=getstartend($seq,\%result,"$1left","$2right",$file,$name,1);}
    if ($x=~m/^([^.-]+)[.]-[.]([^.-]+)$/)  # -r V1.-.V2
      {$region=getstartend($seq,\%result,"$1right","$2left",$file,$name,1);}
    if ($x=~m/^([^.-]+)[.]-([^.-]+)[.]$/)  # -r V1.-V2.
      {$region=getstartend($seq,\%result,"$1right","$2right",$file,$name,1);}

    if ($region!~m/^(\d+)-(\d+)/) {
  #   warn("\n$file:$name:$x: region not found\n");
      next;
    }

    # Write out the fasta sequence, if logically possible:
    if ($1<=$2) {
      my $sname="$x $region ";
      $sname.=($INCLUDEHMM?"include${INCLUDEHMM}hmm":"withouthmm");
      $sname=~s/[\s.-]/_/g; 
      $sname=~s/__/_/g;
      $sname=~s/^_//g;
      print OUT ">${name}_$sname\n";
      print OUT getsubsequence($seq,$region),"\n";
    }
  }
}

# Find the start and end of the requested region:
sub getstartend {
  my ($seq,$bounds,$in1,$in2,$file,$name,$err)=@_;
  my %b=%$bounds;

  # Get information for region 1:
  (my $reg1=$in1)=~s/(left|right)$//;
  my $r=$HMMnames{$reg1};
  if (!defined $r) {die("There is no region named $reg1!  Choose from @REGLIST\n");}

  # Get information for region 2:
  (my $reg2=$in2)=~s/(left|right)$//;
  my $s=$HMMnames{$reg2};
  if (!defined $s) {die("There is no region named $reg2!  Choose from @REGLIST\n");}

  # Get filenames for region 1 and 2:
  my %k=map {basename($_,"[.]HMM")=>$_} grep {/./} split(/ /,$r.$s);

  my $leftlong=$b{$k{"${in1}long"}};
  my $leftshort=$b{$k{"${in1}short"}};
  my ($left,$leftlabel)=checkhmm($leftlong,$leftshort);
  if ($leftlabel ne "") {$leftlabel="${leftlabel}left"}

  my $rightlong=$b{$k{"${in2}long"}};
  my $rightshort=$b{$k{"${in2}short"}};
  my ($right,$rightlabel)=checkhmm($rightlong,$rightshort);
  if ($rightlabel ne "") {$rightlabel="${rightlabel}right"}

  if ($left!~m/^(\d+)-(\d+)/) {return "notfound";}
  my ($l1,$l2)=($1,$2);
  if ($right!~m/^(\d+)-(\d+)/) {return "notfound";}
  my ($r1,$r2)=($1,$2);

# Needs changing: XXX
  #if ($l1>=$l2) {warn("In sequence $name, region $reg1 starts at $l1 but ends at $l2");}
  #if ($l2>=$r1) {warn("In sequence $name, region $reg1 ends at $l2, region $reg2 starts at $r1");}
  #if ($r1>=$r2) {warn("In sequence $name, region $reg1 starts at $r1 but ends at $r2");}
  if ($l2>=$r1) {
    $leftlabel.="-wrongorder";
	if (!defined $err) {
      warn("In sequence $name, region $reg1 starts at $l2 and ends at $r1, which is backwards.\n");
      #die("In sequence $name, region $reg1 starts at $l2 and ends at $r1, which is backwards.\n$SUMMARY\n");
	}
  }

  if ($INCLUDEHMM) {
    return "$l1-$r2".$leftlabel.$rightlabel;
  } else {
    return ($l2+1)."-".($r1-1).$leftlabel.$rightlabel;
  }
}

# This routine checks the long and short hmm for correctness.  The
# function returns either the short or long hmm, plus a diagnostic label.
sub checkhmm {
  my ($hmmlong,$hmmshort)=@_;

  my $hmm;
  my $hmmlabel="";

  # There are five cases:
  # 1. Neither HMM detected:  no extract, no label
  if ($hmmshort eq "notfound" and $hmmlong eq "notfound" )
    {return ("notfound","");}
  # 2. Long HMM is detected, short isn't:  extract long, label "HMM=long"
  if ($hmmshort eq "notfound" and $hmmlong ne "notfound") 
    {return ($hmmlong," HMM=long");}
  # 3. Short HMM is detected, long isn't:  extract short, label "HMM=short"
  if ($hmmlong eq "notfound" and $hmmshort ne "notfound")
    {return ($hmmshort," HMM=short");}
  # 4. Long and short HMM are detected:  extract $INCLUDEHMM, no label
  $hmm=$hmmlong;
  if (defined $INCLUDEHMM and $INCLUDEHMM eq "short") {$hmm=$hmmshort;}
  # 5. Long and short HMM are mismatched:  extract long, label "HMM=mismatch"
  if (!inset($hmmlong,$hmmshort)) { $hmm=$hmmlong; $hmmlabel=" HMM=mismatch"; }
  return ($hmm,$hmmlabel);
}

# Is the short region inside the long region?
# It should be.  If not, it's very strange.
sub inset {
  my ($long,$short)=@_; # Regions
  my ($l1,$l2)=$long=~m/^(\d+)-(\d+)/;
  my ($s1,$s2)=$short=~m/^(\d+)-(\d+)/;
  if (!($l1<=$s1 and $s1<=$l2)) {return 0;}
  if (!($l1<=$s2 and $s2<=$l2)) {return 0;}
  return 1;
}

# Get a region of a sequence.  One hiccup: strings start at 0, but
# sequences are numbered starting at 1.  Another hiccup: the sequence
# starts just before the starting position and ends just after the
# ending position.  So, if start=end, then return the single letter at
# that position.
sub getsubsequence {
  my ($seq,$startend)=@_;
  if (!defined $seq) {die("\nsequence not defined");}
  if (!defined $startend) {return "";}
  if ($startend=~m/notfound/) {return "";}
  my ($s,$e)=($startend=~m/(\d+)-(\d+)/);
  if ($s>=$e) {return "";}
  if (!defined $s or !defined $e) {die("\nstartend bad format");}
  return substr($seq,$s-1,$e-$s+1);
}

# Do an HMM search.  This puts the sequence into a temporary file,
# calls hmmpfam with the appropriate arguments, and checks the answer
# thoroughly.
sub matchhmm { # Hmm=filename, Name=sequencename, Seq=sequence
  my ($hmm,$name,$seq)=@_;

  # Write the sequence to a tempfile:
  open(TMP,">temp.pid.$$") or die("temp.pid.$$: $!");
  print TMP ">$name\n$seq\n";
  close(TMP);

  # Run hmmpfam:
  my $criteria="--domE $HMMevalue{$hmm}";
  if (defined $USEBITSCORE) {
    my $criteria="--domT $HMMscore{$hmm}";
  }
  my $cmd="$HMMSCAN --max $criteria $hmm temp.pid.$$";
    # --domT score = bit score threshold (can be negative)
    # --domE thresh = evalue threshold
  my $output=`$cmd 2>&1`;
  my $oldout=$output;
  unlink("temp.pid.$$");

  # Error checking:
  if (!defined $output) {die("$cmd didn't work");}

#  print "Command:\n$cmd\n";
#  print "Output:\n$output\n";
  if ($output=~m/(FATAL:[^\n]*)/s) {die("\nCMD=$cmd \n$1\n");}
  if ($output=~m/No individual domains that satisfy reporting thresholds/) {return "notfound";}
  if ($output=~m/No hits detected that satisfy reporting thresholds/) {return "notfound";}
  if ($output=~m/No targets detected that satisfy reporting thresholds/) {return "notfound";}

  # Data extraction:
  $output=~s/\r//g; # Delete carriage returns (fine for windows and unix)

  # Hmmer2:
  #$output=~s/.*\nParsed for domains:\nModel[^\n]*\n----[- ]*\n//s; # Cut
  #$output=~s/\nAlignments of top-scoring domains:.*//s; # Cut
  #my ($start,$end)=$output=~m/^\S+\s+\S+\s+(\d+)\s+(\d+)/s; # Match it

  # Hmmer3:
  # Always present, even when no matches found:
  $output=~s/.*\nDomain annotation for each model \(and alignments\):\n//s;
  $output=~s/\nInternal pipeline statistics summary:.*//s;
  # Only present when a match found:
  $output=~s/>>[^\n]*\n   #[^\n]*\n ---[- ]*\n//s;
  $output=~s/\n  Alignments for each domain:.*//s;

  # Take the 'alifrom/ali to' values from the best domain:
  my @outlist=grep {/./} split(/\n/,$output);
  my $bestoutput=$outlist[0];
  my $bestev=10;
  for (my $x=0;$x<@outlist;$x++) {
    my $ev=substr($outlist[$x],31,9);
    $ev=~s/\s+//g;
    $ev=$ev+0;
    if ($ev !~ m/^[\d.eE-]+$/) {
      die("\nMangled output: ev=($ev)\n".join("\n",@outlist)."\n\nCmd: $cmd\nFrom:\n$oldout");
    }
    if ($ev+0<$bestev) {
      $bestev=$ev;$bestoutput=$outlist[$x]
    }
  }
#   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
# ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
#   1 !   31.3   0.0   1.1e-11   1.1e-11       4      38 ..       5      39 ..       2      41 .. 0.93
  my $values=substr($bestoutput,78,15); # Doing envfrom-env to
  my ($start,$end)=$values=~m/^\s*(\d+)\s+(\d+)$/s; # Match it
  if (defined $DEBUG) {
    $values=substr($bestoutput,20,9);
    my ($evalue)=$values=~m/^\s*(\S+)\s*$/s; # Match it
    $values=substr($bestoutput,7,6);
    my ($bitscore)=$values=~m/^\s*(\S+)\s*$/s; # Match it
    open(DEBUG,">>debug.csv");
    print DEBUG "$name,$hmm,$evalue,$bitscore,$start,$end\n".join("\n",@outlist)."\n";
    close(DEBUG);
  }

  # Return something if not found:
  if (!defined $start) {return "notfound";}
  # Shouldn't happen:
  if ($start>=$end) {die("\nhmmpfam reported a match starting at $start but ENDING at $end (wrong order) in $hmm and sequence $name");}
  # Success!
  return "$start-$end"; # = the region that was found
}


# Extract the e-value that's been hardcoded into the HMMs:
sub getevalue {
  my $hmm=$_[0];
  open(HMMIN,"<$hmm") or die("$hmm: $!");
  my @l=grep {/^DESC\s(.*\s)?evalue=[0-9.e-]+/i} <HMMIN>;
  close(HMMIN);
  if (@l) {return $l[0]=~m/\sevalue=([0-9.e-]+)/i;} # First evalue
  return $EVALUE;  # Default value
}

# Extract the score that's been hardcoded into the HMMs:
sub getscore {
  my $hmm=$_[0];
  open(HMMIN,"<$hmm") or die("$hmm: $!");
  my @l=grep {/^DESC\s(.*\s)?score=[0-9.-]+/i} <HMMIN>;
  close(HMMIN);
  if (@l) {return $l[0]=~m/\sscore=([0-9.e-]+)/i;} # First score
  return $SCORE;  # Default value
}

# Remove the directory name from the string:
sub basename {
  (my $aa=$_[0])=~s/.*\///;
  if (defined $_[1]) { $aa=~s/${_[1]}$//i; }
  return $aa;
}

# Progress bar routine:
sub progress {
  my ($pos,$total,$stime,$ltime)=@_;
  if ($pos!=$total and $total!=0) {
    if ($ltime==time() or $stime==$ltime or $total==0) {return time();} # Not ready for a screen update
  }
  $ltime=time();
  my $done=$pos/$total;
  my $ttime=($ltime-$stime)/$done;
  print STDERR "\r";
  print STDERR strftime("%T",localtime($stime))." ";
  print STDERR sprintf("[%-30s]","="x(int($done*30)));
  print STDERR sprintf("%5s",int($done*1000)/10),"% ";
  print STDERR strftime("%T",localtime($stime+$ttime))."  ";
  if ($pos!=$total) {
    print STDERR "Remaining: ".nicetime($stime+$ttime-time());
  } else {
    print STDERR "    Total: ".nicetime($ttime)."\n";
  }
  return $ltime;
}

# Format a number of seconds as hour:minute:second.  Can't handle
# negative numbers, and doesn't know about days.
sub nicetime {
  my $d=$_[0];
  my $hour=int($d/3600); $d-=$hour*3600;
  my $min=int($d/60); $d-=$min*60;
  my $sec=int($d);
  return sprintf("%d:%02d:%02d",$hour,$min,$sec);
}
