#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use File::Basename qw(basename dirname);
use List::Util 'first';
use lib "$FindBin::Bin/./lib";
use Parallel::ForkManager;
use Pod::Usage;
use Env;
####################################################################################
#
# Written by Yi Liao (01/04/2022)
#
# Usage： Genotyping large insertions (equal or larger than 50 bp) for population-scale genome assemblies
#
# Dependencies: UCSC's utilities;stretcher in the EMBOSS package; Perl module, Parallel::ForkManager
#
####################################################################################
my ($svraw,$UCSC,$queries,$singmaf,$stretcher,$refsize,$refname,$speciesname,$refseq,$p,$OUTPUTdir,$help);

GetOptions(
    "ins=s"=>\$svraw,
    "refseq=s"=>\$refseq,
    "refsize=s"=>\$refsize,
    "refname=s"=>\$refname,  
    "query=s"=>\$queries,
    "singmaf=s"=>\$singmaf,
    "stretcher=s"=>\$stretcher,
    "ucsc=s"=>\$UCSC,
    "speciesname=s"=>\$speciesname,
    "t=i"=>\$p,
    "outdir=s"=>\$OUTPUTdir,
    "help"=>\$help,
);

############# Help information ##########
if ($help){
print <<"END.";
Usage: Perl $0 --ins INS.combined.txt -query /path/to/queriesGenome -sing /path/to/sing/ --stretcher /path/to/stretcher --refseq MH63 --refsize MH63.sizes --refname MH63 --t 12

--ins          the file of the merged non-redundant insertions (> 49 bp)               [REQUIRED]
--query        folder kept all query genome sequences                                  [REQUIRED]
--singmaf      folder kept single-coverage maf files (on the reference side)           [REQUIRED]
--refsize      chromosome/contig sizes from the reference                              [REQUIRED]
--refseq       reference genome sequence                                               [REQUIRED]
--stretcher    path to stretcher                                                       [REQUIRED]    
--refname      reference name, e.g. hg38, IRGSP1 [REQUIRED]
--t            how many threads used [default: 4]                                      [Optional]
--refname      reference name, e.g. hg38, IRGSP1                                       [Optional]
--ucsc         path for kent's programs                                                [Optional]
--outdir       set the output dir name                                                 [Optional]
--help         print this help information
END.
exit;
}

my $message = "Try 'perl $0 --help' for more information.";

if (! $svraw){
    pod2usage( {
           -message => $message,
           -verbose => 0,
           -exitval => 1 } );
    }
##########################################
$UCSC ||="$Bin/pub/UCSC/";
#$stretcher ||=".";
my $library_path = "$Bin/pub/Stretcher";
$ENV{'LD_LIBRARY_PATH'} = "$library_path:$ENV{'LD_LIBRARY_PATH'}";
$stretcher ||= "$Bin/pub/Stretcher/";
#$OUTPUTdir ||= ".";
#$OUTPUTdir =~ s/\/$//;
$speciesname ||= "Unknown";
$OUTPUTdir ||= "VCFtemp";
$p ||= 4;
mkdir($OUTPUTdir) unless(-d $OUTPUTdir);
mkdir("$OUTPUTdir/err") unless(-d "$OUTPUTdir/err");

################################## Split Reference by chromosomes and obtain chromosome names
my $ref_file = basename ($refseq);
my $ref_path = dirname ($refseq);
my $byChr = "$ref_path/byChr";

my @chromosomes;

if ( -e $byChr and -d $byChr) {
    opendir (DIR1,"$ref_path/byChr") || die "Error in opening dir\n";
    while ((my $filename = readdir (DIR1))) {
          if ($filename=~/$refname.(.*).fa/) {
               push (@chromosomes,$1);
          }
    }

} else {

    mkdir ("$ref_path/byChr") unless (-d "$ref_path/byChr");
    `${UCSC}faSplit byname $refseq Chr`; # The splited chrmosome sequences will in the current directory
    `mv $refname.*.fa $ref_path/byChr`; # Move the fasta to the directory where reference stored

    opendir (DIR1,"$ref_path/byChr") || die "Error in opening dir\n";
     while ((my $filename = readdir (DIR1))) {
       if ($filename=~/$refname.(.*).fa/) {
            push (@chromosomes,$1);
       }
     }
}

################################## Split maf by chromosomes and obtain query genome names
opendir (DIR2,"$singmaf") || die "Error in opening dir singmaf\n";

my @Genome;
`touch splits.bed`;
  while ((my $filename1 = readdir (DIR2))) {
         if ($filename1=~/$refname.(.*).sing.maf/) {
             push (@Genome,$1);
 
             if ( -e "$singmaf/$1" and -d "$singmaf/$1") {
             } else {
                mkdir("$singmaf/$1"); #  `mkdir $singmaf/$1`;
                `${UCSC}mafSplit splits.bed $refname\_$1\_ $singmaf/$filename1 -byTarget -useFullSequenceName`;
                `mv $refname\_$1\_*.maf $singmaf/$1`;
             }   
         }
}
`rm splits.bed`;

################################## Split queries by chromosomes
foreach my $qg (@Genome) {
`${UCSC}faSplit byname $queries/$qg Chr`;
`mv $qg.*.fa $queries`;
}

################################## Split the svdata by chromosome
my $sv_file = basename($svraw);
my $sv_path = dirname($svraw);

open SV, "$svraw" or die "$!";
while (<SV>) {
chomp;
my @svs = split (/\t/,$_);
my @chrs = split (/\./,$svs[0]);
   foreach my $ctg (@chromosomes) {
     print "$ctg\n";
      if ($chrs[1]=~/^$ctg$/) {
              my $file;
              open ($file, '>>', "$sv_path/$refname.$ctg.INS.bed") or die "Could not open file '$file' $!";
              print $file "$_\n";
      }
    }
}

###############################################################################################

foreach my $ind (@Genome) {
        my $pm = Parallel::ForkManager->new($p);
DATA_LOOP:
       	foreach my $scf (@chromosomes) {
             my $pid = $pm->start and next DATA_LOOP;
	         mkdir("$OUTPUTdir/$ind.VAR.$scf");
	         print "Currently are processing chromosome $scf of sample $ind !!!\n";
	         my $liftover_coord = &liftOver("$singmaf/$ind/$refname\_$ind\_$scf.maf");

################################################################################################
open IN, "$sv_path/$refname.$scf.INS.bed" or die "$!";
open VCF, ">$OUTPUTdir/$refname.$ind.$scf.vcf" or die "$!";
open ERR, ">$OUTPUTdir/err/$refname.$ind.$scf.err.vcf" or die "$!";

while (<IN>) {
chomp;
my $line = $_;
my ($chr,$beg,$end,$len,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
my $ref = (split(/\./,$chr))[0];
my $alt = (split(/\./,$qchr))[0];
my $alt_chr = (split(/\./,$qchr))[1];
my $flank;

if ($len < 1000) {
$flank=100;
} elsif ($len >=1000 and $len < 5000) {
$flank = 500;
} else {
$flank = 800;
}

my $beg0;
my $end0;
my $lenextend;

if ($fiv >= $flank && $thr >= $flank) {
$beg0 = $beg - $flank;
$end0 = $end + $flank;
$lenextend = 2*$flank;
} elsif ($fiv >=$flank && $thr < $flank) {
$beg0 = $beg - $flank;
$end0 = $end + $thr - 5;
$lenextend = $flank + $thr - 5;
} elsif ($fiv < $flank && $thr >= $flank) {
$beg0 = $beg - $fiv + 5;
$end0 = $end + $flank;
$lenextend = $flank + $fiv -5;
} else {
$beg0 = $beg - $fiv + 5;
$end0 = $end + $thr - 5;
$lenextend = $fiv+$thr-10;
}

eval {  ## tag0 begin
####################### tag1 begin
my $altseq;
my ($gap_s,$gap_e);
$altseq = &SubSeq("$queries/$alt.$alt_chr.fa",$qchr,$qs-$flank,$qe+$flank,"plus");
$gap_s = $flank;
$gap_e = $flank;
####################### tag1 end

####################### tag2 begin
my $k;

if (!exists ${$liftover_coord}{$beg0}) {
       for ($k=$beg0;$k<=$beg;$k++) {
        if (exists ${$liftover_coord}{$k}) {
           $beg0 =$k;
           last;
         }
       }
    }

if (!exists ${$liftover_coord}{$end0}) {
       for ($k=$end0;$k>=$end;$k--) {
        if (exists ${$liftover_coord}{$k}) {
           $end0 =$k;
           last;
         }
      }
   }

if (exists ${$liftover_coord}{$beg0} && exists ${$liftover_coord}{$end0}) {
    my @begflank1= split (/_/,${$liftover_coord}{$beg0});
    my @endflank1 = split (/_/,${$liftover_coord}{$end0});

    if ( abs ($begflank1[1]-$endflank1[1]) - $lenextend < 13 ) {
      print VCF "0|0:22\t$line\n";
      next;
    } elsif (abs ($begflank1[1]-$endflank1[1]) - $lenextend > 200000) {
      print VCF "\.|\.:22\t$line\n";
      next;
    }

    my $indseq;
    if ($begflank1[1]<$endflank1[1] && ($begflank1[0] eq $endflank1[0] )) {
    $indseq = &SubSeq("$queries/$ind.$scf.fa",$begflank1[0],$begflank1[1],$endflank1[1],"plus");
    } elsif ($begflank1[1]>$endflank1[1] && ($begflank1[0] eq $endflank1[0])) {
    $indseq = &SubSeq("$queries/$ind/$ind.$scf.fa",$begflank1[0],$endflank1[1],$begflank1[1],"minus");
    } else {
     print VCF "\.|\.:22\t$line\n";
     next;
    }

   my ($seq1,$seq2,$e_value) = &Realign("$alt.ALT.$scf","$ind.VAR.$scf",$altseq,$indseq);
   if (defined ($seq1) && defined ($seq2) && defined ($e_value) ) {
                    if ($e_value>95) {
                             print VCF "1|1:22\t$line\n";
                     } else {
                               my $evaluation = &ParseMaf($seq1,$seq2,$gap_s,$gap_e);
                               if (defined ($evaluation) ) {                                    
                                            print VCF "$evaluation|$evaluation:22\t$line\n";                                     
                               } else {         
                                            print VCF "\.|\.:22\t$line\n";
                                            next;                           
                               }                          
                     }
    } else {
             print VCF "\.|\.:22\t$line\n";
             next;       
    }
  
} else {
   print VCF "\.|\.:22\t$line\n";
   next;
}
####################### tag2 end
##}## tag0 end
}; ## eval
 
if ($@) {
  print ERR "Failed\t$line\n";
}
     
    } ## while 
    $pm->finish;
  } ## chr
    $pm->wait_all_children;
} ## individuals


###############################################################################################
##############################################################################################################
##############################################################################################################

###########################  Start to write to VCF file ######################################################
my %genotypinginf;
opendir (DIR3,"$OUTPUTdir") || die "Error in opening temp vcf directory\n";
while ((my $filename = readdir (DIR3))) {
     if ($filename =~/$refname\.(.*)\.(.*)\.vcf/)  {
       open FILE, "$OUTPUTdir/$filename" or die "$!";
       while (<FILE>) {
            chomp;
            my @string = split (/\t/,$_);
            my $svele = join ("_",($string[1],$string[2],$string[3]));
               $genotypinginf{$1}{$svele} = $string[0];
       }
     }
}

open SVtotal, "$svraw" or die "$!";
open VCF, ">$svraw.vcf" or die "$!";
##############################################################################################################

my $time= strftime "%F", localtime;
print VCF "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=YiLiao_PopASSYSV_version20200901\n\#\#reference=$refname\n";
open Contig, "$refsize" or die "$!";
while (<Contig>) {
chomp;
my @tmp = split(/\s/,$_);
print VCF "\#\#contig=<ID=$tmp[0],length=$tmp[1],assembly=$refname,md5=XXXXXXXX,species=\"$speciesname\",taxonomy=x>\n";
}
close Contig;
print VCF "##INFO=<ID=BKPTID,Number=.,Type=String,Description=\"ID of the assembled alternate allele in the assembly file\">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical micro-homology at event breakpoints\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=CNV,Description=\"Copy number variable region\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Copy number genotype for imprecise events\">
##FORMAT=<ID=CNQ,Number=1,Type=Float,Description=\"Copy number genotype quality for imprecise events\">
";
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
my $n;
for ($n=0;$n<$#Genome;$n++) {
print VCF "$Genome[$n]\t";
}
print VCF "$Genome[-1]\n";

##############################################################################################################
while (<SVtotal>) {
            
            next if ($_=~/\#CHROM/);
            chomp;
            my @string1 = split (/\t/,$_);
            my $svele1 = join ("_",($string1[0],$string1[1],$string1[2]));
            
            print VCF "$string1[0]\t$string1[1]\tSV$string1[1]\t$string1[13]\t$string1[14]\t22\tPASS\tSVTYPE=$string1[4];END=$string1[2];SVLEN=$string1[5]\tGT:GQ\t";
            
            my $m;
            for ($m=0;$m<$#Genome;$m++) {
                       if (exists $genotypinginf{$Genome[$m]}{$svele1}) {
                            print VCF "$genotypinginf{$Genome[$m]}{$svele1}\t";
                       } else {
                            print VCF "\.|\.:22\t";
                       }
            }
            
            if (exists $genotypinginf{$Genome[-1]}{$svele1}) {
                            print VCF "$genotypinginf{$Genome[-1]}{$svele1}\n";
            } else {
                            print VCF "\.|\.:22\n";
            }
   
}

##############################################################################################################
##############################################################################################################

#################
##### Module ####
#################
sub ParseMaf {
my ($sequenceA,$sequenceB,$svbeg,$svend) = @_;
my @seqA = split(//,$sequenceA);
my @seqB = split(//,$sequenceB);

my $i;
my $len_n=0;
my @gap_coord;
my $start;
my $end;
for ($i=0;$i<=$#seqA;$i++) {
 if ($seqA[$i]=~/[ATGCNatgcn]/) {
  $len_n++;
      if ($len_n==$svbeg) {
          $start = $i;
      }
      if ($len_n==$svend) {
          $end = $i;
      }
   }
}

my $j;
my $m=0;
my $n=0;
my $h=0;
my $evalue;

if (defined $start) {
for ($j=0;$j<$start;$j++) {
     if ($seqA[$j] =~/[ATGCNatgcn]/) {
         $n++;
        $seqA[$j] = uc ($seqA[$j]);
       if ($seqB[$j] =~/[ATGCNatgcn]/) {
           $h++;
           $seqB[$j] = uc ($seqB[$j]);
              if ($seqA[$j] eq $seqB[$j]) {
                 $m++;
              }
       }
    }
}

my $gap_identity5 = $m/($h+1);
my $gap_overlap5 = $h/($n+1);

my $q;
my $n1=0;
my $m1=0;
my $h1=0;
for ($q=$end+1;$q<=$#seqA;$q++) {
     if ($seqA[$q] =~/[ATGCNatgcn]/) {
         $n1++;
        $seqA[$q] = uc ($seqA[$q]);
       if ($seqB[$q] =~/[ATGCNatgcn]/) {
              $h1++;
              $seqB[$q] = uc ($seqB[$q]);
              if ($seqA[$q] eq $seqB[$q]) {
                 $m1++;
              }
           }
       }
   }
my $gap_identity3 = $m1/($h1+1);
my $gap_overlap3 = $h1/($n1+1);

my $n2=0;
my $m2=0;
my $y;
my $h2;
 for ($y=$start;$y<=$end;$y++) {
  if ($seqA[$y] =~/[ATGCNatgcn]/) {
    $n2++;
    $seqA[$y] = uc ($seqA[$y]);
    if ($seqB[$y] =~/[ATGCNatgcn]/) {
       $h2++;
       $seqB[$y] = uc ($seqB[$y]);
       if ($seqA[$y] eq $seqB[$y]) {
                    $m2++;
        }
     }
   }
 }
my $gap_region = $m2/($n2+1);
   #my $evalue;
   if ($gap_identity5>=0.75 && $gap_identity3 >=0.75 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region >= 0.4) {
       $evalue=1; ### similar to reference
    } elsif ($gap_identity5>=0.75 && $gap_identity3 >=0.75 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region < 0.1) {
       $evalue=0; ### support a deletion
    } elsif ($gap_region >= 0.9) {
      $evalue=1;
    } else {
      $evalue="\."; ### missing data
    }
} else {
     $evalue="\."; ### missing data
}
return $evalue;
}

#############################################

sub stretcher2maf {
my ($aname,$aseq,$alen,$bname,$bseq,$blen) = @_;
open OUT, ">$OUTPUTdir/$aname.$alen.$bname.$blen.maf" or die "$!";
print OUT "a score=1000000.0\n";
print OUT "s $aname 0 $alen + $alen $aseq\n";
print OUT "s $bname 0 $blen + $blen $bseq\n";
}

##############################################

sub Realign {
my ($t,$q,$seq1,$seq2)=@_;
open SEQ1, ">$OUTPUTdir/$q/$t.fa" or die "$!";
print SEQ1 ">$t\n$seq1\n";
open SEQ2, ">$OUTPUTdir/$q/$q.fa" or die "$!";
print SEQ2 ">$q\n$seq2\n";
`${stretcher}stretcher -asequence $OUTPUTdir/$q/$t.fa -bsequence $OUTPUTdir/$q/$q.fa -outfile $OUTPUTdir/$q/$t.$q.stretcher`;
my ($s1, $s2) = &ParseStretcherOutPut("$OUTPUTdir/$q/$t.$q.stretcher",$t,$q);
my $eValue = &Parse_Stretcher("$OUTPUTdir/$q/$t.$q.stretcher");
`rm $OUTPUTdir/$q/*.fa $OUTPUTdir/$q/*.stretcher`;
return ($s1,$s2,$eValue);
}
##############################################

sub ParseStretcherOutPut {
my ($stretcher,$t,$q) = @_;
local $/="\n";
open EMBOSS, "$stretcher" or die "$!";
my @ref;
my @query;
while (<EMBOSS>) {
chomp;
$_ =~ s/^\s+//;
if ($_=~/\#/) {
}elsif ($_=~/^$t/) {
my @tmp1 = split (/\s/,$_);
push (@ref,$tmp1[1]);
} elsif ($_=~/^$q/) {
my @tmp2 = split (/\s/,$_);
push (@query,$tmp2[1]);
 }
 }
 my $seq1 = join ("",@ref);
 my $seq2 = join ("",@query);
 return ($seq1,$seq2);
 }
#################################################
sub Parse_Stretcher {
my $stretcher = shift;
my $value;
local $/="\n";
open IN2, "$stretcher" or die "$!";
while (<IN2>) {
  if ($_=~/\# Identity:/) {
    $_=~/\# Identity:   (.*) \((.*)%\)/;
    $value = $2;
    last;
  }
}
return $value;
}
################################################## Extract sequence based on coordinate
sub SubSeq {
my ($seq,$chr,$beg,$end,$strand) = @_;
local $/="\n>";
open FA, "$seq" or die "$!";
my $seqrange;
   while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
   if ($id eq $chr) {
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
         my $len = $end - $beg + 1;
         $seqrange = substr $sequence, $beg -1, $len;
         last;
        }
    }

my $sequence;
if ($strand eq "plus") {
$sequence = $seqrange;
} elsif ($strand eq "minus") {
$sequence = &RevCom ($seqrange);
}
return $sequence;
}
################################################# Reverse and complement the sequence
sub RevCom {
my $seq = shift;
my $revcomp = reverse $seq;
$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;
return $revcomp;
}
#################################################
sub Maf2fasta {
my($maf,$ref,$qs,$qe,$extend)=@_;
open MAF, "$maf" or die "Missing maf file!";
$extend ||=0;
my %coord_proxy;
my $beg_index=0;
my $end_index=0;
my $s=1;
my $n=0;
my $ref_end=0;
my $gap_ref = &GetRefCoord("$maf",$ref);
my @species;
$ref = (split(/\./,$ref))[0];
local $/="\n\n"; # record separator as a blank line when reading
my $alt_strand;
my $alt_size;
while (<MAF>) {
chomp;
next if ($_=~/eof/);
my  @temp = split(/\n/,$_);
$s = $s + $ref_end + ${$gap_ref}[$n];
foreach my $line (@temp) {
  if ($line =~/^s/) {
    my ($spec,$start,$extend1,$strand,$qsize,$seq) = (split(/\s+/,$line))[1,2,3,4,5,6];
    push (@species,$spec);
    $spec = (split(/\./,$spec))[0];
    if ($spec !~ /$ref/) {
    $alt_strand = $strand;
    $alt_size = $qsize;
    }
    my $link =&links($s,$start,$strand,$qsize,$seq);
    if (!keys %{$coord_proxy{$spec}}) {
    $coord_proxy{$spec}=$link;
    } else {
    my %linkcombine=(%{$coord_proxy{$spec}},%{$link});
    $coord_proxy{$spec}=\%linkcombine;
     }
    }
   }
  $n++;
  my @array = sort {$a<=>$b} (keys %{$coord_proxy{$ref}});
  $ref_end = $array[-1];
  }

  my @Arr = sort {$a<=>$b} (keys %{$coord_proxy{$ref}});

  if ($alt_strand =~/\+/) {
  my $begin = $qs-1 - $extend;
  my $end2 = $qe + $extend;
  my @all_matches = grep {${$coord_proxy{$ref}}{$_} eq $begin} keys %{$coord_proxy{$ref}};
  if ($#all_matches>-1) {
   $beg_index = $all_matches[0];
  } else {
   $beg_index = $Arr[0];
  }

  my @all_matches1 = grep {${$coord_proxy{$ref}}{$_} eq $end2} keys %{$coord_proxy{$ref}};
  if ($#all_matches1>-1){
    $end_index = $all_matches1[0];
  } else {
    $end_index = $Arr[-1];
  }
} elsif ($alt_strand =~/\-/) {
  my $begin = $qs-$extend;
  my $end2 = $qe+1+$extend;


my @all_matches = grep {${$coord_proxy{$ref}}{$_} eq $begin} keys %{$coord_proxy{$ref}};
 if ($#all_matches>-1) {
 $beg_index = $all_matches[0];
} else {
  $beg_index = $Arr[0];
}

my @all_matches1 = grep {${$coord_proxy{$ref}}{$_} eq $end2} keys %{$coord_proxy{$ref}};
if ($#all_matches1>-1) {
$end_index = $all_matches1[0];
} else {
  $end_index = $Arr[-1];
  }
}
  if ($alt_strand =~/\+/) {
  $alt_strand = "plus";
  } elsif ($alt_strand=~/\-/) {
  $alt_strand = "minus";
  }
  my $ref_hash = \%coord_proxy;
  my $ref_spec = \@species;
  return ($ref_hash,$beg_index,$end_index,$alt_strand,$ref_spec,$alt_size);
}

###############################################################
sub GetRefCoord {
my ($maf,$ref) = @_;
my @gaps;
my $end;
local $/="\n";
open In, "$maf" or die "$!";
  while (<In>) {
  chomp;
  if ($_=~/$ref/ and $_=~/^s/) {
    my ($start,$extend) = (split(/\s+/,$_))[2,3];
    if ($end) {
     my $gap = $start - $end;
     push (@gaps,$gap);
   } else {
     push (@gaps,"0");
   }
    $end = $start + $extend;
   }
  }
  my $gapref = \@gaps;
  return $gapref;
}

###############################################################
sub links {
my ($s,$beg,$strand,$qsize,$seq) = @_;
my %hash;
my @tmp = split(//,$seq);
my $i;
my $j=0;
 for ($i=0;$i<=$#tmp;$i++) {
    if ($strand =~/\+/) {
        if ($tmp[$i]=~/[ATGCNatgcn]/) {
            $hash{$i+$s} = $j + $beg;
            $j++;
        }  else {
            $hash{$i+$s} = $j + $beg;
        }
    } elsif ($strand =~/\-/) {
        if ($tmp[$i]=~/[ATGCNatgcn]/) {
          $hash{$i+$s} = $qsize - $beg - $j;
          $j++;
        }  else {
          $hash{$i+$s} = $qsize - $beg - $j;
        }
    }
 }
my $reflinks = \%hash;
return $reflinks;
}
#################################################################
sub liftOver {
my $maf = shift;
my %liftover;
local $/="\n\n";
open MAF, "<$maf" or die "$!";
my @first = split (/\n/,<MAF>);
my @array;
foreach my $ele (@first) {
  if ($ele=~/^s/) {
   push (@array,$ele);
  }
}

my ($target,$Ts,$Textend,$Tstrand,$Tsize,$Tseq) = (split(/\s+/,$array[0]))[1,2,3,4,5,6];
my ($query,$Qs,$Qextend,$Qstrand,$Qsize,$Qseq) = (split(/\s+/,$array[1]))[1,2,3,4,5,6];
my @targetSEQ = split(//,$Tseq);
my @querySEQ = split (//,$Qseq);
my $i;
my $n=0;
my $m=0;

for ($i=0;$i<=$#targetSEQ;$i++) {
   if ($targetSEQ[$i]!~/\-/) {
     $n++;
   }
   if ($querySEQ[$i]!~/\-/) {
     $m++;
   }

   if ($targetSEQ[$i]!~/\-/ && $querySEQ[$i]!~/\-/) {
     if ($Qstrand =~ /\+/) {
       my $tcoord =  $Ts + $n + 1;
       my $qcoord =  $Qs + $m + 1;
       my $qchrCoord = join  ("_",($query,$qcoord));
       $liftover{$tcoord} = $qchrCoord;
     } elsif ($Qstrand =~/\-/) {
       my $tcoord =  $Ts + $n + 1;
       my $qcoord =  $Qsize - ($Qs + $m);
       my $qchrCoord = join  ("_",($query,$qcoord));
       $liftover{$tcoord} = $qchrCoord;
     }
   }
}

while (<MAF>) {
next if ($_=~/eof/);
my @tmp = split (/\n/,$_);
my ($target1,$Ts1,$Textend1,$Tstrand1,$Tsize1,$Tseq1) = (split(/\s+/,$tmp[1]))[1,2,3,4,5,6];
my ($query1,$Qs1,$Qextend1,$Qstrand1,$Qsize1,$Qseq1) = (split(/\s+/,$tmp[2]))[1,2,3,4,5,6];
my @targetSEQ1 = split(//,$Tseq1);
my @querySEQ1 = split (//,$Qseq1);
my $i1;
my $n1=0;
my $m1=0;
for ($i1=0;$i1<=$#targetSEQ1;$i1++) {
   if ($targetSEQ1[$i1]!~/\-/) {
     $n1++;
   }
   if ($querySEQ1[$i1]!~/\-/) {
     $m1++;
   }
   if ($targetSEQ1[$i1]!~/\-/ && $querySEQ1[$i1]!~/\-/) {
     if ($Qstrand1 =~ /\+/) {
       my $tcoord1 =  $Ts1 + $n1 + 1;
       my $qcoord1 =  $Qs1 + $m1 + 1;
       my $qchrCoord1 = join  ("_",($query1,$qcoord1));
       $liftover{$tcoord1} = $qchrCoord1;;
     } elsif ($Qstrand1 =~/\-/) {
       my $tcoord1 =  $Ts1 + $n1 + 1;
       my $qcoord1 =  $Qsize1 - ($Qs1 + $m1);
       my $qchrCoord1 = join  ("_",($query1,$qcoord1));
       $liftover{$tcoord1} = $qchrCoord1;
     }
   }
 }
}
my $refhash = \%liftover;
return $refhash;
}
#################################################################



