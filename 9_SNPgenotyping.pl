#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use lib "$FindBin::Bin/./lib";
use Parallel::ForkManager;
use Pod::Usage;

####################################################################################
# GenometypingSNPs.pl
#
# Written by Yi Liao (12/25/2022)
#
# Genometype SNPs based on sing-coverage maf files
#
####################################################################################
my ($singmaf,$UCSC,$TBA,$refseq,$refsize,$refname,$speciesname,$OUTPUTdir,$p,$help);

GetOptions(
    "singmaf=s"=>\$singmaf, # Directory for single coverage maf files
    "tba=s"=>\$TBA,
    "ucsc=s"=>\$UCSC, # Path to Kent's utilities
    "refsize=s"=>\$refsize,
    "refseq=s"=>\$refseq,
    "refname=s"=>\$refname,
    "speciesname=s"=>\$speciesname,
    "outdir=s"=>\$OUTPUTdir,
    "t=i"=>\$p,
    "help"=>\$help,
);


#########################################
############# Help information ##########
if ($help){
print <<"END.";
  Usage : perl $0 --refseq ./Ref/MH63 --refsize ./Ref/MH63.sizes --refname MH63 --singmaf ../sing --t 8 --outdir SNPVCF
 
  --refseq       reference sequence                                      [REQUIRED]
  --refsize      reference sizes                                         [REQUIRED]
  --refname      reference name, e.g. hg38, IRGSP1                       [REQUIRED]
  --singmaf      path to single coverage maf file                        [REQUIRED]
  --t            how many processors will invoke                         [default : 4]
  --speciesname  species name, e.g. human, Drosophila, rice , maize      [default: Unknown]
  --outdir       define an output dir                                    [default: VCFtemp]
  --help         print this help information
END.
exit;
}

my $message = "Try 'perl $0 --help' for more information.";

if (! $refseq){
    pod2usage( {
           -message => $message,
           -verbose => 0,
           -exitval => 1 } );
    }

##############################
$UCSC ||="$Bin/pub/UCSC/";
$TBA ||="$Bin/pub/TBA/";
$OUTPUTdir ||= "VCFtemp";
$speciesname ||="Undefined";
$p ||= 4;
mkdir($OUTPUTdir) unless(-d $OUTPUTdir);

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

print "Chr: @chromosomes\n";
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

##################################
###### Start to genotype #########
##################################
foreach my $ind (@Genome) { 
   my $pm = Parallel::ForkManager->new($p);
DATA_LOOP:
  foreach my $scf (@chromosomes) { 
     print "Process 1: Currently are processing chromosome $scf of sample $ind !!!\n";
     my $pid = $pm->start and next DATA_LOOP; 
     my $liftover = &liftOver("$singmaf/$ind/$refname\_$ind\_$scf.maf");
     open VCF, ">$OUTPUTdir/$refname.$ind.$scf.vcf" or die "$!";
     foreach my $key (keys %{$liftover}) {
     print VCF "$key\t${$liftover}{$key}\n";
     } 
     $pm->finish;
   }
     $pm->wait_all_children; 
}

###################################
foreach my $chr (@chromosomes) {
  my %subvcf;
 
  opendir (DIR3,"$OUTPUTdir") || die "Error in opening dir sub vcf directory\n"; 
  while ((my $filename2 = readdir (DIR3))) {
    
    if ($filename2 =~ /$refname.(.*).$chr.vcf/) {
    
    open subVCF, "$OUTPUTdir/$filename2" or die "$!";
     while (<subVCF>) {
       chomp;
       my @line = split (/\t/,$_);
       
       if (exists $subvcf{$line[0]}) {
          push (@{$subvcf{$line[0]}},$line[1]) unless (grep { $line[1] eq $_ } @{$subvcf{$line[0]}});
       } else {
          push (@{$subvcf{$line[0]}},$line[1])
       }
       
     }
    close subVCF;
    
    } 
  }

  
  open VCFall, ">$OUTPUTdir/$refname.$chr.all.vcf" or die "$!";  
  foreach my $locus (keys %subvcf) {
    my $eles = join(",",@{$subvcf{$locus}});
    print VCFall "$locus\t$eles\n";
  }
  close VCFall;  

}

###################################

foreach my $ind1 (@Genome) { 
   my $pm1 = Parallel::ForkManager->new($p);
DATA_LOOP1:
  foreach my $scf1 (@chromosomes) { 
     print "Process 2: Currently are processing chromosome $scf1 of sample $ind1 !!!\n";
     my $pid1 = $pm1->start and next DATA_LOOP1; 
     my $liftover1 = &ReliftOver("$singmaf/$ind1/$refname\_$ind1\_$scf1.maf");
     open SNP, "$OUTPUTdir/$refname.$scf1.all.vcf" or die "$!";
     open VCF1, ">$OUTPUTdir/$refname.$ind1.$scf1.SNP.vcf" or die "$!";
     while (<SNP>) {
     chomp;
     my @tmp = split (/\t/,$_);
     if (exists ${$liftover1}{$tmp[0]}) {
     print VCF1 "$_\t${$liftover1}{$tmp[0]}\n";
     } else {
     print VCF1 "$_\t\-\n";
     }
     }  
     $pm1->finish;
   }
     $pm1->wait_all_children; 
}

#########################################################################################################

my $time= strftime "%F", localtime;
open VCF0,">$OUTPUTdir/$refname.SNP.all.vcf" or die "$!";
print VCF0 "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=SV-GAPs\n\#\#reference=$refname\n";
open Contig, "$refsize" or die "$!";
while (<Contig>) {
chomp;
my @tmp = split(/\s/,$_);
print VCF0 "\#\#contig=<ID=$tmp[0],length=$tmp[1],assembly=$refname,md5=XXXXXXXX,species=\"$speciesname\",taxonomy=x>\n";
}
close Contig;
print VCF0 "##phasing=phased
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
";
print VCF0 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";


my $num;
for ($num=0;$num<$#Genome;$num++) {
print VCF0 "$Genome[$num]\t";
}
print VCF0 "$Genome[-1]\n";

##########

my $i;
for ($i=0;$i<=$#chromosomes;$i++) {

my %genotypinginf;
opendir (DIR5,"$OUTPUTdir") || die "Error in opening temp vcf directory\n";
while ((my $filename = readdir (DIR5))) {
     if ($filename =~/$refname\.(.*)\.$chromosomes[$i]\.SNP\.vcf/)  {
       open FILE, "$OUTPUTdir/$filename" or die "$!";
       while (<FILE>) {
            chomp;
            my @string = split (/\t/,$_);
            $genotypinginf{$1}{$string[0]} = $string[2];
       }
     }
}

open COM, "$OUTPUTdir/$refname.$chromosomes[$i].all.vcf" or die "$!";

while (<COM>) {
 chomp;
 my @ln = split (/\t/,$_);
 my @ln1 = split (/_/,$ln[0]);
 my @Base;
    
   foreach my $ind (@Genome) {
     if (exists ${$genotypinginf{$ind}}{$ln[0]}) {
            push (@Base,${$genotypinginf{$ind}}{$ln[0]});
     } else {
            push (@Base, "\-");
      }
    }
    
    my @alt = split (/,/,$ln[1]);
    my $Qual = 50;
    my $Filter = "PASS";
    my $NS = 0;

    foreach my $ind (@Base) {
       if ($ind ne "\-") {
           $NS++;
       }
    }

     my $DP = 20;
     my @AF;
     foreach my $ele (@alt) {
       my $num = 0;
          foreach my $ele1 (@Base) {
             if ($ele1 eq $ele) {
                 $num++;
             }
          }
     my $fre = $num/$NS;
     my @tmp = split (//,$fre);
     if ($#tmp>3) {
     $fre = substr($fre,0,4);
     }
     push (@AF,$fre);
     }
     
#     my $eles = join(",",@{$subvcf{$locus}});
     my $alt_array = join (",",@alt);
     my $fre_array = join (",",@AF);
     print VCF0 "$ln1[0]\t$ln1[1]\tSNP$ln1[1]\t$ln1[2]\t$alt_array\t$Qual\t$Filter\tNS=$NS;DP=$DP;AF=$fre_array\tGT:GQ:DP:HQ\t";
     
     my $n;
     for ($n=0;$n<=$#Base;$n++) {
      
     if ($Base[$n] eq $ln1[2]) {
        if ($n==$#Base) {
        print VCF0 "0|0:48:10:51,51\n";
        } else {
        print VCF0 "0|0:48:10:51,51\t";
        }
        
     } elsif ($Base[$n] eq "-") {
        if ($n==$#Base) {
        print VCF0 "\.|\.:48:10:51,51\n";
        } else {
        print VCF0 "\.|\.:48:10:51,51\t";
        }
        
     } else {

       my $m;
       for ($m=0;$m<=$#alt;$m++){
       my $t = $m + 1;
          if ($Base[$n] eq $alt[$m]) {
            if ($n==$#Base) {
                print VCF0 "$t|$t:48:10:51,51\n";
            } else {
                print VCF0 "$t|$t:48:10:51,51\t";
            }
         }
       } 
       
    }   
 }   
}
}

`cat $OUTPUTdir/$refname.SNP.all.vcf | grep "#" > $OUTPUTdir/head.txt`;
`cat $OUTPUTdir/$refname.SNP.all.vcf | grep -v "#" | sort -k1,1 -k2,2n > $OUTPUTdir/$refname.SNP.all.tmp.vcf`;
#`cat $OUTPUTdir/head.txt $OUTPUTdir/$refname.SNP.all.tmp.vcf > $OUTPUTdir/$refname.SNP.all.sorted.vcf`;
`cat $OUTPUTdir/head.txt $OUTPUTdir/$refname.SNP.all.tmp.vcf > $OUTPUTdir/All.sorted.SNP.vcf`;
`rm $OUTPUTdir/$refname.*.vcf $OUTPUTdir/head.txt`;

#########################################################################################################

###################################
###################################
###################################
###################################
###################################

sub liftOver {
my $maf = shift;
my %liftover;
local $/="\na score";
open MAF, "<$maf" or die "$!";

while (<MAF>) {
next if ($_=~/##maf version/);
my @tmp = split (/\n/,$_);
my ($target,$Ts,$Textend,$Tstrand,$Tsize,$Tseq) = (split(/\s+/,$tmp[1]))[1,2,3,4,5,6];
my ($query,$Qs,$Qextend,$Qstrand,$Qsize,$Qseq) = (split(/\s+/,$tmp[2]))[1,2,3,4,5,6];
$Tseq = uc ($Tseq);
$Qseq = uc ($Qseq);
my @targetSEQ = split(//,$Tseq);
my @querySEQ = split (//,$Qseq);

my $j;
my $m=0;
for ($j=0;$j<=$#targetSEQ;$j++) {
  if ($targetSEQ[$j]!~/\-/) {
       $m++;
       my $tcoord =  $Ts + $m;
       if (($targetSEQ[$j] ne $querySEQ[$j]) and ($querySEQ[$j]!~/\-/)) {
       my $tinf = join ("_",($target,$tcoord,$targetSEQ[$j]));
       $liftover{$tinf} = $querySEQ[$j];
       }
  }
 }
}

my $refhash = \%liftover;
return $refhash;
}

#################################################################

sub ReliftOver {
my $maf = shift;
my %liftover;
local $/="\na score";
open MAF, "<$maf" or die "$!";

while (<MAF>) {
next if ($_=~/##maf version/);
my @tmp = split (/\n/,$_);
my ($target,$Ts,$Textend,$Tstrand,$Tsize,$Tseq) = (split(/\s+/,$tmp[1]))[1,2,3,4,5,6];
my ($query,$Qs,$Qextend,$Qstrand,$Qsize,$Qseq) = (split(/\s+/,$tmp[2]))[1,2,3,4,5,6];
$Tseq = uc ($Tseq);
$Qseq = uc ($Qseq);
my @targetSEQ = split(//,$Tseq);
my @querySEQ = split (//,$Qseq);

my $j;
my $m=0;
for ($j=0;$j<=$#targetSEQ;$j++) {
  if ($targetSEQ[$j]!~/\-/) {
       $m++;
       my $tcoord =  $Ts + $m;
       if ($querySEQ[$j]!~/\-/) {
       my $tinf = join ("_",($target,$tcoord,$targetSEQ[$j]));
       $liftover{$tinf} = $querySEQ[$j];
       }
  }
 }
}

my $refhash = \%liftover;
return $refhash;
}

#################################################################
