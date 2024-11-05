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
# Written by Yi Liao (12/31/2023)
#
# Usage： Genotyping inversions for population-scale genome assemblies
#
# Dependencies: UCSC's utilities; TBA; Perl modules, Parallel::ForkManager
#
####################################################################################
$/="\n";
my ($svraw,$UCSC,$TBA,$singmaf,$refsize,$refname,$speciesname,$OUTPUTdir,$help);

GetOptions(
    "inv=s"=>\$svraw,
    "refsize=s"=>\$refsize,
    "sing=s"=>\$singmaf,
    "tba=s"=>\$TBA,
    "ucsc=s"=>\$UCSC,
    "refname=s"=>\$refname,  
    "speciesname=s"=>\$speciesname,
    "outdir=s"=>\$OUTPUTdir,
    "help"=>\$help,
);

############# Help information ##########
if ($help){
print <<"END.";
Usage: Perl $0 --ins INV.combined.txt -sing /path/to/sing/ --refsize MH63.sizes --refname MH63 --t 12

--inv          the file of the merged non-redundant insertions (> 49 bp)               [REQUIRED]
--query        folder kept all query genome sequences                                  [REQUIRED]
--singmaf      folder kept single-coverage maf files (on the reference side)           [REQUIRED]
--refsize      chromosome/contig sizes from the reference                              [REQUIRED]
--refseq       reference genome sequence                                               [REQUIRED]
--stretcher    path to stretcher                                                       [REQUIRED]    
--refname      reference name, e.g. hg38, IRGSP1 [REQUIRED]
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
#########################################
$UCSC ||="$Bin/pub/UCSC/";
$TBA ||="$Bin/pub/TBA/";
$speciesname ||= "Unknown";
$OUTPUTdir ||= "VCFtemp";
mkdir($OUTPUTdir) unless(-d $OUTPUTdir);
mkdir("$OUTPUTdir/err") unless(-d "$OUTPUTdir/err");
#########################################

################################## Split maf by chromosomes and obtain query genome names
opendir (DIR1,"$singmaf") || die "Error in opening dir singmaf\n";
my @Genome;
`touch splits.bed`;
  while ((my $filename1 = readdir (DIR1))) {
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
##################################################################################

foreach my $ind (@Genome) {
	         print "Currently are processing chromosome sample $ind !!!\n";
###############################################################################################
open IN, "$svraw" or die "$!";
open VCFtemp, ">$OUTPUTdir/$refname.$ind.vcf" or die "$!";
open ERR, ">$OUTPUTdir/err/$refname.$ind.err.vcf" or die "$!";

while (<IN>) {
chomp;
print "$ind\t$_\n";
next if ($_=~/\#/);
my %exi;
my $line = $_;
my ($chr,$beg,$end,$len,$sample,$qchr,$qs,$qe) = (split (/\t/,$_))[0,1,2,5,7,8,9,10];

my @sampes = split (/,/,$sample);
foreach my $insample (@sampes) {
  my @unitk = split (/\./,$insample);
  $exi{$unitk[0]} = 0;
}
if (exists $exi{$ind}) {
  print VCFtemp "1|1:22\t$line\n";
  next;
}

my @unitChr = split (/\./,$chr);
my $scf = $unitChr[1];
print "chromosome:$scf\n";

eval {  
my $gen = &CheckMAF_ForINV ($chr,$beg,$end,"$singmaf/$ind/$refname\_$ind\_$scf.maf");
print VCFtemp "$gen:22\t$line\n";
}; ## eval
 
if ($@) {
  print ERR "Failed\t$line\n";
}
     
    }  

} 

close IN;
close VCFtemp;
close ERR;


################################################################################
##############################################################################################################

###########################  Start to write to VCF file ######################################################
my %genotypinginf;
opendir (DIR2,"$OUTPUTdir") || die "Error in opening temp vcf directory\n";
while ((my $filename = readdir (DIR2))) {
     if ($filename =~/$refname\.(.*)\.vcf/)  {
       #print "reading $filename\t $1\n";
       open FILE, "$OUTPUTdir/$filename" or die "$!";
       #print "files: $OUTPUTdir/$filename\n";
       #print "reading $filename\t $1\n";

#       my $genty = &readfile ("$OUTPUTdir/$filename");
       
#       print "$filename\n";


      while (<FILE>) {
           # print "$_\n";
            chomp;
            my @string = split (/\t/,$_);
            my $svele = join ("_",($string[1],$string[2],$string[3]));
              # print "$filename\t$svele\n";
               $genotypinginf{$1}{$svele} = $string[0];
              # print "$filename\t$1\t$string[0]\n";
       }
     }
}


=pod
sub readfile {
my $file = shift;
my %genotypes;

open VCFtemp, "$file" or die "$!";
while (<VCFtemp>) {
chomp;
my @string = split (/\t/,$_);
my $svele = join ("_",($string[1],$string[2],$string[3]));
$genotypes{$1}{$svele} = $string[0];
print "$file\t$string[0]\n";
}

my $ref = \ %genotypes;
return $ref;
}
=cut

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
            
            print VCF "$string1[0]\t$string1[1]\tSV$string1[1]\tref\talt\t22\tPASS\tSVTYPE=$string1[4];END=$string1[2];SVLEN=$string1[5]\tGT:GQ\t";
            #print VCF "$string1[0]\t$string1[1]\tSV$string1[1]\t$string1[13]\t$string1[14]\t22\tPASS\tSVTYPE=$string1[4];END=$string1[2];SVLEN=$string1[5]\tGT:GQ\t";
            
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

########################################
########################################
########################################
########################################
sub CheckMAF_ForINV {
my ($target_Chr,$target_start,$target_end,$maf) = @_;
my $genotype; 
open INVMAF, "$maf" or die $!;
local $/="\n\n";
my $invlen = $target_end - $target_start;
my $extend = 0.5*$invlen;

my $five_extend = $target_start - $extend;
my $three_extend = $target_end + $extend;

my %fiveflank;
my %interval;
my %threeflank;

######################################
while (my $line = <INVMAF>) {
my @block = split (/\n/,$line);
my ($target,$Ts,$Textend,$Tstrand) = (split(/\s+/,$block[1]))[1,2,3,4];
my ($query,$Qs,$Qextend,$Qstrand) = (split(/\s+/,$block[2]))[1,2,3,4];

if (exists $fiveflank{$query}){
} else {
${$fiveflank{$query}}[0] = 0;
${$fiveflank{$query}}[1] = 0;
${$interval{$query}}[0] = 0;
${$interval{$query}}[1] = 0;
${$threeflank{$query}}[0] = 0;
${$threeflank{$query}}[1] = 0;
}

if ($target eq $target_Chr) {
    my $C22 = $Ts + $Textend;
    my $ovlap5 = &Overlap($Ts,$C22,$five_extend,$target_start);
    my $ovlap0 = &Overlap($Ts,$C22,$target_start,$target_end);
    my $ovlap3 = &Overlap($Ts,$C22,$target_end,$three_extend);
   
    if ($ovlap5 > 0) {

        print "Five flank: $ovlap5\t$query\t$Qstrand\n";
        if (exists $fiveflank{$query}) {

            if ($Qstrand eq "+") {
                ${$fiveflank{$query}}[0] = ${$fiveflank{$query}}[0] + $ovlap5;
            } elsif ($Qstrand eq "-") {
                ${$fiveflank{$query}}[1] = ${$fiveflank{$query}}[1] + $ovlap5;
            }
        } else {

             if ($Qstrand eq "+") {
                ${$fiveflank{$query}}[0] =  $ovlap5;
             } elsif ($Qstrand eq "-") {
                ${$fiveflank{$query}}[1] =  $ovlap5;
             }
                
       }
    }

   if ($ovlap0 > 0) {
        if (exists $interval{$query}) {
            if ($Qstrand eq "+") {
                ${$interval{$query}}[0] = ${$interval{$query}}[0] + $ovlap0;
            } elsif ($Qstrand eq "-") {
                ${$interval{$query}}[1] = ${$interval{$query}}[1] + $ovlap0;
            }
        } else {
             if ($Qstrand eq "+") {
                ${$interval{$query}}[0] =  $ovlap0;
             } elsif ($Qstrand eq "-") {
                ${$interval{$query}}[1] =  $ovlap0;
             }

       }
    }
    
   if ($ovlap3 > 0) {
        if (exists $threeflank{$query}) {
            if ($Qstrand eq "+") {
                ${$threeflank{$query}}[0] = ${$threeflank{$query}}[0] + $ovlap3;
            } elsif ($Qstrand eq "-") {
                ${$threeflank{$query}}[1] = ${$threeflank{$query}}[1] + $ovlap3;
            }
        } else {

             if ($Qstrand eq "+") {
                ${$threeflank{$query}}[0] =  $ovlap3;
             } elsif ($Qstrand eq "-") {
                ${$threeflank{$query}}[1] =  $ovlap3;
             }

       }
    }
}

}
######################################
my $max = 0;
my $q_chr;

foreach my $key (%interval) {
        my $num = ${$interval{$key}}[0] + ${$interval{$key}}[1];
        if ($num > $max) {
            $max = $num;
            $q_chr = $key;
        }
}

if ((exists $fiveflank{$q_chr}) and (exists $threeflank{$q_chr})) {
    my $fiverate = ${$fiveflank{$q_chr}}[1]/(${$fiveflank{$q_chr}}[0]+1);
    my $threerate = ${$threeflank{$q_chr}}[1]/(${$threeflank{$q_chr}}[0]+1);
    my $intervalrate = ${$interval{$q_chr}}[1]/(${$interval{$q_chr}}[0]+1);
       
    if ($fiverate > 5 && $threerate > 5 && $intervalrate > 5) {
        $genotype = "0|0";
    } elsif ($fiverate > 5 && $threerate > 5 && $intervalrate < 0.2) {
        $genotype = "1|1";
    } elsif ($fiverate < 0.2  && $threerate < 0.2 && $intervalrate < 0.2) {
        $genotype = "0|0";
    } elsif ($fiverate < 0.2 && $threerate < 0.2 && $intervalrate > 5) {
        $genotype = "1|1";
    } else {
        $genotype = ".|.";
    }
} else {
    $genotype = ".|.";
}
##############
return $genotype;
}

########################################################################
sub Overlap {
my ($c1,$c2,$c3,$c4) = @_;
my $overlap=0;

if ($c1>$c4) {
  $overlap =0;
} elsif ($c1>=$c3 and $c1<=$c4 and $c2 <=$c4) {
  $overlap = ($c2-$c1);
} elsif ($c1>=$c3 and $c1<=$c4 and $c2 >$c4) {
  $overlap = ($c4-$c1);
} elsif ($c1<=$c3 and $c2>=$c3 and $c2 <=$c4) {
  $overlap = ($c2-$c3);
} elsif ($c1<=$c3 and $c2>=$c4) {
  $overlap = ($c4-$c3);
}
return $overlap;
}
