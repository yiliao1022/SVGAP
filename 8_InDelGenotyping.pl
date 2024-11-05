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

####################################################################################
#
# Written by Yi Liao (01/02/2023)
#
# Further computationally validated the merge InDels
#
####################################################################################
my ($svraw,$svtype,$refseq,$queries,$singmaf,$UCSC,$refsize,$refname,$speciesname,$OUTPUTdir,$p,$head,$help);

GetOptions(
    "indel=s"=>\$svraw,
    "indeltype=s"=>\$svtype,
    "refseq=s"=>\$refseq,
    "refsize=s"=>\$refsize,
    "refname=s"=>\$refname,   
    "query=s"=>\$queries,
    "singmaf=s"=>\$singmaf,
    "ucsc=s"=>\$UCSC,
    "speciesname=s"=>\$speciesname,
    "t=s"=>\$p,
    'outdir=s'=>\$OUTPUTdir,
    'head'=>\$head,
    'help'=>\$help,
);

############# Help information ##########
if ($help) {
print <<"END.";
  Usage: perl $0 --indel SV.txt --indeltype INS --refseq ./genome/B73 --refsize ./genome/MH63.sizes --refname MH63 --query ./genome/query/ --singmaf /path/to/sing/ --outdir MH63_tmp --t 12
  --indel        raw merged InDels file                                  [Required]
  --indeltype    InDel types [ins|del]                                   [Required]
  --singmaf      path to pairwise genome algnment maf files              [Required]
  --refseq       reference genome sequence                               [Required]        
  --refsize      reference chromosome sizes file                         [Required]
  --refname      reference name, e.g. hg38, IRGSP1                       [Required]
  --query        folder kept all query genomes                           [Required]
  --t            how many threads to use [default:4]                     [Optional]
  --speciesname  species name, e.g. human,Drosophila,rice,maize          [Optional]
  --outdir       output Dir                                              [Optional]
  --head         print VCF header if invoked                             [Optional]
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

############################# Define default parameters
$UCSC ||="$Bin/pub/UCSC/";
#$OUTPUTdir ||= ".";
#$OUTPUTdir =~ s/\/$//;
$speciesname ||= "Unknown";
$OUTPUTdir ||= "VCFtemp";
$p ||= 4;
mkdir($OUTPUTdir) unless(-d $OUTPUTdir);
mkdir("$OUTPUTdir/err") unless(-d "$OUTPUTdir/err");

################################## Split Reference by chromosomes and obtain chromosome names
my $ref_file = basename($refseq);
my $ref_path = dirname($refseq);
my $byChr = "byChr";

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
open ($file, '>>', "$sv_path/$refname.$ctg.$svtype.bed") or die "Could not open file '$file' $!";
print $file "$_\n";
}
}
}

#######################################################
###################### Main ###########################
#######################################################
foreach my $ind (@Genome) {
        my $pm = Parallel::ForkManager->new($p);
DATA_LOOP:
       	foreach my $scf (@chromosomes) {
           my $pid = $pm->start and next DATA_LOOP;
	         mkdir("$OUTPUTdir/$ind.VAR.$scf");
	         print "Currently are processing chromosome $scf of sample $ind !!!\n";
	         my $liftover_coord = &liftOver("$singmaf/$ind/$refname\_$ind\_$scf.maf");
################################################################################################
open IN, "$sv_path/$refname.$scf.$svtype.bed" or die "$!";
open VCF, ">$OUTPUTdir/$refname.$ind.$scf.vcf" or die "$!";
open ERR, ">$OUTPUTdir/err/$refname.$ind.$scf.err.vcf" or die "$!";

while (<IN>) {
chomp;
my $line = $_;
my ($chr,$beg,$end,$len,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,8,9,10,11,12];
my $ref = (split(/\./,$chr))[0];
my $alt = (split(/\./,$qchr))[0];

eval {
####
if ($svtype eq "del") {
   my $beg0 = $beg + 1;
   my $end0 = $end + 2;
   ####
   if (exists ${$liftover_coord}{$beg} && exists ${$liftover_coord}{$end0}) {
          my @tmp1 = split (/_/, ${$liftover_coord}{$beg});
          my @tmp2 = split (/_/, ${$liftover_coord}{$end0});
          my $check = abs ($tmp2[1]-$tmp1[1]);

          if ($check == $len + 2) {
             if ($ind eq $Genome[-1]) {
              print VCF "0|0:22\t$line\n";
              } else {
              print VCF "0|0:22\t$line\n";
              }
           } elsif ($check <= 2) {
             if ($ind eq $Genome[-1]) {
              print VCF "1|1:22\t$line\n";
              } else {
              print VCF "1|1:22\t$line\n";
              }
           } else {
             if ($ind eq $Genome[-1]) {
              print VCF "\.|\.:22\t$line\n";
              } else {
              print VCF "\.|\.:22\t$line\n";
              }
           }
    } else {
       if ($ind eq $Genome[-1]) {
        print VCF "\.|\.:22\t$line\n";
        } else {
        print VCF "\.|\.:22\t$line\n";
        }
     }
   ###
} elsif ($svtype eq "ins") {
my $beg0 = $beg - 2;
my $end0 = $end + 2;
   if (exists ${$liftover_coord}{$beg0} && exists ${$liftover_coord}{$end0}) {
       my @tmp1 = split (/_/, ${$liftover_coord}{$beg0});
       my @tmp2 = split (/_/, ${$liftover_coord}{$end0});
       my $check = abs ($tmp2[1]-$tmp1[1]);
       if ($check == $len + 4) {
          if ($ind eq $Genome[-1]) {
           print VCF "1|1:22\t$line\n";
           } else {
           print VCF "1|1:22\t$line\n";
           }
        } elsif ($check <= 4) {
          if ($ind eq $Genome[-1]) {
           print VCF "0|0:22\t$line\n";
           } else {
           print VCF "0|0:22\t$line\n";
           }
        } else {
          if ($ind eq $Genome[-1]) {
           print VCF "\.|\.:22\t$line\n";
           } else {
           print VCF "\.|\.:22\t$line\n";
           }
        }
  } else {
    if ($ind eq $Genome[-1]) {
     print VCF "\.|\.:22\t$line\n";
     } else {
     print VCF "\.|\.:22\t$line\n";
     }
  }
 }
}; ##eval

if ($@) {
  print ERR "Failed\t$line\n";
}

} ## while
   $pm->finish;
} ## chr
   $pm->wait_all_children;
}

######################################################################################
###########################  Start to write to VCF file ##############################
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

########################################################################################
open SVtotal, "$svraw" or die "$!";
open VCF, ">$svraw.vcf" or die "$!";

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

#########################################################
############################ Subroutines ################
#########################################################

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
################################################################
