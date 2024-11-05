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
# GenometypingDeletions_Sing.pl
#
# Written by Yi Liao (01/04/2023)
#
# Further computationally validated the merged deletions
#
####################################################################################
my ($svraw,$singmaf,$UCSC,$TBA,$refsize,$refname,$speciesname,$refseq,$OUTPUTdir,$extend,$p,$help);

GetOptions(   
    "del=s"=>\$svraw, # Path to the combined non-redundant SV data.
    "singmaf=s"=>\$singmaf, # Directory for single coverage maf files
    "ucsc=s"=>\$UCSC, # Path to Kent's utilities
    "tba=s"=>\$TBA, # Path to TBA utilities
    "refsize=s"=>\$refsize,
    "refseq=s"=>\$refseq,
    "refname=s"=>\$refname,
    "speciesname=s"=>\$speciesname,
    "extend=s"=>\$extend,
    "outdir=s"=>\$OUTPUTdir,
    "t=i"=>\$p,
    "help"=>\$help,
);

############# Help information ##########
if ($help){
print <<"END.";
Usage : perl $0 --del ./Combined.del.txt --refseq ./Ref/MH63 --refsize ./Ref/MH63.sizes --refname MH63 --singmaf ./sing/

--del          the file of the merged non-redundant deletions (> 49 bp)                     [REQUIRED]
--singmaf      folder kept single-coverage maf files (on the reference side)                [REQUIRED]
--refseq       reference sequence                                                           [REQUIRED]
--refsize      chromosome sizes of the reference                                            [REQUIRED]
--refname      reference name, e.g. hg38, IRGSP1                                            [REQUIRED]
--speciesname  species name, e.g. human, Drosophila, rice , maize                           [Default: Undefined]
--ucsc         Kent's utilities                                                             [Default: $Bin/pub/TBA/]
--tba          path for TBA program utilities                                               [Default: $Bin/pub/TBA/]
--extend       the length of sequence around the deletion breakpoints for genotyping        [Default: 1000]
--outdir       set an output dir                                                            [Default: VCFtemp]
--t            how many threads will use                                                    [Default: 4]
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
##################################
$UCSC ||="$Bin/pub/UCSC/";
$TBA ||="$Bin/pub/TBA/";
$extend ||= 1000;
$speciesname ||= "Undefined";
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
              open ($file, '>>', "$sv_path/$refname.$ctg.DEL.bed") or die "Could not open file '$file' $!";
              print $file "$_\n";
      }
    }
}

##################################
###### Start to genotype #########
##################################
foreach my $ind (@Genome) {
  my $pm = Parallel::ForkManager->new($p);

DATA_LOOP:
  foreach my $scf (@chromosomes) {
         my $pid = $pm->start and next DATA_LOOP; 
         print "Currently are processing chromosome $scf of sample $ind !!!\n";
         my $liftover_coord = &liftOver("$singmaf/$ind/$refname\_$ind\_$scf.maf");
############################################################################################################## 
open IN, "$sv_path/$refname.$scf.DEL.bed" or die "$!"; 
open VCF, ">$OUTPUTdir/$refname.$ind.$scf.vcf" or die "$!";
open ERR, ">$OUTPUTdir/err/$refname.$ind.$scf.err.vcf" or die "$!";

while (<IN>) {
chomp;
my %exi;
my $line = $_;
my ($chr,$beg,$end,$len,$sample,$qchr,$qs,$qe,$fiv,$thr) = (split (/\t/,$_))[0,1,2,5,7,8,9,10,11,12];
my @samples = split (/,/,$sample);
foreach my $insample (@samples) {
  my @unitk = split (/\./,$insample);
  $exi{$unitk[0]} = 0;
}
if (exists $exi{$ind}) {
  print VCF "1|1:22\t$line\n";
  next;
}

my $ref = (split(/\./,$chr))[0];
my $alt = (split(/\./,$qchr))[0];
my $flank;

if ($len <= 100) {
  $flank = 100;
} elsif ($len > 100 and $len < 500 ) {
  $flank = 200;
} elsif ($len >= 500 && $len < 2000) {
  $flank = 500;
} elsif ($len >=2000 && $len < 5000) {
  $flank = 2000;
}  else {
  $flank = 5000;
}

if ($fiv > $extend) {
	$fiv = $extend;
} else {
        $fiv = $fiv;
}

if ($thr > 1000) {
        $thr = $extend;
} else {
        $thr = $thr;
}

my $beg0 = $beg - $fiv + 5;
my $end0 = $end + $thr - 5;

eval {

  if (exists ${$liftover_coord}{$beg0} && exists ${$liftover_coord}{$end0}) {
     my @begflank1= split (/_/,${$liftover_coord}{$beg0});
     my @endflank1 = split (/_/,${$liftover_coord}{$end0});
     
     if ( abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10) < 13 ) {
          print VCF "1|1:22\t$_\n";                 
     } elsif ( (abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10)) > $len-13 && (abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10)) < $len + 13) {
          print VCF "0|0:22\t$_\n";  
     } elsif ((abs ($begflank1[1]-$endflank1[1]) - ($fiv + $thr - 10)) > 200000)  {
          print VCF "\.|\.:22\t$_\n";     
     } else {
       my $begin = $beg - $flank;
       my $endin = $end + $flank;

       `${TBA}maf2fasta $ref_path/byChr/$refname.$scf.fa $singmaf/$ind/$refname\_$ind\_$scf.maf $begin $endin > $OUTPUTdir/TBA.$beg\_$end.fasta`;

             if (-z "$OUTPUTdir/TBA.$beg\_$end.fasta") {
                    system "rm $OUTPUTdir/TBA.$beg\_$end.fasta";     
                    print VCF "\.|\.:22\t$_\n";
             } else {
                    my $evaluation = &ParseMaf2fasta("$OUTPUTdir/TBA.$beg\_$end.fasta",$flank,$ind,$len);
                    print VCF "$evaluation|$evaluation:22\t$line\n"; 
                    `rm $OUTPUTdir/TBA.$beg\_$end.fasta`;
	     }
     }       
  } else {
      my $begin = $beg - $flank;
      my $endin = $end + $flank;
     
      `${TBA}maf2fasta $ref_path/byChr/$refname.$scf.fa $singmaf/$ind/$refname\_$ind\_$scf.maf $begin $endin > $OUTPUTdir/TBA.$beg\_$end.fasta`;

        if (-z "$OUTPUTdir/TBA.$beg\_$end.fasta") {
               system "rm $OUTPUTdir/TBA.$beg\_$end.fasta";
               print VCF "\.|\.:22\t$_\n";	       
        } else {
               my $evaluation = &ParseMaf2fasta("$OUTPUTdir/TBA.$beg\_$end.fasta",$flank,$ind,$len);
               print VCF "$evaluation|$evaluation:22\t$line\n";   
               `rm $OUTPUTdir/TBA.$beg\_$end.fasta`;
        }
  }  
    
 }; #eval
  
  if ($@) {
             print ERR "Failed\t$_\n";
    }  
    
}
##############################################################################################################			   
    $pm->finish;
  }
    $pm->wait_all_children;
}

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
print VCF "\#\#fileformat=VCFv4.2\n\#\#fileDate=$time\n\#\#source=SVgapV2023\n\#\#reference=$refname\n";
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
##ALT=<ID=DEL:ME:ALU,Description=\"Deletion of ALU element\">
##ALT=<ID=DEL:ME:L1,Description=\"Deletion of L1 element\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">
##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">
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
##############################################################################################################
##############################################################################################################
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

################################
sub ParseMaf2fasta {
my ($file,$extending,$query,$length) = @_;
my @data;
$/="\n";

open(DATA, "$file") or die "Couldn't open file sub maf2fasta file, $!";
while (<DATA>) {
chomp;
push (@data,$_);
}
close DATA;

my @seq_num = split (/\s/,$data[0]);
my @ref_name = split (/:/,$data[1]);
my $ref = $ref_name[0];
my @ref_seq = split(//,$data[3]);
my @query_seq = split(//,$data[4]);

my $len_n=0;
my $t;
my $start;
my $end;

for ($t=0;$t<$seq_num[1];$t++) {
 if ($ref_seq[$t] =~/[ATGCNatgcn]/) {
  $len_n++;
 if ($len_n==$extending+1) {
     $start = $t;
 }
 if ($len_n==$extending+$length) {
    $end = $t;
   }
  }
}

my $j;
my $m=0;
my $n=0;
my $h=0;
for ($j=0;$j<$start;$j++) {
     if ($ref_seq[$j] =~/[ATGCNatgcn]/) {
         $n++;
        $ref_seq[$j] = uc ($ref_seq[$j]);
       if ($query_seq[$j] =~/[ATGCNatgcn]/) {
           $h++;
           $query_seq[$j] = uc ($query_seq[$j]);
              if ($ref_seq[$j] eq $query_seq[$j]) {
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
for ($q=$end+1;$q<=$#ref_seq;$q++) {
        if ($ref_seq[$q] =~/[ATGCNatgcn]/) {
            $n1++;
           $ref_seq[$q] = uc ($ref_seq[$q]);
          if ($query_seq[$q] =~/[ATGCNatgcn]/) {
                 $h1++;
                 $query_seq[$q] = uc ($query_seq[$q]);
                 if ($ref_seq[$q] eq $query_seq[$q]) {
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
  if ($ref_seq[$y] =~/[ATGCNatgcn]/) {
    $n2++;
    $ref_seq[$y] = uc ($ref_seq[$y]);
    if ($query_seq[$y] =~/[ATGCNatgcn]/) {
       $h2++;
       $query_seq[$y] = uc ($query_seq[$y]);
       if ($ref_seq[$y] eq $query_seq[$y]) {
                    $m2++;
        }
     }
   }
 }

my $gap_region = $m2/($n2+1);

my $evalue;
 #  print "IDENTITY: $gap_identity5\t$gap_identity3\t$gap_region\n";
   if ($gap_identity5>=0.9 && $gap_identity3 >=0.9 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region >= 0.5) {
       $evalue=0; ### similar to reference
    } elsif ($gap_identity5>=0.9 && $gap_identity3 >=0.9 && $gap_overlap5 >=0.2 && $gap_overlap3 >=0.2 && $gap_region < 0.1) {
       $evalue=1; ### support a deletion
    } else {
       $evalue="\."; ### missing data
    }
return $evalue;
}
