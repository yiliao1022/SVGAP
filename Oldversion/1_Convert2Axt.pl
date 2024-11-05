#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin";
use lib "$FindBin::Bin/./pub";
use lib "$FindBin::Bin/./lib";
use Parallel::ForkManager;
use Pod::Usage;
####################################################################################
# Name: 1_Convert2Axt.pl
#
# Authors: Yi Liao (12/15/2022)
#
# Usage: Convert whole genome alignment files from Last (.maf), Lastz (.maf), 
# MuMMer (.delta), minimap2 (.paf), and Anchorwave (.maf) to .axt format which 
# is used for downstream SV discovery and genotyping.
# 
# Importtant Prerequisite:  UCSC utilities;
#
# Usage: perl 1_Convert2Axt.pl --help
#
####################################################################################

my ($Aligner,$WK,$input,$Tname,$Qname,$MUM320,$UCSC,$UTIL,$TBA,$Paftools,$K8,$t,$Help,$verbose);

# read parameters

GetOptions( 'tname=s'=>\$Tname,
            'qname=s'=>\$Qname,
            'ali=s'=>\$Aligner,
            'input=s'=>\$input,
            'wk=s'=>\$WK,
            'mum320=s'=>\$MUM320,
            'ucsc=s'=>\$UCSC,
            'tba=s'=>\$TBA,
            'util=s'=>\$UTIL,
            'k8=s'=>\$K8,
            'paftools=s'=>\$Paftools,
            'query'=>\$verbose,
            't=i'=>\$t,
            'help'=>\$Help
           );

####################### 
 
# pre-defined program path. It will use the ones in $PATH if unspecified.

$K8 ||="$Bin/pub/MINIMAP2/";
$Paftools ||="$Bin/pub/MINIMAP2/";
$TBA ||="$Bin/pub/TBA/";
$UCSC ||="$Bin/pub/UCSC/";
$MUM320 ||="$Bin/pub/MUMMER/";
$UTIL ||="$Bin/util/";
$WK ||='.';
$input ||='.';
$Tname ||="Ref";
$Qname ||="Que";
$t ||= 4;

########################################################################################
if ($Help){
print <<"END.";
Usage:  perl $0 --ali minimap2 --input /path/to/align/ -wk ./output/

Options:    --ali      [last|lastz|mummer|minimap2|anchorwave]                            [Required]
            --input    [path] to the folder store the alignment files                     [Required]
            --wk       [path] to a working folder save the output (default: current)      [Required]
            --query    verbose (produce query-as-reference alignment if invoke)           [Optional]
            --tname    [character] Name of target genome assembly  (default: Ref)         [Optional]
            --qname    [character] Name of query genome assembly   (default: Que)         [Optional]
            --t        [int] Number of theads to run this script (default: 4)             [Optional]
            --help     Print this help information                                        [Optional]
END.
exit;
}

#####
my $message = "Try 'perl $0 --help' for more information.";

if (! $Aligner){
    pod2usage( {
           -message => $message,
           -verbose => 0,
           -exitval => 1 } );
    }




########################################################################################

#### set output directories
`mkdir $WK/Target_$Tname` unless -e "$WK/Target_$Tname" && -d "$WK/Target_$Tname";

if ($verbose) {
`mkdir $WK/Target_$Qname` unless -e "$WK/Target_$Qname" && -d "$WK/Target_$Qname";
}

############################## read the alignment file names to an array. 

opendir (DIR,$input) || die "Error in opening dir $input\n";
my $type;
my @alignment;

if ($Aligner eq "minimap2") {
while ((my $filename = readdir (DIR))) {
  if ($filename =~/(.*).chr.paf/) {
    push (@alignment, $filename);
    $type = "minimap2_chr";
  } elsif ($filename =~/(.*).ctg.paf/) {
    push (@alignment, $filename);
    $type = "minimap2_ctg";
  }
}

} elsif ($Aligner eq "mummer") {
while ((my $filename = readdir (DIR))) {
  if ($filename =~/(.*).delta/) {
    push (@alignment, $filename);
    $type = "mummer";
  }
}

} elsif ($Aligner eq "last") {
while ((my $filename = readdir (DIR))) {
  if ($filename =~/(.*).maf/) {
    push (@alignment, $filename);
    $type = "last";
  }
}

} elsif ($Aligner eq "lastz") {
while ((my $filename = readdir (DIR))) {
  if ($filename =~/(.*).maf/) {
   push (@alignment, $filename);
   $type = "lastz";
  }
}

} elsif ($Aligner eq "anchorwave") {
while ((my $filename = readdir (DIR))) {
 if ($filename =~/(.*).maf/) {
   push (@alignment, $filename);
   $type = "anchorwave";
  }
}

} else {
print "Please provide information regarding the aligner you utilized to generate the genome alignment!!!"
}

############################################################

my $pm = Parallel::ForkManager->new($t);

DATA_LOOP:
  foreach my $file (@alignment) {
  my $pid = $pm->start and next DATA_LOOP;

######### start #########
if ($type eq "minimap2_chr") {
$file =~/(.*)vs(.*).chr.paf/;
`sort -k6,6 -k8,8n $input/$file > $WK/$1.$2.sorted.paf;
${K8}k8 ${Paftools}paftools.js call $WK/$1.$2.sorted.paf > $WK/$1.$2.sorted.var.txt;
${K8}k8 ${Paftools}paftools.js view -f maf $WK/$1.$2.sorted.paf | sed 's/^a /a score=/g; s/_/./g' > $WK/$1.$2.sorted.maf;
${UCSC}mafToAxt $WK/$1.$2.sorted.maf $1 $2 $WK/$1vs$2.axt;
mv $WK/$1vs$2.axt $WK/Target_$Tname;
`;

if ($verbose) {
`${TBA}maf_order $WK/$1.$2.sorted.maf $2 $1 > $WK/$2.$1.sorted.maf;
${UCSC}mafToAxt $WK/$2.$1.sorted.maf $2 $1 $WK/$2vs$1.axt;
mv $WK/$2vs$1.axt $WK/Target_$Qname;
rm $WK/$2.$1.sorted.maf;
`;
}

system "rm $WK/$1.$2.sorted.maf $WK/$1.$2.sorted.paf";

} elsif ($type eq "minimap2_ctg") {  ### Please provide the fasta size files in current folder
$file=~/(.*)vs(.*).ctg.paf/;

`sort -k6,6 -k8,8n $input/$file > $WK/$1.$2.sorted.paf;
${K8}k8 ${Paftools}paftools.js call $WK/$1.$2.sorted.paf > $WK/$1.$2.sorted.var.txt;
${K8}k8 ${Paftools}paftools.js view -f maf $WK/$1.$2.sorted.paf | sed 's/^a /a score=/g; s/_/./g' > $WK/$1.$2.sorted.maf;
perl ${UTIL}/Mafreformat.pl $WK/$1.$2.sorted.maf $WK/$2.sizes; 
${UCSC}mafToAxt $WK/$1.$2.sorted.maf.reformat.maf $1 $2 $WK/$1vs$2.axt;
mv $WK/$1vs$2.axt $WK/Target_$Tname;
rm $WK/$1.$2.sorted.maf $WK/$1.$2.sorted.paf;
`;

if ($verbose) {
`${TBA}maf_order $WK/$1.$2.sorted.maf.reformat.maf $2 $1 > $WK/$2.$1.sorted.maf.reformat.maf;
${UCSC}mafToAxt $WK/$2.$1.sorted.maf.reformat.maf $2 $1 $WK/$2vs$1.axt;
mv $WK/$2vs$1.axt $WK/Target_$Qname;
rm $WK/$2.$1.sorted.maf.reformat.maf;
`;
}
system "rm $WK/$1.$2.sorted.maf.reformat.maf";

} elsif ($type eq "mummer") {
$file =~/(.*)vs(.*).delta/;
`${MUM320}delta2maf -r $input/$file > $WK/$1.$2.sorted.maf; 
${UCSC}mafToAxt $WK/$1.$2.sorted.maf $1 $2 $WK/$1vs$2.axt;
mv $WK/$1vs$2.axt $WK/Target_$Tname;
`;
if ($verbose) {
`${TBA}maf_order $WK/$1.$2.sorted.maf $2 $1 > $WK/$2.$1.sorted.maf;
 ${UCSC}mafToAxt $WK/$2.$1.sorted.maf $2 $1 $WK/$2vs$1.axt;
mv $WK/$2vs$1.axt $WK/Target_$Qname;
rm $WK/$2.$1.sorted.maf;
`;
}
system "rm $WK/$1.$2.sorted.maf";

} elsif ($type eq "last") {

$file=~/(.*)vs(.*).maf/;
`sed -i '1i \#\#maf version=1' $WK/$file; 
${TBA}maf_sort $WK/$file $1 > $WK/$1.$2.sorted.maf;
${UCSC}mafToAxt $WK/$1.$2.sorted.maf $1 $2 $WK/$1vs$2.axt;
mv $WK/$1vs$2.axt $WK/Target_$Tname;
`;

if ($verbose) {
`${TBA}maf_order $WK/$1.$2.sorted.maf $2 $1 > $WK/$2.$1.sorted.maf;
${UCSC}mafToAxt $WK/$2.$1.sorted.maf $2 $1 $WK/$2vs$1.axt;
mv $WK/$2vs$1.axt $WK/Target_$Qname;
rm $WK/$2.$1.sorted.maf;
`;
}
system "rm $WK/$1.$2.sorted.maf"; 

} elsif ($type eq "lastz") {
$file =~/(.*)vs(.*).maf/; 

`sed -i '1i \#\#maf version=1' $file; 
${TBA}maf_sort $WK/$file $1 > $WK/$1.$2.sorted.maf;
${UCSC}mafToAxt $WK/$1.$2.sorted.maf $1 $2 $WK/$1vs$2.axt;
mv $WK/$1vs$2.axt $WK/Target_$Tname;
`;

if ($verbose) {
`${TBA}maf_order $WK/$1.$2.sorted.maf $2 $1 > $WK/$2.$1.sorted.maf;
${UCSC}mafToAxt $WK/$2.$1.sorted.maf $2 $1 $WK/$2vs$1.axt;
mv $WK/$2vs$1.axt $WK/Target_$Qname;
rm $WK/$2.$1.sorted.maf;
`;
}

system "rm $WK/$1.$2.sorted.maf"; 

} elsif ($type eq "anchorwave") {
$file =~/(.*)vs(.*).maf/;
&reformatAnchorWavemaf ("$file");

`${TBA}maf_sort $WK/$1.$2.trans.maf $1 > $WK/$1.$2.sorted.maf;
${UCSC}mafToAxt $WK/$1.$2.sorted.maf $1 $2 $WK/$1vs$2.axt;
mv $WK/$1vs$2.axt $WK/Target_$Tname;
`;

if ($verbose) {
`${TBA}maf_order $WK/$1.$2.sorted.maf $2 $1 > $WK/$2.$1.sorted.maf;
${UCSC}mafToAxt $WK/$2.$1.sorted.maf $2 $1 $WK/$2vs$1.axt;
mv $WK/$2vs$1.axt $WK/Target_$Qname;
rm $WK/$2.$1.sorted.maf;
`;
}

system "rm $WK/$1.$2.sorted.maf";

} else {
print "Please provide information regarding the aligner you utilized to generate the genome alignment!!!";
}

######### end #########
$pm->finish;
}

$pm->wait_all_children;

##############################
##################################################
sub reformatAnchorWavemaf {
my ($file1) = @_;
open MAF, "$input/$file1" or die "$!";
$file1 =~/(.*)vs(.*).maf/;
open TRS,">$WK/$1.$2.trans.maf" or die "$!";
$/="\na";
while (<MAF>) {
if ($_=~/#/) {
print TRS "$_\n";
next;
} else {
my @temp = split (/\n/,$_);
my @unit1 = split (/\s+/,$temp[1]);
my @unit2 = split (/\s+/,$temp[2]);
if ($unit2[3] != 0 and $unit1[3] != 0) {
print TRS "a$temp[0]\n";
print TRS "$unit1[0]\t$1.$unit1[1]\t$unit1[2]\t$unit1[3]\t$unit1[4]\t$unit1[5]\t$unit1[6]\n";
print TRS "$unit2[0]\t$2.$unit2[1]\t$unit2[2]\t$unit2[3]\t$unit2[4]\t$unit2[5]\t$unit2[6]\n";
print TRS "\n";
   }
  }
 }
}

