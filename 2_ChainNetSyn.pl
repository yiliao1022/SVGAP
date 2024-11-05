#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use lib "$FindBin::Bin";
use Pod::Usage;

####################################################################################
# ChainNetSynnet.pl
# Authors: Yi Liao (12/09/2022)
# Constructing Chains, Nets and adding syntenic information based on pairwise alignments
# (.axt format)(LAST, LASTZ, MUMmer, minimap2, et al...)
# Prerequisite: UCSC utilities;
####################################################################################

my ($GenomesDir,$sing,$TBA,$lst,$AxtsDir,$UCSC,$linearGap,$minscore,$SynnetOUTdir,$Help);

GetOptions( 
            'gd=s'=>\$GenomesDir,
            'ad=s'=>\$AxtsDir,
            'lst=s'=>\$lst,
            'ucsc=s'=>\$UCSC,
            'tba=s'=>\$TBA,
            'sing'=>\$sing,
            'lg=s'=>\$linearGap,
            'minscore=i'=>\$minscore,
            'syn=s'=>\$SynnetOUTdir,
            'help'=>\$Help,
            );

############ 
if ($Help){
print << "END.";
Usage:  perl $0 --gd /path/to/genomes/ --ad /path/to/axtfiles/ --lst genome.lst
Options:    --ad      Forder contains the axt alignment files                [REQUIRED]
            --gd      Forder contains all genome sequences                   [REQUIRED]
            --lst     a file list the genome IDs (one line per genome)       [REQUIRED]
            --sing    [verbose] if to report the sing-coverage maf file      [Optional]
            --linearGap linearGap file for ChainNet                          [Optional]
            --minscore  -minScore = 1000 (default)                           [Optional]
            --syn       Output folder for raw synnet files   default:synnet  [Optional]
            --help      print this help information
END.
exit;
}

my $message = "Try 'perl $0 --help' for more information.";

if (! $lst or ! $GenomesDir or !$AxtsDir){
    pod2usage( {
           -message => $message,
           -verbose => 0,
           -exitval => 1 } );
    }

############ set default values
$UCSC ||="$Bin/pub/UCSC/";
$linearGap ||="$Bin/pub/TBA/medium";
$SynnetOUTdir ||="synnet";
$TBA ||="$Bin/pub/TBA/";
$minscore ||=1000;

# set output directories
 my $dir = $ENV{'PWD'};
`mkdir $dir/$SynnetOUTdir` unless -e "$dir/$SynnetOUTdir" && -d "$dir/$SynnetOUTdir";
####################################################################################
open In, "$lst" or die "$!";
my @individual;
while (<In>) {
chomp;
push (@individual, $_);
}

foreach my $filename (@individual) {
  if ($filename !~/\.2bit/ and $filename !~/\.sizes/) {
      unless (-s "$GenomesDir/$filename.sizes") {
        `${UCSC}faSize -detailed $GenomesDir/$filename > $GenomesDir/$filename.sizes`;
         }
      unless (-s "$GenomesDir/$filename.2bit") {
        `${UCSC}faToTwoBit $GenomesDir/$filename $GenomesDir/$filename.2bit`;
         }
  }
}

###################################################################################

opendir (AXT,$AxtsDir) || die "Error in opening dir $AxtsDir\n";
while ((my $filename = readdir (AXT))) {
  if ($filename =~/(.*)vs(.*).axt/) {

`${UCSC}axtChain -linearGap=$linearGap -minScore=$minscore $AxtsDir/$1vs$2.axt $GenomesDir/$1.2bit $GenomesDir/$2.2bit $AxtsDir/$1.$2.chain`;
`${UCSC}chainPreNet $AxtsDir/$1.$2.chain $GenomesDir/$1.sizes $GenomesDir/$2.sizes $AxtsDir/$1.$2.chain.filter`;
`${UCSC}chainNet -minSpace=1 -minScore=$minscore $AxtsDir/$1.$2.chain.filter $GenomesDir/$1.sizes $GenomesDir/$2.sizes $AxtsDir/$1.$2.chain.filter.tnet $AxtsDir/$1.$2.chain.filter.qnet`;
`${UCSC}netSyntenic $AxtsDir/$1.$2.chain.filter.tnet $dir/$SynnetOUTdir/$1vs$2.chain.filter.tnet.synnet`;
`${UCSC}netSyntenic $AxtsDir/$1.$2.chain.filter.qnet $dir/$SynnetOUTdir/$1vs$2.chain.filter.qnet.synnet`;
`rm $AxtsDir/$1.$2.chain $AxtsDir/$1.$2.chain.filter.tnet $AxtsDir/$1.$2.chain.filter.qnet`;

    if ($sing) {
my $singDir = "singRaw";
unless(-d "$dir/$singDir"){`mkdir $dir/$singDir`};
`${UCSC}netToAxt $dir/$SynnetOUTdir/$1vs$2.chain.filter.tnet.synnet  $AxtsDir/$1.$2.chain.filter $GenomesDir/$1.2bit $GenomesDir/$2.2bit $AxtsDir/$1vs$2.chain.filter.tnet.synnet.axt`;
`${UCSC}axtToMaf $AxtsDir/$1vs$2.chain.filter.tnet.synnet.axt $GenomesDir/$1.sizes $GenomesDir/$2.sizes $AxtsDir/$1vs$2.chain.filter.tnet.synnet.axt.maf`;
`${TBA}single_cov2 $AxtsDir/$1vs$2.chain.filter.tnet.synnet.axt.maf > $dir/$singDir/$1.$2.sing.maf`;
`rm $AxtsDir/$1vs$2.chain.filter.tnet.synnet.axt $AxtsDir/$1vs$2.chain.filter.tnet.synnet.axt.maf`;
        }
  } 
}
