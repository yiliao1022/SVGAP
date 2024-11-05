#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(&basename &dirname);

####################################################################################
# SynNetFilter.pl
#
# Authors: Yi Liao (12/16/2022)
#
# Prerequisite: Perl, Kent's Utilities, TBA program.
#
# Usage: Remove non-syntenic fills from *.chain.filter.tnet.synnet, 
# only keep reliable syntenic and conserved orthologous alignments.
#
####################################################################################

my ($synnet,$GenomesDir,$chain,$one2mul,$toplen,$seclen,$nonsynlen,$synlen,$invlen,$alignrate,$TBA,$UCSC,$target,$rawsv,$query,$qFarCutoff,$singfolder,$axtfolder,$filteredSynNet,$homChrpairs,$Help);

# get input parameters
GetOptions(
'synnet=s'=>\$synnet,
'gd=s'=>\$GenomesDir,
'chain=s'=>\$chain,
'one2mul=s'=>\$one2mul,
'alrate=f'=>\$alignrate,
'toplen=i'=>\$toplen,
'nonsynlen=i'=>\$nonsynlen,
'synlen=i'=>\$synlen,
'invlen=i'=>\$invlen,
'qFarCutoff=i'=>\$qFarCutoff,
'ucsc=s'=>\$UCSC,
'rawsv=s'=>\$rawsv,
'tba=s'=>\$TBA,
'target=s'=>\$target,
'query=s'=>\$query,
'sing=s'=>\$singfolder,
'axt=s'=>\$axtfolder,
'cleanSynnet=s'=>\$filteredSynNet,
'homopairs=s'=>\$homChrpairs,
'rawsv=s'=>\$rawsv,
'help'=>\$Help
);

############ Printing help information
if (defined ($Help) || !defined($synnet)){
print <<"END.";
Usage: perl $0 --synnet ./synnet/ --g ./genome/ --chain /path/to/all.chain.filter/ 

Required:

--synnet  <s> : Folder kept the Chain/Net/Synnet files
--gd       <s> : Folder kept genome squences
--chain   <s> : Folder kept the .chain files

Optional:

--one2mul       <s> : List of sequences in query that were aligned to the reference
--rawsv         <s> : output of the syntenic gaps information
--alrate        <i> : the minimum alignable proportion for chains to be processed.
--toplen        <i> : the minimum length for the top syntenic chains to be processed.
--nonsynlen     <i> : the minimum length for the top non-syntenic chains to be processed.
--synlen        <i> : the minimum length for any lower chains that is to be processed.
--invlen        <i> : the minimum length for any chains that is to be considered as inversion.
--qFarCutoff    <i> : the maximum number of the qFar value in the net file when predicting a inversion event.
--help              : print this help information

END.
exit;
}

####################################  Pre-defined parameters
$UCSC ||="$Bin/pub/UCSC/";
$TBA ||="$Bin/pub/TBA/";

$singfolder = 'filtered_sing';
$axtfolder = 'axt';
$filteredSynNet = 'CleanSynNet';
$homChrpairs = 'Hompairs';
$rawsv = 'RawSV';

my $dir = $ENV{'PWD'};
`mkdir $dir/$singfolder` unless -e "$dir/$singfolder" && -d "$dir/$singfolder";
`mkdir $dir/$axtfolder` unless -e "$dir/$axtfolder" && -d "$dir/$axtfolder";
`mkdir $dir/$filteredSynNet` unless -e "$dir/$filteredSynNet" && -d "$dir/$filteredSynNet";
`mkdir $dir/$homChrpairs` unless -e "$dir/$homChrpairs" && -d "$dir/$homChrpairs";
`mkdir $dir/$rawsv` unless -e "$dir/$rawsv" && -d "$dir/$rawsv";

###################################################################################################

=pod
opendir (DIR,$GenomesDir) || die "Error in opening dir $GenomesDir\n";

while ((my $filename = readdir (DIR))) {
  if ($filename !~/\.2bit/ and $filename !~/\.sizes/) {
      unless (-s "$GenomesDir/$filename.sizes") {
        `${UCSC}faSize -detailed $GenomesDir/$filename > $GenomesDir/$filename.sizes`;
         }
      unless (-s "$GenomesDir/$filename.2bit") {
        `${UCSC}faToTwoBit $GenomesDir/$filename $GenomesDir/$filename.2bit`;
         }
  } 
}
=cut

####################### Report single covege of alignments relative to the reference (sing.maf)
opendir (SNY,$synnet) || die "Error in opening dir $synnet\n";

while ((my $filename = readdir (SNY))) {
  if ($filename =~/(.*)vs(.*).chain.filter.tnet.synnet/) {
       print "Processing $filename\n";
       &filter_SynNet("$filename");
       &filter_NET("$filename.refiltered.synnet");
  #####    $singfolder $axtfolder $filteredSynNet $homChrpairs $rawsv 
   `${UCSC}netToAxt $filteredSynNet/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet $chain/$1.$2.chain.filter $GenomesDir/$1.2bit $GenomesDir/$2.2bit $axtfolder/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet.axt`;
   `${UCSC}axtToMaf $axtfolder/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet.axt $GenomesDir/$1.sizes $GenomesDir/$2.sizes $synnet/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet.axt.maf`;
   `${TBA}single_cov2 $synnet/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet.axt.maf > $singfolder/$1.$2.sing.maf`;
   `rm $synnet/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet.axt.maf`;
       
       ####################### Quickly obtain syntenic gaps
       open Gap, "$filteredSynNet/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet" or die "$!";
       open SV, ">$rawsv/$1vs$2.chain.filter.tnet.synnet.refiltered.synnet.rawsv.bed" or die "$!";
       my $Chr;
       while (<Gap>) {
             chomp;
             $_ =~ s/^\s+//;
             my @splits = split (/\s/,$_);
             if ($_=~/^net/) {
                   $Chr = $splits[1];
             } elsif ($_=~/^gap/) {
                my $end = $splits[1]+$splits[2];
                print SV "$Chr\t$splits[1]\t$end\n";
             }
        }
      close Gap;
      close SV;
      #######################
          
  }  elsif ($filename =~/(.*)vs(.*).chain.filter.qnet.synnet/) { 
       print "Processing $filename\n";
       &filter_SynNet("$filename");
       &filter_NET("$filename.refiltered.synnet");
       
        ####################### Quickly obtain syntenic gaps
       open Gap, "$filteredSynNet/$1vs$2.chain.filter.qnet.synnet.refiltered.synnet" or die "$!";
       open SV, ">$rawsv/$1vs$2.chain.filter.qnet.synnet.refiltered.synnet.rawsv.bed" or die "$!";
       my $Chr;
       while (<Gap>) {
             chomp;
             $_ =~ s/^\s+//;
             my @splits = split (/\s/,$_);
             if ($_=~/^net/) {
                   $Chr = $splits[1];
             } elsif ($_=~/^gap/) {
                my $end = $splits[1]+$splits[2];
                print SV "$Chr\t$splits[1]\t$end\n";
             }
        }
      close Gap;
      close SV;
      #######################
  }
}

##########################
#######sub Function ######
##########################
####################################################################

sub array_search {
    my ($arr,$elem) = @_;
    my $idx;
    for my $i (0..$#$arr) {
        if ($arr->[$i] eq $elem) {
            $idx = $i;
            last;
         }
     }
  return $idx;            
}

####################################################################

sub filter_NET {
   my $file = shift;
   local $/="\nnet";
   open NET, "$filteredSynNet/$file" or die "$!";
   open NETout, ">$filteredSynNet/$file.out" or die "$!";
   while (<NET>) {
   chomp;
   my @tmp = split (/\n/,$_);
   if ($#tmp>1) {
   print NETout "net$_\n";
   }
   }
}

####################################################################

sub filter_SynNet {

$toplen ||=500000;
$nonsynlen ||=500000;
$synlen ||=500000;
$invlen ||=10000;
$alignrate ||=0.1;
$qFarCutoff ||=5000000;

my $synnetInd= shift;

open In, "$synnet/$synnetInd" or die "$!";
open Out, ">$filteredSynNet/$synnetInd.refiltered.synnet" or die "$!";
open Out1, ">$homChrpairs/$synnetInd.pairedChr.bed" or die "$!";

my %hash;
my $tag=0;
my $tag1=0;
my $tag2=0;
my $chr;

while (<In>) {
chomp;
if ($_=~/^net|\#/) {
  print Out "$_\n";
     if ($_=~/^net/) {
         my @tmp = split (/\s/,$_);
         my @unit = split (/\./,$tmp[1]);
            $chr = $unit[1];
            print "Processing $chr\n";
      }
} elsif ($_=~/^\s{1}fill/) {
my @temp = split (/\s/, $_);
my @tmp1 = split (/\./,$temp[4]);
my $query_chr = $tmp1[1];
my $ref = \@temp;
my $ali_idx = &array_search ($ref,"ali");
my $topalrate = $temp[$ali_idx+1]/$temp[3];

 if ($temp[3]>$toplen and $topalrate > $alignrate) { # and ($chr eq $query_chr)) {
 print Out "$_\n";
 push (@{$hash{$chr}},$query_chr);
 $tag2=0;
 } else {
 $tag2=1;
 }
} elsif ($_=~/^\s{2}gap/ and $tag2==0) {
 print Out "$_\n";
} elsif ($_=~/^\s{3}fill/ and $_=~/nonSyn/) {

  my @temp = split (/\s+/, $_);
  my @tmp1 = split (/\./,$temp[4]);
  my $query_chr = $tmp1[1];
  my $ref = \@temp;
  my $ali_idx = &array_search ($ref,"ali");
  my $secalrate = $temp[$ali_idx+1]/$temp[3];

   if ($tag2==0 and $temp[3]>$nonsynlen and $secalrate > $alignrate) { # and ($chr eq $query_chr)) {
   print Out "$_\n";
   push (@{$hash{$chr}},$query_chr);
   $tag=0;
   } else {
   $tag=1;
   }

} elsif ($_=~/^\s{3}fill/ and $_=~/syn|inv/) {

     my @temp = split (/\s+/, $_);
     my @tmp1 = split (/\./,$temp[4]);
     my $query_chr = $tmp1[1];
     my $ref = \@temp;
     my $ali_idx = &array_search ($ref,"ali");
     my $qFar_idx = &array_search ($ref,"qFar");
     my $qOver_idx = &array_search ($ref,"qOver");
     my $secalrate = $temp[$ali_idx+1]/$temp[3];  

     if ($_=~/syn/) {
         if ( $tag2==0 and $temp[3]>$synlen and $secalrate > $alignrate) { # or (4*$temp[3]>$temp[$qFar_idx+1]))) {
              print Out "$_\n";
              $tag=0;
         } else {
              $tag=1;
         }
     } elsif ($_=~/inv/) {
       if ( $tag2 ==0 and $temp[3]>$invlen and $secalrate > $alignrate and $temp[$qOver_idx+1] > 0 and $temp[$qFar_idx+1]<$qFarCutoff) {
            print Out "$_\n";
            $tag=0;
       } elsif ( $tag2 ==0 && $temp[3]>5*$invlen && $temp[7] > 5*$invlen && $temp[$qFar_idx+1]<$temp[3]) {
            print Out "$_\n";
            $tag=0;
       } else {
            $tag=1;
       }    
     }
       
} elsif ($_=~/^\s{4}gap/ and $tag==0 and $tag2==0) {
        print Out "$_\n";
} elsif ($_=~/^\s{5}fill/ and $_=~/nonSyn/) {

  my @temp = split (/\s+/, $_);
  my @tmp1 = split (/\./,$temp[4]);
  my $query_chr = $tmp1[1];
  my $ref = \@temp;
  my $ali_idx = &array_search ($ref,"ali");
  my $secalrate = $temp[$ali_idx+1]/$temp[3];
   if ($tag==0 and $tag2==0 and $temp[3]>$nonsynlen and $secalrate > $alignrate) { # and ($chr eq $query_chr)) {
   print Out "$_\n";
   push (@{$hash{$chr}},$query_chr);
   $tag1=0;
   } else {
   $tag1=1;
   }
       
} elsif ($_=~/^\s{5}fill/ and $_=~/syn|inv/) {

  my @temp = split (/\s+/, $_);
  my @tmp1 = split (/\./,$temp[4]);
  my $query_chr = $tmp1[1];
  my $ref = \@temp;
  my $ali_idx = &array_search ($ref,"ali");
  my $qFar_idx = &array_search ($ref,"qFar");
  my $qOver_idx = &array_search ($ref,"qOver");
  my $secalrate = $temp[$ali_idx+1]/$temp[3];  

  if ($_=~/syn/) {
      if ( $tag2==0 and $tag==0 and $temp[3]>$synlen and $secalrate > $alignrate) { #or (4*$temp[3]>$temp[$qFar_idx+1]))) {
           print Out "$_\n";
           $tag1=0;
      } else {
           $tag1=1;
      }
  } elsif ($tag2==0 and $tag==0 and $_=~/inv/) {
    if ($temp[3]>$invlen and $secalrate > $alignrate and $temp[$qOver_idx+1] > 0 and $temp[$qFar_idx+1]<$qFarCutoff) {
         print Out "$_\n";
         $tag1=0;
    } elsif ($temp[3]>5*$invlen && $temp[7] > 5*$invlen && $temp[$qFar_idx+1]<$temp[3]) {
         print Out "$_\n";
         $tag1=0;
    } else {
         $tag1=1;
    }    
  }

} elsif ($_=~/^\s{6}gap/ and $tag1==0 and $tag2==0 and $tag==0) {
     print Out "$_\n";
} else {
  next;
  }
}


foreach my $chrom (sort{$hash{$a}<=>$hash{$b}} keys %hash)  {
  print Out1 "$chrom\t@{$hash{$chrom}}\n";
}

close In;
close Out;
close Out1;

}
