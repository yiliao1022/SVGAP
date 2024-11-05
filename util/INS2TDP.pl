#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


###################
my ($ref,$qry,$ins,$outfile,$Help);

GetOptions(
  "ref:s"=>\$ref,
  "ins:s"=>\$ins,
  "qry:s"=>\$qry,
  "out:s"=>\$outfile,
  "help"=>\$Help
);

####################
if ($Help){
print <<"END.";
  Usage: perl $0 -ref sim.fasta -ins INS.bed -out outfile.bed
  -ref        the reference genome                                  [REQUIRED]
  -qry        the query genome                                      [REQUIRED]
  -ins        the called insertions                                 [REQUIRED]
  -out        out put                                               [REQUIRED]
  -help       print this
END.
exit;
}

####################
my $seqs = &ReadFa("$ref"); ### Read fasta sequences
my $query = &ReadFa("$qry");

open In, "$ins" or die "Can't open the reformated coordinate file";
open Out, ">$outfile" or die "Can't output!!!";

while (<In>) {
chomp;
next if ($_=~/\#/);

my @temp = split (/\t/,$_);
my $len = $temp[2] - $temp[1] + 1;
my $fl5 = $temp[5] - $len;
my $gap1 = $temp[1]-1;
my $fl3 = $temp[6];

my $sequence5 = substr ($$query{$temp[4]},$fl5,$len);
my $sequence0 = substr ($$seqs{$temp[0]},$gap1,$len);
my $suquence3 = substr ($$query{$temp[4]},$fl3,$len);

my ($iden,$orient) = &SeqIdentity ($sequence0,$sequence5,$suquence3);

#open Out, ">$outfile" or die "Can't output!!!";

if ($iden==1) {
  print "$orient\n";
  if ($orient=5) {
     my $beg = $fl5;
     my $end = $temp[5];
     print Out "$temp[0]\t$temp[1]\t$temp[2]\tCNV\t$temp[4]\t$beg\t$end\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[10]\t$temp[11]\t$temp[12]\t$temp[13]\t$temp[14]\t$temp[15]\t$temp[16]\n";

 } elsif ($orient=3) {
     my $beg = $temp[6];
     my $end = $temp[6]+$len;
     print Out "$temp[0]\t$temp[1]\t$temp[2]\tCNV\t$temp[4]\t$beg\t$end\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[10]\t$temp[11]\t$temp[12]\t$temp[13]\t$temp[14]\t$temp[15]\t$temp[16]\n";

 }
}
}


######################################
######################################
##############  Modules ##############
######################################
######################################

sub ReadFa {
local $/="\n>";
open(FA, "<@_") or die("Cannot open FASTA file\n");
my %FA;
while (my $sequence = <FA>) {
         chomp $sequence;
         my ($id) = $sequence =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
         $sequence =~ s/^>*.+\n//;  # remove FASTA header
         $sequence =~ s/\n//g;  # remove endlines
         $FA{$id} = $sequence;
        }
my $seq = \%FA;
return $seq;
}

#####
sub FaSize {
local $/="\n";
open(Size, "<@_") or die ("Cannot open FASTA file\n");
my %Size;
while (my $line = <Size>) {
         chomp $line;
         my @temp = split (/\t/,$line);
         $Size{$temp[0]} = $temp[1];
        }
my $fasize = \%Size;
return $fasize;
}


#####
sub RevCom {
my $origin_seq = shift;
my $revcomp = reverse $origin_seq;
$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;
return $revcomp;
}

#####
sub SeqIdentity {
my ($gap,$flank5,$flank3) = @_;
my @seeds;
my $num;
my $i;
my $ident=0;
my $len = length ($gap);
$num = $len-10+1;

if ($num==0) {
next;
}
  
for ($i=0;$i<$num;$i++) {
  my $seed = substr ($gap,$i,10);
  push (@seeds,$seed);
}

my $m1=0;
my $m2=0;
  
foreach my $sed1 (@seeds) {
  if ($flank5=~/$sed1/) {
    $m1++;
   }
  }
  
foreach my $sed2 (@seeds) {
  if ($flank3=~/$sed2/) {
    $m2++;
   }
  }

my $perc;
my $flank;

my $perc1 = $m1/$num;
my $perc2 = $m2/$num;

if ($perc1>$perc2) {
$perc = $perc1;
$flank=5;
} else {
$perc = $perc2;
$flank=3;
}


if ($perc > 0.5) {
    $ident =1;
}
  return ($ident,$flank);
}

