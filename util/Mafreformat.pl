#!/usr/bin/perl
use strict;
use warnings;

##########################
local $/="\n";
my %sizes;
open Size, "$ARGV[1]" or die $!;;
while (<Size>) {
chomp;
my @temp = split (/\t/,$_);
$sizes{$temp[0]} = $temp[1];
}
close Size;
##########################


$/="\na";
open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].reformat.maf" or die "$!";

while (<In>) {
chomp;
if ($_=~/\#/) {
print Out  "$_\n";
} else {
my @tem = split (/\n/,$_);
my @unit1= split (/\s+/,$tem[1]);
my @unit2= split (/\s+/,$tem[2]);
my @unit = split (/\./,$unit2[1]);
my $qname = join (".",($unit[0],$unit[1]));
my $beg;
my $qsize = $sizes{"$qname"};
if ($unit2[4]=~/\+/) {
$beg = $unit[2] + $unit2[2] -1;
} elsif ($unit2[4]=~/\-/) {
$beg = $qsize - $unit[3] + $unit2[2] ;
}

print Out "a $tem[0]\n$unit1[0]\t$unit1[1]\t$unit1[2]\t$unit1[3]\t$unit1[4]\t$unit1[5]\t$unit1[6]\n$unit2[0]\t$qname\t$beg\t$unit2[3]\t$unit2[4]\t$qsize\t$unit2[6]\n\n";

 }
}

print Out "\#\#eof maf";

