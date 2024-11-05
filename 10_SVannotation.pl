#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

###################
my ($ref,$vcf,$te,$svtype,$cds,$outfile,$help);

GetOptions(
  "ref:s"=>\$ref,
  "vcf:s"=>\$vcf,
  "svtype:s"=>\$svtype,
  "te:s"=>\$te,
  "cds:s"=>\$cds,
  "out:s"=>\$outfile,
  "help"=>\$help
);

####################
if ($help){
print <<"END.";
  Usage: perl $0 --ref ref.fasta --vcf INS.vcf --svtype INS --te TE.fasta --cds Gene.fasta --out annotated.vcf
  --ref        the reference genome                                          [REQUIRED]
  --vcf        the SV vcf file                                               [REQUIRED]
  --svtype     INS or DEL                                                    [REQUIRED]
  --te         a TE liabrary (e.g. output from EDTA pipeline)                [REQUIRED]
  --cds        a cds file from this species or a closely related species     [REQUIRED]
  --out        annotated vcf file                                            [REQUIRED]
  --help       print this help information
END.
exit;
}
####################
my $seqs = &ReadFa("$ref"); ### Read reference sequences

&DUP($vcf,$svtype);
&convert_to_fasta("$vcf", "$vcf.fa", $svtype);
system qq (minimap2 -c $te $vcf.fa > TE.$svtype.paf);
system qq (minimap2 -c -x asm5 $cds $vcf.fa > CDS.$svtype.paf);
system qq (sort -k1,1 -k3,3n TE.$svtype.paf > TE.$svtype.sort.paf);
system qq (sort -k1,1 -k3,3n CDS.$svtype.paf > CDS.$svtype.sort.paf);
&GenomicOverlap("TE.$svtype.sort.paf");
&GenomicOverlap("CDS.$svtype.sort.paf");
&Perc("TE.$svtype.sort.paf.bed","TE");
&Perc("CDS.$svtype.sort.paf.bed","CDS");

############
my %cds;
my %te;
my %td;

open CDS,"CDS.$svtype.sort.paf.bed.perc.bed" or die "$!";
open TE,"TE.$svtype.sort.paf.bed.perc.bed" or die "$!";
open TDP,"$vcf.tdp.bed" or die "$!";

while (<CDS>) {
chomp;
my @temp = split (/\t/,$_);
$cds{$temp[0]} = $temp[1];
}

while (<TE>) {
chomp;
my @temp = split (/\t/,$_);
$te{$temp[0]} = $temp[1];
}

while (<TDP>) {
chomp;
my @temp = split (/\t/,$_);
$td{$temp[0]} = $temp[1];
}
#############
open VCF, "$vcf" or die "$!";
open OUTVCF, ">$outfile" or die "$!";

while (<VCF>) {
chomp;
if ($_=~/^#/) {
print OUTVCF "$_\n";
} else {
my @temp = split (/\t/,$_);
my $key = join ("_",($temp[0],$temp[1],$temp[2]));
my ($cdstmp,$tetmp,$tdtmp);

if (exists $cds{$key}) {
 $cdstmp = $cds{$key};
} else {
 $cdstmp = "NF";
}

if (exists $te{$key}) {
 $tetmp = $te{$key};
} else {
 $tetmp = "NF";
}

if (exists $td{$key}) {
 $tdtmp = $td{$key};
} else {
 $tdtmp = "NF";
}

my $replace = "$temp[7];$tdtmp;$cdstmp;$tetmp";
$temp[7] = $replace;
print OUTVCF join("\t",@temp) . "\n";

 }
}

close (VCF);
close (OUTVCF);

######################################
######################################
##############  Modules ##############
######################################
######################################
sub DUP {
my ($input_line,$svtype) = @_;
open In, "$input_line" or die "Can't open the vcf file";
open Out1, ">$input_line.tdp.bed" or die "$!";

while (<In>) {
chomp;
next if ($_=~/\#/);
my @temp = split (/\t/,$_);
my $sv = join ("_",($temp[0],$temp[1],$temp[2]));
my $len;
my $fl5;
my $fl3;
my $seq0;

if ($svtype eq "DEL") {
$seq0 = $temp[3];
$len = length ($temp[3]);
$fl5 = $temp[1] - $len;
$fl3 = $temp[1] + $len;
} elsif ($svtype eq "INS") {
$seq0 = $temp[4];
$len = length ($temp[4]);
$fl5 = $temp[1] - $len;
$fl3 = $temp[1];
}

my $seq5 = substr ($$seqs{$temp[0]},$fl5,$len);
my $seq3 = substr ($$seqs{$temp[0]},$fl3,$len);
my $seq0rev = &RevCom($seq0);

my ($iden,$orient) = &SeqIdentity ($seq0,$seq5,$seq3);
my ($idenrev,$orientrev) = &SeqIdentity ($seq0rev,$seq5,$seq3);

if ($iden==1) {
  if ($orient==5) {
     print Out1 "$sv\t5DUP\n";
  } elsif ($orient==3) {
     print Out1 "$sv\t3DUP\n";
  } else {
     print Out1 "$sv\tDUP\n";
  }
} elsif ($idenrev==1) {
  if ($orientrev==5) {
     print Out1 "$sv\t5INVDUP\n";
  } elsif ($orientrev==3) {
     print Out1 "$sv\t3INVDUP\n";
  } else {
     print Out1 "$sv\tINVDUP\n";
  }
} else {
    print Out1 "$sv\tNONDUP\n";
}

}

}

#########################
sub Perc {
my ($input,$seqtype) = @_;
my %chromosomes;
my %total_lengths;
open(my $fh, "<", "$input") or die "Cannot open input file: $!";
open OUT2, ">$input.perc.bed" or die "$!";

while (my $line = <$fh>) {
    chomp $line;
    my ($chromosome, $length, $start, $end, $strand, $repeat_name) = split(/\t/, $line);
    my $repeat_type;
    if ($seqtype=~/TE/) {
    ($repeat_type) = $repeat_name =~ /#(.+?)\//;
    } else {
    ($repeat_type) = $repeat_name;
    }
 
    $chromosomes{$chromosome}{$repeat_type} += ($end - $start);
    $total_lengths{$chromosome} = $length;
}
close($fh);

foreach my $chromosome (sort keys %chromosomes) {
    my $total_length = $total_lengths{$chromosome};
    print OUT2 "$chromosome\t(";
    foreach my $repeat_type (sort keys %{$chromosomes{$chromosome}}) {
        my $length = $chromosomes{$chromosome}{$repeat_type};
        my $percentage = $length / $total_length * 100;
        my $formatted_ratio = sprintf("%.1f", $percentage);
         my @remaining_keys = sort keys %{$chromosomes{$chromosome}};
         my $last_key = $remaining_keys[-1];
         if ($repeat_type eq $last_key) {
           print OUT2 "$repeat_type:$formatted_ratio%";
         } else {
           print OUT2 "$repeat_type:$formatted_ratio%;";
         }
    }
    print OUT2 ")";
    print OUT2 "\n";
}

}
##########################
sub GenomicOverlap {
 my ($input_line) = @_;
 my %last_line;
 open Input, "$input_line" or die "$!";
 open OUTfile0, ">$input_line.bed" or die "$!";
 while (my $line = <Input>) {
    my ($chromosome, $start, $end) = (split /\s+/, $line)[0, 2, 3];
    my $overlap = 0;

    if (exists $last_line{$chromosome}) {
        my ($last_start, $last_end) = @{$last_line{$chromosome}}{qw(start end)};
        $overlap = calculate_overlap($start, $end, $last_start, $last_end);
    }

    if ($overlap <= 0.1) {
        print OUTfile0 $line;
        $last_line{$chromosome} = { start => $start, end => $end };
    }
 }
}

###################################
sub calculate_overlap {
    my ($start1, $end1, $start2, $end2) = @_;
    my $overlap = 0;

    if ($start1 <= $end2 && $start2 <= $end1) {
        my $max_start = $start1 > $start2 ? $start1 : $start2;
        my $min_end = $end1 < $end2 ? $end1 : $end2;
        $overlap = ($min_end - $max_start + 1) / ($end1 - $start1 + 1);
    }

    return $overlap;
}

###################################
sub convert_to_fasta {
    my ($input_file, $output_file, $type) = @_;
    open(my $input_fh, '<', $input_file) or die "Unable to open the input file: $!";
    open(my $output_fh, '>', $output_file) or die "Unable to open the output file: $!";

    while (my $line = <$input_fh>) {
        chomp $line;
        next if $line =~ /^#/;
        my ($col1, $col2,$col3, $col5);
        if ($type eq "INS") {
        ($col1, $col2, $col3, $col5) = (split /\t/, $line)[0, 1, 2, 4];
        } elsif ($type eq "DEL") {
        ($col1, $col2, $col3, $col5) = (split /\t/, $line)[0, 1, 2, 3];
        }
        print $output_fh ">$col1\_$col2\_$col3\n$col5\n";
     }
        close $input_fh;                                                                                                                                                          
        close $output_fh;                                                                                                                                                        
}

#####

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
