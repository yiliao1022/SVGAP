#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


####################################################################################
# SVCombine_raw.pl
#
# Authors: Yi Liao (12/20/2022)
#
# Usage: Merge raw SVs from pairwise comparison files
#
####################################################################################
my ($pairsv,$refname,$outdir,$ins_tol,$del_ov,$candidateNum,$Help);

GetOptions(
  "sv=s"=>\$pairsv,
  "refname=s"=>\$refname,
  "outdir=s"=>\$outdir,
  "ins_tol=f"=>\$ins_tol,
  "del_ov=f"=>\$del_ov,
  "cant=i"=>\$candidateNum,
  "help"     =>\$Help
);

$ins_tol ||= 13;
$del_ov ||= 0.90;
$candidateNum ||= 20;
$outdir ||="CombinedSV";
$outdir =~ s/\/$//;
mkdir($outdir) unless(-d $outdir);

##################################################################################################
if ($Help){
print <<"END.";
  Usage: perl $0 --sv ./rice/CleanSV/ --refname MH63
  --sv          folder kept pairwise SVs output             [REQUIRED]
  --refname     reference name [e.g. ISO1,B73, et al...]    [REQUIRED]
  --outdir      output directory [default: CombinedSV]      [Optional]
  --ins_tol     the maximum distance tolerance for two insertion events that are considered or tested if they are the same event [default : 13]
  --del_ov      the minimum overlap proportion for two deletion events that are considered or tested if they are the same event [default : 0.90]
  --cant        How many neighbor SV events used to test together if they belong to the same event, should increase when your sample number is large [default : 20]
  --help        print this help information
END.
exit;
}

my $message = "Try 'perl $0 --help' for more information.";

if (! $pairsv){
    pod2usage( {
           -message => $message,
           -verbose => 0,
           -exitval => 1 } );
    }

###############################################

system ("cat $pairsv/$refname.*.rawSV/*.indel.sort | grep -v '\#' | awk '\$8>49' | awk '\$8<100000' | sort --buffer-size=10G -k1,1 -k2,2n > $outdir/All.DELs.50bplarge.bed");
system ("cat $pairsv/$refname.*.rawSV/*.indel.sort | grep -v '\#' | awk '\$8<50' | sort --buffer-size=10G -k1,1 -k2,2n > $outdir/All.DELs.50bpsmall.bed");
system ("cat $pairsv/*.$refname.rawSV/*.indel.sort | grep -v '\#' | awk '\$8>49' | awk '\$8<100000' | sort --buffer-size=10G -k5,5 -k6,6n > $outdir/All.INTs.50bplarge.bed");
system ("cat $pairsv/*.$refname.rawSV/*.indel.sort | grep -v '\#' | awk '\$8<50' | sort --buffer-size=10G -k5,5 -k6,6n > $outdir/All.INTs.50bpsmall.bed");
system ("cat $pairsv/*.$refname.rawSV/*.cnv.sort | grep -v '\#' | sort --buffer-size=10G -k5,5 -k6,6n > $outdir/All.CNVs.bed");
system ("cat $pairsv/$refname.*.rawSV/*.inv.sort | grep -v '\#' | sort --buffer-size=10G -k1,1 -k2,2n > $outdir/All.INVs.bed");
system ("cat $pairsv/*.$refname.rawSV/*.trl.sort | grep -v '\#' | sort --buffer-size=10G -k5,6 -k2,2n > $outdir/All.TRLs.bed");

&Insertion ("$outdir/All.INTs.50bplarge.bed",$ins_tol);
&Insertion ("$outdir/All.INTs.50bpsmall.bed",$ins_tol);
#&Insertion ("$outdir/All.CNVs.bed",$ins_tol);
#&Insertion ("$outdir/All.TRLs.bed",$ins_tol);
&Deletions ("$outdir/All.DELs.50bplarge.bed",$del_ov);
&Deletions ("$outdir/All.DELs.50bpsmall.bed",$del_ov);
#&Deletions ("$outdir/All.INVs.bed",$del_ov);

############  Combine Deletions

sub Deletions {
my ($del,$overlap_cutoff) = @_;
`sort --buffer-size=10G -k1,1 -k2,2n -k3,3n -k8,8n $del > $del.sort.txt`;
my %hash;
my $log="NA_0_20000000000000_0_0_0";
open DEL, "$del.sort.txt" or die "$!";
open DELOUT, ">$del.sort.merge.txt" or die "$!";

while (<DEL>) {
  chomp;
  my @tmp = split (/_/,$log);
  my ($ref_chr,$ref_s,$ref_e,$alt_chr,$alt_s,$alt_e,$len,$five,$three,$refseq,$altseq) = (split(/\t/,$_))[0,1,2,4,5,6,7,9,10,15,16];
  my $DEL = join ("_",($ref_chr,$ref_s,$ref_e,$len,$five,$three,$refseq,$altseq));
  my $overlap = &Overlap($ref_s,$ref_e,$tmp[1],$tmp[2]);
  if (($ref_chr eq $tmp[0]) and $overlap>$overlap_cutoff and $len/$tmp[3]>0.95 and $len/$tmp[3]<1.05) {
  my $species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$log}},$species;
  } else {
  $log = $DEL;
  my $species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
  push @{$hash{$log}},$species;
  }
}

my $head = "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tMERGE_SAMPLES_NUM\tMERGE_SAMPLES\tQUERY\tQ_POS\tQ_END\tFIVE\tTHREE\tREFseq\tALTseq\n";

foreach my $ind (keys %hash) {
      my @tmp = split (/\_/,$ind);
      my $id = ${$hash{$ind}}[0];
      my ($alt,$alt_s,$alt_e) = (split(/_/,$id))[0,1,2];
      my $num = $#{$hash{$ind}} + 1;
      print DELOUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$id\tDEL\t$tmp[3]\t$num\t";

if ($#{$hash{$ind}} == 0) {
      print DELOUT "${$hash{$ind}}[0]\t";
} else {
  my $last_one = pop @{$hash{$ind}};
  foreach my $svind (@{$hash{$ind}}) {
      print DELOUT "$svind,";
  }
      print DELOUT "$last_one\t";
}
      print DELOUT "$alt\t$alt_s\t$alt_e\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\n";
}

`sort --buffer-size=10G -k1,1 -k2,2n $del.sort.merge.txt > $del.combined.sorted.txt`;
`sed -i '1i $head' $del.combined.sorted.txt`;
`rm $del.sort.txt $del.sort.merge.txt`;
}

###### Combine insertions
sub Insertion {
my ($ins,$tol) = @_;
`sort --buffer-size=10G -k5,5 -k6,6n -k7,7n -k8,8n $ins > $ins.sorted.txt`;
my %hash;
my @content;
my $log="NA_0_0_0_NNNNNN_0_0";
push (@content,$log);
open INS, "$ins.sorted.txt" or die "No insertion file is provided!";
open OUT, ">$ins.sorted.combined.txt" or die "Unable to write  output";

while (<INS>) {
chomp;
my ($ref_chr,$ref_s,$ref_e,$alt_chr,$alt_s,$alt_e,$len,$five,$three,$seq,$altseq,$strand) = (split(/\t/,$_))[4,5,6,0,1,2,7,9,10,16,15,8];
my $sequence; ### = uc ($seq);
if ($strand =~/\-/) {
my $tmpseq = &RevCom($seq);
$sequence = uc ($tmpseq);
} else {
$sequence = uc ($altseq);
}


my $SV = join ("_",($ref_chr,$ref_s,$ref_e,$len,$sequence,$five,$three,$seq,$altseq));

 my $i;
 for ($i=0;$i<=$#content;$i++) {
 my $Comparion = &SeqIdentity($content[$i],$SV,$tol);
    if ($Comparion ==1) {
        my $Species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
        push @{$hash{$content[$i]}},$Species;
        last;
    }

    if ($i == $#content) {
       if ($#content < $candidateNum) {
           unshift (@content,$SV);
       } else {
           pop (@content);
           unshift (@content,$SV);
       }
      my $Species = join ("_",($alt_chr,$alt_s,$alt_e,$len));
      push @{$hash{$SV}},$Species;
      last;
    }

  }

}

my $head = "#CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN\tMERGE_SAMPLES_NUM\tMERGE_SAMPLES\tQUERY\tQ_POS\tQ_END\tFIVE\tTHREE\tREFseq\tALTseq\n";
foreach my $ind (keys %hash) {
      my @tmp = split (/\_/,$ind);
      my $id = ${$hash{$ind}}[0];
      my ($alt,$alt_s,$alt_e) = (split(/_/,$id))[0,1,2];
      my $num = $#{$hash{$ind}} + 1;
      print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$id\tINS\t$tmp[3]\t$num\t";

if ($#{$hash{$ind}} == 0) {
      print OUT "${$hash{$ind}}[0]\t";
} else {
  my $last_one = pop @{$hash{$ind}};
  foreach my $svind (@{$hash{$ind}}) {
      print OUT "$svind,";
  }
      print OUT "$last_one\t";
}
      print OUT "$alt\t$alt_s\t$alt_e\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]\n"
}

`sort --buffer-size=10G -k1,1 -k2,2n $ins.sorted.combined.txt > $ins.combined.sorted.txt`;
`sed -i '1i $head' $ins.combined.sorted.txt`;
`rm $ins.sorted.txt $ins.sorted.combined.txt`;
}

###### Compare identity of the insertion sequences from difference genomes
sub SeqIdentity {
my ($SV1,$SV2,$dis) = @_;
my ($chr1,$beg1,$end1,$len1,$sequence1,$five1,$three1) = (split(/\_/,$SV1))[0,1,2,3,4,5,6];
my ($chr2,$beg2,$end2,$len2,$sequence2,$five2,$three2) = (split(/\_/,$SV2))[0,1,2,3,4,5,6];
my $ident=0;
#if($len1>length $sequence1){
#	print "SV1:".$SV1."\n";
#	print "info:".length($sequence1)."\n";
#        print "info:".$len1."\n";
#	exit;
#}
#print "info:".(length $sequence1)."\n";
#print "info:".$len1."\n";
if (($chr1 eq $chr2) && ($sequence1 eq $sequence2) && (abs($beg1-$beg2) < $dis)) {
  $ident=1;
} elsif (($chr1 eq $chr2) && (abs($beg1-$beg2) < $dis) && ($five1 == $five2 || $three1==$three2)) {
  $ident=1;
} elsif (($chr1 eq $chr2) &&  $len1 > 5 && $len1 < 50 && $len1/$len2 > 0.9 && $len1/$len2<1.1) {

  my @seeds;
  my $num = $len1-4+1;
  my $i;
  for ($i=0;$i<$num;$i++) {
   my $seed = substr ($sequence1,$i,4);
   push (@seeds,$seed);
  }

  my $m=0;
  foreach my $sed (@seeds) {
  if ($sequence2=~/$sed/) {
    $m++;
   }
  }

  my $perc = $m/$num;
   if ($perc > 0.5) {
     $ident =1;
   }

} elsif (($chr1 eq $chr2) && $len1 > 50 && (abs($beg1-$beg2) < $dis)) {
  my @seeds;
  my $num;
  my $i;
  my ($SEQ1,$SEQ2);
  if ($len1 > $len2) {
    $SEQ1 = $sequence2;
    $SEQ2 = $sequence1;
    $num = $len2-10+1;
  } else {
    $SEQ1 = $sequence1;
    $SEQ2 = $sequence2;
    $num = $len1-10+1;
  }

  for ($i=0;$i<$num;$i++) {
  my $seed = substr ($SEQ1,$i,10);
  push (@seeds,$seed);
  }

  my $m=0;
  foreach my $sed (@seeds) {
  if ($SEQ2=~/$sed/) {
    $m++;
   }
  }

  my $perc = $m/$num;
  if ($perc > 0.5) {
    $ident =1;
  }
}
return $ident;
}

############ Calculate the overlap between two coordinate intervals
sub Overlap {
my ($c1,$c2,$c3,$c4) = @_;
my $overlap;

if ($c1>$c4) {
  $overlap =0;
} elsif ($c1>=$c3 and $c1<=$c4 and $c2 <=$c4) {
  $overlap = ($c2-$c1)/($c4-$c3);
} elsif ($c1>=$c3 and $c1<=$c4 and $c2 >$c4) {
  $overlap = ($c4-$c1)/($c4-$c3);
}
return $overlap;
}

############################
sub RevCom {
my ($seq) = @_;
$seq =~ tr/natcguNATCGU/ntagcaNTAGCA/;
my $revseq = reverse $seq;
return $revseq;
}


