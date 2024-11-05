use strict;
use warnings;
$/="\n>";
open In, "$ARGV[0]" or die "$!";
open Out, ">$ARGV[0].subcontigs.fa" or die "$!";

my $bin = $ARGV[1];
my $slide= $ARGV[2];


while (<In>) {
next if ($_=~/^$/);
chomp $_;
my ($id) = $_ =~ /^>*(\S+)/;  # parse ID as first word in FASTA header
   $_ =~ s/^>*.+\n//;  # remove FASTA header
   $_ =~ s/\n//g;  # remove endlines

my $len = length ($_);

my $n = abs ($len/$bin);

my $i;

for ($i=0;$i<=$n;$i++) {
my $beg;
if ($i==0) {
$beg = 0;
} else {
$beg = $i*$bin-$slide;
}
my $readlen = $bin + $slide;
my $seqrange = substr $_, $beg, $readlen;
my $start = $beg + 1;
my $end;
if ($i>$n-1) {
$end = $len;
} else {
$end = $beg + $readlen;
}
print Out ">$id\_$start\_$end\n$seqrange\n";
}

}

