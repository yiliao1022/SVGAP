#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Fcntl qw(SEEK_SET);

my %opt=(vcf=>'',ref=>'',out=>'',allow_symbolic_inv=>1,help=>0);
GetOptions(
  'vcf=s'=>\$opt{vcf}, 'ref|reference=s'=>\$opt{ref}, 'out=s'=>\$opt{out},
  'allow_symbolic_inv|allow-symbolic-inv!'=>\$opt{allow_symbolic_inv},
  'help|h'=>\$opt{help},
) or die usage();
if($opt{help}||!$opt{vcf}||!$opt{ref}){print usage();exit($opt{help}?0:1)}
$opt{out}||="$opt{vcf}.validation.tsv";
die "ERROR: VCF not found: $opt{vcf}\n" unless -s $opt{vcf};
die "ERROR: reference not found: $opt{ref}\n" unless -s $opt{ref};
my $db=VGFasta->new($opt{ref});
open(my $IN,'<',$opt{vcf}) or die $!;
open(my $OUT,'>',$opt{out}) or die $!;
print $OUT join("\t",qw(LINE CHROM POS ID STATUS MESSAGE)),"\n";
my($line_no,$records,$errors,$warnings)=(0,0,0,0);
my($prev_chr,$prev_pos); my %ids;
while(my $line=<$IN>){$line_no++;chomp$line;next if $line=~/^#/||$line eq'';my@f=split(/\t/,$line,-1);if(@f<8){issue('ERROR','fewer than 8 VCF columns','.','.','.');next}$records++;my($chr,$pos,$id,$ref,$alt,$info_s)=@f[0,1,2,3,4,7];
  if(!$db->has_contig($chr)){issue('ERROR','contig absent from reference',$chr,$pos,$id);next}
  if($pos!~/^\d+$/||$pos<1||$pos>$db->length($chr)){issue('ERROR','POS outside contig',$chr,$pos,$id);next}
  if($ids{$id}++){issue('ERROR','duplicate VCF ID',$chr,$pos,$id)}
  if(defined$prev_chr && $chr eq $prev_chr && $pos<$prev_pos){issue('ERROR','VCF is not position sorted within contig',$chr,$pos,$id)}
  ($prev_chr,$prev_pos)=($chr,$pos);
  my%info=parse_info($info_s);my$type=defined$info{SVTYPE}?uc($info{SVTYPE}):$alt=~/^<([^>]+)>$/?uc($1):'';
  if($alt!~/^</){
    if($ref!~/^[ACGTN]+$/i||$alt!~/^[ACGTN]+$/i){issue('ERROR','non-DNA explicit allele',$chr,$pos,$id);next}
    my$actual=$db->get_seq($chr,$pos-1,$pos-1+CORE::length($ref));
    if(!defined$actual||uc($actual)ne uc($ref)){issue('ERROR','REF does not match reference FASTA',$chr,$pos,$id)}
    my$expect_end=$pos+CORE::length($ref)-1;
    if(defined$info{END}&&$info{END}=~/^\d+$/&&$info{END}!=$expect_end){issue('ERROR',"END=$info{END}, expected $expect_end from REF span",$chr,$pos,$id)}
    my$expect_svlen=CORE::length($alt)-CORE::length($ref);
    if(defined$info{SVLEN}&&$info{SVLEN}=~/^-?\d+$/&&$info{SVLEN}!=$expect_svlen){issue('ERROR',"SVLEN=$info{SVLEN}, expected $expect_svlen from ALT-REF",$chr,$pos,$id)}
  }elsif($alt eq '<INV>'&&$opt{allow_symbolic_inv}){
    if($ref!~/^[ACGTN]$/i){issue('ERROR','symbolic INV REF must be one real reference base',$chr,$pos,$id)}else{my$b=$db->get_seq($chr,$pos-1,$pos);issue('ERROR','symbolic INV REF mismatch',$chr,$pos,$id)if uc($b)ne uc($ref)}
    if(!defined$info{END}||$info{END}!~/^\d+$/||$info{END}<$pos){issue('ERROR','symbolic INV lacks valid END',$chr,$pos,$id)}else{my$expect=$info{END}-$pos+1;if(defined$info{SVLEN}&&$info{SVLEN}=~/^-?\d+$/&&$info{SVLEN}!=$expect){issue('ERROR',"INV SVLEN=$info{SVLEN}, expected $expect",$chr,$pos,$id)}}
  }else{issue('ERROR',"unsupported unresolved symbolic allele $alt",$chr,$pos,$id)}
}
close$IN;
print $OUT "#SUMMARY\trecords=$records;errors=$errors;warnings=$warnings\n";
close$OUT;
print STDERR "Validated $records records: errors=$errors warnings=$warnings report=$opt{out}\n";
exit($errors?2:0);
sub issue{my($status,$msg,$chr,$pos,$id)=@_;$errors++ if$status eq'ERROR';$warnings++ if$status eq'WARNING';print $OUT join("\t",$line_no,$chr,$pos,$id,$status,$msg),"\n"}
sub parse_info{my($s)=@_;my%h;for my$x(split(/;/,$s//'')){next if$x eq''||$x eq'.';my($k,$v)=split(/=/,$x,2);$h{$k}=defined$v?$v:undef}return%h}
sub usage{return <<'USAGE';
Usage:
  perl 14_validate_vg_vcf.pl --vcf All.SVGAP.vg.vcf --ref ISO1.fa --out validation.tsv

Checks explicit REF against the reference FASTA, END/SVLEN consistency,
unique IDs, coordinate sorting and supported symbolic <INV> records.
USAGE
}
package VGFasta;
use strict;use warnings;use Fcntl qw(SEEK_SET);
sub new{my($class,$fa)=@_;die"ERROR: uncompressed FASTA required\n"if$fa=~/\.gz$/i;my$self=bless{fa=>$fa,idx=>{}},$class;if(-s"$fa.fai"){$self->_load("$fa.fai")}else{$self->_build}open(my$fh,'<',$fa)or die$!;binmode$fh;$self->{fh}=$fh;return$self}
sub _load{my($s,$f)=@_;open(my$I,'<',$f)or die$!;while(<$I>){chomp;my@x=split/\t/;next unless@x>=5;$s->{idx}{$x[0]}={len=>$x[1]+0,off=>$x[2]+0,lb=>$x[3]+0,lw=>$x[4]+0}}close$I}
sub _build{my($s)=@_;open(my$I,'<',$s->{fa})or die$!;binmode$I;my($n,$len,$off,$lb,$lw);while(1){my$p=tell$I;my$l=<$I>;last unless defined$l;if($l=~/^>(\S+)/){$s->{idx}{$n}={len=>$len,off=>$off,lb=>$lb,lw=>$lw}if defined$n;$n=$1;($len,$off,$lb,$lw)=(0,undef,undef,undef)}else{next unless defined$n;my$raw=CORE::length$l;$l=~s/[\r\n]+$//;my$b=CORE::length$l;next if$b==0;$off=$p unless defined$off;$lb=$b unless defined$lb;$lw=$raw unless defined$lw;$len+=$b}}$s->{idx}{$n}={len=>$len,off=>$off,lb=>$lb,lw=>$lw}if defined$n;close$I}
sub has_contig{exists$_[0]{idx}{$_[1]}}sub length{$_[0]{idx}{$_[1]}{len}}
sub get_seq{my($s,$c,$st,$en)=@_;return undef unless $s->has_contig($c);my$x=$s->{idx}{$c};return undef if $st<0||$en>$x->{len};return '' if $en<=$st;my$need=$en-$st;my$byte=$x->{off}+int($st/$x->{lb})*$x->{lw}+($st%$x->{lb});seek($s->{fh},$byte,SEEK_SET)or die $!;my$q='';while(CORE::length($q)<$need){my$b='';my$n=read($s->{fh},$b,$need-CORE::length($q)+100);last unless$n;$b=~s/[\r\n]//g;$q.=$b}return uc(substr($q,0,$need))}
1;
