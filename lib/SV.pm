package lib::SV;
use strict;
use warnings;

use Exporter qw(import);
our @EXPORT_OK = qw(ParseFill ExtractSeq SynNet);

###################### Pare the fill lines in the SynNet file ##################################
sub ParseFill {
my ($line) = @_;
my %hash;
my @temp = split(/_/,$line);
if ($line=~/inv|syn/) {
$hash{"target_beg"} =$temp[1];
$hash{"target_len"} =$temp[2];
$hash{"query_name"} =$temp[3];
$hash{"query_orient"}=$temp[4];
$hash{"query_beg"}=$temp[5];
$hash{"query_len"}=$temp[6];
$hash{"id"}=$temp[8];
$hash{"score"}=$temp[10];
$hash{"ali"}=$temp[12];
$hash{"qOver"}=$temp[14];
$hash{"qFar"}=$temp[16];
$hash{"qDup"}=$temp[18];
$hash{"type"}=$temp[20];
} elsif ($line=~/nonSyn/) {
$hash{"target_beg"} =$temp[1];
$hash{"target_len"} =$temp[2];
$hash{"query_name"} =$temp[3];
$hash{"query_orient"}=$temp[4];
$hash{"query_beg"}=$temp[5];
$hash{"query_len"}=$temp[6];
$hash{"id"}=$temp[8];
$hash{"score"}=$temp[10];
$hash{"ali"}=$temp[12];
$hash{"qDup"}=$temp[14];
$hash{"type"}=$temp[16];
}
my $ref=\%hash;
return $ref;
}

### End


###################### Load sequence into a Hash ###############################
sub ExtractSeq {
   my ($file) = @_;
   my %hash;
   $/=">";
   open In, "$file" or die "$!";
      while (<In>) {
         next if (length $_ < 2);
         my @unit= split ("\n", $_);
         my $head0=shift @unit;
         my @head1= split (/\s+/,$head0);
         my $head= shift @head1;
         my $seq = join ("", @unit);
         $hash{$head}=$seq;
      }
     my $ref=\%hash;
     return $ref;
 }

### End


################################## Load .synnet file into a Tree topology ######################
sub SynNet {
$/="net";
my $net_file = shift;
open In, "$net_file" || die "Fail open.$net_file";
my $root = Tree::DAG_Node->new; # Root
   $root->name('Root');
   <In>;
while (<In>) {
       if ($_=~/^#/) {
           next;
       } else {
           my @temp = split (/\n/,$_);
           my $head = shift @temp;
           my @unit = split (/\s/,$head);
           my $a1 = Tree::DAG_Node->new;
              $root ->add_daughter($a1); # Chromosome level.
              $a1->name($unit[1]); # Chromosome name.
              print "Processing chromosome $unit[1]\n";
           
       my $i;   
	   my @fill_nodeA;
	   my @gap_nodeA;
	   my @fill_nodeB;
       my @gap_nodeB;
       my @fill_nodeC;
       my @gap_nodeC;
       
       for ($i=0;$i<=$#temp;$i++) {
       
                  if ($temp[$i] =~/^\s{1}fill/) {
                           my @fill=split(/\s/,$temp[$i]);
		                   shift @fill;
                           my $name=join("_",@fill); # Chromosome_Name begin end.
                           my $a2 = Tree::DAG_Node->new; # Fill playerA: Generally, each chromsome may have one or more this level fill, but only one cover most of the chromosome region.
                           $a1->add_daughter($a2);
                           $a2->name($name); # Chromosome_Name begin end.
	                       push(@fill_nodeA,$a2); # record the $a2 which is the ancentor of the .ap?for the main SVs
                  } elsif ($temp[$i] =~/^\s{2}gap/) {
                           my @gap = split(/\s+/,$temp[$i]);
		                   shift @gap; # Remove the first element of @gap, It should be a blank space.
                           my $name=join("_",@gap); # Most simple insertions and deletions will be identified in this player.
                           my $a3 = Tree::DAG_Node->new;
                           $fill_nodeA[$#fill_nodeA]->add_daughter($a3);
                           $a3->name($name); # These names contain the major information of simple insertions and deletions.
			               push (@gap_nodeA,$a3); # record the $a3 (gap) which is the ancentor of complex SVs, including inversions, local/tandem duplications.
                  } elsif ($temp[$i] =~/^\s{3}fill/) {
			               my @fill = split(/\s+/,$temp[$i]);
	                       shift @fill; # Remove the first element of @gap, it should be a blank space.
                           my $name=join("_",@fill); # Most complex SVs will be identified in this player.
			               my $a4 = Tree::DAG_Node->new;
                           $gap_nodeA[$#gap_nodeA]->add_daughter($a4);
                           $a4->name($name); # These names contain the major information of simple insertions and deletions.
			               push (@fill_nodeB,$a4); # record the $a4 (fill) which is the ancentor of complex SVs, including inversions, local/tandem duplications.
		          } elsif ($temp[$i] =~/^\s{4}gap/) {
			               my @gap = split(/\s+/,$temp[$i]);
		                   shift @gap; # Remove the first element of @gap, it should be a blank space.
                           my $name=join("_",@gap); # SVs in duplication regions will be identified in this player.
		                   my $a5 = Tree::DAG_Node->new;
                           $fill_nodeB[$#fill_nodeB]->add_daughter($a5);
                           $a5->name($name); # These names contain the major information of simple insertions and deletions.
			               push (@gap_nodeB,$a5); # record the $a5 (gap) which can be used to infer the mosaic patter of segmental duplications.
		          } elsif ($temp[$i] =~/^\s{5}fill/) {
			               my @fill = split(/\s+/,$temp[$i]);
	                       shift @fill; # Remove the first element of @gap, it should be a blank space.
                           my $name=join("_",@fill); # SVs in duplication regions will be identified in this player.
		                   my $a6 = Tree::DAG_Node->new;
                           $gap_nodeB[$#gap_nodeB]->add_daughter($a6);
                           $a6->name($name); # These names contain the major information of simple insertions and deletions.
		                   push (@fill_nodeC,$a6); # record the $a6 (gap) which can be used to infer the mosaic patter of segmental duplications.	
		          } elsif ($temp[$i] =~/^\s{6}gap/) {
		                   my @gap = split(/\s+/,$temp[$i]);
		                   shift @gap; # Remove the first element of @gap, it should be a blank space.
                           my $name=join("_",@gap); # SVs in duplication regions will be identified in this player.
		                   my $a7 = Tree::DAG_Node->new;
                           $fill_nodeC[$#fill_nodeC]->add_daughter($a7);
                           $a7->name($name); # These names contain the major information of simple insertions and deletions.
			               push (@gap_nodeC,$a7); # record the $a7 (gap) which can be used to infer the mosaic patter of segmental duplications.
		          
		          }
           
              }
         }
    }
return $root;
}
### End

1;
