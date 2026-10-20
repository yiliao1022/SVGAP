#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use Cwd qw(abs_path);

my $root = abs_path(File::Spec->catdir($Bin, '..'));
my @required = (
    [qw(pub UCSC axtChain)],
    [qw(pub UCSC axtToMaf)],
    [qw(pub UCSC chainNet)],
    [qw(pub UCSC chainPreNet)],
    [qw(pub UCSC faSize)],
    [qw(pub UCSC faToTwoBit)],
    [qw(pub UCSC mafToAxt)],
    [qw(pub UCSC netSyntenic)],
    [qw(pub UCSC netToAxt)],
    [qw(pub TBA single_cov2)],
    [qw(pub MINIMAP2 k8)],
    [qw(pub MINIMAP2 paftools.js)],
    [qw(pub MUMMER delta2maf)],
);

my @missing;
for my $parts (@required) {
    my $path = File::Spec->catfile($root, @$parts);
    push @missing, $path unless -f $path && -x $path;
}

if (@missing) {
    die "ERROR: missing or non-executable bundled dependency:\n  "
        . join("\n  ", @missing) . "\n";
}

print "SVGAP dependency check: OK\n";
print "Package root: $root\n";
print "Parallel workers: core Perl fork/wait (no CPAN module required)\n";
exit 0;
