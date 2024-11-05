use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use POSIX qw(strftime);
use File::Basename qw(basename dirname);
use List::Util 'first';
use lib "$FindBin::Bin/./lib";
use Parallel::ForkManager;
use Pod::Usage;
use Env;
#$UCSC ||="$Bin/pub/UCSC/";
#$stretcher ||=".";
my $library_path = "$Bin/pub/Stretcher";
$ENV{'LD_LIBRARY_PATH'} = "$library_path:$ENV{'LD_LIBRARY_PATH'}";
$ENV{'EMBOSS_ACDROOT'} = "$library_path";
print "ENV lib:".$ENV{'LD_LIBRARY_PATH'}."\n";
my $stretcher ||= "$Bin/pub/Stretcher/";
`${stretcher}stretcher -h`
