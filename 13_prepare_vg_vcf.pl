#!/usr/bin/env perl
# Author: Yi Liao @ scau.edu.cn
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use lib "$Bin/lib";
use lib "$Bin/pub";
use lib "$Bin/../lib";
use lib "$Bin/../pub";
use File::Basename qw(basename);
use POSIX qw(strftime);

################################################################################
# 13_prepare_vg_vcf.pl
#
# Convert/merge SVGAP component VCFs into a coordinate-consistent VCF suitable
# for vg construct.  The main rules are:
#   * SVGAP BED/raw-SV intervals are 0-based half-open.
#   * Explicit DEL alleles are rebuilt from the reference with a left anchor.
#   * Explicit INS/CXSV alleles are normalized and checked against the reference.
#   * Symbolic INV is retained with a real REF base and consistent POS/END/SVLEN.
#   * Unresolved symbolic INS/CXSV records are skipped and reported rather than
#     silently producing a graph with unknown sequence.
################################################################################

my %opt = (
    ref                 => '',
    vcfs                => '',
    out                 => '',
    report              => '',
    query_list          => '',
    glist               => '',
    refname             => '',
    refsize             => '',
    missing_gt          => './.',
    max_explicit_allele => 5_000_000,
    keep_symbolic_inv   => 1,
    keep_unresolved     => 0,
    help                => 0,
);

GetOptions(
    'ref|reference|ref_fasta=s' => \$opt{ref},
    'vcfs=s'                    => \$opt{vcfs},
    'out=s'                     => \$opt{out},
    'report=s'                  => \$opt{report},
    'query_list|query-list=s'   => \$opt{query_list},
    'glist=s'                   => \$opt{glist},
    'refname|target=s'          => \$opt{refname},
    'refsize=s'                 => \$opt{refsize},
    'missing_gt=s'              => \$opt{missing_gt},
    'max_explicit_allele|max-explicit-allele=i' => \$opt{max_explicit_allele},
    'keep_symbolic_inv|keep-symbolic-inv!'       => \$opt{keep_symbolic_inv},
    'keep_unresolved|keep-unresolved!'           => \$opt{keep_unresolved},
    'help|h'                    => \$opt{help},
) or die usage();

if ($opt{help} || !$opt{ref} || !$opt{vcfs} || !$opt{out}) {
    print usage();
    exit($opt{help} ? 0 : 1);
}
die "ERROR: reference FASTA does not exist or is empty: $opt{ref}\n" unless -s $opt{ref};
die "ERROR: --max_explicit_allele must be >= 1\n" unless $opt{max_explicit_allele} >= 1;
$opt{report} ||= "$opt{out}.report.tsv";

my $REF = FastFAI->new($opt{ref});
my @samples = read_sample_order($opt{query_list} || $opt{glist}, $opt{refname});
my %chr_rank = read_chr_rank($opt{refsize}, $REF);
my @vcfs = grep { defined($_) && $_ ne '' && -s $_ } split(/,/, $opt{vcfs});
die "ERROR: no readable VCFs in --vcfs\n" unless @vcfs;

my (@records, %meta, %id_seen);
my %stats = map { $_ => 0 } qw(INPUT WRITTEN REBUILT_DEL NORMALIZED_EXPLICIT FIXED_INV SKIPPED UNRESOLVED DUPLICATE_ID);
open(my $REPORT, '>', $opt{report}) or die "ERROR: cannot write $opt{report}: $!\n";
print $REPORT join("\t", qw(SOURCE CHROM OLD_POS ID SVTYPE ACTION NEW_POS NEW_END NEW_SVLEN REASON)), "\n";

for my $vcf (@vcfs) {
    open(my $IN, '<', $vcf) or die "ERROR: cannot read $vcf: $!\n";
    my (@file_samples, %sample_col);
    while (my $line = <$IN>) {
        chomp $line;
        next if $line eq '';
        if ($line =~ /^##/) {
            next if $line =~ /^##fileformat=/ || $line =~ /^##fileDate=/ || $line =~ /^##source=/ || $line =~ /^##contig=/;
            $meta{$line} = 1;
            next;
        }
        if ($line =~ /^#CHROM\t/) {
            my @h = split(/\t/, $line, -1);
            @file_samples = @h > 9 ? @h[9..$#h] : ();
            for my $i (0..$#file_samples) { $sample_col{$file_samples[$i]} = 9 + $i; }
            if (!@samples) { @samples = @file_samples; }
            next;
        }

        my @f = split(/\t/, $line, -1);
        next unless @f >= 8;
        $stats{INPUT}++;
        my ($chrom,$old_pos,$id,$ref,$alt,$qual,$filter,$info_s) = @f[0..7];
        my $format = @f > 8 ? $f[8] : 'GT';
        my %info = parse_info($info_s);
        my $svtype = infer_svtype($alt, \%info, $ref);
        my @out_gt = remap_samples(\@f, \%sample_col, \@samples, $format, $opt{missing_gt});

        my ($ok, $new, $action, $reason) = normalize_record(
            chrom=>$chrom, pos=>$old_pos, id=>$id, ref=>$ref, alt=>$alt,
            qual=>$qual, filter=>$filter, info=>\%info, svtype=>$svtype,
            format=>$format, gt=>\@out_gt, refdb=>$REF, opt=>\%opt,
        );

        if (!$ok) {
            $stats{SKIPPED}++;
            $stats{UNRESOLVED}++ if $action eq 'UNRESOLVED';
            print $REPORT join("\t", $vcf, $chrom, $old_pos, safe($id), $svtype, $action, '.', '.', '.', $reason), "\n";
            next;
        }

        my $base_id = sanitize_id($new->{id}, $new->{chrom}, $new->{pos}, $new->{svtype});
        if ($id_seen{$base_id}++) {
            $stats{DUPLICATE_ID}++;
            $base_id .= '.' . $id_seen{$base_id};
        }
        $new->{id} = $base_id;
        push @records, $new;
        $stats{WRITTEN}++;
        $stats{REBUILT_DEL}++       if $action eq 'REBUILT_DEL';
        $stats{NORMALIZED_EXPLICIT}++ if $action eq 'NORMALIZED_EXPLICIT';
        $stats{FIXED_INV}++         if $action eq 'FIXED_INV';
        print $REPORT join("\t", $vcf, $chrom, $old_pos, safe($id), $new->{svtype}, $action,
            $new->{pos}, $new->{end}, defined($new->{svlen}) ? $new->{svlen} : '.', $reason), "\n";
    }
    close $IN;
}

@records = sort {
    (exists $chr_rank{$a->{chrom}} ? $chr_rank{$a->{chrom}} : 1e12) <=>
    (exists $chr_rank{$b->{chrom}} ? $chr_rank{$b->{chrom}} : 1e12)
    || $a->{chrom} cmp $b->{chrom}
    || $a->{pos} <=> $b->{pos}
    || $a->{end} <=> $b->{end}
    || $a->{id} cmp $b->{id}
} @records;

open(my $OUT, '>', $opt{out}) or die "ERROR: cannot write $opt{out}: $!\n";
print $OUT "##fileformat=VCFv4.2\n";
print $OUT "##fileDate=", strftime('%Y%m%d', localtime), "\n";
print $OUT "##source=SVGAP_13_prepare_vg_vcf\n";
print $OUT "##reference=$opt{ref}\n";
for my $chr ($REF->contigs) {
    print $OUT "##contig=<ID=$chr,length=", $REF->length($chr), ">\n";
}
print $OUT "##ALT=<ID=INV,Description=\"Inversion\">\n";
print $OUT "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n";
print $OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate consistent with REF span or symbolic SV definition\">\n";
print $OUT "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"ALT length minus REF length excluding the shared anchor for explicit alleles\">\n";
print $OUT "##INFO=<ID=SVGAP_ORIG_POS,Number=1,Type=Integer,Description=\"Original SVGAP VCF POS before vg normalization\">\n";
print $OUT "##INFO=<ID=SVGAP_VG_READY,Number=0,Type=Flag,Description=\"Record normalized by 13_prepare_vg_vcf.pl\">\n";
for my $m (sort keys %meta) {
    next if $m =~ /^##ALT=<ID=(?:DEL|INS|CXSV|INV),/;
    next if $m =~ /^##INFO=<ID=(?:SVTYPE|END|SVLEN|SVGAP_ORIG_POS|SVGAP_VG_READY),/;
    print $OUT "$m\n";
}
print $OUT join("\t", '#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',@samples), "\n";
for my $r (@records) {
    print $OUT join("\t", @{$r}{qw(chrom pos id ref alt qual filter info_s format)}, @{$r->{gt}}), "\n";
}
close $OUT;

print $REPORT "#SUMMARY\t", join(';', map { "$_=$stats{$_}" } sort keys %stats), "\n";
close $REPORT;
print STDERR "vg-ready VCF: $opt{out}\n";
print STDERR "Report:       $opt{report}\n";
print STDERR join(' ', map { "$_=$stats{$_}" } sort keys %stats), "\n";
exit 0;

sub normalize_record {
    my %a = @_;
    my ($chrom,$old_pos,$id,$ref,$alt,$qual,$filter,$info,$svtype,$format,$gt,$db,$opt) =
        @a{qw(chrom pos id ref alt qual filter info svtype format gt refdb opt)};
    return (0,undef,'UNRESOLVED','non-numeric POS') unless defined($old_pos) && $old_pos =~ /^\d+$/ && $old_pos >= 1;
    return (0,undef,'UNRESOLVED','reference contig not found') unless $db->has_contig($chrom);
    my $orig_pos = $old_pos + 0;

    if ($alt =~ /^<DEL>$/i || ($svtype eq 'DEL' && $alt =~ /^</)) {
        my ($start0,$end0);
        if (defined($info->{BED_START}) && $info->{BED_START} =~ /^\d+$/ && defined($info->{END}) && $info->{END} =~ /^\d+$/) {
            ($start0,$end0) = ($info->{BED_START}+0, $info->{END}+0);
        } elsif (defined($info->{END}) && $info->{END} =~ /^\d+$/) {
            # Standard symbolic DEL: POS is the 1-based left anchor, so deleted
            # sequence is 0-based [POS, END).
            ($start0,$end0) = ($orig_pos, $info->{END}+0);
        } else {
            return (0,undef,'UNRESOLVED','symbolic DEL lacks END');
        }
        my ($pos,$new_ref,$new_alt,$end,$svlen,$err) = build_explicit_interval_alleles($db,$chrom,$start0,$end0,'');
        return (0,undef,'UNRESOLVED',$err) if $err;
        return (0,undef,'UNRESOLVED','explicit DEL exceeds --max_explicit_allele') if length($new_ref) > $opt->{max_explicit_allele};
        return finalize_new($chrom,$pos,$id,$new_ref,$new_alt,$qual,$filter,$info,'DEL',$end,$svlen,$format,$gt,$orig_pos,'REBUILT_DEL','rebuilt from BED_START/END and reference');
    }

    if ($alt =~ /^<INV>$/i || ($svtype eq 'INV' && $alt =~ /^</)) {
        return (0,undef,'UNRESOLVED','symbolic INV disabled') unless $opt->{keep_symbolic_inv};
        my $pos = (defined($info->{BED_START}) && $info->{BED_START} =~ /^\d+$/) ? $info->{BED_START}+1 : $orig_pos;
        my $end = defined($info->{END}) && $info->{END} =~ /^\d+$/ ? $info->{END}+0 : 0;
        return (0,undef,'UNRESOLVED','INV lacks valid END') unless $end >= $pos;
        my $base = $db->get_seq($chrom,$pos-1,$pos);
        return (0,undef,'UNRESOLVED','cannot fetch INV REF base') unless valid_seq($base);
        my $svlen = $end - $pos + 1;
        return finalize_new($chrom,$pos,$id,$base,'<INV>',$qual,$filter,$info,'INV',$end,$svlen,$format,$gt,$orig_pos,'FIXED_INV','kept symbolic INV with real REF base');
    }

    if ($alt =~ /^</) {
        if ($opt->{keep_unresolved}) {
            my $base = $db->get_seq($chrom,$orig_pos-1,$orig_pos);
            my $end = defined($info->{END}) && $info->{END} =~ /^\d+$/ ? $info->{END}+0 : $orig_pos;
            my $svlen = defined($info->{SVLEN}) && $info->{SVLEN} =~ /^-?\d+$/ ? $info->{SVLEN}+0 : 0;
            return finalize_new($chrom,$orig_pos,$id,$base,$alt,$qual,'UNRESOLVED',$info,$svtype,$end,$svlen,$format,$gt,$orig_pos,'KEPT_UNRESOLVED','kept by --keep_unresolved; may not construct in vg');
        }
        return (0,undef,'UNRESOLVED',"symbolic $alt has no explicit sequence");
    }

    $ref = uc($ref // '');
    $alt = uc($alt // '');
    return (0,undef,'UNRESOLVED','invalid explicit REF/ALT') unless valid_seq($ref) && valid_seq($alt);

    # Canonicalize SVGAP explicit indels from their original 0-based half-open
    # coordinates. This converts legacy right-padded alleles to a left anchor.
    my $raw_start_key = defined($info->{BED_START}) ? 'BED_START' : defined($info->{RAWPOS}) ? 'RAWPOS' : '';
    if (($svtype eq 'DEL' || $svtype eq 'INS') && $raw_start_key ne '' && defined($info->{END}) && $info->{END} =~ /^\d+$/ && $info->{$raw_start_key} =~ /^\d+$/) {
        my ($start0,$end0)=($info->{$raw_start_key}+0,$info->{END}+0);
        my $event_ref=$db->get_seq($chrom,$start0,$end0);
        if (defined $event_ref) {
            my $event_alt;
            if ($svtype eq 'DEL') { $event_alt=''; }
            else { $event_alt=derive_event_alt($ref,$alt,$event_ref); }
            if (defined $event_alt) {
                my ($p,$r,$al,$e,$sl,$err)=build_explicit_interval_alleles($db,$chrom,$start0,$end0,$event_alt);
                return (0,undef,'UNRESOLVED',$err) if $err;
                return (0,undef,'UNRESOLVED','explicit allele exceeds --max_explicit_allele') if length($r)>$opt->{max_explicit_allele} || length($al)>$opt->{max_explicit_allele};
                my $t = length($event_ref)==0 ? 'INS' : length($event_alt)==0 ? 'DEL' : $svtype;
                return finalize_new($chrom,$p,$id,$r,$al,$qual,$filter,$info,$t,$e,$sl,$format,$gt,$orig_pos,$svtype eq 'DEL' ? 'REBUILT_DEL' : 'NORMALIZED_EXPLICIT','rebuilt from raw SVGAP interval with canonical left anchor');
            }
        }
    }
    return (0,undef,'UNRESOLVED','explicit allele exceeds --max_explicit_allele') if length($ref) > $opt->{max_explicit_allele} || length($alt) > $opt->{max_explicit_allele};

    my ($pos,$nref,$nalt) = normalize_explicit($orig_pos,$ref,$alt);
    my $actual = $db->get_seq($chrom,$pos-1,$pos-1+length($nref));

    if (!defined($actual) || uc($actual) ne $nref) {
        # A common SVGAP legacy DEL has explicit or symbolic alleles but carries
        # BED_START/END from a 0-based half-open interval. Rebuild it exactly.
        if ($svtype eq 'DEL' && defined($info->{BED_START}) && $info->{BED_START} =~ /^\d+$/ && defined($info->{END}) && $info->{END} =~ /^\d+$/) {
            my ($p,$r,$al,$e,$sl,$err) = build_explicit_interval_alleles($db,$chrom,$info->{BED_START}+0,$info->{END}+0,'');
            return (0,undef,'UNRESOLVED',$err) if $err;
            return (0,undef,'UNRESOLVED','explicit DEL exceeds --max_explicit_allele') if length($r) > $opt->{max_explicit_allele};
            return finalize_new($chrom,$p,$id,$r,$al,$qual,$filter,$info,'DEL',$e,$sl,$format,$gt,$orig_pos,'REBUILT_DEL','explicit legacy DEL rebuilt from reference');
        }
        return (0,undef,'UNRESOLVED','REF does not match reference after normalization');
    }

    my $end = $pos + length($nref) - 1;
    my $svlen = length($nalt) - length($nref);
    my $new_type = infer_explicit_type($svtype,$nref,$nalt);
    return finalize_new($chrom,$pos,$id,$nref,$nalt,$qual,$filter,$info,$new_type,$end,$svlen,$format,$gt,$orig_pos,'NORMALIZED_EXPLICIT','explicit allele normalized and REF checked');
}

sub finalize_new {
    my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$svtype,$end,$svlen,$format,$gt,$orig_pos,$action,$reason) = @_;
    my %new_info = %$info;
    delete @new_info{qw(END SVLEN SVTYPE BED_START RAWPOS RAWEND)};
    $new_info{SVTYPE} = $svtype if $svtype ne 'SNP';
    $new_info{END} = $end;
    $new_info{SVLEN} = $svlen if $svtype ne 'SNP' || $svlen != 0;
    $new_info{SVGAP_ORIG_POS} = $orig_pos;
    $new_info{SVGAP_VG_READY} = undef;
    my $info_s = format_info(\%new_info);
    my $r = {
        chrom=>$chrom, pos=>$pos, id=>$id, ref=>$ref, alt=>$alt,
        qual=>(defined($qual)&&$qual ne ''?$qual:'.'),
        filter=>(defined($filter)&&$filter ne ''?$filter:'PASS'),
        info_s=>$info_s, format=>(defined($format)&&$format ne ''?$format:'GT'),
        gt=>$gt, end=>$end, svlen=>$svlen, svtype=>$svtype,
    };
    return (1,$r,$action,$reason);
}

sub build_explicit_interval_alleles {
    my ($db,$chrom,$start0,$end0,$alt_event) = @_;
    return (undef,undef,undef,undef,undef,'invalid 0-based interval') unless defined($start0) && defined($end0) && $start0 >= 0 && $end0 >= $start0;
    my $chr_len = $db->length($chrom);
    return (undef,undef,undef,undef,undef,'interval exceeds contig length') if $end0 > $chr_len;
    my $event_ref = $db->get_seq($chrom,$start0,$end0);
    return (undef,undef,undef,undef,undef,'cannot fetch reference interval') unless defined $event_ref;
    $alt_event = '' unless defined $alt_event;
    my ($pos,$ref,$alt);
    if ($start0 > 0) {
        my $anchor = $db->get_seq($chrom,$start0-1,$start0);
        return (undef,undef,undef,undef,undef,'cannot fetch left anchor') unless valid_seq($anchor);
        $pos = $start0; # 1-based position of the left anchor
        $ref = $anchor . $event_ref;
        $alt = $anchor . $alt_event;
    } else {
        return (undef,undef,undef,undef,undef,'cannot right-anchor event spanning entire contig') if $end0 >= $chr_len;
        my $anchor = $db->get_seq($chrom,$end0,$end0+1);
        return (undef,undef,undef,undef,undef,'cannot fetch right anchor') unless valid_seq($anchor);
        $pos = 1;
        $ref = $event_ref . $anchor;
        $alt = $alt_event . $anchor;
    }
    my $end = $pos + length($ref) - 1;
    my $svlen = length($alt) - length($ref);
    return ($pos,uc($ref),uc($alt),$end,$svlen,undef);
}


sub derive_event_alt {
    my ($ref,$alt,$event_ref)=@_;
    $ref=uc($ref); $alt=uc($alt); $event_ref=uc($event_ref // '');
    if (length($ref)==length($event_ref)+1) {
        if (substr($ref,1) eq $event_ref && substr($alt,0,1) eq substr($ref,0,1)) {
            return substr($alt,1);
        }
        if (substr($ref,0,length($ref)-1) eq $event_ref && substr($alt,-1,1) eq substr($ref,-1,1)) {
            return substr($alt,0,length($alt)-1);
        }
    }
    return $alt if $ref eq $event_ref;
    my ($r,$a)=($ref,$alt);
    while (length($r)>0 && length($a)>0 && substr($r,0,1) eq substr($a,0,1)) { $r=substr($r,1); $a=substr($a,1); }
    while (length($r)>0 && length($a)>0 && substr($r,-1,1) eq substr($a,-1,1)) { chop($r); chop($a); }
    return $a if $r eq $event_ref;
    return undef;
}

sub normalize_explicit {
    my ($pos,$ref,$alt) = @_;
    while (length($ref) > 1 && length($alt) > 1 && substr($ref,-1,1) eq substr($alt,-1,1)) {
        chop($ref); chop($alt);
    }
    while (length($ref) > 1 && length($alt) > 1 && substr($ref,0,1) eq substr($alt,0,1)) {
        $ref = substr($ref,1); $alt = substr($alt,1); $pos++;
    }
    return ($pos,$ref,$alt);
}

sub infer_explicit_type {
    my ($old,$ref,$alt) = @_;
    return 'SNP' if length($ref)==1 && length($alt)==1;
    return 'DEL' if length($ref) > length($alt) && length($alt)==1;
    return 'INS' if length($alt) > length($ref) && length($ref)==1;
    return $old && $old ne '' ? $old : 'CXSV';
}

sub infer_svtype {
    my ($alt,$info,$ref) = @_;
    return uc($info->{SVTYPE}) if defined($info->{SVTYPE}) && $info->{SVTYPE} ne '';
    return uc($1) if defined($alt) && $alt =~ /^<([^>:]+)(?::[^>]*)?>$/;
    return 'SNP' if defined($ref) && defined($alt) && $alt !~ /^</ && length($ref)==1 && length($alt)==1;
    return length($ref) > length($alt) ? 'DEL' : length($alt) > length($ref) ? 'INS' : 'CXSV';
}

sub parse_info {
    my ($s) = @_;
    my %h;
    return %h if !defined($s) || $s eq '' || $s eq '.';
    for my $x (split(/;/,$s)) {
        next if $x eq '';
        my ($k,$v) = split(/=/,$x,2);
        $h{$k} = defined($v) ? $v : undef;
    }
    return %h;
}

sub format_info {
    my ($h) = @_;
    my @first = grep { exists $h->{$_} } qw(SVTYPE END SVLEN SVGAP_ORIG_POS SVGAP_VG_READY);
    my %seen = map { $_=>1 } @first;
    my @keys = (@first, sort grep { !$seen{$_} } keys %$h);
    return join(';', map { defined($h->{$_}) ? "$_=$h->{$_}" : $_ } @keys) || '.';
}

sub remap_samples {
    my ($f,$sample_col,$samples,$format,$missing) = @_;
    my @out;
    for my $s (@$samples) {
        if (exists $sample_col->{$s} && defined($f->[$sample_col->{$s}])) {
            push @out, $f->[$sample_col->{$s}];
        } else {
            push @out, normalize_missing_gt($missing,$format);
        }
    }
    return @out;
}

sub normalize_missing_gt {
    my ($gt,$format) = @_;
    my @keys = split(/:/, defined($format)?$format:'GT');
    return $gt if @keys <= 1;
    return join(':', $gt, ('.') x (@keys-1));
}

sub read_sample_order {
    my ($file,$refname) = @_;
    return () unless defined($file) && $file ne '' && -s $file;
    open(my $IN,'<',$file) or die "ERROR: cannot read sample list $file: $!\n";
    my (@s,%seen);
    while (<$IN>) {
        chomp; s/^\s+|\s+$//g; next if $_ eq '' || /^#/;
        next if defined($refname) && $refname ne '' && $_ eq $refname;
        die "ERROR: duplicate sample in $file: $_\n" if $seen{$_}++;
        push @s,$_;
    }
    close $IN;
    return @s;
}

sub read_chr_rank {
    my ($file,$db) = @_;
    my %r; my $i=0;
    if (defined($file) && $file ne '' && -s $file) {
        open(my $IN,'<',$file) or die "ERROR: cannot read $file: $!\n";
        while (<$IN>) { chomp; next if /^\s*$/; my ($c)=split(/\s+/); $r{$c}=++$i unless exists $r{$c}; }
        close $IN;
    }
    for my $c ($db->contigs) { $r{$c}=++$i unless exists $r{$c}; }
    return %r;
}

sub valid_seq { my ($s)=@_; return defined($s) && $s ne '' && $s =~ /^[ACGTNacgtn]+$/; }
sub sanitize_id {
    my ($id,$chr,$pos,$type)=@_;
    $id = "$type.$chr.$pos" if !defined($id) || $id eq '' || $id eq '.';
    $id =~ s/[\s;=,]+/_/g;
    return $id;
}
sub safe { my ($x)=@_; return defined($x)&&$x ne ''?$x:'.'; }

sub usage {
    return <<'USAGE';
Usage:
  perl 13_prepare_vg_vcf.pl \
    --ref ISO1.fa \
    --vcfs All.largeSV.raw.vcf,All.smallIndel.vcf,All.CXSV.vg.vcf,All.sorted.SNP.vcf \
    --query_list query.list \
    --refname ISO1 \
    --refsize ISO1.sizes \
    --out All.SVGAP.vg.vcf

Outputs:
  --out FILE              vg-ready, coordinate-consistent VCF
  --report FILE           per-record normalization/skipping report [OUT.report.tsv]

Important behavior:
  * Symbolic DEL is rebuilt as an explicit left-anchored REF/ALT allele.
  * Explicit alleles are trimmed, reference-checked and assigned END/SVLEN from
    their actual REF/ALT lengths.
  * Symbolic INV is retained because vg construct -S supports <INV>.
  * Symbolic INS/CXSV without a complete sequence is skipped and reported.
  * Very long explicit alleles above --max_explicit_allele are skipped [5 Mb].
USAGE
}

################################################################################
# Lightweight FASTA random access, compatible with samtools .fai.
################################################################################
package FastFAI;
use strict;
use warnings;
use Fcntl qw(SEEK_SET);

sub new {
    my ($class,$fasta)=@_;
    die "ERROR: FASTA must be uncompressed: $fasta\n" if $fasta =~ /\.gz$/i;
    my $self=bless {fasta=>$fasta,idx=>{},order=>[]},$class;
    if (-s "$fasta.fai") { $self->_load_fai("$fasta.fai"); }
    else { $self->_build_index(); }
    open(my $fh,'<',$fasta) or die "ERROR: cannot open FASTA $fasta: $!\n";
    binmode($fh); $self->{fh}=$fh; return $self;
}
sub _load_fai {
    my ($self,$fai)=@_;
    open(my $IN,'<',$fai) or die "ERROR: cannot read $fai: $!\n";
    while (<$IN>) {
        chomp; next if /^\s*$/; my @f=split(/\t/); next unless @f>=5;
        $self->{idx}{$f[0]}={len=>$f[1]+0,offset=>$f[2]+0,line_bases=>$f[3]+0,line_width=>$f[4]+0};
        push @{$self->{order}},$f[0];
    }
    close $IN;
}
sub _build_index {
    my ($self)=@_; my $fa=$self->{fasta};
    open(my $IN,'<',$fa) or die "ERROR: cannot read $fa: $!\n"; binmode($IN);
    my ($name,$len,$offset,$lb,$lw);
    while (1) {
        my $p=tell($IN); my $line=<$IN>; last unless defined $line;
        if ($line =~ /^>(\S+)/) {
            if (defined $name) { $self->{idx}{$name}={len=>$len,offset=>$offset,line_bases=>$lb,line_width=>$lw}; push @{$self->{order}},$name; }
            $name=$1; ($len,$offset,$lb,$lw)=(0,undef,undef,undef);
        } else {
            next unless defined $name; my $raw=length($line); $line =~ s/[\r\n]+$//; my $bases=length($line); next if $bases==0;
            $offset=$p unless defined $offset; $lb=$bases unless defined $lb; $lw=$raw unless defined $lw; $len += $bases;
        }
    }
    if (defined $name) { $self->{idx}{$name}={len=>$len,offset=>$offset,line_bases=>$lb,line_width=>$lw}; push @{$self->{order}},$name; }
    close $IN;
}
sub contigs { return @{$_[0]{order}}; }
sub has_contig { return exists $_[0]{idx}{$_[1]}; }
sub length { return $_[0]{idx}{$_[1]}{len}; }
sub get_seq {
    my ($self,$chr,$start,$end)=@_; return undef unless $self->has_contig($chr);
    my $x=$self->{idx}{$chr}; $start=0 if $start<0; $end=$x->{len} if $end>$x->{len}; return '' if $end<=$start;
    my $need=$end-$start; my $line=$x->{line_bases}; my $width=$x->{line_width};
    my $byte=$x->{offset}+int($start/$line)*$width+($start%$line);
    seek($self->{fh},$byte,SEEK_SET) or die "ERROR: FASTA seek failed for $chr:$start-$end\n";
    my $seq='';
    while (CORE::length($seq)<$need) { my $buf=''; my $n=read($self->{fh},$buf,$need-CORE::length($seq)+int(($need-CORE::length($seq))/$line+2)*($width-$line)); last unless $n; $buf =~ s/[\r\n]//g; $seq.=$buf; }
    return uc(substr($seq,0,$need));
}
1;
