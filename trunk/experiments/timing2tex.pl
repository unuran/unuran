#!/usr/bin/perl
# ----------------------------------------------------------------------------

my $Debug = 0;

# ----------------------------------------------------------------------------

use strict;
use Getopt::Long;

# ----------------------------------------------------------------------------
# Global variables ...

my $distr;
my $method;
my $dxm;
my $mxd;

my $null;

my $tag = \$null;

my $Format = "%-5.2f";

#----------------------------------------------------------------------------
# Read data ...

GetOptions( 'format=s' => \$Format ) || die;

while(<STDIN>) {
    next unless $_=~/[^\s]/;
    if (/^</) {
	$tag = \$distr, next  if /^<\d+\s+distributions>/;
	$tag = \$method, next if /^<\d+\s+methods>/;
	$tag = \$mxd,next     if /^<Method x Distribution>/;
	$tag = \$dxm,next     if /^<Distribution x Method>/;
	$tag = \$null;
    }
    $_ =~ s/_/\\_/g;
    $$tag .= $_;
}

if ($Debug) {
    print STDERR "distributions:\n$distr\n";
    print STDERR "methods:\n$method\n";
    print STDERR "dxm:\n$dxm\n";
    print STDERR "mxd:\n$mxd\n";
    print STDERR "ignore:\n$null\n\n";
}

unless ($distr and $method and $mxd) { die "data missing\n"; }

#----------------------------------------------------------------------------
# Print table headng 
print
    "\\begin{table}\n",
    "  \\centering\n",
    "  \\caption{Timing tests}\n",
    "  \\label{tab:timing}\n";

#----------------------------------------------------------------------------
# Print tables of distributions

print
    "  \\mbox{}\\hfill\n",
    "  \\begin{tabular}[t]{l\@{~~~}l}\n",
    "    \\hlinex\n",
    "     Symbol \& Distribution \\\\\n",
    "    \\xhlinex\n";

my @lines = split /\n/, $distr;
my $n_distr = $#lines + 1;
foreach my $l (@lines) {
    chomp $l;
    $l =~ s/^\s*//;
    (my $symbol, my $name) = split /\s+\.\.\.\s+/, $l;
    $name =~ s/distr=\s*//;
    print "      {$symbol} \& $name \\\\\n";
}

print
    "    \\xhline\n",
    "  \\end{tabular}\n",
    "  \\hfill\\mbox{}\n";

if ($Debug) { print STDERR "number of distributions = $n_distr\n"; }

#----------------------------------------------------------------------------
# Print tables of Methods

print
    "  \\mbox{}\\hfill\n",
    "  \\begin{tabular}[t]{l\@{~~~}l}\n",
    "    \\hlinex\n",
    "      Symbol \& Algorithm \\\\\n",
    "    \\xhlinex\n";

my @lines = split /\n/, $method;
my $n_methods = $#lines + 1;
foreach my $l (@lines) {
    chomp $l;
    (my $symbol, my $name) = split /\s+\.\.\.\s+/, $l;
    $name =~ s/method=\s*//;
    print "      {$symbol} \& $name \\\\\n";
}

print
    "    \\xhline\n",
    "  \\end{tabular}\n",
    "  \\hfill\\mbox{}\n";

if ($Debug) { print STDERR "number of methods = $n_methods\n"; }

#----------------------------------------------------------------------------
# Print timings

print
    "  \\par\\bigskip\n";
printf
    "  \\begin{tabular}{l*{%d}{\@{~~}c}}\n", $n_distr;
print
    "    \\hlinex\n";

my @lines = split /\n/, $mxd;

my $header = shift @lines;
$header =~ s/\s+/ \& /g;
print
    "          $header \\\\\n",
    "    \\xhlinex\n";

foreach my $l (@lines) {
    chomp $l;
    $l =~ s/^\s*//;
    $l =~ s/\s*$//;
    my @results = split /\s+/, $l;
    my $m = shift @results;
    print "    {$m}";
    foreach my $t (@results) {
	my $s = ($t =~ /\d/) ? sprintf($Format,$t) : "--";
	printf " \& %-8s",$s;
    }
    print "\\\\\n";
}

print
    "    \\xhline\n",
    "  \\end{tabular}\n";

#----------------------------------------------------------------------------
# end of table

print
    "\\end{table}\n\n";

#----------------------------------------------------------------------------
# End ...

exit (0);

#----------------------------------------------------------------------------

