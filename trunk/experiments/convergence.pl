#!/usr/bin/perl
#############################################################################
#
# Run experiments with convergence of hat using TDR.
#
#############################################################################

use strict;
$|=1;

#############################################################################
# config

# limits for number of construction points
my $Nmin = 3;   
my $Nmax = 50;   

# sample size for ARS 
my $samplesize = 100000;

# distributions
my @distributions = 
    ( 
#      "beta(5,15)",
      "normal(0,1)",
      "cauchy()", 
#      "normal(0,2)",
#      "normal(0,4)",
#      "normal(0,8)",
#      "normal(0,16)",
#      "normal(0,32)",
#      "normal(0,64)",
#      "normal(0,128)",
#      "normal(0,10000)",
#      "normal(0,0.1)",
#      "normal(0,0.01)",
#      "normal(0,0.001)",
#      "normal(0,0.0001)",
#      "cauchy()", 
      );

#    ( "normal(); domain=(-1e10,1e10)",
#      "cauchy(); domain=(-1e10,1e10)" );
		  
#############################################################################
# global constants

my $rungen = "./rungen -L \"stdout\"";

#############################################################################
# global variables

my $plot_normal;    # plot (normal scale)
my $plot_log;       # plot (logarithmic scale)


#############################################################################
# run

my $method = "method=tdr; variant_gw; usemode=off; usecenter=off; max_sqhratio=1; ";
$method .= "max_intervals=".($Nmax+1);

print "\\documentclass{amsart}\n";
print "\\usepackage{pstricks,pst-plot}\n";
print "\\begin{document}\n";

foreach my $d (@distributions) {

    # -----------------------------------------------------------------------
    # ARS
#    ARS($d);

    # -----------------------------------------------------------------------
    # equiangular rule
#    equiangular($d);

    # -----------------------------------------------------------------------
    # DARS
    DARS1($d);
#    DARS2($d);
#    DARS3($d);

    # -----------------------------------------------------------------------
    # make plot
    make_plot_normal();
    print "\\clearpage\n";
    make_plot_log();
    print "\\clearpage\n";
}

print "\\end{document}\n";
    
#############################################################################
# end

exit (0);

#############################################################################

# ---------------------------------------------------------------------------
# equiangular rule

sub equiangular {
    my $distr = $_[0];

    my $graph_normal = "\t%% equiangular rule\n";
    $graph_normal .= "\t\\psline";
    my $graph_log = $graph_normal;

    for my $i ($Nmin .. $Nmax-1) {
	my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
	my $gen = "$distr & $method; usedars=off; cpoints=$i & $urng";

	my @data = get_data($gen);
	my $n = $data[0];
	my $rho = $data[$n]-1;
	$graph_normal .= sprintf("(%d,%.4f)", $n,$rho);
	$graph_log    .= sprintf("(%d,%.4f)", $n,log10($rho));
    }

    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";
} # end of equiangular()

# ---------------------------------------------------------------------------
# DARS

sub DARS1 {
    my $distr = $_[0];

    my $graph_normal = "\t%% DARS (split at expected point)\n";
    $graph_normal .= "\t\\psline";
    my $graph_log = $graph_normal;

    my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
    my $gen = "$distr & $method; usedars=1; cpoints=$Nmin & $urng";
    my @data = get_data($gen);
    my $nstart = $data[0];
    for my $i ($nstart .. $Nmax) {
	my $rho = $data[$i]-1;
	$graph_normal .= sprintf("(%d,%.4f)", $i,$rho);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($rho));
    }

    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";
} # end of DARS1()

sub DARS2 {
    my $distr = $_[0];

    my $graph_normal = "\t%% DARS (split at arcmean)\n";
    $graph_normal .= "\t\\psline";
    my $graph_log = $graph_normal;

    my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
    my $gen = "$distr & $method; usedars=2; cpoints=$Nmin & $urng";

    my @data = get_data($gen);
    my $nstart = $data[0];
    for my $i ($nstart .. $Nmax) {
	my $rho = $data[$i]-1;
	$graph_normal .= sprintf("(%d,%.4f)", $i,$rho);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($rho));
    }

    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";
} # end of DARS2()

sub DARS3 {
    my $distr = $_[0];

    my $graph_normal = "\t%% DARS (split at mean)\n";
    $graph_normal .= "\t\\psline";
    my $graph_log = $graph_normal;

    my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
    my $gen = "$distr & $method; usedars=3; cpoints=$Nmin & $urng";
    my @data = get_data($gen);
    my $nstart = $data[0];
    for my $i ($nstart .. $Nmax) {
	my $rho = $data[$i]-1;
	$graph_normal .= sprintf("(%d,%.4f)", $i,$rho);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($rho));
    }

    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";
} # end of DARS3()

# ---------------------------------------------------------------------------
# ARS

sub ARS {
    my $distr = $_[0];

    my @data_ARS;
    # run simulations
    for my $s (1 .. $samplesize) {
	unless ($s % 100) {
	    print STDERR "[".(int $s/100)."]";
	}
	my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
	my $gen = "$distr & $method; usedars=off; cpoints=$Nmin & $urng";
	my @data = get_data($gen);
	push @data_ARS, [ @data ];
    }

    # number of intervals at start
    my $nstart = $data_ARS[0][0];

    # analize data
    my $ARS_count;
    my @ARS_min;
    my @ARS_max;
    my @ARS_median;
    my @ARS_quart_up;
    my @ARS_quart_low;
    my @ARS_perz_up;
    my @ARS_perz_low;
    for my $i ($nstart .. $Nmax) {
	# get row ...
	my @data;
	for my $s ( 0 .. $samplesize-1 ) {
	    push @data, $data_ARS[$s][$i]-1;
	}
	# sort data ...
	my @sorted = sort @data;
	# now make statistics ... 
	my $count = $#sorted + 1;
	$ARS_count = $count;
	$ARS_min[$i] = $sorted[0];
	$ARS_max[$i] = $sorted[$#sorted];
	$ARS_median[$i] = $sorted[($count/2)-1];
	$ARS_quart_up[$i] = $sorted[3*$count/4-1];
	$ARS_quart_low[$i] = $sorted[$count/4-1];
	$ARS_perz_up[$i] = $sorted[19*$count/20-1];
	$ARS_perz_low[$i] = $sorted[$count/20-1];
    }

    # print data
    my $graph_normal = "\t%% ARS - range\n";
    $graph_normal .= "\t\\pspolygon[fillcolor=lightgray,fillstyle=solid,linestyle=none]";
    my $graph_log = $graph_normal;
    for my $i ($nstart .. $Nmax) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_max[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_max[$i]));
    }
    for (my $i=$Nmax; $i>=$nstart; --$i) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_min[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_min[$i]));
    }
    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";

    $graph_normal = "\t%% ARS - 5% percentile\n";
    $graph_normal .= "\t\\pspolygon[fillcolor=gray,fillstyle=solid,linestyle=none]";
    $graph_log = $graph_normal;
    for my $i ($nstart .. $Nmax) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_perz_up[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_perz_up[$i]));
    }
    for (my $i=$Nmax; $i>=$nstart; --$i) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_perz_low[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_perz_low[$i]));
    }
    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";

    $graph_normal = "\t%% ARS - quartiles\n";
    $graph_normal .= "\t\\pspolygon[fillcolor=darkgray,fillstyle=solid,linestyle=none]";
    $graph_log = $graph_normal;
    for my $i ($nstart .. $Nmax) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_quart_up[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_quart_up[$i]));
    }
    for (my $i=$Nmax; $i>=$nstart; --$i) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_quart_low[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_quart_low[$i]));
    }
    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";

    $graph_normal = "\t%% ARS - median\n";
    $graph_normal .= "\t\\psline";
    $graph_log = $graph_normal;
    for my $i ($nstart .. $Nmax) {
	$graph_normal .= sprintf("(%d,%.4f)", $i,$ARS_median[$i]);
	$graph_log    .= sprintf("(%d,%.4f)", $i,log10($ARS_median[$i]));
    }
    $plot_normal .= "$graph_normal\n";
    $plot_log    .= "$graph_log\n";


} # end of ARS()

#############################################################################
# make plot for normal scale

sub make_plot_normal {

    my $xmin = 0;
    my $xmax = $Nmax+1;
    my $ymin = 0;
    my $ymax = 4;
    my $xunit = 150/$xmax."mm";
    my $yunit = "10mm";

    print "{\n";
    print "\\newgray{darkgray}{0.70}\n";
    print "\\newgray{gray}{0.80}\n";
    print "\\newgray{lightgray}{0.90}\n";
    print "\\psset{xunit=$xunit,yunit=$yunit}\n";
    print "\\begin{pspicture}($xmin,$ymin)($xmax,$ymax)\n";

    print $plot_normal;

    print "\\psaxes[labels=all,ticks=all,Dx=5,Dy=1]{->}(0,0)($xmin,$ymin)($xmax,$ymax)\n";

    print "\\end{pspicture}\n"; 
    print "}\n";

    undef $plot_normal;
} # end of make_plot_normal()

#############################################################################
# make plot for logarithmic scale

sub make_plot_log {

    my $xmin = 0;
    my $xmax = $Nmax+1;
    my $ymin = -3;
    my $ymax = 1;
    my $xunit = 150/$xmax."mm";
    my $yunit = "10mm";

    print "{\n";
    print "\\newgray{darkgray}{0.70}\n";
    print "\\newgray{gray}{0.80}\n";
    print "\\newgray{lightgray}{0.90}\n";
    print "\\psset{xunit=$xunit,yunit=$yunit}\n";
    print "\\begin{pspicture}($xmin,$ymin)($xmax,$ymax)\n";

    print $plot_log;

    print "\\psaxes[labels=all,ticks=all,Dx=5,Dy=1]{->}(0,0)($xmin,$ymin)($xmax,$ymax)\n";

    print "\\end{pspicture}\n"; 
    print "}\n";

    undef $plot_log;
} # end of make_plot_log()

#############################################################################
# read UNURAN log file

sub get_data {
    my $gen = $_[0];

    # store data
    my @data;
    
    # tmp data
    my $Asqueeze;
    my $Ahat;
    my $rho;
    my $nint = 0;

    # read UNURAN log file
#    print STDERR "$rungen \"$gen\"\n";
    my @unuranlog = `$rungen "$gen"`;

    # scan log file
    foreach my $l (@unuranlog) {

	# stop when no construction points are added
	if ( ($l =~ /^TDR\.\d+\:\s*DARS completed \**\s*/) or
	     ($l =~ /^TDR\.\d+\:\s*GENERATOR destroyed \**\s*/) ) {
	    last;
	}

	# get starting number of intervals
	if ($nint == 0) {
	    if ($l =~ /^TDR\.\d+\:\s*Intervals\:\s*([\d]+)\s*.*$/) {
		$nint = $1;
		# store number in field 0 of data
		$data[0] = $nint;
	    }
	}

	# we are only interested in lines that print area below hat and squeeze
	next unless $l =~ /^TDR\.\d+\:\s*A\((squeeze|total)\)\s*=\s*([-\d\.Ee]+).*$/;

	# area below squeeze (first entry)
	if ($1 eq "squeeze") {
	    $Asqueeze = $2;
	    next;
	}
	
	# area below hat (last entry)
	if ($1 eq "total") {
	    $Ahat = $2;
	    $rho = ($Asqueeze>0) ? $Ahat/$Asqueeze : 999999;
	}

	# store data
	$data[$nint] = $rho;
	++$nint;
    }
    return @data;
} # end of get_data()

#############################################################################
# logarithm to base 10 (common logarithm)

sub log10 {
    my $n = shift;

    if ($n <= 0) { print STDERR "junk\n"; return 999999; }
    return log($n)/log(10);
}

