#!/usr/bin/perl

use strict;
$|=1;

#############################################################################
# config

# limits for number of construction points
my $Nmin = 2;   
my $Nmax = 30;   

# sample size for ARS 
my $samplesize = 3;

# distributions
my @distributions = 
    ( "normal()" );
		  
#############################################################################
# global constants

my $rungen = "./rungen -L \"stdout\"";

#############################################################################
# run

my $method = "method=tdr; variant_gw; usemode=off; usecenter=off; max_intervals=$Nmax+1; max_sqhratio=1";

foreach my $d (@distributions) {

    # -----------------------------------------------------------------------
    # equiangular rule

    my $graph_equi = "\\psline";
    for my $i ($Nmin .. $Nmax) {
	my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
	my $gen = "$d & $method; usedars=off; cpoints=$i & $urng";

	my @data = get_data($gen);
	my $n = $data[0];
	my $rho = $data[$n]-1;
	$graph_equi .= sprintf("(%d,%.4f)", $n,$rho);
    }
    print "$graph_equi\n%\n";

    # -----------------------------------------------------------------------
    # DARS

    my $graph_DARS = "\\psline";
    my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
    my $gen = "$d & $method; usedars=on; cpoints=$Nmin & $urng";
    my @data = get_data($gen);
    my $nstart = $data[0];
    for my $i ($nstart .. $Nmax) {
	my $rho = $data[$i]-1;
	$graph_DARS .= sprintf("(%i,%.4f)", $i,$rho);
    }
    print "$graph_DARS\n%\n";

    # -----------------------------------------------------------------------
    # ARS

    my @data_ARS;
    # run simulations
    for my $s (1 .. $samplesize) {
	my $urng = "urng=mt19937(".(int(rand 12345678)+1).")";
	my $gen = "$d & $method; usedars=off; cpoints=$Nmin & $urng";
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
    my @ARS_quant_up;
    my @ARS_quant_low;
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
	$ARS_quant_up[$i] = $sorted[3*$count/4-1];
	$ARS_quant_low[$i] = $sorted[$count/4-1];
	$ARS_perz_up[$i] = $sorted[19*$count/20-1];
	$ARS_perz_low[$i] = $sorted[$count/20-1];
    }

    # print data
    my $graph_ARS_median = "\\psline";
    for my $i ($nstart .. $Nmax) {
	$graph_ARS_median .= sprintf("(%d,%.4f)", $i,$ARS_median[$i]);
    }
    print "$graph_ARS_median\n%\n";
}
    
#############################################################################
# end

exit (0);

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
	next unless $l =~ /^TDR\.\d+\:\s+A\((squeeze|total)\)\s*=\s*([\d\.]+).*$/;

	# area below squeeze (first entry)
	if ($1 eq "squeeze") {
	    $Asqueeze = $2;
	    next;
	}
	
	# area below hat (last entry)
	if ($1 eq "total") {
	    $Ahat = $2;
	    $rho = $Ahat / $Asqueeze;
	}

	# store data
	$data[$nint] = $rho;
	++$nint;
    }
    return @data;
} # end of get_data()
