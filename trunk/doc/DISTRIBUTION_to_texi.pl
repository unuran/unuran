
############################################################

# valid TAGs for section TAG =DISTRIBUTION
%distr_TAGs =
    (
     "=ROUTINES"    => { "required" => "yes",
			 "scan" => \&scan_ROUTINES },
     
     "=END"         => { "required" => "no",
			 "scan" => \&scan_END },
     );

# store data
$in_DISTRs;

# formated text
$texi_DISTRs;

############################################################

sub scan_DISTR {
    my $file = $_[0];
    my $handle = $_[1];
    my $line = $_[2];

    # get short and long distribution name
    chomp $line;
    $line =~ s/^\s*=DISTRIBUTION\s*//;
    (my $distr, my $order, my $distr_long) = split /\s+/, $line, 3;

    # print distribution name 
    print STDERR "$distr -- $distr_long\n\t" if $VERBOSE;

    # get rank in order (second entry in list: number in square brackets)
    die "no order for distribution" unless ($order =~ /\[(\-?\+?\d+)\]/);
    $order = $1;

    # check for uniqueness of distribution name
    if ($in_DISTRs->{$distr}->{"=NAME"}) {
	die "distribution name $distr not unique";
    }

    # store method name
    $in_DISTRs->{$distr}->{"=NAME"} = $distr_long;

    # store order 
    $in_DISTRs->{$distr}->{"=ORDER"} = $order;

    # store file name
    $in_DISTRs->{$distr}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($distr_TAGs{$this_TAG}) {
		# store subsection text
		$in_DISTRs->{$distr}->{$this_TAG} .= $lines;
	    }
	    else {
		print STDERR "  invalid!!\n\n\t" if $VERBOSE;
	    }
	}
	# anything between =END and next tag is ignored
	$this_TAG = $next_TAG;
    }

    print STDERR "\n" if $VERBOSE;
    
    # analyse entries
    foreach $tag (keys %distr_TAGs) {
	next unless ($distr_TAGs{$tag}{"required"} eq "yes");
	unless ($in_DISTRs->{$distr}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n";
	}
    }

    # scan and format all sections
    foreach $tag (keys %distr_TAGs) {
	&{$distr_TAGs{$tag}{"scan"}}( \($in_DISTRs->{$distr}->{$tag}) );
    }

} # end of scan_DISTR()

############################################################

sub format_DISTR {

    # write texi output
    $texi_DISTRs .= "\@node Distributions\n";
    $texi_DISTRs .= "\@chapter Distributions\n\n";
    
    # sort by distribution type
    foreach my $distr (sort distr_by_order_key keys %$in_DISTRs) {

	# write texi subsection header for distribution type

	# header file name
	$texi_DISTRs .= "\@c -------------------------------------\n";
	$texi_DISTRs .= "\@c ".$in_DISTRs->{$distr}->{"=FILE"}."\n";
	$texi_DISTRs .= "\@c\n\n";

	# node and section
	$texi_DISTRs .= "\@node $type\n";
	$texi_DISTRs .= "\@section ".$in_DISTRs->{$distr}->{"=NAME"}." ($type)\n";

	# function reference
	$texi_DISTRs .= "\n\@subsubsection Function reference\n\n";
	$texi_DISTRs .= $in_DISTRs->{$distr}->{"=ROUTINES"}."\n\n";

	# end of header file
	$texi_DISTRs .= "\@c\n";
	$texi_DISTRs .= "\@c end of ".$in_DISTRs->{$distr}->{"=FILE"}."\n";
	$texi_DISTRs .= "\@c -------------------------------------\n";


    }
    return;

} # end of format_DISTR() 

############################################################

# end of file
1;



