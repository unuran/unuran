############################################################
# $Id$
############################################################

# valid TAGs for section TAG =DISTR
%distr_TAGs =
    (
     "=ROUTINES"    => { "required" => "yes",
			 "scan" => \&scan_ROUTINES },
     
     "=DESCRIPTION" => { "required" => "no",
			 "scan" => \&scan_DESCRIPTION },
     
     "=END"         => { "required" => "no",
			 "scan" => \&scan_END },
     );

############################################################

sub scan_DISTR {
    my $file = $_[0];
    my $handle = $_[1];
    my $line = $_[2];

    # get short and long distribution name
    chomp $line;
    $line =~ s/^\s*=DISTR\s*//;
    (my $distr, my $order, my $distr_long) = split /\s+/, $line, 3;

    # print distribution name 
    print STDERR "$distr -- $distr_long\n\t" if $VERBOSE;

    # get rank in order (second entry in list: number in square brackets)
    die "no order for distribution" unless ($order =~ /\[(\-?\+?\d+)\]/);
    $order = $1;

    # check for uniqueness of distribution name
    if ($in->{"=DISTR"}->{$distr}->{"=NAME"}) {
	die "distribution name $distr not unique";
    }

    # store method name
    $in->{"=DISTR"}->{$distr}->{"=NAME"} = $distr_long;

    # store order 
    $in->{"=DISTR"}->{$distr}->{"=ORDER"} = $order;

    # store file name
    $in->{"=DISTR"}->{$distr}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($distr_TAGs{$this_TAG}) {
		# store subsection text
		$in->{"=DISTR"}->{$distr}->{$this_TAG} .= $lines;
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
	unless ($in->{"=DISTR"}->{$distr}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n" if $VERBOSE;
	}
    }

    # scan and format all sections
    foreach $tag (keys %distr_TAGs) {
	&{$distr_TAGs{$tag}{"scan"}}( \($in->{"=DISTR"}->{$distr}->{$tag}) );
    }

} # end of scan_DISTR()

############################################################

sub format_DISTR {

    # write texi output
    $texi->{"=DISTR"} .= "\@node Distributions\n";
    $texi->{"=DISTR"} .= "\@chapter Distributions\n\n";
    
    # get list of distributions
    my @list_distr = sort distr_by_order_key keys %{$in->{"=DISTR"}};

    # make menu for all distributions
    $texi->{"=DISTR"} .= "\@menu\n";
    foreach my $distr (@list_distr) {
	$texi->{"=DISTR"} .= "* $distr\:: ".$in->{"=DISTR"}->{$distr}->{"=NAME"}."\n";
    }
    $texi->{"=DISTR"} .= "\@end menu\n\n";

    # print subsections for all distribution types
    foreach my $distr (@list_distr) {

	# write texi subsection header for distribution type

	# header file name
	$texi->{"=DISTR"} .= "\@c -------------------------------------\n";
	$texi->{"=DISTR"} .= "\@c ".$in->{"=DISTR"}->{$distr}->{"=FILE"}."\n";
	$texi->{"=DISTR"} .= "\@c\n\n";

	# node and section
	$texi->{"=DISTR"} .= "\@node $distr\n";
	$texi->{"=DISTR"} .= "\@section ".$in->{"=DISTR"}->{$distr}->{"=NAME"}." ($distr)\n\n";

	# description for distribution
	$texi->{"=DISTR"} .= $in->{"=DISTR"}->{$distr}->{"=DESCRIPTION"}."\n\n";

	# function reference
	$texi->{"=DISTR"} .= "\n\@subheading Function reference\n\n";
	$texi->{"=DISTR"} .= $in->{"=DISTR"}->{$distr}->{"=ROUTINES"}."\n\n";

	# end of header file
	$texi->{"=DISTR"} .= "\@c\n";
	$texi->{"=DISTR"} .= "\@c end of ".$in->{"=DISTR"}->{$distr}->{"=FILE"}."\n";
	$texi->{"=DISTR"} .= "\@c -------------------------------------\n";

    }
    return;

} # end of format_DISTR() 

############################################################

sub distr_by_order_key {
    $in->{"=DISTR"}->{$a}->{"=ORDER"} <=> $in->{"=DISTR"}->{$b}->{"=ORDER"};
}

############################################################

# end of file
1;



