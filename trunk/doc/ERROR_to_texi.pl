############################################################
# $Id$
############################################################

# valid TAGs for section TAG =DISTRIBUTION
%error_TAGs =
    (
     "=ROUTINES"    => { "required" => "yes",
			 "scan" => \&scan_ROUTINES },
     
     "=DESCRIPTION" => { "required" => "yes",
			 "scan" => \&scan_DESCRIPTION },
     
     "=END"         => { "required" => "no",
			 "scan" => \&scan_END },
     );

# store data
$in_ERRORs;

# formated text
$texi_ERRORs;

############################################################

sub scan_ERROR {
    my $file = $_[0];
    my $handle = $_[1];
    my $line = $_[2];

    # get title
    chomp $line;
    $line =~ s/^\s*=ERROR\s*//;
    (my $error, my $error_title) = split /\s+/, $line, 2;

    # print title
    print STDERR "$error_title\n\t" if $VERBOSE;

    # store method name
    $in_ERRORs->{$error}->{"=NAME"} = $error_title;

    # store file name
    $in_ERRORs->{$error}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($error_TAGs{$this_TAG}) {
		# store subsection text
		$in_ERRORs->{$error}->{$this_TAG} .= $lines;
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
    foreach $tag (keys %error_TAGs) {
	next unless ($error_TAGs{$tag}{"required"} eq "yes");
	unless ($in_ERRORs->{$error}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n" if $VERBOSE;
	}
    }

    # scan and format all sections
    foreach $tag (keys %error_TAGs) {
	&{$error_TAGs{$tag}{"scan"}}( \($in_ERRORs->{$error}->{$tag}) );
    }

} # end of scan_ERROR()

############################################################

sub format_ERROR {
    my $error = "ERROR";

    # write texi output
    $texi_ERRORs .= "\@node Error and Debuging\n";
    $texi_ERRORs .= "\@chapter Error and Debuging\n";

    # get list of sections
    my @list_error = sort keys %$in_ERRORs;

    # make menu for all distributions
    $texi_ERRORs .= "\@menu\n";
    foreach my $error (@list_error) {
        $texi_ERRORs .= "* $error\:: ".$in_ERRORs->{$error}->{"=NAME"}."\n";
    }
    $texi_ERRORs .= "\@end menu\n\n";

    # print subsections for all distribution types
    foreach my $error (@list_error) {
    
	# write texi subsection header for erroribution type

	# header file name
	$texi_ERRORs .= "\@c -------------------------------------\n";
	$texi_ERRORs .= "\@c ".$in_ERRORs->{$error}->{"=FILE"}."\n";
	$texi_ERRORs .= "\@c\n\n";

        # node and section
        $texi_ERRORs .= "\@node $error\n";
        $texi_ERRORs .= "\@section ".$in_ERRORs->{$error}->{"=NAME"}."\n\n";

	# description for erroribution
	$texi_ERRORs .= $in_ERRORs->{$error}->{"=DESCRIPTION"}."\n\n";

	# function reference
	$texi_ERRORs .= "\n\@heading Function reference\n\n";
	$texi_ERRORs .= $in_ERRORs->{$error}->{"=ROUTINES"}."\n\n";

	# end of header file
	$texi_ERRORs .= "\@c\n";
	$texi_ERRORs .= "\@c end of ".$in_ERRORs->{$error}->{"=FILE"}."\n";
	$texi_ERRORs .= "\@c -------------------------------------\n";
    }

    return;

} # end of format_ERROR() 

############################################################

# end of file
1;



