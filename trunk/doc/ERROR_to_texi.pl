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
    $in->{"=ERROR"}->{$error}->{"=NAME"} = $error_title;

    # store file name
    $in->{"=ERROR"}->{$error}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($error_TAGs{$this_TAG}) {
		# store subsection text
		$in->{"=ERROR"}->{$error}->{$this_TAG} .= $lines;
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
	unless ($in->{"=ERROR"}->{$error}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n" if $VERBOSE;
	}
    }

    # scan and format all sections
    foreach $tag (keys %error_TAGs) {
	&{$error_TAGs{$tag}{"scan"}}( \($in->{"=ERROR"}->{$error}->{$tag}) );
    }

} # end of scan_ERROR()

############################################################

sub format_ERROR {
    my $error = "ERROR";

    # write texi output
    $texi->{"=ERROR"} .= "\@node Error and Debuging\n";
    $texi->{"=ERROR"} .= "\@chapter Error and Debuging\n";

    # get list of sections
    my @list_error = sort keys %{$in->{"=ERROR"}};

    # make menu for all distributions
    $texi->{"=ERROR"} .= "\@menu\n";
    foreach my $error (@list_error) {
        $texi->{"=ERROR"} .= "* $error\:: ".$in->{"=ERROR"}->{$error}->{"=NAME"}."\n";
    }
    $texi->{"=ERROR"} .= "\@end menu\n\n";

    # print subsections for all distribution types
    foreach my $error (@list_error) {
    
	# write texi subsection header for erroribution type

	# header file name
	$texi->{"=ERROR"} .= "\@c -------------------------------------\n";
	$texi->{"=ERROR"} .= "\@c ".$in->{"=ERROR"}->{$error}->{"=FILE"}."\n";
	$texi->{"=ERROR"} .= "\@c\n\n";

        # node and section
        $texi->{"=ERROR"} .= "\@node $error\n";
        $texi->{"=ERROR"} .= "\@section ".$in->{"=ERROR"}->{$error}->{"=NAME"}."\n\n";

	# description for erroribution
	$texi->{"=ERROR"} .= $in->{"=ERROR"}->{$error}->{"=DESCRIPTION"}."\n\n";

	# function reference
	$texi->{"=ERROR"} .= "\n\@heading Function reference\n\n";
	$texi->{"=ERROR"} .= $in->{"=ERROR"}->{$error}->{"=ROUTINES"}."\n\n";

	# end of header file
	$texi->{"=ERROR"} .= "\@c\n";
	$texi->{"=ERROR"} .= "\@c end of ".$in->{"=ERROR"}->{$error}->{"=FILE"}."\n";
	$texi->{"=ERROR"} .= "\@c -------------------------------------\n";
    }

    return;

} # end of format_ERROR() 

############################################################

# end of file
1;
