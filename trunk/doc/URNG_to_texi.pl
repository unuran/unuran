############################################################
# $Id$
############################################################

# valid TAGs for section TAG =DISTRIBUTION
%urng_TAGs =
    (
     "=ROUTINES"    => { "required" => "yes",
			 "scan" => \&scan_ROUTINES },
     
     "=DESCRIPTION" => { "required" => "yes",
			 "scan" => \&scan_DESCRIPTION },
     
     "=END"         => { "required" => "no",
			 "scan" => \&scan_END },
     );

# store data
$in_URNGs;

# formated text
$texi_URNGs;

############################################################

sub scan_URNG {
    my $file = $_[0];
    my $handle = $_[1];
    my $line = $_[2];

    # get title
    chomp $line;
    $line =~ s/^\s*=URNG\s*//;
    (my $urng, my $urng_title) = split /\s+/, $line, 2;

    # print title
    print STDERR "$urng_title\n\t" if $VERBOSE;

    # store method name
    $in_URNGs->{$urng}->{"=NAME"} = $urng_title;

    # store file name
    $in_URNGs->{$urng}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($urng_TAGs{$this_TAG}) {
		# store subsection text
		$in_URNGs->{$urng}->{$this_TAG} .= $lines;
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
    foreach $tag (keys %urng_TAGs) {
	next unless ($urng_TAGs{$tag}{"required"} eq "yes");
	unless ($in_URNGs->{$urng}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n" if $VERBOSE;
	}
    }

    # scan and format all sections
    foreach $tag (keys %urng_TAGs) {
	&{$urng_TAGs{$tag}{"scan"}}( \($in_URNGs->{$urng}->{$tag}) );
    }

} # end of scan_URNG()

############################################################

sub format_URNG {
    my $urng = "URNG";

    # write texi output
    $texi_URNGs .= "\@node URNG\n";
    $texi_URNGs .= "\@chapter ".$in_URNGs->{$urng}->{"=NAME"}."\n\n";
    
    # write texi subsection header for urngibution type

    # header file name
    $texi_URNGs .= "\@c -------------------------------------\n";
    $texi_URNGs .= "\@c ".$in_URNGs->{$urng}->{"=FILE"}."\n";
    $texi_URNGs .= "\@c\n\n";

    # description for urngibution
    $texi_URNGs .= $in_URNGs->{$urng}->{"=DESCRIPTION"}."\n\n";

    # function reference
    $texi_URNGs .= "\n\@heading Function reference\n\n";
    $texi_URNGs .= $in_URNGs->{$urng}->{"=ROUTINES"}."\n\n";

    # end of header file
    $texi_URNGs .= "\@c\n";
    $texi_URNGs .= "\@c end of ".$in_URNGs->{$urng}->{"=FILE"}."\n";
    $texi_URNGs .= "\@c -------------------------------------\n";

    return;

} # end of format_URNG() 

############################################################

# end of file
1;



