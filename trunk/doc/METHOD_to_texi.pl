
############################################################

# valid TAGs for section TAG =METHOD
%method_TAGs =
    (
     
     "=TYPE"        => { "required" => "yes",
			 "scan" => \&scan_METHOD_TYPE },

     "=REQUIRED"    => { "required" => "no", 
			 "scan" => \&scan_METHOD_REQUIRED },
     
     "=OPTIONAL"    => { "required" => "no",
			 "scan" => \&scan_METHOD_OPTIONAL },
     
     "=SEEALSO"     => { "required" => "no",
			 "scan" => \&scan_METHOD_SEEALSO },
     
     "=ABSTRACT"    => { "required" => "no",
			 "scan" => \&scan_METHOD_ABSTRACT },
     
     "=DESCRIPTION" => { "required" => "yes",
			 "scan" => \&scan_METHOD_DESCRIPTION },
     
     "=ROUTINES"    => { "required" => "yes",
			 "scan" => \&scan_ROUTINES },
     
     "=END"         => { "required" => "no",
			 "scan" => \&scan_END },
     );

# store data
$in_METHODs;

# formated text
$texi_METHODs;


############################################################

sub scan_METHOD {
    my $file = $_[0];
    my $handle = $_[1];
    my $line = $_[2];

    # get short and long method name
    chomp $line;
    $line =~ s/^\s*=METHOD\s*//;
    (my $method, my $method_long) = split /\s+/, $line, 2;

    # print method name 
    print STDERR "$method -- $method_long\n\t" if $VERBOSE;

    # check for uniqueness of method name
    if ($in_METHODs->{$method}->{"=NAME"}) {
	die "method name not unique";
    }

    # store method name
    $in_METHODs->{$method}->{"=NAME"} = $method_long;

    # store file name
    $in_METHODs->{$method}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($method_TAGs{$this_TAG}) {
		# store subsection text
		$in_METHODs->{$method}->{$this_TAG} .= $lines;
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
    foreach $tag (keys %method_TAGs) {
	next unless ($method_TAGs{$tag}{"required"} eq "yes");
	unless ($in_METHODs->{$method}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n";
	}
    }

    # scan and format all sections
    foreach $tag (keys %method_TAGs) {
	&{$method_TAGs{$tag}{"scan"}}( \($in_METHODs->{$method}->{$tag}) );
    }
    
} # end of scan_METHOD()

############################################################

sub scan_METHOD_TYPE {
    my $entry = $_[0];    # pointer to TYPE entry

    # we are interested in the very first word only ...
    unless ( $$entry =~ /^\s*([A-Za-z]+)/ ) {
	die "TYPE: missing type";
    }
    $$entry = $1;

    return;
}

############################################################

sub scan_METHOD_REQUIRED {
    return;
}

############################################################

sub scan_METHOD_OPTIONAL {
    return;
}

############################################################

sub scan_METHOD_SEEALSO {
    return;
}

############################################################

sub scan_METHOD_ABSTRACT {
    return;
}

############################################################

sub scan_METHOD_DESCRIPTION {
    return;
}

############################################################

sub format_METHOD {

    # check distribution types
    foreach my $method (keys %$in_METHODs) {
	my $type_ok = 0;
	print STDERR "\t$method -> ".$in_METHODs->{$method}->{"=TYPE"}." ... " if $VERBOSE;
	foreach my $type (keys %$in_DISTRs) {
	    if ($in_METHODs->{$method}->{"=TYPE"} eq $type) {
		$type_ok = 1;
	    }
	}
	print STDERR (($type_ok) ? "ok\n" : "invalid\n")  if $VERBOSE;
	die "invalid distribution type for method $method" unless $type_ok;
    }

    # write texi output
    $texi_METHODs .= "\@node Methods\n";
    $texi_METHODs .= "\@chapter Methods\n\n";
    
    # sort by distribution type
    foreach my $type (sort distr_by_order_key keys %$in_DISTRs) {
	# distributions with non positive order keys are ignored
	next if $in_DISTRs->{$type}->{"=ORDER"} <= 0;

	# write texi subsection header for distribution type
	$texi_METHODs .= "\@node Methods for $type\n";
	$texi_METHODs .= "\@section Methods for ".$in_DISTRs->{$type}->{"=NAME"}." ($type)\n\n";

	# find all methods for distribution type
	foreach my $method (sort keys %$in_METHODs) {
	    next unless $in_METHODs->{$method}->{"=TYPE"} eq $type;
	    
	    # header file name
	    $texi_METHODs .= "\@c -------------------------------------\n";
	    $texi_METHODs .= "\@c ".$in_METHODs->{$method}->{"=FILE"}."\n";
	    $texi_METHODs .= "\@c\n\n";
	    
	    # node and subsection
	    $texi_METHODs .= "\@node $method\n";
	    $texi_METHODs .= "\@subsection $method -- ".$in_METHODs->{$method}->{"=NAME"}."\n\n";
	    
	    # description for method
	    $texi_METHODs .= $in_METHODs->{$method}->{"=DESCRIPTION"}."\n\n";

	    # function reference 
	    $texi_METHODs .= "\n\@subsubsection Function reference\n\n";
	    $texi_METHODs .= $in_METHODs->{$method}->{"=ROUTINES"}."\n\n";

	    # end of header file
	    $texi_METHODs .= "\@c\n";
	    $texi_METHODs .= "\@c end of ".$in_METHODs->{$method}->{"=FILE"}."\n";
	    $texi_METHODs .= "\@c -------------------------------------\n";
	    
	}

    }

    return;
} # end of format_METHOD() 

############################################################

sub distr_by_order_key {
    $in_DISTRs->{$a}->{"=ORDER"} <=> $in_DISTRs->{$b}->{"=ORDER"};
}

############################################################

# end of file
1;

