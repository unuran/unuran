############################################################
# $Id$
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
			 "scan" => \&scan_DESCRIPTION },
     
     "=ROUTINES"    => { "required" => "yes",
			 "scan" => \&scan_ROUTINES },
     
     "=END"         => { "required" => "no",
			 "scan" => \&scan_END },
     );

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
    if ($in->{"=METHOD"}->{$method}->{"=NAME"}) {
	die "method name not unique";
    }

    # store method name
    $in->{"=METHOD"}->{$method}->{"=NAME"} = $method_long;

    # store file name
    $in->{"=METHOD"}->{$method}->{"=FILE"} = $file;

    # scan all subsections 
    my $this_TAG = "=END";  # add =END tag to section tag
    my $lines;

    while (my $next_TAG = scan_subsection($handle,\$lines)) {
	if ($this_TAG ne "=END") {
	    print STDERR "  $this_TAG" if $VERBOSE;
	    if ($method_TAGs{$this_TAG}) {
		# store subsection text
		$in->{"=METHOD"}->{$method}->{$this_TAG} .= $lines;
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
	unless ($in->{"=METHOD"}->{$method}->{$tag}) {
	    print STDERR "\t$tag is missing!!\n\n" if $VERBOSE;
	}
    }

    # scan and format all sections
    foreach $tag (keys %method_TAGs) {
	&{$method_TAGs{$tag}{"scan"}}( \($in->{"=METHOD"}->{$method}->{$tag}) );
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

sub format_METHOD {

    # check distribution types
    foreach my $method (keys %{$in->{"=METHOD"}}) {
	my $type_ok = 0;
	print STDERR "\t$method -> ".$in->{"=METHOD"}->{$method}->{"=TYPE"}." ... " if $VERBOSE > 1;
	foreach my $type (keys %{$in->{"=DISTR"}}) {
	    if ($in->{"=METHOD"}->{$method}->{"=TYPE"} eq $type) {
		$type_ok = 1;
	    }
	}
	print STDERR (($type_ok) ? "ok\n" : "invalid\n")  if $VERBOSE > 1;
	die "invalid distribution type for method $method" unless $type_ok;
    }

    # get list of all methods
    my @list_method = sort keys %{$in->{"=METHOD"}};

    # get list of distributions with existing methods
    my $list_distr_with_method;
    foreach my $distr (sort distr_by_order_key keys %{$in->{"=DISTR"}}) {
	foreach my $method (@list_method) {
	    next unless $in->{"=METHOD"}->{$method}->{"=TYPE"} eq $distr;
	    $list_distr_with_method->{$distr} = 1;
	}
    }
    my @list_distr = sort distr_by_order_key keys %$list_distr_with_method;

    # write texi output
    $texi->{"=METHOD"} .= "\@node Methods\n";
    $texi->{"=METHOD"} .= "\@chapter Methods\n\n";

    # make menu for all distribution types
    $texi->{"=METHOD"} .= "\@menu\n";
    foreach my $type (@list_distr) {
	$texi->{"=METHOD"} .= "* Methods for $type\:: ".$in->{"=DISTR"}->{$type}->{"=NAME"}."\n";
    }
    $texi->{"=METHOD"} .= "\@end menu\n\n";
    
    # print subsections for all distribution types
    foreach my $type (@list_distr) {
	# write texi subsection header for distribution type
	$texi->{"=METHOD"} .= "\@node Methods for $type\n";
	$texi->{"=METHOD"} .= "\@section Methods for ".$in->{"=DISTR"}->{$type}->{"=NAME"}." ($type)\n\n";

	# make menu for all methods for distribution type
	$texi->{"=METHOD"} .= "\@menu\n";
	foreach my $method (@list_method) {
	    next unless $in->{"=METHOD"}->{$method}->{"=TYPE"} eq $type;
	    $texi->{"=METHOD"} .= "* $method\:: ".$in->{"=METHOD"}->{$method}->{"=NAME"}."\n";
	}
	$texi->{"=METHOD"} .= "\@end menu\n\n";

	# print subsubsections for all methods for distribution types
	foreach my $method (@list_method) {
	    next unless $in->{"=METHOD"}->{$method}->{"=TYPE"} eq $type;
	    
	    # header file name
	    $texi->{"=METHOD"} .= "\@c -------------------------------------\n";
	    $texi->{"=METHOD"} .= "\@c ".$in->{"=METHOD"}->{$method}->{"=FILE"}."\n";
	    $texi->{"=METHOD"} .= "\@c\n\n";
	    
	    # node and subsection
	    $texi->{"=METHOD"} .= "\@node $method\n";
	    $texi->{"=METHOD"} .= "\@subsection $method -- ".$in->{"=METHOD"}->{$method}->{"=NAME"}."\n\n";
	    
	    # description for method
	    $texi->{"=METHOD"} .= $in->{"=METHOD"}->{$method}->{"=DESCRIPTION"}."\n\n";

	    # function reference 
	    $texi->{"=METHOD"} .= "\n\@subsubheading Function reference\n\n";
	    $texi->{"=METHOD"} .= $in->{"=METHOD"}->{$method}->{"=ROUTINES"}."\n\n";

	    # end of header file
	    $texi->{"=METHOD"} .= "\@c\n";
	    $texi->{"=METHOD"} .= "\@c end of ".$in->{"=METHOD"}->{$method}->{"=FILE"}."\n";
	    $texi->{"=METHOD"} .= "\@c -------------------------------------\n";
	    
	}

    }

    return;
} # end of format_METHOD() 

############################################################

# end of file
1;

