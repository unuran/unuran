############################################################
# $Id$
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
