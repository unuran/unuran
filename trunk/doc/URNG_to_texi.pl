############################################################
# $Id$
############################################################

############################################################

sub format_URNG {
    my $urng = "URNG";

    # write texi output
    $texi->{"=URNG"} .= "\@node URNG\n";
    $texi->{"=URNG"} .= "\@chapter ".$in->{"=URNG"}->{$urng}->{"=NAME"}."\n\n";
    
    # write texi subsection header for urngibution type

    # header file name
    $texi->{"=URNG"} .= "\@c -------------------------------------\n";
    $texi->{"=URNG"} .= "\@c ".$in->{"=URNG"}->{$urng}->{"=FILE"}."\n";
    $texi->{"=URNG"} .= "\@c\n\n";

    # description for urngibution
    $texi->{"=URNG"} .= $in->{"=URNG"}->{$urng}->{"=DESCRIPTION"}."\n\n";

    # function reference
    $texi->{"=URNG"} .= "\n\@heading Function reference\n\n";
    $texi->{"=URNG"} .= $in->{"=URNG"}->{$urng}->{"=ROUTINES"}."\n\n";

    # end of header file
    $texi->{"=URNG"} .= "\@c\n";
    $texi->{"=URNG"} .= "\@c end of ".$in->{"=URNG"}->{$urng}->{"=FILE"}."\n";
    $texi->{"=URNG"} .= "\@c -------------------------------------\n";

    return;

} # end of format_URNG() 

############################################################

# end of file
1;



