# ----------------------------------------------------------------
# Read config file for testing
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Get list of distributions 

sub get_test_distributions
{
    # config file
    my $test_conf_file = $_[0];

    # list of distributions from PDF files
    # (to verify existance of distributions in config file)
    my $DISTR = $_[1];

    # for storing list of distribution
    my $list_distr;
    my $n_distr = 0;
    
    # scan through config file entries
    foreach my $t (split /\n\s*\n/, read_test_conf($test_conf_file)) {

	# There might be empty entries --> skip
	next unless $t;
	
	# Get name and parameters of distribution
	$t =~ /DISTR:\s*(\w+)\s*\(([^\)]*)\)/ or next;

	# Store data
	my $distr = $1;
	my $params = $2;

	# Check for existing distribution
	die "Unknown distribution: $distr" unless $DISTR->{$distr};
	
	# Process parameters
	my @param_list = split /\,/, $params;
	my $n_params = $#param_list + 1;
	my $fpm;
	if ($n_params > 0) {
	    $fpm = "";
	    foreach my $p (@param_list) {
		if ($p =~ /(.+)\.\.(.+)/) {
		    $fpm .= ($1>0) ? exp(log($1)+rand()*(log($2)-log($1))) : $1+rand()*($2-$1);
		    $fpm .= " ";
		}
		else {
		    $fpm .= "$p ";
		}
	    }
	    $fpm =~ s/,\s*$/ \}/;
	}
	else {
	    $fpm = "";
	}

	# Make distribution object
	$list_distr->{"$distr\_$n_distr"} = $fpm;

	# increment counter for distributions
	++$n_distr;
    }

    # return result
    return $list_distr;

} # end of get_test_distributions()

# ----------------------------------------------------------------
# Read the test.conf file

# content of config file
my $test_conf_content;

sub read_test_conf
{
    my $test_conf_file = $_[0];

    if ($test_conf_content) {
	# file already read --> nothing to do
	return $test_conf_content;
    }

    # open file
    open CONF, $test_conf_file or die "cannot open file $test_conf_file\n";

    # read all lines
    while (<CONF>) {

	# ignore comment lines
	next if /^\#/;

	# concatenate lines that start with blanks
	if (/^\s+[^\s\n]/) {
	    $test_conf_content =~ s/\n$//;
	}

	# store line
	$test_conf_content .= $_;
    }

    # close file handle
    close CONF;
    
    # return result
    $test_conf_content;

} # end of read_test_conf()

# ----------------------------------------------------------------
# End of file
return 1;

# ----------------------------------------------------------------
