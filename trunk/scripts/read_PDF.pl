# ----------------------------------------------------------------
# Read data for PDF from source code
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

use strict;
use File::Find;

# ----------------------------------------------------------------
# Global constants

# Header file for UNURAN standard distributions
my $h_stddistr = "unuran_distributions.h";

# List of distribution types
my %distr_types =
    ( "CONT"  => { "file_prefix" => "c",
		   "PDF_prefix"  => "_unur_pdf_",
		   "PDF_type"    => "double" },

      "DISCR" => { "file_prefix" => "d",
		   "PDF_prefix"  => "_unur_pmf_",
		   "PDF_type"    => "double" } );

# ................................................................
# List of files
my %file_list;

# ................................................................
# List of distributions
my $DISTR;

# Description of data fields (with beta distribution as example)
#
#   $DISTR->{"beta"}                           ... entries for distribution "beta"
#
#   $DISTR->{"beta"}->{"=NAME"}                ... name of distribution
#   $DISTR->{"beta"}->{"=TYPE"}                ... type of distribution (CONT|DISCR)
#   $DISTR->{"beta"}->{"=FILE"}                ... file name + path for C file
#   $DISTR->{"beta"}->{"=ID"}                  ... id of distribution
#
#   $DISTR->{"beta"}->{"=DOC"}                 ... documentation for distribution
#   $DISTR->{"beta"}->{"=DOC"}->{"=PDF"}       ... formula for PDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=CONST"}     ... normalization constant for PDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=CDF"}       ... formula for CDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=DOMAIN"}    ... domain for PDF
#   $DISTR->{"beta"}->{"=DOC"}->{"=FPARAM"}    ... list of parameters with constraints
#      Remark: There exist other fields which are not relevant here
#              (see src/distributions/unuran_distributions.h).
#
#   $DISTR->{"beta"}->{"=PDF"}                 ... PDF of distribution
#   $DISTR->{"beta"}->{"=PDF"}->{"=NAME"}      ... name of PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=RTYPE"}     ... return type of PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=ARGS"}      ... list of arguments for PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=N_PARAMS"}  ... number of parameters for PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=PARAMS"}[2] ... parameter #2 for PDF (starting at 0)
#   $DISTR->{"beta"}->{"=PDF"}->{"=BODY"}      ... function body of PDF
#   $DISTR->{"beta"}->{"=PDF"}->{"=CONST"}     ... macro expansion for (LOG)NORMCONSTANT
#   $DISTR->{"beta"}->{"=PDF"}->{"=DISTR"}     ... macro expansion for DISTR
#

# ----------------------------------------------------------------
# Read data for PDF from source code

sub read_PDFdata {

# Start the search for the files in these directories
    my @Startdirs = @_;

# ................................................................
# List of files
    find (\&find_files, @Startdirs);

# ................................................................
# Scan header file for UNURAN standard distributions
# and get a list of valid distributions
    scan_stddistr( $file_list{$h_stddistr} );

# ................................................................
# End. Return result
    return $DISTR;

} # end of read_PDFdata()

# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Finds recursively all C- or H-files containing the relevant
# informations. 

sub find_files 
{
    if ($File::Find::name =~ /(.*\/)?(.+\.[ch])$/ ) {
	my $file_path =  $File::Find::name;
	my $file_name =  $2;
	$file_list{$file_name} = $file_path;
    }
} # end of find_files()


# ----------------------------------------------------------------
# Scan header file for UNURAN standard distributions

sub scan_stddistr
{
    # header file for standard distributions
    my $file = $_[0];

    # Read header file
    open HFILE, $file or die "cannot open file \"$file\".\n";
    my $file_content;
    while (<HFILE>) {
	chomp;
	$_ =~ s/^\s*(.*)\s*$/$1/;
	$file_content .= "$_\n";
    }
    close HFILE;

    # Split into sections
    my @sections = split /=EON/, $file_content;

    # Scan sections
    foreach my $s (@sections) {
	# Remove trailing part of sections
	(my $dummy, $s) = split /=DISTR/, $s, 2;
	next unless $s;

	# Get name of distribution
	$s =~ s/(\w+)\s+(.+)\n//;
	my $distr = $1;
	$DISTR->{$distr}->{"=NAME"} = $2;

	# Scan distribution section
	scan_distr($distr,$s);

	# Get type of distribution and path of distribution file
	get_distr_file($distr);

	# Read distribution file
	read_distr_file($distr);
    }
    
} # end of scan_h_stddistr()

# ----------------------------------------------------------------
# Scan distribution

sub scan_distr
{
    # distribution
    my $distr = $_[0];

    # description of distribution
    my $distr_text = $_[1];

    # add =END TAG to text (for convenience)
    $distr_text .= "\n=END";

    # split into lines
    my @lines = split /\n/, $distr_text;

    # scan all TAGs (node sections)
    my $this_TAG = "=END";  # add =END tag to node TAG

    foreach my $l (@lines) {
	# next TAG ?
	if ($l =~ /^\s*(\/\*)?\s*(=[A-Z]+)\s*(.*)$/) {
	    # store next TAG
	    $this_TAG = $2;
	    # save rest of line
	    $l = $3;
	}

	# append to stored lines
	# (except for =END TAG)
	unless ($this_TAG eq "=END") {
	    $DISTR->{$distr}->{"=DOC"}->{$this_TAG} .= $l."\n";
	}
    }

} # end of scan_distr() 

# ----------------------------------------------------------------
# Get type of distribution and path of distribution file

sub get_distr_file
{
    # distribution
    my $distr = $_[0];

    my $found = 0;

    foreach my $type (keys %distr_types) {
	next unless $type;
	my $file_name = $distr_types{$type}{"file_prefix"}."\_$distr\.c";
	if ($file_list{$file_name}) {
	    $found = 1;
	    $DISTR->{$distr}->{"=FILE"} = $file_list{$file_name};
	    $DISTR->{$distr}->{"=TYPE"} = $type;
	    last;
	}
    }

    die "Cannot find file for $distr" unless $found;

} # end of get_distr_file()

# ----------------------------------------------------------------
# Read distribution file

sub read_distr_file
{
    # distribution
    my $distr = $_[0];

    # Read file
    my $file = $DISTR->{$distr}->{"=FILE"};
    open CFILE, $file or die "cannot open file $file\n";
    my $file_content;
    while (<CFILE>) {
	$file_content .= $_;
    }
    close CFILE;

    # Check distribution name
    $file_content =~ /static\s+const\s+char\s+distr\_name\s*\[\s*\]\s*=\s*\"$distr\"/
	or die "$distr: distr_name inconsistent";

    # Type of distribution
    my $type = $DISTR->{$distr}->{"=TYPE"};

    # Get PDF source
    my $PDF_name = $distr_types{$type}{"PDF_prefix"}.$distr;
    my $PDF_pattern = 
	"(int|double)\\s+"                   # $1: return type  
	.$PDF_name                           #     name of function
        ."\\s*\\(([^\\)]*)\\)\\s*"           # $2: arguments of function
	."([^;])"                            # $3: first character (to distinguish from prototype) 
	."(.*)"                              # $4: function body
        ."\\/\\*\\s+end\\s+of\\s+$PDF_name"; # end of function marker
    $file_content =~ /$PDF_pattern/s
	or die "cannot find PDF for $distr";

    # Store data
    $DISTR->{$distr}->{"=PDF"}->{"=NAME"}  = $PDF_name;  # name of PDF function
    $DISTR->{$distr}->{"=PDF"}->{"=RTYPE"} = $1;         # return type 
    $DISTR->{$distr}->{"=PDF"}->{"=ARGS"}  = $2;         # arguments for function
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  = $3.$4;      # function body

    # Modify function arguments:
    #   remove DISTR from argument list
    $DISTR->{$distr}->{"=PDF"}->{"=ARGS"}  =~
	s /\,\s*(UNUR_DISTR|struct unur_distr)\s*\*\s*distr\s*//;

    # Modify function body:
    #   remove comments
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s {/\*.*?\*/} []gsx;
    #   remove enclosing brackets
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /^\s*\{(.*)\}\s*$/$1/s;
    #   remove empty lines
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /\n\s*\n/\n/gx;
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /^\s*\n//;
    #   remove declaration of "params"
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ 
	s /.*(register)?\s+double\s*\*\s*params\W.*\n//;
    #   remove all `DISTR.' from body
    $DISTR->{$distr}->{"=PDF"}->{"=BODY"}  =~ s /DISTR\.//g;

    # Get parameters for PDF
    $file_content =~ s {/\*.*?\*/} []gsx;    # remove all comments
    my @lines = split /\n/, $file_content;
    my $n_params = -1;
    foreach my $l (@lines) {
	next unless $l =~ /\#define\s+(\w+)\s+(.*)$/;
	my $macro_name = $1;
	my $macro_body = $2;

	if ($macro_body =~ /params\s*\[(\d)+\]/) {
	    $n_params = ($n_params < $1) ? $1 : $n_params;
	    $DISTR->{$distr}->{"=PDF"}->{"=PARAMS"}[$1] = $macro_name;
	    next;
	}
	
	if ($macro_name =~ /(.*NORMCONSTANT)/) {
	    $DISTR->{$distr}->{"=PDF"}->{"=CONST"} = $macro_body;
	    next;
	}

	if ($macro_name =~ /(DISTR)/) {
	    $DISTR->{$distr}->{"=PDF"}->{"=DISTR"} = $macro_body;
	    next;
	}
    }

    # Number of parameters for PDF
    $DISTR->{$distr}->{"=PDF"}->{"=N_PARAMS"} = $n_params+1;

    # Id of distribution
    $file_content =~ /distr\-\>id\s*=\s*(\w+)/
	or die "cannot find ID for $distr";
    $DISTR->{$distr}->{"=ID"} = $1;
	
} # end of read_distr_file()

# ----------------------------------------------------------------
# Print data on screen (for debugging)

sub print_data
{
    # List of distributions
    my $DISTR = $_[0];

    foreach my $d (keys %{$DISTR}) {
	print "-------------------------------------\n";
	print "distribution = \"$d\"\n";

	print "NAME: ".$DISTR->{$d}->{"=NAME"}."\n";
	print "TYPE: ".$DISTR->{$d}->{"=TYPE"}."\n";
	print "ID  : ".$DISTR->{$d}->{"=ID"}."\n";
	print "FILE: ".$DISTR->{$d}->{"=FILE"}."\n\n";

	print "DOC: PDF    : ".$DISTR->{$d}->{"=DOC"}->{"=PDF"};
	print "DOC: CONST  : ".$DISTR->{$d}->{"=DOC"}->{"=CONST"};
	print "DOC: CDF    : ".$DISTR->{$d}->{"=DOC"}->{"=CDF"};
	print "DOC: DOMAIN : ".$DISTR->{$d}->{"=DOC"}->{"=DOMAIN"};
	print "DOC: FPARAM :\n".$DISTR->{$d}->{"=DOC"}->{"=FPARAM"}."\n";

	print "PDF: NAME  : ".$DISTR->{$d}->{"=PDF"}->{"=NAME"}."\n";
	print "PDF: RTYPE : ".$DISTR->{$d}->{"=PDF"}->{"=RTYPE"}."\n";
	print "PDF: ARGS  : ".$DISTR->{$d}->{"=PDF"}->{"=ARGS"}."\n";

	print "PDF: PARAM : ".$DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"}."\n";
	foreach my $i (0 .. $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"} - 1) {
	    print "\t[$i]: ".$DISTR->{$d}->{"=PDF"}->{"=PARAMS"}[$i]."\n"
	}

	print "PDF: CONST : ".$DISTR->{$d}->{"=PDF"}->{"=CONST"}."\n";
	print "PDF: DISTR : ".$DISTR->{$d}->{"=PDF"}->{"=DISTR"}."\n";

	print "PDF: BODY :\n".$DISTR->{$d}->{"=PDF"}->{"=BODY"}."\n";
    }

} # end of print_data()

# ----------------------------------------------------------------
# End of file
return 1;

# ----------------------------------------------------------------
