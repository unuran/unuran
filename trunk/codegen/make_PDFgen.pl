#!/usr/bin/perl
# ----------------------------------------------------------------
# Make files for code generator
# ----------------------------------------------------------------
# $Id$
# ----------------------------------------------------------------

$DEBUG = 0;

# ----------------------------------------------------------------

require "read_PDF.pl";

# ----------------------------------------------------------------

# C file for code generator  
my $PDFgen_C_file = "PDFgen.c";
my $PDFgen_H_file = "PDFgen_source.h";

# ----------------------------------------------------------------
# List of distributions
my $DISTR = read_PDFdata();

# For description of data fields in this list see file `readPDF.pl'.

# ................................................................

# Print result on screen
    if ($DEBUG) {
	print_data($DISTR);
    }

# ----------------------------------------------------------------

# Make C file with code generator
make_PDFgen_files();

# ----------------------------------------------------------------
# End

exit 0;

# ----------------------------------------------------------------
#
# Subroutines
#
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# Make C file with code generator

sub make_PDFgen_files
{
    my $PDFgen;
    my $PDFgen_static_prototypes;
    my $PDFgen_global_prototypes;
    
    my $empty_line = "\tfprintf (out,\"\\n\");\n\n";

# ................................................................
# Main

    # Mark begin of Main
    $PDFgen .= "/* ------------------------------------- */\n";
    $PDFgen .= "/* Main */\n\n";
    
    # Function name
    my $function = "int \_unurgen\_C\_PDF (UNUR_DISTR *distr, FILE *out, const char *pdf)";
    
    # Function prototype
    $PDFgen_global_prototypes .= "$function;\n";
    
    # Function header
    $PDFgen .= "$function\n{\n";

    # Check for invalid NULL pointer
    $PDFgen .= "\t_unur_check_NULL(\"unurgen\", distr, 0 );\n";

    # Switch
    $PDFgen .= "\tswitch (distr->id) {\n";
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";
	$PDFgen .= "\tcase ".$DISTR->{$d}->{"=ID"}.":\n";
	$PDFgen .= "\t\treturn \_unurgen\_C\_PDF_$d (distr,out,pdf);\n";
    }
    $PDFgen .= "\tdefault:\n";
    $PDFgen .= "\t\treturn 0;\n";
    $PDFgen .= "\t}\n";
    
    # End of function
    $PDFgen .= "}\n\n";

    # Mark end of Main
    $PDFgen .= "/* end of Main */\n";
    $PDFgen .= "/* ------------------------------------- */\n";

# ................................................................
# Continuous distributions
    
    foreach my $d (sort keys %{$DISTR}) {
	next unless $DISTR->{$d}->{"=TYPE"} eq "CONT";

	# Mark begin of distribution
	$PDFgen .= "/* ------------------------------------- */\n";
	$PDFgen .= "/* $d: ".$DISTR->{$d}->{"=NAME"}." */\n\n";

	# Function name
	my $function = "int \_unurgen\_C\_PDF\_$d (UNUR_DISTR *distr, FILE *out, const char *pdf)";

	# Function prototype
	$PDFgen_static_prototypes .= "static $function;\n";

	# Function header
	$PDFgen .= "$function\n{\n";

	# Macro definitions for parameters
	#   (total) number of parameters
	$PDFgen .= "\tfprintf (out,\"\#define  n_params  %d\\n\","
	    .$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".n_params);\n";
	#   parameter list
	foreach my $i (0 .. $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"} - 1) {
	    $PDFgen .= "\tfprintf (out,\"\#define  "
		.$DISTR->{$d}->{"=PDF"}->{"=PARAMS"}[$i]
		."  (%.20g)\\n\","
	        .$DISTR->{$d}->{"=PDF"}->{"=DISTR"}.".params[$i]);\n";
	}
	#   normalization constant
	if ($DISTR->{$d}->{"=PDF"}->{"=CONST"}) {
	    if ($DISTR->{$d}->{"=PDF"}->{"=BODY"} =~ /LOGNORMCONSTANT/) {
		$PDFgen .= "\tfprintf (out,\"\#define  LOGNORMCONSTANT  (%.20g)\\n\","
		    .$DISTR->{$d}->{"=PDF"}->{"=CONST"}.");\n";
	    }
	    elsif ($DISTR->{$d}->{"=PDF"}->{"=BODY"} =~ /NORMCONSTANT/) {
		$PDFgen .= "\tfprintf (out,\"\#define  NORMCONSTANT  (%.20g)\\n\","
		    .$DISTR->{$d}->{"=PDF"}->{"=CONST"}.");\n";
	    }
	}
	$PDFgen .= $empty_line;

	# compose PDF name
	my $PDFname = "\tfprintf (out,\"static ";
	$PDFname .= $DISTR->{$d}->{"=PDF"}->{"=RTYPE"};
	$PDFname .= " %s (";
	$PDFname .= $DISTR->{$d}->{"=PDF"}->{"=ARGS"};
	$PDFname .= ")\\n{\\n\",";
	$PDFname .= " ((pdf) ? pdf : \"pdf_$d\") ";
	$PDFname .= ");\n";

	# Write PDF
	$PDFgen .= $PDFname; 
	foreach my $l (split /\n/, $DISTR->{$d}->{"=PDF"}->{"=BODY"}) {
	    $PDFgen .= "\tfprintf (out,\"$l\\n\");\n";
	}

	$PDFgen .= "\tfprintf (out,\"}\\n\");\n";
	$PDFgen .= $empty_line;

	# Undefine macros
	#   (total) number of parameters
	$PDFgen .= "\tfprintf (out,\"\#undef  n_params\\n\");\n";
	#   parameter list
	foreach my $i (0 .. $DISTR->{$d}->{"=PDF"}->{"=N_PARAMS"} - 1) {
	    $PDFgen .= "\tfprintf (out,\"\#undef  "
		.$DISTR->{$d}->{"=PDF"}->{"=PARAMS"}[$i]."\\n\");\n";
	}
	#   normalization constant
	if ($DISTR->{$d}->{"=PDF"}->{"=CONST"}) {
	    if ($DISTR->{$d}->{"=PDF"}->{"=BODY"} =~ /LOGNORMCONSTANT/) {
		$PDFgen .= "\tfprintf (out,\"\#undef  LOGNORMCONSTANT\\n\");\n";
	    }
	    elsif ($DISTR->{$d}->{"=PDF"}->{"=BODY"} =~ /NORMCONSTANT/) {
		$PDFgen .= "\tfprintf (out,\"\#undef  NORMCONSTANT\\n\");\n";
	    }
	}

	# End of function
	$PDFgen .= "\n";
	$PDFgen .= "\treturn 1;\n";
	$PDFgen .= "}\n\n";

	# Mark end of distribution
	$PDFgen .= "/* end of $d */\n";
	$PDFgen .= "/* ------------------------------------- */\n";
    }

# ................................................................
# Print C code int files

    # Header file
    open HFILE, ">$PDFgen_H_file" or die "cannot open file $PDFgen_H_file\n";
    print HFILE $PDFgen_global_prototypes;
    close HFILE;

    # C file
    open CFILE, ">$PDFgen_C_file" or die "cannot open file $PDFgen_C_file\n";
    print CFILE "\#include <source_unuran.h>\n";
    print CFILE "\#include \"PDFgen_source.h\"\n\n";
    print CFILE $PDFgen_static_prototypes;
    print CFILE $PDFgen;
    close CFILE;

} # end of make_code_generator()

# ----------------------------------------------------------------
