#! /bin/sh
# Run this to generate all the initial makefiles, etc.

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.

ORIGDIR=`pwd`
cd $srcdir
PROJECT=unuran
TEST_TYPE=-f
FILE=src/unuran_config.h          # a file that should exist in the source dir
export WANT_AUTOMAKE=1.7

DIE=0

(autoconf --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have autoconf installed to compile $PROJECT."
	DIE=1
}

(libtool --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have libtool installed to compile $PROJECT."
	DIE=1
}

(automake --version) < /dev/null > /dev/null 2>&1 || {
	echo
	echo "You must have automake installed to compile $PROJECT."
	DIE=1
}

if test "$DIE" -eq 1; then
	exit 1
fi

test $TEST_TYPE $FILE || {
	echo "You must run this script in the top-level $PROJECT directory"
	exit 1
}

if test -z "$*"; then
	echo "I am going to run ./configure with no arguments - if you wish "
        echo "to pass any to it, please specify them on the $0 command line."
fi

case $CC in
*xlc | *xlc\ * | *lcc | *lcc\ *) am_opt=--include-deps;;
esac

#if test -z "$ACLOCAL_FLAGS"; then
#
#	acdir=`aclocal --print-ac-dir`
#        m4list="glib.m4 gettext.m4"
#
#	for file in $m4list
#	do
#		if [ ! -f "$acdir/$file" ]; then
#			echo "WARNING: aclocal's directory is $acdir, but..."
#			echo "         no file $acdir/$file"
#			echo "         You may see fatal macro warnings below."
#			echo "         If these files are installed in /some/dir, set the ACLOCAL_FLAGS "
#			echo "         environment variable to \"-I /some/dir\", or install"
#			echo "         $acdir/$file."
#			echo ""
#		fi
#	done
#fi

#echo "Running gettextize...  Ignore non-fatal messages."
# Hmm, we specify --force here, since otherwise things dont'
# get added reliably, but we don't want to overwrite intl
# while making dist.
#echo "no" | gettextize --copy --force

libtoolize --automake

aclocal $ACLOCAL_FLAGS

# optionally feature autoheader
(autoheader --version)  < /dev/null > /dev/null 2>&1 && autoheader

automake --add-missing $am_opt
autoconf -Wobsolete
cd $ORIGDIR

$srcdir/configure "$@"

echo 
echo "Now type 'make' to compile $PROJECT."





