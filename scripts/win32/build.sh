#! sh

# --- Default building flags ------------------------------------------------

#UNURAN_CONFIGURE_FLAGS="--enable-maintainer-mode --disable-deprecated"

#UNURAN_CONFIGURE_FLAGS="--enable-maintainer-mode --disable-deprecated --enable-logfile --enable-debug-struct"

UNURAN_CONFIGURE_FLAGS="--enable-maintainer-mode --disable-deprecated"

# --- Setup environment -----------------------------------------------------

# Where to find MS Visual Studio 2005 compilers
export VSCYG='/cygdrive/c/Programme/Microsoft Visual Studio 8'
export VCUNIX='c:/Programme/Microsoft Visual Studio 8/VC'
export VCWIN='C:\Programme\Microsoft Visual Studio 8\VC'

# Add path to compiler
export PATH=${VSCYG}/VC/bin:${PATH}
export PATH=${VSCYG}/Common7/IDE:${PATH}

# Add path to MSVC system header files
export INCLUDE="${VCUNIX}/include;${VCWIN}\include;${INCLUDE}"
export INCLUDE="${VCUNIX}/PlatformSDK/Include;${VCWIN}\PlatformSDK\Include;${INCLUDE}"

# Add path to MSVC libraries
export LIB="${VCUNIX}/lib;${VCWIN}\lib;${LIB}"
export LIB="${VCUNIX}/PlatformSDK/Lib;${VCWIN}\PlatformSDK\Lib;${LIB}"

# C prepsocessor, compiler and linker
export CPP="cl -E"
export CC="cl "
export CXXCPP="cl -E"
export CXX="cl "
export LD="cl "

# temporary directory for windows files
export WIN_DIR="win32"

# add path to created DLL (for testing)
export PATH=`pwd`/${WIN_DIR}:${PATH}


# --- Read command line arguments -------------------------------------------

for arg in "$@"; do
	case "${arg}" in
	all)		# create manual
		doc=true
		;;
	examples)	# Compile, link and run all examples
		examples=true
		;;
	check)		# Compile, link and run all tests
		check=true;
		export privatefunctions=true
		;;
	clean)		# Remove all files created by this script
		clean=true;
		;;
	--*)		# ./configure option
		UNURAN_CONFIGURE_FLAGS="${UNURAN_CONFIGURE_FLAGS} ${arg}"
		;;
	*)
		echo "Invalid argment '${arg}'"
		exit 1
		;;
	esac
done


# --- Check working directory -----------------------------------------------

# test whether the script is started from top source directory
if [[ ! -f ./configure.ac || -z `grep unuran ./configure.ac` ]]; then
	echo "You must run this script from UNU.RAN top source directory";
	exit 1;
fi

# test whether script start in cygwin environment
if [[ -z `uname | grep -i cygwin` ]]; then
	echo "This scripts requires CYGWIN";
	exit 1;
fi

# test whether config.h was created using the MSVC compiler
if [[ -f ./config.h && -f ./Makefile && -z `grep "CC = cl" ./Makefile` ]]; then
	echo "config.h created with wrong CC; deleting ..."
	make maintainer-clean
fi

# --- Clear working directory -----------------------------------------------

if [[ "${clean}" ]]; then
	make -f scripts/win32/Makefile.win32 clean;
	echo "working space cleared";
	exit 0;
fi

# --- Prepare UNU.RAN -------------------------------------------------------

# create 'config.h' and Makefiles using autotools
if [[ !( -f ./configure ) ]]; then
	autoreconf -i
fi

if [[ !( -f ./config.h && -f ./Makefile) ]]; then
	./configure ${UNURAN_CONFIGURE_FLAGS}
fi

# create directory for windows files
test -d "${WIN_DIR}" && rm -rf "${WIN_DIR}"
mkdir "${WIN_DIR}"

# create all required UNU.RAN header files 
(cd src/parser; make stringparser_lists.ch)
(cd src; make unuran.h; cp -v unuran.h ../${WIN_DIR})

# create doc
if [[ "${doc}" ]]; then
	(cd doc; make unuran.pdf; cp -v unuran.pdf ../${WIN_DIR});
fi

# --- Create DLL ------------------------------------------------------------

make -f scripts/win32/Makefile.win32

# --- Compile, link and run all examples and test files

if [[ "${examples}" ]]; then
	make -f scripts/win32/Makefile.win32 examples;
fi

if [[ ${check} ]]; then
	make -f scripts/win32/Makefile.win32 check;
fi

# --- Done ------------------------------------------------------------------

exit 0

