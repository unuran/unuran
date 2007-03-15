#! sh

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
	--all)		# create manual
		doc=true;;
	--examples)	# Compile and link all examples and tests
		examples=true;;
	--nix)
		export check=true;
		nix=true;;
	*)
		echo "Invalid argment '${arg}'"
		exit 1;;
	esac
done

# --- Check working directory -----------------------------------------------

# test whether the script is started from top source directory
if [[ ! -f ./configure.ac || -z `grep unuran ./configure.ac` ]]; then
	echo "You must run this script from UNURAN top source directory";
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
	./hmake maintainer-clean
fi

# --- Prepare UNU.RAN -------------------------------------------------------

# create 'config.h' and Makefiles using autotools
if [[ !( -f ./config.h && -f ./Makefile) ]]; then
	./autogen.sh;
fi

# create directory for windows files
test -d "${WIN_DIR}" && rm -rf "${WIN_DIR}"
mkdir "${WIN_DIR}"

# create all required UNURAN header files 
(cd src/parser; ../../hmake stringparser.c)
(cd src; ../hmake unuran.h; cp -v unuran.h ../${WIN_DIR})

# create doc
if [[ "${doc}" ]]; then
	(cd doc; ../hmake unuran.pdf; cp -v unuran.pdf ../${WIN_DIR});
fi

# --- Create DLL ------------------------------------------------------------

./hmake -f scripts/win32/Makefile.win32

# --- Compile, link and run all examples and test files

if [[ "${examples}" ]]; then
	./hmake -f scripts/win32/Makefile.win32 examples;
fi



if [[ ${nix} ]]; then
	./hmake -f scripts/win32/Makefile.win32 nix;
fi





# --- Done ------------------------------------------------------------------

exit 0

