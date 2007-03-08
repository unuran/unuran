#! sh

# Where to find MS Visual Studio 2005 compilers
export VSCYG='/cygdrive/c/Programme/Microsoft Visual Studio 8'
export VCUNIX='c:/Programme/Microsoft Visual Studio 8/VC'
export VCWIN='C:\Programme\Microsoft Visual Studio 8\VC'

# Add path to compiler
export PATH=${VSCYG}/VC/bin:${PATH}
export PATH=${VSCYG}/Common7/IDE:${PATH}

# C prepsocessor, compiler and linker
export CPP="cl -E"
export CC="cl "
export CXXCPP="cl -E"
export CXX="cl "
export LD="cl "

# Add path to MSVC system header files
export INCLUDE="${VCUNIX}/include;${VCWIN}\include;${INCLUDE}"
export INCLUDE="${VCUNIX}/PlatformSDK/Include;${VCWIN}\PlatformSDK\Include;${INCLUDE}"

# Add path to MSVC libraries
export LIB="${VCUNIX}/lib;${VCWIN}\lib;${LIB}"
export LIB="${VCUNIX}/PlatformSDK/Lib;${VCWIN}\PlatformSDK\Lib;${LIB}"

#echo ${INCLUDE}

#./autogen.sh --disable-shared --enable-static
#./hmake

./hmake -f scripts/Makefile.win32



