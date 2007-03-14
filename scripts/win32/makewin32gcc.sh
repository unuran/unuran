#!/bin/bash
# ---------------------------------------------------------------------------
# file for makeing shared libraries using GCC
# ---------------------------------------------------------------------------
# I don't whether it works correctly.
# I have a new script scripts/win32/build.sh that uses CYGWIN and 
# MS Visual Studio that works.
# Thus this script is disabled.

echo "Please use scripts/win32/build.sh"
exit 1

# ---------------------------------------------------------------------------

# mingw 
MINGW=/usr/bin/mingw32

export CC="${MINGW}-gcc"
export STRIP="${MINGW}-strip"
export AR="${MINGW}-ar"
export AS="${MINGW}-as"
export RANLIB="${MINGW}-ranlib"
export NM="${MINGW}-nm"
export LD="${MINGW}-ld"
export DLLTOOL="${MINGW}-dlltool"

# run configure in its own build directory
test -d win || mkdir win;
make distclean
(cd win && ../configure --srcdir=.. --prefix=${HOME}/tmp/win --host=i386-pc-mingw32 --enable-shared)

# compile
(cd win && ../hmake)

exit 0;
