#!/bin/bash
# ---------------------------------------------------------------------------
# $id$
# ---------------------------------------------------------------------------

# home of mingw (for Gentoo Linux)
MINGWHOME="/opt/xmingw"
MINGWPATH="${MINGWHOME}/i386-mingw32msvc/bin"

export CC="${MINGWHOME}/bin/i386-mingw32msvc-gcc"
export STRIP="${MINGWPATH}/strip"
export AR="${MINGWPATH}/ar"
export RANLIB="${MINGWPATH}/ranlib"
export NM="${MINGWPATH}/nm"

# run configure in its own build directory
test -d win || mkdir win;
make distclean
(cd win && ../configure --srcdir=.. --prefix=/home/leydold/tmp/win --host=i386-pc-mingw32msvc --enable-shared)

# compile
(cd win && make)

exit 0;