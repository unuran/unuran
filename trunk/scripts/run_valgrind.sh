#!/bin/bash

# check arguments 
test -z $1 && echo "Argument missing" && exit 1

test -f suppressions-valgrind && SUPPRESS="--suppressions=suppressions-valgrind"

VALGRINDOPTIONS="-v --tool=memcheck --leak-check=yes --leak-resolution=high --num-callers=40 --show-reachable=yes $SUPPRESS"

PROG=`echo $1 | sed -e "s#\./##"`

LOGFILENAME="valgrind-${PROG}"

echo "run valgrind on ${PROG} ..."

valgrind --logfile=${LOGFILENAME} ${VALGRINDOPTIONS} ./${PROG}

exit 0
