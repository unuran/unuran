#!/bin/bash

# check arguments 
test -z $1 && echo "Argument missing" && exit 1

VALGRINDOPTIONS="-v --tool=memcheck --leak-check=yes --leak-resolution=high --num-callers=40 --show-reachable=yes"

PROG=`echo $1 | sed -e "s#\./##"`

LOGFILENAME="valgrind-${PROG}"

echo "run valgrind on ${PROG} ..."

valgrind --logfile=${LOGFILENAME} ${VALGRINDOPTIONS} ./${PROG}

exit 0
