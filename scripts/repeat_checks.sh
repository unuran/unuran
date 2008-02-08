#!/bin/bash
# -----------------------------------

# directory for storing results and log files
RESULTS=checkresults

# -----------------------------------

usage() {
cat << "USAGE_END"
Usage: repeat_checks.sh [options] <rep> 

Repeat checks with different seeds <rep> times

Options:
  -f             run in fullcheck mode
  -h, --help     help (print this page)

USAGE_END
        exit $1
}

# -----------------------------------
# Read command line arguments

for x in "$@" ; do
    case "${x}" in
        -f)        # fullcheck mode
	    export UNURANFULLCHECK=true
	    ;;
        -h|--help) # helpfullcheck mode
	    usage 0
	    ;;
        -*)  # invalid option
	    echo "$0: Invalid switch \"${x}\"!."
	    usage 1
	    ;;
        *)   # this should be the number of repetitions
             rep="${x}"
             ;;
    esac
done

# -----------------------------------
# number of repetitions

if [[ -z "$rep" ]] ; then rep=10; fi

echo "number of repetitions = $rep"

# -----------------------------------
# make

if [[ -x ../hmake ]] ; then
    MAKE=../hmake
else
    MAKE=make
fi

# -----------------------------------
# make directory with test results

if [[ -d "$RESULTS" ]] ; then rm -rf "$RESULTS"; fi
mkdir -v "$RESULTS"

# -----------------------------------
# run ...

for (( i=0; i<$rep; i++)); do
    SEED=`echo $RANDOM$RANDOM | cut -b1-8`
    echo $i:$SEED
    mkdir -v "$RESULTS/$i"
    echo "SEED = ${SEED}" > $RESULTS/$i/ENV
    if [[ -n "${UNURANFULLCHECK}" ]] ; then
	echo "check mode: fullcheck" >> $RESULTS/$i/ENV
    else
	echo "check mode: installation" >> $RESULTS/$i/ENV
    fi
    $MAKE check > $RESULTS/$i/RESULTS
    echo "copy log file"
    mv *.log $RESULTS/$i
done

# -----------------------------------
# end

exit 0

# -----------------------------------
