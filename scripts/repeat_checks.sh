#!/bin/bash
# -----------------------------------

# directory for storing results and log files
RESULTS=RESULTS

RULER="-------------------------------"

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

MAKE="make -j2"

# -----------------------------------
# make directory with test results

if [[ -d "$RESULTS" ]] ; then rm -rf "$RESULTS"; fi
mkdir -v "$RESULTS"

# -----------------------------------
# start

date
echo "${RULER}" >> $RESULTS/RESULTS
echo -n "started: " >> $RESULTS/RESULTS
date >> $RESULTS/RESULTS
echo "${RULER}" >> $RESULTS/RESULTS

# -----------------------------------
# run ...

for (( i=0; i<$rep; i++)); do
    SEED=`echo $RANDOM$RANDOM | cut -b1-8`
    echo $i:$SEED
    date
    mkdir -v "$RESULTS/$i"
    echo "SEED = ${SEED}" > $RESULTS/$i/ENV
    if [[ -n "${UNURANFULLCHECK}" ]] ; then
	echo "check mode: fullcheck" >> $RESULTS/$i/ENV
    else
	echo "check mode: installation" >> $RESULTS/$i/ENV
    fi
    export SEED="${SEED}"
    export UNURANSTOPWATCH=true
    time ($MAKE check > $RESULTS/$i/RESULTS) > $RESULTS/$i/TIMING 2>&1
    echo "move log files"
    mv *.log $RESULTS/$i
    echo "#${i}: SEED=${SEED}" >> $RESULTS/RESULTS
    grep "real" $RESULTS/$i/TIMING  | sed -e 's/real/time =/' >> $RESULTS/RESULTS
    tail -n 10 $RESULTS/$i/RESULTS | grep "[0-9]\+ of [0-9]\+ tests failed" >> $RESULTS/RESULTS
    tail -n 10 $RESULTS/$i/RESULTS | grep "All [0-9]\+ tests passed" >> $RESULTS/RESULTS
    grep "FAIL: " $RESULTS/$i/RESULTS | sed -e 's/FAIL:/   FAIL:/' >> $RESULTS/RESULTS
    echo "${RULER}" >> $RESULTS/RESULTS
    ## (cd $RESULTS/$i; grep "PASS: " RESULTS | sed -e 's/PASS: //' | sed -e 's/$/_test.log/' | xargs -i= rm "=")
    (cd $RESULTS/$i; grep "PASS: " RESULTS | sed -e 's/PASS: //' | sed -e 's/$/_unuran.log/' | xargs -i= rm "=")
done

# -----------------------------------
# stop

echo -n "stopped: " >> $RESULTS/RESULTS
date >> $RESULTS/RESULTS
echo "${RULER}" >> $RESULTS/RESULTS
date

# -----------------------------------
# end

exit 0

# -----------------------------------
