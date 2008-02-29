#!/bin/bash
# -----------------------------------

# directory for storing results and log files
RESULTS=RESULTS

RULER="-------------------------------"

MAKE="make -j2"

# -----------------------------------

usage() {
cat << "USAGE_END"
Usage: repeat_checks.sh [options]

Options:
  -n=N           repeat N times [default: N=10]
  -t=TEST        run particular program TEST [default: run all tests]
  -f             run in fullcheck mode
  -h, --help     help (print this page)

USAGE_END
        exit $1
}

# -----------------------------------
# Read command line arguments

for x in "$@" ; do
    case "${x}" in
	-n=*) # number of repetitions
	    rep=`echo "${x}" | sed -e "s/-n=//"`
	    ;;
	-t=*) # name of test program
	    testname=`echo "${x}" | sed -e "s/-t=//"`
	    ;;
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
        *)   # invalid argument
	    echo "$0: Invalid argument \"${x}\"!."
	    usage 1
	    ;;
    esac
done

# -----------------------------------
# number of repetitions

if [[ -z "$rep" ]] ; then rep=10; fi

echo "number of repetitions = $rep"

# -----------------------------------
# make directory with test results

if [[ -d "$RESULTS" ]] ; then rm -rf "$RESULTS"; fi
mkdir -v "$RESULTS"

# -----------------------------------
# check test program

if [[ -n "$testname" ]] ; then 
    $MAKE $testname >> $RESULTS/RESULTS 2>&1
    if [[ ! -x "$testname" ]] ; then 
	echo "invalid test: \"$testname\""
	usage 1
    fi
fi

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
    ## create directory for results
    mkdir -v "$RESULTS/$i"
    ## create random seed
    SEED=`echo $RANDOM$RANDOM | cut -b1-8`
    echo $i:$SEED
    date
    echo "SEED = ${SEED}" > $RESULTS/$i/ENV
    ## check mode
    if [[ -n "${UNURANFULLCHECK}" ]] ; then
	echo "check mode: fullcheck" >> $RESULTS/$i/ENV
    else
	echo "check mode: installation" >> $RESULTS/$i/ENV
    fi
    ## set environment variables for tests
    export SEED="${SEED}"
    export UNURANSTOPWATCH=true
    ## run tests
    if [[ -z "$testname" ]] ; then 
	time ($MAKE check >> $RESULTS/$i/RESULTS) >> $RESULTS/$i/TIMING 2>&1
    else
	time ("./$testname" >> $RESULTS/$i/RESULTS) >> $RESULTS/$i/TIMING 2>&1
	if [ $? != 0 -a $? != 77 ] ; then
	    echo "FAIL: $testname" >> $RESULTS/$i/RESULTS;
	fi
    fi
    ## move log files
    mv *.log $RESULTS/$i
    #### (cd $RESULTS/$i; grep "PASS: " RESULTS | sed -e 's/PASS: //' | sed -e 's/$/_test.log/' | xargs -i= rm "=")
    (cd $RESULTS/$i; grep "PASS: " RESULTS | sed -e 's/PASS: //' | sed -e 's/$/_unuran.log/' | xargs -i= rm "=")
    ## store test results 
    echo "#${i}: SEED=${SEED}" >> $RESULTS/RESULTS
    grep "real" $RESULTS/$i/TIMING  | sed -e 's/real/time =/' >> $RESULTS/RESULTS
    tail -n 10 $RESULTS/$i/RESULTS | grep "[0-9]\+ of [0-9]\+ tests failed" >> $RESULTS/RESULTS
    tail -n 10 $RESULTS/$i/RESULTS | grep "All [0-9]\+ tests passed" >> $RESULTS/RESULTS
    grep "FAIL: " $RESULTS/$i/RESULTS | sed -e 's/FAIL:/   FAIL:/' >> $RESULTS/RESULTS
    echo "${RULER}" >> $RESULTS/RESULTS
    ## print info on screen
    echo `grep "FAIL: " $RESULTS/$i/RESULTS | sed -e 's/FAIL:/   FAIL:/'`
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
