
## -lrngstream

all: run rungen test

run: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O3 -o run run.c -lunuran -lm -lprng

rundmalloc: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O3 -o run run.c -lunuran -lm -lprng -ldmalloc

rungen: rungen.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o rungen rungen.c -lunuran -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o test test.c -lunuran -lm -lprng 

clean:
	rm -f *~ *.o run run_acg rungen test unuran.log
