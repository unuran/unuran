
## -lrngstream

all: run rungen test

run: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o run run.c ../src/libunuran.la -lm -lprng

rundmalloc: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o run run.c ../src/libunuran.la -lm -lprng -ldmalloc

rungen: rungen.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o rungen rungen.c ../src/libunuran.la -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o test test.c ../src/libunuran.la -lm -lprng 

clean:
	rm -f *~ *.o run run_acg rungen test unuran.log
	rm -rf .libs
