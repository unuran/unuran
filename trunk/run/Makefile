
## -lrngstream

all: run rungen test

run: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o run run.c ../src/libunuran.la -lm -lprng

mcrun: mcrun.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o mcrun mcrun.c ../src/libunuran.la -lm -lprng

rungen: rungen.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o rungen rungen.c ../src/libunuran.la -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o test test.c ../src/libunuran.la -lm -lprng 

clean:
	rm -f *~ *.o run rungen test unuran.log
	rm -rf .libs
