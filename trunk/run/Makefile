
## -lrngstream

all: run rungen test

something: something.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o something something.c ../src/libunuran.la -lm -lprng

run: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o run run.c ../src/libunuran.la -lm -lprng

mcrun: mcrun.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o mcrun mcrun.c ../src/libunuran.la -lm -lprng

rungen: rungen.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o rungen rungen.c ../src/libunuran.la -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o test test.c ../src/libunuran.la -lm -lprng 

something: something.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o something something.c ../src/libunuran.la -lm -lprng 

f: ../src/utils/fminmax.c fmain.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -ansi -pedantic -o f fmain.c ../src/libunuran.la -lm -lprng 

try: try.c try_source.h 
	gcc -lm -ansi -pedantic -DHAVE_IEEE try.c -o try 

clean:
	rm -f *~ *.o run rungen mcrun test something unuran.log
	rm -rf .libs
