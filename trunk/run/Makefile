
## -lrngstream

all: run rungen test something

run: run.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o run run.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm 

mcrun: mcrun.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O3 -o mcrun mcrun.c ../src/libunuran.la -lm -lprng

rungen: rungen.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o rungen rungen.c ../src/libunuran.la -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -o test test.c ../src/libunuran.la -lm -lprng 

q: q.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o q q.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm 

something: something.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o something something.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm 

matr: matr.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la ../tests/testdistributions/libtestdistributions.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o matr matr.c ../src/experimental/libunurexperimental.la ../tests/testdistributions/libtestdistributions.la ../src/libunuran.la -lm 

spheres: spheres.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la ../tests/testdistributions/libtestdistributions.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o spheres spheres.c ../src/experimental/libunurexperimental.la ../tests/testdistributions/libtestdistributions.la ../src/libunuran.la -lm 

lotto: lotto.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o lotto lotto.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm -lprng 

steps: steps.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o steps steps.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm 

other: other.c ../src/unuran.h ../src/libunuran.la ../src/experimental/libunurexperimental.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o other other.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm 
		
hit-shared: hit.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -g -O2 -shared hit.c -o hit.so ../src/.libs/libunuran.so -lm 

hit: hit.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -Wall -DHAVE_CONFIG_H -g -O2 -m32 -o hit hit.c ../src/experimental/libunurexperimental.la ../src/libunuran.la -lm 

f: ../src/utils/fminmax.c fmain.c ../src/unuran.h ../src/libunuran.la
	../libtool gcc -I. -I.. -I../src -I../src/tests -ansi -pedantic -o f fmain.c ../src/libunuran.la -lm -lprng 

try: try.c try_source.h 
	gcc -lm -ansi -pedantic -DHAVE_IEEE try.c -o try 

clean:
	rm -f *~ *.o run rungen mcrun test hit something f matr lotto try unuran.log
	rm -rf .libs
