
## -lrngstream

all: run run_acg rungen test

run: run.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O3 -o run run.c -lunuran -lm -lprng -ldmalloc

run_acg: run_acg.c ../src/unuran.h ../codegen/libunuracg.a ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -I../codegen -L../src -L../codegen -Wall -g -O2 -o run_acg run_acg.c -lunuracg -lunuran -lm -lprng

rungen: rungen.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o rungen rungen.c -lunuran -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o test test.c -lunuran -lm -lprng 

clean:
	rm -f *~ *.o run run_acg rungen test unuran.log
