 

all: run run_acg test

run: run.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O3 -o run run.c -lunuran -lm -lprng -lrngstream

run_timing: run_timing.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O3 -o run_timing run_timing.c -lunuran -lm -lprng -lrngstream

run_acg: run_acg.c ../src/unuran.h ../codegen/libunuracg.a ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -I../codegen -L../src -L../codegen -Wall -g -O2 -o run_acg run_acg.c -lunuracg -lunuran -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o test test.c -lunuran -lm -lprng 

clean:
	rm -f *~ *.o run test run_timing run_acg unuran.log acg
