 

all: run run_acg test parsetest

run: run.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o run run.c -lunuran -lm -lprng

run_acg: run_acg.c ../src/unuran.h ../codegen/libunuracg.a ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -I../codegen -L../src -L../codegen -Wall -g -O2 -o run_acg run_acg.c -lunuracg -lunuran -lm -lprng

fstrmain: fstrmain.c fstr.c 
	gcc -DHAVE_CONFIG_H -I. -I.. -I../src -I../src/methods -I../src/distributions -I../src/parser -I../src -Wall -g -O2 -o fstrmain fstrmain.c fstr.c -lm

test: test.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o test test.c -lunuran -lm -lprng 

getopttst: getopttst.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o getopttst getopttst.c -lunuracg -lunuran -lm -lprng 

acg: acg.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o acg acg.c -lunuracg -lunuran -lm -lprng 

parsetest: parsetest.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o parsetest parsetest.c -lunuran -lm -lprng

clean:
	rm -f *~ *.o run test run_acg fstrmain unuran.log acg
	rm -f uniftest_prng.c uniftest_prng
	rm -f uniftest_c.c uniftest_c
	rm -f uniftest_f.f uniftest_f
	rm -f uniftest_java.java Uniftest.class Urand.class
	rm -f parsetest stringpars.c
