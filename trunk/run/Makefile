 

all: run run_acg test parsetest

run: run.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o run run.c -lunuran -lm -lprng

run_acg: run_acg.c ../src/unuran.h ../codegen/libunuracg.a ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -I../codegen -L../src -L../codegen -Wall -g -O2 -o run_acg run_acg.c -lunuracg -lunuran -lm -lprng

test: test.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o test test.c -lunuran -lm -lprng 

getopttst: getopttst.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o getopttst getopttst.c -lunuracg -lunuran -lm -lprng 

acg: acg.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o acg acg.c -lunuracg -lunuran -lm -lprng 

parsetest: parsetest.c stringpars.c.base
	./gen_stringpars.pl
	gcc -DHAVE_CONFIG_H -Wall -I. -I.. -I../src -I../src/methods -I../src/distributions -L../src -o parsetest parsetest.c stringpars.c -lunuran -lm -lprng

clean:
	rm -f *~ *.o run test run_acg getopttst unuran.log acg
	rm -f uniftest_prng.c uniftest_prng
	rm -f uniftest_c.c uniftest_c
	rm -f uniftest_f.f uniftest_f
	rm -f uniftest_java.java Uniftest.class Urand.class
	rm -f stringpars.c
