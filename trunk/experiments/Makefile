 

all: rungen

rungen: rungen.c ../src/unuran.h ../src/libunuran.a
	gcc -I. -I.. -I../src -I../src/tests -L../src -Wall -g -O2 -o rungen rungen.c -lunuran -lm -lprng

clean:
	rm -f *~ *.o rungen unuran.log
