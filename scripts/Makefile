OPTIONS=-g -pg -O3
all:	tmp.x

tmp.x:	Franck.o Neighbours.o main.o global.o MC.o Parameters.o
	g++ ${OPTIONS}  -o tmp.x MC.o Franck.o Neighbours.o global.o main.o Parameters.o

clear: 
	rm -f *.o

Parameters.o: Parameters.cc
	g++ ${OPTIONS}   -c Parameters.cc


Franck.o: code.cc
	g++ ${OPTIONS}  -o Franck.o -c code.cc

MC.o: MC.cc
	g++ ${OPTIONS}  -c MC.cc

Neighbours.o: Neighbours.cc
	g++ ${OPTIONS}  -c Neighbours.cc

main.o: main.cc
	g++ ${OPTIONS}  -c main.cc
 
global.o: global.cc
	g++ ${OPTIONS}  -c global.cc
