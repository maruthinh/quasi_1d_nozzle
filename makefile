#declare the variable
CC=g++

CFLAGS=-c -Wall

all: nozzle

nozzle:  memory2d.o iniflow.o bcond.o rgrid.o inigrid.o tstep.o lrroe.o lrmovers.o fluxroe.o fluxmovers.o solver.o srcterm.o minmodlim.o conver.o wsolut.o main.o 
	$(CC) memory2d.o iniflow.o bcond.o rgrid.o inigrid.o tstep.o lrroe.o lrmovers.o fluxroe.o fluxmovers.o solver.o srcterm.o minmodlim.o conver.o wsolut.o main.o -o nozzle

memory2d.o: memory2d.cpp
	$(CC) $(CFLAGS) memory2d.cpp

rgrid.o: rgrid.cpp
	$(CC) $(CFLAGS) rgrid.cpp

iniflow.o: iniflow.cpp
	$(CC) $(CFLAGS) iniflow.cpp

bcond.o: bcond.cpp
	$(CC) $(CFLAGS) bcond.cpp

inigrid.o: inigrid.cpp
	$(CC) $(CFLAGS) inigrid.cpp

tstep.o: tstep.cpp
	$(CC) $(CFLAGS) tstep.cpp

lrroe.o: lrroe.cpp
	$(CC) $(CFLAGS) lrroe.cpp

lrmovers.o: lrmovers.cpp
	$(CC) $(CFLAGS) lrmovers.cpp

fluxroe.o: fluxroe.cpp
	$(CC) $(CFLAGS) fluxroe.cpp

fluxmovers.o: fluxmovers.cpp
	$(CC) $(CFLAGS) fluxmovers.cpp

solver.o: solver.cpp
	$(CC) $(CFLAGS) solver.cpp

srcterm.o: srcterm.cpp
	$(CC) $(CFLAGS) srcterm.cpp

minmodlim.o: minmodlim.cpp
	$(CC) $(CFLAGS) minmodlim.cpp

conver.o: conver.cpp
	$(CC) $(CFLAGS) conver.cpp

wsolut.o: wsolut.cpp
	$(CC) $(CFLAGS) wsolut.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp


clean: 
	rm -rf *o nozzle
