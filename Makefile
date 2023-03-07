# Compiler and complilation flag are included in Makefile.inc
include ./Makefile.inc

all:
	cd PaScaL_TDMA; make lib
	cd src; mkdir -p obj; make all

lib:
	cd PaScaL_TDMA; make lib

exe:
#	cd src; mkdir -p obj; make all
	cd run; mpirun -np 8 ./PaScaL_TCS.ex
clean:
	cd PaScaL_TDMA; make clean
	cd src; make clean

execlean:
	cd src; make clean

rm:
	rm -r ./run/data