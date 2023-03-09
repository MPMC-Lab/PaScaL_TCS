HOME = ..
include $(HOME)/Makefile.inc

EXE = PaScaL_TCS.ex 

OBJDIR = ./obj

global      = $(OBJDIR)/module_global.o
module_mpi  = $(OBJDIR)/module_mpi_topology.o $(OBJDIR)/module_mpi_subdomain.o
solver	    = $(OBJDIR)/module_solve_thermal.o $(OBJDIR)/module_solve_momentum.o $(OBJDIR)/module_solve_pressure.o 
analysis    = $(OBJDIR)/module_post.o $(OBJDIR)/module_mpi_statistics.o
main 	    = $(OBJDIR)/main.o 

OBJS   	    = $(global) $(module_mpi) $(solver) $(analysis) $(main)

MKL_HOME    = /opt/intel/oneapi/mkl/latest
PaScaL_HOME = ../PaScaL_TDMA

PaScaLINC = -I$(PaScaL_HOME)/include
PaScaLLIB = -L$(PaScaL_HOME)/lib -lpascal_tdma

#FFTWINC = -I$(MKL_HOME)/include/fftw
FFTWINC = -I/scratch/x2434a03/fftw3/fftw3-library/include

FFTWFLAGS 	=  -qmkl=sequential

$(OBJDIR)/%.o: %.f90
	$(FC) $(FLAG) -c $^ -o $@ -module $(OBJDIR) -I$(OBJDIR) $(FFTWINC) $(FFTWLIB) ${PaScaLINC} $(FFTWFLAGS)

all:$(OBJS)
	$(FC) $(FLAG) $(FFTWINC) $(FFTWLIB) $(OBJS) -o $(EXE) $(PaScaLINC) $(PaScaLLIB) $(FFTWFLAGS)
	mv $(EXE) ../run

clean:
	rm -rf $(OBJDIR); rm -f ../run/$(EXE); rm -rf ../run/data 
