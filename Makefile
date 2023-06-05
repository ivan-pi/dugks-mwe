


FC = gfortran
FCFLAGS = -Wall -fopenmp
LDFLAGS =
LDLIBS =

all: main_taylor_green

OBJECTS = precision.o gnuplot_io.o lattice.o periodic_dugks.o taylor_green.o

main_taylor_green: main_taylor_green.f90 $(OBJECTS)
	$(FC) $(LDFLAGS) $(FCFLAGS) -o $@ $^ $(LDLIBS)

#
# Modules
#
precision.o precision.mod: precision.F90
	$(FC) $(FCFLAGS) -c $<

gnuplot_io.o gnuplot_io.mod: gnuplot_io.F90
	$(FC) $(FCFLAGS) -c $<

lattice.o lattice.mod: lattice.F90
	$(FC) $(FCFLAGS) -c $<

periodic_dugks.o periodic_dugks.mod: periodic_dugks.F90
	$(FC) $(FCFLAGS) -c $<

taylor_green.o taylor_green.mod: taylor_green.f90
	$(FC) $(FCFLAGS) -c $<

#
# Module dependencies
#
gnuplot_io.mod: precision.mod
lattice.mod: precision.mod gnuplot_io.mod
periodic_dugks.mod: precision.mod lattice.mod
taylor_green.mod: precision.mod lattice.mod

.phony: clean

clean:
	rm -rf *.o *.mod main_taylor_green



