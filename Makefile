


FC = gfortran
FCFLAGS = -Wall -fopenmp

all: main_taylor_green

OBJECTS = precision.o output_gnuplot.o lattice.o periodic_dugks.o taylor_green.o

main_taylor_green: main_taylor_green.f90 $(OBJECTS)
	$(FC) $(FCFLAGS) -o $@ $^

#
# Modules
#
precision.o precision.mod: precision.F90
	$(FC) $(FCFLAGS) -c $<

output_gnuplot.o output_gnuplot.mod: output_gnuplot.F90
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
output_gnuplot.mod: precision.mod
lattice.mod: precision.mod output_gnuplot.mod
periodic_dugks.mod: precision.mod lattice.mod
taylor_green.mod: precision.mod lattice.mod

.phony: clean

clean:
	rm -rf *.o *.mod main_taylor_green



