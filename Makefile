


FC = gfortran
FCFLAGS = -Wall -fopenmp
LDFLAGS = -L /opt/homebrew/Cellar/flang/21.1.8/lib
LDLIBS = -lflang_rt.runtime

all: main_taylor_green

OBJECTS = precision.o gnuplot_io.o lattice.o periodic_dugks.o periodic_lw.o taylor_green.o

TUNE= -D NCOLLAPSE=$(NCOLLAPSE) -D TILE_1=$(TILE_1) -D TILE_2=$(TILE_2)




main_taylor_green: main_taylor_green.F90 $(OBJECTS) flang_l2_norm.o
	$(FC) $(TUNE) $(LDFLAGS) $(FCFLAGS)  -o $@ $^ $(LDLIBS)

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

periodic_lw.o periodic_lw.mod: periodic_lw.F90
	$(FC) $(FCFLAGS) -c $<

taylor_green.o taylor_green.mod: taylor_green.f90
	$(FC) $(FCFLAGS) -c $<

flang_l2_norm.o: flang_l2_norm.f90
	flang -O3 -static-libflangrt -c $<

#
# Module dependencies
#
gnuplot_io.mod: precision.mod
lattice.mod: precision.mod gnuplot_io.mod
periodic_dugks.mod: precision.mod lattice.mod
periodic_lw.mod: precision.mod lattice.mod
taylor_green.mod: precision.mod lattice.mod

.phony: clean

clean:
	rm -rf *.o *.mod main_taylor_green



