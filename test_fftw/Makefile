# Compiler settings
FC=gfortran
FFLAGS=-Wall -Wextra -pedantic -Wno-unused-parameter -O2
CXX=g++
CXXFLAGS=-std=c++17 -Wall -Wextra -pedantic -Og

# Default directory if not provided on command line
FFTW3_DIR ?= /usr/local

INC      := -I$(FFTW3_DIR)/include
LDFLAGS  := -L$(FFTW3_DIR)/lib -Wl,-rpath,$(FFTW3_DIR)/lib
LDLIBS   := -lfftw3

# Targets
TARGETS=test_fftw test_suite

.PHONY: all
all: $(TARGETS)

test_fftw: test_fftw.f90
	$(FC) $(FFLAGS) $(INC) $(LDFLAGS) -o $@ $< $(LDLIBS)

test_suite: test_suite.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(LDFLAGS) -o $@ $< $(LDLIBS)

run: all
	./test_fftw
	./test_suite

# Package Metadata
PACKAGE_NAME := test_fftw
DIST_DIR     := $(PACKAGE_NAME)

# --- Distribution Rule ---
.PHONY: dist
dist: clean
	@echo "Creating distribution tarball..."
	@mkdir -p $(DIST_DIR)
	@# Copy only source files, headers, and the Makefile
	cp README.txt Makefile *.cpp *.f90 $(DIST_DIR)/
	@# Create the tarball and remove the temporary directory
	tar -czvf $(DIST_DIR).tar.gz $(DIST_DIR)
	@rm -rf $(DIST_DIR)
	@echo "Distribution ready: $(DIST_DIR).tar.gz"

.PHONY: clean
clean:
	$(RM) $(TARGETS) *.o *.mod *.tar.gz
	$(RM) -rf $(DIST_DIR)