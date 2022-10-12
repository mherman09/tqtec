#####
#	TQTEC MAKEFILE
#	1-D THERMAL MODELING AND THERMOCHRONOLOGY
#####

# Fortran compiler
FC = gfortran

# Fortran compiler options
FOPT = -Wall -O3 -J ./obj

# C compiler
CC = gcc

# Directory to install compiled programs
BIN = ./bin



##### YOU SHOULD NOT NEED TO EDIT BELOW HERE #####

.PHONY: clean \
        tqtec \
        readtqtec


# Programs to compile
all: tqtec \
     readtqtec


# Compilation rules
tqtec: $(BIN)/tqtec
$(BIN)/tqtec: src/tqtec.f90 obj/error_exit.o
	test -d ./obj || mkdir ./obj
	$(FC) $^ -o $@ $(FOPT)

readtqtec: $(BIN)/readtqtec
$(BIN)/readtqtec: src/readtqtec.f90
	test -d ./obj || mkdir ./obj
	$(FC) $< -o $@ $(FOPT)

obj/error_exit.o: src/error_exit.c
	test -d ./obj || mkdir ./obj
	$(CC) -c src/error_exit.c -o obj/error_exit.o


# Old executables, no longer compiled by default
tqtec_old: $(BIN)/tqtec_old
$(BIN)/tqtec_old: src/TQTec.f
	test -d ./obj || mkdir ./obj
	$(FC) $< -o $@ $(FOPT)

readtqtec_old: $(BIN)/readtqtec_old
$(BIN)/readtqtec_old: src/readTQTec02.f
	test -d ./obj || mkdir ./obj
	$(FC) $< -o $@ $(FOPT)


# Clean up
clean:
	-rm $(BIN)/*
	-rm -rf ./obj
	-rm *.mod
