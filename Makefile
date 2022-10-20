#####
#	TQTEC MAKEFILE
#	1-D THERMAL MODELING AND THERMOCHRONOLOGY
#####

# Fortran compiler
FC = gfortran

# Fortran compiler options
FOPT = -Wall -Wextra -pedantic -fcheck=all -J ./obj -g -ffpe-trap=invalid,zero,overflow -fbacktrace
# FOPT = -Wall -Wextra -J ./obj

# C compiler
CC = gcc

# Directory to install compiled programs
BIN = ./bin



##### YOU SHOULD NOT NEED TO EDIT BELOW HERE #####

.PHONY: clean \
        tqtec \
        readtqtec \
        fission_track_module \
        minage


# Programs to compile
all: tqtec \
     readtqtec \
     minage


# Compilation rules
# Executables
tqtec: $(BIN)/tqtec
$(BIN)/tqtec: src/tqtec.f90 ./obj/error_exit.o
	$(FC) src/tqtec.f90 ./obj/error_exit.o   -o $(BIN)/tqtec   $(FOPT)

readtqtec: $(BIN)/readtqtec
$(BIN)/readtqtec: src/readtqtec.f90
	$(FC) src/readtqtec.f90   -o $(BIN)/readtqtec   $(FOPT)

minage: $(BIN)/minage
$(BIN)/minage: src/minage.f90    ./obj/fission_track_module.o ./obj/error_exit.o
	$(FC) src/minage.f90 ./obj/fission_track_module.o ./obj/error_exit.o   -o $(BIN)/minage   $(FOPT)

# Object/module files
fission_track_module: ./obj/fission_track_module.o
./obj/fission_track_module.o: src/fission_track_module.f90
	$(FC) src/fission_track_module.f90    -c -o ./obj/fission_track_module.o   $(FOPT)

./obj/error_exit.o: src/error_exit.c
	$(CC) -c src/error_exit.c -o obj/error_exit.o








# Old executables, no longer compiled by default
tqtec_old: $(BIN)/tqtec_old
$(BIN)/tqtec_old: src/TQTec.f
	$(FC) $< -o $@ $(FOPT)

readtqtec_old: $(BIN)/readtqtec_old
$(BIN)/readtqtec_old: src/readTQTec02.f
	$(FC) $< -o $@ $(FOPT)

ftage_old: $(BIN)/ftage_old
$(BIN)/ftage_old: src/ftage.f90
	$(FC) $< -o $@ $(FOPT)


# Clean up
clean:
	-rm $(BIN)/*
	-rm obj/*
	-rm *.mod
