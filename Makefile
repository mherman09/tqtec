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
        readtqtec #\
        #fission_track_module \
        #ftage


# Programs to compile
all: tqtec \
     readtqtec #\
     #fission_track_module \
    #  #ftage


# Compilation rules
# Executables
tqtec: $(BIN)/tqtec
$(BIN)/tqtec: ./obj src/tqtec.f90 ./obj/error_exit.o
	$(FC) src/tqtec.f90 ./obj/error_exit.o   -o $(BIN)/tqtec   $(FOPT)

readtqtec: $(BIN)/readtqtec
$(BIN)/readtqtec: ./obj src/readtqtec.f90
	$(FC) src/readtqtec.f90   -o $(BIN)/readtqtec   $(FOPT)

# # Object/module files
# fission_track_module: ./obj/fission_track_module.o
# ./obj/fission_track_module.o: ./obj src/fission_track_module.f90
# 	$(FC) src/fission_track_module.f90    -c -o ./obj/fission_track_module.o   $(FOPT)

./obj/error_exit.o: ./obj src/error_exit.c
	$(CC) -c src/error_exit.c -o obj/error_exit.o

./obj:
	test -d ./obj || mkdir ./obj







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
