all: bin/tqtec \
     bin/tqtec_2021 \
     bin/readtqtec \
     bin/readtqtec_2021

FC = gfortran
FOPT = -Wall

bin/tqtec: src/TQTec.f
	$(FC) $< -o $@ $(FOPT)

bin/tqtec_2021: src/tqtec.f90
	$(FC) $< -o $@ $(FOPT)

bin/readtqtec: src/readTQTec02.f
	$(FC) $< -o $@ $(FOPT)

bin/readtqtec_2021: src/readtqtec.f90
	$(FC) $< -o $@ $(FOPT)
