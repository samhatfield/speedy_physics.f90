PYFC = f2py
FC = gfortran

physics: utils physical_constants humidity convection condensation radiation
	$(PYFC) -c -I. utils.o humidity.o convection.o condensation.o radiation.o -m physics physics.f90

physical_constants:
	$(FC) -c physical_constants.f90

humidity:
	$(FC) -c -fPIC humidity.f90

convection: physical_constants utils
	$(FC) -c -fPIC convection.f90

condensation: physical_constants utils
	$(FC) -c -fPIC condensation.f90

radiation:
	$(FC) -c -fPIC radiation.f90

utils: physical_constants
	$(FC) -c -fPIC utils.f90

.PHONY: clean
clean:
	rm -f *.o *.mod *.so *.pyf

