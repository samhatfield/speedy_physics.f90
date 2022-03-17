PYFC = f2py
FC = gfortran

physics: physical_constants humidity convection
	$(PYFC) -c -I. humidity.o convection.o -m physics physics.f90

physical_constants:
	$(FC) -c physical_constants.f90

humidity:
	$(FC) -c -fPIC humidity.f90

convection:
	$(FC) -c -fPIC convection.f90

.PHONY: clean
clean:
	rm -f *.o *.mod *.so *.pyf

