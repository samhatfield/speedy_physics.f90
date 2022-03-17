PYFC = f2py
FC = gfortran

physics: physical_constants humidity
	$(PYFC) -c -I. humidity.o -m physics physics.f90

humidity:
	$(FC) -c -fPIC humidity.f90

physical_constants:
	$(FC) -c physical_constants.f90

.PHONY: clean
clean:
	rm -f *.o *.mod *.so *.pyf

