PYFC = f2py
FC = gfortran

all: physics

physics:
	$(PYFC) --overwrite-signature -h physics.pyf $^ physics.f90 -m $@
	$(PYFC) -c physics.pyf $^ physics.f90
	mv $@.*.so physics.so

humidity:
	$(PYFC) -c $^ humidity.f90 -m $@
	mv $@.*.so humidity.so


.PHONY: clean
clean:
	rm -f *.o *.mod *.so *.pyf

