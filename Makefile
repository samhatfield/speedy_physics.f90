PYFC = f2py
FC = gfortran

all: physics params

physics: params.f90
	$(PYFC) --overwrite-signature -h physics.pyf $^ physics.f90 -m $@
	$(PYFC) -c physics.pyf $^ physics.f90
	mv $@.*.so physics.so

params:
	$(PYFC) -c -m params params.f90
	mv $@.*.so params.so

.PHONY: clean
clean:
	rm -f *.o *.mod *.so *.pyf

