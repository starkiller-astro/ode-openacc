FC := gfortran
#FC := pgf95

ACC :=

include ../GMake.defs

programs += test_react

ALL: $(pnames)


SRC_DIRS := . ../lapack/ ../blas/ ../boxlib/

include ../GMake.common


# default target -- this is the executable
%.$(suf).exe: $(OBJECTS)
	@echo " "
	@echo "Linking..."
	$(link) -o $@ $(OBJECTS)

