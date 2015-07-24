FC := pgf95
#FC := pgf95_local
#FC := gfortran
#FC := Cray

ACC :=

include ../GMake.defs

programs += test_react

ALL: $(pnames)



SRC_DIRS := . ../lapack/ ../blas/

include ../GMake.common


# default target -- this is the executable
%.$(suf).exe: $(OBJECTS)
	@echo " "
	@echo "Linking..."
	$(link) -o $@  $(OBJECTS)

