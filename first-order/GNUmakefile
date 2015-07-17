#FC := pgf95
FC := pgf95_local
#FC := gfortran
#FC := Cray

ALL: test_react

SRC_DIRS := . ../lapack/ ../blas/

include ../GMake.common


# default target -- this is the executable
test_react: $(OBJECTS)
	@echo " "
	@echo "Linking..."
	$(link) -o test_react $(OBJECTS)

