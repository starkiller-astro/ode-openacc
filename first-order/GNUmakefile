#FC := pgf95
FC := gfortran
#FC := Cray

ALL: test_react

SRC_DIRS := . ../lapack/ ../blas/

include ../GMake.common


# default target -- this is the executable
test_react: $(OBJECTS)
	@echo " "
	@echo "Linking..."
	$(wrapper) -o test_react $(OBJECTS)

