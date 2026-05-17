FC = gfortran
FFLAGS = -O3

SRC = node.f90 tree.f90

TARGET ?= test

all:
	$(FC) $(FFLAGS) $(SRC) $(TARGET).f90 -o $(TARGET)

clean:
	rm -rf test tree_time forces_time *.o *.mod