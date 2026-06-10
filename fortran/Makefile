FC = gfortran
FFLAGS = -O3 -march=native -fopenmp

SRC = octree.f90 morton.f90
# SRC = node_pointer.f90 octree_pointer.f90

TARGET ?= test

all:
	$(FC) $(FFLAGS) $(SRC) $(TARGET).f90 -o $(TARGET)

clean:
	rm -rf test tree_time forces_time *.o *.mod