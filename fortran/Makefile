FC = gfortran
FFLAGS = -O3 -march=native -fopenmp
# FFLAGS = -O0 -march=native -fopenmp -fsanitize=address

# SRC = octree.f90 morton.f90
# SRC = node_pointer.f90 octree_pointer.f90
SRC = octree.f90 morton/morton_octree.f90

TARGET ?= test

all:
	$(FC) $(FFLAGS) $(SRC) $(TARGET).f90 -o $(TARGET)

clean:
	rm -rf test tree_time forces_time *.o *.mod