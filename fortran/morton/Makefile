FC = gfortran
FFLAGS = -O3 -march=native -fopenmp

SRC = morton_octree.f90 ../octree.f90 test_morton.f90

all:
	$(FC) $(FFLAGS) $(SRC) -o test_morton

clean:
	rm -rf test tree_time forces_time *.o *.mod