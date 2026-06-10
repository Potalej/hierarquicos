FC = gfortran
FFLAGS = -O3 -march=native -fopenmp

# SRC = ../fortran/octree.f90 ../fortran/api3d.f90 parareal.f90
SRC = octree.f90 api3d.f90 parareal.f90

all:
	$(FC) $(FFLAGS) $(SRC) main.f90 -o main

clean:
	rm -rf test tree_time forces_time *.o *.mod