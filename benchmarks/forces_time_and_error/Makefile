FC = gfortran
FFLAGS = -O3 -march=native -fopenmp
# FFLAGS = -O0 -march=native -fopenmp -fsanitize=address

ROOT_DIR = ../..
FORTRAN_DIR = $(ROOT_DIR)/fortran

# libs
IV_LIB = $(ROOT_DIR)/initial_values/initial_values.f90
DF_LIB = $(FORTRAN_DIR)/direct_forces.f90

SRC = ${FORTRAN_DIR}/octree_default.f90 ${DF_LIB} $(IV_LIB)

TARGET = forces_time_and_error

all:
	mkdir -p build
	$(FC) -Jbuild $(FFLAGS) $(SRC) $(TARGET).f90 -o $(TARGET)

clean:
	rm -rf forces_time_and_error *.o *.mod