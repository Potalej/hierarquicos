FC = gfortran
FFLAGS = -O3 -march=native -fopenmp
# FFLAGS = -O0 -march=native -fopenmp -fsanitize=address

ROOT_DIR = ../..
FORTRAN_DIR = $(ROOT_DIR)/fortran

# libs
FILES_IV := $(shell cat $(ROOT_DIR)/initial_values/module.txt)
SOURCES_IV := $(addprefix $(ROOT_DIR)/initial_values/, $(FILES_IV))
FORCES_LIB = $(FORTRAN_DIR)/direct_forces.f90
SRC = ../test_forces.f90 ${FORTRAN_DIR}/octree_default.f90 ${FORCES_LIB} $(SOURCES_IV)

all:
	mkdir -p build

time_and_error:
	$(FC) -Jbuild $(FFLAGS) $(SRC) forces_time_and_error.f90 -o forces_time_and_error

only_time:
	$(FC) -Jbuild $(FFLAGS) $(SRC) forces_only_time.f90 -o forces_only_time

clean:
	rm -rf forces_time_and_error *.o *.mod