# compile the module with f2py
python -m numpy.f2py -c -m api octree.f90 api3d.f90 --f90flags="-fopenmp" -lgomp