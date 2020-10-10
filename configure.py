#!/usr/bin/python
if __name__ == '__main__':
 im:port sys
 import os
 sys.path.insert(0, os.path.abspath('config'))
 import configure
 configure_options = [
   '--download-fblaslapack=1',
   '--known-bits-per-byte=8',
   '--known-level1-dcache-assoc=0',
   '--known-level1-dcache-linesize=32',
   '--known-level1-dcache-size=32768',
   '--known-memcmp-ok=1',
   '--known-mpi-c-double-complex=0',
   '--known-mpi-long-double=1',
   '--known-mpi-shared-libraries=0',
   '--known-mpi-shared=0',
  # '--download-mpi=1',
   '--with-mpi-dir=/home/<USERNAME>/<MPI_INSTALL_DIR>',
   '--known-sizeof-MPI_Comm=8',
   '--known-sizeof-MPI_Fint=4',
   '--known-sizeof-char=1',
   '--known-sizeof-double=8',
   '--known-sizeof-float=4',
   '--known-sizeof-int=4',
   '--known-sizeof-long-long=8',
   '--known-sizeof-long=8',
   '--known-sizeof-short=2',
   '--known-sizeof-size_t=8',
   '--known-sizeof-void-p=8',
   '--with-debugging=0',
   '--with-fortran-kernels=1',
   '--with-gcov=0',
   '--with-x=0',
   '--with-64-bit-indices=1',
   '--download-elemental=1',
   '--with-cxx-dialect=C++11',
   '--download-metis=1',
   '--download-parmetis=1'
#   'PETSC_ARCH=try5',
 ]
 configure.petsc_configure(configure_options)

