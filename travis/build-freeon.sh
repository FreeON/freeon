#!/bin/bash -x

./autogen.sh
./configure F77=gfortran FC=gfortran CC=gfortran LIBS="-l:liblapack.so.3gf -l:libblas.so.3gf" --disable-openmp
make -j3 -C Modules
make -j3 -C $COMPONENT
#make test
