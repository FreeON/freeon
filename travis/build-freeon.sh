#!/bin/bash

./autogen.sh
./configure F77=gfortran FC=gfortran CC=gfortran LIBS="-l:liblapack.so.3gf -l:libblas.so.3gf" --disable-openmp
make -C $COMPONENT
#make test
