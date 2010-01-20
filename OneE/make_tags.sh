#!/bin/bash

ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.f90 ../Modules/*.{c,f,f90}
