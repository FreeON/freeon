#!/bin/bash

ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{c,f90} ../Modules/*.{c,f,f90}
