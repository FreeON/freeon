#!/bin/bash

ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{f,f90} ../Modules/*.{c,f,f90}
