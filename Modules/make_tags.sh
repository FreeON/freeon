#!/bin/bash

ctags --Fortran-kinds=+i+L --exclude="*-preprocessed*" *.{c,f,f90}
