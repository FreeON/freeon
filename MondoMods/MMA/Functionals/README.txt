The package Format.m extends Mathematica's formatting rules.
Assignments to expressions and formatting of lists is now possible.
Code translation includes C, FORTRAN77, TeX and a limited subset
of Maple.  Optimized computational sequences are possible utilising
the auxiliary package Optimize.m, item number 0206-592 (for more
details, refer to the documentation for that package).

LaTeX documentation for the package Format.m is available in the
file Format.tex.  Typeset the file twice to generate the table of
contents and the reference list.  This file uses the LaTeX style
file format.sty. The page format will need to be adjusted for
printing in the US (appropriate settings are given in the file
commented out using the character %).

Various files with .f and .mf extensions provide examples of the
use of Splice to process template files. These contain generalised
formulations to problems containing a mixture of FORTRAN77 and
Mathematica code. The files are:

driver.f - a main driver routine for the files func.f and sub.f
(resulting from the application of Splice to func.mf and sub.mf).

example.mf - a simple FORTRAN77 main routine template file.

func.mf - a FORTRAN77 function routine template file.

linsolv.f - the LAPACK LU decomposition, linear solution and related
routines for use with the file newton.mf.

newton.mf - a multi-dimensional FORTRAN77 template file for newton's
method.  This file uses the lapack linear solution routine linsolv.f.

rksub.mf - a FORTRAN77 function routine template file for the
Runge-Kutta method for a Hamiltonian system.

sub.mf - a FORTRAN77 subroutine template file.

