!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE ONXParameters
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
   REAL(DOUBLE), PARAMETER :: re0      = 0.0D0
   REAL(DOUBLE), PARAMETER :: re1      = 1.0D0
   REAL(DOUBLE), PARAMETER :: re2      = 2.0D0
   REAL(DOUBLE), PARAMETER :: rpi      = 3.1415926535898D0       ! Pi
   REAL(DOUBLE), PARAMETER :: PSch     = 4.9738746919472D0
   REAL(DOUBLE), PARAMETER :: Prev     = 34.986836655254D0       ! 2 Pi^(5/2)
   REAL(DOUBLE), PARAMETER :: Prev1    = 5.914967172796D0        ! Sqrt[2 Pi^(5/2)]
   REAL(DOUBLE), PARAMETER :: Prev2    = Prev1*Prev1
   REAL(DOUBLE), PARAMETER :: Small    = 1.0D-10
   INTEGER, PARAMETER      :: IDmn(0:4)= (/1,4,17,24,36/)

   REAL(DOUBLE), SAVE      :: DisRange = 0.0D0
   REAL(DOUBLE), SAVE      :: DenRange = 0.0D0
   REAL(DOUBLE), SAVE      :: ONXRange = 0.0D0
   REAL(DOUBLE), SAVE      :: MatRange = 0.0D0

   REAL(DOUBLE), SAVE      :: xNERIs   = 0.0D0

   INTEGER, SAVE           :: NRows=0,NCols=0,NElem=0
   INTEGER, SAVE           :: MaxN2=0

   LOGICAL, SAVE           :: Gradient=.FALSE.

   LOGICAL, SAVE           :: NoSym=.FALSE.

END MODULE ONXParameters
