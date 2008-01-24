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
MODULE ShellPairStruct
  USE DerivedTypes
  INTEGER,PARAMETER :: PairLngth=5
  INTEGER,PARAMETER :: MaxSPairs=100
  !
  TYPE SmallAtomInfo
     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
  END TYPE SmallAtomInfo
  !
  TYPE ShellPairG
     INTEGER :: IntType
     INTEGER :: L
     LOGICAL :: Switch
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(10,MaxSPairs) :: Cst
  END TYPE ShellPairG
  !
  TYPE AtomInfo
     INTEGER      :: K1,K2,NFPair
     REAL(DOUBLE) :: Atm1X,Atm1Y,Atm1Z
     REAL(DOUBLE) :: Atm2X,Atm2Y,Atm2Z
     REAL(DOUBLE) :: R12
     REAL(DOUBLE) :: Atm12X,Atm12Y,Atm12Z
  END TYPE AtomInfo

  TYPE ShellPair
     INTEGER :: IntType
     INTEGER :: L
     LOGICAL :: Switch
     TYPE(SmallAtomInfo) :: AtmInfo
     REAL(DOUBLE), DIMENSION(10,100) :: Cst
  END TYPE ShellPair

  TYPE AtomPr
     TYPE(ShellPair) :: SP
  END TYPE AtomPr

  TYPE AtomPrG
     TYPE(ShellPairG) :: SP
  END TYPE AtomPrG

  TYPE ONX2OffSt
     INTEGER :: A,B,C,D
  END TYPE ONX2OffSt
  !
END MODULE ShellPairStruct
