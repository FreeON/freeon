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
PROGRAM MakeM
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE MBlock
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                :: M,M2
#else
  TYPE(BCSR)                 :: M,M2
#endif
  INTEGER                    :: NCA,NCB
  REAL(DOUBLE), DIMENSION(3) :: B,A
  TYPE(AtomPair)             :: Pair
!
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(DBL_RNK4)             :: MD
!
  TYPE(ARGMT)                :: Args
  TYPE(DBL_VECT)             :: COrig
  INTEGER                    :: P,R,AtA,AtB,NN
  INTEGER                    :: IXYZ
  CHARACTER(LEN=*)              , PARAMETER :: Prog='MakeM'
  CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: Cart=(/'X','Y','Z'/)
!
!---------------------------------------
! Start up macro
!
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!
! Get Multipole origine. TODO
  CALL New(COrig,3)
  CALL SetEQ(COrig,0.0d0)
  !COrig%D(:)=GM%PBC%CellCenter%D(:)
!----------------------------------------------
! Allocations
!
  CALL New(MD,(/3,BS%NASym+2,BS%NASym+2,2*BS%NASym+4/),(/1,-1,-1,-1/))
  CALL New(M)
!-----------------------------------------------
! Run over cartisian componants
!
  DO IXYZ=1,3
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
     P=1; R=1; M%RowPt%I(1)=1
     CALL SetEq(M%MTrix,Zero)
!-----------------------------------------------
! Main loops
!
#ifdef PARALLEL
     M%NAtms=0
     DO AtA=Beg%I(MyId),End%I(MyId)
        M%NAtms=M%NAtms+1
#else
     M%NAtms=NAtoms
     DO AtA=1,NAtoms
#endif
        DO AtB=1,NAtoms
           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
              NN = Pair%NA*Pair%NB
!--------------------------------------------------------
              A(:) = Pair%A(:)
              B(:) = Pair%B(:)
              DO NCA = 1,CS_OUT%NCells
                 Pair%A(:) = A(:)+CS_OUT%CellCarts%D(:,NCA)
                 Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                           + (Pair%A(2)-Pair%B(2))**2 &
                           + (Pair%A(3)-Pair%B(3))**2
!--------------------------------------------------------
                 IF(TestAtomPair(Pair)) THEN
                    M%MTrix%D(R:R+NN-1)=M%MTrix%D(R:R+NN-1)+DBlok(BS,MD,Pair,IXYZ,COrig)
                 ENDIF
              ENDDO
              M%ColPt%I(P)=AtB
              M%BlkPt%I(P)=R

              R=R+NN
              P=P+1
#ifdef PARALLEL
              M%RowPt%I(M%NAtms+1)=P
              IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                 WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                 WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                 CALL Halt(' DBCSR dimensions blown in MakeM ')
              ENDIF
#else
              M%RowPt%I(AtA+1)=P
              IF(R>MaxNon0.OR.P>MaxBlks) &
                   CALL Halt(' BCSR dimensions blown in MakeM ')
#endif
           ENDIF
        ENDDO
     ENDDO
     M%NBlks=P-1
     M%NNon0=R-1
!------------------------------------------------------------
! Put D to disk
!
     CALL Filter(M2,M)
     CALL Put(M2,TrixFile('Dipole'//Cart(IXYZ),Args))
!------------------------------------------------------------
! Printing
!
     CALL PChkSum(M2,'Dipole'//Cart(IXYZ),Prog)
     CALL PPrint( M2,'Dipole'//Cart(IXYZ))
     CALL Plot(   M2,'Dipole'//Cart(IXYZ))
!     IF(Cart(IXYZ)=='X') CALL Print_BCSR(M2,'Dipole X',Unit_O=6)
!     IF(Cart(IXYZ)=='Y') CALL Print_BCSR(M2,'Dipole Y',Unit_O=6)
!     IF(Cart(IXYZ)=='Z') CALL Print_BCSR(M2,'Dipole Z',Unit_O=6)
!------------------------------------------------------------
  ENDDO ! End Loop over Cartesian Componants
!------------------------------------------------------------
! Tidy up
!

  CALL Delete(M )
  CALL Delete(M2)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL Delete(COrig)
!
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero
!
  CALL ShutDown(Prog)
!
END PROGRAM MakeM
