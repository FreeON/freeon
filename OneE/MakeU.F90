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
!    COMPUTE THE EFFECTIVE CORE POTENTIAL MATRIX U
!    Author: Matt Challacombe
!------------------------------------------------------------
PROGRAM MakeU
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
  USE ECPBlock
  USE BraBloks
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                :: U,T1
#else
  TYPE(BCSR)                 :: U,T1
#endif
  INTEGER                    :: NCB,NCC
  REAL(DOUBLE),DIMENSION(3)  :: B
  TYPE(AtomPair)             :: Pair
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(DBL_RNK4)             :: MD
  TYPE(ARGMT)                :: Args
  REAL(DOUBLE)               :: Cx,Cy,Cz
  INTEGER                    :: P,R,AtA,AtB,AtC,KC,NN
  CHARACTER(LEN=5),PARAMETER :: Prog='MakeU'
  !---------------------------------------
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  !----------------------------------------------
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  !----------------------------------------------
  ! Allocations
  CALL NewBraBlok(BS)
  CALL New(U)
  !-----------------------------------------------
  ! Initialize the matrix and associated indecies
  P=1; R=1; U%RowPt%I(1)=1
  CALL SetEq(U%MTrix,Zero)
  !-----------------------------------------------
  ! Main loops
#ifdef PARALLEL
  U%NAtms=0
  DO AtA=Beg%I(MyId),End%I(MyId)
     U%NAtms=U%NAtms+1
#else
     U%NAtms=NAtoms
     DO AtA=1,NAtoms
#endif
        DO AtB=1,NAtoms
           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
              NN=Pair%NA*Pair%NB
              B=Pair%B
              DO NCB=1,CS_OUT%NCells
                 Pair%B=B+CS_OUT%CellCarts%D(:,NCB)
                 Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                         +(Pair%A(2)-Pair%B(2))**2 &
                         +(Pair%A(3)-Pair%B(3))**2
                 IF(TestAtomPair(Pair)) THEN
                    ! Go over ECP centers
                    DO AtC=1,NAtoms
                       KC=GM%AtTyp%I(AtC)
                       IF(BS%NTyp1PF%I(KC)>0)THEN
                   DO NCC=1,CS_OUT%NCells
                      Cx=GM%Carts%D(1,AtC)+CS_OUT%CellCarts%D(1,NCC)
                      Cy=GM%Carts%D(2,AtC)+CS_OUT%CellCarts%D(2,NCC)
                            Cz=GM%Carts%D(3,AtC)+CS_OUT%CellCarts%D(3,NCC)
                            U%MTrix%D(R:R+NN-1)=U%MTrix%D(R:R+NN-1)+UBlock(BS,Pair,KC,Cx,Cy,Cz)
                         ENDDO
                       ENDIF
                    ENDDO
                 ENDIF
              ENDDO
              U%ColPt%I(P)=AtB
              U%BlkPt%I(P)=R
              R=R+NN
              P=P+1
#ifdef PARALLEL
              U%RowPt%I(U%NAtms+1)=P
              IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                 WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                 WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                 CALL Halt(' DBCSR dimensions blown in MakeS ')
              ENDIF
#else
              U%RowPt%I(AtA+1)=P
              IF(R>MaxNon0.OR.P>MaxBlks) &
                   CALL Halt(' BCSR dimensions blown in MakeS ')
#endif
           ENDIF
        ENDDO
     ENDDO
     U%NBlks=P-1
     U%NNon0=R-1
     !------------------------------------------------------------
     ! Put U to disk
     CALL Filter(T1,U)
     CALL Put(T1,TrixFile('U',Args))
!     CALL PPrint( T1,'U',Unit_O=6)
     !-----------------------------------------------------------
     ! Printing
     CALL PChkSum(T1,'U',Prog)
     CALL PPrint( T1,'U')
     CALL Plot(   T1,'U')
     !------------------------------------------------------------
     ! Tidy up
     CALL Delete(U)
     CALL Delete(T1)
     CALL Delete(BS)
     CALL Delete(GM)
     CALL DeleteBraBlok()
     ! didn't count flops, any accumulation is residual
     ! from matrix routines
     PerfMon%FLOP=Zero
     CALL ShutDown(Prog)
   END PROGRAM MakeU
