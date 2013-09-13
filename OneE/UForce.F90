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
!    COMPUTE THE FORCE CORESPONDING TO THE DERIVATIVE OF THE
!    KINETIC ENERGY MATRIX, TForce=2*Tr{P.dT}
!    Authors: Matt Challacombe
!----------------------------------------------------------------------------------
PROGRAM UForce
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
#ifdef PARALLEL
  USE MondoMPI
  TYPE(BCSR)                  :: P
  TYPE(DBL_VECT)              :: TotUFrc
  INTEGER                     :: IErr,TotFrcComp
#else
  TYPE(BCSR)                  :: P
#endif
  INTEGER                     :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)   :: B,F
  TYPE(AtomPair)              :: Pair
  TYPE(BSET)                  :: BS
  TYPE(CRDS)                  :: GM
  TYPE(ARGMT)                 :: Args
  INTEGER                     :: Q,R,AtA,AtB,AtC,KC,JP,MB,MA,NB,MN1,MN,A1,A2,C1,C2,NCB,NCC
  TYPE(HGRho)                 :: Rho
  TYPE(DBL_VECT)              :: UFrc,Frc
  REAL(DOUBLE)                :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,UFrcChk,F1a,F2a,F3a,f11,f10

  CHARACTER(LEN=6),PARAMETER  :: Prog='UForce'
!-------------------------------------------------------------------------------------
! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
  CALL New(UFrc,3*NAtoms)
  UFrc%D=Zero
  ! First compute
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     A1=3*(AtA-1)+1
     A2=3*AtA
     MA=BSiz%I(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
           MN=MN1+1
           !Quick and dirty!
           SELECT CASE(P%NSMat)
           CASE(1)
              !We don't need to do anything!
           CASE(2)
              !We add up the two density martices!
              CALL DAXPY(MN,1D0,P%MTrix%D(Q+MN),1,P%MTrix%D(Q),1)
           CASE(4)
              !We add up the diagonal density martices!
              CALL DAXPY(MN,1D0,P%MTrix%D(Q+3*MN),1,P%MTrix%D(Q),1)
           CASE DEFAULT;CALL Halt(' UForce: P%NSMat doesn''t have an expected value! ')
           END SELECT
           B=Pair%B
           DO NCB=1,CS_OUT%NCells
              Pair%B=B +CS_OUT%CellCarts%D(:,NCB)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 ! Go over ECP centers
                 DO AtC=1,NAtoms
                    C1=3*(AtC-1)+1
                    C2=3*AtC
                    KC=GM%AtTyp%I(AtC)
                    IF(BS%NTyp1PF%I(KC)>0)THEN
                       DO NCC=1,CS_OUT%NCells
                          Cx=GM%Carts%D(1,AtC)+CS_OUT%CellCarts%D(1,NCC)
                          Cy=GM%Carts%D(2,AtC)+CS_OUT%CellCarts%D(2,NCC)
                          Cz=GM%Carts%D(3,AtC)+CS_OUT%CellCarts%D(3,NCC)
                          F=Two*TrPdU(BS,Pair,KC,Cx,Cy,Cz,P%MTrix%D(Q:Q+MN1))
                          ! Derivative contribution from bra <A| and ket |B>
                          UFrc%D(A1:A2)=UFrc%D(A1:A2)+F
                          ! then add in the ECP center |C| part
                          UFrc%D(C1:C2)=UFrc%D(C1:C2)-F
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO


!!$!=====================================================================================
!!$     f10=Zero
!!$     DO AtA=1,NAtoms
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     MA=BSiz%I(AtA)
!!$     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
!!$        AtB=P%ColPt%I(JP)
!!$           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
!!$              Q=P%BlkPt%I(JP)
!!$              B=Pair%B
!!$          NB=BSiz%I(AtB)
!!$           MN1=MA*NB-1
!!$              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
!!$                   +(Pair%A(2)-Pair%B(2))**2 &
!!$                   +(Pair%A(3)-Pair%B(3))**2
!!$              IF(TestAtomPair(Pair)) THEN
!!$                 ! Go over ECP centers
!!$                 DO AtC=1,NAtoms
!!$                    KC=GM%AtTyp%I(AtC)
!!$                    IF(BS%NTyp1PF%I(KC)>0)THEN
!!$                       Cx=GM%Carts%D(1,AtC)
!!$                       Cy=GM%Carts%D(2,AtC)
!!$                       Cz=GM%Carts%D(3,AtC)
!!$                       Ax = Pair%A(1)
!!$                       Ay = Pair%A(2)
!!$                       Az = Pair%A(3)
!!$                       Bx = Pair%B(1)
!!$                       By = Pair%B(2)
!!$                       Bz = Pair%B(3)
!!$                       f10=f10+DUBlock(BS,Pair,Ax,Ay,Az,Bx,By,Bz,KC,Cx,Cy,Cz,P%MTrix%D(Q:Q+MN1))
!!$                    ENDIF
!!$                 ENDDO
!!$              ENDIF
!!$           ENDIF
!!$        ENDDO
!!$     ENDDO
!!$!=====================================================================================
!!$     GM%Carts%D(1,1)=GM%Carts%D(1,1)+1D-5
!!$     f11=Zero
!!$     DO AtA=1,NAtoms
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     MA=BSiz%I(AtA)
!!$     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
!!$        AtB=P%ColPt%I(JP)
!!$           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
!!$           Q=P%BlkPt%I(JP)
!!$          NB=BSiz%I(AtB)
!!$           MN1=MA*NB-1
!!$              B=Pair%B
!!$              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
!!$                   +(Pair%A(2)-Pair%B(2))**2 &
!!$                   +(Pair%A(3)-Pair%B(3))**2
!!$              IF(TestAtomPair(Pair)) THEN
!!$                 ! Go over ECP centers
!!$                 DO AtC=1,NAtoms
!!$                    KC=GM%AtTyp%I(AtC)
!!$                    IF(BS%NTyp1PF%I(KC)>0)THEN
!!$                       Cx=GM%Carts%D(1,AtC)
!!$                       Cy=GM%Carts%D(2,AtC)
!!$                       Cz=GM%Carts%D(3,AtC)
!!$                       Ax = Pair%A(1)
!!$                       Ay = Pair%A(2)
!!$                       Az = Pair%A(3)
!!$                       Bx = Pair%B(1)
!!$                       By = Pair%B(2)
!!$                       Bz = Pair%B(3)
!!$                       f11=f11+DUBlock(BS,Pair,Ax,Ay,Az,Bx,By,Bz,KC,Cx,Cy,Cz,P%MTrix%D(Q:Q+MN1))
!!$                    ENDIF
!!$                 ENDDO
!!$              ENDIF
!!$           ENDIF
!!$        ENDDO
!!$     ENDDO
UFrc%D=Two*UFrc%D
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotUFrc,TotFrcComp)
  CALL MPI_Reduce(UFrc%D(1),TotUFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    UFrc%D(1:TotFrcComp)=TotUFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotUFrc)
#endif
! Rescale the Forces if needed.
  IF(P%NSMat.GT.1) CALL DSCAL(3*NAtoms,0.5D0,    UFrc%D(1  ),1)
  !Dont forget to rescale the lattice forces for unrestricted theory!
  !IF(P%NSMat.GT.1) CALL DSCAL(       9,0.5D0,LatFrc_U%D(1,1),1)

!  CALL PPrint(UFrc,'dU',Unit_O=6)
! Do some checksumming, resumming and IO
  CALL PChkSum(UFrc,'dU/dR',Proc_O=Prog)
! Sum in contribution to total force
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
!     WRITE(*,11)AtA,UFrc%D(A1:A2)
     11 format(I3," ",3(F10.5," "))
     GM%Gradients%D(1:3,AtA)=GM%Gradients%D(1:3,AtA)+UFrc%D(A1:A2)
  ENDDO
  CALL Put(GM,Tag_O=CurGeom)
!------------------------------------------------------------------------------
! Tidy up
  CALL Delete(UFrc)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL ShutDown(Prog)
END PROGRAM UForce
