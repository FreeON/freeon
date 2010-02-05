!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
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
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------
!    COMPUTE THE FORCE CORESPONDING TO THE DERIVATIVE OF THE 
!    KINETIC ENERGY MATRIX, TForce=2*Tr{P.dT}
!    Authors: Matt Challacombe and CJ Tymczak
!----------------------------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM TForce
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
  USE BlokTrPdT

#ifdef PARALLEL
  USE MondoMPI
  TYPE(BCSR)                  :: P
  TYPE(DBL_VECT)              :: TotTFrc
  TYPE(DBL_RNK2)              :: TmpLatFrc_T
  INTEGER                     :: IErr,TotFrcComp
#else
  TYPE(BCSR)                  :: P
#endif
  INTEGER                     :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)   :: A,B,F_nlm,nlm
  TYPE(DBL_RNK2)              :: LatFrc_T
  TYPE(AtomPair)              :: Pair
  TYPE(BSET)                  :: BS
  TYPE(CRDS)                  :: GM
  TYPE(ARGMT)                 :: Args
  INTEGER                     :: Q,R,AtA,AtB,JP,MB,MA,NB,MN1,MN,A1,A2
  TYPE(HGRho)                 :: Rho
  TYPE(DBL_VECT)              :: TFrc,Frc
  REAL(DOUBLE)                :: TFrcChk
 
  CHARACTER(LEN=6),PARAMETER  :: Prog='TForce'
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
!
  CALL New(TFrc,3*NAtoms)
  TFrc%D   = Zero
  CALL New(LatFrc_T,(/3,3/))
  LatFrc_T%D = Zero
!--------------------------------------------------------------------------------
! TForce=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
!--------------------------------------------------------------------------------
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
           CASE DEFAULT;CALL Halt(' TForce: P%NSMat doesn''t have an expected value! ')
           END SELECT
           A=Pair%A
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%A=A
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:3)    = TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
                 TFrc%D(A1:A2) = TFrc%D(A1:A2) + Two*F_nlm(1:3)
!                Lattice Forces
                 nlm        = AtomToFrac(GM,Pair%A)
                 LatFrc_T%D = LatFrc_T%D + Two*LaticeForce(GM,nlm,F_nlm)
              ENDIF
              Pair%A=A+CS_OUT%CellCarts%D(:,NC)
              Pair%B=B
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:3)    = TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
                 TFrc%D(A1:A2) = TFrc%D(A1:A2) + Two*F_nlm(1:3)
!                Lattice Forces
                 nlm        = AtomToFrac(GM,Pair%A)
                 LatFrc_T%D = LatFrc_T%D + Two*LaticeForce(GM,nlm,F_nlm)
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
#ifdef PARALLEL
! Collect the Forces
  TotFrcComp = 3*NAtoms
  CALL New(TotTFrc,TotFrcComp)
  CALL MPI_Reduce(TFrc%D(1),TotTFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
    TFrc%D(1:TotFrcComp) = TotTFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotTFrc)
! Collect the Lattice Forces
  CALL New(TmpLatFrc_T,(/3,3/))
  CALL DBL_VECT_EQ_DBL_SCLR(9,TmpLatFrc_T%D(1,1),0.0d0)
  CALL MPI_REDUCE(LatFrc_T%D(1,1),TmpLatFrc_T%D(1,1),9,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,ROOT,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     LatFrc_T%D = TmpLatFrc_T%D
  ENDIF
  CALL Delete(TmpLatFrc_T)
#endif
#ifdef PARALLEL  
  IF(MyID == ROOT) THEN
#endif
! Rescale the Forces if needed.
     IF(P%NSMat.GT.1) CALL DSCAL(3*NAtoms,0.5D0,    TFrc%D(1  ),1)
     IF(P%NSMat.GT.1) CALL DSCAL(       9,0.5D0,LatFrc_T%D(1,1),1)
!    Zero the Lower Triange
     DO I=1,3
        DO J=1,I-1
           LatFrc_T%D(I,J)   = Zero
        ENDDO
     ENDDO
!    Sum in the T contribution to total force
     DO AtA=1,NAtoms
        A1=3*(AtA-1)+1
        A2=3*AtA
        GM%Gradients%D(1:3,AtA) =  GM%Gradients%D(1:3,AtA)+TFrc%D(A1:A2)
     ENDDO
!    Sum in the T contribution to total Lattice force
     GM%PBC%LatFrc%D = GM%PBC%LatFrc%D+LatFrc_T%D
#ifdef PARALLEL  
  ENDIF
#endif
! Do some printing
  CALL Print_Force(GM,TFrc,'T Force')
  CALL Print_Force(GM,TFrc,'T Force',Unit_O=6)
  CALL Print_LatForce(GM,LatFrc_T%D,'T Lattice Force')
  CALL Print_LatForce(GM,LatFrc_T%D,'T Lattice Force',Unit_O=6)
! Do some checksumming, resumming and IO 
  CALL PChkSum(TFrc,    'dT/dR',Proc_O=Prog)  
  CALL PChkSum(LatFrc_T,'LFrcT',Proc_O=Prog)  
! Save Forces to Disk
#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL MondoLog(DEBUG_NONE, Prog, "writing forces to hdf", "Clone "//TRIM(IntToChar(MyClone)))
    CALL Put(GM,Tag_O=CurGeom)
  ELSE
    CALL MondoLog(DEBUG_NONE, Prog, "sending forces to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
  ENDIF
#else
  CALL Put(GM,Tag_O=CurGeom)
#endif
! Tidy up 
  CALL Delete(TFrc)
  CALL Delete(LatFrc_T)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL ShutDown(Prog)
END PROGRAM TForce
