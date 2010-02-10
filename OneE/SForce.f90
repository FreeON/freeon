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
!    COMPUTE THE FORCE DUE TO CHANGES IN THE DENSITY MATRIX:
!    dP.(2T+J+K)=-2*dS.P.F.P (Early work by McWeeny, need a cite...)
!    Authors: Matt Challacombe and CJ Tymczak
!------------------------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM SForce
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
  USE BraBloks
  USE BlokTrWdS
  USE MondoLogger

#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)                :: T1,F,P_DBCSR
  TYPE(BCSR)                 :: P
  TYPE(DBL_VECT)             :: TotSFrc
  TYPE(DBL_RNK2)             :: TmpLatFrc_S
  INTEGER                    :: IErr,TotFrcComp
#else
  TYPE(BCSR)                 :: F,P,S,Z,ZT,T1,T2,T3
#endif
  TYPE(AtomPair)             :: Pair
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(ARGMT)                :: Args
  INTEGER                    :: Q,R,AtA,AtB,NN,JP,MB,MA,NB,MN1,MN,A1,A2
  TYPE(HGRho)                :: Rho
  TYPE(DBL_VECT)             :: SFrc
  REAL(DOUBLE)               :: SFrcChk

  INTEGER                        :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)      :: A,B,F_nlm,nlm
  TYPE(DBL_RNK2)                 :: LatFrc_S
  CHARACTER(LEN=6),PARAMETER     :: Prog='SForce'
  LOGICAL                        :: Present

#if defined(PARALLEL_CLONES)
  INTEGER :: oldClone, rank
#endif

  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)

  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)

  ! Allocations
#if ! defined(PARALLEL)
  CALL New(P,OnAll_O=.TRUE.)
  CALL New(F)
  CALL New(S)
  CALL New(Z)
  CALL New(ZT)
  CALL New(T1)
  CALL New(T2)
  CALL New(T3)
#endif

  ! Compute W=P.F.P
#ifdef PARALLEL
  CALL Get(P_DBCSR,TrixFile('D',Args,1))
#else
  CALL Get(P,TrixFile('OrthoD',Args,1))
#endif
  ! Is this a bug if we don't use the extrapolated Fockian?
  INQUIRE(FILE=TrixFile('F_DIIS',Args,0),EXIST=Present)
  IF(Present) THEN
!    CALL MondoLog(DEBUG_NONE, "SForce", "getting F from "//TRIM(TrixFile('F_DIIS',Args,0)))
    CALL Get(F,TrixFile('F_DIIS',Args,0))
  ELSE
!    CALL MondoLog(DEBUG_NONE, "SForce", "getting F from "//TRIM(TrixFile('OrthoF',Args,0)))
    CALL Get(F,TrixFile('OrthoF',Args,0))
  ENDIF

  INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
  IF(Present)THEN
    CALL Get(Z,TrixFile("X",Args))
    CALL Get(ZT,TrixFile("X",Args))
  ELSE
    CALL Get(Z,TrixFile("Z",Args))
    CALL Get(ZT,TrixFile("ZT",Args))
  ENDIF

!!$  ! This tests the commutator of P and F. This should be zero.
!!$  CALL Multiply(P, F, T1)
!!$  CALL Multiply(F, P, T2)
!!$  CALL Multiply(T2, -1D0)
!!$  CALL Add(T1, T2, T3)
!!$  CALL MondoLog(DEBUG_NONE, "SForce", "DIIS error = "//TRIM(DblToChar(SQRT(Dot(T3, T3))/DBLE(NBasF))))

#ifdef PARALLEL
  CALL Multiply(P_DBCSR,F,T1)       ! T1:=P.F
  CALL Multiply(T1,P_DBCSR,F)       ! F:=P.F.P
  CALL Filter(P_DBCSR,F)            ! P=Filter[P.F.P]
  CALL SetEq(P,P_DBCSR)
  CALL BcastBCSR(P)
  CALL Delete(P_DBCSR)
#else

!!  old
!!  CALL Multiply(P,F,T1)       ! T1:=P.F
!!  CALL Multiply(T1,P,F)       ! F:=P.F.P
!!  CALL Filter(P,F)            ! P=Filter[P.F.P]
!!  CALL Multiply(Z,P,T1)       ! P:= Z.P.ZT
!!  CALL Multiply(T1,ZT,P)      ! Back AO basis.
!! Begin Change
  CALL Multiply(F,P,T1)       ! T1:=F.P
  CALL Filter(P,T1)            ! P=Filter[P.F]
  CALL Multiply(Z,P,T1)       ! P:= Z.P.ZT
  CALL Multiply(T1,ZT,P)      ! Back AO basis.
!! End Change

#endif

  CALL Delete(F)
  CALL Delete(S)
  CALL Delete(Z)
  CALL Delete(T1)
  CALL Delete(T2)
  CALL Delete(T3)

  CALL New(SFrc,3*NAtoms)
  SFrc%D = Zero
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(LatFrc_S,(/3,3/))
  LatFrc_S%D = Zero

  !--------------------------------------------------------------------------------
  ! SForce=-2*Tr{P.F.P.dS} (Extra 2 to account for symmetry of S in the trace)
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
        ! Quick and dirty!
        SELECT CASE(P%NSMat)
        CASE(1)
          ! We don't need to do anything!
        CASE(2)
          ! We add up the two density martices!
          CALL DAXPY(MN,1D0,P%MTrix%D(Q+MN),1,P%MTrix%D(Q),1)
        CASE(4)
          ! We add up the diagonal density martices!
          CALL DAXPY(MN,1D0,P%MTrix%D(Q+3*MN),1,P%MTrix%D(Q),1)
        CASE DEFAULT
          CALL Halt(' SForce: P%NSMat doesn''t have an expected value! ')
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
            F_nlm(1:3)    = TrWdS(BS,Pair,P%MTrix%D(Q:Q+MN1))
            SFrc%D(A1:A2) = SFrc%D(A1:A2) - Two*F_nlm(1:3)

            ! Lattice Forces
            nlm        = AtomToFrac(GM,Pair%A)
            LatFrc_S%D = LatFrc_S%D - Two*LaticeForce(GM,nlm,F_nlm)
          ENDIF
          Pair%A=A+CS_OUT%CellCarts%D(:,NC)
          Pair%B=B
          Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
               +(Pair%A(2)-Pair%B(2))**2  &
               +(Pair%A(3)-Pair%B(3))**2
          IF(TestAtomPair(Pair)) THEN
            F_nlm(1:3)    = TrWdS(BS,Pair,P%MTrix%D(Q:Q+MN1))
            SFrc%D(A1:A2) = SFrc%D(A1:A2) - Two*F_nlm(1:3)
            ! Lattice Forces
            nlm        = AtomToFrac(GM,Pair%A)
            LatFrc_S%D = LatFrc_S%D - Two*LaticeForce(GM,nlm,F_nlm)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO

#ifdef PARALLEL
  ! Collect the Forces
  TotFrcComp = 3*NAtoms
  CALL New(TotSFrc,TotFrcComp)
  CALL MPI_Reduce(SFrc%D(1),TotSFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
    SFrc%D(1:TotFrcComp) = TotSFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotSFrc)
  ! Collect the Lattice Forces
  CALL New(TmpLatFrc_S,(/3,3/))
  CALL DBL_VECT_EQ_DBL_SCLR(9,TmpLatFrc_S%D(1,1),0.0d0)
  CALL MPI_REDUCE(LatFrc_S%D(1,1),TmpLatFrc_S%D(1,1),9,MPI_DOUBLE_PRECISION,MPI_SUM,ROOT,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
    LatFrc_S%D = TmpLatFrc_S%D
  ENDIF
  CALL Delete(TmpLatFrc_S)
#endif
#ifdef PARALLEL
  IF(MyID == ROOT) THEN
#endif
    ! Rescale the Forces if needed.
    IF(P%NSMat.GT.1) CALL DSCAL(3*NAtoms,0.5D0,    SFrc%D(1  ),1)
    IF(P%NSMat.GT.1) CALL DSCAL(       9,0.5D0,LatFrc_S%D(1,1),1)
    ! Zero the Lower Triange
    DO I=1,3
      DO J=1,I-1
        LatFrc_S%D(I,J) = Zero
      ENDDO
    ENDDO
    ! Sum in the S contribution to total force
    DO AtA=1,NAtoms
      A1=3*(AtA-1)+1
      A2=3*AtA
      GM%Gradients%D(1:3,AtA) =  SFrc%D(A1:A2)
    ENDDO
    ! Sum in the S contribution to total lattice force
    GM%PBC%LatFrc%D = LatFrc_S%D
#ifdef PARALLEL
  ENDIF
#endif

  ! Do some printing
  CALL Print_Force(GM,SFrc,'S Force')
  CALL Print_LatForce(GM,LatFrc_S%D,'S Lattice Force')

  ! Do some checksumming and IO
  CALL PChkSum(SFrc,    'dS/dR',Proc_O=Prog)
  CALL PChkSum(LatFrc_S,'LFrcS',Proc_O=Prog)

  ! Save Forces to Disk
#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL Put(GM, Tag_O = CurGeom)

    oldClone = MyClone
    DO rank = 1, MSize(MPI_COMM_WORLD)-1
      CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(GM, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

      ! Put to correct HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
      CALL Put(GM, Tag_O = CurGeom)
    ENDDO
    MyClone = oldClone

    ! Reopen old HDFGroup.
    CALL CloseHDFGroup(H5GroupID)
    H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID = H5GroupID
  ELSE
    !CALL MondoLog(DEBUG_MAXIMUM, Prog, "sending density and multipoles to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
    CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(GM, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
  ENDIF
#else
  CALL Put(GM, Tag_O = CurGeom)
#endif

  ! Tidy up
  CALL Delete(SFrc)
  CALL Delete(LatFrc_S)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  CALL ShutDown(Prog)

END PROGRAM SForce
