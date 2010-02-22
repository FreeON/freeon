! vim: tw=0
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
!    COMPUTE THE OVERLAP MATRIX S
!    Authors: Matt Challacombe and C.J. Tymczak
!------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM MakeS
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
  USE OverlapBlock
  USE BraBloks
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                :: S,T1
#else
  TYPE(BCSR)                 :: S,T1
#endif
  INTEGER                    :: NC
  REAL(DOUBLE),DIMENSION(3)  :: B
  TYPE(AtomPair)             :: Pair
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(INT_VECT)             :: Stat
  TYPE(ARGMT)                :: Args
  INTEGER                    :: P,R,AtA,AtB,NN,OldFileID                        
  CHARACTER(LEN=5),PARAMETER :: Prog='MakeS'

! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry
  IF(SCFActn=='RestartGeomSwitch') THEN
     CALL Get(BS,Tag_O=CurBase)
     ! Close current HDF group 
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! Open the old group and HDF
     HDF_CurrentID=OpenHDF(Restart)
     OldFileID=HDF_CurrentID 
     CALL New(Stat,3)
     CALL Get(Stat,'current_state')
     HDF_CurrentID=OpenHDFGroup(HDF_CurrentID,"Clone #"//TRIM(IntToChar(MyClone)))
     ! Get old basis set stuff
     CurGeom=TRIM(IntToChar(Stat%I(3)))
     CurBase=TRIM(IntToChar(Stat%I(2)))
     CALL Get(GM,Tag_O=CurGeom)
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
  ELSE
     CALL Get(BS,Tag_O=CurBase)
     CALL Get(GM,Tag_O=CurGeom)
  ENDIF
!---------------------------------------------- 
! Allocations 
  CALL NewBraBlok(BS)
  CALL New(S)
!-----------------------------------------------
! Initialize the matrix and associated indecies

  P=1; R=1; S%RowPt%I(1)=1
  CALL SetEq(S%MTrix,Zero)
!-----------------------------------------------
! Main loops
#ifdef PARALLEL
  S%NAtms=0
  DO AtA=Beg%I(MyId),End%I(MyId)
     S%NAtms=S%NAtms+1  
#else
  S%NAtms=NAtoms
  DO AtA=1,NAtoms
#endif
     DO AtB=1,NAtoms
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
           NN = Pair%NA*Pair%NB
           B = Pair%B
           DO NC = 1,CS_OUT%NCells
              Pair%B = B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                        + (Pair%A(2)-Pair%B(2))**2 &
                        + (Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 S%MTrix%D(R:R+NN-1)=S%MTrix%D(R:R+NN-1)+SBlok(BS,Pair)
              ENDIF
           ENDDO
           S%ColPt%I(P)=AtB
           S%BlkPt%I(P)=R
           R=R+NN 
           P=P+1 
#ifdef PARALLEL           
           S%RowPt%I(S%NAtms+1)=P
           IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
              WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
              WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
              CALL Halt(' DBCSR dimensions blown in MakeS ')
           ENDIF
#else  
           S%RowPt%I(AtA+1)=P        
           IF(R>MaxNon0.OR.P>MaxBlks) THEN
              CALL MondoLog(DEBUG_NONE, Prog, "MaxNon0 = "//TRIM(IntToChar(MaxNon0)))
              CALL MondoLog(DEBUG_NONE, Prog, "R = "//TRIM(IntToChar(R)))
              CALL MondoLog(DEBUG_NONE, Prog, "MaxBlks = "//TRIM(IntToChar(MaxBlks)))
              CALL MondoLog(DEBUG_NONE, Prog, "P = "//TRIM(IntToChar(P)))
              CALL Halt(' BCSR dimensions blown in MakeS ')
           ENDIF
#endif
        ENDIF
     ENDDO
  ENDDO
  S%NBlks=P-1
  S%NNon0=R-1
!------------------------------------------------------------
! Put S to disk
  Thresholds%Trix = Thresholds%Trix*1.D-2
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "S (before filter) checksum = "//TRIM(DblToChar(CheckSum(S))), "Clone "//TRIM(IntToChar(MyClone)), serialize_O = .TRUE.)
  CALL Filter(T1,S)
  Thresholds%Trix = Thresholds%Trix*1.D2
  IF(SCFActn=='RestartGeomSwitch') THEN
     CALL Put(T1,TrixFile('S',Args,Stats_O=(/Current(1),Current(2),Current(3)-1/)))
  ELSE
     CALL Put(T1,TrixFile('S',Args))
  ENDIF
!-----------------------------------------------------------
! Printing
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "S (after filter) checksum = "//TRIM(DblToChar(CheckSum(T1))), "Clone "//TRIM(IntToChar(MyClone)), serialize_O = .TRUE.)
  CALL PChkSum(T1,'S Clone '//TRIM(IntToChar(MyClone)),Prog)
  CALL PPrint( T1,'S')
  CALL Plot(   T1,'S')
!------------------------------------------------------------
! Tidy up
  CALL Delete(S)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL DeleteBraBlok()
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM MakeS
