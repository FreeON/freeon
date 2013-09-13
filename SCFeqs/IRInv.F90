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
!    COMPUTE THE S^(-1/2) of the overlap Matrix
!    Author: Anders Niklasson and C. J. Tymczak
!--------------------------------------------------------------
PROGRAM IRInv
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
  USE SetXYZ
  USE DenMatMethods
  IMPLICIT NONE
!
  TYPE(BCSR)                     :: S,SInvL,SInvR,Tmp1,Tmp2,Tmp3
!
  TYPE(BSET)                     :: BS
  TYPE(CRDS)                     :: GM
  TYPE(ARGMT)                    :: Args
  TYPE(AtomPair)                 :: Pair
  INTEGER                        :: I,J,K,L,IC
  REAL(DOUBLE)                   :: Error,Smax,Smin,Factor
  LOGICAL                        :: ErrorTrue
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=8),PARAMETER     :: Prog='IRInv'
!----------------------------------------------------------------------------
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Get the Overlap Matrix
  CALL Get(S,TrixFile('S',Args))
!----------------------------------------------------------------------------
! Allocate
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(Tmp3)
  CALL New(SInvL)
  CALL New(SInvR)
!-------------------------------------------------
! Calculate Eigen Bounds of S
!  CALL EigenBounds(S,Smin,Smax)
  Smin=5D-3
  Smax=5.D0
!-------------------------------------------------
! Iterate Using Anders Algorithm
! Rescale I as the Initial Guess
  Factor = SQRT(1.9D0)/SQRT(ABS(Smax))
  CALL SetToI(SInvL)
  CALL SetToI(SInvR)
  CALL Multiply(SInvL,Factor)
  CALL Multiply(SInvR,Factor)
! Calculate the Inverse of S
  IC        = 0
  ErrorTrue = .FALSE.
  DO I=1,20
!    Calculate X
     CALL Multiply(SInvL,S,Tmp1)
     CALL Filter(Tmp3,Tmp1)
     WRITE(*,'(F8.2,1x,F8.2)') 100.D0*DBLE(Tmp3%NNon0)/DBLE(NBasF*NBasF),100.D0*DBLE(Tmp1%NNon0)/DBLE(NBasF*NBasF)
     CALL Multiply(Tmp3,SInvR,Tmp2)
     CALL Filter(Tmp1,Tmp2)
     WRITE(*,'(F8.2,1x,F8.2)') 100.D0*DBLE(Tmp1%NNon0)/DBLE(NBasF*NBasF),100.D0*DBLE(Tmp2%NNon0)/DBLE(NBasF*NBasF)
!    Compute the Error
     CALL Add(Tmp2,-One)
     Error = MAX(Tmp2)
     WRITE(*,'(i4,1x,D14.6,1x,F8.2,1x,F8.2)') I,Error,                                   &
                                           100.D0*DBLE(SInvL%NNon0)/DBLE(NBasF*NBasF),   &
                                           100.D0*DBLE(Tmp1%NNon0)/DBLE(NBasF*NBasF)
     IF(Error < Thresholds%Trix .AND. .NOT. ErrorTrue) ErrorTrue = .TRUE.
     IF(ErrorTrue) IC = IC + 1
     IF(IC > 1) EXIT
!    Calculate X2
     CALL Multiply(Tmp1,Tmp1,Tmp3)
     CALL Filter(Tmp2,Tmp3)
!    Calculate 1.875*I-1.250*X+0.375*X*X
     CALL Multiply(Tmp1,-1.250D0)
     CALL Multiply(Tmp2, 0.375D0)
     CALL Add(Tmp1,Tmp2,Tmp3)
     CALL SetEq(Tmp1,Tmp3)
     CALL Add(Tmp1,1.875D0)
!    Calculate SInvL = (1.875*I-1.250*X+0.375*X*X)*SinvL
     CALL Multiply(Tmp1,SInvL,Tmp2)
     CALL Filter(SInvL,Tmp2)
!    Calculate SInvR = SInvR*(1.875*I-1.250*X+0.375*X*X)
     CALL Multiply(SInvR,Tmp1,Tmp2)
     CALL Filter(SInvR,Tmp2)
  ENDDO
!---------------------------------------------------
  IF(.TRUE.) STOP
!---------------------------------------------------
!!$    IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
!!$       CALL PPrint(PerfMon,Prog)
!!$       CALL PPrint(PerfMon,Prog,Unit_O=6)
!!$    ENDIF
!!$    ! Consistency check
!!$    IF(TEST_AINV)THEN
!!$       CALL New(T1)
!!$       CALL New(T2)
!!$       CALL Multiply(Zt,A,T1)
!!$       CALL Multiply(T1,Z,T2)
!!$       !
!!$       CALL SetToI(T1)
!!$       CALL Multiply(T1,-One)
!!$       CALL Add(T1,T2,A)
!!$       Mx0=Max(A)
!!$       !
!!$       IF(AInvDistanceThresh/=Zero)THEN
!!$          Mssg='Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
!!$               //', DistThrsh='//TRIM(DblToShrtChar(AInvDistanceThresh))  &
!!$               //', TrixThrsh='//TRIM(DblToShrtChar(Thresholds%Trix))
!!$       ELSE
!!$          Mssg='Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
!!$               //', TrixThrsh='//TRIM(DblToShrtChar(Thresholds%Trix))
!!$       ENDIF
!!$       IF(Mx0>1.D2*Thresholds%Trix)THEN
!!$          CALL Warn('In BlokAInv, failed test: '//TRIM(Mssg))
!!$       ENDIF
!!$       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!!$          Mssg=ProcessName(Prog)//TRIM(Mssg)
!!$          WRITE(*,*)TRIM(Mssg)
!!$          CALL OpenASCII(OutFile,Out)
!!$          CALL PrintProtectL(Out)
!!$          WRITE(Out,*)TRIM(Mssg)
!!$          CALL PrintProtectR(Out)
!!$          CLOSE(UNIT=Out,STATUS='KEEP')
!!$       ENDIF
!!$    ENDIF
!!$    !  Put Z and ZT to disk
!!$    CALL Put(Z,TrixFile('Z',Args))
!!$    CALL Put(Zt,TrixFile('ZT',Args))
!!$    !  Debug
!!$    CALL PChkSum(Z,'Z',Proc_O=Prog)
!!$    CALL PPrint(Z,'Z')
!!$    CALL Plot(Z,'Z')
!!$    !  Tidy up
!!$    CALL Delete(Z)
!!$    CALL Delete(ZT)
    CALL ShutDown(Prog)
!
END PROGRAM IRInv
