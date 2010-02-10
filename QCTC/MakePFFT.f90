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
PROGRAM MakePFFT
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE BraBloks
  USE PFFTen
  USE DivPFFTen
  USE AtomPairs

  IMPLICIT NONE

  TYPE(TIME)                     :: TimePFFT
  TYPE(DBL_VECT)                 :: TenC,TenS
  TYPE(DBL_RNK2)                 :: BoxShape
  TYPE(DBL_RNK3)                 :: dTenC,dTenS
  TYPE(CellSet)                  :: CSS
  INTEGER                        :: MaxEll,I,J,K
  REAL(DOUBLE)                   :: DDelta,Rad,AtoAU,SUM
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=8),PARAMETER     :: Prog='MakePFFT'

  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
  ! Get the geometry
  CALL Get(GM,Tag_O=CurGeom)
  ! If Periodic, Make Tensors
  IF(GM%PBC%Dimen>0)THEN
    IF(GM%PBC%PFFMaxEll>FFELL) THEN
      CALL MondoHalt(0,'MaxELL > FFELL Halting in CalculatePFFT')
    ENDIF
    !    Set up the multipole arrays
    CALL MultipoleSetUp()
    !    Allocate the tensors
    MaxEll=GM%PBC%PFFMaxEll
    CALL New(TenC,LSP(2*MaxEll),0)
    CALL New(TenS,LSP(2*MaxEll),0)
    TenC%D=Zero
    TenS%D=Zero
    !    Calculate the tensors ...
    CALL CalculatePFFT(MaxEll,GM,Args,GM%InCells,TenC,TenS)
    !    Put them to HDF
    CALL Put(TenC,'PFFTensorC')
    CALL Put(TenS,'PFFTensorS')
    !    Do the Derivaltive Arrays
    !    First, Allocate the derivative tensors
    CALL New(dTenC,(/LSP(2*MaxEll),3,3/),(/0,1,1/))
    CALL New(dTenS,(/LSP(2*MaxEll),3,3/),(/0,1,1/))
    dTenC%D=Zero
    dTenS%D=Zero

    IF(GM%PBC%Dimen/=2)THEN
       !    Calculate the derivative Tensors
       CALL CalculateDivPFFT(MaxEll,GM,Args,GM%InCells,dTenC,dTenS)
    ELSE
       CALL New(BoxShape,(/3,3/))
       CALL New(CSS,CS_IN%NCells)
       DDelta = 1.D-5
       !       Initialize CS
       BoxShape%D=GM%PBC%BoxShape%D
       dTenC%D=Zero
       dTenS%D=Zero
       DO K=1,CS_IN%NCells
          CS_IN%CellCarts%D(:,K) = AtomToFrac(GM,CS_IN%CellCarts%D(:,K))
       ENDDO
!       Calculate the tensors numerically
       DO I=1,3
          DO J=1,3
             IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
                WRITE(*,*) 'I,J = ',I,J
                GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J) + DDelta
                GM%PBC%InvBoxSh%D      = InverseMatrix3x3(GM%PBC%BoxShape%D)
                GM%PBC%CellVolume      = ABS(CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I))
!
                DO K=1,CS_IN%NCells
                   CSS%CellCarts%D(:,K) = FracToAtom(GM,CS_IN%CellCarts%D(:,K))
                ENDDO
!
                CALL CalculatePFFT(MaxEll,GM,Args,CSS,TenC,TenS)
                dTenC%D(:,I,J) = TenC%D(:)
                dTenS%D(:,I,J) = TenS%D(:)
                !
                GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J) - DDelta
                GM%PBC%InvBoxSh%D      = InverseMatrix3x3(GM%PBC%BoxShape%D)
                GM%PBC%CellVolume      = ABS(CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I))
                DO K=1,CS_IN%NCells
                   CSS%CellCarts%D(:,K) = FracToAtom(GM,CS_IN%CellCarts%D(:,K))
                ENDDO
                !
                CALL CalculatePFFT(MaxEll,GM,Args,CSS,TenC,TenS)
                dTenC%D(:,I,J) = (TenC%D(:)-dTenC%D(:,I,J))/(Two*DDelta)
                dTenS%D(:,I,J) = (TenS%D(:)-dTenS%D(:,I,J))/(Two*DDelta)

                GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J)
                GM%PBC%InvBoxSh%D      = InverseMatrix3x3(GM%PBC%BoxShape%D)
                GM%PBC%CellVolume      = ABS(CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I))
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !    Put them to HDF
    CALL Put(dTenC,'dPFFTensorC')
    CALL Put(dTenS,'dPFFTensorS')
    !    Print sum Checksums
    CALL PChkSum(TenC,'TenC',Proc_O=Prog)
    CALL PChkSum(TenS,'TenS',Proc_O=Prog)
!!$
!!$     CALL PChkSum(TenC,'TenC',Proc_O=Prog,Unit_O=6)
!!$     CALL PChkSum(TenS,'TenS',Proc_O=Prog,Unit_O=6)
    !    Delete
    CALL Delete(TenC)
    CALL Delete(TenS)
    CALL Delete(dTenC)
    CALL Delete(dTenS)
  ENDIF
  ! Delete
  CALL Delete(GM)
  CALL Delete(Args)
  ! didn't count flops, any accumulation is residual
  ! from matrix routines
  PerfMon%FLOP=Zero
  ! Shutdown
  CALL ShutDown(Prog)
  !
END PROGRAM MakePFFT
!!$     LM = LTD(2)
!!$     IF(.TRUE.) THEN
!!$        CALL CalculateDivPFFT(MaxEll,GM,Args,CS_IN,dTenC,dTenS)
!!$        CALL Print_SP(4,dTenC%D(:,1,1),dTenS%D(:,1,1),Tag_O='dTen(1,1)',Pre_O='Long')
!!$        WRITE(*,*) 'dTenC(L=2,M=2)'
!!$        DO I=1,3
!!$           WRITE(*,*) (dTenC%D(LM,I,J),J=1,3)
!!$        ENDDO
!!$        WRITE(*,*) 'dTenS(L=2,M=2)'
!!$        DO I=1,3
!!$           WRITE(*,*) (dTenS%D(LM,I,J),J=1,3)
!!$        ENDDO
!!$     ENDIF
!!$!
!!$     IF(.TRUE.) THEN
!!$        CALL New(BoxShape,(/3,3/))
!!$        CALL New(CSS,CS_IN%NCells)
!!$        DDelta = 1.D-6
!!$!       Initialize CS
!!$        BoxShape%D=GM%PBC%BoxShape%D
!!$        dTenC%D=Zero
!!$        dTenS%D=Zero
!!$        DO K=1,CS_IN%NCells
!!$           CS_IN%CellCarts%D(:,K) = AtomToFrac(GM,CS_IN%CellCarts%D(:,K))
!!$        ENDDO
!!$!       Calculate the tensors numerically
!!$
!!$        DO I=1,3
!!$           DO J=1,3
!!$              IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
!!$                 WRITE(*,*) 'I,J = ',I,J
!!$                 GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J) + DDelta
!!$                 GM%PBC%InvBoxSh%D      = InverseMatrix(GM%PBC%BoxShape%D)
!!$                 GM%PBC%CellVolume      = ABS(CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I))
!!$!
!!$                 DO K=1,CS_IN%NCells
!!$                    CSS%CellCarts%D(:,K) = FracToAtom(GM,CS_IN%CellCarts%D(:,K))
!!$                 ENDDO
!!$!
!!$                 CALL CalculatePFFT(MaxEll,GM,Args,CSS,TenC,TenS)
!!$                 dTenC%D(:,I,J) = TenC%D(:)
!!$                 dTenS%D(:,I,J) = TenS%D(:)
!!$!
!!$                 GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J) - DDelta
!!$                 GM%PBC%InvBoxSh%D      = InverseMatrix(GM%PBC%BoxShape%D)
!!$                 GM%PBC%CellVolume      = ABS(CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I))
!!$                 DO K=1,CS_IN%NCells
!!$                    CSS%CellCarts%D(:,K) = FracToAtom(GM,CS_IN%CellCarts%D(:,K))
!!$                 ENDDO
!!$!
!!$                 CALL CalculatePFFT(MaxEll,GM,Args,CSS,TenC,TenS)
!!$                 dTenC%D(:,I,J) = (TenC%D(:)-dTenC%D(:,I,J))/(Two*DDelta)
!!$                 dTenS%D(:,I,J) = (TenS%D(:)-dTenS%D(:,I,J))/(Two*DDelta)
!!$
!!$                 GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J)
!!$                 GM%PBC%InvBoxSh%D      = InverseMatrix(GM%PBC%BoxShape%D)
!!$                 GM%PBC%CellVolume      = ABS(CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I))
!!$
!!$              ENDIF
!!$           ENDDO
!!$        ENDDO
!!$        CALL Print_SP(4,dTenC%D(:,1,1),dTenS%D(:,1,1),Tag_O='dTen(1,1)',Pre_O='Long')
!!$        WRITE(*,*) 'dTenC(L=2,M=2)'
!!$        DO I=1,3
!!$           WRITE(*,*) (dTenC%D(LM,I,J),J=1,3)
!!$        ENDDO
!!$        WRITE(*,*) 'dTenS(L=2,M=2)'
!!$        DO I=1,3
!!$           WRITE(*,*) (dTenS%D(LM,I,J),J=1,3)
!!$        ENDDO
!!$     ENDIF
!!$     IF(.TRUE.) STOP
