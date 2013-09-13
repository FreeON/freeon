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
PROGRAM DDIIS
!H=================================================================================
!H PROGRAM DDIIS
!H
!H  OPTIONS:
!H  DEBUGING: Use -DDDIIS_DBUG to print some stuff.
!H  INFO    : Use -DDDIIS_INFO to print some stuff.
!H
!H Comment:
!H
!H Ref:
!H  V. Weber, C. Daul Chem. Phys. Let. 370, 99-105, 2003.
!H
!H=================================================================================
  !
#ifdef DDIIS_DBUG
#define DDIIS_INFO
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE SetXYZ
  USE LinAlg
  USE MatFunk
#ifdef PARALLEL
  USE MondoMPI
#endif
  !
  IMPLICIT NONE
  !-------------------------------------------------------------------
#ifdef PARALLEL
  TYPE(DBCSR)                      :: F,P,EPrm,Tmp1,TmpFPrm
  TYPE(DBCSR)                      :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
  TYPE(DBCSR)                      :: FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#else
  TYPE(BCSR )                      :: F,P,EPrm,Tmp1,TmpFPrm
  TYPE(BCSR )                      :: PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
  TYPE(BCSR )                      :: FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#endif
  TYPE(ARGMT)                      :: Args
  TYPE(INT_VECT)                   :: Idx,SCFOff
  TYPE(DBL_VECT)                   :: V,DIISCo,AbsDIISCo
  TYPE(DBL_RNK2)                   :: B,BInv
  !-------------------------------------------------------------------
  REAL(DOUBLE)                     :: DIISErr,Damp,EigThresh
  REAL(DOUBLE)                     :: CondA,DDIISMaxCond
  INTEGER                          :: I,J,I0,J0,N,M,BMax,DoDIIS,iOffSet
  INTEGER                          :: CPSCFCycl,DDIISStart
  INTEGER                          :: LastSCFCycle,LastCPSCFCycle
  INTEGER                          :: DDIISBeg,DDIISEnd,DDIISCurDim
  INTEGER                          :: RespOrder
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=1)                 :: Chr1,Chr2,Chr3
  CHARACTER(LEN=*), PARAMETER      :: Prog='DDIIS'
  LOGICAL                          :: IsPresent
  !-------------------------------------------------------------------
  REAL(DOUBLE), PARAMETER          :: DEFAULT_EIGTHRESH = 1.00D-10
  REAL(DOUBLE), PARAMETER          :: DEFAULT_DAMP      = 1.00D-01 ! DEFAULT_DAMP_MIN < damp < DEFAULT_DAMP_MAX
  REAL(DOUBLE), PARAMETER          :: DEFAULT_DAMP_MAX  = 5.00D-01
  REAL(DOUBLE), PARAMETER          :: DEFAULT_DAMP_MIN  = 1.00D-03
  REAL(DOUBLE), PARAMETER          :: DEFAULT_MAXCOND   = 1.00D+08
  INTEGER     , PARAMETER          :: DEFAULF_BMAX      = 15
  INTEGER     , PARAMETER          :: DEFAULF_START     = 1
  !-------------------------------------------------------------------
  !type(DBL_RNK2) :: BTmp
  integer :: ierr
  !
  ! Initial setup.
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  !
  CPSCFCycl=Args%I%I(1)
  !
  ! Get the response order.
  RespOrder=LEN(TRIM(Args%C%C(3)))
  !
  ! Get Last SCF cycle.
  LastSCFCycle=0
  LastCPSCFCycle=0
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  ! Get the directions.
  Chr1=TRIM(Args%C%C(3)(1:1))
  SELECT CASE(RespOrder)
  CASE(1);Chr2=' '                   ;Chr3=' '
  CASE(2);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=' '
  CASE(3);Chr2=TRIM(Args%C%C(3)(2:2));Chr3=TRIM(Args%C%C(3)(3:3))
  CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Parse for DDIIS options.
  !-------------------------------------------------------------------
  !
  ! Open input
  CALL OpenASCII(InpFile,Inp)
  !
  ! Threshold for projection of small eigenvalues.
  IF(.NOT.OptDblQ(Inp,'DDIISThresh',EigThresh)) EigThresh=DEFAULT_EIGTHRESH ! EigThresh=1.0D-10
  IF(EigThresh.LT.Zero) &
         CALL Halt(' DDIISThresh cannot be smaller than zero, DDIISThresh=' &
                   //TRIM(DblToShrtChar(EigThresh))//'.')
  !
  ! Damping coefficient for first cycle
  IF(.NOT.OptDblQ(Inp,'DDIISDamp',Damp)) Damp=DEFAULT_DAMP !Damp=1.0D-1
  IF(Damp.LT.Zero) &
         CALL Halt(' DDIISDamp cannot be smaller than zero, DDIISDamp=' &
                   //TRIM(DblToShrtChar(Damp))//'.')
  ! Dont allow damping below DEFAULT_DAMP_MIN, as this can cause false convergence.
  Damp=MAX(Damp,DEFAULT_DAMP_MIN)
  ! Dont allow damping above DEFAULT_DAMP_MAX, as this is silly (DIIS would certainly work better).
  Damp=MIN(Damp,DEFAULT_DAMP_MAX)
  !
  ! Max number of equations to keep in DIIS.
  IF(.NOT.OptIntQ(Inp,'DDIISDimension',BMax)) BMax=DEFAULF_BMAX !BMax=15
  IF(BMax.LT.0) &
         CALL Halt(' DDIISDimension cannot be smaller than zero, DDIISDimension=' &
                   //TRIM(IntToChar(BMax))//'.')
  !
  ! Condition number of B.
  IF(.NOT.OptDblQ(Inp,'DDIISMaxCond',DDIISMaxCond)) DDIISMaxCond=DEFAULT_MAXCOND !DDIISMaxCond=1.0D+8
  IF(DDIISMaxCond.LT.Zero) &
         CALL Halt(' DDIISMaxCond cannot be smaller than zero, DDIISMaxCond=' &
                   //TRIM(DblToShrtChar(DDIISMaxCond))//'.')
  !
  ! The iteration where we want to start the DDIIS.
  IF(.NOT.OptIntQ(Inp,'DDIISStart',DDIISStart)) DDIISStart=DEFAULF_START !DDIISStart=1
  IF(DDIISStart.LT.1) &
         CALL Halt(' DDIISStart cannot be smaller than one, DDIISStart=' &
                   //TRIM(IntToChar(DDIISStart))//'.')
  !
  ! Close the input.
  CLOSE(Inp)
  !
  !-------------------------------------------------------------------
  ! Do we need to do the DDIIS?
  !-------------------------------------------------------------------
  !
  IF(CPSCFCycl<=DDIISStart)THEN
     ! No DIIS yet, but damp non-extrapolated Fock matrices
     DoDIIS=-1
  ELSEIF(BMax/=0)THEN
     ! We are doing DIIS, extrapolating non-extrapolated Fock matrices
     DoDIIS=1
  ELSEIF(BMax==0)THEN
     ! We are purely damping, using previously extrapolated Fock matrices
     DoDIIS=0
  ENDIF
  !
  !-------------------------------------------------------------------
  ! Initialize or get Beg and End DDIIS variables.
  !-------------------------------------------------------------------
  !
  !CALL New(BTmp,(/BMax+1,BMax+1/))
  IF(CPSCFCycl.LE.1) THEN
     DDIISBeg=1
     DDIISEnd=1
     CALL Put(DDIISBeg,'DDIISBeg'//IntToChar(CPSCFCycl))
     CALL Put(DDIISEnd,'DDIISEnd'//IntToChar(CPSCFCycl))
     !BTmp%D=Zero
     !CALL Put(BTmp,'DDIISBMtrix')
  ELSE
     CALL Get(DDIISBeg,'DDIISBeg'//IntToChar(CPSCFCycl))
     CALL Get(DDIISEnd,'DDIISEnd'//IntToChar(CPSCFCycl))
     !CALL Get(BTmp,'DDIISBMtrix')
  ENDIF
#ifdef DDIIS_DBUG
  WRITE(*,*) 'DDIISBeg ',DDIISBeg
  WRITE(*,*) 'DDIISEnd ',DDIISEnd
  WRITE(*,*) 'CPSCFCycl',CPSCFCycl
#endif
  !CALL PrintMatrix(BTmp%D,BMax+1,BMax+1,2)
  !
  !-------------------------------------------------------------------
  ! Allocations.
  !-------------------------------------------------------------------
  !
  CALL New(Tmp1)
  CALL New(EPrm)
  CALL New(F)
  CALL New(P)
  SELECT CASE(RespOrder)
  CASE(1)
     CALL New(FPrm1_1)
     CALL New(PPrm1_1)
  CASE(2)
     CALL New(FPrm1_1)
     CALL New(PPrm1_1)
     CALL New(FPrm2_1)
     CALL New(PPrm2_1)
     IF(Chr1.EQ.Chr2) THEN
     ELSE
        CALL New(FPrm1_2)
        CALL New(PPrm1_2)
     ENDIF
  CASE(3)
     CALL New(FPrm1_1)
     CALL New(PPrm1_1)
     CALL New(FPrm2_1)
     CALL New(PPrm2_1)
     CALL New(FPrm3_1)
     CALL New(PPrm3_1)
     IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
     ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
        CALL New(FPrm1_3)
        CALL New(PPrm1_3)
        CALL New(FPrm2_2)
        CALL New(PPrm2_2)
     ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
        CALL New(FPrm1_2)
        CALL New(PPrm1_2)
        CALL New(FPrm2_3)
        CALL New(PPrm2_3)
     ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
        CALL New(FPrm1_2)
        CALL New(PPrm1_2)
        CALL New(FPrm1_3)
        CALL New(PPrm1_3)
        CALL New(FPrm2_2)
        CALL New(PPrm2_2)
        CALL New(FPrm2_3)
        CALL New(PPrm2_3)
     ELSE
        CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
     ENDIF
  CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Loading matrices.
  !-------------------------------------------------------------------
  !
  ! Load the orthogonal GS Fock matrix.
  CALL Get(F,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)))
  !
  ! Load the orthogonal GS Density matrix.
  CALL Get(P,TrixFile('OrthoD',Args,LastSCFCycle-Args%I%I(1)))
  !
  SELECT CASE(RespOrder)
  CASE(1)
     !
     CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)),Args,0))
     !
  CASE(2)
     !
     IF(Chr1.EQ.Chr2) THEN
        ! PPrm2_1 <-> aa
        ! PPrm1_1 <-> a
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
        CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
     ELSE
        ! PPrm2_1 <-> ab
        ! PPrm1_1 <-> a
        ! PPrm1_2 <-> b
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
        CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:2)))
        CALL Get(FPrm1_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_2,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(2:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
     ENDIF
     !
     CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL Get(PPrm2_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)),Args,0))
     !
  CASE(3)
     !
     IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
        ! A.EQ.B.EQ.C
        ! PPrm3_1 <-> aaa
        ! PPrm2_1 <-> aa
        ! PPrm1_1 <-> a
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
        CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
        CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
     ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
        ! A.EQ.B.NE.C
        ! PPrm3_1 <-> aac
        ! PPrm2_1 <-> aa
        ! PPrm2_2 <-> ac
        ! PPrm1_1 <-> a
        ! PPrm1_3 <-> c
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
        CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(3:3)))
        CALL Get(FPrm1_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(3:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_3,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(3:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
        CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:3)))
        CALL Get(FPrm2_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_2,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(2:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
     ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
        ! A.NE.B.EQ.C
        ! PPrm3_1 <-> abb
        ! PPrm2_1 <-> ab
        ! PPrm2_3 <-> bb
        ! PPrm1_1 <-> a
        ! PPrm1_2 <-> b
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
        CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:2)))
        CALL Get(FPrm1_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_2,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(2:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
        CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:3)))
        CALL Get(FPrm2_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_3,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(2:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
     ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
        ! A.NE.B.NE.C
        ! PPrm3_1 <-> abc
        ! PPrm2_1 <-> ab
        ! PPrm2_2 <-> ac
        ! PPrm2_3 <-> bc
        ! PPrm1_1 <-> a
        ! PPrm1_2 <-> b
        ! PPrm1_3 <-> c
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1)))
        CALL Get(FPrm1_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:2)))
        CALL Get(FPrm1_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_2,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(2:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(3:3)))
        CALL Get(FPrm1_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(3:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm1_3,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(3:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:2)))
        CALL Get(FPrm2_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:2)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(1:1))//TRIM(Args%C%C(3)(3:3)))
        CALL Get(FPrm2_2,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(1:1))//TRIM(Args%C%C(3)(3:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_2,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(1:1))//TRIM(Args%C%C(3)(3:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        !
        CALL Get(LastCPSCFCycle,'lastcpscfcycle'//TRIM(Args%C%C(3)(2:3)))
        CALL Get(FPrm2_3,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)(2:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
        CALL Get(PPrm2_3,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)(2:3)), &
                 Args,LastCPSCFCycle-Args%I%I(1)))
     ELSE
        CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
     ENDIF
     !
     ! PPrm3_1 <-> abc
     CALL Get(FPrm3_1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL Get(PPrm3_1,TrixFile('OrthoDPrime'//TRIM(Args%C%C(3)),Args,0))
     !
  CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Build up new DDIIS error.
  !-------------------------------------------------------------------
  !
  SELECT CASE(RespOrder)
  CASE(1)
     ! Create a new error vector E1=[F_(i+1),P_i]1     (Could be done in a better way!)
     ! PPrm1_1 <-> a
     CALL Multiply(PPrm1_1,F      ,EPrm     )  !E'=P'F
     CALL Multiply(P      ,FPrm1_1,EPrm, One)  !E'=PF'+P'F
     CALL Multiply(FPrm1_1,P      ,EPrm,-One)  !E'=F'P-(PF'+P'F)
     CALL Multiply(F      ,PPrm1_1,EPrm, One)  !E'=FP'+F'P-(PF'+P'F)
  CASE(2)
     ! Create a new error vector E2=[F_(i+1),P_i]2   (Could be done in a better way!)
     ! PPrm2_1 <-> ab
     ! PPrm1_1 <-> a
     ! PPrm1_2 <-> b
     IF(Chr1.EQ.Chr2) THEN
        CALL Multiply(PPrm2_1,F      ,EPrm     )  !Eaa=Paa*F
        CALL Multiply(PPrm1_1,FPrm1_1,EPrm, One)  !Eaa=Pa*Fa
        CALL Multiply(P      ,FPrm2_1,EPrm, One)  !Eaa=P*Faa
        !
        CALL Multiply(FPrm2_1,P      ,EPrm,-One)  !Eaa=Faa*P
        CALL Multiply(FPrm1_1,PPrm1_1,EPrm, One)  !Eaa=Fa*Pa
        CALL Multiply(F      ,PPrm2_1,EPrm, One)  !Eaa=F*Paa
     ELSE
        CALL Multiply(PPrm1_1,0.5d0)
        CALL Multiply(PPrm1_2,0.5d0)
        !
        CALL Multiply(PPrm2_1,F      ,EPrm     )  !Eab=Pab*F
        CALL Multiply(PPrm1_1,FPrm1_2,EPrm, One)  !Eab=0.5*Pa*Fb
        CALL Multiply(PPrm1_2,FPrm1_1,EPrm, One)  !Eab=0.5*Pb*Fa
        CALL Multiply(P      ,FPrm2_1,EPrm, One)  !Eab=P*Fab
        !
        CALL Multiply(FPrm2_1,P      ,EPrm,-One)  !Eab=Fab*P
        CALL Multiply(FPrm1_1,PPrm1_2,EPrm, One)  !Eab=0.5*Fa*Pb
        CALL Multiply(FPrm1_2,PPrm1_1,EPrm, One)  !Eab=0.5*Fb*Pa
        CALL Multiply(F      ,PPrm2_1,EPrm, One)  !Eab=F*Pab
     ENDIF
     !
  CASE(3)
     ! Create a new error vector E3=[F_(i+1),P_i]3 (Could be done in a better way!)
     ! PPrm3_1 <-> abc
     ! PPrm2_1 <-> ab
     ! PPrm2_2 <-> ac
     ! PPrm2_3 <-> bc
     ! PPrm1_1 <-> a
     ! PPrm1_2 <-> b
     ! PPrm1_3 <-> c
     IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
        ! A.EQ.B.EQ.C
        ! PPrm3_1 <-> aaa
        ! PPrm2_1 <-> aa
        ! PPrm1_1 <-> a
        CALL Multiply(PPrm3_1,F      ,EPrm     )  !Eaaa=Paaa*F
        CALL Multiply(PPrm2_1,FPrm1_1,EPrm, One)  !Eaaa=Paa*Fa
        CALL Multiply(PPrm1_1,FPrm2_1,EPrm, One)  !Eaaa=Pa*Faa
        CALL Multiply(P      ,FPrm3_1,EPrm, One)  !Eaaa=P*Faaa
        !
        CALL Multiply(FPrm3_1,P      ,EPrm,-One)  !Eaaa=Faaa*P
        CALL Multiply(FPrm2_1,PPrm1_1,EPrm, One)  !Eaaa=Faa*Pa
        CALL Multiply(FPrm1_1,PPrm2_1,EPrm, One)  !Eaaa=Fa*Paa
        CALL Multiply(F      ,PPrm3_1,EPrm, One)  !Eaaa=F*Paaa
     ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
        ! A.EQ.B.NE.C
        ! PPrm3_1 <-> aac
        ! PPrm2_1 <-> aa
        ! PPrm2_2 <-> ac
        ! PPrm1_1 <-> a
        ! PPrm1_3 <-> c
        CALL Multiply(PPrm1_1,2.0d0/3.0d0)
        CALL Multiply(PPrm1_3,1.0d0/3.0d0)
        CALL Multiply(FPrm1_1,2.0d0/3.0d0)
        CALL Multiply(FPrm1_3,1.0d0/3.0d0)
        !
        CALL Multiply(PPrm3_1,F      ,EPrm     )  !Eaac=Paac*F
        CALL Multiply(PPrm2_1,FPrm1_3,EPrm, One)  !Eaac=Paa*Fc
        CALL Multiply(PPrm2_2,FPrm1_1,EPrm, One)  !Eaac=Pac*Fa
        CALL Multiply(PPrm1_1,FPrm2_2,EPrm, One)  !Eaac=Pa*Fac
        CALL Multiply(PPrm1_3,FPrm2_1,EPrm, One)  !Eaac=Pc*Faa
        CALL Multiply(P      ,FPrm3_1,EPrm, One)  !Eaac=P*Faac
        !
        CALL Multiply(FPrm3_1,P      ,EPrm,-One)  !Eaac=Faac*P
        CALL Multiply(FPrm2_1,PPrm1_3,EPrm, One)  !Eaac=Faa*Pc
        CALL Multiply(FPrm2_2,PPrm1_1,EPrm, One)  !Eaac=Fac*Pa
        CALL Multiply(FPrm1_1,PPrm2_2,EPrm, One)  !Eaac=Fa*Pac
        CALL Multiply(FPrm1_3,PPrm2_1,EPrm, One)  !Eaac=Fc*Paa
        CALL Multiply(F      ,PPrm3_1,EPrm, One)  !Eaac=F*Paac
     ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
        ! A.NE.B.EQ.C
        ! PPrm3_1 <-> abb
        ! PPrm2_1 <-> ab
        ! PPrm2_3 <-> bb
        ! PPrm1_1 <-> a
        ! PPrm1_2 <-> b
        CALL Multiply(PPrm1_1,1.0d0/3.0d0)
        CALL Multiply(PPrm1_2,2.0d0/3.0d0)
        CALL Multiply(FPrm1_1,1.0d0/3.0d0)
        CALL Multiply(FPrm1_2,2.0d0/3.0d0)
        !
        CALL Multiply(PPrm3_1,F      ,EPrm     )  !Eabb=Pabb*F
        CALL Multiply(PPrm2_1,FPrm1_2,EPrm, One)  !Eabb=Pab*Fb
        CALL Multiply(PPrm2_3,FPrm1_1,EPrm, One)  !Eabb=Pbb*Fa
        CALL Multiply(PPrm1_1,FPrm2_3,EPrm, One)  !Eabb=Pa*Fbb
        CALL Multiply(PPrm1_2,FPrm2_1,EPrm, One)  !Eabb=Pb*Fab
        CALL Multiply(P      ,FPrm3_1,EPrm, One)  !Eabb=P*Fabb
        !
        CALL Multiply(FPrm3_1,P      ,EPrm,-One)  !Eabb=Fabb*P
        CALL Multiply(FPrm2_1,PPrm1_2,EPrm, One)  !Eabb=Fab *Pb
        CALL Multiply(FPrm2_3,PPrm1_1,EPrm, One)  !Eabb=Fbb *Pa
        CALL Multiply(FPrm1_1,PPrm2_3,EPrm, One)  !Eabb=Fa  *Pbb
        CALL Multiply(FPrm1_2,PPrm2_1,EPrm, One)  !Eabb=Fb  *Pab
        CALL Multiply(F      ,PPrm3_1,EPrm, One)  !Eabb=F   *Pabb
     ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
        ! A.NE.B.NE.C
        ! PPrm3_1 <-> abc
        ! PPrm2_1 <-> ab
        ! PPrm2_2 <-> ac
        ! PPrm2_3 <-> bc
        ! PPrm1_1 <-> a
        ! PPrm1_2 <-> b
        ! PPrm1_3 <-> c
        CALL Multiply(PPrm1_1,1.0d0/3.0d0)
        CALL Multiply(PPrm1_2,1.0d0/3.0d0)
        CALL Multiply(PPrm1_3,1.0d0/3.0d0)
        CALL Multiply(FPrm1_1,1.0d0/3.0d0)
        CALL Multiply(FPrm1_2,1.0d0/3.0d0)
        CALL Multiply(FPrm1_3,1.0d0/3.0d0)
        !
        CALL Multiply(PPrm3_1,F      ,EPrm     )  !Eabc=Pabc*F
        CALL Multiply(PPrm2_1,FPrm1_3,EPrm, One)  !Eabc=Pab*Fc
        CALL Multiply(PPrm2_2,FPrm1_2,EPrm, One)  !Eabc=Pac*Fb
        CALL Multiply(PPrm2_3,FPrm1_1,EPrm, One)  !Eabc=Pbc*Fa
        CALL Multiply(PPrm1_1,FPrm2_3,EPrm, One)  !Eabc=Pa*Fbc
        CALL Multiply(PPrm1_2,FPrm2_2,EPrm, One)  !Eabc=Pb*Fac
        CALL Multiply(PPrm1_3,FPrm2_1,EPrm, One)  !Eabc=Pc*Fab
        CALL Multiply(P      ,FPrm3_1,EPrm, One)  !Eabc=P*Fabc
        !
        CALL Multiply(FPrm3_1,P      ,EPrm,-One)  !Eabc=Fabc*P
        CALL Multiply(FPrm2_1,PPrm1_3,EPrm, One)  !Eabc=Fab*Pc
        CALL Multiply(FPrm2_2,PPrm1_2,EPrm, One)  !Eabc=Fac*Pb
        CALL Multiply(FPrm2_3,PPrm1_1,EPrm, One)  !Eabc=Fbc*Pa
        CALL Multiply(FPrm1_1,PPrm2_3,EPrm, One)  !Eabc=Fa*Pbc
        CALL Multiply(FPrm1_2,PPrm2_2,EPrm, One)  !Eabc=Fb*Pac
        CALL Multiply(FPrm1_3,PPrm2_1,EPrm, One)  !Eabc=Fc*Pab
        CALL Multiply(F      ,PPrm3_1,EPrm, One)  !Eabc=F*Pabc
     ELSE
        CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
     ENDIF
     !
  CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
  END SELECT
  !
  ! Deallocate local arrays.
  CALL Delete(F)
  CALL Delete(P)
  SELECT CASE(RespOrder)
  CASE(1)
     CALL Delete(FPrm1_1)
     CALL Delete(PPrm1_1)
  CASE(2)
     CALL Delete(FPrm1_1)
     CALL Delete(PPrm1_1)
     CALL Delete(FPrm2_1)
     CALL Delete(PPrm2_1)
     IF(Chr1.NE.Chr2) THEN
        CALL Delete(FPrm1_2)
        CALL Delete(PPrm1_2)
     ENDIF
  CASE(3)
     CALL Delete(FPrm1_1)
     CALL Delete(PPrm1_1)
     CALL Delete(FPrm2_1)
     CALL Delete(PPrm2_1)
     CALL Delete(FPrm3_1)
     CALL Delete(PPrm3_1)
     IF(Chr1.EQ.Chr2.AND.Chr1.EQ.Chr3) THEN
     ELSEIF(Chr1.EQ.Chr2.AND.Chr1.NE.Chr3) THEN
        CALL Delete(FPrm1_3)
        CALL Delete(PPrm1_3)
        CALL Delete(FPrm2_2)
        CALL Delete(PPrm2_2)
     ELSEIF(Chr1.NE.Chr2.AND.Chr2.EQ.Chr3) THEN
        CALL Delete(FPrm1_2)
        CALL Delete(PPrm1_2)
        CALL Delete(FPrm2_3)
        CALL Delete(PPrm2_3)
     ELSEIF(Chr1.NE.Chr2.AND.Chr1.NE.Chr3.AND.Chr2.NE.Chr3) THEN
        CALL Delete(FPrm1_2)
        CALL Delete(PPrm1_2)
        CALL Delete(FPrm1_3)
        CALL Delete(PPrm1_3)
        CALL Delete(FPrm2_2)
        CALL Delete(PPrm2_2)
        CALL Delete(FPrm2_3)
        CALL Delete(PPrm2_3)
     ELSE
        CALL Halt('Response: unknown symmetry <'//Chr1//Chr2//Chr3//'>.')
     ENDIF
  CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
  END SELECT
  !
  ! We dont filter E' for obvious reasons .
  CALL Put(EPrm,TrixFile('EPrime'//TRIM(Args%C%C(3)),Args,0))
  !
#ifdef DDIIS_DBUG
  if(myid.eq.0) WRITE(*,*) 'Save E''=<'//TRIM(TrixFile('EPrime'//TRIM(Args%C%C(3)),Args,0))//'>'
#endif
  !
  ! Compute the DDIIS error.
  DIISErr=SQRT(Dot(EPrm,EPrm))/DBLE(NBasF)
  !write(*,*) 'DIISErr',DIISErr
  !
  ! IO save DDIIS error.
  CALL Put(DIISErr,'ddiiserr')
  !
  !-------------------------------------------------------------------
  ! Build the B matrix and solve the linear problem.
  !-------------------------------------------------------------------
  !
  ! Build the B matrix if on second SCF cycle (starting from 0)
  ! and if pure damping flag is not on (DIISDimension=0)
  SELECT CASE(DoDIIS)
  CASE(1)
     ! Do DDIIS
     !
     DDIISCurDim=DDIISEnd-DDIISBeg+1
     N=DDIISCurDim+1
     M=DDIISBeg
     CALL New(B,(/N,N/))
     !
     ! Build the B matrix.
     I0=DDIISBeg-CPSCFCycl   ! -1!I added that.
     DO I=1,N-1
#ifdef DDIIS_DEBUG
        WRITE(*,*) 'Load <ei|=<'//TRIM(TrixFile('EPrime'//TRIM(Args%C%C(3)),Args,I0))//'>'
#endif
        CALL Get(Tmp1,TrixFile('EPrime'//TRIM(Args%C%C(3)),Args,I0))
        J0=I0
        DO J=I,N-1
#ifdef DDIIS_DEBUG
           WRITE(*,*) 'Load |ej>=<'//TRIM(TrixFile('EPrime'//TRIM(Args%C%C(3)),Args,J0))//'>'
#endif
           CALL Get(EPrm,TrixFile('EPrime'//TRIM(Args%C%C(3)),Args,J0))
           B%D(I,J)=Dot(Tmp1,EPrm)
           B%D(J,I)=B%D(I,J)
           J0=J0+1
        ENDDO
        I0=I0+1
     ENDDO
     B%D(N,1:N)=One
     B%D(1:N,N)=One
     B%D(N,N)=Zero
     !
     ! Solve the least squares problem to obtain new DIIS coeficients.
     CALL New(DIISCo,N)
     CALL SetEq(DIISCo,Zero)
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
     !
     CALL New(BInv,(/N,N/))
     CALL New(V,N)
     !
     ! Solve the linear system and remove
     ! linear dependancies if needed.
     DO
        V%D=Zero
        V%D(N)=One
        CALL SetDSYEVWork(N)
        BInv%D=Zero
        IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                            EigenThresh_O=EigThresh,PrintCond_O=.TRUE.,Prog_O=Prog, &
                            CoNo_O=CondA)
        ELSE
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                            EigenThresh_O=EigThresh, &
                            CoNo_O=CondA)
        ENDIF
        CALL UnSetDSYEVWork()
        CALL DGEMV('N',N,N,One,BInv%D(1,1),N,V%D(1),1,Zero,DIISCo%D(1),1)
        !
#ifdef DDIIS_INFO
        !CALL PrintMatrix(B%D,N,N,2)
        WRITE(*,*) 'CondA=',CondA
#endif
        !
        ! Do we need to reduce the size?
        IF(.NOT.(CondA.GT.DDIISMaxCond.AND.DDIISCurDim.GT.1)) EXIT
        !
        ! We need to remove the oldest element.
        DDIISBeg=DDIISBeg+1
        !
        ! New current dimension.
        DDIISCurDim=DDIISEnd-DDIISBeg+1
        !
        ! Set the new size.
        N=DDIISCurDim+1
        M=DDIISBeg
        !
        ! Copy temporarly the B matrix.
        BInv%D=B%D
        CALL Delete(B)
        CALL New(B,(/N,N/))
        B%D=Zero
        !
        ! Copy back the B matrix.
        B%D(1:N,1:N)=BInv%D(2:N+1,2:N+1)
        !
        ! Reallocate and set old arrays.
        CALL Delete(DIISCo)
        CALL New(DIISCo,N)
        DIISCo%D=Zero
        CALL Delete(BInv)
        CALL New(BInv,(/N,N/))
        BInv%D=Zero
#ifdef DDIIS_INFO
        WRITE(*,*) 'Reduce DDIIS Space from '//TRIM(IntToChar(DDIISCurDim+1)) &
                   //' to '//TRIM(IntToChar(DDIISCurDim))//'.'
#endif
     ENDDO
     !
     CALL Delete(B   )
     CALL Delete(V   )
     CALL Delete(BInv)
     !
     !BTmp%D=Zero
     !BTmp%D(1:N,1:N)=B%D(1:N,1:N)
     !CALL PrintMatrix(BTmp%D,BMax+1,BMax+1,2)
     !
#ifdef PARALLEL
     ENDIF
     CALL BCast(DIISCo)
     CALL BCast(N)
#endif
     Mssg=ProcessName(Prog,'Pulay C1')//'DIISCo = '
  CASE DEFAULT
     ! Do Damping.
     !
     Mssg=ProcessName(Prog,'Damping')//'Co = '
     N=3
     M=CPSCFCycl-N+2
     CALL New(DIISCo,2)
     !
     ! Damping on the second cycle
     DIISCo%D(1)=One-Damp
     DIISCo%D(2)=Damp
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Printing.
  !-------------------------------------------------------------------
  !
  IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
     DO I=1,N-2
        IF(MOD(I,4)==0)THEN
           Mssg=TRIM(Mssg)//RTRN//ProcessName() &
                  //'          '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
        ELSE
           Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
        ENDIF
     ENDDO
     Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(DIISCo%D(N-1)))
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(Out)
        ENDIF
#ifdef PARALLEL
     ENDIF
#endif
  ENDIF
  !
  !-------------------------------------------------------------------
  ! Reordering.
  !-------------------------------------------------------------------
  !
  ! Allocate some indecies for re-ordering
  CALL New(Idx      ,N)
  CALL New(SCFOff   ,N)
  CALL New(AbsDIISCo,N)
  !
  ! Reorder the DIIS, starting with smallest values and summing to the largest
  IOffSet=M-CPSCFCycl
  DO I=1,N-1
     Idx%I(I)=I
     SCFOff%I(I)=IOffSet
     AbsDIISCo%D(I)=ABS(DIISCo%D(I))
     IOffSet=IOffSet+1
  ENDDO
  CALL Sort(AbsDIISCo,Idx,N-1,1)
  !
  ! Start with a matrix of diagonal zeros...
  CALL New(TmpFPrm)
  CALL SetToI(TmpFPrm)
  CALL Multiply(TmpFPrm,Zero)
  !
  ! And do the summation
  DO I=1,N-1
     CALL Get(Tmp1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(3)),Args,SCFOff%I(Idx%I(I))))
     CALL Multiply(Tmp1,DIISCo%D(Idx%I(I)))
     CALL Add(TmpFPrm,Tmp1,EPrm)
     CALL SetEq(TmpFPrm,EPrm)
  ENDDO
  !
  ! Deallocate local arrays.
  CALL Delete(Idx      )
  CALL Delete(SCFOff   )
  CALL Delete(AbsDIISCo)
  !
  !-------------------------------------------------------------------
  ! Increment DDIIS and IO.
  !-------------------------------------------------------------------
  !
  DDIISEnd=DDIISEnd+1
  DDIISCurDim=DDIISEnd-DDIISBeg+1
  IF(DDIISCurDim.GT.BMax) DDIISBeg=DDIISBeg+1
  !
  ! Put in HDF Beg and End DDIIS variables.
  CALL Put(DDIISBeg,'DDIISBeg'//IntToChar(CPSCFCycl+1))
  CALL Put(DDIISEnd,'DDIISEnd'//IntToChar(CPSCFCycl+1))
  !CALL Put(BTmp,'DDIISBMtrix')
  !
  !-------------------------------------------------------------------
  ! IO for the orthogonal, extrapolated FPrim
  !-------------------------------------------------------------------
  !
  CALL Put(TmpFPrm,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(3)),Args,0))
  CALL PChkSum(TmpFPrm,'FPrime_DDIIS'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( TmpFPrm,'FPrime_DDIIS'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
  CALL Plot(   TmpFPrm,'FPrime_DDIIS'//TRIM(Args%C%C(3))//'_'//TRIM(SCFCycl))
  !
  !-------------------------------------------------------------------
  ! Tidy up
  !-------------------------------------------------------------------
  !
  CALL Delete(Tmp1     )
  CALL Delete(TmpFPrm  )
  CALL Delete(EPrm     )
  CALL Delete(DIISCo   )
  !
  !CALL Delete(BTmp)
  !
  CALL ShutDown(Prog)
  !
CONTAINS
  !
  SUBROUTINE PrintMatrix(V,M,N,IOpt,IOut_O,SHFTM_O,SHFTN_O,TEXT_O)
    IMPLICIT NONE
!    ++++ PRINT OUT A RECTANGULAR MATRIX +++++
!
!    V    : MATRIX MxN
!    M    : NUMBER OF ROW
!    N    : NUMBER OF COLUMN
!
!    IOpt = 0; MAX COLUMN = 10
!    IOpt = 1; MAX COLUMN =  7
!    IOpt = 2; MAX COLUMN =  5
!    IOpt = 3; MAX COLUMN =  3
!
    INTEGER                                  :: IMax,IMin,NMAX,I,J
    INTEGER                                  :: IOut,SHFTM,SHFTN
    INTEGER                     , INTENT(IN) :: M,N,IOpt
    INTEGER         , OPTIONAL  , INTENT(IN) :: IOut_O,SHFTM_O,SHFTN_O
    REAL(DOUBLE), DIMENSION(:,:), INTENT(IN) :: V
    CHARACTER(LEN=*), OPTIONAL  , INTENT(IN) :: TEXT_O
!
    IOut = 6
    IF(PRESENT(IOut_O)) IOut = IOut_O
!
    SHFTM = 0
    SHFTN = 0
    IF(PRESENT(SHFTM_O)) SHFTM = SHFTM_O
    IF(PRESENT(SHFTN_O)) SHFTN = SHFTN_O
!
    WRITE(IOut,100)
!
    SELECT CASE(IOPT)
    CASE(0); NMAX = 10
    CASE(1); NMAX =  7
    CASE(2); NMAX =  5
    CASE(3); NMAX =  3
    CASE DEFAULT
       ! WRITE ERROR MESSAGE OR PUT SOMETHING THERE
       RETURN
    END SELECT
!
    IF(N.EQ.0) RETURN
!
    IF(PRESENT(TEXT_O)) WRITE(IOUT,*) TEXT_O
    IMAX = 0
!
    DO WHILE (IMAX.LT.M)
       IMIN = IMAX+1
       IMAX = IMAX+NMAX
       IF(IMAX .GT. M) IMAX = M
       SELECT CASE(IOPT)
       CASE(0); WRITE(IOUT,1000) (I+SHFTN,I=IMIN,IMAX)
       CASE(1); WRITE(IOUT,2000) (I+SHFTN,I=IMIN,IMAX)
       CASE(2); WRITE(IOUT,3000) (I+SHFTN,I=IMIN,IMAX)
       CASE(3); WRITE(IOUT,4000) (I+SHFTN,I=IMIN,IMAX)
       END SELECT
       DO J = 1,N
          SELECT CASE(IOPT)
          CASE(0); WRITE(IOUT,1100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(1); WRITE(IOUT,2100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(2); WRITE(IOUT,3100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(3); WRITE(IOUT,4100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          END SELECT
       ENDDO
       WRITE (IOUT,100)
    ENDDO
    !
    WRITE (IOUT,100)
    !
100 FORMAT('')
1000 FORMAT(6X,10(4X,I3,4X))
1100 FORMAT(I5,1X,10F11.5  )
2000 FORMAT(6X,7(6X,I3,6X ))
2100 FORMAT(I5,1X,7F15.10  )
3000 FORMAT(6X,7(7X,I3,6X ))
3100 FORMAT(I5,1X,7E16.8   )
4000 FORMAT(6X,7(7X,I3,6X ))
4100 FORMAT(I5,1X,7E16.8   )
  END SUBROUTINE PrintMatrix

END PROGRAM DDIIS
