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
!  MODULE FOR SPARSE BLOCKED LINEAR ALGEBRA
!  Author:  Matt Challacombe
!-------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE COMMON_DEBUG
  USE DerivedTypes
  IMPLICIT NONE
  INTEGER :: JPrc,FFrom
  REAL(DOUBLE) :: T1
  TYPE(INT_VECT) :: Pt,NActual,NPredct
  TYPE(DBL_VECT) :: Tim
END MODULE COMMON_DEBUG

MODULE LinAlg
  USE GlobalCharacters
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE PrettyPrint
  USE MemMan
  USE SetXYZ
  USE InOut
#ifdef NAG
  USE F90_UNIX_PROC
#endif
  USE Order
  USE Clock
#ifdef PARALLEL
  USE MondoMPI
#endif
  USE Thresholding
  IMPLICIT NONE
  !===============================================================================
  !  Interface blocks for generic linear algebra routines
  !===============================================================================
  INTERFACE Multiply
#ifdef PARALLEL
    MODULE PROCEDURE MultiplyM_BCSR, MultiplyM_DBCSR,           &
         MultiplyM_BCSR_SCLR, MultiplyM_DBCSR_SCLR, &
         !MultiplyM_BCSR_VECT, MultiplyM_DBCSR_VECT
         MultiplyM_BCSR_VECT
#else
    MODULE PROCEDURE MultiplyM_BCSR, MultiplyM_BCSR_SCLR,       &
         MultiplyM_BCSR_VECT
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE SymbolikMM
#ifdef PARALLEL
    MODULE PROCEDURE SymbolikMM_BCSR, SymbolikMM_DBCSR
#else
    MODULE PROCEDURE SymbolikMM_BCSR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE NumerikMM
#ifdef PARALLEL
    MODULE PROCEDURE NumerikMM_BCSR, NumerikMM_DBCSR
#else
    MODULE PROCEDURE NumerikMM_BCSR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE Add
#ifdef PARALLEL
    MODULE PROCEDURE Add_BCSR, Add_DBCSR,  &
         Add_BCSR_SCLR, Add_DBCSR_SCLR
#else
    MODULE PROCEDURE Add_BCSR, Add_BCSR_SCLR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE NumerikADD
#ifdef PARALLEL
    MODULE PROCEDURE NumerikADD_BCSR, NumerikADD_DBCSR
#else
    MODULE PROCEDURE NumerikADD_BCSR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE Trace
#ifdef PARALLEL
    MODULE PROCEDURE TraceMM_BCSR, TraceMM_DBCSR, &
         Trace_BCSR, Trace_DBCSR
#else
    MODULE PROCEDURE TraceMM_BCSR, Trace_BCSR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE Dot
#ifdef PARALLEL
    MODULE PROCEDURE Dot_BCSR, Dot_DBCSR, Dot_DBL_VECT
#else
    MODULE PROCEDURE Dot_BCSR, Dot_DBL_VECT
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE Filter
#ifdef PARALLEL
    MODULE PROCEDURE FilterM_BCSR,FilterM_DBCSR, &
         FilterM_InPlace_BCSR,FilterM_InPlace_DBCSR
#else
    MODULE PROCEDURE FilterM_BCSR,FilterM_InPlace_BCSR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE XPose
    MODULE PROCEDURE XPose_BCSR,XPose_Simple
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE Max
#ifdef PARALLEL
    MODULE PROCEDURE Max_BCSR, Max_DBCSR
#else
    MODULE PROCEDURE Max_BCSR
#endif
  END INTERFACE
  !-------------------------------------------------------------------------------
  INTERFACE Plot
    MODULE PROCEDURE Plot_BCSR
#ifdef PARALLEL
    MODULE PROCEDURE Plot_DBCSR
#endif
  END INTERFACE

  INTERFACE SetToI
    MODULE PROCEDURE SetToI_BCSR
#ifdef PARALLEL
    MODULE PROCEDURE SetToI_DBCSR
#endif
  END INTERFACE

  INTERFACE FNorm
    MODULE PROCEDURE FNorm_BCSR
#ifdef PARALLEL
    MODULE PROCEDURE FNorm_DBCSR
#endif
  END INTERFACE

  !  Global stuff for matrix algebra
  TYPE(INT_VECT) :: GlobalRowPtA,GlobalRowPtB

  !===============================================================================
CONTAINS
  !  MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY
  !  PLY MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTI
  !  LTIPLY MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MULTIPLY, MATRIX MU
#ifdef PARALLEL
  !===============================================================================
  !     F90 wrapper for F77 BCSR numeric matrix multiply: C=C+Beta*A.B
  !===============================================================================
  SUBROUTINE SymbolikMM_DBCSR(A,B,C,Flag,UpDate)
    IMPLICIT NONE
    TYPE(DBCSR), INTENT(IN)         :: A,B
    LOGICAL,     INTENT(IN)         :: UpDate
    TYPE(DBCSR), INTENT(INOUT)      :: C
    TYPE(INT_VECT), INTENT(INOUT)   :: Flag
    INTEGER                         :: CNNon0_Old
    TYPE(INT_VECT)                  :: DRowPt,DColPt,DBlkPt
    INTEGER                         :: Status,DNAtms,DNBlks
    INTEGER, EXTERNAL               :: SymbolikMM_GENERIC77
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: SzCpt,SzMat

    CALL Initialize(DRowPt)
    CALL Initialize(DColPt)
    CALL Initialize(DBlkPt)

    IF(UpDate)THEN
      DNAtms=C%NAtms
      DNBlks=C%NBlks
      CALL New(DRowPt,C%NAtms+1)
      CALL New(DColPt,C%NBlks)
      CALL New(DBlkPt,C%NBlks)
      CNNon0_Old=C%NNon0
    ELSE
      DNAtms=1
      DNBlks=1
      CALL New(DRowPt,1+1)
      CALL New(DColPt,1)
      CALL New(DBlkPt,1)
    ENDIF
    CALL Load(A,GlobalRowPtA)
    CALL Load(B,GlobalRowPtB)
    C%NAtms=A%NAtms
    Status=SymbolikMM_GENERIC_77(C%NSMat,               & !<<<SPIN
         NAtoms,A%NAtms,A%NBlks,NAtoms,                 &
         B%NAtms,B%NBlks,NAtoms,                        &
         DNAtms,DNBlks,SIZE(C%ColPt%I),SIZE(C%MTrix%D), &
         C%NAtms,C%NBlks,C%NNon0,UpDate,OffSt%I(MyID),  &
         GlobalRowPtA%I,A%ColPt%I,                      &
         GlobalRowPtB%I,B%ColPt%I,                      &
         C%RowPt%I,C%ColPt%I,C%BlkPt%I,                 &
         DRowPt%I,  DColPt%I, DBlkPt%I,                 &
         BSiz%I,Flag%I)
    IF(Status==FAIL)THEN
      SzCpt='  SIZE(C%ColPt) = '//IntToChar(SIZE(C%ColPt%I))
      SzMat='  SIZE(C%MTrix) = '//IntToChar(SIZE(C%MTrix%D))
      CALL Halt(' Dimensions in SymbolikMM_DBCSR:' &
           //Rtrn//SzCpt//Rtrn//SzMat)
    ENDIF
    CALL Clear(A,GlobalRowPtA)
    CALL Clear(B,GlobalRowPtB)
    CALL Delete(DRowPt)
    CALL Delete(DColPt)
    CALL Delete(DBlkPt)
    IF(UpDate)THEN
      C%MTrix%D(CNNon0_Old+1:C%NNon0)=Zero
    ELSE
      C%MTrix%D(1:C%NNon0)=Zero
    ENDIF
  END SUBROUTINE SymbolikMM_DBCSR
#endif
  !===============================================================================
  !     BCSR wrapper for generic BCSR symbolic matrix multiply: C=Beta*C+A.B
  !===============================================================================
  SUBROUTINE SymbolikMM_BCSR(A,B,C,Flag,UpDate)
    IMPLICIT NONE
    TYPE(BCSR),  INTENT(IN)       :: A,B
    LOGICAL,     INTENT(IN)       :: UpDate
    TYPE(BCSR),  INTENT(INOUT)    :: C
    TYPE(INT_VECT), INTENT(INOUT) :: Flag
    INTEGER                       :: CNNon0_Old
    TYPE(INT_VECT)                :: DRowPt,DColPt,DBlkPt
    INTEGER, EXTERNAL             :: SymbolikMM_GENERIC77
    INTEGER                       :: Status,DNAtms,DNBlks

    CALL Initialize(DRowPt)
    CALL Initialize(DColPt)
    CALL Initialize(DBlkPt)

    IF(UpDate)THEN
      DNAtms=C%NAtms
      DNBlks=C%NBlks
      CALL New(DRowPt,C%NAtms+1)
      CALL New(DColPt,C%NBlks)
      CALL New(DBlkPt,C%NBlks)
      CNNon0_Old=C%NNon0
    ELSE
      DNAtms=1
      DNBlks=1
      CALL New(DRowPt,1)
      CALL New(DColPt,1)
      CALL New(DBlkPt,1)
    ENDIF

    Status=SymbolikMM_GENERIC_77(C%NSMat,               & !<<<SPIN
         NAtoms,A%NAtms,A%NBlks,A%NAtms,                &
         B%NAtms,B%NBlks,B%NAtms,                       &
         DNAtms,DNBlks,SIZE(C%ColPt%I),SIZE(C%MTrix%D), &
         C%NAtms,C%NBlks,C%NNon0,UpDate,0,              &
         A%RowPt%I,A%ColPt%I,                           &
         B%RowPt%I,B%ColPt%I,                           &
         C%RowPt%I,C%ColPt%I,C%BlkPt%I,                 &
         DRowPt%I,  DColPt%I, DBlkPt%I,                 &
         BSiz%I,Flag%I)
    IF(Status==FAIL) THEN
      CALL Halt('Dimensions in SymbolikMM_BCSR')
    ENDIF

    CALL Delete(DRowPt)
    CALL Delete(DColPt)
    CALL Delete(DBlkPt)

    IF(UpDate)THEN
      C%MTrix%D(CNNon0_Old+1:C%NNon0)=Zero
    ELSE
      C%MTrix%D(1:C%NNon0)=Zero
    ENDIF
  END SUBROUTINE SymbolikMM_BCSR

  !C=========================================================================
  !C     Generic F77 style (D,B)CSR symbolic matrix multiply: C=C+Beta*A.B
  !C=========================================================================
  !C23456789012345678901234567890123456789012345678901234567890123456789012
  INTEGER FUNCTION SymbolikMM_GENERIC_77(NSMat,NAtoms77,               & !<<<SPIN
       ANAtms,ANBlks,ANAtoms,        &
       BNAtms,BNBlks,BNAtoms,        &
       DNAtms,DNBlks,MxBlks,MxNon0,  &
       CNAtms,CNBlks,CNNon0,UpDate,  &
       AOffSt,ARowPt,AColPt,         &
       BRowPt,BColPt,                &
       CRowPt,CColPt,CBlkPt,         &
       DRowPt,DColPt,DBlkPt,         &
       BSiz77,Flag77)
    IMPLICIT NONE
    INTEGER :: NSMat,                                                & !<<<SPIN
         NAtoms77,ANAtms,ANBlks,ANAtoms,BNAtms,BNBlks,BNAtoms, &
         DNAtms,DNBlks,MxBlks,MxNon0,CNAtms,CNBlks,CNNon0,     &
         AOffSt,JP,KP,S,T,IL,IG,JG,KG,MA,MB,MM,                &
         IStrtA,IStopA,IStrtB,IStopB,IStrtD,IStopD
    INTEGER, DIMENSION(ANAtoms+1) :: ARowPt
    INTEGER, DIMENSION(ANBlks)    :: AColPt
    INTEGER, DIMENSION(BNAtoms+1) :: BRowPt
    INTEGER, DIMENSION(BNBlks)    :: BColPt
    INTEGER, DIMENSION(CNAtms+1)  :: CRowPt
    INTEGER, DIMENSION(MxBlks)    :: CColPt
    INTEGER, DIMENSION(MxBlks)    :: CBlkPt
    INTEGER, DIMENSION(DNAtms+1)  :: DRowPt
    INTEGER, DIMENSION(DNBlks)    :: DColPt
    INTEGER, DIMENSION(DNBlks)    :: DBlkPt
    INTEGER, DIMENSION(NAtoms77)  :: BSiz77
    INTEGER, DIMENSION(NAtoms77)  :: Flag77
    LOGICAL UpDate
    !C-----------------------------------------------------------------------------------
    T=1
    IF(UpDate)THEN
      S=CNNon0+1
      CALL INT_VECT_EQ_INT_VECT(CNAtms+1,DRowPt,CRowPt)
      CALL INT_VECT_EQ_INT_VECT(CNBlks  ,DColPt,CColPt)
      CALL INT_VECT_EQ_INT_VECT(CNBlks  ,DBlkPt,CBlkPt)
    ELSE
      S=1
      CNAtms=ANAtms
    ENDIF
    DO IL=1,ANAtms
      IG=IL+AOffSt
      CRowPt(IL)=T
      MA=BSiz77(IG)
      IStrtA=ARowPt(IG)
      IStopA=ARowPt(IG+1)-1
      IF(UpDate)THEN
        IStrtD=DRowPt(IL)
        IStopD=DRowPt(IL+1)-1
        DO JP=IStrtD,IStopD
          JG=DColPt(JP)
          IF(Flag77(JG).EQ.0)THEN
            Flag77(JG)=1
            CColPt(T)=JG
            CBlkPt(T)=DBlkPt(JP)
            T=T+1
          ENDIF
        ENDDO
      ENDIF
      DO JP=IStrtA,IStopA
        JG=AColPt(JP)
        IStrtB=BRowPt(JG)
        IStopB=BRowPt(JG+1)-1
        IF(IStrtB.NE.0.AND.IStopB.NE.0)THEN
          MB=BSiz77(JG)
          MM=MA*MB
          DO KP=IStrtB,IStopB
            KG=BColPt(KP)
            IF(Flag77(KG).EQ.0)THEN
              CColPt(T)=KG
              CBlkPt(T)=S
              Flag77(KG) =1
              T=T+1
              S=S+MA*BSiz77(KG)*NSMat !<<<SPIN
              IF(T.GT.MxBlks.OR.S.GT.MxNon0)THEN
                SymbolikMM_Generic_77=-1
                RETURN
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DO KP=CRowPt(IL),T-1
        Flag77(CColPt(KP))=0
      ENDDO
    ENDDO
    CNNon0=S-1
    CNBlks=T-1
    CRowPt(CNAtms+1)=T
    SymbolikMM_Generic_77=0
    DO IG=1,NAtoms77
      Flag77(IG)=0
    ENDDO
    RETURN
  END FUNCTION SymbolikMM_Generic_77
#ifdef PARALLEL
  !===============================================================================
  !     Wrapper for generic F77 style DBCSR numeric matrix multiply: C=C+A.B
  !===============================================================================
  SUBROUTINE NumerikMM_DBCSR(A,B,BMTrix,C,ASMat,BSMat,CSMat,Flag,Perf_O) !<<<SPIN
    IMPLICIT NONE
    TYPE(DBCSR),              INTENT(IN)    :: A,B
    REAL(DOUBLE),DIMENSION(:),INTENT(IN)    :: BMTrix
    TYPE(DBCSR),              INTENT(INOUT) :: C
    INTEGER, INTENT(IN)                     :: ASMat,BSMat,CSMat !<<<SPIN
    TYPE(INT_VECT), INTENT(INOUT)           :: Flag
    TYPE(TIME),OPTIONAL,      INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                            :: FlOp

    FlOp=Zero
    CALL Load(A,GlobalRowPtA)
    CALL Load(B,GlobalRowPtB)
    CALL NumerikMM_GENERIC(ASMat,BSMat,CSMat,A%NAtms,OffSt%I(MyId),      & !<<<SPIN
         GlobalRowPtA%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         GlobalRowPtB%I,B%ColPt%I,B%BlkPt%I,BMTrix,    &
         C%RowPt%I,C%ColPt%I,C%BlkPt%I,C%MTrix%D,      &
         BSiz%I,Flag%I,FlOp)
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O)) &
         Perf_O%FLOP=Perf_O%FLOP+FlOp
    CALL Clear(A,GlobalRowPtA)
    CALL Clear(B,GlobalRowPtB)
  END SUBROUTINE NumerikMM_DBCSR
#endif

  !===============================================================================
  !     Wrapper for generic F77 style BCSR numeric matrix multiply: C=C+A.B
  !===============================================================================
  SUBROUTINE NumerikMM_BCSR(A,B,C,ASMat,BSMat,CSMat,Flag,Perf_O) !<<<SPIN
    IMPLICIT NONE
    TYPE(BCSR),         INTENT(IN)    :: A,B
    TYPE(BCSR),         INTENT(INOUT) :: C
    INTEGER, INTENT(IN)               :: ASMat,BSMat,CSMat !<<<SPIN
    TYPE(INT_VECT), INTENT(INOUT)     :: Flag
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: FlOp

    FlOp=Zero
    CALL NumerikMM_GENERIC(ASMat,BSMat,CSMat,A%NAtms,0,             & !<<<SPIN
         A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
         C%RowPt%I,C%ColPt%I,C%BlkPt%I,C%MTrix%D, &
         BSiz%I,Flag%I,FlOp)
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE NumerikMM_BCSR

  !===============================================================================
  !     Generic F77 style (D,B)CSR numeric matrix multiply: C=C+A.B
  !===============================================================================
  SUBROUTINE NumerikMM_GENERIC(ASMat,BSMat,CSMat, & !<<<SPIN
       ANAtms,AOffSt,                             &
       ARowPt,AColPt,ABlkPt,AMTrix,               &
       BRowPt,BColPt,BBlkPt,BMTrix,               &
       CRowPt,CColPt,CBlkPt,CMTrix,               &
       BSiz,Flag,FlOp)

    IMPLICIT NONE
    INTEGER,                  INTENT(IN)    :: ASMat,BSMat,CSMat,ANAtms,AOffSt !<<<SPIN
    INTEGER,     DIMENSION(:),INTENT(IN)    :: ARowPt,AColPt,ABlkPt, &
                                               BRowPt,BColPt,BBlkPt,BSiz
    REAL(DOUBLE),DIMENSION(:),INTENT(IN)    :: AMTrix,BMTrix
    INTEGER,     DIMENSION(:),INTENT(IN)    :: CRowPt,CColPt,CBlkPt
    REAL(DOUBLE),DIMENSION(:),INTENT(INOUT) :: CMTrix
    INTEGER,     DIMENSION(:),INTENT(INOUT) :: Flag
    REAL(DOUBLE),             INTENT(INOUT) :: FlOp
    REAL(DOUBLE) :: Op
    INTEGER      :: K,JP,KP,P,Q,R,IL,IG,JG,KG, &
                    MA,MB,NB,MNA,IStrtA,IStopA, &
                    IStrtB,IStopB,IStrtC,IStopC

    Op=Zero
    DO IL=1,ANAtms
      IG=IL+AOffSt
      MA=BSiz(IG)
      IStrtA=ARowPt(IG)
      IStopA=ARowPt(IG+1)-1
      IStrtC=CRowPt(IL)
      IStopC=CRowPt(IL+1)-1
      DO JP=IStrtC,IStopC
        Flag(CColPt(JP))=CBlkPt(JP)
      ENDDO
      DO JP=IStrtA,IStopA
        JG=AColPt(JP)
        IStrtB=BRowPt(JG)
        IStopB=BRowPt(JG+1)-1
        IF(IStrtB/=0.AND.IStopB/=0)THEN
          MB=BSiz(JG)
          MNA=MA*MB
          P=ABlkPt(JP)+(ASMat-1)*MNA !<<<SPIN
          DO KP=IStrtB,IStopB
            KG=BColPt(KP)
            NB=BSiz(KG)
            Q=BBlkPt(KP)+(BSMat-1)*MB*NB !<<<SPIN
            R=Flag(KG)+(CSMat-1)*MA*NB !<<<SPIN
            Op=Op+DBLE(MNA*NB)
            CALL DGEMM_NN(MA,MB,NB,One,AMTrix(P),BMTrix(Q),CMTrix(R))
            !                     CALL DGEMM_NN2(MA,MB,NB,One,AMTrix(P),BMTrix(Q),CMTrix(R))
          ENDDO
        ENDIF
      ENDDO
      DO K=IStrtC,IStopC
        Flag(CColPt(K))=0
      ENDDO
    ENDDO
    FlOp=FlOp+Two*Op
  END SUBROUTINE NumerikMM_GENERIC

  !===============================================================================
  !     Double GEneral Matrix Multiply: C_MxN=C_MxN+Beta*(A_MxK).(B_KxN)
  !===============================================================================
  SUBROUTINE DGEMM_NN2(M,K,N,Beta,A,B,C)
    IMPLICIT NONE
    INTEGER,                   INTENT(IN)  :: M,K,N
    REAL(DOUBLE),              INTENT(IN)  :: Beta
    REAL(DOUBLE),DIMENSION(:) ,INTENT(IN)  :: A,B
    REAL(DOUBLE),DIMENSION(:) ,INTENT(OUT) :: C
    REAL(DOUBLE)                           :: BB
    INTEGER                                :: I,J,J1,L,IMJ, &
                                              K_J_Off,M_J_Off,M_L_Off

    DO J=1,N
      J1=J-1
      K_J_Off=J1*K
      M_J_Off=J1*M
      DO L=1,K
        M_L_Off=(L-1)*M
        BB=Beta*B(L+K_J_Off)! BB=Beta*B(L,J)
        DO I=1,M
          IMJ=I+M_J_Off
          C(IMJ)=C(IMJ)+BB*A(I+M_L_Off)  ! C(I,J)=C(I,J)+BB*A(I,L)
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE DGEMM_NN2

  !===============================================================================
  !     Matrix multiply for BCSR matrices
  !===============================================================================
  SUBROUTINE MultiplyM_BCSR(A,B,C,Beta_O,Perf_O)
    TYPE(BCSR),           INTENT(IN)    :: A,B
    REAL(DOUBLE),OPTIONAL,INTENT(IN)    :: Beta_O
    TYPE(BCSR),           INTENT(INOUT) :: C
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    LOGICAL                             :: UpDate
    INTEGER                             :: NSMat,N1,N2
    TYPE(INT_VECT)                      :: Flag

#ifdef CHECK_UNRESTRICTED

    TYPE(DBL_RNK2)                      :: AD,BD,CD,DD

    CALL Initialize(AD)
    CALL Initialize(BD)
    CALL Initialize(CD)
    CALL Initialize(DD)
#endif
    CALL Initialize(Flag)

    NSMat=MAX(A%NSMat,B%NSMat)
    IF(PRESENT(Beta_O))THEN
      IF(Beta_O/=Zero)THEN

        ! Check allocation
        IF(.NOT.AllocQ(C%Alloc)) THEN
          CALL Halt('Asking for update of non-allocated BCSR matrix in MultiplyM_BCSR')
        ENDIF

        ! If Beta /=1 then rescale
        IF(Beta_O/=One)THEN
          C%MTrix%D(1:C%NNon0)=Beta_O*C%MTrix%D(1:C%NNon0)
          PerfMon%FLOP=PerfMon%FLOP+DBLE(C%NNon0)
          IF(PRESENT(Perf_O)) THEN
            Perf_O%FLOP=Perf_O%FLOP+DBLE(C%NNon0)
          ENDIF
        ENDIF
        UpDate=.TRUE.
      ELSE
        IF(.NOT.AllocQ(C%Alloc)) THEN
          CALL New(C,NSMat_O=NSMat)
        ENDIF
        UpDate=.FALSE.
      ENDIF
    ELSE
      IF(.NOT.AllocQ(C%Alloc)) THEN
        CALL New(C,NSMat_O=NSMat)
      ENDIF
      UpDate=.FALSE.
    ENDIF

    IF(C%NSMat.NE.NSMat) THEN
      CALL MondoLog(DEBUG_MAXIMUM, "MultiplyM_BCSR", "C%NSMat.NE.NSMat! Deallocate and reallocate.")
      CALL Delete(C)
      CALL New(C,NSMat_O=NSMat)
    ENDIF

#ifndef CHECK_UNRESTRICTED

    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
    CALL SymbolikMM(A,B,C,Flag,UpDate)
    IF(A%NSMat.EQ.1.AND.B%NSMat.EQ.1)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
    ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,1,2,2,Flag,Perf_O)
    ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,2,1,2,Flag,Perf_O)
    ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,2,2,2,Flag,Perf_O)
    ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,1,2,2,Flag,Perf_O)
      CALL NumerikMM(A,B,C,1,3,3,Flag,Perf_O)
      CALL NumerikMM(A,B,C,1,4,4,Flag,Perf_O)
    ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,2,1,2,Flag,Perf_O)
      CALL NumerikMM(A,B,C,3,1,3,Flag,Perf_O)
      CALL NumerikMM(A,B,C,4,1,4,Flag,Perf_O)
    ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
      CALL NumerikMM(A,B,C,1,1,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,2,3,1,Flag,Perf_O)
      CALL NumerikMM(A,B,C,1,2,2,Flag,Perf_O)
      CALL NumerikMM(A,B,C,2,4,2,Flag,Perf_O)
      CALL NumerikMM(A,B,C,3,1,3,Flag,Perf_O)
      CALL NumerikMM(A,B,C,4,3,3,Flag,Perf_O)
      CALL NumerikMM(A,B,C,3,2,4,Flag,Perf_O)
      CALL NumerikMM(A,B,C,4,4,4,Flag,Perf_O)
    ELSE
      CALL Halt('MultiplyM_BCSR: Error with NSMat!')
    ENDIF
    CALL Delete(Flag)

#else

    IF(NSMat.GT.1)THEN
      CALL New(Flag,NAtoms)
      CALL SetEq(Flag,0)
      CALL SymbolikMM(A,B,C,UpDate)
      N1=NBasF
      N2=NBasF*2
      write(*,'(A,3I3)') 'MultiplyM_BCSR: DGEMM are going to be done in dense format!',A%NSMat,B%NSMat,C%NSMat
      call seteq(AD,A)
      call seteq(BD,B)
      call seteq(CD,C)
      if(.NOT.UpDate)CD%D=0d0
      IF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
        CD%d(1:N1,   1:N1)=matmul(AD%d,BD%d(1:N1,   1:N1))+CD%d(1:N1,   1:N1)
        CD%d(1:N1,N1+1:N2)=matmul(AD%d,BD%d(1:N1,N1+1:N2))+CD%d(1:N1,N1+1:N2)
        !Check BCSR!
        CALL NumerikMM(A,B,C,1,1,1,Perf_O)
        CALL NumerikMM(A,B,C,1,2,2,Perf_O)
        call SetEq(DD,C)
        IF(SUM(ABS(CD%D-DD%D)).GT.1D-10) THEN
          WRITE(*,*) 'MMM:1-2 SUM(ABS(CD%D-DD%D)',SUM(ABS(CD%D-DD%D))
          STOP '56789'
        ENDIF
        CALL Delete(DD)
      ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
        CD%d(1:N1,   1:N1)=matmul(AD%d(1:N1,   1:N1),BD%d)+CD%d(1:N1,   1:N1)
        CD%d(1:N1,N1+1:N2)=matmul(AD%d(1:N1,N1+1:N2),BD%d)+CD%d(1:N1,N1+1:N2)
        !Check BCSR!
        CALL NumerikMM(A,B,C,1,1,1,Perf_O)
        CALL NumerikMM(A,B,C,2,1,2,Perf_O)
        call SetEq(DD,C)
        IF(SUM(ABS(CD%D-DD%D)).GT.1D-10) THEN
          WRITE(*,*) 'MMM:1-2 SUM(ABS(CD%D-DD%D)',SUM(ABS(CD%D-DD%D))
          STOP '56789'
        ENDIF
        CALL Delete(DD)
      ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
        CD%d(1:N1,   1:N1)=matmul(AD%d(1:N1,   1:N1),BD%d(1:N1,   1:N1))+CD%d(1:N1,   1:N1)
        CD%d(1:N1,N1+1:N2)=matmul(AD%d(1:N1,N1+1:N2),BD%d(1:N1,N1+1:N2))+CD%d(1:N1,N1+1:N2)
        !Check BCSR!
        CALL NumerikMM(A,B,C,1,1,1,Perf_O)
        CALL NumerikMM(A,B,C,2,2,2,Perf_O)
        call SetEq(DD,C)
        IF(SUM(ABS(CD%D-DD%D)).GT.1D-10) THEN
          WRITE(*,*) 'MMM:1-2 SUM(ABS(CD%D-DD%D)',SUM(ABS(CD%D-DD%D))
          STOP '56789'
        ENDIF
        CALL Delete(DD)
      ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
        CD%d(   1:N1,   1:N1)=matmul(AD%d,BD%d(   1:N1,   1:N1))+CD%d(   1:N1,   1:N1)
        CD%d(N1+1:N2,   1:N1)=matmul(AD%d,BD%d(N1+1:N2,   1:N1))+CD%d(N1+1:N2,   1:N1)
        CD%d(   1:N1,N1+1:N2)=matmul(AD%d,BD%d(   1:N1,N1+1:N2))+CD%d(   1:N1,N1+1:N2)
        CD%d(N1+1:N2,N1+1:N2)=matmul(AD%d,BD%d(N1+1:N2,N1+1:N2))+CD%d(N1+1:N2,N1+1:N2)
        !Check BCSR!
        CALL NumerikMM(A,B,C,1,1,1,Perf_O)
        CALL NumerikMM(A,B,C,1,2,2,Perf_O)
        CALL NumerikMM(A,B,C,1,3,3,Perf_O)
        CALL NumerikMM(A,B,C,1,4,4,Perf_O)
        call SetEq(DD,C)
        IF(SUM(ABS(CD%D-DD%D)).GT.1D-10) THEN
          WRITE(*,*) 'MMM:1-2 SUM(ABS(CD%D-DD%D)',SUM(ABS(CD%D-DD%D))
          STOP '56789'
        ENDIF
        CALL Delete(DD)
      ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
        CD%d(   1:N1,   1:N1)=matmul(AD%d(   1:N1,   1:N1),BD%d)+CD%d(   1:N1,   1:N1)
        CD%d(N1+1:N2,   1:N1)=matmul(AD%d(N1+1:N2,   1:N1),BD%d)+CD%d(N1+1:N2,   1:N1)
        CD%d(   1:N1,N1+1:N2)=matmul(AD%d(   1:N1,N1+1:N2),BD%d)+CD%d(   1:N1,N1+1:N2)
        CD%d(N1+1:N2,N1+1:N2)=matmul(AD%d(N1+1:N2,N1+1:N2),BD%d)+CD%d(N1+1:N2,N1+1:N2)
        !Check BCSR!
        CALL NumerikMM(A,B,C,1,1,1,Perf_O)
        CALL NumerikMM(A,B,C,2,1,2,Perf_O)
        CALL NumerikMM(A,B,C,3,1,3,Perf_O)
        CALL NumerikMM(A,B,C,4,1,4,Perf_O)
        call SetEq(DD,C)
        IF(SUM(ABS(CD%D-DD%D)).GT.1D-10) THEN
          WRITE(*,*) 'MMM:1-2 SUM(ABS(CD%D-DD%D)',SUM(ABS(CD%D-DD%D))
          STOP '56789'
        ENDIF
        CALL Delete(DD)
      ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
        CD%d=matmul(AD%d,BD%d)+CD%d
        !Check BCSR!
        CALL NumerikMM(A,B,C,1,1,1,Perf_O)
        CALL NumerikMM(A,B,C,2,3,1,Perf_O)
        CALL NumerikMM(A,B,C,1,2,2,Perf_O)
        CALL NumerikMM(A,B,C,2,4,2,Perf_O)
        CALL NumerikMM(A,B,C,3,1,3,Perf_O)
        CALL NumerikMM(A,B,C,4,3,3,Perf_O)
        CALL NumerikMM(A,B,C,3,2,4,Perf_O)
        CALL NumerikMM(A,B,C,4,4,4,Perf_O)
        call SetEq(DD,C)
        IF(SUM(ABS(CD%D-DD%D)).GT.1D-10) THEN
          WRITE(*,*) 'MMM:1-2 SUM(ABS(CD%D-DD%D)',SUM(ABS(CD%D-DD%D))
          STOP '56789'
        ENDIF
        CALL Delete(DD)
      ELSE
        CALL Halt('MultiplyM_BCSR: Error with NSMat!')
      ENDIF
      call seteq(C,CD,nsmat_o=nsmat)
      call delete(AD)
      call delete(BD)
      call delete(CD)
      CALL Delete(Flag)
    ELSE
      CALL New(Flag,NAtoms)
      CALL SetEq(Flag,0)
      CALL SymbolikMM(A,B,C,UpDate)
      CALL NumerikMM( A,B,C,1,1,1,Perf_O)
      CALL Delete(Flag)
    ENDIF

#endif
  END SUBROUTINE MultiplyM_BCSR

#ifdef PARALLEL
  !===============================================================================
  !     Matrix multiply for DBCSR matrices: C=Beta*C+A.B
  !===============================================================================
  SUBROUTINE MultiplyM_DBCSR(A,B,C,Beta_O,Perf_O,NoGlobal_O,Stop_O)
    USE COMMON_DEBUG
    REAL(DOUBLE),OPTIONAL,INTENT(IN)      :: Beta_O
    TYPE(DBCSR),          INTENT(INOUT)   :: A,B
    TYPE(DBCSR),          INTENT(INOUT)   :: C
    TYPE(TIME),  OPTIONAL,INTENT(INOUT)   :: Perf_O
    TYPE(DBCSR),    DIMENSION(0:NPrc-1)   :: U
    !        Pointer allocations
    TYPE(INT_VECT)                        :: SendReqst,RecvReqst, &
         SendSched,RecvSched, &
         SendPrior,RecvPrior, &
         ToDo
    TYPE(INT_RNK2)                        :: Stats
    !        Contiguous memory arrays
    REAL(DOUBLE),ALLOCATABLE,DIMENSION(:) :: UMTrix,BMTrix
    INTEGER,ALLOCATABLE,DIMENSION(:)      :: VBlks,VDisp
    !        small stuff
    INTEGER,        DIMENSION(0:NPrc-1)   :: V,VType
    INTEGER,EXTERNAL                      :: PostRecv77,PostSend77
    TYPE(TIME)                            :: Time0,Time1,Time2,Time3,  &
         Time4,Time5,Time6,Time7,Time8
    INTEGER,        DIMENSION(0:NPrc-1)   :: UMatPt,VMatPt,VNon0
    INTEGER                               :: I,J,K,L,N,Q,To,From,Tag,UTotal,VTotal, &
         ANBlks,BNBlks,Status,StartMem,StopMem,NSMat
    LOGICAL                               :: UpDate
    LOGICAL,OPTIONAL, INTENT(IN)          :: Stop_O,NoGlobal_O
    REAL(DOUBLE)                          :: FlOp
    CHARACTER(LEN=15)                     :: Prog='MultiplyM_DBCSR'
    TYPE(INT_VECT)                        :: Flag

    CALL Initialize(Flag)

    CALL AlignNodes()
    !         CALL Elapsed_TIME(Time0,'Init')
    !         CALL Elapsed_TIME(Time1,'Init')
    !-------------------------------------------------------------------------------
    !         IF(PRESENT(Stop_O))
    !         CALL PPrint(MemStats,'Begin MultiplyM_DBCSR')
    NSMat=MAX(A%NSMat,B%NSMat)
    IF(PRESENT(Beta_O))THEN
      IF(Beta_O/=Zero)THEN
        !              Check allocation
        IF(.NOT.AllocQ(C%Alloc)) &
             CALL Halt(' Asking for update of non-allocated' &
             //' DBCSR matrix in MultiplyM_DBCSR')
        !              If Beta /=1 then rescale
        IF(Beta_O/=One)THEN
          C%MTrix%D(1:C%NNon0)=Beta_O*C%MTrix%D(1:C%NNon0)
          PerfMon%FLOP=PerfMon%FLOP+DBLE(C%NNon0)
          IF(PRESENT(Perf_O)) &
               Perf_O%FLOP=Perf_O%FLOP+DBLE(C%NNon0)
        ENDIF
        UpDate=.TRUE.
      ELSE
        IF(.NOT.AllocQ(C%Alloc))CALL New(C,NSMat_O=NSMat)
        UpDate=.FALSE.
      ENDIF
    ELSE
      IF(.NOT.AllocQ(C%Alloc))CALL New(C,NSMat_O=NSMat)
      UpDate=.FALSE.
    ENDIF

    IF(C%NSMat.NE.NSMat) THEN
      write(*,*) 'MultiplyM_DBCSR: C%NSMat.NE.NSMat! Deallocate and reallocate.'
      CALL Delete(C)
      CALL New(C,NSMat_O=NSMat)
    ENDIF



    StartMem=MemStats%MemTab
    !-------------------------------------------------------------------------------
    !        Find global symbolic matrix structure if not alread known.

    CALL AlignNodes()
    CALL LocalToGlobal(A)
    CALL LocalToGlobal(B)
    C%GUpDate=STATUS_FALSE
    !-------------------------------------------------------------------------------
    !        Set up flag and global row pointers

    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
    CALL New(GlobalRowPtA,NAtoms+1)
    CALL New(GlobalRowPtB,NAtoms+1)
    CALL SetEq(GlobalRowPtA,0)
    CALL SetEq(GlobalRowPtB,0)
    !        Scheduling arrays
    CALL New(SendSched,NPrc-1)
    CALL New(RecvSched,NPrc-1)
    CALL New(SendPrior,NPrc-1)
    CALL New(RecvPrior,NPrc-1)
    CALL New(SendReqst,NPrc-1)
    CALL New(RecvReqst,NPrc-1)
    !-------------------------------------------------------------------------------
    K=0
    V=0
    U%NNon0=0
    DO N=0,NPrc-1
      IF(N/=MyId)THEN
        K=K+1
        RecvPrior%I(K)=RecvStruct(N,A,B,U(N))
        SendPrior%I(K)=SendStruct(N,A,B,V(N))
        RecvSched%I(K)=N
        SendSched%I(K)=N
      ENDIF
    ENDDO
    CALL Sort(RecvPrior,RecvSched)
    CALL Sort(SendPrior,SendSched)
    !-------------------------------------------------------------------------------
    !        Allocate contigous recieve buffers

    U(:)%NNon0=U(:)%NNon0
    UTotal=MAX(1,SUM(U(:)%NNon0))
    IF(UTotal/=0)THEN
      ALLOCATE(UMTrix(1:UTotal),STAT=Status)
      CALL IncMem(Status,0,UTotal)
    ENDIF
    !        Compute recieve segmentation pointers and allocate symbolic matrices
    Q=1
    DO N=1,NPrc-1
      IF(RecvPrior%I(N)/=FAIL)THEN
        From=RecvSched%I(N)
        UMatPt(From)=Q
        Q=Q+U(From)%NNon0
        CALL New(U(From),(/U(From)%NAtms,U(From)%NBlks,1/), &
             Node_O=From,NoGlobals_O=.TRUE.)
      ENDIF
    ENDDO
    !-------------------------------------------------------------------------------
    !        Allocate contigous send buffers

    ALLOCATE(BMTrix(1:B%NNon0),STAT=Status)
    CALL IncMem(Status,0,B%NNon0)
    BMTrix(1:B%NNon0)=B%MTrix%D(1:B%NNon0)

    VType=MPI_DATATYPE_NULL
    VTotal=MAX(1,SUM(V(:)))
    IF(VTotal/=0)THEN
      ALLOCATE(VBlks(1:VTotal),STAT=Status)
      CALL IncMem(Status,VTotal,0)
      ALLOCATE(VDisp(1:VTotal),STAT=Status)
      CALL IncMem(Status,VTotal,0)
    ENDIF
    !        -----------------------------------------------
    !         CALL Elapsed_TIME(Time1,'Accum') ! Preliminarys
    !         CALL Elapsed_TIME(Time2,'Init')  ! Recvs
    !-------------------------------------------------------------------------------
    !        Post non-blocking recieves

    DO N=1,NPrc-1
      IF(RecvPrior%I(N)/=FAIL)THEN
        From=RecvSched%I(N)
        Tag=From*MaxProc+MyId
        Q=UMatPt(From)
        BNBlks=B%GRwPt%I(NAtoms+1)-1
        !              Stupididy check #2
        IF(UTotal<Q+U(From)%NNon0-1) &
             CALL Halt(' Stupid error 1 in MultiplyMM_DBCSR')

        RecvReqst%I(N)=PostRecv_77(NAtoms,MyId,From,Tag,MONDO_COMM,            &
                                  A%NAtms,A%NBlks,B%NSMat,BNBlks,              &
                                  U(From)%NAtms,U(From)%NBlks,U(From)%NNon0,   &
                                  A%RowPt%I(1:A%NAtms+1),A%ColPt%I(1:A%NBlks), &
                                  B%GRwPt%I(1:NAtoms+1),B%GClPt%I(1:BNBlks),   &
                                  U(From)%RowPt%I(1:U(From)%NAtms+1),          &
                                  U(From)%ColPt%I(1:U(From)%NBlks),            &
                                  U(From)%BlkPt%I(1:U(From)%NBlks),            &
                                  UMTrix(Q:Q+U(From)%NNon0-1),                 &
                                  Flag%I(1:NAtoms),BSiz%I(1:NAtoms),           &
                                  Beg%I(From),End%I(From),                     &
                                  OffSt%I(MyId),OffSt%I(From))

        !RecvReqst%I(N)=PostRecv77(NAtoms,MyId,From,Tag,MONDO_COMM,             &
        !     A%NAtms,A%NBlks,B%NSMat,BNBlks,              &
        !     U(From)%NAtms,U(From)%NBlks,U(From)%NNon0,   &
        !     A%RowPt%I(1),A%ColPt%I(1),                   &
        !     B%GRwPt%I(1),B%GClPt%I(1),                   &
        !     U(From)%RowPt%I(1),U(From)%ColPt%I(1),       &
        !     U(From)%BlkPt%I(1),UMTrix(Q),                &
        !     Flag%I(1),BSiz%I(1),                         &
        !     Beg%I(From),End%I(From),                     &
        !     OffSt%I(MyId),OffSt%I(From))

      ELSE
        RecvReqst%I(N)=MPI_REQUEST_NULL
      ENDIF
    ENDDO
    !-------------------------------------------------------------------------------
    !        Make sure all recieves have been posted before comencing the sends
    CALL AlignNodes()
    !        ----------------------------------------
    !         CALL Elapsed_TIME(Time2,'Accum') ! Recvs
    !         CALL Elapsed_TIME(Time3,'Init')  ! Sends
    !-------------------------------------------------------------------------------
    !        Post sends

    Q=1
    DO N=1,NPrc-1
      IF(SendPrior%I(N)/=FAIL)THEN
        To=SendSched%I(N)
        Tag=MyId*MaxProc+To
        ANBlks=A%GRwPt%I(NAtoms+1)-1

        SendReqst%I(N)=PostSend_77(                                       &
                       NAtoms,ANBlks,B%NSMat,B%NAtms,B%NBlks,             &
                       B%NNon0,V(To),To,Tag,MyId,VType(To),MONDO_COMM,    &
                       A%GrwPt%I(1:NAtoms+1),A%GClPt%I(1:ANBlks),         &
                       B%RowPt%I(1:B%NAtms+1),B%ColPt%I(1:B%NBlks),       &
                       B%BlkPt%I(1:B%NBlks),BMTrix(1:B%NNon0),            &
                       VBlks(Q:Q+V(To)-1),VDisp(Q:Q+V(To)-1),             &
                       Flag%I(1:NAtoms),BSiz%I(1:NAtoms),OffSt%I(MyId),   &
                       Beg%I(MyId),End%I(MyId),Beg%I(To),End%I(To)        )

        !SendReqst%I(N)=PostSend77(                                            &
        !     NAtoms,ANBlks,B%NSMat,B%NAtms,B%NBlks,                 &
        !     B%NNon0,V(To),To,Tag,MyId,VType(To),MONDO_COMM,        &
        !     A%GrwPt%I(1),A%GClPt%I(1),B%RowPt%I(1),B%ColPt%I(1),   &
        !     B%BlkPt%I(1),BMTrix(1),VBlks(Q),VDisp(Q),              &
        !     Flag%I(1),BSiz%I(1),OffSt%I(MyId),                     &
        !     Beg%I(MyId),End%I(MyId),Beg%I(To),End%I(To)        )

        Q=Q+V(To)
      ELSE
        SendReqst%I(N)=MPI_REQUEST_NULL
      ENDIF
    ENDDO
#ifdef NONBLOCKING
#else
    !-------------------------------------------------------------------------------
    !        This call can really speed things up with blocking sends
    CALL AlignNodes()
#endif
    !        --------------------------------------------
    !         CALL Elapsed_TIME(Time3,'Accum') ! Sends
    !         CALL Elapsed_TIME(Time4,'Init')  ! Symbolic
    !-------------------------------------------------------------------------------
    !        Symbolic matrix-matrix multiply

    CALL SymbolikMM(A,B,C,UpDate) ! Local
    DO N=1,NPrc-1
      IF(RecvPrior%I(N)/=FAIL)THEN
        From=RecvSched%I(N)
        CALL SymbolikMM(A,U(From),C,.TRUE.) ! Non-local
      ENDIF
    ENDDO
    !        --------------------------------------------
    !         CALL Elapsed_TIME(Time4,'Accum') ! Symbolic
    !         CALL Elapsed_TIME(Time5,'Init')  ! Numeric
    !-------------------------------------------------------------------------------
    !        Local numerical matrix-matrix multiply
!!$         FlOp=Zero
!!$         CALL Load(A,GlobalRowPtA)
!!$         CALL Load(B,GlobalRowPtB)
!!$         CALL NumerikMM_GENERIC77(A%NAtms,OffSt%I(MyId),                        &
!!$                                  GlobalRowPtA%I(1),A%ColPt%I(1),A%BlkPt%I(1),A%MTrix%D(1), &
!!$                                  GlobalRowPtB%I(1),B%ColPt%I(1),B%BlkPt%I(1),BMTrix(1),    &
!!$                                  C%RowPt%I(1),C%ColPt%I(1),C%BlkPt%I(1),C%MTrix%D(1),      &
!!$                                  BSiz%I(1),Flag%I(1),FlOp)
!!$         PerfMon%FLOP=PerfMon%FLOP+FlOp
!!$         IF(PRESENT(Perf_O)) &
!!$           Perf_O%FLOP=Perf_O%FLOP+FlOp
!!$         CALL Clear(A,GlobalRowPtA)
!!$         CALL Clear(B,GlobalRowPtB)
!!$         !CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O) ! Local MM

    IF(    A%NSMat.EQ.1.AND.B%NSMat.EQ.1)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
    ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,1,2,2,Perf_O)
    ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,2,1,2,Perf_O)
    ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,2,2,2,Perf_O)
    ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,1,2,2,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,1,3,3,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,1,4,4,Perf_O)
    ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,2,1,2,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,3,1,3,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,4,1,4,Perf_O)
    ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
      CALL NumerikMM(A,B,BMTrix,C,1,1,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,2,3,1,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,1,2,2,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,2,4,2,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,3,1,3,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,4,3,3,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,3,2,4,Perf_O)
      CALL NumerikMM(A,B,BMTrix,C,4,4,4,Perf_O)
    ELSE
      CALL Halt('MultiplyM_DBCSR: Error with NSMat!')
    ENDIF

    !        --------------------------------------------
    !         CALL Elapsed_TIME(Time5,'Accum') ! Numeric
    !         CALL Elapsed_TIME(Time6,'Init')  ! WaitSome
    !        --------------------------------------------
    DO
      !           --------------------------------------------
      !            CALL Elapsed_TIME(Time6,'Start') ! WaitSome
      !           --------------------------------------------
      !           Wait for some buffers to fill
      CALL WaitSome(RecvReqst,ToDo)
      !           --------------------------------------------
      !            CALL Elapsed_TIME(Time6,'Accum') ! WaitSome
      !           --------------------------------------------
      !           If all done, exit
      IF(ToDo%I(1)==FAIL)EXIT
      !           Go over filled buffers
      DO N=1,SIZE(ToDo%I)
        From=RecvSched%I(ToDo%I(N))
        Q=UMatPt(From)
        !              ---------------------------------------------
        !               CALL Elapsed_TIME(Time5,'Start') ! Numeric
        !              ---------------------------------------------
        !              Non-local MM
!!$               FlOp=Zero
!!$               CALL Load(A,GlobalRowPtA)
!!$               CALL Load(U(From),GlobalRowPtB)
!!$               CALL NumerikMM_GENERIC77(A%NAtms,OffSt%I(MyId),                        &
!!$                                        GlobalRowPtA%I(1),A%ColPt%I(1),A%BlkPt%I(1),A%MTrix%D(1), &
!!$                                        GlobalRowPtB%I(1),U(From)%ColPt%I(1),                     &
!!$                                        U(From)%BlkPt%I(1),UMTrix(Q),                             &
!!$                                        C%RowPt%I(1),C%ColPt%I(1),C%BlkPt%I(1),C%MTrix%D(1),      &
!!$                                        BSiz%I(1),Flag%I(1),FlOp)
!!$               PerfMon%FLOP=PerfMon%FLOP+FlOp
!!$               IF(PRESENT(Perf_O)) &
!!$                  Perf_O%FLOP=Perf_O%FLOP+FlOp
!!$               CALL Clear(A,GlobalRowPtA)
!!$               CALL Clear(U(From),GlobalRowPtB)
!!$               !CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)

        IF(    A%NSMat.EQ.1.AND.B%NSMat.EQ.1)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
        ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,2,2,Perf_O)
        ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,2,1,2,Perf_O)
        ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,2,2,2,Perf_O)
        ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,2,2,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,3,3,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,4,4,Perf_O)
        ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,2,1,2,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,3,1,3,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,4,1,4,Perf_O)
        ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,1,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,2,3,1,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,1,2,2,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,2,4,2,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,3,1,3,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,4,3,3,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,3,2,4,Perf_O)
          CALL NumerikMM(A,U(From),UMTrix(Q:),C,4,4,4,Perf_O)
        ELSE
          CALL Halt('MultiplyM_DBCSR: Error with NSMat!')
        ENDIF

        !              ---------------------------------------------
        !               CALL Elapsed_TIME(Time5,'Accum') ! Numeric
        !              ---------------------------------------------
      ENDDO
    ENDDO
    !        ---------------------------------------
    !         CALL Elapsed_TIME(Time8,'Init')
    CALL AlignNodes()
    !         CALL Elapsed_TIME(Time8,'Accum')! Anomolous imbalance
    !         CALL Elapsed_TIME(Time7,'Init') ! Tidy
    !---------------------------------------------------------------------------
    !        Tidy up ...

    DO N=0,NPrc-1
      IF(AllocQ(U(N)%Alloc))  &
           CALL Delete(U(N))
    ENDDO
    IF(UTotal/=0)THEN
      DEALLOCATE(UMTrix,STAT=Status)
      CALL DecMem(Status,0,UTotal)
    ENDIF
    CALL Delete(RecvReqst)
    CALL Delete(SendSched)
    CALL Delete(RecvSched)
    CALL Delete(SendPrior)
    CALL Delete(RecvPrior)
    IF(AllocQ(ToDo%Alloc)) &
         CALL Delete(ToDo)
    CALL Delete(GlobalRowPtA)
    CALL Delete(GlobalRowPtB)
    CALL Delete(Flag)
#ifdef NONBLOCKING
    !-------------------------------------------------------------------------------

    !        Make sure sends complete before deallocating send buffers
    CALL WaitAll(SendReqst)
#endif
    CALL Delete(SendReqst)
    DEALLOCATE(BMTrix,STAT=Status)
    CALL DecMem(Status,0,B%NNon0)
    IF(VTotal/=0)THEN
      DEALLOCATE(VBlks,STAT=Status)
      CALL DecMem(Status,VTotal,0)
      DEALLOCATE(VDisp,STAT=Status)
      CALL DecMem(Status,VTotal,0)
    ENDIF
    DO N=0,NPrc-1
      IF(VType(N)/=MPI_DATATYPE_NULL)  &
           CALL FreeType(VType(N))
    ENDDO
    !-------------------------------------------------------------------------------
    !        Check for leaks ....

    StopMem=MemStats%MemTab
    IF(StartMem/=StopMem)&
         CALL Halt('Memory Leak in MultiplyM_DBCSR')
    !-------------------------------------------------------------------------------
    CALL AlignNodes()
    !        -------------------------------------------
    !         CALL Elapsed_TIME(Time7,'Accum') ! Tidy
    !-------------------------------------------------------------------------------
    !        Print preformance statistics
    !
    !         CALL Elapsed_TIME(Time0,'Accum')
    !         CALL Elapsed_TIME(Time1,Mssg_O='PrelimStff',Proc_O=Prog)
    !         CALL Elapsed_TIME(Time2,Mssg_O='Post Recvs',Proc_O=Prog)
    !         CALL Elapsed_TIME(Time2,Mssg_O='Post Sends',Proc_O=Prog)
    !         CALL Elapsed_TIME(Time4,Mssg_O='SymbolikMM',Proc_O=Prog)
    !         CALL Elapsed_TIME(Time5,Mssg_O='NumerikMM', Proc_O=Prog)
    !         CALL Elapsed_TIME(Time6,Mssg_O='WaitSome',  Proc_O=Prog)
    !         CALL Elapsed_TIME(Time8,Mssg_O='Imbalance', Proc_O=Prog)
    !         CALL Elapsed_TIME(Time7,Mssg_O='TidyUp',    Proc_O=Prog)
    !         CALL Elapsed_TIME(Time0,Mssg_O='Total',     Proc_O=Prog)

  END SUBROUTINE  MultiplyM_DBCSR

  !===============================================================================
  !     Determine structure of the recieve buffer DBCSR matrix U
  !===============================================================================
  FUNCTION RecvStruct(From,A,B,U)
    INTEGER                       :: RecvStruct
    INTEGER,       INTENT(INOUT)  :: From
    TYPE(DBCSR),   INTENT(INOUT)  :: A,B,U
    INTEGER                       :: S,T,JP,KP,IL,IG,JH,JL,JG,KG,MA,MB,MN_A,MN_B
    REAL(DOUBLE)                  :: FlOp
    !-------------------------------------------------------------------------------
    !        Check for a quick return

    IF(MyId==From.OR.From==MPI_PROC_NULL)THEN
      RecvStruct=FAIL
      RETURN
    ENDIF

    Flag%I=0
    !-------------------------------------------------------------------------------
    !        Check for overlaping row-col blocks and compute DBCSR dimensions

    T=0
    S=1
    FlOp=Zero
    DO JH=Beg%I(From),End%I(From)
      DO IL=1,A%NAtms
        IG=IL+OffSt%I(MyId)
        MA=BSiz%I(IG)
        DO JP=A%RowPt%I(IL),A%RowPt%I(IL+1)-1
          JG=A%ColPt%I(JP)
          IF(JG==JH)THEN
            JL=JG-OffSt%I(From)
            MB=BSiz%I(JG)
            MN_A=MA*MB                       *A%NSMat   !<<<<< add spin here?
            IF(Flag%I(JL)==0)THEN
              Flag%I(JL)=1
              DO KP=B%GRwPt%I(JG),B%GRwPt%I(JG+1)-1
                KG=B%GClPt%I(KP)
                T=T+1
                MN_B=MB*BSiz%I(KG)         *B%NSMat   !<<<<< add spin here?
                S=S+MN_B
                FlOp=FlOp+MA*MN_B
              ENDDO
            ELSE
              DO KP=B%GRwPt%I(JG),B%GRwPt%I(JG+1)-1
                KG=B%GClPt%I(KP)
                FlOp=FlOp+MN_A*BSiz%I(KG)
              ENDDO
            ENDIF
            GOTO 11 ! Sometimes, ya got to have one ...
          ENDIF
        ENDDO
      ENDDO
11    CONTINUE
    ENDDO
    CALL SetEq(Flag,0)
    !        Return if no overlap ...
    IF(T==0)THEN
      RecvStruct=FAIL
      RETURN
    ENDIF
    !        Assign the matrix dimensions ...
    U%NAtms=End%I(From)-Beg%I(From)+1
    U%NNon0=S
    U%NBlks=T

    RecvStruct=FLOOR(FlOp)/DBLE(S)
    !         RecvStruct=S
    !         RecvStruct=RANDOM_INT((/1,1000/))

  END FUNCTION RecvStruct
  !===============================================================================
  !     Determine structure of the send buffer DBCSR matrix V
  !===============================================================================
  FUNCTION SendStruct(To,A,B,V)
    INTEGER                       :: SendStruct
    INTEGER,       INTENT(INOUT)  :: To,V
    TYPE(DBCSR),   INTENT(INOUT)  :: A,B
    INTEGER                       :: S,T,JH,JP,KP,IG,JG,JL,KG, &
         MA,MB,MN_A,MN_B
    REAL(DOUBLE)                  :: FlOp
    !-------------------------------------------------------------------------------
    !        Check for a quick return

    IF(MyId==To.OR.To==MPI_PROC_NULL)THEN
      SendStruct=FAIL
      RETURN
    ENDIF
    !-------------------------------------------------------------------------------
    !        Overkill....

    Flag%I=0
    !-------------------------------------------------------------------------------
    !        Find overlaping row-col blocks in B%MTrix and commark for SEND

    T=0
    S=1
    FlOp=Zero
    DO JH=Beg%I(MyId),End%I(MyId)
      DO IG=Beg%I(To),End%I(To)
        MA=BSiz%I(IG)
        DO JP=A%GRwPt%I(IG),A%GRwPt%I(IG+1)-1
          JG=A%GClPt%I(JP)
          IF(JG==JH)THEN
            JL=JG-OffSt%I(MyId)
            IF(Flag%I(JL)==0)THEN
              MB=BSiz%I(JG)
              MN_A=MA*MB                          *A%NSMat  !<<<<< add spin here?
              Flag%I(JL)=1
              DO KP=B%RowPt%I(JL),B%RowPt%I(JL+1)-1
                T=T+1
                KG=B%ColPt%I(KP)
                MN_B=MB*BSiz%I(KG)               *B%NSMat  !<<<<< add spin here?
                S=S+MN_B
                FlOp=FlOp+MA*MN_B
              ENDDO
            ELSE
              DO KP=B%RowPt%I(JL),B%RowPt%I(JL+1)-1
                KG=B%ColPt%I(KP)
                FlOp=FlOp+MN_A*BSiz%I(KG)
              ENDDO
            ENDIF
            GOTO 11 ! Sometimes, ya got to have one ...
          ENDIF
        ENDDO
      ENDDO
11    CONTINUE
    ENDDO
    CALL SetEq(Flag,0)
    !        Return if no computed overlap
    IF(T==0)THEN
      SendStruct=FAIL
      RETURN
    ENDIF
    !        Assign dimensions of the buffer DBCSR matrix
    V=T
    SendStruct=FLOOR(FlOp/DBLE(S))
    !         SendStruct=S
    !         SendStruct=RANDOM_INT((/1,1000/))
  END FUNCTION SendStruct
  !===============================================================================
  !     Post non-blocking ISend and IRecv for DBCSR data exchange
  !===============================================================================
  INTEGER FUNCTION PostRecv_77(NAtoms,MyId,From,Tag,COMM,ANAtms,           &
       ANBlks,BNSMat,BNBlks,UNAtms,UNBlks,UNNon0,  &
       ARowPt,AColPt,BGRwPt,BGClPt,                &
       URowPt,UColPt,UBlkPt,UMTrix,                &
       Flag77,BSiz77,FromBeg,FromEnd,              &
       MyOff,FromOff)
    IMPLICIT NONE
    INTEGER :: I,S,T,IG,IL,KG,KP,JL,JP,JG,JH,                 &
         NAtoms,MyId,From,Tag,COMM,                         &
         ANAtms,ANBlks,BNSMat,BNBlks,UNAtms,UNBlks,UNNon0,  &
         FromBeg,FromEnd,MyOff,FromOff,IErr
    REAL(DOUBLE), DIMENSION(1:UNNon0) :: UMTrix
    INTEGER,      DIMENSION(1:ANAtms+1)    :: ARowPt
    INTEGER,      DIMENSION(1:ANBlks)      :: AColPt
    INTEGER,      DIMENSION(1:NAtoms+1)    :: BGRwPt
    INTEGER,      DIMENSION(1:BNBlks)      :: BGClPt
    INTEGER,      DIMENSION(1:UNAtms+1)    :: URowPt
    INTEGER,      DIMENSION(1:UNBlks)      :: UColPt
    INTEGER,      DIMENSION(1:UNBlks)      :: UBlkPt
    INTEGER,      DIMENSION(1:NAtoms)      :: Flag77
    INTEGER,      DIMENSION(1:NAtoms)      :: BSiz77

    !         INTEGER BIG_INT
    !         PARAMETER(BIG_INT=100000)
    !         INCLUDE 'MONDO_MPI_INCLUDE.Inc'
    !C
    DO I=1,NAtoms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Flag77(I)=0
    ENDDO
    !-------------------------------------------------------------------------------
    !        Find overlaping row-col blocks and compute pointers

    T=1
    S=1
    CALL INT_VECT_EQ_INT_SCLR(UNAtms,URowPt,BIG_INT)
    DO JH=FromBeg,FromEnd
      DO IL=1,ANAtms
        IG=IL+MyOff
        DO JP=ARowPt(IL),ARowPt(IL+1)-1
          JG=AColPt(JP)
          IF(JG==JH)THEN
            JL=JG-FromOff
            IF(Flag77(JL)==0)THEN
              Flag77(JL)=1
              URowPt(JL)=T
              DO KP=BGRwPt(JG),BGRwPt(JG+1)-1
                KG=BGClPt(KP)
                UColPt(T)=KG
                UBlkPt(T)=S
                T=T+1
                S=S+BSiz77(JG)*BSiz77(KG)          *BNSMat       !<<<<< add spin here?
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    CALL INT_VECT_EQ_INT_SCLR(NAtoms,Flag77,0)
    !-------------------------------------------------------------------------------
    URowPt(UNAtms+1)=T
    DO I=UNAtms+1,2,-1
      IF(URowPt(I-1)==BIG_INT)THEN
        URowPt(I-1)=URowPt(I)
      ENDIF
    ENDDO
    !-------------------------------------------------------------------------------
    CALL MPI_IRECV(UMTrix,UNNon0,MPI_DOUBLE_PRECISION,From,Tag, &
         COMM,PostRecv_77,IErr)
    CALL ErrChk(IErr,'PostRecv_77')
    RETURN
  END FUNCTION PostRecv_77


  INTEGER FUNCTION PostSend_77(NAtoms,ANBlks,BNSMat,BNAtms,BNBlks,  &
       BNNon0,V,To,Tag,MyId,VType,COMM,       &
       AGrwPt,AGClPt,BRowPt,BColPt,           &
       BBlkPt,BMTrix,VBlks,VDisp,Flag77,      &
       BSiz77,MyOff,MyBeg,MyEnd,ToBeg,ToEnd)
    IMPLICIT NONE
    INTEGER :: NAtoms,ANBlks,BNSMat,BNAtms,BNBlks,BNNon0,V,To,Tag,MyId, &
         VType,COMM,MyOff,MyBeg,MyEnd,ToBeg,ToEnd
    INTEGER :: IErr,NAtms,I,K,MN,S,T,TT,JP,KP,IG,JG,JL,JH,KG

    REAL(DOUBLE), DIMENSION(BNNon0)   :: BMTrix
    INTEGER,      DIMENSION(NAtoms)   :: BSiz77
    INTEGER,      DIMENSION(NAtoms)   :: Flag77
    INTEGER,      DIMENSION(V)        :: VBlks
    INTEGER,      DIMENSION(V)        :: VDisp
    INTEGER,      DIMENSION(NAtoms+1) :: AGrwPt
    INTEGER,      DIMENSION(ANBlks)   :: AGClPt
    INTEGER,      DIMENSION(BNAtms+1) :: BRowPt
    INTEGER,      DIMENSION(BNBlks)   :: BColPt
    INTEGER,      DIMENSION(BNBlks)   :: BBlkPt
    !         INCLUDE 'MONDO_MPI_INCLUDE.Inc'

    DO I=1,NAtoms !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Flag77(I)=0
    ENDDO
    !-------------------------------------------------------------------------------
    !        Find overlaping row-col blocks in BMTrix and mark for SEND

    T=0
    S=1
    NAtms=MyEnd-MyBeg+1
    DO JH=MyBeg,MyEnd
      DO IG=ToBeg,ToEnd
        DO JP=AGRwPt(IG),AGRwPt(IG+1)-1
          JG=AGClPt(JP)
          IF(JG.EQ.JH)THEN
            JL=JG-MyOff
            IF(Flag77(JL).EQ.0)THEN
              Flag77(JL)=1
              DO KP=BRowPt(JL),BRowPt(JL+1)-1
                KG=BColPt(KP)
                MN=BSiz77(JG)*BSiz77(KG)                 *BNSMat       !<<<<< add spin here?
                !C                          Check for contigous blocks ..
                !C                          IF(T.NE.0.AND.VDisp(T)+MN.EQ.BBlkPt(KP)-1)THEN
                !C                              VBlks(T)=VBlks(T)+MN
                !C                           ELSE
                T=T+1
                VDisp(T)=BBlkPt(KP)-1
                VBlks(T)=MN
                !C                           ENDIF
                S=S+MN
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CALL INT_VECT_EQ_INT_SCLR(NAtoms,Flag77,0)
    CALL MPI_TYPE_INDEXED(T,VBlks,VDisp, &
         MPI_DOUBLE_PRECISION,VType,IErr)
    CALL ErrChk(IErr,'PostSend #1')
    CALL MPI_TYPE_COMMIT(VType,IErr)
    CALL ErrChk(IErr,'PostSend #2')
#ifdef NONBLOCKING
    !        Post a non-blocking send
    CALL MPI_ISEND(BMTrix,1,VType,To,Tag,COMM,PostSend_77,IErr)
    CALL ErrChk(IErr,'Nonblocking PostSend #3')
#else
    !        Post a blocking send
    CALL MPI_SEND(BMTrix,1,VType,To,Tag,COMM,IErr)
    CALL ErrChk(IErr,'Blocking PostSend #3')
    PostSend_77=MPI_REQUEST_NULL
#endif
    RETURN
  END FUNCTION PostSend_77
  !===============================================================================
  !     Matrix-scalar multiply for DCSR matrices
  !===============================================================================
  SUBROUTINE MultiplyM_DBCSR_SCLR(A,B,Expert_O,Perf_O)
    REAL(DOUBLE),         INTENT(IN)    :: B
    TYPE(DBCSR),          INTENT(INOUT) :: A
    INTEGER   ,  OPTIONAL,INTENT(IN)    :: Expert_O
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                             :: Expert
    INTEGER                             :: AtA,AtB,MA,MB,P,R
    !-------------------------------------------------------------------------------
    Expert=0
    IF(PRESENT(Expert_O))Expert=Expert_O
    IF(Expert.EQ.0)THEN
      CALL DBL_Scale(A%NNon0,A%MTrix%D,B)
      !CALL DSCAL(A%NNon0,B,A%MTrix%D,1)
      PerfMon%FLOP=PerfMon%FLOP+DBLE(A%NNon0)
    ELSEIF(Expert.GT.0.AND.Expert.LT.5)THEN
      DO AtA=1,A%NAtms
        MA=BSiz%I(AtA)
        DO P=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1
          AtB =A%ColPt%I(P)
          R   =A%BlkPt%I(P)
          MB  =BSiz%I(AtB)
          CALL DBL_Scale(MA*MB,A%MTrix%D(R+(Expert-1)*MA*MB),B)
          !CALL DSCAL(MA*MB,B,A%MTrix%D(R+(Expert-1)*MA*MB),1)
        ENDDO
      ENDDO
      PerfMon%FLOP=PerfMon%FLOP+DBLE(A%NNon0)/DBLE(A%NSMat)
    ELSE
      CALL Halt('MultiplyM_DBCSR_SCLR: wrong value for Expert!')
    ENDIF
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+DBLE(A%NNon0)
  END SUBROUTINE MultiplyM_DBCSR_SCLR
#endif
  !===============================================================================
  !     Matrix-scalar multiply for BCSR matrices
  !===============================================================================
  SUBROUTINE MultiplyM_BCSR_SCLR(A,B,Expert_O,Perf_O)
    REAL(DOUBLE),         INTENT(IN)    :: B
    TYPE(BCSR),           INTENT(INOUT) :: A
    INTEGER   ,  OPTIONAL,INTENT(IN)    :: Expert_O
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                             :: Expert
    INTEGER                             :: AtA,AtB,MA,MB,P,R
    !-------------------------------------------------------------------------------
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      Expert=0
      IF(PRESENT(Expert_O))Expert=Expert_O
      IF(Expert.EQ.0)THEN
        CALL DBL_Scale(A%NNon0,A%MTrix%D,B)
        !CALL DSCAL(A%NNon0,B,A%MTrix%D,1)
        PerfMon%FLOP=PerfMon%FLOP+DBLE(A%NNon0)
      ELSEIF(Expert.GT.0.AND.Expert.LT.5)THEN
        DO AtA=1,A%NAtms
          MA=BSiz%I(AtA)
          DO P=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1
            AtB =A%ColPt%I(P)
            R   =A%BlkPt%I(P)
            MB  =BSiz%I(AtB)
            CALL DBL_Scale(MA*MB,A%MTrix%D(R+(Expert-1)*MA*MB),B)
            !CALL DSCAL(MA*MB,B,A%MTrix%D(R+(Expert-1)*MA*MB),1)
          ENDDO
        ENDDO
        PerfMon%FLOP=PerfMon%FLOP+DBLE(A%NNon0)/DBLE(A%NSMat)
      ELSE
        CALL Halt('MultiplyM_BCSR_SCLR: wrong value for Expert!')
      ENDIF
#ifdef PARALLEL
    ENDIF
#endif
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+DBLE(A%NNon0)
  END SUBROUTINE MultiplyM_BCSR_SCLR
  !===============================================================================
  !     BCSR: Matrix Vector Multiply
  !===============================================================================
  SUBROUTINE  MultiplyM_BCSR_VECT(A,Vin,Vout)
    TYPE(BCSR)       :: A
    TYPE(DBL_VECT)   :: Vin,Vout
    INTEGER          :: AtA,AtB,Abeg,Aend,P,R,MA,MB,IJ,I,J,IA,IB
    !
    DO I=1,NBasF
      Vout%D(I)=Zero
    ENDDO
    !
    IA = 0
    DO AtA=1,A%NAtms
      MA   = BSiz%I(AtA)
      Abeg = A%RowPt%I(AtA)
      Aend = A%RowPt%I(AtA+1)-1
      DO P=Abeg,Aend
        AtB  = A%ColPt%I(P)
        R    = A%BlkPt%I(P)
        MB   = BSiz%I(AtB)
        IB   = 0
        DO I=1,AtB-1
          IB = IB + BSiz%I(I)
        ENDDO
        DO J=1,MB
          DO I=1,MA
            Vout%D(IA+I) = Vout%D(IA+I) + A%MTrix%D(R)*Vin%D(IB+J)
            R = R + 1
          ENDDO
        ENDDO
      ENDDO
      IA = IA + MA
    ENDDO
    !
  END SUBROUTINE MultiplyM_BCSR_VECT
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !     MATRIX ADDITION, MATRIX ADDITION, MATRIX ADDITION, MATRIX ADDITION, MATRIX ADDITION
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef PARALLEL
  !===============================================================================
  !     DBCSR wrapper for generic matrix addition
  !===============================================================================
  SUBROUTINE Add_DBCSR(A,B,C,Perf_O,NoGlobal_O)
    IMPLICIT NONE
    TYPE(DBCSR),        INTENT(INOUT) :: A,B
    TYPE(DBCSR),        INTENT(INOUT) :: C
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                           :: Status,NSMat
    REAL(DOUBLE)                      :: FlOp
    LOGICAL,OPTIONAL, INTENT(IN)      :: NoGlobal_O
    !-------------------------------------------------------------------------------
    NSMat=MAX(A%NSMat,B%NSMat)
    IF(.NOT.AllocQ(C%Alloc)   ) CALL New(C,NSMat_O=NSMat)
    !
    IF(NSMat.NE.C%NSMat) THEN
      write(*,*) 'Add_DBCSR: NSMat.NE.C%NSMat! Deallocate-reallocate'
      CALL Delete(C)
      CALL New(C,NSMat_O=NSMat)
    ENDIF
    C%GUpDate=STATUS_FALSE
    !
    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
    IF(    A%NSMat.EQ.1.AND.B%NSMat.EQ.1)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
    ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
      CALL NumerikADD(A,B,C,1,2,2,Perf_O)
    ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
      CALL NumerikADD(A,B,C,2,1,2,Perf_O)
    ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
      CALL NumerikADD(A,B,C,2,2,2,Perf_O)
    ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
      CALL NumerikADD(A,B,C,1,4,4,Perf_O)
    ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
      CALL NumerikADD(A,B,C,4,1,4,Perf_O)
    ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
      CALL NumerikADD(A,B,C,1,1,1,Perf_O)
      CALL NumerikADD(A,B,C,2,2,2,Perf_O)
      CALL NumerikADD(A,B,C,3,3,3,Perf_O)
      CALL NumerikADD(A,B,C,4,4,4,Perf_O)
    ELSE
      CALL Halt('Add_BCSR: Error with NSMat!')
    ENDIF
    CALL Delete(Flag)
    !
    IF(Status==FAIL) THEN
      WRITE(*,*) "Status = ",Status
      WRITE(*,*) A%NAtms,A%NBlks,A%NNon0,SIZE(A%MTrix%D)
      WRITE(*,*) B%NAtms,B%NBlks,B%NNon0,SIZE(B%MTrix%D)
      WRITE(*,*) C%NAtms,C%NBlks,C%NNon0,SIZE(C%MTrix%D)
      CALL Halt('Dimensions in Add_BCSR')
    ENDIF
  END SUBROUTINE Add_DBCSR

  SUBROUTINE NumerikADD_DBCSR(A,B,C,ASMat,BSMat,CSMat,Perf_O) !<<<SPIN
    IMPLICIT NONE
    TYPE(DBCSR),        INTENT(IN)    :: A,B
    TYPE(DBCSR),        INTENT(INOUT) :: C
    INTEGER            ,INTENT(IN)    :: ASMat,BSMat,CSMat !<<<SPIN
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                           :: Status
    REAL(DOUBLE)                      :: FlOp
    !-------------------------------------------------------------------------------
    CALL SetEq(Flag,0)
    FlOp=Zero
    Status=Add_GENERIC(ASMat,BSMat,CSMat,C%NSMat,         &
         SIZE(C%ColPt%I),SIZE(C%MTrix%D),         &
         A%NAtms,OffSt%I(MyId),                   &
         C%NAtms,C%NBlks,C%NNon0,                 &
         A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
         C%RowPt%I,C%ColPt%I,C%BlkPt%I,C%MTrix%D, &
         BSiz%I,Flag%I,Flop)
    IF(Status==FAIL) CALL Halt('Dimensions in NumerikADD_DBCSR')
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE NumerikADD_DBCSR
#endif

  !===============================================================================
  !     Wrapper for generic F77 style BCSR numeric matrix add
  !===============================================================================
  SUBROUTINE NumerikADD_BCSR(A,B,C,ASMat,BSMat,CSMat,Flag,Perf_O) !<<<SPIN
    IMPLICIT NONE
    TYPE(BCSR),         INTENT(IN)    :: A,B
    TYPE(BCSR),         INTENT(INOUT) :: C
    INTEGER, INTENT(IN)               :: ASMat,BSMat,CSMat !<<<SPIN
    TYPE(INT_VECT), INTENT(INOUT)     :: Flag
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                           :: Status
    REAL(DOUBLE)                      :: FlOp

    CALL SetEq(Flag,0)
    FlOp=Zero
    Status=Add_GENERIC(ASMat,BSMat,CSMat,C%NSMat,               &
         SIZE(C%ColPt%I),SIZE(C%MTrix%D),         &
         A%NAtms,0,                               &
         C%NAtms,C%NBlks,C%NNon0,                 &
         A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
         C%RowPt%I,C%ColPt%I,C%BlkPt%I,C%MTrix%D, &
         BSiz%I,Flag%I,Flop)
    IF(Status==FAIL) CALL Halt('Dimensions in NumerikADD_BCSR')
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE NumerikADD_BCSR

  !===============================================================================
  !     BCSR wrapper for generic matrix addition
  !===============================================================================
  SUBROUTINE Add_BCSR(A,B,C,Perf_O)
    IMPLICIT NONE
    TYPE(BCSR),         INTENT(INOUT) :: A,B
    TYPE(BCSR),         INTENT(INOUT) :: C
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                           :: Status
    INTEGER                           :: NSMat,N2
    TYPE(INT_VECT)                    :: Flag
#ifdef CHECK_UNRESTRICTED
    REAL(DOUBLE)                      :: FlOp
    TYPE(DBL_RNK2)                    :: AD,BD,CD,DD
    CALL Initialize(AD)
    CALL Initialize(BD)
    CALL Initialize(CD)
    CALL Initialize(DD)
#endif
    CALL Initialize(Flag)

    NSMat=MAX(A%NSMat,B%NSMat)
    IF(.NOT.AllocQ(C%Alloc)) THEN
      !CALL MondoLog(DEBUG_MAXIMUM, "Add_BCSR", "C not allocated, allocating")
      CALL New(C,NSMat_O=NSMat)
      CALL DSCAL(SIZE(C%MTrix%D),Zero,C%MTrix%D,1)
    ELSE
      !CALL MondoLog(DEBUG_MAXIMUM, "Add_BCSR", "C already allocated, reallocating")
      CALL Delete(C)
      CALL New(C,NSMat_O=NSMat)
      CALL DSCAL(SIZE(C%MTrix%D),Zero,C%MTrix%D,1)
    ENDIF

    IF(NSMat.NE.C%NSMat) THEN
      !CALL MondoLog(DEBUG_MAXIMUM, "Add_BCSR", "NSMat.NE.C%NSMat! Deallocate-reallocate")
      CALL Delete(C)
      CALL New(C,NSMat_O=NSMat)
      CALL DSCAL(SIZE(C%MTrix%D),Zero,C%MTrix%D,1)
    ENDIF

#ifndef CHECK_UNRESTRICTED
    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(A%NSMat.EQ.1.AND.B%NSMat.EQ.1)THEN
        CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
      ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
        CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
        CALL NumerikADD(A,B,C,1,2,2,Flag,Perf_O)
      ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
        CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
        CALL NumerikADD(A,B,C,2,1,2,Flag,Perf_O)
      ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
        CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
        CALL NumerikADD(A,B,C,2,2,2,Flag,Perf_O)
      ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
         ! Fill in off diagonal spin blocks ...
         CALL SetEq(C,B)
         ! ... then add in diagonal blocks
         CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
         CALL NumerikADD(A,B,C,1,4,4,Flag,Perf_O)
      ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
         ! Fill in off diagonal spin blocks ...
         CALL SetEq(C,A)
         ! ... then add in diagonal blocks
        CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
        CALL NumerikADD(A,B,C,4,1,4,Flag,Perf_O)
      ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
        CALL NumerikADD(A,B,C,1,1,1,Flag,Perf_O)
        CALL NumerikADD(A,B,C,2,2,2,Flag,Perf_O)
        CALL NumerikADD(A,B,C,3,3,3,Flag,Perf_O)
        CALL NumerikADD(A,B,C,4,4,4,Flag,Perf_O)
      ELSE
        CALL Halt('Add_BCSR: Error with NSMat!')
      ENDIF
#ifdef PARALLEL
    ENDIF
#endif
    CALL Delete(Flag)
#else
    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(A%NSMat.EQ.1.AND.B%NSMat.EQ.1) THEN
        FlOp=Zero
        Status=Add_GENERIC(1,1,1,1,                                 &
             SIZE(C%ColPt%I),SIZE(C%MTrix%D),         &
             A%NAtms,0,                               &
             C%NAtms,C%NBlks,C%NNon0,                 &
             A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
             B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
             C%RowPt%I,C%ColPt%I,C%BlkPt%I,C%MTrix%D, &
             BSiz%I,Flag%I,Flop)
        IF(Status==FAIL) THEN
          WRITE(*,*) "Status = ",Status
          WRITE(*,*) A%NAtms,A%NBlks,A%NNon0,SIZE(A%MTrix%D)
          WRITE(*,*) B%NAtms,B%NBlks,B%NNon0,SIZE(B%MTrix%D)
          WRITE(*,*) C%NAtms,C%NBlks,C%NNon0,SIZE(C%MTrix%D)
          CALL Halt('Dimensions in Add_BCSR')
        ENDIF
      ELSE
        write(*,'(A,3I3)') 'Add_BCSR: The add will be done in a dense way!',A%NSMat,B%NSMat,C%NSMat
        N1=NBasF
        N2=2*NBasF
        call seteq(AD,A)
        call seteq(BD,B)
        call seteq(CD,C)
        CD%D=0d0
        IF(    A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
          CD%D=AD%D+BD%D
          !Check ADD!
          CALL NumerikADD(A,B,C,1,1,1,Perf_O)
          CALL NumerikADD(A,B,C,2,2,2,Perf_O)
          CALL SetEq(dd,c)
          IF(SUM(ABS(cd%d-dd%d)).GT.1D-10) THEN
            write(*,*) 'Add 2-2 cd-dd=',SUM(ABS(cd%d-dd%d)),SUM(ABS(cd%d)),SUM(ABS(dd%d))
            STOP '123456'
          ENDIF
        ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1)THEN
          CD%D(1:N1,   1:N1)=AD%D(1:N1,   1:N1)+BD%D
          CD%D(1:N1,N1+1:N2)=AD%D(1:N1,N1+1:N2)+BD%D
          !Check ADD!
          CALL NumerikADD(A,B,C,1,1,1,Perf_O)
          CALL NumerikADD(A,B,C,2,1,2,Perf_O)
          CALL SetEq(dd,c)
          IF(SUM(ABS(cd%d-dd%d)).GT.1D-10) THEN
            write(*,*) 'Add 2-1 cd-dd=',SUM(ABS(cd%d-dd%d)),SUM(ABS(cd%d)),SUM(ABS(dd%d))
            STOP '123456'
          ENDIF
        ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2)THEN
          CD%D(1:N1,   1:N1)=AD%D+BD%D(1:N1,   1:N1)
          CD%D(1:N1,N1+1:N2)=AD%D+BD%D(1:N1,N1+1:N2)
          !Check ADD!
          CALL NumerikADD(A,B,C,1,1,1,Perf_O)
          CALL NumerikADD(A,B,C,1,2,2,Perf_O)
          CALL SetEq(dd,c)
          IF(SUM(ABS(cd%d-dd%d)).GT.1D-10) THEN
            write(*,*) 'Add 1-2 cd-dd=',SUM(ABS(cd%d-dd%d)),SUM(ABS(cd%d)),SUM(ABS(dd%d))
            STOP '123456'
          ENDIF
        ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1)THEN
          CD%D(1:N1   ,   1:N1)=AD%D(1:N1   ,   1:N1)+BD%D
          CD%D(N1+1:N2,N1+1:N2)=AD%D(N1+1:N2,N1+1:N2)+BD%D
          CD%D(1:N1   ,N1+1:N2)=AD%D(1:N1   ,N1+1:N2)
          CD%D(N1+1:N2,   1:N1)=AD%D(N1+1:N2,   1:N1)
          write(*,*) 'Warning: Add 4-1 isn''t implemented in BCSR form yet'
          call seteq(dd,a)
        ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4)THEN
          CD%D(1:N1   ,   1:N1)=AD%D+BD%D(1:N1   ,   1:N1)
          CD%D(N1+1:N2,N1+1:N2)=AD%D+BD%D(N1+1:N2,N1+1:N2)
          CD%D(1:N1   ,N1+1:N2)=     BD%D(1:N1   ,N1+1:N2)
          CD%D(N1+1:N2,   1:N1)=     BD%D(N1+1:N2,   1:N1)
          write(*,*) 'Warning: Add 1-4 isn''t implemented in BCSR form yet'
          call seteq(dd,a)
        ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
          CD%D=AD%D+BD%D
          !Check ADD!
          CALL NumerikADD(A,B,C,1,1,1,Perf_O)
          CALL NumerikADD(A,B,C,2,2,2,Perf_O)
          CALL NumerikADD(A,B,C,3,3,3,Perf_O)
          CALL NumerikADD(A,B,C,4,4,4,Perf_O)
          CALL SetEq(dd,c)
          IF(SUM(ABS(cd%d-dd%d)).GT.1D-10) THEN
            write(*,*) 'Add 4-4 cd-dd=',SUM(ABS(cd%d-dd%d)),SUM(ABS(cd%d)),SUM(ABS(dd%d))
            STOP '123456'
          ENDIF
        ELSE
          call Halt('Add_BCSR: Something wrong!')
        ENDIF
        call seteq(C,CD,nsmat_o=nsmat)

        CALL Delete(dd)
        call delete(AD)
        call delete(BD)
        call delete(CD)
      ENDIF
#ifdef PARALLEL
    ENDIF
#endif
    CALL Delete(Flag)
#endif
  END SUBROUTINE Add_BCSR

  !===============================================================================
  !     Generic F77 style (D)BCSR matrix addition
  !===============================================================================
  FUNCTION Add_GENERIC(ASMat,BSMat,CSMat,NSMat,       & !<<<SPIN
       MxBlks,MxNon0,ANAtms,AOffSt,   &
       CNAtms,CNBlks,CNNon0,          &
       ARowPt,AColPt,ABlkPt,AMTrix,   &
       BRowPt,BColPt,BBlkPt,BMTrix,   &
       CRowPt,CColPt,CBlkPt,CMTrix,   &
       BSiz,Flag,Flop)
    IMPLICIT NONE
    INTEGER                                 :: Add_GENERIC
    INTEGER,                  INTENT(IN)    :: ASMat,BSMat,CSMat,NSMat,MxBlks,MxNon0,ANAtms,AOffSt !<<<SPIN
    REAL(DOUBLE),DIMENSION(:),INTENT(IN)    :: AMTrix,BMTrix
    INTEGER,     DIMENSION(:),INTENT(IN)    :: ARowPt,AColPt,ABlkPt, &
         BRowPt,BColPt,BBlkPt, &
         BSiz
    INTEGER,                  INTENT(INOUT)   :: CNAtms,CNBlks,CNNon0
    REAL(DOUBLE),DIMENSION(:),INTENT(INOUT)   :: CMTrix
    INTEGER,     DIMENSION(:),INTENT(INOUT)   :: CRowPt,CColPt,CBlkPt
    INTEGER,     DIMENSION(:),INTENT(INOUT) :: Flag
    REAL(DOUBLE),             INTENT(INOUT) :: FlOp
    REAL(DOUBLE)                            :: Op
    INTEGER                                 :: JP,P,Q,Q1,R,S,IL,KL,IG,JG,M,MN,MN1,  &
         IStrtA,IStopA,IStrtB,IStopB

    Q=1
    R=1
    Op=Zero
    CNAtms=ANAtms
    DO IL=1,ANAtms
      IG=IL+AOffSt
      CRowPt(IL)=R
      M=BSiz(IG)
      IStrtA=ARowPt(IL)
      IStopA=ARowPt(IL+1)-1
      IStrtB=BRowPt(IL)
      IStopB=BRowPt(IL+1)-1
      DO JP=IStrtA,IStopA
        JG=AColPt(JP)
        MN=M*BSiz(JG)
        MN1=MN-1
        P=ABlkPt(JP)+(ASMat-1)*MN !<<<SPIN
        Q1=Q+(CSMat-1)*MN !<<<SPIN
        CMTrix(Q1:Q1+MN1)=AMTrix(P:P+MN1)
        CColPt(R)=JG
        CBlkPt(R)=Q
        Flag(JG)=Q
        Q=Q+MN*NSMat !<<<SPIN
        R=R+1
      ENDDO
      DO JP=IStrtB,IStopB
        JG=BColPt(JP)
        MN=M*BSiz(JG)
        MN1=MN-1
        P=BBlkPt(JP)+(BSMat-1)*MN !<<<SPIN
        S=Flag(JG)
        IF(S/=0)THEN
          S=S+(CSMat-1)*MN !<<<SPIN
          CMTrix(S:S+MN1)=CMTrix(S:S+MN1)+BMTrix(P:P+MN1)
          Op=Op+MN
        ELSE
          IF(Q+MN>MxNon0.OR.R>MxBlks)THEN
            Add_GENERIC=-1
            RETURN
          ENDIF
          Q1=Q+(CSMat-1)*MN !<<<SPIN
          CMTrix(Q1:Q1+MN1)=BMTrix(P:P+MN1)
          CColPt(R)=JG
          CBlkPt(R)=Q
          Flag(JG)=Q
          Q=Q+MN*NSMat !<<<SPIN
          R=R+1
        ENDIF
      ENDDO
      DO KL=CRowPt(IL),R-1
        Flag(CColPt(KL))=0
      ENDDO
    ENDDO
    Add_GENERIC=SUCCEED
    CNBlks=R-1
    ! minus one
    CNNon0=Q-1
    CRowPt(CNAtms+1)=R
    FlOp=FlOp+Op
  END FUNCTION Add_GENERIC
#ifdef PARALLEL
  !===============================================================================
  !     DBCSR wrapper for generic matrix-scalar addition
  !===============================================================================
  SUBROUTINE Add_DBCSR_SCLR(A,B,Perf_O)
    IMPLICIT NONE
    TYPE(DBCSR),        INTENT(INOUT) :: A
    REAL(DOUBLE),       INTENT(IN)    :: B
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    !-------------------------------------------------------------------------------
    CALL Add_GENERIC_SCLR(A%NAtms,OffSt%I(MyId),A%RowPt%I,A%ColPt%I,  &
         A%BlkPt%I,A%MTrix%D,B,BSiz%I)
    PerfMon%FLOP=PerfMon%FLOP+DBLE(NBasF)
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+DBLE(NBasF)
  END SUBROUTINE Add_DBCSR_SCLR
#endif
  !===============================================================================
  !     BCSR wrapper for generic matrix-scalar addition
  !===============================================================================
  SUBROUTINE Add_BCSR_SCLR(A,B,Perf_O)
    IMPLICIT NONE
    TYPE(BCSR),         INTENT(INOUT) :: A
    REAL(DOUBLE),       INTENT(IN)    :: B
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    !-------------------------------------------------------------------------------
    CALL Add_GENERIC_SCLR(A%NAtms,0,A%RowPt%I,A%ColPt%I,  &
         A%BlkPt%I,A%MTrix%D,B,BSiz%I)
    PerfMon%FLOP=PerfMon%FLOP+DBLE(NBasF)
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+DBLE(NBasF)
  END SUBROUTINE Add_BCSR_SCLR
  !===============================================================================
  !     Generic F77 style (D)BCSR matrix-scalar addition
  !===============================================================================
  SUBROUTINE Add_GENERIC_SCLR(ANAtms,AOffSt,ARowPt,AColPt,ABlkPt,AMTrix,B,BSiz)
    IMPLICIT NONE
    INTEGER,                  INTENT(IN)    :: ANAtms,AOffSt
    INTEGER,     DIMENSION(:),INTENT(IN)    :: ARowPt,AColPt,ABlkPt,BSiz
    REAL(DOUBLE),             INTENT(IN)    :: B
    REAL(DOUBLE),DIMENSION(:),INTENT(INOUT) :: AMTrix
    INTEGER                                 :: I,IG,JP,Q,MA,IStrtA,IStopA
    !-------------------------------------------------------------------------------
    DO I=1,ANAtms
      IG=I+AOffSt
      MA=BSiz(IG)
      IStrtA=ARowPt(I)
      IStopA=ARowPt(I+1)-1
      DO JP=IStrtA,IStopA
        IF(AColPt(JP)==IG)THEN
          Q=ABlkPt(JP)
          CALL AddToDiag(MA,AMTrix(Q),B)
          EXIT
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE Add_GENERIC_SCLR
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !     MATRIX TRACE MATRIX TRACE MATRIX TRACE MATRIX TRACE MATRIX TRACE MATRIX TRACE MATRIX
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef PARALLEL
  !===============================================================================
  !     DBCSR trace of a matrix product Tr{A.B}=Sum_{ki} A_ki*B_ik
  !===============================================================================
  FUNCTION TraceMM_DBCSR(A,B,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: TraceMM_DBCSR
    TYPE(DBCSR),        INTENT(IN)    :: A,B
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: Tr,Op
    REAL(DOUBLE), EXTERNAL            :: BlkTrace_2
    INTEGER                           :: I,MA,J,JP,K,KP,P,Q,MB, &
         IStrtA,IStopA,IStrtB,IStopB
    !-------------------------------------------------------------------------------


    CALL New(GlobalRowPtA,NAtoms+1)
    CALL New(GlobalRowPtB,NAtoms+1)
    CALL SetEq(GlobalRowPtA,0)
    CALL SetEq(GlobalRowPtB,0)
    CALL Load(A,GlobalRowPtA)
    CALL Load(B,GlobalRowPtB)
    !-------------------------------------------------------------------------------
    Tr=Zero
    Op=Zero
    DO I=Beg%I(MyId),End%I(MyId)
      MA=BSiz%I(I)
      IStrtA=GlobalRowPtA%I(I)
      IStopA=GlobalRowPtA%I(I+1)-1
      DO JP=IStrtA,IStopA
        J=A%ColPt%I(JP)
        P=A%BlkPt%I(JP)
        MB=BSiz%I(J)
        IStrtB=GlobalRowPtB%I(J)
        IStopB=GlobalRowPtB%I(J+1)-1
        IF(IStrtB/=0.AND.IStopB/=0)THEN
          DO KP=IStrtB,IStopB
            K=B%ColPt%I(KP)
            IF(K.EQ.I)THEN
              Q=B%BlkPt%I(KP)
              Tr=Tr+BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              Op=Op+DBLE(MA*MB)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    TraceMM_DBCSR=AllReduce(Tr)
    !-----------------------------------------
    !        Count the muscle

    PerfMon%FLOP=PerfMon%FLOP+Two*Op
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+Two*Op
    !-------------------------------------------------------------------------------
    !        Tidy up ...

    CALL Delete(GlobalRowPtA)
    CALL Delete(GlobalRowPtB)

  END FUNCTION TraceMM_DBCSR
#endif
  !===============================================================================
  !     BCSR trace of a matrix product Tr{A.B}=Sum_{ki} A_ki*B_ik
  !===============================================================================
  FUNCTION TraceMM_BCSR(A,B,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: TraceMM_BCSR
    TYPE(BCSR),         INTENT(IN)    :: A,B
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: Op
    REAL(DOUBLE), EXTERNAL            :: BlkTrace_2
    INTEGER                           :: IL,MA,J,JP,K,KP,P,Q,MB, &
         IStrtA,IStopA,IStrtB,IStopB
    !-------------------------------------------------------------------------------
    Op=Zero
    TraceMM_BCSR=Zero
    DO IL=1,A%NAtms
      MA=BSiz%I(IL)
      IStrtA=A%RowPt%I(IL)
      IStopA=A%RowPt%I(IL+1)-1
      DO JP=IStrtA,IStopA
        J=A%ColPt%I(JP)
        P=A%BlkPt%I(JP)
        MB=BSiz%I(J)
        IStrtB=B%RowPt%I(J)
        IStopB=B%RowPt%I(J+1)-1
        DO KP=IStrtB,IStopB
          K=B%ColPt%I(KP)
          IF(K.EQ.IL)THEN
            Q=B%BlkPt%I(KP)
            IF(    A%NSMat.EQ.1.AND.B%NSMat.EQ.1) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              Op=Op+DBLE(MA*MB)
            ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.2) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q+MB*MA))
              Op=Op+DBLE(2*MA*MB)
            ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.1) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P+MB*MA),B%MTrix%D(Q))
              Op=Op+DBLE(2*MA*MB)
            ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P+MB*MA),B%MTrix%D(Q+MB*MA))
              Op=Op+DBLE(2*MA*MB)
            ELSEIF(A%NSMat.EQ.1.AND.B%NSMat.EQ.4) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q+3*MB*MA))
              Op=Op+DBLE(2*MA*MB)
            ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.1) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P+3*MB*MA),B%MTrix%D(Q))
              Op=Op+DBLE(2*MA*MB)
            ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4) THEN
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P),B%MTrix%D(Q))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P+MA*MB),B%MTrix%D(Q+2*MA*MB))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P+3*MA*MB),B%MTrix%D(Q+MA*MB))
              TraceMM_BCSR=TraceMM_BCSR  &
                   +BlkTrace_2(MA,MB,A%MTrix%D(P+3*MA*MB),B%MTrix%D(Q+3*MA*MB))
              Op=Op+DBLE(4*MA*MB)
            ELSE
              CALL Halt('TraceMM_BCSR 2: Something wrong there!')
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    PerfMon%FLOP=PerfMon%FLOP+Two*Op
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+Two*Op
  END FUNCTION TraceMM_BCSR
#ifdef PARALLEL
  !===============================================================================
  !     DBCSR trace of a matrix Tr{A}=Sum_i A_ii
  !===============================================================================
  FUNCTION Trace_DBCSR(A,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: Trace_DBCSR
    TYPE(DBCSR),        INTENT(INOUT)    :: A
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: Op,Tr
    REAL(DOUBLE),EXTERNAL             :: BlkTrace_1
    INTEGER                           :: IL,IG,MA,J,JP,P,IStrtA,IStopA
    !-------------------------------------------------------------------------------
    Op=Zero
    Tr=Zero
    DO IL=1,A%NAtms
      IG=IL+OffSt%I(MyId)
      MA=BSiz%I(IG)
      IStrtA=A%RowPt%I(IL)
      IStopA=A%RowPt%I(IL+1)-1
      DO JP=IStrtA,IStopA
        J=A%ColPt%I(JP)
        IF(J==IG)THEN
          P=A%BlkPt%I(JP)
          IF(    A%NSMat.EQ.1) THEN
            Tr=Tr+BlkTrace_1(MA,A%MTrix%D(P))
            Op=Op+DBLE(MA)
          ELSEIF(A%NSMat.EQ.2) THEN
            Tr=Tr+BlkTrace_1(MA,A%MTrix%D(P))
            Tr=Tr+BlkTrace_1(MA,A%MTrix%D(P+MA*MA))
            Op=Op+2D0*DBLE(MA)
          ELSEIF(A%NSMat.EQ.4) THEN
            Tr=Tr+BlkTrace_1(MA,A%MTrix%D(P))
            Tr=Tr+BlkTrace_1(MA,A%MTrix%D(P+3*MA*MA))
            Op=Op+2D0*DBLE(MA)
          ELSE
            CALL Halt('Trace_DBCSR: Something wrong!')
          ENDIF
          !P=A%BlkPt%I(JP)
          !Tr=Tr+BlkTrace_1(MA,A%MTrix%D(P))
          !Op=Op+DBLE(MA)
          EXIT
        ENDIF
      ENDDO
    ENDDO
    Trace_DBCSR=AllReduce(Tr)
    PerfMon%FLOP=PerfMon%FLOP+Op
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+Op
  END FUNCTION Trace_DBCSR
#endif
  !===============================================================================
  !     BCSR trace of a matrix Tr{A}=Sum_i A_ii
  !===============================================================================
  FUNCTION Trace_BCSR(A,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: Trace_BCSR
    TYPE(BCSR),         INTENT(IN)    :: A
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: Op
    REAL(DOUBLE), EXTERNAL            :: BlkTrace_1
    INTEGER                           :: IL,MA,J,JP,P,IStrtA,IStopA
    !-------------------------------------------------------------------------------
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      Op=Zero
      Trace_BCSR=Zero
      DO IL=1,A%NAtms
        MA=BSiz%I(IL)
        IStrtA=A%RowPt%I(IL)
        IStopA=A%RowPt%I(IL+1)-1
        DO JP=IStrtA,IStopA
          J=A%ColPt%I(JP)
          IF(J==IL)THEN
            P=A%BlkPt%I(JP)
            IF(    A%NSMat.EQ.1) THEN
              Trace_BCSR=Trace_BCSR+BlkTrace_1(MA,A%MTrix%D(P))
              Op=Op+DBLE(MA)
            ELSEIF(A%NSMat.EQ.2) THEN
              Trace_BCSR=Trace_BCSR+BlkTrace_1(MA,A%MTrix%D(P))
              Trace_BCSR=Trace_BCSR+BlkTrace_1(MA,A%MTrix%D(P+MA*MA))
              Op=Op+2D0*DBLE(MA)
            ELSEIF(A%NSMat.EQ.4) THEN
              Trace_BCSR=Trace_BCSR+BlkTrace_1(MA,A%MTrix%D(P))
              Trace_BCSR=Trace_BCSR+BlkTrace_1(MA,A%MTrix%D(P+3*MA*MA))
              Op=Op+2D0*DBLE(MA)
            ELSE
              CALL Halt('Trace_BCSR: Something wrong!')
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ENDDO
      PerfMon%FLOP=PerfMon%FLOP+Op
      IF(PRESENT(Perf_O))  &
           Perf_O%FLOP=Perf_O%FLOP+Op
#ifdef PARALLEL
    ENDIF
    !! write(*,*) 'warning from parallel: trace_bcsr not broadcasted!'
    !! CALL Bcast(Trace_BCSR)
#endif
  END FUNCTION Trace_BCSR

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !     MATRIX EQUALS IDENTITY MATRIX EQUALS IDENTITY MATRIX EQUALS IDENTITY MATRIX E
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE SetToI_BCSR(A,B,Expert_O)
    TYPE(BCSR)              :: A
    TYPE(DBL_RNK2),OPTIONAL :: B
    INTEGER       ,OPTIONAL :: Expert_O
    INTEGER                 :: I,N,N2,Q,R,Expert
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      Expert=0
      IF(PRESENT(Expert_O))Expert=Expert_O
      IF(.NOT.AllocQ(A%Alloc))CALL New(A)
      Q=1
      R=1
      A%NAtms=0
      A%RowPt%I(1)=1
      DO I=1,NAtoms
        A%NAtms=A%NAtms+1
        N=BSiz%I(I)
        N2=N*N
        A%ColPt%I(Q)=I
        A%BlkPt%I(Q)=R
        IF(Expert.EQ.0)THEN
          IF(A%NSMat.EQ.1)THEN
            IF(PRESENT(B))THEN
              A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
            ELSE
              CALL DiagI(N,A%MTrix%D(R:R+N2-1))
            ENDIF
          ELSEIF(A%NSMat.EQ.2) THEN
            IF(PRESENT(B))THEN
              A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
              R=R+N2
              A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
            ELSE
              CALL DiagI(N,A%MTrix%D(R:R+N2-1))
              R=R+N2
              CALL DiagI(N,A%MTrix%D(R:R+N2-1))
            ENDIF
          ELSEIF(A%NSMat.EQ.4) THEN
            IF(PRESENT(B))THEN
              A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
              R=R+N2
              A%MTrix%D(R:R+N2-1)=0D0
              R=R+N2
              A%MTrix%D(R:R+N2-1)=0D0
              R=R+N2
              A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
            ELSE
              CALL DiagI(N,A%MTrix%D(R:R+N2-1))
              R=R+N2
              A%MTrix%D(R:R+N2-1)=0D0
              R=R+N2
              A%MTrix%D(R:R+N2-1)=0D0
              R=R+N2
              CALL DiagI(N,A%MTrix%D(R:R+N2-1))
            ENDIF
          ELSE
            CALL Halt('SetToI_BCSR: wrong value for A%NSMat')
          ENDIF
          R=R+N2
        ELSEIF(Expert.GT.0.AND.Expert.LT.5)THEN
          !The expert SetToI
          IF(PRESENT(B))THEN
            R=R+(Expert-1)*N2
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
          ELSE
            R=R+(Expert-1)*N2
            CALL DiagI(N,A%MTrix%D(R:R+N2-1))
          ENDIF
          R=R+(A%NSMat-Expert+1)*N2
        ELSE
          CALL Halt('SetToI_BCSR: wrong value for Expert')
        ENDIF
        Q=Q+1
        A%RowPt%I(A%NAtms+1)=Q
      ENDDO
      A%NBlks=Q-1
      A%NNon0=R-1
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE SetToI_BCSR
#ifdef PARALLEL
  SUBROUTINE SetToI_DBCSR(A,B,Expert_O)
    TYPE(DBCSR)             :: A
    TYPE(DBL_RNK2),OPTIONAL :: B
    INTEGER       ,OPTIONAL :: Expert_O
    INTEGER                 :: I,N,N2,Q,R,Expert
    Expert=0
    IF(PRESENT(Expert_O))Expert=Expert_O
    IF(.NOT.AllocQ(A%Alloc))CALL New(A)
    Q=1
    R=1
    A%NAtms=0
    A%RowPt%I(1)=1
    DO I=Beg%I(MyId),End%I(MyId)
      A%NAtms=A%NAtms+1
      N=BSiz%I(I)
      N2=N*N
      A%ColPt%I(Q)=I
      A%BlkPt%I(Q)=R
      !<<<
      IF(Expert.EQ.0)THEN
        IF(A%NSMat.EQ.1)THEN
          IF(PRESENT(B))THEN
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
          ELSE
            CALL DiagI(N,A%MTrix%D(R:R+N2-1))
          ENDIF
        ELSEIF(A%NSMat.EQ.2) THEN
          IF(PRESENT(B))THEN
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
            R=R+N2
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
          ELSE
            CALL DiagI(N,A%MTrix%D(R:R+N2-1))
            R=R+N2
            CALL DiagI(N,A%MTrix%D(R:R+N2-1))
          ENDIF
        ELSEIF(A%NSMat.EQ.4) THEN
          IF(PRESENT(B))THEN
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
            R=R+N2
            A%MTrix%D(R:R+N2-1)=0D0
            R=R+N2
            A%MTrix%D(R:R+N2-1)=0D0
            R=R+N2
            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
          ELSE
            CALL DiagI(N,A%MTrix%D(R:R+N2-1))
            R=R+N2
            A%MTrix%D(R:R+N2-1)=0D0
            R=R+N2
            A%MTrix%D(R:R+N2-1)=0D0
            R=R+N2
            CALL DiagI(N,A%MTrix%D(R:R+N2-1))
          ENDIF
        ELSE
          CALL Halt('SetToI_BCSR: wrong value for A%NSMat')
        ENDIF
        R=R+N2
      ELSEIF(Expert.GT.0.AND.Expert.LT.5)THEN
        !The expert SetToI
        IF(PRESENT(B))THEN
          R=R+(Expert-1)*N2
          A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
        ELSE
          R=R+(Expert-1)*N2
          CALL DiagI(N,A%MTrix%D(R:R+N2-1))
        ENDIF
        R=R+(A%NSMat-Expert+1)*N2
      ELSE
        CALL Halt('SetToI_BCSR: wrong value for Expert')
      ENDIF
      Q=Q+1
      !>>>
!!$         IF(PRESENT(B))THEN
!!$            A%MTrix%D(R:R+N2-1)=B%D(1:N2,I)
!!$         ELSE
!!$           CALL DiagI(N,A%MTrix%D(R:R+N2-1))
!!$         ENDIF
!!$         Q=Q+1
!!$         R=R+N2
      !<<<
      A%RowPt%I(A%NAtms+1)=Q
    ENDDO
    A%NBlks=Q-1
    A%NNon0=R-1
    CALL LocalToGlobal(A)
  END SUBROUTINE SetToI_DBCSR
#endif
  SUBROUTINE DiagI(N,A)
    INTEGER                      :: I,N
    REAL(DOUBLE), DIMENSION(N,N) :: A
    A=Zero
    DO I=1,N
      A(I,I)=One
    ENDDO
  END SUBROUTINE DiagI
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !     MATRIX DOT PRODUCT MATRIX DOT PRODUCT MATRIX DOT PRODUCT MATRIX DOT PRODUCT
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef PARALLEL
  !===============================================================================
  !     Wrapper for generic F77 style DBCSR matrix inner product
  !===============================================================================
  FUNCTION Dot_DBCSR(A,B,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: Dot_DBCSR
    TYPE(DBCSR),        INTENT(INOUT)    :: A,B
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: Tmp,FlOp

    FlOp=Zero
    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
    Tmp=Dot_GENERIC(A%NSMat,A%NAtms,OffSt%I(MyId),BSiz%I,Flag%I,FlOp,  &
         A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D,   &
         B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D)
    Dot_DBCSR=AllReduce(Tmp)
    CALL Delete(Flag)
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+FlOp
  END FUNCTION Dot_DBCSR
#endif
  !===============================================================================
  !     Wrapper for generic F77 style BCSR matrix inner product
  !===============================================================================
  FUNCTION Dot_BCSR(A,B,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: Dot_BCSR
    TYPE(BCSR),         INTENT(IN)    :: A,B
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: FlOp
    TYPE(INT_VECT)                    :: Flag

    CALL Initialize(Flag)

    !write(*,*) 'Dot_BCSR: ',A%NSMat,B%NSMat
    IF(A%NSMat.NE.B%NSMat)CALL Halt('Dot_BCSR: A%NSMat.NE.B%NSMat!')
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      FlOp=Zero
      CALL New(Flag,NAtoms)
      CALL SetEq(Flag,0)
      Dot_BCSR=                                            &
           Dot_GENERIC(A%NSMat,A%NAtms,0,BSiz%I,Flag%I,FlOp,    &
           A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
           B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D)
      CALL Delete(Flag)
      PerfMon%FLOP=PerfMon%FLOP+FlOp
      IF(PRESENT(Perf_O))  &
           Perf_O%FLOP=Perf_O%FLOP+FlOp
#ifdef PARALLEL
    ENDIF
#endif
  END FUNCTION Dot_BCSR

  !===============================================================================
  !     Generic F77 style (D,B)CSR inner product of two matrices
  !     MatrixDot = (A,B)
  !===============================================================================
  FUNCTION Dot_GENERIC(NSMat,ANAtms,OffStA,BSiz,Flag,FlOp, &
       ARowPt,AColPt,ABlkPt,AMTrix,  &
       BRowPt,BColPt,BBlkPt,BMTrix)
    IMPLICIT NONE
    REAL(DOUBLE)                            :: Dot_GENERIC
    INTEGER,                  INTENT(IN)    :: NSMat,ANAtms,OffStA
    INTEGER,     DIMENSION(:),INTENT(IN)    :: ARowPt,AColPt,ABlkPt, &
         BRowPt,BColPt,BBlkPt,BSiz
    REAL(DOUBLE),DIMENSION(:),INTENT(IN)    :: AMTrix,BMTrix
    INTEGER,     DIMENSION(:),INTENT(INOUT) :: Flag
    REAL(DOUBLE),             INTENT(INOUT) :: FlOp
    REAL(DOUBLE)                            :: Op
    REAL(DOUBLE), EXTERNAL                  :: DBL_Dot
    INTEGER                                 :: I,IG,J,L,JP,P,Q,MA,MN, &
         IStrtA,IStopA,IStrtB,IStopB

    Op=Zero
    Dot_GENERIC=Zero
    DO I=1,ANAtms
      IG=I+OffStA
      MA=BSiz(IG)
      IStrtA=ARowPt(I)
      IStopA=ARowPt(I+1)-1
      IStrtB=BRowPt(I)
      IStopB=BRowPt(I+1)-1
      DO JP=IStrtB,IStopB
        Flag(BColPt(JP))=JP
      ENDDO
      DO JP=IStrtA,IStopA
        J=AColPt(JP)
        L=Flag(J )
        IF(L/=0)THEN
          P=ABlkPt(JP)
          Q=BBlkPt(L )
          MN=MA*BSiz(J)*NSMat!<<< SPIN
          Dot_GENERIC=Dot_GENERIC &
               +DBL_DOT(MN,AMTrix(P),BMTrix(Q))
          Op=Op+DBLE(MN)
        ENDIF
      ENDDO
      DO JP=IStrtB,IStopB
        Flag(BColPt(JP))=0
      ENDDO
    ENDDO
    FlOp=FlOp+Two*Op
  END FUNCTION Dot_GENERIC
  !===============================================================================
  !     Take the dot (inner) product of DBL_VECTs
  !===============================================================================
  FUNCTION Dot_DBL_VECT(A,B,N_O,Perf_O)
    REAL(DOUBLE)                      :: Dot_DBL_VECT
    TYPE(DBL_VECT),     INTENT(IN)    :: A,B
    INTEGER,   OPTIONAL,INTENT(IN)    :: N_O
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE), EXTERNAL            :: DBL_Dot
    INTEGER                           :: N
    N=SIZE(A%D); IF(PRESENT(N_O))N=N_O
    IF(SIZE(B%D)<N)CALL Halt('Dimensions in Dot_DBL_VECT')
    Dot_DBL_VECT=DBL_DOT(N,A%D(1),B%D(1))
    PerfMon%FLOP=PerfMon%FLOP+Two*DBLE(N)
    IF(PRESENT(Perf_O))  &
         Perf_O%FLOP=Perf_O%FLOP+Two*DBLE(N)
  END FUNCTION Dot_DBL_VECT
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !     MATRIX FILTRATION, MATRIX FILTRATION, MATRIX FILTRATION, MATRIX FILTRATION
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#ifdef PARALLEL
  !===============================================================================
  !     Wrapper for generic F77 style DBCSR matrix filtration
  !===============================================================================
  SUBROUTINE FilterM_DBCSR(A,B,Tol_O,Perf_O,SetEq_O)
    REAL(DOUBLE),OPTIONAL,INTENT(IN)    :: Tol_O
    LOGICAL,     OPTIONAL,INTENT(IN)    :: SetEq_O
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    TYPE(DBCSR),          INTENT(INOUT) :: A
    TYPE(DBCSR),          INTENT(INOUT) :: B
    REAL(DOUBLE)                        :: Tol,FlOp
    !-------------------------------------------------------------------------------
    IF(.NOT.AllocQ(A%Alloc))CALL New(A,NSMat_O=B%NSMat)
    FlOp=Zero
    IF(PRESENT(Tol_O))THEN
      Tol=Tol_O
    ELSE
      Tol=Thresholds%Trix
    ENDIF
    IF(PRESENT(SetEq_O))THEN
      IF(SetEq_O)THEN
        CALL SetEq(A,B)
        RETURN
      ENDIF
    ENDIF
    CALL FilterM_GENERIC(B%NSMat,B%NAtms,OffSt%I(MyId),           & !<<< SPIN
         A%NAtms,A%NBlks,A%NNon0,                 &
         A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
         BSiz%I,Tol,FlOp)
    A%Node=MyId
    A%GUpDate=STATUS_FALSE
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE FilterM_DBCSR
  !===============================================================================
  !     Gather local indecies to global ones, bcast global indecies
  !===============================================================================
  SUBROUTINE LocalToGlobal2(A)
    TYPE(DBCSR), INTENT(INOUT) :: A
    INTEGER                    :: I,K,NAtms
    TYPE(INT_VECT)             :: MA,MB,MN,NA,NB,NN
    !        Allocate offset indecies
    IF(MyId==ROOT)THEN
      CALL New(MA,NPrc,M_O=0)
      CALL New(MB,NPrc,M_O=0)
      CALL New(MN,NPrc,M_O=0)
      CALL New(NA,NPrc,M_O=0)
      CALL New(NB,NPrc,M_O=0)
      CALL New(NN,NPrc,M_O=0)
    ENDIF
    !        Number of atoms per node
    NAtms=End%I(MyId)-Beg%I(MyId)+1
    IF(MyID==NPrc-1)NAtms=NAtms+1
    !        Gather the indecies to root
    CALL Gather(NAtms  ,MA)
    CALL Gather(A%NBlks,MB)
    CALL Gather(A%NNon0,MN)
    !        Calculate index offsets (displacements)
    IF(MyID==ROOT)THEN
      NA%I(0)=  NAtms
      NB%I(0)=A%NBlks
      NN%I(0)=A%NNon0
      DO I=1,NPrc-1
        NA%I(I)=MA%I(I)+NA%I(I-1)
        NB%I(I)=MB%I(I)+NB%I(I-1)
        NN%I(I)=MN%I(I)+NN%I(I-1)
      ENDDO
      DO I=NPrc,1,-1
        NA%I(I)=NA%I(I-1)
        NB%I(I)=NB%I(I-1)
        NN%I(I)=NN%I(I-1)
      ENDDO
      NA%I(0)=0
      NB%I(0)=0
      NN%I(0)=0
    ENDIF
    !        Gather the global row-col pointers to ROOT
    CALL Gather(A%RowPt,A%GRwPt,  NAtms,MA,NA)
    CALL Gather(A%ColPt,A%GClPt,A%NBlks,MB,NB)
    !        Add offsets to achieve correct indexing
    IF(MyId==ROOT)THEN
      DO I=1,NPrc-1
        DO K=NA%I(I)+1,NA%I(I+1)
          A%GRwPt%I(K)=A%GRwPt%I(K)+NB%I(I)
        ENDDO
      ENDDO
    ENDIF
    CALL Bcast(A%GRwPt)
    CALL Bcast(A%GClPt)
    !-------------------------------------------------------------------------------
    !        Tidy up

    IF(MyId==ROOT)THEN
      CALL Delete(MA); CALL Delete(MB); CALL Delete(MN)
      CALL Delete(NA); CALL Delete(NB); CALL Delete(NN)
    ENDIF
  END SUBROUTINE LocalToGlobal2

  SUBROUTINE LocalToGlobal(A)
    TYPE(DBCSR),     INTENT(INOUT) :: A
    INTEGER                        :: I,K,NAtms,IErr
    TYPE(INT_VECT)                 :: MA,MB,MN,NA,NB,NN

    IF(A%GUpDate==STATUS_TRUE)THEN
      RETURN
    ELSE
      A%GUpDate=STATUS_TRUE
    ENDIF

    !         CALL LocalToGlobal2(A)
    !         RETURN

    CALL New(MA,NPrc,M_O=0)
    CALL New(MB,NPrc,M_O=0)
    CALL New(MN,NPrc,M_O=0)
    CALL New(NA,NPrc,M_O=0)
    CALL New(NB,NPrc,M_O=0)
    CALL New(NN,NPrc,M_O=0)
    !        Number of atoms per node
    NAtms=End%I(MyId)-Beg%I(MyId)+1
    IF(MyID==NPrc-1)NAtms=NAtms+1
    !        Gather the indecies to root
    MA%I=0; MB%I=0; MN%I=0
    CALL MPI_ALLGATHER(NAtms,1,MPI_INTEGER, &
         MA%I, 1,MPI_INTEGER, &
         MONDO_COMM,IErr)
    CALL ErrChk(IErr,'LocalToGlobal#1')
    CALL MPI_ALLGATHER(A%NBlks,1,MPI_INTEGER, &
         MB%I,   1,MPI_INTEGER, &
         MONDO_COMM,IErr)
    CALL ErrChk(IErr,'LocalToGlobal#2')
    CALL MPI_ALLGATHER(A%NNon0,1,MPI_INTEGER, &
         MN%I,   1,MPI_INTEGER, &
         MONDO_COMM,IErr)
    CALL ErrChk(IErr,'LocalToGlobal#3')
    !         CALL PPrint(MA,' 1 MA I',M_O=0,N_O=NPrc-1,Unit_O=6)
    !         CALL PPrint(MB,' 1 MB I',M_O=0,N_O=NPrc-1,Unit_O=6)
    !         CALL PPrint(MN,' 1 MN I',M_O=0,N_O=NPrc-1,Unit_O=6)
    !        Calculate index offsets (displacements)
    NA%I(0)=MA%I(0)
    NB%I(0)=MB%I(0)
    NN%I(0)=MN%I(0)
    DO I=1,NPrc-1
      NA%I(I)=MA%I(I)+NA%I(I-1)
      NB%I(I)=MB%I(I)+NB%I(I-1)
      NN%I(I)=MN%I(I)+NN%I(I-1)
    ENDDO
    DO I=NPrc,1,-1
      NA%I(I)=NA%I(I-1)
      NB%I(I)=NB%I(I-1)
      NN%I(I)=NN%I(I-1)
    ENDDO
    NA%I(0)=0
    NB%I(0)=0
    NN%I(0)=0
    !         CALL PPrint(NA,' 1 NA II',M_O=0,N_O=NPrc,Unit_O=6)
    !         CALL PPrint(NB,' 1 NB II',M_O=0,N_O=NPrc,Unit_O=6)
    !         CALL PPrint(NN,' 1 NN II',M_O=0,N_O=NPrc,Unit_O=6)
    !         CALL Halt(' Localto global ')
    CALL MPI_ALLGATHERV(A%RowPt%I,NAtms,    MPI_INTEGER,  &
         A%GRwPt%I,MA%I,NA%I,MPI_INTEGER,  &
         MONDO_COMM,IErr)
    DO I=1,NPrc-1
      DO K=NA%I(I)+1,NA%I(I+1)
        A%GRwPt%I(K)=A%GRwPt%I(K)+NB%I(I)
      ENDDO
    ENDDO
    CALL ErrChk(IErr,'LocalToGlobal#4')
    CALL MPI_ALLGATHERV(A%ColPt%I,A%NBlks,  MPI_INTEGER,  &
         A%GClPt%I,MB%I,NB%I,MPI_INTEGER,  &
         MONDO_COMM,IErr)
    CALL ErrChk(IErr,'LocalToGlobal#5')
    !        Tidy up
    CALL Delete(MA); CALL Delete(MB); CALL Delete(MN)
    CALL Delete(NA); CALL Delete(NB); CALL Delete(NN)
  END SUBROUTINE LocalToGlobal
#endif
  !===============================================================================
  !     Wrapper for generic F77 style BCSR matrix filtration
  !===============================================================================
  SUBROUTINE FilterM_BCSR(A,B,Tol_O,Perf_O,SetEq_O)
    TYPE(BCSR),           INTENT(INOUT) :: A,B
    REAL(DOUBLE),OPTIONAL,INTENT(IN)    :: Tol_O
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                        :: Tol,FlOp
    LOGICAL,     OPTIONAL,INTENT(IN)    :: SetEq_O
    !-------------------------------------------------------------------------------
    IF(.NOT.AllocQ(A%Alloc))CALL New(A,NSMat_O=B%NSMat)
    FlOp=Zero
    IF(PRESENT(Tol_O))THEN
      Tol=Tol_O
    ELSE
      Tol=Thresholds%Trix
    ENDIF
    IF(PRESENT(SetEq_O))THEN
      IF(SetEq_O)THEN
        CALL SetEq(A,B)
        RETURN
      ENDIF
    ENDIF

    !write(*,*) 'LinAlg',A%NSMat,B%NSMat

#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      CALL FilterM_GENERIC(B%NSMat,B%NAtms,0,A%NAtms,A%NBlks,A%NNon0, &
           A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
           B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
           BSiz%I,Tol,FlOp)
#ifdef PARALLEL
    ENDIF
#endif
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE FilterM_BCSR
  !===============================================================================
  !     Wrapper for generic F77 style BCSR matrix filtration
  !===============================================================================
  SUBROUTINE FilterM_InPlace_BCSR(A,Tol_O,Perf_O)
    TYPE(BCSR),           INTENT(INOUT) :: A
    REAL(DOUBLE),OPTIONAL,INTENT(IN)    :: Tol_O
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                        :: Tol,FlOp
    !-------------------------------------------------------------------------------
    FlOp=Zero
    IF(PRESENT(Tol_O))THEN
      Tol=Tol_O
    ELSE
      Tol=Thresholds%Trix
    ENDIF
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      CALL FilterM_InPlace_GENERIC(0,A%NAtms,A%NBlks,A%NNon0, &
           A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
           BSiz%I,Tol,FlOp)
#ifdef PARALLEL
    ENDIF
#endif
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE FilterM_InPlace_BCSR
#ifdef PARALLEL
  !===============================================================================
  !     Wrapper for generic F77 style DBCSR matrix filtration
  !===============================================================================
  SUBROUTINE FilterM_InPlace_DBCSR(A,Tol_O,Perf_O)
    REAL(DOUBLE),OPTIONAL,INTENT(IN)    :: Tol_O
    TYPE(TIME),  OPTIONAL,INTENT(INOUT) :: Perf_O
    TYPE(DBCSR),          INTENT(INOUT) :: A
    REAL(DOUBLE)                        :: Tol,FlOp
    !-------------------------------------------------------------------------------
    FlOp=Zero
    IF(PRESENT(Tol_O))THEN
      Tol=Tol_O
    ELSE
      Tol=Thresholds%Trix
    ENDIF
    CALL FilterM_InPlace_GENERIC(OffSt%I(MyId),A%NAtms,A%NBlks,A%NNon0, &
         A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         BSiz%I,Tol,FlOp)
    A%Node=MyId
    A%GUpDate=STATUS_FALSE
    PerfMon%FLOP=PerfMon%FLOP+FlOp
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+FlOp
  END SUBROUTINE FilterM_InPlace_DBCSR
#endif
  !===============================================================================
  !     Generic F77 style routine for filtering (D/B)CSR matrices:
  !     A=FILTER(B,Tol)
  !===============================================================================
  SUBROUTINE FilterM_InPlace_GENERIC(BOffSt,ANAtms,ANBlks,ANNon0, &
       ARowPt,AColPt,ABlkPt,AMTrix, &
       MSiz,Tol,FlOp)
    INTEGER,                  INTENT(IN)  :: BOffSt
    INTEGER,     DIMENSION(:),INTENT(IN)  :: MSiz
    REAL(DOUBLE),             INTENT(IN)  :: Tol
    INTEGER,                  INTENT(INOUT) :: ANAtms,ANBlks,ANNon0
    INTEGER,     DIMENSION(:),INTENT(INOUT) :: ARowPt,AColPt,ABlkPt
    REAL(DOUBLE),DIMENSION(:),INTENT(INOUT) :: AMTrix
    REAL(DOUBLE),             INTENT(INOUT) :: FlOp
    REAL(DOUBLE),EXTERNAL                 :: DBL_Dot
    REAL(DOUBLE)                          :: Op
    INTEGER                               :: IL,IG,J,JP,K,P,Q,MA,NA,MN,IStrtA,IStopA,iMvAtm,iMvBlk,iS
    !-------------------------------------------------------------------------------
    iMvAtm=0;iMvBlk=0;K=0;Q=0;
    Op=Zero
    DO IL=1,ANAtms
      IG=IL+BOffSt
      MA=MSiz(IG)
      IStrtA=ARowPt(IL)
      IStopA=ARowPt(IL+1)-1
      ARowPt(IL)=ARowPt(IL)-iMvAtm
      IF(IStrtA.NE.0.AND.IStopA.NE.0)THEN
        ! Reorder the ABlkPt
        ! Here we suppose that ABlkPt(Row)<ABlkPt(Row+1)!
        CALL IntIntSort77(IStopA-IStrtA+1,ABlkPt(IStrtA),AColPt(IStrtA),1)
        DO JP=IStrtA,IStopA
          J=AColPt(JP)
          P=ABlkPt(JP)
          NA=MSiz(J)
          MN=MA*NA
          IF(SQRT(DBL_Dot(MN,AMTrix(P),AMTrix(P))).LT.Tol)THEN
            AColPt(JP)=-9999
            ABlkPt(JP)=-9999
            CALL DBL_VECT_EQ_DBL_SCLR(MN,AMTrix(P),BIG_DBL)
            iMvAtm=iMvAtm+1
            iMvBlk=iMvBlk+MN
          ELSE
            IF(iMvAtm.NE.0) THEN
              AColPt(JP-iMvAtm)=J
              ABlkPt(JP-iMvAtm)=P-iMvBlk
              CALL DBL_VECT_EQ_DBL_VECT(MN,AMTrix(P-iMvBlk),AMTrix(P))
            ENDIF
            K=K+1
            Q=Q+MN
          ENDIF
          Op=Op+DBLE(MN)
        ENDDO
      ENDIF
    ENDDO
    ARowPt(ANAtms+1)=ARowPt(ANAtms+1)-iMvAtm
    IF(K.GT.ANBlks.OR.Q.GT.ANNon0) THEN
      WRITE(*,*) 'WRONG LOGIC IN Filter_InPlace 1!'
      WRITE(*,*) 'ANBlks',ANBlks,' K',K
      WRITE(*,*) 'ANNon0',ANNon0,' Q',Q
      STOP
    ENDIF
    ANBlks=K
    ANNon0=Q
    !
    iS=SIZE(AMTrix)
    CALL DBL_VECT_EQ_DBL_SCLR(iS-ANNon0,AMTrix(ANNon0+1),-BIG_DBL)
    iS=SIZE(AColPt)
    CALL INT_VECT_EQ_INT_SCLR(iS-ANBlks,AColPt(ANBlks+1),BIG_INT)
    CALL INT_VECT_EQ_INT_SCLR(iS-ANBlks,ABlkPt(ANBlks+1),BIG_INT)
    !
    DO JP=1,ANBlks
      IF(AColPt(JP).EQ.-9999.OR.ABlkPt(JP).EQ.-9999) THEN
        WRITE(*,*) 'WRONG LOGIC IN Filter_InPlace 2!'
        WRITE(*,*) 'MyID',MyID
        WRITE(*,*) 'ANBlks',ANBlks
        WRITE(*,*) 'AColPt(',JP,')=',AColPt(JP)
        WRITE(*,*) 'ABlkPt(',JP,')=',ABlkPt(JP)
        STOP
      ENDIF
    ENDDO
    DO P=1,ANNon0
      IF(ABS(AMTrix(P)).GT.1D15) THEN
        WRITE(*,*) 'WRONG LOGIC IN Filter_InPlace 3!'
        WRITE(*,*) 'MyID',MyID
        WRITE(*,*) 'ANNon0',ANNon0
        WRITE(*,*) 'AMTrix(',P,')=',AMTrix(P)
        STOP
      ENDIF
    ENDDO
    FlOp=FlOp+Two*Op
    !WRITE(*,*) 'iMvAtm=',iMvAtm,' iMvBlk=',iMvBlk
    !
  END SUBROUTINE FilterM_InPlace_GENERIC



  !===============================================================================
  !     Generic F77 style routine for filtering (D/B)CSR matrices:
  !     A=FILTER(B,Tol)
  !===============================================================================
  SUBROUTINE FilterM_GENERIC(NSMat,BNAtms,BOffSt,         &
       ANAtms,ANBlks,ANNon0,        &
       ARowPt,AColPt,ABlkPt,AMTrix, &
       BRowPt,BColPt,BBlkPt,BMTrix, &
       MSiz,Tol,FlOp)
    INTEGER,                  INTENT(IN)  :: NSMat,BNAtms,BOffSt
    INTEGER,     DIMENSION(:),INTENT(IN)  :: BRowPt,BColPt,BBlkPt,MSiz
    REAL(DOUBLE),DIMENSION(:),INTENT(IN)  :: BMTrix
    REAL(DOUBLE),             INTENT(IN)  :: Tol
    INTEGER,                  INTENT(INOUT) :: ANAtms,ANBlks,ANNon0
    INTEGER,     DIMENSION(:),INTENT(INOUT) :: ARowPt,AColPt,ABlkPt
    REAL(DOUBLE),DIMENSION(:),INTENT(INOUT) :: AMTrix
    REAL(DOUBLE),             INTENT(INOUT) :: FlOp
    REAL(DOUBLE),EXTERNAL                 :: DBL_Dot
    REAL(DOUBLE)                          :: Op
    INTEGER                               :: IL,IG,J,JP,K,P,Q,MA,NA,MN,MN1, &
         IStrtB,IStopB
    !-------------------------------------------------------------------------------
    K=1
    Q=1
    Op=Zero
    ANAtms=BNAtms
    DO IL=1,BNAtms
      ARowPt(IL)=K
      IG=IL+BOffSt
      MA=MSiz(IG)
      IStrtB=BRowPt(IL)
      IStopB=BRowPt(IL+1)-1
      IF(IStrtB/=0.AND.IStopB/=0)THEN
        DO JP=IStrtB,IStopB
          J=BColPt(JP)
          P=BBlkPt(JP)
          NA=MSiz(J)
          MN=MA*NA*NSMat!<< SPIN
          MN1=MN-1
          IF(SQRT(DBL_Dot(MN,BMTrix(P),BMTrix(P)))>Tol)THEN
            !                  IF(DSQRT(DDot(MN,BMTrix(P),1,BMTrix(P),1))>Tol)THEN
            !                  IF(FNorm(MN,BMTrix(P:))>Tol)THEN
            CALL DBL_VECT_EQ_DBL_VECT(MN,AMTrix(Q),BMTrix(P))
            AColPt(K)=J
            ABlkPt(K)=Q
            K=K+1
            Q=Q+MN
          ENDIF
          Op=Op+DBLE(MN)
        ENDDO
      ENDIF
    ENDDO
    ARowPt(ANAtms+1)=K
    ANBlks=K-1
    ANNon0=Q-1
    FlOp=FlOp+Two*Op
  END SUBROUTINE FilterM_GENERIC

  !===============================================================================
  SUBROUTINE XPose_Simple(M,N,A,AT)
    INTEGER                      :: I,J,M,N,IDex,JDex
    REAL(DOUBLE), DIMENSION(M*N) :: A,AT
    DO I=1,N
      DO J=1,M
        IDex=(I-1)*M+J
        JDex=(J-1)*N+I
        AT(JDex)=A(IDex)
      ENDDO
    ENDDO
  END SUBROUTINE XPose_Simple

  !===============================================================================
  !     Max block of a BCSR matrix
  !===============================================================================
  SUBROUTINE XPose_BCSR(A,B)
    TYPE(BCSR)           :: A
    TYPE(BCSR), OPTIONAL :: B
    TYPE(INT_VECT)       :: RowPt,ColPt,BlkPt

    IF(PRESENT(B))THEN
      IF(.NOT.AllocQ(B%Alloc))  &
           CALL New(B,NSMat_O=A%NSMat)
      IF(A%NSMat.NE.B%NSMat) THEN
        CALL Delete(B)
        CALL New(B,NSMat_O=A%NSMat)
      ENDIF
      CALL XPose_GENERIC(A%NSMat,A%RowPt%I,A%ColPt%I,A%BlkPt%I, &
           B%RowPt%I,B%ColPt%I,B%BlkPt%I,A%MTrix%D,B%MTrix%D)
      B%NAtms=A%NAtms
      B%NBlks=A%NBlks
      B%NNon0=A%NNon0
    ELSE
      CALL New(RowPt,NAtoms+1)
      CALL New(ColPt,A%NBlks)
      CALL New(BlkPt,A%NBlks)
      CALL XPose_GENERIC(A%NSMat,A%RowPt%I,A%ColPt%I,A%BlkPt%I, &
           RowPt%I,ColPt%I,BlkPt%I)
      A%RowPt%I(1:NAtoms+1)=RowPt%I(1:NAtoms+1)
      A%ColPt%I(1:A%NBlks)=ColPt%I(1:A%NBlks)
      A%BlkPt%I(1:A%NBlks)=BlkPt%I(1:A%NBlks)
      CALL Delete(RowPt)
      CALL Delete(ColPt)
      CALL Delete(BlkPt)
    ENDIF

  END SUBROUTINE XPose_BCSR

  SUBROUTINE XPose_GENERIC(ANSMat,ARowPt,AColPt,ABlkPt,BRowPt,BColPt,BBlkPt,AMTrix,BMTrix)
    REAL(DOUBLE), DIMENSION(:), OPTIONAL :: AMTrix,BMTrix
    INTEGER,      DIMENSION(:)           :: ARowPt,AColPt,ABlkPt, &
         BRowPt,BColPt,BBlkPt
    INTEGER                              :: ANSMat,I,J,K,MI,NJ,MN,MN1,P,Q,NEXT
    LOGICAL                              :: SymbolicOnly

    IF(PRESENT(AMTrix).AND.PRESENT(BMTrix))THEN
      SymbolicOnly=.FALSE.
    ELSE
      SymbolicOnly=.TRUE.
    ENDIF

    DO I=1,NAtoms
      BRowPt(I)=0
    ENDDO
    DO I=1,NAtoms
      IF((ARowPt(I).NE.0).AND.(ARowPt(I+1)-1.NE.0))THEN
        DO K=ARowPt(I),ARowPt(I+1)-1
          J=AColPt(K)+1
          BRowPt(J)=BRowPt(J)+1
        ENDDO
      ENDIF
    ENDDO
    BRowPt(1)=1
    DO I=1,NAtoms
      BRowPt(I+1)=BRowPt(I)+BRowPt(I+1)
    ENDDO
    Q=1
    IF(SymbolicOnly)THEN
      DO I=1,NAtoms
        MI=BSiz%I(I)
        IF((ARowPt(I).NE.0).AND.(ARowPt(I+1)-1.NE.0))THEN
          DO K=ARowPt(I),ARowPt(I+1)-1
            J=AColPt(K)
            P=ABlkPt(K)
            NJ=BSiz%I(J)
            NEXT=BRowPt(J)
            MN=MI*NJ
            BBlkPt(NEXT)=P
            BColPt(NEXT)=I
            BRowPt(J)=NEXT+1
            Q=Q+MN
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO I=1,NAtoms
        MI=BSiz%I(I)
        IF((ARowPt(I).NE.0).AND.(ARowPt(I+1)-1.NE.0))THEN
          DO K=ARowPt(I),ARowPt(I+1)-1
            J=AColPt(K)
            P=ABlkPt(K)
            NJ=BSiz%I(J)
            NEXT=BRowPt(J)
            MN=MI*NJ
            MN1=MN-1
            IF(ANSMat.EQ.1)THEN
              CALL XPoseSqMat(MI,NJ,AMTrix(P:P+MN1),BMTrix(Q:Q+MN1))
            ELSEIF(ANSMat.EQ.2)THEN
              CALL XPoseSqMat(MI,NJ,AMTrix(P:P+MN1),BMTrix(Q:Q+MN1))
              CALL XPoseSqMat(MI,NJ,AMTrix(P+MN:P+MN+MN1),BMTrix(Q+MN:Q+MN+MN1))
              MN=2*MN
            ELSE
              CALL Halt('Error in XPose_GENERIC')
            ENDIF

            BBlkPt(NEXT)=Q
            BColPt(NEXT)=I
            BRowPt(J)=NEXT+1
            Q=Q+MN
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    DO I=NAtoms,1,-1
      BRowPt(I+1)=BRowPt(I)
    ENDDO
    BRowPt(1)=1
  END SUBROUTINE XPose_GENERIC

  SUBROUTINE XPoseSqMat(M,N,A,AT)
    INTEGER,                     INTENT(IN)  :: M,N
    REAL(DOUBLE), DIMENSION(M,N),INTENT(IN)  :: A
    REAL(DOUBLE), DIMENSION(N,M),INTENT(OUT) :: AT
    AT=TRANSPOSE(A)
  END SUBROUTINE XPoseSqMat

#ifdef PARALLEL
  !===============================================================================
  !     Max block of a DBCSR matrix
  !===============================================================================
  FUNCTION Max_DBCSR(A,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: Max_DBCSR
    TYPE(DBCSR),        INTENT(INOUT) :: A
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    REAL(DOUBLE)                      :: Mx
    INTEGER                           :: IL,IG,MA,NA,J,JP,P,IStrtA,IStopA
    REAL(DOUBLE), EXTERNAL            :: DBL_Dot
    !         TYPE(BCSR) :: B
    !-------------------------------------------------------------------------------
    !         CALL SetEq(B,A)
    !         Max_DBCSR=Max_BCSR(B)
    !         CALL delete(B)
    !         RETURN

    Mx=Zero
    DO IL=1,A%NAtms
      IG=IL+OffSt%I(MyId)
      MA=BSiz%I(IG)
      IStrtA=A%RowPt%I(IL)
      IStopA=A%RowPt%I(IL+1)-1
      DO JP=IStrtA,IStopA
        J=A%ColPt%I(JP)
        P=A%BlkPt%I(JP)
        NA=BSiz%I(J)
        Mx=MAX(Mx,SQRT(DBL_Dot(MA*NA,A%MTrix%D(P),A%MTrix%D(P))))
        !               Mx=MAX(Mx,DSQRT(DDot(MA*NA,A%MTrix%D(P),1,A%MTrix%D(P),1)))
        !               Mx=MAX(Mx,FNorm(MA*NA,A%MTrix%D(P:)))
      ENDDO
    ENDDO
    Max_DBCSR=AllReduce(Mx,Op_O=MPI_MAX)
    PerfMon%FLOP=PerfMon%FLOP+DBLE(2*A%NNon0)
    IF(PRESENT(Perf_O)) &
         Perf_O%FLOP=Perf_O%FLOP+DBLE(2*A%NNon0)
  END FUNCTION Max_DBCSR
#endif
  !===============================================================================
  !     Max block of a BCSR matrix
  !===============================================================================
  FUNCTION Max_BCSR(A,Perf_O)
    IMPLICIT NONE
    REAL(DOUBLE)                      :: Max_BCSR
    TYPE(BCSR),         INTENT(IN)    :: A
    TYPE(TIME),OPTIONAL,INTENT(INOUT) :: Perf_O
    INTEGER                           :: IL,MA,NA,MNA,J,JP,P,IStrtA,IStopA,iSMat
    REAL(DOUBLE), EXTERNAL            :: DBL_Dot
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      IF(A%NSMat.NE.1.AND.A%NSMat.NE.2.AND.A%NSMat.NE.4) &
           CALL Halt('Max_BCSR: NSMat doesn''t have an expected value!')

      Max_BCSR=Zero
      DO IL=1,A%NAtms
        MA=BSiz%I(IL)
        IStrtA=A%RowPt%I(IL)
        IStopA=A%RowPt%I(IL+1)-1
        DO JP=IStrtA,IStopA
          J=A%ColPt%I(JP)
          P=A%BlkPt%I(JP)
          NA=BSiz%I(J)
          MNA=MA*NA
          DO iSMat=1,A%NSMat
            Max_BCSR=MAX(Max_BCSR,SQRT(DBL_Dot(MNA,A%MTrix%D(P),A%MTrix%D(P))))
            P=P+MNA
          ENDDO
        ENDDO
      ENDDO
      PerfMon%FLOP=PerfMon%FLOP+DBLE(2*A%NNon0*A%NSMat)
      IF(PRESENT(Perf_O)) &
           Perf_O%FLOP=Perf_O%FLOP+DBLE(2*A%NNon0*A%NSMat)
#ifdef PARALLEL
    ENDIF
    IF(InParallel) &
         CALL BCast(Max_BCSR)
#endif
  END FUNCTION Max_BCSR
  !===================================================================

  !===================================================================

  SUBROUTINE RowColChk(N,RowPt,ColPt,Mssg)
    INTEGER :: I,N
    TYPE(INT_VECT) :: RowPt,ColPt
    CHARACTER(LEN=*) :: MSSG
    DO I=1,RowPt%I(N+1)-1
      IF(ColPt%I(I)<1.OR.ColPt%I(I)>NAtoms)THEN
        WRITE(*,*)' I = ',I,' ColPt = ',ColPt%I(I)
        CALL Halt(Mssg)
      ENDIF
    ENDDO
  END SUBROUTINE RowColChk

  !===================================================================

  !===================================================================
  SUBROUTINE MatrixCheck(A,Mssg)
    TYPE(BCSR),       INTENT(IN) :: A
    CHARACTER(LEN=*), INTENT(IN) :: Mssg
    INTEGER                      :: I,J,SzRow,SzCol,SzMat, &
         IStrtA,IStopA
    !-------------------------------------------------------------------------------
    SzRow=SIZE(A%RowPt%I)
    SzCol=MIN(SIZE(A%ColPt%I),SIZE(A%BlkPt%I))
    SzMat=SIZE(A%MTrix%D)
    DO I=1,NAtoms+1
      IF(A%RowPt%I(I)>SzCol) &
           CALL Halt(TRIM(Mssg)//': RowPt error in MatrixCheck, I = '//TRIM(IntToChar(I)))
    ENDDO
    DO I=1,NAtoms
      IStrtA=A%RowPt%I(I)
      IStopA=A%RowPt%I(I+1)-1
      DO J=IStrtA,IStopA
        IF(A%ColPt%I(J)>NAtoms) &
             CALL Halt(TRIM(Mssg)//': ColPt error 1 in MatrixCheck, J = '//TRIM(IntToChar(J)))
        IF(A%ColPt%I(J)<=0) &
             CALL Halt(TRIM(Mssg)//': ColPt error 2 in MatrixCheck, J = '//TRIM(IntToChar(J)))
        IF(A%BlkPt%I(J)>SzMat) &
             CALL Halt(TRIM(Mssg)//': BlkPt error 1 in MatrixCheck, J = '//TRIM(IntToChar(J)))
        IF(A%BlkPt%I(J)<=0) &
             CALL Halt(TRIM(Mssg)//': BlkPt error 2 in MatrixCheck, J = '//TRIM(IntToChar(J)))
      ENDDO
    ENDDO
  END SUBROUTINE MatrixCheck
#ifdef PARALLEL
  !===============================================================================
  !     More intelegent domain decomposition based on a matrices symbolic structure
  !==================================================================================
  SUBROUTINE RePart(MatrixFile)
    CHARACTER(LEN=*), INTENT(IN)        :: MatrixFile
    TYPE(INT_VECT)                      :: RowPt,ColPt,Count,NAv
    TYPE(CHR_VECT)                      :: Chr,CBeg,CEnd
    INTEGER                             :: NNon0Av,I,J,K,IPrc, &
         NAtms,NNon0,NBlks,IOS, &
         MA,IStrt,IStop,JP,KKK
    INTEGER, PARAMETER                  :: Mns=1,Pls=2
    INTEGER, DIMENSION(Mns:Pls)         :: Beg3,End2
    INTEGER, DIMENSION(Mns:Pls,Mns:Pls) :: End3,Dv
    LOGICAL                             :: Exists
    !-------------------------------------------------------------------------------
    !        Allocate domain limits if nessesary

    IF(.NOT.AllocQ(Beg%Alloc))  CALL New(Beg,NPrc-1,0)
    IF(.NOT.AllocQ(End%Alloc))  CALL New(End,NPrc-1,0)
    IF(.NOT.AllocQ(OffSt%Alloc))CALL New(OffSt,NPrc-1,0)
    IF(MyId==ROOT)THEN
      !-------------------------------------------------------------------------------
      !           Get global structure from BCSR matrix

      INQUIRE(FILE=MatrixFile,EXIST=Exists)
      IF(Exists)THEN
        OPEN(UNIT=Seq,FILE=MatrixFile,STATUS='OLD', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ELSE
        CALL Halt(' Matrix File '//TRIM(MatrixFile)//' not found')
      ENDIF
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)NAtms,NNon0,NBlks
      CALL New(RowPt,NAtoms+1)
      CALL New(ColPt,NBlks)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(RowPt%I(I),I=1,NAtoms+1)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(ColPt%I(I),I=1,NBlks)
      CLOSE(Seq)
      !-------------------------------------------------------------------------------
      !           Calculate the number of NNon0s per row

      CALL New(Count,NAtoms)
      CALL New(NAv,NPrc-1,0)

      DO I=1,NAtoms
        MA=BSiz%I(I)
        Count%I(I)=0
        IStrt=RowPt%I(I)
        IStop=RowPt%I(I+1)-1
        IF(IStrt/=0.AND.IStop/=0)THEN
          DO JP=IStrt,IStop
            J=ColPt%I(JP)
            Count%I(I)=Count%I(I)+MA*BSiz%I(J)
          ENDDO
        ENDIF
        WRITE(55,*)I,Count%I(I)
      ENDDO
      !-------------------------------------------------------------------------------
      !           Look ahead 1 node algorithm for decomposition
      !           If NPrc==2^p, should use bisection instead.
      !           -- See ORB.F90 for work in progress...

      Beg%I(0)=1
      End%I(NPrc-1)=NAtoms

      !            NNon0Av=SUM(Count%I(1:NAtoms))/DBLE(NPrc)
      !            WRITE(*,*)' NNon0Av = ',NNon0Av

      DO IPrc=0,NPrc-2
        !-------------------------------------------------------------------------------
        !             Running averages to section

        NNon0Av=SUM(Count%I(Beg%I(IPrc):NAtoms))/DBLE(NPrc-IPrc)
        NAv%I(IPrc)=NNon0Av
        !-------------------------------------------------------------------------------
        !              Forcast Beg and End for (de/inc)rements of (+/- 1)
        !
        DO K=Beg%I(IPrc),NAtoms
          IF(SUM(Count%I(Beg%I(IPrc):K))>=NNon0Av)THEN
            End2(Pls)=K
            EXIT
          ENDIF
        ENDDO
        End2(Mns)=End2(Pls)-1
        Beg3(Pls)=End2(Pls)+1
        Beg3(Mns)=End2(Mns)+1
        End3(Pls,Pls)=NAtoms
        DO K=Beg3(Pls),NAtoms
          IF(SUM(Count%I(Beg3(Pls):K))>=NNon0Av)THEN
            End3(Pls,Pls)=K
            EXIT
          ENDIF
        ENDDO
        End3(Pls,Mns)=End3(Pls,Pls)-1
        DO K=Beg3(Mns),NAtoms
          IF(SUM(Count%I(Beg3(Mns):K))>=NNon0Av)THEN
            End3(Mns,Pls)=K
            EXIT
          ENDIF
        ENDDO
        End3(Mns,Mns)=End3(Mns,Pls)-1
        !-------------------------------------------------------------------------------
        !              These are deviations from the running average for each choice
        !              of a place to section and the possilbe sectioning in the next
        !              iteration
        !
        Dv(Pls,Pls)= &
             ABS(NNon0Av-SUM(Count%I(Beg%I(IPrc):End2(Pls)))) &
             +ABS(NNon0Av-SUM(Count%I(Beg3(Pls):End3(Pls,Pls))))
        Dv(Pls,Mns)= &
             ABS(NNon0Av-SUM(Count%I(Beg%I(IPrc):End2(Pls))))   &
             +ABS(NNon0Av-SUM(Count%I(Beg3(Pls):End3(Pls,Mns))))
        Dv(Mns,Pls)= &
             ABS(NNon0Av-SUM(Count%I(Beg%I(IPrc):End2(Mns))))   &
             +ABS(NNon0Av-SUM(Count%I(Beg3(Mns):End3(Mns,Pls))))
        Dv(Mns,Mns)= &
             ABS(NNon0Av-SUM(Count%I(Beg%I(IPrc):End2(Mns))))   &
             +ABS(NNon0Av-SUM(Count%I(Beg3(Mns):End3(Mns,Mns))))

        !write(*,*)' 0th Deviation +:',TRIM(IntToChar(ABS(NNon0Av-SUM(Count%I(Beg%I(IPrc):End2(Pls)))))),   &
        !                         '-:',TRIM(IntToChar(ABS(NNon0Av-SUM(Count%I(Beg%I(IPrc):End2(Mns))))))
        !write(*,*)' 1st Deviation ++:',TRIM(IntToChar(ABS(NNon0Av-SUM(Count%I(Beg3(Pls):End3(Pls,Pls)))))) ,&
        !                         '+-:',TRIM(IntToChar(ABS(NNon0Av-SUM(Count%I(Beg3(Pls):End3(Pls,Mns)))))) ,&
        !                         '-+:',TRIM(IntToChar(ABS(NNon0Av-SUM(Count%I(Beg3(Mns):End3(Mns,Pls)))))) ,&
        !                         '--:',TRIM(IntToChar(ABS(NNon0Av-SUM(Count%I(Beg3(Mns):End3(Mns,Mns))))))
        !-------------------------------------------------------------------------------
        !              Pick the best section based on minimizing the I and I+1 deviation
        !              from the average
        !
        IF(MIN(Dv(Pls,Pls),Dv(Pls,Mns))< &
             MIN(Dv(Mns,Pls),Dv(Mns,Mns)))THEN
          End%I(IPrc)=End2(Pls)
          Beg%I(IPrc+1)=End%I(IPrc)+1

          !                  WRITE(*,12)Dv(Pls,Pls),Dv(Pls,Mns),Dv(Mns,Pls),Dv(Mns,Mns), &
          !                             IPrc,Beg%I(IPrc),End%I(IPrc),                    &
          !                             SUM(Count%I(Beg%I(IPrc):End%I(IPrc)))

        ELSE
          End%I(IPrc)=MAX(Beg%I(IPrc),End2(Mns))
          Beg%I(IPrc+1)=End%I(IPrc)+1
          !                  WRITE(*,13)Dv(Pls,Pls),Dv(Pls,Mns),Dv(Mns,Pls),Dv(Mns,Mns), &
          !                             IPrc,Beg%I(IPrc),End%I(IPrc) ,                   &
          !                             SUM(Count%I(Beg%I(IPrc):End%I(IPrc)))
        ENDIF

        !12 format('I  ++ = ',I5,', +- = ',I5,', -+ = ',I5,', -- = ',I5, &
        !          ' IPrc = ',I5,' B = ',I5,' E = ',I5,' Sum = ',I5)
        !13 format('II ++ = ',I5,', +- = ',I5,', -+ = ',I8,', -- = ',I5, &
        !          ' IPrc = ',I5,' B = ',I5,' E = ',I5,'Sum = ',I5)
        !         DO I=0,NPrc-1
        !         WRITE(*,*)I,Beg%I(I),End%I(I),SUM(Count%I(Beg%I(I):End%I(I))),NAv%I(I), &
        !                   ' Diff = ',SUM(Count%I(Beg%I(I):End%I(I)))-NAv%I(I)
        !         ENDDO


      ENDDO
      !---------------------------------------------------------
      !           Calculate DBCSR matrix RowPt -> GRwPt off-sets

      DO I=0,NPrc-1
        OffSt%I(I)=End%I(I)-Beg%I(I)+1
      ENDDO
      DO I=1,NPrc-1
        OffSt%I(I)=OffSt%I(I)+OffSt%I(I-1)
      ENDDO
      DO I=NPrc-1,1,-1
        OffSt%I(I)=OffSt%I(I-1)
      ENDDO
      OffSt%I(0)=0
      !-------------------------------------------------------------------------------


      MaxAtmsNode=0
      DO K=0,NPrc-1
        MaxAtmsNode=MAX(MaxAtmsNode,End%I(K)-Beg%I(K)+1)
      ENDDO
      MaxAtmsNode=1+MaxAtmsNode
      !-------------------------------------------------------------------------------
      !           Put the domain boundaries to disk

      !            CALL OpenHDF(Ctrl%Info(Ctrl%ISet))
      !            CALL Put(Beg,'beg')
      !            CALL Put(End,'end')
      !            CALL Put(OffSt,'dbcsroffsets')
      !            CALL CloseHDF()
      !-------------------------------------------------------------------------------
      !           Debug if asked

      IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)'New decomposition based on: ',TRIM(MatrixFile)
        CALL New(CBeg,NPrc-1,0)
        CALL New(CEnd,NPrc-1,0)
        DO I=0,NPrc-1
          CBeg%C(I)=IntToChar(Beg%I(I))
          CEnd%C(I)=IntToChar(End%I(I))
        ENDDO
        WRITE(Out,*)'Atomic re-partitioning = ',                              &
             ('['//TRIM(CBeg%C(K))//'-'//TRIM(CEnd%C(K))//'], ',K=0,NPrc-2),  &
             '['//TRIM(CBeg%C(NPrc-1))//'-'//TRIM(CEnd%C(NPrc-1)),']'
        DO K=0,NPrc-1
          WRITE(*,*)' K = ',K,'Beg = ', Beg%I(K),' End = ',End%I(K)
          CBeg%C(K)=IntToChar(SUM(Count%I(Beg%I(K):End%I(K))))
        ENDDO
        WRITE(Out,*)'Number of Non0 elements per node = ',                &
             (TRIM(CBeg%C(K))//', ',K=0,NPrc-2),TRIM(CBeg%C(NPrc-1))
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
        CALL Delete(CBeg)
        CALL Delete(CEnd)
      ENDIF
      !-------------------------------------------------------------------------------
      !           Tidy up

      CALL Delete(RowPt)
      CALL Delete(ColPt)
      CALL Delete(Count)
      CALL Delete(NAv)
    ENDIF
    CALL Bcast(MaxAtmsNode)
    CALL Bcast(Beg)
    CALL Bcast(End)
    CALL Bcast(OffSt)
    RETURN
1   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in RePart.')
  END SUBROUTINE RePart
  !===============================================================================
  !     Simulated Annealing Domain Decomposition
  !===============================================================================
  SUBROUTINE SADD(MatrixFile)
    CHARACTER(LEN=*), INTENT(IN)        :: MatrixFile
    TYPE(INT_VECT)                      :: RowPt,ColPt,Counts,NAv
    TYPE(CHR_VECT)                      :: Chr,CBeg,CEnd
    INTEGER                             :: NNon0Av,I,J,K,M,IPrc,    &
         NAtms,NNon0,NBlks,IOS,   &
         MA,IStrt,IStop,JP,       &
         F1,F2,B,BOld,NOvr,NSuc,NLim
    LOGICAL                             :: Exists
    REAL(DOUBLE)                        :: AtsPerPrc,Ats,Temp,TFct,R,DF

    !-------------------------------------------------------------------------------
    !        Allocate domain limits if nessesary

    IF(.NOT.AllocQ(Beg%Alloc))  CALL New(Beg,NPrc-1,0)
    IF(.NOT.AllocQ(End%Alloc))  CALL New(End,NPrc-1,0)
    IF(.NOT.AllocQ(OffSt%Alloc))CALL New(OffSt,NPrc-1,0)
    IF(MyId==ROOT)THEN
      !-------------------------------------------------------------------------------
      !           Get global structure from BCSR matrix

      INQUIRE(FILE=MatrixFile,EXIST=Exists)
      IF(Exists)THEN
        OPEN(UNIT=Seq,FILE=MatrixFile,STATUS='OLD', &
             FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      ELSE
        CALL Halt(' Matrix File '//TRIM(MatrixFile)//' not found')
      ENDIF
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)NAtms,NNon0,NBlks
      CALL New(RowPt,NAtoms+1)
      CALL New(ColPt,NBlks)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(RowPt%I(I),I=1,NAtoms+1)
      READ(UNIT=Seq,Err=1,IOSTAT=IOS)(ColPt%I(I),I=1,NBlks)
      CLOSE(Seq)
      !-------------------------------------------------------------------------------
      !           Calculate the number of NNon0s per row

      CALL New(Counts,NAtoms)
      CALL New(NAv,NPrc-1,0)

      DO I=1,NAtoms
        MA=BSiz%I(I)
        Counts%I(I)=0
        IStrt=RowPt%I(I)
        IStop=RowPt%I(I+1)-1
        IF(IStrt/=0.AND.IStop/=0)THEN
          DO JP=IStrt,IStop
            J=ColPt%I(JP)
            Counts%I(I)=Counts%I(I)+MA*BSiz%I(J)
          ENDDO
        ENDIF
        WRITE(55,*)I,Counts%I(I)
      ENDDO
      !-------------------------------------------------------------------------------
      AtsPerPrc=DBLE(NAtoms)/DBLE(NPrc)
      Beg%I(0)=1
      Ats=0
      DO I=1,NPrc-1
        Ats=Ats+AtsPerPrc
        Beg%I(I)=Ats
      ENDDO
      F1=ObjectiveF(Beg%I,Counts%I)
      TFct=0.90D0
      Temp=F1*8.0D-2
      NLim=10*NPrc
      NOvr=1000*NPrc
      DO J=1,100
        NSuc=0
        DO K=1,NOvr
          M=RANDOM_INT((/1,NPrc-1/))
          IF(M==NPrc-1)THEN
            B=RANDOM_INT((/Beg%I(M-1)+2,NAtoms-2/))
          ELSEIF(M==0)THEN
            B=RANDOM_INT((/0,Beg%I(M+1)-2/))
          ELSE
            B=RANDOM_INT((/Beg%I(M-1)+2,Beg%I(M+1)-2/))
          ENDIF
          BOld=Beg%I(M)
          Beg%I(M)=B
          F2=ObjectiveF(Beg%I,Counts%I)
          DF=DBLE(F2-F1)
          R=RANDOM_DBL((/Zero,One/))
          IF(DF<Zero.OR.R<EXP(-DF/Temp))THEN
            F1=F2
            NSuc=NSuc+1
          ELSE
            Beg%I(M)=BOld
          ENDIF
          IF(NSuc>NLim)EXIT
        ENDDO
        !               WRITE(*,11)Temp,NSuc,F1
        !            11 FORMAT('Temp = ',F14.6,' NSuc = ',I3,' F = ',I14)
        Temp=Temp*TFct
      ENDDO

      DO I=0,NPrc-2
        End%I(I)=Beg%I(I+1)-1
      ENDDO
      End%I(NPrc-1)=NAtoms
      !---------------------------------------------------------
      !           Calculate DBCSR matrix RowPt -> GRwPt off-sets

      DO I=0,NPrc-1
        OffSt%I(I)=End%I(I)-Beg%I(I)+1
      ENDDO
      DO I=1,NPrc-1
        OffSt%I(I)=OffSt%I(I)+OffSt%I(I-1)
      ENDDO
      DO I=NPrc-1,1,-1
        OffSt%I(I)=OffSt%I(I-1)
      ENDDO
      OffSt%I(0)=0
      !-------------------------------------------------------------------------------


      MaxAtmsNode=0
      DO K=0,NPrc-1
        MaxAtmsNode=MAX(MaxAtmsNode,End%I(K)-Beg%I(K)+1)
      ENDDO
      MaxAtmsNode=1+MaxAtmsNode
      !-------------------------------------------------------------------------------
      !           Debug if asked

      IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)'New decomposition based on: ',TRIM(MatrixFile)
        CALL New(CBeg,NPrc-1,0)
        CALL New(CEnd,NPrc-1,0)
        DO I=0,NPrc-1
          CBeg%C(I)=IntToChar(Beg%I(I))
          CEnd%C(I)=IntToChar(End%I(I))
        ENDDO
        WRITE(Out,*)'Atomic re-partitioning = ',                              &
             ('['//TRIM(CBeg%C(K))//'-'//TRIM(CEnd%C(K))//'], ',K=0,NPrc-2),  &
             '['//TRIM(CBeg%C(NPrc-1))//'-'//TRIM(CEnd%C(NPrc-1)),']'
        DO K=0,NPrc-1
          !                  WRITE(*,*)' K = ',K,'Beg = ', Beg%I(K),' End = ',End%I(K)
          CBeg%C(K)=IntToChar(SUM(Counts%I(Beg%I(K):End%I(K))))
        ENDDO
        WRITE(Out,*)'Number of Non0 elements per node = ',                &
             (TRIM(CBeg%C(K))//', ',K=0,NPrc-2),TRIM(CBeg%C(NPrc-1))
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
        CALL Delete(CBeg)
        CALL Delete(CEnd)
      ENDIF
      !-------------------------------------------------------------------------------
      !           Tidy up

      CALL Delete(RowPt)
      CALL Delete(ColPt)
      CALL Delete(Counts)
      CALL Delete(NAv)
    ENDIF
    CALL Bcast(MaxAtmsNode)
    CALL Bcast(Beg)
    CALL Bcast(End)
    CALL Bcast(OffSt)
    RETURN
1   CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in GreedySCF.')
  END SUBROUTINE SADD
  !===============================================================================
  !     Cumulative deviation of NNon0s per proc from the average
  !===============================================================================
  INTEGER FUNCTION ObjectiveF(Beg,Counts)
    INTEGER,INTENT(IN),DIMENSION(:) :: Counts,Beg
    INTEGER                         :: I,NP,NA,B,E,C,AverageNon0
    NA=SIZE(Counts)
    NP=SIZE(Beg)
    AverageNon0=INT(DBLE(SUM(Counts(1:NA)))/DBLE(NP))
    ObjectiveF=0
    DO I=1,NP
      B=Beg(I)
      IF(I==NP)THEN
        E=NA
      ELSE
        E=Beg(I+1)-1
      ENDIF
      C=SUM(Counts(B:E))
      ObjectiveF=MAX(ObjectiveF,MAX(0,C-AverageNon0))
    ENDDO
  END FUNCTION ObjectiveF
  !===============================================================================

  !===============================================================================
  SUBROUTINE Plot_DBCSR(A,Name)
    TYPE(DBCSR) :: A
    CHARACTER(LEN=*) :: Name
    INTEGER :: I,K
    IF(PrintFlags%Mat/=DEBUG_PLOT_MATRICES)RETURN
    CALL LocalToGlobal(A)
    IF(MyID==ROOT)THEN
      CALL OpenASCII('PlotFile_1',Plt,NewFile_O=.TRUE.)
      DO I=1,NAtoms
        DO K=A%GRwPt%I(I),A%GRwPt%I(I+1)-1
          WRITE(Plt,1)I,NAtoms-A%GClPt%I(K)
        ENDDO
      ENDDO
      CLOSE(Plt)
      CALL OpenASCII('PlotFile_2',Plt,NewFile_O=.TRUE.)
      WRITE(Plt,2)
      WRITE(Plt,3)TRIM(Name)//'.eps'
      WRITE(Plt,4)NAtoms+1
      WRITE(Plt,5)NAtoms+1
      WRITE(Plt,6); WRITE(Plt,7); WRITE(Plt,8)
      !          IF(NPrc>1)THEN
      !             WRITE(Plt,9)NAtoms+1,End%I(0),BakSlash
      !             DO I=1,NPrc-2
      !                WRITE(Plt,10)End%I(I),BakSlash
      !             ENDDO
      !             WRITE(Plt,11)
      !          ELSE
      WRITE(Plt,12)
      !          ENDIF
      CLOSE(Plt)
      !          CALL SYSTEM(' gnuplot < PlotFile_2 ')
      !          CALL SYSTEM(' rm -f PlotFile_1 PlotFile_2 ')
    ENDIF
2   FORMAT('set term  postscript eps  "Times-Roman" 18')
3   FORMAT('set output "',A,'"')
4   FORMAT('set xrange [0 :',I12,']')
5   FORMAT('set yrange [0 :',I12,']')
6   FORMAT('set size 0.8, 1.0 ')
7   FORMAT('set noxtics ')
8   FORMAT('set noytics ')
11  FORMAT("     'PlotFile_1' using 1:2 with points 4 ")
12  FORMAT("plot 'PlotFile_1' using 1:2 with points 4 ")
13  FORMAT("plot 'PlotFile_1' using 1:2 with points 4 ")
  END SUBROUTINE Plot_DBCSR
#endif

#ifdef USE_METIS
  !===============================================================================
  !     Performs P.A, where P is a permutation opperator
  !===============================================================================
  SUBROUTINE MetisReorder(A)
    TYPE(BCSR)     :: A
    TYPE(INT_VECT) :: Perm,IPerm
    TYPE(INT_VECT) :: RowPt,ColPt,NSiz
    INTEGER        :: N,I,J,K,I1,I2,J1,J2,K1,K2,R,S
    !                                                0 1 2 3 4 5 6 7
    INTEGER,DIMENSION(8),PARAMETER :: Opt=(/1,3,1,1,0,1,60,2/)
    !----------------------------------------------------------------------------
    !        Allocate CSR arrays for graph
    CALL New(RowPt,A%NAtms+1)
    CALL New(ColPt,A%NBlks-A%NAtms)
    !        Allocate permutation vectors
    CALL New(Perm,A%NAtms)
    CALL New(IPerm,A%NAtms)
    !----------------------------------------------------------------------------
    !        Compute graph
    N=0
    DO I=1, A%NAtms
      RowPt%I(I)=N+1
      DO K=A%RowPt%I(I),A%RowPt%I(I+1)-1
        J=A%ColPt%I(K)
        IF(I/=J)THEN
          N=N+1
          ColPt%I(N)=J
        ENDIF
      ENDDO
    ENDDO
    RowPt%I(A%NAtms+1)=N+1
    !----------------------------------------------------------------------------
    !        Compute Metis ordering for this graph
    CALL METIS_NodeND(A%NAtms,RowPt%I,ColPt%I,1,Opt,Perm%I,IPerm%I)
    !        tidy
    CALL Delete(RowPt)
    CALL Delete(ColPt)
    !----------------------------------------------------------------------------
    !        Perform the symmetric permutation on A defined by Perm
    CALL SymPe(A,Perm)
    !        Recompute the new BCSR blocking
    CALL New(NSiz,A%NAtms)
    NSiz%I=BSiz%I
    DO I=1,A%NAtms
      BSiz%I(Perm%I(I))=NSiz%I(I)
    ENDDO
    OffS%I(1)=1
    DO I=2,NAtoms
      OffS%I(I)=OffS%I(I-1)+BSiz%I(I-1)
    ENDDO
    CALL Delete(NSiz)
    CALL Delete(Perm)
    CALL Delete(IPerm)
  END SUBROUTINE MetisReorder
#endif

  !===============================================================================
  !     Performs the symetric permutation P^T.A.P, where P is in vector form
  !===============================================================================
  SUBROUTINE SymPe(A,P)
    TYPE(BCSR)     :: A
    TYPE(INT_VECT) :: P
    TYPE(INT_VECT) :: RowPt,ColPt,BlkPt
    INTEGER        :: I,J,K
    !----------------------------------------------------------------------------
    !        Allocations
    CALL New(RowPt,A%NAtms+1)
    CALL New(ColPt,A%NBlks)
    CALL New(BlkPt,A%NBlks)
    !----------------------------------------------------------------------------
    !        Permute rows
    DO I=1,A%NAtms
      RowPt%I(P%I(I)+1)=A%RowPt%I(I+1)-A%RowPt%I(I)
    ENDDO
    RowPt%I(1)=1
    DO I=1,A%NAtms
      RowPt%I(I+1)=RowPt%I(I+1)+RowPt%I(I)
    ENDDO
    !        Copy col and blk pointers
    J=1
    DO I=1,A%NAtms
      J=RowPt%I(P%I(I))
      DO K=A%RowPt%I(I),A%RowPt%I(I+1)-1
        ColPt%I(J)=A%ColPt%I(K)
        BlkPt%I(J)=A%BlkPt%I(K)
        J=J+1
      ENDDO
    ENDDO
    !----------------------------------------------------------------------------
    !        Permute cols
    DO I=1,A%NBlks
      A%ColPt%I(I)=P%I(ColPt%I(I))
    ENDDO
    !----------------------------------------------------------------------------
    !        Remaining back copy
    A%RowPt%I(1:A%NAtms+1)=RowPt%I(1:A%NAtms+1)
    A%BlkPt%I(1:A%NBlks)=BlkPt%I(1:A%NBlks)
    !----------------------------------------------------------------------------
    !        Clean up
    CALL Delete(RowPt)
    CALL Delete(ColPt)
    CALL Delete(BlkPt)

  END SUBROUTINE SymPe

  SUBROUTINE PlotDecay(A,GM,Name)
    TYPE(BCSR)     :: A
    TYPE(CRDS)     :: GM
    CHARACTER(LEN=*) :: Name
    TYPE(DBL_VECT) :: Distance,Magnitude
    TYPE(INT_VECT) :: ISort
    INTEGER :: I,J,K,L,M,N,P
    REAL(DOUBLE),EXTERNAL   :: DBL_Dot
    CALL New(ISort,A%NBlks)
    CALL New(Distance,A%NBlks)
    CALL New(Magnitude,A%NBlks)
    L=0
    DO I=1,NAtoms
      M=BSiz%I(I)
      DO K=A%RowPt%I(I),A%RowPt%I(I+1)-1
        J=A%ColPt%I(K)
        P=A%BlkPt%I(K)
        N=BSiz%I(J)
        L=L+1
        ISort%I(L)=L
        Distance%D(L)=SQRT((GM%Carts%D(1,I)-GM%Carts%D(1,J))**2 &
             +(GM%Carts%D(2,I)-GM%Carts%D(2,J))**2 &
             +(GM%Carts%D(3,I)-GM%Carts%D(3,J))**2 )
        Magnitude%D(L)=SQRT(DBL_Dot(M*N,A%MTrix%D(P),A%MTrix%D(P)))
      ENDDO
    ENDDO
    CALL Sort(Distance,ISort,A%NBlks)
    !
    CALL OpenASCII(TRIM(Name)//'_MatrixDecayData',Plt,NewFile_O=.TRUE.)
    DO I=1,A%NBlks
      WRITE(Plt,*)Distance%D(I),Magnitude%D(ISort%I(I))
    ENDDO
    CLOSE(Plt)

    CALL OpenASCII(TRIM(Name)//'_MatrixDecayPlotMe',Plt,NewFile_O=.TRUE.)
    WRITE(Plt,2)
    WRITE(Plt,3)TRIM(Name)//'.eps'
    WRITE(Plt,6)
    WRITE(Plt,*)'set pointsize 0.1'
    WRITE(Plt,*)'set logscale y'
    WRITE(Plt,*)"plot '"//TRIM(Name)//"_MatrixDecayData' using 1:2 notitle with points 1 "
    CLOSE(Plt)

    CALL Delete(ISort)
    CALL Delete(Distance)
    CALL Delete(Magnitude)

2   FORMAT('set term  postscript eps  "Times-Roman" 18')
    !   2   FORMAT('set term  jpeg transparent')
3   FORMAT('set output "',A,'"')
6   FORMAT('set size square ')
  END SUBROUTINE PlotDecay

  SUBROUTINE Plot_BCSR(A,Name)
    TYPE(BCSR) :: A
    CHARACTER(LEN=*) :: Name
    INTEGER :: I,K
    IF(PrintFlags%Mat/=DEBUG_PLOT_MATRICES)RETURN
    CALL OpenASCII(TRIM(Name)//'_PlotFile_1',Plt,NewFile_O=.TRUE.)
    DO I=1,NAtoms
      DO K=A%RowPt%I(I),A%RowPt%I(I+1)-1
        WRITE(Plt,1)I,NAtoms-A%ColPt%I(K)
      ENDDO
    ENDDO
    CLOSE(Plt)

    CALL OpenASCII(TRIM(Name)//'_GnuPlotMe',Plt,NewFile_O=.TRUE.)
    WRITE(Plt,2)
    WRITE(Plt,3)TRIM(Name)//'.eps'
    WRITE(Plt,4)DBLE(NAtoms+1)
    WRITE(Plt,5)DBLE(NAtoms)
    WRITE(Plt,6); WRITE(Plt,7); WRITE(Plt,8)
    WRITE(Plt,*)'set pointsize '//TRIM(FltToChar(50.D0/DBLE(NAtoms)))
    WRITE(Plt,*)"plot '"//TRIM(Name)//"_PlotFile_1' using 1:2 notitle with points 1 "
    CLOSE(Plt)
1   FORMAT(2(1x,I16))
2   FORMAT('set term  postscript eps  "Times-Roman" 18')
    !   2   FORMAT('set term  jpeg transparent')
3   FORMAT('set output "',A,'"')
4   FORMAT('set xrange [0:',F12.3,']')
5   FORMAT('set yrange [-1:',F12.3,']')
6   FORMAT('set size square ')
7   FORMAT('set noxtics ')
8   FORMAT('set noytics ')
  END SUBROUTINE Plot_BCSR

  SUBROUTINE MStats(A,FileName)
    TYPE(BCSR),           INTENT(INOUT) :: A
    CHARACTER(LEN=*)                    :: FileName
    TYPE(CRDS)                          :: GM
    REAL(DOUBLE),EXTERNAL               :: DBL_Dot
    REAL(DOUBLE)                        :: R,E,Op
    INTEGER                             :: I,J,JP,K,P,Q,MA,NA,MN,MN1, &
         IStrtA,IStopA

    TYPE(INT_VECT)            :: Stat

    !-------------------------------------------------------------------------------


    CALL New(Stat,3)
    CALL Get(Stat,'current')

    CALL Get(GM,Tag_O=CurGeom)

    CALL OpenASCII(FileName,Tmp,NewFile_O=.TRUE.)
    K=1
    Q=1
    Op=Zero
    DO I=1,NAtoms
      MA=BSiz%I(I)
      IStrtA=A%RowPt%I(I)
      IStopA=A%RowPt%I(I+1)-1
      IF(IStrtA/=0.AND.IStopA/=0)THEN
        DO JP=IStrtA,IStopA
          J=A%ColPt%I(JP)
          P=A%BlkPt%I(JP)
          NA=BSiz%I(J)
          MN=MA*NA
          MN1=MN-1
          E=SQRT(DBL_Dot(MN,A%MTrix%D(P),A%MTrix%D(P)))
          R=SQRT( (GM%Carts%D(1,I)-GM%Carts%D(1,J))**2 &
               +(GM%Carts%D(2,I)-GM%Carts%D(2,J))**2 &
               +(GM%Carts%D(3,I)-GM%Carts%D(3,J))**2 )
          WRITE(Tmp,*)R,E
        ENDDO
      ENDIF
    ENDDO
    CALL Delete(GM)
    CLOSE(Tmp)
  END SUBROUTINE MStats
  !-------------------------------------------------------------------------------
  !      Frobenius Norm
  !-------------------------------------------------------------------------------
#ifdef PARALLEL
  FUNCTION FNorm_DBCSR(M)
    TYPE(DBCSR)       :: M
    INTEGER          :: I
    REAL(DOUBLE)     :: FNorm_DBCSR,Frob
    FNorm_DBCSR = Zero
    DO I = 1, M%NNon0
      FNorm_DBCSR = FNorm_DBCSR + (M%MTrix%D(I))**2
    ENDDO
    Frob=AllReduce(FNorm_DBCSR)
    FNorm_DBCSR=SQRT(Frob)
  END FUNCTION FNorm_DBCSR
#endif
  FUNCTION FNorm_BCSR(M)
    TYPE(BCSR)       :: M
    INTEGER          :: I
    REAL(DOUBLE)     :: FNorm_BCSR
    FNorm_BCSR = Zero
    DO I = 1, M%NNon0
      FNorm_BCSR = FNorm_BCSR + (M%MTrix%D(I))**2
    ENDDO
    FNorm_BCSR = SQRT(FNorm_BCSR)
  END FUNCTION FNorm_BCSR
  !
  !------------------------------------------------------------------
  !


      FUNCTION CROSS_PRODUCT(V1,V2)
        REAL(DOUBLE),DIMENSION(3) :: V1,V2,CROSS_PRODUCT
        !
        CROSS_PRODUCT(1)=V1(2)*V2(3)-V1(3)*V2(2)
        CROSS_PRODUCT(2)=V1(3)*V2(1)-V1(1)*V2(3)
        CROSS_PRODUCT(3)=V1(1)*V2(2)-V1(2)*V2(1)

      END FUNCTION CROSS_PRODUCT

      FUNCTION VABS(V)
        REAL(DOUBLE),DIMENSION(3) :: V
        REAL(DOUBLE) :: VABS
        !
        VABS=SQRT(V(1)**2+V(2)**2+V(3)**2)
      END FUNCTION VABS

END MODULE LinAlg
