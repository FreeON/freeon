!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
PROGRAM SliceAndDIIS
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
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR) &
#else
  TYPE(BCSR)  & 
#endif
                                 :: F,P,E,Tmp1,Tmp2
  TYPE(ARGMT)                    :: Args
  TYPE(INT_VECT)                 :: IWork
  TYPE(DBL_VECT)                 :: V,C,EigenV,Work
  TYPE(DBL_RNK2)                 :: B,BOld,BInv,BTmp
  INTEGER                        :: I,J,K,N,ISCF,JSCF,KSCF,LgN,LWORK,LIWORK,Info
  CHARACTER(LEN=2)               :: Cycl,NxtC
  CHARACTER(LEN=12),PARAMETER    :: Prog='SliceAndDIIS'
  LOGICAL                        :: Present

  INTERFACE DSYEV
     SUBROUTINE DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
         USE GlobalScalars
         CHARACTER(LEN=1), INTENT(IN)    :: JOBZ, UPLO
         INTEGER,          INTENT(IN)    :: LDA,  LWORK, N
         INTEGER,          INTENT(OUT)   :: INFO
         REAL(DOUBLE),     INTENT(INOUT) :: A(LDA,*)
         REAL(DOUBLE),     INTENT(OUT)   :: W(*)
         REAL(DOUBLE),     INTENT(OUT)   :: WORK(*)
     END SUBROUTINE DSYEV
  END INTERFACE
!------------------------------------------------------------------ 
!
!
  CALL StartUp(Args,Prog)
  ISCF=Args%I%I(1)
  Cycl=IntToChar(ISCF)
  NxtC=IntToChar(ISCF+1)
!------------------------------------------
!
!
  CALL New(P)
  CALL New(F)
  CALL New(E)
  CALL New(Tmp1)
!-----------------------------------------------------------------------------
!
!
  CALL Get(P,TrixFile('OrthoD',Args,0))    ! the orthogonalized Fock matrix
  CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonalized Density matrix
  CALL Multiply(F,P,E)  
  CALL Multiply(P,F,E,-One)
! This filter is dodgey, as it can cause premature zeros in B, stalling DIIS.
!  CALL Filter(Tmp1,E)
! Copy for now.
  CALL SetEq(Tmp1,E)
  CALL Put(Tmp1,TrixFile('E',Args,0))
!-------------------------------------------------------------
  CALL New(B,(/ISCF+1,ISCF+1/))
  IF(ISCF>1)THEN
     CALL New(BOld,(/ISCF,ISCF/))
     CALL Get(BOld,'bmat_diis')
!     CALL PPrint(BOld,'BOld',Unit_O=6)
     DO I=1,ISCF
        DO J=1,ISCF
           B%D(I,J)=BOld%D(I,J)
        ENDDO
     ENDDO      
     CALL Delete(BOld)
  ENDIF
  DO JSCF=ISCF-1,1,-1
     KSCF=JSCF-ISCF
     CALL Get(E,TrixFile('E',Args,KSCF))
     B%D(ISCF,JSCF)=Dot(Tmp1,E)
     B%D(JSCF,ISCF)=B%D(ISCF,JSCF)
  ENDDO       
  DO I=1,ISCF+1 
     B%D(ISCF+1,I)=One
     B%D(I,ISCF+1)=One
  ENDDO
  B%D(ISCF,ISCF)=Dot(Tmp1,Tmp1)
  B%D(ISCF+1,ISCF+1)=Zero
  CALL Put(B,'bmat_diis',UnLimit_O=.TRUE.)
!  CALL PPrint(B,'B',Unit_O=6)
!---------------------------------------------------------------
  N=ISCF+1
  CALL New(C,N)
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
     CALL New(V,N)
     CALL New(BInv,(/N,N/))
     CALL New(BTmp,(/N,N/))
     CALL New(EigenV,N)
     LWORK=MAX(1,3*N)
     CALL New(Work,LWORK)
     CALL DSYEV('V','U',N,B%D,N,EigenV%D,Work%D,LWORK,Info)
     IF(Info/=SUCCEED) &
       CALL Halt('DSYEV hosed in DIIS. INFO='//TRIM(IntToChar(Info)))
     CALL Delete(Work)
     V%D=Zero
     V%D(ISCF+1)=One
     C%D=Zero
     BInv%D=Zero
     BTmp%D=Zero
     DO I=1,N
        IF(DABS(EigenV%D(I))<1.D-8)THEN
           BInv%D(I,I)=Zero 
        ELSE
           BInv%D(I,I)=One/EigenV%D(I)
        ENDIF
     ENDDO     
     CALL DGEMM('N','N',N,N,N,One,B%D,N,BInv%D,N,Zero,BTmp%D,N)
     CALL DGEMM('N','T',N,N,N,One,BTmp%D,N,B%D,N,Zero,BInv%D,N)
     CALL DGEMV('N',N,N,One,BInv%D,N,V%D,1,Zero,C%D,1)
     CALL Delete(V)
     CALL Delete(BInv)
     CALL Delete(BTmp)
     CALL Delete(EigenV)
     CALL Put(C%D(ISCF+1),'diiserr',Tag_O='_'//TRIM(CurGeom) &
                                        //'_'//TRIM(CurBase) &
                                        //'_'//TRIM(SCFCycl))
#ifdef PARALLEL
  ENDIF
  CALL BCast(C)
#else
  IF(PrintFlags%Key>=DEBUG_MEDIUM) &
     CALL PPrint(C,'DIIS Coefficients')
#endif
!--------------------------------------------------------------------------------
  CALL Multiply(F,C%D(ISCF))     
  DO JSCF=iSCF-1,1,-1
     KSCF=JSCF-ISCF
     CALL Get(Tmp1,TrixFile('OrthoF',Args,KSCF))
     CALL Multiply(Tmp1,C%D(JSCF))
     CALL Add(F,Tmp1,E)
     CALL Filter(F,E)
  ENDDO
!--------------------------------------------------------------------
!  IO for the non-orthogonal P 
!
   CALL Put(F,TrixFile('F_DIIS',Args,0))
   CALL PChkSum(F,'F_DIIS['//TRIM(Cycl)//']',Prog)
   CALL PPrint(F,'F_DIIS['//TRIM(Cycl)//']')
   CALL Plot(F,'F_DIIS_'//TRIM(Cycl))
!------------------------------
!  Tidy up 
!
!   CALL Delete(F)
!   CALL Delete(P)
!   CALL Delete(E)
!   CALL Delete(Tmp1)
   CALL ShutDown(Prog)   
!------------------------------
END PROGRAM ! DIIS




