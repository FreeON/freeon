MODULE DenMatMethods
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  IMPLICIT NONE

  INTERFACE NormTrace
    MODULE PROCEDURE NormTrace_BCSR
#ifdef PARALLEL
    MODULE PROCEDURE NormTrace_DBCSR
#endif
  END INTERFACE

  INTERFACE FockGuess
    MODULE PROCEDURE FockGuess_BCSR
#ifdef PARALLEL
    MODULE PROCEDURE FockGuess_DBCSR
#endif
  END INTERFACE

  INTERFACE CnvrgChck
#ifdef PARALLEL
    MODULE PROCEDURE CnvrgChck_DBCSR
#endif
    MODULE PROCEDURE CnvrgChck_BCSR
  END INTERFACE

  INTERFACE PutXForm
#ifdef PARALLEL
    MODULE PROCEDURE PutXForm_DBCSR
#endif
    MODULE PROCEDURE PutXForm_BCSR
  END INTERFACE

  REAL(DOUBLE)            :: TrP,TrP2,TrP3,TrP4
  REAL(DOUBLE) :: CurThresh
CONTAINS
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  SUBROUTINE SussTrix(TrixName,Prog)
    CHARACTER(LEN=*) :: TrixName,Prog
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
    REAL(DOUBLE)     :: Trix
!-------------------------------------------------------------------------------
    ! Check for thresholds overide 
    CALL OpenASCII(InpFile,Inp)         
    IF(OptDblQ(Inp,TrixName,Trix))THEN
       Thresholds%Trix=Trix
       Mssg=TRIM(ProcessName(Prog))//' Trix = '  &
          //TRIM(DblToShrtChar(Thresholds%Trix))
       CALL OpenASCII(OutFile,Out)         
       WRITE(Out,*)TRIM(Mssg)
       CLOSE(Out)
    ENDIF
    CLOSE(Inp)
    CALL SetVarThresh()
  END SUBROUTINE SussTrix
!-------------------------------------------------------------------------------

#ifdef PARALLEL
!-------------------------------------------------------------------------------
  SUBROUTINE PutXForm_DBCSR(Prog,Args,P,Z,Tmp1)
    TYPE(DBCSR)       :: P,Z,Tmp1
    CHARACTER(LEN=*) :: Prog
    TYPE(ARGMT)      :: Args
    LOGICAL          :: Present
!-------------------------------------------------------------------------------
    ! IO for the orthogonal P
    CALL Put(P,TrixFile('OrthoD',Args,1))
    ! Archive the orthogonal DM 
    CALL Put(P,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL PChkSum(P,'OrthoP['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( P,'OrthoP['//TRIM(NxtCycl)//']')
    CALL Plot(   P,'OrthoP_'//TRIM(NxtCycl))
    ! Convert to AO representation
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
       CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
       CALL Multiply(Z,P,Tmp1)
       CALL Multiply(Tmp1,Z,P)
    ELSE
       CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
       CALL Multiply(Z,P,Tmp1)
       CALL Get(Z,TrixFile('ZT',Args))
       CALL Multiply(Tmp1,Z,P)
    ENDIF
    CALL Filter(Tmp1,P)     ! Thresholding
    ! IO for the non-orthogonal P
    CALL Put(Tmp1,TrixFile('D',Args,1))
    CALL Put(Zero,'homolumogap')
    CALL PChkSum(Tmp1,'P['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint(Tmp1,'P['//TRIM(NxtCycl)//']')
    CALL Plot(Tmp1,'P_'//TRIM(NxtCycl))
  END SUBROUTINE PutXForm_DBCSR
#endif
!-------------------------------------------------------------------------------
  SUBROUTINE PutXForm_BCSR(Prog,Args,P,Z,Tmp1)
    TYPE(BCSR)       :: P,Z,Tmp1
    CHARACTER(LEN=*) :: Prog
    TYPE(ARGMT)      :: Args
    LOGICAL          :: Present
!-------------------------------------------------------------------------------
    ! IO for the orthogonal P
!    CALL Put(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)
    CALL Put(P,TrixFile('OrthoD',Args,1))
    CALL PChkSum(P,'OrthoP['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( P,'OrthoP['//TRIM(NxtCycl)//']')
    CALL Plot(   P,'OrthoP_'//TRIM(NxtCycl))
    ! Convert to AO representation
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
       CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
       CALL Multiply(Z,P,Tmp1)
       CALL Multiply(Tmp1,Z,P)
    ELSE
       CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
       CALL Multiply(Z,P,Tmp1)
       CALL Get(Z,TrixFile('ZT',Args))
       CALL Multiply(Tmp1,Z,P)
    ENDIF
    CALL Filter(Tmp1,P)     ! Thresholding
    ! IO for the non-orthogonal P
    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,1))
    CALL Put(Zero,'homolumogap')
    CALL PChkSum(Tmp1,'P['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint(Tmp1,'P['//TRIM(NxtCycl)//']')
    CALL Plot(Tmp1,'P_'//TRIM(NxtCycl))
  END SUBROUTINE PutXForm_BCSR
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  SUBROUTINE SetVarThresh(MM_O)
    INTEGER,OPTIONAL  :: MM_O
    INTEGER           :: MM
    REAL(DOUBLE),SAVE :: OldThresh=0D0
    RETURN
    IF(PRESENT(MM_O))THEN
       MM=MM_O
    ELSE
       MM=0
    ENDIF
    IF(OldThresh==0D0)OldThresh=Thresholds%Trix
    Thresholds%Trix=MIN(0.05D0*1.093D0**MM,One)*OldThresh
    CurThresh=OldThresh
  END SUBROUTINE SetVarThresh
!-------------------------------------------------------------------------------

#ifdef PARALLEL
!-------------------------------------------------------------------------------
  FUNCTION CnvrgChck_DBCSR(Prog,NPur,Ne,MM,F,P,POld,Tmp1,Tmp2)
    LOGICAL              :: CnvrgChck_DBCSR
    TYPE(DBCSR)          :: F,P,POld,Tmp1,Tmp2
    REAL(DOUBLE)         :: Ne,Energy,AbsErrP,FNormErrP,TwoNP,N2F,  &
         AbsErrE,RelErrE,AveErrE,MaxCommErr,FNormCommErr,PNon0,TraceP

    REAL(DOUBLE),DIMENSION(2) :: CErr
    REAL(DOUBLE),SAVE    :: OldE,OldAEP
    INTEGER              :: MM,NPur
    CHARACTER(LEN=*)     :: Prog
    CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg,CnvrgCmmnt
#ifdef PRINT_PURE_EVALS
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
        TYPE(DBL_RNK2)                 :: dP
        TYPE(DBL_VECT)                 :: EigenV,Work
        TYPE(INT_VECT)                 :: IWork
        INTEGER                        :: LWORK,LIWORK,Info

        CALL New(dP,(/NBasF,NBasF/))
        CALL SetEq(dP,P)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dP%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in FockGuess. INFO='//TRIM(IntToChar(Info)))
        PrintFlags%Fmt=DEBUG_MMASTYLE
        CALL Print_DBL_VECT(EigenV,'Values['//TRIM(IntToChar(NPur))//']',Unit_O=6)
        PrintFlags%Fmt=DEBUG_DBLSTYLE
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dP)
#endif

     IF(NPur==0)THEN
        OldE=BIG_DBL
        OldAEP=BIG_DBL
     ENDIF
     PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
     ! Density matrix errors
     CALL Multiply(Pold,-One)
     CALL Add(Pold,P,Tmp1)
     AbsErrP=ABS(Max(Tmp1)+1.D-20)
     FNormErrP=FNorm(Tmp1)
     ! Energy errors

     CALL Multiply(P,F,Tmp1)
     Energy=Trace(Tmp1)    

     AbsErrE=ABS(OldE-Energy)
     RelErrE=AbsErrE/ABS(Energy)
     ! Convergence check
     CnvrgChck_DBCSR=.FALSE.
     ! Absolute convergence test
     IF(RelErrE<Thresholds%ETol*1D-2.AND. &
        AbsErrP<Thresholds%DTol*1D-1)THEN
        CnvrgChck_DBCSR=.TRUE.
        CnvrgCmmnt='Met dE/dP goals'
     ENDIF
     ! Test in the asymptotic regime for stall out
     IF(RelErrE<1D2*Thresholds%Trix**2)THEN
!    IF(RelErrE<Thresholds%ETol)THEN
        ! Check for increasing /P
        IF(AbsErrP>OldAEP)THEN
           CnvrgChck_DBCSR=.TRUE.
           CnvrgCmmnt='Hit dP increase'
        ENDIF
        ! Check for an increasing energy
        IF(Energy>OldE)THEN
           CnvrgChck_DBCSR=.TRUE.
           CnvrgCmmnt='Hit dE increase'
        ENDIF
     ENDIF


     CALL BCast(CnvrgChck_DBCSR)

!     IF(NPur<35)CnvrgChck_DBCSR=.FALSE.
!     CALL OpenASCII('CommErr_'//TRIM(DblToShrtChar(Thresholds%Trix))//'.dat',77)
!     CALL OpenASCII('CommErr_'//TRIM(DblToShrtChar(CurThresh))//'.dat',77)
!     CErr=CommutatorErrors(F,P)
!     WRITE(77,22)NPur,Thresholds%Trix,AbsErrP,FNormErrP,CErr(1),CErr(2), &
!                 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
!     22 FORMAT(I3,8(1x,F20.14))
!     CLOSE(77)

     ! Updtate previous cycle values
     OldE=Energy
     OldAEP=AbsErrP
     ! Print convergence stats
     Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))      &
          //'dE='//TRIM(DblToShrtChar(RelErrE))                 &
          //', dP='//TRIM(DblToShrtChar(AbsErrP))               &
          //', %Non0='//TRIM(DblToShrtChar(PNon0))              

     IF(MyId==ROOT)THEN

        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(*,*)TRIM(Mssg)
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF

     ENDIF

     ! Set thresholding for next cycle
     CALL SetVarThresh(MM)
     ! Look for convergence
     IF(.NOT.CnvrgChck_DBCSR)THEN
        CALL SetEq(Pold,P)
        RETURN
     ENDIF
     ! Normalize Trace
     CALL NormTrace(P,Tmp2,Tmp1,Ne,1)
     MM=MM+1
#ifdef COMPUTE_COMMUTATOR
     ! Commutator [F,P] 
     N2F=FNorm(F)
     CALL Multiply(F,P,Tmp1)
     CALL Multiply(P,F,POld)
     CALL Multiply(POld,-One)
     CALL Add(Tmp1,POld,F) 
     MM=MM+2
     FNormCommErr=FNorm(F)
     MaxCommErr=Max(F)
#endif
     ! Print summary stats
     TraceP = Trace(P)

     IF(MyId==ROOT)THEN

        IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           Mssg=ProcessName(Prog,CnvrgCmmnt) &
                //'Tr{FP}='//TRIM(DblToChar(Energy)) &
                //', dNel = '//TRIM(DblToShrtChar(Two*ABS(TraceP-Ne))) 
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)

           Mssg=ProcessName(Prog)//TRIM(IntToChar(NPur))//' purification steps, ' &
                //TRIM(IntToChar(MM))//' matrix multiplies'
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog)//'Fractional occupation = '              &
                //TRIM(DblToShrtChar(Half*DBLE(NEl)/DBLE(NBasF)))          &
                //', ThrX='//TRIM(DblToShrtChar(Thresholds%Trix))          &
                //', %Non0s = '//TRIM(DblToShrtChar(PNon0))              
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog,'Max abs errors') &
                //'dE='//TRIM(DblToShrtChar(AbsErrE))//', '                 &
                //'dP='//TRIM(DblToShrtChar(AbsErrP))
#ifdef COMPUTE_COMMUTATORS
           Mssg=TRIM(Mssg)//', '//'[F,P]='//TRIM(DblToShrtChar(MaxCommErr))
#endif
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog) &
                //'Rel dE='//TRIM(DblToShrtChar(RelErrE))//', '                &
                //'||dP||_F='//TRIM(DblToShrtChar(FNormErrP))
#ifdef COMPUTE_COMMUTATORS
           Mssg=TRIM(Mssg)//', '//'||[F,P]||_F='//TRIM(DblToShrtChar(FNormCommErr))
#endif
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF

     ENDIF

   END FUNCTION CnvrgChck_DBCSR
#endif

!-------------------------------------------------------------------------------
  FUNCTION CnvrgChck_BCSR(Prog,NPur,Ne,MM,F,P,POld,Tmp1,Tmp2)
    LOGICAL              :: CnvrgChck_BCSR
    TYPE(BCSR)           :: F,P,POld,Tmp1,Tmp2
    REAL(DOUBLE)         :: Ne,Energy,AbsErrP,FNormErrP,TwoNP,N2F,  &
         AbsErrE,RelErrE,AveErrE,MaxCommErr,FNormCommErr,PNon0,TraceP

    REAL(DOUBLE),DIMENSION(2) :: CErr
    REAL(DOUBLE),SAVE    :: OldE,OldAEP
    INTEGER              :: MM,NPur
    CHARACTER(LEN=*)     :: Prog
    CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg,CnvrgCmmnt
#ifdef PRINT_PURE_EVALS
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
        TYPE(DBL_RNK2)                 :: dP
        TYPE(DBL_VECT)                 :: EigenV,Work
        TYPE(INT_VECT)                 :: IWork
        INTEGER                        :: LWORK,LIWORK,Info

        CALL New(dP,(/NBasF,NBasF/))
        CALL SetEq(dP,P)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dP%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in FockGuess. INFO='//TRIM(IntToChar(Info)))
        PrintFlags%Fmt=DEBUG_MMASTYLE
        CALL Print_DBL_VECT(EigenV,'Values['//TRIM(IntToChar(NPur))//']',Unit_O=6)
        PrintFlags%Fmt=DEBUG_DBLSTYLE
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dP)
#endif

     IF(NPur==0)THEN
        OldE=BIG_DBL
        OldAEP=BIG_DBL
     ENDIF
     PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
     ! Density matrix errors
     CALL Multiply(Pold,-One)
     CALL Add(Pold,P,Tmp1)
     AbsErrP=ABS(Max(Tmp1)+1.D-20)
     FNormErrP=FNorm(Tmp1)
     ! Energy errors

     Energy=Trace(P,F)

     AbsErrE=ABS(OldE-Energy)
     RelErrE=AbsErrE/ABS(Energy)
     ! Convergence check
     CnvrgChck_BCSR=.FALSE.
     ! Absolute convergence test
     IF(RelErrE<Thresholds%ETol*1D-2.AND. &
        AbsErrP<Thresholds%DTol*1D-1)THEN
        CnvrgChck_BCSR=.TRUE.
        CnvrgCmmnt='Met dE/dP goals'
     ENDIF
     ! Test in the asymptotic regime for stall out
     IF(RelErrE<1D2*Thresholds%Trix**2)THEN
!    IF(RelErrE<Thresholds%ETol)THEN
        ! Check for increasing /P
        IF(AbsErrP>OldAEP)THEN
           CnvrgChck_BCSR=.TRUE.
           CnvrgCmmnt='Hit dP increase'
        ENDIF
        ! Check for an increasing energy
        IF(Energy>OldE)THEN
           CnvrgChck_BCSR=.TRUE.
           CnvrgCmmnt='Hit dE increase'
        ENDIF
     ENDIF

!     IF(NPur<35)CnvrgChck_BCSR=.FALSE.
!     CALL OpenASCII('CommErr_'//TRIM(DblToShrtChar(Thresholds%Trix))//'.dat',77)
!     CALL OpenASCII('CommErr_'//TRIM(DblToShrtChar(CurThresh))//'.dat',77)
!     CErr=CommutatorErrors(F,P)
!     WRITE(77,22)NPur,Thresholds%Trix,AbsErrP,FNormErrP,CErr(1),CErr(2), &
!                 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
!     22 FORMAT(I3,8(1x,F20.14))
!     CLOSE(77)

     ! Updtate previous cycle values
     OldE=Energy
     OldAEP=AbsErrP
     ! Print convergence stats
     Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))      &
          //'dE='//TRIM(DblToShrtChar(RelErrE))                 &
          //', dP='//TRIM(DblToShrtChar(AbsErrP))               &
          //', %Non0='//TRIM(DblToShrtChar(PNon0))              
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(*,*)TRIM(Mssg)
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
     ! Set thresholding for next cycle
     CALL SetVarThresh(MM)
     ! Look for convergence
     IF(.NOT.CnvrgChck_BCSR)THEN
        CALL SetEq(Pold,P)
        RETURN
     ENDIF
     ! Normalize Trace
     CALL NormTrace(P,Tmp2,Tmp1,Ne,1)
     MM=MM+1
#ifdef COMPUTE_COMMUTATOR
     ! Commutator [F,P] 
     N2F=FNorm(F)
     CALL Multiply(F,P,Tmp1)
     CALL Multiply(P,F,POld)
     CALL Multiply(POld,-One)
     CALL Add(Tmp1,POld,F) 
     MM=MM+2
     FNormCommErr=FNorm(F)
     MaxCommErr=Max(F)
#endif
     ! Print summary stats
     TraceP = Trace(P)
        IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           Mssg=ProcessName(Prog,CnvrgCmmnt) &
                //'Tr{FP}='//TRIM(DblToChar(Energy)) &
                //', dNel = '//TRIM(DblToShrtChar(Two*ABS(TraceP-Ne))) 
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)

           Mssg=ProcessName(Prog)//TRIM(IntToChar(NPur))//' purification steps, ' &
                //TRIM(IntToChar(MM))//' matrix multiplies'
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog)//'Fractional occupation = '              &
                //TRIM(DblToShrtChar(Half*DBLE(NEl)/DBLE(NBasF)))          &
                //', ThrX='//TRIM(DblToShrtChar(Thresholds%Trix))          &
                //', %Non0s = '//TRIM(DblToShrtChar(PNon0))              
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog,'Max abs errors') &
                //'dE='//TRIM(DblToShrtChar(AbsErrE))//', '                 &
                //'dP='//TRIM(DblToShrtChar(AbsErrP))
#ifdef COMPUTE_COMMUTATORS
           Mssg=TRIM(Mssg)//', '//'[F,P]='//TRIM(DblToShrtChar(MaxCommErr))
#endif
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog) &
                //'Rel dE='//TRIM(DblToShrtChar(RelErrE))//', '                &
                //'||dP||_F='//TRIM(DblToShrtChar(FNormErrP))
#ifdef COMPUTE_COMMUTATORS
           Mssg=TRIM(Mssg)//', '//'||[F,P]||_F='//TRIM(DblToShrtChar(FNormCommErr))
#endif
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
   END FUNCTION CnvrgChck_BCSR

   FUNCTION CommutatorErrors(F,P) RESULT(CErrs)
     REAL,DIMENSION(2) :: CErrs
     TYPE(BCSR)        :: F,P,T1,T2,T3
     CALL New(T1)
     CALL New(T2)
     CALL New(T3)
     CALL Multiply(F,P,T1)
     CALL Multiply(P,F,T2)
     CALL Multiply(T2,-One)
     CALL Add(T1,T2,T3)
     CErrs(1)=MAX(T3)
     CErrs(2)=FNorm(T3)
     CALL Delete(T1)
     CALL Delete(T2)
     CALL Delete(T3)
   END FUNCTION CommutatorErrors
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
    SUBROUTINE SP2(P,P2,Tmp1,Norm,MMMs)
#ifdef PARALLEL
      TYPE(DBCSR)   :: P,P2,Tmp1
#else
      TYPE(BCSR)   :: P,P2,Tmp1
#endif
      REAL(DOUBLE) :: Norm,CR
      INTEGER      :: MMMs
!-------------------------------------------------------------------------------
      CALL Multiply(P,P,P2)             ! The only multiplication is a square
      MMMs=MMMs+1
      TrP=Trace(P)
      CR=TrP-Norm                       ! CR = Occupation error criteria
      IF (CR >  0) THEN                 ! Too many states
         CALL Filter(P,P2)              ! P = P^2
      ELSE                              ! Too few states 
         CALL Multiply(P ,Two) 
         CALL Multiply(P2,-One)
         CALL Add(P,P2,Tmp1)           ! P = 2P-P^2
         CALL Filter(P,Tmp1)
      ENDIF
    END SUBROUTINE SP2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
    SUBROUTINE SP4(P,P2,Tmp1,Tmp2,Norm,MMMs)
      TYPE(BCSR)                     :: P,P2,Tmp1,Tmp2
      REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Thresh_old
      INTEGER                        :: MMMs
      REAL(DOUBLE)                   :: EPS = 1.D-12
!-------------------------------------------------------------------------------
      CALL Multiply(P,P,P2)
      MMMs = MMMs+1
      TrP  = Trace(P)
      TrP2 = Trace(P2)
      TrP3 = Trace(P,P2)
      TrP4 = Trace(P2,P2)
      Gt   = (Norm-Four*TrP3+3*TrP4)
      Gn   = (TrP2-Two*TrP3+TrP4)
      CR   = TrP - Norm                   ! Trace correction measure
      IF(ABS(Gn).LT.EPS)THEN              ! Close to idempotency
         G=Three                          ! To avoid numerical errors
      ELSE
         G=Gt/Gn                          ! Boundary measure
      ENDIF
      IF((G>Zero).AND.(G<Six))THEN        ! Check the bounds
         IF (CR  > Zero) THEN             ! Too many states
           CALL Filter(Tmp2,P2)
           CALL Multiply(P,4.d0)
           CALL Multiply(P2,-3.d0)
           CALL Add(P,P2,Tmp1)
           CALL Filter(P,Tmp1)
           CALL Multiply(Tmp2,P,Tmp1)     ! P = P^2*(4P-3P^2)
           MMMs = MMMs+1
        ELSE                              ! Too few states
           CALL Filter(Tmp2,P2)
           CALL Multiply(P ,-8.d0)
           CALL Multiply(P2, 3.d0)
           CALL Add(P,P2,Tmp1)
           CALL Add(Tmp1,Six)
           CALL Filter(P,Tmp1)
           CALL Multiply(Tmp2,P,Tmp1)     ! P = P^2*(6I-8P+3P^2)
           MMMs = MMMs+1
        ENDIF
        CALL Filter(P,Tmp1)

     ELSE 
        IF(G<0) THEN                      ! Too many states
           CALL Filter(P,P2)              ! P = P^2
        ELSE                              ! Too few states 
           CALL Multiply(P ,Two) 
           CALL Multiply(P2,-One)
           CALL Add(P ,P2,Tmp1)           ! P = 2P-P^2
           CALL Filter(P,Tmp1)
        ENDIF
     ENDIF
   END SUBROUTINE SP4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE TS4(P,P2,Tmp1,Tmp2,Norm,MMMs)
     TYPE(BCSR)                     :: P,P2,Tmp1,Tmp2
     REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Coeff
     INTEGER                        :: MMMs
     REAL(DOUBLE)                   :: EPS = 1.D-12
!-------------------------------------------------------------------------------
     CALL Multiply(P,P,P2)    
     MMMs = MMMs+1
     TrP  = Trace(P)
     TrP2 = Trace(P2)
     TrP3 = Trace(P,P2)
     TrP4 = Trace(P2,P2)
     Gt   = (Norm-Four*TrP3+3*TrP4)
     Gn   = (TrP2-Two*TrP3+TrP4)
     CR   = TrP - Norm     
     IF(ABS(Gn).LT.EPS)THEN 	              ! Close to idempotency
        G=Three 	        	      ! To avoid numerical errors
     ELSE
        G=Gt/Gn	                       	      ! Boundary measure
     ENDIF
     IF ((G>Zero).AND.(G<Six))THEN 	      ! Check the bounds
        CALL Filter(Tmp2,P2)
        Coeff = Four-Two*G
        CALL Multiply(P,Coeff)
        Coeff = G-Three
        CALL Multiply(P2,Coeff)
        CALL Add(P,P2,Tmp1)
        CALL Add(Tmp1,G)
        CALL Filter(P,Tmp1)
        CALL Multiply(Tmp2,P,Tmp1)         ! P^2*(G+(4-2*G)*P+(G-3)*P^2)
        CALL Filter(P,Tmp1)
        MMMs = MMMs+1
     ELSE
        IF (G < 0) THEN                      ! Too many states
           CALL Filter(P,P2)                 ! P = P^2
        ELSE                                 ! Too few states 
           CALL Multiply(P ,Two) 
           CALL Multiply(P2,-One)
           CALL Add(P,P2,Tmp1)             ! P = 2P-P^2
           CALL Filter(P,Tmp1)
        ENDIF
     ENDIF
   END SUBROUTINE TS4
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE PM1(P,P2,P3,Tmp1,Tmp2,Norm,MMMs)
     TYPE(BCSR)                     :: P,P2,P3,Tmp1,Tmp2
     REAL(DOUBLE)                   :: CR,Norm,Coeff
     INTEGER                        :: MMMs
     REAL(DOUBLE)                   :: EPS = 1.D-8
!-------------------------------------------------------------------------------
     CALL Multiply(P,P,P2)                   
     MMMs=MMMs+1
     CALL Filter(Tmp1,P2)
     CALL Multiply(P,Tmp1,P3)               
     MMMs=MMMs+1
     TrP=Trace(P)
     TrP2=Trace(P2)
     TrP3=Trace(P3)
     IF(ABS(TrP-TrP2)<EPS)THEN                ! Close to idempotency
       CALL Multiply(P2, Three)
       CALL Multiply(P3,-Two)
       CALL Add(P2,P3,Tmp1)                   ! T = 3P^2-2P^2 (McWeeny close to convergence)
     ELSE
        CR = (TrP2 - TrP3)/(TrP-TrP2)
        IF (CR <= Half) THEN
           Coeff = (One-Two*CR)/(One-CR)
           CALL Multiply(P,Coeff)
           Coeff = (One+CR)/(One-CR)
           CALL Multiply(P2,Coeff)
           Coeff = -One/(One-CR)
           CALL Multiply(P3,Coeff)
           CALL Add(P,P2,Tmp2)
           CALL Add(Tmp2,P3,Tmp1)             ! T = [(1-2CR)*P+(1+CR)*P^2-P^3]/(1-CR)
        ELSE
           Coeff = (One+CR)/CR
           CALL Multiply(P2,Coeff)
           Coeff = -One/CR
           CALL Multiply(P3,Coeff)
           CALL Add(P2,P3,Tmp1)               ! T = [(1+CR)*P^2-P^3]/CR
        ENDIF
     ENDIF
     CALL Filter(P,Tmp1)
   END SUBROUTINE PM1

   SUBROUTINE PM2(P,P2,P3,Tmp1,MMMs)
     TYPE(BCSR)    :: P,P2,P3,Tmp1
     REAL(DOUBLE)  :: c,u,v,w,TrP1
     INTEGER       :: MMMs
     LOGICAL,SAVE  :: FixedUVW=.FALSE.
!-------------------------------------------------------------------------------
      CALL Multiply(P,P,Tmp1) 
      MMMs=MMMs+1
      CALL Filter(P2,Tmp1)    
      CALL Multiply(P2,P,Tmp1) 
      MMMs=MMMs+1
      CALL Filter(P3,Tmp1)    
      TrP1=Trace(P)
      TrP2=Trace(P2)
      TrP3=Trace(P3)
      c=(TrP2-TrP3)/(TrP1-TrP2+1D-30)
WRITE(*,*)' C = ',C
      IF(ABS(c-Half)<1.D-5.OR.FixedUVW)THEN
         c=Half
         u=Zero
         v=Three
         w=-Two
         FixedUVW=.TRUE.
      ELSEIF(c<=Half)THEN
         u=(One-Two*c)/(One-c+1D-30) 
         v=(One+c)/(One-c+1D-30) 
         w=-One/(One-c+1D-30)
      ELSE
         u=Zero
         v=(One+c)/(c+1D-30)     
         w=-One/(c+1D-30)
      ENDIF
!     Assemble purified P
      CALL Multiply(P ,u)
      CALL Multiply(P2,v)
      CALL Multiply(P3,w)
      CALL Add(P,P2,Tmp1)
      CALL Add(Tmp1,P3,P2) !     P[J+1] = u*P[J] + v*P[J].P[J] + w*P[J].P[J].P[J]
      CALL Filter(P,P2)    !     P=Filter[P[N+1,I+1]]
   END SUBROUTINE PM2
!-------------------------------------------------------------------------------


#ifdef PARALLEL
!-------------------------------------------------------------------------------
    SUBROUTINE FockGuess_DBCSR(F,P,Norm,Order)
      TYPE(DBCSR)                    :: P
      TYPE(BCSR)                     :: F
      REAL(DOUBLE)                   :: Fmin,Fmax,DF,Coeff,Mu,Lmbd1,  &
                                        Lmbd2,Norm
      INTEGER                        :: Order
#ifdef PRINT_GUESS_EVALS
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
        TYPE(DBL_RNK2)                 :: dP
        TYPE(DBL_VECT)                 :: EigenV,Work
        TYPE(INT_VECT)                 :: IWork
        INTEGER                        :: LWORK,LIWORK,Info
#endif
!-------------------------------------------------------------------------------
!     Estimate spectral bounds 
      IF(MyId==ROOT) &
        CALL SpectralBounds(F,Fmin,Fmax)

      CALL BCast(FMin)
      CALL BCast(FMax)

!     Set up the Density Matrix
     IF(Order==1) THEN
         Call SetEq(P,F)
         DF = (Fmax - Fmin)
         CALL Add(P,-Fmax)
         Coeff = -One/DF
         CALL Multiply(P,Coeff)          ! P = (I*F_max-F)/DF
#ifdef PRINT_GUESS_EVALS
        CALL New(dP,(/NBasF,NBasF/))
        CALL SetEq(dP,P)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dP%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in FockGuess. INFO='//TRIM(IntToChar(Info)))
        CALL Print_DBL_VECT(EigenV,'Guess1Values',Unit_O=6)
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dP)
#endif
      ELSEIF(Order==2) THEN
         Mu    = Trace(F)/DBLE(NBasF)
         Lmbd1 = Norm/(Fmax-Mu)
         Lmbd2 = (DBLE(NBasF)-Norm)/(Mu-Fmin)
         IF(Lmbd1 < Lmbd2) THEN
            CALL SetEq(P,F)
            CALL Add(P,-Mu)
            Coeff = -Lmbd1/DBLE(NBasF)
            CALL Multiply(P,Coeff) 
            Coeff = Norm/DBLE(NBasF)
            CALL Add(P,Coeff) 
         ELSE
            CALL SetEq(P,F)
            CALL Add(P,-Mu)
            Coeff = -Lmbd2/DBLE(NBasF)
            CALL Multiply(P,Coeff) 
            Coeff = Norm/DBLE(NBasF)
            CALL Add(P,Coeff)
         ENDIF
#ifdef PRINT_GUESS_EVALS
        CALL New(dP,(/NBasF,NBasF/))
        CALL SetEq(dP,P)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dP%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in FockGuess. INFO='//TRIM(IntToChar(Info)))
        CALL Print_DBL_VECT(EigenV,'Guess2Values',Unit_O=6)
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dP)
#endif
      ELSE
         CALL MondoHalt(99,'Wrong Order in FockGuess')
      ENDIF
    END SUBROUTINE FockGuess_DBCSR
#endif

!-------------------------------------------------------------------------------
    SUBROUTINE FockGuess_BCSR(F,P,Norm,Order)
      TYPE(BCSR)                     :: F,P
      REAL(DOUBLE)                   :: Fmin,Fmax,DF,Coeff,Mu,Lmbd1,  &
                                        Lmbd2,Norm
      INTEGER                        :: Order
#ifdef PRINT_GUESS_EVALS
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
        TYPE(DBL_RNK2)                 :: dP
        TYPE(DBL_VECT)                 :: EigenV,Work
        TYPE(INT_VECT)                 :: IWork
        INTEGER                        :: LWORK,LIWORK,Info
#endif
!-------------------------------------------------------------------------------
!     Estimate spectral bounds 
      CALL SpectralBounds(F,Fmin,Fmax)

!     Set up the Density Matrix
     IF(Order==1) THEN
         Call SetEq(P,F)
         DF = (Fmax - Fmin)
         CALL Add(P,-Fmax)
         Coeff = -One/DF
         CALL Multiply(P,Coeff)          ! P = (I*F_max-F)/DF
#ifdef PRINT_GUESS_EVALS
        CALL New(dP,(/NBasF,NBasF/))
        CALL SetEq(dP,P)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dP%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in FockGuess. INFO='//TRIM(IntToChar(Info)))
        CALL Print_DBL_VECT(EigenV,'Guess1Values',Unit_O=6)
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dP)
#endif
      ELSEIF(Order==2) THEN
         Mu    = Trace(F)/DBLE(NBasF)
         Lmbd1 = Norm/(Fmax-Mu)
         Lmbd2 = (DBLE(NBasF)-Norm)/(Mu-Fmin)
         IF(Lmbd1 < Lmbd2) THEN
            CALL SetEq(P,F)
            CALL Add(P,-Mu)
            Coeff = -Lmbd1/DBLE(NBasF)
            CALL Multiply(P,Coeff) 
            Coeff = Norm/DBLE(NBasF)
            CALL Add(P,Coeff) 
         ELSE
            CALL SetEq(P,F)
            CALL Add(P,-Mu)
            Coeff = -Lmbd2/DBLE(NBasF)
            CALL Multiply(P,Coeff) 
            Coeff = Norm/DBLE(NBasF)
            CALL Add(P,Coeff)
         ENDIF
#ifdef PRINT_GUESS_EVALS
        CALL New(dP,(/NBasF,NBasF/))
        CALL SetEq(dP,P)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dP%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in FockGuess. INFO='//TRIM(IntToChar(Info)))
        CALL Print_DBL_VECT(EigenV,'Guess2Values',Unit_O=6)
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dP)
#endif
      ELSE
         CALL MondoHalt(99,'Wrong Order in FockGuess')
      ENDIF
    END SUBROUTINE FockGuess_BCSR

!-------------------------------------------------------------------------------
! Estimate spectral bounds via Gersgorin approximation, [F_min-F_max].
!-------------------------------------------------------------------------------
!   S(R) = Sum_R,C { ABS(F(R,C) }
!   F_max = Max_R { F(R,R) + S(R) - ABS(F(R,R)) }
!   F_min = Min_R { F(R,R) - S(R) + ABS(F(R,R)) }
!-------------------------------------------------------------------------------
      SUBROUTINE SpectralBounds(F,F_min,F_max)
        TYPE(BCSR)       :: F
        INTEGER          :: I,R,J,M,Col,Blk,N,C,Check
        REAL(DOUBLE)     :: F_min,F_max,Tmp_max,Tmp_min,Diag,Sum
#ifdef EXACT_EIGEN_VALUES
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
        TYPE(DBL_RNK2)                 :: dF
        TYPE(DBL_VECT)                 :: EigenV,Work
        TYPE(INT_VECT)                 :: IWork
        INTEGER                        :: LWORK,LIWORK,Info
!-------------------------------------------------------------------------------
        CALL New(dF,(/NBasF,NBasF/))
        CALL SetEq(dF,F)
        CALL New(EigenV,NBasF)
        CALL SetEq(EigenV,Zero)
        LWORK=MAX(1,3*NBasF+10)
        CALL New(Work,LWork)
        CALL DSYEV('V','U',NBasF,dF%D,NBasF,EigenV%D,Work%D,LWORK,Info)
        IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in RHEqs. INFO='//TRIM(IntToChar(Info)))
        CALL PPrint(EigenV,'FockValues',Unit_O=6)
        F_Min=EigenV%D(1)
        F_Max=EigenV%D(NBasF)
        CALL Delete(EigenV)
        CALL Delete(Work)
        CALL Delete(dF)
#else
        F_min =  1.0D10
        F_max = -1.0D10 
        DO I = 1,F%NAtms                                   ! Step over row blocks
           M = BSiz%I(I)
           DO R = 0,M-1                                    ! Step over rows in each block
              Sum = Zero
              Check = 0
              DO J = F%RowPt%I(I),F%RowPt%I(I+1)-1         ! Step over column blocks
                 Col = F%ColPt%I(J)
                 Blk = F%BlkPt%I(J)
                 N = BSiz%I(Col)
                 DO C = 0,N-1                              ! Step over colums in each block
                    Sum = Sum + ABS(F%MTrix%D(Blk+C*M+R))  ! Assume col by col storage?
                    IF(I.EQ.Col) THEN
                       IF(R.EQ.C) THEN
                          Diag = F%MTrix%D(Blk+C*M+R)      ! Diagonal elements
                          Check = 1
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO
              IF (Check.EQ.0) THEN                         ! In case there wasn't a diagonal element
                 Diag = Zero
              ENDIF
              Tmp_max = Diag + Sum - ABS(Diag)             ! Upper bound of Grsg. circle
              Tmp_min = Diag - Sum + ABS(Diag)             ! Lower bound of Grsg. circle
              IF(F_max .LE.Tmp_max) THEN                   ! Check for the highest bound
                 F_max = Tmp_max
              ENDIF
              IF (F_min.GE.Tmp_min) THEN                   ! Check for the lowest bound
                 F_min = Tmp_min
              ENDIF
           ENDDO
        ENDDO
#endif
      END SUBROUTINE SpectralBounds
!-------------------------------------------------------------------------------

#ifdef PARALLEL
!-------------------------------------------------------------------------------
   SUBROUTINE NormTrace_DBCSR(P,P2,Tmp1,Norm,Order)
     TYPE(DBCSR)    :: P,P2,Tmp1
     REAL(DOUBLE)  :: Norm,CR,C1,C2
     INTEGER       :: Order
!-------------------------------------------------------------------------------
     IF(Order==0) THEN
        TrP  = Trace(P)
        CR   = Norm/TrP
        CALL Multiply(P,CR)
     ELSEIF(Order==1) THEN
        CALL Multiply(P,P,P2)
        TrP  = Trace(P)
        TrP2 = Trace(P2)
        CR   = TrP - TrP2
        IF(CR==Zero) RETURN
        C1   = (Norm - TrP2)/CR
        C2   = (TrP- Norm)/CR
        CALL Multiply(P  ,C1)
        CALL Multiply(P2 ,C2)
        CALL Add(P,P2,Tmp1)
        CALL Filter(P,Tmp1)
     ELSE
        CALL MondoHalt(99,'Wrong Order in NormTrace')
     ENDIF
   END SUBROUTINE NormTrace_DBCSR
!-------------------------------------------------------------------------------
#endif

!-------------------------------------------------------------------------------
   SUBROUTINE NormTrace_BCSR(P,P2,Tmp1,Norm,Order)
     TYPE(BCSR)    :: P,P2,Tmp1
     REAL(DOUBLE)  :: Norm,CR,C1,C2
     INTEGER       :: Order
!-------------------------------------------------------------------------------
     IF(Order==0) THEN
        TrP  = Trace(P)
        CR   = Norm/TrP
        CALL Multiply(P,CR)
     ELSEIF(Order==1) THEN
        CALL Multiply(P,P,P2)
        TrP  = Trace(P)
        TrP2 = Trace(P2)
        CR   = TrP - TrP2
        IF(CR==Zero) RETURN
        C1   = (Norm - TrP2)/CR
        C2   = (TrP- Norm)/CR
        CALL Multiply(P  ,C1)
        CALL Multiply(P2 ,C2)
        CALL Add(P,P2,Tmp1)
        CALL Filter(P,Tmp1)
     ELSE
        CALL MondoHalt(99,'Wrong Order in NormTrace')
     ENDIF
   END SUBROUTINE NormTrace_BCSR
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   SUBROUTINE CalculateDegen(Norm,Degen,Occpan)
     REAL(DOUBLE)                 :: Norm,Degen,Occpan

     Degen     = ((Norm-TrP2)**3)/((TrP2-TrP3)*(Norm-Two*TrP2+TrP3))
     Occpan    = Two*(TrP2-TrP3)/(Norm-TrP2)  

   END SUBROUTINE CalculateDegen
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
#ifdef BROKEN_GAP
   SUBROUTINE CalculateGap(F,P,Norm,lumo_occ,Gap)
     TYPE(BCSR)                   :: F,P
     REAL(DOUBLE)                 :: Norm,Gap,lumo_occ
     REAL(DOUBLE)                 :: TrFP,TrFP2,TrFP3

     CALL Multiply(P,P,P2)
     CALL Multiply(P,P2,P3)
     TrFP  = Trace(F,P)
     TrFP2 = Trace(F,P2)
     TrFP3 = Trace(F,P3)    
! 
     lumo_occ = Half + Half*Sqrt(ABS(One-Two*Norm+Two*TrP2)) 

     Gap      = -(TrFP - Three*TrFP2 + Two*TrFP2)/(lumo_occ*(lumo_occ-One)*(Two*lumo_occ-One))

   END SUBROUTINE CalculateGap
#endif
END MODULE DenMatMethods
