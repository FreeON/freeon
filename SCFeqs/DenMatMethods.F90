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
  REAL(DOUBLE)            :: TrP,TrP2,TrP3,TrP4
  TYPE(BCSR)              :: P2,P3,Ptmp1,Ptmp2
  REAL(DOUBLE),PARAMETER  :: ThFAC1 = 0.125D0
  REAL(DOUBLE),PARAMETER  :: ThFAC2 = 1.D0/ThFAC1
!
!  WE NEED TO DEAL WITH THE ALLOCATION OF TEMPORARY MATRICES, EITHER DECLARE THEM
!  GLOBALLY (HERE FOR EXAMPLE), OR PASS THE MEMORY IN.  IMPLICIT ALLOCATION EXPLICIT
!  DEALLOCATION SHOULD BE AVOIDED...
!
  CONTAINS
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE Commute(F,P,T)
      TYPE(BCSR)                     :: F,P,T
!
      CALL Multiply(F,P,Ptmp1)
      CALL Multiply(P,F,Ptmp2)
      CALL Multiply(Ptmp2,-One)
      CALL Add(Ptmp1,Ptmp2,T)
!
    END SUBROUTINE Commute
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE FockGuess(F,P,Norm,Order)
      TYPE(BCSR)                     :: F,P
      REAL(DOUBLE)                   :: Fmin,Fmax,DF,Coeff,Mu,Lmbd1,  &
                                        Lmbd2,Norm
      INTEGER                        :: Order
!----------------------------------------------------------------------------
!     Estimate spectral bounds 
!----------------------------------------------------------------------------
      CALL SpectralBounds(F,Fmin,Fmax)
!----------------------------------------------------------------------------
!     Set up the Density Matrix
!----------------------------------------------------------------------------
      IF(Order==1) THEN
         Call SetEq(P,F)
         DF = (Fmax - Fmin)
         CALL Add(P,-Fmax)
         Coeff = -One/DF
         CALL Multiply(P,Coeff)          ! P = (I*F_max-F)/DF
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
      ELSE
         CALL MondoHalt(99,'Wrong Order in FockGuess')
      ENDIF
!
    END SUBROUTINE FockGuess
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE NormTrace(P,Norm,Order)
     TYPE(BCSR)                     :: P
     REAL(DOUBLE)                   :: Norm,CR,C1,C2
     INTEGER                        :: Order
! 
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
        CALL Add(P,P2,Ptmp1)
        CALL Filter(P,Ptmp1)
     ELSE
        CALL MondoHalt(99,'Wrong Order in NormTrace')
     ENDIF
!
   END SUBROUTINE NormTrace
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE CalculateDegen(Norm,Degen,Occpan)
     REAL(DOUBLE)                 :: Norm,Degen,Occpan
!
     Degen     = ((Norm-TrP2)**3)/((TrP2-TrP3)*(Norm-Two*TrP2+TrP3))
     Occpan    = Two*(TrP2-TrP3)/(Norm-TrP2)  
!
   END SUBROUTINE CalculateDegen
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE CalculateGap(F,P,Norm,lumo_occ,Gap)
     TYPE(BCSR)                   :: F,P
     REAL(DOUBLE)                 :: Norm,Gap,lumo_occ
     REAL(DOUBLE)                 :: TrFP,TrFP2,TrFP3
!
     CALL Multiply(P,P,P2)
     CALL Multiply(P,P2,P3)
     TrFP  = Trace(F,P)
     TrFP2 = Trace(F,P2)
     TrFP3 = Trace(F,P3)    
! 
     lumo_occ = Half + Half*Sqrt(ABS(One-Two*Norm+Two*TrP2)) 
!
     Gap      = -(TrFP - Three*TrFP2 + Two*TrFP2)/(lumo_occ*(lumo_occ-One)*(Two*lumo_occ-One))
!
   END SUBROUTINE CalculateGap
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE ShrtOut(ProgName,Iter,Energy,PNon0,ErrorE)
     REAL(DOUBLE)                   :: Energy,ErrorE
     INTEGER                        :: PNon0,Iter
     CHARACTER(LEN=*)               :: ProgName
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!
     IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        Mssg=TRIM(ProgName)//'   '//TRIM(IntToChar(Iter))          &
             //' Tr(P.F) = '//TRIM(DblToMedmChar(Energy))          &
             //', %Non0  = '//TRIM(IntToChar(PNon0))               &
             //', Error(E) = '//TRIM(DblToShrtChar(ErrorE))         
        WRITE(*,*)   TRIM(Mssg)
        CALL OpenASCII(OutFile,Out)
        WRITE(Out,*) TRIM(Mssg)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
!
   END SUBROUTINE ShrtOut
   
   FUNCTION CnvrgChck(Prog,NPur,MM,F,P,POld)
     LOGICAL              :: CnvrgChck
     TYPE(BCSR)           :: F,P,POld
     REAL(DOUBLE)         :: Energy,AbsErrP,AveErrP,MedErrP,  &
                             AbsErrE,RelErrE,AveErrE,MaxCommErr
     REAL(DOUBLE),SAVE    :: OldE,OldAEP
     INTEGER              :: MM,PNon0,NPur
     CHARACTER(LEN=*)     :: Prog
     CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
!---------------------------------------------------------------------
     IF(NPur==0)THEN
        OldE=BIG_DBL
        OldAEP=BIG_DBL
     ENDIF
     PNon0   = 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
     ! Density matrix errors
     CALL Multiply(Pold,-One)
     CALL Add(Pold,P,PTmp1)
     AbsErrP=ABS(Max(PTmp1)+1.D-20)
     AveErrP=OneNorm(PTmp1)/DBLE(PNon0)
     MedErrP=TwoNorm(PTmp1)/TwoNorm(P)
     ! Energy errors
#ifdef PARALLEL
     CALL Multiply(P,F,PTmp1)
     Energy=Trace(PTmp1)    
#else
     Energy=Trace(P,F)
#endif     
     AbsErrE=ABS(OldE-Energy)
     RelErrE=AbsErrE/ABS(Energy)
     ! Convergence check
     CnvrgChck=.FALSE.
     ! Absolute convergence test
     IF(RelErrE<Thresholds%ETol*1D-2.AND. &
        AbsErrP<Thresholds%DTol*1D-1)CnvrgChck=.TRUE.           
     ! Test in the asymptotic regime for stall out
     IF(RelErrE<Thresholds%ETol)THEN
        ! Check for increasing /P
        IF(AbsErrP>OldAEP)CnvrgChck=.TRUE.
        ! Check for an increasing energy
        IF(Energy>OldE)CnvrgChck=.TRUE.
     ENDIF
     ! Updtate previous cycle values
     OldE=Energy
     OldAEP=AbsErrP
     ! Print convergence stats
     Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))      &
          //'dE='//TRIM(DblToShrtChar(RelErrE))                 &
          //', dP='//TRIM(DblToShrtChar(AbsErrP))                &
          //', %Non0='//TRIM(IntToChar(PNon0))              
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(*,*)TRIM(Mssg)
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
#ifdef PARALLEL
     ENDIF
#endif
     IF(.NOT.CnvrgChck)RETURN
     CALL Commute(F,P,PTmp1)      
     MaxCommErr=Max(PTmp1)
     ! Print summary stats
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           Mssg=ProcessName(Prog)//TRIM(IntToChar(NPur))//' purification steps, ' &
                //TRIM(IntToChar(MM))//' matrix multiplies, %Non0s = '//TRIM(IntToChar(PNon0))              
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog)//'Tr{FP}='//TRIM(DblToMedmChar(Energy))//', '                    &
                //'dE='//TRIM(DblToShrtChar(AbsErrP))//', '                      &
                //'dP='//TRIM(DblToShrtChar(AbsErrP))//', '                 &
                //'[F,P]='//TRIM(DblToShrtChar(MaxCommErr))
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
#ifdef PARALLEL
     ENDIF
#endif

   END FUNCTION CnvrgChck
!----------------------------------------------------------------------------
   SUBROUTINE MednOut(ProgName,Iter,Energy,PNon0,ErrorE,ErrorP)
     REAL(DOUBLE)                   :: Energy,ErrorE,ErrorP
     INTEGER                        :: PNon0,Iter
     CHARACTER(LEN=*)               :: ProgName
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!


     IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        Mssg=TRIM(ProgName)//'   '//TRIM(IntToChar(Iter))          &
             //' Tr(P.F) = '//TRIM(DblToMedmChar(Energy))          &
             //', %Non0  = '//TRIM(IntToChar(PNon0))               &
             //', Error(E) = '//TRIM(DblToShrtChar(ErrorE))        &
             //', Error(P) = '//TRIM(DblToShrtChar(ErrorP))    
        WRITE(*,*)   TRIM(Mssg)
        CALL OpenASCII(OutFile,Out)
        WRITE(Out,*) TRIM(Mssg)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
   END SUBROUTINE MednOut
!----------------------------------------------------------------------------
   SUBROUTINE LongOut(ProgName,Iter,Energy,PNon0,ErrorE,ErrorP,ErrorFP)
     REAL(DOUBLE)                   :: Energy,ErrorE,ErrorP,ErrorFP
     INTEGER                        :: PNon0,Iter
     CHARACTER(LEN=*)               :: ProgName
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!
     IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        Mssg=TRIM(ProgName)//'   '//TRIM(IntToChar(Iter))          &
             //' Tr(P.F) = '//TRIM(DblToMedmChar(Energy))          &
             //', %Non0  = '//TRIM(IntToChar(PNon0))               &
             //', Error(E) = '//TRIM(DblToShrtChar(ErrorE))        &
             //', Error(P) = '//TRIM(DblToShrtChar(ErrorP))        &         
             //', Error(FP) = '//TRIM(DblToShrtChar(ErrorFP))    
        WRITE(*,*)   TRIM(Mssg)
        CALL OpenASCII(OutFile,Out)
        WRITE(Out,*) TRIM(Mssg)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
   END SUBROUTINE LongOut
!----------------------------------------------------------------------------
   SUBROUTINE FinalOut(ProgName,Energy,ErrorE,ErrorN,ErrorP,ErrorFP,PNon0,Count)
     REAL(DOUBLE)                   :: Energy,ErrorN,ErrorE,ErrorP,ErrorFP
     INTEGER                        :: PNon0,Iter,Count
     CHARACTER(LEN=*)               :: ProgName
     CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!
     IF(PrintFlags%Key .NE. DEBUG_NONE) THEN
        CALL OpenASCII(OutFile,Out)
        Mssg = TRIM(ProgName) // '               : Energy             = '//TRIM(DblToMedmChar(Energy))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        Mssg = TRIM(ProgName) // '               : Error(E(n+1)-E(n)) = '//TRIM(DblToMedmChar(ErrorE))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        Mssg = TRIM(ProgName) // '               : Error(Tr[P]-N)     = '//TRIM(DblToMedmChar(ErrorN))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        Mssg = TRIM(ProgName) // '               : Error(P(n)-P(n+1)) = '//TRIM(DblToMedmChar(ErrorP))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        Mssg = TRIM(ProgName) // '               : Error([F,P])       = '//TRIM(DblToMedmChar(ErrorFP))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        Mssg = TRIM(ProgName) // '               : Sparcity of Matrix = '//TRIM(IntToChar(PNon0))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        Mssg = TRIM(ProgName) // '               : No. of Matrix Mult = '//TRIM(IntToChar(Count))
        WRITE(*,*)   TRIM(Mssg)
        WRITE(Out,*) TRIM(Mssg)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
!
   END SUBROUTINE FinalOut
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE SP2(P,Norm,Count,DoThresh)
      TYPE(BCSR)                     :: P
      REAL(DOUBLE)                   :: Norm,CR
      INTEGER                        :: Count
      LOGICAL                        :: DoThresh
!
      Count = Count+1
      CALL Multiply(P,P,P2)             ! The only multiplication is a square
      TrP  = Trace(P)
      CR   = TrP - Norm                 ! CR = Occupation error criteria
      IF (CR >  0) THEN                 ! Too many states
         CALL SetEq(Ptmp1,P2)           ! P = P^2
      ELSE                              ! Too few states 
         CALL Multiply(P ,Two) 
         CALL Multiply(P2,-One)
         CALL Add(P,P2,Ptmp1)         ! P = 2P-P^2
      ENDIF
!
      IF(DoThresh) THEN 
         CALL Filter(P,Ptmp1)
      ELSE
         CALL  SetEq(P,Ptmp1)
      ENDIF
!
    END SUBROUTINE SP2
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE SP4(P,Norm,Count,DoThresh)
      TYPE(BCSR)                     :: P
      REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Thresh_old
      INTEGER                        :: Count
      LOGICAL                        :: DoThresh
      REAL(DOUBLE)                   :: EPS = 1.D-12
!
      Count = Count+1
      CALL Multiply(P,P,P2)
!
      TrP  = Trace(P)
      TrP2 = Trace(P2)
      TrP3 = Trace(P,P2)
      TrP4 = Trace(P2,P2)
      Gt   = (Norm-Four*TrP3+3*TrP4)
      Gn   = (TrP2-Two*TrP3+TrP4)
      CR   = TrP - Norm                      ! Trace correction measure
!
      IF (ABS(Gn).LT.EPS) THEN               ! Close to idempotency
         G = Three                           ! To avoid numerical errors
      ELSE
         G = Gt/Gn                           ! Boundary measure
      ENDIF
!
      IF ( (G > Zero) .AND. (G < Six) ) THEN ! Check the bounds
         IF (CR  > Zero) THEN                ! Too many states
           Count = Count+1
           IF(DoThresh) THEN
              Thresholds%Trix = Thresholds%Trix*ThFAC1
              CALL Filter(Ptmp2,P2)
              Thresholds%Trix = Thresholds%Trix*ThFAC2
           ELSE
              CALL  SetEq(Ptmp2,P2)
           ENDIF
           CALL Multiply(P , 4.d0)
           CALL Multiply(P2,-3.d0)
           CALL Add(P,P2,Ptmp1)
           IF(DoThresh) THEN
              CALL Filter(P,Ptmp1)
           ELSE
              CALL  SetEq(P,Ptmp1)
           ENDIF
           CALL Multiply(Ptmp2,P,Ptmp1)      ! P = P^2*(4P-3P^2)
        ELSE                                 ! Too few states
           Count = Count+1
           IF(DoThresh) THEN
              Thresholds%Trix = Thresholds%Trix*ThFAC1
              CALL Filter(Ptmp2,P2)
              Thresholds%Trix = Thresholds%Trix*ThFAC2
           ELSE
              CALL  SetEq(Ptmp2,P2)
           ENDIF
           CALL Multiply(P ,-8.d0)
           CALL Multiply(P2, 3.d0)
           CALL Add(P,P2,Ptmp1)
           CALL Add(Ptmp1,Six)
           IF(DoThresh) THEN
              CALL Filter(P,Ptmp1)
           ELSE
              CALL  SetEq(P,Ptmp1)
           ENDIF
           CALL Multiply(Ptmp2,P,Ptmp1)      ! P = P^2*(6I-8P+3P^2)
        ENDIF
     ELSE 
        IF (G < 0) THEN                      ! Too many states
           CALL SetEq(Ptmp1,P2)              ! P = P^2
        ELSE                                 ! Too few states 
           CALL Multiply(P ,Two) 
           CALL Multiply(P2,-One)
           CALL Add(P ,P2,Ptmp1)            ! P = 2P-P^2
        ENDIF
     ENDIF
!
     IF(DoThresh) THEN 
        CALL Filter(P,Ptmp1)
     ELSE
        CALL  SetEq(P,Ptmp1)
     ENDIF
!
   END SUBROUTINE SP4
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE NT4(P,Norm,Count,DoThresh)
     TYPE(BCSR)                     :: P
     REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Coeff
     INTEGER                        :: Count
     LOGICAL                        :: DoThresh   
     REAL(DOUBLE)                   :: EPS = 1.D-12
!
     Count = Count+1
     CALL Multiply(P,P,P2)    
     TrP  = Trace(P)
     TrP2 = Trace(P2)
     TrP3 = Trace(P,P2)
     TrP4 = Trace(P2,P2)
     Gt   = (Norm-Four*TrP3+3*TrP4)
     Gn   = (TrP2-Two*TrP3+TrP4)
     CR   = TrP - Norm     
!
     IF (ABS(Gn).LT.EPS) THEN               ! Close to idempotency
        G = Three                           ! To avoid numerical errors
     ELSE
        G = Gt/Gn                           ! Boundary measure
     ENDIF
!
     IF ( (G > Zero) .AND. (G < Six) ) THEN ! Check the bounds
        Count = Count+1
        IF(DoThresh) THEN
!          Thresholds%Trix = Thresholds%Trix*ThFAC1
           CALL Filter(Ptmp2,P2)
!           Thresholds%Trix = Thresholds%Trix*ThFAC2
        ELSE
           CALL  SetEq(Ptmp2,P2)
        ENDIF
        Coeff = Four-Two*G
        CALL Multiply(P ,Coeff)
        Coeff = G-Three
        CALL Multiply(P2,Coeff)
        CALL Add(P,P2,Ptmp1)
        CALL Add(Ptmp1,G)
        IF(DoThresh) THEN
!           Thresholds%Trix = Thresholds%Trix*ThFAC1
           CALL Filter(P,Ptmp1)
!           Thresholds%Trix = Thresholds%Trix*ThFAC2
        ELSE
           CALL  SetEq(P,Ptmp1)
        ENDIF
        CALL Multiply(Ptmp2,P,Ptmp1)         ! P^2*(G+(4-2*G)*P+(G-3)*P^2)
     ELSE
        IF (G < 0) THEN                      ! Too many states
           CALL SetEq(Ptmp1,P2)              ! P = P^2
        ELSE                                 ! Too few states 
           CALL Multiply(P ,Two) 
           CALL Multiply(P2,-One)
           CALL Add(P,P2,Ptmp1)             ! P = 2P-P^2
        ENDIF
     ENDIF
!
     IF(DoThresh) THEN 
        CALL Filter(P,Ptmp1)
     ELSE
        CALL  SetEq(P,Ptmp1)
     ENDIF
!
   END SUBROUTINE NT4
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE PM(P,Norm,Count,DoThresh)
     TYPE(BCSR)                     :: P
     REAL(DOUBLE)                   :: CR,Norm,Coeff
     INTEGER                        :: Count
     LOGICAL                        :: DoThresh
     REAL(DOUBLE)                   :: EPS = 1.D-8
!
     Count = Count+2
     CALL Multiply(P,P,P2)                        ! The only multiplication is a square
!
     IF(DoThresh) THEN 
        Thresholds%Trix = Thresholds%Trix*ThFAC1
        CALL Filter(Ptmp1,P2)
        Thresholds%Trix = Thresholds%Trix*ThFAC2
        CALL Multiply(P,Ptmp1,P3)                 ! The only multiplication is a square
     ELSE
        CALL Multiply(P,P2,P3)     
     ENDIF
!
     TrP  = Trace(P)
     TrP2 = Trace(P2)
     TrP3 = Trace(P3)
!
     IF (ABS(TrP-TrP2) < EPS) THEN                ! Close to idempotency
       CALL Multiply(P2, Three)
       CALL Multiply(P3,-Two)
       CALL Add(P2,P3,Ptmp1)                      ! T = 3P^2-2P^2 (McWeeny close to convergence)
     ELSE
        CR = (TrP2 - TrP3)/(TrP-TrP2)
        IF (CR <= Half) THEN
           Coeff = (One-Two*CR)/(One-CR)
           CALL Multiply(P,Coeff)
           Coeff = (One+CR)/(One-CR)
           CALL Multiply(P2,Coeff)
           Coeff = -One/(One-CR)
           CALL Multiply(P3,Coeff)
           CALL Add(P,P2,Ptmp2)
           CALL Add(Ptmp2,P3,Ptmp1)              ! T = [(1-2CR)*P+(1+CR)*P^2-P^3]/(1-CR)
        ELSE
           Coeff = (One+CR)/CR
           CALL Multiply(P2,Coeff)
           Coeff = -One/CR
           CALL Multiply(P3,Coeff)
           CALL Add(P2,P3,Ptmp1)                 ! T = [(1+CR)*P^2-P^3]/CR
        ENDIF
     ENDIF
!
     IF(DoThresh) THEN 
        CALL Filter(P,Ptmp1)
     ELSE
        CALL  SetEq(P,Ptmp1)
     ENDIF
!
   END SUBROUTINE PM
!
END MODULE DenMatMethods
