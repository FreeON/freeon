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
      TYPE(BCSR)                     :: F,P,FP,PF,T
!
      CALL Multiply(F,P,FP)
      CALL Multiply(P,F,PF)
      CALL Multiply(PF,-One)
      CALL Add(FP,PF,T)
!
      CALL Delete(FP)
      CALL Delete(PF)
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
   SUBROUTINE NormTrace(Pin,Norm,Order)
     TYPE(BCSR)                     :: Pin,P2,Ptmp1
     REAL(DOUBLE)                   :: Norm,CR,C1,C2
     INTEGER                        :: Order
! 
     IF(Order==1) THEN
        TrP  = Trace(Pin)
        CALL Multiply(Pin,Pin,Ptmp1)
        CALL Filter(P2,Ptmp1)
        TrP2 = Trace(P2)
        CR   = TrP - TrP2
        IF(CR==Zero)RETURN
        C1   = (Norm - TrP2)/CR
        C2   = (TrP- Norm)/CR
        CALL Multiply(Pin ,C1)
        CALL Multiply(P2  ,C2)
        CALL Add(Pin,P2,Ptmp1)
        CALL Filter(Pin,Ptmp1)
     ELSE
        CALL MondoHalt(99,'Wrong Order in NormTrace')
     ENDIF
     CALL Delete(P2)
     CALL Delete(Ptmp1)
!
   END SUBROUTINE NormTrace
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE CalculateDegen(Norm,Degen,Occpan,lumo_occ)
     REAL(DOUBLE)                 :: Norm,Degen,Occpan,lumo_occ
!
     Degen     = ((Norm-TrP2)**3)/((TrP2-TrP3)*(Norm-Two*TrP2+TrP3))
     Occpan    = Two*(TrP2-TrP3)/(Norm-TrP2)  
     lumo_occ  = Half + Half*Sqrt(ABS(One-Two*Norm+Two*TrP2)) 
!
   END SUBROUTINE CalculateDegen
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE CalculateGap(F,P,Norm,lumo_occ,Gap)
     TYPE(BCSR)                   :: F,P,P2,P3
     REAL(DOUBLE)                 :: Norm,Gap,lumo_occ
     REAL(DOUBLE)                 :: TrFP,TrFP2,TrFP3
!
     CALL Multiply(P,P,P2)
     CALL Multiply(P,P2,P3)
     TrFP  = Trace(F,P)
     TrFP2 = Trace(F,P2)
     TrFP3 = Trace(F,P3)    
! 
     Gap       = -(TrFP - Three*TrFP2 + Two*TrFP2)/(lumo_occ*(lumo_occ-One)*(Two*lumo_occ-One))
!
     CALL Delete(P2)
     CALL Delete(P3)
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
    SUBROUTINE SP2(Pin,Pout,Norm,Count,DoThresh)
      TYPE(BCSR)                     :: Pin,Pout,P2,Ptmp1
      REAL(DOUBLE)                   :: Norm,CR
      INTEGER                        :: Count
      LOGICAL                        :: DoThresh
!
      Count = Count+1
      CALL Multiply(Pin,Pin,P2)         ! The only multiplication is a square
      TrP  = Trace(Pin)
      CR   = TrP - Norm                 ! CR = Occupation error criteria
      IF (CR >  0) THEN                 ! Too many states
         CALL SetEq(Ptmp1,P2)           ! P = P^2
      ELSE                              ! Too few states 
         CALL Multiply(Pin,Two) 
         CALL Multiply(P2,-One)
         CALL Add(Pin,P2,Ptmp1)         ! P = 2P-P^2
         CALL Multiply(Pin,Half) 
      ENDIF
!
      IF(DoThresh) THEN 
         CALL Filter(Pout,Ptmp1)
      ELSE
         CALL  SetEq(Pout,Ptmp1)
      ENDIF
      CALL Delete(P2)
      CALL Delete(Ptmp1)
!
    END SUBROUTINE SP2
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE SP4(Pin,Pout,Norm,Count,DoThresh)
      TYPE(BCSR)                     :: Pin,Pout,P2,Ptmp1,Ptmp2
      REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Thresh_old
      INTEGER                        :: Count
      LOGICAL                        :: DoThresh
      REAL(DOUBLE)                   :: EPS = 1.D-12
!
      Count = Count+1
      CALL Multiply(Pin,Pin,P2)   
!
      IF(DoThresh) THEN 
         Thresholds%Trix = 0.1D0*Thresholds%Trix
         CALL Filter(Ptmp1,P2)
         CALL  SetEq(P2,Ptmp1)
         Thresholds%Trix = 10.D0*Thresholds%Trix
      ENDIF
!
      TrP  = Trace(Pin)
      TrP2 = Trace(P2)
      TrP3 = Trace(Pin,P2)
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
         IF (CR  > Zero) THEN              ! Too many states
           Count = Count+1
           CALL SetEq(Ptmp1,Pin)
           CALL SetEq(Ptmp2,P2)
           CALL Multiply(Ptmp1, 4.d0)
           CALL Multiply(Ptmp2,-3.d0)
           CALL Add(Ptmp1,Ptmp2,Pout)
           CALL Multiply(P2,Pout,Ptmp1)      ! P = P^2*(4P-3P^2)
        ELSE                                 ! Too few states
           Count = Count+1
           CALL SetEq(Ptmp1,Pin)
           CALL SetEq(Ptmp2,P2)
           CALL Multiply(Ptmp1,-8.d0)
           CALL Multiply(Ptmp2, 3.d0)
           CALL Add(Ptmp1,Ptmp2,Pout)
           CALL Add(Pout,Six)
           CALL Multiply(P2, Pout,Ptmp1)     ! P = P^2*(6I-8P+3P^2)
        ENDIF
     ELSE 
        IF (G < 0) THEN                      ! Too many states
           CALL SetEq(Ptmp1,P2)              ! P = P^2
        ELSE                                 ! Too few states 
           CALL Multiply(Pin,Two) 
           CALL Multiply(P2,-One)
           CALL Add(Pin,P2,Ptmp1)            ! P = 2P-P^2
           CALL Multiply(Pin,Half) 
        ENDIF
     ENDIF
!
     IF(DoThresh) THEN 
        CALL Filter(Pout,Ptmp1)
     ELSE
        CALL  SetEq(Pout,Ptmp1)
     ENDIF

     CALL Delete(P2)
     IF(AllocQ(Ptmp1%Alloc))&
     CALL Delete(Ptmp1)
     IF(AllocQ(Ptmp2%Alloc))&
     CALL Delete(Ptmp2)
   END SUBROUTINE SP4
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE NT4(Pin,Pout,Norm,Count,DoThresh)
     TYPE(BCSR)                     :: Pin,Pout,P2,Ptmp1,Ptmp2
     REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Coeff
     INTEGER                        :: Count
     LOGICAL                        :: DoThresh   
     REAL(DOUBLE)                   :: EPS = 1.D-12
!
     Count = Count+1
     CALL Multiply(Pin,Pin,P2)    
!
     TrP  = Trace(Pin)
     TrP2 = Trace(P2)
     TrP3 = Trace(Pin,P2)
     TrP4 = Trace(P2,P2)
!
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
     IF (G > 6) THEN                        ! Out of bound, increase eigenvalues
        CALL Multiply(Pin,Two)
        CALL Multiply(P2,-One)
        CALL Add(Pin,P2,Ptmp1)              ! P = 2P-P^2
        CALL Multiply(Pin,Half)
     ELSEIF (G < Zero) THEN                 ! Out of bound, decrease eigenvalues
        CALL SetEq(Ptmp1,P2)                ! P = P^2
     ELSE     
        Count = Count+1
        CALL SetEq(Ptmp1,Pin)
        CALL SetEq(Ptmp2,P2)
        Coeff = G-Three
        CALL Multiply(Ptmp2,Coeff)
        Coeff = Four-Two*G
        CALL Multiply(Ptmp1,Coeff)
        CALL Add(Ptmp1,Ptmp2,Pout)
        CALL Add(Pout,G)
!
        IF(DoThresh) THEN 
           Thresholds%Trix = 0.1D0*Thresholds%Trix
           CALL Filter(Ptmp1,P2)
           CALL SetEq(P2,Ptmp1)
           CALL Filter(Ptmp1,Pout)
           CALL SetEq(Pout,Ptmp1)
           Thresholds%Trix = 10.D0*Thresholds%Trix
        ENDIF
!
        CALL Multiply(P2,Pout,Ptmp1) 
     ENDIF
!
     IF(DoThresh) THEN 
        CALL Filter(Pout,Ptmp1)
     ELSE
        CALL  SetEq(Pout,Ptmp1)
     ENDIF

     CALL Delete(P2)
     IF(AllocQ(Ptmp1%Alloc))&
     CALL Delete(Ptmp1)
     IF(AllocQ(Ptmp2%Alloc))&
     CALL Delete(Ptmp2)
!
   END SUBROUTINE NT4
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE PM(Pin,Pout,Norm,Count,DoThresh)
     TYPE(BCSR)                     :: Pin,Pout,P2,P3,Ptmp1,Ptmp2
     REAL(DOUBLE)                   :: CR,Norm,Coeff
     INTEGER                        :: Count
     LOGICAL                        :: DoThresh
     REAL(DOUBLE)                   :: EPS = 1.D-8
!
     Count = Count+2
     CALL Multiply(Pin,Pin,P2)                   ! The only multiplication is a square
!
     IF(DoThresh) THEN 
        Thresholds%Trix = 0.1D0*Thresholds%Trix
        CALL Filter(Ptmp1,P2)
        CALL Multiply(Pin,Ptmp1,P3)              ! The only multiplication is a square
        Thresholds%Trix = 10.D0*Thresholds%Trix
     ELSE
        CALL Multiply(Pin,P2,P3)     
     ENDIF
!
     TrP  = Trace(Pin)
     TrP2 = Trace(P2)
     TrP3 = Trace(P3)
!
     IF (ABS(TrP-TrP2) < EPS) THEN        ! Close to idempotency
       CALL Multiply(P2, Three)
       CALL Multiply(P3,-Two)
       CALL Add(P2,P3,Ptmp1)              ! T = 3P^2-2P^2 (McWeeny close to convergence)
     ELSE
        CR = (TrP2 - TrP3)/(TrP-TrP2)
        IF (CR <= Half) THEN
           CALL SetEq(Ptmp1,Pin)
           Coeff = (One-Two*CR)/(One-CR)
           CALL Multiply(Ptmp1,Coeff)
           Coeff = (One+CR)/(One-CR)
           CALL Multiply(P2,Coeff)
           Coeff = -One/(One-CR)
           CALL Multiply(P3,Coeff)
           CALL Add(Ptmp1,P2,Ptmp2)
           CALL Add(Ptmp2,P3,Ptmp1)       ! T = [(1-2CR)*P+(1+CR)*P^2-P^3]/(1-CR)
        ELSE
           Coeff = (One+CR)/CR
           CALL Multiply(P2,Coeff)
           Coeff = -One/CR
           CALL Multiply(P3,Coeff)
           CALL Add(P2,P3,Ptmp1)          ! T = [(1+CR)*P^2-P^3]/CR
        ENDIF
     ENDIF
!
     IF(DoThresh) THEN 
        CALL Filter(Pout,Ptmp1)
     ELSE
        CALL  SetEq(Pout,Ptmp1)
     ENDIF

     CALL Delete(P2)
     CALL Delete(P3)
     CALL Delete(Ptmp1)
     CALL Delete(Ptmp2)
!
   END SUBROUTINE PM
!
END MODULE DenMatMethods
