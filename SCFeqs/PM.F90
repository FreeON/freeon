!----------------------------------------------------------------------
!                         May 15th 2002
!                        Anders Niklasson
!           Palser Manoulopulos canonical purification scheme
!----------------------------------------------------------------------
PROGRAM PM
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
  TYPE(BCSR)                     :: F,P,P2,P3,Pold,T,T2,Z
  TYPE(ARGMT)                    :: Args
!
  REAL(DOUBLE)                   :: EigH,EigL,E_old,E_new, Error_Energy,  &
                                    Tr,Tol_Energy,Diag,Tmp_min,Tmp_max,F_min,F_max, &
                                    Sum,Energy,DF,Coeff,CR,EPS,Mu,Lmbd,Lmbd1, &
                                    Ne,TrP,TrP2,TrP3,Nbsf,ErrorP,Fill,MinErrorP
  INTEGER                        :: I,J,K,MM,M,N,Nr_Max_It,R,C,col,Blk,Check,PNon0
  LOGICAL                        :: Present,Converged
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=2),PARAMETER     :: Prog='PM'
!------------------------------------------------------------------
! Start
  CALL StartUp(Args,Prog)
!-------------------------------------------------------------------------
  CALL New(F)
!-----------------------------------------------------------------------------
! Get The Fock Matrix
!-----------------------------------------------------------------------------
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF
!=============================================================================
!  PM-Solver 
!=============================================================================
!
!  Set up the Density Matrix from the Fock Matrix
!
!=============================================================================
  Tol_Energy = Thresholds%ETol*1.D-2      ! Convergence criterea
  Ne         = Half*DBLE(NEl)             ! 0.5 because of spin degeneracy
  Fill       = Ne/DBLE(NBasF) 
!
  CALL New(P)
  CALL New(P2)
  CALL New(P3)
  CALL New(T)
  CALL New(T2)
!-----------------------------------------------------
! Estimate spectral bounds 
!-----------------------------------------------------
  CALL SpectralBounds(F,F_min,F_max)
!-----------------------------------------------------
! Set up the starting guess P
!-----------------------------------------------------
!   
  CALL SetEq(T,F)
  CALL Multiply(T,Zero)
  CALL Add(T,One)
! 
  Nbsf = DBLE(NBasF)
  Mu = Trace(F)/DBLE(Nbsf)
  Lmbd1 = Ne/(F_max-Mu)
  Lmbd = (Nbsf-Ne)/(Mu-F_min)
  IF (Lmbd1.LT.Lmbd) THEN
     Lmbd = Lmbd1
  ENDIF
  CALL SetEq(P,F)
  CALL Add(P,-Mu)
  Coeff = -Lmbd/Nbsf
  CALL Multiply(P,Coeff) 
  Coeff = Ne/Nbsf 
  CALL Add(P,Coeff)    ! P = (Lmbd/Nbsf)*(I*mu-F) + I*Ne/Nbsf
!==============================================================
! Main Loop: Iterate until convergence
!==============================================================
  EPS       = 1.D-8
  MM        = 0                  ! Nr of matrix multiplies
  E_old     = Zero
  MinErrorP = 1.0D10
  Nr_Max_It = 100                ! Maximal number of allowed iterations
  Converged = .FALSE.
!
  DO I = 1,Nr_Max_It
     CALL SetEq(Pold,P)
     CALL Multiply(P,P,P2)       ! The only multiplication is a square
     CALL Multiply(P2,P,P3)      ! The only multiplication is a square
     MM   = MM + 2                 ! Two multiplications
     TrP  = Trace(P)
     TrP2 = Trace(P2)
     TrP3 = Trace(P3)
     IF (ABS(TrP-TrP2).LT.EPS) THEN  ! Cloese to idempotency
       Coeff = Three
       CALL Multiply(P2,Coeff)
       Coeff = Two
       CALL Multiply(P3,-Coeff)
       CALL ADD(P2,P3,T)          ! T = 3P^2-2P^2 (McWeeny close to convergence)
     ELSE
        CR = (TrP2 - TrP3)/(TrP-TrP2)
        IF (CR.LE.Half) THEN
           Coeff = -One+Two*CR
           CALL Multiply(P,Coeff)
           Coeff = -One-CR
           CALL Multiply(P2,Coeff)
           CALL Add(P,P2,T2)
           CALL Add(T2,P3,T)
           Coeff = -One/(One-CR)
           CALL Multiply(T,Coeff)  ! T = [(1-2CR)*P+(1+CR)*P^2-P^3]/(1-CR)
        ELSE
           Coeff = One+CR
           CALL Multiply(P2,-Coeff)
           CALL Add(P2,P3,T)
           Coeff = One/CR
           CALL Multiply(T,-Coeff) ! T = [(1+CR)*P^2-P^3]/CR
        ENDIF
     ENDIF
     CALL Filter(P,T)
!
     TrP            = Trace(P)
     Energy         = Trace(P,F)
     Error_Energy   = ABS(Energy - E_old)
     E_old          = Energy
     PNon0          = 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF) 
!
     IF(I .GT. 20) THEN
        CALL Multiply(Pold,-One)
        CALL Add(Pold,P,T) 
        ErrorP    = MAX(T)
        MinErrorP = MIN(MinErrorP,ErrorP)
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           Mssg=TRIM(Prog)//'   '//TRIM(IntToChar(I))                 &
                //' Tr(P.F) = '//TRIM(DblToMedmChar(Energy))          &
                //', %Non0 = '//TRIM(IntToChar(PNon0))                &
                //', ErrorP = '//TRIM(DblToShrtChar(ErrorP))         
           WRITE(*,*)   TRIM(Mssg)
           CALL OpenASCII(OutFile,Out)
           WRITE(Out,*) TRIM(Mssg)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
        IF(ErrorP/MinErrorP .GT. 1.5D0) THEN
           IF(Error_Energy .LT. 1.D4*Tol_Energy) THEN
               Converged=.TRUE.
           ENDIF
           IF(ErrorP/MinErrorP .GT. 4.0D0) THEN
              Converged=.TRUE.
           ENDIF
        ENDIF
        IF(Error_Energy .LT. Tol_Energy) Converged=.TRUE.
     ELSE        
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           Mssg=TRIM(Prog)//'   '//TRIM(IntToChar(I))                 &
                //' Tr(P.F) = '//TRIM(DblToMedmChar(Energy))          &
                //', %Non0 = '//TRIM(IntToChar(PNon0))                    
           WRITE(*,*)   TRIM(Mssg)
           CALL OpenASCII(OutFile,Out)
           WRITE(Out,*) TRIM(Mssg)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
        IF(Error_Energy .LT. Tol_Energy) Converged=.TRUE.
     ENDIF
!
     IF(Converged) EXIT
  ENDDO 
! 
!--------------------------------------------------------------------------
! P now contains the Fock density matrix, Energy is
! the sum of the occupied eigenvalues, and MM is
! the number of necessary matrix multiplications.
!--------------------------------------------------------------------------
! Write Out Statisitcs
!--------------------------------------------------------------------------
  TrP = Trace(P)
  CALL Multiply(P,F,P2)
  CALL Multiply(F,P,T)
  CALL Multiply(T,-One)
  CALL Add(T,P2,T)
!        
  CALL OpenASCII(OutFile,Out)
  Mssg = TRIM(Prog) // '               : |Trace(P)-NEl|= '//TRIM(DblToMedmChar(ABS(Two*TrP-NEl)))
  WRITE(*,*)   TRIM(Mssg)
  WRITE(Out,*) TRIM(Mssg)
  Mssg = TRIM(Prog) // '               : Tol_Energy    = '//TRIM(DblToMedmChar(Tol_Energy))
  WRITE(*,*)   TRIM(Mssg)
  WRITE(Out,*) TRIM(Mssg)
  Mssg = TRIM(Prog) // '               : |E_n-E_{n-1}| = '//TRIM(DblToMedmChar(Error_Energy))
  WRITE(*,*)   TRIM(Mssg)
  WRITE(Out,*) TRIM(Mssg)
  Mssg = TRIM(Prog) // '               : Max|PF-FP|    = '//TRIM(DblToMedmChar( MAX(T)/MAX(P2) ))
  WRITE(*,*)   TRIM(Mssg)
  WRITE(Out,*) TRIM(Mssg)
  Mssg = TRIM(Prog) // '               : Nr Mtrx Mlt   = '//TRIM(IntToChar(MM))
  WRITE(*,*)   TRIM(Mssg)
  WRITE(Out,*) TRIM(Mssg)
  Mssg = TRIM(Prog) // '               : Filling       = '//TRIM(DblToMedmChar(Fill))
  WRITE(*,*)   TRIM(Mssg)
  WRITE(Out,*) TRIM(Mssg)
  CLOSE(UNIT=Out,STATUS='KEEP')
!
  CALL Delete(F)
  CALL Delete(P2)
  CALL Delete(P3)
  CALL Delete(Pold)
!=============================================================================
!
!  TRANSFORMATION TO AN AO REPRESENTATION AND IO
!
!=============================================================================
!-----------------------------------------------------------------------------
!  IO for the orthogonal P
!
   CALL Put(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)
   CALL Put(P,TrixFile('OrthoD',Args,1))
   CALL PChkSum(P,'OrthoP['//TRIM(NxtCycl)//']',Prog)
   CALL PPrint( P,'OrthoP['//TRIM(NxtCycl)//']')
   CALL Plot(   P,'OrthoP_'//TRIM(NxtCycl))
!-----------------------------------------------------------------------------
!  Convert to AO representation
!
   INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
   IF(Present)THEN
      CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Z,P,T)
      CALL Multiply(T,Z,P)
   ELSE
      CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Z,P,T)
      CALL Get(Z,TrixFile('ZT',Args))
      CALL Multiply(T,Z,P)
   ENDIF
   CALL Filter(T,P)
!----------------------------------------------------------------------------0-
!  IO for the non-orthogonal P
!
   CALL Put(T,TrixFile('D',Args,1))
   CALL Put(Zero,'homolumogap')
   CALL PChkSum(T,'P['//TRIM(NxtCycl)//']',Prog)
   CALL PPrint(T,'P['//TRIM(NxtCycl)//']')
   CALL Plot(T,'P_'//TRIM(NxtCycl))
!-----------------------------------------------------------------------------
!  Tidy up
!
   CALL Delete(P)
   CALL Delete(T)
   CALL Delete(Z)
   CALL ShutDown(Prog)
!
 END PROGRAM PM
