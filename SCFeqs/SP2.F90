!----------------------------------------------------------------------
!                         May 13th 2002
! Anders M. N. Niklasson: "Expansion Algorithm for the Density Matrix".
! Constructs the density matrix from the Hamiltonian in terms of a
! trace correcting purification expansion with 2nd order purifications.
!----------------------------------------------------------------------
PROGRAM SP2
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
  TYPE(BCSR)                     :: F,P,Pold,P2,T,Z
!-------------------------------------------------------------------------------------
! Note, with Add(A,B,A) and Thresh(P) possible T is not necessary and less memory is
! required! Only P^2 has to be calculated which possibly could be optimized!
!-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
!
  REAL(DOUBLE)                   :: E_old,E_new,Error_Energy,Tol_Energy,Diag, &
                                    Tmp_min,Tmp_max,F_min,F_max,Sum,Energy,DF, &
                                    TrP,TrP2,Coeff,CR,Ne,Slask,ErrorP,MinErrorP,Fill,C1,C2,Norm
  INTEGER                        :: I,J,K,MM,M,N,Nr_Max_It,R,C,col,Blk,Check,PNon0
  LOGICAL                        :: Present,Converged
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='SP2'
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
!=====================================================
!  "N"-Solver of 2nd order                           
!=====================================================  
  Tol_Energy = Thresholds%ETol*1D-2    ! Convergence criterea
  Ne         = Half*DBLE(NEl)          ! 0.5 because of spin degeneracy?
  Fill       = Ne/DBLE(NBasF) 
  CALL New(P)
  CALL New(P2)
  CALL New(Pold)
  CALL New(T)                          ! In principle not necessary for the algorithm
!-----------------------------------------------------
! Estimate spectral bounds 
!-----------------------------------------------------
  CALL SpectralBounds(F,F_min,F_max)
!-----------------------------------------------------------
!  Set up the starting Density Matrix from the Fock Matrix
!-----------------------------------------------------------
  Call SetEq(P,F)
  DF = (F_max - F_min)
  CALL Add(P,-F_max)
  Coeff = -One/DF
  CALL Multiply(P,Coeff)               ! P = (I*F_max-F)/DF
!-----------------------------------------------------------
! Main Loop: Iterate until convergence              
!-----------------------------------------------------------
  MM         = 0                        ! Nr of matrix multiplies
  E_old      = Zero
  MinErrorP  = 1.0D10
  Nr_Max_It  = 200                      ! Maximal number of allowed iterations
  Converged  = .FALSE.
!
  DO I = 1,Nr_Max_It
     CALL SetEq(Pold,P)
     CALL Multiply(P,P,P2)             ! The only multiplication is a square
     MM = MM + 1                       ! One multiplication
     CR = Trace(P) - Ne                ! CR = Occupation error criteria
     IF (CR.GT.0) THEN                 ! Too many states
        CALL SetEq(T,P2)               ! P = P^2
     ELSE                              ! Too few states 
        CALL Multiply(P,Two) 
        CALL Multiply(P2,-One)
        CALL Add(P,P2,T)               ! P = 2P-P^2
     ENDIF
     CALL Filter(P,T)                  ! Thresholding
!
     Energy       = Trace(P,F)
     Error_Energy = ABS(Energy - E_old)
     E_old        = Energy
     PNon0        = 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF) 
!
     IF(I .GT. 20) THEN
        CALL Multiply(Pold,-One)
        CALL Add(Pold,P,T) 
        ErrorP    = MAX(T)
        MinErrorP = MIN(MinErrorP,ErrorP)
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           Mssg=TRIM(Prog)//'   '//TRIM(IntToChar(I))                 &
                //' Tr(P.F) = '//TRIM(DblToMedmChar(Energy))          &
                //', Tr(P)-N = '//TRIM(DblToShrtChar(CR))             &
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
                //', Tr(P)-N = '//TRIM(DblToShrtChar(CR))             &
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
!--------------------------------------------------------------------------
! Normalize Trace
!--------------------------------------------------------------------------
  TrP  = Trace(P)
  CALL Multiply(P,P,P2) 
  TrP2 = Trace(P2)
  Norm = TrP - TrP2
  C1   = (Ne - TrP2)/Norm 
  C2   = (TrP- Ne)/Norm
  CALL Multiply(P ,C1)
  CALL Multiply(P2,C2)
  CALL Add(P,P2,T)
  CALL Filter(P,T)
!--------------------------------------------------------------------------
! Write Out Statisitcs
!--------------------------------------------------------------------------
  TrP     = Trace(P)
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
  CALL Delete(Pold)
!=============================================================================
!  TRANSFORMATION TO AN AO REPRESENTATION AND IO
!=============================================================================
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
   CALL Filter(T,P)     ! Thresholding
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
 END PROGRAM SP2




