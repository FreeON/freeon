!----------------------------------------------------------------------
!                         May 13th 2002
! Anders M. N. Niklasson: "Expansion Algorithm for the Density Matrix".
! Constructs the density matrix from the Hamiltonian in terms of a
! trace correcting purification expansion with 2nd order purifications.
!----------------------------------------------------------------------
PROGRAM NT4p
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
  USE DenMatMethods
  USE MatFunk
  IMPLICIT NONE
  TYPE(BCSR)                     :: F,P,Pold,T,Z
  TYPE(DBL_RNK2)                 :: B,C
  !-------------------------------------------------------------------------------------
  ! Trace purserving NT4
  !-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  !
  REAL(DOUBLE)                   :: Energy,Energy_old,Thresh_old
  REAL(DOUBLE)                   :: ErrorE,ErrorN,ErrorP,ErrorFP,Ne
  REAL(DOUBLE)                   :: Degen,Occpan,lumo_occ,Gap
  !
  INTEGER                        :: I,Nr_Max_It,PNon0,MM
  LOGICAL                        :: Present,Converged
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='NT4'
  !-------------------------------------------------------------------------
  ! Start
  !-------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  !-------------------------------------------------------------------------
  CALL New(F)
  !-------------------------------------------------------------------------
  ! Get The Fock Matrix
  !-------------------------------------------------------------------------
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))   
  ENDIF
  !-------------------------------------------------------------------------
  ! Initialize                     
  !-------------------------------------------------------------------------
  MM         = 0                        
  Energy     = Zero
  Nr_Max_It  = 50                   
  Converged  = .FALSE.  
  Ne         = Half*DBLE(NEl)    
  Thresh_old = Thresholds%Trix
  !
  CALL New(P)
  CALL New(Pold)
  CALL New(T)    
  CALL SetEq(Pold,P)    
  !-------------------------------------------------------------------------
  !  Set up the starting Density Matrix from the Fock Matrix
  !-------------------------------------------------------------------------
  CALL FockGuess(F,P,Ne,1)
  !-------------------------------------------------------------------------
  ! Main Loop: Iterate until convergence              
  !-------------------------------------------------------------------------
  DO I = 1,Nr_Max_It
     CALL SetEq(Pold,P)
     Energy_old = Energy
     !--------------------------------------------------------------------------
     !    Set the Threshold for reduction of the error
     !--------------------------------------------------------------------------
     IF(I==1) THEN
        Thresholds%Trix = Thresh_old*1.D-2
     ELSEIF(I >= 2) THEN
        Thresholds%Trix = MIN(1.5D0*Thresholds%Trix,Thresh_old)
     ENDIF
     !--------------------------------------------------------------------------
     !    One Step of the Algorithm
     !--------------------------------------------------------------------------
     CALL NT4(P,P,Ne,MM,.TRUE.)
     !--------------------------------------------------------------------------
     !    Output Convergence Infomation
     !--------------------------------------------------------------------------
     Energy  = Trace(P,F)
     PNon0   = 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
     ErrorE  = ABS(Energy-Energy_old)/ABS(Energy)
     CALL Multiply(Pold,-One)
     CALL Add(Pold,P,T)
     ErrorP  = TwoNorm(T)/TwoNorm(P)
     CALL MednOut(Prog,I,Energy,PNon0,ErrorE,ErrorP)
     !     CALL CalculateGap(F,P,Ne,lumo_occ,Gap)
     !--------------------------------------------------------------------------
     !    Output Convergence Infomation
     !--------------------------------------------------------------------------
     IF(ErrorP < 5.0D0*Thresh_old) EXIT
  ENDDO
  !--------------------------------------------------------------------------
  ! Normalize Trace
  !--------------------------------------------------------------------------
  CALL NormTrace(P,Ne,1)
#ifdef WORK_ON_DEGENERACY
  !--------------------------------------------------------------------------
  !    Calculate Degeneracy, Occupancy, lumo_occ, and Gap
  !--------------------------------------------------------------------------
  CALL CalculateDegen(Ne,Degen,Occpan,lumo_occ)
  WRITE(*,*)' Degen    = ',Degen
  WRITE(*,*)' Occ      = ',Occpan
  WRITE(*,*)' LumoOcc= = ',Lumo_occ
  CALL New(B,(/NBasF,NBasF/))
  CALL New(C,(/NBasF,NBasF/))
  CALL SetEq(B,P)
  CALL SetDSYEVWork(NBasF)
  CALL FunkOnSqMat(NBasF,Inverse,B%D,C%D,PrintValues_O=.TRUE.,Unit_O=6)
  CALL Delete(B)
  CALL Delete(C)
  CALL UnSetDSYEVWork()
#endif
  !--------------------------------------------------------------------------
  ! Write Out Statisitcs
  !--------------------------------------------------------------------------
  Energy  = Trace(P,F)
  PNon0   = 100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
  ErrorE  = ABS(Energy-Energy_old)/ABS(Energy)
  ErrorN  = Trace(P)-Ne
  CALL Add(Pold,P,T)
  ErrorP  = TwoNorm(T)/TwoNorm(P)
  CALL Commute(F,P,T)
  ErrorFP = TwoNorm(T)/TwoNorm(F)
  CALL FinalOut(Prog,Energy,ErrorE,ErrorN,ErrorP,ErrorFP,PNon0,MM)
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
  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(Pold)
  CALL Delete(T)
  CALL Delete(Z)
  CALL ShutDown(Prog)
  !
END PROGRAM NT4p
