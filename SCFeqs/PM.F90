!----------------------------------------------------------------------
!                         May 13th 2002
! Anders M. N. Niklasson: "Expansion Algorithm for the Density Matrix".
! Constructs the density matrix from the Hamiltonian in terms of a
! trace correcting purification expansion with 2nd order purifications.
!----------------------------------------------------------------------
PROGRAM PMp
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
  IMPLICIT NONE
  TYPE(BCSR)                     :: F,P,Pold,T,Z
  !-------------------------------------------------------------------------------------
  ! Trace perserving Palser-Manopolus
  !-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  !
  REAL(DOUBLE)                   :: Energy,Energy_old,Thresh_old
  REAL(DOUBLE)                   :: ErrorE,ErrorN,ErrorP,ErrorFP,Ne
  !
  INTEGER                        :: I,Nr_Max_It,PNon0,MM
  LOGICAL                        :: Present,Converged
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=2),PARAMETER     :: Prog='PM'
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
  CALL FockGuess(F,P,Ne,2)
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
        Thresholds%Trix = MIN(1.25D0*Thresholds%Trix,Thresh_old)
     ENDIF
     !--------------------------------------------------------------------------
     !    One Step of the Algorithm
     !--------------------------------------------------------------------------
     CALL PM(P,P,Ne,MM,.TRUE.)
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
     !--------------------------------------------------------------------------
     !    Output Convergence Infomation
     !--------------------------------------------------------------------------
     IF(ErrorP < 5.0D0*Thresh_old) EXIT
  ENDDO
  !--------------------------------------------------------------------------
  ! Normalize Trace
  !--------------------------------------------------------------------------
  CALL NormTrace(P,Ne,1)
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
  CALL ShutDown(Prog)
  IF(.TRUE.) STOP
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
END PROGRAM PMp
