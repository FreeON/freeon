PROGRAM DMP_NT4 ! Density matrix purification, NT4 variation
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
! Trace purserving NT4
!-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Energy,Energy_old,Thresh_old
  REAL(DOUBLE)                   :: ErrorE,ErrorN,ErrorP,ErrorFP,Ne
  REAL(DOUBLE)                   :: Degen,Occpan,lumo_occ,Gap
  INTEGER                        :: I,Nr_Max_It,PNon0,MM
  LOGICAL                        :: Present,Converged
  REAL(DOUBLE),PARAMETER         :: FirstIter = 1.0D-2
  REAL(DOUBLE),PARAMETER         :: GrowFac   = 1.5D0
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='NT4'
!-------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
! Get The Fock Matrix
  CALL New(F)
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))   
  ENDIF
  ! Init
  CALL New(P)
  CALL New(Pold)
  CALL New(T)    
  CALL New(P2)
  CALL New(P3)
  CALL New(Ptmp1)
  CALL New(Ptmp2)
  ! Guess P from P
  CALL FockGuess(F,P,Ne,1)
  MM=0                        
  Ne=Half*DBLE(NEl)    
  CALL SetEq(Pold,P)    
  DO I=1,100 ! Main Loop
     CALL SetEq(Pold,P)
     CALL NT4(P,Ne,MM,.TRUE.)
     IF(CnvrgChck(Prog,I,MM,F,P,POld))EXIT
  ENDDO
  ! Normalize Trace
  CALL NormTrace(P,Ne,1)
  ! IO for the orthogonal P
  CALL Put(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)
  CALL Put(P,TrixFile('OrthoD',Args,1))
  CALL PChkSum(P,'OrthoP['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint( P,'OrthoP['//TRIM(NxtCycl)//']')
  CALL Plot(   P,'OrthoP_'//TRIM(NxtCycl))
  ! Convert to AO representation
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
  ! IO for the non-orthogonal P
  CALL Put(T,'CurrentDM',CheckPoint_O=.TRUE.)
  CALL Put(T,TrixFile('D',Args,1))
  CALL Put(Zero,'homolumogap')
  CALL PChkSum(T,'P['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(T,'P['//TRIM(NxtCycl)//']')
  CALL Plot(T,'P_'//TRIM(NxtCycl))
  ! Tidy up
  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(Pold)
  CALL Delete(T)
  CALL Delete(Z)
  CALL Delete(P2)
  CALL Delete(P3)
  CALL Delete(Ptmp1)
  CALL Delete(Ptmp2)
!
  CALL ShutDown(Prog)
END PROGRAM DMP_NT4




