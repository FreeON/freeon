PROGRAM DMP_PM ! Density matrix purification, PM variation
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
  TYPE(ARGMT)                    :: Args
  TYPE(BCSR)                     :: F,P,POld,Tmp1,Tmp2,Tmp3,Tmp4
  REAL(DOUBLE)                   :: Ne
  INTEGER                        :: I,MM
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=2),PARAMETER     :: Prog='PM'
!-------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  ! Suss for matrix threshold overide
  CALL SussTrix('PMTrix',Prog)  
  ! Get the Fock matrix
  CALL New(F)
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))   
  ENDIF
  ! Allocate some more matrices
  CALL New(P)
  CALL New(Pold)
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(Tmp3)
  CALL New(Tmp4)
  MM=0                        
  Ne=Half*DBLE(NEl)    
  ! Guess P from F
  CALL FockGuess(F,P,Ne,2)
  CALL SetEq(Pold,P)    
  ! Do PM iterations
  DO I=1,100 
      CALL PM2(P,Tmp1,Tmp2,Tmp3,MM)
!     CALL PM1(P,Tmp1,Tmp2,Tmp3,Tmp4,Ne,MM)
     IF(CnvrgChck(Prog,I,Ne,MM,F,P,POld,Tmp1,Tmp2))EXIT
  ENDDO
  ! Delete some obsolete matrices
  CALL Delete(F)
  ! Orthogonal put and xform to AO rep and put
  CALL PutXForm(Prog,Args,P,POld,Tmp1)
  ! Tidy up
  CALL Delete(P)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL Delete(Tmp3)
  CALL ShutDown(Prog)
END PROGRAM DMP_PM




