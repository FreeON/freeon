!----------------------------------------------------------------------
PROGRAM DMP_SP4 ! Density matrix purification, SP4 variation
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
  TYPE(BCSR)                     :: F,P,Pold,Tmp1,Tmp2,Tmp3
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Ne
  INTEGER                        :: I,MM
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='SP4'
!-------------------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  ! Suss for matrix threshold overide
  CALL SussTrix('SPTwoTrix',Prog)  
  CALL New(F)
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))   
  ENDIF
! Initialize                     
  CALL New(P)
  CALL New(Pold)
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(Tmp3)
  MM=0                        
  Ne=Half*DBLE(NEl)    
  ! Guess P from F
  CALL FockGuess(F,P,Ne,1)
  CALL SetEq(Pold,P)    
  ! Do SP2 iterations
  DO I=1,100
     CALL SP4(P,Tmp1,Tmp2,Tmp3,Ne,MM)
     IF(CnvrgChck(Prog,I,Ne,MM,F,P,POld,Tmp1,Tmp2))EXIT
  ENDDO
  ! Orthogonal put and xform to AO rep and put
  CALL PutXForm(Prog,Args,P,POld,Tmp1)
  ! Tidy up
  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(Pold)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL Delete(Tmp3)
  CALL ShutDown(Prog)
END PROGRAM DMP_SP4




