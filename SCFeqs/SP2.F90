!-------------------------------------------------------------------------------
!                         May 13th 2002
! Anders M. N. Niklasson: "Expansion Algorithm for the Density Matrix".
! Constructs the density matrix from the Hamiltonian in terms of a
! trace correcting purification expansion with 2nd order purifications.
!-------------------------------------------------------------------------------
PROGRAM DMP_SP2 ! Density matrix purification, SP2 variation
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
#ifdef PARALLEL
  TYPE(DBCSR)                     :: F,P,Pold,Tmp1,Tmp2
  TYPE(BCSR)                      :: F_BCSR
#else
  TYPE(BCSR)                     :: F,P,Pold,Tmp1,Tmp2
#endif
!-------------------------------------------------------------------------------
! Trace Setting SP2
!-------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Ne,Lambda
  INTEGER                        :: I,MM
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='SP2'
!-------------------------------------------------------------------------------
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
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
  MM=0                        
  Ne=Half*DBLE(NEl)    
  ! Guess P from F
#ifdef PARALLEL
  CALL SetEq(F_BCSR,F)
  CALL FockGuess(F_BCSR,P,Ne,1)
  CALL Delete(F_BCSR)
#else
  CALL FockGuess(F,P,Ne,1)
#endif
  CALL SetEq(Pold,P)    
  ! Do SP2 iterations
  DO I=1,100
     CALL SP2(P,Tmp1,Tmp2,Ne,MM)
     IF(CnvrgChck(Prog,I,Ne,MM,F,P,POld,Tmp1,Tmp2))EXIT
  ENDDO
  CALL OpenASCII(InpFile,Inp)         
  ! If we are called without the DIIS Fockian, consider a levelshift
  IF(.NOT.PRESENT.AND.OptDblQ(Inp,'LevelShift',Lambda))THEN
     ! Get the Fock matrix back ... 
     CALL Get(F,TrixFile('OrthoF',Args,0))   
     ! Construct the virtual projector, Q=I-P
     CALL SetEq(POld,P)
     CALL Multiply(POld,-One)   
     CALL Add(POld,One)     
     ! The shifted Fockian is F[Lambda] = P.F + (1+Lambda)*Q.F
     Lambda=One+ABS(Lambda)
     CALL Multiply(P,F,Tmp1)
     CALL Multiply(POld,F,Tmp2)
     CALL Multiply(Tmp2,Lambda)
     CALL Add(Tmp1,Tmp2,POld)
     ! Put the Fock matrix back all tidy like
     CALL Put(POld,TrixFile('OrthoF',Args,0))   
     Lambda=ABS(Lambda)-One
     Mssg=TRIM(ProcessName(Prog))//' LevlShift = '  &
          //TRIM(DblToMedmChar(Lambda))
     CALL OpenASCII(OutFile,Out)         
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
  ENDIF
  CLOSE(Inp)
  ! Orthogonal put and xform to AO rep and put
  CALL PutXForm(Prog,Args,P,POld,Tmp1)
  ! Tidy up
  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(Pold)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL ShutDown(Prog)
END PROGRAM DMP_SP2
