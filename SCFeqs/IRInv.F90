!    COMPUTE THE S^(-1/2) of the overlap Matrix
!    Author: Anders Niklasson and C. J. Tymczak
!--------------------------------------------------------------
PROGRAM IRInv
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE SetXYZ
  IMPLICIT NONE
!
  INTERFACE DSYEV
     SUBROUTINE DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
       USE GlobalScalars
       CHARACTER(LEN=1), INTENT(IN)    :: JOBZ, UPLO
       INTEGER,          INTENT(IN)    :: LDA,  LWORK, N
       INTEGER,          INTENT(OUT)   :: INFO
       REAL(DOUBLE),     INTENT(INOUT) :: A(LDA,*)
       REAL(DOUBLE),     INTENT(OUT)   :: W(*)
       REAL(DOUBLE),     INTENT(OUT)   :: WORK(*)
     END SUBROUTINE DSYEV
  END INTERFACE
!
  TYPE(BCSR)                     :: S,S0,SDg,SNotDg,SInvH,Tmp1,Tmp2,Tmp3
!
  TYPE(BSET)                     :: BS
  TYPE(CRDS)                     :: GM
  TYPE(ARGMT)                    :: Args
  TYPE(AtomPair)                 :: Pair 
  TYPE(DBL_VECT)                 :: Values,Work
  TYPE(DBL_RNK2)                 :: SBlk,InvSBlk,Vectors
  INTEGER                        :: I,J,K,AtA,AtB,NN,P,R,Sbeg,Send,   &
                                    NA,NB,NANBmax,LWork,Info
  REAL(DOUBLE)                   :: F2N,Beta,Factor,OneOvF2N,Error,Error_old     
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=8),PARAMETER     :: Prog='IRInv'
!--------------------------------------------------------------------------------------------------------
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Get the Overlap Matrix
  CALL Get(S,TrixFile('S',Args))
!-------------------------------------------------
! Find the Max and Min Eignevalue of S
  CALL New(Vectors,(/NBasF,NBasF/))
  CALL New(Values,NBasF)
  LWork=MAX(1,3*NBasF)
  CALL New(Work,LWork)
  CALL SetEq(Vectors,S)
  CALL DSYEV('V','U',NBasF,Vectors%D,NBasF,Values%D,Work%D,LWork,Info)
  WRITE(*,*) "S:"
  WRITE(*,*) "Smallest Eigenvalue = ",Values%D(1)
  WRITE(*,*) "Largest  Eigenvalue = ",Values%D(NBasF)
  Beta = Values%D(NBasF)*0.01D0
  CALL New(S0)
  CALL SetEq(S0,S)
  CALL Add(S0,Beta)
!-------------------------------------------------
! Iterate Using Anders Algorithm
!
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(Tmp3)
! Rescale S0 as the Initial Guess
  Factor = One/SQRT(Two*Values%D(NBasF)**3)
  CALL SetEq(SInvH,S0)
  CALL Multiply(SInvH,Factor)
! Calculate the Inverse of S
  DO J=1,2
     Error_old=1D8
     DO I=1,20
!       Calculate X
        CALL Multiply(SInvH,S0,Tmp1)
        CALL Multiply(Tmp1,SInvH,Tmp2)
        CALL SetEq(Tmp1,Tmp2)
!       Compute the Error
        CALL Add(Tmp2,-One)
        Error = MAX(Tmp2)
        WRITE(*,*) I,Error,100.D0*DBLE(SInvH%NNon0)/DBLE(NBasF*NBasF)
        IF(Error > Error_old) THEN
           CALL SetEq(SInvH,Tmp3)
           EXIT
        ENDIF
        IF(Error < 1.D-8) EXIT
        Error_old = Error
!       Calculate X2
        CALL Multiply(Tmp1,Tmp1,Tmp3)
        CALL SetEq(Tmp2,Tmp3)
!       Calculate 1.875*I-1.250*X+0.375*X*X
        CALL Multiply(Tmp1,-1.250D0)
        CALL Multiply(Tmp2, 0.375D0)
        CALL Add(Tmp1,1.875D0)
        CALL Add(Tmp1,Tmp2,Tmp3)
        CALL SetEq(Tmp1,Tmp3)
!       Calculate SInvH*(1.875*I-1.250*X+0.375*X*X)
        CALL Multiply(SInvH,Tmp1,Tmp2)
        CALL Multiply(Tmp1,SInvH,Tmp3)
        CALL Add(Tmp2,Tmp3,Tmp1)
        CALL Multiply(Tmp1,Half)
        CALL SetEq(Tmp3,SInvH)
        CALL Filter(SInvH,Tmp1)
     ENDDO
!    Rescale delta+S0
     WRITE(*,*) 
     CALL Add(S0,-Beta)
  ENDDO
! Calculte the Error In the Inverse
  CALL Multiply(SInvH,S,Tmp1)
  CALL Multiply(Tmp1,SInvH,Tmp2)
  CALL Add(Tmp2,-One)
  WRITE(*,*) "MaxError = ",MAX(Tmp2)
!
  IF(.TRUE.) STOP
!---------------------------------------------------
!!$    IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
!!$       CALL PPrint(PerfMon,Prog)  
!!$       CALL PPrint(PerfMon,Prog,Unit_O=6)  
!!$    ENDIF
!!$    ! Consistency check
!!$    IF(TEST_AINV)THEN
!!$       CALL New(T1)
!!$       CALL New(T2)
!!$       CALL Multiply(Zt,A,T1)
!!$       CALL Multiply(T1,Z,T2)
!!$       !
!!$       CALL SetToI(T1)
!!$       CALL Multiply(T1,-One)
!!$       CALL Add(T1,T2,A)
!!$       Mx0=Max(A)
!!$       !
!!$       IF(AInvDistanceThresh/=Zero)THEN     
!!$          Mssg='Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
!!$               //', DistThrsh='//TRIM(DblToShrtChar(AInvDistanceThresh))  &
!!$               //', TrixThrsh='//TRIM(DblToShrtChar(Thresholds%Trix))
!!$       ELSE
!!$          Mssg='Max(Z^t.A.Z-I)='//TRIM(DblToShrtChar(Mx0))                 &
!!$               //', TrixThrsh='//TRIM(DblToShrtChar(Thresholds%Trix))
!!$       ENDIF
!!$       IF(Mx0>1.D2*Thresholds%Trix)THEN
!!$          CALL Warn('In BlokAInv, failed test: '//TRIM(Mssg))
!!$       ENDIF
!!$       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!!$          Mssg=ProcessName(Prog)//TRIM(Mssg)
!!$          WRITE(*,*)TRIM(Mssg)
!!$          CALL OpenASCII(OutFile,Out)
!!$          CALL PrintProtectL(Out)
!!$          WRITE(Out,*)TRIM(Mssg)
!!$          CALL PrintProtectR(Out)
!!$          CLOSE(UNIT=Out,STATUS='KEEP')
!!$       ENDIF
!!$    ENDIF
!!$    !  Put Z and ZT to disk
!!$    CALL Put(Z,TrixFile('Z',Args))
!!$    CALL Put(Zt,TrixFile('ZT',Args))
!!$    !  Debug
!!$    CALL PChkSum(Z,'Z',Proc_O=Prog)
!!$    CALL PPrint(Z,'Z')
!!$    CALL Plot(Z,'Z')
!!$    !  Tidy up
!!$    CALL Delete(Z)
!!$    CALL Delete(ZT)
    CALL ShutDown(Prog) 
!
END PROGRAM IRInv



