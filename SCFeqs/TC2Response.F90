PROGRAM TC2R
!H=================================================================================
!H PROGRAM TC2Response
!H
!H  OPTIONS:
!H  DEBUGING: Use -DTC2R_DBUG.
!H  INFO    : Use -DTC2R_INFO.
!H  EIGENVAL: Use -DTC2R_EIGENVAL.
!H
!H Comment:
!H
!H
!H=================================================================================
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE Parse
  USE Macros
  USE LinAlg
  USE DenMatMethods
  !
  IMPLICIT NONE
  !
#ifdef TC2R_EIGENVAL
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
#endif
  !
  !------------------------------------------------------------------
  TYPE RINFO
     INTEGER :: Order
  END TYPE RINFO
  !
  !------------------------------------------------------------------
  !old TYPE(BCSR )                    :: F,FPrim
  TYPE(BCSR )                    :: F,FPrim,FPrm2,FPrm3
#ifdef PARALLEL
  !old TYPE(DBCSR)                    :: T,P,PPrim
  TYPE(DBCSR)                    :: T,P,PPrim_1,PPrim_2,PPrm2_1,PPrm2_2,PPrm3_1
  TYPE(DBCSR)                    :: PPrmOld,Tmp1,Tmp2,Tmp3
#else
  !old TYPE(BCSR )                    :: T,P,PPrim
  TYPE(BCSR )                    :: T,P,PPrim_1,PPrim_2,PPrm2_1,PPrm2_2,PPrm3_1
  TYPE(BCSR )                    :: PPrmOld,Tmp1,Tmp2,Tmp3
#endif
  TYPE(ARGMT)                    :: Args
  !-------------------------------------------------------------------
  INTEGER                        :: MM,I,LastSCFCycle,RespOrder
  REAL(DOUBLE)                   :: Ne
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: FFile
  CHARACTER(LEN=*), PARAMETER    :: Prog='TC2R'
  LOGICAL                        :: EXIST,IsPresent
  !-------------------------------------------------------------------
#ifdef TC2R_EIGENVAL
  INTEGER :: LWORK,Info
  TYPE(DBL_VECT)                 :: EigenV,Work
  TYPE(DBL_RNK2)                 :: PPrimA
#endif
  !-------------------------------------------------------------------
!!$  type(bcsr) :: work2


  !
  ! Initial setup.
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
  !
  ! Get Last SCF cycle.
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  RespOrder=1
  !
  !-------------------------------------------------------------------
  ! Allocations.
  !-------------------------------------------------------------------
  !
  CALL New(F       )
  CALL New(T       )
  CALL New(FPrim   )
  !
  IF(RespOrder.GE.2) CALL New(FPrm2)
  IF(RespOrder.GE.3) CALL New(FPrm3)
  !
  CALL New(P       )
  CALL New(Tmp1    )
  CALL New(Tmp2    )
  CALL New(Tmp3    )
  CALL New(PPrmOld)
  !
  CALL New(PPrim_1)
  IF(RespOrder.GE.2) CALL New(PPrim_2)
  IF(RespOrder.GE.2) CALL New(PPrm2_1)
  IF(RespOrder.GE.3) CALL New(PPrm2_2)
  IF(RespOrder.GE.3) CALL New(PPrm3_1)
  !
  !-------------------------------------------------------------------
  ! Loading matrices.
  !-------------------------------------------------------------------
  !
  ! Load Fock Matrix.
  CALL Get(F,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
  !
  ! Load FockPrim Matrices.
  SELECT CASE(RespOrder)
  CASE(1)
     ! Load FockPrime Matrix.
     FFile=TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0)
     INQUIRE(FILE=FFile,EXIST=IsPresent)
     IF(IsPresent) THEN
        CALL Get(FPrim,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     ELSE
        CALL Get(FPrim,TrixFile('OrthoFPrime' //TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     ENDIF
  CASE(2)
     ! Load FockPrime Matrix.
     CALL Get(F,TrixFile('OrthoFPrime',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.) !to change
     ! Load FockPrim2 Matrix.
     FFile=TrixFile('FPrim2_DDIIS'//TRIM(Args%C%C(4)),Args,0)
     INQUIRE(FILE=FFile,EXIST=IsPresent)
     IF(IsPresent) THEN
        CALL Get(FPrim,TrixFile('FPrim2_DDIIS'//TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     ELSE
        CALL Get(FPrim,TrixFile('OrthoFPrim2' //TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     ENDIF
  CASE(3)
     ! Load FockPrime Matrix.
     CALL Get(F,TrixFile('OrthoFPrime',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.) !to change
     ! Load FockPrim2 Matrix.
     CALL Get(F,TrixFile('OrthoFPrim2',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.) !to change
     ! Load FockPrim3 Matrix.
     FFile=TrixFile('FPrim3_DDIIS'//TRIM(Args%C%C(4)),Args,0)
     INQUIRE(FILE=FFile,EXIST=IsPresent)
     IF(IsPresent) THEN
        CALL Get(FPrim,TrixFile('FPrim3_DDIIS'//TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     ELSE
        CALL Get(FPrim,TrixFile('OrthoFPrim3' //TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     ENDIF
  CASE DEFAULT; STOP 'Problem with RespOrder'
  END SELECT


  !
  ! Load perturbation (dipole moment at this moment).
  FFile=TRIM(Args%C%C(3))//TRIM(Args%C%C(4))
  CALL Get(T,TrixFile(FFile,Args))
  !
  ! Initialize some variables.
  MM=0
  Ne=Half*DBLE(NEl)
  !
  !-------------------------------------------------------------------
  ! Get Density matrix guess.
  !-------------------------------------------------------------------
  !
  ! Guess P, PPrim_1 from F, FPrim.
  !old CALL FockPrimGuess(F,FPrim,P,PPrim_1)
  CALL FockPrimGuess(F,FPrim,FPrm2,FPrm3,P,PPrim_1,PPrm2_1,PPrm3_1,RespOrder)
  !
  ! Delete no more use buffers.
  CALL Delete(F    )
  CALL Delete(FPrim)
  IF(RespOrder.GE.2) CALL Delete(FPrm2)
  IF(RespOrder.GE.3) CALL Delete(FPrm3)
  !
  ! Set PPrmOld.
  SELECT CASE(RespOrder)
  CASE(1); CALL SetEq(PPrmOld,PPrim_1)
  CASE(2); CALL SetEq(PPrmOld,PPrm2_1)
  CASE(3); CALL SetEq(PPrmOld,PPrm3_1)
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Let's TC2R.
  !-------------------------------------------------------------------
  !
  ! Do TC2R iterations.

!!$  call new(work2)

  DO I=1,100
     !old CALL TC2R_DMP(P,PPrim_1,Tmp1,Tmp2,Tmp3,Ne,MM)
     CALL TC2R_DMP(P,PPrim_1,PPrim_2,PPrm2_1,PPrm2_2,PPrm3_1,Tmp1,Tmp2,Tmp3,Ne,MM,RespOrder)
     !old IF(CnvrgChckPrim(Prog,I,Ne,MM,T,PPrim_1,PPrmOld,Tmp1,Tmp2)) EXIT
     IF(CnvrgChckPrim(Prog,I,Ne,MM,T,PPrim_1,PPrm2_1,PPrm3_1,PPrmOld,Tmp1,Tmp2,RespOrder)) EXIT
  ENDDO
!!$  write(*,*) 'beta=',-2.0d0*Trace(work2,T)
!!$  call delete(work2)

  !
#ifdef TC2R_EIGENVAL
  CALL New(PPrimA,(/NBasF,NBasF/))
  CALL SetEq(PPrimA,PPrim_1)
  CALL New(EigenV,NBasF)
  CALL SetEq(EigenV,Zero)
  LWORK=MAX(1,3*NBasF+10)
  CALL New(Work,LWork)
  CALL DSYEV('V','U',NBasF,PPrimA%D(1,1),NBasF,EigenV%D(1),Work%D(1),LWORK,Info)
  IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in TC2R. INFO='//TRIM(IntToChar(Info)))
  CALL PLOT_ERROR_WITH_BOXES(EigenV,41,'test'//TRIM(IntToChar(Args%I%I(1))))
  CALL Delete(PPrimA)
  CAll Delete(EigenV)
  CALL Delete(Work)
#endif
  !
  !-------------------------------------------------------------------
  ! Save the matrices.
  !-------------------------------------------------------------------
  !
  ! Orthogonal put and xform to AO rep and put.
  !old CALL PutXFormPrim(Prog,Args,PPrim,Tmp1,Tmp2)
  SELECT CASE(RespOrder)
  CASE(1); CALL PutXFormPrim(Prog,Args,PPrim_1,Tmp1,Tmp2,'DPrime')
  CASE(2); CALL PutXFormPrim(Prog,Args,PPrm2_1,Tmp1,Tmp2,'DPrim2')
  CASE(3); CALL PutXFormPrim(Prog,Args,PPrm3_1,Tmp1,Tmp2,'DPrim3')
  CASE DEFAULT; STOP 'Problem with RespOrder'
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Deallocations the matrices.
  !-------------------------------------------------------------------
  !
  CALL Delete(P       )
  CALL Delete(T       )
  CALL Delete(Tmp1    )
  CALL Delete(Tmp2    )
  CALL Delete(Tmp3    )
  CALL Delete(PPrim_1 )
  IF(RespOrder.GE.2) CALL Delete(PPrim_2)
  IF(RespOrder.GE.2) CALL Delete(PPrm2_1)
  IF(RespOrder.GE.3) CALL Delete(PPrm2_2)
  IF(RespOrder.GE.3) CALL Delete(PPrm3_1)
  CALL Delete(PPrmOld)
  !
  CALL ShutDown(Prog)
  !
CONTAINS
  !
  !
  !old SUBROUTINE TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
  SUBROUTINE TC2R_DMP(P,PPrim_1,PPrim_2,PPrm2_1,PPrm2_2,PPrm3_1,Tmp1,Tmp2,Tmp3,Ne,MM,RespOrder)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
!H  This routine does the linear TC2Response scheme.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    !old TYPE(DBCSR) , INTENT(INOUT) :: P,PPrim,Tmp1,Tmp2,Tmp3
    TYPE(DBCSR) , INTENT(INOUT) :: P,PPrim_1,PPrim_2,PPrm2_1,PPrm2_2,PPrm3_1,Tmp1,Tmp2,Tmp3
#else
    !old TYPE(BCSR ) , INTENT(INOUT) :: P,PPrim,Tmp1,Tmp2,Tmp3
    TYPE(BCSR ) , INTENT(INOUT) :: P,PPrim_1,PPrim_2,PPrm2_1,PPrm2_2,PPrm3_1,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE), INTENT(IN   ) :: Ne
    INTEGER     , INTENT(IN   ) :: RespOrder
    INTEGER     , INTENT(INOUT) :: MM
    !-------------------------------------------------------------------
    REAL(DOUBLE)                :: N
    !-------------------------------------------------------------------
    !
    MM=MM+1
    N=Trace(P)
    IF(N>=Ne) THEN
       !----------------------
       ! Daaa = Daaa*X0 + X0*Daaa + Daa*Da + Da*Daa
       ! Dabc = Dabc*X0 + X0*Dabc + Dac*Db + Db*Dac + Dbc*Da + Da*Dbc
       ! PPrm3_1 <-> abc
       ! PPrm2_1 <-> ac
       ! PPrm2_2 <-> bc
       ! PPrim_1 <-> a
       ! PPrim_2 <-> b
       !----------------------
!!$       IF(RespOrder.GE.3) THEN
!!$          CALL Multiply(P,PPrm3_1,Tmp1)           !Tmp1=P*PPrm3
!!$          CALL Multiply(PPrm3_1,P,Tmp2)           !Tmp2=PPrm3*P
!!$          CALL Add(Tmp1,Tmp2,Tmp3)                !Tmp3=Tmp1+Tmp2
!!$          CALL Multiply(PPrm2_1,PPrm1_2,Tmp2)     !Tmp2=PPrm1*PPrm2
!!$          CALL Add(Tmp2,Tmp3,Tmp1)                !Tmp1=Tmp2+Tmp3
!!$          CALL Multiply(PPrm1_2,PPrm2_1,Tmp2)     !Tmp2=PPrm2*PPrm1
!!$          CALL Add(Tmp1,Tmp2,Tmp3)                !Tmp3=Tmp1+Tmp2
!!$          IF(A.EQ.B.AND.A.EQ.C) THEN
!!$             !----------------------
!!$             CALL Filter(PPrm3_1,Tmp3)            !PPrm3=Tmp3
!!$          ELSE
!!$             CALL Multiply(PPrm2_2,PPrm1_1,Tmp2)  !Tmp2=PPrm1*PPrm2
!!$             CALL Add(Tmp2,Tmp3,Tmp1)             !Tmp1=Tmp2+Tmp3
!!$
!!$             CALL Multiply(PPrm1_1,PPrm2_2,Tmp2)  !Tmp2=PPrm2*PPrm1
!!$             CALL Add(Tmp1,Tmp2,Tmp3)             !Tmp3=Tmp1+Tmp2
!!$             !----------------------
!!$             CALL Filter(PPrm3_1,Tmp3)            !PPrm3=Tmp3
!!$          ENDIF
!!$       ENDIF
       !
       !----------------------
       ! Daa = X0*Daa + Daa*X0 + Da*Da;
       ! Dab = X0*Dab + Dab*X0 + Da*Db + Db*Da;
       ! PPrm2_1 <-> ab
       ! PPrim_1 <-> a
       ! PPrim_2 <-> b
       !----------------------
!!$       IF(RespOrder.GE.2) THEN
!!$          CALL Multiply(P,PPrm2_1,Tmp1)         !Tmp1=P*PPrm2
!!$          CALL Multiply(PPrm2_1,P,Tmp2)         !Tmp2=PPrm2*P
!!$          CALL Add(Tmp1,Tmp2,Tmp3)              !Tmp3=Tmp1+Tmp2
!!$          CALL Multiply(PPrim_1,PPrim_2,Tmp1)   !
!!$          CALL Add(Tmp1,Tmp3,Tmp2)              !Tmp2=Tmp1+Tmp3
!!$          IF(A.EQ.B) THEN
!!$             !----------------------
!!$             CALL Filter(PPrm2_1,Tmp2)          !PPrm2=Tmp2
!!$          ELSE
!!$             CALL Multiply(PPrim_2,PPrim_1,Tmp1)!
!!$             CALL Add(Tmp1,Tmp2,Tmp3)           !Tmp3=Tmp1+Tmp2
!!$             !----------------------
!!$             CALL Filter(PPrm2_1,Tmp3)          !PPrm2=Tmp2
!!$          ENDIF
!!$       ENDIF

!!$       ! Daa = X0*Daa + Daa*X0 + Da*Da;
!!$       if(MM.le.1) then
!!$          CALL Multiply(PPrim_1,PPrim_1,Tmp1)   !
!!$          CALL Filter(work2,Tmp1)          !PPrm2=Tmp2
!!$       else
!!$          CALL Multiply(P,work2,Tmp1)         !Tmp1=P*PPrm2
!!$          CALL Multiply(work2,P,Tmp2)         !Tmp2=PPrm2*P
!!$          CALL Add(Tmp1,Tmp2,Tmp3)              !Tmp3=Tmp1+Tmp2
!!$          CALL Multiply(PPrim_1,PPrim_1,Tmp1)   !
!!$          CALL Add(Tmp1,Tmp3,Tmp2)              !Tmp2=Tmp1+Tmp3
!!$          CALL Filter(work2,Tmp2)          !PPrm2=Tmp2
!!$       endif

       !
       !----------------------
       ! D1 = X0*D1 + D1*X0;
       ! PPrim_1 <-> a
       !----------------------
       CALL Multiply(P,PPrim_1,Tmp1)    !Tmp1=P*PPrim
       CALL Multiply(PPrim_1,P,Tmp2)    !Tmp2=PPrim*P
       CALL Add(Tmp1,Tmp2,Tmp3)         !Tmp3=Tmp1+Tmp2 !PPrim=P*PPrim+PPrim*P
       !----------------------
       CALL Filter(PPrim_1,Tmp3)        !PPrim=Tmp3
       !----------------------
       !
       !----------------------
       ! X0 = X0^2;
       !----------------------
       CALL Multiply(P,P,Tmp3)        !Tmp3=P*P
       !----------------------
       CALL Filter(P,Tmp3)            !P=Tmp3
       !----------------------
    ELSE
       !----------------------
       ! Daaa = 2*Daaa - (Daaa*X0 + X0*Daaa + Daa*Da + Da*Daa)
       ! Dabc = 2*Dabc - (Dabc*X0 + X0*Dabc + Dac*Db + Db*Dac + Dbc*Da + Da*Dbc)
       ! PPrm3_1 <-> abc
       ! PPrm2_1 <-> ac
       ! PPrm2_2 <-> bc
       ! PPrim_1 <-> a
       ! PPrim_2 <-> b
       !----------------------
!!$       IF(RespOrder.GE.3) THEN
!!$          CALL Multiply(P,PPrm3_1,Tmp1)           !Tmp1=P*PPrm3
!!$          CALL Multiply(PPrm3_1,P,Tmp2)           !Tmp2=PPrm3*P
!!$          CALL Add(Tmp1,Tmp2,Tmp3)                !Tmp3=Tmp1+Tmp2
!!$          CALL Multiply(PPrm2_1,PPrm1_2,Tmp2)     !Tmp2=PPrm1*PPrm2
!!$          CALL Add(Tmp2,Tmp3,Tmp1)                !Tmp1=Tmp2+Tmp3
!!$          CALL Multiply(PPrm1_2,PPrm2_1,Tmp2)     !Tmp2=PPrm2*PPrm1
!!$          CALL Add(Tmp1,Tmp2,Tmp3)                !Tmp3=Tmp1+Tmp2
!!$          IF(A.EQ.B.AND.A.EQ.C) THEN
!!$             CALL Multiply(Tmp3,-One)       !Tmp3=-Tmp3
!!$             CALL Multiply(PPrm3_1,Two)       !PPrm3=2*PPrm3
!!$             CALL Add(PPrm3_1,Tmp3,Tmp1)      !Tmp1=PPrm3+Tmp3
!!$             !----------------------
!!$             CALL Filter(PPrm3_1,Tmp1)            !PPrm3=Tmp3
!!$          ELSE
!!$             CALL Multiply(PPrm2_2,PPrm1_1,Tmp2)  !Tmp2=PPrm1*PPrm2
!!$             CALL Add(Tmp2,Tmp3,Tmp1)             !Tmp1=Tmp2+Tmp3
!!$             CALL Multiply(PPrm1_1,PPrm2_2,Tmp2)  !Tmp2=PPrm2*PPrm1
!!$             CALL Add(Tmp1,Tmp2,Tmp3)             !Tmp3=Tmp1+Tmp2
!!$             CALL Multiply(Tmp3,-One)       !Tmp3=-Tmp3
!!$             CALL Multiply(PPrm3_1,Two)       !PPrm3=2*PPrm3
!!$             CALL Add(PPrm3_1,Tmp3,Tmp1)      !Tmp1=PPrm3+Tmp3
!!$             !----------------------
!!$             CALL Filter(PPrm3_1,Tmp1)            !PPrm3=Tmp3
!!$          ENDIF
!!$       ENDIF
       !
       !----------------------
       ! Daa = 2*Daa - (X0*Daa + Daa*X0 + Da*Da);
       ! Dab = 2*Dab - (X0*Dab + Dab*X0 + Da*Db + Db*Da);
       ! PPrm2_1 <-> ab
       ! PPrim_1 <-> a
       ! PPrim_2 <-> b
       !----------------------
!!$       IF(RespOrder.GE.2) THEN
!!$          CALL Multiply(P,PPrm2_1,Tmp1)         !Tmp1=P*PPrm2
!!$          CALL Multiply(PPrm2_1,P,Tmp2)         !Tmp2=PPrm2*P
!!$          CALL Add(Tmp1,Tmp2,Tmp3)              !Tmp3=Tmp1+Tmp2
!!$          CALL Multiply(PPrim_1,PPrim_2,Tmp1)   !
!!$          CALL Add(Tmp1,Tmp3,Tmp2)              !Tmp2=Tmp1+Tmp3
!!$          IF(A.EQ.B) THEN
!!$             CALL Multiply(Tmp2,-One)           !Tmp2=-Tmp2
!!$             CALL Multiply(PPrm2_1,Two)         !PPrm2=2*PPrm2
!!$             CALL Add(PPrm2_1,Tmp2,Tmp1)        !Tmp1=PPrm2+Tmp2
!!$             !----------------------
!!$             CALL Filter(PPrm2_1,Tmp1)          !PPrm2=Tmp2
!!$          ELSE
!!$             CALL Multiply(PPrim_2,PPrim_1,Tmp1)!
!!$             CALL Add(Tmp1,Tmp2,Tmp3)           !Tmp3=Tmp1+Tmp2
!!$             CALL Multiply(Tmp3,-One)           !Tmp2=-Tmp2
!!$             CALL Multiply(PPrm2_1,Two)         !PPrm2=2*PPrm2
!!$             CALL Add(PPrm2_1,Tmp3,Tmp1)        !Tmp1=PPrm2+Tmp2
!!$             !----------------------
!!$             CALL Filter(PPrm2_1,Tmp3)          !PPrm2=Tmp2
!!$          ENDIF
!!$       ENDIF


!!$       ! Daa = 2*Daa - (X0*Daa + Daa*X0 + Da*Da);
!!$       if(MM.le.1) then
!!$          CALL Multiply(PPrim_1,PPrim_1,Tmp1)   !
!!$          CALL Multiply(Tmp1,-One)              !Tmp2=-Tmp2
!!$          CALL Filter(work2,Tmp1)               !PPrm2=Tmp2
!!$       else
!!$          CALL Multiply(P,work2,Tmp1)         !Tmp1=P*PPrm2
!!$          CALL Multiply(work2,P,Tmp2)         !Tmp2=PPrm2*P
!!$          CALL Add(Tmp1,Tmp2,Tmp3)              !Tmp3=Tmp1+Tmp2
!!$          CALL Multiply(PPrim_1,PPrim_1,Tmp1)   !
!!$          CALL Add(Tmp1,Tmp3,Tmp2)              !Tmp2=Tmp1+Tmp3
!!$          CALL Multiply(Tmp2,-One)           !Tmp2=-Tmp2
!!$          CALL Multiply(work2,Two)         !PPrm2=2*PPrm2
!!$          CALL Add(work2,Tmp2,Tmp1)        !Tmp1=PPrm2+Tmp2
!!$          CALL Filter(work2,Tmp1)          !PPrm2=Tmp2
!!$       endif

       !
       !----------------------
       ! Da = 2*Da - (X0*Da+Da*X0);
       ! PPrim_1 <-> a
       !----------------------
       CALL Multiply(P,PPrim_1,Tmp1)    !Tmp1=P*PPrim
       CALL Multiply(PPrim_1,P,Tmp2)    !Tmp2=PPrim*P
       CALL Add(Tmp1,Tmp2,Tmp3)         !Tmp3=Tmp1+Tmp2
       CALL Multiply(Tmp3,-One)         !Tmp3=-Tmp3
       CALL Multiply(PPrim_1,Two)       !PPrim=2*PPrim
       CALL Add(PPrim_1,Tmp3,Tmp1)      !Tmp1=PPrim+Tmp3 
       !----------------------
       CALL Filter(PPrim_1,Tmp1)        !D1=Tmp1
       !
       !----------------------
       ! X0 = 2*X0-X02;
       !----------------------
       CALL Multiply(P,P,Tmp1)        !Tmp1=P*P
       CALL Multiply(Tmp1,-One)       !Tmp1=-Tmp1
       CALL Multiply(P,Two)           !P=2*P
       CALL Add(P,Tmp1,Tmp3)          !Tmp3=P+Tmp1
       !----------------------
       CALL Filter(P,Tmp3)            !X0=Tmp3
       !----------------------
    ENDIF
    !
  END SUBROUTINE TC2R_DMP
  !
  !
  !SUBROUTINE PutXFormPrim(Prog,Args,PPrim,Z,Tmp)
  SUBROUTINE PutXFormPrim(Prog,Args,PPrm,Z,Tmp,DName)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PutXFormPrim(Prog,Args,PPrim,Z,Tmp)
!H  Save onto the disc the OrthoPPrim and the PPrim matrices.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    !old TYPE(DBCSR)     , INTENT(INOUT) :: PPrim,Z,Tmp
    TYPE(DBCSR)     , INTENT(INOUT) :: PPrm,Z,Tmp
#else
    !old TYPE(BCSR )     , INTENT(INOUT) :: PPrim,Z,Tmp
    TYPE(BCSR )     , INTENT(INOUT) :: PPrm,Z,Tmp
#endif
    TYPE(ARGMT)     , INTENT(IN   ) :: Args
    CHARACTER(LEN=*), INTENT(IN   ) :: Prog,DName
    !-------------------------------------------------------------------
    LOGICAL                         :: IsPresent
    !-------------------------------------------------------------------
    !
    ! IO for the orthogonal PPrim
    !old CALL Put(PPrim,TrixFile('OrthoDPrime'//TRIM(Args%C%C(4)),Args,1))
    !old CALL PChkSum(PPrim,'OrthoDPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']',Prog)
    !old CALL PPrint( PPrim,'OrthoDPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']')
    !old CALL Plot(   PPrim,'OrthoDPrime'//TRIM(Args%C%C(4))//'_'//TRIM(NxtCycl))

    CALL Put(PPrm,TrixFile('Ortho'//TRIM(DName)//TRIM(Args%C%C(4)),Args,1))
    CALL PChkSum(PPrm,'Ortho'//TRIM(DName)//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( PPrm,'Ortho'//TRIM(DName)//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']')
    CALL Plot(   PPrm,'Ortho'//TRIM(DName)//TRIM(Args%C%C(4))//'_'//TRIM(NxtCycl))
     
    !
    ! Convert to AO representation
    INQUIRE(FILE=TrixFile('X',Args),EXIST=IsPresent)
    IF(IsPresent)THEN
       CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
       !old CALL Multiply(Z,PPrim,Tmp)
       !old CALL Multiply(Tmp,Z,PPrim)

       CALL Multiply(Z,PPrm,Tmp)
       CALL Multiply(Tmp,Z,PPrm)
    ELSE
       CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
       !old CALL Multiply(Z,PPrim,Tmp)
       CALL Multiply(Z,PPrm,Tmp)
       CALL Get(Z,TrixFile('ZT',Args))
       !old CALL Multiply(Tmp,Z,PPrim)
       CALL Multiply(Tmp,Z,PPrm)
    ENDIF
    !old CALL Filter(Tmp,PPrim)              ! Thresholding
    CALL Filter(Tmp,PPrm)              ! Thresholding
    !
    ! IO for the non-orthogonal PPrim
    !old CALL Put(Tmp,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,1))
    !old CALL PChkSum(Tmp,'DPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']',Prog)
    !old CALL PPrint( Tmp,'DPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']')
    !old CALL Plot(   Tmp,'DPrime'//TRIM(Args%C%C(4))//'_'//TRIM(NxtCycl))
    !
    CALL Put(Tmp,TrixFile(TRIM(DName)//TRIM(Args%C%C(4)),Args,1))
    CALL PChkSum(Tmp,TRIM(DName)//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( Tmp,TRIM(DName)//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']')
    CALL Plot(   Tmp,TRIM(DName)//TRIM(Args%C%C(4))//'_'//TRIM(NxtCycl))
    !
  END SUBROUTINE PutXFormPrim
  !
  !
  !old SUBROUTINE FockPrimGuess(F,FPrim,P,PPrim)
  SUBROUTINE FockPrimGuess(F,FPrim,FPrm2,FPrm3,P,PPrim,PPrm2,PPrm3,RespOrder)
!H---------------------------------------------------------------------------------
!H SUBROUTINE FockPrimGuess(F,FPrim,P,PPrim)
!H  Set the FockPrim guess for TC2R.
!H  This routine can be improved.
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    !old TYPE(BCSR ), INTENT(INOUT) :: F,FPrim
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrim,FPrm2,FPrm3
#ifdef PARALLEL
    !old TYPE(DBCSR), INTENT(  OUT) :: P,PPrim
    TYPE(DBCSR), INTENT(  OUT) :: P,PPrim,PPrm2,PPrm3
    TYPE(DBCSR)                :: Tmp
#else
    !old TYPE(BCSR ), INTENT(  OUT) :: P,PPrim
    TYPE(BCSR ), INTENT(  OUT) :: P,PPrim,PPrm2,PPrm3
    TYPE(BCSR )                :: Tmp
#endif
    INTEGER    , INTENT(IN   ) :: RespOrder
    !-------------------------------------------------------------------
    REAL(DOUBLE)               :: Fmin,Fmax,DF,Coeff
    !-------------------------------------------------------------------
    !
    ! Estimate spectral bounds.
#ifdef PARALLEL
    IF(MyId==ROOT) CALL SpectralBounds(F,Fmin,Fmax)
    CALL BCast(FMin)
    CALL BCast(FMax)
#else
    CALL SpectralBounds(F,Fmin,Fmax)
#endif
    !
    DF=(Fmax-Fmin)
    Coeff=-One/DF
    !
    ! Set up P.
    Call SetEq(P,F)
    CALL Add(P,-Fmax)
    CALL Multiply(P,Coeff)          ! P = (I*F_max-F)/DF = X0
    !
    ! Set up PPrim.
    CALL SetEq(PPrim,FPrim)
    CALL Add(PPrim,-Fmax)
    CALL Multiply(PPrim,-Coeff)     ! PPrim = -(I*F_max-FPrim)/DF = -X1
    !
    CALL New(Tmp)
    CALL Add(P,PPrim,Tmp)           ! Tmp = PPrim+P = X0-X1
    CALL SetEq(PPrim,Tmp)           ! PPrim = Tmp = X0-X1
    CALL Delete(Tmp)
    !
    CALL Multiply(PPrim,-One)       ! PPrim = -PPrim = X1-X0
    !
    IF(RespOrder.GE.2) THEN
       CALL SetEq(PPrm2,FPrm2)
       CALL Add(PPrm2,-Fmax)
       CALL Multiply(PPrm2,-Coeff)  ! PPrm2 = -(I*F_max-FPrm2)/DF = -X1
       !
       CALL New(Tmp)
       CALL Add(P,PPrm2,Tmp)        ! Tmp = PPrm2+P = X0-X1
       CALL SetEq(PPrm2,Tmp)        ! PPrm2 = Tmp = X0-X1
       CALL Delete(Tmp)
       !
       CALL Multiply(PPrm2,-One)    ! PPrm2 = -PPrm2 = X1-X0
    ENDIF
    !
    IF(RespOrder.GE.3) THEN
       CALL SetEq(PPrm3,FPrm3)
       CALL Add(PPrm3,-Fmax)
       CALL Multiply(PPrm3,-Coeff)  ! PPrm3 = -(I*F_max-FPrm3)/DF = -X1
       !
       CALL New(Tmp)
       CALL Add(P,PPrm3,Tmp)        ! Tmp = PPrm3+P = X0-X1
       CALL SetEq(PPrm3,Tmp)        ! PPrm3 = Tmp = X0-X1
       CALL Delete(Tmp)
       !
       CALL Multiply(PPrm3,-One)    ! PPrm3 = -PPrm3 = X1-X0
    ENDIF
    !
  END SUBROUTINE FockPrimGuess
  !
  !
  !old LOGICAL FUNCTION CnvrgChckPrim(Prog,NPur,Ne,MM,T,PPrim,PPrmOld,Tmp1,Tmp2)
  LOGICAL FUNCTION CnvrgChckPrim(Prog,NPur,Ne,MM,T,PPrim,PPrm2,PPrm3,PPrmOld,Tmp1,Tmp2,RespOrder)
!H---------------------------------------------------------------------------------
!H FUNCTION CnvrgChckPrim(Prog,NPur,Ne,MM,T,PPrim,PPrmOld,Tmp1,Tmp2)
!H  Check the convergence during TC2R steps.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    !old TYPE(DBCSR)                     :: T,PPrim,PPrmOld,Tmp1,Tmp2
    TYPE(DBCSR)                     :: T,PPrim,PPrm2,PPrm3,PPrmOld,Tmp1,Tmp2
#else
    !old TYPE(BCSR )                     :: T,PPrim,PPrmOld,Tmp1,Tmp2
    TYPE(BCSR )                     :: T,PPrim,PPrm2,PPrm3,PPrmOld,Tmp1,Tmp2
#endif
    REAL(DOUBLE)    , INTENT(IN   ) :: Ne
    CHARACTER(LEN=*), INTENT(IN   ) :: Prog
    INTEGER         , INTENT(IN   ) :: NPur,RespOrder
    INTEGER         , INTENT(INOUT) :: MM
    !-------------------------------------------------------------------
    REAL(DOUBLE)                    :: Prop,AbsErrProp,RelErrProp
    REAL(DOUBLE)                    :: AveErrProp,PNon0,PPrmNon0
    REAL(DOUBLE)                    :: AbsErrPPrim,FNormErrPPrim,TRacePPrim
    REAL(DOUBLE)    ,SAVE           :: OldProp,OldAEP
    CHARACTER(LEN=2*DEFAULT_CHR_LEN):: Mssg,CnvrgCmmnt
    !-------------------------------------------------------------------
    !
    IF(NPur==0)THEN
       OldProp=BIG_DBL
       OldAEP=BIG_DBL
    ENDIF
#ifdef PARALLEL
    PNon0=DBLE(Reduce(P%NNon0))
    PNon0=100.D0*PNon0/DBLE(NBasF*NBasF)
    SELECT CASE(RespOrder)
    CASE(1)
       PPrmNon0=DBLE(Reduce(PPrim%NNon0))
       PPrmNon0=100.D0*DBLE(PPrmNon0)/DBLE(NBasF*NBasF)
    CASE(2)
       PPrmNon0=DBLE(Reduce(PPrm2%NNon0))
       PPrmNon0=100.D0*DBLE(PPrmNon0)/DBLE(NBasF*NBasF)
    CASE(3)
       PPrmNon0=DBLE(Reduce(PPrm3%NNon0))
       PPrmNon0=100.D0*DBLE(PPrmNon0)/DBLE(NBasF*NBasF)
    END SELECT
#else
    PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
    SELECT CASE(RespOrder)
    CASE(1); PPrmNon0=100.D0*DBLE(PPrim%NNon0)/DBLE(NBasF*NBasF)
    CASE(2); PPrmNon0=100.D0*DBLE(PPrm2%NNon0)/DBLE(NBasF*NBasF)
    CASE(3); PPrmNon0=100.D0*DBLE(PPrm3%NNon0)/DBLE(NBasF*NBasF)
    END SELECT
#endif
    !
    ! Density matrix errors.
    CALL Multiply(PPrmOld,-One)

    !old CALL Add(PPrmOld,PPrim,Tmp1)
    SELECT CASE(RespOrder)
    CASE(1); CALL Add(PPrmOld,PPrim,Tmp1)
    CASE(2); CALL Add(PPrmOld,PPrm2,Tmp1)
    CASE(3); CALL Add(PPrmOld,PPrm3,Tmp1)
    END SELECT

    AbsErrPPrim=ABS(Max(Tmp1)+1.D-20)
    !
    ! Check for divergence steps.
    IF(AbsErrPPrim.GT.1D+10) &
         & CALL Halt('The absolute derivative density matrix error is too big! AbsErrPPrim='&
         &            //TRIM(DblToShrtChar(AbsErrPPrim)))
    !
    ! Compute Frobenus norm.
    FNormErrPPrim=FNorm(Tmp1)
    !
    ! Energy errors
#ifdef PARALLEL
    CALL Multiply(PPrim,T,Tmp1)
    Prop=-Two*Trace(Tmp1)
#else
    !old Prop=-Two*Trace(PPrim,T)
    SELECT CASE(RespOrder)
    CASE(1); Prop=-Two*Trace(PPrim,T)
    CASE(2); Prop=-Two*Trace(PPrm2,T) !to check
    CASE(3); Prop=-Two*Trace(PPrm3,T) !to check
    END SELECT
#endif
    !write(*,*) 'Prop=',Prop
    !
    ! Compute Abs and Rel error.
    AbsErrProp=ABS(OldProp-Prop)
    ! should check if it is zero
    RelErrProp=AbsErrProp/ABS(Prop)
    ! write(*,*) 'AbsErrProp',AbsErrProp,'AbsErrPPrim',AbsErrPPrim
    ! write(*,*) 'Thresholds%ETol',Thresholds%ETol,'Thresholds%DTol',Thresholds%DTol
    !
    ! Set convergence check.
    CnvrgChckPrim=.FALSE.
    !
    ! Absolute convergence test.
    !IF(RelErrProp<Thresholds%ETol*1D-2.AND.AbsErrPPrim<Thresholds%DTol*1D-1)THEN
    !IF(RelErrProp<Thresholds%ETol*1D-0.AND.AbsErrPPrim<Thresholds%DTol*1D+1)THEN
    !IF(AbsErrProp<1.0d-5.AND.AbsErrPPrim<1.0d-4)THEN
    !IF(RelErrProp<1.0d-5.AND.AbsErrPPrim<1.0d-4)THEN
    !IF(RelErrProp<SQRT(Thresholds%ETol*1D-2).AND.AbsErrPPrim<Thresholds%DTol*1D-1)THEN
    IF(RelErrProp<SQRT(Thresholds%ETol).AND.AbsErrPPrim<Thresholds%DTol) THEN
       CnvrgChckPrim=.TRUE.
       !old CnvrgCmmnt='Met dE'' goals'
       SELECT CASE(RespOrder)
       CASE(1); CnvrgCmmnt='Met dE1 goals'
       CASE(2); CnvrgCmmnt='Met dE2 goals'
       CASE(3); CnvrgCmmnt='Met dE3 goals'
       END SELECT
    ENDIF
    !
    ! Test in the asymptotic regime for stall out
    !IF(RelErrProp<1D2*Thresholds%Trix**2)THEN
    IF(RelErrProp<Thresholds%ETol)THEN
       ! Check for increasing /P
       IF(AbsErrPPrim>OldAEP)THEN
          CnvrgChckPrim=.TRUE.
          !old CnvrgCmmnt='dP'' increase'
          SELECT CASE(RespOrder)
          CASE(1); CnvrgCmmnt='dP1 increase'
          CASE(2); CnvrgCmmnt='dP2 increase'
          CASE(3); CnvrgCmmnt='dP3 increase'
          END SELECT
       ENDIF
       ! Check for an increasing energy
       !   IF(Prop>OldProp)THEN
       !       CnvrgChckPrim=.TRUE.
       !       CnvrgCmmnt='dE'' increase'
       !   ENDIF
    ENDIF
    !
    ! Updtate previous cycle values
    OldProp=Prop
    OldAEP=AbsErrPPrim
    !
    ! Print convergence stats
    !old Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))        &
    !old      //'dE''='//TRIM(DblToShrtChar(RelErrProp))              &
    !old      //', dP''='//TRIM(DblToShrtChar(AbsErrPPrim))           &
    !old      //', %Non0='//TRIM(DblToShrtChar(PNon0))              

    SELECT CASE(RespOrder)
    CASE(1)
       Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))       &
            //'dE1='//TRIM(DblToShrtChar(RelErrProp))              &   
            //', dP1='//TRIM(DblToShrtChar(AbsErrPPrim))           &
            //', %Non0='//TRIM(DblToShrtChar(PPrmNon0))              
    CASE(2)
       Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))       &
            //'dE2='//TRIM(DblToShrtChar(RelErrProp))              &   
            //', dP2='//TRIM(DblToShrtChar(AbsErrPPrim))           &
            //', %Non0='//TRIM(DblToShrtChar(PPrmNon0))              
    CASE(3)
       Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))       &
            //'dE3='//TRIM(DblToShrtChar(RelErrProp))              &   
            //', dP3='//TRIM(DblToShrtChar(AbsErrPPrim))           &
            //', %Non0='//TRIM(DblToShrtChar(PPrmNon0))              
    END SELECT

    IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       CALL OpenASCII(OutFile,Out)
       CALL PrintProtectL(Out)
       WRITE(*,*)TRIM(Mssg)
       WRITE(Out,*)TRIM(Mssg)
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL
    ENDIF
#endif
    ENDIF
    !
    ! Set thresholding for next cycle
    CALL SetVarThresh(MM)
    !
    ! Look for convergence
    IF(.NOT.CnvrgChckPrim)THEN
       !old CALL SetEq(PPrmOld,PPrim)
       SELECT CASE(RespOrder)
       CASE(1); CALL SetEq(PPrmOld,PPrim)
       CASE(2); CALL SetEq(PPrmOld,PPrm2)
       CASE(3); CALL SetEq(PPrmOld,PPrm3)
       END SELECT
       RETURN
    ENDIF
    !
    ! Increment the TC2R counter.
    MM=MM+1
    !
    ! Print summary stats.
    !old TracePPrim = Trace(PPrim)
    SELECT CASE(RespOrder)
    CASE(1); TracePPrim = Trace(PPrim)
    CASE(2); TracePPrim = Trace(PPrm2)
    CASE(3); TracePPrim = Trace(PPrm3)
    END SELECT
    !
    IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       CALL OpenASCII(OutFile,Out)
       CALL PrintProtectL(Out)
       Mssg=ProcessName(Prog,CnvrgCmmnt)                                             &
            &                //'Tr{PT}='//TRIM(DblToChar(Prop))
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       !
       Mssg=ProcessName(Prog)//TRIM(IntToChar(NPur))//' purification steps, '        &
            &                //TRIM(IntToChar(MM))//' matrix multiplies'
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       Mssg=ProcessName(Prog)//'ThrX='//TRIM(DblToShrtChar(Thresholds%Trix))         &
            &                //', %Non0s = '//TRIM(DblToShrtChar(PPrmNon0))              
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       !old Mssg=ProcessName(Prog,'Max abs errors') &
       !old      &                //'dE''='//TRIM(DblToShrtChar(RelErrProp))//', '        &
       !old      &                //'dP''='//TRIM(DblToShrtChar(AbsErrPPrim))
       SELECT CASE(RespOrder)
       CASE(1) 
          Mssg=ProcessName(Prog,'Max abs errors') &
               &                //'dE1='//TRIM(DblToShrtChar(RelErrProp))//', '        &
               &                //'dP1='//TRIM(DblToShrtChar(AbsErrPPrim))
       CASE(2)
          Mssg=ProcessName(Prog,'Max abs errors') &
               &                //'dE2='//TRIM(DblToShrtChar(RelErrProp))//', '        &
               &                //'dP2='//TRIM(DblToShrtChar(AbsErrPPrim))
       CASE(3)
          Mssg=ProcessName(Prog,'Max abs errors') &
               &                //'dE3='//TRIM(DblToShrtChar(RelErrProp))//', '        &
               &                //'dP3='//TRIM(DblToShrtChar(AbsErrPPrim))
       END SELECT
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)

       !old Mssg=ProcessName(Prog)//'Rel dE''='//TRIM(DblToShrtChar(RelErrProp))//', '    &
       !old      &                //'||dP''||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       SELECT CASE(RespOrder)
       CASE(1)
          Mssg=ProcessName(Prog)//'Rel dE1='//TRIM(DblToShrtChar(RelErrProp))//', '    &
               &                //'||dP1||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       CASE(2)
          Mssg=ProcessName(Prog)//'Rel dE2='//TRIM(DblToShrtChar(RelErrProp))//', '    &
               &                //'||dP2||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       CASE(3)
          Mssg=ProcessName(Prog)//'Rel dE3='//TRIM(DblToShrtChar(RelErrProp))//', '    &
               &                //'||dP3||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       END SELECT
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL
    ENDIF
#endif
    ENDIF
    !
  END FUNCTION CnvrgChckPrim
  !  
  !
#ifdef TC2R_EIGENVAL
  SUBROUTINE PLOT_ERROR_WITH_BOXES(Eig,NBrBox,Name)
    TYPE(DBL_VECT)               :: Eig
    INTEGER         , INTENT(IN) :: NBrBox
    CHARACTER(LEN=*), INTENT(IN) :: Name
    TYPE(INT_VECT)               :: Count
    REAL(DOUBLE)                 :: MaxE,MinE,DeltaE
    INTEGER                      :: ICount,IBox,IBas,Plt
    !
    CALL New(Count,NBrBox)
    !
    MaxE=MAXVAL(Eig%D(1:NBasF))
    MinE=MINVAL(Eig%D(1:NBasF))
    !
    DeltaE=(MaxE-MinE)/DBLE(NBrBox)
    !
    DO IBox=1,NBrBox
       ICount=0
       DO IBas=1,NBasF
          IF(MinE+DBLE((IBox-1))*DeltaE.LE.Eig%D(IBas).AND.Eig%D(IBas).LE.MinE+DBLE(IBox)*DeltaE) ICount=ICount+1
       ENDDO
       Count%I(IBox)=ICount
    ENDDO
    !
!    IF(SUM(Count%I).NE.NBasF) STOP 11111
    !
    CALL OpenASCII(TRIM(Name)//'_PLOT_BAR_data',Plt,NewFile_O=.TRUE.)
    DO IBox=1,NBrBox
       WRITE(Plt,*) MinE+(DBLE(IBox-1)+0.5d0)*DeltaE ,Count%I(IBox)
    ENDDO
    CLOSE(Plt)
    !
    CALL OpenASCII(TRIM(Name)//'_PLOT_BAR_PlotMe',Plt,NewFile_O=.TRUE.)
    WRITE(Plt,2)
    WRITE(Plt,3)TRIM(Name)//'.eps'
    WRITE(Plt,6)
    WRITE(Plt,*)'set pointsize 0.1'
!    WRITE(Plt,*)'set logscale y'
    WRITE(Plt,*)"plot '"//TRIM(Name)//"_PLOT_BAR_data' using 1:2 notitle with boxes "
    CLOSE(Plt)
    !
    CALL Delete(Count)
    !
   1   FORMAT(2(1x,I16))
   2   FORMAT('set term  postscript eps  "Times-Roman" 18')
!   2   FORMAT('set term  jpeg transparent')
   3   FORMAT('set output "',A,'"')
   6   FORMAT('set size square ')
   9   FORMAT('plot [0 : ',I12,' ] ',I12,', \\')
  10   FORMAT('                    ',I12,', \\')
  END SUBROUTINE PLOT_ERROR_WITH_BOXES
#endif
  !
  !
END PROGRAM TC2R





