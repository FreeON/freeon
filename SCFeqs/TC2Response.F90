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
  USE DenMatMethods, ONLY: SpectralBounds,SetVarThresh
  USE DenMatResponse
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
  TYPE(BCSR )                    :: F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1
#ifdef PARALLEL
  TYPE(DBCSR)                    :: T,P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
  TYPE(DBCSR)                    :: PPrmOld,Tmp1,Tmp2,Tmp3
#else
  TYPE(BCSR )                    :: T,P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1
  TYPE(BCSR )                    :: PPrmOld,Tmp1,Tmp2,Tmp3
#endif
  TYPE(ARGMT)                    :: Args
  !-------------------------------------------------------------------
  INTEGER                        :: MM,I,RespOrder
  INTEGER                        :: LastSCFCycle,LastCPSCFCycle
  REAL(DOUBLE)                   :: Ne,SwitchThresh
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: FFile
  CHARACTER(LEN=1)               :: Chr1,Chr2,Chr3
  CHARACTER(LEN=*), PARAMETER    :: Prog='TC2R'
  LOGICAL                        :: EXIST,IsPresent
  !-------------------------------------------------------------------
  REAL(DOUBLE)    , PARAMETER    :: DEFAULT_SWITCHTHRESH=1.00D-3
  !-------------------------------------------------------------------
  !
  ! Initial setup.
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
  !
  ! Open input
  CALL OpenASCII(InpFile,Inp)  
  !
  ! Switching.
  IF(.NOT.OptDblQ(Inp,'TC2RSwitchingThresh',SwitchThresh)) SwitchThresh=DEFAULT_SWITCHTHRESH
  !
  ! Close the input.
  CLOSE(Inp)
  !
  ! Get the response order.
  RespOrder=LEN(TRIM(Args%C%C(3)))
  !
  !-------------------------------------------------------------------
  ! Allocations.
  !-------------------------------------------------------------------
  !
  CALL AllocArray(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
       &          P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
       &          PPrmOld,T,Tmp1,Tmp2,Tmp3,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Loading matrices.
  !-------------------------------------------------------------------
  !
  CALL LoadMatrices(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
       &            P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
       &            T,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Get Density matrix guess.
  !-------------------------------------------------------------------
  !
  CALL FockPrimGuess(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
       &             P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
       &             RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Set up the old densiry derivative.
  !-------------------------------------------------------------------
  !
  CALL SetPPrmOld(PPrmOld,PPrm1_1,PPrm2_1,PPrm3_1,RespOrder)
  !
  !-------------------------------------------------------------------
  ! Let's TC2R.
  !-------------------------------------------------------------------
  !
  CALL DoTC2R(P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1,PPrmOld, &
       &      T,Tmp1,Tmp2,Tmp3,RespOrder,SwitchThresh,Args)
  !
  !-------------------------------------------------------------------
  ! Save the matrices.
  !-------------------------------------------------------------------
  !
  CALL SavePPrm(PPrm1_1,PPrm2_1,PPrm3_1,Tmp1,Tmp2,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! Deallocations the matrices.
  !-------------------------------------------------------------------
  !
  CALL DeAllocArray(F,FPrm1_1,FPrm1_2,FPrm1_3,FPrm2_1,FPrm2_2,FPrm2_3,FPrm3_1, &
       &            P,PPrm1_1,PPrm1_2,PPrm1_3,PPrm2_1,PPrm2_2,PPrm2_3,PPrm3_1, &
       &            PPrmOld,T,Tmp1,Tmp2,Tmp3,RespOrder,Args)
  !
  !-------------------------------------------------------------------
  ! We are done.
  !-------------------------------------------------------------------
  !
  CALL ShutDown(Prog)
  !
END PROGRAM TC2R
