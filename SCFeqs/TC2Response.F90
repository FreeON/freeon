PROGRAM TC2R
!H=================================================================================
!H PROGRAM TC2Response
!H
!H  OPTIONS:
!H  DEBUGING: Use -DTC2R_DBUG.
!H  INFO    : Use -DTC2R_INFO.
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
  !------------------------------------------------------------------
  TYPE(BCSR )                    :: F,FPrim
#ifdef PARALLEL
  TYPE(DBCSR)                    :: T,P,PPrim
  TYPE(DBCSR)                    :: PPrimOld,Tmp1,Tmp2,Tmp3
#else
  TYPE(BCSR )                    :: T,P,PPrim
  TYPE(BCSR )                    :: PPrimOld,Tmp1,Tmp2,Tmp3
#endif
  TYPE(ARGMT)                    :: Args
  !-------------------------------------------------------------------
  INTEGER                        :: MM,I,LastSCFCycle
  REAL(DOUBLE)                   :: Ne
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: FFile
  CHARACTER(LEN=4), PARAMETER    :: Prog='TC2R'
  LOGICAL                        :: EXIST,IsPresent
  !-------------------------------------------------------------------
  !
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  !
  ! Get Last SCF cycle.
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  ! Load Fock Matrix.
  CALL New(F)
  CALL Get(F,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
  !CALL Get(F,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)))
  !
  ! Load FockPrim Matrix.
  CALL New(FPrim)
  FFile=TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0)
  INQUIRE(FILE=FFile,EXIST=IsPresent)
  IF(IsPresent) THEN
     CALL Get(FPrim,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     !CALL Get(FPrim,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0))
  ELSE
     CALL Get(FPrim,TrixFile('OrthoFPrime' //TRIM(Args%C%C(4)),Args,0),BCast_O=.FALSE.)
     !CALL Get(FPrim,TrixFile('OrthoFPrime' //TRIM(Args%C%C(4)),Args,0))
  ENDIF
  !
  ! Load perturbation (dipole moment at this moment).
  CALL New(T)
  FFile=''
  DO I=3,4!Args%NC
     FFile=TRIM(FFile)//TRIM(Args%C%C(I))
  ENDDO
  CALL Get(T,TrixFile(FFile,Args))
  !
  ! Allocate some arrays.
  CALL New(P       )
  CALL New(Tmp1    )
  CALL New(Tmp2    )
  CALL New(Tmp3    )
  CALL New(PPrim   )
  CALL New(PPrimOld)
  !
  ! Initialize some variables.
  MM=0
  Ne=Half*DBLE(NEl)
  !
  ! Guess P, PPrim from F, FPrim.
  CALL FockPrimGuess(F,FPrim,P,PPrim)
  !
  ! Delete no more use buffers.
  CALL Delete(F    )
  CALL Delete(FPrim)
  !
  ! Set PPrimOld.
  CALL SetEq(PPrimOld,PPrim)
  !
  ! Do TC2R iterations.
  DO I=1,100
     CALL TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
     IF(CnvrgChckPrim(Prog,I,Ne,MM,T,PPrim,PPrimOld,Tmp1,Tmp2)) EXIT
  ENDDO
  !
  ! Orthogonal put and xform to AO rep and put.
  CALL PutXFormPrim(Prog,Args,PPrim,Tmp1,Tmp2)
  !
  ! Tidy up.
  CALL Delete(P       )
  CALL Delete(T       )
  CALL Delete(Tmp1    )
  CALL Delete(Tmp2    )
  CALL Delete(Tmp3    )
  CALL Delete(PPrim   )
  CALL Delete(PPrimOld)
  !
  CALL ShutDown(Prog)
  !
CONTAINS
  !
  !
  SUBROUTINE TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
!H  This routine does the linear TC2Response scheme.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR) , INTENT(INOUT) :: P,PPrim,Tmp1,Tmp2,Tmp3
#else
    TYPE(BCSR ) , INTENT(INOUT) :: P,PPrim,Tmp1,Tmp2,Tmp3
#endif
    REAL(DOUBLE), INTENT(IN   ) :: Ne
    INTEGER     , INTENT(INOUT) :: MM
    !-------------------------------------------------------------------
    REAL(DOUBLE)                :: N
    !-------------------------------------------------------------------
    !
    MM=MM+1
    N=Trace(P)
    IF(N>=Ne) THEN
       ! D3 = X0*D3 + D3*X0 + D1*D2+D2*D1;
       ! D2 = X0*D2 + D2*X0 + D1*D1;
       !
       ! D1 = X0*D1 + D1*X0;
       CALL Multiply(P,PPrim,Tmp1)    !Tmp1=P*PPrim
       CALL Multiply(PPrim,P,Tmp2)    !Tmp2=PPrim*P
       CALL Add(Tmp1,Tmp2,Tmp3)       !Tmp3=Tmp1+Tmp2 !PPrim=P*PPrim+PPrim*P
       !----------------------
       !IF(MOD(MM,2).EQ.0) THEN
       !CALL SetEq(PPrim,Tmp3)
       !ELSE
       CALL Filter(PPrim,Tmp3)        !PPrim=Tmp3
       !ENDIF
       !----------------------
       ! X0 = X02;
       CALL Multiply(P,P,Tmp3)        !Tmp3=P*P
       !----------------------
       !CALL SetEq(P,Tmp3)
       CALL Filter(P,Tmp3)            !P=Tmp3
       !----------------------
    ELSE
       ! D3 = 2*D3 - (X0*D3+D3*X0) - (D1*D2+D2*D1);
       ! D2 = 2*D2 - (X0*D2+D2*X0) - (D1*D1);
       !
       ! D1 = 2*D1 - (X0*D1+D1*X0);      
       CALL Multiply(P,PPrim,Tmp1)    !Tmp1=P*PPrim
       CALL Multiply(PPrim,P,Tmp2)    !Tmp2=PPrim*P
       CALL Add(Tmp1,Tmp2,Tmp3)       !Tmp3=Tmp1+Tmp2
       CALL Multiply(Tmp3,-One)       !Tmp3=-Tmp3
       CALL Multiply(PPrim,Two)       !PPrim=2*PPrim
       CALL Add(PPrim,Tmp3,Tmp1)      !Tmp1=PPrim+Tmp3 
       !----------------------
       !IF(MOD(MM,2).EQ.0) THEN
       !CALL SetEq(PPrim,Tmp1)
       !ELSE
       CALL Filter(PPrim,Tmp1)        !D1=Tmp1
       !ENDIF
       !----------------------
       ! X0 = 2*X0-X02;
       CALL Multiply(P,P,Tmp1)        !Tmp1=P*P
       CALL Multiply(Tmp1,-One)       !Tmp1=-Tmp1
       CALL Multiply(P,Two)           !P=2*P
       CALL Add(P,Tmp1,Tmp3)          !Tmp3=P+Tmp1
       !----------------------
       !CALL SetEq(P,Tmp3)
       CALL Filter(P,Tmp3)            !X0=Tmp3
       !----------------------
    ENDIF
    !
  END SUBROUTINE TC2R_DMP
  !
  !
  SUBROUTINE PutXFormPrim(Prog,Args,PPrim,Z,Tmp)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PutXFormPrim(Prog,Args,PPrim,Z,Tmp)
!H  Save onto the disc the OrthoPPrim and the PPrim matrices.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR)     , INTENT(INOUT) :: PPrim,Z,Tmp
#else
    TYPE(BCSR )     , INTENT(INOUT) :: PPrim,Z,Tmp
#endif
    TYPE(ARGMT)     , INTENT(IN   ) :: Args
    CHARACTER(LEN=*), INTENT(IN   ) :: Prog
    !-------------------------------------------------------------------
    LOGICAL                         :: IsPresent
    !-------------------------------------------------------------------
    !
    ! IO for the orthogonal PPrim
    CALL Put(PPrim,TrixFile('OrthoDPrime'//TRIM(Args%C%C(4)),Args,1))
    CALL PChkSum(PPrim,'OrthoDPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( PPrim,'OrthoDPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']')
    CALL Plot(   PPrim,'OrthoDPrime'//TRIM(Args%C%C(4))//'_'//TRIM(NxtCycl))
    !
    ! Convert to AO representation
    INQUIRE(FILE=TrixFile('X',Args),EXIST=IsPresent)
    IF(IsPresent)THEN
       CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
       CALL Multiply(Z,PPrim,Tmp)
       CALL Multiply(Tmp,Z,PPrim)
    ELSE
       CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
       CALL Multiply(Z,PPrim,Tmp)
       CALL Get(Z,TrixFile('ZT',Args))
       CALL Multiply(Tmp,Z,PPrim)
    ENDIF
    CALL Filter(Tmp,PPrim)              ! Thresholding
    !
    ! IO for the non-orthogonal PPrim
    CALL Put(Tmp,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,1))
    CALL PChkSum(Tmp,'DPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( Tmp,'DPrime'//TRIM(Args%C%C(4))//'['//TRIM(NxtCycl)//']')
    CALL Plot(   Tmp,'DPrime'//TRIM(Args%C%C(4))//'_'//TRIM(NxtCycl))
    !
  END SUBROUTINE PutXFormPrim
  !
  !
  SUBROUTINE FockPrimGuess(F,FPrim,P,PPrim)
!H---------------------------------------------------------------------------------
!H SUBROUTINE FockPrimGuess(F,FPrim,P,PPrim)
!H  Set the FockPrim guess for TC2R.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrim
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(  OUT) :: P,PPrim
#else
    TYPE(BCSR ), INTENT(  OUT) :: P,PPrim
#endif
    !-------------------------------------------------------------------
    REAL(DOUBLE)             :: Fmin,Fmax,DF,Coeff
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
    CALL Multiply(P,Coeff)       ! P = (I*F_max-F)/DF = X0
    !
    ! Set up PPrim.
    Call SetEq(PPrim,FPrim)
    CALL Add(PPrim,-Fmax)
    CALL Multiply(PPrim,-Coeff)  ! PPrim = -(I*F_max-FPrim)/DF = -X1
    CALL Add(PPrim,P,PPrim)      ! PPrim = PPrim+P = X0-X1
    CALL Multiply(PPrim,-One)    ! PPrim = -PPrim = X1-X0
    !
  END SUBROUTINE FockPrimGuess
  !
  !
  LOGICAL FUNCTION CnvrgChckPrim(Prog,NPur,Ne,MM,T,PPrim,PPrimOld,Tmp1,Tmp2)
!H---------------------------------------------------------------------------------
!H FUNCTION CnvrgChckPrim(Prog,NPur,Ne,MM,T,PPrim,PPrimOld,Tmp1,Tmp2)
!H  Check the convergence during TC2R steps.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL
    TYPE(DBCSR)                     :: T,PPrim,PPrimOld,Tmp1,Tmp2
#else
    TYPE(BCSR )                     :: T,PPrim,PPrimOld,Tmp1,Tmp2
#endif
    REAL(DOUBLE)    , INTENT(IN   ) :: Ne
    CHARACTER(LEN=*), INTENT(IN   ) :: Prog
    INTEGER         , INTENT(IN   ) :: NPur
    INTEGER         , INTENT(INOUT) :: MM
    !-------------------------------------------------------------------
    REAL(DOUBLE)                    :: Prop,AbsErrProp,RelErrProp
    REAL(DOUBLE)                    :: AveErrProp,PNon0
    REAL(DOUBLE)                    :: AbsErrPPrim,FNormErrPPrim,TRacePPrim
    REAL(DOUBLE)    ,SAVE           :: OldProp,OldAEP
    CHARACTER(LEN=2*DEFAULT_CHR_LEN):: Mssg,CnvrgCmmnt
    !-------------------------------------------------------------------
    !
    IF(NPur==0)THEN
       OldProp=BIG_DBL
       OldAEP=BIG_DBL
    ENDIF
    PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
    !
    ! Density matrix errors.
    CALL Multiply(PPrimOld,-One)
    CALL Add(PPrimOld,PPrim,Tmp1)
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
    Prop=-Two*Trace(PPrim,T)
#endif
    !
    ! Compute Abs and Rel error.
    AbsErrProp=ABS(OldProp-Prop)
    RelErrProp=AbsErrProp/ABS(Prop)
    !
    ! Set convergence check.
    CnvrgChckPrim=.FALSE.
    !
    ! Absolute convergence test.
    IF(RelErrProp<Thresholds%ETol*1D-2.AND.AbsErrPPrim<Thresholds%DTol*1D-1)THEN
       CnvrgChckPrim=.TRUE.
       CnvrgCmmnt='Met dE'' goals'
    ENDIF
    !
    ! Test in the asymptotic regime for stall out
    IF(RelErrProp<1D2*Thresholds%Trix**2)THEN
       !    IF(RelErrE<Thresholds%ETol)THEN
       ! Check for increasing /P
       IF(AbsErrPPrim>OldAEP)THEN
          CnvrgChckPrim=.TRUE.
          CnvrgCmmnt='dP'' increase'
       ENDIF
       ! Check for an increasing energy
       IF(Prop>OldProp)THEN
          CnvrgChckPrim=.TRUE.
          CnvrgCmmnt='dE'' increase'
       ENDIF
    ENDIF
    !
    ! Updtate previous cycle values
    OldProp=Prop
    OldAEP=AbsErrPPrim
    !
    ! Print convergence stats
    Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))      &
         //'dE''='//TRIM(DblToShrtChar(RelErrProp))              &
         //', dP''='//TRIM(DblToShrtChar(AbsErrPPrim))           &
         //', %Non0='//TRIM(DblToShrtChar(PNon0))              
    IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
       CALL OpenASCII(OutFile,Out)
       CALL PrintProtectL(Out)
       WRITE(*,*)TRIM(Mssg)
       WRITE(Out,*)TRIM(Mssg)
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
    ENDIF
    !
    ! Set thresholding for next cycle
    CALL SetVarThresh(MM)
    !
    ! Look for convergence
    IF(.NOT.CnvrgChckPrim)THEN
       CALL SetEq(PPrimOld,PPrim)
       RETURN
    ENDIF
    ! Normalize Trace                  ! i do not need that
    ! CALL NormTrace(P,Tmp2,Tmp1,Ne,1) ! i do not need that
    !
    ! Increment the TC2R counter.
    MM=MM+1
    !
    ! Print summary stats.
    TracePPrim = Trace(PPrim)
    IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
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
            &                //', %Non0s = '//TRIM(DblToShrtChar(PNon0))              
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       Mssg=ProcessName(Prog,'Max abs errors') &
            &                //'dE''='//TRIM(DblToShrtChar(AbsErrProp))//', '        &
            &                //'dP''='//TRIM(DblToShrtChar(AbsErrPPrim))
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       Mssg=ProcessName(Prog)//'Rel dE''='//TRIM(DblToShrtChar(RelErrProp))//', '    &
            &                //'||dP''||_F='//TRIM(DblToShrtChar(FNormErrPPrim))
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
    ENDIF
    !
  END FUNCTION CnvrgChckPrim
  !  
  !
END PROGRAM TC2R





