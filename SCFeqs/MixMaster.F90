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
#ifdef TC2R_EIGENVAL
  INTEGER :: LWORK,Info
  TYPE(DBL_VECT)                 :: EigenV,Work
  TYPE(DBL_RNK2)                 :: PPrimA
#endif
  !-------------------------------------------------------------------
  !
  ! Initial setup.
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
  !
  ! Get Last SCF cycle.
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  !-------------------------------------------------------------------
  ! Allocations.
  !-------------------------------------------------------------------
  !
  CALL New(F       )
  CALL New(T       )
  CALL New(FPrim   )
  CALL New(P       )
  CALL New(Tmp1    )
  CALL New(Tmp2    )
  CALL New(Tmp3    )
  CALL New(PPrim   )
  CALL New(PPrimOld)
  !
  !-------------------------------------------------------------------
  ! F_(n-2)
  CALL Get(FnMns2,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)-1),BCast_O=.FALSE.)
  ! F_(n-1)
  CALL Get(FnMns1,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
  ! First expand P_n=Theta[L*F_(n-2)+(1-L)*F_(n-1)] about L=0
  ! via P_n(L=0) ~ E_n(0)+Tr{ F_(n-2) PPrime(0) }
  CALL FockPrimGuess(FnMns2,FnMns1,P,PPrim)
  CALL SetEq(PPrimOld,PPrim)
  DO I=1,100
     CALL TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
     IF(CnvrgChckPrim(Prog,I,Ne,MM,T,PPrim,PPrimOld,Tmp1,Tmp2)) EXIT
  ENDDO
  e0=
  e0p=
  ! Now switch around and expand P_n=Theta[L*F_(n-2)+(1-L)*F_(n-1)] about L=1
  ! as P_n(L=1) ~ E_n(1)+Tr{ F_(n-1) PPrime(1) }
  CALL FockPrimGuess(FnMns1,FnMns2,P,PPrim)
  CALL SetEq(PPrimOld,PPrim)
  DO I=1,100
     CALL TC2R_DMP(P,PPrim,Tmp1,Tmp2,Tmp3,Ne,MM)
     IF(CnvrgChckPrim(Prog,I,Ne,MM,T,PPrim,PPrimOld,Tmp1,Tmp2)) EXIT
  ENDDO
  e1=
  e1p=
  ! Minimize the cubic E[L]=a+b*L+c*L^2+d*L^3 WRT L
  a=e0
  b=e0p
  c=-3D0*e0-2D0*e0p+3D0*e1-e1p
  d=2D0*e0+e0p-2D0*e1+e1p
  ! Two extrema for a cubic
  LMns=(-c-SQRT(c*c-3*b*d])/(3*d)
  LPls=(-c+SQRT(c*c-3*b*d])/(3*d)
  EMns=a+b*LMns+c*LMns**2+d*LMns**3
  EPls=a+b*LPls+c*LPls**2+d*LPls**3
  IF(EMns<EPls)THEN
     Lambda=LMns
     OneMnL=One-LMns
  ELSE
     Lambda=LPls
     OneMnL=One-LPls
  ENDIF
  !  F(L) = Lambda*F_(n-2)+(1-Lambda)*F_(n-1)
  CALL Multiply(FnMns2,Lambda)
  CALL Multiply(FnMns1,OneMnL)
  ! Tmp1 is the new, mixed Fock matrix that hopefully 
  ! minimizes the energy.        
  CALL Add(FnMns2,FnMns1,Tmp1)
  ! Compute the corresponding density matrix

  MM=0                        
  Ne=Half*DBLE(NEl)    
  ! Guess P from F
  CALL FockGuess(Tmp1,P,Ne,1)
#endif
  CALL SetEq(Pold,P)    
  ! Do SP2 iterations
  DO I=1,100
     CALL SP2(P,Tmp1,Tmp2,Ne,MM)
     IF(CnvrgChck(Prog,I,Ne,MM,F,P,POld,Tmp1,Tmp2))EXIT
  ENDDO
  ! Orthogonal put and xform to AO rep and put
  CALL PutXForm(Prog,Args,P,POld,Tmp1)
  ! Tidy up

  
  !
  CALL ShutDown(Prog)

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
    !IF(MOD(MM,10).EQ.0) THEN
    !IF(MM==10) THEN
    !   write(*,*) 'loading new D',MM
    !   CALL Delete(P)
    !   CALL New(P)
    !   CALL Get(P,TrixFile('OrthoD',Args,LastSCFCycle-Args%I%I(1)),BCast_O=.FALSE.)
    !ENDIF
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
!       IF(MOD(MM,3).NE.0) THEN
!       IF(MM.LT.3) THEN
!          CALL SetEq(PPrim,Tmp3)
!       ELSE
          CALL Filter(PPrim,Tmp3)        !PPrim=Tmp3
!       ENDIF
       !----------------------
       ! X0 = X02;
       CALL Multiply(P,P,Tmp3)        !Tmp3=P*P
       !----------------------
       !CALL SetEq(P,Tmp3)
!       IF(MM.LT.1) THEN
!          CALL SetEq(P,Tmp3)
!       ELSE
          CALL Filter(P,Tmp3)            !P=Tmp3
!       ENDIF
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
!       IF(MOD(MM,3).NE.0) THEN
!       IF(MM.LT.3) THEN
!          CALL SetEq(PPrim,Tmp1)
!       ELSE
          CALL Filter(PPrim,Tmp1)        !D1=Tmp1
!       ENDIF
       !----------------------
       ! X0 = 2*X0-X02;
       CALL Multiply(P,P,Tmp1)        !Tmp1=P*P
       CALL Multiply(Tmp1,-One)       !Tmp1=-Tmp1
       CALL Multiply(P,Two)           !P=2*P
       CALL Add(P,Tmp1,Tmp3)          !Tmp3=P+Tmp1
       !----------------------
       !CALL SetEq(P,Tmp3)
!       IF(MM.LT.1) THEN
!          CALL SetEq(P,Tmp3)
!       ELSE
          CALL Filter(P,Tmp3)            !X0=Tmp3
!       ENDIF
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
!H  This routine can be improved.
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BCSR ), INTENT(INOUT) :: F,FPrim
#ifdef PARALLEL
    TYPE(DBCSR), INTENT(  OUT) :: P,PPrim
    TYPE(DBCSR)                :: Tmp
#else
    TYPE(BCSR ), INTENT(  OUT) :: P,PPrim
    TYPE(BCSR )                :: Tmp
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
    !
    CALL New(Tmp)
    CALL Add(P,PPrim,Tmp)        ! Tmp = PPrim+P = X0-X1
    CALL SetEq(PPrim,Tmp)        ! PPrim = Tmp = X0-X1
    CALL Delete(Tmp)
    !
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
#ifdef PARALLEL
    PNon0=DBLE(Reduce(P%NNon0))
    PNon0=100.D0*PNon0/DBLE(NBasF*NBasF)
#else
    PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
#endif
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
       CnvrgCmmnt='Met dE'' goals'
    ENDIF
    !
    ! Test in the asymptotic regime for stall out
    !IF(RelErrProp<1D2*Thresholds%Trix**2)THEN
    IF(RelErrProp<Thresholds%ETol)THEN
       ! Check for increasing /P
       IF(AbsErrPPrim>OldAEP)THEN
          CnvrgChckPrim=.TRUE.
          CnvrgCmmnt='dP'' increase'
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
    Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))        &
         //'dE''='//TRIM(DblToShrtChar(RelErrProp))              &
         //', dP''='//TRIM(DblToShrtChar(AbsErrPPrim))           &
         //', %Non0='//TRIM(DblToShrtChar(PNon0))              
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
            &                //', %Non0s = '//TRIM(DblToShrtChar(PNon0))              
       IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
          WRITE(*,*)TRIM(Mssg)
       ENDIF
       WRITE(Out,*)TRIM(Mssg)
       Mssg=ProcessName(Prog,'Max abs errors') &
            &                //'dE''='//TRIM(DblToShrtChar(RelErrProp))//', '        &
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





