MODULE DenMatMethods
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
  IMPLICIT NONE
  REAL(DOUBLE)            :: TrP,TrP2,TrP3,TrP4
CONTAINS
!--------------------------------------------------------------
!
!--------------------------------------------------------------
  SUBROUTINE SussTrix(TrixName,Prog)
    CHARACTER(LEN=*) :: TrixName,Prog
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
    REAL(DOUBLE)     :: Trix
!--------------------------------------------------------------
    ! Check for thresholds overide 
    CALL OpenASCII(InpFile,Inp)         
    IF(OptDblQ(Inp,TrixName,Trix))THEN
       Thresholds%Trix=Trix
       Mssg=TRIM(ProcessName(Prog))//' Trix = '  &
          //TRIM(DblToShrtChar(Thresholds%Trix))
       CALL OpenASCII(OutFile,Out)         
       WRITE(Out,*)TRIM(Mssg)
       CLOSE(Out)
    ENDIF
    CLOSE(Inp)
  END SUBROUTINE SussTrix
!--------------------------------------------------------------
!
!--------------------------------------------------------------
  SUBROUTINE PutXForm(Prog,Args,P,Z,Tmp1)
    CHARACTER(LEN=*) :: Prog
    TYPE(ARGMT)      :: Args
    TYPE(BCSR)       :: P,Z,Tmp1
    LOGICAL          :: Present
!----------------------------------------------------------
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
       CALL Multiply(Z,P,Tmp1)
       CALL Multiply(Tmp1,Z,P)
    ELSE
       CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
       CALL Multiply(Z,P,Tmp1)
       CALL Get(Z,TrixFile('ZT',Args))
       CALL Multiply(Tmp1,Z,P)
    ENDIF
    CALL Filter(Tmp1,P)     ! Thresholding
    ! IO for the non-orthogonal P
    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,1))
    CALL Put(Zero,'homolumogap')
    CALL PChkSum(Tmp1,'P['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint(Tmp1,'P['//TRIM(NxtCycl)//']')
    CALL Plot(Tmp1,'P_'//TRIM(NxtCycl))
  END SUBROUTINE PutXForm
!--------------------------------------------------------------
!
!--------------------------------------------------------------
  SUBROUTINE SetVarThresh(MM)
    INTEGER :: MM
    REAL(DOUBLE),SAVE :: OldThresh=0D0
    RETURN
    IF(OldThresh==0D0)THEN
       OldThresh=Thresholds%Trix
       Thresholds%Trix=1D-2*Thresholds%Trix
    ElSE
       Thresholds%Trix=MIN(15D-1**(MM)*Thresholds%Trix,OldThresh)
    ENDIF
  END SUBROUTINE SetVarThresh
!--------------------------------------------------------------
!
!--------------------------------------------------------------
  FUNCTION CnvrgChck(Prog,NPur,Ne,MM,F,P,POld,Tmp1,Tmp2)
    LOGICAL              :: CnvrgChck
    TYPE(BCSR)           :: F,P,POld,Tmp1,Tmp2
    REAL(DOUBLE)         :: Ne,Energy,AbsErrP,AveErrP,L2ErrP,TwoNP,N2F,  &
         AbsErrE,RelErrE,AveErrE,MaxCommErr,L2CommErr
    REAL(DOUBLE),SAVE    :: OldE,OldAEP
    INTEGER              :: MM,PNon0,NPur
    CHARACTER(LEN=*)     :: Prog
    CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg,CnvrgCmmnt
!---------------------------------------------------------------------
     IF(NPur==0)THEN
        OldE=BIG_DBL
        OldAEP=BIG_DBL
     ENDIF
     PNon0=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
     ! Density matrix errors
     CALL Multiply(Pold,-One)
     CALL Add(Pold,P,Tmp1)
     AbsErrP=ABS(Max(Tmp1)+1.D-20)
     AveErrP=OneNorm(Tmp1)/DBLE(PNon0)
     L2ErrP=TwoNorm(Tmp1)/TwoNorm(P)
     ! Energy errors
#ifdef PARALLEL
     CALL Multiply(P,F,Tmp1)
     Energy=Trace(Tmp1)    
#else
     Energy=Trace(P,F)
#endif     
     AbsErrE=ABS(OldE-Energy)
     RelErrE=AbsErrE/ABS(Energy)
     ! Convergence check
     CnvrgChck=.FALSE.
     ! Absolute convergence test
     IF(RelErrE<Thresholds%ETol*1D-2.AND. &
        AbsErrP<Thresholds%DTol*1D-1)THEN
        CnvrgChck=.TRUE.
        CnvrgCmmnt='Met dE/dP goals'
     ENDIF
     ! Test in the asymptotic regime for stall out
     IF(RelErrE<Thresholds%ETol)THEN
        ! Check for increasing /P
        IF(AbsErrP>OldAEP)THEN
           CnvrgChck=.TRUE.
           CnvrgCmmnt='Hit dP increase'
        ENDIF
        ! Check for an increasing energy
        IF(Energy>OldE)THEN
           CnvrgChck=.TRUE.
           CnvrgCmmnt='Hit dE increase'
        ENDIF
     ENDIF
!     IF(NPur<99)CnvrgChck=.FALSE.

     ! Updtate previous cycle values
     OldE=Energy
     OldAEP=AbsErrP
     ! Print convergence stats
     Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))      &
          //'dE='//TRIM(DblToShrtChar(RelErrE))                 &
          //', dP='//TRIM(DblToShrtChar(AbsErrP))               &
          //', ThrX='//TRIM(DblToShrtChar(Thresholds%Trix))     &
          //', %Non0='//TRIM(IntToChar(PNon0))              
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(*,*)TRIM(Mssg)
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
#ifdef PARALLEL
     ENDIF
#endif
     IF(.NOT.CnvrgChck)THEN
        CALL SetEq(Pold,P)
        RETURN
     ENDIF
     ! Normalize Trace
     CALL NormTrace(P,Tmp2,Tmp1,Ne,1)
     MM=MM+1
     ! Commutator [F,P] 
     N2F=TwoNorm(F)
     CALL Multiply(F,P,Tmp1)
     CALL Multiply(P,F,POld)
     CALL Multiply(POld,-One)
     CALL Add(Tmp1,POld,F) 
     MM=MM+2
     L2CommErr=TwoNorm(F)/N2F
     MaxCommErr=Max(F)
     ! Print summary stats
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)

           Mssg=ProcessName(Prog,CnvrgCmmnt) &
                //'Tr{FP}='//TRIM(DblToChar(Energy)) &
                //', dNel = '//TRIM(DblToShrtChar(Two*ABS(Trace(P)-Ne))) 
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)

           Mssg=ProcessName(Prog)//TRIM(IntToChar(NPur))//' purification steps, ' &
                //TRIM(IntToChar(MM))//' matrix multiplies, %Non0s = '//TRIM(IntToChar(PNon0))              
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog,'Max abs errors') &
                //'dE='//TRIM(DblToShrtChar(AbsErrE))//', '                 &
                //'dP='//TRIM(DblToShrtChar(AbsErrP))//', '                 &
                //'[F,P]='//TRIM(DblToShrtChar(MaxCommErr))
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           Mssg=ProcessName(Prog,'Rel L2  errors') &
                //'dE='//TRIM(DblToShrtChar(RelErrE))//', '                &
                //'dP='//TRIM(DblToShrtChar(L2ErrP))//', '                 &
                //'[F,P]='//TRIM(DblToShrtChar(L2CommErr))
           IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              WRITE(*,*)TRIM(Mssg)
           ENDIF
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
#ifdef PARALLEL
     ENDIF
#endif
   END FUNCTION CnvrgChck
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE SP2(P,P2,Tmp1,Norm,MMMs)
      TYPE(BCSR)   :: P,P2,Tmp1
      REAL(DOUBLE) :: Norm,CR
      INTEGER      :: MMMs
!----------------------------------------------------------------------------
      CALL Multiply(P,P,P2)             ! The only multiplication is a square
      MMMs=MMMs+1
      TrP=Trace(P)
      CR=TrP-Norm                       ! CR = Occupation error criteria
      IF (CR >  0) THEN                 ! Too many states
         CALL Filter(P,P2)              ! P = P^2
      ELSE                              ! Too few states 
         CALL Multiply(P ,Two) 
         CALL Multiply(P2,-One)
         CALL Add(P,P2,Tmp1)           ! P = 2P-P^2
         CALL Filter(P,Tmp1)
      ENDIF
    END SUBROUTINE SP2
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
    SUBROUTINE SP4(P,P2,Tmp1,Tmp2,Norm,MMMs)
      TYPE(BCSR)                     :: P,P2,Tmp1,Tmp2
      REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Thresh_old
      INTEGER                        :: MMMs
      REAL(DOUBLE)                   :: EPS = 1.D-12
!----------------------------------------------------------------------------
      CALL Multiply(P,P,P2)
      MMMs = MMMs+1
      TrP  = Trace(P)
      TrP2 = Trace(P2)
      TrP3 = Trace(P,P2)
      TrP4 = Trace(P2,P2)
      Gt   = (Norm-Four*TrP3+3*TrP4)
      Gn   = (TrP2-Two*TrP3+TrP4)
      CR   = TrP - Norm                   ! Trace correction measure
      IF(ABS(Gn).LT.EPS)THEN              ! Close to idempotency
         G=Three                          ! To avoid numerical errors
      ELSE
         G=Gt/Gn                          ! Boundary measure
      ENDIF
      IF((G>Zero).AND.(G<Six))THEN        ! Check the bounds
         IF (CR  > Zero) THEN             ! Too many states
           CALL Filter(Tmp2,P2)
           CALL Multiply(P,4.d0)
           CALL Multiply(P2,-3.d0)
           CALL Add(P,P2,Tmp1)
           CALL Filter(P,Tmp1)
           CALL Multiply(Tmp2,P,Tmp1)     ! P = P^2*(4P-3P^2)
           MMMs = MMMs+1
        ELSE                              ! Too few states
           CALL Filter(Tmp2,P2)
           CALL Multiply(P ,-8.d0)
           CALL Multiply(P2, 3.d0)
           CALL Add(P,P2,Tmp1)
           CALL Add(Tmp1,Six)
           CALL Filter(P,Tmp1)
           CALL Multiply(Tmp2,P,Tmp1)     ! P = P^2*(6I-8P+3P^2)
           MMMs = MMMs+1
        ENDIF
        CALL Filter(P,Tmp1)

     ELSE 
        IF(G<0) THEN                      ! Too many states
           CALL Filter(P,P2)              ! P = P^2
        ELSE                              ! Too few states 
           CALL Multiply(P ,Two) 
           CALL Multiply(P2,-One)
           CALL Add(P ,P2,Tmp1)           ! P = 2P-P^2
           CALL Filter(P,Tmp1)
        ENDIF
     ENDIF
   END SUBROUTINE SP4
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
   SUBROUTINE NT4(P,P2,Tmp1,Tmp2,Norm,MMMs)
     TYPE(BCSR)                     :: P,P2,Tmp1,Tmp2
     REAL(DOUBLE)                   :: CR,Norm,Gt,Gn,G,Coeff
     INTEGER                        :: MMMs
     REAL(DOUBLE)                   :: EPS = 1.D-12
!----------------------------------------------------------------------------
     CALL Multiply(P,P,P2)    
     MMMs = MMMs+1
     TrP  = Trace(P)
     TrP2 = Trace(P2)
     TrP3 = Trace(P,P2)
     TrP4 = Trace(P2,P2)
     Gt   = (Norm-Four*TrP3+3*TrP4)
     Gn   = (TrP2-Two*TrP3+TrP4)
     CR   = TrP - Norm     
     IF(ABS(Gn).LT.EPS)THEN 	              ! Close to idempotency
        G=Three 	        	      ! To avoid numerical errors
     ELSE
        G=Gt/Gn	                       	      ! Boundary measure
     ENDIF
     IF ((G>Zero).AND.(G<Six))THEN 	      ! Check the bounds
        CALL Filter(Tmp2,P2)
        Coeff = Four-Two*G
        CALL Multiply(P,Coeff)
        Coeff = G-Three
        CALL Multiply(P2,Coeff)
        CALL Add(P,P2,Tmp1)
        CALL Add(Tmp1,G)
        CALL Filter(P,Tmp1)
        CALL Multiply(Tmp2,P,Tmp1)         ! P^2*(G+(4-2*G)*P+(G-3)*P^2)
        MMMs = MMMs+1
     ELSE
        IF (G < 0) THEN                      ! Too many states
           CALL SetEq(Tmp1,P2)              ! P = P^2
        ELSE                                 ! Too few states 
           CALL Multiply(P ,Two) 
           CALL Multiply(P2,-One)
           CALL Add(P,P2,Tmp1)             ! P = 2P-P^2
        ENDIF
     ENDIF
     CALL Filter(P,Tmp1)
   END SUBROUTINE NT4
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
   SUBROUTINE PM1(P,P2,P3,Tmp1,Tmp2,Norm,MMMs)
     TYPE(BCSR)                     :: P,P2,P3,Tmp1,Tmp2
     REAL(DOUBLE)                   :: CR,Norm,Coeff
     INTEGER                        :: MMMs
     REAL(DOUBLE)                   :: EPS = 1.D-8
!----------------------------------------------------------------------------
     CALL Multiply(P,P,P2)                   
     MMMs=MMMs+1
     CALL Filter(Tmp1,P2)
     CALL Multiply(P,Tmp1,P3)               
     MMMs=MMMs+1
     TrP=Trace(P)
     TrP2=Trace(P2)
     TrP3=Trace(P3)
     IF(ABS(TrP-TrP2)<EPS)THEN                ! Close to idempotency
       CALL Multiply(P2, Three)
       CALL Multiply(P3,-Two)
       CALL Add(P2,P3,Tmp1)                   ! T = 3P^2-2P^2 (McWeeny close to convergence)
     ELSE
        CR = (TrP2 - TrP3)/(TrP-TrP2)
        IF (CR <= Half) THEN
           Coeff = (One-Two*CR)/(One-CR)
           CALL Multiply(P,Coeff)
           Coeff = (One+CR)/(One-CR)
           CALL Multiply(P2,Coeff)
           Coeff = -One/(One-CR)
           CALL Multiply(P3,Coeff)
           CALL Add(P,P2,Tmp2)
           CALL Add(Tmp2,P3,Tmp1)             ! T = [(1-2CR)*P+(1+CR)*P^2-P^3]/(1-CR)
        ELSE
           Coeff = (One+CR)/CR
           CALL Multiply(P2,Coeff)
           Coeff = -One/CR
           CALL Multiply(P3,Coeff)
           CALL Add(P2,P3,Tmp1)               ! T = [(1+CR)*P^2-P^3]/CR
        ENDIF
     ENDIF
     CALL Filter(P,Tmp1)
   END SUBROUTINE PM1

   SUBROUTINE PM2(P,P2,P3,Tmp1,MMMs)
     TYPE(BCSR)    :: P,P2,P3,Tmp1
     REAL(DOUBLE)  :: c,u,v,w,TrP1
     INTEGER       :: MMMs
     LOGICAL,SAVE  :: FixedUVW=.FALSE.
!----------------------------------------------------------------------------
      CALL Multiply(P,P,Tmp1) 
      MMMs=MMMs+1
      CALL Filter(P2,Tmp1)    
      CALL Multiply(P2,P,Tmp1) 
      MMMs=MMMs+1
      CALL Filter(P3,Tmp1)    
      TrP1=Trace(P)
      TrP2=Trace(P2)
      TrP3=Trace(P3)
      c=(TrP2-TrP3)/(TrP1-TrP2+1D-30)
      IF(ABS(c-Half)<1.D-5.OR.FixedUVW)THEN
         c=Half
         u=Zero
         v=Three
         w=-Two
         FixedUVW=.TRUE.
      ELSEIF(c<=Half)THEN
         u=(One-Two*c)/(One-c+1D-30) 
         v=(One+c)/(One-c+1D-30) 
         w=-One/(One-c+1D-30)
      ELSE
         u=Zero
         v=(One+c)/(c+1D-30)     
         w=-One/(c+1D-30)
      ENDIF
!     Assemble purified P
      CALL Multiply(P ,u)
      CALL Multiply(P2,v)
      CALL Multiply(P3,w)
      CALL Add(P,P2,Tmp1)
      CALL Add(Tmp1,P3,P2) !     P[J+1] = u*P[J] + v*P[J].P[J] + w*P[J].P[J].P[J]
      CALL Filter(P,P2)    !     P=Filter[P[N+1,I+1]]
   END SUBROUTINE PM2
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
    SUBROUTINE FockGuess(F,P,Norm,Order)
      TYPE(BCSR)                     :: F,P
      REAL(DOUBLE)                   :: Fmin,Fmax,DF,Coeff,Mu,Lmbd1,  &
                                        Lmbd2,Norm
      INTEGER                        :: Order
!----------------------------------------------------------------------------
!     Estimate spectral bounds 
!----------------------------------------------------------------------------
      CALL SpectralBounds(F,Fmin,Fmax)
!----------------------------------------------------------------------------
!     Set up the Density Matrix
!----------------------------------------------------------------------------
      IF(Order==1) THEN
         Call SetEq(P,F)
         DF = (Fmax - Fmin)
         CALL Add(P,-Fmax)
         Coeff = -One/DF
         CALL Multiply(P,Coeff)          ! P = (I*F_max-F)/DF
      ELSEIF(Order==2) THEN
         Mu    = Trace(F)/DBLE(NBasF)
         Lmbd1 = Norm/(Fmax-Mu)
         Lmbd2 = (DBLE(NBasF)-Norm)/(Mu-Fmin)
         IF(Lmbd1 < Lmbd2) THEN
            CALL SetEq(P,F)
            CALL Add(P,-Mu)
            Coeff = -Lmbd1/DBLE(NBasF)
            CALL Multiply(P,Coeff) 
            Coeff = Norm/DBLE(NBasF)
            CALL Add(P,Coeff) 
         ELSE
            CALL SetEq(P,F)
            CALL Add(P,-Mu)
            Coeff = -Lmbd2/DBLE(NBasF)
            CALL Multiply(P,Coeff) 
            Coeff = Norm/DBLE(NBasF)
            CALL Add(P,Coeff)
         ENDIF
      ELSE
         CALL MondoHalt(99,'Wrong Order in FockGuess')
      ENDIF
    END SUBROUTINE FockGuess
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
   SUBROUTINE NormTrace(P,P2,Tmp1,Norm,Order)
     TYPE(BCSR)    :: P,P2,Tmp1
     REAL(DOUBLE)  :: Norm,CR,C1,C2
     INTEGER       :: Order
!----------------------------------------------------------------------------
     IF(Order==0) THEN
        TrP  = Trace(P)
        CR   = Norm/TrP
        CALL Multiply(P,CR)
     ELSEIF(Order==1) THEN
        CALL Multiply(P,P,P2)
        TrP  = Trace(P)
        TrP2 = Trace(P2)
        CR   = TrP - TrP2
        IF(CR==Zero) RETURN
        C1   = (Norm - TrP2)/CR
        C2   = (TrP- Norm)/CR
        CALL Multiply(P  ,C1)
        CALL Multiply(P2 ,C2)
        CALL Add(P,P2,Tmp1)
        CALL Filter(P,Tmp1)
     ELSE
        CALL MondoHalt(99,'Wrong Order in NormTrace')
     ENDIF
   END SUBROUTINE NormTrace
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
   SUBROUTINE CalculateDegen(Norm,Degen,Occpan)
     REAL(DOUBLE)                 :: Norm,Degen,Occpan
!
     Degen     = ((Norm-TrP2)**3)/((TrP2-TrP3)*(Norm-Two*TrP2+TrP3))
     Occpan    = Two*(TrP2-TrP3)/(Norm-TrP2)  
!
   END SUBROUTINE CalculateDegen
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
#ifdef BROKEN_GAP
   SUBROUTINE CalculateGap(F,P,Norm,lumo_occ,Gap)
     TYPE(BCSR)                   :: F,P
     REAL(DOUBLE)                 :: Norm,Gap,lumo_occ
     REAL(DOUBLE)                 :: TrFP,TrFP2,TrFP3
!
     CALL Multiply(P,P,P2)
     CALL Multiply(P,P2,P3)
     TrFP  = Trace(F,P)
     TrFP2 = Trace(F,P2)
     TrFP3 = Trace(F,P3)    
! 
     lumo_occ = Half + Half*Sqrt(ABS(One-Two*Norm+Two*TrP2)) 
!
     Gap      = -(TrFP - Three*TrFP2 + Two*TrFP2)/(lumo_occ*(lumo_occ-One)*(Two*lumo_occ-One))
!
   END SUBROUTINE CalculateGap
#endif
END MODULE DenMatMethods
