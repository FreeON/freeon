!    Matt Challacombe,  "A simplified density matrix minimization for linear 
!    scaling SCF theory", Journal of Chemical Physics,  110, 2332 (1999) 
!    Major hack by MatCha on 12/10/00:
!    Modified to use the Parser-Manopulis purification scheme.
!    Suss of CG convergence with P-M value c
!---------------------------------------------------------------------
PROGRAM SDMM
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
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)  & 
#else
  TYPE(BCSR)   & 
#endif
                                 :: P,P2,P3,F,G,H,Z,T1,T2,T3
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Fract,Factr,ChemP,Gamma,         &
                                    Numerator,Denominator,Powell,    & 
                                    B,C,D,StepL,Discrim,SqDis,Root1, &
                                    Root2,Ar1,Ar2,Ar,DeltaN,u,v,w,TrP1,TrP2,TrP3, &
                                    OldDeltaN,OldDeltaP,DeltaP,TrP,Fact,Denom,&
                                    Diff,LShift,ENew,CnvgQ,DeltaNQ,DeltaPQ,MinDeltaPQ,   &
                                    FixedPoint,NumPot,DenPot,OldPoint,NewE,OldE,DeltaE,DeltaEQ, &
                                    OldDeltaEQ,OldDeltaPQ,OldDeltaE,MxCom
  REAL(DOUBLE),PARAMETER         :: ThreshAmp=1.D1
  INTEGER                        :: I,J,ICG,NCG,IPur,NPur,OnExit,MinPQCount,PcntPNon0,CGCycles
  LOGICAL                        :: Present,FixedUVW,ToExit
   CHARACTER(LEN=20)             :: ItAnounce
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=4),PARAMETER     :: Prog='SDMM'
!------------------------------------------------------------------ 
! Start
  CALL StartUp(Args,Prog)
! Look for overide of default number of CG cycles
  CALL OpenASCII(InpFile,Inp)  
  IF(.NOT.OptIntQ(Inp,Prog,CGCycles))CGCycles=8
  CLOSE(Inp)
#ifdef PARALLEL
!-----------------------------------------------------------------------
! Repartion based on previous matrices structure
!
  IF(Current(1)>0)THEN
     CALL SADD(TrixFile('OrthoD',Args,0))  ! previous P if possible
  ELSE
     CALL SADD(TrixFile('OrthoF',Args,0))  ! current F
  ENDIF
#endif
!-------------------------------------------------------------------------
  CALL New(P)
  CALL New(T1)
  CALL New(F)
  CALL New(T2)
  CALL New(G)
  CALL New(H)
  CALL New(T3)
!-----------------------------------------------------------------------------
!
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
    CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF
!-----------------------------------------------------------------------------
! Compute the 0th order SDMM gradient 
!
  Fract=DBLE(NEl)/DBLE(2*NBasF)
  Factr=Six*(Fract-One)*Fract
  CALL SetEq(G,F)                 ! G=F
  CALL Multiply(G,Factr)          ! G=[6(N_El/N_BasF-1)*(N_El/N_BasF)]*F
  ChemP=-Trace(G)/DBLE(NBasF)     ! mu=-Trace[G]/NBasF
  CALL Add(G,ChemP)               ! G[0]=G[0]+mu*I
  CALL SetEq(H,G)                 ! H[0]=G[0]
!=============================================================================  
!
! CONGUGATE GRADIENT MINIMIZATION OF Tr{(3P^2-2P^3).F}
!
!=============================================================================  
!
  NCG=0
  OnExit=0
  OldPoint=BIG_DBL
  OldE=Zero
  DO ICG=0,CGCycles-1
     NCG=NCG+1
!
#ifdef PARALLEL
     CALL AlignNodes()
     IF(MyId==ROOT)THEN
#endif
!-----------------------------------------------------------------------------
!    LINE MINIMIZATION: compute coeficients
!
!    Can do these more cheaply, use loose thresholding 
     Thresholds%Trix=ThreshAmp*Thresholds%Trix
!
#ifdef PARALLEL
     CALL Multiply(H,G,T1)
     B=-Trace(T1)                         ! B=-Tr{H.G}=6*Tr{H.(I-P).P.F}
#else
     B=-Trace(H,G)                        ! B=-Tr{H.G}=6*Tr{H.(I-P).P.F}
#endif
     CALL Multiply(H,F,T1)                ! T1=H.F
     CALL Filter(T2,T1)                   ! T2=Filter[H.F]
     CALL Multiply(H,T2,T1)               ! T1=H.HF
#ifdef PARALLEL
     CALL Multiply(H,T1,T2)
     D=-Two*Trace(T2)                     ! D=-2*Trace[H.HHF]      
#else
     D=-Two*Trace(H,T1)                   ! D=-2*Trace[H.HHF]      
#endif

     IF(ICG==0)THEN
        Fract=DBLE(NEl)/DBLE(2*NBasF)     ! P[0]=(N_El/2*N_BasF)*I
        C=Three*(One-Two*Fract)*Trace(T1) ! C=3*Trace{(I-2*P).HHF} 
     ELSE
        CALL SetEq(T2,P)                  ! T2=P
        CALL Multiply(T2,-Two)            ! T2=-2*P
        CALL Add(T2,One)                  ! T2=(I-2*P) 
#ifdef PARALLEL
        CALL Multiply(T2,T1,T3)
        C=Three*Trace(T3)                 ! C=3*Trace{(I-2*P).HHF}        
#else
        C=Three*Trace(T2,T1)              ! C=3*Trace{(I-2*P).HHF}        
#endif
     ENDIF         
!-----------------------------------------------------------------------------
!    LINE MINIMIZATION: solve quadratic
!
     IF(C==Zero.AND.D==Zero)THEN
        IF(ICG<5)THEN 
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(Out,*)' Warning: exiting CG loop in SDMM with c=0, d=0 '
           WRITE(*  ,*)' Warning: exiting CG loop in SDMM with c=0, d=0 '
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
!       Reset thresholding before exit
        Thresholds%Trix=Thresholds%Trix/ThreshAmp
        EXIT
     ELSEIF(d==Zero)THEN
        StepL=-Half*B/C
     ELSE
        Discrim=C*C-Three*B*D
        Denom=Three*D
        IF(Discrim<Zero)THEN 
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(Out,*)'Warning: Complex solution in SDMM Line Minimization '
           WRITE(*  ,*)'Warning: Complex solution in SDMM Line Minimization '
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
!          Reset thresholding before exit
           Thresholds%Trix=Thresholds%Trix/ThreshAmp
           EXIT
        ENDIF
        SqDis=DSQRT(Discrim)
        Root1=-(C+SqDis)/Denom
        Root2=(SqDis-c)/Denom
        AR1=B*Root1+C*Root1**2+D*Root1**3
        AR2=B*Root2+C*Root2**2+D*Root2**3
        AR=MIN(AR1,AR2)
        IF(AR1==AR)THEN
           StepL=Root1
        ELSE
           StepL=Root2
        ENDIF
     ENDIF
!----------------------------------------------------------------------------
!    End cheap thresholding 
     Thresholds%Trix=Thresholds%Trix/ThreshAmp
!-----------------------------------------------------------------------------
!    DENSITY UPTDATE
!
     IF(ICG==0)THEN
        Fract=DBLE(NEl)/DBLE(2*NBasF)
        CALL SetEq(T1,H)                ! T1=H[0]
        CALL Multiply(T1,StepL)         ! T1=StepL*H[0]
        CALL Add(T1,Fract)              ! P[n+1,1]=(N_El/2*N_BasF)*I+StepL[0]*H[0]        
     ELSE       
        CALL Multiply(H,StepL/Two)      ! H=StepL*H 
!       Symmetrization to keep [P,F] small 
        CALL XPose(H,T2)
        CALL Add(H,T2,T3)
        CALL SetEq(H,T3)
        CALL Add(H,P,T1)                ! T1=P[N+1,I+1]=P[N+1,I]+StepL[I]*H[I]   
        CALL Multiply(H,One/StepL)      ! H=H/StepL 
     ENDIF
!    Begin expensive thresholding 
     Thresholds%Trix=Thresholds%Trix/ThreshAmp        
     CALL Filter(P,T1)                  ! P=Filter[P[N+1,I+1]]
!    End expensive thresholding 
     Thresholds%Trix=Thresholds%Trix*ThreshAmp        
!       CALL SetEq(P,T1)  
     PcntPNon0=1D2*DBLE(P%NNon0)/DBLE(NBasF*NBasF)        
!-----------------------------------------------------------------------------
!    COMPUTE CONVERGENCE STATS
!
     CALL Multiply(P,P,T2)
     CALL Multiply(T2,-One)
     CALL Add(P,T2,T3)
     DenPot=Trace(T3)
     CALL Multiply(T3,P,T2) 
     NumPot=Trace(T2)
     c=NumPot/DenPot
     FixedPoint=Half-c
!
     CALL Multiply(P,F,T2)              ! T2=P.F    
     CALL Filter(T3,T2)                 ! T3=Filter[P.F]

!    Add symetrization step to maintain small [F,P]
     CALL XPose(T3,T1)
     CALL Multiply(T1,-One)
     CALL ADD(T3,T1,T2)
     MxCom=MAX(T2)
!
     OldE=NewE
     NewE=Trace(T3)                     ! Tr{P.F}
!-----------------------------------------------------------------------------
!    PRINT CONVERGENCE STATS IF REQUESTED
!
     IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
        Mssg=ProcessName(Prog,'CG '//TRIM(IntToChar(NCG))) &
           //'Tr(P.F) = '//TRIM(DblToMedmChar(NewE))       &
           //', %Non0= '//TRIM(IntToChar(PcntPNon0)) &
           //', MAX[P,F]= '//TRIM(DblToShrtChar(MxCom))    
        WRITE(*,*)TRIM(Mssg)
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)TRIM(Mssg)
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
!-----------------------------------------------------------------------------
!    CHECK MINIMIZATION CONVERGENCE
!
!    Detect convergence inflection
     IF(NINT(OldPoint/ABS(OldPoint))/=NINT(FixedPoint/ABS(FixedPoint)) &
        .AND.NCG/=1.AND.OnExit==0)THEN
!       Do a few more cycles for good measure
        OnExit=NCG+3
     ENDIF
!     IF(OnExit==NCG)EXIT
     OldPoint=FixedPoint
!-----------------------------------------------------------------------------
!    GRADIENT EVALUATION
!
     CALL SetEq(T1,P)
     CALL Multiply(T1,Six)                 ! T1=6*P
     CALL Add(T1,-Six)                     ! T1=6*(P-I)
     CALL Multiply(T1,T3,T2)               ! T2=6*(P-I).PF
     ChemP=-Trace(T2)/DBLE(NBasF)          ! mu=-Trace[GradE]/NBasF
     CALL Add(T2,ChemP)                    ! T2 = G[i+1]=G[i+1]+mu*I     
     Denominator=Dot(G,G)                  ! D=(G[i],G[i])
     CALL Multiply(G,-One)                 ! G=-G[i]          
     CALL Add(T2,G,T1)                     ! T1=(G[i+1]-G[i])
     CALL SetEq(G,T2)                      ! G=G[i+1]
     Numerator=Dot(T1,G)                   ! N=(G[i+1]-G[i]).G[i+1]
     Gamma=MAX(Zero,Numerator/Denominator) 
     CALL Multiply(H,Gamma)                ! H=Gamma*H[i]
     CALL Add(G,H,T1)                      ! T1=G[i+1]+Gamma*H[i]
     CALL Filter(H,T1)                     ! H[i+1]=Filter[G[i+1]+Gamma*H[i]]
   ENDDO
!-----------------------------------------------------------------------------
!  Tidy up a bit ...
!
   CALL Delete(G)
   CALL Delete(H) 
   CALL Delete(T3)
!=============================================================================
!
!  PARSER-MANOPULIS PURIFICATION CYCLES
!
!=============================================================================
!
   CALL New(P2)
   CALL New(P3)
   CALL SetEq(T2,P)          
   OldE=1.D20
   OldDeltaN=1.D20
   OldDeltaP=1.D20
   MinDeltaPQ=1.D20
   MinPQCount=1
   NPur=0
   FixedUVW=.FALSE.
   DO J=0,40
      NPur=NPur+1
      CALL Multiply(P,P,T1) 
!      CALL SetEq(P2,T1)   
      CALL Filter(P2,T1)    
      CALL Multiply(P2,P,T1) 
!      CALL SetEq(P3,T1)   
      CALL Filter(P3,T1)    
      TrP1=Trace(p)
      TrP2=Trace(p2)
      TrP3=Trace(P3)
      c=(TrP2-TrP3)/(TrP1-TrP2)
      IF(ABS(c-Half)<1.D-5.OR.FixedUVW)THEN
         c=Half
         u=Zero
         v=Three
         w=-Two
         FixedUVW=.TRUE.
      ELSEIF(c<=Half)THEN
         u=(One-Two*c)/(One-c) 
         v=(One+c)/(One-c) 
         w=-One/(One-c)
      ELSE
         u=Zero
         v=(One+c)/c     
         w=-One/c
      ENDIF
!------------------------------------------------------------------------------
!     ASSEMBLE PURIFIED P
!
      CALL Multiply(P ,u)
      CALL Multiply(P2,v)
      CALL Multiply(P3,w)
      CALL Add(P,P2,T1)
      CALL Add(T1,P3,P2)     !     P[J+1] = u*P[J] + v*P[J].P[J] + w*P[J].P[J].P[J]
      CALL Filter(P,P2)      !     P=Filter[P[N+1,I+1]]
      PcntPNon0=1D2*DBLE(P%NNon0)/DBLE(NBasF*NBasF)        
!-----------------------------------------------------------------------------
!     COMPUTE PURIFICATION STATS
!
      TrP=Two*Trace(P)
      CALL Multiply(T2,-One)  
      CALL Add(T2,P,T1)                          ! T1=P[J+1]-P[J]
      DeltaP=ABS(Max(T1)+1.D-20)                 ! DeltaP=MAX(P[J+1]-P[J])
      DeltaPQ=ABS(DeltaP-OldDeltaP)/DeltaP

#ifdef PARALLEL
      CALL Multiply(P,F,T1)
      NewE=Trace(T1)                             ! E=Tr{P.F}
#else
      NewE=Trace(P,F)
#endif     
      DeltaE=ABS(OldE-NewE)
      DeltaEQ=ABS((NewE-OldE)/NewE)
!-----------------------------------------------------------------------------
!     PURIFICATION CONVERGENCE CHECK 
!
      ToExit=.FALSE.
!     Test when in the asymptotic regime
      IF(NPur>6)THEN
         CALL OpenASCII(OutFile,Out)
         IF(DeltaEQ<Thresholds%ETol*5.D-4)THEN
!           Check for non-decreasing /P
            IF(DeltaP>OldDeltaP)ToExit=.TRUE.
        ENDIF
!       Check for absolute convergence in total energy
        IF(DeltaEQ<1.D-11)ToExit=.TRUE.
        CLOSE(Out)
      ENDIF
!     Updtate previous cycle values
      OldE=NewE
      OldDeltaE=DeltaE
      OldDeltaP=DeltaP
      OldDeltaEQ=DeltaEQ
      OldDeltaPQ=DeltaPQ
      CALL SetEq(T2,P)
!-----------------------------------------------------------------------------
!     PRINT CONVERGENCE INFO IF REQUESTED
!
#ifdef PARALLEL
      IF(MyId==ROOT)THEN
#endif
         IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
            CALL OpenASCII(OutFile,Out)
            CALL PrintProtectL(Out)
            Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))      &
               //'Tr(P.F) = '//TRIM(DblToMedmChar(NewE))               &
               //', %Non0= '//TRIM(IntToChar(PcntPNon0))               &
!               //', c = '//TRIM(DblToShrtChar(c))                     &
               //', MAX(/P) = '//TRIM(DblToShrtChar(DeltaP)) 
            WRITE(*,*)TRIM(Mssg)
            WRITE(Out,*)TRIM(Mssg)
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
         ENDIF
#ifdef PARALLEL
      ENDIF
#endif
      IF(ToExit)EXIT
   ENDDO
!  Clean up a bit ...
   CALL Delete(F)
   CALL Delete(P2)
   CALL Delete(P3)
!-----------------------------------------------------------------------------
!  Report SDMM statistics
!
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
          CALL OpenASCII(OutFile,Out)
          CALL PrintProtectL(Out)
          Mssg=ProcessName(Prog)//TRIM(IntToChar(NCG))//' CG, and ' &
                //TRIM(IntToChar(NPur))//' purification steps taken.'
          IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
             WRITE(*,*)TRIM(Mssg)
          ENDIF
          WRITE(Out,*)TRIM(Mssg)
          Mssg=ProcessName(Prog)//'MAX(/\P) = '//TRIM(DblToShrtChar(DeltaP))//', ' &
                //'|Tr(P)-NEl| = '//TRIM(DblToShrtChar(DeltaN))//' .'
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
!=============================================================================
!
!  TRANSFORMATION TO AN AO REPRESENTATION AND IO
!
!=============================================================================
!  Renormalization
!
   Fact=DBLE(NEl)/TrP
   CALL Multiply(P,Fact)
!-----------------------------------------------------------------------------
!  IO for the orthogonal P 
!
   CALL Put(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)   
   CALL Put(P,TrixFile('OrthoD',Args,1))
   CALL PChkSum(P,'OrthoP['//TRIM(NxtCycl)//']',Prog)
   CALL PPrint( P,'OrthoP['//TRIM(NxtCycl)//']')
   CALL Plot(   P,'OrthoP_'//TRIM(NxtCycl))
!-----------------------------------------------------------------------------
!  Convert to AO representation
!
   CALL Get(Z,TrixFile('Z',Args))
   CALL Multiply(Z,P,T1)
   CALL Get(Z,TrixFile('ZT',Args))
   CALL Multiply(T1,Z,P)
   CALL Filter(T1,P) 
!-----------------------------------------------------------------------------
!  IO for the non-orthogonal P 
!
   CALL Put(T1,TrixFile('D',Args,1))
   CALL Put(Zero,'homolumogap')
   CALL PChkSum(T1,'P['//TRIM(NxtCycl)//']',Prog)
   CALL PPrint(T1,'P['//TRIM(NxtCycl)//']')
   CALL Plot(T1,'P_'//TRIM(NxtCycl))
!-----------------------------------------------------------------------------
!  Tidy up 
!
   CALL Delete(P)
   CALL Delete(T1)
   CALL Delete(Z)
   CALL ShutDown(Prog)   
END PROGRAM

