!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    Based on
!    Matt Challacombe,  "A simplified density matrix minimization for linear 
!    scaling SCF theory", Journal of Chemical Physics,  110, 2332 (1999) 
!
!    Major hack by MatCha on 11/7/00:
!    Carefull handling of sensitive numerics by avoiding over filtration
!    Suss of CG convergence with Tr(P.P-P)
!    All around improved convergence
!
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
                             :: P,F,G,H,Z,T1,T2,T3
  TYPE(ARGMT)                :: Args
  REAL(DOUBLE)               :: Fract,Factr,ChemP,Gamma,         &
                                Numerator,Denominator,Powell,    & 
                                B,C,D,StepL,Discrim,SqDis,Root1, &
                                Root2,Ar1,Ar2,Ar,DeltaN,         &
                                OldDeltaN,OldDeltaP,DeltaP,TrP,Fact,Denom,&
                                Diff,LShift,ENew,CnvgQ,DeltaNQ,DeltaPQ,   &
                                FixedPoint,NumPot,DenPot,OldPoint
  INTEGER                        :: I,J,ICG,NCG,IPur,NPur
  LOGICAL                        :: Present,LevelShift=.FALSE.
  CHARACTER(LEN=2)               :: Cycl,NxtC
  CHARACTER(LEN=20)              :: ItAnounce
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=4),PARAMETER     :: Prog='SDMM'
!------------------------------------------------------------------ 
! Set up ...
!
  CALL StartUp(Args,Prog)
  Cycl=IntToChar(Args%I%I(1))
  NxtC=IntToChar(Args%I%I(1)+1)
#ifdef PARALLEL
!-----------------------------------------------------------------------
! Repartion according to 
!
!  CALL RePart('/users/mchalla/Estane_1.3-21G.Cyc0.OrthoF')
!  CALL SADD('/users/mchalla/Estane_1.3-21G.Cyc0.OrthoF')
! CALL GreedySCF('/users/mchalla/Estane_1.3-21G.Cyc0.OrthoF')

  FFile=TrixFile('OrthoD',Args,-1)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL SADD(TrixFile('OrthoD',Args,-1)) ! (A) previous P 
  ELSE
     CALL SADD(TrixFile('OrthoF',Args,0))  ! (B) orthognl F
  ENDIF
#endif
!-------------------------------------------------------------------------
  CALL New(F)
  CALL New(P)
  CALL New(G)
  CALL New(H)
  CALL New(T1)
  CALL New(T2)
  CALL New(T3)
!-----------------------------------------------------------
!
!
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
!     CALL Get(F,'/zinc/mchalla/Estane_1.3-21G.Cyc0.OrthoF')
!     CALL Get(F,'/users/mchalla/Estane_1.3-21G.Cyc0.OrthoF')
    CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF


!---------------------------------------------
! Compute the 0th order SDMM gradient 
!
  Fract=DBLE(NEl)/DBLE(2*NBasF)
  Factr=Six*(Fract-One)*Fract
  CALL SetEq(G,F)                 ! G=F
  CALL Multiply(G,Factr)          ! G=[6(N_El/N_BasF-1)*(N_El/N_BasF)]*F
  ChemP=-Trace(G)/DBLE(NBasF)     ! mu=-Trace[G]/NBasF
  CALL Add(G,ChemP)               ! G[0]=G[0]+mu*I
  CALL SetEq(H,G)                 ! H[0]=G[0]
!
!=============================================================  
! CONGUGATE GRADIENT MINIMIZATION OF Tr{(3P^2-2P^3).F}
!=============================================================  
!
  NCG=0
  OldPoint=BIG_DBL
!  Thresholds%Trix=1.D-20

!
  DO ICG=0,19
     NCG=NCG+1
!
#ifdef PARALLEL
     CALL AlignNodes()
     IF(MyId==ROOT)THEN
#endif
!--------------------------------------------
!    LINE MINIMIZATION: compute coeficients
!
!    Can do these more cheaply
!     Thresholds%Trix=1.D2*Thresholds%Trix
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
!--------------------------------------------
!    LINE MINIMIZATION: solve quadratic
!
     IF(C==Zero.AND.D==Zero)THEN
        IF(ICG<5)THEN 
           WRITE(*,*)' Warning: exiting CG loop in SDMM with c=0, d=0 '
        ENDIF
        EXIT
     ELSEIF(d==Zero)THEN
        StepL=-Half*B/C
     ELSE
        Discrim=C*C-Three*B*D
        Denom=Three*D
        IF(Discrim<Zero)THEN 
           WRITE(*,*)'Warning: Complex solution in SDMM Line Minimization '
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
!    Reset thresholds
!     Thresholds%Trix=1.D-2*Thresholds%Trix
!------------------------------------------------
!    DENSITY UPTDATE
!
     IF(ICG==0)THEN
        Fract=DBLE(NEl)/DBLE(2*NBasF)
        CALL SetEq(T1,H)                ! T1=H[0]
        CALL Multiply(T1,StepL)         ! T1=StepL*H[0]
        CALL Add(T1,Fract)              ! P[n+1,1]=(N_El/2*N_BasF)*I+StepL[0]*H[0]        
     ELSE       
        CALL Multiply(H,StepL)          ! H=StepL*H 
        CALL Add(H,P,T1)                ! T1=P[N+1,I+1]=P[N+1,I]+StepL[I]*H[I]   
        CALL Multiply(H,One/StepL)      ! H=H/StepL 
     ENDIF
!    This filter leads to instability in early cycles
     IF(NCG>3)THEN
       CALL Filter(P,T1)                ! P=Filter[P[N+1,I+1]]
     ELSE
       CALL SetEq(P,T1)                 ! P=P[N+1,I+1]
     ENDIF
!------------------------------------------------
!    COMPUTE CONVERGENCE STATS
!
     CALL Multiply(P,P,T2)
     CALL Multiply(T2,-One)
     CALL Add(P,T2,T3)
     DenPot=Trace(T3)
     CALL Multiply(T3,P,T2) 
     NumPot=Trace(T2)
     FixedPoint=NumPot/DenPot
     FixedPoint=Half-FixedPoint     
!
     CALL Multiply(P,F,T2)              ! T2=P.F    
     CALL Filter(T3,T2)                 ! T3=Filter[P.F]
!
     ENew=Trace(T3)                     ! Tr{P.F}
!------------------------------------------------
!    PRINT CONVERGENCE STATS IF REQUESTED
!
     IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        Mssg=ProcessName(Prog,'CG '//TRIM(IntToChar(NCG))) &
           //'FixedPoint = '//TRIM(DblToShrtChar(FixedPoint))     &
           //', Tr(P.F) = '//TRIM(DblToMedmChar(ENew)) 
        WRITE(*,*)TRIM(Mssg)
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)TRIM(Mssg)
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
!------------------------------------------------
!    CHECK FOR MINIMIZATION CONVERGENCE
!
!    Convergence inflection
!          
     IF(NINT(OldPoint/ABS(OldPoint))/=NINT(FixedPoint/ABS(FixedPoint)).AND.NCG/=1)EXIT
     OldPoint=FixedPoint
!--------------------------------------------------------
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
!
   ENDDO
!-------------------------
!  Tidy up a bit ...
!
   CALL Delete(G)
   CALL Delete(H) 
   CALL Delete(T3)
!
!=============================================================
!  PURIFICATION CYCLES
!=============================================================
!
   NPur=0
   CALL SetEq(T1,P)          ! T1=P
   OldDeltaN=BIG_DBL
   OldDeltaP=BIG_DBL
!
   DO J=0,20
      NPur=NPur+1
!--------------------------------------------
!     P[J+1]=(2*P.P.P-3*P.P)[J]
!
      CALL Multiply(T1,-Two) ! T1=-2P[J]
      CALL Add(T1,Three)     ! T1=(3*I-2*P[J])
      CALL Multiply(T1,P,T2) ! T2=(3*I-2*P[J]).P[J]
      CALL Filter(T1,T2)     ! T1=Filter[(3*I-2*P[J]).P[J]]     
      CALL Multiply(T1,P,T2) ! T2=P[J+1]=(3*I-2*P[J]).P[J].P[J] 
      CALL Filter(T1,T2)     ! T1=Filter[(3*I-2*P[J]).P[J].P[J]]
!--------------------------------------------------------------------
!     COMPUTE PURIFICATION STATS
!
      TrP=Trace(T1)                              ! N[J+1]=Tr{P[J+1]}
      CALL Multiply(P,-One)                      ! P=-P[J]
      CALL Add(T1,P,T2)                          ! T2=P[J+1]-P[J]
      DeltaP=Max(T2)+1.D-20                      ! DeltaP=MAX(P[J+1]-P[J])
      DeltaN=DABS(Two*TrP-DBLE(NEl))+1.D-20      ! DeltaN=2*Tr{P}-N_El
      DeltaPQ=ABS(DeltaP-OldDeltaP)/DeltaP
      DeltaNQ=ABS(DeltaN-OldDeltaN)/DeltaN
!
      IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
#ifdef PARALLEL
         CALL Multiply(P,F,T3)
         ENew=Trace(T3)                          ! Tr{P.F}
#else
         ENew=Trace(P,F)
#endif     
      ENDIF
!------------------------------------------------------------------
!     PRINT CONVERGENCE INFO IF REQUESTED
!
#ifdef PARALLEL
      IF(MyId==ROOT)THEN
#endif
         IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
            CALL OpenASCII(OutFile,Out)
            CALL PrintProtectL(Out)
            Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur))) &
             //'    Tr(P)-Nel = '//TRIM(DblToShrtChar(DeltaN))         &
             //', Tr(P.F) = '//TRIM(DblToMedmChar(ENew))
            WRITE(*,*)TRIM(Mssg)
            WRITE(Out,*)TRIM(Mssg)
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
         ENDIF
#ifdef PARALLEL
      ENDIF
#endif
!--------------------------------------------------------------------
!     CHECK CONVERGENCE
!
!     Test only in asymptotic regime
!      WRITE(*,*)' DeltaN = ',DeltaN,' DeltaP = ',DeltaP
      IF(DeltaN<1.0D-2.AND.DetlaP<1.D-1.AND.NPur>3)THEN
!        Check for low digit rebound
         IF(DeltaP>OldDeltaP)THEN
!            WRITE(*,*)' exit 1 '
!            WRITE(*,*)' DeltaPs = ',DeltaP,OldDeltaP
           EXIT
         ENDIF
!        Check for stagnation stallout
!         WRITE(*,*)' Delta NQ = ',DeltaNQ,' DeltaPQ = ',DeltaPQ
         IF(DeltaNQ<2.D-1.OR.DeltaPQ<2.D-1)THEN
!            WRITE(*,*)' exit 2'
            EXIT
         ENDIF
      ENDIF  
      OldDeltaP=DeltaP
      OldDeltaN=DeltaN
      CALL SetEq(P,T1)                           ! P=P[J+1]
!
   ENDDO
!  Check for failure
   IF(DeltaN>1.0D0)CALL Warn(' Convergence failure in SDMM, lost /N = ' &
                            //TRIM(DblToMedmChar(DeltaN))//' electrons.')
!  Clean up a bit ...
   CALL Delete(F)
   CALL Delete(P)
!-----------------------------------------------
!  Report SDMM statistics
!
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       CALL OpenASCII(OutFile,Out)
       CALL PrintProtectL(Out)
       Mssg=ProcessName(Prog)//TRIM(IntToChar(NCG))//' CG, and ' &
             //TRIM(IntToChar(NPur))//' purification steps taken.'
       WRITE(*,*)TRIM(Mssg)
       WRITE(Out,*)TRIM(Mssg)
       Mssg=ProcessName(Prog)//'MAX(/\P) = '//TRIM(DblToShrtChar(DeltaP))//', ' &
             //'|TrP-NEl| = '//TRIM(DblToShrtChar(DeltaN))//' .'
       WRITE(*,*)TRIM(Mssg)
       WRITE(Out,*)TRIM(Mssg)
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL
   ENDIF
#endif
!
!--------------------------------------------------------------------
!  Renormalization
!         
   Fact=DBLE(NEl)/(Two*TrP)
   CALL Multiply(T1,Fact)
!--------------------------------------------------------------------
!  IO for the orthogonal P 
!
   CALL Put(T1,TrixFile('OrthoD',Args,1))
   CALL PChkSum(T1,'OrthoP['//TRIM(NxtC)//']',Prog)
   CALL PPrint( T1,'OrthoP['//TRIM(NxtC)//']')
   CALL Plot(   T1,'OrthoP_'//TRIM(NxtC))
!---------------------------------------------
!  Convert back to an AO representation
!
   CALL Get(Z,TrixFile('Z',Args))
   CALL Multiply(Z,T1,T2)
   CALL Get(Z,TrixFile('ZT',Args))
   CALL Multiply(T2,Z,T1)
   CALL Filter(T2,T1) 
!--------------------------------------------------------------------
!  IO for the non-orthogonal P 
!
   CALL Put(T2,TrixFile('D',Args,1),      &
            BlksName_O='ndi'//TRIM(NxtC), &
            Non0Name_O='ndm'//TRIM(NxtC)  )
   IF(PrintFlags%Key>DEBUG_MEDIUM) &
      CALL PChkSum(T2,'P['//TRIM(NxtC)//']',Prog)
   IF(PrintFlags%Mat==DEBUG_MATRICES)THEN
      CALL PPrint(T2,'P['//TRIM(NxtC)//']')
   ELSEIF(PrintFlags%Mat==PLOT_MATRICES)THEN
      CALL Plot(T2,'P_'//TRIM(NxtC))
   ENDIF
!------------------------------
!
   CALL Put(0.123456D0,'homolumogap')
!------------------------------
!  Tidy up 
!
   CALL Delete(T1)
   CALL Delete(T2)
   CALL Delete(Z)
!
   CALL ShutDown(Prog)   
!----------------------------------------------------------------------

!
END PROGRAM ! SDMM

