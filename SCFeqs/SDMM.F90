!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    Matt Challacombe,  "A simplified density matrix minimization for linear 
!    scaling SCF theory", Journal of Chemical Physics,  110, 2332 (1999) 
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
                                DeltaN_Old,DeltaP,TrP,Fact,Denom,&
                                Diff,LShift
  INTEGER                        :: I,J,ICG,NCG,IPur,NPur
  LOGICAL                        :: Present,LevelShift=.FALSE.
  CHARACTER(LEN=2)               :: Cycl,NxtC
  CHARACTER(LEN=4),PARAMETER     :: Prog='SDMM'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
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
!
!
!  IF(Args%NI>=2)THEN
!     NCG=Args%I%I(2)
!  ELSE
     NCG=20
!  ENDIF
!-----------------------------------------
! Allocations 
!
  CALL New(F)
  CALL New(P)
  CALL New(G)
  CALL New(H)
  CALL New(T1)
  CALL New(T2)
#ifdef PARALLEL
  CALL New(T3)
#endif
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
  ChemP=Trace(G)/DBLE(NBasF)      ! mu=Trace[G]/NBasF
  CALL Add(G,ChemP)               ! G[0]=G[0]-mu*I
  CALL SetEq(H,G)                 ! H[0]=G[0]
!=============================================================  
! CONGUGATE GRADIENT MINIMIZATION
!
  DO ICG=0,NCG
#ifdef PARALLEL
     CALL AlignNodes()
     IF(MyId==ROOT)THEN
#endif
        Mssg='SDMM: CG Step # '//TRIM(IntToChar(ICG))
        WRITE(*,*)TRIM(Mssg)
        IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(Out,*)TRIM(Mssg)
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
#ifdef PARALLEL
     ENDIF
!     CALL AlignNodes()
!      CALL PPrint(PerfMon,Prog)
!     IF(ICG==1)CALL HALT(Prog)   
!     IF(ICG==2)CALL ShutDown(Prog)   
#endif
!--------------------------------------------
!    LINE MINIMIZATION: compute coeficients
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
!#ifdef PARALLEL
!     IF(MyId==ROOT)  &
!#endif
!     WRITE(*,*)' I = ',ICG,' B = ',B,' C = ',C,'D = ',D
!--------------------------------------------
!    LINE MINIMIZATION: solve quadratic
!
     IF(C==Zero.AND.D==Zero)THEN
        IF(ICG<5)THEN 
           WRITE(*,*)' Warning: exiting CG loop in SDMM '
           WRITE(*,*)' Warning: with ICG = ',ICG
        ENDIF
        EXIT
     ELSEIF(d==Zero)THEN
        StepL=-Half*B/C
     ELSE
        Discrim=C*C-Three*B*D
        Denom=Three*D
        IF(Discrim<Zero)  &
           CALL Halt(' Complex solution in SDMM Line Minimization ')
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
     CALL Filter(P,T1)                  ! P=Filter[P[N+1,I+1]]
!--------------------------------------------------------
!    GRADIENT EVALUATION
!
     CALL Multiply(T1,Six)              ! T1=6*P
     CALL Add(T1,-Six)                  ! T1=6*(P-I)
     CALL Multiply(T1,P,T2)             ! T2=6*(P-I).P
     CALL Filter(T1,T2)                 ! T1=Filter[6*(P-I).P]    
     CALL Multiply(T1,F,T2)             ! T2=G[i+1]=6*(P-I).P.F
     ChemP=Trace(T2)/DBLE(NBasF)        ! mu=Trace[GradE]/NBasF
     CALL Add(T2,ChemP)                 ! T2 = G[i+1]=G[i+1]-mu*I     
     Denominator=Dot(G,G)               ! D=(G[i],G[i])
     CALL Multiply(G,-One)              ! G=-G[i]          
     CALL Add(T2,G,T1)                  ! T1=(G[i+1]-G[i])
     CALL SetEq(G,T2)                   ! G=G[i+1]
     Numerator=Dot(T1,G)                ! N=(G[i+1]-G[i]).G[i+1]
     Gamma=MAX(Zero,Numerator/Denominator) 
     CALL Multiply(H,Gamma)             ! H=Gamma*H[i]
     CALL Add(G,H,T1)                   ! T1=G[i+1]+Gamma*H[i]
     CALL Filter(H,T1)                  ! H[i+1]=Filter[G[i+1]+Gamma*H[i]]

!     WRITE(*,*)' Beta[',ICG+1,'] = ',Gamma
!     WRITE(*,*)' HDot = ',Dot(H,H)
!     WRITE(*,*)' GDot = ',Dot(G,G)     
!--------
!      IF(ICG==1)CALL Halt(' ICG')
   ENDDO
!-------------------------
!  Tidy up a bit ...
!
   CALL Delete(G)
   CALL Delete(H) 
   CALL Delete(F)
#ifdef PARALLEL
   CALL Delete(T3)
#endif
!=============================================================
!  PURIFICATION CYCLE
!
   DeltaN_Old=1.D10
   CALL SetEq(T1,P)          ! T1=P
   NPur=0
   DO J=0,30
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
!     Stagnation Check
!
      CALL Multiply(P,-One)                ! P=-P[J]
      CALL Add(T1,P,T2)                    ! T2=P[J+1]-P[J]
      DeltaP=Max(T2)                       ! DeltaP=MAX(P[J+1]-P[J])
!      IF(DeltaP<Thresholds%Trix)EXIT       ! Convergence check
      TrP=Trace(T1)                        ! N[J+1]=Tr{P[J+1]}
      DeltaN=DABS(Two*TrP-DBLE(NEl))       ! DeltaN=2*Tr{P}-N_El
      IF(DeltaN/=Zero)THEN
         Diff=DABS(DLOG10(DeltaN/DeltaN_Old))
      ELSE
         Diff=Zero
      ENDIF
      IF(Diff<0.1D0.AND.DeltaP<1.D-3.AND.J>4)EXIT  ! Stagnation ?
      CALL SetEq(P,T1)                     ! P=P[J+1]
      DeltaN_Old=DeltaN
#ifdef PARALLEL
      IF(MyId==ROOT)THEN
#endif
         IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
            CALL OpenASCII(OutFile,Out)
            CALL PrintProtectL(Out)
            WRITE(Out,111)J,Two*TrP,DeltaP,Two*DeltaN
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
         ENDIF
         WRITE(*,111)J,Two*TrP,DeltaP,Two*DeltaN
#ifdef PARALLEL
      ENDIF
#endif
   ENDDO
   CALL Delete(P)
!   CALL PPrint(PerfMon,File_O=TRIM(SCFName)//'.SCALING_RESULTS', &
!               BareBones_O=.TRUE.)
!-----------------------------------------------
!  Report SDMM statistics
!
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       CALL OpenASCII(OutFile,Out)
       CALL PrintProtectL(Out)
       IF(LevelShift)THEN
          Mssg=ProcessName(Prog)//'Level shift = '//TRIM(DblToShrtChar(LShift))//' .'
          WRITE(Out,*)TRIM(Mssg)
       ENDIF
       Mssg=ProcessName(Prog)//TRIM(IntToChar(NCG))//' CG, and ' &
             //TRIM(IntToChar(NPur))//' purification steps taken.'
       WRITE(Out,*)TRIM(Mssg)
       Mssg=ProcessName(Prog)//'MAX(/\P) = '//TRIM(DblToShrtChar(DeltaP))//', ' &
             //'||Tr{P}-N_El|| = '//TRIM(DblToShrtChar(DeltaN))//' .'
       WRITE(Out,*)TRIM(Mssg)
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL
   ENDIF
#endif
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
111 FORMAT(' ',I2,' Tr{P} = ',D16.8,' /P = ',D10.3,' /NEl = ',D10.3)
!
END PROGRAM ! SDMM

