!    COMPUTE AND SET THRESHOLDS AND INTERMEDIATE VALUES USED IN THRESHOLDING
!    Author: Matt Challacombe and C.J. Tymczak
!------------------------------------------------------------------------------
MODULE Thresholding
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  USE SpecFun
  USE Parse
  USE McMurchie
#ifdef MMech
  USE Mechanics
#endif
!-------------------------------------------------  
!  Primary thresholds
   TYPE(TOLS), SAVE :: Thresholds
!-------------------------------------------------------------------------
! Intermediate thresholds
  REAL(DOUBLE), SAVE  :: MinZab,MinXab
  REAL(DOUBLE), SAVE  :: AtomPairDistanceThreshold  ! Atom pairs
  REAL(DOUBLE), SAVE  :: PrimPairDistanceThreshold  ! Prim pairs
  CONTAINS
!====================================================================================================
!    Set and load global threholding parameters
!====================================================================================================
     SUBROUTINE SetThresholds(Base)
         INTEGER          :: NExpt,Lndx
         TYPE(DBL_VECT)   :: Expts
         CHARACTER(LEN=*) :: Base
!        Get the primary thresholds
         CALL Get(Thresholds,Tag_O=Base)
#ifdef MMech
         IF(HasQM())THEN
#endif
!        Get distribution exponents
         CALL Get(NExpt,'nexpt',Tag_O=Base)
         CALL Get(Lndx ,'lndex',Tag_O=Base)
         CALL New(Expts,NExpt)
         CALL Get(Expts,'dexpt',Tag_O=Base)
!        MinZab=MinZa+MinZb, MinZa=MinZb
         MinZab=Expts%D(1)
!        MinXab=MinZa*MinZb/(MinZa+MinZb)=MinZab/4
         MinXab=MinZab/Four
!        Delete Exponents
         CALL Delete(Expts)
!        Set Atom-Atom thresholds
         CALL SetAtomPairThresh(Thresholds%Dist)
!        Set Prim-Prim thresholds
         CALL SetPrimPairThresh(Thresholds%Dist)
#ifdef MMech
         ENDIF
#endif
     END SUBROUTINE SetThresholds
!====================================================================================================
!    Preliminary worst case thresholding at level of atom pairs
!    Using Exp[-MinXab*|A-B|^2] = Tau 
!====================================================================================================
     SUBROUTINE SetAtomPairThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        AtomPairDistanceThreshold=-LOG(Tau)/MinXab
     END SUBROUTINE SetAtomPairThresh
!
     FUNCTION TestAtomPair(Pair,Box_O)
        LOGICAL                   :: TestAtomPair
        TYPE(AtomPair)            :: Pair
        TYPE(BBox),OPTIONAL       :: Box_O
        IF(Pair%AB2>AtomPairDistanceThreshold) THEN
           TestAtomPair = .FALSE.
        ELSE
           TestAtomPair = .TRUE.
           IF(PRESENT(Box_O)) &
              TestAtomPair=TestBoxPairOverlap(Pair,Box_O)
        ENDIF
     END FUNCTION TestAtomPair

     FUNCTION TestBoxPairOverlap(Pair,Box)
        LOGICAL                   :: TestBoxPairOverlap
        TYPE(AtomPair)            :: Pair
        TYPE(BBox),OPTIONAL       :: Box
        REAL(DOUBLE)              :: PrBxSeprtn,PairExtent,PairHlfWdt
        REAL(DOUBLE),DIMENSION(3) :: PairMidPnt,PairUntVct, &
                                     PairBoxVct,PrBxUntVct
!----------------------------------------------------------------
!        TestBoxPairOverlap=.TRUE.
!        RETURN

!       Midpoint vector of AB
        PairMidPnt=(Pair%A+Pair%B)*Half
!       Distance from Box center to line midpoint
        PairBoxVct=Box%Center-PairMidPnt
!       If pair and box are not sepertated, return with overlap=true
        PrBxSeprtn=SQRT(PairBoxVct(1)**2+PairBoxVct(2)**2+PairBoxVct(3)**2)
        IF(PrBxSeprtn==Zero)THEN
           TestBoxPairOverlap=.TRUE.
           RETURN
        ENDIF
!       Pair box unit vector
        PrBxUntVct=PairBoxVct/PrBxSeprtn
!       Max extent of atom pair from the line AB 
        IF(ABS(Pair%AB2)<1D-20)THEN
!          Line is a point: Exp[-MinZab*R^2]<Tau
           PairExtent=SQRT(PrimPairDistanceThreshold/MinZab)
           PairHlfWdt=Zero
           PairUntVct=Zero
        ELSE
!          Exp[-MinXab*|A-B|^2]*Exp[-MinZab*R^2]<Tau
           PairExtent=SQRT(MAX(Zero,(PrimPairDistanceThreshold-MinXab*Pair%AB2/MinZab)))
           PairUntVct=(Pair%A-Pair%B)*Half
!          Half length (width) of AB
           PairHlfWdt=SQRT(PairUntVct(1)**2+PairUntVct(2)**2+PairUntVct(3)**2)
!          AB Unit vector
           PairUntVct=PairUntVct/PairHlfWdt
        ENDIF
!       Check for line box overlap
        TestBoxPairOverlap=.FALSE.
        IF(ABS(PairBoxVct(1))-PairExtent*ABS(PrBxUntVct(1)) > &
          +Box%Half(1)+PairHlfWdt*ABS(PairUntVct(1)))RETURN
        IF(ABS(PairBoxVct(2))-PairExtent*ABS(PrBxUntVct(2)) > &
          +Box%Half(2)+PairHlfWdt*ABS(PairUntVct(2)))RETURN
        IF(ABS(PairBoxVct(3))-PairExtent*ABS(PrBxUntVct(3)) > &
          +Box%Half(3)+PairHlfWdt*ABS(PairUntVct(3)))RETURN
!       Could use some tests invovling cross products here...
        TestBoxPairOverlap=.TRUE.
     END FUNCTION TestBoxPairOverlap
!====================================================================================================
!    Secondary thresholding of primitive pairs using current
!    value of Xab and Exp[-Xab |A-B|^2] = Tau
!====================================================================================================
     SUBROUTINE SetPrimPairThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        PrimPairDistanceThreshold=-LOG(Tau)
     END SUBROUTINE SetPrimPairThresh
!
     FUNCTION TestPrimPair(Xi,Dist)
        LOGICAL                :: TestPrimPair
        TYPE(AtomPair)         :: Pair
        REAL(DOUBLE)           :: Xi,Dist
        IF(Xi*Dist>PrimPairDistanceThreshold) THEN
           TestPrimPair = .FALSE.
        ELSE
           TestPrimPair = .TRUE.
        ENDIF
     END FUNCTION TestPrimPair
!
     FUNCTION Extent(Ell,Zeta,HGTF,Tau_O,ExtraEll_O,Potential_O) RESULT (R)
       INTEGER                         :: Ell
       REAL(DOUBLE)                    :: Zeta
       REAL(DOUBLE),DIMENSION(:)       :: HGTF
       REAL(DOUBLE),OPTIONAL           :: Tau_O
       INTEGER,OPTIONAL                :: ExtraEll_O
       LOGICAL,OPTIONAL                :: Potential_O
       REAL(DOUBLE)                    :: R,R0,R1,R2
       IF(Zeta>=(NuclearExpnt-1D1))THEN
          ! Quick turn around for nuclei
          R=1.D-10
       ELSE
          R=Extent0(Ell,Zeta,HGTF,Tau_O=Tau_O,ExtraEll_O=ExtraEll_O,Potential_O=Potential_O)
       ENDIF
    END FUNCTION Extent
!===================================================================================================
!     Simple expressions to determine largest extent R for a distribution rho_LMN(R)
!     outside of which its value at a point is less than Tau (default) or outside 
!     of which the error made using the classical potential is less than Tau (Potential option) 
!===================================================================================================
     FUNCTION Extent0(Ell,Zeta,HGTF,Tau_O,ExtraEll_O,Potential_O) RESULT (R)
       INTEGER                         :: Ell
       REAL(DOUBLE)                    :: Zeta
       REAL(DOUBLE),DIMENSION(:)       :: HGTF
       REAL(DOUBLE),OPTIONAL           :: Tau_O
       INTEGER,OPTIONAL                :: ExtraEll_O
       LOGICAL,OPTIONAL                :: Potential_O
       INTEGER                         :: L,M,N,Lp,Mp,Np,LMN,ExtraEll       
       LOGICAL                         :: Potential
       REAL(DOUBLE)                    :: Tau,T,R,CramCo,MixMax,ScaledTau,ZetaHalf,HGInEq
       REAL(DOUBLE),PARAMETER          :: K3=1.09D0**3
       REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
                                                        720D0,5040D0,40320D0,362880D0,   &
                                                        3628800D0,39916800D0,479001600D0/)
!------------------------------------------------------------------------------------------------------------------
       IF(PRESENT(ExtraEll_O))THEN 
          ExtraEll=ExtraEll_O
       ELSE 
          ExtraEll=0
       ENDIF
       IF(PRESENT(Tau_O)) THEN
          Tau=Tau_O
       ELSE
          Tau=Thresholds%Dist
       ENDIF
       IF(PRESENT(Potential_O)) THEN
          Potential=Potential_O
       ELSE
          Potential=.FALSE.
       ENDIF
       IF(Ell+ExtraEll==0)THEN
          ! For S functions we can use tighter bounds (no halving of exponents)
          ScaledTau=Tau/(ABS(HGTF(1))+SMALL_DBL)
          IF(Potential)THEN
             ! R is the boundary of the quantum/classical potential approximation
             ! HGTF(1)*Int dr [(Pi/Zeta)^3/2 delta(r)-Exp(-Zeta r^2)]/|r-R| < Tau
             R=PFunk(Zeta,ScaledTau)
          ELSE
             ! Gaussian solution gives HGTF(1)*Exp[-Zeta*R^2] <= Tau
             R=SQRT(MAX(SMALL_DBL,-LOG(ScaledTau)/Zeta))
          ENDIF
       ELSE
          ! Universal prefactor based on Cramers inequality:
          ! H_n(t) < K 2^(n/2) SQRT(n!) EXP(t^2/2), with K=1.09
          CramCo=SMALL_DBL
          DO L=0,Ell
             DO M=0,Ell-L
                DO N=0,Ell-L-M
                   LMN=LMNDex(L,M,N)
                   MixMax=Fact(L+ExtraEll)*Fact(M)*Fact(N)
                   MixMax=MAX(MixMax,Fact(L)*Fact(M+ExtraEll)*Fact(N))
                   MixMax=MAX(MixMax,Fact(L)*Fact(M)*Fact(N+ExtraEll))
                   HGInEq=SQRT(MixMax*(Two*Zeta)**(L+M+N+ExtraEll))*HGTF(LMN)
                   CramCo=MAX(CramCo,ABS(HGInEq))
                   ! CramCo=CramCo+ABS(HGInEq)
                ENDDO
             ENDDO       
          ENDDO
          ! Now we just use expresions based on spherical symmetry but with half the exponent ...
          ZetaHalf=Half*Zeta
          ! and the threshold rescaled by the Cramer coefficient:
          ScaledTau=Tau/CramCo
          IF(Potential)THEN
             ! R is the boundary of the quantum/classical potential approximation
             ! CCo*Int dr [(2 Pi/Zeta)^3/2 delta(r)-Exp(-Zeta/2 r^2)]/|r-R| < Tau
             R=PFunk(ZetaHalf,ScaledTau)
          ELSE
             ! Gaussian solution gives CCo*Exp[-Zeta*R^2/2] <= Tau
             R=SQRT(MAX(SMALL_DBL,-LOG(ScaledTau)/ZetaHalf))
          ENDIF
       ENDIF
     END FUNCTION Extent0
!====================================================================================================
!    COMPUTE THE R THAT SATISFIES (Pi/z)^(3/2) Erfc[Sqrt[z]*R]/R < Tau 
!====================================================================================================
     FUNCTION PFunk(Zeta,Tau) RESULT(R)
        REAL(DOUBLE)  :: Tau,Zeta,SqZ,NewTau,Ec,R,BisR,DelR,X,CTest
        INTEGER       :: J,K              
!---------------------------------------------------------------------- 
        SqZ=SQRT(Zeta)
        NewTau=(Pi/Zeta)**(-1.5D0)*Tau
        DelR=Erf_Switch/(Two*SqZ)
        BisR=Zero
        ! Root finding 
        DO K=1,1000
           ! Half the step size
           DelR=Half*DelR
           ! New midpoint
           R=BisR+DelR
           X=SqZ*R
           ! Compute Erfc[Sqrt(Zeta)*R]
           IF(X>Erf_Switch)THEN
              Ec=Zero
           ELSE
              J=AINT(X*Erf_Grid)
              Ec=One-(Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*(Erf_3(J)+X*Erf_4(J)))))
           ENDIF
           Val=Ec/R
           CTest=(Val-NewTau)/Tau
           IF(R<1.D-30)THEN
              R=SMALL_DBL
              RETURN
           ELSEIF(ABS(CTest)<1D-8)THEN
              RETURN
           ELSEIF(DelR<Two*SMALL_DBL)THEN
              WRITE(*,*)' Tau= ',Tau
              WRITE(*,*)' Zeta = ',Zeta
              WRITE(*,*)' R    = ',R
              WRITE(*,*)' DelR = ',DelR
              CALL Halt(' Faild in Extent of Overlap ')
           ENDIF
           ! If still to the left, increment bisection point
           IF(CTest>Zero)BisR=R
        ENDDO
     END FUNCTION PFunk    
!===================================================================================================
!    Recursive bisection to determine largest extent for this distribution
!    outside of which its contribution to the density and gradient is less than Tau
!===================================================================================================
     FUNCTION Extent1(Ell,Zeta,HGTF,Tau_O,ExtraEll_O,Potential_O) RESULT (R)
       INTEGER                         :: Ell,ExtraEll
       REAL(DOUBLE)                    :: Zeta
       REAL(DOUBLE),DIMENSION(:)       :: HGTF
       REAL(DOUBLE),OPTIONAL           :: Tau_O
       INTEGER,OPTIONAL                :: ExtraEll_O
       LOGICAL,OPTIONAL                :: Potential_O
       REAL(DOUBLE),DIMENSION(0:HGEll) :: Co,HGPot
       REAL(DOUBLE)                    :: Tau,FUN,F0,F1
       INTEGER                         :: J,L,K,M,N,LMN,SN,LL
       REAL(DOUBLE)                    :: ConvergeTo,RMIN,RMAX,R,RErr
       LOGICAL                         :: Potential
!
       ConvergeTo=1D-8
!
       IF(PRESENT(ExtraEll_O))THEN
          ExtraEll = ExtraEll_O
       ELSE 
          ExtraEll = 0
       ENDIF
!
       IF(PRESENT(Potential_O)) THEN
          Potential = Potential_O
       ELSE
          Potential = .FALSE.
       ENDIF
!
       IF(PRESENT(Tau_O)) THEN
          Tau=Tau_O
       ELSE
          Tau=Thresholds%Dist
       ENDIF
!      Take the spherical average of HGTF coefficients      

       DO L=0,Ell
          Co(L)=Zero
          DO LMN=LBegin(L),LEnd(L)
!             Co(L)=MAX(Co(L),ABS(HGTF(LMN)))
             Co(L)=Co(L)+ABS(HGTF(LMN))
          ENDDO
       ENDDO
!
!      Compute extent of a Hermite Gaussian overlap
!
       IF(.NOT. Potential )THEN
          RMIN = Zero
          RMAX = SQRT(EXP_SWITCH/Zeta)
!
          F0 = Zero
          DO L=ExtraELL,Ell+ExtraEll
             F0 = F0 + Co(L-ExtraEll)*AsymHGTF(L,Zeta,RMIN)
          ENDDO
          IF(F0 < Zero) THEN
             CALL MondoHalt(-100,'F0 < 0')
          ENDIF
          IF(F0 < Tau) THEN
             R = Zero
             RETURN
          ENDIF
!
          F1 = Zero
          DO L=ExtraELL,Ell+ExtraEll
             F1 = F1 + Co(L-ExtraEll)*AsymHGTF(L,Zeta,RMAX)
          ENDDO
          IF(F1>Tau) THEN
             R = RMAX
             RETURN
          ENDIF
!
          R = Half*(RMIN+RMAX)
          DO J=1,200
             FUN  = Zero
             DO L=ExtraEll,Ell+ExtraEll             
                FUN  = FUN +Co(L-ExtraEll)*AsymHGTF(L,Zeta,R)
             ENDDO
             FUN = FUN-Tau
             IF(FUN < Zero) THEN
                RMAX = R
             ELSEIF(FUN > Zero) THEN
                RMIN = R
             ENDIF
             RErr = ABS(R-Half*(RMIN+RMAX))
             R = Half*(RMIN+RMAX)
             IF(RErr < ConvergeTo) GOTO 100
          ENDDO
          CALL MondoHalt(-100,'Overlap did not converge in 200 iterations')
100       CONTINUE
!
!      Do a Potential overlap
!
       ELSEIF( Potential ) THEN
          RMIN = 1.D-14
          RMAX = SQRT(GAMMA_SWITCH/Zeta)
!
          F0    = Zero
          HGPot = AsymPot(Ell+ExtraEll,Zeta,RMIN)
          DO L=ExtraELL,Ell+ExtraEll
             F0 = F0 + Co(L-ExtraEll)*ABS(HGPot(L))
          ENDDO
          IF(F0 < Zero) THEN
             CALL MondoHalt(-100,'F0 < 0')
          ENDIF
          IF(F0 < Tau) THEN
             R = Zero
             RETURN
          ENDIF
!
          F1 = Zero
          HGPot = AsymPot(Ell+ExtraEll,Zeta,RMAX)
          DO L=ExtraELL,Ell+ExtraEll
             F1 = F1 + Co(L-ExtraEll)*ABS(HGPot(L))
          ENDDO
          IF(F1>Tau) THEN
             R = RMAX
             RETURN
          ENDIF
!
          R = Half*(RMIN+RMAX)
          DO J=1,200
             FUN  = Zero
             HGPot = AsymPot(Ell+ExtraEll+1,Zeta,R)
             DO L=ExtraELL,Ell+ExtraEll 
                FUN  = FUN +Co(L-ExtraEll)*ABS(HGPot(L))
             ENDDO
             FUN = FUN-Tau
             IF(FUN < Zero) THEN
                RMAX = R
             ELSEIF(FUN > Zero) THEN
                RMIN = R
             ENDIF
             RErr = ABS(R-Half*(RMIN+RMAX))
             R = Half*(RMIN+RMAX)
             IF(RErr < ConvergeTo) GOTO 200
          ENDDO
          CALL MondoHalt(-100,'Potential Overlap did not converge in 200 iterations')
200       CONTINUE                     
       ENDIF
     END FUNCTION Extent1
!===================================================================================
!     Recursive bisection to determine largest extent for this distribution, outside
!     outside of which its contribution to the density and gradient is less than Tau
!     Has the possible drawback of finding zeros (roots) in HG functions.  
!===================================================================================
      FUNCTION Extent2(Ell,Zeta,HGTF,Tau_O,ExtraEll_O,Potential_O) RESULT (R)
         INTEGER                         :: Ell
         REAL(DOUBLE)                    :: Zeta
         REAL(DOUBLE),OPTIONAL           :: Tau_O
         REAL(DOUBLE),DIMENSION(:)       :: HGTF
         INTEGER,OPTIONAL                :: ExtraEll_O
         LOGICAL,OPTIONAL                :: Potential_O
         REAL(DOUBLE),DIMENSION(0:HGEll) :: Co,ErrR
         REAL(DOUBLE),DIMENSION(0:HGEll, &
                                0:HGEll) :: HGErr
         REAL(DOUBLE),DIMENSION(0:20)    :: LambdaR
         REAL(DOUBLE)                    :: R
         INTEGER                         :: J,L,K,M,N,LMN,ExtraEll,LTot,n1,n2,j1
         REAL(DOUBLE)                    :: ConvergeTo,Tau,R2,DelR,BisR,CTest, &
                                            RhoR,dRho,MidRho,Xpt,TwoZ, &
                                            Omega,RPE,RTE,T,upq
         LOGICAL                         :: PotentialQ
!----------------------------------------------------------------------------------
         ConvergeTo=1D-8
!
         IF(PRESENT(ExtraEll_O))THEN
            ExtraEll=ExtraEll_O
         ELSE
            ExtraEll=1
         ENDIF

         IF(PRESENT(Tau_O)) THEN
            Tau=Tau_O
         ELSE
            Tau=Thresholds%Dist
         ENDIF
!        Take the spherical max of HGTF coefficients
         DO L=0,Ell
            Co(L)=Zero
            DO LMN=LBegin(L),LEnd(L)
               Co(L)=Co(L)+ABS(HGTF(LMN))
            ENDDO
         ENDDO
!        Zero extent check
         IF(SUM(Co(0:L))==Zero)THEN
            R=Zero
            RETURN
         ENDIF
!        Do straight overlap type extent
         IF(.NOT.PRESENT(Potential_O))THEN
            DelR=SQRT(EXP_SWITCH/Zeta)*4D0
            BisR=Zero
            DO J=1,200
!              Half the step size
               DelR=Half*DelR
!              New midpoint
               R=BisR+DelR
!              Compute radial HGTF[R]
               RhoR=Zero
               dRho=Zero
               R2=R*R
               Xpt=EXP(-Zeta*R2)
               TwoZ=Two*Zeta
               LambdaR(0)=Xpt
               LambdaR(1)=TwoZ*R*Xpt
               DO L=2,Ell+ExtraEll
                  LambdaR(L)=TwoZ*(R*LambdaR(L-1)-DBLE(L-1)*LambdaR(L-2))
               ENDDO
               DO L=0,Ell
                  RhoR=RhoR+Co(L)*LambdaR(L)
                  dRho=dRho+Co(L)*LambdaR(L+ExtraEll)
               ENDDO
               MidRho=MAX(ABS(dRho),ABS(RhoR))
!              Convergence test
               CTest=(MidRho-Tau)/Tau
               IF(R<1.D-30)THEN
                  R=Zero
                  RETURN
               ELSEIF(ABS(CTest)<ConvergeTo)THEN
                  RETURN
               ELSEIF(DelR<1.D-32)THEN
                  RETURN
               ENDIF
!              If still to the left, increment bisection point
               IF(CTest>Zero)BisR=R
            ENDDO
            CALL Halt(' Faild in Extent of Overlap ')
          ELSE
            RTE=Zeta*NuclearExpnt
            RPE=Zeta+NuclearExpnt
            Omega=RTE/RPE
            Upq=TwoPi5x2/(RTE*SQRT(RPE)) &
               *(NuclearExpnt/Pi)**(ThreeHalves) ! add on moment for delta function...
            DelR=SQRT(GAMMA_SWITCH/Zeta)
            LTot=Ell+ExtraEll
            BisR=Zero
            DO K=1,200
!              Half the step size
               DelR=Half*DelR
!              New midpoint
               R=BisR+DelR
               T=Omega*R*R
               CALL ErrInts(HGEll,LTot,ErrR,Omega,T)
               DO J=0,LTot
                  HGErr(0,J)=Upq*ErrR(J)
               ENDDO
               DO J=0,LTot-1
                  J1=J+1
                  HGErr(1,J)=HGErr(0,J1)*R
               ENDDO
               DO N=2,LTot
                  N1=N-1
                  N2=N-2
                  DO J=0,LTot-N
                     J1=J+1
                     HGErr(N,J)=HGErr(N1,J1)*R+HGErr(N2,J1)*DBLE(N1)
                  ENDDO
               ENDDO
               RhoR=Zero
               dRho=Zero
               DO L=0,Ell
                  RhoR=RhoR+Co(L)*HGErr(L,0)
                  dRho=dRho+Co(L)*HGErr(L+ExtraEll,0)
               ENDDO
               MidRho=MAX(ABS(dRho),ABS(RhoR))
!              Convergence test
               CTest=(MidRho-Tau)/Tau
               IF(R<1.D-30)THEN
                  R=Zero
                  RETURN
               ELSEIF(ABS(CTest)<ConvergeTo)THEN
                  RETURN
               ELSEIF(DelR<1.D-32)THEN
                  RETURN
               ENDIF
!              If still to the left, increment bisection point
               IF(CTest>Zero)BisR=R
            ENDDO
            CALL Halt(' Faild in Extent of Potential ')
         ENDIF
       END FUNCTION Extent2
!===================================================================================
!    Norm*Gamma[L/2+3/2,Zeta*R*R]
!===================================================================================
     FUNCTION AsymHGTF(L,Zeta,R)
       INTEGER                      :: L,N,LL
       REAL(DOUBLE)                 :: Zeta,R,AsymHGTF,Norm
       REAL(DOUBLE),DIMENSION(0:16) :: NFactor = (/ 1.00000000000000000D0,&
                              1.12837916709551257D0,2.00000000000000000D0,&
                              4.51351666838205030D0,1.20000000000000000D1,&
                              3.61081333470564024D1,1.20000000000000000D2,&
                              4.33297600164676828D2,1.68000000000000000D3,&
                              6.93276160263482925D3,3.02400000000000000D4,&
                              1.38655232052696585D5,6.65280000000000000D5,&
                              3.32772556926471804D6,1.72972800000000000D7,&
                              9.31763159394121052D7,5.18918400000000000D8 /)
       N    = (-1)**L
       Norm = (Zeta**(Half*DBLE(L)))*NFactor(L)
       LL   = (L+3)/2
       IF(L == 0) THEN
          AsymHGTF = EXP(-Zeta*R*R)
       ELSE
          IF(N == 1) THEN
             AsymHGTF = Norm*GammaHalf(LL,Zeta*R*R)
          ELSE
             AsymHGTF = Norm*GammaOne(LL,Zeta*R*R)
          ENDIF
       ENDIF
     END FUNCTION AsymHGTF
!===================================================================================
!
!===================================================================================
     FUNCTION AsymPot(L,Zeta,R)
       INTEGER                                 :: L,J,N1,N2,N,J1
       REAL(DOUBLE)                            :: Zeta,R
       REAL(DOUBLE)                            :: Upq,T
       REAL(DOUBLE),DIMENSION(0:HGEll)         :: ErrR,AsymPot
       REAL(DOUBLE),DIMENSION(0:HGEll,0:HGEll) :: HGErr
!---------------------------------------------------------------------------------- 
       Upq=Two*Pi/Zeta
       T=Zeta*R*R
       CALL ErrInts(HGEll,L,ErrR,Zeta,T)
       DO J=0,L
          HGErr(0,J)=Upq*ErrR(J)
       ENDDO
       DO J=0,L-1
          J1=J+1
          HGErr(1,J)=HGErr(0,J1)*R
       ENDDO
       DO N=2,L
          N1=N-1
          N2=N-2
          DO J=0,L-N
             J1=J+1
             HGErr(N,J)=HGErr(N1,J1)*R+HGErr(N2,J1)*DBLE(N1)
          ENDDO
       ENDDO
       AsymPot=Zero
       DO J=0,L
          AsymPot(J) = HGErr(J,0)
       ENDDO
     END FUNCTION AsymPot

END MODULE Thresholding
