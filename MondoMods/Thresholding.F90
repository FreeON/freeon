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
              TestAtomPair=BoxPairOverlap(Pair,Box_O)
        ENDIF
     END FUNCTION TestAtomPair
     !==================================================================================================
     ! Modified Miguel Gomez' Ray-AABB test (Gamasutra '99) to work for Cylinder-AABB
     ! for Ray=Point, uses A Simple Method for Box-Sphere Intersection Testing, Graphics Gems, pp. 247-250
     !==================================================================================================
     FUNCTION BoxPairOverlap(Pair,Box)
       LOGICAL                   :: BoxPairOverlap
       INTEGER                   :: I
       TYPE(AtomPair)            :: Pair
       TYPE(BBox),OPTIONAL       :: Box
       REAL(DOUBLE)              :: Ext,R,LMag,HL,d,s
       REAL(DOUBLE),DIMENSION(3) :: RayAB,LHat,T,THat,TMag
       !----------------------------------------------------------------
       BoxPairOverlap=.FALSE.
       IF(Pair%AB2>1.D-20)THEN ! This is the Cylinder-Box test
          ! Solve  Exp[-MinXab*|A-B|^2]*Exp[-MinZab*Ext^2]<Tau
          Ext=SQRT(MAX(Zero,(PrimPairDistanceThreshold-MinXab*Pair%AB2)/MinZab))
          ! Ray from A to B
          RayAB=Pair%A-Pair%B
          ! Length of RayAB
          LMag=SQRT(RayAB(1)**2+RayAB(2)**2+RayAB(3))
          ! AB unit vector
          LHat=RayAB/LMag
          ! Half length
          HL=Half*LMag
          ! Distance from Box center to RayAB midpoint
          T=Box%Center-(Pair%A+Pair%B)*Half
          TMag=SQRT(T(1)**2+T(2)**2+T(3)**2)
          ! Unit vector in direction T
          THat=T/TMag
          ! Recompute T to change test from Ray-Box to Cylinder-Box
          T=THat*MAX(Zero,TMag-Ext)
          ! Seperating axis tests
          DO I=1,3
             IF(ABS(T(I))>Box%Half(I)+HL*ABS(LHat(I)))RETURN
          ENDDO
          ! Cross product tests
          r=Box%Half(2)*ABS(LHat(3))+Box%Half(3)*ABS(LHat(2))
          IF(ABS(T(2)*LHat(3)-T(3)*LHat(2))>r)RETURN
          r=Box%Half(1)*ABS(LHat(3))+Box%Half(3)*ABS(LHat(1))
          IF(ABS(T(3)*LHat(1)-T(1)*LHat(3))>r)RETURN
          r=Box%Half(1)*ABS(LHat(2))+Box%Half(2)*ABS(LHat(1))
          IF(ABS(T(1)*LHat(2)-T(2)*LHat(1))>r)RETURN
       ELSE ! Sphere-AABB from Graphics Gems
          d=Zero
          DO I=1,3
             IF(Pair%A(I)<Box%BndBox(I,1))d=d+(Pair%A(I)-Box%BndBox(I,1))**2
             IF(Pair%A(I)>Box%BndBox(I,2))d=d+(Pair%A(I)-Box%BndBox(I,2))**2
          ENDDO
          ! PrimPairDistanceThreshold is the squared radius of the sphere
          IF(d>PrimPairDistanceThreshold)RETURN
       ENDIF
       BoxPairOverlap=.TRUE.
     END FUNCTION BoxPairOverlap
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
!===================================================================================================
!     Simple expressions to determine largest extent R for a distribution rho_LMN(R)
!     outside of which its value at a point is less than Tau (default) or outside 
!     of which the error made using the classical potential is less than Tau (Potential option) 
!===================================================================================================
     FUNCTION Extent(Ell,Zeta,HGTF,Tau_O,ExtraEll_O,Potential_O) RESULT (R)
       INTEGER                         :: Ell
       REAL(DOUBLE)                    :: Zeta
       REAL(DOUBLE),DIMENSION(:)       :: HGTF
       REAL(DOUBLE),OPTIONAL           :: Tau_O
       INTEGER,OPTIONAL                :: ExtraEll_O
       LOGICAL,OPTIONAL                :: Potential_O
       INTEGER                         :: L,M,N,Lp,Mp,Np,LMN,ExtraEll       
       LOGICAL                         :: Potential
       REAL(DOUBLE)                    :: Tau,T,R,CramCo,MixMax,ScaledTau,ZetaHalf,HGInEq,TMP
       REAL(DOUBLE),PARAMETER          :: K3=1.09D0**3
       REAL(DOUBLE),DIMENSION(0:12),PARAMETER :: Fact=(/1D0,1D0,2D0,6D0,24D0,120D0,      &
                                                        720D0,5040D0,40320D0,362880D0,   &
                                                        3628800D0,39916800D0,479001600D0/)
!------------------------------------------------------------------------------------------------------------------
       ! Quick turn around for nuclei
       IF(Zeta>=(NuclearExpnt-1D1))THEN
          R=1.D-10
          RETURN
       ENDIF
       ! Misc options....
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
       ! Spherical symmetry check
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
                   CramCo=CramCo+ABS(HGInEq)
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
     END FUNCTION Extent
!====================================================================================================
!    COMPUTE THE R THAT SATISFIES (Pi/z)^(3/2) Erfc[Sqrt[z]*R]/R < Tau 
!====================================================================================================
     FUNCTION PFunk(Zeta,Tau) RESULT(R)
        REAL(DOUBLE)  :: Tau,Zeta,SqZ,NewTau,Val,Ec,R,BisR,DelR,X,CTest
        INTEGER       :: J,K              
        LOGICAL :: pp
!---------------------------------------------------------------------- 
        SqZ=SQRT(Zeta)
        NewTau=Tau*(Zeta/Pi)**(1.5D0)
        ! Quick check for max resolution of Erfc approx
        IF(1D-13*SqZ/Erf_Switch>NewTau)THEN
           R=Erf_Switch/SqZ
           RETURN
        ENDIF
!        WRITE(*,*)'=-=========================================='
        ! Ok, within resolution--do root finding...
        DelR=Erf_Switch/SqZ
        BisR=Zero
        DO K=1,100
           ! New midpoint
           R=BisR+DelR
           X=SqZ*R
           ! Compute Erfc[Sqrt(Zeta)*R]
           IF(X>=Erf_Switch)THEN
              Ec=Zero
           ELSE
              J=AINT(X*Erf_Grid)
              Ec=One-(Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*(Erf_3(J)+X*Erf_4(J)))))
           ENDIF           
           Val=Ec/R
           CTest=(Val-NewTau)/Tau
!           WRITE(*,33)R,DelR,X,CTest; 33 format(5(2x,D12.6))
           ! Go for relative error to get smoothness, but bail if 
           ! absolute accuracy of erf interpolation is exceeded         
           IF(ABS(CTest)<1D-4.OR.Tau*ABS(CTest)<1.D-12)THEN
              ! Converged
              RETURN           
           ELSEIF(DelR<1D-40)THEN
              ! This is unacceptable
              EXIT
           ENDIF
           ! If still to the left, increment bisection point
           IF(CTest>Zero)BisR=R
           DelR=Half*DelR
        ENDDO
        CALL Halt(' Failed to converge in PFunk: Tau = '//TRIM(DblToShrtChar(NewTau))//RTRN &
                                           //' CTest = '//TRIM(DblToShrtChar(CTest))//RTRN &
                                           //' Zeta = '//TRIM(DblToShrtChar(Zeta))//RTRN &
                                           //' R = '//TRIM(DblToShrtChar(R))//RTRN &
                                           //' Ec = '//TRIM(DblToShrtChar(Ec))//RTRN &
                                           //' dR = '//TRIM(DblToShrtChar(DelR))//RTRN &
                                           //' SqZ*R = '//TRIM(DblToMedmChar(X)))
     END FUNCTION PFunk    
END MODULE Thresholding
