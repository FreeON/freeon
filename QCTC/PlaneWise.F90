!------------------------------------------------------------------------------
!
!   Written by Matt Challacombe (2007)
!   A major pain in the ass, this was
!
MODULE PlaneWise
  USE Derivedtypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Macros
  USE SpecFun
  USE MemMan
  USE MondoPoles
  IMPLICIT NONE
  ! Globals
CONTAINS
!========================================================================================
!
!========================================================================================
  SUBROUTINE InnPlane(Ell,A1,A2,B1,B2,Ylm,WellSepCells_O)
    INTEGER                                :: Ell,L,M,LDex,LMDex
    INTEGER                                :: Mu1,Mu2,Lm1,Lm2,MuMu,LmLm,NC,NCell
    REAL(DOUBLE)                           :: Beta,Oa
    REAL(DOUBLE),DIMENSION(3)              :: A1,A2,B1,B2
    TYPE(CellSet),OPTIONAL                 :: WellSepCells_O
    LOGICAL                                :: WellSepQ,DoCell
    REAL(DOUBLE),DIMENSION(2)              :: S,H
    REAL(DOUBLE)                           :: Ph,Rh,Kh,Ps,Rs,Ks
    REAL(DOUBLE),DIMENSION(-2*Ell:2*Ell+1) :: Gamma
    COMPLEX(DOUBLE),DIMENSION(0:Ell,0:Ell) :: YRecp,YReal,Ylm
    REAL(DOUBLE),DIMENSION(0:Ell,0:Ell)    :: ALP
    COMPLEX(DOUBLE)                        :: Eye=(0D0,1D0)
!    REAL(DOUBLE),EXTERNAL                  :: DGamma
    REAL(DOUBLE),DIMENSION(3)              :: A1xA2
    REAL(DOUBLE)                           :: KMax,SMax,RMax,R,RMin,RToThEllPlsOne
    TYPE(DBL_RNK2)                         :: HPlane
    TYPE(DBL_VECT)                         :: HSort
    TYPE(INT_VECT)                         :: ISort
    INTEGER, PARAMETER                     :: EllSwitch=11
    REAL(DOUBLE),PARAMETER                 :: Accuracy=1D-30
    !-----------------------------------------------------------------------------------
    IF(PRESENT(WellSepCells_O))THEN
       WellSepQ=.TRUE.
    ELSE
       WellSepQ=.FALSE.
    ENDIF
    !-----------------------------------------------------------------------------------
    RMin=1D10
    R=SQRT(DOT_PRODUCT(A1,A1))
    RMin=MIN(RMin,R)
    R=SQRT(DOT_PRODUCT(A2,A2))
    RMin=MIN(RMin,R)
    Beta=SQRT(Pi)/RMin
    SMax=(One/Beta)*SQRT(ABS(LOG(Accuracy)))
    KMax=    Beta *SQRT(ABS(LOG(Accuracy)))/Pi
    RMax=(One/Accuracy)**(One/DBLE(EllSwitch))

    Ylm=(0D0,0D0)
    YRecp=(0D0,0D0)
    YReal=(0D0,0D0)
    !-----------------------------------------------------------------------------------
    A1xA2=CROSS_PRODUCT(A1,A2)
    Oa=SQRT(DOT_PRODUCT(A1xA2,A1xA2))
    !-----------------------------------------------------------------------------------
    ! ALP==LegendreP[L,M,0]/Gamma[(L+M+1)/2]
    !-----------------------------------------------------------------------------------
    DO M=0,Ell
       ALP(M,M)=FactMlm0(M)
    ENDDO
    DO M=0,Ell-1
       ALP(M+1,M)=0
    ENDDO
    DO L=2,Ell
       LDex=LTD(L)
       DO M=0,L-2
          ALP(L,M)=-ALP(L-2,M)*FactMlm2(LDex+M)
       ENDDO
    ENDDO
    DO L=0,Ell
       DO M=0,L
          ALP(L,M)=ALP(L,M)/Factorial(L-M)
       ENDDO
    ENDDO
    DO L=0,Ell
       DO M=0,L
          ALP(L,M)=ALP(L,M)/(DGamma(DBLE(L+M+1)*5D-1))
       ENDDO
    ENDDO
    !-----------------------------------------------------------------------------------
    ! RECIPROCAL SPACE PART
    !-----------------------------------------------------------------------------------
    NCell=0
    DO Mu1=-30,30
       DO Mu2=-30,30
          IF(.NOT.(Mu1==0.AND.Mu2==0))THEN
             H=Mu1*B1(1:2)+Mu2*B2(1:2)
             Rh=SQRT(H(1)*H(1)+H(2)*H(2))
             IF(Rh<KMax)THEN
                NCell=NCell+1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    CALL New(ISort,NCell)
    CALL New(HSort,NCell)
    CALL New(HPlane,(/2,NCell/))
    !
    NCell=0
    DO Mu1=-30,30
       DO Mu2=-30,30
          IF(.NOT.(Mu1==0.AND.Mu2==0))THEN
             H=Mu1*B1(1:2)+Mu2*B2(1:2)
             Rh=SQRT(H(1)*H(1)+H(2)*H(2))
             IF(Rh<KMax)THEN
                NCell=NCell+1
                HPlane%D(:,NCell)=H
                HSort%D(NCell)=Rh
                ISort%I(NCell)=NCell
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    CALL Sort(HSort,ISort,NCell,2)
    !
    DO MuMu=1,NCell
       H=HPlane%D(:,ISort%I(MuMu))
       Ph=ATAN2(H(2),H(1))
       Rh=HSort%D(MuMu)
       Kh=(Rh*Pi/Beta)**2
       CALL InPlaneNegativeGammas(Ell,Kh,Gamma(-2*Ell:1))
       DO L=0,MIN(EllSwitch,Ell)
          DO M=0,L
             YRecp(L,M)=YRecp(L,M)+Eye**M*EXP(-Eye*DBLE(M)*Ph)*Rh**(L-1)*Gamma(M-L+1)*(Pi**L)
          ENDDO
       ENDDO
    ENDDO
    !
    DO L=2,MIN(EllSwitch,Ell)
       LDex=LTD(L)
       YRecp(L,0)=YRecp(L,0)+Two*Pi*Beta**(L-1)/DBLE(L-1) & ! This is the Mu1,Mu2==0 term
            -Two*Oa*Beta**(L+1)/DBLE(L+1)   ! This is subtracting out the Lam1,Lam2==0 term
    ENDDO
    !
    DO L=0,MIN(EllSwitch,Ell)
       DO M=0,L
          YRecp(L,M)=YRecp(L,M)/Oa
       ENDDO
    ENDDO
    !
    CALL Delete(ISort)
    CALL Delete(HSort)
    CALL Delete(HPlane)
    !---------------------------------------------------------------------------------
    ! SUBTRACT THE RECIPROCAL SPACE PART OF THE NEAR FIELD FROM THE PLANEWISE CELLS
    !---------------------------------------------------------------------------------
    IF(WellSepQ)THEN
       DO NC = 1,WellSepCells_O%NCells-1
          S=WellSepCells_O%CellCarts%D(1:2,NC)
          Ps=ATAN2(S(2),S(1))
          Rs=SQRT(S(1)**2+S(2)**2)
          Ks=(Rs*Beta)**2
          DO L=0,MIN(EllSwitch,Ell)
             DO M=0,L
                YRecp(L,M)=YRecp(L,M)-EXP(-Eye*DBLE(M)*Ps)*DGamma(DBLE(L+M+1)*5D-1)     &
                                     *RegularizedGammaP(5D-1*DBLE(L+M+1),Ks)/Rs**(L+1)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !-----------------------------------------------------------------------------------
    ! REAL SPACE PART EWALD-LIKE PART (FOR LOWER ELL<=ELLSWITCH)
    !-----------------------------------------------------------------------------------
    !
    NCell=0
    DO Lm1=-30,30
       DO Lm2=-30,30
          IF(.NOT.(Lm1==0.AND.Lm2==0))THEN
             S=Lm1*A1(1:2)+Lm2*A2(1:2)
             Rs=SQRT(S(1)*S(1)+S(2)*S(2))
             IF(Rs<SMax)THEN
                NCell=NCell+1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    CALL New(ISort,NCell)
    CALL New(HSort,NCell)
    CALL New(HPlane,(/2,NCell/))
    !
    NCell=0
    DO Lm1=-30,30
       DO Lm2=-30,30
          IF(.NOT.(Lm1==0.AND.Lm2==0))THEN
             S=Lm1*A1(1:2)+Lm2*A2(1:2)
             Rs=SQRT(S(1)*S(1)+S(2)*S(2))
             IF(Rs<SMax)THEN
                NCell=NCell+1
                HPlane%D(:,NCell)=S
                HSort%D(NCell)=Rs
                ISort%I(NCell)=NCell
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    CALL Sort(HSort,ISort,NCell,2)
    !
    DO LmLm=1,NCell
       S=HPlane%D(:,ISort%I(LmLm))
       IF(WellSepQ)THEN
          DoCell=.NOT.InCell_CellSet(WellSepCells_O,S(1),S(2),0D0)
       ELSE
          DoCell=.TRUE.
       ENDIF
       IF(DoCell)THEN
          Rs=HSort%D(LmLm)
          Ps=ATAN2(S(2),S(1))
          Ks=(Rs*Beta)**2
          CALL InPlanePositiveGammas(Ell,Ks,Gamma(0:2*Ell+1))
          DO L=0,MIN(EllSwitch,Ell)
             DO M=0,L
                YReal(L,M)=YReal(L,M)+EXP(-Eye*DBLE(M)*Ps)*Gamma(L+M+1)/Rs**(L+1)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !
    CALL Delete(ISort)
    CALL Delete(HSort)
    CALL Delete(HPlane)
    !-----------------------------------------------------------------------------------
    ! REAL SPACE PART (FOR HIGH ELL)
    !-----------------------------------------------------------------------------------
    !
    NCell=0
    DO Lm1=-100,100
       DO Lm2=-100,100
          IF(.NOT.(Lm1==0.AND.Lm2==0))THEN
             S=Lm1*A1(1:2)+Lm2*A2(1:2)
             Rs=SQRT(S(1)*S(1)+S(2)*S(2))
             IF(Rs<RMax)THEN
                NCell=NCell+1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    CALL New(ISort,NCell)
    CALL New(HSort,NCell)
    CALL New(HPlane,(/2,NCell/))
    !
    NCell=0
    DO Lm1=-100,100
       DO Lm2=-100,100
          IF(.NOT.(Lm1==0.AND.Lm2==0))THEN
             S=Lm1*A1(1:2)+Lm2*A2(1:2)
             Rs=SQRT(S(1)*S(1)+S(2)*S(2))
             IF(Rs<RMax)THEN
                NCell=NCell+1
                HPlane%D(:,NCell)=S
                HSort%D(NCell)=Rs
                ISort%I(NCell)=NCell
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    !
    CALL Sort(HSort,ISort,NCell,2)
    !
    DO LmLm=1,NCell
       S=HPlane%D(:,ISort%I(LmLm))
       IF(WellSepQ)THEN
          DoCell=.NOT.InCell_CellSet(WellSepCells_O,S(1),S(2),0D0)
       ELSE
          DoCell=.TRUE.
       ENDIF
       IF(DoCell)THEN
          Rs=HSort%D(LmLm)
          Ps=ATAN2(S(2),S(1))
          DO L=EllSwitch+1,Ell
             DO M=0,L
                YReal(L,M)=YReal(L,M)+EXP(-Eye*DBLE(M)*Ps)/Rs**(L+1)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
!
    CALL Delete(ISort)
    CALL Delete(HSort)
    CALL Delete(HPlane)
    !---------------------------------------------------------------------------------
    ! COMBINE THE EWALD AND REAL SPACE COMPONENTS AND NORMALIZE THE HARMONICS
    !---------------------------------------------------------------------------------
    DO L=0,Ell
       DO M=0,L
          Ylm(L,M)=Ylm(L,M)+ALP(L,M)*(YRecp(L,M)+YReal(L,M))*(0D0,1D0)**(M)
       ENDDO
    ENDDO
    !
  END SUBROUTINE InnPlane

   SUBROUTINE OffPlane(Ell,A1,A2,A3,B1,B2,Ylm)
     INTEGER                                :: Ell,L,M,LDex,LMDex
     INTEGER                                :: Mu1,Mu2,Lm1,Lm2,NCell,MuMu
     REAL(DOUBLE)                           :: Beta,Oa,HMax
     REAL(DOUBLE)                           :: XiOne,XiTwo,Zeta,ZPhz,SgnP,SgnM
     REAL(DOUBLE),DIMENSION(3)              :: A1,A2,A3,B1,B2
     REAL(DOUBLE),DIMENSION(2)              :: S,H
     REAL(DOUBLE)                           :: Ph,Rh,Kh,Ps,Rs,Ks
     REAL(DOUBLE),DIMENSION(-2*Ell:2*Ell)   :: Gamma
     COMPLEX(DOUBLE),DIMENSION(0:Ell,0:Ell) :: YOff,Ylm
     COMPLEX(DOUBLE)                        :: XYPhz
     COMPLEX(DOUBLE),PARAMETER              :: Eye=(0D0,1D0)
     REAL(DOUBLE),DIMENSION(3)              :: A1xA2
     TYPE(DBL_RNK2)                         :: HPlane
     TYPE(DBL_VECT)                         :: HSort
     TYPE(INT_VECT)                         :: ISort
     !
     A1xA2=CROSS_PRODUCT(A1,A2)
     Oa=SQRT(DOT_PRODUCT(A1xA2,A1xA2))
     !
     XiOne=A3(1)
     XiTwo=A3(2)
     Zeta =A3(3)
     !
     HMax=(One/(Pi*Zeta))*ABS(LOG(1D-20))
     !
     NCell=0
     DO Mu1=-100,100
        DO Mu2=-100,100
           IF(.NOT.(Mu1==0.AND.Mu2==0))THEN
              H=Mu1*B1(1:2)+Mu2*B2(1:2)
              Rh=SQRT(H(1)*H(1)+H(2)*H(2))
              IF(Rh<HMax)THEN
                 NCell=NCell+1
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !
     CALL New(ISort,NCell)
     CALL New(HSort,NCell)
     CALL New(HPlane,(/2,NCell/))
     !
     NCell=0
     DO Mu1=-100,100
        DO Mu2=-100,100
           IF(.NOT.(Mu1==0.AND.Mu2==0))THEN
              H=Mu1*B1(1:2)+Mu2*B2(1:2)
              Rh=SQRT(H(1)*H(1)+H(2)*H(2))
              IF(Rh<HMax)THEN
                 NCell=NCell+1
                 HPlane%D(:,NCell)=H
                 HSort%D(NCell)=Rh
                 ISort%I(NCell)=NCell
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !
     CALL Sort(HSort,ISort,NCell,2)
     !
     YOff=0D0
     !
     SgnP=SIGN(One,Zeta)
     SgnM=SIGN(One,-Zeta)
     !
     DO MuMu=1,NCell
        H=HPlane%D(:,ISort%I(MuMu))
        Ph=ATAN2(H(2),H(1))
        Rh=HSort%D(MuMu)
        !
        ZPhz=EXP(Two*Pi*Rh*Zeta)
        XYPhz=EXP(Two*Pi*Eye*Rh*(XiOne*COS(Ph)+XiTwo*SIN(Ph)))
        !
        DO L=0,Ell
           DO M=0,L
              YOff(L,M)=YOff(L,M)+EXP(-Eye*DBLE(M)*Ph)*Rh**(L-1) &
                   *(One/(ZPhz*XYPhz-One)*SgnP**(L+M)      &
                    +One/(ZPhz/XYPhz-One)*SgnM**(L+M))

           ENDDO
        ENDDO
     ENDDO
     !
     DO L=0,Ell
        DO M=0,L
           Ylm(L,M)=Ylm(L,M)+YOff(L,M)*(Two*Pi)**(L)/(Oa*Factorial(L-M))
        ENDDO
     ENDDO
     !
     CALL Delete(ISort)
     CALL Delete(HSort)
     CALL Delete(HPlane)
     !
   END SUBROUTINE OffPlane
   !---------------------------------------------------------------------
   !  This code computes Gamma[n,x], with -L<n<1/2 by halves.
   !  Indexing is by negative integers, with -2*L<M<1.
   !---------------------------------------------------------------------
   SUBROUTINE InPlaneNegativeGammas(Ell,X,G)
     INTEGER :: Ell,M,TwoM
     REAL(DOUBLE) :: X,XM,SQX,EX
     REAL(DOUBLE),DIMENSION(-2*Ell:1) :: G
     REAL(DOUBLE),EXTERNAL :: DEI,DERFC
     !-------------------------------------------------------------------
     EX=EXP(-X)
     XM=SQRT(X)
     !
!     G(1)=SQRT(Pi)*DERFC(XM)
     G(1)=SQRT(Pi)*ERFC(XM)
     G(0)=-DEI(-X)
     !
     XM=One/XM
     SQX=XM
     ! Downward recurrence
     DO M=-1,-Ell,-1
        TwoM=2*M
        G(TwoM+1)=(G(TwoM+3)-XM*EX)/(M+5D-1)
        XM=XM*SQX
        G(TwoM)=(G(TwoM+2)-XM*EX)/DBLE(M)
        XM=XM*SQX
     ENDDO
     !
   END SUBROUTINE InPlaneNegativeGammas
   !---------------------------------------------------------------------
   !  This code computes Gamma[n,x], with -L<n<1/2 by halves.
   !  Indexing is by negative integers, with -2*L<M<1.
   !---------------------------------------------------------------------
   SUBROUTINE InPlanePositiveGammas(Ell,X,G)
     INTEGER :: Ell,M,TwoM
     REAL(DOUBLE) :: X,XM,SQX,EX
     REAL(DOUBLE),DIMENSION(0:2*Ell+1) :: G
     REAL(DOUBLE),EXTERNAL :: DEI,DERFC
     !-------------------------------------------------------------------
     EX=EXP(-X)
     SQX=SQRT(X)
     !
!     G(1)=SQRT(Pi)*DERFC(SQX)
     G(1)=SQRT(Pi)*ERFC(SQX)
     G(0)=-DEI(-X)
     !
     XM=One
     ! Upward recurrence
     DO M=1,Ell
        TwoM=2*M
        G(TwoM)=XM*EX+5D-1*DBLE(TwoM-2)*G(TwoM-2)
        XM=XM*SQX
        G(TwoM+1)=XM*EX+5D-1*DBLE(TwoM-1)*G(TwoM-1)
        XM=XM*SQX
     ENDDO
     !
   END SUBROUTINE InPlanePositiveGammas
   !---------------------------------------------------------------------
   !  This code computes Gamma[n,x,0], with -L<n<1/2 by halves.
   !  Indexing is by negative integers, with -2*L<M<1.
   !---------------------------------------------------------------------
   SUBROUTINE InPlanePositiveGCompliment(Ell,X,G)
     INTEGER :: Ell,M,TwoM
     REAL(DOUBLE) :: X,XM,SQX,EX,Gamma
     REAL(DOUBLE),DIMENSION(0:2*Ell+1) :: G
     REAL(DOUBLE),EXTERNAL :: DEI,DERFC
     !-------------------------------------------------------------------
     XM=One
     EX=EXP(-X)
     SQX=SQRT(X)
     !
     G(0)=0D0
!     G(1)=SQRT(Pi)*(DERFC(SQX)-1D0)
     G(1)=SQRT(Pi)*(ERFC(SQX)-1D0)
     G(2)=EX-1D0
     XM=XM*SQX
     G(3)=XM*EX+5D-1*G(1)
     XM=XM*SQX
!     G(4)=XM*EX+G(2)
!!$
!!$     !
!!$     WRITE(*,*)0,G(0)
!!$     WRITE(*,*)1,G(1)
!!$     WRITE(*,*)2,G(2)
!!$     WRITE(*,*)3,G(3)
!!$     WRITE(*,*)4,G(4)
     !
     ! Upward recurrence
     DO M=2,Ell
        TwoM=2*M
        G(TwoM)=XM*EX+5D-1*DBLE(TwoM-2)*G(TwoM-2)
!        WRITE(*,*)TwoM,G(TwoM)
        XM=XM*SQX
        G(TwoM+1)=XM*EX+5D-1*DBLE(TwoM-1)*G(TwoM-1)
!        WRITE(*,*)(TwoM+1),G(TwoM+1)
        XM=XM*SQX
     ENDDO
!     STOP
     !
   END SUBROUTINE InPlanePositiveGCompliment

  FUNCTION RegularizedGammaQ(A,X)
    REAL(DOUBLE) :: A,X,Gamma,RegularizedGammaQ,GGammaI,XAM,S,R,T0,Fct,DK
    INTEGER      :: IA,K
    XAM=-X+A*LOG(X)
    IF(XAM.GT.7D2.OR.A>170D0) THEN
       WRITE(*,*)'a and/or x too large'
       STOP
    ENDIF
    Gamma=DGamma(A)
    IF(X.EQ.0D0) THEN
       RegularizedGammaQ=1D0
    ELSEIF(X<1D0+A) THEN
       S=1D0/A
       R=S
       DO K=1,60
          R=R*X/(A+DBLE(K))
          S=S+R
          IF(ABS(R/S)<1D-17)EXIT
       ENDDO
       RegularizedGammaQ=(Gamma-EXP(XAM)*S)/Gamma
    ELSEIF (X.GT.1.0+A) THEN
       T0=0D0
       DO K=60,1,-1
          T0=(K-A)/(1.0D0+K/(X+T0))
       ENDDO
       RegularizedGammaQ=(EXP(XAM)/(X+T0))/Gamma
    ENDIF
  END FUNCTION RegularizedGammaQ

  FUNCTION RegularizedGammaP(A,X)
    REAL(DOUBLE) :: A,X,Gamma,RegularizedGammaP,GGammaI,XAM,S,R,T0,Fct,DK
    INTEGER      :: IA,K
    XAM=-X+A*LOG(X)
    IF(XAM.GT.7D2.OR.A>170D0) THEN
       WRITE(*,*)'a and/or x too large'
       STOP
    ENDIF
    Gamma=DGamma(A)
    IF(X.EQ.0D0) THEN
       RegularizedGammaP=1D0
    ELSEIF(X<1D0+A) THEN
       S=1D0/A
       R=S
       DO K=1,60
          R=R*X/(A+DBLE(K))
          S=S+R
          IF(ABS(R/S)<1D-17)EXIT
       ENDDO
       RegularizedGammaP=EXP(XAM)*S/Gamma
    ELSEIF(X.GT.1.0+A) THEN
       T0=0D0
       DO K=60,1,-1
          T0=(K-A)/(1.0D0+K/(X+T0))
       ENDDO
       RegularizedGammaP=(Gamma-EXP(XAM)/(X+T0))/Gamma
    ENDIF
  END FUNCTION RegularizedGammaP

  FUNCTION dRegularizedGammaQ(A,X)
    REAL(DOUBLE) :: A,X,Gamma,dRegularizedGammaQ
    dRegularizedGammaQ=-EXP(-X)*X**(A-HALF)/DGamma(A+Half)
  END FUNCTION DRegularizedGammaQ


  function dgamma ( x )

    !*****************************************************************************80
    !
    !! DGAMMA evaluates Gamma(X) for a real argument.
    !
    !  Discussion:
    !
    !    This routine calculates the GAMMA function for a real argument X.
    !    Computation is based on an algorithm outlined in reference 1.
    !    The program uses rational functions that approximate the GAMMA
    !    function to at least 20 significant decimal digits.  Coefficients
    !    for the approximation over the interval (1,2) are unpublished.
    !    Those for the approximation for 12 <= X are from reference 2.
    !
    !  Modified:
    !
    !    03 April 2007
    !
    !  Author:
    !
    !    William Cody, Laura Stoltz
    !
    !  Reference:
    !
    !    William Cody,
    !    An Overview of Software Development for Special Functions,
    !    in Numerical Analysis Dundee, 1975,
    !    edited by GA Watson,
    !    Lecture Notes in Mathematics 506,
    !    Springer, 1976.
    !
    !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    !    Charles Mesztenyi, John Rice, Henry Thatcher,
    !    Christoph Witzgall,
    !    Computer Approximations,
    !    Wiley, 1968,
    !    LC: QA297.C64.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, the argument of the function.
    !
    !    Output, real ( kind = 8 ) DGAMMA, the value of the function.
    !
    implicit none

    real ( double ) c(7)
    real ( double ) dgamma
    real ( double ) eps
    real ( double ) fact
    real ( double ) half
    integer i
    integer n
    real ( double ) one
    real ( double ) p(8)
    logical parity
    real ( double ) pi
    real ( double ) q(8)
    real ( double ) res
    real ( double ) sqrtpi
    real ( double ) sum
    real ( double ) twelve
    real ( double ) two
    real ( double ) x
    real ( double ) xbig
    real ( double ) xden
    real ( double ) xinf
    real ( double ) xminin
    real ( double ) xnum
    real ( double ) y
    real ( double ) y1
    real ( double ) ysq
    real ( double ) z
    real ( double ) zero
    !
    !  Mathematical constants
    !
    data one /1.0d0 /
    data half /0.5d0/
    data twelve /12.0d0/
    data two /2.0d0 /
    data zero /0.0d0/
    data sqrtpi /0.9189385332046727417803297d0/
    data pi /3.1415926535897932384626434d0/
    !
    !  Machine dependent parameters
    !
    data xbig / 171.624d0 /
    data xminin / 2.23d-308 /
    data eps /2.22d-16/
    data xinf /1.79d308/
    !
    !  Numerator and denominator coefficients for rational minimax
    !  approximation over (1,2).
    !
    data p/-1.71618513886549492533811d+0,2.47656508055759199108314d+1, &
         -3.79804256470945635097577d+2,6.29331155312818442661052d+2, &
         8.66966202790413211295064d+2,-3.14512729688483675254357d+4, &
         -3.61444134186911729807069d+4,6.64561438202405440627855d+4/
    data q/-3.08402300119738975254353d+1,3.15350626979604161529144d+2, &
         -1.01515636749021914166146d+3,-3.10777167157231109440444d+3, &
         2.25381184209801510330112d+4,4.75584627752788110767815d+3, &
         -1.34659959864969306392456d+5,-1.15132259675553483497211d+5/
    !
    !  Coefficients for minimax approximation over (12, INF).
    !
    data c/ &
         -1.910444077728d-03, &
         8.4171387781295d-04, &
         -5.952379913043012d-04, &
         7.93650793500350248d-04, &
         -2.777777777777681622553d-03, &
         8.333333333333333331554247d-02, &
         5.7083835261d-03/

    parity = .false.
    fact = one
    n = 0
    y = x
    !
    !  Argument is negative.
    !
    if ( y <= zero ) then

       y = -x
       y1 = aint ( y )
       res = y - y1

       if ( res .ne. zero ) then

          if ( y1 .ne. aint ( y1 * half ) * two ) then
             parity = .true.
          end if

          fact = -pi / sin ( pi * res )
          y = y + one

       else

          res = xinf
          dgamma = res
          return

       end if

    end if
    !
    !  Argument is positive.
    !
    if ( y < eps ) then
       !
       !  Argument < EPS.
       !
       if ( y .ge. xminin ) then
          res = one / y
       else
          res = xinf
          dgamma = res
          return
       end if

    else if ( y < twelve ) then

       y1 = y
       !
       !  0.0 < argument < 1.0.
       !
       if ( y < one ) then

          z = y
          y = y + one
          !
          !  1.0 < argument < 12.0.
          !  Reduce argument if necessary.
          !
       else

          n = int ( y ) - 1
          y = y - real ( n, double )
          z = y - one

       end if
       !
       !  Evaluate approximation for 1.0 < argument < 2.0.
       !
       xnum = zero
       xden = one
       do i = 1, 8
          xnum = ( xnum + p(i) ) * z
          xden = xden * z + q(i)
       end do

       res = xnum / xden + one
       !
       !  Adjust result for case  0.0 < argument < 1.0.
       !
       if ( y1 < y ) then

          res = res / y1
          !
          !  Adjust result for case 2.0 < argument < 12.0.
          !
       else if ( y < y1 ) then

          do i = 1, n
             res = res * y
             y = y + one
          end do

       end if

    else
       !
       !  Evaluate for 12.0 <= argument.
       !
       if ( y <= xbig ) then

          ysq = y * y
          sum = c(7)
          do i = 1, 6
             sum = sum / ysq + c(i)
          end do
          sum = sum / y - y + sqrtpi
          sum = sum + ( y - half ) * log ( y )
          res = exp ( sum )

       else

          res = xinf
          dgamma = res
          return

       end if

    end if
    !
    !  Final adjustments and return.
    !
    if ( parity ) then
       res = -res
    end if

    if ( fact .ne. one ) then
       res = fact / res
    end if

    dgamma = res

    !  return
  end function dgamma


END MODULE PlaneWise
