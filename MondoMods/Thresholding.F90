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
  USE McMurchie
!-------------------------------------------------  
!  Primary thresholds
!
   TYPE(TOLS), SAVE :: Thresholds
!-------------------------------------------------------------------------
! Intermediate thresholds
!
  INTEGER,SAVE        :: MinRadialAngSym
  REAL(DOUBLE), SAVE  :: MinRadialExponent
  REAL(DOUBLE), SAVE  :: MinDensityExponent
  REAL(DOUBLE), SAVE  :: PrimPairDistanceThreshold  ! Prim pairs
  REAL(DOUBLE), SAVE  :: PenetratDistanceThreshold  ! Penetration threshold
!
  REAL(DOUBLE),SAVE   :: AtomPairDistanceThreshold  ! General Atom Pairs
  TYPE(DBL_RNK2),SAVE :: ABDistanceThreshold        ! AB Atom Pairs
  CONTAINS
!====================================================================================================
!    Set and load global threholding parameters
!====================================================================================================
     SUBROUTINE SetThresholds(CurBase) 
         INTEGER          :: NExpt,Lndx
         TYPE(DBL_VECT)   :: Expts
         CHARACTER(LEN=*) :: CurBase
!        Get the primary thresholds
         CALL Get(Thresholds,Tag_O=CurBase)
!        Get distribution exponents
         CALL Get(NExpt,'nexpt',Tag_O=CurBase)
         CALL Get(Lndx ,'lndex',Tag_O=CurBase)
         CALL New(Expts,NExpt)
         CALL Get(Expts,'dexpt',Tag_O=CurBase)
!        Xi=Zeta*Zeta/(Zeta+Zeta)=Zeta_Min/2
         MinDensityExponent=Half*Expts%D(1)
!        Xi_Min=Zeta*Zeta/(Zeta+Zeta)=Zeta_Min/2
         MinRadialExponent=Half*Expts%D(1)
!        L_min = Lnsx/2
         MinRadialAngSym = Lndx
!        Dlete Exponents
         CALL Delete(Expts)
!        Set Atom-Atom thresholds
         CALL SetAtomPairThresh(Thresholds%Dist,CurBase)
!        Set Prim-Prim thresholds
         CALL SetPrimPairThresh(Thresholds%Dist)
     END SUBROUTINE SetThresholds
!====================================================================================================
!    Set the Atom Pair Distance Threshhold: 
!====================================================================================================
     SUBROUTINE SetAtomPairThresh(Tau,CurBase)
       REAL(DOUBLE),INTENT(IN)         :: Tau
       CHARACTER(LEN=*)                :: CurBase
!
       TYPE(BSET)                      :: BS
       INTEGER                         :: I,J,CFA,CFB,PFA,PFB,LA,LB,ELL
       REAL(DOUBLE)                    :: ZA,ZB,Zeta
!
       CALL Get(BS,Tag_O=CurBase)
       CALL New(ABDistanceThreshold,(/BS%NKind,BS%NKind/))
       DO I=1,BS%NKind
          CFA = BS%NCFnc%I(I)
          PFA = BS%NPFnc%I(CFA,I)
          LA  = BS%ASymm%I(2,CFA,I)
          ZA  = BS%Expnt%D(PFA,CFA,I)
          DO J=1,BS%NKind
             CFB = BS%NCFnc%I(J)
             PFB = BS%NPFnc%I(CFB,J)
             LB  = BS%ASymm%I(2,CFB,J)
             ZB  = BS%Expnt%D(PFB,CFB,J)
             Ell  = LA+LB
             Zeta = ZA*ZB/(ZA+ZB)
             ABDistanceThreshold%D(I,J) = AtomPairExtent(Ell,Zeta,Tau)
             ABDistanceThreshold%D(I,J) = ABDistanceThreshold%D(I,J)**2
             AtomPairDistanceThreshold  = MAX(AtomPairDistanceThreshold,ABDistanceThreshold%D(I,J))
          ENDDO
       ENDDO
!
     END SUBROUTINE SetAtomPairThresh
!====================================================================================================
!    Test the Distance for atom pair
!====================================================================================================
     FUNCTION TestAtomPair(Pair)
        LOGICAL                   :: TestAtomPair
        TYPE(AtomPair)            :: Pair
!        IF(Pair%AB2 > AtomPairDistanceThreshold) THEN
        IF(Pair%AB2 > ABDistanceThreshold%D(Pair%KA,Pair%KB)) THEN
           TestAtomPair = .FALSE.
        ELSE
           TestAtomPair = .TRUE.
        ENDIF
     END FUNCTION TestAtomPair
!====================================================================================================
!
!====================================================================================================
     FUNCTION AtomPairExtent(L,Zeta,Tau)
       INTEGER                         :: L
       REAL(DOUBLE)                    :: Tau,Zeta,Coef 
       REAL(DOUBLE)                    :: AtomPairExtent
!
       Coef = SQRT(DBLE(1+L))
       AtomPairExtent = SQPLog(L,Coef,Tau,Zeta)
       RETURN
!
     END FUNCTION AtomPairExtent
!====================================================================================================
!    Set the Primitive Pair Distance Threshhold: Xi_ab*Min*|A-B|^2 > -Log(Tau)
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
!====================================================================================================
!    Compute the extent of a Gaussian with exponent Zeta and amplitude Amp:
!    Amp*Exp[-Zeta*Extent^2] > Tau  
!====================================================================================================
     SUBROUTINE SetPenetrationThresh(Tau,ExpSwitch_O)
        REAL(DOUBLE),         INTENT(IN) :: Tau
        REAL(DOUBLE),OPTIONAL,INTENT(IN) :: ExpSwitch_O
        PenetratDistanceThreshold=-LOG(Tau)
        IF(PRESENT(ExpSwitch_O)) &
           PenetratDistanceThreshold=MIN(PenetratDistanceThreshold,ExpSwitch_O)
     END SUBROUTINE SetPenetrationThresh
!====================================================================================================
!    
!====================================================================================================
     FUNCTION GaussianExtent(Zeta,Amp)
        REAL(DOUBLE),INTENT(IN) :: Zeta,Amp
        REAL(DOUBLE)            :: GaussianExtent
        GaussianExtent=SQRT(MAX(1.D-10,PenetratDistanceThreshold+Amp)/Zeta)
     END FUNCTION GaussianExtent
!===================================================================================================
!     Recursive bisection to determine largest extent for this distribution
!     outside of which its contribution to the density and gradient is less than Tau
!===================================================================================================
     FUNCTION Extent(Ell,Zeta,HGTF,Tau,ExtraEll_O,Potential_O,ConvergeTo_O) RESULT (R)
       INTEGER                         :: Ell,ExtraEll
       REAL(DOUBLE)                    :: Zeta,Tau,Coef
       REAL(DOUBLE),DIMENSION(:)       :: HGTF
       INTEGER,OPTIONAL                :: ExtraEll_O
       LOGICAL,OPTIONAL                :: Potential_O
       REAL(DOUBLE),OPTIONAL           :: ConvergeTo_O
       REAL(DOUBLE),DIMENSION(0:HGEll) :: Co,HGPot
       REAL(DOUBLE)                    :: FUN,F0,F1
       INTEGER                         :: J,L,K,M,N,LMN,SN,LL
       REAL(DOUBLE)                    :: ConvergeTo,RMIN,RMAX,R,RErr,X
       LOGICAL                         :: Potential
!
       IF(PRESENT(ConvergeTo_O))THEN
          ConvergeTo=ConvergeTo_O
       ELSE 
          ConvergeTo=1.0D-8
       ENDIF
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
!      Take the spherical average of HGTF coefficients    
       DO L=0,Ell
          Co(L)=Zero
          DO LMN=LBegin(L),LEnd(L)
             Co(L)=Co(L)+HGTF(LMN)**2
          ENDDO
          Co(L)=SQRT(Co(L))
       ENDDO
!
!      Compute extent of a Hermite Gaussian overlap
!
       IF(.NOT. Potential )THEN
          Coef = Zero
          DO L=0,Ell
             Coef = Coef +Co(L)**2
          ENDDO
          Coef = SQRT(Coef)
          L = Ell+ExtraEll
          R = SQPLog(L,Coef,Tau,Zeta)
          RETURN
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
     END FUNCTION Extent
!===================================================================================
!    Solution to (2*Zeta*X)^L Exp[-Zeta*X*X] == tau 
!===================================================================================
     FUNCTION SQPlog(L,Coef,tau,Zeta)
       INTEGER           :: L
       REAL(DOUBLE)           :: Coef,tau,Zeta,X,OZ,OZZ,SQPlog
       REAL(DOUBLE),PARAMETER :: Xmax = 0.36787944117144232D0
!
       IF(Coef == Zero) THEN
          SQPLog = Zero
          RETURN
       ENDIF
!
       SELECT CASE(L)
       CASE(0)
          X = tau/Coef
          IF(X > One) THEN
             SQPLog = Zero
          ELSE
             SQPLog = SQRT(-LOG(X)/Zeta)
          ENDIF
       CASE(1)
          OZ  = 0.50D0/Zeta
          X   = tau/Coef
          X   = OZ*X*X
          IF(X > Xmax) THEN
             SQPLog = Zero
          ELSE
             SQPLog = SQRT(OZ*ProductLog1(X))
          ENDIF
       CASE(2)
          OZ  = 1.00D0/Zeta
          OZZ = 0.25D0*OZ
          X   = tau/Coef
          X   = OZZ*X
          IF(X > Xmax) THEN
             SQPLog = Zero
          ELSE
             SQPLog = SQRT(OZ*ProductLog1(X))
          ENDIF
       CASE(3)
          OZ  = 1.5D0/Zeta
          OZZ = 0.11111111111111111D0*OZ
          X   = tau/Coef
          X   = OZZ*(X**(0.6666666666666666D0))
          IF(X > Xmax) THEN
             SQPLog = Zero
          ELSE
             SQPLog = SQRT(OZ*ProductLog1(X))
          ENDIF
       CASE(4)
          OZ  = 2.0D0/Zeta
          OZZ = 0.0625D0*OZ
          X   = tau/Coef
          X   = OZZ*SQRT(X)
          IF(X > Xmax) THEN
             SQPLog = Zero
          ELSE
             SQPLog = SQRT(OZ*ProductLog1(X))
          ENDIF
       CASE(5:)
          OZ  = Half*DBLE(L)/Zeta
          OZZ = Half/(DBLE(L)*Zeta)
          X   = tau/Coef
          X   = OZZ*(X**(Two/DBLE(L)))
          IF(X > Xmax) THEN
             SQPLog = Zero
          ELSE
             SQPLog = SQRT(OZ*ProductLog1(X))
          ENDIF
       END SELECT
!
     END FUNCTION SQPlog
!===================================================================================
!    Norm*Gamma[L/2+3/2,Zeta*R*R]
!===================================================================================
     FUNCTION AsymHGTF(L,Zeta,R)
       INTEGER                      :: L,N,LL
       REAL(DOUBLE)                 :: Zeta,R,R0,TwoZR,AsymHGTF,ZRR
       REAL(DOUBLE),PARAMETER       :: A0=1.00000000000000000D0
       REAL(DOUBLE),PARAMETER       :: A1=1.00000000000000000D0
       REAL(DOUBLE),PARAMETER       :: A2=1.35914091422952261D0
       REAL(DOUBLE),PARAMETER       :: A3=1.19035469743248456D0
       REAL(DOUBLE),PARAMETER       :: A4=1.38544801854949691D0
!
       R0  = MAX(R,SQRT(Half*DBLE(L)/Zeta))
       ZR  = Zeta*R0
       ZRR = ZR*R0
!
       SELECT CASE(L)
       CASE(0)
          AsymHGTF = EXP(-ZRR)
       CASE(1)
          AsymHGTF = Two*ZR*EXP(-ZRR)
       CASE(2)
          TwoZR    = Two*ZR
          AsymHGTF = A2*TwoZR*TwoZR*EXP(-ZRR)
       CASE(3)
          TwoZR    = Two*ZR
          AsymHGTF = A3*TwoZR*TwoZR*TwoZR*EXP(-ZRR)
       CASE(4:)
          TwoZR    = Two*ZR
          AsymHGTF = (TwoZR**L)*EXP(-ZRR)
       END SELECT
! 
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
!====================================================================================================
!     COMPUTE FUNCTIONS THAT RETURN THE ARGUMENT T TO THE GAMMA FUNCTIONS F[m,T]
!     THAT RESULT FROM USING THE THE MULTIPOLE APPROXIMATION TO WITHIN A SPECIFIED 
!     ERROR:  THESE FUNCTIONS GIVE T (THE PAC) FOR A GIVEN PENETRATION ERROR.
!====================================================================================================
     FUNCTION PFunk(Ell,Ack)
        INTEGER                    :: Ell
        REAL(DOUBLE)               :: PFunk,X,Ack,MinAcc,MaxAcc
        REAL(DOUBLE),DIMENSION(30) :: W
        INCLUDE 'GammaDimensions.Inc'
        X=-LOG(Ack)
        INCLUDE 'PFunk.Inc'
        PFunk=MIN(PFunk,Gamma_Switch)
     END FUNCTION PFunk    
END MODULE
!
!!$     FUNCTION AsymHGTF(L,Zeta,R)
!!$       INTEGER                      :: L,N,LL
!!$       REAL(DOUBLE)                 :: Zeta,R,AsymHGTF,Norm
!!$       REAL(DOUBLE),DIMENSION(0:16) :: NFactor = (/ 1.00000000000000000D0,&
!!$                              1.12837916709551257D0,2.00000000000000000D0,&
!!$                              4.51351666838205030D0,1.20000000000000000D1,&
!!$                              3.61081333470564024D1,1.20000000000000000D2,&
!!$                              4.33297600164676828D2,1.68000000000000000D3,&
!!$                              6.93276160263482925D3,3.02400000000000000D4,&
!!$                              1.38655232052696585D5,6.65280000000000000D5,&
!!$                              3.32772556926471804D6,1.72972800000000000D7,&
!!$                              9.31763159394121052D7,5.18918400000000000D8 /)
!!$!
!!$       N    = (-1)**L
!!$       Norm = (Zeta**(Half*DBLE(L)))*NFactor(L)
!!$       LL   = (L+3)/2
!!$       IF(L == 0) THEN
!!$          AsymHGTF = EXP(-Zeta*R*R)
!!$       ELSE
!!$          IF(N == 1) THEN
!!$             AsymHGTF = Norm*GammaHalf(LL,Zeta*R*R)
!!$          ELSE
!!$             AsymHGTF = Norm*GammaOne(LL,Zeta*R*R)
!!$          ENDIF
!!$       ENDIF
!!$     END FUNCTION AsymHGTF
!
!!$          RMIN = Zero
!!$          RMAX = SQRT(EXP_SWITCH/Zeta)
!!$!
!!$          F0 = Zero
!!$          DO L=ExtraELL,Ell+ExtraEll
!!$             F0 = F0 + Co(L-ExtraEll)*AsymHGTF(L,Zeta,RMIN)
!!$          ENDDO
!!$          IF(F0 < Zero) THEN
!!$             CALL MondoHalt(-100,'F0 < 0')
!!$          ENDIF
!!$          IF(F0 < Tau) THEN
!!$             R = Zero
!!$             GOTO 100
!!$          ENDIF
!!$!
!!$          F1 = Zero
!!$          DO L=ExtraELL,Ell+ExtraEll
!!$             F1 = F1 + Co(L-ExtraEll)*AsymHGTF(L,Zeta,RMAX)
!!$          ENDDO
!!$          IF(F1>Tau) THEN
!!$             R = RMAX
!!$             GOTO 100
!!$          ENDIF
!!$!
!!$          R = Half*(RMIN+RMAX)
!!$          DO J=1,200
!!$             FUN  = Zero
!!$             DO L=ExtraEll,Ell+ExtraEll             
!!$                FUN  = FUN +Co(L-ExtraEll)*AsymHGTF(L,Zeta,R)
!!$             ENDDO
!!$             FUN = FUN-Tau
!!$             IF(FUN < Zero) THEN
!!$                RMAX = R
!!$             ELSEIF(FUN > Zero) THEN
!!$                RMIN = R
!!$             ENDIF
!!$             RErr = ABS(R-Half*(RMIN+RMAX))
!!$             R = Half*(RMIN+RMAX)
!!$             IF(RErr < ConvergeTo) GOTO 100
!!$          ENDDO
!!$          CALL MondoHalt(-100,'Overlap did not converge in 200 iterations')
!!$100       CONTINUE
!
!!$       RMIN = SQRT(Half*DBLE(L)/Zeta)
!!$       RMAX = SQRT(EXP_SWITCH/Zeta)
!!$!
!!$       F1 = AsymHGTF(L,Zeta,RMAX)
!!$       IF(F1>Tau) THEN
!!$          AtomPairExtent = RMAX
!!$          GOTO 100
!!$       ENDIF
!!$       F0 = AsymHGTF(L,Zeta,RMIN)
!!$       IF(F0<Tau) THEN
!!$          AtomPairExtent = Zero
!!$          GOTO 100
!!$       ENDIF
!!$!
!!$       AtomPairExtent = Half*(RMIN+RMAX)
!!$       DO J=1,200
!!$          FUN = AsymHGTF(L,Zeta,AtomPairExtent)-Tau
!!$          IF(FUN < Zero) THEN
!!$             RMAX = AtomPairExtent
!!$          ELSEIF(FUN > Zero) THEN
!!$             RMIN = AtomPairExtent
!!$          ENDIF
!!$          RErr = ABS(AtomPairExtent - Half*(RMIN+RMAX))
!!$          AtomPairExtent = Half*(RMIN+RMAX) 
!!$          IF(RErr < 1.D-8) GOTO 100
!!$       ENDDO
!!$       CALL MondoHalt(-100,'Overlap did not converge in 200 iterations')
!!$100    CONTINUE
