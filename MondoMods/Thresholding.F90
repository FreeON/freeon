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
  REAL(DOUBLE), SAVE  :: AtomPairDistanceThreshold  ! Atom pairs
  REAL(DOUBLE), SAVE  :: PrimPairDistanceThreshold  ! Prim pairs
  REAL(DOUBLE), SAVE  :: PenetratDistanceThreshold  ! Penetration threshold
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
         MinDensityExponent=Half*Expts%D(1)
!        Xi_Min=Zeta*Zeta/(Zeta+Zeta)=Zeta_Min/2
         MinRadialExponent=Half*Expts%D(1)
!        L_min = Lnsx/2
         MinRadialAngSym = Lndx
!        Dlete Exponents
         CALL Delete(Expts)
!        Set Atom-Atom thresholds
         CALL SetAtomPairThresh(Thresholds%Dist)
!        Set Prim-Prim thresholds
         CALL SetPrimPairThresh(Thresholds%Dist)
     END SUBROUTINE SetThresholds
!====================================================================================================
!    Set the Atom Pair Distance Threshhold: 
!====================================================================================================
     SUBROUTINE SetAtomPairThresh(Tau)
       INTEGER                         :: J,L
       REAL(DOUBLE),INTENT(IN)         :: Tau
       REAL(DOUBLE)                    :: RMIN,RMAX,F0,F1,R,FUN,RErr
!
       RMIN = SQRT(-LOG(Tau)/MinRadialExponent)
       RMAX = 10.D0*RMIN
!
       F0 = AsymHGTF(MinRadialAngSym+1,MinRadialExponent,RMIN)
       F1 = AsymHGTF(MinRadialAngSym+1,MinRadialExponent,RMAX)
       IF(F1>Tau) THEN
          R = RMAX
          GOTO 100
       ENDIF
!
       R = Half*(RMIN+RMAX)
       DO J=1,200
          FUN = AsymHGTF(MinRadialAngSym+1,MinRadialExponent,R)-Tau
          IF(FUN < Zero) THEN
             RMAX = R
          ELSEIF(FUN > Zero) THEN
             RMIN = R
          ENDIF
          RErr = ABS(R-Half*(RMIN+RMAX))
          R = Half*(RMIN+RMAX) 
          IF(RErr < 1.D-8) GOTO 100
       ENDDO
       CALL MondoHalt(-100,'Overlap did not converge in 200 iterations')
100    CONTINUE
       AtomPairDistanceThreshold=R*R
!
     END SUBROUTINE SetAtomPairThresh
!====================================================================================================
!    Test the Distance for atom pair
!====================================================================================================
     FUNCTION TestAtomPair(Pair)
        LOGICAL                   :: TestAtomPair
        TYPE(AtomPair)            :: Pair
        IF(Pair%AB2>AtomPairDistanceThreshold) THEN
           TestAtomPair = .FALSE.
        ELSE
           TestAtomPair = .TRUE.
        ENDIF
     END FUNCTION TestAtomPair
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
     FUNCTION Extent(Ell,Zeta,HGTF,Tau,Digits_O,ExtraEll_O,Potential_O) RESULT (R)
       INTEGER                         :: Ell,ExtraEll
       REAL(DOUBLE)                    :: Zeta,Tau
       REAL(DOUBLE),DIMENSION(:)       :: HGTF
       INTEGER,OPTIONAL                :: Digits_O,ExtraEll_O
       LOGICAL,OPTIONAL                :: Potential_O
       REAL(DOUBLE),DIMENSION(0:HGEll) :: Co,HGPot
       REAL(DOUBLE)                    :: FUN,F0,F1
       INTEGER                         :: J,L,K,M,N,LMN,SN,LL
       REAL(DOUBLE)                    :: ConvergeTo,RMIN,RMAX,R,RErr
       LOGICAL                         :: Potential
!
       IF(PRESENT(Digits_O))THEN
          ConvergeTo=10.D0**(-Digits_O)
       ELSE 
          ConvergeTo=10.D0**(-8)
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
     END FUNCTION Extent
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
!!$!===================================================================================
!!$!     Recursive bisection to determine largest extent for this distribution
!!$!     outside of which its contribution to the density and gradient is less than Tau
!!$!===================================================================================
!!$     FUNCTION Extent(LLow,LHig,LShift,Zeta,HGTF,Tau,Digits_O,Potential_O) RESULT (R)
!!$       INTEGER                          :: LLow,LHig,LSft
!!$       REAL(DOUBLE)                     :: Zeta,Tau
!!$       REAL(DOUBLE),DIMENSION(:)        :: HGTF
!!$       INTEGER,OPTIONAL                 :: Digits_O
!!$       LOGICAL,OPTIONAL                 :: Potential_O
!!$       REAL(DOUBLE),DIMENSION(0:HGEll):: Co,HGPot
!!$       REAL(DOUBLE)                     :: FUN,F0,F1
!!$       INTEGER                          :: J,L,K,M,N,LMN,SN,LL
!!$       REAL(DOUBLE)                     :: ConvergeTo,RMIN,RMAX,R,RErr
!!$       LOGICAL                          :: Potential
!!$!
!!$       IF(PRESENT(Digits_O))THEN
!!$          ConvergeTo=10.D0**(-Digits_O)
!!$       ELSE 
!!$          ConvergeTo=10.D0**(-8)
!!$       ENDIF
!!$!
!!$       IF(PRESENT(Potential_O)) THEN
!!$          Potential = Potential_O
!!$       ELSE
!!$          Potential = .FALSE.
!!$       ENDIF
!!$!      Take the spherical average of HGTF coefficients      
!!$       DO L = LLow,LHig
!!$          Co(L)=Zero
!!$          DO LMN=LBegin(L),LEnd(L)
!!$             Co(L)=Co(L)+HGTF(LMN)**2
!!$          ENDDO
!!$          Co(L)=SQRT(Co(L))
!!$       ENDDO
!!$!
!!$!      Compute extent of a Hermite Gaussian overlap
!!$!
!!$       IF(.NOT. Potential )THEN
!!$          RMIN = Zero
!!$          RMAX = SQRT(EXP_SWITCH/Zeta)
!!$!
!!$          F0 = Zero
!!$          DO L = LLow,LHig
!!$             F0 = F0 + Co(L)*AsymHGTF(L+LSft,Zeta,RMIN)
!!$          ENDDO
!!$          IF(F0 < Zero) THEN
!!$             CALL MondoHalt(-100,'F0 < 0')
!!$          ENDIF
!!$          IF(F0 < Tau) THEN
!!$             R = Zero
!!$             RETURN
!!$          ENDIF
!!$!
!!$          F1 = Zero
!!$          DO L = LLow,LHig
!!$             F1 = F1 + Co(L)*AsymHGTF(L+LSft,Zeta,RMAX)
!!$          ENDDO
!!$          IF(F1>Tau) THEN
!!$             R = RMAX
!!$             RETURN
!!$          ENDIF
!!$!
!!$          R = Half*(RMIN+RMAX)
!!$          DO J=1,200
!!$             FUN  = Zero
!!$             DO L = LLow,LHig       
!!$                FUN  = FUN +Co(L)*AsymHGTF(L+LSft,Zeta,R)
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
!!$!
!!$!      Do a Potential overlap
!!$!
!!$       ELSEIF( Potential ) THEN
!!$          RMIN = 1.D-14
!!$          RMAX = SQRT(GAMMA_SWITCH/Zeta)
!!$!
!!$          F0 = Zero
!!$          HGPot = AsymPot(LHig+LSft,Zeta,RMIN)
!!$          DO L = LLow,LHig
!!$             F0 = F0 + Co(L)*ABS(HGPot(L+LShift))
!!$          ENDDO
!!$          IF(F0 < Zero) THEN
!!$             CALL MondoHalt(-100,'F0 < 0')
!!$          ENDIF
!!$          IF(F0 < Tau) THEN
!!$             R = Zero
!!$             RETURN
!!$          ENDIF
!!$!
!!$          F1 = Zero
!!$          HGPot = AsymPot(LHig+LSft,Zeta,RMAX)
!!$          DO L = LLow,LHig
!!$             F1 = F1 + Co(L)*ABS(HGPot(L+LSft))
!!$          ENDDO
!!$          IF(F1>Tau) THEN
!!$             R = RMAX
!!$             RETURN
!!$          ENDIF
!!$!
!!$          R = Half*(RMIN+RMAX)
!!$          DO J=1,200
!!$             FUN  = Zero
!!$             HGPot = AsymPot(LHig+LSft,Zeta,R)
!!$             DO L = LLow,LHig
!!$                FUN  = FUN +Co(L)*ABS(HGPot(L+LSft))
!!$             ENDDO
!!$             FUN = FUN-Tau
!!$             IF(FUN < Zero) THEN
!!$                RMAX = R
!!$             ELSEIF(FUN > Zero) THEN
!!$                RMIN = R
!!$             ENDIF
!!$             RErr = ABS(R-Half*(RMIN+RMAX))
!!$             R = Half*(RMIN+RMAX)
!!$             IF(RErr < ConvergeTo) GOTO 200
!!$          ENDDO
!!$          CALL MondoHalt(-100,'Potential Overlap did not converge in 200 iterations')
!!$200       CONTINUE                     
!!$       ENDIF
!!$     END FUNCTION Extent
