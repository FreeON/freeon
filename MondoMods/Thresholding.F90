!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
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
         INTEGER          :: NExpt
         TYPE(DBL_VECT)   :: Expts
         CHARACTER(LEN=*) :: CurBase
!        Get the primary thresholds
         CALL Get(Thresholds,Tag_O=CurBase)
!        Get distribution exponents
         CALL Get(NExpt,'nexpt',Tag_O=CurBase)
         CALL New(Expts,NExpt)
         CALL Get(Expts,'dexpt',Tag_O=CurBase)
         MinDensityExponent=Half*Expts%D(1)
!        Xi_Min=Zeta*Zeta/(Zeta+Zeta)=Zeta_Min/2
         MinRadialExponent=Half*Expts%D(1)
         CALL Delete(Expts)
!        Set Atom-Atom thresholds
         CALL SetAtomPairThresh(Thresholds%Dist)
!        Set Prim-Prim thresholds
         CALL SetPrimPairThresh(Thresholds%Dist)
     END SUBROUTINE SetThresholds
!====================================================================================================
!    Set the Atom Pair Distance Threshhold: Zeta_Min*|A-B|^2 > -Log(Tau)
!====================================================================================================
     SUBROUTINE SetAtomPairThresh(Tau)
        REAL(DOUBLE),INTENT(IN) :: Tau
        AtomPairDistanceThreshold=-LOG(Tau)/MinRadialExponent
     END SUBROUTINE SetAtomPairThresh
!
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
!
     FUNCTION GaussianExtent(Zeta,Amp)
        REAL(DOUBLE),INTENT(IN) :: Zeta,Amp
        REAL(DOUBLE)            :: GaussianExtent
        GaussianExtent=SQRT(MAX(1.D-10,PenetratDistanceThreshold+Amp)/Zeta)
     END FUNCTION GaussianExtent


!===================================================================================
!     Recursive bisection to determine largest extent for this distribution
!     outside of which its contribution to the density and gradient is less than Tau
!===================================================================================
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
          ConvergeTo=10.D0**(-4)
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
!
!      Take the spherical average of HGTF coefficients      
       DO L=0,Ell
          Co(L)=Zero
          DO LMN=LBegin(L),LEnd(L)
             Co(L)=Co(L)+HGTF(LMN)**2
          ENDDO
          Co(L)=SQRT(Co(L))
       ENDDO
!      Compute extent of a Hermite Gaussian
       IF(.NOT. Potential )THEN
          RMIN = Zero
          RMAX = SQRT(EXP_SWITCH/Zeta)
          R = Half*(RMIN+RMAX)
          F0 = Zero
          F1 = Zero
          DO L=ExtraELL,Ell+ExtraEll
             F0 = F0 + Co(L-ExtraEll)*AsymHGTF(L,Zeta,RMIN)
             F1 = F1 + Co(L-ExtraEll)*AsymHGTF(L,Zeta,RMAX)
          ENDDO
          IF(F0 < Zero) THEN
             CALL MondoHalt(-100,'F0 < 0')
          ENDIF
          IF(F0 < Tau) THEN
             R = Zero
          ELSEIF((F1-Tau) > Zero) THEN
             R = RMAX
          ELSE 
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
100          CONTINUE
          ENDIF
!      Do a Potential overlap
       ELSEIF( Potential ) THEN
          RMIN = 1.D-8
          RMAX = SQRT(GAMMA_SWITCH/Zeta)
          R = Half*(RMIN+RMAX)
          F0 = Zero
          HGPot = AsymPot(Ell+ExtraEll,Zeta,RMIN)
          DO L=ExtraELL,Ell+ExtraEll
             F0 = F0 + Co(L-ExtraEll)*ABS(HGPot(L))
          ENDDO
          F1 = Zero
          HGPot = AsymPot(Ell+ExtraEll,Zeta,RMAX)
          DO L=ExtraELL,Ell+ExtraEll
             F1 = F1 + Co(L-ExtraEll)*ABS(HGPot(L))
          ENDDO
          IF(F0 < Zero) THEN
             CALL MondoHalt(-100,'F0 < 0')
          ENDIF
          IF(F0 < Tau) THEN
             R = Zero
          ELSEIF((F1-Tau) > Zero) THEN
             R = RMAX
          ELSE 
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
200          CONTINUE            
          ENDIF             
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
!
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
!===================================================================================
!     Recursive bisection to determine largest extent for this distribution, outside
!     outside of which its contribution to the density and gradient is less than Tau
!===================================================================================
      FUNCTION Extent_Old(Ell,Zeta,HGTF,Tau,Digits_O,ExtraEll_O,Potential_O) RESULT (R)
         INTEGER                         :: Ell
         REAL(DOUBLE)                    :: Zeta,Tau
         REAL(DOUBLE),DIMENSION(:)       :: HGTF
         INTEGER,OPTIONAL                :: Digits_O,ExtraEll_O
         LOGICAL,OPTIONAL                :: Potential_O
         REAL(DOUBLE),DIMENSION(0:HGEll) :: Co,ErrR
         REAL(DOUBLE),DIMENSION(0:HGEll, &
                                0:HGEll) :: HGErr
         REAL(DOUBLE),DIMENSION(0:20)    :: LambdaR
         REAL(DOUBLE)                    :: R
         INTEGER                         :: J,L,K,M,N,LMN,ExtraEll,LTot
         REAL(DOUBLE)                    :: ConvergeTo,R2,DelR,BisR,CTest, &
                                            RhoR,dRho,MidRho,Xpt,TwoZ, &
                                            Omega,RPE,RTE,T
         LOGICAL                         :: PotentialQ
!----------------------------------------------------------------------------------
         IF(PRESENT(Digits_O))THEN
            ConvergeTo=10.D0**(-Digits_O)
         ELSE 
            ConvergeTo=10.D0**(-5)
         ENDIF
         IF(PRESENT(ExtraEll_O))THEN
            ExtraEll=ExtraEll_O
         ELSE
            ExtraEll=1
         ENDIF
!        Take the spherical average of HGTF coefficients      
         DO L=0,Ell
            Co(L)=Zero
            DO LMN=LBegin(L),LEnd(L)
               Co(L)=Co(L)+HGTF(LMN)**2
            ENDDO
            Co(L)=SQRT(Co(L))
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
!            WRITE(*,*)' Zeta  = ',Zeta
!            WRITE(*,*)' Omega = ',Omega
!            WRITE(*,*)' UPQ   = ',Upq
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
!               WRITE(*,*)BisR,' MidRho = ',MidRho
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
!            WRITE(*,*)' Zeta  = ',Zeta
!            WRITE(*,*)' Omega = ',Omega
!            WRITE(*,*)' UPQ   = ',Upq
!            WRITE(*,*)' Ell   = ',Ell
!            WRITE(*,*)' Co    = ',Co(0)
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
      END FUNCTION Extent_old
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
