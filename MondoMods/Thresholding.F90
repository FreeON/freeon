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
!    Author: Matt Challacombe
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
!     Recursive bisection to determine largest extent for this distribution, outside
!     outside of which its contribution to the density and gradient is less than Tau
!===================================================================================
      FUNCTION Extent(Ell,Zeta,HGTF,Tau,Digits_O,ExtraEll_O,Potential_O) RESULT (R)
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
            ConvergeTo=10.D0**(-4)
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
!        Do straight overlap type extent 
         IF(.NOT.PRESENT(Potential_O))THEN
            DelR=SQRT(EXP_SWITCH/Zeta)*1.5D0
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
               ENDIF
!              If still to the left, increment bisection point
               IF(CTest>Zero)BisR=R
            ENDDO
            IF(MidRho>Tau)CALL Halt(' Faild in Extent of Overlap ')
          ELSE
            RTE=Zeta*NuclearExpnt
            RPE=Zeta+NuclearExpnt
            Omega=RTE/RPE 
            Upq=TwoPi5x2/(RTE*SQRT(RPE)) &
               *(NuclearExpnt/Pi)**(ThreeHalves) ! add on moment for delta function...
            DelR=SQRT(GAMMA_SWITCH/Zeta)*Two
            LTot=Ell+ExtraEll
            BisR=Zero
!            WRITE(*,*)' Zeta  = ',Zeta
!            WRITE(*,*)' Omega = ',Omega
!            WRITE(*,*)' UPQ   = ',Upq
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
!               WRITE(*,*)R,MidRho
!              Convergence test
               CTest=(MidRho-Tau)/Tau
               IF(R<1.D-30)THEN
                  R=Zero
                  RETURN
               ELSEIF(ABS(CTest)<ConvergeTo)THEN
                  RETURN
               ENDIF
!              If still to the left, increment bisection point
               IF(CTest>Zero)BisR=R
            ENDDO
            IF(MidRho>Tau)CALL Halt(' Faild in Extent of Potential ')
         ENDIF
      END FUNCTION Extent
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
