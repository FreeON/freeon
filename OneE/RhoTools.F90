!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE RhoTools
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  USE AtomPairs
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
!---------------------------------------------------------------------------------------
!  Density in a Hermite Gaussian basis
!---------------------------------------------------------------------------------------
  TYPE HGRho_new
     INTEGER          :: NSDen   !-- Number of Spin densities
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: NDist   !-- Number of distributions
     INTEGER          :: NCoef   !-- Number of coefficients
     TYPE(INT_VECT)   :: OffCo   !-- Coefficient offset     (NDist)
     TYPE(INT_VECT)   :: Ell     !-- Angular symmetry       (NDist)
     TYPE(DBL_VECT)   :: Zeta    !-- Distriution Exponent   (NDist)
     TYPE(DBL_VECT)   :: Qx      !-- x-coordinates          (NDist)
     TYPE(DBL_VECT)   :: Qy      !-- y-coordinates          (NDist)
     TYPE(DBL_VECT)   :: Qz      !-- z-coordinates          (NDist)
     TYPE(DBL_VECT)   :: Co      !-- Density coefficients   (NCoef)
  END TYPE HGRho_new
!
  CONTAINS
!---------------------------------------------------------------------------------------
! Allocate the Density
!---------------------------------------------------------------------------------------
  SUBROUTINE New_HGRho_new(A,N_O)
    TYPE(HGRho_new)                 :: A
    INTEGER,OPTIONAL,DIMENSION(3)   :: N_O
!
    IF(AllocQ(A%Alloc)) THEN
       CALL MondoHalt(0,'Attemp to Allocate HGRho_new already allocated')
    ELSE
       A%Alloc=ALLOCATED_TRUE
       IF(PRESENT(N_O)) THEN
          A%NDist=N_O(1)
          A%NCoef=N_O(2)
          A%NSDen=N_O(3)
       ELSE
          A%NDist=0
          A%NCoef=0
          A%NSDen=0
       ENDIF
       CALL New(A%OffCo,A%NDist)
       CALL New(A%Ell  ,A%NDist)
       CALL New(A%Zeta ,A%NDist)
       CALL New(A%Qx   ,A%NDist)
       CALL New(A%Qy   ,A%NDist)
       CALL New(A%Qz   ,A%NDist)
       CALL New(A%Co   ,A%NCoef*A%NSDen)
!old       CALL New(A%Co   ,A%NCoef)
    ENDIF
!
  END SUBROUTINE New_HGRho_new
!---------------------------------------------------------------------------------------
! Delete the density
!---------------------------------------------------------------------------------------
  SUBROUTINE Delete_HGRho_new(A)
    TYPE(HGRho_new)            :: A
!
    IF(AllocQ(A%Alloc)) THEN
       CALL Delete(A%OffCo)
       CALL Delete(A%Ell)
       CALL Delete(A%Zeta)
       CALL Delete(A%Qx)
       CALL Delete(A%Qy)
       CALL Delete(A%Qz)
       CALL Delete(A%Co)
       A%Alloc=ALLOCATED_FALSE
    ELSE
       CALL MondoHalt(0,'Attemp to Deallocate HGRho_new not allocated')
    ENDIF
!
  END SUBROUTINE Delete_HGRho_new
!---------------------------------------------------------------------------------------
! Copy One density to Another
!---------------------------------------------------------------------------------------
  SUBROUTINE Copy_HGRho_new(A,B)
    TYPE(HGRho_new)                  :: A,B
    INTEGER                          :: I
!
    IF(AllocQ(A%Alloc)) THEN
       CALL Delete_HGRho_new(A)
    ENDIF
    CALL New_HGRho_new(A,(/B%NDist,B%NCoef,B%NSDen/))
    DO I=1,B%NDist
       A%OffCo%I(I) = B%OffCo%I(I)
       A%Ell%I(I)   = B%Ell%I(I)
       A%Zeta%D(I)  = B%Zeta%D(I)
       A%Qx%D(I)    = B%Qx%D(I)
       A%Qy%D(I)    = B%Qy%D(I)
       A%Qz%D(I)    = B%Qz%D(I)
    ENDDO

    DO I=1,B%NCoef*B%NSDen
       A%Co%D(I)   = B%Co%D(I)
    ENDDO
!
  END SUBROUTINE Copy_HGRho_new
!---------------------------------------------------------------------------------------
!  Integrate The Distibutions from NLow to NHig
!---------------------------------------------------------------------------------------
  FUNCTION Integrate_HGRho_new(Rho,iSDen,NLow,NHig) RESULT(RhoSum)
    TYPE(HGRho_new)                  :: Rho
    INTEGER                          :: iSDen,NLow,NHig,nd,OffCo,I0
    REAL(DOUBLE)                     :: Zeta,RhoSum
!
    I0=(iSDen-1)*Rho%NCoef !<<< SPIN
    RhoSum = Zero
    DO nd = NLow,NHig
       Zeta   = Rho%Zeta%D(nd)
       OffCo  = Rho%OffCo%I(nd)+I0 !<<< SPIN
       RhoSum = RhoSum + Rho%Co%D(OffCo)*((Pi/Zeta)**ThreeHalves)
    ENDDO
!
  END FUNCTION Integrate_HGRho_new
!---------------------------------------------------------------------------------------
!  Add Distributions of s type and with exponent zeta
!---------------------------------------------------------------------------------------
 SUBROUTINE AddDist(Rho,GM,Zeta,BNumDist,ENumDist)
    TYPE(CRDS)                    :: GM
    TYPE(HGRho_new)               :: Rho,RhoTmp
    INTEGER                       :: BNumDist,ENumDist,NumDist,I,iq,ir,AtA,iSMat,I0,I1
    REAL(DOUBLE)                  :: DDelta,Zeta
!
    CALL Copy_HGRho_new(RhoTmp,Rho)
    CALL Delete_HGRho_new(Rho)
!
    DDelta = Half*(Zeta/Pi)**(ThreeHalves)
    NumDist = ENumDist-BNumDist+1
    CALL New_HGRho_new(Rho,(/RhoTmp%NDist+NumDist,RhoTmp%NCoef+NumDist,RhoTmp%NSDen/))
!
    DO I=1,RhoTmp%NDist
       Rho%OffCo%I(I) = RhoTmp%OffCo%I(I)
       Rho%Ell%I(I)   = RhoTmp%Ell%I(I)
       Rho%Zeta%D(I)  = RhoTmp%Zeta%D(I)
       Rho%Qx%D(I)    = RhoTmp%Qx%D(I)
       Rho%Qy%D(I)    = RhoTmp%Qy%D(I)
       Rho%Qz%D(I)    = RhoTmp%Qz%D(I)
    ENDDO

    I0=1
    I1=1
    DO iSMat=1,RhoTmp%NSDen
       DO I=1,RhoTmp%NCoef
!old       DO I=1,RhoTmp%NCoef
!old          Rho%Co%D(I)    = RhoTmp%Co%D(I) !<<< SPIN
          Rho%Co%D(I1)    = RhoTmp%Co%D(I0) !<<< SPIN
          I0=I0+1
          I1=I1+1
       ENDDO
       I1=I1+NumDist
    ENDDO
!
    DO I = 1,NumDist
       AtA = I+BNumDist-1
       Iq = RhoTmp%NDist+I
       Ir = RhoTmp%NCoef+I
       Rho%OffCo%I(Iq) = Ir
       Rho%Ell%I(Iq)   = 0
       Rho%Zeta%D(Iq)  = Zeta
       Rho%Qx%D(Iq)    = GM%Carts%D(1,AtA)
       Rho%Qy%D(Iq)    = GM%Carts%D(2,AtA)
       Rho%Qz%D(Iq)    = GM%Carts%D(3,AtA)
       IF(GM%AtNum%D(AtA) < 105.D0) THEN
          Rho%Co%D(Ir)    =-GM%AtNum%D(AtA)*DDelta !<<< SPIN
       ELSE
          Rho%Co%D(Ir)    = Zero !<<< SPIN
       ENDIF
    ENDDO
    CALL Delete_HGRho_new(RhoTmp)
!
  END SUBROUTINE AddDist
!---------------------------------------------------------------------------------------
! Subroutine to Remove unnessesary distributions in rho
!---------------------------------------------------------------------------------------
  SUBROUTINE Prune_Rho_new(TOL,Rho)
    TYPE(HGRho_new)     :: Rho,RhoTmp
    INTEGER             :: I,NDist,NCoef,Ell,LenKet,OffCo1,OffCo2,iSMat,I1,I2
    REAL(DOUBLE)        :: Zeta,Mag,TOL
!
    CALL Copy_HGRho_new(RhoTmp,Rho)
    CALL Delete_HGRho_new(Rho)
!
!   Count the Number of Distribution whose Magnitude is above TOL
!
    NDist = 0
    NCoef = 0
    DO I=1,RhoTmp%NDist
       Ell    = RhoTmp%Ell%I(I)
       Zeta   = RhoTmp%Zeta%D(I)
       LenKet = LHGTF(Ell)
       OffCo1 = RhoTmp%OffCo%I(I)
       Mag    = MagDist(Ell,Zeta,RhoTmp%Co%D(OffCo1:OffCo1+LenKet-1)) !<<< SPIN
       IF(Mag > TOL) THEN
          NDist = NDist+1
          NCoef = NCoef+LenKet
       ENDIF
    ENDDO
!
!   Create New Rho
!
    CALL New_HGRho_new(Rho,(/NDist,NCoef,RhoTmp%NSDen/))
!
    NDist = 0
    NCoef = 0
    DO I=1,RhoTmp%NDist
       Ell    = RhoTmp%Ell%I(I)
       Zeta   = RhoTmp%Zeta%D(I)
       LenKet = LHGTF(Ell)
       OffCo1 = RhoTmp%OffCo%I(I)
       Mag    = MagDist(Ell,Zeta,RhoTmp%Co%D(OffCo1:OffCo1+LenKet-1))!<<< SPIN
       IF(Mag > TOL) THEN
          NDist = NDist+1
          NCoef = NCoef+LenKet
          OffCo2= NCoef-LenKet+1
          Rho%Ell%I(NDist)   = RhoTmp%Ell%I(I)
          Rho%Zeta%D(NDist)  = RhoTmp%Zeta%D(I)
          Rho%Qx%D(NDist)    = RhoTmp%Qx%D(I)
          Rho%Qy%D(NDist)    = RhoTmp%Qy%D(I)
          Rho%Qz%D(NDist)    = RhoTmp%Qz%D(I)
          Rho%OffCo%I(NDist) = OffCo2
!old          Rho%Co%D(OffCo2:OffCo2+LenKet-1) = RhoTmp%Co%D(OffCo1:OffCo1+LenKet-1)!<<< SPIN
          DO iSMat=1,RhoTmp%NSDen
             I1=OffCo1+(iSMat-1)*RhoTmp%NCoef
             I2=OffCo2+(iSMat-1)*Rho%NCoef
             CALL DCOPY(LenKet,RhoTmp%Co%D(I1),1,Rho%Co%D(I2),1)
          ENDDO
       ENDIF
    ENDDO
    CALL Delete_HGRho_new(RhoTmp)
!
  END SUBROUTINE Prune_Rho_new
!---------------------------------------------------------------------------------------
! Collect Distributuins of a perticulate Ell type
!---------------------------------------------------------------------------------------
  SUBROUTINE PruneEll_Rho_new(EllKeep,Rho)
    TYPE(HGRho_new)     :: Rho,RhoTmp
    INTEGER             :: EllKeep,I,NDist,NCoef,Ell,LenKet,OffCo1,OffCo2
    REAL(DOUBLE)        :: Zeta,Mag
!
    CALL Copy_HGRho_new(RhoTmp,Rho)
    CALL Delete_HGRho_new(Rho)
!
!   Count the Number of Distribution with Ell = EllKeep
!
    NDist = 0
    NCoef = 0
    DO I=1,RhoTmp%NDist
       Ell    = RhoTmp%Ell%I(I)
       LenKet = LHGTF(Ell)
       IF(Ell==EllKeep) THEN
          NDist = NDist+1
          NCoef = NCoef+LenKet
       ENDIF
    ENDDO
!
!   Create New Rho
!
    CALL New_HGRho_new(Rho,(/NDist,NCoef,RhoTmp%NSDen/))
!
    NDist = 0
    NCoef = 0
    DO I=1,RhoTmp%NDist
       Ell    = RhoTmp%Ell%I(I)
       Zeta   = RhoTmp%Zeta%D(I)
       LenKet = LHGTF(Ell)
       OffCo1 = RhoTmp%OffCo%I(I)
       IF(Ell==EllKeep) THEN
          NDist = NDist+1
          NCoef = NCoef+LenKet
          OffCo2= NCoef-LenKet+1
          Rho%Ell%I(NDist)   = RhoTmp%Ell%I(I)
          Rho%Zeta%D(NDist)  = RhoTmp%Zeta%D(I)
          Rho%Qx%D(NDist)    = RhoTmp%Qx%D(I)
          Rho%Qy%D(NDist)    = RhoTmp%Qy%D(I)
          Rho%Qz%D(NDist)    = RhoTmp%Qz%D(I)
          Rho%OffCo%I(NDist) = OffCo2
          Rho%Co%D(OffCo2:OffCo2+LenKet-1) = RhoTmp%Co%D(OffCo1:OffCo1+LenKet-1)!<<< SPIN
       ENDIF
    ENDDO
    CALL Delete_HGRho_new(RhoTmp)
!
  END SUBROUTINE PruneEll_Rho_new
!---------------------------------------------------------------------------------------
! Fold the Distributions Back into the Box
!---------------------------------------------------------------------------------------
  SUBROUTINE Fold_Rho_new(GM,Rho)
    TYPE(HGRho_new)           :: Rho
    TYPE(CRDS)                :: GM
    INTEGER                   :: I
    REAL(DOUBLE),DIMENSION(3) :: Q
!
    DO I = 1,Rho%NDist
       Q(1) = Rho%Qx%D(I)
       Q(2) = Rho%Qy%D(I)
       Q(3) = Rho%Qz%D(I)
       CALL AtomCyclic(GM,Q)
       Rho%Qx%D(I) = Q(1)
       Rho%Qy%D(I) = Q(2)
       Rho%Qz%D(I) = Q(3)
    ENDDO
!
  END SUBROUTINE Fold_Rho_new
!---------------------------------------------------------------------------------------
! Calculate the total Dipole of Rho
!---------------------------------------------------------------------------------------
  SUBROUTINE CalRhoPoles_new(MP,Rho,GMLoc)
    TYPE(CMPoles)             :: MP
    TYPE(HGRho_new)           :: Rho
    TYPE(CRDS)                :: GMLoc
    INTEGER                   :: I,Ell,LenKet,OffCo,iadd,jadd
    REAL(DOUBLE)              :: RX,RY,RZ,Zeta
    REAL(DOUBLE),DIMENSION(3) :: Center
!
    Center(:) = GMLoc%PBC%CellCenter%D(:)
!
    MP%DPole%D = Zero
    MP%QPole%D = Zero
!
    DO I=1,Rho%NDist
       Zeta   = Rho%Zeta%D(I)
       Ell    = Rho%Ell%I(I)
       OffCo  = Rho%OffCo%I(I)
       LenKet = LHGTF(Ell)
       RX   = Rho%Qx%D(I)-Center(1)
       RY   = Rho%Qy%D(I)-Center(2)
       RZ   = Rho%Qz%D(I)-Center(3)
       MP%DPole%D=MP%DPole%D   +   CalculateDiPole(Ell,Zeta,RX,RY,RZ,Rho%Co%D(OffCo:OffCo+LenKet-1))!<<< SPIN
       MP%QPole%D=MP%QPole%D + CalculateQuadruPole(Ell,Zeta,RX,RY,RZ,Rho%Co%D(OffCo:OffCo+LenKet-1)) !<<< SPIN
    ENDDO
!
  END SUBROUTINE CalRhoPoles_new
!---------------------------------------------------------------------------------------
! Function MagnitudeDistribution
!---------------------------------------------------------------------------------------
  FUNCTION MagDist(LL,Expt,Coefs)
    INTEGER                      :: LL,L,M,N,LP,MP,NP,LMN,LMNP,LSUM
    REAL(DOUBLE)                 :: MagDist,Expt,Factor
    REAL(DOUBLE),DIMENSION(1:)    :: Coefs
!
    MagDist = Zero
    DO L=0,LL
       DO M=0,LL-L
          DO N=0,LL-L-M
             LMN=LMNDex(L,M,N)
             DO LP=0,LL
                DO MP=0,LL-LP
                   DO NP=0,LL-LP-MP
                      LMNP=LMNDex(LP,MP,NP)
                      LSUM = LP+MP+NP
                      Factor  =  ((-1)**LSUM)*DoubleFac(Expt,L+LP,M+MP,N+NP)
                      MagDist = MagDist + Factor*Coefs(LMN)*Coefs(LMNP)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    MagDist = SQRT(ABS(MagDist))*((Pi/(Expt))**ThreeHalves)
!
  END FUNCTION MagDist
!---------------------------------------------------------------------------------------
! Function DoubleFac
!---------------------------------------------------------------------------------------
  FUNCTION DoubleFac(Expt,L,M,N)
    INTEGER                       :: L,M,N
    REAL(DOUBLE)                  :: DoubleFac,Expt
    REAL(DOUBLE),DIMENSION(-1:11) :: DFac = (/1.D0,0.D0,1.D0,0.D0,3.D0,0.D0,15.D0,0.D0, &
                                              105.D0,0.D0,945.D0,0.D0,10395.D0/)
!
    DoubleFac = DFac(L-1)*DFac(M-1)*DFac(N-1)
    IF(DoubleFac == Zero) RETURN
    DoubleFac =  (Expt**((L+M+N)/2))*DoubleFac
!
  END FUNCTION DoubleFac
!#######################################################################################################################
!---------------------------------------------------------------------------------------
! Converts New Rho to Old Rho
!---------------------------------------------------------------------------------------
  SUBROUTINE ConvertToOldRho(Rho_old,Rho_new)
    TYPE(HGRho)           :: Rho_old
    TYPE(HGRho_new)       :: Rho_new
!
    INTEGER               :: NExpt,Ell,I,J,K,IE,OffCo,LenKet,iq,ir,NCoef,NDist,LKet1,LKet2, &
                             iSMat
    TYPE(DBL_VECT)        :: Expt
    TYPE(INT_VECT)        :: Lndx,NQ
    REAL(DOUBLE)          :: Zeta
!
!   Count the number of Exponents and MAX angular symetry
!
    IE=0
    NExpt = 0
    CALL New(Expt,Rho_new%NDist)
    CALL New(Lndx,Rho_new%NDist)
    DO I=1,Rho_new%NDist
       Zeta = Rho_new%Zeta%D(I)
       Ell  = Rho_new%Ell%I(I)
       DO J=1,NExpt
          IF(ABS(Expt%D(J)-Zeta) < 1.D-10) THEN
             IE = J
             GOTO 100
          ENDIF
       ENDDO
       NExpt         = NExpt+1
#ifdef PARALLEL
         IE = NExpt
#endif
       Expt%D(NExpt) = Zeta
       Lndx%I(NExpt) = Ell
100    CONTINUE
       IF(IE/=0)Lndx%I(IE) = MAX(Lndx%I(IE),Ell)
    ENDDO
!
!   Count the number of distributions per  Exponents, the number of distibutions and coefs
!
    NDist = 0
    NCoef = 0
    CALL New(NQ,NExpt)
    NQ%I      = 0
    DO I=1,Rho_new%NDist
       Zeta = Rho_new%Zeta%D(I)
       DO J=1,NExpt
          IF(ABS(Expt%D(J)-Zeta) < 1.D-10) THEN
             IE=J
             GOTO 200
          ENDIF
       ENDDO
       CALL Halt('MakeRho:ConvertToOldRho:Exponent Not Found')
200    CONTINUE
       NQ%I(IE) = NQ%I(IE)+1
       NDist = NDist + 1
       NCoef = NCoef + LHGTF(Lndx%I(IE))
    ENDDO
!
!   Allocate Rho_old
!
    CALL New(Rho_old,(/NExpt,NDist,NCoef,Rho_new%NSDen/))
!
!   Store the Indexing Info to Rho
!
    Rho_old%NQ%I(1:NExpt)   = NQ%I(1:NExpt)
    Rho_old%Expt%D(1:NExpt) = Expt%D(1:NExpt)
    Rho_old%Lndx%I(1:NExpt) = Lndx%I(1:NExpt)
!
!   Calculate OffQ and OffR
!
    Rho_old%OffQ%I = CalOffQ(Rho_old)
    Rho_old%OffR%I = CalOffR(Rho_old)
!
!   Store the Coef and Positions
!
    NQ%I            = 0
    DO I=1,Rho_new%NDist
!
       Zeta   = Rho_new%Zeta%D(I)
       Ell    = Rho_new%Ell%I(I)
       LKet1  = LHGTF(Ell)
!
       IE = 0
       DO J=1,NExpt
          IF(ABS(Rho_old%Expt%D(J)-Zeta) < 1.D-10) THEN
             IE=J
             GOTO 300
          ENDIF
       ENDDO
       CALL Halt('MakeRho:ConvertToOldRho:Exponent Not Found')
300    CONTINUE
!
       Ell    = Rho_old%Lndx%I(IE)
       LKet2  = LHGTF(Ell)
       NQ%I(IE) = NQ%I(IE)+1
       iq  = Rho_old%OffQ%I(IE) + NQ%I(IE)
!
       Rho_old%Qx%D(iq) = Rho_new%Qx%D(I)
       Rho_old%Qy%D(iq) = Rho_new%Qy%D(I)
       Rho_old%Qz%D(iq) = Rho_new%Qz%D(I)
!
       DO iSMat=1,Rho_old%NSDen !<<< SPIN
          OffCo= Rho_new%OffCo%I(I)                       +(iSMat-1)*Rho_new%NCoef !<<< SPIN
          ir   = Rho_old%OffR%I(IE)+(NQ%I(IE)-1)*LKet2+1  +(iSMat-1)*Rho_old%NCoef !<<< SPIN
          DO J=0,LKet2-1
             Rho_old%Co%D(ir+J) = Zero
          ENDDO
          DO J=0,MIN(LKet1,LKet2)-1
             Rho_old%Co%D(ir+J) = Rho_new%Co%D(OffCo+J)
          ENDDO
       ENDDO
!
    ENDDO
!
!   Clean Up
!
    CALL Delete(NQ)
    CALL Delete(Expt)
    CALL Delete(Lndx)
!
  END SUBROUTINE ConvertToOldRho
!--------------------------------------------------------------
! Calculate OffQ from NQ
!--------------------------------------------------------------
  FUNCTION CalOffQ(Rho)
    TYPE(HGRho)                  :: Rho
    INTEGER,DIMENSION(Rho%NExpt) :: CalOffQ
    INTEGER                      :: I,J,Iq
!
    CalOffQ = 0
    DO I = 1,Rho%NExpt
       Iq = 0
       DO J = 1,I-1
          Iq = Iq + Rho%NQ%I(J)
       ENDDO
       CalOffQ(I) = Iq
    ENDDO
!
  END FUNCTION CalOffQ
!--------------------------------------------------------------
! Calculate OffR from NQ and Lndx
!--------------------------------------------------------------
  FUNCTION CalOffR(Rho)
    TYPE(HGRho)                  :: Rho
    INTEGER,DIMENSION(Rho%NExpt) :: CalOffR
    INTEGER                      :: I,J,Ir
!
    CalOffR = 0
    DO I = 1,Rho%NExpt
       Ir = 0
       DO J = 1,I-1
          Ir = Ir + Rho%NQ%I(J)*LHGTF(Rho%Lndx%I(J))
       ENDDO
       CalOffR(I) = Ir
    ENDDO
!
  END FUNCTION CalOffR
!--------------------------------------------------------------
! Calculate the Number of Distributions from NQ
!--------------------------------------------------------------
  FUNCTION  CalNDist(Rho)
    TYPE(HGRho)                :: Rho
    INTEGER                    :: I,CalNDist
    CalNDist = 0
    DO I=1,Rho%NExpt
       CalNDist = CalNDist+Rho%NQ%I(I)
    ENDDO
  END FUNCTION CalNDist
!--------------------------------------------------------------
! Calculate the Number of Coefficients  from NQ and Lndx
!--------------------------------------------------------------
  FUNCTION CalNCoef(Rho)
    TYPE(HGRho)                :: Rho
    INTEGER                    :: I,CalNCoef
    CalNCoef = 0
    DO I=1,Rho%NExpt
       CalNCoef = CalNCoef+Rho%NQ%I(I)*LHGTF(Rho%Lndx%I(I))
    ENDDO
  END FUNCTION CalNCoef
!--------------------------------------------------------------------------------------------
!
! OLD STUFF
!
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
! Subroutine to Remove unnessesary distributions in rho
!--------------------------------------------------------------
  SUBROUTINE Prune_Rho(TOL,Rho_in,Rho_out)
    TYPE(HGRho)         :: Rho_in,Rho_out
    INTEGER             :: zq,iq,iqq,NQ,iadd,iiadd,jadd,jjadd
    INTEGER             :: LL,LenKet,NExpt,NDist,NCoef
    INTEGER             :: L,M,N,LMN,LP,MP,NP,LMNP
    REAL(DOUBLE)        :: Expt,Mag,TOL,Factor
    LOGICAL             :: HasQMI,HasMMI
!
    NExpt = Rho_in%NExpt
    CALL Delete_HGRho(Rho_out)
    CALL New_HGRho(Rho_out,(/NExpt,0,0,0/))
    Rho_out%Expt%D = Rho_in%Expt%D
    Rho_out%Lndx%I = Rho_in%Lndx%I
    Rho_out%NQ%I = 0
!
!   Count the Number of Distribution whose Magnitude is above TOL
!
    DO zq = 1,NExpt
       Expt   = Rho_in%Expt%D(zq)
       NQ     = Rho_in%NQ%I(zq)
       LL     = Rho_in%Lndx%I(zq)
       LenKet = LHGTF(LL)
       DO iq = 1,NQ
          iadd  = Rho_in%OffQ%I(zq) + iq
          jadd  = Rho_in%OffR%I(zq) + (iq-1)*LenKet+1
!         Calculate Magnitude
          Mag = MagDist(LL,Expt,Rho_in%Co%D(jadd:jadd+LenKet-1))
          IF(Mag>TOL) THEN
             Rho_out%NQ%I(zq) = Rho_out%NQ%I(zq) + 1
          ENDIF
       ENDDO
    ENDDO
!
!   Calculate NDist and NCoef from NQ and Lndx
!
    NDist = CalNDist(Rho_out)
    NCoef = CalNCoef(Rho_out)
!
!   Initailize  OffQ,OffR and RhoCo
!
    Rho_out%OffQ%I=CalOffQ(Rho_out)
    Rho_out%OffR%I=CalOffR(Rho_out)
!
!   Re-allocate the density
!
    CALL New_HGRho(Rho_out,(/NExpt,NDist,NCoef,Rho_out%NSDen/))
!
!   Add the Distributions to Rho_out
!
    DO zq = 1,NExpt
       Expt   = Rho_in%Expt%D(zq)
       NQ     = Rho_in%NQ%I(zq)
       LL     = Rho_in%Lndx%I(zq)
       LenKet = LHGTF(LL)
       iqq    = 1
       DO iq = 1,NQ
          iadd  = Rho_in%OffQ%I(zq) + iq
          jadd  = Rho_in%OffR%I(zq) + (iq-1)*LenKet+1
!         Calculate Magnitude
          Mag = MagDist(LL,Expt,Rho_in%Co%D(jadd:jadd+LenKet-1))
          IF(Mag>TOL)THEN
             iiadd = Rho_out%OffQ%I(zq) + iqq
             jjadd = Rho_out%OffR%I(zq) + (iqq-1)*LenKet+1
             Rho_out%Qx%D(iiadd) = Rho_in%Qx%D(iadd)
             Rho_out%Qy%D(iiadd) = Rho_in%Qy%D(iadd)
             Rho_out%Qz%D(iiadd) = Rho_in%Qz%D(iadd)
             DO LMN=0,LenKet-1
                Rho_out%Co%D(jjadd+LMN) = Rho_in%Co%D(jadd+LMN)
             ENDDO
             iqq = iqq + 1
          ENDIF
       ENDDO
    ENDDO
!
  END SUBROUTINE Prune_Rho
!---------------------------------------------------------------------
! Fold the Distributions Back into the Box
!---------------------------------------------------------------------
  SUBROUTINE Fold_Rho(GM,Rho)
    TYPE(HGRho)               :: Rho
    TYPE(CRDS)                :: GM
    INTEGER                   :: zq,iq,NQ,iadd
    INTEGER                   :: LL,LenKet
    REAL(DOUBLE)              :: Expt
    REAL(DOUBLE),DIMENSION(3) :: Q
!
    DO zq = 1,Rho%NExpt
       Expt   = Rho%Expt%D(zq)
       NQ     = Rho%NQ%I(zq)
       LL     = Rho%Lndx%I(zq)
       LenKet = LHGTF(LL)
       DO iq = 1,NQ
          iadd = Rho%OffQ%I(zq) + iq
          Q(1) = Rho%Qx%D(iadd)
          Q(2) = Rho%Qy%D(iadd)
          Q(3) = Rho%Qz%D(iadd)
          CALL AtomCyclic(GM,Q)
          Rho%Qx%D(iadd) = Q(1)
          Rho%Qy%D(iadd) = Q(2)
          Rho%Qz%D(iadd) = Q(3)
       ENDDO
    ENDDO
!
  END SUBROUTINE Fold_Rho
!---------------------------------------------------------------------
!  Integrate Rho
!---------------------------------------------------------------------
  SUBROUTINE Integrate_HGRho(Rho,RhoSumE,RhoSumN,RhoSumMM,MM_NATOMS,QM_NATOMS)
    TYPE(HGRho)                      :: Rho
    INTEGER                          :: zq,iq,oq,orr,NQ,jadd,LenKet,LL
    INTEGER                          :: MM_NATOMS,QM_NATOMS
    REAL(DOUBLE)                     :: RhoSumE,RhoSumN,RhoSumMM,Expt,Weig
    LOGICAL                          :: HasQMI,HasMMI
!
    HasMMI=.FALSE.
    HasQMI=.FALSE.
    IF(MM_NATOMS/=0) HasMMI=.TRUE.
    IF(QM_NATOMS/=0) HasQMI=.TRUE.
!
    RhoSumE = Zero
    RhoSumN = Zero
    RhoSumMM= Zero
!
    IF(HasQMI) THEN
!
       DO zq = 1,Rho%NExpt-1
          Expt   = Rho%Expt%D(zq)
          NQ     = Rho%NQ%I(zq)
          oq     = Rho%OffQ%I(zq)
          orr    = Rho%OffR%I(zq)
          LL     = Rho%Lndx%I(zq)
          LenKet = LHGTF(LL)
          DO iq = 1,NQ
             jadd = orr+(iq-1)*LenKet+1
!             WRITE(*,*) Expt,NQ,oq,orr,LL,LenKet,jadd
!             IF(.TRUE.) STOP
             Weig = Rho%Co%D(jadd)*((Pi/Expt)**1.5D0)
             RhoSumE = RhoSumE + Weig
          ENDDO
       ENDDO
!
!      Nuclear integration
!
       zq     = Rho%NExpt
       Expt   = Rho%Expt%D(zq)
       IF(HasMMI) THEN
!         MM charges may have been pruned, however this is very unlikely for QM nuclear charges
          NQ=QM_NATOMS
       ELSE
          NQ=Rho%NQ%I(zq)
       ENDIF
       oq     = Rho%OffQ%I(zq)
       orr    = Rho%OffR%I(zq)
       LL     = Rho%Lndx%I(zq)
       LenKet = LHGTF(LL)
       DO iq = 1,NQ
          jadd = orr+(iq-1)*LenKet+1
          Weig = Rho%Co%D(jadd)*((Pi/Expt)**ThreeHalves)
          RhoSumN = RhoSumN + Weig
       ENDDO
!
    ENDIF
!
    IF(HasMMI) THEN
!
!   MM integration
!
    zq     = Rho%NExpt
    Expt   = Rho%Expt%D(zq)
!     NQ=MM_NATOMS
      NQ=Rho%NQ%I(zq)-QM_NATOMS
      oq     = Rho%OffQ%I(zq)
      orr    = Rho%OffR%I(zq)
    IF(HasQMI) THEN
      oq     = oq+QM_NATOMS
      orr    = orr+QM_NATOMS
    ENDIF
    LL     = Rho%Lndx%I(zq)
    LenKet = LHGTF(LL)
    DO iq = 1,NQ
       jadd = orr+(iq-1)*LenKet+1
       Weig = Rho%Co%D(jadd)*((Pi/Expt)**ThreeHalves)
       RhoSumMM= RhoSumMM+ Weig
    ENDDO
!
    ENDIF
!
  END SUBROUTINE Integrate_HGRho

#ifdef PARALLEL
!---------------------------------------------------------------------
!  Parallel Integrate Rho
!---------------------------------------------------------------------
  SUBROUTINE ParaIntegrate_HGRho(Rho,RhoSumE,RhoSumN,RhoSumMM,MM_NATOMS,QM_NATOMS)
    TYPE(HGRho)                      :: Rho
    INTEGER                          :: zq,iq,oq,orr,NQ,jadd,LenKet,LL
    INTEGER                          :: MM_NATOMS,QM_NATOMS
    REAL(DOUBLE)                     :: RhoSumE,RhoSumN,RhoSumMM,Expt,Weig
    LOGICAL                          :: HasQMI,HasMMI
!
    HasMMI=.FALSE.
    HasQMI=.FALSE.
    IF(MM_NATOMS/=0) HasMMI=.TRUE.
    IF(QM_NATOMS/=0) HasQMI=.TRUE.
!
    RhoSumE = Zero
    RhoSumN = Zero
    RhoSumMM= Zero
!
    IF(HasQMI) THEN
!
    DO zq = 1,Rho%NExpt-1
       Expt   = Rho%Expt%D(zq)
       NQ     = Rho%NQ%I(zq)
       oq     = Rho%OffQ%I(zq)
       orr    = Rho%OffR%I(zq)
       LL     = Rho%Lndx%I(zq)
       LenKet = LHGTF(LL)
       DO iq = 1,NQ
          jadd = orr+(iq-1)*LenKet+1
          Weig = Rho%Co%D(jadd)*((Pi/Expt)**ThreeHalves)
          RhoSumE = RhoSumE + Weig
       ENDDO
    ENDDO
!
! Nuclear integration
!
    zq     = Rho%NExpt
    Expt   = Rho%Expt%D(zq)
    IF(HasMMI) THEN
      NQ=QM_NATOMS !!!! MM charges may have been pruned, however this is very unlikely for QM nuclear charges
!     NQ=NQ-MM_NATOMS
    ELSE
      NQ=Rho%NQ%I(zq)
    ENDIF
    oq     = Rho%OffQ%I(zq)
    orr    = Rho%OffR%I(zq)
    LL     = Rho%Lndx%I(zq)
    LenKet = LHGTF(LL)
#ifdef PARALLEL
    DO iq = Beg%I(MyID),End%I(MyID)
#else
    DO iq = 1,NQ
#endif
       jadd = orr+(iq-1)*LenKet+1
       Weig = Rho%Co%D(jadd)*((Pi/Expt)**ThreeHalves)
       RhoSumN = RhoSumN + Weig
    ENDDO
!
    ENDIF
!
    IF(HasMMI) THEN
!
! MM integration
!
    zq     = Rho%NExpt
    Expt   = Rho%Expt%D(zq)
!     NQ=MM_NATOMS
      NQ=Rho%NQ%I(zq)-QM_NATOMS
      oq     = Rho%OffQ%I(zq)
      orr    = Rho%OffR%I(zq)
    IF(HasQMI) THEN
      oq     = oq+QM_NATOMS
      orr    = orr+QM_NATOMS
    ENDIF
    LL     = Rho%Lndx%I(zq)
    LenKet = LHGTF(LL)
    DO iq = 1,NQ
       jadd = orr+(iq-1)*LenKet+1
       Weig = Rho%Co%D(jadd)*((Pi/Expt)**ThreeHalves)
       RhoSumMM= RhoSumMM+ Weig
    ENDDO
!
    ENDIF
!
  END SUBROUTINE ParaIntegrate_HGRho
#endif
!========================================================================================
! Calculate the total Dipole of Rho
!========================================================================================
  SUBROUTINE CalRhoPoles(MP,Center,Rho,GMLoc)
    TYPE(CMPoles)             :: MP
    TYPE(HGRho)               :: Rho
    INTEGER                   :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LQ,LenQ
    REAL(DOUBLE)              :: RX,RY,RZ,R2,Expt
    REAL(DOUBLE),DIMENSION(3) :: Center
    TYPE(CRDS)                :: GMLoc
!
    MP%DPole%D = Zero
    MP%QPole%D = Zero
    DO zq=1,Rho%NExpt
       NQ     = Rho%NQ%I(zq)
       Expt   = Rho%Expt%D(zq)
       OffQ   = Rho%OffQ%I(zq)
       OffR   = Rho%OffR%I(zq)
       LQ     = Rho%Lndx%I(zq)
       LenQ   = LHGTF(LQ)
       IF(NQ > 0) THEN
          DO iq = 1,NQ
             iadd = Rho%OffQ%I(zq)+iq
             jadd = Rho%OffR%I(zq)+(iq-1)*LenQ+1
             RX   = Rho%Qx%D(iadd)-Center(1)
             RY   = Rho%Qy%D(iadd)-Center(2)
             RZ   = Rho%Qz%D(iadd)-Center(3)
             MP%DPole%D=MP%DPole%D+CalculateDiPole(LQ,Expt,RX,RY,RZ,Rho%Co%D(jadd:jadd+LenQ-1))
             MP%QPole%D=MP%QPole%D+CalculateQuadruPole(LQ,Expt,RX,RY,RZ,Rho%Co%D(jadd:jadd+LenQ-1))
          ENDDO
       ENDIF
    ENDDO
!
  END SUBROUTINE CalRhoPoles
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
  SUBROUTINE RhoFill(Rho2,Rho,IExpt,IDist,ICoef)
!
    INTEGER :: IExpt,IDist,ICoef
    TYPE(HGRho) :: Rho,Rho2
!
!   Fill the content of Rho2 into Rho upto the given numbers
!
    Rho%NQ%I(1:IExpt)=Rho2%NQ%I(1:IExpt)
    Rho%Lndx%I(1:IExpt)=Rho2%Lndx%I(1:IExpt)
    Rho%OffQ%I(1:IExpt)=Rho2%OffQ%I(1:IExpt)
    Rho%OffR%I(1:IExpt)=Rho2%OffR%I(1:IExpt)
    Rho%Expt%D(1:IExpt)=Rho2%Expt%D(1:IExpt)
    Rho%Qx%D(1:IDist)=Rho2%Qx%D(1:IDist)
    Rho%Qy%D(1:IDist)=Rho2%Qy%D(1:IDist)
    Rho%Qz%D(1:IDist)=Rho2%Qz%D(1:IDist)
    Rho%Co%D(1:ICoef)=Rho2%Co%D(1:ICoef)
!
  END SUBROUTINE RhoFill
!
END MODULE RhoTools
!
