MODULE RhoTools
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  USE RhoBlok
  USE AtomPairs
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
CONTAINS
!--------------------------------------------------------------
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
    CALL New_HGRho(Rho_out,(/NExpt,0,0/))
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
    CALL New_HGRho(Rho_out,(/NExpt,NDist,NCoef/))
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
#ifdef PERIODIC
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
#endif
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
!#ifdef PERIODIC
!    Center(:) = GMLoc%PBC%CellCenter(:)
!#else
!    Center(:) = Half*(GMLoc%BndBox%D(:,2)+GMLoc%BndBox%D(:,1))
!#endif
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
    MagDist = SQRT(ABS(MagDist)*((Pi/(Two*Expt))**ThreeHalves))
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
!
!-----------------------------------------------------------
!
    SUBROUTINE RhoFill(Rho2,Rho,IExpt,IDist,ICoef)
!
    INTEGER :: IExpt,IDist,ICoef
    TYPE(HGRho) :: Rho,Rho2
!
! Fill the content of Rho2 into Rho upto the given numbers
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



END MODULE RhoTools
!
