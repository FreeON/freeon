MODULE RhoTools
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  USE RhoBlok
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
          IF(Mag  .GT. TOL) THEN
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
          IF(Mag .GT. TOL) THEN
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
  SUBROUTINE Integrate_HGRho(Rho)
    TYPE(HGRho)                      :: Rho
    INTEGER                          :: zq,iq,oq,orr,NQ,jadd,LenKet,LL
    REAL(DOUBLE)                     :: RhoSumE,RhoSumN,Expt,Weig
!
    RhoSumE = Zero
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
    RhoSumN = Zero
    zq     = Rho%NExpt
    Expt   = Rho%Expt%D(zq)
    NQ     = Rho%NQ%I(zq)
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
    IF(ABS(RhoSumE+RhoSumN)>1.D-2) &
       CALL Halt(' Density hosed! Rho_e = '//TRIM(DblToMedmChar(Two*RhoSumE)) &
                              //',Rho_n = '//TRIM(DblToMedmChar(Two*RhoSumN)))

!    WRITE(*,*) ' Int[Rho_E] = ',RhoSumE
!    WRITE(*,*) ' Int[Rho_N] = ',RhoSumN   
!    WRITE(*,*) ' Int[Rho_T] = ',RhoSumN+RhoSumE
  END SUBROUTINE Integrate_HGRho

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
END MODULE RhoTools
!
!
!!$    INTEGER             :: QAdd,RAdd,SAdd,TAdd,Or,Oq,zq,Negl,K,LMN,Ell,LenKet
!!$    REAL(DOUBLE)        :: NEps,Chg
!!$    TYPE(INT_VECT)      :: IOrdr
!!$    TYPE(DBL_VECT)      :: ROrdr
!!$    Negl=1
!!$    SAdd=1
!!$    TAdd=0
!!$    Rho_out%NExpt=Rho_in%NExpt
!!$    Rho_out%Expt%D=Rho_in%Expt%D
!!$    DO Zq=1,Rho_in%NExpt
!!$       Chg=(Pi/Rho_in%Expt%D(Zq))**ThreeHalves
!!$       Oq=Rho_in%OffQ%I(Zq)
!!$       Or=Rho_in%OffR%I(Zq)
!!$       Ell=Rho_in%Lndx%I(Zq) 
!!$       Rho_out%Lndx%I(Zq)=Ell 
!!$       LenKet=LHGTF(Ell)
!!$       CALL New(IOrdr,Rho_in%Nq%I(Zq))
!!$       CALL New(ROrdr,Rho_in%Nq%I(Zq))
!!$       DO Iq=1,Rho_in%Nq%I(Zq)
!!$          QAdd=Oq+Iq
!!$          IOrdr%I(IQ)=IQ
!!$          ROrdr%D(IQ)=Rho_in%Co%D(QAdd)*Chg
!!$       ENDDO
!!$       IF(Zq<Rho_in%NExpt.AND.Rho_in%Nq%I(ZQ)>1) THEN
!!$          CALL DblIntSort77(Rho_in%NQ%I(Zq),ROrdr%D,IOrdr%I,+2)
!!$       ENDIF
!!$       Rho_out%NQ%I(Zq)=0
!!$       DO K=1,Rho_in%Nq%I(Zq)
!!$          IF(ROrdr%D(K)>NEps/DBLE(Negl).OR.Zq==Rho_in%NExpt)THEN
!!$             Negl=Negl+1
!!$             Rho_out%NQ%I(Zq)=Rho_out%NQ%I(Zq)+1
!!$             Iq=IOrdr%I(K)
!!$             QAdd=Oq+Iq
!!$             RAdd=Or+(Iq-1)*LenKet
!!$             Rho_out%Qx%D(SAdd)=Rho_in%Qx%D(QAdd)
!!$             Rho_out%Qy%D(SAdd)=Rho_in%Qy%D(QAdd)
!!$             Rho_out%Qz%D(SAdd)=Rho_in%Qz%D(QAdd)
!!$             DO LMN=1,LenKet
!!$                Rho_out%Co%D(TAdd+LMN)=Rho_in%Co%D(RAdd+LMN)
!!$             ENDDO
!!$             SAdd=SAdd+1
!!$             TAdd=TAdd+LenKet
!!$          ENDIF
!!$       ENDDO
!!$       CALL Delete(IOrdr)
!!$       CALL Delete(ROrdr)
!!$    ENDDO
