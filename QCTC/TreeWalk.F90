!    PERFORM O(Lg N) WALKS ON THE POLETREE REPRESENTATION OF THE TOTAL 
!    ELECTRON DENSITY, USING THE K-D TREE DATA STRUCTURE TO APPLY THE 
!    PENETRATION ACCESABILITY CRITERIA
!
!    Author: Matt Challacombe 
!------------------------------------------------------------------------------
MODULE TreeWalk
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Globals
  USE PoleTree
  USE MondoPoles 
  USE GammaFunctions
  IMPLICIT NONE
!---------------------------------------------------------------------
! GLOBAL ..
  TYPE(PrimPair)                  :: Prim
  TYPE(BBox)                      :: PBox
  REAL(DOUBLE)                    :: DP2
  REAL(DOUBLE)                    :: PoleSwitch
  INTEGER                         :: At
  REAL(DOUBLE),DIMENSION(1:HGLen) :: HGKet
  REAL(DOUBLE),DIMENSION(0:SPLen) :: SPKetC
  REAL(DOUBLE),DIMENSION(0:SPLen) :: SPKetS
! WORKSPACE FOR ERI INTERMEDIATES
  REAL(DOUBLE),DIMENSION(500)      :: R
  REAL(DOUBLE),DIMENSION(20)       :: W
  REAL(DOUBLE),DIMENSION(0:20)     :: AuxR,G
!-----------
  CONTAINS  !
     SUBROUTINE SetKet(P,E)
       TYPE(PrimPair) :: P
       REAL(DOUBLE)   :: E
!---------------------------------------------
!      Zero <KET| accumulators
       HGKet=Zero
       SPKetC=Zero
       SPKetS=Zero
!      Set global BBox for this primitive
       PBox%BndBox(:,1)=P%P
       PBox%BndBox(:,2)=P%P
!      Expand this BBox from zero to the correct extent
       PBox=ExpandBox(PBox,E)
    END SUBROUTINE SetKet
!===================================================================================================
!
!===================================================================================================
     RECURSIVE SUBROUTINE JWalk(Q)
       TYPE(PoleNode)                   :: Q
       REAL(DOUBLE)                     :: PQ2,PQ,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
                                           CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
                                           TT,TwoC,COne,SOne,CTwo,STwo,RTE,RPE,T,Upq
       INTEGER                          :: J,Ell,LCode
#ifdef EXPLICIT_SOURCE
       REAL(DOUBLE)                     :: o1,o2,ET,TwoT
#else
       INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex
       REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
#endif
!-----------------------------------------------------------------------------------------------
!      PAC: 
       PQx=Prim%P(1)-Q%Box%Center(1)
       PQy=Prim%P(2)-Q%Box%Center(2)
       PQz=Prim%P(3)-Q%Box%Center(3)
       PQ2=PQx*PQx+PQy*PQy+PQz*PQz
       RTE=Prim%Zeta*Q%Zeta
       RPE=Prim%Zeta+Q%Zeta
       Omega=RTE/RPE 
       T=Omega*PQ2
       IF(ABS(PQx)>PBox%Half(1)+Q%Box%Half(1).OR.  &
          ABS(PQy)>PBox%Half(2)+Q%Box%Half(2).OR.  &
          ABS(PQz)>PBox%Half(3)+Q%Box%Half(3).OR.  &
          T>Gamma_Switch)THEN
!         MAC: 
          IF(PQ2>(Q%Strength*DP2+Q%DMax2).OR.Q%Leaf)THEN 
!            Evaluate multipoles
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE
             INCLUDE "IrRegulars.Inc"
             INCLUDE "CTraX.Inc"
#else
             CALL CTrax(Prim,Q,SPKetC,SPKetS)
#endif
          ELSE
!            Keep on truckin
             CALL JWalk(Q%Descend)
             CALL JWalk(Q%Descend%Travrse)
          ENDIF
       ELSEIF(Q%Leaf)THEN
          PQx=-PQx; PQy=-PQy; PQz=-PQz
          Upq=TwoPi5x2/(RTE*SQRT(RPE))
!         Compute a Hermite Gaussian Electron Repulsion Integral (HGERI)
          Ell=Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE
          LCode=Prim%Ell*100+Q%Ell
          INCLUDE 'HGTraX.Inc'
#else
          CALL AuxInts(2*HGEll,Ell,AuxR,Omega,T)
          CALL MD3TRR(2*HGEll,Ell,MDR,AuxR,Upq,PQx,PQy,PQz)
          DO LP=0,Prim%Ell
             DO MP=0,Prim%Ell-LP
                DO NP=0,Prim%Ell-LP-MP 
                   PDex=LMNDex(LP,MP,NP)
                   DO LQ=0,Q%Ell
                      DO MQ=0,Q%Ell-LQ
                         DO NQ=0,Q%Ell-LQ-MQ
                            QDex=LMNDex(LQ,MQ,NQ)
                            HGKet(PDex)=HGKet(PDex)+MDR(LP+LQ,MP+MQ,NP+NQ,0)*Q%Co(QDex)
                         ENDDO 
                      ENDDO 
                   ENDDO
                ENDDO 
             ENDDO 
          ENDDO
#endif
       ELSE
!         Keep on truckin 
          CALL JWalk(Q%Descend)
          CALL JWalk(Q%Descend%Travrse)
       ENDIF
     END SUBROUTINE JWalk
!===================================================================================================
!
!===================================================================================================
     RECURSIVE SUBROUTINE VWalk(Q)
       TYPE(PoleNode)                   :: Q
       REAL(DOUBLE)                     :: PQ2,PQ,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
                                           CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
                                           TT,TwoC,COne,SOne,CTwo,STwo,RTE,RPE,T,Upq
       INTEGER                          :: J,Ell,LCode
#ifdef EXPLICIT_SOURCE
       REAL(DOUBLE)                     :: o1,o2,ET,TwoT
#else
       INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex
       REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
#endif
       REAL(DOUBLE),PARAMETER           :: VTol = 1.D-8
!---------------------------------------------------------------------------------------------------
!      PAC:
       PQx=Prim%P(1)-Q%Box%Center(1)
       PQy=Prim%P(2)-Q%Box%Center(2)
       PQz=Prim%P(3)-Q%Box%Center(3)
       PQ2=PQx*PQx+PQy*PQy+PQz*PQz
       RTE=Prim%Zeta*Q%Zeta
       RPE=Prim%Zeta+Q%Zeta
       Omega=RTE/RPE 
       T=Omega*PQ2
       IF(ABS(PQx)>PBox%Half(1)+Q%Box%Half(1).OR.  &
          ABS(PQy)>PBox%Half(2)+Q%Box%Half(2).OR.  &
          ABS(PQz)>PBox%Half(3)+Q%Box%Half(3).OR.  &
              T>Gamma_Switch) THEN
!         MAC:
          IF((PQ2>Q%Strength*DP2+Q%DMax2).OR.Q%Leaf)THEN
!            Evaluate multipoles
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE
             INCLUDE "IrRegulars.Inc"
             INCLUDE "CTraX.Inc"
#else
            CALL CTrax(Prim,Q,SPKetC,SPKetS)
#endif
          ELSE
!            Keep on truckin
             CALL VWalk(Q%Descend)
             CALL VWalk(Q%Descend%Travrse)
          ENDIF
       ELSEIF(Q%Leaf)THEN
!         Check for self-interaction
          IF(Q%Zeta==NuclearExpnt .AND. PQ2<VTol) RETURN
!         Flip sign...
          PQx=-PQx; PQy=-PQy; PQz=-PQz
          Upq=TwoPi5x2/(RTE*SQRT(RPE))
!         Compute a Hermite Gaussian Electron Repulsion Integral (HGERI)
          Ell=Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE
          LCode=Prim%Ell*100+Q%Ell
          INCLUDE 'HGTraX.Inc'
#else
          CALL AuxInts(2*HGEll,Ell,AuxR,Omega,T)
          CALL MD3TRR(2*HGEll,Ell,MDR,AuxR,Upq,PQx,PQy,PQz)
          DO LP=0,Prim%Ell
             DO MP=0,Prim%Ell-LP
                DO NP=0,Prim%Ell-LP-MP 
                   PDex=LMNDex(LP,MP,NP)
                   DO LQ=0,Q%Ell
                      DO MQ=0,Q%Ell-LQ
                         DO NQ=0,Q%Ell-LQ-MQ
                            QDex=LMNDex(LQ,MQ,NQ)
                            HGKet(PDex)=HGKet(PDex)+MDR(LP+LQ,MP+MQ,NP+NQ,0)*Q%Co(QDex)
                         ENDDO 
                      ENDDO 
                   ENDDO
                ENDDO 
             ENDDO 
          ENDDO
#endif
       ELSE
!         Keep on trunkin
          CALL VWalk(Q%Descend)
          CALL VWalk(Q%Descend%Travrse)
       ENDIF
     END SUBROUTINE VWalk
END MODULE
                                                                   
