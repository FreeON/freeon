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
!    PERFORM O(Lg N) WALKS ON THE POLETREE ACCUMULATING COULOMB SUMS
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
!===================================================================================================
!
!===================================================================================================
     RECURSIVE SUBROUTINE JWalk(Q)
       TYPE(PoleNode)                   :: Q
       REAL(DOUBLE)                     :: PQ2,PQ,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
                                           CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
                                           TT,TwoC,COne,SOne,CTwo,STwo,RTE,RPE,T,Upq
       INTEGER                          :: J,Ell,LCode
#ifdef EXPLICIT_SOURCE_DIRECT
       REAL(DOUBLE)                     :: o1,o2,ET,TwoT
#else
       INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex
       REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
#endif
!-----------------------------------------------------------------------------------------------
!      PAC: Exp[-W_pq PQ^2]/2 < Tau/Amp
!
       PQx=Prim%P(1)-Q%Box%Center(1)
       PQy=Prim%P(2)-Q%Box%Center(2)
       PQz=Prim%P(3)-Q%Box%Center(3)
       PQ2=PQx*PQx+PQy*PQy+PQz*PQz
       RTE=Prim%Zeta*Q%Zeta
       RPE=Prim%Zeta+Q%Zeta
       Omega=RTE/RPE 
       IF(Omega*PQ2>PoleSwitch)THEN
!         MAC: (d_q/|PQ|)^(p+1) < Tau/Amp
          IF(PQ2>Q%D2*DP2 .OR. Q%Leaf)THEN 
!            Evaluate multipoles
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE_POLES
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
!         Set up
          T=Omega*PQ2          
          Upq=TwoPi5x2/(RTE*SQRT(RPE))
!         Sign inversion
          PQx=-PQx; PQy=-PQy; PQz=-PQz
!         Compute a Hermite Gaussian Electron Repulsion Integral (HGERI)
          Ell=Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE_DIRECT
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
#ifdef EXPLICIT_SOURCE_DIRECT
       REAL(DOUBLE)                     :: o1,o2,ET,TwoT
#else
       INTEGER                          :: LP,MP,NP,LQ,MQ,NQ,PDex,QDex
       REAL(DOUBLE),DIMENSION(0:2*HGEll,0:2*HGEll,0:2*HGEll,0:2*HGEll) :: MDR
#endif
       REAL(DOUBLE),PARAMETER           :: VTol = 1.D-14
!---------------------------------------------------------------------------------------------------
!      PAC: Exp[-W_Min PQ^2]/2 < Tau/Amp
       PQx=Prim%P(1)-Q%Box%Center(1)
       PQy=Prim%P(2)-Q%Box%Center(2)
       PQz=Prim%P(3)-Q%Box%Center(3)
       PQ2=PQx*PQx+PQy*PQy+PQz*PQz
       RTE=Prim%Zeta*Q%Zeta
       RPE=Prim%Zeta+Q%Zeta
       Omega=RTE/RPE
       IF(Omega*PQ2>PoleSwitch)THEN
!         MAC: (d_q/|PQ|)^(p+1) < Tau/Amp
          IF(PQ2>Q%D2*DP2.OR.Q%Leaf)THEN
             IF(PQ2==Zero)RETURN
!            Evaluate multipoles
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE_POLES
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
          IF(Q%BDex>0 .AND. PQ2 < VTol)RETURN
!         Set up
          T=Omega*PQ2          
          Upq=TwoPi5x2/(RTE*SQRT(RPE))
!         Sign inversion
          PQx=-PQx; PQy=-PQy; PQz=-PQz
!         Compute a Hermite Gaussian Electron Repulsion Integral (HGERI)
          Ell=Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE_DIRECT
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
                                                                   
