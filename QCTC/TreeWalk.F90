!    PERFORM O(Lg N) WALKS ON THE POLETREE REPRESENTATION OF THE TOTAL 
!    ELECTRON DENSITY, USING THE K-D TREE DATA STRUCTURE TO APPLY THE 
!    PENETRATION ACCESABILITY CRITERIA
!
!    Author: Matt Challacombe 
!------------------------------------------------------------------------------
MODULE TreeWalk
  USE DerivedTypes
  USE PoleNodeType
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Globals
  USE PoleTree
  USE MondoPoles 
  USE GammaFunctions
  USE PoleGlobals
  USE ERIGlobals
  IMPLICIT NONE
!---------------------------------------------------------------------
  CONTAINS  
!===================================================================================================
!
!===================================================================================================
     SUBROUTINE SetKet(P,E)
       TYPE(PrimPair) :: P
       REAL(DOUBLE)   :: E
!---------------------------------------------
!      Zero <KET| accumulators
       ! HGKet=Zero
       ! SPKetC=Zero
       ! SPKetS=Zero
       CALL DBL_VECT_EQ_DBL_SCLR(HGLen,HGKet(1),Zero)
       CALL DBL_VECT_EQ_DBL_SCLR(SPLen+1,SPKetC(0),Zero)
       CALL DBL_VECT_EQ_DBL_SCLR(SPLen+1,SPKetS(0),Zero)
!      Set global BBox for this primitive
       PBox%BndBox(1,1)=P%P(1)
       PBox%BndBox(2,1)=P%P(2)
       PBox%BndBox(3,1)=P%P(3)
       PBox%BndBox(1,2)=P%P(1)
       PBox%BndBox(2,2)=P%P(2)
       PBox%BndBox(3,2)=P%P(3)
!      Expand this BBox from zero to the correct extent
#ifdef NewPAC
       E=1.D-16
#endif
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
#ifdef NewPAC
       REAL(DOUBLE)                     :: SqrtW,SqrtT,LeftS,RightS
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
#ifdef NewPAC
       SqrtW  = SQRT(PrimBeta*Q%Beta/(PrimBeta+Q%Beta))
       SqrtT  = SqrtW*(SQRT(PQ2)-SQRT(Q%DMax2))
       LeftS  = PrimWCoef*SqrtW*Erfcc(SqrtT)
       RightS = Q%PACStr*SqrtT
       IF(SqrtT > Zero .AND. LeftS < RightS .OR. T>Gamma_Switch) THEN
#else
       IF(ABS(PQx)>PBox%Half(1)+Q%Box%Half(1).OR.  &
          ABS(PQy)>PBox%Half(2)+Q%Box%Half(2).OR.  &
          ABS(PQz)>PBox%Half(3)+Q%Box%Half(3).OR.  &
          T>Gamma_Switch)THEN
#endif
!         MAC: 
          IF(PQ2>(Q%Strength*DP2+Q%DMax2).OR.Q%Leaf)THEN 
!            Evaluate multipoles
#ifdef EXPLICIT_SOURCE
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
!
             SELECT CASE(Ell)
                INCLUDE "IrRegulars.Inc"
             CASE DEFAULT
                CALL Halt('No explicit code for case Ell   = '//TRIM(IntToChar(Ell))//' in IrRegulars.Inc')
             END SELECT
!
             SELECT CASE(LCode)
                INCLUDE "CTraX.Inc"
             CASE DEFAULT
                CALL Halt('No explicit code for case LCode = '//TRIM(IntToChar(LCode))//' in CTraX.Inc')
             END SELECT
!
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
          Ell=Prim%Ell+Q%Ell
          LCode=Prim%Ell*100+Q%Ell
          J=AINT(T*Gamma_Grid)
!
          SELECT CASE(Ell)
             INCLUDE 'AuxInt.Inc'
          CASE DEFAULT
             CALL Halt('No explicit code for case Ell   = '//TRIM(IntToChar(Ell))//' in AuxInt.Inc')
          END SELECT
!
          SELECT CASE(LCode)
             INCLUDE 'HGTraX.Inc'
          CASE DEFAULT
             CALL Halt('No explicit code for case LCode = '//TRIM(IntToChar(LCode))//' in HGTraX.Inc')
          END SELECT
!
#else
          Ell=Prim%Ell+Q%Ell
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
       REAL(DOUBLE),PARAMETER           :: VTol = 1.D-12
#ifdef NewPAC
       REAL(DOUBLE)                     :: SqrtW,SqrtT,LeftS,RightS
#endif
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
#ifdef NewPAC
       SqrtW  = SQRT(GFactor*Omega)
       SqrtT  = SqrtW*(SQRT(PQ2)-SQRT(Q%DMax2))
       LeftS  = PrimWCoef*SqrtW*Erfcc(SqrtT)
       RightS = Q%PACStr*SqrtT
       IF(LeftS < RightS) THEN
#else
       IF(ABS(PQx)>PBox%Half(1)+Q%Box%Half(1).OR.  &
          ABS(PQy)>PBox%Half(2)+Q%Box%Half(2).OR.  &
          ABS(PQz)>PBox%Half(3)+Q%Box%Half(3).OR.  &
          T>Gamma_Switch)THEN
#endif
!         MAC:
          IF((PQ2>Q%Strength*DP2+Q%DMax2).OR.Q%Leaf)THEN
             IF(Q%Zeta==NuclearExpnt .AND. PQ2<VTol) RETURN
!            Evaluate multipoles
#ifdef EXPLICIT_SOURCE
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
!
             SELECT CASE(Ell)
                INCLUDE "IrRegulars.Inc"
             CASE DEFAULT
                CALL Halt('No explicit code for case Ell   = '//TRIM(IntToChar(Ell))//' in IrRegulars.Inc')
             END SELECT
!
             SELECT CASE(LCode)
                INCLUDE "CTraX.Inc"
             CASE DEFAULT
                CALL Halt('No explicit code for case LCode = '//TRIM(IntToChar(LCode))//' in CTraX.Inc')
             END SELECT
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
#ifdef EXPLICIT_SOURCE
          Ell=Prim%Ell+Q%Ell
          LCode=Prim%Ell*100+Q%Ell
          J=AINT(T*Gamma_Grid)
!
          SELECT CASE(Ell)
             INCLUDE 'AuxInt.Inc'
          CASE DEFAULT
             CALL Halt('No explicit code for case Ell   = '//TRIM(IntToChar(Ell))//' in AuxInt.Inc')
          END SELECT
!
          SELECT CASE(LCode)
             INCLUDE 'HGTraX.Inc'
          CASE DEFAULT
             CALL Halt('No explicit code for case LCode = '//TRIM(IntToChar(LCode))//' in HGTraX.Inc')
          END SELECT
!
#else
          Ell=Prim%Ell+Q%Ell
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
                                                                   
