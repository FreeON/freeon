!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and Michael R. Salazar
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!--  COMPUTE THE NUCLEAR-TOTAL ELECTROSTATIC ENERGY IN O(N Lg N) CPU TIME
!==============================================================================
MODULE GradNukE
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE Thresholding
  USE Globals
  USE AtomPairs
  USE BraKetBloks
  USE PoleTree 
  USE MondoPoles
  IMPLICIT NONE
!---------------------------------------
  TYPE(PrimPair)  :: Prim
  REAL(DOUBLE)    :: DP2
  REAL(DOUBLE)    :: PoleSwitch
  INTEGER         :: At
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    SUBROUTINE GradNukE(dE_Nuc) 
       REAL(DOUBLE) :: NukE,Tau,NukeCo,NukePole
       TYPE(DBL_RNK2)      :: dE_Nuc
!---------------------------------------------------------------------------------------------
!      Set MAC and PAC parameters
       DO At=1,NAtoms
          NukeCo=-GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          NukePole=-GM%AtNum%I(At)
          Tau=Thresholds%TwoE
          DP2=(GM%AtNum%I(At)/Tau)**(Two/DBLE(SPEll+1))
          PoleSwitch=MAX(1.D-12,-LOG(Two*Tau/GM%AtNum%I(At)))
!         Set atomic "primitive"
          Prim%P=GM%Carts%D(:,At) 
          Prim%Zeta=NuclearExpnt
          Prim%HGKet(1)=Zero
          Prim%HGKet(2)=Zero
          Prim%HGKet(3)=Zero
          Prim%HGKet(4)=Zero
          Prim%SPKetC(0)=Zero
          Prim%SPKetC(1)=Zero
          Prim%SPKetC(2)=Zero
          Prim%SPKetS(2)=Zero

          Prim%Ell=1
!         Walk the walk
          CALL VWalk(PoleRoot)

!         Accumulate the atomic contribution
!          NukE=NukE+NukeCo*Prim%HGKet(1)+NukePole*Prim%SPKetC(0)

!         See HGToSP.Inc to understand the following:
          dE_Nuc%D(1,At)=-NukeCo*Prim%HGKet(2)-Half*NukePole*Prim%SPKetS(2)
          dE_Nuc%D(2,At)=-NukeCo*Prim%HGKet(3)-Half*NukePole*Prim%SPKetC(2)
          dE_Nuc%D(3,At)=-NukeCo*Prim%HGKet(4)+NukePole*Prim%SPKetC(1)

       ENDDO
     END SUBROUTINE GradNukE
!=============================================================================================
!
!=============================================================================================
     RECURSIVE SUBROUTINE VWalk(Q)
       TYPE(PoleNode)                   :: Q
       REAL(DOUBLE)                     :: PQ2,PQ,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
                                           CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
                                           TT,TwoC,COne,SOne,CTwo,STwo  
       REAL(DOUBLE),DIMENSION(0:2*SPEll):: CoFact
       REAL(DOUBLE),DIMENSION(0:SPLen2) :: ALegendreP
       REAL(DOUBLE),DIMENSION(50)       :: W
       INTEGER                          :: Ell,LCode
!---------------------------------------------------------------------------------------------------
!      PAC: Exp[-W_Min PQ^2]/2 < Tau/Amp
       PQx=Prim%P(1)-Q%Box%Center(1)
       PQy=Prim%P(2)-Q%Box%Center(2)
       PQz=Prim%P(3)-Q%Box%Center(3)
       PQ2=PQx*PQx+PQy*PQy+PQz*PQz
       Omega=Prim%Zeta*Q%Zeta/(Prim%Zeta+Q%Zeta)
       IF(Omega*PQ2>PoleSwitch)THEN
!         MAC: (d_q/|PQ|)^(p+1) < Tau/Amp
          IF(PQ2>Q%D2*DP2.OR.Q%Leaf)THEN
!          IF(Q%Leaf)THEN
             IF(PQ2==Zero)RETURN
!            Evaluate multipoles
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell
#ifdef EXPLICIT_SOURCE
             INCLUDE "IrRegulars.Inc"
             INCLUDE "CTraX.Inc"
#else
             CALL CTrax(Prim,Q)
#endif
          ELSE
!            Keep walking
             CALL VWalk(Q%Descend)
             CALL VWalk(Q%Descend%Travrse)
          ENDIF
       ELSEIF(Q%Leaf)THEN
          IF(Q%BDex>0.AND.PQ2==Zero)RETURN
          CALL MDERI(Prim,Q)
       ELSE
!         Keep walking
          CALL VWalk(Q%Descend)
          CALL VWalk(Q%Descend%Travrse)
       ENDIF
     END SUBROUTINE VWalk
END MODULE
