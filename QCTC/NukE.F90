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
!    Authors: Matt Challacombe and C.J. Tymczak
!    COMPUTE THE NUCLEAR-TOTAL ELECTROSTATIC ENERGY IN O(N Lg N) CPU TIME
!==============================================================================
MODULE NuklarE
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE QCTCThresholds
  USE Thresholding
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree 
  USE TreeWalk
#ifdef PERIODIC
  USE PBCFarField
#endif
  IMPLICIT NONE
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    FUNCTION NukE()
       REAL(DOUBLE)                   :: NukE
       REAL(DOUBLE),DIMENSION(1)      :: NukeCo,NukePole
       REAL(DOUBLE),DIMENSION(1:1)    :: HGBra
       REAL(DOUBLE),DIMENSION(0:0)    :: SPBraC,SPBraS
#ifdef PERIODIC
       INTEGER                        :: NC
       REAL(DOUBLE),DIMENSION(3)      :: QC
#endif
!---------------------------------------------------------------------------------------------
       NukE=Zero 
       DO At=1,NAtoms
          NukeCo=-GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          NukePole=-GM%AtNum%I(At)
          DP2=(GM%AtNum%I(At)/TauMAC)**(Two/DBLE(SPEll+2))
!         Set atomic "primitive"  
          Prim%P=GM%Carts%D(:,At)
          Prim%Zeta=NuclearExpnt
          PBox%BndBox(:,1)=Prim%P
          PBox%BndBox(:,2)=Prim%P
          PBox=ExpandBox(PBox,Extent(0,NuclearExpnt,NukeCo,TauPAC))
          HGKet(1)=Zero
          SPKetC(0)=Zero       
#ifdef PERIODIC
          QC(:)=GM%Carts%D(:,At)
          DO NC=1,CSMM1%NCells
!            Set Atomic Coordinates
             Prim%P(:)=QC(:)+CSMM1%CellCarts%D(:,NC)
             PBox%Center=Prim%P
!            Walk the walk
             CALL VWalk(PoleRoot)
          ENDDO
!         Reset the Atomic Coordinates
          Prim%P(:) = QC(:)
          PBox%Center=Prim%P
#else
!         Walk the walk
          CALL VWalk(PoleRoot)
#endif
!         Accumulate the atomic contribution
          NukE=NukE+NukeCo(1)*HGKet(1)+NukePole(1)*SPKetC(0) 
       ENDDO
     END FUNCTION NukE
!
END MODULE
