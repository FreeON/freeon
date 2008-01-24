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
!    COMPUTE THE NUCLEAR-TOTAL ELECTROSTATIC ENERGY IN O(N Lg N) CPU TIME
!    Authors: Matt Challacombe and C.J. Tymczak
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
  USE PBCFarField
#ifdef PARALLEL
  USE ParallelQCTC
#endif
  IMPLICIT NONE
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    FUNCTION NukE(GMLoc)
      TYPE(CRDS),  INTENT(IN)         :: GMLoc
      REAL(DOUBLE)                    :: NukE,NukeCo,NukePole,PExtent
      REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
      REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
      INTEGER                         :: NC
      REAL(DOUBLE),DIMENSION(3)       :: PTmp
!---------------------------------------------------------------------------------------------
      NukE=Zero 
      DO At=1,GMLoc%Natms 
#ifdef PARALLEL
      IF(At >= BegAtInd%I(MyID) .AND. At <= EndAtInd%I(MyID)) THEN
#endif
         IF(GMLoc%AtNum%D(At)<105.D0)THEN
!           Initialize |BRA>
            HGBra(1) =-GMLoc%AtNum%D(At)*(NuclearExpnt/Pi)**(ThreeHalves)
            SPBraC(0)=-GMLoc%AtNum%D(At)
            Prim%Ell = 0
            Prim%P   = GMLoc%Carts%D(:,At)
            Prim%Zeta= NuclearExpnt
!           Set the MAC
            DP2 = (FudgeFactorial(0,SPELL+1)*ABS(GMLoc%AtNum%D(At))/TauMAC)**(Two/DBLE(SPEll+2))
            DP2 = MIN(1.D10,DP2)
!           Set the PAC
            PrimWCoef = ABS(GMLoc%AtNum%D(At))
!           Initialize <KET|
            CALL SetKet(Prim,PExtent)
            PTmp=GMLoc%Carts%D(:,At)
            DO NC=1,CS_IN%NCells
!              Set Atomic Coordinates
               Prim%P=PTmp+CS_IN%CellCarts%D(:,NC) 
               PBox%Center= Prim%P
!              Walk the walk
               CALL VWalk(PoleRoot)
            ENDDO
!           Reset the Atomic Coordinates
            Prim%P=PTmp
!           Accumulate the atomic contribution
            NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!           Add in the Far Field, Dipole and Quadripole  Correction 
            IF(GMLoc%PBC%Dimen > 0) THEN
               NukE = NukE + CTraxFF(Prim,HGBra,GMLoc)
            ENDIF
         ENDIF
#ifdef PARALLEL
      ENDIF
#endif
      ENDDO
!
    END FUNCTION NukE
!
END MODULE
