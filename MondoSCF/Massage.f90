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
!============================================================================
!                 ** THE MYSTERIES OF MONDO MASSAGE **
!
! WHY DO WE NEED TO TRANSLATE TO GET THE RIGHT NUMBERS??????
! SHOULDNT WRAPPING TAKE CARE OF IT ALL ?????
! ALSO, WHY DOES WRAPPING APPARENTLY NOT BRING THE ATOMS BACK INTO THE UC????
! AND WHY IS ALL THE PBC STUFF LIKE WRAPPING AND TRANSLATION IN ATOM PAIRS INSTEAD
! OF PBC?  AND DO WE REALLY NEED FOUR COORDINATE ARRAYS IN CRDS?
!============================================================================
MODULE Massage
  USE Parse
  USE InOut
  USE AtomPairs
  USE OptionKeys
  USE DynamicsKeys
  USE GeometryKeys
  USE ControlStructures
  USE MondoLogger

  IMPLICIT NONE

CONTAINS
  !============================================================================
  ! ALL REORDERING, RESCALING, WRAPPING AND TRANSLATING OF COORDINATES OCCURS
  ! HERE AND NO WHERE ELSE!
  !============================================================================
  SUBROUTINE MassageCoordinates(O,G,P)
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
    INTEGER          :: I,J,GBeg,GEnd

    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
      GBeg=0
      GEnd=G%Clones+1
    ELSE
      GBeg=1
      GEnd=G%Clones
    ENDIF
    DO I=GBeg,GEnd
      CALL ToAtomicUnits(G%Clone(I))
      CALL PeriodicXLate(G%Clone(I))
    ENDDO
  END SUBROUTINE MassageCoordinates
  !============================================================================
  ! RESCALE VALUES IF ORIGINALLY IN ANGSTROMS
  !============================================================================
  SUBROUTINE ToAtomicUnits(G)
    TYPE(CRDS) :: G

    IF(G%InAU) THEN
      CALL MondoLog(DEBUG_NONE, "ToAtomicUnits", "coordinates are already in atomic units")
      RETURN
    ENDIF

    CALL MondoLog(DEBUG_NONE, "ToAtomicUnits", "converting coordinates into atomic units")
    G%Carts%D    = AngstromsToAU*G%Carts%D
    G%Velocity%D = AngstromsToAU*G%Velocity%D

    G%PBC%CellCenter%D = G%PBC%CellCenter%D*AngstromsToAU
    G%PBC%BoxShape%D   = AngstromsToAU*G%PBC%BoxShape%D
    G%PBC%InvBoxSh%D   = G%PBC%InvBoxSh%D/AngstromsToAU
    G%PBC%CellVolume   = G%PBC%CellVolume*AngstromsToAU**G%PBC%Dimen
    G%PBC%DipoleFAC    = G%PBC%DipoleFAC/(AngstromsToAU**G%PBC%Dimen)
    G%PBC%QupoleFAC    = G%PBC%QupoleFAC/(AngstromsToAU**G%PBC%Dimen)

  END SUBROUTINE ToAtomicUnits
  !============================================================================
  !
  !============================================================================
  SUBROUTINE PeriodicXLate(G)
    TYPE(CRDS)                  :: G
    REAL(DOUBLE),DIMENSION(1:3) :: CMVec
    INTEGER                     :: I
    !-------------------------------------------------------------------------!
    IF(G%PBC%Translate) THEN
      CMVec(:)=Zero
      DO I=1,G%NAtms
        CMVec(:)=CMVec(:)+G%BoxCarts%D(:,I)
      ENDDO
      CMVec(:)=Half-CMVec(:)/DBLE(G%NAtms)
      G%PBC%TransVec%D(:)=FracToAtom(G,CMVec(:))
      DO I=1,3
        IF(G%PBC%AutoW%I(I)==0) G%PBC%TransVec%D(I)=Zero
      ENDDO
      CALL Translate(G,G%PBC%TransVec%D)
    ELSE
      G%PBC%TransVec%D(:)=Zero
    ENDIF
    !
  END SUBROUTINE PeriodicXLate
END MODULE Massage

