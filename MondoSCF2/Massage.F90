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
    !-------------------------------------------------------------------------!
    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
       GBeg=0
       GEnd=G%Clones+1
    ELSE	
       GBeg=1
       GEnd=G%Clones
    ENDIF
    DO I=GBeg,GEnd
       CALL ToAtomicUnits(G%Clone(I))
#ifdef PERIODIC
       CALL PeriodicXLate(G%Clone(I))
#endif
    ENDDO
  END SUBROUTINE MassageCoordinates
  !============================================================================
  ! RESCALE VALUES IF ORIGINALLY IN ANGSTROMS
  !============================================================================
  SUBROUTINE ToAtomicUnits(G)
    TYPE(CRDS) :: G
    IF(G%InAU)RETURN
    G%InAU=.TRUE.
    G%Carts%D=AngstromsToAU*G%Carts%D               
    G%AbCarts%D=AngstromsToAU*G%AbCarts%D               
#ifdef PERIODIC
    G%PBC%BoxShape=AngstromsToAU*G%PBC%BoxShape
    G%PBC%InvBoxSh=G%PBC%InvBoxSh/AngstromsToAU
    G%PBC%CellVolume=G%PBC%CellVolume*AngstromsToAU**G%PBC%Dimen
    G%PBC%CellCenter=G%PBC%CellCenter*AngstromsToAU
    G%PBC%DipoleFac=G%PBC%DipoleFac/AngstromsToAU**G%PBC%Dimen
    G%PBC%QupoleFac=G%PBC%QupoleFac/AngstromsToAU**G%PBC%Dimen
#endif
  END SUBROUTINE ToAtomicUnits
  !============================================================================
  ! 
  !============================================================================
  SUBROUTINE PeriodicXLate(G)
    TYPE(CRDS)                  :: G
    REAL(DOUBLE),DIMENSION(1:3) :: CMVec
    INTEGER                     :: I
    !-------------------------------------------------------------------------!
    CMVec(:)=Zero
    DO I=1,G%NAtms
       CMVec(:)=CMVec(:)+G%BoxCarts%D(:,I)
    ENDDO
    CMVec(:)=Half-CMVec(:)/DBLE(G%NAtms)
    G%PBC%TransVec(:)=FracToAtom(G,CMVec(:))
    DO I=1,3
       IF(.NOT.G%PBC%AutoW(I))G%PBC%TransVec(I)=Zero
    ENDDO
    CALL Translate(G,G%PBC%TransVec)
  END SUBROUTINE PeriodicXLate
END MODULE Massage

