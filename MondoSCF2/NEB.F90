MODULE NEB
  !===============================================================================
  ! Module for calculating reaction (minimum energy) paths between known 
  ! reactant and product states.  This module impliments the climbing image
  ! NEB so that the highest energy image will converge to a saddle point.
  !
  ! Module written by Graeme Henkelman
  ! Email: graeme@lanl.gov
  !
  ! NEB References:
  !  H. Jonsson, G. Mills, and K.W. Jacobsen, "Nudged elastic band method 
  !    for finding minimum energy paths of transitions," in Classical and
  !    Quantum Dynamics in Condensed Phase Simulations, edited by B.J.Berne,
  !    G. Ciccotti, and D. F. Coker (World Scientific, Singapore, 1998), p.385
  !  G. Henkelman, B.P. Uberuaga, and H. Jonsson, "A climbing-image NEB method
  !     for finding saddle points and minimum energy paths", J. Chem. Phys,
  !     v113, 9901 (2000).
  !        
  !===============================================================================
  USE InOut
  USE PrettyPrint
  USE ControlStructures
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC :: NEBInit,NEBForce
  TYPE(Crds) :: NEBReac,NEBProd      ! Reactant and Product structures for the NEB
CONTAINS
  !===============================================================================
  ! Initialize the NEB by generating an linear interpolation between intial
  ! and final states. 
  !===============================================================================
  SUBROUTINE NEBInit(G)
    TYPE(Geometries) :: G
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms) :: ReactionVector
    REAL(DOUBLE)     :: ImageFraction
    INTEGER          :: iCLONE,I
    !----------------------------------------------------------------------------
    !Initialize each clone to initial state then interpolate Cartesian coordinates
    ReactionVector=G%Clone(0)%AbCarts%D-G%Clone(G%Clones+1)%AbCarts%D
    DO iCLONE=1,G%Clones
       ImageFraction=DBLE(iCLONE)/DBLE(G%Clones)
       CALL SetEq_CRDS(G%Clone(0),G%Clone(iCLONE))
       G%Clone(iCLONE)%AbCarts%D=G%Clone(0)%AbCarts%D+ImageFraction*ReactionVector
    ENDDO
  END SUBROUTINE NEBInit

  SUBROUTINE SetEq_CRDS(G1,G2)
    TYPE(CRDS) :: G1,G2
    G2%NElec=G1%NElec
    G2%Ordrd=G1%Ordrd
    G2%Multp=G1%Multp
    G2%TotCh=G1%TotCh
    G2%NAlph=G1%NAlph
    G2%NBeta=G1%NBeta
    G2%Carts%D=G1%Carts%D
    G2%AbCarts%D=G2%AbCarts%D
    G2%NAtms=G1%NAtms      
    G2%Nkind=G1%Nkind 
    G2%AtNum%D=G1%AtNum%D 
    G2%AtMss%D=G1%AtMss%D 
    G2%AtNam%C=G1%AtNam%C 
    G2%AtTyp%I=G1%AtTyp%I      
    G2%CConstrain%I=G1%CConstrain%I 
!    CALL SetEq_PBCInfo(G1%PBC,G2%PBC)
  END SUBROUTINE SetEq_CRDS
  !===============================================================================
  ! Project out the force along the band and add spring forces along the band.
  !===============================================================================
  SUBROUTINE NEBForce(G,O)
    TYPE(Geometries) :: G
    TYPE(Options)    :: O
    INTEGER          :: I,UMaxI,NAtms
    LOGICAL          :: UPm,UPp
    REAl(DOUBLE)     :: UMin,UMax,Um,Up,Rm,Rp
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms) :: N
    !----------------------------------------------------------------------------
    write(*,*)'Into NEBForce'
    ! Find the image with the maximum total energy
    UMax=G%Clone(1)%ETotal
    UMaxI=1
    DO I=2,G%Clones
       IF(G%Clone(I)%ETotal>UMax) THEN
          UMaxI=I
          UMax=G%Clone(I)%ETotal
       ENDIF
    ENDDO
    ! Find the tangent to the path at each image
    ! Project out potential forces along the band
    ! Add spring forces along the band
    DO I=1,G%Clones
       ! Are the neighboring images higher in energy?
       IF(I==1)THEN
          UPm=.FALSE.
       ELSE
          UPm=G%Clone(I-1)%ETotal>G%Clone(I)%ETotal
       ENDIF
       IF(I==G%Clones)THEN
          UPp=.FALSE.
       ELSEIF(I==G%Clones)THEN
          UPp=G%Clone(I+1)%ETotal>G%Clone(I)%ETotal
       ENDIF
       IF((UPm.AND..NOT.UPp).OR.(.NOT.UPm.AND.UPp))THEN   ! IF(.NOT.(UPm==UPp))THEN 
          ! If we are not at an extrema of energy, 
          ! the tangent is the vector to the lower energy neighbour
          IF(UPm)THEN
             N=G%Clone(I)%AbCarts%D-G%Clone(I-1)%AbCarts%D
          ELSE
             N=G%Clone(I+1)%AbCarts%D-G%Clone(I)%AbCarts%D
          ENDIF
       ELSE
          ! At an extrema of energy, 
          ! interpolate the tangent linearly with the energy
          Um=G%Clone(I-1)%ETotal-G%Clone(I)%ETotal
          Up=G%Clone(I+1)%ETotal-G%Clone(I)%ETotal
          UMin=MIN(ABS(Up),ABS(Um))
          UMax=MAX(ABS(Up),ABS(Um))
          IF(Um>Up)THEN
             N=(G%Clone(I+1)%AbCarts%D-G%Clone(I)%AbCarts%D)*UMin
             N=N+(G%Clone(I)%AbCarts%D-G%Clone(I-1)%AbCarts%D)*UMax
          ELSE
             N=(G%Clone(I+1)%AbCarts%D-G%Clone(I)%AbCarts%D)*UMax
             N=N+(G%Clone(I+1)%AbCarts%D-G%Clone(I)%AbCarts%D)*UMin
          ENDIF
       ENDIF
       !????????????????????????
       N=N/SQRT(SUM(N**2))
       ! Project out the force along the tangent
       G%Clone(I)%AbCarts%D=G%Clone(I)%AbCarts%D-G%Clone(I)%AbCarts%D*SUM(G%Clone(I)%AbCarts%D*N)
       ! Add spring forces along the band
       IF(I==1)THEN
          Rm=SQRT(SUM((G%Clone(I)%AbCarts%D-G%Clone(0)%AbCarts%D)**2))
       ELSE
          Rm=SQRT(SUM((G%Clone(I)%AbCarts%D-G%Clone(I-1)%AbCarts%D)**2))
       ENDIF
       IF(I==G%Clones)THEN
          Rp=SQRT(SUM((G%Clone(I)%AbCarts%D-G%Clone(G%Clones+1)%AbCarts%D)**2))
       ELSE
          Rp=SQRT(SUM((G%Clone(I)%AbCarts%D-G%Clone(I+1)%AbCarts%D)**2))
       ENDIF
       G%Clone(I)%Vects%D=G%Clone(I)%Vects%D+O%NEBSpring*N*(Rp-Rm)
       ! If climbing image, zero forces along the band for that image
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          G%Clone(I)%Vects%D=G%Clone(I)%Vects%D-N*SUM(N*G%Clone(I)%Vects%D)
       ENDIF
    ENDDO
  END SUBROUTINE NEBForce
  !===============================================================================
  ! Generate a cubic spline along the band and interpolate to find extrema
  !===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
