module NEB
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
  SUBROUTINE NEBInit(G,Reac,Prod)
    TYPE(Geometries) :: G
    TYPE(CRDS)       :: Reac,Prod
    TYPE(DBL_RNK2)   :: ReactionVector
    REAL(DOUBLE)     :: ImageFraction
    INTEGER          :: I
    !----------------------------------------------------------------------------

    NEBReac=Reac
    NEBProd=Prod

    write(*,*)'Into NEBInit'

    CALL New(ReactionVector,(/3,NEBReac%NAtms/))

    !Initialize each clone to initial state then interpolate Cartesian coordinates
    ReactionVector%D=NEBReac%Carts%D-NEBProd%Carts%D
    write(*,*)'Reacs'
    write(*,*) NEBReac%Carts%D
    write(*,*)'Prods'
    write(*,*) NEBProd%Carts%D
    write(*,*)'ReactionVector'
    write(*,*) ReactionVector%D

    DO I=1,G%Clones
       ImageFraction=REAL(I)/REAL(G%Clones)
       G%Clone(I)%Carts%D=Reac%Carts%D+ImageFraction*ReactionVector%D
       write(*,*)'Clone ',I
       write(*,*) G%Clone(I)%Carts%D
    ENDDO

    write(*,*)'Done NEBInit'
  END SUBROUTINE NEBInit


!===============================================================================
! Project out the force along the band and add spring forces along the band.
!===============================================================================
  SUBROUTINE NEBForce(G,O)
    TYPE(Geometries) :: G
    TYPE(Options)    :: O
    INTEGER          :: I,UMaxI,NAtms
    LOGICAL          :: UPm,UPp
    REAl(DOUBLE)     :: UMin,UMax,Um,Up,Rm,Rp
    TYPE(DBL_RNK2)   :: N
    !----------------------------------------------------------------------------

    write(*,*)'Into NEBForce'

    ! Allocate local vectors
    CALL New(N,(/3,NEBReac%NAtms/))

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
       IF(.NOT.(UPm==UPp))THEN
          ! If we are not at an extrema of energy, 
          ! the tangent is the vector to the lower energy neighbour
          IF(UPm)THEN
             N%D=G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D
          ELSE
             N%D=G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D
          ENDIF
       ELSE
          ! At an extrema of energy, 
          ! interpolate the tangent linearly with the energy
          Um=G%Clone(I-1)%ETotal-G%Clone(I)%ETotal
          Up=G%Clone(I+1)%ETotal-G%Clone(I)%ETotal
          UMin=MIN(ABS(Up),ABS(Um))
          UMax=MAX(ABS(Up),ABS(Um))
          IF(Um>Up)THEN
             N%D=(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMin
             N%D=N%D+(G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)*UMax
          ELSE
             N%D=(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMax
             N%D=N%D+(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMin
          ENDIF
       ENDIF
       N%D=N%D/SQRT(SUM(N%D**2))

       ! Project out the force along the tangent
       G%Clone(I)%Carts%D=G%Clone(I)%Carts%D-G%Clone(I)%Carts%D*SUM(G%Clone(I)%Carts%D*N%D)

       ! Add spring forces along the band
       IF(I==1)THEN
          Rm=SQRT(SUM((G%Clone(I)%Carts%D-NEBReac%Carts%D)**2))
       ELSE
          Rm=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)**2))
       ENDIF
       IF(I==G%Clones)THEN
          Rp=SQRT(SUM((G%Clone(I)%Carts%D-NEBProd%Carts%D)**2))
       ELSE
          Rp=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I+1)%Carts%D)**2))
       ENDIF
       G%Clone(I)%Vects%D=G%Clone(I)%Vects%D+O%NEBSpring*N%D*(Rp-Rm)
       
       ! If climbing image, zero forces along the band for that image
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          G%Clone(I)%Vects%D=G%Clone(I)%Vects%D-N%D*SUM(N%D*G%Clone(I)%Vects%D)
       ENDIF
    ENDDO
    write(*,*)'Done NEBForce'

  END SUBROUTINE NEBForce

!===============================================================================
! Generate a cubic spline along the band and interpolate to find extrema
!===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
