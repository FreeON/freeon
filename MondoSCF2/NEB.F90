MODULE NEB
  !===============================================================================
  ! Module for calculating reaction (minimum energy) paths between known 
  ! reactant and product states.  This module impliments the climbing image
  ! NEB so that the highest energy image will converge to a saddle point.
  !
  ! Module written by Graeme Henkelman and Matt Challacombe
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
  CONTAINS
  !===============================================================================
  ! Initialize the NEB by generating an linear interpolation between intial
  ! and final states. 
  !===============================================================================
  SUBROUTINE NEBInit(G)
    TYPE(Geometries) :: G
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms) :: ReactionVector
    REAL(DOUBLE)     :: ImageFraction
    INTEGER          :: iCLONE,j
    !----------------------------------------------------------------------------
    !Initialize each clone to initial state then interpolate Cartesian coordinates
    write(*,*)'NEB: Into NEBInit'
    ReactionVector=G%Clone(G%Clones+1)%AbCarts%D-G%Clone(0)%AbCarts%D
    iClone=0
    write(*,'(1A7,I3)')'Image',iClone
    write(*,'(3F13.5)') (G%Clone(iCLONE)%AbCarts%D(:,j),j=1,G%Clone(0)%NAtms)
    DO iCLONE=1,G%Clones
       ImageFraction=DBLE(iCLONE)/DBLE(G%Clones+1)
       CALL SetEq_CRDS(G%Clone(0),G%Clone(iCLONE))
       G%Clone(iCLONE)%AbCarts%D=G%Clone(0)%AbCarts%D+ImageFraction*ReactionVector
       write(*,'(1A7,I3)')'Image',iClone
       write(*,'(3F13.5)') (G%Clone(iCLONE)%AbCarts%D(:,j),j=1,G%Clone(0)%NAtms)
    ENDDO
    iClone=G%Clones+1
    write(*,'(1A7,I3)')'Image',iClone
    write(*,'(3F13.5)') (G%Clone(iCLONE)%AbCarts%D(:,j),j=1,G%Clone(0)%NAtms)
    write(*,*)'NEB: Done NEBInit'
  END SUBROUTINE NEBInit

  !===============================================================================
  ! Make a deep copy of the CRDS structure
  ! (This should move.  Also figure out PBC issue.) 
  !===============================================================================
  SUBROUTINE SetEq_CRDS(G1,G2)
    TYPE(CRDS) :: G1,G2
    G2%InAU=G1%InAU
    G2%NElec=G1%NElec
    G2%Ordrd=G1%Ordrd
    G2%Multp=G1%Multp
    G2%TotCh=G1%TotCh
    G2%NAlph=G1%NAlph
    G2%NBeta=G1%NBeta
    G2%Carts%D=G1%Carts%D
    G2%AbCarts%D=G1%AbCarts%D
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
    INTEGER          :: I,j,UMaxI,NAtms
    LOGICAL          :: UPm,UPp
    REAl(DOUBLE)     :: UMin,UMax,Um,Up,Rm,Rp,Dist,FProj
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms) :: N
    !----------------------------------------------------------------------------
    write(*,*)'NEB: Into NEBForce'
    Dist=0
    ! Find the image with the maximum total energy
    UMax=G%Clone(1)%ETotal
    UMaxI=1
    DO I=2,G%Clones
       IF(G%Clone(I)%ETotal>UMax) THEN
          UMaxI=I
          UMax=G%Clone(I)%ETotal
       ENDIF
    ENDDO
    write(*,*)'NEB: Found max energy image, ',UMaxI
    ! Find the tangent to the path at each image
    ! Project out potential forces along the band
    ! Add spring forces along the band

!GH    write(*,*)'React Crds'
!GH    write(*,'(3F13.5)') (G%Clone(0)%AbCarts%D(:,j),j=1,G%Clone(0)%NAtms)

!GH    write(*,*)'Prod Crds'
!GH    write(*,'(3F13.5)') (G%Clone(G%Clones+1)%AbCarts%D(:,j),j=1,G%Clone(0)%NAtms)
    write(*,*)'NEB: Distance, Energies, and Forces'


    write(*,'(A,I5,3F13.5)') 'NEB: ',0,Dist,G%Clone(0)%ETotal
    DO I=1,G%Clones
       ! Are the neighboring images higher in energy?
       IF(I==1)THEN
          UPm=.FALSE.
       ELSE
          UPm=G%Clone(I-1)%ETotal>G%Clone(I)%ETotal
       ENDIF
       IF(I==G%Clones)THEN
          UPp=.FALSE.
       ELSE
          UPp=G%Clone(I+1)%ETotal>G%Clone(I)%ETotal
       ENDIF
       IF(UPm.NEQV.UPp)THEN 
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
       N=N/SQRT(SUM(N**2))

!GH       write(*,*)'Crds'
!GH       write(*,'(3F13.5)') (G%Clone(I)%AbCarts%D(:,j),j=1,G%Clone(0)%NAtms)
!GH       write(*,*)'Normal'
!GH       write(*,'(3F13.5)') (N(:,j),j=1,G%Clone(0)%NAtms)
!GH       write(*,*)'In Force'
!GH       write(*,'(3F13.5)') (G%Clone(I)%Gradients%D(:,j),j=1,G%Clone(0)%NAtms)

       ! Project out the force along the tangent, unless this is the
       ! climbing image, for which the force along the tangent is inverted
       FProj=SUM(G%Clone(I)%Gradients%D*N)
!GH       write(*,*)'Prj Force'
       IF(O%NEBClimb.AND.I==UMaxI)THEN
!GH          write(*,*)'Climbing Image',I
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D-2.0*N*FProj
       ELSE
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D-N*FProj
       ENDIF
!GH       write(*,'(3F13.5)') (G%Clone(I)%Gradients%D(:,j),j=1,G%Clone(0)%NAtms)

       ! Add spring forces along the band (if we are not the climbing image)
       Rm=SQRT(SUM((G%Clone(I)%AbCarts%D-G%Clone(I-1)%AbCarts%D)**2))
       Rp=SQRT(SUM((G%Clone(I)%AbCarts%D-G%Clone(I+1)%AbCarts%D)**2))
!GH       write(*,*)'Rm,Rp',Rm,Rp
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          ! Do nothing (no springs for the climbing image)
       ELSE
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D+O%NEBSpring*N*(Rp-Rm)
       ENDIF

!GH       write(*,*)'Out Force'
!GH       write(*,'(3F13.5)') (G%Clone(I)%Gradients%D(:,j),j=1,G%Clone(0)%NAtms)
       ! Write distance, energies and forces
       Dist=Dist+Rm
       write(*,'(A,I5,3F13.5)') 'NEB: ',I,Dist,G%Clone(I)%ETotal,FProj
    ENDDO
    Rm=SQRT(SUM((G%Clone(G%Clones+1)%AbCarts%D-G%Clone(G%Clones)%AbCarts%D)**2))
    Dist=Dist+Rm
    write(*,'(A,I5,3F13.5)') 'NEB: ',G%Clones+1,Dist,G%Clone(G%Clones+1)%ETotal
    write(*,*)'NEB: Done NEBForce'

  END SUBROUTINE NEBForce

  !===============================================================================
  ! Generate a cubic spline along the band and interpolate to find extrema
  !===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
