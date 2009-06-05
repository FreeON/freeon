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

#define NEB_DEBUG

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
  USE ls_rmsd
  USE PrettyPrint
  USE ControlStructures
  USE Order
  USE MondoLogger

  IMPLICIT NONE

  SAVE

CONTAINS
  !===============================================================================
  ! Initialize the NEB by generating an linear interpolation between intial
  ! and final states.
  !===============================================================================
  SUBROUTINE NEBInit(G)
    TYPE(Geometries)                            :: G
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms)  :: ReactionVector
    REAL(DOUBLE)                                :: ImageFraction
    INTEGER                                     :: iCLONE,j
    CHARACTER(LEN=DEFAULT_CHR_LEN)              :: Message

    !----------------------------------------------------------------------------
    !Initialize each clone to initial state then interpolate Cartesian coordinates
#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "starting...")
#endif

    ReactionVector=G%Clone(G%Clones+1)%Carts%D-G%Clone(0)%Carts%D
    iClone=0
#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Image "//TRIM(IntToChar(iCLONE)))
    WRITE(Message,'(3F13.5)') (G%Clone(iCLONE)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
    CALL MondoLog(DEBUG_NONE, "NEBInit", Message)
#endif
    DO iCLONE=1,G%Clones
       ImageFraction=DBLE(iCLONE)/DBLE(G%Clones+1)
       CALL SetEq_CRDS(G%Clone(0),G%Clone(iCLONE))
       G%Clone(iCLONE)%Carts%D=G%Clone(0)%Carts%D+ImageFraction*ReactionVector
#ifdef NEB_DEBUG
       CALL MondoLog(DEBUG_NONE, "NEBInit", "Image "//TRIM(IntToChar(iCLONE)))
       WRITE(Message,'(3F13.5)') (G%Clone(iCLONE)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
       CALL MondoLog(DEBUG_NONE, "NEBInit", Message)
#endif
    ENDDO
    iClone=G%Clones+1
#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Image "//TRIM(IntToChar(iCLONE)))
    write(Message,'(3F13.5)') (G%Clone(iCLONE)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
    CALL MondoLog(DEBUG_NONE, "NEBInit", Message)
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Done NEBInit")
#endif
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
    G2%Carts%D=G1%Carts%D
    G2%NAtms=G1%NAtms
    G2%Nkind=G1%Nkind
    G2%AtNum%D=G1%AtNum%D
    G2%AtMss%D=G1%AtMss%D
    G2%AtNam%C=G1%AtNam%C
    G2%AtTyp%I=G1%AtTyp%I
    G2%CConstrain%I=G1%CConstrain%I
    !    CALL SetEq_PBCInfo(G1%PBC,G2%PBC)
  END SUBROUTINE SetEq_CRDS

  SUBROUTINE NEBPurify(G,Init_O,Print_O)
    TYPE(Geometries)                      :: G
    LOGICAL,OPTIONAL                      :: Init_O,Print_O
    LOGICAL                               :: Init
    INTEGER                               :: I,iCLONE,bCLONE,eCLONE,nCLONE,J
    TYPE(DBL_RNK2), DIMENSION(G%Clones+1) :: GTmp
    INTEGER, DIMENSION(G%Clones+1)        :: I2
    REAL(DOUBLE),DIMENSION(G%Clones+1)    :: R2
    REAL(DOUBLE),DIMENSION(3,3)           :: U
    REAL(DOUBLE),DIMENSION(3)             :: Center1,Center2
    REAL(DOUBLE)                          :: Error
    CHARACTER(LEN=4*DCL)                  :: Mssg

    IF(PRESENT(Init_O))THEN
       Init=.TRUE.
    ELSE
       Init=.FALSE.
    ENDIF

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEB", "Init = "//TRIM(LogicalToChar(Init)))
#endif

    IF(Init)THEN
       bCLONE=G%Clones+1
       eCLONE=G%Clones+1
       ! Check for stupid input
       DO I=1,G%Clone(0)%NAtms
          IF(G%Clone(0)%AtNum%D(I).NE.G%Clone(G%Clones+1)%AtNum%D(I))THEN
             CALL MondoHalt(NEBS_ERROR,'Ordering of Reactant and Product is different!')
          ENDIF
       ENDDO
    ELSE
       bCLONE=1
       eCLONE=G%Clones+1
    ENDIF

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEB", "bCLONE = "//TRIM(IntToChar(bCLONE)))
    CALL MondoLog(DEBUG_NONE, "NEB", "eCLONE = "//TRIM(IntToChar(eCLONE)))
    CALL MondoLog(DEBUG_NONE, "NEB", "CLONE zero before anything = ")
    CALL PPrint(G%Clone(0),Unit_O=6,PrintGeom_O='XYZ')
#endif

!!$    ! Scale the coordinates by Z
!!$    DO I=1,G%Clone(0)%NAtms
!!$       G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)*G%Clone(0)%AtNum%D(I)
!!$    ENDDO

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEB", "CLONE ZERO AFTER SCALING = ")
    CALL PPrint(G%Clone(0),Unit_O=6,PrintGeom_O='XYZ')
#endif

    ! Constraints over-ride purification
    DO I=1,G%Clone(0)%NAtms
       IF(G%Clone(0)%CConstrain%I(I)/=0)THEN
          ! No RMSD alignment with constraints
          GOTO 101 ! can still re-order based on RMSD
       ENDIF
    ENDDO

    ! Translate and rotate each clone to minimize the rmsd relative to clone zero
    DO iCLONE=bCLONE,eCLONE
#ifdef NEB_DEBUG
      CALL MondoLog(DEBUG_NONE, "NEB", "purifying clone "//TRIM(IntToChar(iclone)))
#endif

!!$       ! Scale the coordinates by Z
!!$       DO I=1,G%Clone(0)%NAtms
!!$          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)*G%Clone(iCLONE)%AtNum%D(I)
!!$       ENDDO

       ! Find the transformation that minimizes the RMS deviation between the
       ! reactants (clone 0), the clones (1-N) and the products (N+1)
       CALL RMSD(G%Clone(0)%NAtms,G%Clone(iCLONE)%Carts%D,G%Clone(0)%Carts%D,  &
            1, U, Center2, Center1, error )! , calc_g, grad)
#ifdef NEB_DEBUG
       WRITE(*,333)1,Center1
       WRITE(*,333)2,Center2
       WRITE(*,333)-1,-(Center2-Center1)
333    FORMAT('Center',I2,' = ',3(F10.5,', '))
#endif
       IF(Init)THEN
          ! Translate the reactants JUST ONCE to C1
          DO I=1,G%Clone(0)%NAtms
             G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)-Center1
          ENDDO
       ENDIF
#ifdef NEB_DEBUG
       WRITE(*,*)' CLONE = ',0,' AFTER TRANSLATION '
       CALL PPrint(G%Clone(0),Unit_O=6,PrintGeom_O='XYZ')
       WRITE(*,*)' CLONE = ',iCLONE,' BEFORE TRANSLATION'
       CALL PPrint(G%Clone(iCLONE),Unit_O=6,PrintGeom_O='XYZ')
#endif
       ! Translation to C2 ...
       DO I=1,G%Clone(0)%NAtms
          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)-Center2
       ENDDO
       ! ... and rotation
       DO I=1,G%Clone(0)%NAtms
          G%Clone(iCLONE)%Carts%D(:,I)=MATMUL(U,G%Clone(iCLONE)%Carts%D(:,I))
       ENDDO
#ifdef NEB_DEBUG
       WRITE(*,*)' CLONE = ',iCLONE,' AFTER TRANSLATION AND ROTATION '
       CALL PPrint(G%Clone(iCLONE),Unit_O=6,PrintGeom_O='XYZ')
#endif
    ENDDO

!!$    ! Un-scale the coordinates by Z
!!$    DO I=1,G%Clone(0)%NAtms
!!$       G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)/G%Clone(0)%AtNum%D(I)
!!$    ENDDO
!!$    DO iCLONE=bCLONE,eCLONE
!!$       DO I=1,G%Clone(0)%NAtms
!!$          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)/G%Clone(iCLONE)%AtNum%D(I)
!!$       ENDDO
!!$    ENDDO
!    IF(PRESENT(Print_O))THEN

101 IF(Init) RETURN

    ! Compute RMSD from first clone
    R2(:)=Zero
    DO iCLONE=bCLONE,eCLONE
       R2(iCLONE)=Zero
       DO I=1,G%Clone(0)%NAtms
          DO J=1,3
             R2(iCLONE)=R2(iCLONE)+(G%Clone(iCLONE)%Carts%D(J,I)-G%Clone(0)%Carts%D(J,I))**2
          ENDDO
       ENDDO
       R2(iCLONE)=SQRT(R2(iCLONE))/G%Clone(0)%NAtms
    ENDDO

    ! Order based on RMSD
    nCLONE=eCLONE-bCLONE+1
    DO I=1,nCLONE
       I2(I)=I
       CALL New(GTmp(I),(/3,G%Clone(0)%NAtms/))
       GTmp(I)%D=G%Clone(I)%Carts%D
    ENDDO
    CALL DblIntSort77(nCLONE,R2,I2,2)
    WRITE(*,*)R2
    WRITE(*,*)I2

    DO I=1,nCLONE
       G%Clone(I)%Carts%D=GTmp(I2(I))%D
    ENDDO
    DO I=1,nCLONE
       CALL Delete(GTmp(I))
    ENDDO

    Mssg='RMSDs = '
    DO I=1,G%Clones
       Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(R2(I)))//','
    ENDDO
    Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(R2(G%Clones+1)))
    !CALL MondoLog(DEBUG_NONE, "NEBPurify", TRIM(Mssg))
!   ENDIF

    CALL MondoLog(DEBUG_NONE, "FreeON", Mssg, "NEBPurify("//TRIM(IntToChar(G%Clone(1)%Confg))//')')
  END SUBROUTINE NEBPurify
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
    CHARACTER(LEN=DCL) :: Mssg
    !----------------------------------------------------------------------------
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
    ! Find the tangent to the path at each image
    ! Project out potential forces along the band
    ! Add spring forces along the band

    !GH    write(*,*)'React Crds'
    !GH    write(*,'(3F13.5)') (G%Clone(0)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)

    !GH    write(*,*)'Prod Crds'
    !GH    write(*,'(3F13.5)') (G%Clone(G%Clones+1)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)

    CALL MondoLog(DEBUG_NONE, "NEBForce", &
          "Dist = "//TRIM(FltToShrtChar(Dist)) &
      //', E = '//TRIM(DblToMedmChar(G%Clone(0)%ETotal)), "Reactant")

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
             N=G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D
          ELSE
             N=G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D
          ENDIF
       ELSE
          ! At an extrema of energy,
          ! interpolate the tangent linearly with the energy
          Um=G%Clone(I-1)%ETotal-G%Clone(I)%ETotal
          Up=G%Clone(I+1)%ETotal-G%Clone(I)%ETotal
          UMin=MIN(ABS(Up),ABS(Um))
          UMax=MAX(ABS(Up),ABS(Um))
          IF(Um>Up)THEN
             N=(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMin
             N=N+(G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)*UMax
          ELSE
             N=(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMax
             N=N+(G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)*UMin
          ENDIF
       ENDIF
       N=N/SQRT(SUM(N**2))

       !GH       write(*,*)'Crds'
       !GH       write(*,'(3F13.5)') (G%Clone(I)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
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
       Rm=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)**2))
       Rp=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I+1)%Carts%D)**2))
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

       CALL MondoLog(DEBUG_NONE, "NEBForce", &
                "Rm = "//TRIM(FltToShrtChar(Rm)) &
            //", Rp = "//TRIM(FltToShrtChar(Rp)) &
            //", Dist = "//TRIM(FltToShrtChar(Dist)) &
            //', E = '//TRIM(DblToMedmChar(G%Clone(I)%ETotal)) &
            //', F = '//TRIM(DblToMedmChar(FProj)), "Image "//TRIM(IntToChar(I)))
    ENDDO

    Rm=SQRT(SUM((G%Clone(G%Clones+1)%Carts%D-G%Clone(G%Clones)%Carts%D)**2))
    Dist=Dist+Rm
    CALL MondoLog(DEBUG_NONE, "NEBForce", &
             "Dist = "//TRIM(FltToShrtChar(Dist)) &
         //', E = '//TRIM(DblToMedmChar(G%Clone(G%Clones+1)%ETotal)), "Product")
  END SUBROUTINE NEBForce

  !===============================================================================
  ! Generate a cubic spline along the band and interpolate to find extrema
  !===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
