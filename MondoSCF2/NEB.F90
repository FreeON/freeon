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
  IMPLICIT NONE  
  SAVE
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
    !write(*,*)'NEB: Into NEBInit'
    ReactionVector=G%Clone(G%Clones+1)%Carts%D-G%Clone(0)%Carts%D
    iClone=0
    !write(*,'(1A7,I3)')'Image',iClone
    !write(*,'(3F13.5)') (G%Clone(iCLONE)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
    DO iCLONE=1,G%Clones
       ImageFraction=DBLE(iCLONE)/DBLE(G%Clones+1)
       CALL SetEq_CRDS(G%Clone(0),G%Clone(iCLONE))
       G%Clone(iCLONE)%Carts%D=G%Clone(0)%Carts%D+ImageFraction*ReactionVector
       !write(*,'(1A7,I3)')'Image',iClone
       !write(*,'(3F13.5)') (G%Clone(iCLONE)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
    ENDDO
    iClone=G%Clones+1
    !write(*,'(1A7,I3)')'Image',iClone
    !write(*,'(3F13.5)') (G%Clone(iCLONE)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
    !write(*,*)'NEB: Done NEBInit'
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
    TYPE(Geometries) :: G
    LOGICAL,OPTIONAL :: Init_O,Print_O
    LOGICAL          :: Init
    INTEGER          :: I,iCLONE,bCLONE,eCLONE,J
    REAL(DOUBLE),DIMENSION(G%Clones+1) :: R2
    REAL(DOUBLE),DIMENSION(3,3) :: U
    REAL(DOUBLE),DIMENSION(3)   :: Center1,Center2
    REAL(DOUBLE) :: Error
    CHARACTER(LEN=4*DCL) :: Mssg
    IF(PRESENT(INIT_O))THEN
       INIT=.TRUE.
    ELSE
       INIT=.FALSE.
    ENDIF
    IF(INIT)THEN
       bCLONE=G%Clones+1
       eCLONE=G%Clones+1
       ! Check for stupid input
       DO I=1,G%Clone(0)%NAtms
          IF(G%Clone(0)%AtNum%D(I).NE.G%Clone(G%Clones+1)%AtNum%D(I))THEN
             CALL MondoHalt(NEBS_ERROR,'Ordering of Reactant and Product is different! ')
          ENDIF
       ENDDO
    ELSE
       bCLONE=1
       eCLONE=G%Clones+1
    ENDIF
    !
!!$    ! Scale the coordinates by Z 
!!$    DO I=1,G%Clone(0)%NAtms
!!$       G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)*G%Clone(0)%AtNum%D(I)
!!$    ENDDO


#ifdef NEB_DEBUG       
    WRITE(*,*)' bCLONE = ',bCLONE
    WRITE(*,*)' eCLONE = ',eCLONE
    WRITE(*,*)' CLONE ZERO AFTER SCALING = '
    CALL PPrint(G%Clone(0),Unit_O=6,PrintGeom_O='XYZ')
#endif
    !
    DO iCLONE=bCLONE,eCLONE
#ifdef NEB_DEBUG       
       WRITE(*,*)'==========',iclone,'============='
#endif



!!$       ! Scale the coordinates by Z
!!$       DO I=1,G%Clone(0)%NAtms
!!$          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)*G%Clone(iCLONE)%AtNum%D(I)
!!$       ENDDO
       ! Find the transformation that minimizes the RMS deviation between the 
       ! reactants (clone 0), the clones (1-N) and the products (N+1)
       CALL RMSD(G%Clone(0)%NAtms,G%Clone(iCLONE)%Carts%D,G%Clone(0)%Carts%D,  &
            1, U, center2, center1, error )! , calc_g, grad)
#ifdef NEB_DEBUG       
       WRITE(*,333)1,center1
       WRITE(*,333)2,center2
       WRITE(*,333)-1,-(center2-center1)
333    FORMAT('CENTER',I2,' = ',3(F10.5,', '))
#endif
       IF(INIT)THEN
          ! Translate the reactants JUST ONCE to C1
          DO I=1,G%Clone(0)%NAtms
             G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)-center1
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
          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)-center2
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
       Mssg=ProcessName('MondoSCF','NEBPurify')//'RMSDs = '
       DO I=1,G%Clones
          IF(MOD(I,4)==0)THEN
             Mssg=TRIM(Mssg)//RTRN//ProcessName() &
                  //'          '//TRIM(DblToShrtChar(R2(I)))//','
          ELSE
             Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(R2(I)))//','
          ENDIF
       ENDDO
       Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(R2(G%Clones+1)))
       WRITE(*,*)TRIM(Mssg)
!    ENDIF
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

    Mssg=ProcessName('NEBForce','Reactant')//' D = '//TRIM(FltToShrtChar(Dist)) &
         //', E = '//TRIM(DblToMedmChar(G%Clone(0)%ETotal))
    WRITE(*,*)TRIM(Mssg) 
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

       Mssg=ProcessName('NEBForce','Image '//TRIM(IntToChar(I))) &
            //' D = '//TRIM(FltToShrtChar(Dist))                 &
            //', E = '//TRIM(DblToMedmChar(G%Clone(I)%ETotal))   &
            //', F = '//TRIM(DblToMedmChar(FProj))
       WRITE(*,*)TRIM(Mssg) 
    ENDDO
    Rm=SQRT(SUM((G%Clone(G%Clones+1)%Carts%D-G%Clone(G%Clones)%Carts%D)**2))
    Dist=Dist+Rm
    Mssg=ProcessName('NEBForce','Product')   &
         //' D = '//TRIM(FltToShrtChar(Dist)) &
         //', E = '//TRIM(DblToMedmChar(G%Clone(G%Clones+1)%ETotal)) 
    WRITE(*,*)TRIM(Mssg) 
  END SUBROUTINE NEBForce

  !===============================================================================
  ! Generate a cubic spline along the band and interpolate to find extrema
  !===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
