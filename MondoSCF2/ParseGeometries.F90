MODULE ParseGeometries
  USE Parse
  USE InOut
  USE OptionKeys
  USE DynamicsKeys
  USE GeometryKeys
  USE ControlStructures
CONTAINS
  !
  ! <<README GRAME, HUGH AND KAROLY>> <<README GRAME, HUGH AND KAROLY>> <<README GRAME, HUGH AND KAROLY>> <<README GRAME, HUGH AND KAROLY>>
  !
  ! >>>>>>>>>>>>>> THIS IS WHERE ALL COORDINATES FOR TS SEARCH, DYNAMICS, AND QMMM SHOULD BE PARSED IN AND SORTED OUT <<<<<<<<<<<<<<<<<<<<
  !
  !                IF WE ARE RESTARTING, THEN GEOMETRIES SHOULD BE OBTAINED FROM O%RestartHDF, INCLUDING FOR DYNAMICS
  ! 
  !                THIS MEANS THAT VELOCITIES, FORCES, AND GEOMETRIES SHOULD BE ARCHIVED INTO THE HDF FILE. 
  !
  !                THIS SHOULD NOT BE UNDERTAKEN UNTILL WE HAVE HDF GROUPS ENABLED FOR EACH KLONE
  !
  ! <<README GRAME, HUGH AND KAROLY>> <<README GRAME, HUGH AND KAROLY>> <<README GRAME, HUGH AND KAROLY>> <<README GRAME, HUGH AND KAROLY>>
  !
  !================================================================================================================
  !
  !================================================================================================================
  SUBROUTINE LoadCoordinates(N,O,D,G)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Dynamics)   :: D
    TYPE(Geometries) :: G
    !--------------------------------------------------------------------------------------------------------------!
    IF(O%GradOpt==GRAD_TS_SEARCH_NEB)THEN
       CALL MondoHalt(PRSE_ERROR,' To here, G%Klones has not been parsed, needs some more work ')
       ALLOCATE(G%Klone(1:G%Klones))
       IF(O%Guess==GUESS_EQ_RESTART)THEN
       ELSE
          CALL ParseCoordinates(REACTANTS_BEGIN,REACTANTS_END,G%Klone(1))        ! Read in the reactants geometry
          CALL ParseCoordinates(PRODUCTS_BEGIN,PRODUCTS_END,G%Klone(G%Klones))   ! Read in the products geometry
       ENDIF
       !
       ! GRAEME PUTS A SUBROUTINE HERE TO ALLOCATE AND SET THE REMAINING (G%Klones-2) KLONES ?
       !
    ELSEIF(O%GradOpt==GRAD_DO_DYNAMICS.AND.D%MDAlgorithm==MD_PARALLEL_REP)THEN
#ifdef !defined(PARALLEL)
       CALL MondoHalt(PRSE_ERROR,' MondoSCF must be compiled in parallel for replica exchange to be active.')
#endif
       ALLOCATE(G%Klone(1:G%Klones))
       IF(O%Guess==GUESS_EQ_RESTART)THEN
       ELSE       
          ! HUGH PUTS A ROUTINE TO GENERATE (G%Klones) PERTURBED TRAJECTORIES AND VELOCITIES?
       ENDIF
    ELSE
       G%Klones=1
       ALLOCATE(G%Klone(1))
       IF(O%Guess==GUESS_EQ_RESTART)THEN
          HDFFileID=OpenHDFFile(N%RFile)
          CALL Get(G%Klone(1),TAG_O=GeoTag(O%RestartState))
          CALL CloseHDF(HDFFileID)
          HDFFileID=N%NewFileID          
       ELSE       
          CALL ParseCoordinates(GEOMETRY_BEGIN,GEOMETRY_END,G%Klone(1))
       ENDIF
    ENDIF

  END SUBROUTINE LoadCoordinates

  FUNCTION GeoTag(State) RESULT(Tag)
    TYPE(INT_VECT)   :: State
    CHARACTER(LEN=4) :: Tag
    Tag=IntToChar(State%I(3))
  END FUNCTION GeoTag

  SUBROUTINE ParseCoordinates(BeginDelimiter,EndDelimiter,G)
    CHARACTER(LEN=*)          :: BeginDelimiter,EndDelimiter
    TYPE(CRDS)                :: G
    REAL(DOUBLE),DIMENSION(3) :: Carts(3)
    INTEGER                   :: J,N
    CHARACTER(LEN=2)          :: At
    CHARACTER(LEN=DCL)        :: Line,LineLowCase
    !------------------------------------------------------------------------!
    ! Determine the number of atoms in this block
    N=0
    CALL AlignLowCase(BeginDelimiter,Inp)
    DO 
       READ(Inp,DEFAULT_CHR_FMT,END=1)Line
       LineLowCase = Line
       Call LowCase(LineLowCase)
       IF(INDEX(LineLowCase,EndDelimiter)/=0)EXIT
       N=N+1
    ENDDO
    G%NAtms=N
    CALL New(G) ! Allocate a geometry
    ! Parse the coordinates
    CALL AlignLowCase(BeginDelimiter,Inp)
    DO 
       READ(Inp,DEFAULT_CHR_FMT,END=1)Line
       LineLowCase = Line
       Call LowCase(LineLowCase)
       IF(INDEX(LineLowCase,EndDelimiter)/=0)EXIT
       N=N+1
       CALL LineToGeom(Line,At,Carts)
       G%Carts%D(:,N)=Carts(1:3) 
       !       G%Vects%D(:,N)=Carts(4:6)
       ! Find the atom number (element 105 is a ghost function) 
       DO J=1,105
          IF(At==Ats(J))THEN
             G%AtNum%D(N)=J
             G%AtMss%D(N)=AtsMss(J)
             EXIT
          ENDIF
       ENDDO
    ENDDO
    IF(N/=G%NAtms) &
         CALL MondoHalt(PRSE_ERROR,'Atom number mismatch in ParseCoordinates')
!
!    ULTIMATELY, THE FOLLOWING ITEMS SHOULD BE ASSOCIATED WITH THE Geometries TYPE
!    RATHER THAN THE CRDS TYPE
!
    ! Parse coordinates determining spin coordinates
    IF(.NOT.OptDblQ(Inp,TOTAL_CHARGE,G%TotCh))  &
         CALL MondoHalt(PRSE_ERROR,TOTAL_Charge//' not found in input.')
    IF(.NOT.OptIntQ(Inp,MULTIPLICITY,G%Multp))  &
         CALL MondoHalt(PRSE_ERROR,MULTIPLICITY//' not found in input.')
    CALL SpinCoords(G) ! Compute closed shell spin coordinates
    ! Parsing of misc geometry info
    G%InAU=OptKeyQ(Inp,GEOMETRY,IN_AU)
    IF(OptKeyQ(Inp,GEOMETRY,Z_ORDER))THEN
       G%Ordrd=SFC_PEANO
    ELSEIF(OptKeyQ(Inp,GEOMETRY,RANDOM_ORDER))THEN
       G%Ordrd=SFC_RANDOM
    ELSEIF(OptKeyQ(Inp,GEOMETRY,H_Order))THEN
       G%Ordrd=SFC_HILBERT 
    ELSEIF(OptKeyQ(Inp,GEOMETRY,TRAVEL_Order))THEN
       G%Ordrd=SFC_TRAVEL 
    ELSEIF(OptKeyQ(Inp,GEOMETRY,TABLETRAV_ORDER))THEN
       G%Ordrd=SFC_TABLETRAV
    ELSE
       G%Ordrd=SFC_NONE
    ENDIF
    RETURN
    ! End of file error message
1   CALL Halt('While parsing, failed to find '//EndDelimiter)
  END SUBROUTINE ParseCoordinates

  SUBROUTINE SpinCoords(G) 
    TYPE(CRDS)     :: G
    INTEGER        :: I,J,NUnPEl
    !----------------------------------------------------------------------------
    ! Calculate the electronic coordinates
    G%NElec=0
    DO I=1,G%NAtms
       G%NElec=G%NElec+G%AtNum%D(I)
    ENDDO
    G%NElec=G%NElec-G%TotCh
    NUnPEl=G%Multp-1
    IF(NUnPEl.NE.0) &
         CALL MondoHalt(PRSE_ERROR,'Open shell not supported yet.'   &
         //' NElectrons = '//TRIM(IntToChar(G%NElec)) &
         //' NUnPairedE = '//TRIM(IntToChar(NUnPEl)))
    G%NAlph=DBLE(G%NElec+NUnPEl)*Half
    G%NBeta=DBLE(G%NElec-NUnPEl)*Half
  END SUBROUTINE SpinCoords

!  << C.J. READMEM >> << C.J. READMEM >> << C.J. READMEM >> << C.J. READMEM >> << C.J. READMEM >>
!
!  >>>>>>> ALL RESCALING, TRANLATIONS, TRANSFORMATIONS ETC ETC SHOULD ALL BE DONE HERE <<<<<<<<<
!
!  << C.J. READMEM >> << C.J. READMEM >> << C.J. READMEM >> << C.J. READMEM >> << C.J. READMEM >>

  SUBROUTINE MassageCoordinates(G,P)
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
    INTEGER          :: I
    !-------------------------------------------------------------------------!
       DO I=1,G%Klones
#ifdef PERIODIC
          IF(P%InAtomCrd)THEN
             ! Convert coordinates and unit cell to AU
             IF(.NOT.G%Klone(I)%InAU) &
                G%Klone(I)%Carts%D=AngstromsToAU*G%Klone(I)%Carts%D           
             G%Klone(I)%AbCarts%D=G%Klone(I)%Carts%D
             ! Don't really compute fractional coordinates as they are never actually used!
!             CALL CalFracCarts( G%Klone(I) )
          ELSE
             CALL MondoHalt(PRSE_ERROR,' Need to get rid of box carts etc and just overwrite carts in AtomCarts ' &
                                       //' since they are only used in input' )
             ! CALL CalAtomCarts(GM)
          ENDIF
          IF(P%Trans_COM) &
             CALL CalTransVec(G%Klone(I))
          CALL Translate(G%Klone(I),G%Klone(I)%PBC%TransVec)
          CALL WrapAtoms(G%Klone(I))
#else
          ! Convert to AU
          IF(.NOT.G%Klone(I)%InAU) &
              G%Klone(I)%Carts%D=AngstromsToAU*G%Klone(I)%Carts%D           
          ! Find the box bounding the atomic positions after conversion to AU
          G%Klone(I)%BndBox%D(:,1)=+1.D8
          G%Klone(I)%BndBox%D(:,2)=-1.D8
          DO J=1,G%Klone(I)%NAtms
             G%Klone(I)%BndBox%D(:,1)=MIN(G%Klone(I)%BndBox%D(:,1),G%Klone(I)%Carts%D(:,J) )
             G%Klone(I)%BndBox%D(:,2)=MAX(G%Klone(I)%BndBox%D(:,2),G%Klone(I)%Carts%D(:,J))
          ENDDO
       ENDIF
#endif
    ENDDO
  END SUBROUTINE MassageCoordinates

END MODULE ParseGeometries
