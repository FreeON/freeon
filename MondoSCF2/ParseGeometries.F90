MODULE ParseGeometries
  USE Parse
  USE InOut
  USE OptionKeys
  USE PrettyPrint
  USE DynamicsKeys
  USE GeometryKeys
  USE PrettyPrint
  USE ControlStructures
  USE NEB
  USE Conflicted
CONTAINS
  !================================================================================================================
  !
  !================================================================================================================
  SUBROUTINE LoadCoordinates(N,O,D,G)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Dynamics)   :: D
    TYPE(Geometries) :: G
    TYPE(Int_Vect)   :: CurrentState
    INTEGER          :: I,iCLONE,HDFFileID
    !--------------------------------------------------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)
!   NEB
    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
       ! Parse for the number of clones to use in NEB
       IF(.NOT.OptIntQ(Inp,CLONES,G%Clones))THEN
          CALL MondoHalt(PRSE_ERROR,'Doing NEB, did not find input number of beads (clones) to use ')
       ENDIF
       IF(O%Guess==GUESS_EQ_RESTART)THEN
          ! Overide any potential change in the number of clones ...
          HDFFileID=OpenHDF(N%RFile)
          HDF_CurrentID=HDFFileID
          CALL Get(G%Clones,'clones')
          CALL CloseHDF(HDFFileID)          
       ENDIF
       ! Allocate Clones+2 geometries for NEB
       ALLOCATE(G%Clone(0:G%Clones+1))
       ! Get left and right endpoints
       IF(O%EndPts==ENDPOINTS_FROM_HDF)THEN            
          CALL New(CurrentState,3)
          ! Get the left endpoint from reactants HDF
          HDFFileID=OpenHDF(N%ReactantsFile)
          HDF_CurrentID=HDFFileID
          ! Get the current reactants state 
          CALL Get(CurrentState,'current_state')
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
          CALL Get(G%Clone(0),TAG_O=GTag(CurrentState))
          CALL CloseHDFGroup(HDF_CurrentID)
          CALL CloseHDF(HDFFileID)          
          ! Get the right endpoint from products HDF
          HDFFileID=OpenHDF(N%ProductsFile)
          HDF_CurrentID=HDFFileID
          ! Get the current reactants state 
          CALL Get(CurrentState,'current_state')
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(1)))
          CALL Get(G%Clone(G%Clones+1),TAG_O=GTag(CurrentState))
          CALL CloseHDFGroup(HDF_CurrentID)
          CALL CloseHDF(HDFFileID)          
          CALL Delete(CurrentState)
          ! CALL PPrint(G%Clone(0),Unit_O=6)
          ! CALL PPrint(G%Clone(G%Clones+1),Unit_O=6)
       ELSE 
          ! Read in the reactants geometry from input
          CALL ParseCoordinates(REACTANTS_BEGIN,REACTANTS_END,G%Clone(0),O%Coordinates)          
          ! Read in the products geometry from input
          CALL ParseCoordinates(PRODUCTS_BEGIN,PRODUCTS_END,G%Clone(G%Clones+1),O%Coordinates)   
       ENDIF
       IF(O%Guess==GUESS_EQ_RESTART)THEN
          ! Get midpoints from HDF
          HDFFileID=OpenHDF(N%RFile)
          HDF_CurrentID=HDFFileID
          DO iCLONE=1,G%Clones
             HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
             CALL Get(G%Clone(iCLONE),TAG_O=GTag(O%RestartState))
             CALL CloseHDFGroup(HDF_CurrentID)
          ENDDO
          CALL CloseHDF(HDFFileID)          
          ! DO I=0,G%Clones+1
          !    CALL PPrint(G%Clone(I),Unit_O=6)
          ! ENDDO
       ELSE
          ! Initialize midpoints via interpolation
          DO iCLONE=1,G%Clones
             G%Clone(iCLONE)%NAtms=G%Clone(0)%NAtms
             G%Clone(iCLONE)%NKind=G%Clone(0)%NKind
             CALL New(G%Clone(iCLONE))
          ENDDO
          CALL NEBInit(G)
       ENDIF
!   Parrelel Rep
    ELSEIF(O%Grad==GRAD_DO_DYNAMICS .AND. D%Parallel_Rep)THEN
       CALL MondoHalt(PRSE_ERROR,'Parralel Rep not implimented')
       ! Parse for the number of clones to use in Parallel Rep
       IF(.NOT.OptIntQ(Inp,CLONES,G%Clones))THEN
          CALL MondoHalt(PRSE_ERROR,'Doing Parralel Rep, input number of beads (clones) to use ')
       ENDIF
       ALLOCATE(G%Clone(1:G%Clones))
       IF(O%Guess==GUESS_EQ_RESTART)THEN
       ELSE         
       ENDIF 
!   Not doing any fancy shmancy with clones ...
    ELSE 
       IF(O%Guess==GUESS_EQ_RESTART .AND. .NOT. OptKeyQ(Inp,RESTART_OPTION,RESTART_NEWGEOM))THEN      
          HDFFileID=OpenHDF(N%RFile)
          HDF_CurrentID=HDFFileID
          CALL Get(G%Clones,'clones')
          ALLOCATE(G%Clone(1:G%Clones))          
          DO iCLONE=1,G%Clones
             HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
             CALL Get(G%Clone(iCLONE),TAG_O=GTag(O%RestartState))
             CALL CloseHDFGroup(HDF_CurrentID)
          ENDDO
          CALL CloseHDF(HDFFileID)          
       ELSE       
          G%Clones=1
          ALLOCATE(G%Clone(1))
          CALL ParseCoordinates(GEOMETRY_BEGIN,GEOMETRY_END,G%Clone(1),O%Coordinates)
          IF(O%Guess==GUESS_EQ_RESTART .AND. OptKeyQ(Inp,RESTART_OPTION,RESTART_NEWGEOM)) THEN
!             WRITE(*,*)' REPARSING GEOMETRY ON RESTART'
          ENDIF
       ENDIF
    ENDIF
    CLOSE(UNIT=Inp,STATUS='KEEP')
!
  END SUBROUTINE LoadCoordinates
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  FUNCTION GTag(State) RESULT(Tag)
    TYPE(INT_VECT)   :: State
    CHARACTER(LEN=4) :: Tag
    Tag=IntToChar(State%I(3))
  END FUNCTION GTag
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  SUBROUTINE ParseCoordinates(BeginDelimiter,EndDelimiter,G,Coordinates)
    CHARACTER(LEN=*)          :: BeginDelimiter,EndDelimiter
    TYPE(CRDS)                :: G
    TYPE(CHR_VECT)            :: C
    REAL(DOUBLE),DIMENSION(3) :: Carts(3)
    INTEGER                   :: J,N,Coordinates
    CHARACTER(LEN=2)          :: At
    CHARACTER(LEN=DCL)        :: Line,LineLowCase
!------------------------------------------------------------------------!
!   Determine the number of atoms in this block
    N=0
    CALL Align(BeginDelimiter,Inp)
    DO 
       READ(Inp,DEFAULT_CHR_FMT,END=1)Line
       LineLowCase = Line
       IF(INDEX(LineLowCase,EndDelimiter)/=0)EXIT
       N=N+1
    ENDDO
    G%NAtms=N
    IF(Coordinates==GRAD_CART_OPT) THEN
      IntCMax=1
    ELSE
      IntCMax=30*G%NAtms
    ENDIF
    G%IntCs%N=IntCMax
    N=0
    CALL New(G) ! Allocate a geometry
!   Parse the coordinates
    CALL Align(BeginDelimiter,Inp)
    DO 
       READ(Inp,DEFAULT_CHR_FMT,END=1)Line
       LineLowCase = Line
       IF(INDEX(LineLowCase,EndDelimiter)/=0)EXIT
       Call LowCase(LineLowCase)
       N=N+1
       CALL LineToChars(LineLowCase,C)
       IF(SIZE(C%C)<4)CALL MondoHalt(PRSE_ERROR,' bad data on parsing goemetry at line = <<' &
                                      //TRIM(LineLowCase)//'>>')
       IF(SIZE(C%C)==4) THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%AbCarts%D(1,N)=CharToDbl(C%C(2))
          G%AbCarts%D(2,N)=CharToDbl(C%C(3))
          G%AbCarts%D(3,N)=CharToDbl(C%C(4))
          G%CConstrain%I(N)=0
       ELSEIF(SIZE(C%C)==5)THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%AbCarts%D(1,N)=CharToDbl(C%C(2))
          G%AbCarts%D(2,N)=CharToDbl(C%C(3))
          G%AbCarts%D(3,N)=CharToDbl(C%C(4))
          IF(TRIM(C%C(5))=='u')THEN
             G%CConstrain%I(N)=0
          ELSE IF(TRIM(C%C(5))=='c')THEN
             G%CConstrain%I(N)=1
          ELSE IF(TRIM(C%C(5))=='r')THEN
             G%CConstrain%I(N)=2
          ELSE
             G%CConstrain%I(N)=0
          ENDIF
       ELSEIF(SIZE(C%C)==7)THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%AbCarts%D(1,N) =CharToDbl(C%C(2))
          G%AbCarts%D(2,N) =CharToDbl(C%C(3))
          G%AbCarts%D(3,N) =CharToDbl(C%C(4))
          G%Velocity%D(1,N)=CharToDbl(C%C(5))
          G%Velocity%D(2,N)=CharToDbl(C%C(6))
          G%Velocity%D(3,N)=CharToDbl(C%C(7))
          G%CConstrain%I(N)=0
       ELSEIF(SIZE(C%C)==8)THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%AbCarts%D(1,N) =CharToDbl(C%C(2))
          G%AbCarts%D(2,N) =CharToDbl(C%C(3))
          G%AbCarts%D(3,N) =CharToDbl(C%C(4))
          G%Velocity%D(1,N)=CharToDbl(C%C(5))
          G%Velocity%D(2,N)=CharToDbl(C%C(6))
          G%Velocity%D(3,N)=CharToDbl(C%C(7))
          IF(TRIM(C%C(5))=='u')THEN
             G%CConstrain%I(N)=0
          ELSE IF(TRIM(C%C(5))=='c')THEN
             G%CConstrain%I(N)=1
          ELSE IF(TRIM(C%C(5))=='r')THEN
             G%CConstrain%I(N)=2
          ELSE
             G%CConstrain%I(N)=0
          ENDIF
       ENDIF
!
       CALL Delete(C)
!      Find the atom number (elements >= 105 are ghost functions) 
       DO J=1,107
          IF(At==Ats(J))THEN
             G%AtNum%D(N)=J
             G%AtNam%C(N)=Ats(J)
             G%AtMMTyp%C(N)='UNK' 
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
       IF(G%AtNum%D(i) < 105.D0) THEN
          G%NElec=G%NElec+G%AtNum%D(I)
       ENDIF
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

END MODULE ParseGeometries
