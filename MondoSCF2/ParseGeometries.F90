MODULE ParseGeometries
  USE Parse
  USE InOut
  USE OptionKeys
  USE PrettyPrint
  USE DynamicsKeys
  USE GeometryKeys
  USE PrettyPrint
  USE ControlStructures
  USE ls_rmsd
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
    INTEGER          :: I,iCLONE,HDFFileID,J,IGeo
    REAL(DOUBLE),DIMENSION(3,3) :: U
    REAL(DOUBLE),DIMENSION(3)   :: Center1,Center2
    REAL(DOUBLE) :: Error,R2
    !--------------------------------------------------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)
    !   NEB
    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
       ! Parse for the number of clones to use in NEB
       IF(.NOT.OptIntQ(Inp,CLONES,G%Clones))THEN
          CALL MondoHalt(PRSE_ERROR,'Doing NEB, did not find input number of beads (clones) to use ')
       ENDIF
       IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS)THEN
          ! Overide any potential change in the number of clones ...
          HDFFileID=OpenHDF(N%RFile)
          HDF_CurrentID=HDFFileID
          CALL Get(G%Clones,'clones')
          CALL CloseHDF(HDFFileID)          
       ENDIF
       ! Allocate Clones+2 geometries for NEB
       ALLOCATE(G%Clone(0:G%Clones+1))
       IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS)THEN
          ! Get all the images, including endpoints from HDF
          HDFFileID=OpenHDF(N%RFile)
          HDF_CurrentID=HDFFileID
          ! Get and print the begining state
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #0")
          CALL Get(G%Clone(0))
          CALL CloseHDFGroup(HDF_CurrentID)
          ! Reset the atomic number in the case that we had/have atomic numbers
          ! cooresponding to an effective core potential basis set
          CALL ReSetAtNum(G%Clone(0))
          CALL PPrint(G%Clone(0),FileName_O=N%GFile,Unit_O=Geo, &
               PrintGeom_O=O%GeomPrint,Clone_O=0)
          ! Get and print the end state
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(G%Clones+1)))
          CALL Get(G%Clone(G%Clones+1))
          CALL CloseHDFGroup(HDF_CurrentID)
          ! Reset the atomic number in the case that we had/have atomic numbers
          ! cooresponding to an effective core potential basis set
          CALL ReSetAtNum(G%Clone(G%Clones+1))
          CALL PPrint(G%Clone(G%Clones+1),FileName_O=N%GFile,Unit_O=Geo, &
               PrintGeom_O=O%GeomPrint,Clone_O=G%Clones+1)
          ! Get and print the midpoints, past and present
          CALL New(CurrentState,3)
          CurrentState%I=O%RestartState%I
          DO IGeo=MAX(O%RestartState%I(3)-101,1),O%RestartState%I(3)
             CurrentState%I(3)=IGeo
             DO iCLONE=1,G%Clones
                HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
                CALL Get(G%Clone(iCLONE),TAG_O=GTag(CurrentState))
                ! Reset the atomic number in the case that we had/have atomic numbers
                ! cooresponding to an effective core potential basis set
                CALL ReSetAtNum(G%Clone(iCLONE))
                IF(iGEO/=O%RestartState%I(3))THEN
                   CALL PPrint(G%Clone(iCLONE),FileName_O=N%GFile,Unit_O=Geo, &
                        PrintGeom_O=O%GeomPrint,Clone_O=iCLONE)
                ENDIF
                CALL CloseHDFGroup(HDF_CurrentID)
             ENDDO
          ENDDO
          CALL Delete(CurrentState)
          CALL CloseHDF(HDFFileID)          
       ELSE ! Get some endpoints 
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
             ! Reset the atomic number in the case that we had/have atomic numbers
             ! cooresponding to an effective core potential basis set
             CALL ReSetAtNum(G%Clone(0))
             CALL ReSetAtNum(G%Clone(G%Clones+1))
          ELSE 
             ! Read in the reactants geometry from input
             CALL ParseCoordinates(REACTANTS_BEGIN,REACTANTS_END,G%Clone(0),O%Coordinates)          
             ! Read in the products geometry from input
             CALL ParseCoordinates(PRODUCTS_BEGIN,PRODUCTS_END,G%Clone(G%Clones+1),O%Coordinates)  
          ENDIF
          CALL PPrint(G%Clone(0),FileName_O=N%GFile,Unit_O=Geo, &
               PrintGeom_O=O%GeomPrint,Clone_O=0,CrdInAng_O=.TRUE.)
          CALL PPrint(G%Clone(G%Clones+1),FileName_O=N%GFile,Unit_O=Geo, &
               PrintGeom_O=O%GeomPrint,Clone_O=G%Clones+1,CrdInAng_O=.TRUE.)
          ! Purify R and P images ...
          CALL NEBPurify(G,Init_O=.TRUE.)
          ! ... then interpolate ...
          DO iCLONE=1,G%Clones
             G%Clone(iCLONE)%NAtms=G%Clone(0)%NAtms
             G%Clone(iCLONE)%NKind=G%Clone(0)%NKind
             CALL New(G%Clone(iCLONE))
          ENDDO
          CALL NEBInit(G)
          ! .. and purify the NEB images ...
          CALL NEBPurify(G)
       ENDIF
       !  Parrelel Rep
    ELSEIF(O%Grad==GRAD_DO_DYNAMICS.AND.D%Parallel_Rep)THEN
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
       IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS)THEN
          HDFFileID=OpenHDF(N%RFile)
          HDF_CurrentID=HDFFileID
          CALL Get(G%Clones,'clones')
          ALLOCATE(G%Clone(1:G%Clones))          
          CALL New(CurrentState,3)
          CurrentState%I=O%RestartState%I
          DO IGeo=MAX(O%RestartState%I(3)-101,1),O%RestartState%I(3)
             CurrentState%I(3)=IGeo
             DO iCLONE=1,G%Clones
                HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
                CALL Get(G%Clone(iCLONE),TAG_O=GTag(CurrentState))
                ! Reset the atomic number in the case that we had/have atomic numbers
                ! cooresponding to an effective core potential basis set
                CALL ReSetAtNum(G%Clone(iCLONE))
                IF(iGEO/=O%RestartState%I(3))THEN
                   CALL PPrint(G%Clone(iCLONE),FileName_O=N%GFile,Unit_O=Geo, &
                        PrintGeom_O=O%GeomPrint,Clone_O=iCLONE)
                   ! CALL PPrint(G%Clone(iCLONE),Unit_O=6,PrintGeom_O=O%GeomPrint,Clone_O=iCLONE)
                ENDIF
                CALL CloseHDFGroup(HDF_CurrentID)
             ENDDO
          ENDDO
          CALL Delete(CurrentState)
          CALL CloseHDF(HDFFileID)          
       ELSE       
          G%Clones=1
          ALLOCATE(G%Clone(1))
          CALL ParseCoordinates(GEOMETRY_BEGIN,GEOMETRY_END,G%Clone(1),O%Coordinates)
          !          CALL PPrint(G%Clone(iCLONE),FileName_O=N%GFile,Unit_O=Geo,PrintGeom_O=O%GeomPrint)
       ENDIF
    ENDIF
    CLOSE(UNIT=Inp,STATUS='KEEP')
  END SUBROUTINE LoadCoordinates
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  SUBROUTINE ParseCoordinates(BeginDelimiter,EndDelimiter,G,Coordinates)
    CHARACTER(LEN=*)          :: BeginDelimiter,EndDelimiter
    TYPE(CRDS)                :: G
    TYPE(CHR_VECT)            :: C
    REAL(DOUBLE),DIMENSION(3) :: Carts(3)
    INTEGER                   :: J,N,Coordinates,L
    CHARACTER(LEN=3)          :: At,AtTmp
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
       ! Set the atom for freq calculation
       G%DoFreq%I(N)=0
       IF(SIZE(C%C)==4) THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%Carts%D(1,N)=CharToDbl(C%C(2))
          G%Carts%D(2,N)=CharToDbl(C%C(3))
          G%Carts%D(3,N)=CharToDbl(C%C(4))
          G%CConstrain%I(N)=0
          G%Velocity%D(1:3,N)=0.0D0
       ELSEIF(SIZE(C%C)==5)THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%Carts%D(1,N)=CharToDbl(C%C(2))
          G%Carts%D(2,N)=CharToDbl(C%C(3))
          G%Carts%D(3,N)=CharToDbl(C%C(4))
          IF(TRIM(C%C(5))=='u')THEN
             G%CConstrain%I(N)=0
          ELSE IF(TRIM(C%C(5))=='c')THEN
             G%CConstrain%I(N)=1
          ELSE IF(TRIM(C%C(5))=='r')THEN
             G%CConstrain%I(N)=2
          ELSE IF(TRIM(C%C(5))=='f')THEN
             ! We found an atom for freq calculation
             G%DoFreq%I(N)=1
             ! We set the constrain for any case
             G%CConstrain%I(N)=0
          ELSE
             G%CConstrain%I(N)=0
          ENDIF
          G%Velocity%D(1:3,N)=0.0D0
       ELSEIF(SIZE(C%C)==7)THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%Carts%D(1,N) =CharToDbl(C%C(2))
          G%Carts%D(2,N) =CharToDbl(C%C(3))
          G%Carts%D(3,N) =CharToDbl(C%C(4))
          G%Velocity%D(1,N)=CharToDbl(C%C(5))
          G%Velocity%D(2,N)=CharToDbl(C%C(6))
          G%Velocity%D(3,N)=CharToDbl(C%C(7))
          G%CConstrain%I(N)=0
       ELSEIF(SIZE(C%C)==8)THEN
          At=TRIM(ADJUSTL(C%C(1)))
          G%Carts%D(1,N) =CharToDbl(C%C(2))
          G%Carts%D(2,N) =CharToDbl(C%C(3))
          G%Carts%D(3,N) =CharToDbl(C%C(4))
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
!VW    Check for multiple definition of the same atom.
       L=SCAN(At,Numbers)
       SELECT CASE(L)
       CASE(0);AtTmp=At
       CASE(2);AtTmp(1:2)=At(1:1)//' '! we need this blank to be compatible with the definition of Ats.
       CASE(3);AtTmp(1:2)=At(1:2)
       CASE DEFAULT;CALL MondoHalt(PRSE_ERROR,'Cannot regonize this atom At=<'//At//'>')
       END SELECT
!      Find the atom number (elements >= 105 are ghost functions) 
       DO J=1,107
          !vwIF(At==Ats(J))THEN
          IF(AtTmp(1:2)==Ats(J)(1:2))THEN
             G%AtNum%D(N)=J
             G%AtNam%C(N)=At!vwAts(J)
             G%AtMMTyp%C(N)='UNK' 
             G%AtMss%D(N)=AtsMss(J)
             EXIT
          ENDIF
       ENDDO
    ENDDO
    !
    ! If no freqs have been explicitly set, then do all of them.
    IF(SUM(G%DoFreq%I).EQ.0)G%DoFreq%I=1
    write(*,*) 'ParseGeometies: G%DoFreq%I=',G%DoFreq%I
    !
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
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  SUBROUTINE SpinCoords(G) 
    TYPE(CRDS)     :: G
    INTEGER        :: I,J,NUnPEl
    ! Calculate the electronic coordinates
    G%NElec=0
    DO I=1,G%NAtms
       IF(G%AtNum%D(i) < 105.D0) THEN
          G%NElec=G%NElec+INT(G%AtNum%D(I))
       ENDIF
    ENDDO
    G%NElec=G%NElec-G%TotCh
    NUnPEl=G%Multp-1
!vw    IF(NUnPEl.NE.0) &
!vw         CALL MondoHalt(PRSE_ERROR,'Open shell not supported yet.'   &
!vw         //' NElectrons = '//TRIM(IntToChar(G%NElec)) &
!vw         //' NUnPairedE = '//TRIM(IntToChar(NUnPEl)))
    G%NAlph=(G%NElec+NUnPEl)/2
    G%NBeta=(G%NElec-NUnPEl)/2
    IF(G%NAlph+G%NBeta.NE.G%NElec)CALL Halt('SpinCoords: Did you give the right charge/multiplicity!')
  END SUBROUTINE SpinCoords
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  SUBROUTINE ReSetAtNum(G)
    TYPE(CRDS) :: G
    INTEGER :: I,J,IAtomNo
    DO I=1,G%NAtms
       IAtomNo=G%AtNum%D(I)
       IF(TRIM(Ats(I))/=TRIM(G%AtNam%C(I)))THEN
          DO J=1,105
             IF(TRIM(Ats(J))==TRIM(G%AtNam%C(I)))THEN
                G%AtNum%D(I)=J
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE ReSetAtNum
!------------------------------------------------------------------------!
!
!------------------------------------------------------------------------!
  FUNCTION GTag(State) RESULT(Tag)
    TYPE(INT_VECT)   :: State
    CHARACTER(LEN=4) :: Tag
    Tag=IntToChar(State%I(3))
  END FUNCTION GTag

END MODULE ParseGeometries
