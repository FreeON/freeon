MODULE ParseExtraCoords
   USE ControlStructures
   USE DerivedTypes
   USE GlobalScalars
!  USE Macros
   USE OptionKeys
   USE GeometryKeys
   USE ParseGeometries
   USE Massage
   USE InCoords
   USE SetXYZ
   USE MemMan 
   IMPLICIT NONE
   CONTAINS
   !
   SUBROUTINE LoadExtraCoords(GOpt,Opts,Nams,Geos)
     !
     ! This subroutine parses the inPut file for
     ! additional internal coordinate definitions
     ! and constraints.
     ! WARNING! Constraint values must be given in 
     ! Angstroems and Degrees
     !
     TYPE(FileNames)             :: Nams
     TYPE(Options)               :: Opts
     TYPE(Geometries)            :: Geos
     TYPE(GeomOpt)               :: GOpt
     CHARACTER(LEN=DCL)          :: Line,LineLowCase,Atomname,chGEO
     CHARACTER(LEN=5)            :: CHAR
     INTEGER                     :: I1,I2,J,NIntCs,SerNum,NConstr
     INTEGER                     :: NatmsLoc
     INTEGER                     :: NCartConstr,iCLONE
     REAL(DOUBLE)                :: V,Value,DegToRad
     TYPE(DBL_RNK2)              :: XYZ 
     TYPE(INT_VECT)              :: CConstrain
     TYPE(CRDS)                  :: GMLoc
     INTEGER                     :: HDFFileID,LastGeom
     !
     CALL OpenASCII(Nams%IFile,Inp)
     !
     NatmsLoc=Geos%Clone(1)%Natms
     CALL New(XYZ,(/3,NatmsLoc/))
     CALL New(CConstrain,NatmsLoc)
     !
     IF(Opts%Guess==GUESS_EQ_RESTART) THEN
       HDFFileID=OpenHDF(Nams%RFile)
       HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                     "Clone #"//TRIM(IntToChar(1)))
         LastGeom=Opts%RestartState%I(3)
         chGEO=IntToChar(LastGeom)
         CALL Get(GMLoc,chGEO)
       CALL CloseHDFGroup(HDF_CurrentID)
       CALL CloseHDF(HDFFileID)
     ! CALL ToAtomicUnits(GMLoc)
       XYZ%D=GMLoc%AbCarts%D
       CConstrain%I=GMLoc%CConstrain%I
       CALL Delete(GMLoc)
     ELSE
       XYZ%D=Geos%Clone(1)%AbCarts%D
       CConstrain%I=Geos%Clone(1)%CConstrain%I
     ENDIF
     !
     DegToRad=PI/180.D0
     !
     ! Find extra internal coordinates and constraints
     !
     NIntCs=0
     NConstr=0
     NCartConstr=0
     !
     IF(FindMixedCaseKey('BEGIN_ADD_INTERNALS',Inp)) THEN
       !
       CALL AlignLowCase('begin_add_internals',Inp)
       DO 
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
         LineLowCase = Line
         Call LowCase(LineLowCase)
       !
       IF(INDEX(LineLowCase,'stre')/=0.OR.&
          INDEX(LineLowCase,'bend')/=0.OR.&
          INDEX(LineLowCase,'alpha')/=0.OR.&
          INDEX(LineLowCase,'beta')/=0.OR.&
          INDEX(LineLowCase,'gamma')/=0.OR.&
          INDEX(LineLowCase,'tors')/=0.OR.&
          INDEX(LineLowCase,'outp')/=0.OR.&
          INDEX(LineLowCase,'linb1')/=0.OR.&
          INDEX(LineLowCase,'linb2')/=0) THEN
          NIntCs=NIntCs+1 
       ELSE IF(INDEX(LineLowCase,'cart')/=0) THEN
          NIntCs=NIntCs+3 
       ENDIF
       !
         IF(INDEX(LineLowCase,'end_add_internals')/=0) GO TO 2
       ENDDO
       1  CALL MondoHalt(PRSE_ERROR, &
          ' Found no <end_add_internals> in inPut file '//TRIM(InpFile))
       2  CONTINUE
       !
       ! Generate an addition to the IntC set!
       !
       CALL New(GOpt%ExtIntCs,NIntCs)
       GOpt%ExtIntCs%Active%L=.TRUE.
       !
       ! Parse again and fill GOpt%ExtIntCs!
       !
       NIntCs=0
       CALL AlignLowCase('begin_add_internals',Inp)
       DO 
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
         LineLowCase = Line
         Call LowCase(LineLowCase)
         !
         IF(INDEX(LineLowCase,'stre_a')/=0) THEN
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:8)='STRE_A  ' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:6)=(/0,0,0,1,0,0/)
                GOpt%ExtIntCs%Active%L(NIntCs)=.FALSE.
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAu
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'stre_b')/=0) THEN
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen==1) CALL Halt('Extra coord stre_b while PBC dimension is 1')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:8)='STRE_B  ' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:6)=(/0,0,0,0,1,0/)
                GOpt%ExtIntCs%Active%L(NIntCs)=.FALSE.
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAu
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'stre_c')/=0) THEN
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen<3) CALL Halt('Extra coord stre_c while PBC dimension is < 3')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:8)='STRE_C  ' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:6)=(/0,0,0,0,0,1/)
                GOpt%ExtIntCs%Active%L(NIntCs)=.FALSE.
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAu
                  NConstr=NConstr+1
                ENDIF
         ELSE IF(INDEX(LineLowCase,'stre')/=0) THEN
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:5)='STRE ' 
         !--------------------
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:2),Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAu
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'bend')/=0) THEN 
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:5)='BEND ' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3),Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'gamma')/=0) THEN 
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen<2) CALL Halt('Extra coord gamma while PBC dimension is 1')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:8)='GAMMA   ' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/1,0,0,0,0,0,0,1,0/)
                GOpt%ExtIntCs%Active%L(NIntCs)=.FALSE.
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'area_l')/=0) THEN 
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen/=2) CALL Halt('Extra coord area_l while PBC dimension is 1')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:6)='AREA_L' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/0,0,0,1,0,0,0,1,0/)
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)= &
                                             Value*AngstromsToAu**2
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'volm_l')/=0) THEN 
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen/=3) CALL Halt('Extra coord volm_l while PBC dimension is /= 3')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:6)='VOLM_L' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:12)=(/0,0,0,1,0,0,0,1,0,0,0,1/)
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAu**3
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'beta')/=0) THEN 
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen/=3) CALL Halt('Extra coord beta while PBC dimension is /= 3')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:8)='BETA    ' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/1,0,0,0,0,0,0,0,1/)
                GOpt%ExtIntCs%Active%L(NIntCs)=.FALSE.
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'alpha')/=0) THEN 
         !--------------------
                IF(Geos%Clone(1)%PBC%Dimen/=3) CALL Halt('Extra coord alpha while PBC dimension is /= 3')
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:8)='ALPHA   ' 
                GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
                GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/0,1,0,0,0,0,0,0,1/)
                GOpt%ExtIntCs%Active%L(NIntCs)=.FALSE.
                IF(INDEX(LineLowCase,'.')==0) THEN
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'tors')/=0) THEN 
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:5)='TORS ' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4),Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'outp')/=0) THEN 
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:5)='OUTP ' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4),Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'linb1')/=0) THEN 
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:5)='LINB1' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3),Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'linb2')/=0) THEN 
         !--------------------
                NIntCs=NIntCs+1 
                GOpt%ExtIntCs%Def%C(NIntCs)(1:5)='LINB2'
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3),Value
                  GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
                  GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'cart')/=0) THEN
         !--------------------
                GOpt%ExtIntCs%Def%C(NIntCs+1)(1:5)='CARTX' 
                GOpt%ExtIntCs%Def%C(NIntCs+2)(1:5)='CARTY' 
                GOpt%ExtIntCs%Def%C(NIntCs+3)(1:5)='CARTZ' 
         !
                READ(LineLowCase,*) CHAR,SerNum  
         !
                GOpt%ExtIntCs%Atoms%I(NIntCs+1,1)=SerNum
                GOpt%ExtIntCs%Atoms%I(NIntCs+2,1)=SerNum
                GOpt%ExtIntCs%Atoms%I(NIntCs+3,1)=SerNum
         !
                GOpt%ExtIntCs%Constraint%L(NIntCs+1)=.TRUE.
                GOpt%ExtIntCs%Constraint%L(NIntCs+2)=.TRUE.
                GOpt%ExtIntCs%Constraint%L(NIntCs+3)=.TRUE.
                  NConstr=NConstr+3
         !!!! supposing that constraints are the same for all clones
                GOpt%ExtIntCs%ConstrValue%D(NIntCs+1)=XYZ%D(1,SerNum)
                GOpt%ExtIntCs%ConstrValue%D(NIntCs+2)=XYZ%D(2,SerNum)
                GOpt%ExtIntCs%ConstrValue%D(NIntCs+3)=XYZ%D(3,SerNum)
         !
                NIntCs=NIntCs+3 
                NCartConstr=NCartConstr+3
         !
         !               IF(INDEX(LineLowCase,'.')==0) THEN
         ! READ(LineLowCase,*) CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1)
         !               ELSE
         !               ENDIF
         !--------------------
         ELSE 
         ! dont do anything with non-conform lines
         !
         ENDIF
         !          
         IF(INDEX(LineLowCase,'end_add_internals')/=0) GO TO 12
       ENDDO
       11  CALL MondoHalt(PRSE_ERROR, &
           ' Found no <end_add_internals> in inPut file '&
           //TRIM(InpFile))
       12  CONTINUE
       !
     ELSE
       NIntCs=0
       CALL New(GOpt%ExtIntCs,NIntCs)
     ENDIF !!! key for extra internals found
     ! 
     ! Fill in Cartesian constraints stored 
     ! in Geos%Clone(1)%CConstrain%I
     ! Again, supposing that constraints are the same for all clones.
     ! Then, renumber atoms in IntCs by exclusion of Rigid atoms
     !
     CALL MergeConstr(GOpt%ExtIntCs,XYZ%D,CConstrain%I, &
                      NIntCs,NConstr,NCartConstr)
     CALL ReNumbIntC(GOpt%ExtIntCs,CConstrain%I) 
     CALL ReOrdIntC(GOpt%ExtIntCs,NIntCs)
     !
     GOpt%ExtIntCs%N=NIntCs
     GOpt%CoordCtrl%NExtra=NIntCs
     GOpt%Constr%NConstr=NConstr
     GOpt%Constr%NCartConstr=NCartConstr
     !
     CLOSE(Inp,STATUS='KEEP')
     CALL Delete(XYZ)
     CALL Delete(CConstrain)
    !CALL PrtIntCoords(GOPt%ExtIntCs,GOpt%ExtIntCs%Value%D,&
    !                  'chk extra ')
   END SUBROUTINE LoadExtraCoords
!
!------------------------------------------------------------------
!
   SUBROUTINE ReOrdIntC(IntCs,NIntC)
     TYPE(INTC)           :: IntCs
     INTEGER              :: NIntC,I,J
     INTEGER,DIMENSION(4) :: Atoms
     !
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         Atoms(1:2)=IntCs%Atoms%I(I,1:2)
         IF(Atoms(1)>Atoms(2)) THEN
           DO J=1,2 ; IntCs%Atoms%I(I,J)=Atoms(3-J) ; ENDDO
         ENDIF
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND'.OR. &
               IntCs%Def%C(I)(1:4)=='LINB') THEN
         Atoms(1:3)=IntCs%Atoms%I(I,1:3)
         IF(Atoms(1)>Atoms(3)) THEN
           DO J=1,3 ; IntCs%Atoms%I(I,J)=Atoms(4-J) ; ENDDO
         ENDIF
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
          Atoms(1:4)=IntCs%Atoms%I(I,1:4)
          IF(Atoms(2)>Atoms(3)) THEN
            DO J=1,4 ; IntCs%Atoms%I(I,J)=Atoms(5-J) ; ENDDO
          ENDIF
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
          Atoms(1:4)=IntCs%Atoms%I(I,1:4)
          IF(Atoms(3)>Atoms(4)) THEN
            IntCs%Atoms%I(I,3)=Atoms(4)
            IntCs%Atoms%I(I,4)=Atoms(3)
          ENDIF
       ENDIF
     ENDDO 
   END SUBROUTINE ReOrdIntC
!
!------------------------------------------------------------------
!
   SUBROUTINE MergeConstr(IntCs,XYZ,CConstrain,&
                          NIntCs,NConstr,NCartConstr)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER,DIMENSION(:)        :: CConstrain
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER                     :: NIntCs,NConstr,NCartConstr
     INTEGER                     :: I,NNewC,NatmsLoc
     !
     NatmsLoc=SIZE(XYZ,2)
     NNewC=0
     DO I=1,NatmsLoc    
       IF(CConstrain(I)==1) NNewC=NNewC+3
     ENDDO
     IF(NNewC==0) RETURN
     NConstr=NConstr+NNewC
     NCartConstr=NCartConstr+NNewC
     ! fill in old intcs
     CALL New(IntC_New,NIntCs+NNewC)
     CALL Set_INTC_EQ_INTC(IntCs,IntC_New,1,NIntCs,1)
     ! fill in intcs related to Cartesian constraints
     NNewC=NIntCs
     DO I=1,NatmsLoc    
       IF(CConstrain(I)==1) THEN
         IntC_New%Def%C(NNewC+1)(1:5)='CARTX' 
         IntC_New%Def%C(NNewC+2)(1:5)='CARTY' 
         IntC_New%Def%C(NNewC+3)(1:5)='CARTZ' 
         !
         IntC_New%Atoms%I(NNewC+1:NNewC+3,1:4)=0
         IntC_New%Atoms%I(NNewC+1,1)=I
         IntC_New%Atoms%I(NNewC+2,1)=I
         IntC_New%Atoms%I(NNewC+3,1)=I
         !
         IntC_New%Constraint%L(NNewC+1)=.TRUE.
         IntC_New%Constraint%L(NNewC+2)=.TRUE.
         IntC_New%Constraint%L(NNewC+3)=.TRUE.
         !
         IntC_New%ConstrValue%D(NNewC+1)=XYZ(1,I)
         IntC_New%ConstrValue%D(NNewC+2)=XYZ(2,I)
         IntC_New%ConstrValue%D(NNewC+3)=XYZ(3,I)
         NNewC=NNewC+3
       ENDIF
     ENDDO
     !
     CALL Delete(IntCs)
     NIntCs=NNewC
     CALL New(IntCs,NIntCs)
     CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntCs,1)
     CALL Delete(IntC_New)
   END SUBROUTINE MergeConstr
!
!----------------------------------------------------------------
!
   SUBROUTINE ReNumbIntC(IntCs,CConstr) 
     TYPE(INTC)           :: IntCs
     INTEGER,DIMENSION(:) :: CConstr
     INTEGER              :: I,J,K,NatmsOld,NatmsNew
     TYPE(INT_VECT)       :: Map
     !
     NatmsOld=SIZE(CConstr)
     NatmsNew=0
     CALL New(Map,NatmsOld)
     DO I=1,NatmsOld
       IF(CConstr(I)/=2) THEN
         NatmsNew=NatmsNew+1
         Map%I(I)=NatmsNew
       ELSE
         Map%I(I)=0
       ENDIF
     ENDDO
     !
     DO I=1,IntCs%N
       DO J=1,4
         K=IntCs%Atoms%I(I,J)
         IF(K==0) EXIT
         IntCs%Atoms%I(I,J)=Map%I(K)
       ENDDO
     ENDDO
     !
     CALL Delete(Map)
   END SUBROUTINE ReNumbIntC
END MODULE ParseExtraCoords
