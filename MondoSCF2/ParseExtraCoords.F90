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
     TYPE(Geometries)            :: Geos
     TYPE(GeomOpt)               :: GOpt
     CHARACTER(LEN=DCL)          :: Line,LineLowCase,Atomname
     TYPE(INTC)                  :: IntC_Extra
     CHARACTER(LEN=5)            :: CHAR
     INTEGER                     :: I1,I2,J,NIntC_Extra,SerNum,NConstr
     INTEGER                     :: NatmsLoc
     INTEGER                     :: NCartConstr,iCLONE
     REAL(DOUBLE)                :: V,Value,DegToRad
     TYPE(DBL_RNK2)              :: XYZ 
     TYPE(INT_VECT)              :: CConstrain
     TYPE(CRDS)                  :: GMLoc
     TYPE(Options)               :: Opts
     !
     CALL OpenASCII(Nams%IFile,Inp)
     !
     NatmsLoc=Geos%Clone(1)%Natms
     CALL New(XYZ,(/3,NatmsLoc/))
     CALL New(CConstrain,NatmsLoc)
     !
     IF(Opts%Guess==GUESS_EQ_RESTART) THEN
       ! reparse for values of constraints which have not been stored
       ! in CRDS or which might have been modified
       CALL ParseCoordinates(GEOMETRY_BEGIN,GEOMETRY_END,GMLoc)
       CALL ToAtomicUnits(GMLoc)
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
     NIntC_Extra=0
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
          INDEX(LineLowCase,'tors')/=0.OR.&
          INDEX(LineLowCase,'outp')/=0.OR.&
          INDEX(LineLowCase,'linb1')/=0.OR.&
          INDEX(LineLowCase,'linb2')/=0) THEN
          NIntC_Extra=NIntC_Extra+1 
       ELSE IF(INDEX(LineLowCase,'cart')/=0) THEN
          NIntC_Extra=NIntC_Extra+3 
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
       CALL New(IntC_Extra,NIntC_Extra)
       IntC_Extra%Atoms=0
       IntC_Extra%Value=Zero
       IntC_Extra%Constraint=.FALSE.
       IntC_Extra%ConstrValue=Zero   
       !
       ! Parse again and fill IntC_Extra!
       !
       NIntC_Extra=0
       CALL AlignLowCase('begin_add_internals',Inp)
       DO 
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
         LineLowCase = Line
         Call LowCase(LineLowCase)
         !
         IF(INDEX(LineLowCase,'stre')/=0) THEN

                NIntC_Extra=NIntC_Extra+1 
                IntC_Extra%DEF(NIntC_Extra)='STRE ' 
         !--------------------
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:2)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:2),Value
                  IntC_Extra%Constraint(NIntC_Extra)=.TRUE.
                  IntC_Extra%ConstrValue(NIntC_Extra)=Value*AngstromsToAu
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'bend')/=0) THEN 
         !--------------------
                NIntC_Extra=NIntC_Extra+1 
                IntC_Extra%DEF(NIntC_Extra)='BEND ' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:3)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:3),Value
                  IntC_Extra%Constraint(NIntC_Extra)=.TRUE.
                  IntC_Extra%ConstrValue(NIntC_Extra)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'tors')/=0) THEN 
         !--------------------
                NIntC_Extra=NIntC_Extra+1 
                IntC_Extra%DEF(NIntC_Extra)='TORS ' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:4)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:4),Value
                  IntC_Extra%Constraint(NIntC_Extra)=.TRUE.
                  IntC_Extra%ConstrValue(NIntC_Extra)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'outp')/=0) THEN 
         !--------------------
                NIntC_Extra=NIntC_Extra+1 
                IntC_Extra%DEF(NIntC_Extra)='OUTP ' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:4)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:4),Value
                  IntC_Extra%Constraint(NIntC_Extra)=.TRUE.
                  IntC_Extra%ConstrValue(NIntC_Extra)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'linb1')/=0) THEN 
         !--------------------
                NIntC_Extra=NIntC_Extra+1 
                IntC_Extra%DEF(NIntC_Extra)='LINB1' 
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:3)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:3),Value
                  IntC_Extra%Constraint(NIntC_Extra)=.TRUE.
                  IntC_Extra%ConstrValue(NIntC_Extra)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'linb2')/=0) THEN 
         !--------------------
                NIntC_Extra=NIntC_Extra+1 
                IntC_Extra%DEF(NIntC_Extra)='LINB2'
                IF(INDEX(LineLowCase,'.')==0) THEN
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:3)
                ELSE
                  READ(LineLowCase,*) &
                  CHAR,IntC_Extra%ATOMS(NIntC_Extra,1:3),Value
                  IntC_Extra%Constraint(NIntC_Extra)=.TRUE.
                  IntC_Extra%ConstrValue(NIntC_Extra)=Value*DegToRad 
                  NConstr=NConstr+1
                ENDIF
         !--------------------
         ELSE IF(INDEX(LineLowCase,'cart')/=0) THEN
         !--------------------
                IntC_Extra%DEF(NIntC_Extra+1)='CARTX' 
                IntC_Extra%DEF(NIntC_Extra+2)='CARTY' 
                IntC_Extra%DEF(NIntC_Extra+3)='CARTZ' 
         !
                  READ(LineLowCase,*) CHAR,SerNum  
         !
                IntC_Extra%ATOMS(NIntC_Extra+1,1)=SerNum
                IntC_Extra%ATOMS(NIntC_Extra+2,1)=SerNum
                IntC_Extra%ATOMS(NIntC_Extra+3,1)=SerNum
         !
                IntC_Extra%Constraint(NIntC_Extra+1)=.TRUE.
                IntC_Extra%Constraint(NIntC_Extra+2)=.TRUE.
                IntC_Extra%Constraint(NIntC_Extra+3)=.TRUE.
                  NConstr=NConstr+3
         !!!! supposing that constraints are the same for all clones
         !      DO iCLONE=1,Geos%Clones 
                  iCLONE=1
                  IntC_Extra%ConstrValue(NIntC_Extra+1)=&
                       XYZ%D(1,SerNum)
                  IntC_Extra%ConstrValue(NIntC_Extra+2)=&
                       XYZ%D(2,SerNum)
                  IntC_Extra%ConstrValue(NIntC_Extra+3)=&
                       XYZ%D(3,SerNum)
         !      ENDDO
         !
                  NIntC_Extra=NIntC_Extra+3 
                  NCartConstr=NCartConstr+3
         !
         !               IF(INDEX(LineLowCase,'.')==0) THEN
         ! READ(LineLowCase,*) CHAR,IntC_Extra%ATOMS(NIntC_Extra,1)
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
       NIntC_Extra=0
       CALL New(IntC_Extra,NIntC_Extra)
     ENDIF !!! key for extra internals found
     ! 
     ! Fill in Cartesian constraints stored 
     ! in Geos%Clone(iCLONE)%CConstrain%I
     ! Again, supposing that constraints are the same for all clones.
     ! Then, renumber atoms in IntC_Extra by exclusion of Rigid atoms
     !
     !DO iCLONE=1,Geos%Clones
        iCLONE=1
        CALL MergeConstr(IntC_Extra,XYZ%D,CConstrain%I, &
                         NIntC_Extra,NConstr,NCartConstr)
        CALL ReNumbIntC(IntC_Extra,CConstrain%I) 
     !ENDDO
     !
     ! Redefine size of Lagrangian arrays and save IntC_Extra to disk 
     !
     DO iCLONE=1,Geos%Clones
       Geos%Clone(iCLONE)%NLagr=NConstr
       CALL Delete(Geos%Clone(iCLONE)%LagrMult)
       CALL Delete(Geos%Clone(iCLONE)%GradMult)
       CALL Delete(Geos%Clone(iCLONE)%LagrDispl)
       CALL New(Geos%Clone(iCLONE)%LagrMult,NConstr)
       CALL New(Geos%Clone(iCLONE)%GradMult,NConstr)
       CALL New(Geos%Clone(iCLONE)%LagrDispl,NConstr)
       CALL WriteIntCs(IntC_Extra,TRIM(Nams%M_SCRATCH)//&
         TRIM(Nams%SCF_NAME)//'.'//TRIM(IntToChar(iCLONE))//'IntC_Extra')
     ENDDO
     !
     GOpt%CoordCtrl%NExtra=NIntC_Extra
     GOpt%Constr%NConstr=NConstr
     GOpt%Constr%NCartConstr=NCartConstr
     !
     !IF(NIntC_Extra>0) THEN
     !  CALL OPENAscii(Nams%OFile,Out)
     !  CALL PrtIntCoords(IntC_Extra,IntC_Extra%Value, &
     !  'Extra Internals & Constraints')
     !  CLOSE(Out,STATUS='KEEP')
     !ENDIF
     !
     CALL Delete(IntC_Extra)
     !
     CLOSE(Inp,STATUS='KEEP')
     !
     CALL Delete(XYZ)
     CALL Delete(CConstrain)
   END SUBROUTINE LoadExtraCoords
!
!------------------------------------------------------------------
!
   SUBROUTINE MergeConstr(IntC_Extra,XYZ,CConstrain,&
                          NIntC_Extra,NConstr,NCartConstr)
     TYPE(INTC)                  :: IntC_Extra,IntC_New
     INTEGER,DIMENSION(:)        :: CConstrain
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER                     :: NIntC_Extra,NConstr,NCartConstr
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
     CALL New(IntC_New,NIntC_Extra+NNewC)
     CALL Set_INTC_EQ_INTC(IntC_Extra,IntC_New,1,NIntC_Extra,1)
     ! fill in intcs related to Cartesian constraints
     NNewC=NIntC_Extra
     DO I=1,NatmsLoc    
       IF(CConstrain(I)==1) THEN
         IntC_New%Def(NNewC+1)='CARTX' 
         IntC_New%Def(NNewC+2)='CARTY' 
         IntC_New%Def(NNewC+3)='CARTZ' 
         !
         IntC_New%Atoms(NNewC+1:NNewC+3,1:4)=0
         IntC_New%ATOMS(NNewC+1,1)=I
         IntC_New%ATOMS(NNewC+2,1)=I
         IntC_New%ATOMS(NNewC+3,1)=I
         !
         IntC_New%Constraint(NNewC+1)=.TRUE.
         IntC_New%Constraint(NNewC+2)=.TRUE.
         IntC_New%Constraint(NNewC+3)=.TRUE.
         !
         IntC_New%ConstrValue(NNewC+1)=XYZ(1,I)
         IntC_New%ConstrValue(NNewC+2)=XYZ(2,I)
         IntC_New%ConstrValue(NNewC+3)=XYZ(3,I)
         NNewC=NNewC+3
       ENDIF
     ENDDO
     !
     CALL Delete(IntC_Extra)
     NIntC_Extra=NNewC
     CALL New(IntC_Extra,NIntC_Extra)
     CALL Set_INTC_EQ_INTC(IntC_New,IntC_Extra,1,NIntC_Extra,1)
     CALL Delete(IntC_New)
   END SUBROUTINE MergeConstr
!
!----------------------------------------------------------------
!
   SUBROUTINE ReNumbIntC(IntCs,CConstr) 
     TYPE(INTC)           :: IntCs
     INTEGER,DIMENSION(:) :: CConstr
     INTEGER              :: I,J,K,NIntC,NatmsOld,NatmsNew
     TYPE(INT_VECT)       :: Map
     !
     NIntC=SIZE(IntCs%Def)
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
     DO I=1,NIntC
       DO J=1,4
         K=IntCs%Atoms(I,J)
         IF(K==0) EXIT
         IntCs%Atoms(I,J)=Map%I(K)
       ENDDO
     ENDDO
     !
     CALL Delete(Map)
   END SUBROUTINE ReNumbIntC
END MODULE ParseExtraCoords
