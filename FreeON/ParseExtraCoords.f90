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

MODULE ParseExtraCoords
  USE ControlStructures
  USE DerivedTypes
  USE GlobalScalars
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
    CHARACTER(LEN=4)            :: CharAux4
    INTEGER,DIMENSION(4)        :: AtomAux4
    INTEGER                     :: I1,I2,J,NIntCs,SerNum,NConstr
    INTEGER                     :: NatmsLoc
    INTEGER                     :: NCartConstr,iCLONE
    REAL(DOUBLE)                :: V,Value,DegToRad
    REAL(DOUBLE)                :: InvBoxSh(3,3),Vect1(3),Vect2(3)
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
      XYZ%D=GMLoc%Carts%D
      CConstrain%I=GMLoc%CConstrain%I
      CALL Delete(GMLoc)
    ELSE
      XYZ%D=Geos%Clone(1)%Carts%D
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
      CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "reading internals")
      CALL AlignLowCase('begin_add_internals',Inp)
      DO
        READ(Inp,DEFAULT_CHR_FMT,END=1)Line
        CALL RemoveComments(Line)
        LineLowCase = Line
        Call LowCase(LineLowCase)
        !
        IF(INDEX(LineLowCase,'stre')/=0.OR.&
           INDEX(LineLowCase,'bend')/=0.OR.&
           INDEX(LineLowCase,'alpha')/=0.OR.&
           INDEX(LineLowCase,'beta')/=0.OR.&
           INDEX(LineLowCase,'gamma')/=0.OR.&
           INDEX(LineLowCase,'volm_l')/=0.OR.&
           INDEX(LineLowCase,'area_l')/=0.OR.&
           INDEX(LineLowCase,'tors')/=0.OR.&
           INDEX(LineLowCase,'outp')/=0.OR.&
           INDEX(LineLowCase,'linb1')/=0.OR.&
           INDEX(LineLowCase,'linb2')/=0.OR.&
           INDEX(LineLowCase,'cartx')/=0.OR.&
           INDEX(LineLowCase,'carty')/=0.OR.&
           INDEX(LineLowCase,'cartz')/=0) THEN
          NIntCs=NIntCs+1
        ENDIF

        IF(INDEX(LineLowCase,'end_add_internals')/=0) GO TO 2
      ENDDO
1     CALL MondoHalt(PRSE_ERROR, ' Found no <end_add_internals> in inPut file '//TRIM(InpFile))
2     CONTINUE
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
        READ(Inp,DEFAULT_CHR_FMT,END=11)Line
        CALL RemoveComments(Line)
        IF(LEN(TRIM(Line)) == 0) CYCLE
        LineLowCase = Line
        Call LowCase(LineLowCase)

        IF(INDEX(LineLowCase,'stre_a')/=0) THEN
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='STRE_A    '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:6)=(/0,0,0,1,0,0/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAU
          NConstr=NConstr+1
          CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "STRE_A = "//TRIM(FltToChar(Value))//" A")
        ELSE IF(INDEX(LineLowCase,'stre_b')/=0) THEN
          IF(Geos%Clone(1)%PBC%Dimen==1) CALL Halt('Extra coord stre_b while PBC dimension is 1')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='STRE_B    '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:6)=(/0,0,0,0,1,0/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAU
          NConstr=NConstr+1
          CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "STRE_B = "//TRIM(FltToChar(Value))//" A")
        ELSE IF(INDEX(LineLowCase,'stre_c')/=0) THEN
          IF(Geos%Clone(1)%PBC%Dimen<3) CALL Halt('Extra coord stre_c while PBC dimension is < 3')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='STRE_C    '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:6)=(/0,0,0,0,0,1/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAU
          NConstr=NConstr+1
          CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "STRE_C = "//TRIM(FltToChar(Value))//" A")
        ELSE IF(INDEX(LineLowCase,'stre')/=0) THEN
          !--------------------
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='STRE      '
          !--------------------
          IF(INDEX(LineLowCase,'.')==0) THEN
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:2)
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:2), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:6)
            ENDIF
          ELSE
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:2),Value
              GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
              GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAU
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:2), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:6),Value
              GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
              GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAU
            ENDIF
            NConstr=NConstr+1
          ENDIF
          !--------------------
        ELSE IF(INDEX(LineLowCase,'bend')/=0) THEN
          !--------------------
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='BEND      '
          IF(INDEX(LineLowCase,'.')==0) THEN
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3),&
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:9)
            ENDIF
          ELSE
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3),Value
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:3), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:9),Value
            ENDIF
            GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
            GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
            NConstr=NConstr+1
          ENDIF
          !--------------------
        ELSE IF(INDEX(LineLowCase,'gamma')/=0) THEN
          IF(Geos%Clone(1)%PBC%Dimen<2) CALL Halt('Extra coord gamma while PBC dimension is 1')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='GAMMA     '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/1,0,0,0,0,0,0,1,0/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
          NConstr=NConstr+1
          CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "Gamma = "//TRIM(FltToChar(Value))//" degree")
        ELSE IF(INDEX(LineLowCase,'area_l')/=0) THEN
          !--------------------
          IF(Geos%Clone(1)%PBC%Dimen/=2) CALL Halt('Extra coord area_l while PBC dimension is 1')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='AREA_L    '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/0,0,0,1,0,0,0,1,0/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)= Value*AngstromsToAU**2
          NConstr=NConstr+1
          !--------------------
        ELSE IF(INDEX(LineLowCase,'volm_l')/=0) THEN
          !--------------------
          IF(Geos%Clone(1)%PBC%Dimen/=3) CALL Halt('Extra coord volm_l while PBC dimension is /= 3')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='VOLM_L    '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:12)=(/0,0,0,1,0,0,0,1,0,0,0,1/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*AngstromsToAU**3
          NConstr=NConstr+1
          !--------------------
        ELSE IF(INDEX(LineLowCase,'beta')/=0) THEN
          IF(Geos%Clone(1)%PBC%Dimen/=3) CALL Halt('Extra coord beta while PBC dimension is /= 3')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='BETA      '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/1,0,0,0,0,0,0,0,1/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
          NConstr=NConstr+1
          CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "Beta = "//TRIM(FltToChar(Value))//" degree")
        ELSE IF(INDEX(LineLowCase,'alpha')/=0) THEN
          IF(Geos%Clone(1)%PBC%Dimen/=3) CALL Halt('Extra coord alpha while PBC dimension is /= 3')
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='ALPHA     '
          GOpt%ExtIntCs%Atoms%I(NIntCs,1:3)=1
          GOpt%ExtIntCs%Cells%I(NIntCs,1:9)=(/0,1,0,0,0,0,0,0,1/)
          GOpt%ExtIntCs%Active%L(NIntCs)=.TRUE.
          READ(LineLowCase,*) CHAR,Value
          GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
          GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
          NConstr=NConstr+1
          CALL MondoLog(DEBUG_MAXIMUM, "LoadExtraCoords", "Alpha = "//TRIM(FltToChar(Value))//" degree")
        ELSE IF(INDEX(LineLowCase,'tors')/=0) THEN
          !--------------------
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='TORS      '
          IF(INDEX(LineLowCase,'.')==0) THEN
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:12)
            ENDIF
          ELSE
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4),Value
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:12),Value
            ENDIF
            GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
            GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
            NConstr=NConstr+1
          ENDIF
          !--------------------
        ELSE IF(INDEX(LineLowCase,'outp')/=0) THEN
          !--------------------
          NIntCs=NIntCs+1
          GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='OUTP      '
          IF(INDEX(LineLowCase,'.')==0) THEN
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:12)
            ENDIF
          ELSE
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4),Value
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:12),Value
            ENDIF
            GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
            GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
            NConstr=NConstr+1
          ENDIF
          !--------------------
        ELSE IF(INDEX(LineLowCase,'linb1')/=0.OR. &
             INDEX(LineLowCase,'linb2')/=0) THEN
          !--------------------
          NIntCs=NIntCs+1
          IF(INDEX(LineLowCase,'linb1')/=0) THEN
            GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='LINB1     '
          ELSE
            GOpt%ExtIntCs%Def%C(NIntCs)(1:10)='LINB2     '
          ENDIF
          IF(INDEX(LineLowCase,'.')==0) THEN
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4)
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:12)
            ENDIF
          ELSE
            IF(INDEX(LineLowCase,'cell')==0) THEN
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4),Value
            ELSE
              READ(LineLowCase,*) &
                   CHAR,GOpt%ExtIntCs%Atoms%I(NIntCs,1:4), &
                   CharAux4,GOpt%ExtIntCs%Cells%I(NIntCs,1:12),Value
            ENDIF
            GOpt%ExtIntCs%Constraint%L(NIntCs)=.TRUE.
            GOpt%ExtIntCs%ConstrValue%D(NIntCs)=Value*DegToRad
            NConstr=NConstr+1
          ENDIF
          !--------------------
        ELSE IF(INDEX(LineLowCase,'cartx')/=0) THEN
          !--------------------
          GOpt%ExtIntCs%Def%C(NIntCs+1)(1:10)='CARTX     '
          READ(LineLowCase,*) CHAR,SerNum
          GOpt%ExtIntCs%Atoms%I(NIntCs+1,1)=SerNum
          GOpt%ExtIntCs%Constraint%L(NIntCs+1)=.TRUE.
          NConstr=NConstr+1
          Vect1=XYZ%D(1:3,SerNum)
          IF(Geos%Clone(1)%PBC%Dimen>0) THEN
            InvBoxSh=InverseBoxShape(Geos%Clone(1)%PBC%BoxShape%D,Geos%Clone(1)%PBC%Dimen)
            CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
            Vect1=Vect2
          ENDIF
          GOpt%ExtIntCs%ConstrValue%D(NIntCs+1)=Vect1(1)
          NIntCs=NIntCs+1
          NCartConstr=NCartConstr+1
        ELSE IF(INDEX(LineLowCase,'carty')/=0) THEN
          GOpt%ExtIntCs%Def%C(NIntCs+1)(1:10)='CARTY     '
          READ(LineLowCase,*) CHAR,SerNum
          GOpt%ExtIntCs%Atoms%I(NIntCs+1,1)=SerNum
          GOpt%ExtIntCs%Constraint%L(NIntCs+1)=.TRUE.
          NConstr=NConstr+1
          Vect1=XYZ%D(1:3,SerNum)
          IF(Geos%Clone(1)%PBC%Dimen>0) THEN
            InvBoxSh=InverseBoxShape(Geos%Clone(1)%PBC%BoxShape%D,Geos%Clone(1)%PBC%Dimen)
            CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
            Vect1=Vect2
          ENDIF
          GOpt%ExtIntCs%ConstrValue%D(NIntCs+1)=Vect1(2)
          NIntCs=NIntCs+1
          NCartConstr=NCartConstr+1
        ELSE IF(INDEX(LineLowCase,'cartz')/=0) THEN
          GOpt%ExtIntCs%Def%C(NIntCs+1)(1:10)='CARTZ     '
          READ(LineLowCase,*) CHAR,SerNum
          GOpt%ExtIntCs%Atoms%I(NIntCs+1,1)=SerNum
          GOpt%ExtIntCs%Constraint%L(NIntCs+1)=.TRUE.
          NConstr=NConstr+1
          Vect1=XYZ%D(1:3,SerNum)
          IF(Geos%Clone(1)%PBC%Dimen>0) THEN
            InvBoxSh=InverseBoxShape(Geos%Clone(1)%PBC%BoxShape%D,Geos%Clone(1)%PBC%Dimen)
            CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
            Vect1=Vect2
          ENDIF
          GOpt%ExtIntCs%ConstrValue%D(NIntCs+1)=Vect1(3)
          NIntCs=NIntCs+1
          NCartConstr=NCartConstr+1
        ELSE IF(INDEX(LineLowCase,'end_add_internals')/=0) THEN
          GO TO 12
        ELSE
          ! dont do anything with non-conform lines
          Call Warn("illegal input: "//TRIM(Line))
        ENDIF
      ENDDO
11    CALL MondoHalt(PRSE_ERROR, ' Found no <end_add_internals> in inPut file '//TRIM(InpFile))
12    CONTINUE
    ELSE
      NIntCs=0
      CALL New(GOpt%ExtIntCs,NIntCs)
    ENDIF !!! key for extra internals found
    ! CALL PrtIntCoords(GOPt%ExtIntCs,GOpt%ExtIntCs%Value%D,&
    !                   'chk extra 1',PBCDim_O=1)
    !
    ! Fill in Cartesian constraints stored
    ! in Geos%Clone(1)%CConstrain%I
    ! Again, supposing that constraints are the same for all clones.
    ! Then, renumber atoms in IntCs by exclusion of Rigid atoms
    !
    CALL MergeConstr(GOpt%ExtIntCs,XYZ%D,CConstrain%I, &
         Geos%Clone(1)%PBC%BoxShape%D, &
         Geos%Clone(1)%PBC%Dimen, &
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
  END SUBROUTINE LoadExtraCoords
  !
  !------------------------------------------------------------------
  !
  SUBROUTINE ReOrdIntC(IntCs,NIntC)
    TYPE(INTC)           :: IntCs
    INTEGER              :: NIntC,I,J,K,L,M,N
    INTEGER              :: Atoms(4),Cells(12)
    !
    ! Use similar ordering as in InCoords
    !
    DO I=1,NIntC
      Atoms(1:4)=IntCs%Atoms%I(I,1:4)
      Cells(1:12)=IntCs%Cells%I(I,1:12)
      IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
        IF(Atoms(1)>Atoms(2)) THEN
          DO J=1,2
            IntCs%Atoms%I(I,J)=Atoms(3-J)
          ENDDO
          DO J=1,2
            K=3*(4-J-1)+1
            L=K+2
            M=3*(J-1)+1
            N=M+2
            IntCs%Cells%I(I,M:N)=Cells(K:L)
          ENDDO
        ENDIF
      ELSE IF(IntCs%Def%C(I)(1:4)=='BEND'.OR. &
           IntCs%Def%C(I)(1:4)=='LINB') THEN
        IF(Atoms(1)>Atoms(3)) THEN
          DO J=1,3
            IntCs%Atoms%I(I,J)=Atoms(4-J)
          ENDDO
          DO J=1,3
            K=3*(4-J-1)+1
            L=K+2
            M=3*(J-1)+1
            N=M+2
            IntCs%Cells%I(I,M:N)=Cells(K:L)
          ENDDO
        ENDIF
      ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
        IF(Atoms(1)>Atoms(4)) THEN
          DO J=1,4
            IntCs%Atoms%I(I,J)=Atoms(5-J)
          ENDDO
          DO J=1,4
            K=3*(5-J-1)+1
            L=K+2
            M=3*(J-1)+1
            N=M+2
            IntCs%Cells%I(I,M:N)=Cells(K:L)
          ENDDO
        ENDIF
      ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
        IF(Atoms(3)>Atoms(4)) THEN
          IntCs%Atoms%I(I,3)=Atoms(4)
          IntCs%Atoms%I(I,4)=Atoms(3)
          IntCs%Cells%I(I,7:9)=Cells(10:12)
          IntCs%Cells%I(I,10:12)=Cells(7:9)
        ENDIF
      ENDIF
    ENDDO
  END SUBROUTINE ReOrdIntC
  !
  !------------------------------------------------------------------
  !
  SUBROUTINE MergeConstr(IntCs,XYZ,CConstrain,BoxShape,PBCDim, &
       NIntCs,NConstr,NCartConstr)
    TYPE(INTC)                  :: IntCs,IntC_New
    INTEGER,DIMENSION(:)        :: CConstrain
    REAL(DOUBLE),DIMENSION(:,:) :: XYZ
    REAL(DOUBLE)                :: BoxShape(3,3),InvBoxSh(3,3)
    REAL(DOUBLE)                :: Vect1(3),Vect2(3)
    INTEGER                     :: NIntCs,NConstr,NCartConstr
    INTEGER                     :: I,NNewC,NatmsLoc,PBCDim
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
        IntC_New%Def%C(NNewC+1)(1:10)='CARTX     '
        IntC_New%Def%C(NNewC+2)(1:10)='CARTY     '
        IntC_New%Def%C(NNewC+3)(1:10)='CARTZ     '
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
        Vect1=XYZ(1:3,I)
        IF(PBCDim>0) THEN
          InvBoxSh=InverseBoxShape(BoxShape,PBCDim)
          CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
          Vect1=Vect2
        ENDIF
        IntC_New%ConstrValue%D(NNewC+1)=Vect1(1)
        IntC_New%ConstrValue%D(NNewC+2)=Vect1(2)
        IntC_New%ConstrValue%D(NNewC+3)=Vect1(3)
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
    IF(IntCs%N==0)RETURN
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
