#ifdef MMech
      SUBROUTINE ParseMM(Ctrl)
        Implicit None
        TYPE(SCFControls)          :: Ctrl
        TYPE(CRDS)                 :: GM_MM,GM
        INTEGER :: II
        REAL(DOUBLE)               :: Sum,SumO,dx,dy,dz,dd
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: OPLS_DATABASE,MM_COORDS,MM_SEQ
        CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line, LineLowCase
        CHARACTER(LEN=DEFAULT_CHR_LEN)  :: QMMM_DATA
        TYPE(INT_VECT) :: ATMMARK,Active_Bond,Active_Angle
        TYPE(INT_VECT) :: MMAtNum
        TYPE(CHR_VECT) :: MMAtNam
        TYPE(INT_VECT) :: Active_Torsion,Active_OutOfPlane
        TYPE(INT_VECT) :: GlobalQMNum
        TYPE(DBL_VECT) :: AuxVect     
        INTEGER :: I,J,K,L,N1,N2,N3,N4
        TYPE(DBL_VECT) :: Charge14,LJEps14,LJEps,LJRAD
        TYPE(DBL_VECT) :: BondEQ,BondFC
        TYPE(DBL_VECT) :: AngleEQ,AngleFC
        TYPE(DBL_VECT) :: TorsionEQ,TorsionFC
        TYPE(DBL_VECT) :: OutOfPlaneEQ,OutOfPlaneFC
        TYPE(INT_RNK2) :: BondIJ
        TYPE(INT_VECT) :: OutOfPlanePeriod,TorsionPeriod
        TYPE(INT_VECT) :: Stat
        TYPE(DBL_VECT) :: GrdMM
        CHARACTER(LEN=6) :: CurG
        TYPE(INT_RNK2) :: AngleIJK,TorsionIJKL,OutOfPlaneIJKL
        REAL(DOUBLE) :: SUM_OLD,XTrans,YTrans,ZTrans
        INTEGER :: IA,IB,IC,IAtTrans,IAtOrig
!
        CALL OpenASCII(InpFile,Inp)
        CALL OpenASCII(OutFile,Out)
!
          CALL New(Stat,3)
          Stat%I=Ctrl%Current
          CurG=TRIM(IntToChar(Stat%I(3)))
          CALL Delete(Stat)
!
        IF(HasQM()) CALL Get(GM,Tag_O=CurG)
!  
! Find names of datafiles for MM calculation
!
         CALL AlignLowCase('begin_mm_filenames',Inp)
           READ(Inp,DEFAULT_CHR_FMT) OPLS_DATABASE
           READ(Inp,DEFAULT_CHR_FMT) MM_SEQ
           READ(Inp,DEFAULT_CHR_FMT) MM_COORDS
           IF(HasQM()) READ(Inp,DEFAULT_CHR_FMT) QMMM_DATA
         DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1)Line
           LineLowCase = Line
           Call LowCase(LineLowCase)
           IF(INDEX(LineLowCase,'end_mm_filenames')/=0) GO TO 2
         ENDDO
      1  CALL MondoHalt(PRSE_ERROR,' Found no <end_mm_filenames> in inPut file '//TRIM(InpFile))
      2  CONTINUE
!
! Check whether filenames are given 
!
         N1 = LEN_TRIM(OPLS_DATABASE) 
         N2 = LEN_TRIM(MM_SEQ)
         N3 = LEN_TRIM(MM_COORDS)
         N4 = LEN_TRIM(QMMM_DATA)
         IF(N1 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <opls_database> in inPut file '//TRIM(InpFile))
         IF(N2 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <mm_seq> in inPut file '//TRIM(InpFile))
         IF(N3 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <mm_coords> in inPut file '//TRIM(InpFile))
         IF(HasQM()) THEN
         IF(N4 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <qmmmdata> in inPut file '//TRIM(InpFile))
         ENDIF
!
! Now process MM data files
!
      CALL MM_FILE_PROCESS ( "tmp_opls_bin", OPLS_DATABASE(1:N1))
      CALL MM_SYSTEM_CONSTRUCT ( "tmp_opls_bin", MM_SEQ(1:N2) )
      CALL COORDINATES_READ ( MM_COORDS(1:N3) )
!
      CALL New(ATMMARK,MM_NATOMS)
      CALL New(Active_Bond,NBonds)
      CALL New(Active_Angle,NAngles)
      CALL New(Active_Torsion,NDIHEDRALS)
      CALL New(Active_OutOfPlane,NImpropers)
!
      IF(HasQM()) THEN
!
! Calling QMMMDATA_READ also modifies ATMCHG (MM Charges)
!
        CALL QMMMDATA_READ(QMMM_DATA,ATMMARK,Active_Bond,Active_Angle, &
        Active_Torsion,Active_OutOfPlane)
!
      ELSE
!
        ATMMARK%I(1:MM_NATOMS)=0
        Active_Bond%I(1:NBonds)=1
        Active_Angle%I(1:NAngles)=1
        Active_Torsion%I(1:NDIHEDRALS)=1
        Active_OutOfPlane%I(1:NImpropers)=1
!
      ENDIF
!
! SAVE 14Charges and other parameters
!
        CALL Put(NBonds,'MM_NBond')
        CALL Put(NAngles,'MM_NAngle')
        CALL Put(NDihedrals,'MM_NTorsion')
        CALL Put(NImpropers,'MM_NOutOfPlane')
!
      IF(NBONDS/=0) THEN
      CALL New(BondIJ,(/2,NBONDS/))
      CALL New(BondEQ,NBONDS)
      CALL New(BondFC,NBONDS)
        BondIJ%I(1,:)=BONDS(:)%I
        BondIJ%I(2,:)=BONDS(:)%J
        BondEQ%D(:)=BONDS(:)%EQ
        BondFC%D(:)=BONDS(:)%FC
        CALL Put(BondIJ,'MM_BondIJ')
        CALL Put(BondFC,'BondFC')
        CALL Put(BondEQ,'BondEQ')
      CALL Delete(BondIJ)
      CALL Delete(BondFC)
      CALL Delete(BondEQ)
      ENDIF
!
      IF(NAngles/=0) THEN
      CALL New(AngleIJK,(/3,NAngles/))
      CALL New(AngleEQ,NAngles)
      CALL New(AngleFC,NAngles)
        AngleIJK%I(1,:)=ANGLES(:)%I
        AngleIJK%I(2,:)=ANGLES(:)%J
        AngleIJK%I(3,:)=ANGLES(:)%K
        AngleEQ%D(:)=ANGLES(:)%EQ
        AngleFC%D(:)=ANGLES(:)%FC
        CALL Put(AngleIJK,'MM_AngleIJK')
        CALL Put(AngleEQ,'MM_AngleEQ')
        CALL Put(AngleFC,'MM_AngleFC')
      CALL Delete(AngleIJK)
      CALL Delete(AngleEQ)
      CALL Delete(AngleFC)
      ENDIF
!
      IF(NDihedrals/=0) THEN
      CALL New(TorsionIJKL,(/4,NDihedrals/))
      CALL New(TorsionEQ,NDihedrals)
      CALL New(TorsionFC,NDihedrals)
      CALL New(TorsionPeriod,NDihedrals)
        TorsionIJKL%I(1,:)=DIHEDRALS(:)%I
        TorsionIJKL%I(2,:)=DIHEDRALS(:)%J
        TorsionIJKL%I(3,:)=DIHEDRALS(:)%K
        TorsionIJKL%I(4,:)=DIHEDRALS(:)%L
        TorsionEQ%D(:)=DIHEDRALS(:)%PHASE
        TorsionFC%D(:)=DIHEDRALS(:)%FC
        TorsionPeriod%I(:)=DIHEDRALS(:)%PERIOD
        CALL Put(TorsionIJKL,'MM_TorsionIJKL')
        CALL Put(TorsionEQ,'TorsionEQ')
        CALL Put(TorsionFC,'TorsionFC')
        CALL Put(TorsionPeriod,'TorsionPeriod')
      CALL Delete(TorsionIJKL)
      CALL Delete(TorsionEQ)
      CALL Delete(TorsionFC)
      CALL Delete(TorsionPeriod)
      ENDIF
!
      IF(NImpropers/=0) THEN
      CALL New(OutOfPlaneIJKL,(/4,NImpropers/))
      CALL New(OutOfPlaneEQ,NImpropers)
      CALL New(OutOfPlaneFC,NImpropers)
      CALL New(OutOfPlanePeriod,NImpropers)
        OutOfPlaneIJKL%I(1,:)=Impropers(:)%I
        OutOfPlaneIJKL%I(2,:)=Impropers(:)%J
        OutOfPlaneIJKL%I(3,:)=Impropers(:)%K
        OutOfPlaneIJKL%I(4,:)=Impropers(:)%L
        OutOfPlaneEQ%D(:)=Impropers(:)%PHASE
        OutOfPlaneFC%D(:)=Impropers(:)%FC
        OutOfPlanePeriod%I(:)=Impropers(:)%PERIOD
        CALL Put(OutOfPlaneIJKL,'MM_OutOfPlaneIJKL')
        CALL Put(OutOfPlaneEQ,'OutOfPlaneEQ')
        CALL Put(OutOfPlaneFC,'OutOfPlaneFC')
        CALL Put(OutOfPlanePeriod,'OutOfPlanePeriod')
      CALL Delete(OutOfPlaneIJKL)
      CALL Delete(OutOfPlaneEQ)
      CALL Delete(OutOfPlaneFC)
      CALL Delete(OutOfPlanePeriod)
      ENDIF
!
      CALL New(Charge14,MM_NATOMS)
      CALL New(LJEps14,MM_NATOMS)
      CALL New(LJEps,MM_NATOMS)
      CALL New(LJRAD,MM_NATOMS)
        Charge14%D(:)=ATMCHG14(:)
        LJEps14%D(:)=ATMEps14(:)
        LJEps%D(:)=ATMEps(:)
        LJRAD%D(:)=ATMSIG(:)
        CALL Put(Charge14,'Charge14')
        CALL Put(LJEps14,'LJEps14')
        CALL Put(LJEps,'LJEps')
        CALL Put(LJRAD,'LJRAD')
      CALL Delete(Charge14)
      CALL Delete(LJEps14)
      CALL Delete(LJEps)
      CALL Delete(LJRAD)
!
! Save atomic numbers of MM atoms
!
        CALL New(MMAtNum,MM_NATOMS)
        MMAtNum%I=ATMNUM
        CALL Put(MMAtNum,'MMAtNum')
        CALL Delete(MMAtNum)
!
! Save MM atomic types of MM atoms
!
        CALL New(MMAtNam,MM_NATOMS)
        MMAtNam%C=ATMNam
        CALL Put(MMAtNam,'MMAtNam')
        CALL Delete(MMAtNam)
!
! Calculate topology matrices, which will be needed later
! for exclusion energy calculations. These matrices
! are based only on the topology given by the force-field.
! They do not involve any extra bonds, angles, etc.
!
      CALL TOPOLOGIES_MM(MM_NATOMS,NBonds, &
          Bonds(:)%I,Bonds(:)%J,Ctrl%Info)
!
! Calculate total MM Charge
!
        SUM=Zero
      DO I=1,MM_NATOMS
        SUM=SUM+ATMCHG(I)
      ENDDO
        GM_MM%TotCh=SUM
!
      GM_MM%NAtms = MM_Natoms
        CALL New(GM_MM)
!
!-----------------------------------------------------------
! Parse for periodic options, here logical options only
!-----------------------------------------------------------
!
!     GM_MM%InAu = .FALSE. !!! for MM
      GM_MM%Carts%D(:,:) = ATMCRD(:,:)
!
#ifdef PERIODIC
        CALL ParsePeriodic(Ctrl,GM_MM)
#endif
!
!     GM_MM%InAu = .TRUE. !!! for MM
!     IF(.NOT.PBC_On) GM_MM%Carts%D = AngstromsToAu*GM_MM%Carts%D
      GM_MM%Nkind = NTYPES
      GM_MM%AtMss%D(:) = ATMMAS(:)
      GM_MM%AtNum%D(:) = ATMCHG(:)
!
! Set GM_MM atomkinds
! WARNING! In the present version, for the QM/MM case
! kinds for the QM part are all the same!
!
      CALL FindKind(GM_MM)
!
! Fractional coordinates handling
!
        CALL ConvertCoords(GM_MM)
!
! Print out MM coordinates into outPut file
!
        WRITE(Out,*) 
        WRITE(Out,*) 'MM system Cartesian Coordinates in Angstroems:'
        WRITE(Out,*) 
        DO I=1,GM_MM%Natms
          WRITE(Out,120) I,ATMNAM(I),GM_MM%Carts%D(1:3,I)/AngstromsToAu
!         WRITE(Out,110) I,ATMNAM(I),GM_MM%Carts%D(1:3,I)/AngstromsToAu
        ENDDO
!       WRITE(Out,*) 
!       WRITE(Out,*) 'MM system Cartesian Coordinates in Bohrs:'
!       WRITE(Out,*) 
!       DO I=1,GM_MM%Natms
!         WRITE(Out,*) I,ATMNAM(I),GM_MM%Carts%D(1:3,I)
!       ENDDO
110   FORMAT(I7,2X,A8,3F30.16)
120   FORMAT(I7,2X,A8,3F12.6)
!
! Zero GrdMM (initial gradients) for the case HasMM 
!
      IF(HasMM().AND.Ctrl%Current(3)==1) THEN
        CALL New(GrdMM,3*GM_MM%Natms)
        GrdMM%D(:)=Zero
          CALL Put(GrdMM,'GradEMM',Tag_O=CurGeom)
        CALL Delete(GrdMM)
      ENDIF
!
! Calculate connection between QMMM and QM atom numbers
!
      IF(HasQM()) THEN
        II=0
            CALL New(GlobalQMNum,GM%NAtms)
            GlobalQMNum%I(:)=0 
!
        DO I=1,MM_Natoms
          IF(ATMMARK%I(I)==1) THEN
            II=II+1
            GlobalQMNum%I(II)=I
            CALL New(AuxVect,3)
            AuxVect%D=GM%Carts%D(1:3,II)-GM_MM%Carts%D(1:3,I)
            IF(DOT_PRODUCT(AuxVect%D,AuxVect%D)>0.001D0) THEN
              CALL Halt('Order of QM atoms is not the same as the one in QMMM atomlist')
            ENDIF
            CALL Delete(AuxVect)
          ENDIF
        ENDDO
!
! Link atoms have a GlobalQMNum of zero.
!
          CALL Put(GlobalQMNum,'GlobalQMNum',N_O=GM%NAtms)
          CALL Delete(GlobalQMNum)
      ENDIF
!
        CALL Put(GM_MM,Tag_O='GM_MM'//CurG)
          IF(MM_NATOMS/=0) THEN
            CALL Put(ATMMARK,'ATMMARK',N_O=MM_NATOMS)
          ENDIF
          IF(NBonds/=0) CALL Put(Active_Bond,'Active_Bond',N_O=NBonds)
          IF(NAngles/=0) CALL Put(Active_Angle,'Active_Angle',N_O=NAngles)
          IF(NDIHEDRALS/=0) CALL Put(Active_Torsion,'Active_Torsion',N_O=NDIHEDRALS)
          IF(NImpropers/=0) CALL Put(Active_OutOfPlane,'Active_OutOfPlane',N_O=NImpropers)
!
          CALL Delete(ATMMARK)
          CALL Delete(Active_Bond)
          CALL Delete(Active_Angle)
          CALL Delete(Active_Torsion)
          CALL Delete(Active_OutOfPlane)
!
          CALL Delete(GM_MM)
!
        CLOSE(Out,STATUS='KEEP')
        CLOSE(Inp,STATUS='KEEP')
!
      END SUBROUTINE ParseMM
#endif
!
!------------------------------------------------------------------------
!
#ifdef MMech
!
      SUBROUTINE ParseMech(Ctrl)
         TYPE(SCFControls)          :: Ctrl
!----------------------------------------------------------------------------
         CALL OpenASCII(InpFile,Inp)
         CALL OpenASCII(OutFile,Out)
!
!        Parse <OPTIONS> for QM_MM
!
         Ctrl%Mechanics(1) = .false.
         Ctrl%Mechanics(2) = .true.  !!!! default is QM 
         IF(    OptKeyQ(Inp,mechanics_option,pureMM))THEN
            Ctrl%Mechanics(1) = .true.
            Ctrl%Mechanics(2) = .false.
         ELSEIF(OptKeyQ(Inp,mechanics_option,pureQM))THEN
            Ctrl%Mechanics(1) = .false.
            Ctrl%Mechanics(2) = .true.
         ELSEIF(OptKeyQ(Inp,mechanics_option,QMandMM))THEN
            Ctrl%Mechanics(1) = .true.
            Ctrl%Mechanics(2) = .true.
         ENDIF
!
        CALL Put(Ctrl%Mechanics(1),'Ctrl_Mechanics1')
        CALL Put(Ctrl%Mechanics(2),'Ctrl_Mechanics2')
!
        CLOSE(Out,STATUS='KEEP')
        CLOSE(Inp,STATUS='KEEP')
!
      END SUBROUTINE ParseMech
#endif


!
!----------------------------------------------------------------------------
#ifdef MMech
!
      SUBROUTINE QMMMDATA_READ(FILE,ATMMARK,Active_Bond,Active_Angle, &
        Active_Torsion,Active_OutOfPlane)
      IMPLICIT NONE
      INTEGER III,K,L,M,N,N4,ISUB,IRES,II,J
      CHARACTER(LEN=DEFAULT_CHR_LEN)  :: FILE
      TYPE(INT_VECT) :: ATMMARK,Active_Bond,Active_Angle
      TYPE(INT_VECT) :: Active_Torsion,Active_OutOfPlane
!
      N4 = LEN_TRIM(FILE)
!
! Read QM-MM distinquishing list and modified Charges
! which must be zero on all QM atoms and 'scaled' for the MM ones.
!
      OPEN(UNIT=60,FILE=FILE(1:N4),STATUS='OLD')
      DO ISUB=1,NSUBSYS
        DO IRES=SUBIND(ISUB)+1,SUBIND(ISUB+1)
            READ(60,750) 
          DO I=RESIND(IRES)+1,RESIND(IRES+1)
            READ(60,700) ATMNAM(I),ATMMARK%I(I),ATMCHG(I)
          ENDDO
        ENDDO
      ENDDO
!
600   FORMAT(I10,3F15.7)
700   FORMAT(10X,3X,A8,I4,2F20.10)
750   FORMAT(I10,A8,F20.10)
!
! Bonds
!
         READ(60,*) 
      DO I=1,NBonds
       READ(60,61) II,K,L,Active_Bond%I(I)
      ENDDO
61    Format(I10,3X,2I5,3X,I5)
!
! Angles
!
         READ(60,*) 
      DO II=1,NAngles
        READ(60,62) III,I,J,K,Active_Angle%I(II)
      ENDDO
62    Format(I10,3X,3I5,3X,I5)
!
! dihedrals
!
         READ(60,*) 
      DO II=1,NDIHEDRALS
        READ(60,81) III,I,J,K,L,Active_Torsion%I(II)
      ENDDO
81  FORMAT(I10,3X,4I5,3X,I5,F12.4)
!
! Impropers
!
        READ(60,*) 
      DO II=1,NImpropers
        READ(60,81) III,I,J,K,L,Active_OutOfPlane%I(II)
      ENDDO
!
      CLOSE(60)
!
      END SUBROUTINE QMMMDATA_READ
#endif
