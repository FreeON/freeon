MODULE IntCoo
!
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE InOut
   USE MemMan
   USE SetXYZ
   USE ProcessControl
   USE PrettyPrint
   USE ParsingConstants
#ifdef MMech
!
IMPLICIT NONE
!
CONTAINS
!--------------------------------------------------------------
!
SUBROUTINE Topology_12(NAtoms_Loc,NBond,BondI,BondJ,Top12,InfFile)
! Set up a table which shows the atom numbers of Atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2) :: Top12_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBond,NAtoms_Loc,NMax12
INTEGER,DIMENSION(1:NBond) :: BondI,BondJ
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
!
    NMax12=5
    CALL New(Top12,(/NAtoms_Loc,NMax12+1/))
    Top12%I(:,:)=0 
!
    DO I=1,NBond
!
      II=BondI(I)
      JJ=BondJ(I)
      NI=Top12%I(II,1)
      NJ=Top12%I(JJ,1)
!
! check matrix Size, increase Size if necessary
      IF(NI>=NMax12 .OR. NJ>=NMax12) THEN
        NMax12=NMax12+5
        CALL New(Top12_2,(/NAtoms_Loc,NMax12+1/))
        Top12_2%I(1:NAtoms_Loc,1:NMax12+1)=0 
        Top12_2%I(1:NAtoms_Loc,1:NMax12+1-5)=Top12%I(1:NAtoms_Loc,1:NMax12+1-5)
        CALL Delete(Top12)
        CALL New(Top12,(/NAtoms_Loc,NMax12+1/))
        Top12%I(:,:)=Top12_2%I(:,:)
        CALL Delete(Top12_2)
      ENDIF
!
      IF(NI/=0) THEN
      DO M=1,NI         
        IF(Top12%I(II,M+1)==JJ) THEN
          EXIT
        ELSE
          Top12%I(II,1)=NI+1
          Top12%I(II,1+(NI+1))=JJ
          EXIT
        ENDIF
      ENDDO 
      ELSE
          Top12%I(II,1)=NI+1
          Top12%I(II,1+(NI+1))=JJ
      ENDIF
!
      IF(NJ/=0) THEN
      DO M=1,NJ         
        IF(Top12%I(JJ,M+1)==II) THEN
          EXIT
        ELSE
          Top12%I(JJ,1)=NJ+1
          Top12%I(JJ,1+(NJ+1))=II
          EXIT
        ENDIF
      ENDDO 
      ELSE
          Top12%I(JJ,1)=NJ+1
          Top12%I(JJ,1+(NJ+1))=II
      ENDIF
!
    ENDDO
!
    IF(PRESENT(InfFile)) THEN
      CALL Put(NMax12,'NMax12')
      CALL Put(Top12,'Top12')
    ENDIF
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
!
END SUBROUTINE Topology_12 
!--------------------------------------------------------------
!
SUBROUTINE Topology_13(NAtoms_Loc,Top12,Top13,InfFile)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2),OPTIONAL :: Top13
TYPE(INT_RNK2) :: Top13_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NAtoms_Loc,NMax13,NMax12,KK,IN12,JN12
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
!
    IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NAtoms_Loc,K/))
        CALL Get(Top12,'Top12')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix')
      ENDIF
    ELSE
      NMax12=Size(Top12%I,2)-1
    ENDIF
!
    NMax13=10
    K=NMax13+1
    CALL New(Top13,(/NAtoms_Loc,K/))
    Top13%I(1:NAtoms_Loc,1:NMax13+1)=0 
!
    DO II=1,NAtoms_Loc
      IN12=Top12%I(II,1)
      DO J=1,IN12
        JJ=Top12%I(II,J+1)
        JN12=Top12%I(JJ,1)
          DO K=1,JN12
          KK=Top12%I(JJ,K+1)
          IF(II/=KK) THEN
!
      NI=Top13%I(II,1)
!
! check matrix Size, increase Size if necessary
!
      IF(NI>=NMax13) THEN
        NMax13=NMax13+10
        CALL New(Top13_2,(/NAtoms_Loc,NMax13+1/))
        Top13_2%I(1:NAtoms_Loc,1:NMax13+1)=0 
        Top13_2%I(1:NAtoms_Loc,1:NMax13+1-10)=Top13%I(1:NAtoms_Loc,1:NMax13+1-10)
        CALL Delete(Top13)
        CALL New(Top13,(/NAtoms_Loc,NMax13+1/))
        Top13%I(1:NAtoms_Loc,1:NMax13+1)=Top13_2%I(1:NAtoms_Loc,1:NMax13+1)
        CALL Delete(Top13_2)
      ENDIF
!
      IF(NI/=0) THEN
        IF(ANY(Top13%I(II,2:NI+1)==KK)) THEN
          CYCLE
        ELSE
          Top13%I(II,1)=NI+1
          Top13%I(II,1+(NI+1))=KK
        ENDIF
      ELSE
          Top13%I(II,1)=NI+1
          Top13%I(II,1+(NI+1))=KK
      ENDIF
!
        ENDIF !!! II/=KK
        ENDDO !!!! KK
      ENDDO !!!! JJ
    ENDDO !!!! II
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
!
    IF(PRESENT(InfFile)) THEN
      CALL Put(NMax13,'NMax13')
      CALL Put(Top13,'Top13')
    ENDIF
!
    IF(.NOT.PRESENT(Top13)) CALL Delete(Top13)
!
END SUBROUTINE Topology_13 
!--------------------------------------------------------------
!
SUBROUTINE Topology_14(NAtoms_Loc,Top12,Top14,InfFile)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL:: Top12
TYPE(INT_RNK2),OPTIONAL:: Top14
TYPE(INT_RNK2) :: Top14_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
INTEGER :: NAtoms_Loc,NMax14,NMax12,IN12,JN12,KN12
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
!
     IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NAtoms_Loc,K/))
        CALL Get(Top12,'Top12')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix')
      ENDIF
    ELSE
      NMax12=Size(Top12%I,2)-1
    ENDIF
!
    NMax14=10
    K=NMax14+1
    CALL New(Top14,(/NAtoms_Loc,K/))
    Top14%I(1:NAtoms_Loc,1:NMax14+1)=0 
!
    DO II=1,NAtoms_Loc
      IN12=Top12%I(II,1)
      DO J=1,IN12
        JJ=Top12%I(II,J+1)
        JN12=Top12%I(JJ,1)
          DO K=1,JN12
          KK=Top12%I(JJ,K+1)
          IF(II/=KK) THEN
            KN12=Top12%I(KK,1)
              DO L=1,KN12
              LL=Top12%I(KK,L+1)
          IF(JJ/=LL.AND.II/=LL) THEN
!
      NI=Top14%I(II,1)
!
! check matrix Size, increase Size if necessary
!
      IF(NI>=NMax14) THEN
        NMax14=NMax14+10
        CALL New(Top14_2,(/NAtoms_Loc,NMax14+1/))
        Top14_2%I(1:NAtoms_Loc,1:NMax14+1)=0 
        Top14_2%I(1:NAtoms_Loc,1:NMax14+1-10)=Top14%I(1:NAtoms_Loc,1:NMax14+1-10)
        CALL Delete(Top14)
        CALL New(Top14,(/NAtoms_Loc,NMax14+1/))
        Top14%I(1:NAtoms_Loc,1:NMax14+1)=Top14_2%I(1:NAtoms_Loc,1:NMax14+1)
        CALL Delete(Top14_2)
      ENDIF
!
      IF(NI/=0) THEN
        IF(ANY(Top14%I(II,2:NI+1)==LL)) THEN
          CYCLE
        ELSE
          Top14%I(II,1)=NI+1
          Top14%I(II,1+(NI+1))=LL
        ENDIF
      ELSE
          Top14%I(II,1)=NI+1
          Top14%I(II,1+(NI+1))=LL
      ENDIF
!
          ENDIF !!! II/=LL and JJ/=LL
          ENDDO !!! LL
        ENDIF !!! II/=KK
        ENDDO !!!! KK
      ENDDO !!!! JJ
    ENDDO !!!! II
!
    IF(PRESENT(InfFile)) THEN
      CALL Put(NMax14,'NMax14')
      CALL Put(Top14,'Top14')
    ENDIF
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
    IF(.NOT.PRESENT(Top14)) CALL Delete(Top14)
!
END SUBROUTINE Topology_14 
!--------------------------------------------------------------
!
SUBROUTINE SORT_INTO_Box1(BoxSize,C,NAtoms_Loc,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
! sort the Atoms of a molecule into Boxes
!
! BoxSize: linear Box Size
!
IMPLICIT NONE
INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NAtoms_Loc
REAL(DOUBLE) :: BoxSize,VBIG,C(1:3,NAtoms_Loc),BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
SAVE VBIG
DATA VBIG/1.D+90/ 
!
! First count the number of Atoms in the individual Boxes 
!
!find borders of the global Box
!
    BXMIN= VBIG
    BXMax=-VBIG
    BYMIN= VBIG
    BYMax=-VBIG
    BZMIN= VBIG
    BZMax=-VBIG
  DO I=1,NAtoms_Loc
    IF(C(1,I)<BXMIN) BXMIN=C(1,I)
    IF(C(1,I)>BXMax) BXMax=C(1,I)
    IF(C(2,I)<BYMIN) BYMIN=C(2,I)
    IF(C(2,I)>BYMax) BYMax=C(2,I)
    IF(C(3,I)<BZMIN) BZMIN=C(3,I)
    IF(C(3,I)>BZMax) BZMax=C(3,I)
  ENDDO
!
  NX=INT((BXMax-BXMIN)/BoxSize)+1
  NY=INT((BYMax-BYMIN)/BoxSize)+1
  NZ=INT((BZMax-BZMIN)/BoxSize)+1
  NBox=NX*NY*NZ
!
END SUBROUTINE SORT_INTO_Box1
!
!--------------------------------------------------------------
SUBROUTINE SORT_INTO_Box2(BoxSize,C,NAtoms_Loc,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,InfFile,ISet)
!
! sort the Atoms of a molecule into Boxes
!
! BoxI(I) : contains the ordering number of the first atom of the I-th Box (like in sparse row-wise)
! BoxJ(J) : gives the original serial number of the atom desribed by the J-th ordering number
! C: contains Cartesian coordinates of Atoms
! BoxSize: linear Box Size
!
IMPLICIT NONE
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER,OPTIONAL :: ISET 
TYPE(INT_VECT),OPTIONAL :: BoxI1,BoxJ1
TYPE(INT_VECT) :: BoxI,BoxJ
INTEGER,ALLOCATABLE,DIMENSION(:) :: ISIGN
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: BoxCOUNTER 
INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NAtoms_Loc
REAL(DOUBLE) :: BoxSize,VBIG,C(1:3,NAtoms_Loc),BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
SAVE VBIG
DATA VBIG/1.D+90/ 
!
  NBox=NX*NY*NZ
  CALL New(BoxI,NBox+1)
  CALL New(BoxJ,NAtoms_Loc)
!
  ALLOCATE(ISIGN(1:NAtoms_Loc))
!
  ALLOCATE(BoxCOUNTER(1:NX,1:NY,1:NZ))
  BoxCOUNTER(1:NX,1:NY,1:NZ)=0
  BoxI%I(1:NBox+1)=0
!
! COUNT NUMBER OF Atoms IN THE Box
!
  DO I=1,NAtoms_Loc
!
! identify Box
!     
    IX=INT((C(1,I)-BXMIN)/BoxSize)+1
    IY=INT((C(2,I)-BYMIN)/BoxSize)+1
    IZ=INT((C(3,I)-BZMIN)/BoxSize)+1
!     
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY !order parameter: ZXY
    BoxCOUNTER(IX,IY,IZ)=BoxCOUNTER(IX,IY,IZ)+1
    ISIGN(I)=IORD*(NAtoms_Loc+1)+BoxCOUNTER(IX,IY,IZ) !!! shows both Box and index within Box
!     
  ENDDO
!     
! Count pointers to individual Boxes within array A
!     
     BoxI%I(1)=1
  DO IZ=1,NZ
  DO IX=1,NX
  DO IY=1,NY
!
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY
    BoxI%I(IORD+1)=BoxI%I(IORD)+BoxCOUNTER(IX,IY,IZ)
!     
  ENDDO
  ENDDO
  ENDDO
!
! Set up contents of Boxes as represented in A, sparse row-wise
!
  DO I=1,NAtoms_Loc
!
    IORD=INT(ISIGN(I)/(NAtoms_Loc+1))
    IADD=ISIGN(I)-IORD*(NAtoms_Loc+1)
    BoxJ%I(BoxI%I(IORD)-1+IADD)=I
!     
  ENDDO
!
! Now check ordering algorithm
!
! DO I=1,NBox  
! DO J=BoxI%I(I),BoxI%I(I+1)-1 !!! j is ordering index
!   JJ=BoxJ%I(J) !!! jj is original atom number
!   write(*,100) i,C(1:3,JJ)
! ENDDO
! ENDDO
100 format(I8,3F12.8)
!
  IF(PRESENT(BoxI1)) THEN
    BoxI1%I(:)=BoxI%I(:)
  ENDIF
  IF(PRESENT(BoxJ1)) THEN
    BoxJ1%I(:)=BoxJ%I(:)
  ENDIF
!
  IF(PRESENT(InfFile)) THEN
    CALL Put(NBox,'NBox'//TRIM(IntToChar(ISet)))
    CALL Put(BoxI,'BoxI'//TRIM(IntToChar(ISet)))
    CALL Put(BoxJ,'BoxJ'//TRIM(IntToChar(ISet)))
  ENDIF
!
  DEALLOCATE(BoxCOUNTER)
  CALL Delete(BoxI)
  CALL Delete(BoxJ)
!
END SUBROUTINE SORT_INTO_Box2
!
!----------------------------------------------------------------
!
SUBROUTINE SORT_INTO_Box(BoxSize,C,NAtoms_Loc,InfFile,ISet)
IMPLICIT NONE
INTEGER,OPTIONAL :: ISet
INTEGER :: NAtoms_Loc,NX,NY,NZ,NBox
REAL(DOUBLE) :: BoxSize,BXMIN,BYMIN,BZMIN
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
TYPE(INT_VECT) :: BoxI1,BoxJ1
REAL(DOUBLE),DIMENSION(1:3,1:NAtoms_Loc) :: C
!
CALL SORT_INTO_Box1(BoxSize,C,NAtoms_Loc,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
NBox=NX*NY*NZ
CALL New(BoxI1,NBox+1)
CALL New(BoxJ1,NAtoms_Loc)
!
CALL SORT_INTO_Box2(BoxSize,C,NAtoms_Loc,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,InfFile,ISet)
!
CALL Delete(BoxI1)
CALL Delete(BoxJ1)
!
END SUBROUTINE SORT_INTO_Box
!----------------------------------------------------------------
!
SUBROUTINE Topologies_MM(NAtoms_Loc,NBond,BondI,BondJ,InfFile,Top12OUT)
! Set up a table which shows the atom numbers of Atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12OUT
TYPE(INT_RNK2) :: Top12OUT_2
TYPE(INT_RNK2) :: Top12,Top13,Top14
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBond,NAtoms_Loc
INTEGER,DIMENSION(1:NBond) :: BondI,BondJ
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
!
   CALL Topology_12(NAtoms_Loc,NBond,BondI,BondJ,Top12,InfFile)
   CALL Topology_13(NAtoms_Loc,Top12,Top13,InfFile)
   CALL Topology_14(NAtoms_Loc,Top12,Top14,InfFile)
   CALL Excl_List(NAtoms_Loc,Top12,Top13,Top14,InfFile=InfFile)
   CALL Excl_List14(NAtoms_Loc,Top12,Top13,Top14,InfFile=InfFile)
!
   IF(PRESENT(Top12OUT)) THEN
     IF(AllocQ(Top12OUT%Alloc)) CALL Delete(Top12OUT)
     N=Size(Top12%I,2)
     CALL New(Top12OUT_2,(/NAtoms_Loc,N/))
     Top12OUT_2=Top12
     CALL Delete(Top12)
     CALL Delete(Top13)
     CALL Delete(Top14)
     CALL New(Top12OUT,(/NAtoms_Loc,N/))
     Top12OUT=Top12OUT_2
     CALL Delete(Top12OUT_2)
   ELSE
     CALL Delete(Top12)
     CALL Delete(Top13)
     CALL Delete(Top14)
   ENDIF
!
END SUBROUTINE Topologies_MM
!----------------------------------------------------------------
!
SUBROUTINE Excl_List(NAtoms_Loc,Top12,Top13,Top14,Top_Excl_Out,InfFile)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER :: NAtoms_Loc,NMax12,NMax13,NMax14,NMax_Excl,NNew,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ
!
    IF(PRESENT(Top12)) THEN
      NMax12=Size(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NAtoms_Loc,NMax12+1/))
        CALL Get(Top12,'Top12')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top13)) THEN
      NMax13=Size(Top13%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax13,'NMax13')
        CALL New(Top13,(/NAtoms_Loc,NMax13+1/))
        CALL Get(Top13,'Top13')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_13 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top14)) THEN
      NMax14=Size(Top14%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax14,'NMax14')
        CALL New(Top14,(/NAtoms_Loc,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl
!
    NMax_Excl=NMax12+NMax13+NMax14
    CALL New(Top_Excl,(/NAtoms_Loc,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
    Top_Excl%I(1:NAtoms_Loc,1:NMax12+1)=Top12%I(1:NAtoms_Loc,1:NMax12+1)
!
! Now merge Topologies, in order to avoid double counting in 
! Exclusion energies
!
      NMax_Excl_Out=0
    DO I=1,NAtoms_Loc
!
      NNew=0
      NOLD=Top_Excl%I(I,1)
      DO J=1,Top13%I(I,1)
        JJ=Top13%I(I,J+1)
        IF(ANY(Top_Excl%I(I,2:NOLD+1)==JJ)) THEN
          CYCLE
        ELSE
          NNew=NNew+1
          Top_Excl%I(I,NOLD+1+NNew)=JJ
        ENDIF
      ENDDO 
      NNew=NOLD+NNew
      Top_Excl%I(I,1)=NNew
      IF(NMax_Excl_Out<NNew) NMax_Excl_Out=NNew
!
      NNew=0
      NOLD=Top_Excl%I(I,1)
      DO J=1,Top14%I(I,1)
        JJ=Top14%I(I,J+1)
        IF(ANY(Top_Excl%I(I,2:NOLD+1)==JJ)) THEN
          CYCLE
        ELSE
          NNew=NNew+1
          Top_Excl%I(I,NOLD+1+NNew)=JJ
        ENDIF
      ENDDO 
      NNew=NOLD+NNew
      Top_Excl%I(I,1)=NNew
      IF(NMax_Excl_Out<NNew) NMax_Excl_Out=NNew
!
    ENDDO 
!
! STORE RESULT
!
    IF(PRESENT(Top_Excl_Out)) THEN
      CALL New(Top_Excl_Out,(/NAtoms_Loc,NMax_Excl_Out+1/))
      Top_Excl_Out%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)=&
      Top_Excl%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMax_Excl_Out,'NMax_Excl')
        CALL Put(Top_Excl_Out,'TOP_Excl')
      ENDIF
    ELSE IF(PRESENT(InfFile)) THEN
      CALL New(Top_New,(/NAtoms_Loc,NMax_Excl_Out+1/))
      Top_New%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)=&
      Top_Excl%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)
        CALL Put(NMax_Excl_Out,'NMax_Excl')
        CALL Put(Top_New,'TOP_Excl')
        CALL Delete(Top_New)
    ELSE
        CALL Halt('Do not know where to put exclusion topology')
    ENDIF
!
    CALL Delete(Top_Excl)
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
    IF(.NOT.PRESENT(Top13)) CALL Delete(Top13)
    IF(.NOT.PRESENT(Top14)) CALL Delete(Top14)
!
END SUBROUTINE Excl_LIST 
!
!----------------------------------------------------------------
!
SUBROUTINE Excl_List14(NAtoms_Loc,Top12,Top13,Top14,Top_Excl_Out,InfFile)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
! of Atoms in 14 distance. From the Top14 list
! Top13 and Top12 occurences must be filtered out
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER :: NAtoms_Loc,NMax12,NMax13,NMax14,NMax_Excl,NNew,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ,NExcl,N12,N13
!
    IF(PRESENT(Top12)) THEN
      NMax12=Size(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NAtoms_Loc,NMax12+1/))
        CALL Get(Top12,'Top12')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top13)) THEN
      NMax13=Size(Top13%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax13,'NMax13')
        CALL New(Top13,(/NAtoms_Loc,NMax13+1/))
        CALL Get(Top13,'Top13')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_13 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top14)) THEN
      NMax14=Size(Top14%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax14,'NMax14')
        CALL New(Top14,(/NAtoms_Loc,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl to the Size of Top14 and to zero 
!
    NMax_Excl=NMax14
    CALL New(Top_Excl,(/NAtoms_Loc,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
!
! Now merge Topologies, in order to avoid double counting in 
! Exclusion energies
!
      NMax_Excl_Out=0
    DO I=1,NAtoms_Loc
!
      DO J=1,Top14%I(I,1)
          NExcl=Top_Excl%I(I,1)
          N12=Top12%I(I,1)
          N13=Top13%I(I,1)
          JJ=Top14%I(I,J+1)
        IF(ANY(Top_Excl%I(I,2:NExcl+1)==JJ).OR. &
           ANY(Top12%I(I,2:N12+1)==JJ).OR. & 
           ANY(Top13%I(I,2:N13+1)==JJ)  ) THEN
          CYCLE
        ELSE
          Top_Excl%I(I,1)=NExcl+1
          Top_Excl%I(I,NExcl+2)=JJ
        ENDIF
      ENDDO 
          NExcl=Top_Excl%I(I,1)
      IF(NMax_Excl_Out<NExcl) NMax_Excl_Out=NExcl
!
    ENDDO 
!
! STORE RESULT
!
    IF(PRESENT(Top_Excl_Out)) THEN
      CALL New(Top_Excl_Out,(/NAtoms_Loc,NMax_Excl_Out+1/))
      Top_Excl_Out%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)=Top_Excl%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMax_Excl_Out,'NMax_Excl14')
        CALL Put(Top_Excl_Out,'TOP_Excl14')
      ENDIF
    ELSE
      CALL New(Top_New,(/NAtoms_Loc,NMax_Excl_Out+1/))
      Top_New%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)=Top_Excl%I(1:NAtoms_Loc,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMax_Excl_Out,'NMax_Excl14')
        CALL Put(Top_New,'TOP_Excl14')
      ENDIF
    ENDIF
!
    CALL Delete(Top_Excl)
    CALL Delete(Top_New)
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
    IF(.NOT.PRESENT(Top13)) CALL Delete(Top13)
    IF(.NOT.PRESENT(Top14)) CALL Delete(Top14)
!
END SUBROUTINE Excl_LIST14 
!----------------------------------------------------------------
      SUBROUTINE BMatrix(NAtoms_Loc,XYZ,NIntC,IntCs,B)
!
! This subroutine calculates the sparse, TYPE(BMATR) type 
! representation of the B matrix 
! In current version Cartesian coordinates 
! must be passed in in Angstroems.
! Later use au representation of everything, including B-matrix.
! Linear bendings _must_ always appear in pairs, Defd as LINB1 and LINB2
!
      IMPLICIT NONE
      INTEGER :: NIntC,NIntC2,I,J,K,L,NAtoms_Loc,IntCoo
      REAL(DOUBLE) :: Zero,AngstromsToAU,thresh_B
      REAL(DOUBLE),DIMENSION(1:3,1:NAtoms_Loc) :: XYZ
      TYPE(INTC) :: IntCs
      TYPE(BMATR):: B
!
      thresh_B=1.d-8
!
! allocate B matrix
!
      CALL NEW(B,NIntC)
        B%IB=0
        B%B=Zero
!
      DO IntCoo=1,NIntC     
!
        IF(IntCs%Def(IntCoo)(1:4)=='STRE') THEN
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2)
          CALL STRE(I,J,XYZ(1:3,I),XYZ(1:3,J),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
          write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
        IF(IntCs%Def(IntCoo)(1:4)=='BEND') THEN
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) !!! central atom
          K=IntCs%Atoms(IntCoo,3) 
       CALL BEND(I,J,K,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
          write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
! out of plane of i from the plane of jkl, with center at j
        IF(IntCs%Def(IntCoo)(1:4)=='OUTP') THEN
          I=IntCs%Atoms(IntCoo,1) !!! end atom
          J=IntCs%Atoms(IntCoo,2) !!! central atom
          K=IntCs%Atoms(IntCoo,3) !!! Def plane
          L=IntCs%Atoms(IntCoo,4) !!! Def plane
       CALL OUTP(I,J,K,L,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ(1:3,L),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
          write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
! torsion of i-j-k-l
        IF(IntCs%Def(IntCoo)(1:4)=='TORS') THEN
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) 
          K=IntCs%Atoms(IntCoo,3) 
          L=IntCs%Atoms(IntCoo,4) 
       CALL TORS(I,J,K,L,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ(1:3,L),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
          write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
! linear bendig of i-j-k
        IF(IntCs%Def(IntCoo)(1:5)=='LINB1') THEN
          IF(IntCs%Def(IntCoo+1)(1:5)/='LINB2') CALL Halt('LINB Definitions are not paired!')
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) 
          K=IntCs%Atoms(IntCoo,3) 
       CALL LINB(I,J,K,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),B%IB(IntCoo+1,1:12),B%B(IntCoo+1,1:12),thresh_B)
          write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
          write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo+1,1:12),B%B(IntCoo+1,1:12)
        ENDIF
        IF(IntCs%Def(IntCoo)(1:5)=='LINB2') CYCLE
!
      ENDDO !!!! loop over internal coords
100 format(A5,12I3,/,6F12.6,/,6F12.6)
!
      END SUBROUTINE BMatrix
!--------------------------------------------------------------------------------------
!
      SUBROUTINE  Stre (i,j,XI,XJ,IBB,BB,thresh_B)
!...
!...  Calculate the B matrix elements for a bond stretch as Defined
!...  by Wilson.
!...

      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,RIJ
      REAL(DOUBLE) :: dijsq,t,dij,thresh_B   
      INTEGER :: I,J,K,L,M,N
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)

      dijsq = 0.0d0

      DO m =1,3
         t = xj(m) - xi(m)
         rij(m) = t
         dijsq = dijsq + t * t
      end DO

      dij = sqrt(dijsq)

      DO m = 1, 3
         t = rij(m)
         if (dabs(t) > thresh_B) then
            IBB(m) = (3*(i-1)+m) 
            IBB(3+m) = (3*(j-1)+m) 
            t = t / dij
            BB(m) = -t
            BB(3+m) = t
         end if
      end DO

      RETURN
      END SUBROUTINE STRE
!
      SUBROUTINE Bend(I,J,K,XI,XJ,XK,IBB,BB,thresh_B)
!
!  this subroutine computes the b matrix elements of a valence
!  angle bending coordinate as Defined by wilson.
!  i and k are the numbers of the end Atoms.
!  j is the number of the central atom.
!  noat: number of Atoms
!  ic: serial number of internal coordinate
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,RJI,RJK,EJI,EJK
      REAL(DOUBLE) :: djisq,tp,t,djksq,dxsq,dx,DOtj,djk,dji
      REAL(DOUBLE) :: sinj,SMI,SMK,SUM,thresh_B
      INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
   10 djisq=0.d0
      djksq=0.d0
      DO 20 m=1,3
        tp=xj(m)
        t=xi(m)-tp
        rji(m)=t
        djisq=djisq+t*t
        t=xk(m)-tp
        rjk(m)=t
        djksq=djksq+t*t
   20 Continue        
      dji=dsqrt(djisq)
      djk=dsqrt(djksq)
      dx=1.d0
      DOtj=0.d0
      DO 30 m=1,3
      t=rji(m)/dji
      eji(m)=t
      tp=rjk(m)/djk
      ejk(m)=tp
   30 DOtj=DOtj+t*tp
      if(dabs(DOtj).gt.0.99995d0)go to 60
      sinj=dsqrt(1.d0-DOtj*DOtj)
      ii=3*(i-1)
      jj=3*(j-1)
      kk=3*(k-1)
      DO 40 m=1,3
      smi=dx*(DOtj*eji(m)-ejk(m))/(dji*sinj)
      IF(dabs(smi)>thresh_B) THEN
        IBB(m) = ii+m
        BB(m) = smi
      ENDIF
      smk=dx*(DOtj*ejk(m)-eji(m))/(djk*sinj)
      IF(dabs(smk)>thresh_B) THEN
        IBB(6+m) = kk+m
        BB(6+m) = smk
      ENDIF
      sum=smi+smk
      IF(dabs(sum)>thresh_B) THEN
        IBB(3+m) = jj+m
        BB(3+m) = -sum
      ENDIF
   40 CONTINUE
      RETURN
   50 ier=1
      RETURN
   60 ier=-1
      WRITE(6,1000)
      RETURN
!
 1000 format(' i-j-k is colinear - use linear bend')
!
      END SUBROUTINE BEND
!
      SUBROUTINE OutP(I,J,K,L,XI,XJ,XK,XL,IBB,BB,thresh_B)
!
!  this subroutine computes the b matrix elements for an out of
!  plane wagging coordinate as Defined by decius, mcintosh, michaelian
!  and peterson.  subroutine coded by m peterson, univ of toronto.
!  i is the end atom (atom wagged with respect to j-k-l plane).
!  j is the apex atom (Atoms i, k and l are attached to j).
!  k and l are the anchor Atoms (Define the j-k-l plane).
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,XL,RJI,RJK,RJL,EJI,EJK,EJL,C1
      REAL(DOUBLE) :: djisq,tp,t,djksq,djlsq,dx,DOtj
      REAL(DOUBLE) :: djk,dji,djl,cosi,cosk,cosl,tpp
      REAL(DOUBLE) :: sinsin,sini,dot,sint,cost,tant,cossin,sml
      REAL(DOUBLE) :: sinj,SMI,SMK,SUM,thresh_B
      INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
   10 djisq=0.d0
      djksq=0.d0
      djlsq=0.d0
      do 20 m=1,3
      tp=xj(m)
      t=xi(m)-tp
      rji(m)=t
      djisq=djisq+t*t
      t=xk(m)-tp
      rjk(m)=t
      djksq=djksq+t*t
      t=xl(m)-xj(m)
      rjl(m)=t
      djlsq=djlsq+t*t
   20 CONTINUE      
      dji=dsqrt(djisq)
      djk=dsqrt(djksq)
      djl=dsqrt(djlsq)
      dx=1.d0
      cosi=0.d0
      cosk=0.d0
      cosl=0.d0
      do 30 m=1,3
      t=rji(m)/dji
      eji(m)=t
      tp=rjk(m)/djk
      ejk(m)=tp
      tpp=rjl(m)/djl
      ejl(m)=tpp
      cosi=cosi+tp*tpp
      cosk=cosk+t*tpp
   30 cosl=cosl+t*tp
      if(dabs(cosi).gt.0.99995d0)go to 70
      sinsin=1.d0-cosi*cosi
      sini=dsqrt(sinsin)
      c1(1)=ejk(2)*ejl(3)-ejk(3)*ejl(2)
      c1(2)=ejk(3)*ejl(1)-ejk(1)*ejl(3)
      c1(3)=ejk(1)*ejl(2)-ejk(2)*ejl(1)
      dot=eji(1)*c1(1)+eji(2)*c1(2)+eji(3)*c1(3)
      sint=dot/sini
!     if(dabs(sint).gt.0.00005d0) then
!        write(6,1020)nob
!        write(6,'(15x,a36)')intch(nob)
!     end if
      if(dabs(sint).gt.0.99995d0)go to 80
      cost=dsqrt(1.d0-sint*sint)
      tant=sint/cost
      ii=3*(i-1)
      jj=3*(j-1)
      kk=3*(k-1)
      ll=3*(l-1)
      cossin=cost*sini
      do 50 m=1,3
      t=c1(m)/cossin
      smi=(t-tant*eji(m))/dji
      IF(dabs(smi)>thresh_B) THEN
        IBB(m) = ii+m
        BB(m) = dx*smi
      ENDIF
      smk=t*(cosi*cosk-cosl)/(sinsin*djk)
      IF(dabs(smk)>thresh_B) THEN
        IBB(6+m) = kk+m
        BB(6+m) = dx*smk
      ENDIF
      sml=t*(cosi*cosl-cosk)/(sinsin*djl)
      IF(dabs(sml)>thresh_B) THEN
        IBB(9+m) = ll+m
        BB(9+m) = dx*sml
      ENDIF
      sum=smi+smk+sml
      IF(dabs(sum)>thresh_B) THEN
        IBB(6+m) = jj+m
        BB(6+m) = -dx*sum
      ENDIF
   50 CONTINUE
      return
   60 ier=1
      return
   70 ier=-1
      write(6,1000)
      return
   80 ier=-1
      write(6,1010)
      return
!
 1000 format(/,2x,'<!> k-j-l is colinear (no plane Defined for wag of i)')
 1010 format(/,2x,'<!> i is perpendicular to j-k-l plane - use valence angle bends')
 1020 format(/,2x,'<!> warning: wag of a non-planar system at internal',/,15x,'coordinate no.',i5,':')
!
      END SUBROUTINE OUTP
!
      SUBROUTINE Tors(I,J,K,L,XI,XJ,XK,XL,IBB,BB,thresh_B)

!...
!...  Calculate the B matrix elements for torsion as Defined by
!...  R.L. Hilderbrandt, J. Mol. Spec., 44 (1972) 599.
!...
!...  i-j-k-l : torsion around j-k bond
!... 
!...  Two data cards are read :
!...   (1) contains ni atom numbers for the i-type Atoms
!...   (2) contains nl atom numbers for the l-type Atoms
!...
!...  iatom, latom : atom numbers for the i- and l-type Atoms (size: 5)
!...
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,XL,rij,rjk,rlk,eij,ejk,elk,cr,sj,sk
      REAL(DOUBLE) :: SMI,SMK,thresh_B,dij,dijsq,dx,djk,t,dxsq,djksq
      REAL(DOUBLE) :: cosj,sin2j,smj,dlksq,dlk,cosk,sin2k,sml 
      INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
      djksq = 0.0d0
      dxsq = 0.0d0
!
      do m = 1, 3
         sj(m) = 0.0d0
         sk(m) = 0.0d0
         t = xk(m) - xj(m)
         rjk(m) = t
         djksq = djksq + t * t
      end do
!
      djk = 1.0d0 / sqrt(djksq)
      dx = 1.0d0
!
      do m = 1, 3
        ejk(m) = rjk(m) * djk
      end do
!
      jj = 3 * (j - 1)
      kk = 3 * (k - 1)
!
!...
        dijsq=0.d0
      do 40 m=1,3
        t=xj(m)-xi(m)
        rij(m)=t
        dijsq=dijsq+t*t
   40 CONTINUE        
        dij=1.d0/dsqrt(dijsq)
        cosj=0.d0
      do 50 m=1,3
        t=rij(m)*dij
        eij(m)=t
        cosj=cosj-t*ejk(m)
   50 CONTINUE
      if(dabs(cosj).gt.0.99995d0)go to 120
      sin2j=(1.d0-cosj*cosj)
      ii=3*(i-1)
      cr(1)=eij(2)*ejk(3)-eij(3)*ejk(2)
      cr(2)=eij(3)*ejk(1)-eij(1)*ejk(3)
      cr(3)=eij(1)*ejk(2)-eij(2)*ejk(1)
      do 60 m=1,3
        t=cr(m)/sin2j
        smi=t*dij
      IF(dabs(smi)>thresh_B) THEN
        IBB(m) = ii+m
        BB(m) = -dx*smi
      ENDIF
        smk=t*cosj*djk
        sk(m)=sk(m)+smk
        smj=smi-smk
        sj(m)=sj(m)+smj
   60 CONTINUE
!
      dlksq=0.d0
      do 70 m=1,3
      t=xk(m)-xl(m)
      rlk(m)=t
   70 dlksq=dlksq+t*t
      dlk=1.d0/dsqrt(dlksq)
      cosk=0.d0
      do 80 m=1,3
        t=rlk(m)*dlk
        elk(m)=t
      cosk=cosk+ejk(m)*t
   80 CONTINUE
      if(dabs(cosk).gt.0.99995d0)go to 120
      sin2k=(1.d0-cosk*cosk)
      ll=3*(l-1)
      cr(1)=elk(3)*ejk(2)-elk(2)*ejk(3)
      cr(2)=elk(1)*ejk(3)-elk(3)*ejk(1)
      cr(3)=elk(2)*ejk(1)-elk(1)*ejk(2)
      do 90 m=1,3
        t=cr(m)/sin2k
        sml=t*dlk
        IF(dabs(sml)>thresh_B) THEN
          IBB(9+m) = ll+m
          BB(9+m) = -dx*sml
        ENDIF
        smj=t*cosk*djk
        sj(m)=sj(m)+smj
        smk=sml-smj
        sk(m)=sk(m)+smk
   90 CONTINUE
      do 100 m=1,3
        smj=sj(m)
      IF(dabs(smj)>thresh_B) THEN
        IBB(3+m) = jj+m
        BB(3+m) = dx*smj
      ENDIF
        smk=sk(m)
      IF(dabs(smk)>thresh_B) THEN
        IBB(6+m) = kk+m
        BB(6+m) = dx*smk
      ENDIF
  100 CONTINUE
      return
  110 ier=1
      return
  120 ier=-1
      write(6,1030)
      return
!
 1030 format(' i-j-k or j-k-l is colinear (no torsion possIBBle)')
!
      END SUBROUTINE TORS
!
      SUBROUTINE LinB(I,J,K,XI,XJ,XK,IBB1,BB1,IBB2,BB2,thresh_B)
!
!  this subroutine computes the b matrix elements for a linear bend
!  or for a pair of perpendicular linear bends.
!  i and k are the end Atoms.
!  j is the central atom.
!
!  a gives the cartesian coordinates of a point in space, such
!  that the vector from atom j to point a is perpendicular to
!  the line i-j-k and serves to orient the coordinates in space.
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,RJI,RJK,EJK,UP,UN,UNIT,A
      REAL(DOUBLE) :: t,SUM,djisq,djksq,djasq,tp,dji,djk,dja,dx
      REAL(DOUBLE) :: dotj,dotp,test,smi,smk
      REAL(DOUBLE) :: thresh_B
      INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER,IntCoo
      INTEGER :: IBB1(1:12),IBB2(1:12)
      REAL(DOUBLE) :: BB1(1:12),BB2(1:12)
!
! Set up coordinates of point a
!
        SUM=0.D0
      DO m=1,3
        RJI(M)=XI(M)-XJ(M)
        RJK(M)=XK(M)-XJ(M)
        SUM=RJI(M)**2+RJK(M)**2
      ENDDO 
        SUM=SQRT(SUM)
! form arbitrary vector rjk which is no more parallel with jk bond
        RJK(1)=RJK(1)+SUM 
        RJK(2)=RJK(2)+SUM 
        RJK(3)=RJK(3)-SUM 
! form vector a, perpendicular to the plane of rjk and rji
        A(1)= RJI(2)*RJK(3)-RJI(3)*RJK(2)
        A(2)=-RJI(1)*RJK(3)+RJI(3)*RJK(1)
        A(3)= RJI(1)*RJK(2)-RJI(2)*RJK(1)
        SUM=A(1)**2+A(2)**2+A(3)**2
        SUM=SQRT(SUM)
        a(1:3)=a(1:3)/SUM
! generate point a
        a(1:3)=xj(1:3)+a(1:3)        
!
      djisq=0.d0
      djksq=0.d0
      djasq=0.d0
      do 20 m=1,3
        tp=xj(m)
        t=xi(m)-tp
        rji(m)=t
        djisq=djisq+t*t
        t=xk(m)-tp
        rjk(m)=t
        djksq=djksq+t*t
        t=a(m)-tp
        un(m)=t
        djasq=djasq+t*t
   20 CONTINUE
      dji=dsqrt(djisq)
      djk=dsqrt(djksq)
      dja=dsqrt(djasq)
      dx=1.d0
      dotj=0.d0
      dotp=0.d0
      do 30 m=1,3
        t=rji(m)/dji
        tp=rjk(m)/djk
        ejk(m)=tp
        dotj=dotj+t*tp
        tp=un(m)/dja
        unit(m)=tp
        dotp=dotp+t*tp
   30 CONTINUE       
      test=dabs(dotj)-1.d0
      if(dabs(test).gt.0.00005d0)go to 70
      if(dabs(dotp).gt.0.00005d0)go to 80
      ii=3*(i-1)
      jj=3*(j-1)
      kk=3*(k-1)
      do 40 m=1,3
        t=unit(m)
        if(dabs(t)<thresh_B)go to 40
        t=-dx*t
        smi=t/dji
        IBB1(m) = ii+m
        BB1(m) = smi
        smk=t/djk
        IBB1(3+m) = jj+m
        BB1(3+m) = -smi-smk
        IBB1(6+m) = kk+m
        BB1(6+m) = smk
   40 continue
      up(1)=ejk(2)*unit(3)-ejk(3)*unit(2)
      up(2)=ejk(3)*unit(1)-ejk(1)*unit(3)
      up(3)=ejk(1)*unit(2)-ejk(2)*unit(1)
      do 50 m=1,3
        t=up(m)
        if(dabs(t)<thresh_B)go to 50
        t=-dx*t
        smi=t/dji
        IBB2(m) = ii+m
        BB2(m) = smi
        smk=t/djk
        IBB2(3+m) = jj+m
        BB2(3+m) = -smi-smk
        IBB2(6+m) = kk+m
        BB2(6+m) = smk
   50 continue
      return
   60 ier=1
      return
   70 ier=-1
      write(6,1020)
      return
   80 ier=-1
      write(6,1030)
      return
!
 1000 format(3g12.6)
 1010 format('+',86x,'a = (',2(f11.7,','),f11.7,')')
 1020 format(' i-j-k not colinear - use valence angle bend')
 1030 format(' atom a not perpendicular to i-j-k at j')
!
      END SUBROUTINE LinB
!
!----------------------------------------------------------------
!
      SUBROUTINE DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,IntSet,IntCs,NIntC)
!
! This routine defines internal coordinates
! being used in geometry manipulations.
! They may be STRE, BEND, TORS, OUTP and LINB1&LINB2 -s.
! right now the complete set of primitive internals 
! is defd. only.
! IntSet=1 : defines internal coords which will be present
!            for a long (entire) course of the optimization
! IntSet=2 : defines more temporary coordinates, coming
!            from more temporary VdW interactions
!
      IMPLICIT NONE
      INTEGER :: I,J,N,Natoms_Loc,NBond,NIntC,IntSet
      INTEGER :: NMax_Excl,NMax12,ILast
      REAL(DOUBLE),DIMENSION(1:3,1:Natoms_Loc) :: XYZ
      TYPE(INT_RNK2) :: BondIJ 
      TYPE(INT_RNK2) :: AngleIJK
      TYPE(INT_RNK2) :: TorsionIJKL
      INTEGER :: NX,NY,NZ,NBOX,NAngle,NTorsion
      REAL(DOUBLE) :: BXMIN,BYMIN,BZMIN,BoxSize,Fact
      TYPE(DBL_VECT) :: CritRad   
      TYPE(INT_RNK2) :: Top12,Top13,Top14,Top_Excl
      TYPE(INTC) :: IntCs
      TYPE(INT_VECT) :: BoxI,BoxJ,MMAtNum
      CHARACTER(LEN=DefAULT_CHR_LEN) :: InfFile
!
      NIntC=0
!
! first sort atoms into boxes      
!
   BoxSize=5.D0*AngstromsToAU !in A
   CALL SORT_INTO_Box1(BoxSize,XYZ,Natoms_Loc,&
                   NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
   NBox=NX*NY*NZ
   CALL New(BoxI,NBox+1)
   CALL New(BoxJ,Natoms_Loc)
!
   CALL SORT_INTO_Box2(BoxSize,XYZ,Natoms_Loc,&
                   NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI,BoxJ)
!
! now define _covalent_ bonding scheme, based on Slater-radii
!
      IF(IntSet==1) THEN
        N=SIZE(SLRadii,1)
        CALL New(CritRad,N)
        Fact=1.3D0 !!! Scaling factor for Slater Radii
        CritRad%D=Fact*SLRadii*AngstromsToAU
      ELSE
        N=SIZE(VDWRadii,1)
        CALL New(CritRad,N)
        Fact=1.0D0 !!! Scaling factor for VDW Radii
        CritRad%D=Fact*VDWRadii*AngstromsToAU
      ENDIF
!
      CALL BondList(Natoms_Loc,XYZ,NBond,MMAtNum, &
           BoxI,BoxJ,NBox,NX,NY,NZ,CritRad)
      CALL New(BondIJ,(/2,NBond/))
      CALL BondList(Natoms_Loc,XYZ,NBond,MMAtNum, &
           BoxI,BoxJ,NBox,NX,NY,NZ,CritRad,BondIJ)
!
      IF(IntSet==1) THEN
!
! Now define covalent topology matrices
!
        CALL Topology_12(NAtoms_Loc,NBond,BondIJ%I(1,1:NBond),BondIJ%I(2,1:NBond),Top12,InfFile=InfFile)
        CALL Topology_13(NAtoms_Loc,Top12,Top13)
        CALL Topology_14(NAtoms_Loc,Top12,Top14)
        CALL Excl_List(NAtoms_Loc,Top12,Top13,Top14,Top_Excl,InfFile=InfFile)
!
      write(*,*) 'Top12'
      do i=1,Natoms_Loc
        write(*,100) i,Top12%I(i,1),(Top12%I(i,j+1),j=1,Top12%I(i,1))
      enddo
      write(*,*) 'Top13'
      do i=1,Natoms_Loc
        write(*,100) i,Top13%I(i,1),(Top13%I(i,j+1),j=1,Top13%I(i,1))
      enddo
      write(*,*) 'Top14'
      do i=1,Natoms_Loc
        write(*,100) i,Top14%I(i,1),(Top14%I(i,j+1),j=1,Top14%I(i,1))
      enddo
      write(*,*) 'Top_Excl'
      do i=1,Natoms_Loc
        write(*,100) i,Top_Excl%I(i,1),(Top_Excl%I(i,j+1),j=1,Top_Excl%I(i,1))
      enddo
!
        CALL Delete(Top_Excl)
        CALL Delete(Top14)
        CALL Delete(Top13)
!
      ELSE
!
! Filter out covalent 12, 13 and 14 connections from VWD bondlist
!
write(*,*) 'bef get excl  '  
        CALL Get(NMax_Excl,'NMax_Excl')
        CALL New(Top_Excl,(/Natoms_Loc,NMax_Excl+1/))
        CALL Get(Top_Excl,'TOP_EXCL')
write(*,*) 'bef vdw filter'  
          CALL VDWFilter(BondIJ,NBond,Top_Excl)
write(*,*) 'aft vdw filter'  
        CALL Delete(Top_Excl)
!     
      ENDIF
!
! Now define bond angles and torsions
!
      IF(IntSet==1) THEN
!
        CALL AngleList(NAtoms_Loc,Top12,NAngle=NAngle)
        CALL New(AngleIJK,(/3,NAngle/))
        CALL AngleList(NAtoms_Loc,Top12,AngleIJK,NAngle)
!
        CALL TorsionList(NAtoms_Loc,Top12,NTorsion=NTorsion)
        CALL New(TorsionIJKL,(/4,NTorsion/))
        CALL TorsionList(NAtoms_Loc,Top12,TorsionIJKL,NTorsion)
!
        CALL Delete(Top12)
!
      ELSE
!
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/Natoms_Loc,NMax12+1/))
        CALL Get(Top12,'TOP12')
!
write(*,*) 'bef vdw angles'  
        CALL VDWAngleList(NAtoms_Loc,Top12,AngleIJK,NAngle,BondIJ,NBond)
!
write(*,*) 'bef vdw torsions'  
        CALL VDWTorsionList(NAtoms_Loc,Top12,TorsionIJKL,NTorsion,BondIJ,NBond)
!
        CALL Delete(Top12)
!
      ENDIF
!
      write(*,*) 'NBond= ',NBond 
      do i=1,NBond      
        write(*,110) i,(BondIJ%I(j,i),j=1,2)
      enddo
      write(*,*) 'NAngle= ',NAngle
      do i=1,NAngle     
        write(*,110) i,(AngleIJK%I(j,i),j=1,3)
      enddo
      write(*,*) 'NTorsion= ',NTorsion
      do i=1,NTorsion   
        write(*,110) i,(TorsionIJKL%I(j,i),j=1,4)
      enddo
100   format(I4,2X,I4,2X,20I4)
110   format(I4,2X,20I4)
!
! Fill Data into IntCs
!
      NIntC=NBond+NAngle+NTorsion
      CALL New(IntCs,NIntC)
      IntCs%Def='     '
      IntCs%Atoms(1:NIntC,1:4)=0   
      IntCs%Value=Zero
      IntCs%Constraint=.FALSE.
!
      DO I=1,NBond
        IntCs%Def(I)='STRE '
        IntCs%Atoms(I,1:2)=BondIJ%I(1:2,I)
      ENDDO
        ILast=NBond
      DO I=1,NAngle
        IntCs%Def(ILast+I)='BEND '
        IntCs%Atoms(ILast+I,1:3)=AngleIJK%I(1:3,I)
      ENDDO
        ILast=NBond+NAngle
      DO I=1,NTorsion
        IntCs%Def(ILast+I)='TORS '
        IntCs%Atoms(ILast+I,1:4)=TorsionIJKL%I(1:4,I)
      ENDDO
!
write(*,*) 'after intcs defined:'
    do i=1,NIntC
    write(*,210) i,IntCs%Def(i),IntCs%Atoms(i,1:4)
    enddo
210 format(i3,' ',A5,4(I3,' '),12F8.3)
!
! tidy up
!
      CALL Delete(TorsionIJKL)
      CALL Delete(AngleIJK)
      CALL Delete(CritRad)
      CALL Delete(BondIJ)
      CALL Delete(BoxI)
      CALL Delete(BoxJ)
!
      END SUBROUTINE DefineIntCoos
!
!--------------------------------------------------------
      SUBROUTINE BondList(NAtoms_Loc,XYZ,NBond,MMAtNum,&
             BoxI,BoxJ,NBox,NX,NY,NZ,CritRad,BondIJ)
!
! Define number of bonds
!
      IMPLICIT NONE
      INTEGER :: I,J,Natoms_Loc,NBond
      REAL(DOUBLE),DIMENSION(1:3,1:Natoms_Loc) :: XYZ
      TYPE(DBL_VECT) :: CritRad !!! VdW or Slater Radii
      REAL(DOUBLE) :: R12,R12_2
      INTEGER :: IZ,IX,IY,I1,I2,JJ1,JJ2,NBOX,IORD,IORDD
      INTEGER :: IZD,IXD,IYD,NX,NY,NZ,NJJ1,NJJ2
      LOGICAL :: FillBondIJ
      TYPE(INT_RNK2),OPTIONAL :: BondIJ
      TYPE(DBL_VECT) :: DVect
      TYPE(INT_VECT) :: BoxI,BoxJ,MMATNUM
!
      FillBondIJ=.FALSE.
      IF(PRESENT(BondIJ)) FillBondIJ=.TRUE.
!
      CALL New(DVect,3)
!
! Go through all boxes 
!
   NBond=0
   DO IZ=1,NZ
   DO IX=1,NX
   DO IY=1,NY
!
! indices of central and neighboring Boxes
!
       IOrd=NX*NY*(IZ-1)+NY*(IX-1)+IY
!
     DO I1=BoxI%I(IOrd),BoxI%I(IOrd+1)-1
       JJ1=BoxJ%I(I1) !!! atom in central box
       NJJ1=MMAtNum%I(JJ1)
! second atom may come from central or neigbouring Boxes
! and must be an MM atom, LJ is not calculated for QM-QM pairs
     DO IZD=-1,1
       IF(IZ+IZD>0 .AND. IZ+IZD<=NZ) THEN
     DO IXD=-1,1
       IF(IX+IXD>0 .AND. IX+IXD<=NX) THEN
     DO IYD=-1,1
       IF(IY+IYD>0 .AND. IY+IYD<=NY) THEN
!
       IOrdD=NX*NY*(IZ-1+IZD)+NY*(IX-1+IXD)+IY+IYD
         DO I2=BoxI%I(IOrdD),BoxI%I(IOrdD+1)-1
           JJ2=BoxJ%I(I2) !!! second atom
           NJJ2=MMAtNum%I(JJ2)
       IF(JJ2<=JJ1) CYCLE !!!avoid double counting
           DVect%D(:)=XYZ(:,JJ1)-XYZ(:,JJ2)
           r12_2=DOT_PRODUCT(DVect%D,DVect%D)
           R12=SQRT(r12_2)
           IF(R12<CritRad%D(NJJ1)+CritRad%D(NJJ2)) THEN
             NBond=NBond+1
             IF(FillBondIJ) THEN
               BondIJ%I(1,NBond)=JJ1
               BondIJ%I(2,NBond)=JJ2
               write(*,*) 'bond= ',NBond,' atom1= ',jj1,' atom2= ',jj2
             ENDIF 
           ENDIF
         ENDDO
!
! ends on neighboring boxes
!
       ENDIF
     ENDDO
       ENDIF
     ENDDO
       ENDIF
     ENDDO
!
     ENDDO !!! central box atoms
!
   ENDDO
   ENDDO
   ENDDO !!! ends on central box indices
!
! tidy up
!
   CALL Delete(DVect)      
!
      END SUBROUTINE BondList
!
!----------------------------------------------------------------
SUBROUTINE AngleList(NAtoms_Loc,Top12,AngleIJK,NAngle)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2) :: Top12
TYPE(INT_RNK2),OPTIONAL :: AngleIJK 
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ
INTEGER :: NAtoms_Loc,KK,IN12,JN12,NAngle
LOGICAL :: AngleFill
!
            AngleFill=.FALSE.
            IF(PRESENT(AngleIJK)) AngleFill=.TRUE.
            NAngle=0
    DO II=1,NAtoms_Loc
      IN12=Top12%I(II,1)
      DO J=1,IN12
        JJ=Top12%I(II,J+1)
        JN12=Top12%I(JJ,1)
          DO K=1,JN12
          KK=Top12%I(JJ,K+1)
          IF(II<KK) THEN
!
            NAngle=NAngle+1
              IF(AngleFill) THEN
                AngleIJK%I(1,NAngle)=II
                AngleIJK%I(2,NAngle)=JJ
                AngleIJK%I(3,NAngle)=KK
              ENDIF
!
          ENDIF !!! II/=KK
        ENDDO !!!! KK
      ENDDO !!!! JJ
    ENDDO !!!! II
!
END SUBROUTINE AngleList   
!------------------------------------------------------
SUBROUTINE TorsionList(NAtoms_Loc,Top12,TorsionIJKL,NTorsion)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2):: Top12
TYPE(INT_RNK2),OPTIONAL :: TorsionIJKL
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL,NTorsion
INTEGER :: NAtoms_Loc,IN12,JN12,KN12
LOGICAL :: TorsionFill
!
            TorsionFill=.FALSE.
            IF(PRESENT(TorsionIJKL)) TorsionFill=.TRUE.
            NTorsion=0
    DO II=1,NAtoms_Loc
      IN12=Top12%I(II,1)
      DO J=1,IN12
        JJ=Top12%I(II,J+1)
        JN12=Top12%I(JJ,1)
          DO K=1,JN12
          KK=Top12%I(JJ,K+1)
          IF(II/=KK) THEN
            KN12=Top12%I(KK,1)
              DO L=1,KN12
              LL=Top12%I(KK,L+1)
          IF(JJ/=LL.AND.II<LL) THEN
!
            NTorsion=NTorsion+1
            IF(TorsionFill) THEN
              TorsionIJKL%I(1,NTorsion)=II
              TorsionIJKL%I(2,NTorsion)=JJ
              TorsionIJKL%I(3,NTorsion)=KK
              TorsionIJKL%I(4,NTorsion)=LL
            ENDIF
!
          ENDIF !!! II/=LL and JJ/=LL
          ENDDO !!! LL
        ENDIF !!! II/=KK
        ENDDO !!!! KK
      ENDDO !!!! JJ
    ENDDO !!!! II
!
END SUBROUTINE TorsionList 
!-------------------------------------------------
SUBROUTINE VDWFilter(BondIJ,NBond,Top_Excl)
!
    IMPLICIT NONE
    INTEGER :: I,N,I1,I2,J,K,L,NBond,NBond2
    TYPE(INT_RNK2) :: BondIJ,BondIJ2,Top_Excl
!
    CALL New(BondIJ2,(/2,NBond/)) 
!
      NBond2=0
    DO I=1,NBond
      I1=BondIJ%I(1,I)
      I2=BondIJ%I(2,I)
      N=Top_Excl%I(I1,1)
      IF(ANY(Top_Excl%I(I1,2:N+1)==I2)) CYCLE
      NBond2=NBond2+1
      BondIJ2%I(1,NBond2)=I1
      BondIJ2%I(2,NBond2)=I2
    ENDDO
!
    CALL Delete(BondIJ) 
    NBond=NBond2
    CALL New(BondIJ,(/2,NBond/))
    BONDIJ%I(1:2,1:NBond)=BONDIJ2%I(1:2,1:NBond)
!
    CALL Delete(BondIJ2) 
!
END SUBROUTINE VDWFilter
!------------------------------------------------------------
SUBROUTINE VDWAngleList(NAtoms_Loc,Top12,AngleIJK,NAngle,BondIJ,NBond)
! this routine generates bond-angles associated with WDV bonds
!
    IMPLICIT NONE
    INTEGER :: NAtoms_Loc,I,I1,I2,N1,N2,J,J1,J2,NBond,NAngle
    TYPE(INT_RNK2) :: Top12   
    TYPE(INT_RNK2) :: AngleIJK
    TYPE(INT_RNK2) :: BondIJ
!
      NAngle=0
    DO I=1,NBond
      I1=BondIJ%I(1,I)
      I2=BondIJ%I(2,I)
      N1=Top12%I(I1,1)
      N2=Top12%I(I2,1)
      NAngle=NAngle+N1+N2
    ENDDO
!
    CALL New(AngleIJK,(/3,NAngle/))
!
      NAngle=0
    DO I=1,NBond
      I1=BondIJ%I(1,I)
      I2=BondIJ%I(2,I)
      N1=Top12%I(I1,1)
      N2=Top12%I(I2,1)
      DO J=1,N1
        NAngle=NAngle+1
        AngleIJK%I(1,NAngle)=Top12%I(I1,J+1)
        AngleIJK%I(2,NAngle)=I1
        AngleIJK%I(3,NAngle)=I2
      ENDDO
      DO J=1,N2
        NAngle=NAngle+1
        AngleIJK%I(1,NAngle)=I1
        AngleIJK%I(2,NAngle)=I2
        AngleIJK%I(3,NAngle)=Top12%I(I2,J+1)
      ENDDO
    ENDDO
!
END SUBROUTINE VDWAngleList
!-------------------------------------------------
SUBROUTINE VDWTorsionList(NAtoms_Loc,Top12,TorsionIJKL,NTorsion,BondIJ,NBond)
! this routine generates bond-angles associated with WDV bonds
!
    IMPLICIT NONE
    INTEGER :: NAtoms_Loc,I,I1,I2,J,J1,J2,N1,N2,NBond,NTorsion
    TYPE(INT_RNK2) :: Top12
    TYPE(INT_RNK2) :: TorsionIJKL
    TYPE(INT_RNK2) :: BondIJ
!
      NTorsion=0
    DO I=1,NBond
      I1=BondIJ%I(1,I)
      I2=BondIJ%I(2,I)
      N1=Top12%I(I1,1)
      N2=Top12%I(I2,1)
      NTorsion=NTorsion+N1*N2
    ENDDO
!
    CALL New(TorsionIJKL,(/4,NTorsion/))
!
      NTorsion=0
    DO I=1,NBond
      I1=BondIJ%I(1,I)
      I2=BondIJ%I(2,I)
      N1=Top12%I(I1,1)
      N2=Top12%I(I2,1)
      DO J1=1,N1
      DO J2=1,N2
        NTorsion=NTorsion+1
        TorsionIJKL%I(1,NTorsion)=Top12%I(I1,J1+1)
        TorsionIJKL%I(2,NTorsion)=I1
        TorsionIJKL%I(3,NTorsion)=I2
        TorsionIJKL%I(4,NTorsion)=Top12%I(I2,J2+1)
      ENDDO
      ENDDO
    ENDDO
!
END SUBROUTINE VDWTorsionList
!--------------------------------------------------------
SUBROUTINE GetIntCs(XYZ,Natoms_Loc,InfFile,IntCs,NIntC,Refresh)
!
! This subroutine constructs the IntCs array, which holds
! definitions of internal coordinates to be used in the 
! forthcoming geometry manipulation procedure.
! Refresh=1 : Refresh all definitions
!        =2 : Refresh only definitions based on VDW interaction
!        =3 : Do not refresh definitions, use the one from HDF
! 
      IMPLICIT NONE
      TYPE(INTC) :: IntCs,IntC_Cov,IntC_VDW,IntC_Extra,IntC_New
      INTEGER :: NIntC,NIntC_Cov,NIntC_VDW,NIntC_Extra,NNew
      INTEGER :: I,J,Refresh,Natoms_Loc,II,ILast
      TYPE(INT_VECT) :: MMAtNum
      TYPE(CRDS) :: GMLoc
      CHARACTER(LEN=DefAULT_CHR_LEN) :: InfFile
      REAL(DOUBLE),DIMENSION(1:3,1:Natoms_Loc) :: XYZ
!
! Get atomnames (numbers) from HDF 
!
      CALL New(MMAtNum,Natoms_Loc)
      CALL Get(MMAtNum,'MMATNUM')
!
      IF(Refresh==1) Then !!! Total refresh
!
!define covalent bonding scheme
write(*,*) 'int c 1'
        CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,1,IntC_Cov,NIntC_Cov)
!define Van der Waals bonding scheme
write(*,*) 'int c 2'
        CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,2,IntC_VDW,NIntC_VDW)
!get extra internal coordinates and constraints
          CALL Get(NIntC_Extra,'NIntC_Extra')
          IF(NIntC_Extra/=0) THEN
            CALL New(IntC_Extra,NIntC_Extra)
            CALL Get(IntC_Extra,'IntC_Extra')
          ENDIF
          CALL Put(NIntC_Cov,'NIntC_Cov')
          IF(NIntC_Cov/=0) CALL Put(IntC_Cov,'IntC_Cov')
          CALL Put(NIntC_VDW,'NIntC_VDW')
          IF(NIntC_VDW/=0) CALL Put(IntC_VDW,'IntC_VDW')
!
      ELSE IF(Refresh==2) THEN !!! refresh VDW terms
!
! Refresh only the VDW terms
!
          CALL Get(NIntC_Extra,'NIntC_Extra')
          IF(NIntC_Extra/=0) THEN
            CALL New(IntC_Extra,NIntC_Extra)
            CALL Get(IntC_Extra,'IntC_Extra')
          ENDIF
          CALL Get(NIntC_Cov,'NIntC_Cov')
          CALL New(IntC_Cov,NIntC_Cov)
          CALL Get(IntC_Cov,'IntC_Cov')
!       define Van der Waals bonding scheme
        CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,2,IntC_VDW,NIntC_VDW)
          CALL Put(NIntC_VDW,'NIntC_VDW')
          IF(NIntC_VDW/=0) CALL Put(IntC_VDW,'IntC_VDW')
!
      ELSE IF(Refresh==3) THEN !!! no refresh, get everything from HDF
!
          CALL Get(NIntC_Extra,'NIntC_Extra')
          IF(NIntC_Extra/=0) THEN
            CALL New(IntC_Extra,NIntC_Extra)
            CALL Get(IntC_Extra,'IntC_Extra')
          ENDIF
write(*,*) 'after get extra'
          CALL Get(NIntC_Cov,'NIntC_Cov')
          CALL New(IntC_Cov,NIntC_Cov)
          CALL Get(IntC_Cov,'IntC_Cov')
          CALL Get(NIntC_VDW,'NIntC_VDW')
          CALL New(IntC_VDW,NIntC_VDW)
          CALL Get(IntC_VDW,'IntC_VDW')
!
      ENDIF
write(*,*) 'INTC_Cov '
    do i=1,NIntC
    write(*,110) i,IntC_Cov%Def(i),IntC_Cov%Atoms(i,1:4)
    enddo
write(*,*) 'INTC_VDW '
    do i=1,NIntC
    write(*,110) i,IntC_VDW%Def(i),IntC_VDW%Atoms(i,1:4)
    enddo
110 format(i3,' ',A5,4(I3,' '),12F8.3)
!
! Merge INTC arrays
!
        NIntC=NIntC_Cov+NIntC_VDW+NIntC_Extra
write(*,*) 'NIntC= ',NIntC
        CALL New(IntCs,NIntC)
!
          ILast=0
        IntCs%Def(ILast+1:ILast+NIntC_Cov)=IntC_Cov%Def(1:NIntC_Cov)
   IntCs%Atoms(ILast+1:ILast+NIntC_Cov,1:4)=IntC_Cov%Atoms(1:NIntC_Cov,1:4)
        IntCs%Value(ILast+1:ILast+NIntC_Cov)=IntC_Cov%Value(1:NIntC_Cov)
        IntCs%Constraint(ILast+1:ILast+NIntC_Cov)=IntC_Cov%Constraint(1:NIntC_Cov)
          ILast=NIntC_Cov
        IntCs%Def(ILast+1:ILast+NIntC_VDW)=IntC_VDW%Def(1:NIntC_VDW)
   IntCs%Atoms(ILast+1:ILast+NIntC_VDW,1:4)=IntC_VDW%Atoms(1:NIntC_VDW,1:4)
        IntCs%Value(ILast+1:ILast+NIntC_VDW)=IntC_VDW%Value(1:NIntC_VDW)
        IntCs%Constraint(ILast+1:ILast+NIntC_VDW)=IntC_VDW%Constraint(1:NIntC_VDW)
!
       IF(NIntC_Extra/=0) THEN
!
          ILast=ILast+NIntC_VDW
       IntCs%Def(ILast+1:ILast+NIntC_Extra)=IntC_Extra%Def(1:NIntC_Extra)
     IntCs%Atoms(ILast+1:ILast+NIntC_Extra,1:4)=IntC_Extra%Atoms(1:NIntC_Extra,1:4)
     IntCs%Value(ILast+1:ILast+NIntC_Extra)=IntC_Extra%Value(1:NIntC_Extra)
IntCs%Constraint(ILast+1:ILast+NIntC_Extra)=IntC_Extra%Constraint(1:NIntC_Extra)
!
       ENDIF
!
    write(*,*) 'merged intcs'
    do i=1,NIntC
    write(*,110) i,IntCs%Def(i),IntCs%Atoms(i,1:4)
    enddo
!
! tidy up
!
        CALL Delete(IntC_Cov)
        CALL Delete(IntC_VDW)
        IF(NIntC_Extra/=0) CALL Delete(IntC_Extra)
!
! Now filter out repeated definitions, 
! which may have occured in INTC_Extra
!
        ILast=NIntC_Cov+NIntC_VDW
100     II=0
        DO I=ILast+1,ILast+NIntC_Extra
          DO J=1,ILast
            IF(IntCs%Atoms(J,1)==IntCs%Atoms(I,1).AND.&
               IntCs%Atoms(J,2)==IntCs%Atoms(I,2).AND.&
               IntCs%Atoms(J,3)==IntCs%Atoms(I,3).AND.&
               IntCs%Atoms(J,4)==IntCs%Atoms(I,4)) THEN
               IntCs%Def(I)='BLANK'
               IntCs%Atoms(I,1:4)=0      
            ENDIF
          ENDDO
        ENDDO
        IF(II/=0) GO TO 100
!
! Compress IntCs array, get rid of BLANK-s
!
        IF(ANY(IntCs%Def(:)=='BLANK')) THEN
          NNew=NIntC
          DO I=1,NIntc
            IF(IntCs%Def(I)=='BLANK') NNew=NNew-1
          ENDDO
          CALL New(IntC_New,NNew)
	  NNew=0    
          DO I=1,NIntc
            IF(IntCs%Def(I)/='BLANK') THEN
              NNew=NNew+1
              IntC_New%Def(NNew)=IntCs%Def(I)
              IntC_New%Atoms(NNew,1:4)=IntCs%Atoms(I,1:4)
              IntC_New%Value(NNew)=IntCs%Value(I)
              IntC_New%Constraint(NNew)=IntCs%Constraint(I)
            ENDIF
          ENDDO
          CALL Delete(IntCs)
          NIntC=NNew
          CALL New(IntCs,NIntC)
          IntCs=IntC_New
          CALL Delete(IntC_New)
        ENDIF
!
      CALL Delete(MMAtNum)
    write(*,*) 'filtered intcs'
    do i=1,NIntC
    write(*,110) i,IntCs%Def(i),IntCs%Atoms(i,1:4)
    enddo
!
END SUBROUTINE GetIntCs   
!
!-------------------------------------------------------
!
    SUBROUTINE CoordTrf(GMLoc,Refresh,VectCart,VectInt,TrfTyp)
!
! Routine to carry out coordinate transformations
! Refresh controls refreshing of internal coordinate set
! In case of Int-> Cartesian refresh is not allowed,
! any refresh request will be overwritten
! TrfTyp=1 : Cartesian -> Internal
! TrfTyp=2 : Internal -> Cartesian
!
    TYPE(CRDS) :: GMLoc
    TYPE(DBL_VECT) :: VectCart,VectInt
    TYPE(BCSR) :: SpVectCart,SpVectInt,SpB
    INTEGER :: NCart,I,J,Refresh,TrfTyp,NIntC
    INTEGER :: NVBlocksB,NHBlocksB
    TYPE(INTC) :: IntCs
    TYPE(BMATR):: B
!
! Check refresh
!
    CALL OpenASCII(OutFile,Out)
    IF(TrfTyp==2.AND.Refresh/=3) THEN
      Refresh=3
      WRITE(Out,*) 'Internal Coordinate refresh refused in Internal -> Cartesian transformation'
    ENDIF
    CLOSE(Out)
!
    NCart=3*GMLoc%Natms
!
! Now, generate/refresh internal coordinates
! also put them into HDF
!
    CALL GetIntCs(GMLoc%Carts%D,GMLoc%Natms,InfFile,IntCs,NIntC,Refresh)
!
! Generate vibrational B matrix in quasi-sparse TYPE(BMATR) representation 
!
      GMLoc%Carts%D=GMLoc%Carts%D/AngstromsToAU 
    CALL BMatrix(GMLoc%Natms,GMLoc%Carts%D,NIntC,IntCs,B)
      GMLoc%Carts%D=GMLoc%Carts%D*AngstromsToAU 
write(*,*) 'aft bmatrix'
    do i=1,NIntC
    write(*,100) i,IntCs%Def(i),IntCs%Atoms(i,1:4),B%IB(i,1:12),B%B(i,1:12)
    enddo
100 format(i3,' ',A5,4(I3,' '),/,12(I3,' '),12F8.3)
!
! Now, turn B matrix into sparse blocked representation
!
! first, determine sparse matr. params.
    NVBlocksB=(MAX(NIntC,GMLoc%Natms)-1)/3+1 !!! vertical matrix dimension in blocks of size 3x3
    NHBlocksB=(NCart-1)/3+1 !!! horizontal matrix dimension in blocks of size 3x3
 write(*,*) 'dimensions of B= ',NVBlocksB,NHBlocksB
    Natoms=NVBlocksB !!!! temporary
    MaxAtms=Natoms+1
    MaxBlks=NVBlocksB*NHBlocksB
    MaxNon0=9*MaxBlks !!! including zeros of blocks
    CALL New(BSiz,Natoms)
    CALL New(OffS,Natoms)
    BSiz%I=3
    OffS%I(1)=1
    DO I=2,Natoms
      OffS%I(I)=OffS%I(I-1)+3
    ENDDO
!
    CALL Set_BCSR_EQ_BMATR(SpB,B)
!
! flag for pprint
    NBasF=3*NVBlocksB
    PrintFlags%Mat=DEBUG_MATRICES
    CALL PPrint(SpB,'SpB',Unit_O=6)
!
! Now, turn input vector into sparse matrix form
!
!         CALL PPrint(VectCart,'VectCart',Unit_O=6)
    CALL Set_BCSR_EQ_VECT(SpVectCart,VectCart)
          CALL PPrint(SpVectCart,'SpVectCart',Unit_O=6)
!
    CALL Delete(SpB)
    CALL Delete(B)
    CALL Delete(IntCs)
    CALL Delete(SpVectCart)
    CALL Delete(BSiz)
    CALL Delete(OffS)
!
END SUBROUTINE CoordTrf
#endif
!-------------------------------------------------------
END MODULE IntCoo

