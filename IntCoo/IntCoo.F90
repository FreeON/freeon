MODULE IntCoo
!
   USE DerivedTypes
   USE GlobalScalars
   USE InOut
   Use MemMan
   Use ProcessControl
!
IMPLICIT NONE
!
CONTAINS
!--------------------------------------------------------------
!
SUBROUTINE Topology_12(NatomsMM,NBonds,BondI,BondJ,Top12,InfFile)
! Set up a table which shows the atom numbers of atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2) :: Top12_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBonds,NatomsMM,NMax12
INTEGER,DIMENSION(1:NBonds) :: BondI,BondJ
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    NMax12=5
    CALL New(Top12,(/NatomsMM,NMax12+1/))
    Top12%I(:,:)=0 
!
    DO I=1,NBonds
!
      II=BondI(I)
      JJ=BondJ(I)
      NI=Top12%I(II,1)
      NJ=Top12%I(JJ,1)
!
! check matrix Size, increase Size if necessary
      IF(NI>=NMax12 .OR. NJ>=NMax12) THEN
        NMax12=NMax12+5
        CALL New(Top12_2,(/NatomsMM,NMax12+1/))
        Top12_2%I(1:NatomsMM,1:NMax12+1)=0 
        Top12_2%I(1:NatomsMM,1:NMax12+1-5)=Top12%I(1:NatomsMM,1:NMax12+1-5)
        CALL Delete(Top12)
        CALL New(Top12,(/NatomsMM,NMax12+1/))
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
      CALL OpenHDF(InfFile)
      CALL Put(NMax12,'NMax12')
      CALL Put(Top12,'Top12')
      CALL CloseHDF()
    ENDIF
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
!
END SUBROUTINE Topology_12 
!--------------------------------------------------------------
!
SUBROUTINE Topology_13(NatomsMM,Top12,Top13,InfFile)
! Set up a table which shows the atom numbers of atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2),OPTIONAL :: Top13
TYPE(INT_RNK2) :: Top13_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NatomsMM,NMax13,NMax12,KK,IN12,JN12
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL OpenHDF(InfFile)
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NatomsMM,K/))
        CALL Get(Top12,'Top12')
        CALL CloseHDF()
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix')
      ENDIF
    ELSE
      NMax12=Size(Top12%I,2)
    ENDIF
!
    NMax13=10
    K=NMax13+1
    CALL New(Top13,(/NatomsMM,K/))
    Top13%I(1:NatomsMM,1:NMax13+1)=0 
!
    DO II=1,NatomsMM
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
        CALL New(Top13_2,(/NatomsMM,NMax13+1/))
        Top13_2%I(1:NatomsMM,1:NMax13+1)=0 
        Top13_2%I(1:NatomsMM,1:NMax13+1-10)=Top13%I(1:NatomsMM,1:NMax13+1-10)
        CALL Delete(Top13)
        CALL New(Top13,(/NatomsMM,NMax13+1/))
        Top13%I(1:NatomsMM,1:NMax13+1)=Top13_2%I(1:NatomsMM,1:NMax13+1)
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
      CALL OpenHDF(InfFile)
      CALL Put(NMax13,'NMax13')
      CALL Put(Top13,'Top13')
      CALL CloseHDF()
    ENDIF
!
    IF(.NOT.PRESENT(Top13)) CALL Delete(Top13)
!
END SUBROUTINE Topology_13 
!--------------------------------------------------------------
!
SUBROUTINE Topology_14(NatomsMM,Top12,Top14,InfFile)
! Set up a table which shows the atom numbers of atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL:: Top12
TYPE(INT_RNK2),OPTIONAL:: Top14
TYPE(INT_RNK2) :: Top14_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
INTEGER :: NatomsMM,NMax14,NMax12,IN12,JN12,KN12
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
     IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL OpenHDF(InfFile)
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NatomsMM,K/))
        CALL Get(Top12,'Top12')
        CALL CloseHDF()
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix')
      ENDIF
    ELSE
      NMax12=Size(Top12%I,2)
    ENDIF
!
    NMax14=10
    K=NMax14+1
    CALL New(Top14,(/NatomsMM,K/))
    Top14%I(1:NatomsMM,1:NMax14+1)=0 
!
    DO II=1,NatomsMM
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
        CALL New(Top14_2,(/NatomsMM,NMax14+1/))
        Top14_2%I(1:NatomsMM,1:NMax14+1)=0 
        Top14_2%I(1:NatomsMM,1:NMax14+1-10)=Top14%I(1:NatomsMM,1:NMax14+1-10)
        CALL Delete(Top14)
        CALL New(Top14,(/NatomsMM,NMax14+1/))
        Top14%I(1:NatomsMM,1:NMax14+1)=Top14_2%I(1:NatomsMM,1:NMax14+1)
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
      CALL OpenHDF(InfFile)
      CALL Put(NMax14,'NMax14')
      CALL Put(Top14,'Top14')
      CALL CloseHDF()
    ENDIF
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
    IF(.NOT.PRESENT(Top14)) CALL Delete(Top14)
!
END SUBROUTINE Topology_14 
!--------------------------------------------------------------
!
SUBROUTINE SORT_INTO_Box1(BoxSize,C,NatomsMM,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
! sort the atoms of a molecule into Boxes
!
! BoxSize: linear Box Size
!
IMPLICIT NONE
REAL(DOUBLE) :: BoxSize,VBIG,C(1:3,NatomsMM),BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NatomsMM
SAVE VBIG
DATA VBIG/1.D+90/ 
!
! First count the number of atoms in the individual Boxes 
!
!find borders of the global Box
!
    BXMIN= VBIG
    BXMax=-VBIG
    BYMIN= VBIG
    BYMax=-VBIG
    BZMIN= VBIG
    BZMax=-VBIG
  DO I=1,NatomsMM
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
SUBROUTINE SORT_INTO_Box2(BoxSize,C,NatomsMM,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,InfFile,ISet)
!
! sort the atoms of a molecule into Boxes
!
! BoxI(I) : contains the ordering number of the first atom of the I-th Box (like in sparse row-wise)
! BoxJ(J) : gives the original serial number of the atom desribed by the J-th ordering number
! C: contains Cartesian coordinates of atoms
! BoxSize: linear Box Size
!
IMPLICIT NONE
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER,OPTIONAL :: ISET 
REAL(DOUBLE) :: BoxSize,VBIG,C(1:3,NatomsMM),BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
TYPE(INT_VECT),OPTIONAL :: BoxI1,BoxJ1
TYPE(INT_VECT) :: BoxI,BoxJ
INTEGER,ALLOCATABLE,DIMENSION(:) :: ISIGN
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: BoxCOUNTER 
INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NatomsMM
SAVE VBIG
DATA VBIG/1.D+90/ 
!
  NBox=NX*NY*NZ
  CALL New(BoxI,NBox+1)
  CALL New(BoxJ,NatomsMM)
!
  ALLOCATE(ISIGN(1:NatomsMM))
!
  ALLOCATE(BoxCOUNTER(1:NX,1:NY,1:NZ))
  BoxCOUNTER(1:NX,1:NY,1:NZ)=0
  BoxI%I(1:NBox+1)=0
!
! COUNT NUMBER OF ATOMS IN THE Box
!
  DO I=1,NatomsMM
!
! identify Box
!     
    IX=INT((C(1,I)-BXMIN)/BoxSize)+1
    IY=INT((C(2,I)-BYMIN)/BoxSize)+1
    IZ=INT((C(3,I)-BZMIN)/BoxSize)+1
!     
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY !order parameter: ZXY
    BoxCOUNTER(IX,IY,IZ)=BoxCOUNTER(IX,IY,IZ)+1
    ISIGN(I)=IORD*(NatomsMM+1)+BoxCOUNTER(IX,IY,IZ) !!! shows both Box and index within Box
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
  DO I=1,NatomsMM
!
    IORD=INT(ISIGN(I)/(NatomsMM+1))
    IADD=ISIGN(I)-IORD*(NatomsMM+1)
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
    CALL OpenHDF(InfFile)
    CALL Put(NBox,'NBox'//TRIM(IntToChar(ISet)))
    CALL Put(BoxI,'BoxI'//TRIM(IntToChar(ISet)))
    CALL Put(BoxJ,'BoxJ'//TRIM(IntToChar(ISet)))
    CALL CloseHDF()
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
SUBROUTINE SORT_INTO_Box(BoxSize,C,NatomsMM,InfFile,ISet)
IMPLICIT NONE
INTEGER,OPTIONAL :: ISet
INTEGER :: NatomsMM,NX,NY,NZ,NBox
REAL(DOUBLE) :: BoxSize,BXMIN,BYMIN,BZMIN
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
TYPE(INT_VECT) :: BoxI1,BoxJ1
REAL(DOUBLE),DIMENSION(1:3,1:NatomsMM) :: C
!
CALL SORT_INTO_Box1(BoxSize,C,NatomsMM,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
NBox=NX*NY*NZ
CALL New(BoxI1,NBox+1)
CALL New(BoxJ1,NatomsMM)
!
CALL SORT_INTO_Box2(BoxSize,C,NatomsMM,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,InfFile,ISet)
!
CALL Delete(BoxI1)
CALL Delete(BoxJ1)
!
END SUBROUTINE SORT_INTO_Box
!----------------------------------------------------------------
!
SUBROUTINE Topologies_MM(NatomsMM,NBonds,BondI,BondJ,InfFile,Top12OUT)
! Set up a table which shows the atom numbers of atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12OUT
TYPE(INT_RNK2) :: Top12OUT_2
TYPE(INT_RNK2) :: Top12,Top13,Top14
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBonds,NatomsMM
INTEGER,DIMENSION(1:NBonds) :: BondI,BondJ
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
   CALL Topology_12(NatomsMM,NBonds,BondI,BondJ,Top12,InfFile)
   CALL Topology_13(NatomsMM,Top12,Top13,InfFile)
   CALL Topology_14(NatomsMM,Top12,Top14,InfFile)
   CALL Excl_List(NatomsMM,InfFile,Top12,Top13,Top14)
   CALL Excl_List14(NatomsMM,InfFile,Top12,Top13,Top14)
!
   IF(PRESENT(Top12OUT)) THEN
     IF(AllocQ(Top12OUT%Alloc)) CALL Delete(Top12OUT)
     N=Size(Top12%I,2)
     CALL New(Top12OUT_2,(/NatomsMM,N/))
     Top12OUT_2=Top12
     CALL Delete(Top12)
     CALL Delete(Top13)
     CALL Delete(Top14)
     CALL New(Top12OUT,(/NatomsMM,N/))
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
SUBROUTINE Excl_List(NatomsMM,InfFile,Top12,Top13,Top14,Top_Excl_Out)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER :: NatomsMM,NMax12,NMax13,NMax14,NMax_Excl,NNew,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ
!
      IF(PRESENT(InfFile)) CALL OpenHDF(InfFile)
!
    IF(PRESENT(Top12)) THEN
      NMax12=Size(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NatomsMM,NMax12+1/))
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
        CALL New(Top13,(/NatomsMM,NMax13+1/))
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
        CALL New(Top14,(/NatomsMM,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl
!
    NMax_Excl=NMax12+NMax13+NMax14
    CALL New(Top_Excl,(/NatomsMM,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
    Top_Excl%I(1:NatomsMM,1:NMax12+1)=Top12%I(1:NatomsMM,1:NMax12+1)
!
! Now merge Topologies, in order to avoid double counting in 
! Exclusion energies
!
      NMax_Excl_Out=0
    DO I=1,NatomsMM
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
      CALL New(Top_Excl_Out,(/NatomsMM,NMax_Excl_Out+1/))
      Top_Excl_Out%I(1:NatomsMM,1:NMax_Excl_Out+1)=Top_Excl%I(1:NatomsMM,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMax_Excl_Out,'NMax_Excl')
        CALL Put(Top_Excl_Out,'TOP_Excl')
      ENDIF
    ELSE
      CALL New(Top_New,(/NatomsMM,NMax_Excl_Out+1/))
      Top_New%I(1:NatomsMM,1:NMax_Excl_Out+1)=Top_Excl%I(1:NatomsMM,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMax_Excl_Out,'NMax_Excl')
        CALL Put(Top_New,'TOP_Excl')
      ENDIF
    ENDIF
!
    CALL Delete(Top_Excl)
    CALL Delete(Top_New)
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
    IF(.NOT.PRESENT(Top13)) CALL Delete(Top13)
    IF(.NOT.PRESENT(Top14)) CALL Delete(Top14)
!
    IF(PRESENT(InfFile)) CALL CloseHDF()
!
END SUBROUTINE Excl_LIST 
!
!----------------------------------------------------------------
!
SUBROUTINE Excl_List14(NatomsMM,InfFile,Top12,Top13,Top14,Top_Excl_Out)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
! of atoms in 14 distance. From the Top14 list
! Top13 and Top12 occurences must be filtered out
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER :: NatomsMM,NMax12,NMax13,NMax14,NMax_Excl,NNew,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ,NExcl,N12,N13
!
      IF(PRESENT(InfFile)) CALL OpenHDF(InfFile)
!
    IF(PRESENT(Top12)) THEN
      NMax12=Size(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NatomsMM,NMax12+1/))
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
        CALL New(Top13,(/NatomsMM,NMax13+1/))
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
        CALL New(Top14,(/NatomsMM,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl to the Size of Top14 and to zero 
!
    NMax_Excl=NMax14
    CALL New(Top_Excl,(/NatomsMM,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
!
! Now merge Topologies, in order to avoid double counting in 
! Exclusion energies
!
      NMax_Excl_Out=0
    DO I=1,NatomsMM
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
      CALL New(Top_Excl_Out,(/NatomsMM,NMax_Excl_Out+1/))
      Top_Excl_Out%I(1:NatomsMM,1:NMax_Excl_Out+1)=Top_Excl%I(1:NatomsMM,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMax_Excl_Out,'NMax_Excl14')
        CALL Put(Top_Excl_Out,'TOP_Excl14')
      ENDIF
    ELSE
      CALL New(Top_New,(/NatomsMM,NMax_Excl_Out+1/))
      Top_New%I(1:NatomsMM,1:NMax_Excl_Out+1)=Top_Excl%I(1:NatomsMM,1:NMax_Excl_Out+1)
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
    IF(PRESENT(InfFile)) CALL CloseHDF()
!
END SUBROUTINE Excl_LIST14 
!----------------------------------------------------------------
END MODULE IntCoo
