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
SUBROUTINE Topology_12(NATOMS,NBONDS,BONDI,BONDJ,Top12,InfFile)
! Set up a table which shows the atom numbers of atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2) :: Top12_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBONDS,NATOMS,NMax12
INTEGER,DIMENSION(1:NBONDS) :: BONDI,BONDJ
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    NMax12=5
    CALL New(Top12,(/NATOMS,NMax12+1/))
    Top12%I(:,:)=0 
!
    DO I=1,NBONDS
!
      II=BONDI(I)
      JJ=BONDJ(I)
      NI=Top12%I(II,1)
      NJ=Top12%I(JJ,1)
!
! check matrix size, increase size if necessary
      IF(NI>=NMax12 .OR. NJ>=NMax12) THEN
        NMax12=NMax12+5
        CALL New(Top12_2,(/NATOMS,NMax12+1/))
        Top12_2%I(1:NATOMS,1:NMax12+1)=0 
        Top12_2%I(1:NATOMS,1:NMax12+1-5)=Top12%I(1:NATOMS,1:NMax12+1-5)
        CALL DELETE(Top12)
        CALL New(Top12,(/NATOMS,NMax12+1/))
        Top12%I(:,:)=Top12_2%I(:,:)
        CALL DELETE(Top12_2)
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
    IF(.NOT.PRESENT(Top12)) CALL DELETE(Top12)
!
END SUBROUTINE Topology_12 
!--------------------------------------------------------------
!
SUBROUTINE Topology_13(NATOMS,Top12,Top13,InfFile)
! Set up a table which shows the atom numbers of atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2),OPTIONAL :: Top13
TYPE(INT_RNK2) :: Top13_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NATOMS,NMax13,NMax12,KK,IN12,JN12
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL OpenHDF(InfFile)
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NATOMS,K/))
        CALL Get(Top12,'Top12')
        CALL CloseHDF()
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix')
      ENDIF
    ELSE
      NMax12=SIZE(Top12%I,2)
    ENDIF
!
    NMax13=10
    K=NMax13+1
    CALL New(Top13,(/NATOMS,K/))
    Top13%I(1:NATOMS,1:NMax13+1)=0 
!
    DO II=1,NATOMS
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
! check matrix size, increase size if necessary
!
      IF(NI>=NMax13) THEN
        NMax13=NMax13+10
        CALL New(Top13_2,(/NATOMS,NMax13+1/))
        Top13_2%I(1:NATOMS,1:NMax13+1)=0 
        Top13_2%I(1:NATOMS,1:NMax13+1-10)=Top13%I(1:NATOMS,1:NMax13+1-10)
        CALL DELETE(Top13)
        CALL New(Top13,(/NATOMS,NMax13+1/))
        Top13%I(1:NATOMS,1:NMax13+1)=Top13_2%I(1:NATOMS,1:NMax13+1)
        CALL DELETE(Top13_2)
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
    IF(.NOT.PRESENT(Top12)) CALL DELETE(Top12)
!
    IF(PRESENT(InfFile)) THEN
      CALL OpenHDF(InfFile)
      CALL Put(NMax13,'NMax13')
      CALL Put(Top13,'Top13')
      CALL CloseHDF()
    ENDIF
!
    IF(.NOT.PRESENT(Top13)) CALL DELETE(Top13)
!
END SUBROUTINE Topology_13 
!--------------------------------------------------------------
!
SUBROUTINE Topology_14(NATOMS,Top12,Top14,InfFile)
! Set up a table which shows the atom numbers of atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL:: Top12
TYPE(INT_RNK2),OPTIONAL:: Top14
TYPE(INT_RNK2) :: Top14_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
INTEGER :: NATOMS,NMax14,NMax12,IN12,JN12,KN12
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
     IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL OpenHDF(InfFile)
        CALL Get(NMax12,'NMAX12')
        K=NMax12+1
        CALL New(Top12,(/NATOMS,K/))
        CALL Get(Top12,'Top12')
        CALL CloseHDF()
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix')
      ENDIF
    ELSE
      NMax12=SIZE(Top12%I,2)
    ENDIF
!
    NMax14=10
    K=NMax14+1
    CALL New(Top14,(/NATOMS,K/))
    Top14%I(1:NATOMS,1:NMax14+1)=0 
!
    DO II=1,NATOMS
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
! check matrix size, increase size if necessary
!
      IF(NI>=NMax14) THEN
        NMax14=NMax14+10
        CALL New(Top14_2,(/NATOMS,NMax14+1/))
        Top14_2%I(1:NATOMS,1:NMax14+1)=0 
        Top14_2%I(1:NATOMS,1:NMax14+1-10)=Top14%I(1:NATOMS,1:NMax14+1-10)
        CALL DELETE(Top14)
        CALL New(Top14,(/NATOMS,NMax14+1/))
        Top14%I(1:NATOMS,1:NMax14+1)=Top14_2%I(1:NATOMS,1:NMax14+1)
        CALL DELETE(Top14_2)
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
    IF(.NOT.PRESENT(Top12)) CALL DELETE(Top12)
    IF(.NOT.PRESENT(Top14)) CALL DELETE(Top14)
!
END SUBROUTINE Topology_14 
!--------------------------------------------------------------
!
SUBROUTINE SORT_INTO_BOX1(BOXSIZE,C,NATOMS,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
! sort the atoms of a molecule into boxes
!
! BOXSIZE: linear box size
!
IMPLICIT NONE
REAL(DOUBLE) :: BOXSIZE,VBIG,C(1:3,NATOMS),BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
INTEGER :: I,J,JJ,NX,NY,NZ,NBOX,IX,IY,IZ,IORD,IADD,NATOMS
SAVE VBIG
DATA VBIG/1.D+90/ 
!
! First count the number of atoms in the individual boxes 
!
!find borders of the global box
!
    BXMIN= VBIG
    BXMAX=-VBIG
    BYMIN= VBIG
    BYMAX=-VBIG
    BZMIN= VBIG
    BZMAX=-VBIG
  DO I=1,NATOMS
    IF(C(1,I)<BXMIN) BXMIN=C(1,I)
    IF(C(1,I)>BXMAX) BXMAX=C(1,I)
    IF(C(2,I)<BYMIN) BYMIN=C(2,I)
    IF(C(2,I)>BYMAX) BYMAX=C(2,I)
    IF(C(3,I)<BZMIN) BZMIN=C(3,I)
    IF(C(3,I)>BZMAX) BZMAX=C(3,I)
  ENDDO
!
  NX=INT((BXMAX-BXMIN)/BOXSIZE)+1
  NY=INT((BYMAX-BYMIN)/BOXSIZE)+1
  NZ=INT((BZMAX-BZMIN)/BOXSIZE)+1
  NBOX=NX*NY*NZ
!
END SUBROUTINE SORT_INTO_BOX1
!
!--------------------------------------------------------------
SUBROUTINE SORT_INTO_BOX2(BOXSIZE,C,NATOMS,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BOXI1,BOXJ1,InfFile,ISet)
!
! sort the atoms of a molecule into boxes
!
! BOXI(I) : contains the ordering number of the first atom of the I-th box (like in sparse row-wise)
! BOXJ(J) : gives the original serial number of the atom desribed by the J-th ordering number
! C: contains Cartesian coordinates of atoms
! BOXSIZE: linear box size
!
IMPLICIT NONE
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER,OPTIONAL :: ISET 
REAL(DOUBLE) :: BOXSIZE,VBIG,C(1:3,NATOMS),BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
TYPE(INT_VECT),OPTIONAL :: BOXI1,BOXJ1
TYPE(INT_VECT) :: BOXI,BOXJ
INTEGER,ALLOCATABLE,DIMENSION(:) :: ISIGN
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: BOXCOUNTER 
INTEGER :: I,J,JJ,NX,NY,NZ,NBOX,IX,IY,IZ,IORD,IADD,NATOMS
SAVE VBIG
DATA VBIG/1.D+90/ 
!
  NBOX=NX*NY*NZ
  CALL New(BOXI,NBOX+1)
  CALL New(BOXJ,NATOMS)
!
  ALLOCATE(ISIGN(1:NATOMS))
!
  ALLOCATE(BOXCOUNTER(1:NX,1:NY,1:NZ))
  BOXCOUNTER(1:NX,1:NY,1:NZ)=0
  BOXI%I(1:NBOX+1)=0
!
! COUNT NUMBER OF ATOMS IN THE BOX
!
  DO I=1,NATOMS
!
! identify box
!     
    IX=INT((C(1,I)-BXMIN)/BOXSIZE)+1
    IY=INT((C(2,I)-BYMIN)/BOXSIZE)+1
    IZ=INT((C(3,I)-BZMIN)/BOXSIZE)+1
!     
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY !order parameter: ZXY
    BOXCOUNTER(IX,IY,IZ)=BOXCOUNTER(IX,IY,IZ)+1
    ISIGN(I)=IORD*(NATOMS+1)+BOXCOUNTER(IX,IY,IZ) !!! shows both box and index within box
!     
  ENDDO
!     
! Count pointers to individual boxes within array A
!     
     BOXI%I(1)=1
  DO IZ=1,NZ
  DO IX=1,NX
  DO IY=1,NY
!
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY
    BOXI%I(IORD+1)=BOXI%I(IORD)+BOXCOUNTER(IX,IY,IZ)
!     
  ENDDO
  ENDDO
  ENDDO
!
! Set up contents of boxes as represented in A, sparse row-wise
!
  DO I=1,NATOMS
!
    IORD=INT(ISIGN(I)/(NATOMS+1))
    IADD=ISIGN(I)-IORD*(NATOMS+1)
    BOXJ%I(BOXI%I(IORD)-1+IADD)=I
!     
  ENDDO
!
! Now check ordering algorithm
!
! DO I=1,NBOX  
! DO J=BOXI%I(I),BOXI%I(I+1)-1 !!! j is ordering index
!   JJ=BOXJ%I(J) !!! jj is original atom number
!   write(*,100) i,C(1:3,JJ)
! ENDDO
! ENDDO
100 format(I8,3F12.8)
!
  IF(PRESENT(BOXI1)) THEN
    BOXI1%I(:)=BOXI%I(:)
  ENDIF
  IF(PRESENT(BOXJ1)) THEN
    BOXJ1%I(:)=BOXJ%I(:)
  ENDIF
!
  IF(PRESENT(InfFile)) THEN
    CALL OpenHDF(InfFile)
    CALL Put(NBOX,'NBOX'//TRIM(IntToChar(ISet)))
    CALL Put(BOXI,'BOXI'//TRIM(IntToChar(ISet)))
    CALL Put(BOXJ,'BOXJ'//TRIM(IntToChar(ISet)))
    CALL CloseHDF()
  ENDIF
!
  DEALLOCATE(BOXCOUNTER)
  CALL DELETE(BOXI)
  CALL DELETE(BOXJ)
!
END SUBROUTINE SORT_INTO_BOX2
!
!----------------------------------------------------------------
!
SUBROUTINE SORT_INTO_BOX(BOXSIZE,C,NATOMS,InfFile,ISet)
IMPLICIT NONE
INTEGER,OPTIONAL :: ISet
INTEGER :: NATOMS,NX,NY,NZ,NBOX
REAL(DOUBLE) :: BOXSIZE,BXMIN,BYMIN,BZMIN
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
TYPE(INT_VECT) :: BOXI1,BOXJ1
REAL(DOUBLE),DIMENSION(1:3,1:NATOMS) :: C
!
CALL SORT_INTO_BOX1(BOXSIZE,C,NATOMS,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
CALL New(BOXI1,NBOX+1)
CALL New(BOXJ1,NATOMS)
!
CALL SORT_INTO_BOX2(BOXSIZE,C,NATOMS,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BOXI1,BOXJ1,InfFile,ISet)
!
CALL DELETE(BOXI1)
CALL DELETE(BOXJ1)
!
END SUBROUTINE SORT_INTO_BOX
!----------------------------------------------------------------
!
SUBROUTINE Topologies_MM(NATOMS,NBONDS,BONDI,BONDJ,InfFile,Top12OUT)
! Set up a table which shows the atom numbers of atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12OUT
TYPE(INT_RNK2) :: Top12OUT_2
TYPE(INT_RNK2) :: Top12,Top13,Top14
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBONDS,NATOMS
INTEGER,DIMENSION(1:NBONDS) :: BONDI,BONDJ
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
   CALL Topology_12(NATOMS,NBONDS,BONDI,BONDJ,Top12,InfFile)
   CALL Topology_13(NATOMS,Top12,Top13,InfFile)
   CALL Topology_14(NATOMS,Top12,Top14,InfFile)
   CALL Excl_List(Natoms,InfFile,Top12,Top13,Top14)
   CALL Excl_List14(Natoms,InfFile,Top12,Top13,Top14)
!
   IF(PRESENT(Top12OUT)) THEN
     IF(AllocQ(Top12OUT%Alloc)) CALL DELETE(Top12OUT)
     N=SIZE(Top12%I,2)
     CALL New(Top12OUT_2,(/NATOMS,N/))
     Top12OUT_2=Top12
     CALL DELETE(Top12)
     CALL DELETE(Top13)
     CALL DELETE(Top14)
     CALL New(Top12OUT,(/NATOMS,N/))
     Top12OUT=Top12OUT_2
     CALL DELETE(Top12OUT_2)
   ELSE
     CALL DELETE(Top12)
     CALL DELETE(Top13)
     CALL DELETE(Top14)
   ENDIF
!
END SUBROUTINE Topologies_MM
!----------------------------------------------------------------
!
SUBROUTINE Excl_List(Natoms,InfFile,Top12,Top13,Top14,Top_Excl_Out)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER :: Natoms,NMax12,NMax13,NMax14,NMax_Excl,NNEW,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ
!
      IF(PRESENT(InfFile)) CALL OpenHDF(InfFile)
!
    IF(PRESENT(Top12)) THEN
      NMax12=SIZE(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/Natoms,NMax12+1/))
        CALL Get(Top12,'Top12')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top13)) THEN
      NMax13=SIZE(Top13%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax13,'NMax13')
        CALL New(Top13,(/Natoms,NMax13+1/))
        CALL Get(Top13,'Top13')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_13 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top14)) THEN
      NMax14=SIZE(Top14%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax14,'NMax14')
        CALL New(Top14,(/Natoms,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl
!
    NMax_Excl=NMax12+NMax13+NMax14
    CALL New(Top_Excl,(/Natoms,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
    Top_Excl%I(1:Natoms,1:NMAX12+1)=Top12%I(1:Natoms,1:NMAX12+1)
!
! Now merge Topologies, in order to avoid double counting in 
! exclusion energies
!
      NMAX_Excl_Out=0
    DO I=1,Natoms
!
      NNEW=0
      NOLD=Top_Excl%I(I,1)
      DO J=1,Top13%I(I,1)
        JJ=Top13%I(I,J+1)
        IF(ANY(Top_Excl%I(I,2:NOLD+1)==JJ)) THEN
          CYCLE
        ELSE
          NNEW=NNEW+1
          Top_Excl%I(I,NOLD+1+NNEW)=JJ
        ENDIF
      ENDDO 
      NNEW=NOLD+NNEW
      Top_Excl%I(I,1)=NNEW
      IF(NMAX_Excl_Out<NNEW) NMAX_Excl_Out=NNEW
!
      NNEW=0
      NOLD=Top_Excl%I(I,1)
      DO J=1,Top14%I(I,1)
        JJ=Top14%I(I,J+1)
        IF(ANY(Top_Excl%I(I,2:NOLD+1)==JJ)) THEN
          CYCLE
        ELSE
          NNEW=NNEW+1
          Top_Excl%I(I,NOLD+1+NNEW)=JJ
        ENDIF
      ENDDO 
      NNEW=NOLD+NNEW
      Top_Excl%I(I,1)=NNEW
      IF(NMAX_Excl_Out<NNEW) NMAX_Excl_Out=NNEW
!
    ENDDO 
!
! STORE RESULT
!
    IF(PRESENT(Top_Excl_Out)) THEN
      CALL New(Top_Excl_Out,(/Natoms,NMAX_Excl_Out+1/))
      Top_Excl_Out%I(1:Natoms,1:NMAX_Excl_Out+1)=Top_Excl%I(1:Natoms,1:NMAX_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMAX_Excl_Out,'NMAX_EXCL')
        CALL Put(Top_Excl_Out,'TOP_EXCL')
      ENDIF
    ELSE
      CALL New(Top_New,(/Natoms,NMAX_Excl_Out+1/))
      Top_New%I(1:Natoms,1:NMAX_Excl_Out+1)=Top_Excl%I(1:Natoms,1:NMAX_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMAX_Excl_Out,'NMAX_EXCL')
        CALL Put(Top_New,'TOP_EXCL')
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
SUBROUTINE Excl_List14(Natoms,InfFile,Top12,Top13,Top14,Top_Excl_Out)
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
INTEGER :: Natoms,NMax12,NMax13,NMax14,NMax_Excl,NNEW,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ,NEXCL,N12,N13
!
      IF(PRESENT(InfFile)) CALL OpenHDF(InfFile)
!
    IF(PRESENT(Top12)) THEN
      NMax12=SIZE(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/Natoms,NMax12+1/))
        CALL Get(Top12,'Top12')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_12 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top13)) THEN
      NMax13=SIZE(Top13%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax13,'NMax13')
        CALL New(Top13,(/Natoms,NMax13+1/))
        CALL Get(Top13,'Top13')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_13 matrix in Excl_list')
      ENDIF
    ENDIF
!
    IF(PRESENT(Top14)) THEN
      NMax14=SIZE(Top14%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax14,'NMax14')
        CALL New(Top14,(/Natoms,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl to the size of Top14 and to zero 
!
    NMax_Excl=NMax14
    CALL New(Top_Excl,(/Natoms,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
!
! Now merge Topologies, in order to avoid double counting in 
! exclusion energies
!
      NMAX_Excl_Out=0
    DO I=1,Natoms
!
      DO J=1,Top14%I(I,1)
          NEXCL=Top_Excl%I(I,1)
          N12=Top12%I(I,1)
          N13=Top13%I(I,1)
          JJ=Top14%I(I,J+1)
        IF(ANY(Top_Excl%I(I,2:NEXCL+1)==JJ).OR. &
           ANY(Top12%I(I,2:N12+1)==JJ).OR. & 
           ANY(Top13%I(I,2:N13+1)==JJ)  ) THEN
          CYCLE
        ELSE
          Top_Excl%I(I,1)=NEXCL+1
          Top_Excl%I(I,NEXCL+2)=JJ
        ENDIF
      ENDDO 
          NEXCL=Top_Excl%I(I,1)
      IF(NMAX_Excl_Out<NEXCL) NMAX_Excl_Out=NEXCL
!
    ENDDO 
!
! STORE RESULT
!
    IF(PRESENT(Top_Excl_Out)) THEN
      CALL New(Top_Excl_Out,(/Natoms,NMAX_Excl_Out+1/))
      Top_Excl_Out%I(1:Natoms,1:NMAX_Excl_Out+1)=Top_Excl%I(1:Natoms,1:NMAX_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMAX_Excl_Out,'NMAX_EXCL14')
        CALL Put(Top_Excl_Out,'TOP_EXCL14')
      ENDIF
    ELSE
      CALL New(Top_New,(/Natoms,NMAX_Excl_Out+1/))
      Top_New%I(1:Natoms,1:NMAX_Excl_Out+1)=Top_Excl%I(1:Natoms,1:NMAX_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        CALL Put(NMAX_Excl_Out,'NMAX_EXCL14')
        CALL Put(Top_New,'TOP_EXCL14')
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
