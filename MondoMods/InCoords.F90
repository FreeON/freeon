MODULE InCoords
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
   USE LinAlg
   USE AInv   
!
!  USE GlobalCharacters
!  USE Parse
!  USE Macros
#ifdef NAG
    USE F90_UNIX_ENV
#endif
!
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
! BoXI(I) : contains the ordering number of the first atom of the I-th Box (like in sparse row-wise)
! BoXJ(J) : gives the original serial number of the atom desribed by the J-th ordering number
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
! must be passed in in Angstroems. !!!! Changed, pass in Au-s!!!
! Later use au representation of everything, including B-matrix.
! Linear bendings _must_ always appear in pairs, Defd as LINB1 and LINB2
!
      IMPLICIT NONE
      INTEGER :: NIntC,NIntC2,I,J,K,L,NAtoms_Loc,IntCoo
      REAL(DOUBLE) :: Zero,thresh_B
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
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
        IF(IntCs%Def(IntCoo)(1:4)=='BEND') THEN
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) !!! central atom
          K=IntCs%Atoms(IntCoo,3) 
       CALL BEND(I,J,K,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
! out of plane of i from the plane of jkl, with center at j
        IF(IntCs%Def(IntCoo)(1:4)=='OutP') THEN
          I=IntCs%Atoms(IntCoo,1) !!! end atom
          J=IntCs%Atoms(IntCoo,2) !!! central atom
          K=IntCs%Atoms(IntCoo,3) !!! Def plane
          L=IntCs%Atoms(IntCoo,4) !!! Def plane
       CALL OutP(I,J,K,L,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ(1:3,L),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
! torsion of i-j-k-l
        IF(IntCs%Def(IntCoo)(1:4)=='TORS') THEN
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) 
          K=IntCs%Atoms(IntCoo,3) 
          L=IntCs%Atoms(IntCoo,4) 
       CALL TORS(I,J,K,L,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ(1:3,L),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
!        write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ENDIF
! linear bendig of i-j-k
        IF(IntCs%Def(IntCoo)(1:5)=='LINB1') THEN
          IF(IntCs%Def(IntCoo+1)(1:5)/='LINB2') CALL Halt('LINB Definitions are not paired!')
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) 
          K=IntCs%Atoms(IntCoo,3) 
       CALL LINB(I,J,K,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),B%IB(IntCoo+1,1:12),B%B(IntCoo+1,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo+1,1:12),B%B(IntCoo+1,1:12)
        ENDIF
        IF(IntCs%Def(IntCoo)(1:5)=='LINB2') CYCLE
        IF(IntCs%Def(IntCoo)(1:4)=='CART' ) CALL Halt('CART is not yet coded')
!
      ENDDO !!!! loop over internal coords
100 format(A5,12I3,/,6F12.6,/,6F12.6)
!
      END SUBROUTINE BMatrix
!--------------------------------------------------------------------------------------
!
      SUBROUTINE Stre(i,j,XI,XJ,IBB,BB,thresh_B)
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
         t = XJ(m) - XI(m)
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
!----------------------------------------------------------------
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
      REAL(DOUBLE) :: DIJSq,Tp,t,DJKSq,dxsq,dx,DOtj,DJK,DJI
      REAL(DOUBLE) :: sinj,SMI,SMK,SUM,thresh_B
      INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
   10 DIJSq=0.d0
      DJKSq=0.d0
      DO 20 m=1,3
        Tp=XJ(m)
        t=XI(m)-Tp
        RJI(m)=t
        DIJSq=DIJSq+t*t
        t=xk(m)-Tp
        rjk(m)=t
        DJKSq=DJKSq+t*t
   20 Continue        
      DJI=dsqrt(DIJSq)
      DJK=dsqrt(DJKSq)
      dx=1.d0
      DOtj=0.d0
      DO 30 m=1,3
      t=RJI(m)/DJI
      eji(m)=t
      Tp=rjk(m)/DJK
      ejk(m)=Tp
   30 DOtj=DOtj+t*Tp
      if(dabs(DOtj).gt.0.99995d0)go to 60  !!!! colinearity
      sinj=dsqrt(1.d0-DOtj*DOtj)
      ii=3*(i-1)
      jj=3*(j-1)
      kk=3*(k-1)
      DO 40 m=1,3
      smi=dx*(DOtj*eji(m)-ejk(m))/(DJI*sinj)
      IF(dabs(smi)>thresh_B) THEN
        IBB(m) = ii+m
        BB(m) = smi
      ENDIF
      smk=dx*(DOtj*ejk(m)-eji(m))/(DJK*sinj)
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
        IBB(:)=0
        BB(:)=Zero 
      RETURN
!
 1000 format(' i-j-k is colinear - use linear bend! For this bending the row of the B matrix will be deleted!')
!
      END SUBROUTINE BEND
!
!---------------------------------------------------------------------
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
      REAL(DOUBLE) :: DIJSq,Tp,t,DJKSq,djlsq,dx,DOtj
      REAL(DOUBLE) :: DJK,DJI,djl,cosi,cosk,cosl,Tpp
      REAL(DOUBLE) :: sinsin,sini,dot,sint,cost,tant,cossin,sml
      REAL(DOUBLE) :: sinj,SMI,SMK,SUM,thresh_B
      INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
   10 DIJSq=0.d0
      DJKSq=0.d0
      djlsq=0.d0
      do 20 m=1,3
      Tp=XJ(m)
      t=XI(m)-Tp
      RJI(m)=t
      DIJSq=DIJSq+t*t
      t=xk(m)-Tp
      rjk(m)=t
      DJKSq=DJKSq+t*t
      t=xl(m)-XJ(m)
      rjl(m)=t
      djlsq=djlsq+t*t
   20 CONTINUE      
      DJI=dsqrt(DIJSq)
      DJK=dsqrt(DJKSq)
      djl=dsqrt(djlsq)
      dx=1.d0
      cosi=0.d0
      cosk=0.d0
      cosl=0.d0
      do 30 m=1,3
      t=RJI(m)/DJI
      eji(m)=t
      Tp=rjk(m)/DJK
      ejk(m)=Tp
      Tpp=rjl(m)/djl
      ejl(m)=Tpp
      cosi=cosi+Tp*Tpp
      cosk=cosk+t*Tpp
   30 cosl=cosl+t*Tp
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
      if(dabs(sint).gt.0.99995d0)go to 80 !!!! sint show OutP angle
      cost=dsqrt(1.d0-sint*sint)
      tant=sint/cost
      ii=3*(i-1)
      jj=3*(j-1)
      kk=3*(k-1)
      ll=3*(l-1)
      cossin=cost*sini
      do 50 m=1,3
      t=c1(m)/cossin
      smi=(t-tant*eji(m))/DJI
      IF(dabs(smi)>thresh_B) THEN
        IBB(m) = ii+m
        BB(m) = dx*smi
      ENDIF
      smk=t*(cosi*cosk-cosl)/(sinsin*DJK)
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
      END SUBROUTINE OutP
!
!----------------------------------------------------------------------
!
      SUBROUTINE Tors(I,J,K,L,XI,XJ,XK,XL,IBB,BB,thresh_B)
!...
!...  Calculate the B matrix elements for torsion as Defined by
!...  R.L. Hilderbrandt, J. Mol. Spec., 44 (1972) 599.
!...
!...  i-j-k-l : torsion around j-k bond
!... 
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,XL,rij,rjk,rlk,eij,ejk,elk,cr,sj,sk
      REAL(DOUBLE) :: SMI,SMK,thresh_B,dij,dijsq,dx,DJK,t,dxsq,DJKSq
      REAL(DOUBLE) :: cosj,sin2j,smj,dlksq,dlk,cosk,sin2k,sml 
      INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
      DJKSq = 0.0d0
      dxsq = 0.0d0
!
      do m = 1, 3
         sj(m) = 0.0d0
         sk(m) = 0.0d0
         t = xk(m) - XJ(m)
         rjk(m) = t
         DJKSq = DJKSq + t * t
      end do
!
      DJK = 1.0d0 / sqrt(DJKSq)
      dx = 1.0d0
!
      do m = 1, 3
        ejk(m) = rjk(m) * DJK
      end do
!
      jj = 3 * (j - 1)
      kk = 3 * (k - 1)
!
!...
        dijsq=0.d0
      do 40 m=1,3
        t=XJ(m)-XI(m)
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
        smk=t*cosj*DJK
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
      if(dabs(cosk).gt.0.99995d0) go to 120
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
        smj=t*cosk*DJK
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
! colinear bonds in torsion, delete the corresponding row of the B matrix!
        IBB(:)=0    
        BB(:)=Zero
      return
!
1030 format(' i-j-k or j-k-l is colinear (no torsion possIBBle), this row of the B matrix will be deleted.')
!
      END SUBROUTINE TORS
!
!------------------------------------------------------------------
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
      REAL(DOUBLE) :: t,SUM,DIJSq,DJKSq,djasq,Tp,DJI,DJK,dja,dx
      REAL(DOUBLE) :: dotj,doTp,test,smi,smk
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
! form vector a, perpendicular to the plane of rjk and RJI
        A(1)= RJI(2)*RJK(3)-RJI(3)*RJK(2)
        A(2)=-RJI(1)*RJK(3)+RJI(3)*RJK(1)
        A(3)= RJI(1)*RJK(2)-RJI(2)*RJK(1)
        SUM=A(1)**2+A(2)**2+A(3)**2
        SUM=SQRT(SUM)
        a(1:3)=a(1:3)/SUM
! generate point a
        a(1:3)=XJ(1:3)+a(1:3)        
!
      DIJSq=0.d0
      DJKSq=0.d0
      djasq=0.d0
      do 20 m=1,3
        Tp=XJ(m)
        t=XI(m)-Tp
        RJI(m)=t
        DIJSq=DIJSq+t*t
        t=xk(m)-Tp
        rjk(m)=t
        DJKSq=DJKSq+t*t
        t=a(m)-Tp
        un(m)=t
        djasq=djasq+t*t
   20 CONTINUE
      DJI=dsqrt(DIJSq)
      DJK=dsqrt(DJKSq)
      dja=dsqrt(djasq)
      dx=1.d0
      dotj=0.d0
      doTp=0.d0
      do 30 m=1,3
        t=RJI(m)/DJI
        Tp=rjk(m)/DJK
        ejk(m)=Tp
        dotj=dotj+t*Tp
        Tp=un(m)/dja
        unit(m)=Tp
        doTp=doTp+t*Tp
   30 CONTINUE       
      test=dabs(dotj)-1.d0
      if(dabs(test).gt.0.00005d0)go to 70
      if(dabs(doTp).gt.0.00005d0)go to 80
      ii=3*(i-1)
      jj=3*(j-1)
      kk=3*(k-1)
      do 40 m=1,3
        t=unit(m)
        if(dabs(t)<thresh_B)go to 40
        t=-dx*t
        smi=t/DJI
        IBB1(m) = ii+m
        BB1(m) = smi
        smk=t/DJK
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
        smi=t/DJI
        IBB2(m) = ii+m
        BB2(m) = smi
        smk=t/DJK
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
! They may be STRE, BEND, TORS, OutP and LINB1&LINB2 -s.
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
      INTEGER :: NX,NY,NZ,NBOX,NAngle,NTorsion,NLinb
      REAL(DOUBLE) :: BXMIN,BYMIN,BZMIN,BoxSize,Fact,Value
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
! now define bonding scheme, based on Slater or Van der Waals radii
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
        CALL Delete(Top_Excl)
        CALL Delete(Top14)
        CALL Delete(Top13)
!
      ELSE
!
! Filter out covalent 12, 13 and 14 connections from VWD bondlist
!
        CALL Get(NMax_Excl,'NMax_Excl')
        CALL New(Top_Excl,(/Natoms_Loc,NMax_Excl+1/))
        CALL Get(Top_Excl,'TOP_EXCL')
          CALL VDWFilter(BondIJ,NBond,Top_Excl)
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
        CALL VDWAngleList(NAtoms_Loc,Top12,AngleIJK,NAngle,BondIJ,NBond)
!
        CALL VDWTorsionList(NAtoms_Loc,Top12,TorsionIJKL,NTorsion,BondIJ,NBond)
!
        CALL Delete(Top12)
!
      ENDIF
!
! Calculate number of linear bendings
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
          IntCs%Def(ILast+NLinB+I)='BEND '
          IntCs%Atoms(ILast+NLinB+I,1:3)=AngleIJK%I(1:3,I)
      ENDDO
        ILast=NBond+NAngle
      DO I=1,NTorsion
        IntCs%Def(ILast+I)='TORS '
        IntCs%Atoms(ILast+I,1:4)=TorsionIJKL%I(1:4,I)
      ENDDO
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
    INTEGER :: NAtoms_Loc,I,I1,I2,J,J1,J2,N1,N2,NBond,NTorsion,K1,K2
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
        K1=Top12%I(I1,J1+1)
        K2=Top12%I(I2,J2+1)
        IF(K1/=K2.AND.K1/=I2.AND.K2/=I1) THEN 
          NTorsion=NTorsion+1
          TorsionIJKL%I(1,NTorsion)=K1
          TorsionIJKL%I(2,NTorsion)=I1
          TorsionIJKL%I(3,NTorsion)=I2
          TorsionIJKL%I(4,NTorsion)=K2
        ENDIF
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
! WARNING! In the present version bending -> linear bending transitions are 
! always checked and refreshed
! 
      IMPLICIT NONE
      TYPE(INTC) :: IntCs,IntC_Cov,IntC_VDW,IntC_Extra,IntC_New
      INTEGER :: NIntC,NIntC_Cov,NIntC_VDW,NIntC_Extra,NNew,Nintc_New
      INTEGER :: I,J,K,Refresh,Natoms_Loc,II,ILast
      INTEGER :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
      TYPE(INT_VECT) :: MMAtNum,LinAtom,MarkLinb
      TYPE(INT_RNK2) :: LinBBridge,Top12
      TYPE(CRDS) :: GMLoc
      CHARACTER(LEN=DefAULT_CHR_LEN) :: InfFile
      REAL(DOUBLE),DIMENSION(1:3,1:Natoms_Loc) :: XYZ
      REAL(DOUBLE) :: Value
!
! Get atomnames (numbers) from HDF 
!
      CALL New(MMAtNum,Natoms_Loc)
      CALL Get(MMAtNum,'MMATNUM')
!
      IF(Refresh==1) Then !!! Total refresh
!
!define covalent bonding scheme
        CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,1,IntC_Cov,NIntC_Cov)
!define Van der Waals bonding scheme
        CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,2,IntC_VDW,NIntC_VDW)
!!!!    CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,2,IntC_VDW,NIntC_VDW)
!!!NIntC_VDW=0 !!!! temporary, see same comment also a few lines below
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
!!!!!   CALL DefineIntCoos(NAtoms_Loc,XYZ,MMAtNum,InfFile,2,IntC_VDW,NIntC_VDW)
!!!!NIntC_VDW=0 !!!!! temporary !!!!
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
          CALL Get(NIntC_Cov,'NIntC_Cov')
          CALL New(IntC_Cov,NIntC_Cov)
          CALL Get(IntC_Cov,'IntC_Cov')
          CALL Get(NIntC_VDW,'NIntC_VDW')
          CALL New(IntC_VDW,NIntC_VDW)
          CALL Get(IntC_VDW,'IntC_VDW')
!
      ENDIF
!
! Merge INTC arrays
!
        NIntC=NIntC_Cov+NIntC_VDW+NIntC_Extra
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
! tidy up
!
        IF(NIntC_Cov/=0) CALL Delete(IntC_Cov)
        IF(NIntC_VDW/=0) CALL Delete(IntC_VDW)
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
!
! Now check for bending -> linear bending transitions
!
      CALL New(MarkLinB,NIntC)
!
      MarkLinB%I=0
          NLinB=0
      DO I=1,NIntC 
        IF(IntCs%Def(I)(1:4)=='BEND') THEN
          I1=IntCs%Atoms(I,1) 
          I2=IntCs%Atoms(I,2) 
          I3=IntCs%Atoms(I,3) 
          CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value)
          IF(ABS(Value-PI)*180.D0/PI < LinCrit) THEN  !!!linear within 1 degree
            NLinB=NLinB+1
            MarkLinB%I(I)=1
          ENDIF
        ENDIF
      ENDDO
!
    IF(NLinB/=0) THEN
!
      NIntc_New=NIntc+NLinB
      CALL NEW(IntC_New,NIntc_New)
!
        NLinB=0
      DO I=1,NIntC
        IF(MarkLinB%I(I)==0) THEN
          IntC_New%Def(NLinB+I)=Intcs%Def(I)
          IntC_New%Atoms(NLinB+I,1:4)=Intcs%Atoms(I,1:4)
          IntC_New%Value(NLinB+I)=Intcs%Value(I)
          IntC_New%Constraint(NLinB+I)=Intcs%Constraint(I)
        ELSE
          IntC_New%Def(NLinB+I)='LINB1'
          IntC_New%Atoms(NLinB+I,1:4)=Intcs%Atoms(I,1:4)
          IntC_New%Value(NLinB+I)=Intcs%Value(I)
          IntC_New%Constraint(NLinB+I)=Intcs%Constraint(I)
          NLinB=NLinB+1
          IntC_New%Def(NLinB+I)='LINB2'
          IntC_New%Atoms(NLinB+I,1:4)=Intcs%Atoms(I,1:4)
          IntC_New%Value(NLinB+I)=Intcs%Value(I)
          IntC_New%Constraint(NLinB+I)=Intcs%Constraint(I)
        ENDIF
      ENDDO
!
      CALL Delete(MarkLinB)
!
      CALL Delete(IntCs)
      NIntC=NIntC_New
      CALL New(IntCs,NIntC)
          Intcs%Def(:)=IntC_New%Def(:)
          Intcs%Atoms(:,:)=IntC_New%Atoms(:,:)
          Intcs%Value(:)=IntC_New%Value(:)      
          Intcs%Constraint(:)=IntC_New%Constraint(:)
      CALL Delete(IntC_New)
!
! Now recognize colinear atoms of the molecule and 
! introduce torsions. Introduction of extra bends is not done,
! since it is supposed to be done in previous steps
!
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/Natoms_Loc,NMax12+1/))
        CALL Get(Top12,'TOP12')
!
        CALL NEW(LinBBridge,(/2,NIntc/))
        CALL NEW(LinAtom,NIntc)
        LinBBridge%I(1:2,:)=0
        LinAtom%I(:)=0
        NLinB=0
      DO I=1,NIntc
        IF(IntCs%Def(I)(1:5)=='LINB1'.AND.LinAtom%I(IntCs%Atoms(I,1))==0) THEN 
          NLinB=NLinB+1 
          LinAtom%I(IntCs%Atoms(I,1))=1
          LinAtom%I(IntCs%Atoms(I,2))=1
          LinAtom%I(IntCs%Atoms(I,3))=1
          LinBBridge%I(1,NLinB)=IntCs%Atoms(I,1)
          LinBBridge%I(2,NLinB)=IntCs%Atoms(I,2)
! now, go on left side, use Top12, then go on right, then define torsions
            I1=LinBBridge%I(2,NLinB)
            I2=LinBBridge%I(1,NLinB)
180         CONTINUE
          DO J=1,Top12%I(I2,1)
            I3=Top12%I(I2,1+J) 
            CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value)
            IF(ABS(Value-PI)*180.D0/Pi<LinCrit) THEN !!!!linear within 1 degree
              LinBBridge%I(1,NLinB)=I3
              LinAtom%I(I3)=1
              I1=I2
              I2=I3
              GO TO 180
            ENDIF
          ENDDO
!
            I1=LinBBridge%I(1,NLinB)
            I2=LinBBridge%I(2,NLinB)
120         CONTINUE
          DO J=1,Top12%I(I2,1)
            I3=Top12%I(I2,1+J) 
            CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value)
            IF(ABS(Value-PI)*180.D0/Pi<LinCrit) THEN !!!!linear within 1 degree
              LinBBridge%I(2,NLinB)=I3
              LinAtom%I(I3)=1
              I1=I2
              I2=I3
              GO TO 120
            ENDIF
          ENDDO
!
        ENDIF
      ENDDO 
!
! bridges are set now, add torsions
! first, count number of torsions, to be added
!
      NTorsLinB=0
      DO I=1,NLinB
        I1=LinBBridge%I(1,I)
        I2=LinBBridge%I(2,I)
        NTorsLinB=NTorsLinB+Top12%I(I1,1)*Top12%I(I2,1)
      ENDDO
!
! Now generate the INTCs for the new torsions
!
      NIntc_New=NIntC+NTorsLinB
      CALL New(IntC_New,NIntc_New)
        IntC_New%Def(1:NIntC)=Intcs%Def(1:NIntC)
        IntC_New%Atoms(1:NIntC,1:4)=Intcs%Atoms(1:NIntC,1:4)
        IntC_New%Value(1:NIntC)=Intcs%Value(1:NIntC)      
        IntC_New%Constraint(1:NIntC)=Intcs%Constraint(1:NIntC)
!
          NTorsLinB=0
        DO I=1,NLinB
          I2=LinBBridge%I(1,I)
          I3=LinBBridge%I(2,I)
          DO J=1,Top12%I(I2,1)
            I1=Top12%I(I2,1+J)
            IF(LinAtom%I(I1)==0) THEN
          DO K=1,Top12%I(I3,1)
            I4=Top12%I(I3,1+K)
            IF(LinAtom%I(I4)==0) THEN
            IF(I1/=I4) THEN
              NTorsLinB=NTorsLinB+1
              IntC_New%Def(NIntC+NTorsLinB)='TORS '
              IntC_New%Atoms(NIntC+NTorsLinB,1)=I1
              IntC_New%Atoms(NIntC+NTorsLinB,2)=I2
              IntC_New%Atoms(NIntC+NTorsLinB,3)=I3
              IntC_New%Atoms(NIntC+NTorsLinB,4)=I4
              IntC_New%Value(1:NIntC)=Zero
              IntC_New%Constraint(1:NIntC)=.FALSE.
            ENDIF
            ENDIF
          ENDDO
            ENDIF
          ENDDO
        ENDDO
!
        CALL Delete(Intcs)
        NIntC=NIntC+NTorsLinB
        CALL New(Intcs,NIntC)
        Intcs%Def(1:NIntC)=IntC_New%Def(1:NIntC)
        Intcs%Atoms(1:NIntC,1:4)=IntC_New%Atoms(1:NIntC,1:4)
        Intcs%Value(1:NIntC)=IntC_New%Value(1:NIntC)      
        Intcs%Constraint(1:NIntC)=IntC_New%Constraint(1:NIntC)
        CALL Delete(IntC_New)
!
        CALL Delete(LinAtom)
        CALL Delete(LinBBridge)
        CALL Delete(Top12)
!
! Put IntCs onto HDF
!
        CALL Put(NIntC,'NIntC')
        CALL Put(IntCs,'IntCs')
!
    ENDIF
!
END SUBROUTINE GetIntCs   
!
!-------------------------------------------------------
!
SUBROUTINE INTCValue(IntCs,XYZ,Natoms_Loc)
!
! Determine value of internal coordinates.
! Input coordintes are now in atomic units!
! 
    IMPLICIT NONE
    TYPE(INTC) :: IntCs
    INTEGER :: NIntCs,I,J,K,L,I1,I2,I3,I4,Natoms_Loc
    REAL(DOUBLE),DIMENSION(1:3,1:NAtoms_Loc) :: XYZ
!
    NIntCs=SIZE(IntCs%Def(:))
!
    DO I=1,NIntCs
      IF(IntCs%Def(I)(1:4)=='STRE') THEN
        I1=IntCs%Atoms(I,1)
        I2=IntCs%Atoms(I,2)
        CALL STREValue(XYZ(1:3,I1),XYZ(1:3,I2),IntCs%Value(I))
!       IntCs%Value(I)=IntCs%Value(I)/AngstromsToAU 
      ENDIF
      IF(IntCs%Def(I)(1:4)=='BEND'.OR.IntCs%Def(I)=='LINB') THEN
        I1=IntCs%Atoms(I,1)
        I2=IntCs%Atoms(I,2)
        I3=IntCs%Atoms(I,3)
        CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),IntCs%Value(I))
      ENDIF
      IF(IntCs%Def(I)(1:4)=='OutP'.OR.IntCs%Def(I)=='TORS') THEN
        I1=IntCs%Atoms(I,1)
        I2=IntCs%Atoms(I,2)
        I3=IntCs%Atoms(I,3)
        I4=IntCs%Atoms(I,4)
        CALL TORSValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),XYZ(1:3,I4),IntCs%Value(I),IntCs%Def(I))
      ENDIF
    ENDDO 
!
END SUBROUTINE INTCValue
!
!------------------------------------------------------------
!
SUBROUTINE STREValue(XI,XJ,RIJ)
!    
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,DV
      REAL(DOUBLE) :: RIJ,RIJ2
      INTEGER :: I,J,K,L,M,N
!
      DV=XI-XJ
      RIJ2=DOT_PRODUCT(DV,DV)
      RIJ=SQRT(RIJ2)
!
END SUBROUTINE STREValue
!
!-------------------------------------------------------
!
SUBROUTINE BENDValue(XI,XJ,XK,AIJK)
!
! Calculate bond angle in radians
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,RJI,RJK,EJI,EJK
      REAL(DOUBLE) :: DIJSq,Tp,t,DJKSq,dxsq,dx,DOtj,DJK,DJI
      REAL(DOUBLE) :: sinj,SMI,SMK,SUM,thresh_B,AIJK
      INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER
!
   10 DIJSq=Zero
      DJKSq=Zero
   DO M=1,3
        Tp=XJ(M)
        T=XI(M)-Tp
        RJI(M)=T
        DIJSq=DIJSq+T*T
        T=XK(M)-Tp
        RJK(M)=T
        DJKSq=DJKSq+T*T
   ENDDO              
        DJI=SQRT(DIJSq)
        DJK=SQRT(DJKSq)
        DX=One  
        DotJ=Zero 
   DO M=1,3
        T=RJI(M)/DJI
        EJI(M)=T
        Tp=RJK(M)/DJK
        EJK(M)=Tp
        Dotj=Dotj+T*Tp
   ENDDO
!
      AIJK=ACOS(DotJ)
!     AIJK=AIJK-DBLE(INT(AIJK/PI))*PI
!
END SUBROUTINE BENDValue
!
!-------------------------------------------------------
!
SUBROUTINE TORSValue(XI,XJ,XK,XL,VTors,Def)
!
   IMPLICIT NONE
   INTEGER :: NCoord     
   INTEGER            :: I,IFAC,J,JFAC,K,KFAC,L,LFAC
   REAL(DOUBLE) :: CosNPhi,CosPhi,CosPhi2,DCos,DForce
   REAL(DOUBLE) :: IJCOSJK,KLCOSJK,D2Cos,Fact1,Fact2,NFJKL,NFIJK
   REAL(DOUBLE) :: DJK,DIJ,DKL,VTors
   REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XL,XK
   TYPE(DBL_VECT) :: FIJK,DVectIJ,DVectJK,DVectKL,FJKL
   TYPE(DBL_VECT) :: DGradI,DGradJ,DGradK,DGradL
   CHARACTER(LEN=5) :: Def
!
   CALL New(DVectIJ,3)
   CALL New(DVectJK,3)
   CALL New(DVectKL,3)
   CALL New(FIJK,3)
   CALL New(FJKL,3)
!
      DVectIJ%D=XI-XJ
      DVectJK%D=XK-XJ
      DVectKL%D=XL-XK
!
      DIJ=SQRT(DOT_PRODUCT(DVectIJ%D,DVectIJ%D))
      DJK=SQRT(DOT_PRODUCT(DVectJK%D,DVectJK%D))
      DKL=SQRT(DOT_PRODUCT(DVectKL%D,DVectKL%D))
!
      DVectJK%D=DVectJK%D/DJK
!
      IJCOSJK=DOT_PRODUCT(DVectIJ%D,DVectJK%D)
      KLCOSJK=DOT_PRODUCT(DVectKL%D,DVectJK%D)
!
! Check for linearity of three atoms within the torsion
! In case linearity appears, set Value of torsion to zero
! Current criterium is 1 degree.
!
      IF((ABS(ABS(IJCOSJK/DIJ)-One))*180.D0/PI<LinCrit .OR. &
         (ABS(ABS(KLCOSJK/DKL)-One))*180.D0/PI<LinCrit) THEN
         WRITE(*,*) 'WARNING, LINEARITY IN TORSION OR OUTP'
         WRITE(*,*) 'COORDINATES: '
         WRITE(*,*) XI(1:3)
         WRITE(*,*) XJ(1:3)
         WRITE(*,*) XK(1:3)
         WRITE(*,*) XL(1:3)
         CALL Delete(DVectIJ)
         CALL Delete(DVectJK)
         CALL Delete(DVectKL)
         CALL Delete(FIJK)
         CALL Delete(FJKL)
!        Def='BLANC' !!! do not change it, since in the next geometry it may be non-linear
         VTors=Zero
         RETURN
      ENDIF
!
      FIJK%D=DVectIJ%D-IJCOSJK*DVectJK%D
      FJKL%D=DVectKL%D-KLCOSJK*DVectJK%D
!
      NFIJK=SQRT(DOT_PRODUCT(FIJK%D,FIJK%D))
      NFJKL=SQRT(DOT_PRODUCT(FJKL%D,FJKL%D))
!
      FIJK%D=FIJK%D/NFIJK
      FJKL%D=FJKL%D/NFJKL
!
      CosPhi=DOT_PRODUCT(FIJK%D,FJKL%D)
!
      VTors=ACOS(CosPhi)
!
   CALL Delete(DVectIJ)
   CALL Delete(DVectJK)
   CALL Delete(DVectKL)
   CALL Delete(FIJK)
   CALL Delete(FJKL)
!
END SUBROUTINE TorsValue
!
!-------------------------------------------------------
!
    SUBROUTINE CoordTrf(GMLoc,TrixThresh,AInvDistanceThresh,IntCs,NIntC,VectCart,VectInt,TrfType)
!
! Routine to carry out coordinate transformations
! Refresh controls refreshing of internal coordinate set
! In case of Int-> Cartesian refresh is not allowed,
! any refresh request will be overwritten
! TrfType=1 : Cartesian -> Internal
! TrfType=2 : Internal -> Cartesian
!
    TYPE(CRDS) :: GMLoc
    TYPE(DBL_VECT) :: VectCart
    TYPE(DBL_VECT) :: VectInt,VectAux
    TYPE(BCSR) :: SpVectCart,SpVectInt,SpB,SpBt,Gc,Z,Zt,SpVectActInt
    TYPE(BCSR) :: SpVectDiff,SpVectAux
    REAL(DOUBLE) :: TrixThresh,AInvDistanceThresh,DiffMax,DiffNorm
    REAL(DOUBLE) :: GrdTrfCrit,CooTrfCrit,ScaleTo,MaxGradDiff,Sum
    INTEGER :: NCart,I,II,J,TrfType,NIntC,DiffLength
    INTEGER :: NVBlocksB,NHBlocksB,BlockGeomSize
    INTEGER :: MaxIt_GrdTrf,MaxIt_CooTrf
    TYPE(INTC) :: IntCs
    TYPE(BMATR):: B
    LOGICAL :: RefreshB,RefreshAct
!
!   CALL OpenASCII(OutFile,Out)
    GrdTrfCrit=TrixThresh*1.D0 !!! For Gradient transformation
!   CooTrfCrit=TrixThresh*1.D0 !!! For Backtransformation
    CooTrfCrit=1.D-4 !!! For Backtransformation
    MaxIt_GrdTrf=20
    MaxIt_CooTrf=20
    NCart=3*GMLoc%Natms
    MaxGradDiff=0.5D0
!
! Refresh B matrix during iterative back-trf, if displacements are too large?
!
    RefreshB=.TRUE.
    RefreshAct=.TRUE.
!
! Set some block-matrix params
!
    CALL Get(BlockGeomSize,'BlkGeomSize')
    MaxBlks=Natoms*4*BlockGeomSize !!! for primitive internals, suppose, that Cartesians are ordered by atoms (atom blocking is efficient)
    MaxNon0=(BlockGeomSize*BlockGeomSize)*MaxBlks !!! including zeros of blocks
!
! Turn input vector into sparse matrix form
! And initialize Transformed vectors
!
      CALL New(VectAux,NBasF)
!
    IF(TrfType==1) THEN
!
! Print List of internal coordinates with actual values
!
      CALL INTCValue(IntCs,GMLoc%Carts%D,GMLoc%Natms)
      WRITE(*,*) 'INTERNAL COORDINATES'
      WRITE(*,*) 'DEFINITION  ATOMS_INVOLVED  VALUE'
      WRITE(Out,*) 'INTERNAL COORDINATES'
      WRITE(Out,*) 'DEFINITION  ATOMS_INVOLVED  VALUE'
      DO I=1,NIntC
        IF(IntCs%Def(I)(1:4)=='STRE') THEN
          SUM=IntCs%Value(I)/AngstromsToAu
        ELSE
          SUM=IntCs%Value(I)*180.D0/PI
        ENDIF
        WRITE(*,111) I,IntCs%Def(I),IntCs%Atoms(I,1:4),SUM
        WRITE(Out,111) I,IntCs%Def(I),IntCs%Atoms(I,1:4),SUM
      ENDDO
 111 FORMAT(I4,2X,A5,2X,4I3,2X,F12.6)
!
      VectAux%D(1:NCart)=VectCart%D(1:NCart)
      IF(NCart<NBasF) VectAux%D(NCart+1:NBasF)=Zero
    CALL Set_BCSR_EQ_VECT(SpVectCart,VectAux)
!initialization
      VectAux%D=Zero     
    CALL Set_BCSR_EQ_VECT(SpVectInt,VectAux)
!
    ELSE IF(TrfType==2) THEN !---------------------------
!
! Calc values of internals in atomic unit, and add displacement,
! which is stored in VectInt. Then convert into Sparse matrx repr.
!
      CALL INTCValue(IntCs,GMLoc%Carts%D,GMLoc%Natms)
!
      VectAux%D(1:NIntC)=VectInt%D(1:NIntC)+IntCs%Value(1:NIntC)
      CALL MapBackAngle(IntCs,NIntC,VectAux%D) 
      IF(NIntC<NBasF) VectAux%D(NIntC+1:NBasF)=Zero
    CALL Set_BCSR_EQ_VECT(SpVectInt,VectAux)
!
!initialization of new Cartesians
!
      DO I=1,GMLoc%Natms 
        J=3*(I-1)
        VectAux%D(J+1)=GMLoc%Carts%D(1,I)
        VectAux%D(J+2)=GMLoc%Carts%D(2,I)
        VectAux%D(J+3)=GMLoc%Carts%D(3,I)
      ENDDO
      IF(NCart<NBasF) VectAux%D(NCart+1:NBasF)=Zero
    CALL Set_BCSR_EQ_VECT(SpVectCart,VectAux)
!
    ENDIF
!-------------------------------------------------------------------
!
! Perform Coordinate transformation on the input vectors
!
    IF(TrfType==1) THEN
!
      WRITE(*,*) 'Gradient transformation, No. Int. Coords= ',NIntC
      WRITE(Out,*) 'Gradient transformation, No. Int. Coords= ',NIntC
!
! Cartesian --> Internal transformation
!
! Get B matrix in sparse representations
! and get Inverse factors Z and Zt
!
    CALL GetSpBMatr(GMLoc,IntCs,NIntC,B,SpB,SpBt) 
    CALL GetCoordTrfZ(GMLoc,SpB,SpBt,TrixThresh,AInvDistanceThresh,Z,Zt)
!
      II=0
100   CONTINUE      
!
      CALL Multiply(SpBt,SpVectInt,SpVectAux)
!       CALL PPrint(SpVectAux,'Bt*gi',Unit_O=6)
      CALL Multiply(SpVectAux,-One)
      CALL Add(SpVectCart,SpVectAux,SpVectDiff)
!       CALL PPrint(SpVectDiff,'gc-Bt*gi',Unit_O=6)
      CALL Delete(SpVectAux)
      CALL Multiply(Zt,SpVectDiff,SpVectAux)
!       CALL PPrint(SpVectAux,'Zt*[gc-Bt*gi]',Unit_O=6)
      CALL Delete(SpVectDiff)
      CALL Multiply(Z,SpVectAux,SpVectDiff)
!       CALL PPrint(SpVectDiff,'Z*Zt*[gc-Bt*gi]',Unit_O=6)
      CALL Delete(SpVectAux)
      CALL Multiply(SpB,SpVectDiff,SpVectAux)
!       CALL PPrint(SpVectAux,'B*Z*Zt*[gc-Bt*gi]',Unit_O=6)
      CALL Delete(SpVectDiff)
!
! Check convergence
!
      II=II+1
        DiffLength=SpVectAux%NNon0
        DiffMax=Zero
      DO I=1,DiffLength     
        DiffMax=MAX(DiffMax,ABS(SpVectAux%MTrix%D(I)))
      ENDDO
        DiffNorm=DOT_PRODUCT(SpVectAux%MTrix%D(1:DiffLength),SpVectAux%MTrix%D(1:DiffLength))
        DiffNorm=DiffNorm/DBLE(NIntC)
!
! IF DiffMax is too large, eg. due to the 'bad' quality 
! of the preconditioner, rescale gradients
!
      IF(DiffMax>MaxGradDiff) THEN
        WRITE(*,*) 'Rescale Step from ',DiffMax,' to ',MaxGradDiff
        WRITE(Out,*) 'Rescale Step from ',DiffMax,' to ',MaxGradDiff
        SUM=MaxGradDiff/DiffMax
        SpVectAux%MTrix%D(:)=SUM*SpVectAux%MTrix%D(:)
        DiffMax=MaxGradDiff
      ENDIF
!
      CALL Add(SpVectInt,SpVectAux,SpVectDiff)
      CALL Delete(SpVectAux)
!       CALL PPrint(SpVectDiff,'gi+B*Z*Zt*[gc-Bt*gi]',Unit_O=6)
      CALL Delete(SpVectInt)
      CALL Set_BCSR_EQ_BCSR(SpVectInt,SpVectDiff)
      CALL Delete(SpVectDiff)
!      
! Review iteration
!
      WRITE(*,110) II,DiffMax,DiffNorm
      WRITE(Out,110) II,DiffMax,DiffNorm
110   FORMAT('Grad Trf, step= ',I3,' MaxChange= ',F12.6,' ChangeNorm= ',F12.6)
!      
      IF(DiffMax>GrdTrfCrit.AND.II<=MaxIt_GrdTrf) THEN
          GO TO 100      
      ELSE
          IF(II>MaxIt_GrdTrf) THEN
            WRITE(*,*) 'Stop Gradient Trf, max. number of Iterations exceeded!'
            WRITE(*,*) 'Use current gradient vector!'
            WRITE(Out,*) 'Stop Gradient Trf, max. number of Iterations exceeded!'
            WRITE(Out,*) 'Use current gradient vector!'
          ENDIF
          CALL Put(SpB,'SpB')
          CALL Put(SpBt,'SpBt')
          CALL Put(Z,'Z-fact')
          CALL Put(Zt,'Z-fact-t')
          CALL Delete(SpB)
          CALL Delete(SpBt)
          CALL Delete(Z)
          CALL Delete(Zt)
      ENDIF
!
      WRITE(*,120) II
      WRITE(Out,120) II
120   FORMAT('Gradient transformation converged in ',I3,' steps')
!
! Fill result into VectInt
!
      CALL Set_DBL_VECT_EQ_BCSRColVect(VectAux,SpVectInt)
      VectInt%D(1:NIntC)=VectAux%D(1:NIntC)
!
      CALL Delete(SpVectCart)
      CALL Delete(SpVectInt)
!
    ELSE !----------------------------------------------------------------
!
! Internal --> Cartesian transformation
!
      WRITE(*,*) 'Iterative back-transformation, No. Int. Coords= ',NIntC
      WRITE(Out,*) 'Iterative back-transformation, No. Int. Coords= ',NIntC
      II=0
200   CONTINUE
!
! Get B matrix in sparse representations
! and get Inverse factors Z and Zt
!
      IF(II==0) THEN
        CALL Get(SpB,'SpB')
        CALL Get(SpBt,'SpBt')
        CALL Get(Z,'Z-fact')
        CALL Get(Zt,'Z-fact-t')
      ELSE IF(RefreshB.AND.RefreshAct) THEN
        CALL GetSpBMatr(GMLoc,IntCs,NIntC,B,SpB,SpBt) 
        CALL GetCoordTrfZ(GMLoc,SpB,SpBt,TrixThresh,AInvDistanceThresh,Z,Zt)
      ENDIF
!
! Set actual internals at present Cartesian coordinates
!
      IF(II>0) CALL INTCValue(IntCs,GMLoc%Carts%D,GMLoc%Natms)
!
      VectAux%D(1:NIntC)=IntCs%Value(1:NIntC)
      IF(NIntC<NBasF) VectAux%D(NIntC+1:NBasF)=Zero
    CALL Set_BCSR_EQ_VECT(SpVectActInt,VectAux)
!
! Calculate difference between required and actual internals
!
      CALL Multiply(SpVectActInt,-One)
      CALL Add(SpVectActInt,SpVectInt,SpVectDiff)
!
! Do transformation
!
      CALL Multiply(SpBt,SpVectDiff,SpVectAux)
      CALL Delete(SpVectDiff)
      CALL Multiply(Zt,SpVectAux,SpVectDiff)
      CALL Delete(SpVectAux)
      CALL Multiply(Z,SpVectDiff,SpVectAux)
      CALL Delete(SpVectDiff)
!
! Check convergence
!
      II=II+1
        DiffLength=SpVectAux%NNon0
        DiffMax=Zero
      DO I=1,DiffLength     
        DiffMax=MAX(DiffMax,ABS(SpVectAux%MTrix%D(I)))
      ENDDO
        DiffNorm=DOT_PRODUCT(SpVectAux%MTrix%D(1:DiffLength),SpVectAux%MTrix%D(1:DiffLength))
        DiffNorm=DiffNorm/DBLE(NCart)
!
! Scale Displacements
!
      ScaleTo=0.6D0  !!! In atomic units. No larger cartesian displacements in a back-trf step are allowed
      IF(DiffMax>ScaleTo) THEN
        CALL Multiply(SpVectAux,ScaleTo/DiffMax)
        DiffMax=ScaleTo
      ENDIF
      IF(DiffMax>ScaleTo/Two) THEN
        RefreshAct=.TRUE.
      ELSE
        RefreshAct=.FALSE.
      ENDIF
!
! Modify Cartesians
!
      CALL Add(SpVectCart,SpVectAux,SpVectDiff)
      CALL Delete(SpVectCart)
      CALL Set_BCSR_EQ_BCSR(SpVectCart,SpVectDiff)
      CALL Delete(SpVectDiff)
!
! tidy up
!
      CALL Delete(SpVectAux)
      CALL Delete(SpVectActInt) 
!
! Fill new Cartesians into GMLoc
!
      CALL Set_DBL_VECT_EQ_BCSRColVect(VectAux,SpVectCart)
      DO I=1,GMLoc%Natms 
        J=3*(I-1)
        GMLoc%Carts%D(1,I)=VectAux%D(J+1)
        GMLoc%Carts%D(2,I)=VectAux%D(J+2)
        GMLoc%Carts%D(3,I)=VectAux%D(J+3)
      ENDDO
!
! Review iteration
!
      WRITE(*,210) II,DiffMax,DiffNorm
      WRITE(Out,210) II,DiffMax,DiffNorm
210   FORMAT('Coord Back-Trf, step= ',I3,' MaxChange= ',F12.6,' ChangeNorm= ',F12.6)
!      
      IF(DiffMax>CooTrfCrit .AND. II<=MaxIt_CooTrf) THEN 
        IF(RefreshB.AND.RefreshAct) THEN
          CALL Delete(SpB)
          CALL Delete(SpBt)
          CALL Delete(Z)
          CALL Delete(Zt)
        ENDIF
        GO TO 200
      ELSE
          IF(II>MaxIt_CooTrf) THEN
            WRITE(*,*) 'Stop Coord Back-Trf, max. number of Iterations exceeded!'
            WRITE(*,*) 'Use current gradient vector!'
            WRITE(Out,*) 'Stop Coord Back-Trf, max. number of Iterations exceeded!'
            WRITE(Out,*) 'Use current Cartesian coordinates!'
          ENDIF
          CALL Delete(SpB)
          CALL Delete(SpBt)
          CALL Delete(Z)
          CALL Delete(Zt)
      ENDIF
!
      WRITE(*,220) II
      WRITE(Out,220) II
220   FORMAT('Coordinate back-transformation converged in ',I3,' steps')
!
    ENDIF !-------------------------------------------------------------
!
    CALL Delete(VectAux)
!
!   Close(Out)
!
END SUBROUTINE CoordTrf
!-------------------------------------------------------
    SUBROUTINE SpectralShift(Matr,Shift)
!
    IMPLICIT NONE
    TYPE(BCSR) :: Matr,SpAux,Identity
    REAL(DOUBLE) :: Shift,Buf
!
    CALL SetToI(Identity)
    CALL Multiply(Identity,Shift)
    CALL SetEq(SpAux,Matr)
    CALL Delete(Matr)
    CALL Add(SpAux,Identity,Matr)
    CALL Delete(SpAux)
    CALL Delete(Identity)
!
    END SUBROUTINE SpectralShift
!
!-------------------------------------------------------
!
    SUBROUTINE GetCoordTrfZ(GMLoc,SpB,SpBt,TrixThresh,AInvDistanceThresh,Z,Zt)
!
! Calculate the overlap, Gc=Bt*B and its inverse factors: Gc^-1=Zt*Z
!
    TYPE(BCSR) :: SpB,SpBt,Gc,Z,Zt,MtrxAux,MtrxAux2
    REAL(DOUBLE) :: TrixThresh,AInvDistanceThresh 
    INTEGER :: I,J,BlockGeomSize
    TYPE(CRDS) :: GMLoc
!
! redefine some block data
!
    CALL Get(BlockGeomSize,'BlkGeomSize')
!   MaxBlks=MAX(MaxBlks,Natoms*Natoms/4) 
    MaxBlks=MAX(MaxBlks,Natoms*Natoms) 
    MaxNon0=(BlockGeomSize*BlockGeomSize)*MaxBlks !!! including zeros of blocks
!
! Calculate Gc=Bt*B matrix!
!
    CALL Multiply(SpBt,SpB,Gc)
!
! Apply spectral shift to Gc!
!
    CALL SpectralShift(Gc,1.D-3) !!! large spectral shift may result in greater sparsity
!
! Factorize Gc inverse! Gc^-1=Z*Zt
!
!   CALL BlockedAInv(Gc,TrixThresh,GMLoc,AInvDistanceThresh,Z,Zt,PerfMon)
    CALL BlockedAInv(Gc,TrixThresh,GMLoc,1.d6,Z,Zt,PerfMon)
!
!! check inverse
!!
!    CALL Multiply(Zt,Gc,MtrxAux)
!    CALL Multiply(Z,MtrxAux,MtrxAux2)
!    CALL PPrint(MtrxAux2,'Unity? ',Unit_O=6)
!    CALL Delete(MtrxAux)
!    CALL Delete(MtrxAux2)
!
! Tidy up
!
    CALL Delete(Gc)
!
    END SUBROUTINE GetCoordTrfZ
!
!------------------------------------------------------------------------
!
    SUBROUTINE GetSpBMatr(GMLoc,IntCs,NIntC,B,SpB,SpBt) 
!
! Generate vibrational B matrix in quasi-sparse TYPE(BMATR) representation 
! and transform into sparse one
!
    TYPE(CRDS) :: GMLoc
    TYPE(BCSR) :: SpB,SpBt
    TYPE(INTC) :: IntCs
    INTEGER    :: NIntC,NCart,I
    TYPE(BMATR):: B
!
! Calculate B matrix in Atomic Units
!
    CALL BMatrix(GMLoc%Natms,GMLoc%Carts%D,NIntC,IntCs,B)
!
! Calculate 'folded' B matrix for composit 'B matrix coordinates'
!
!   CALL FoldBMatrix(GMLoc%Natms,GMLoc%Carts%D,NIntC,IntCs,B)
!
    NCart=3*GMLoc%Natms
!
! Now, turn B matrix into sparse blocked representation
!
    CALL Set_BCSR_EQ_BMATR(SpB,B,NCart=NCart)
!     CALL PPrint(SpB,'SpB',Unit_O=6)
!
! Generate sparse transpose B
!
    CALL Set_BCSR_EQ_BMATR(SpBt,B,NCart=NCart,Transpose=.True.)
!     CALL PPrint(SpBt,'SpBt',Unit_O=6)
!
! Tidy up
!
    CALL Delete(B)
!
    END SUBROUTINE GetSpBMatr
!
!---------------------------------------------------------
!
    SUBROUTINE MapBackAngle(IntCs,NIntC,VectAux)
    IMPLICIT NONE
    TYPE(INTC) :: IntCs
    INTEGER :: I,J,NIntC
    REAL(DOUBLE),DIMENSION(:) :: VectAux
    REAL(DOUBLE) :: SUM,TwoPi
!
! Map back angles into the ranges the iterative back-transformation
! can cope with.
!
    TwoPi=Two*PI
!
    DO I=1,NIntC
      IF(IntCs%Def(I)(1:4)=='BEND'.OR.&
         IntCs%Def(I)(1:4)=='LINB') THEN
        SUM=VectAux(I)
          J=INT(SUM/TwoPi)
        SUM=SUM-J*TwoPi
        IF(SUM<Zero) SUM=-SUM
        IF(SUM>PI) THEN
          VectAux(I)=TwoPi-SUM
        ELSE
          VectAux(I)=SUM
        ENDIF
      ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
        SUM=VectAux(I)
          J=INT(SUM/TwoPi)
        SUM=SUM-J*TwoPi
        IF(SUM<Zero) SUM=TwoPi+SUM
          VectAux(I)=SUM
      ENDIF
    ENDDO
!
    END SUBROUTINE MapBackAngle
!
!-------------------------------------------------------
!
    SUBROUTINE GDIIS(CurStat,GDIISMemory,GMLoc)
    IMPLICIT NONE
    TYPE(CRDS) :: GMLoc
    INTEGER :: I,J,K,L,GDIISMemory,DimGDIIS,IGeom,DimOverl
    INTEGER :: II,ILow
    INTEGER,DIMENSION(3) :: CurStat
    TYPE(DBL_RNK2) :: Displacements,OverlapOfDispl,Structures,ActCarts
    REAL(DOUBLE) :: Sum
!
! Current implementation runs on Cartesian displacements only
!
    DimGDIIS=3*GMLoc%Natms
!
! Get Displacements from geometries, stored in HDF
! and fill into Displacements columns.
!
    II=CurStat(3) ! number of geometries stored in HDF
    CALL New(Displacements,(/DimGDIIS,II+1/))
    CALL New(Structures,(/DimGDIIS,II+1/))
    CALL New(ActCarts,(/3,GMLoc%Natms/))
!
! Later, tag geometries for basis sets, as well.
!
    ILow=II-GDIISMemory
    IF(ILow<0) ILow=0
    DO IGeom=ILow+1,II
      CALL Get(ActCarts,'cartesians',Tag_O=TRIM(IntToChar(IGeom)))
        K=IGeom-ILow
      DO I=1,GMLoc%Natms
        J=3*(I-1)
        Structures%D(J+1,K)=ActCarts%D(1,I)
        Structures%D(J+2,K)=ActCarts%D(2,I)
        Structures%D(J+3,K)=ActCarts%D(3,I)
      ENDDO
    ENDDO
        K=II-ILow
      DO I=1,GMLoc%Natms
        J=3*(I-1)
        Structures%D(J+1,K+1)=GMLoc%Carts%D(1,I)
        Structures%D(J+2,K+1)=GMLoc%Carts%D(2,I)
        Structures%D(J+3,K+1)=GMLoc%Carts%D(3,I)
      ENDDO
!
! get displacements
!
    DO I=1,II-ILow
      IGeom=ILow+I
      Displacements%D(1:DimGDIIS,I)=Structures%D(1:DimGDIIS,IGeom+1)-Structures%D(1:DimGDIIS,IGeom)
    ENDDO
!
! Now calculate overlap. Presently only non-sparse representation is available
! dimension of overlap is (II-Ilow)+1
!
    DimOverl=II-Ilow
    CALL New(OverlapOfDispl,(/DimOverl+1,DimOverl+1/))
!
    DO I=1,DimOverl
        SUM=DOT_PRODUCT(Displacements%D(1:DimGDIIS,I),Displacements%D(1:DimGDIIS,I))
        OverlapOfDispl%D(I,I)=SUM
      DO J=I+1,DimOverl
        SUM=DOT_PRODUCT(Displacements%D(1:DimGDIIS,I),Displacements%D(1:DimGDIIS,J))
        OverlapOfDispl%D(I,J)=SUM
        OverlapOfDispl%D(J,I)=SUM
      ENDDO
    ENDDO
!
! Add 'surface' terms to overlap
!
    OverlapOfDispl%D(1:DimOverl,DimOverl+1)=One
    OverlapOfDispl%D(DimOverl+1,1:DimOverl)=One
    OverlapOfDispl%D(DimOverl+1,DimOverl+1)=Zero 
!
! Now, Calculate SVD inverse of the Overlap
!
!
! Tidy up
!
    CALL Delete(OverlapOfDispl)
    CALL Delete(ActCarts)
    CALL Delete(Structures)
    CALL Delete(Displacements)
!
    END SUBROUTINE GDIIS
!
!--------------------------------------------------------------------
!
!    SUBROUTINE FoldBMatrix(NatLoc,XYZ,NIntC,IntCs,B,BFolded)
!!
!    TYPE(INTC) :: IntCs
!    TYPE(BMATR) :: B,BFolded
!    INTEGER :: I,J,K,L,NIntC,NatLoc
!!
!! Find three atoms which are non-linear,
!! possibly find the ones of having the larges distance and angles
!! close to 60 degree
!! Inclusion of the Cartesian coordinates of three atoms would also 
!! help to construct
!!
!!
!!
!    END SUBROUTINE FoldBMatrix
!
#endif
!-------------------------------------------------------
END MODULE InCoords

