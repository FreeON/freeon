!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form <<These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics>>, and given appropriate citation.  
!------------------------------------------------------------------------------
MODULE InCoords
!
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE GlobalCharacters
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

!
IMPLICIT NONE
!
CONTAINS
!--------------------------------------------------------------
!
SUBROUTINE Topology_12(NatmsLoc,NBond,BondI,BondJ,Top12,InfFile,Tag_O)
!
! Set up a table which shows the atom numbers of Atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based on input list
! of bonds
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2) :: Top12_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBond,NatmsLoc,NMax12
INTEGER,DIMENSION(1:NBond) :: BondI,BondJ
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile,Tag_O
!
    NMax12=5
    CALL New(Top12,(/NatmsLoc,NMax12+1/))
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
        CALL New(Top12_2,(/NatmsLoc,NMax12+1/))
        Top12_2%I(1:NatmsLoc,1:NMax12+1)=0 
        Top12_2%I(1:NatmsLoc,1:NMax12+1-5)=Top12%I(1:NatmsLoc,1:NMax12+1-5)
        CALL Delete(Top12)
        CALL New(Top12,(/NatmsLoc,NMax12+1/))
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
      IF(PRESENT(Tag_O)) THEN
        CALL Put(NMax12,'NMax12'//TRIM(Tag_O))
        CALL Put(Top12,'Top12'//TRIM(Tag_O))
      ELSE
        CALL Put(NMax12,'NMax12')
        CALL Put(Top12,'Top12')
      ENDIF
    ENDIF
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
!
END SUBROUTINE Topology_12 
!--------------------------------------------------------------
!
SUBROUTINE Topology_13(NatmsLoc,Top12,Top13,InfFile,Tag_O)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12
TYPE(INT_RNK2),OPTIONAL :: Top13
TYPE(INT_RNK2) :: Top13_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NatmsLoc,NMax13,NMax12,KK,IN12,JN12
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile,Tag_O
!
    IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NatmsLoc,K/))
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
    CALL New(Top13,(/NatmsLoc,K/))
    Top13%I(1:NatmsLoc,1:NMax13+1)=0 
!
    DO II=1,NatmsLoc
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
        CALL New(Top13_2,(/NatmsLoc,NMax13+1/))
        Top13_2%I(1:NatmsLoc,1:NMax13+1)=0 
        Top13_2%I(1:NatmsLoc,1:NMax13+1-10)=Top13%I(1:NatmsLoc,1:NMax13+1-10)
        CALL Delete(Top13)
        CALL New(Top13,(/NatmsLoc,NMax13+1/))
        Top13%I(1:NatmsLoc,1:NMax13+1)=Top13_2%I(1:NatmsLoc,1:NMax13+1)
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
      IF(PRESENT(Tag_O)) THEN
        CALL Put(NMax13,'NMax13'//TRIM(Tag_O))
        CALL Put(Top13,'Top13'//TRIM(Tag_O))
      ELSE
        CALL Put(NMax13,'NMax13')
        CALL Put(Top13,'Top13')
      ENDIF
    ENDIF
!
    IF(.NOT.PRESENT(Top13)) CALL Delete(Top13)
!
END SUBROUTINE Topology_13 
!--------------------------------------------------------------
!
SUBROUTINE Topology_14(NatmsLoc,Top12,Top14,InfFile,Tag_O)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL:: Top12
TYPE(INT_RNK2),OPTIONAL:: Top14
TYPE(INT_RNK2) :: Top14_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
INTEGER :: NatmsLoc,NMax14,NMax12,IN12,JN12,KN12
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile,Tag_O
!
     IF(.NOT.PRESENT(Top12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        K=NMax12+1
        CALL New(Top12,(/NatmsLoc,K/))
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
    CALL New(Top14,(/NatmsLoc,K/))
    Top14%I(1:NatmsLoc,1:NMax14+1)=0 
!
    DO II=1,NatmsLoc
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
        CALL New(Top14_2,(/NatmsLoc,NMax14+1/))
        Top14_2%I(1:NatmsLoc,1:NMax14+1)=0 
        Top14_2%I(1:NatmsLoc,1:NMax14+1-10)=Top14%I(1:NatmsLoc,1:NMax14+1-10)
        CALL Delete(Top14)
        CALL New(Top14,(/NatmsLoc,NMax14+1/))
        Top14%I(1:NatmsLoc,1:NMax14+1)=Top14_2%I(1:NatmsLoc,1:NMax14+1)
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
      IF(PRESENT(Tag_O)) THEN
        CALL Put(NMax14,'NMax14'//TRIM(Tag_O))
        CALL Put(Top14,'Top14'//TRIM(Tag_O))
      ELSE
        CALL Put(NMax14,'NMax14')
        CALL Put(Top14,'Top14')
      ENDIF
    ENDIF
!
    IF(.NOT.PRESENT(Top12)) CALL Delete(Top12)
    IF(.NOT.PRESENT(Top14)) CALL Delete(Top14)
!
END SUBROUTINE Topology_14 
!--------------------------------------------------------------
!
SUBROUTINE SORT_INTO_Box1(BoxSize,C,NatmsLoc,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
! sort the Atoms of a molecule into Boxes
!
! BoxSize: linear Box Size
!
IMPLICIT NONE
INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NatmsLoc
REAL(DOUBLE) :: BoxSize,VBIG,C(1:3,NatmsLoc),BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
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
    DO I=1,NatmsLoc
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
SUBROUTINE SORT_INTO_Box2(BoxSize,C,NatmsLoc,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,InfFile,ISet)
!
! sort the Atoms of a molecule into Boxes
!
! BoxI(I) : contains the ordering (box) number of the first atom of the I-th Box (like in sparse row-wise)
! BoxJ(J) : gives the original serial number of the atom desribed by the J-th ordering number
! C: contains Cartesian coordinates of Atoms
! BoxSize: linear Box Size
!
IMPLICIT NONE
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
INTEGER,OPTIONAL :: ISET 
TYPE(INT_VECT),OPTIONAL :: BoxI1,BoxJ1
TYPE(INT_VECT) :: BoxI,BoxJ
TYPE(INT_RNK2) :: ISign
TYPE(INT_RNK3) :: BoxCounter 
INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NatmsLoc
REAL(DOUBLE) :: BoxSize,VBIG,C(1:3,NatmsLoc),BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
SAVE VBIG
DATA VBIG/1.D+90/ 
!
  NBox=NX*NY*NZ
  CALL New(BoxI,NBox+1)
  CALL New(BoxJ,NatmsLoc)
!
  CALL New(ISign,(/2,NatmsLoc/))
!
  CALL New(BoxCounter,(/NX,NY,NZ/))
  BoxCounter%I(1:NX,1:NY,1:NZ)=0
  BoxI%I(1:NBox+1)=0
!
! COUNT NUMBER OF Atoms IN THE Box
!
  DO I=1,NatmsLoc
!
! identify Box
!     
    IX=INT((C(1,I)-BXMIN)/BoxSize)+1
    IY=INT((C(2,I)-BYMIN)/BoxSize)+1
    IZ=INT((C(3,I)-BZMIN)/BoxSize)+1
!     
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY !order parameter: ZXY
!
    BoxCounter%I(IX,IY,IZ)=BoxCounter%I(IX,IY,IZ)+1
    ISign%I(1,I)=IORD
    ISign%I(2,I)=BoxCounter%I(IX,IY,IZ) !!! shows both Box and index within Box
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
    BoxI%I(IORD+1)=BoxI%I(IORD)+BoxCounter%I(IX,IY,IZ)
!     
  ENDDO
  ENDDO
  ENDDO
!
! Set up contents of Boxes as represented in A, sparse row-wise
!
  DO I=1,NatmsLoc
    IORD=ISign%I(1,I)
    IADD=ISign%I(2,I)
    BoxJ%I(BoxI%I(IORD)-1+IADD)=I
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
  CALL Delete(ISign)
  CALL Delete(BoxCounter)
  CALL Delete(BoxI)
  CALL Delete(BoxJ)
!
END SUBROUTINE SORT_INTO_Box2
!
!----------------------------------------------------------------
!
SUBROUTINE SORT_INTO_Box(BoxSize,C,NatmsLoc,InfFile,ISet)
IMPLICIT NONE
INTEGER,OPTIONAL :: ISet
INTEGER :: NatmsLoc,NX,NY,NZ,NBox
REAL(DOUBLE) :: BoxSize,BXMIN,BYMIN,BZMIN
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
TYPE(INT_VECT) :: BoxI1,BoxJ1
REAL(DOUBLE),DIMENSION(1:3,1:NatmsLoc) :: C
!
CALL SORT_INTO_Box1(BoxSize,C,NatmsLoc,NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
NBox=NX*NY*NZ
CALL New(BoxI1,NBox+1)
CALL New(BoxJ1,NatmsLoc)
!
CALL SORT_INTO_Box2(BoxSize,C,NatmsLoc,NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,InfFile,ISet)
!
CALL Delete(BoxI1)
CALL Delete(BoxJ1)
!
END SUBROUTINE SORT_INTO_Box
!----------------------------------------------------------------
!
SUBROUTINE Topologies_MM(NatmsLoc,NBond,BondI,BondJ,InfFile,Top12OUT)
!
! Set up a table which shows the atom numbers of Atoms 
! connected to a certain atom by the input bonds (Topology mtr)
! Here, the generation of the Topology mtr is based 
! EXCLUSIVELY on input list of bonds
! WARNING! This subroutine is mainly used to set up the
! topology files necessary for the calculation of 
! MM exclusion energies!!!
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top12OUT
TYPE(INT_RNK2) :: Top12OUT_2
TYPE(INT_RNK2) :: Top12,Top13,Top14
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBond,NatmsLoc
INTEGER,DIMENSION(1:NBond) :: BondI,BondJ
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile
CHARACTER(LEN=DefAULT_CHR_LEN) :: Tag_O
!
   Tag_O='MMTop'
!
   CALL Topology_12(NatmsLoc,NBond,BondI,BondJ,Top12=Top12,InfFile=InfFile,Tag_O=Tag_O)
   CALL Topology_13(NatmsLoc,Top12,Top13=Top13,InfFile=InfFile,Tag_O=Tag_O)
   CALL Topology_14(NatmsLoc,Top12,Top14=Top14,InfFile=InfFile,Tag_O=Tag_O)
   CALL Excl_List(NatmsLoc,Top12,Top13,Top14=Top14,InfFile=InfFile,Tag_O=Tag_O)
   CALL Excl_List14(NatmsLoc,Top12,Top13,Top14=Top14,InfFile=InfFile,Tag_O=Tag_O)
!
   IF(PRESENT(Top12OUT)) THEN
     IF(AllocQ(Top12OUT%Alloc)) CALL Delete(Top12OUT)
     N=Size(Top12%I,2)
     CALL New(Top12OUT_2,(/NatmsLoc,N/))
     Top12OUT_2=Top12
     CALL Delete(Top12)
     CALL Delete(Top13)
     CALL Delete(Top14)
     CALL New(Top12OUT,(/NatmsLoc,N/))
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
SUBROUTINE Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl_Out,InfFile,Tag_O)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile,Tag_O
INTEGER :: NatmsLoc,NMax12,NMax13,NMax14,NMax_Excl,NNew,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ
!
    IF(PRESENT(Top12)) THEN
      NMax12=Size(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NatmsLoc,NMax12+1/))
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
        CALL New(Top13,(/NatmsLoc,NMax13+1/))
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
        CALL New(Top14,(/NatmsLoc,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl
!
    NMax_Excl=NMax12+NMax13+NMax14
    CALL New(Top_Excl,(/NatmsLoc,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
    Top_Excl%I(1:NatmsLoc,1:NMax12+1)=Top12%I(1:NatmsLoc,1:NMax12+1)
!
! Now merge Topologies, in order to avoid double counting in 
! Exclusion energies
!
      NMax_Excl_Out=0
    DO I=1,NatmsLoc
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
      CALL New(Top_Excl_Out,(/NatmsLoc,NMax_Excl_Out+1/))
      Top_Excl_Out%I(1:NatmsLoc,1:NMax_Excl_Out+1)=&
      Top_Excl%I(1:NatmsLoc,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        IF(PRESENT(Tag_O)) THEN
          CALL Put(NMax_Excl_Out,'NMax_Excl'//TRIM(Tag_O))
          CALL Put(Top_Excl_Out,'TOP_Excl'//TRIM(Tag_O))
        ELSE
          CALL Put(NMax_Excl_Out,'NMax_Excl')
          CALL Put(Top_Excl_Out,'TOP_Excl')
        ENDIF
      ENDIF
    ELSE IF(PRESENT(InfFile)) THEN
      CALL New(Top_New,(/NatmsLoc,NMax_Excl_Out+1/))
      Top_New%I(1:NatmsLoc,1:NMax_Excl_Out+1)=&
      Top_Excl%I(1:NatmsLoc,1:NMax_Excl_Out+1)
        IF(PRESENT(Tag_O)) THEN
          CALL Put(NMax_Excl_Out,'NMax_Excl'//TRIM(Tag_O))
          CALL Put(Top_New,'TOP_Excl'//TRIM(Tag_O))
        ELSE
          CALL Put(NMax_Excl_Out,'NMax_Excl')
          CALL Put(Top_New,'TOP_Excl')
        ENDIF
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
SUBROUTINE Excl_List14(NatmsLoc,Top12,Top13,Top14,Top_Excl_Out,InfFile,Tag_O)
!
! This subroutine merges Topological information
! to get the list for Exclusion energy calculation
! of Atoms in 14 distance. From the Top14 list
! Top13 and Top12 occurences must be filtered out
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: Top_Excl_Out,Top12,Top13,Top14
TYPE(INT_RNK2) :: Top_Excl,Top_New
CHARACTER(LEN=DefAULT_CHR_LEN),OPTIONAL :: InfFile,Tag_O
INTEGER :: NatmsLoc,NMax12,NMax13,NMax14,NMax_Excl,NNew,NOLD
INTEGER :: NMax_Excl_Out,I,J,K,KK,JJ,NExcl,N12,N13
!
    IF(PRESENT(Top12)) THEN
      NMax12=Size(Top12%I,2)-1
    ELSE
      IF(PRESENT(InfFile)) THEN
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NatmsLoc,NMax12+1/))
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
        CALL New(Top13,(/NatmsLoc,NMax13+1/))
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
        CALL New(Top14,(/NatmsLoc,NMax14+1/))
        CALL Get(Top14,'Top14')
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing Top_14 matrix in Excl_list')
      ENDIF
    ENDIF
!
! Initialize Top_Excl to the Size of Top14 and to zero 
!
    NMax_Excl=NMax14
    CALL New(Top_Excl,(/NatmsLoc,NMax_Excl+1/))
    Top_Excl%I(:,:)=0
!
! Now merge Topologies, in order to avoid double counting in 
! Exclusion energies
!
      NMax_Excl_Out=0
    DO I=1,NatmsLoc
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
      CALL New(Top_Excl_Out,(/NatmsLoc,NMax_Excl_Out+1/))
      Top_Excl_Out%I(1:NatmsLoc,1:NMax_Excl_Out+1)=Top_Excl%I(1:NatmsLoc,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        IF(PRESENT(Tag_O)) THEN
          CALL Put(NMax_Excl_Out,'NMax_Excl14'//TRIM(Tag_O))
          CALL Put(Top_Excl_Out,'TOP_Excl14'//TRIM(Tag_O))
        ELSE
          CALL Put(NMax_Excl_Out,'NMax_Excl14')
          CALL Put(Top_Excl_Out,'TOP_Excl14')
        ENDIF
      ENDIF
    ELSE
      CALL New(Top_New,(/NatmsLoc,NMax_Excl_Out+1/))
      Top_New%I(1:NatmsLoc,1:NMax_Excl_Out+1)=Top_Excl%I(1:NatmsLoc,1:NMax_Excl_Out+1)
      IF(PRESENT(InfFile)) THEN
        IF(PRESENT(Tag_O)) THEN
          CALL Put(NMax_Excl_Out,'NMax_Excl14'//TRIM(Tag_O))
          CALL Put(Top_New,'TOP_Excl14'//TRIM(Tag_O))
        ELSE
          CALL Put(NMax_Excl_Out,'NMax_Excl14')
          CALL Put(Top_New,'TOP_Excl14')
        ENDIF
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
      SUBROUTINE BMatrix(NatmsLoc,XYZ,NIntC,IntCs,B)
!
! This subroutine calculates the sparse, TYPE(BMATR) type 
! representation of the B matrix 
! In current version Cartesian coordinates 
! must be passed in in Angstroems. !!!! Changed, pass in Au-s!!!
! Later use au representation of everything, including B-matrix.
! Linear bendings _must_ always appear in pairs, Defd as LINB1 and LINB2
!
      IMPLICIT NONE
      INTEGER :: NIntC,NIntC2,I,J,K,L,NatmsLoc,IntCoo
      REAL(DOUBLE) :: thresh_B
      REAL(DOUBLE),DIMENSION(1:3,1:NatmsLoc) :: XYZ
      TYPE(INTC) :: IntCs
      TYPE(BMATR):: B
!
      thresh_B=GeOpCtrl%BMatThrsh
!
!     IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
!       CALL INTCValue(IntCs,XYZ)
!       CALL PrtIntCoords(IntCs,IntCs%Value,'Internals in Bmatr gen')
!     ENDIF
!
! allocate B matrix
!
      CALL NEW(B,NIntC)
        B%IB=0
        B%B=Zero
!
      DO IntCoo=1,NIntC     
!
        IF(.NOT.IntCs%Active(IntCoo)) CYCLE
!
        IF(IntCs%Def(IntCoo)(1:4)=='STRE') THEN
! stre
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2)
          CALL STRE(I,J,XYZ(1:3,I),XYZ(1:3,J), &
                    B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5), &
!                      B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
!
        ELSE IF(IntCs%Def(IntCoo)(1:4)=='BEND') THEN
! bend
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) !!! central atom
          K=IntCs%Atoms(IntCoo,3) 
          CALL BEND(I,J,K,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
                    B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5), &
!                      B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
        ELSE IF(IntCs%Def(IntCoo)(1:4)=='OutP') THEN
! out of plane
          I=IntCs%Atoms(IntCoo,1) !!! end atom
          J=IntCs%Atoms(IntCoo,2) !!! central atom
          K=IntCs%Atoms(IntCoo,3) !!! Def plane
          L=IntCs%Atoms(IntCoo,4) !!! Def plane
          CALL OutP(I,J,K,L,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
            XYZ(1:3,L),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12), &
!                      B%B(IntCoo,1:12)
        ELSE IF(IntCs%Def(IntCoo)(1:4)=='TORS') THEN
! torsion of i-j-k-l
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) 
          K=IntCs%Atoms(IntCoo,3) 
          L=IntCs%Atoms(IntCoo,4) 
          CALL TORS(I,J,K,L,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
            XYZ(1:3,L),B%IB(IntCoo,1:12),B%B(IntCoo,1:12),&
            thresh_B,IntCs%Active(IntCoo))
!         write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:12), &
!         B%B(IntCoo,1:12)
        ELSE IF(IntCs%Def(IntCoo)(1:5)=='LINB1') THEN
! linear bendig of i-j-k
          IF(IntCs%Def(IntCoo+1)(1:5)/='LINB2') &
             CALL Halt('LINB Definitions are not paired!')
          I=IntCs%Atoms(IntCoo,1)
          J=IntCs%Atoms(IntCoo,2) 
          K=IntCs%Atoms(IntCoo,3) 
          CALL LINB(I,J,K,XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
            B%IB(IntCoo,1:12),B%B(IntCoo,1:12), &
            B%IB(IntCoo+1,1:12),B%B(IntCoo+1,1:12),thresh_B)
!         write(*,100) IntCs%Def(IntCoo)(1:5), &
!           B%IB(IntCoo,1:12),B%B(IntCoo,1:12)
!         write(*,100) IntCs%Def(IntCoo)(1:5), &
!           B%IB(IntCoo+1,1:12),B%B(IntCoo+1,1:12)
        ELSE IF(IntCs%Def(IntCoo)(1:5)=='LINB2') THEN
          CYCLE
        ELSE IF(IntCs%Def(IntCoo)(1:5)=='CARTX' ) THEN
          I=IntCs%Atoms(IntCoo,1)
          CALL BCART(I,'X',B%IB(IntCoo,1:12),B%B(IntCoo,1:12))
        ELSE IF(IntCs%Def(IntCoo)(1:5)=='CARTY' ) THEN
          I=IntCs%Atoms(IntCoo,1)
          CALL BCART(I,'Y',B%IB(IntCoo,1:12),B%B(IntCoo,1:12))
        ELSE IF(IntCs%Def(IntCoo)(1:5)=='CARTZ' ) THEN
          I=IntCs%Atoms(IntCoo,1)
          CALL BCART(I,'Z',B%IB(IntCoo,1:12),B%B(IntCoo,1:12))
        ENDIF
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
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,RIJ
      REAL(DOUBLE) :: dijsq,t,dij,thresh_B   
      INTEGER :: I,J,K,L,M,N
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
      DijSq = 0.0d0
!  
      DO M=1,3
         T=XJ(M)-XI(M)
         Rij(M)=T
         DijSq=DijSq+T*T
      ENDDO 
!
      Dij = SQRT(DijSq)
!
      DO M=1,3
         T=Rij(M)
         IF(ABS(T)>Thresh_B) THEN
            IBB(M)=(3*(I-1)+M) 
            IBB(3+M)=(3*(J-1)+M) 
            T=T/Dij
            BB(M)=-T
            BB(3+M)=T
         ENDIF   
      ENDDO 
!
      END SUBROUTINE STRE
!
!----------------------------------------------------------------
!
      SUBROUTINE Bend(I,J,K,XI,XJ,XK,IBB,BB,thresh_B)
!
!  This subroutine computes the b matrix elements of a valence
!  angle bending coordinate as Defined by wilson.
!  i and k are the numbers of the end Atoms.
!  j is the number of the central atom.
!  noat: number of Atoms
!  ic: serial number of internal coordinate
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,RJI,RJK,EJI,EJK
      REAL(DOUBLE) :: DIJSq,Tp,T,DJKSq,DxSq,Dx,DOtj,DJK,DJI
      REAL(DOUBLE) :: SinJ,SMI,SMK,SUM,thresh_B
      INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
!
      DijSq=Zero 
      DjkSq=Zero 
      DO 20 M=1,3
        Tp=XJ(M)
        T=XI(M)-Tp
        Rji(M)=T
        DijSq=DijSq+T*T
        T=XK(M)-Tp
        Rjk(m)=T
        DjkSq=DjkSq+T*T
   20 Continue
        Dji=SQRT(DIJSq)
        Djk=SQRT(DJKSq)
        DX=One 
        DotJ=Zero 
        DO 30 M=1,3
        T=RJI(M)/DJI
        Eji(M)=T
        Tp=Rjk(M)/Djk
        Ejk(M)=Tp
        DotJ=DotJ+T*Tp
   30 CONTINUE      
      IF(ABS(DotJ) > 0.99995D0) GO TO 60  !!!! colinearity
      SinJ=SQRT(One-DotJ*DotJ)
      II=3*(I-1)
      JJ=3*(J-1)
      KK=3*(K-1)
      DO 40 M=1,3
      SMI=Dx*(DotJ*Eji(M)-Ejk(M))/(Dji*SinJ)
      IF(ABS(SMI)>Thresh_B) THEN
        IBB(M) = II+M
        BB(M) = SMI
      ENDIF
      SMK=DX*(DotJ*Ejk(M)-Eji(M))/(DJK*SinJ)
      IF(ABS(SMK)>Thresh_B) THEN
        IBB(6+M) = KK+M
        BB(6+M) = SMK
      ENDIF
      Sum=SMI+SMK
      IF(ABS(Sum)>Thresh_B) THEN
        IBB(3+M) = JJ+M
        BB(3+M) = -Sum
      ENDIF
   40 CONTINUE
      RETURN
   50 IER=1
      RETURN
   60 IER=-1
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
      DijSq=Zero
      DjkSq=Zero
      DjlSq=Zero
      DO 20 M=1,3
        Tp=XJ(M)
        T=XI(M)-Tp
        Rji(M)=T
        DijSq=DijSq+T*T
        T=XK(M)-Tp
        Rjk(M)=T
        DjkSq=DjkSq+T*T
        T=XK(M)-XJ(M)
        Rjl(M)=T
        DjlSq=DjlSq+T*T
   20 CONTINUE      
      Dji=SQRT(DijSq)
      Djk=SQRT(DjkSq)
      Djl=SQRT(DjlSq)
      Dx=One 
      CosI=Zero 
      CosK=Zero 
      CosL=Zero 
      DO 30 M=1,3
        T=Rji(M)/Dji
        Eji(M)=T
        Tp=Rjk(M)/Djk
        Ejk(M)=Tp
        Tpp=Rjl(M)/Djl
        Ejl(M)=Tpp
        CosI=CosI+Tp*Tpp
        CosK=CosK+T*Tpp
        CosL=CosL+T*Tp
   30 CONTINUE         
      IF(ABS(CosI) > 0.99995D0) GO TO 70
      SinSin=One-CosI*CosI
      SinI=SQRT(SinSin)
      C1(1)=Ejk(2)*Ejl(3)-Ejk(3)*Ejl(2)
      C1(2)=Ejk(3)*Ejl(1)-Ejk(1)*Ejl(3)
      C1(3)=Ejk(1)*Ejl(2)-Ejk(2)*Ejl(1)
      Dot=Eji(1)*C1(1)+Eji(2)*C1(2)+Eji(3)*C1(3)
      SinT=Dot/SinI
!     if(dabs(sint).gt.0.00005d0) then
!        write(6,1020)nob
!        write(6,'(15x,a36)')intch(nob)
!     end if
      IF(ABS(SinT)>0.99995D0) GO TO 80 !!!! sint show OutP angle
      CosT=SQRT(One-SinT*SinT)
      TanT=SinT/CosT
      II=3*(I-1)
      JJ=3*(J-1)
      KK=3*(K-1)
      LL=3*(L-1)
      CosSin=CosT*SinI
      DO 50 M=1,3
        T=C1(M)/CosSin
        SMI=(T-TanT*Eji(M))/Dji
        IF(ABS(SMI)>Thresh_B) THEN
          IBB(M) = II+M
          BB(M) = Dx*SMI
        ENDIF
        SMK=T*(CosI*CosK-CosL)/(SinSin*Djk)
        IF(ABS(SMK)>Thresh_B) THEN
          IBB(6+M) = KK+M
          BB(6+M) = Dx*SMK
        ENDIF
        SML=T*(CosI*CosL-CosK)/(SinSin*Djl)
        IF(ABS(SML)>Thresh_B) THEN
          IBB(9+M) = LL+M
          BB(9+M) = Dx*SML
        ENDIF
        Sum=SMI+SMK+SML
        IF(ABS(Sum)>Thresh_B) THEN
          IBB(6+M) = JJ+M
          BB(6+M) = -Dx*Sum
        ENDIF
   50 CONTINUE
        RETURN
   60 CONTINUE
        IEr=1
        RETURN
   70 CONTINUE
        IEr=-1
        WRITE(6,1000)
        RETURN
   80 CONTINUE
      IEr=-1
      WRITE(6,1010)
      RETURN
!
 1000 format(/,2x,'<!> k-j-l is colinear (no plane Defined for wag of i)')
 1010 format(/,2x,'<!> i is perpendicular to j-k-l plane - use valence angle bends')
 1020 format(/,2x,'<!> warning: wag of a non-planar system at internal',/,15x,'coordinate no.',i5,':')
!
      END SUBROUTINE OutP
!
!----------------------------------------------------------------------
!
      SUBROUTINE Tors(I,J,K,L,XI,XJ,XK,XL,IBB,BB,thresh_B,Active)
!...
!...  Calculate the B matrix elements for torsion as Defined by
!...  R.L. Hildebrandt, J. Mol. Spec., 44 (1972) 599.
!...
!...  i-j-k-l : torsion around j-k bond
!... 
!
      IMPLICIT NONE
      REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,XL,rij,rjk,rlk
      REAL(DOUBLE),DIMENSION(1:3) :: eij,ejk,elk,cr,sj,sk
      REAL(DOUBLE) :: SMI,SMK,Thresh_B,dij,dijsq,dx,DJK,t,dxsq,DJKSq
      REAL(DOUBLE) :: cosj,sin2j,smj,dlksq,dlk,cosk,sin2k,sml 
      INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
      INTEGER :: IBB(1:12)
      REAL(DOUBLE) :: BB(1:12)
      LOGICAL      :: Active
!
      DjkSq=Zero 
      DxSq=Zero 
!
      DO M=1,3
         Sj(M) = Zero 
         Sk(M) = Zero 
         T=XK(M)-XJ(M)
         Rjk(M)=T
         DjkSq=DjkSq+T*T
      ENDDO 
!
      DJK=One/SQRT(DjkSq)
      Dx=One  
!
      DO M=1,3
        Ejk(M)=Rjk(M)*Djk
      ENDDO 
!
      JJ=3*(J-1)
      KK=3*(K-1)
!
!...
        DijSq=Zero 
      DO 40 M=1,3
        T=XJ(M)-XI(M)
        Rij(M)=T
        DijSq=DijSq+T*T
   40 CONTINUE        
        Dij=One/SQRT(DijSq)
        CosJ=Zero
      DO 50 M=1,3
        T=Rij(M)*Dij
        Eij(M)=T
        CosJ=CosJ-T*Ejk(M)
   50 CONTINUE
      IF(ABS(CosJ)>0.99995d0) GO TO 120
      Sin2J=(One-CosJ*CosJ)
      II=3*(I-1)
      Cr(1)=Eij(2)*Ejk(3)-Eij(3)*Ejk(2)
      Cr(2)=Eij(3)*Ejk(1)-Eij(1)*Ejk(3)
      Cr(3)=Eij(1)*Ejk(2)-Eij(2)*Ejk(1)
      DO 60 M=1,3
        T=Cr(M)/Sin2J
        SMI=T*Dij
      IF(ABS(SMI)>Thresh_B) THEN
        IBB(M) = II+M
        BB(M) = -Dx*SMI
      ENDIF
        SMK=T*CosJ*DJK
        SK(M)=SK(M)+SMK
        SMJ=SMI-SMK
        SJ(M)=SJ(M)+SMJ
   60 CONTINUE
!
      DlkSq=Zero
      DO 70 M=1,3
        T=XK(M)-XL(M)
        Rlk(M)=T
        DlkSq=DlkSq+T*T
   70 CONTINUE          
      Dlk=One/SQRT(DlkSq)
      CosK=Zero
      DO 80 M=1,3
        T=Rlk(M)*Dlk
        Elk(m)=T
      CosK=CosK+Ejk(M)*T
   80 CONTINUE
      IF(ABS(CosK)>0.99995d0) GO TO 120
      Sin2K=(One-CosK*CosK)
      LL=3*(L-1)
      Cr(1)=Elk(3)*Ejk(2)-Elk(2)*Ejk(3)
      Cr(2)=Elk(1)*Ejk(3)-Elk(3)*Ejk(1)
      Cr(3)=Elk(2)*Ejk(1)-Elk(1)*Ejk(2)
      DO 90 M=1,3
        T=Cr(M)/Sin2K
        SML=T*Dlk
        IF(ABS(SML)>Thresh_B) THEN
          IBB(9+M) = LL+M
          BB(9+M) = -Dx*SML
        ENDIF
        SMJ=T*CosK*DJK
        SJ(M)=SJ(M)+SMJ
        SMK=SML-SMJ
        SK(M)=SK(M)+SMK
   90 CONTINUE
      DO 100 M=1,3
        SMJ=SJ(M)
      IF(ABS(SMJ)>Thresh_B) THEN
        IBB(3+M) = JJ+M
        BB(3+M) = Dx*SMJ
      ENDIF
        SMK=SK(M)
      IF(ABS(SMK)>Thresh_B) THEN
        IBB(6+M) = KK+M
        BB(6+M) = Dx*SMK
      ENDIF
  100 CONTINUE
      RETURN
  110 CONTINUE
      IEr=1
      RETURN
  120 CONTINUE
      IEr=-1
      WRITE(6,1030)
! colinear bonds in torsion, delete the corresponding row of the B matrix!
        IBB(:)=0    
        BB(:)=Zero
        Active=.FALSE.
        RETURN
!
1030 format(' i-j-k or j-k-l is colinear (no torsion possibble), this row of the B matrix will be deleted.')
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
      REAL(DOUBLE),DIMENSION(1:3) :: RJKN
      REAL(DOUBLE) :: t,SUM,DIJSq,DJKSq,djasq,Tp,DJI,DJK,dja,dx
      REAL(DOUBLE) :: DotJ,DotP,Test,SMI,SMK,SumX
      REAL(DOUBLE) :: thresh_B
      INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER,IntCoo
      INTEGER :: IBB1(1:12),IBB2(1:12)
      REAL(DOUBLE) :: BB1(1:12),BB2(1:12)
!
! Set up coordinates of point a
!
        SUM=0.D0
        RJI=XI-XJ
        RJK=XK-XJ
        Sum=DOT_PRODUCT(RJI,RJI)+DOT_PRODUCT(RJK,RJK)
        Sum=SQRT(Sum)
! construct orthogonal vector to the plane of IJK 
        CALL CROSS_PRODUCT(RJI,RJK,RJKN)
        SumX=DOT_PRODUCT(RJKN,RJKN)
! form arbitrary vector rjk which is no more parallel with jk bond
        IF(SumX<1.D-5) THEN
          RJKN(1)=RJK(1)+SUM 
          RJKN(2)=RJK(2)+SUM 
          RJKN(3)=RJK(3)-SUM 
        ENDIF
! form vector a, perpendicular to the plane of rjkn and RJI
        CALL CROSS_PRODUCT(RJKN,RJI,A)
        Sum=DOT_PRODUCT(A,A)        
        Sum=SQRT(Sum)
        A(1:3)=A(1:3)/Sum
! generate point a
        A(1:3)=XJ(1:3)+A(1:3)        
!
      DIJSq=Zero
      DJKSq=Zero
      DjaSq=Zero
      DO 20 M=1,3
        Tp=XJ(M)
        T=XI(M)-Tp
        RJI(M)=T
        DijSq=DijSq+T*T
        T=XK(M)-Tp
        Rjk(M)=T
        DjkSq=DjkSq+T*T
        T=A(M)-Tp
        UN(M)=T
        DjaSq=DjaSq+T*T
   20 CONTINUE
      Dji=SQRT(DijSq)
      Djk=SQRT(DjkSq)
      dja=SQRT(DjaSq)
      Dx=One 
      DotJ=Zero 
      DotP=Zero 
      DO 30 M=1,3
        T=Rji(M)/Dji
        Tp=Rjk(M)/Djk
        Ejk(M)=Tp
        DotJ=DotJ+T*Tp
        Tp=UN(M)/Dja
        Unit(M)=Tp
        DotP=DotP+T*Tp
   30 CONTINUE       
      Test=ABS(DotJ)-One 
      IF(ABS(Test)>0.00005d0) GO TO 70
      IF(ABS(DotP)>0.00005d0) GO TO 80
      II=3*(I-1)
      JJ=3*(J-1)
      KK=3*(K-1)
      DO 40 M=1,3
        T=Unit(M)
        IF(ABS(T)<Thresh_B) GO TO 40
        T=-Dx*T
        SMI=T/Dji
        IBB1(M)=II+M
        BB1(M)=SMI
        SMK=T/Djk
        IBB1(3+M) = jj+M
        BB1(3+M) = -SMI-SMK
        IBB1(6+M) = KK+M
        BB1(6+M) = SMK
   40 CONTINUE
      Up(1)=Ejk(2)*Unit(3)-Ejk(3)*Unit(2)
      Up(2)=Ejk(3)*Unit(1)-Ejk(1)*Unit(3)
      Up(3)=Ejk(1)*Unit(2)-Ejk(2)*Unit(1)
      DO 50 M=1,3
        T=Up(M)
        IF(ABS(T)<Thresh_B) GO TO 50
        T=-Dx*T
        SMI=T/Dji
        IBB2(M)=II+M
        BB2(M)=SMI
        SMK=T/Djk
        IBB2(3+M)=JJ+M
        BB2(3+M)=-SMI-SMK
        IBB2(6+M)=KK+M
        BB2(6+M)=SMK
   50 CONTINUE
      RETURN
   60 CONTINUE
      IEr=1
      RETURN
   70 CONTINUE
      IEr=-1
      WRITE(6,1020)
      RETURN
   80 CONTINUE
      IEr=-1
      WRITE(6,1030)
      RETURN
!
 1010 FORMAT('+',86x,'a = (',2(f11.7,','),f11.7,')')
 1020 FORMAT(' i-j-k not colinear - use valence angle bend')
 1030 FORMAT(' atom a not perpendicular to i-j-k at j')
!
      END SUBROUTINE LinB
!
!----------------------------------------------------------------
!
      SUBROUTINE BCART(I,CHAR,IB,B)
!
      INTEGER :: I,J
      INTEGER,DIMENSION(1:12)      :: IB
      REAL(DOUBLE),DIMENSION(1:12) :: B
      CHARACTER :: CHAR
!
        J=(I-1)*3
	IB=0
	 B=Zero
!
        IF(CHAR=='X') THEN
	  IB(1)=J+1
	   B(1)=One
	ENDIF
!
        IF(CHAR=='Y') THEN
	  IB(2)=J+2
	   B(2)=One
	ENDIF
!
        IF(CHAR=='Z') THEN
	  IB(3)=J+3
	   B(3)=One
	ENDIF
!
      END SUBROUTINE BCART
!
!----------------------------------------------------------------
!
      SUBROUTINE DefineIntCoos(NatmsLoc,XYZ,MMAtNum,InfFile,IntSet,IntCs,NIntC)
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
      INTEGER :: I,J,N,NatmsLoc,NBond,NIntC,IntSet
      INTEGER :: NMax_Excl,NMax12,ILast
      REAL(DOUBLE),DIMENSION(1:3,1:NatmsLoc) :: XYZ
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
   CALL SORT_INTO_Box1(BoxSize,XYZ,NatmsLoc,&
                   NX,NY,NZ,BXMIN,BYMIN,BZMIN)
!
   NBox=NX*NY*NZ
   CALL New(BoxI,NBox+1)
   CALL New(BoxJ,NatmsLoc)
!
   CALL SORT_INTO_Box2(BoxSize,XYZ,NatmsLoc,&
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
      CALL BondList(NatmsLoc,XYZ,NBond,MMAtNum, &
           BoxI,BoxJ,NBox,NX,NY,NZ,CritRad)
      CALL New(BondIJ,(/2,NBond/))
      CALL BondList(NatmsLoc,XYZ,NBond,MMAtNum, &
           BoxI,BoxJ,NBox,NX,NY,NZ,CritRad,BondIJ)
!
      IF(IntSet==1) THEN
!
! Now define covalent topology matrices
!
        CALL Topology_12(NatmsLoc,NBond,BondIJ%I(1,1:NBond),BondIJ%I(2,1:NBond),Top12,InfFile=InfFile)
        CALL Topology_13(NatmsLoc,Top12,Top13)
        CALL Topology_14(NatmsLoc,Top12,Top14)
        CALL Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl,InfFile=InfFile)
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
        CALL New(Top_Excl,(/NatmsLoc,NMax_Excl+1/))
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
        CALL AngleList(NatmsLoc,Top12,NAngle=NAngle)
        CALL New(AngleIJK,(/3,NAngle/))
        CALL AngleList(NatmsLoc,Top12,AngleIJK,NAngle)
!
        CALL TorsionList(NatmsLoc,Top12,NTorsion=NTorsion)
        CALL New(TorsionIJKL,(/4,NTorsion/))
        CALL TorsionList(NatmsLoc,Top12,TorsionIJKL,NTorsion)
!
        CALL Delete(Top12)
!
      ELSE
!
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NatmsLoc,NMax12+1/))
        CALL Get(Top12,'TOP12')
!
        CALL VDWAngleList(NatmsLoc,Top12,AngleIJK,NAngle,BondIJ,NBond)
!
        CALL VDWTorsionList(NatmsLoc,Top12,TorsionIJKL,NTorsion,BondIJ,NBond)
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
      IF(NIntC==0) GO TO 300
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
300   CONTINUE
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
      SUBROUTINE BondList(NatmsLoc,XYZ,NBond,MMAtNum,&
             BoxI,BoxJ,NBox,NX,NY,NZ,CritRad,BondIJ)
!
! Define number of bonds
!
      IMPLICIT NONE
      INTEGER :: I,J,NatmsLoc,NBond
      REAL(DOUBLE),DIMENSION(1:3,1:NatmsLoc) :: XYZ
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
SUBROUTINE AngleList(NatmsLoc,Top12,AngleIJK,NAngle)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2) :: Top12
TYPE(INT_RNK2),OPTIONAL :: AngleIJK 
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ
INTEGER :: NatmsLoc,KK,IN12,JN12,NAngle
LOGICAL :: AngleFill
!
            AngleFill=.FALSE.
            IF(PRESENT(AngleIJK)) AngleFill=.TRUE.
            NAngle=0
    DO II=1,NatmsLoc
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
SUBROUTINE TorsionList(NatmsLoc,Top12,TorsionIJKL,NTorsion)
! Set up a table which shows the atom numbers of Atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2):: Top12
TYPE(INT_RNK2),OPTIONAL :: TorsionIJKL
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL,NTorsion
INTEGER :: NatmsLoc,IN12,JN12,KN12
LOGICAL :: TorsionFill
!
            TorsionFill=.FALSE.
            IF(PRESENT(TorsionIJKL)) TorsionFill=.TRUE.
            NTorsion=0
    DO II=1,NatmsLoc
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
SUBROUTINE VDWAngleList(NatmsLoc,Top12,AngleIJK,NAngle,BondIJ,NBond)
! this routine generates bond-angles associated with WDV bonds
!
    IMPLICIT NONE
    INTEGER :: NatmsLoc,I,I1,I2,N1,N2,J,J1,J2,NBond,NAngle
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
SUBROUTINE VDWTorsionList(NatmsLoc,Top12,TorsionIJKL,NTorsion,BondIJ,NBond)
! this routine generates bond-angles associated with WDV bonds
!
    IMPLICIT NONE
    INTEGER :: NatmsLoc,I,I1,I2,J,J1,J2,N1,N2,NBond,NTorsion,K1,K2
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
SUBROUTINE GetIntCs(XYZ,NatmsLoc,InfFile,IntCs,NIntC,Refresh)
!
! This subroutine constructs the IntCs array, which holds
! definitions of internal coordinates to be used in the 
! forthcoming geometry manipulation procedure.
! Refresh=1 : Refresh all definitions
!        =2 : Refresh only definitions based on VDW interaction
!        =3 : Do not refresh definitions, use the one from HDF
!        =4 : Refresh/generate only the covalent coordinates
!        =5 : only the extra coordinates from input
! WARNING! In the present version bending -> linear bending transitions are 
! always checked and refreshed
! 
      IMPLICIT NONE
      TYPE(INTC) :: IntCs,IntC_Cov,IntC_VDW,IntC_Extra,IntC_New
      INTEGER :: NIntC,NIntC_Cov,NIntC_VDW,NIntC_Extra,NNew,Nintc_New
      INTEGER :: I,J,K,Refresh,NatmsLoc,II,ILast
      INTEGER :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
      INTEGER :: NStreGeOp,NBendGeOp,NLinBGeOp,NOutPGeOp,NTorsGeOp
      TYPE(INT_VECT) :: MMAtNum,LinAtom,MarkLinb
      TYPE(DBL_VECT) :: NuclCharge,AuxVect
      TYPE(INT_RNK2) :: LinBBridge,Top12
      TYPE(CRDS) :: GMLoc
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: InfFile
      REAL(DOUBLE),DIMENSION(1:3,1:NatmsLoc) :: XYZ
      REAL(DOUBLE) :: Value
!
! Get atomnames (numbers) from HDF 
!
      CALL New(MMAtNum,NatmsLoc)
#ifdef MMech
      IF(HasMM()) THEN
        CALL Get(MMAtNum,'MMATNUM')
      ELSE
        CALL New(NuclCharge,NatmsLoc)
        CALL Get(NuclCharge,  'atomicnumbers',Tag_O=CurGeom)
        DO I=1,NatmsLoc ; MMAtNum%I(I)=INT(NuclCharge%D(I)) ; ENDDO
        CALL Delete(NuclCharge)
      ENDIF
#else
        CALL New(NuclCharge,NatmsLoc)
        CALL Get(NuclCharge,  'atomicnumbers',Tag_O=CurGeom)
        DO I=1,NatmsLoc ; MMAtNum%I(I)=INT(NuclCharge%D(I)) ; ENDDO
        CALL Delete(NuclCharge)
#endif

!
      IF(Refresh==1) Then !!! Total refresh
!define covalent bonding scheme
          CALL DefineIntCoos(NatmsLoc,XYZ,MMAtNum,InfFile,1, &
                             IntC_Cov,NIntC_Cov)
!define Van der Waals bonding scheme
          CALL DefineIntCoos(NatmsLoc,XYZ,MMAtNum,InfFile,2, &
                             IntC_VDW,NIntC_VDW)
          CALL Put(NIntC_Cov,'NIntC_Cov')
          IF(NIntC_Cov/=0) CALL Put(IntC_Cov,'IntC_Cov')
          CALL Put(NIntC_VDW,'NIntC_VDW')
          IF(NIntC_VDW/=0) CALL Put(IntC_VDW,'IntC_VDW')
!
      ELSE IF(Refresh==2) THEN !!! refresh VDW terms
!
! Refresh only the VDW terms
!
          CALL Get(NIntC_Cov,'NIntC_Cov')
          CALL New(IntC_Cov,NIntC_Cov)
          CALL Get(IntC_Cov,'IntC_Cov')
          CALL DefineIntCoos(NatmsLoc,XYZ,MMAtNum,InfFile,2, &
                             IntC_VDW,NIntC_VDW)
          CALL Put(NIntC_VDW,'NIntC_VDW')
          IF(NIntC_VDW/=0) CALL Put(IntC_VDW,'IntC_VDW')
!
      ELSE IF(Refresh==3) THEN !!! no refresh, get everything from HDF
          CALL Get(NIntC_Cov,'NIntC_Cov')
          CALL New(IntC_Cov,NIntC_Cov)
          CALL Get(IntC_Cov,'IntC_Cov')
          CALL Get(NIntC_VDW,'NIntC_VDW')
          CALL New(IntC_VDW,NIntC_VDW)
          CALL Get(IntC_VDW,'IntC_VDW')
!
      ELSE IF(Refresh==4) THEN !!! refresh covalent bonds only
        CALL DefineIntCoos(NatmsLoc,XYZ,MMAtNum,InfFile,1, &
                           IntC_Cov,NIntC_Cov)
        NIntC_VDW=0 
!
      ELSE IF(Refresh==5) THEN !!! use only extra coords from input
        NIntC_Cov=0 
        NIntC_VDW=0 
      ENDIF
!
! Put Covalent terms
!
        CALL Put(NIntC_Cov,'NIntC_Cov')
        IF(NIntC_Cov/=0) CALL Put(IntC_Cov,'IntC_Cov')
        CALL Put(NIntC_VDW,'NIntC_VDW')
        IF(NIntC_VDW/=0) CALL Put(IntC_VDW,'IntC_VDW')
!
! Get extra internal coordinates and constraints
!
        CALL Get(NIntC_Extra,'NIntC_Extra')
        IF(NIntC_Extra/=0) THEN
          CALL New(IntC_Extra,NIntC_Extra)
          CALL Get(IntC_Extra,'IntC_Extra')
            IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
              CALL PrtIntCoords(IntC_Extra, &
                 IntC_Extra%Value,'Extra Internals after Read')
            ENDIF
        ENDIF
!
! Merge INTC arrays
!
        NIntC=NIntC_Cov+NIntC_VDW+NIntC_Extra
!
        IF(AllocQ(IntCs%Alloc)) CALL Delete(IntCs)
        CALL New(IntCs,NIntC)
!
          ILast=0
        IF(NIntC_Cov/=0) THEN
          CALL Set_INTC_EQ_INTC(IntC_Cov,IntCs,1,NIntC_Cov,ILast+1)
        ENDIF
          ILast=NIntC_Cov
        IF(NIntC_VDW/=0) THEN 
          CALL Set_INTC_EQ_INTC(IntC_VDW,IntCs,1,NIntC_VDW,ILast+1)
        ENDIF
!
        IF(NIntC_Extra/=0) THEN
          ILast=ILast+NIntC_VDW
          CALL Set_INTC_EQ_INTC(IntC_Extra,IntCs,1,NIntC_Extra,ILast+1)
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
! keep those values, which were defd. in extras.
! Do not check for Cartesians, since they
! are not supposed to be generated by anything, but
! EXTRA-s.
!
        ILast=NIntC_Cov+NIntC_VDW
100     II=0
        DO I=ILast+1,ILast+NIntC_Extra
            IF(IntCs%Def(J)=='CART ') CYCLE
          DO J=1,ILast
            IF(IntCs%Def(J)=='BLANK') CYCLE
            IF(IntCs%Atoms(J,1)==IntCs%Atoms(I,1).AND.&
               IntCs%Atoms(J,2)==IntCs%Atoms(I,2).AND.&
               IntCs%Atoms(J,3)==IntCs%Atoms(I,3).AND.&
               IntCs%Atoms(J,4)==IntCs%Atoms(I,4)) THEN
               IntCs%Def(J)='BLANK'
               IntCs%Atoms(J,1:4)=0      
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
              CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NNew)
            ENDIF
          ENDDO
          CALL Delete(IntCs)
          NIntC=NNew
          CALL New(IntCs,NIntC)
            CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntC,1)
          CALL Delete(IntC_New)
        ENDIF
!
      CALL Delete(MMAtNum)
!
! Check for bending - lin.bending transions
! also for long range torsions
!
      CALL ChkBendToLinB(IntCs,NIntC,XYZ)
!
! Print actual set of internals for testing
!
      IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
!       CALL INTCValue(IntCs,XYZ)
        CALL PrtIntCoords(IntCs, &
          IntCs%Value,'internals after ChkBendToLinB')
      ENDIF 
!
! Set active all internal coords defd so far.
! 'Linear torsions will be deactivated when their value is calculated
! by some subroutine (eg. INTCValue).
! The set of active coordinates may vary during the process of optimization.
!
      IntCs%Active=.TRUE.
!
! Save values of constraints into HDF
!
      CALL New(AuxVect,NIntC)
        AuxVect%D=IntCs%Value
        CALL Put(Current(3),'LastIntcGeom')
        CALL Put(AuxVect,'Constraints'//TRIM(CurGeom))
      CALL Delete(AuxVect)
!
! Count number of different internal coord types
!
      NStreGeOp=0;NBendGeOp=0;NLinBGeOp=0;NOutPGeOp=0;NTorsGeOp=0
        DO I=1,NIntC
           IF(IntCs%Def(I)(1:4)=='STRE') THEN
             NStreGeOp=NStreGeOp+1
           ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
             NBendGeOp=NBendGeOp+1
           ELSE IF(IntCs%Def(I)(1:4)=='LINB') THEN
             NLinBGeOp=NLinBGeOp+1
           ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
             NOutPGeOp=NOutPGeOp+1
           ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
             NTorsGeOp=NTorsGeOp+1
           ENDIF
        ENDDO
      CALL Put(NStreGeOp,'NStreGeOp')
      CALL Put(NBendGeOp,'NBendGeOp')
      CALL Put(NLinBGeOp,'NLinBGeOp')
      CALL Put(NOutPGeOp,'NOutPGeOp')
      CALL Put(NTorsGeOp,'NTorsGeOp')
!
! Save current internals to HDF
!
!     CALL Put(NIntC,'NIntC')
!     CALL Put(IntCs,'IntCs')
!
END SUBROUTINE GetIntCs
!-------------------------------------------------------
!
SUBROUTINE INTCValue(IntCs,XYZ)
!
! Determine value of internal coordinates.
! Input coordintes are now in atomic units!
! 
    IMPLICIT NONE
    TYPE(INTC) :: IntCs
    INTEGER :: NIntCs,I,J,K,L,I1,I2,I3,I4,NatmsLoc
    REAL(DOUBLE),DIMENSION(:,:) :: XYZ
!
    NIntCs=SIZE(IntCs%Def(:))
    NatmsLoc=SIZE(XYZ,2)
!
    DO I=1,NIntCs
      IF(IntCs%Def(I)(1:4)=='STRE') THEN
        I1=IntCs%Atoms(I,1)
        I2=IntCs%Atoms(I,2)
        CALL STREValue(XYZ(1:3,I1),XYZ(1:3,I2),IntCs%Value(I))
!       IntCs%Value(I)=IntCs%Value(I)/AngstromsToAU 
!
      ELSE IF(IntCs%Def(I)(1:4)=='BEND'.OR.IntCs%Def(I)(1:4)=='LINB') THEN
        I1=IntCs%Atoms(I,1)
        I2=IntCs%Atoms(I,2)
        I3=IntCs%Atoms(I,3)
        CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),IntCs%Value(I))
!
      ELSE IF(IntCs%Def(I)(1:4)=='OutP'.OR.IntCs%Def(I)(1:4)=='TORS') THEN
        I1=IntCs%Atoms(I,1)
        I2=IntCs%Atoms(I,2)
        I3=IntCs%Atoms(I,3)
        I4=IntCs%Atoms(I,4)
        CALL TORSValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),XYZ(1:3,I4),IntCs%Value(I),IntCs%Def(I),IntCs%Active(I))
!
      ELSE IF(IntCs%Def(I)(1:5)=='CARTX') THEN
        I1=IntCs%Atoms(I,1)
        IntCs%Value(I)=XYZ(1,I1)
      ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
        I1=IntCs%Atoms(I,1)
        IntCs%Value(I)=XYZ(2,I1)
      ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
        I1=IntCs%Atoms(I,1)
        IntCs%Value(I)=XYZ(3,I1)
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
      IF(DABS(DotJ)>0.999999D0) THEN
        DotJ=SIGN(DABS(DotJ),DotJ)
      ENDIF
      IF(ABS(DotJ-One)<1.D-6) THEN
        AIJK=0.D0
      ELSE IF(ABS(DotJ+One)<1.D-6) THEN
        AIJK=PI  
      ELSE
        AIJK=ACOS(DotJ)
      ENDIF
!
END SUBROUTINE BENDValue
!
!-------------------------------------------------------------
!
      SUBROUTINE TORSValue(XIin,XJin,XKin,XLin,VTors,Def,Active)
!
      REAL(DOUBLE),DIMENSION(1:3) :: XIin,XJin,XKin,XLin
      TYPE(DBL_VECT)              :: XI,XJ,XK,XL,V1,V2
      TYPE(DBL_RNK2)              :: Rot
      REAL(DOUBLE)                :: VTors,CosPhi,Sum,V1ABS,V2ABS
      CHARACTER(LEN=5)            :: Def
      INTEGER                     :: I,J,K,L
      LOGICAL                     :: Active
!
      CALL New(XI,3)
      CALL New(XJ,3)
      CALL New(XK,3)
      CALL New(XL,3)
      CALL New(Rot,(/3,3/))
      CALL New(V1,3)
      CALL New(V2,3)
!
      XI%D=XIin
      XJ%D=XJin
      XK%D=XKin
      XL%D=XLin
!
! Check for linearity
!
      V1%D=XI%D-XJ%D
      V2%D=XJ%D-XK%D
      V1ABS=SQRT(DOT_PRODUCT(V1%D,V1%D))
      V2ABS=SQRT(DOT_PRODUCT(V2%D,V2%D))
      IF(ABS(V1ABS)<1.D-3.OR.ABS(V2ABS)<1.D-3) THEN
        Active=.FALSE. 
        GO TO 1000
      ENDIF
      SUM=ABS(DOT_PRODUCT(V1%D,V2%D))/V1ABS/V2ABS
      SUM=ACOS(SUM)*180.D0/PI
      IF(Sum<LinCrit) THEN
        Active=.FALSE. 
        GO TO 1000
      ENDIF
!
      V1%D=XL%D-XK%D
      V2%D=XK%D-XJ%D
      V1ABS=SQRT(DOT_PRODUCT(V1%D,V1%D))
      V2ABS=SQRT(DOT_PRODUCT(V2%D,V2%D))
      IF(ABS(V1ABS)<1.D-3.OR.ABS(V2ABS)<1.D-3) THEN
        Active=.FALSE. 
        GO TO 1000
      ENDIF
      SUM=ABS(DOT_PRODUCT(V1%D,V2%D))/V1ABS/V2ABS
      SUM=ACOS(SUM)*180.D0/PI
      IF(Sum<LinCrit) THEN
        Active=.FALSE. 
        GO TO 1000
      ENDIF
!
! Translate, so that XK be in the origin
!
      XI%D=XI%D-XK%D 
      XJ%D=XJ%D-XK%D 
      XL%D=XL%D-XK%D 
      XK%D=Zero        
!
! Rotate system, such that J be on Z axis
!
      V1%D=Zero
      V1%D(3)=One 
      V2%D=XJ%D
      CALL Rotate(V1%D,V2%D,Rot%D)
!
      CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,XI%D,V2%D)
      XI%D=V2%D
      CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,XJ%D,V2%D)
      XJ%D=V2%D
      V2%D=Zero
      CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,XL%D,V2%D)
      XL%D=V2%D
!
! Now, rotate system around z axis such, that IJ bond be 
! parallel with x axis, I pointing to positive X
!
      V1%D=Zero
      V1%D(1)=One 
      V2%D=Zero
      V2%D(1)=XI%D(1)
      V2%D(2)=XI%D(2)
      CALL Rotate(V1%D,V2%D,Rot%D)
!
      CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,XI%D,V2%D)
      XI%D=V2%D
      CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,XJ%D,V2%D)
      XJ%D=V2%D
      CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,XL%D,V2%D)
      XL%D=V2%D
!
! Now, calculate torsional angle
!
      V1%D=Zero
      V1%D(1)=One 
      V2%D(1)=XL%D(1)
      V2%D(2)=XL%D(2)
      V2%D(3)=Zero   
      Sum=SQRT(DOT_PRODUCT(V2%D,V2%D))
      V2%D=V2%D/Sum
      CosPhi=DOT_PRODUCT(V1%D,V2%D)
!
      IF(DABS(CosPhi)>0.999999D0) THEN
        CosPhi=SIGN(DABS(CosPhi),CosPhi)
      ENDIF
      IF(ABS(CosPhi-One)<1.D-6) THEN
        VTors=0.D0
      ELSE IF(ABS(CosPhi+One)<1.D-6) THEN
        VTors=PI  
      ELSE
        VTors=ACOS(CosPhi)
      ENDIF
!
! Take orientation of torsional chain into account
!
      IF(XL%D(2)>Zero) VTors=Two*PI-VTors
!
1000  CONTINUE
      CALL Delete(V2)
      CALL Delete(V1)
      CALL Delete(Rot)
      CALL Delete(XI)
      CALL Delete(XJ)
      CALL Delete(XK)
      CALL Delete(XL)
!
      END SUBROUTINE TORSValue
!-------------------------------------------------------
!
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
    REAL(DOUBLE) :: TrixThresh,AInvDistanceThresh,DiffMax,RMSD
    REAL(DOUBLE) :: GrdTrfCrit,CooTrfCrit,ScaleTo,MaxGradDiff,Sum
    INTEGER :: NCart,I,II,J,TrfType,NIntC,DiffLength
    INTEGER :: NVBlocksB,NHBlocksB,BlockGeomSize
    INTEGER :: MaxIt_GrdTrf,MaxIt_CooTrf
    TYPE(INTC) :: IntCs
    TYPE(BMATR):: B
    LOGICAL :: RefreshB,RefreshAct
!
    CALL OpenASCII(OutFile,Out)
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
      CALL INTCValue(IntCs,GMLoc%Carts%D)
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
      CALL INTCValue(IntCs,GMLoc%Carts%D)
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
      IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
        WRITE(*,*) 'Gradient transformation, No. Int. Coords= ',NIntC
        WRITE(Out,*) 'Gradient transformation, No. Int. Coords= ',NIntC
      ENDIF
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
        RMSD=DOT_PRODUCT(SpVectAux%MTrix%D(1:DiffLength),SpVectAux%MTrix%D(1:DiffLength))
        RMSD=SQRT(RMSD/DBLE(NIntC))
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
      WRITE(*,110) II,DiffMax,RMSD
      WRITE(Out,110) II,DiffMax,RMSD
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
      IF(II>0) CALL INTCValue(IntCs,GMLoc%Carts%D)
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
        RMSD=DOT_PRODUCT(SpVectAux%MTrix%D(1:DiffLength),SpVectAux%MTrix%D(1:DiffLength))
        RMSD=SQRT(RMSD/DBLE(NCart))
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
      WRITE(*,210) II,DiffMax,RMSD
      WRITE(Out,210) II,DiffMax,RMSD
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
    Close(Out)
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
!       IF(SUM>PI) SUM=SUM-TwoPi
          VectAux(I)=SUM
      ENDIF
    ENDDO
!
    END SUBROUTINE MapBackAngle
!
!-------------------------------------------------------
!
    SUBROUTINE GDIIS(GMLoc)
    IMPLICIT NONE
    TYPE(CRDS)     :: GMLoc
    INTEGER        :: I,J,K,L,GDIISMemoryIn,DimGDIIS,IGeom,DimOverl
    INTEGER        :: SRMemory,RefMemory,CartGradMemory
    INTEGER        :: LastStruct,ILow,NCart,GDIISMemory
    INTEGER        :: INFO,NewDim 
    TYPE(DBL_RNK2) :: SRDispl,AMat,SRStruct,ActCarts
    TYPE(DBL_RNK2) :: TrfDispl,TrfGrad,EigVect
    TYPE(DBL_RNK2) :: RefGrad
    TYPE(DBL_RNK2) :: AuxStruct,RefStruct
    TYPE(DBL_VECT) :: Coeffs,AuxVect,AuxVect2,EigVal
    TYPE(DBL_VECT) :: GrdDotGrd,GrdDotDX,DXDotDX,MetricA
    TYPE(INT_VECT) :: Selection
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
    REAL(DOUBLE)                   :: Sum,S1,S2
    REAL(DOUBLE)                   :: Eig,MFact
    REAL(DOUBLE)                   :: SumX,SumY,SumZ
    REAL(DOUBLE)                   :: MinGrdDotGrd,MinGrdDotDX
    REAL(DOUBLE)                   :: MinDXDotDX
    REAL(DOUBLE)                   :: EDI,ED0,EDT 
!
! Current implementation runs on Cartesian displacements only
! Input GMLoc contains the actual geometry
! CurGeom is set to the last (simple relaxation) structure
!
      CALL Get(SRMemory,'SRMemory')
      CALL Get(RefMemory,'RefMemory')
      CALL Get(CartGradMemory,'CartGradMemory')
      NCart=3*GMLoc%Natms
      DimGDIIS=NCart   !!! later dimension may become NIntC
      IF(SRMemory/=RefMemory) CALL Halt('SRMemory/=RefMemory in GDIIS')
      IF(SRMemory/=CartGradMemory) CALL Halt('SRMemory/=CartGradMemory in GDIIS')
      GDIISMemory=SRMemory
!
      IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
        CALL OpenAscii(OutFile,Out)
        WRITE(Out,*) 'GDIISMemory= ',GDIISMemory
        WRITE(*,*) 'GDIISMemory= ',GDIISMemory
        CLOSE(Unit=Out,STATUS='KEEP')
      ENDIF
!
    GMTag=''
#ifdef MMech
    IF(HasMM()) GMTag='GM_MM'
#endif
!
! Get Displacements from geometries, stored in HDF
! and fill them into SRDispl columns.
!
      CALL New(SRDispl,(/DimGDIIS,GDIISMemory/))
      CALL New(SRStruct,(/DimGDIIS,GDIISMemory/))
      CALL New(RefStruct,(/DimGDIIS,GDIISMemory/))
      CALL New(RefGrad,(/DimGDIIS,GDIISMemory/))
      SRStruct%D=Zero
      CALL New(AuxVect,DimGDIIS)
      CALL New(GrdDotDX,GDIISMemory)
      CALL New(DXDotDX,GDIISMemory)
      CALL New(GrdDotGrd,GDIISMemory)
      CALL New(Selection,GDIISMemory)
      Selection%I=1
      CALL New(MetricA,GDIISMemory)
      CALL New(Coeffs,GDIISMemory)
!
! Get recent Ref and SR structures, as well as the Gradients
!
      DO IGeom=1,GDIISMemory
        CALL Get(AuxVect,'SR'//TRIM(IntToChar(IGeom)))
        SRStruct%D(:,IGeom)=AuxVect%D
        CALL Get(AuxVect,'Ref'//TRIM(IntToChar(IGeom)))
        RefStruct%D(:,IGeom)=AuxVect%D
        CALL Get(AuxVect,'CartGrad'//TRIM(IntToChar(IGeom)))
        RefGrad%D(:,IGeom)=AuxVect%D
      ENDDO
!
! Get simple relaxation displacements ('error vectors')
! and normalize them!
!
      SRDispl%D=SRStruct%D-RefStruct%D
!
! Now calculate overlap ('A' matrix). 
! Presently only non-sparse representation is available
!
      CALL New(AMat,(/GDIISMemory,GDIISMemory/))
!
!     CALL DGEMM_TNc(GDIISMemory,DimGDIIS,GDIISMemory,One,Zero,  &
!                  SRDispl%D,SRDispl%D,AMat%D)
      CALL DGEMM_TNc(GDIISMemory,DimGDIIS,GDIISMemory,One,Zero,  &
                   RefGrad%D,RefGrad%D,AMat%D)
!
! Now, calculate eigenvalues and eigenvectors of Overlap
!
      CALL New(EigVal,GDIISMemory)
      EigVal%D=Zero
      CALL New(TrfDispl,(/GDIISMemory,GDIISMemory/))
      CALL New(TrfGrad,(/GDIISMemory,GDIISMemory/))
!
        CALL SetDSYEVWork(GDIISMemory)
!
        BLKVECT%D=AMat%D
        CALL DSYEV('V','U',GDIISMemory,BLKVECT%D,BIGBLOK,BLKVALS%D, &
          BLKWORK%D,BLKLWORK,INFO)
        IF(INFO/=SUCCEED) &
        CALL Halt('DSYEV hosed in RotationsOff. INFO='&
                   //TRIM(IntToChar(INFO)))
!
! Construct transformation matrices
!
        DO I=1,GDIISMemory
            Sum=Zero
          DO J=1,GDIISMemory
            Sum=Sum+BLKVECT%D(J,I)
          ENDDO
          IF(ABS(Sum)>1.D-8) THEN
            SumX=One/Sum
            TrfDispl%D(:,I)=SumX*BLKVECT%D(:,I)
            TrfGrad%D(:,I)=Sum*BLKVECT%D(:,I)
          ELSE
            TrfDispl%D(:,I)=Zero 
            TrfGrad%D(:,I)=Zero
            Selection%I(I)=0
          ENDIF
        ENDDO
!
! Now, transform SR and Ref structures and RefGrad
! to get the unselected basis for new GDIIS steps
!
         CALL New(AuxStruct,(/DimGDIIS,GDIISMemory/))
!!
!          CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero,  &
!            SRStruct%D,TrfDispl%D,AuxStruct%D)
!          SRStruct%D=AuxStruct%D
!!
!          CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero,  &
!            RefStruct%D,TrfDispl%D,AuxStruct%D)
!          RefStruct%D=AuxStruct%D
!!
!! Gradient transformation
!!
!          CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero,  &
!            RefGrad%D,TrfGrad%D,AuxStruct%D)
!          RefGrad%D=AuxStruct%D
!!
!          SRDispl%D=SRStruct%D-RefStruct%D
!!
!! Calculate scalar product of gradients and displacement in the 
!! orthogonal-displacements' basis
!
          MinDXDotDX  =1.D99
          MinGrdDotDX =1.D99
          MinGrdDotGrd=1.D99
        DO I=1,GDIISMemory
          IF(Selection%I(I)/=1) CYCLE
          SumX=DOT_PRODUCT(RefGrad%D(:,I),SRDispl%D(:,I))
          SumY=DOT_PRODUCT(SRDispl%D(:,I),SRDispl%D(:,I))
          SumZ=DOT_PRODUCT(RefGrad%D(:,I),RefGrad%D(:,I))
          GrdDotDX%D(I) =SumX
          DXDotDX%D(I)  =SumY
          GrdDotGrd%D(I)=SumZ
          MinGrdDotDX =MIN(MinGrdDotDX,ABS(SumX))
          MinDXDotDX  =MIN(MinDXDotDX      ,SumY)
          MinGrdDotGrd=MIN(MinGrdDotGrd    ,SumZ)
        ENDDO  
!
! Now, redefine metric of coordinates, and gradients
! Currently only metric of coordinates determines all other metrics.
!
        IF(GeOpCtrl%GDIISMetricOn) THEN
          MFact=One/SQRT(MinDXDotDX)
!         MFact=One/SQRT(MinGrdDotGrd)
          GeOpCtrl%GDIISMetric=MFact*GeOpCtrl%GDIISMetric
          RefStruct%D         =MFact*RefStruct%D 
          SRStruct%D          =MFact*SRStruct%D 
          SRDispl%D           =MFact*SRDispl%D 
          Sum=One/MFact
          RefGrad%D           =Sum*RefGrad%D 
          MFact=MFact*MFact
          GrdDotDX%D          =GrdDotDX%D
          DXDotDX%D           =MFact*DXDotDX%D
          GrdDotGrd%D         =GrdDotGrd%D/MFact
          DO I=1,GDIISMemory
            IF(Selection%I(I)/=1) THEN
              RefStruct%D(:,I)    =Zero
              SRStruct%D(:,I)     =Zero
              SRDispl%D(:,I)      =Zero
              GrdDotDX%D(I)       =Zero
              DXDotDX%D(I)        =Zero
              GrdDotGrd%D(I)      =Zero
            ENDIF
          ENDDO
        ENDIF !!!! Metric On
!
! Check for 'uphill' trajectories and 'reverse' them.
!
!        DO I=1,GDIISMemory
!          Sum=GrdDotDX%D(I)
!          IF(Sum>Zero) THEN
!            SRStruct%D(:,I)=RefStruct%D(:,I)-SRDispl%D(:,I)
!            GrdDotDX%D(I)=-GrdDotDX%D(I)
!          ENDIF
!        ENDDO
!
! Calculate total energy lowering over images
!
          ED0=GrdDotDX%D(1)
          EDT=ED0
        DO I=2,GDIISMemory
	  Sum=GrdDotDX%D(I)
	  EDT=EDT+Sum
	  ED0=MAX(ED0,Sum)
	ENDDO
!
! Define unnormalized, unselected set of coefficients
!
        DO I=1,GDIISMemory
            Sum=Zero
          DO J=1,GDIISMemory
            Sum=Sum+BLKVECT%D(J,I)
          ENDDO
          IF(ABS(BLKVALS%D(I))>1.D-8) THEN
            Coeffs%D(I)=Sum/BLKVALS%D(I)
!           Coeffs%D(I)=One/DXDotDX%D(I)
          ELSE
            Coeffs%D(I)=Zero   
          ENDIF
        ENDDO
!
! Select out subspace for final GDIIS
!
        CALL DGEMM_NNc(GDIISMemory,GDIISMemory,1,One,Zero,  &
          BLKVECT%D,Coeffs%D,DXDotDX%D)
          Coeffs%D=DXDotDX%D
!
        CALL GDIISSelect(Coeffs%D,Selection%I)
selection%i=1
!
      CALL UnSetDSYEVWork()
!
! Rescale coeffs to get a sum of One.
!
          Sum=Zero
        DO I=1,GDIISMemory 
          IF(Selection%I(I)==1) THEN 
            Sum=Sum+Coeffs%D(I)  
          ELSE
            Coeffs%D(I)=Zero
          ENDIF
        ENDDO
        Sum=One/Sum
        Coeffs%D=Sum*Coeffs%D
!
! Calculate new geometry, linearcombine previous steps 
!
          AuxVect%D=Zero
        DO I=1,GDIISMemory
          IF(Selection%I(I)/=1) CYCLE
          AuxVect%D=AuxVect%D+Coeffs%D(I)*SRStruct%D(:,I)
        ENDDO
!
! Restore original metric on resulting GDIIS geometry
!
        IF(GeOpCtrl%GDIISMetricOn) THEN
          MFact=One/GeOpCtrl%GDIISMetric
          AuxVect%D=MFact*AuxVect%D
        ENDIF
!
! Fill new geometry into GM array
!
        CALL CartRNK1ToCartRNK2(AuxVect%D,GMLoc%Carts%D)
!
! Save unitary transformed structures and gradients
!
        CALL Put(0,'SrMemory')
        CALL Put(0,'RefMemory')
        CALL Put(0,'CartGradMemory')
!       CALL Put(0,'IntGradMemory')
        DO I=1,GDIISMemory
          IF(Selection%I(I)==1) THEN
            AuxVect%D=SRStruct%D(:,I)
            CALL PutSRStep(Vect_O=AuxVect,Tag_O='SR',Metric_O=.FALSE.)
            AuxVect%D=RefStruct%D(:,I)
            CALL PutSRStep(Vect_O=AuxVect,Tag_O='Ref',Metric_O=.FALSE.)
            AuxVect%D=RefGrad%D(:,I)
            CALL PutSRStep(Vect_O=AuxVect,Tag_O='CartGrad',&
                           Metric_O=.FALSE.)
          ENDIF
        ENDDO
!
!call prtxyz(GMLoc%Carts%D,Title_O='Final Carts.')
!
! Tidy up
!
       CALL Delete(Selection)
       CALL Delete(AuxStruct)
       CALL Delete(GrdDotGrd)
       CALL Delete(GrdDotDX)
       CALL Delete(DXDotDX)
       CALL Delete(AuxVect)
       CALL Delete(Coeffs)
       CALL Delete(EigVal)
       CALL Delete(TrfDispl)
       CALL Delete(TrfGrad)
       CALL Delete(MetricA)
       CALL Delete(AMat)
       CALL Delete(RefGrad)
       CALL Delete(RefStruct)
       CALL Delete(SRStruct)
       CALL Delete(SRDispl)
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
!-------------------------------------------------------
!

!
     SUBROUTINE LJCell(GMLoc,LJCutOff,AtmMark,LJEps,LJRad, &
       XYZLJCell,AtmMarkLJCell,LJEpsLJCell,LJRadLJCell,NAtomsLJCell)
!
! WARNING! Pass in coordinates in Angstroms!
!
! Calculate coordinates of multiples of the elementary cell
! including the central cell, and store them in XYZLJCell.
! Storage is saved by filtering out replica atoms, which are
! far from the LJ sphere of the central cell. It's worth
! for the if statements, since later we can save at least 
! the same amount of if-s.
!
! WARNING! Pass in coordinates in wrapped representation,
! so that all coordinates be greater than zero for central cell !!!
!
! WARNING! Current implementation is optimal for cells with rectangular
! lattice vectors! This is the case for large MD simulations.
! For small elementary cells the LJ sum should be fast enough,
! even with this technique.
!
     TYPE(CRDS) :: GMLoc
     REAL(DOUBLE) :: LJCutOff,MaxCart,DX,DY,DZ
     REAL(DOUBLE) :: XTrans,YTrans,ZTrans
     INTEGER :: I,J,K,L,NX,NZ,NY,NAtomsLJCell
     INTEGER :: II,IA,IB,IC
     TYPE(DBL_VECT) :: LJEps,LJRad
     TYPE(INT_VECT) :: AtmMark
     TYPE(DBL_RNK2) :: XYZLJCell
     TYPE(DBL_RNK2) :: XYZLJCell2
     TYPE(DBL_VECT) :: LJEpsLJCell,LJRadLJCell
     TYPE(INT_VECT) :: AtmMarkLJCell
     TYPE(DBL_VECT) :: LJEpsLJCell2,LJRadLJCell2
     TYPE(INT_VECT) :: AtmMarkLJCell2
!
! First, calculate maximum Cartesian extension of the elementary cell
!
     DX=MAX(GMLoc%PBC%BoxShape(1,1),GMLoc%PBC%BoxShape(2,1),GMLoc%PBC%BoxShape(3,1))
     DY=MAX(GMLoc%PBC%BoxShape(1,2),GMLoc%PBC%BoxShape(2,2),GMLoc%PBC%BoxShape(3,2))
     DZ=MAX(GMLoc%PBC%BoxShape(1,3),GMLoc%PBC%BoxShape(2,3),GMLoc%PBC%BoxShape(3,3))
     MaxCart=MAX(DX,DY,DZ) 
!
     IF(LJCutOff<DX) THEN
       NX=1
     ELSE
       NX=INT(LJCutOff/DX)+1
     ENDIF
!
     IF(LJCutOff<DY) THEN
       NY=1
     ELSE
       NY=INT(LJCutOff/DY)+1
     ENDIF
!
     IF(LJCutOff<DZ) THEN
       NZ=1
     ELSE
       NZ=INT(LJCutOff/DZ)+1
     ENDIF
!
! Calc. Number of atoms in the additional system
! and allocate XYZLJCell
!
     NAtomsLJCell=(2*NX+1)*(2*NY+1)*(2*NZ+1)*GMLoc%Natms
     CALL New(XYZLJCell2,(/3,NAtomsLJCell/))
     CALL New(AtmMarkLJCell2,NAtomsLJCell)
     CALL New(LJEpsLJCell2,NAtomsLJCell)
     CALL New(LJRadLJCell2,NAtomsLJCell)
!
! First, copy central cell coordinates into new array
!
      DO I=1,GMLoc%Natms
        XYZLJCell2%D(1:3,I)=GMLoc%Carts%D(1:3,I)
        AtmMarkLJCell2%I(I)=AtmMark%I(I)
        LJEpsLJCell2%D(I)=LJEps%D(I)
        LJRadLJCell2%D(I)=LJRad%D(I)
      ENDDO
!
! Now, copy atoms from LJ region of central cell
!
      II=GMLoc%Natms
      DO IA=-NX,NX
      DO IB=-NY,NY
      DO IC=-NZ,NZ
  IF(IA==0.AND.IB==0.AND.IC==0) CYCLE
       DO I=1,GMLoc%Natms
         XTrans=GMLoc%Carts%D(1,I)+IA*GMLoc%PBC%BoxShape(1,1)+&
                IB*GMLoc%PBC%BoxShape(2,1)+IC*GMLoc%PBC%BoxShape(3,1)
                IF(XTrans<-LJCutOff .OR. XTrans>DX+LJCutOff) CYCLE
         YTrans=GMLoc%Carts%D(2,I)+IA*GMLoc%PBC%BoxShape(1,2)+&
                IB*GMLoc%PBC%BoxShape(2,2)+IC*GMLoc%PBC%BoxShape(3,2)
                IF(YTrans<-LJCutOff .OR. YTrans>DY+LJCutOff) CYCLE
         ZTrans=GMLoc%Carts%D(3,I)+IA*GMLoc%PBC%BoxShape(1,3)+&
                IB*GMLoc%PBC%BoxShape(2,3)+IC*GMLoc%PBC%BoxShape(3,3)
                IF(ZTrans<-LJCutOff .OR. ZTrans>DZ+LJCutOff) CYCLE
         II=II+1
         XYZLJCell2%D(1,II)=XTrans
         XYZLJCell2%D(2,II)=YTrans
         XYZLJCell2%D(3,II)=ZTrans
         AtmMarkLJCell2%I(II)=AtmMark%I(I)
         LJEpsLJCell2%D(II)=LJEps%D(I)
         LJRadLJCell2%D(II)=LJRad%D(I)
       ENDDO
     ENDDO
     ENDDO
     ENDDO
!
! Compress
!
     NAtomsLJCell=II 
     CALL New(XYZLJCell,(/3,NAtomsLJCell/))
     CALL New(AtmMarkLJCell,NAtomsLJCell)
     CALL New(LJEpsLJCell,NAtomsLJCell)
     CALL New(LJRadLJCell,NAtomsLJCell)
     XYZLJCell%D(1:3,1:NAtomsLJCell)=XYZLJCell2%D(1:3,1:NAtomsLJCell)
     AtmMarkLJCell%I(1:NAtomsLJCell)=AtmMarkLJCell2%I(1:NAtomsLJCell)
     LJEpsLJCell%D(1:NAtomsLJCell)=LJEpsLJCell2%D(1:NAtomsLJCell)
     LJRadLJCell%D(1:NAtomsLJCell)=LJRadLJCell2%D(1:NAtomsLJCell)
!
! Tidy up
!
     CALL Delete(XYZLJCell2)
     CALL Delete(AtmMarkLJCell2)
     CALL Delete(LJEpsLJCell2)
     CALL Delete(LJRadLJCell2)
!
     END SUBROUTINE LJCell
!

!
!-------------------------------------------------------------------
!
     SUBROUTINE SetOneLJCell(GMLoc,AtmMark,LJEps,LJRad, &
         XYZLJCell,AtmMarkLJCell,LJEpsLJCell,LJRadLJCell,NAtomsLJCell)
!
     TYPE(CRDS)     :: GMLoc
     INTEGER        :: NAtomsLJCell
     TYPE(DBL_VECT) :: LJEps,LJRad,LJEpsLJCell,LJRadLJCell
     TYPE(INT_VECT) :: AtmMark,AtmMarkLJCell
     TYPE(DBL_RNK2) :: XYZLJCell
     REAL(DOUBLE)   :: LJCutOff
!
     NAtomsLJCell=GMLoc%Natms
     CALL New(XYZLJCell,(/3,NAtomsLJCell/))
     CALL New(AtmMarkLJCell,NAtomsLJCell)
     CALL New(LJEpsLJCell,NAtomsLJCell)
     CALL New(LJRadLJCell,NAtomsLJCell)
     XYZLJCell%D=GMLoc%Carts%D
     AtmMarkLJCell%I=AtmMark%I
     LJEpsLJCell%D=LJEps%D
     LJRadLJCell%D=LJRad%D
!
     END SUBROUTINE SetOneLJCell
!
!---------------------------------------------------------------------
    SUBROUTINE GetFullBMatr(NatmsLoc,CartsLoc,IntCs,NIntC,FullB,FullBt)
!
! Generate vibrational B matrix in quasi-sparse TYPE(BMATR) representation
! and transform into sparse one
!
    TYPE(CRDS) :: GMLoc
    TYPE(INTC) :: IntCs
    INTEGER    :: NIntC,NCart,I,J,K,NatmsLoc
    TYPE(BMATR):: B
    TYPE(DBL_RNK2) :: FullB,FullBt
    REAL(DOUBLE),DIMENSION(:,:) :: CartsLoc
    REAL(DOUBLE) :: SUM
!
! Calculate B matrix in Atomic Units
!
    CALL BMatrix(NatmsLoc,CartsLoc,NIntC,IntCs,B)
!
    NCart=3*NatmsLoc
!
! Now, turn B matrix into full representation
!
        FullB%D=Zero
        FullBt%D=Zero
    DO I=1,NIntC
      DO J=1,12
        K=B%IB(I,J) 
        IF(K==0) CYCLE
        SUM=B%B(I,J)
        FullB%D(I,K)=SUM
        FullBt%D(K,I)=SUM
      ENDDO
!write(*,200) i,b%ib(i,1:12)
!write(*,201) i,b%b(i,1:12)
!write(*,100) i,(fullb%d(i,k),k=1,ncart)
    ENDDO
!100 format(I4,100F7.3)
!200 format(I4,100I4)
!201 format(I4,100F7.3)
!stop
!
! Tidy up
!
    CALL Delete(B)
!
    END SUBROUTINE GetFullBMatr
!
!-------------------------------------------------------
!
    SUBROUTINE GetFullGcInv(FullB,FullBt,FullGcInv,NIntC,NCart)
!
    TYPE(DBL_RNK2) :: FullB,FullBt,FullGcInv
    TYPE(DBL_RNK2) :: FullGc,TestMat
    INTEGER :: NIntC,NCart,I,J,K
    REAL(DOUBLE) :: SUM
!
    CALL New(FullGc,(/NCart,NCart/))
    CALL New(TestMat,(/NCart,NCart/))
!
    FullGc%D=Zero
    CALL DGEMM_NNc(NCart,NIntC,NCart,One,Zero,FullBt%D,FullB%D,FullGc%D)
!
! Get SVD inverse
!
    CALL SetDSYEVWork(NCart)
    CALL FunkOnSqMat(NCart,Inverse,FullGc%D,FullGcInv%D)
    CALL UnSetDSYEVWork()
!!
!! Test inverse
!!
!    CALL DGEMM_NNc(NCart,NCart,NCart,One,Zero,FullGcInv%D,FullGc%D,TestMat%D)
!    DO I=1,NCart
!      write(*,100) (TestMat%D(I,J),J=1,NCart)
!    ENDDO
!100 FORMAT(10F8.4)
!call pprint(FullGcInv,' FullGcInv ',unit_o=6)
!
    CALL Delete(TestMat)
    CALL Delete(FullGc)
!
    END SUBROUTINE GetFullGcInv
!-------------------------------------------------------
!
    SUBROUTINE CartToInternal(XYZ,IntCs,VectCart,VectInt)
!
      REAL(DOUBLE),DIMENSION(:,:)  :: XYZ   
      REAL(DOUBLE),DIMENSION(:)    :: VectCart,VectInt
      TYPE(DBL_VECT) :: VectCartAux,VectIntAux
      TYPE(DBL_VECT) :: VectCartAux2,VectIntAux2
      TYPE(DBL_RNK2) :: FullB,FullBt,FullGcInv
      REAL(DOUBLE)   :: DiffMax,RMSD
      REAL(DOUBLE)   :: GrdTrfCrit,MaxGradDiff,Sum
      INTEGER        :: NCart,I,II,J,NIntC
      INTEGER        :: MaxIt_GrdTrf,NatmsLoc
      TYPE(INTC)     :: IntCs
!
      CALL OpenASCII(OutFile,Out)
!
! Iteration control parameters
!
      GrdTrfCrit  =GeOpCtrl%GrdTrfCrit !!! For Gradient transformation
      MaxIt_GrdTrf=GeOpCtrl%MaxIt_GrdTrf
      MaxGradDiff =GeOpCtrl%MaxGradDiff
      NatmsLoc=SIZE(XYZ,2)
      NIntC=SIZE(IntCs%Def)
      NCart=3*NatmsLoc        
!
      CALL New(FullB,(/NIntC,NCart/))
      CALL New(FullBt,(/NCart,NIntC/))
      CALL New(FullGcInv,(/NCart,NCart/))
!
      CALL New(VectCartAux,NCart)
      CALL New(VectCartAux2,NCart)
      CALL New(VectIntAux,NIntC)
      CALL New(VectIntAux2,NIntC)
!
      VectInt=Zero
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,*) 'Gradient transformation, No. Int. Coords= ',NIntC
        WRITE(Out,*) 'Gradient transformation, No. Int. Coords= ',NIntC
      ENDIF
!
! Cartesian --> Internal transformation
!
! Get B matrix and Bt*B inverse
!
      CALL GetFullBMatr(NatmsLoc,XYZ,IntCs,NIntC,FullB,FullBt)
      CALL GetFullGcInv(FullB,FullBt,FullGcInv,NIntC,NCart)
!
      II=0
100   CONTINUE      
!
      VectCartAux%D=Zero
!
! gc-Bt*gi
!
      CALL DGEMM_NNc(NCart,NIntC,1,One,Zero,FullBt%D,VectInt,VectCartAux%D)
      VectCartAux%D=VectCart-VectCartAux%D
!
! GcInv*[gc-Bt*gi]
!
      CALL DGEMM_NNc(NCart,NCart,1,One,Zero,&
           FullGcInv%D,VectCartAux%D,VectCartAux2%D)
!
! B*GcInv*[gc-Bt*gi]
!
      CALL DGEMM_NNc(NIntC,NCart,1,One,Zero, &
           FullB%D,VectCartAux2%D,VectIntAux%D)
!
! Check convergence
!
      II=II+1
        DiffMax=Zero
        DO I=1,NIntC ; DiffMax=MAX(DiffMax,ABS(VectIntAux%D(I))) ; ENDDO
        RMSD=DOT_PRODUCT(VectIntAux%D,VectIntAux%D)
        RMSD=SQRT(RMSD/DBLE(NIntC))
!
! IF DiffMax is too large, eg. due to the 'bad' quality 
! of the preconditioner, rescale gradients
!
      IF(DiffMax>MaxGradDiff) THEN
        WRITE(*,*) 'Rescale Step from ',DiffMax,' to ',MaxGradDiff
        WRITE(Out,*) 'Rescale Step from ',DiffMax,' to ',MaxGradDiff
        SUM=MaxGradDiff/DiffMax
        VectIntAux%D(:)=SUM*VectIntAux%D(:)
        DiffMax=MaxGradDiff
      ENDIF
!
! gi+B*GcInv*[gc-Bt*gi]
!
      VectInt=VectInt+VectIntAux%D
!
! Review iteration
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,110) II,DiffMax,RMSD
        WRITE(Out,110) II,DiffMax,RMSD
      ENDIF
110   FORMAT('Grad Trf, step= ',I3,' MaxChange= ',F12.6,' ChangeNorm= ',F12.6)
!      
      IF(DiffMax>GrdTrfCrit.AND.II<=MaxIt_GrdTrf) THEN
        GO TO 100      
      ELSE
        IF(II>MaxIt_GrdTrf) THEN
          IF(PrintFlags%Geop==DEBUG_GEOP) THEN
            WRITE(*,*) 'Stop Gradient Trf, max. '//&
                       'number of Iterations exceeded!'
            WRITE(*,*) 'Use current gradient vector!'
            WRITE(Out,*) 'Stop Gradient Trf, max. number '//&
                         'of Iterations exceeded!'
            WRITE(Out,*) 'Use current gradient vector!'
          ENDIF
        ENDIF
          CALL Put(FullB,'FullB')
          CALL Put(FullBt,'FullBt')
          CALL Put(FullGcInv,'FullGcInv')
      ENDIF
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,120) II
        WRITE(Out,120) II
      ENDIF
120   FORMAT('Gradient transformation converged in ',I3,' steps')
!
      CLOSE(Out,STATUS='KEEP')
!
!
! Tidy up
!
      CALL Delete(VectIntAux2)
      CALL Delete(VectIntAux)
      CALL Delete(VectCartAux2)
      CALL Delete(VectCartAux)
      CALL Delete(FullGcInv)
      CALL Delete(FullBt)
      CALL Delete(FullB)
!
    END SUBROUTINE CartToInternal
!
!---------------------------------------------------------------------
!
    SUBROUTINE InternalToCart(XYZ,IntCs,VectInt)
!
    REAL(DOUBLE),DIMENSION(:,:) :: XYZ
    REAL(DOUBLE),DIMENSION(:) :: VectInt
    TYPE(DBL_VECT)            :: VectCart
    TYPE(DBL_VECT)            :: VectCartAux,VectIntAux
    TYPE(DBL_VECT)            :: VectCartAux2,VectIntAux2
    TYPE(DBL_VECT)            :: VectIntReq
    TYPE(DBL_RNK2)            :: FullB,FullBt,FullGcInv,ActCarts
    REAL(DOUBLE)              :: DiffMax,RMSD,RMSDOld,TrixThresh
    REAL(DOUBLE)              :: CooTrfCrit,MaxCartDiff,Sum
    REAL(DOUBLE)              :: DistRefresh
    REAL(DOUBLE)              :: SumX,SumY,SumZ
    REAL(DOUBLE)              :: ConstrMax,ConstrRMS
    REAL(DOUBLE)              :: ConstrRMSOld
    REAL(DOUBLE)              :: ConstrMaxCrit,RMSCrit
    INTEGER                   :: NCart,I,IStep,J,NIntC,NConstr
    INTEGER                   :: MaxIt_CooTrf,NatmsLoc
    INTEGER                   :: NCartConstr
    TYPE(INTC)                :: IntCs
    TYPE(BMATR)               :: B
    LOGICAL                   :: RefreshB,RefreshAct
    LOGICAL                   :: DoIterate
    TYPE(INT_VECT)            :: MMAtNum
!
      NatmsLoc=SIZE(XYZ,2)
      NCart=3*NatmsLoc   
      NIntC=SIZE(IntCs%Def)
!
! iteration control parameters
!
      CooTrfCrit   = GeOpCtrl%CooTrfCrit 
      MaxIt_CooTrf = GeOpCtrl%MaxIt_CooTrf
      MaxCartDiff  = GeOpCtrl%MaxCartDiff
      DistRefresh  = GeOpCtrl%DistRefresh
      ConstrMaxCrit= GeOpCtrl%ConstrMaxCrit
      RMSCrit      = GeOpCtrl%RMSCrit
      NCartConstr  = GeOpCtrl%NCartConstr
!
! Refresh B matrix during iterative back-trf, if displacements are too large?
!
      RefreshB=.TRUE.
      RefreshAct=.TRUE.
!
! Auxiliary arrays
!
      CALL New(ActCarts,(/3,NatmsLoc/))
      CALL New(FullB,(/NIntC,NCart/))
      CALL New(FullBt,(/NCart,NIntC/))
      CALL New(FullGcInv,(/NCart,NCart/))
!
      CALL New(VectCart,NCart)
      CALL New(VectCartAux,NCart)
      CALL New(VectCartAux2,NCart)
      CALL New(VectIntAux,NIntC)
      CALL New(VectIntAux2,NIntC)
      CALL New(VectIntReq,NIntC)
!
! Calc values of internals in atomic unit, and add displacement,
! which is stored in VectInt. Then convert into Sparse matrx repr.
!
      CALL INTCValue(IntCs,XYZ)
!
! The required new value of internal coordinates
!
      VectIntReq%D=VectInt+IntCs%Value
      CALL MapBackAngle(IntCs,NIntC,VectIntReq%D) 
!
!initialization of new Cartesians
!
        ActCarts%D=XYZ            
      DO I=1,NatmsLoc
        J=3*(I-1)
        VectCart%D(J+1)=XYZ(1,I)
        VectCart%D(J+2)=XYZ(2,I)
        VectCart%D(J+3)=XYZ(3,I)
      ENDDO
!
! Internal --> Cartesian transformation
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,*) 'Iterative back-transformation, '//&
                   'No. Int. Coords= ',NIntC
        CALL OpenAscii(OutFile,Out)
        WRITE(Out,*) 'Iterative back-transformation,'//&
                     ' No. Int. Coords= ',NIntC
        CLOSE(Out,STATUS='KEEP')
      ENDIF
!
      ConstrMax=ConstrMaxCrit*10.D0
      ConstrRMS=1.D0
      ConstrRMSOld=2.D0
      RMSD=1.D+9
!
      IStep=0
200   CONTINUE
!
! Get B and GcInv
!
       IF(IStep==0) THEN
         CALL Get(FullB,'FullB')
         CALL Get(FullBt,'FullBt')
         CALL Get(FullGcInv,'FullGcInv')
       ELSE IF(RefreshB.AND.RefreshAct) THEN
         CALL GetFullBMatr(NatmsLoc,ActCarts%D,IntCs,NIntC,FullB,FullBt)
         CALL GetFullGcInv(FullB,FullBt,FullGcInv,NIntC,NCart)
       ENDIF
!
! Compute actual value of internals 
!
      IF(IStep>0) CALL INTCValue(IntCs,ActCarts%D)
!
! Calculate difference between required and actual internals
! Calc [phi_r-phi_a]
!
      VectIntAux%D=VectIntReq%D-IntCs%Value
      CALL MapDisplTors(IntCs,NIntC,VectIntAux%D) 
!
! Check convergence on constraints
!
      IF(NConstr/=0) THEN
        ConstrRMSOld=ConstrRMS
        CALL ConstrConv(IntCs,VectIntAux%D,ConstrMax,ConstrRMS)
      ENDIF
!
! Do transformation
! 
! Bt*[phi_r-phi_a]
!
      CALL DGEMM_NNc(NCart,NIntC,1,One,Zero,FullBt%D,VectIntAux%D,VectCartAux%D)
!
! GcInv*Bt*[phi_r-phi_a]
!
      CALL DGEMM_NNc(NCart,NCart,1,One,Zero,FullGcInv%D,VectCartAux%D,VectCartAux2%D)
!
! Project out translations and rotations
!
!     CALL TranslsOff(VectCartAux2%D)
!     CALL RotationsOff(VectCartAux2%D,ActCarts%D)
!
! Check convergence
!
      IStep=IStep+1
!
      RMSDOld=RMSD
      CALL CartConv(VectCartAux2%D,MaxCartDiff,DiffMax,RMSD, &
                    IntCs,NCartConstr)
!
! Refresh B matrix?  
!
      IF(DiffMax>DistRefresh) THEN
        RefreshAct=.TRUE.
      ELSE
        RefreshAct=.FALSE.
      ENDIF
!
! Modify Cartesians
!
      VectCart%D=VectCart%D+VectCartAux2%D
      CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
!
! Review iteration
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,210) IStep,DiffMax,RMSD
        CALL OpenAscii(OutFile,Out)
        WRITE(Out,210) IStep,DiffMax,RMSD
        CLOSE(Out,STATUS='KEEP')
      ENDIF
210   FORMAT('Step= ',I3,'   Max_DX= ',F12.6,'  X_RMSD= ',F12.6)
!      
      CALL BackTrfConvg(DoIterate,DiffMax,CooTrfCrit, &
        RMSD,RMSDOld,RMSCrit,ConstrMax,ConstrMaxCrit, &
        ConstrRMS,ConstrRMSOld,NConstr,MaxIt_CooTrf,IStep)
!
      IF(DoIterate) THEN
!
!       IF(RefreshB.AND.RefreshAct) THEN
!         CALL Delete(FullB)
!         CALL Delete(FullBt)
!         CALL Delete(FullGcInv)
!       ENDIF
        GO TO 200
      ELSE
          IF(IStep>MaxIt_CooTrf) THEN
            IF(RMSD>0.01D0) THEN
              CALL Halt('Iterative backtransformation did not converge')
            ENDIF
            IF(PrintFlags%Geop==DEBUG_GEOP) THEN
              CALL OpenAscii(OutFile,Out)
              WRITE(*,180) 
              WRITE(Out,180) 
              WRITE(*,190) 
              WRITE(Out,190) 
              CLOSE(Out,STATUS='KEEP')
            ENDIF
            GO TO 300
          ENDIF
      ENDIF
180   FORMAT('Stop Coord Back-Trf, max. number of Iterations exceeded!')
190   FORMAT('Use Current Geometry!')
!
! At this point, convergence has been reached
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        CALL OpenAscii(OutFile,Out)
        WRITE(*,220) IStep
        WRITE(Out,220) IStep
        CLOSE(Out,STATUS='KEEP')
      ENDIF
220   FORMAT('Coordinate back-transformation converged in ',I3,' steps')
!
300   CONTINUE
!
! Fill new Cartesians into XYZ  
!
      XYZ=ActCarts%D
!
! Final internal coordinates
!
!     CALL INTCValue(IntCs,XYZ)
!     CALL PrtIntCoords(IntCs,IntCs%Value,'Final Internals')
!
! Tidy up
!
      CALL Delete(VectIntReq)
      CALL Delete(VectIntAux2)
      CALL Delete(VectIntAux)
      CALL Delete(VectCartAux2)
      CALL Delete(VectCartAux)
      CALL Delete(VectCart)
      CALL Delete(FullGcInv)
      CALL Delete(FullB)
      CALL Delete(FullBt)
      CALL Delete(ActCarts)
!
      END SUBROUTINE InternalToCart
!
!----------------------------------------------------------
!
      SUBROUTINE PrtIntCoords(IntCs,Value,CHAR)
!
      TYPE(INTC) :: IntCs
      INTEGER    :: I,NIntC,J
      REAL(DOUBLE),DIMENSION(:) :: Value
      REAL(DOUBLE) :: SUM
      CHARACTER(LEN=*) :: CHAR   
!
      NIntC=SIZE(IntCs%Def)
!
      CALL OpenAscii(OutFile,Out)
!
      WRITE(*,*) TRIM(CHAR)
      WRITE(Out,*) TRIM(CHAR)
      WRITE(*,*) '         INTERNAL COORDINATES'
      WRITE(*,*) 'DEFINITION  ATOMS_INVOLVED     VALUE'
      WRITE(Out,*) 'INTERNAL COORDINATES'
      WRITE(Out,*) 'DEFINITION  ATOMS_INVOLVED  VALUE'
      DO I=1,NIntC
        IF(IntCs%Def(I)(1:4)=='STRE') THEN
          SUM=Value(I)/AngstromsToAu
        ELSE IF(IntCs%Def(I)(1:4)=='BEND'.OR. &
                IntCs%Def(I)(1:4)=='TORS'.OR. &
                IntCs%Def(I)(1:4)=='LINB'.OR. &
                IntCs%Def(I)(1:4)=='OUTP') THEN
          SUM=Value(I)*180.D0/PI
        ELSE IF(IntCs%Def(I)(1:4)=='CART') THEN
          SUM=Value(I)
        ENDIF
          WRITE(*,111) I,IntCs%Def(I),IntCs%Atoms(I,1:4),SUM,IntCs%Constraint(I),IntCs%Active(I)
          WRITE(Out,111) I,IntCs%Def(I),IntCs%Atoms(I,1:4),SUM,IntCs%Constraint(I),IntCs%Active(I)
      ENDDO
!      
 111 FORMAT(I4,2X,A5,2X,4I3,2X,F12.6,L5,L5)
 222 FORMAT(I4,2X,A5,2X,4I3,2X,3F12.6,L5,L5)
!
      CLOSE(Out,STATUS='KEEP')
!
      END SUBROUTINE PrtIntCoords
!
!----------------------------------------------------------
!
      SUBROUTINE RedundancyOffFull(IntDispl,NCart,FullB_O,FullBt_O,FullGcInv_O)
!
      REAL(DOUBLE),DIMENSION(:)  :: IntDispl
      TYPE(DBL_RNK2),OPTIONAL :: FullB_O,FullBt_O,FullGcInv_O
      TYPE(DBL_RNK2)          :: FullB,FullBt,FullGcInv
      TYPE(DBL_VECT)          :: IntAux
      TYPE(DBL_VECT)          :: CartAux1,CartAux2
      INTEGER                 :: NIntC,NCart
      REAL(DOUBLE)            :: SUM1,SUM2
!
! Project out first order redundancies
!
      NIntC=SIZE(IntDispl)
      CALL New(IntAux,NIntC)
      CALL New(CartAux1,NCart)
      CALL New(CartAux2,NCart)
        CALL New(FullB,(/NIntC,NCart/))
        CALL New(FullBt,(/NCart,NIntC/))
        CALL New(FullGcInv,(/NCart,NCart/))
!
      IntAux%D=IntDispl
!
      IF(.NOT.PRESENT(FullB_O).OR..NOT.PRESENT(FullBt_O).OR. &
         .NOT.PRESENT(FullGcInv_O)) THEN
        CALL Get(FullB,'FullB')
        CALL Get(FullBt,'FullBt')
        CALL Get(FullGcInv,'FullGcInv')
      ELSE
        FullB%D=FullB_O%D
        FullBt%D=FullBt_O%D
        FullGcInv%D=FullGcInv_O%D
      ENDIF
!
! Carry out projection intdispl'=B*GcInv*Bt*intdispl
!
! Bt*intdispl
!
      CALL DGEMM_NNc(NCart,NIntC,1,One,Zero,FullBt%D,IntDispl,CartAux1%D)
!
! GcInv*Bt*intdispl
!
      CALL DGEMM_NNc(NCart,NCart,1,One,Zero,FullGcInv%D,CartAux1%D,CartAux2%D)
!
! intdispl'=B*GcInv*Bt*intdispl
!
      CALL DGEMM_NNc(NIntC,NCart,1,One,Zero,FullB%D,CartAux2%D,IntDispl)
!
! Calculate rate of redundancy projected out
!
      SUM1=SQRT(DOT_PRODUCT(IntAux%D,IntDispl))
      SUM2=SQRT(DOT_PRODUCT(IntAux%D,IntAux%D))
      SUM1=(SUM2-SUM1)/SUM2*100.D0 !!! percentage, projected out
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,100) SUM1
100     FORMAT('PERCENT OF REDUNDANCY PROJECTED OUT= ',F7.3,'%')
      ENDIF
!
      CALL Delete(FullB)
      CALL Delete(FullBt)
      CALL Delete(FullGcInv)
      CALL Delete(IntAux)
      CALL Delete(CartAux1)
      CALL Delete(CartAux2)
!
      END SUBROUTINE RedundancyOffFull
!
!----------------------------------------------------------
!
      SUBROUTINE RotationsOff(CartVect,XYZ)
      REAL(DOUBLE),DIMENSION(:) :: CartVect
      REAL(DOUBLE),DIMENSION(:,:) :: XYZ
      REAL(DOUBLE)                :: X,Y,Z,XX,YY,ZZ,XY,YZ,ZX
      REAL(DOUBLE)                :: CMX,CMY,CMZ
      REAL(DOUBLE)                :: SUM,SUM1,SUM2,SUM3
      REAL(DOUBLE)                :: V1X,V1Y,V1Z
      REAL(DOUBLE)                :: V2X,V2Y,V2Z
      REAL(DOUBLE)                :: V3X,V3Y,V3Z
      TYPE(DBL_RNK2)              :: Theta,Theta2,CMCarts
      TYPE(DBL_Vect)              :: Rot1,Rot2,Rot3,Vect
      INTEGER :: NCart,NatmsLoc,I,J,INFO
!
      NCart=SIZE(CartVect)
      NatmsLoc=NCart/3
      CALL New(Theta,(/3,3/))
      CALL New(Theta2,(/3,3/))
      CALL New(CMCarts,(/3,NatmsLoc/))
      CALL New(Rot1,NCart)
      CALL New(Rot2,NCart)
      CALL New(Rot3,NCart)
      CALL New(Vect,3)
!
! Calculate center of Mass (Masses are unit here)
!
      CALL CenterOfMass(CMX,CMY,CMZ,XYZ_O=XYZ)
!
! Mass centered Cartesians
!
      DO I=1,NatmsLoc
        CMCarts%D(1,I)=XYZ(1,I)-CMX
        CMCarts%D(2,I)=XYZ(2,I)-CMY
        CMCarts%D(3,I)=XYZ(3,I)-CMZ
      ENDDO
!
! Build inertial momentum tensor
!
      CALL PrincMomTens(CMCarts%D,Theta%D)
!
! Get eigenvectors of inertial momentum tensor
!
      CALL SetDSYEVWork(3)
        BLKVECT%D(1:3,1:3)=Theta%D(1:3,1:3)
        CALL DSYEV('V','U',3,BLKVECT%D,BIGBLOK,BLKVALS%D, &
        BLKWORK%D,BLKLWORK,INFO)
        IF(INFO/=SUCCEED) &
	CALL Halt('DSYEV hosed in RotationsOff. INFO='&
                   //TRIM(IntToChar(INFO)))
        Theta2%D(1:3,1:3)=BLKVECT%D(1:3,1:3)
      CALL UnSetDSYEVWork()
!
! Calculate displacements resulting from rotations
! around principal axis by a unit rotation vector pointing along
! principal axis. 
!
      DO I=1,NatmsLoc
        J=(I-1)*3
! r1
        CALL CROSS_PRODUCT(Theta2%D(:,1),CMCarts%D(:,I),Rot1%D(J+1:J+3))
! r2
        CALL CROSS_PRODUCT(Theta2%D(:,2),CMCarts%D(:,I),Rot2%D(J+1:J+3))
! r3
        CALL CROSS_PRODUCT(Theta2%D(:,3),CMCarts%D(:,I),Rot3%D(J+1:J+3))
!
      ENDDO
!
! Normalize Rot vectors!
!
      SUM=SQRT(DOT_PRODUCT(Rot1%D,Rot1%D))
      IF(Sum>1.D-6) THEN
        Sum=One/Sum
        Rot1%D=Sum*Rot1%D
      ENDIF
      SUM=SQRT(DOT_PRODUCT(Rot2%D,Rot2%D))
      IF(Sum>1.D-6) THEN
        Sum=One/Sum
        Rot2%D=Sum*Rot2%D
      ENDIF
      SUM=SQRT(DOT_PRODUCT(Rot3%D,Rot3%D))
      IF(Sum>1.D-6) THEN
        Sum=One/Sum
        Rot3%D=Sum*Rot3%D
      ENDIF
!!
!! test orthogonalities
!!
!      CALL TranslsOff(Rot1%D)
!      CALL TranslsOff(Rot2%D)
!      CALL TranslsOff(Rot3%D)
!      SUM=DOT_PRODUCT(Rot1%D,Rot1%D)
!      write(*,*) 'rot test 1 1 ',sum
!      SUM=DOT_PRODUCT(Rot1%D,Rot2%D)
!      write(*,*) 'rot test 1 2 ',sum
!      SUM=DOT_PRODUCT(Rot1%D,Rot3%D)
!      write(*,*) 'rot test 1 3 ',sum
!      SUM=DOT_PRODUCT(Rot2%D,Rot2%D)
!      write(*,*) 'rot test 2 2 ',sum
!      SUM=DOT_PRODUCT(Rot2%D,Rot3%D)
!      write(*,*) 'rot test 2 3 ',sum
!      SUM=DOT_PRODUCT(Rot3%D,Rot3%D)
!      write(*,*) 'rot test 3 3 ',sum
!
! Now, project out rotations from Cartesian displacement vector
!
      SUM=DOT_PRODUCT(CartVect(1:NCart),CartVect(1:NCart))
      SUM1=DOT_PRODUCT(Rot1%D(1:NCart),CartVect(1:NCart))
      SUM2=DOT_PRODUCT(Rot2%D(1:NCart),CartVect(1:NCart))
      SUM3=DOT_PRODUCT(Rot3%D(1:NCart),CartVect(1:NCart))
      CartVect(1:NCart)=CartVect(1:NCart)-SUM1*Rot1%D(1:NCart)
      CartVect(1:NCart)=CartVect(1:NCart)-SUM2*Rot2%D(1:NCart)
      CartVect(1:NCart)=CartVect(1:NCart)-SUM3*Rot3%D(1:NCart)
!
! Percentage of rotations
!
      SUM1=SUM1*SUM1/SUM*100.D0
      SUM2=SUM2*SUM2/SUM*100.D0
      SUM3=SUM3*SUM3/SUM*100.D0
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,100) SUM1,SUM2,SUM3
      ENDIF
100   FORMAT('Rot1= ',F7.3,'%    Rot2= ',F7.3,'%     Rot3= ',F7.3,'% ')
!
      CALL Delete(Vect)
      CALL Delete(Rot3)
      CALL Delete(Rot2)
      CALL Delete(Rot1)
      CALL Delete(CMCarts)
      CALL Delete(Theta2)
      CALL Delete(Theta)
!
      END SUBROUTINE RotationsOff
!
!----------------------------------------------------------
!
      SUBROUTINE TranslsOff(CartVect)
!
      REAL(DOUBLE),DIMENSION(:) :: CartVect
      REAL(DOUBLE)              :: SUM,SUM1,SUM2,SUM3
      TYPE(DBL_VECT)            :: Tr1,Tr2,Tr3
      INTEGER                   :: I,J,NCart,NatmsLoc
!
      NCart=SIZE(CartVect)
      NatmsLoc=NCart/3
      CALL New(Tr1,NCart) 
      CALL New(Tr2,NCart) 
      CALL New(Tr3,NCart) 
!
        Tr1%D=Zero
        Tr2%D=Zero
        Tr3%D=Zero
	Sum=One/DBLE(NatmsLoc)
      DO I=1,NatmsLoc
        J=(I-1)*3
        Tr1%D(J+1)=Sum
        Tr2%D(J+2)=Sum
        Tr3%D(J+3)=Sum
      ENDDO
!
! Now, project out translations from CartVect
!
      SUM =DOT_PRODUCT(CartVect,CartVect)
      IF(SUM > 1.D-6) THEN
        SUM1=DOT_PRODUCT(Tr1%D,CartVect)
          CartVect=CartVect-SUM1*Tr1%D
        SUM2=DOT_PRODUCT(Tr2%D,CartVect)
          CartVect=CartVect-SUM2*Tr2%D
        SUM3=DOT_PRODUCT(Tr3%D,CartVect)
          CartVect=CartVect-SUM3*Tr3%D
        SUM1=SUM1*SUM1/SUM*100.D0
        SUM2=SUM2*SUM2/SUM*100.D0
        SUM3=SUM3*SUM3/SUM*100.D0
      ELSE
        SUM1=Zero
        SUM2=Zero
        SUM3=Zero
      ENDIF
!
! Percentage of translations
!
      IF(PrintFlags%Geop==DEBUG_GEOP) THEN
        WRITE(*,100) SUM1,SUM2,SUM3
      ENDIF
100   FORMAT(' Tr1= ',F7.3,'%     Tr2= ',F7.3,'%      Tr3= ',F7.3,'% ')
!
      CALL Delete(Tr3)
      CALL Delete(Tr2)
      CALL Delete(Tr1)
!
      END SUBROUTINE TranslsOff
!-------------------------------------------------------
!
      SUBROUTINE Rotate(V1,V2,Rot)
!
! Compute matrix Rot, which rotates vector V2 into vector V1
!
      REAL(DOUBLE),DIMENSION(1:3)      :: V1,V2
      REAL(DOUBLE),DIMENSION(1:3,1:3)  :: Rot  
      REAL(DOUBLE)                 :: Sum,SumM,Sum1,Sum2,Sum3,V1N,V2N
      REAL(DOUBLE) :: PnX,PnY,PnZ
      REAL(DOUBLE) :: CosPhi,SinPhi
      INTEGER      :: I,J,Step
      TYPE(DBL_VECT) :: Vect,CrossProd
!          
      CALL New(Vect,3)
      CALL New(CrossProd,3)
!
      V1N=DSQRT(DOT_PRODUCT(V1,V1))
      V2N=DSQRT(DOT_PRODUCT(V2,V2))
      V1=V1/V1N
      V2=V2/V2N
      Step=1
!
      CrossProd%D(1)=V1(2)*V2(3)-V1(3)*V2(2)
      CrossProd%D(2)=V1(3)*V2(1)-V1(1)*V2(3) 
      CrossProd%D(3)=V1(1)*V2(2)-V1(2)*V2(1) 
!
      Sum1=DOT_PRODUCT(CrossProd%D(1:3),CrossProd%D(1:3))
!
      IF(DABS(Sum1) < 1.D-6) THEN  !!! V1 & V2 are parallel
        Rot=Zero
        Sum=DOT_PRODUCT(V1,V2)
        DO I=1,3 ; Rot(I,I)=Sum ; ENDDO
      ELSE      
        Sum=One/SQRT(Sum1)
        CrossProd%D=SUM*CrossProd%D
        CosPhi=DOT_PRODUCT(V1,V2)     
        SinPhi=SQRT(SUM1)
        Sum=One-CosPhi
!
77      CONTINUE
!
        Rot(1,1)=CosPhi + Sum*CrossProd%D(1)**2
        Rot(2,2)=CosPhi + Sum*CrossProd%D(2)**2
        Rot(3,3)=CosPhi + Sum*CrossProd%D(3)**2
!
        Sum1=Sum*CrossProd%D(1)*CrossProd%D(2)
        Sum2=SinPhi*CrossProd%D(3)
        Rot(1,2)=Sum1-Sum2
        Rot(2,1)=Sum1+Sum2
!
        Sum1=Sum*CrossProd%D(1)*CrossProd%D(3)
        Sum2=SinPhi*CrossProd%D(2)
        Rot(1,3)=Sum1+Sum2
        Rot(3,1)=Sum1-Sum2
!
        Sum1=Sum*CrossProd%D(2)*CrossProd%D(3)
        Sum2=SinPhi*CrossProd%D(1)
        Rot(2,3)=Sum1-Sum2
        Rot(3,2)=Sum1+Sum2
!
! Test rotation
!
        CALL DGEMM_NNc(3,3,1,One,Zero,Rot,V2,Vect%D(1:3))
!
        SumM=DOT_PRODUCT((V1-Vect%D),(V1-Vect%D))
        IF(SumM>1.D-6) GO TO 889
!       IF(SumM<1.D-6) write(*,*) 'rotation succesful'
        GO TO 100
889     SinPhi=-SinPhi
        Step=Step+1
        IF(Step > 2) THEN
          Call Halt('Rotation unsuccesful for torsvalue')
        ENDIF
        GOTO 77        
      ENDIF
!
100   CONTINUE
      CALL Delete(CrossProd)
      CALL Delete(Vect)
!
      END SUBROUTINE Rotate
!-------------------------------------------------------
      SUBROUTINE MapDisplTors(IntCs,NIntC,VectInt) 
!
! This routine maps displacements of torsions:
! if  0 < DisplTors < PI      => DisplTors=DisplTors
! if  0 > DisplTors <-PI      => DisplTors=DisplTors
! if PI < DisplTors < 2*PI    => DisplTors=DisplTors-2*PI
! if-PI > DisplTors >-2*PI    => DisplTors=DisplTors+2*PI
!
      TYPE(INTC)                :: IntCs
      INTEGER                   :: NIntC,I,II
      REAL(DOUBLE),DIMENSION(:) :: VectInt
      REAL(DOUBLE)              :: Sum,ASum,TwoPi
!
      TwoPi=Two*PI
!
      DO I=1,NIntC
        IF(IntCs%Def(I)(1:4)/='STRE'.OR.IntCs%Def(I)(1:4)/='CART') THEN
	  Sum=VectInt(I)
	  II=INT(Sum/TwoPI)
	  Sum=Sum-DBLE(II)*TwoPi
	  VectInt(I)=Sum
	  ASum=ABS(Sum)
	  IF(ASum>PI) THEN
            IF(Sum>Zero) VectInt(I)=VectInt(I)-TwoPi
            IF(Sum<Zero) VectInt(I)=VectInt(I)+TwoPi 
	  ENDIF
	ENDIF
      ENDDO
!
      END SUBROUTINE MapDisplTors
!
!--------------------------------------------------------------------
!
      SUBROUTINE CartRNK1ToCartRNK2(VectCart,ActCarts,Add_O)
!
      REAL(DOUBLE),DIMENSION(:)   :: VectCart
      REAL(DOUBLE),DIMENSION(:,:) :: ActCarts(:,:)
      INTEGER :: I,J,NatmsLoc
      LOGICAL,OPTIONAL :: Add_O
!
      NatmsLoc=SIZE(ActCarts,2)
!
      IF(PRESENT(Add_O)) THEN
        IF(Add_O) THEN
          DO I=1,NatmsLoc 
            J=3*(I-1)
            ActCarts(1,I)=ActCarts(1,I)+VectCart(J+1)
            ActCarts(2,I)=ActCarts(2,I)+VectCart(J+2)
            ActCarts(3,I)=ActCarts(3,I)+VectCart(J+3)
          ENDDO
        ENDIF
      ELSE
        DO I=1,NatmsLoc 
          J=3*(I-1)
          ActCarts(1,I)=VectCart(J+1)
          ActCarts(2,I)=VectCart(J+2)
          ActCarts(3,I)=VectCart(J+3)
        ENDDO
      ENDIF
!
      END SUBROUTINE CartRNK1ToCartRNK2
!
!---------------------------------------------------------------
!
      SUBROUTINE CartRNK2ToCartRNK1(VectCart,ActCarts,Add_O)
!
      REAL(DOUBLE),DIMENSION(:)   :: VectCart
      REAL(DOUBLE),DIMENSION(:,:) :: ActCarts(:,:)
      INTEGER :: I,J,NatmsLoc
      LOGICAL,OPTIONAL :: Add_O
!
      NatmsLoc=SIZE(ActCarts,2)
!
      IF(PRESENT(Add_O)) THEN
        IF(Add_O) THEN
          DO I=1,NatmsLoc 
            J=3*(I-1)
            VectCart(J+1)=VectCart(J+1)+ActCarts(1,I)
            VectCart(J+2)=VectCart(J+2)+ActCarts(2,I)
            VectCart(J+3)=VectCart(J+3)+ActCarts(3,I)
          ENDDO
        ENDIF
      ELSE
        DO I=1,NatmsLoc 
          J=3*(I-1)
          VectCart(J+1)=ActCarts(1,I)
          VectCart(J+2)=ActCarts(2,I)
          VectCart(J+3)=ActCarts(3,I)
        ENDDO
      ENDIF
!
      END SUBROUTINE CartRNK2ToCartRNK1
!-------------------------------------------------------
!
      SUBROUTINE ChkBendToLinB(IntCs,NIntC,XYZ)
!
      TYPE(INTC) :: IntCs,IntC_New
      INTEGER    :: NIntC,Nintc_New
      TYPE(INT_VECT) :: LinAtom,MarkLinb
      REAL(DOUBLE),DIMENSION(:,:) :: XYZ
      REAL(DOUBLE)                :: Value
      INTEGER :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
      INTEGER :: I,J,K,NatmsLoc
      TYPE(INT_RNK2) :: LinBBridge,Top12
!
! Now check for bending -> linear bending transitions
!
      NatmsLoc=SIZE(XYZ,2)
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
      IF(NLinB==0) THEN
        CALL Delete(MarkLinB)
        RETURN
      ENDIF
!     CALL Halt('Bend - LinB transition. Subroutine must be tested')
!
! Modify IntCs by adding new coords (LinB-s).
!
      NIntc_New=NIntc+NLinB
      CALL NEW(IntC_New,NIntc_New)
!
        NLinB=0
      DO I=1,NIntC
        IF(MarkLinB%I(I)==0) THEN
          CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NLinB+I)
        ELSE
          CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NLinB+I)
          IntC_New%Def(NLinB+I)='LINB1'
          NLinB=NLinB+1
          CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NLinB+I)
          IntC_New%Def(NLinB+I)='LINB2'
        ENDIF
      ENDDO
!
      CALL Delete(IntCs)
      NIntC=NIntC_New
      CALL New(IntCs,NIntC)
        CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntC,1)
      CALL Delete(IntC_New)
!
! Now recognize colinear atoms of the molecule and 
! introduce long-range torsions. 
!
        CALL Get(NMax12,'NMax12')
        CALL New(Top12,(/NatmsLoc,NMax12+1/))
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
        CALL Set_INTC_EQ_INTC(IntCs,IntC_New,1,NIntC,1)
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
              NIntC=NIntC+1
              IntC_New%Def(NIntC)='TORS '
              IntC_New%Atoms(NIntC,1)=I1
              IntC_New%Atoms(NIntC,2)=I2
              IntC_New%Atoms(NIntC,3)=I3
              IntC_New%Atoms(NIntC,4)=I4
              IntC_New%Value(NIntC)=Zero
              IntC_New%Constraint(NIntC)=.FALSE.
            ENDIF
            ENDIF
          ENDDO
            ENDIF
          ENDDO
        ENDDO
!
        CALL Delete(Intcs)
        CALL New(Intcs,NIntC)
          CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntC,1)
        CALL Delete(IntC_New)
!
        CALL Delete(LinAtom)
        CALL Delete(LinBBridge)
        CALL Delete(Top12)
!
! Tidy up
!
        CALL Delete(MarkLinB)
!
END SUBROUTINE ChkBendToLinB    
!
!-------------------------------------------------------
!
      SUBROUTINE SetConstraint(IntCs,XYZ,Displ)
!
      TYPE(INTC)     :: IntCs
      INTEGER        :: I,J,NIntC,LastIntcGeom,NDim,JJ
      TYPE(DBL_VECT) :: Displ
      REAL(DOUBLE),DIMENSION(:,:) :: XYZ  
      TYPE(DBL_VECT) :: Constraints
      LOGICAL        :: DoInternals
!
        NIntC=SIZE(IntCs%Def)
        NDim=SIZE(Displ%D)
        DoInternals=GeOpCtrl%DoInternals
        IF(NDim/=NIntC.AND.DoInternals) &
            Call Halt('Dimensionality error in SetConstraint')
!
! Get values of constraints from HDF
!
        CALL New(Constraints,NIntC)
        CALL Get(LastIntcGeom,'LastIntcGeom')
        CALL Get(Constraints,'Constraints'//TRIM(IntToChar(LastIntcGeom)))
!
        IF(DoInternals) THEN
          CALL INTCValue(IntCs,XYZ)
          DO I=1,NIntC
            IF(IntCs%Constraint(I)) THEN
              Displ%D(I)=Constraints%D(I)-IntCs%Value(I)
!write(*,100) i,Constraints%D(I),IntCs%Value(I),Displ%D(I)
!100 format('int= ',i4,3F12.6)
  	    ENDIF
          ENDDO	
        ELSE
          DO I=1,NIntC
            IF(IntCs%Constraint(I)) THEN
              IF(IntCs%Def(I)(1:4)/='CART') CALL Halt('Internal '// &
                 'coord constraints in Cartesian optimization')
                J=IntCs%Atoms(I,1)
              IF(IntCs%Def(I)(1:5)/='CARTX') JJ=(J-1)*3+1
              IF(IntCs%Def(I)(1:5)/='CARTY') JJ=(J-1)*3+2
              IF(IntCs%Def(I)(1:5)/='CARTZ') JJ=(J-1)*3+3
              Displ%D(JJ)=Constraints%D(I)-IntCs%Value(I)
  	    ENDIF
          ENDDO	
        ENDIF
!
! Tidy up
!
        CALL Delete(Constraints)
!
      END SUBROUTINE SetConstraint
!
!-------------------------------------------------------
!
      SUBROUTINE ConstrConv(IntCs,IntDiff,ConstrMax,ConstrRMS)
      TYPE(INTC) :: IntCs
      REAL(DOUBLE),DIMENSION(:) :: IntDiff
      REAL(DOUBLE) :: ConstrMax,ConstrRMS,Sum
      INTEGER      :: I,J,NIntC,NConstr
!
      NIntC=SIZE(IntCs%Def)
!
      ConstrMax=Zero
      ConstrRMS=Zero
      NConstr=0
      DO I=1,NIntC 
        IF(IntCs%Constraint(I)) THEN
          NConstr=NConstr+1
          Sum=ABS(IntDiff(I))
          IF(Sum>ConstrMax) ConstrMax=Sum
          ConstrRMS=ConstrRMS+Sum*Sum
        ENDIF
      ENDDO
      ConstrRMS=SQRT(ConstrRMS/DBLE(NConstr))
!
      END SUBROUTINE ConstrConv
!
!-------------------------------------------------------------
!
      SUBROUTINE CartConv(CartDispl,MaxCartDiff,DiffMax,RMSD,&
                          IntCs,NConstr)
!
      REAL(DOUBLE),DIMENSION(:) :: CartDispl 
      REAL(DOUBLE)              :: MaxCartDiff,DiffMax,RMSD,Sum
      INTEGER                   :: I,J,JJ,NCart,NIntC,NConstr
      TYPE(INTC)                :: IntCs
!
      NCart=SIZE(CartDispl)
      NIntC=SIZE(IntCs%Def)
!
! Make constraints on Cartesians 'hard'
!
      IF(NConstr>0) THEN
        DO I=1,NIntC
          IF(IntCs%Def(I)(1:4)=='CART'.AND.IntCs%Constraint(I)) THEN
            J=IntCs%Atoms(I,1) 
            JJ=(J-1)*3
            IF(IntCs%Def(I)(5:5)=='X') THEN
              JJ=JJ+1
            ELSE IF(IntCs%Def(I)(5:5)=='Y') THEN
              JJ=JJ+2
            ELSE IF(IntCs%Def(I)(5:5)=='Z') THEN
              JJ=JJ+3
            ENDIF
            CartDispl(JJ)=Zero
          ENDIF
        ENDDO
      ENDIF      
!
      DiffMax=Zero
      DO I=1,NCart     
        DiffMax=MAX(DiffMax,ABS(CartDispl(I)))
      ENDDO
!
! Scale Displacements
!
      IF(DiffMax>MaxCartDiff) THEN
        Sum=MaxCartDiff/DiffMax
        CartDispl=Sum*CartDispl
        DiffMax=MaxCartDiff
      ENDIF
!
      RMSD=DOT_PRODUCT(CartDispl,CartDispl)
      RMSD=SQRT(RMSD/DBLE(NCart))
!
      END SUBROUTINE CartConv
!
!-------------------------------------------------------
      SUBROUTINE BackTrfConvg(DoIterate,DiffMax,CooTrfCrit, &
        RMSD,RMSDOld,RMSCrit,ConstrMax,ConstrMaxCrit, &
        ConstrRMS,ConstrRMSOld,NConstr,MaxIt_CooTrf,IStep)
!
        REAL(DOUBLE) :: DiffMax,CooTrfCrit,RMSD,RMSDOld,RMSCrit
        REAL(DOUBLE) :: ConstrMax,ConstrMaxCrit,ConstrRMS
        REAL(DOUBLE) :: ConstrRMSOld            
        INTEGER      :: NConstr,MaxIt_CooTrf,IStep
        LOGICAL      :: DoIterate
        LOGICAL      :: ConvConstr
!
        DoIterate=.TRUE.
        ConvConstr=.TRUE.
!
        IF(NConstr==0) THEN
          ConvConstr=(ConstrMax>ConstrMaxCrit.AND. &
          ConstrRMS<ConstrRMSOld*RMSCrit) 
        ENDIF
!
        DoIterate=(DiffMax>CooTrfCrit.AND.IStep<=MaxIt_CooTrf)
!       DoIterate=(DiffMax>CooTrfCrit.AND.IStep<=MaxIt_CooTrf.AND. &
!                  RMSD<RMSDOld*RMSCrit)
        DoIterate=DoIterate.AND.ConvConstr
!
      END SUBROUTINE BackTrfConvg
!
!----------------------------------------------------------------------
      SUBROUTINE PrtXYZ(XYZ,Title_O,PrtU_O,Convert_O)
        REAL(DOUBLE),DIMENSION(:,:) :: XYZ
        TYPE(DBL_RNK2)              :: XYZPrint
        INTEGER,OPTIONAL            :: PrtU_O
        INTEGER                     :: I,J,NatmsLoc,PrtU
        TYPE(INT_VECT)              :: MMAtNum
        TYPE(DBL_VECT)              :: NuclCharge
        REAL(DOUBLE)                :: Sum
        CHARACTER(LEN=*),OPTIONAL   :: Title_O
        LOGICAL                     :: Opened,Exist
        LOGICAL,OPTIONAL            :: Convert_O    
!
        NatmsLoc=SIZE(XYZ,2)
        PrtU=6 
        IF(PRESENT(PrtU_O)) THEN
          PrtU=PrtU_O
!         CALL OpenAscii('coords',PrtU)
        ENDIF
        CALL New(MMAtNum,NatmsLoc) 
!
#ifdef MMech
        IF(HasMM()) THEN
          CALL Get(MMAtNum,'MMAtNum')
        ELSE
          CALL New(NuclCharge,NatmsLoc)
          CALL Get(NuclCharge,  'atomicnumbers',Tag_O=CurGeom)
          DO I=1,NatmsLoc ; MMAtNum%I(I)=INT(NuclCharge%D(I)) ; ENDDO
          CALL Delete(NuclCharge)
        ENDIF
#else
          CALL New(NuclCharge,NatmsLoc)
          CALL Get(NuclCharge,  'atomicnumbers',Tag_O=CurGeom)
          DO I=1,NatmsLoc ; MMAtNum%I(I)=INT(NuclCharge%D(I)) ; ENDDO
          CALL Delete(NuclCharge)
#endif

!
        CALL New(XYZPrint,(/3,NatmsLoc/))
          XYZPrint%D=XYZ
          IF(PRESENT(Convert_O)) THEN
            IF(Convert_O) THEN
              Sum=One/AngstromsToAu
              XYZPrint%D=Sum*XYZPrint%D
            ENDIF
          ENDIF
!
        WRITE(PrtU,*) NatmsLoc
          IF(PRESENT(Title_O)) WRITE(PrtU,*) Title_O
        DO I=1,NatmsLoc
          WRITE(PrtU,100) MMAtNum%I(I),XYZPrint%D(1:3,I)
        ENDDO
100     FORMAT(I4,3X,3F12.6)
!
        CALL Delete(XYZPrint) 
        CALL Delete(MMAtNum) 
!       CLOSE(PrtU,STATUS='KEEP')
!
      END SUBROUTINE PrtXYZ
!----------------------------------------------------------------------
      SUBROUTINE PutSRStep(XYZ_O,Vect_O,Tag_O,Metric_O)
        INTEGER                              :: IGeom,NCart
        INTEGER                              :: NatmsLoc,J,I
        REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: XYZ_O
        TYPE(DBL_VECT),OPTIONAL              :: Vect_O
        REAL(DOUBLE)                         :: CMX,CMY,CMZ
        REAL(DOUBLE)                         :: MFact
        TYPE(DBL_VECT)                       :: AuxVect
        CHARACTER(LEN=*),OPTIONAL            :: Tag_O
        INTEGER                              :: GDIISMemory
        INTEGER                              :: SRMemory
        INTEGER                              :: RefMemory
        INTEGER                              :: CartGradMemory
        INTEGER                              :: IntGradMemory
        TYPE(DBL_VECT)                       :: Vect
        TYPE(DBL_RNK2)                       :: XYZ 
        LOGICAL                              :: Metric
        LOGICAL,OPTIONAL                     :: Metric_O
!
! Metric of the coordinates
!
          Metric=GeOpCtrl%GDIISMetricOn
          MFact=One
          IF(Metric) MFact=GeOpCtrl%GDIISMetric
          IF(PRESENT(Metric_O)) THEN
              Metric=Metric_O
            IF(Metric) THEN
              MFact=GeOpCtrl%GDIISMetric
            ELSE
              MFact=One
            ENDIF
          ENDIF
!
          IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='CartGrad') THEN
            MFact=One/MFact
          ENDIF
!
! Increment GDIISMemory
!
        IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='SR') THEN
          CALL Get(SrMemory,'SRMemory')
          SRMemory=SRMemory+1
          IGeom=SRMemory
          CALL Put(SRMemory,'SRMemory')
        ENDIF
!
        IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='Ref') THEN
          CALL Get(RefMemory,'RefMemory')
          RefMemory=RefMemory+1
          IGeom=RefMemory
          CALL Put(RefMemory,'RefMemory')
        ENDIF
!
        IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='CartGrad') THEN
          CALL Get(CartGradMemory,'CartGradMemory')
          CartGradMemory=CartGradMemory+1
          IGeom=CartGradMemory
          CALL Put(CartGradMemory,'CartGradMemory')
        ENDIF
!
!       IF(PRESENT(Tag_O).AND.TRIM(Tag_O)=='IntGrad') THEN
!         CALL Get(IntGradMemory,'IntGradMemory')
!         IntGradMemory=IntGradMemory+1
!         IGeom=IntGradMemory
!         CALL Put(IntGradMemory,'IntGradMemory')
!       ENDIF
!
        IF(PRESENT(XYZ_O)) THEN
!
          NatmsLoc=SIZE(XYZ_O,2)
          NCart=3*NatmsLoc
          CALL New(XYZ,(/3,NatmsLoc/))
          IF(Metric) THEN
            XYZ%D=MFact*XYZ_O
          ELSE
            XYZ%D=XYZ_O
          ENDIF
!!!!        CALL CenterOfMass(CMX,CMY,CMZ,XYZ_O=XYZ%D,Move_O=.TRUE.)
!
! Put Geometry of simple relaxation set into HDF, for GDIIS.
!
            CALL New(AuxVect,NCart)
              CALL CartRNK2ToCartRNK1(AuxVect%D,XYZ%D)
              CALL Put(AuxVect,TRIM(Tag_O)//TRIM(IntToChar(IGeom)))
            CALL Delete(AuxVect)
          CALL Delete(XYZ)
!
        ELSE IF(PRESENT(Vect_O)) THEN
!
          NCart=SIZE(Vect_O%D)
          CALL New(Vect,NCart)
          IF(Metric) THEN
            Vect%D=MFact*Vect_O%D
          ELSE
            Vect%D=Vect_O%D
          ENDIF
!!!!      CALL CenterOfMass(CMX,CMY,CMZ,Vect_O=Vect%D,Move_O=.TRUE.)
            CALL Put(Vect,TRIM(Tag_O)//TRIM(IntToChar(IGeom)))
          CALL Delete(Vect)
!
        ELSE
!
          CALL Halt('No input coordinates in PutSRStep')
!
        ENDIF
!
      END SUBROUTINE PutSRStep
!----------------------------------------------------------------------
!
      SUBROUTINE CenterOfMass(CX,CY,CZ,XYZ_O,Vect_O,Move_O)
        REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: XYZ_O
        REAL(DOUBLE),DIMENSION(:)  ,OPTIONAL :: Vect_O
        REAL(DOUBLE)                :: CX,CY,CZ
        INTEGER                     :: I,J,NCart,NatmsLoc
        LOGICAL,OPTIONAL            :: Move_O
!
        IF((PRESENT(XYZ_O).AND.PRESENT(Vect_O)) .OR. &
           (.NOT.PRESENT(XYZ_O).AND..NOT.PRESENT(Vect_O))) THEN
          CALL Halt('Inappropriate input in CenterOfMass')
        ENDIF
!
            CX=Zero
            CY=Zero
            CZ=Zero
        IF(PRESENT(XYZ_O)) THEN
          NatmsLoc=SIZE(XYZ_O,2)
          DO I=1,NatmsLoc
            CX=CX+XYZ_O(1,I) 
            CY=CY+XYZ_O(2,I) 
            CZ=CZ+XYZ_O(3,I) 
          ENDDO
        ELSE IF(PRESENT(Vect_O)) THEN
          NatmsLoc=SIZE(Vect_O)/3
          NCart=3*NatmsLoc
          DO I=1,NatmsLoc
            J=(I-1)*3
            CX=CX+Vect_O(J+1) 
            CY=CY+Vect_O(J+2) 
            CZ=CZ+Vect_O(J+3) 
          ENDDO
        ENDIF
          CX=CX/DBLE(NatmsLoc)
          CY=CY/DBLE(NatmsLoc)
          CZ=CZ/DBLE(NatmsLoc)
!
! Move coordinates into Center of Mass?
!
        IF(PRESENT(Move_O)) THEN
          IF(Move_O) THEN
            IF(PRESENT(XYZ_O)) THEN
              XYZ_O(1,:)=XYZ_O(1,:)-CX
              XYZ_O(2,:)=XYZ_O(2,:)-CY
              XYZ_O(3,:)=XYZ_O(3,:)-CZ
            ELSE IF(PRESENT(Vect_O)) THEN
              DO I=1,NatmsLoc
                J=(I-1)*3
                Vect_O(J+1)=Vect_O(J+1)-CX
                Vect_O(J+2)=Vect_O(J+2)-CY
                Vect_O(J+3)=Vect_O(J+3)-CZ
              ENDDO
            ENDIF
          ENDIF
        ENDIF
!
      END SUBROUTINE CenterOfMass
!
!-------------------------------------------------------------------
!
      SUBROUTINE PrincMomTens(CMCarts,Theta)
!
! Calculates inertial momentum tensor
! Pass in mass-centered Cartesians
!
        REAL(DOUBLE),DIMENSION(:,:) :: Theta,CMCarts
        INTEGER                     :: I,NatmsLoc
        REAL(DOUBLE)                :: X,Y,Z,XX,YY,ZZ,XY,YZ,ZX
!                
      NatmsLoc=SIZE(CMCarts,2)
!                
        Theta=Zero
      DO I=1,NatmsLoc
        X=CMCarts(1,I)
        Y=CMCarts(2,I)
        Z=CMCarts(3,I)
	XX=X*X
	YY=Y*Y
	ZZ=Z*Z
	XY=X*Y
	YZ=Y*Z
	ZX=Z*X
        Theta(1,1)=Theta(1,1)+YY+ZZ
        Theta(2,2)=Theta(2,2)+ZZ+XX
        Theta(3,3)=Theta(3,3)+XX+YY
        Theta(1,2)=Theta(1,2)-XY    
        Theta(2,1)=Theta(2,1)-XY    
        Theta(1,3)=Theta(1,3)-ZX    
        Theta(3,1)=Theta(3,1)-ZX    
        Theta(2,3)=Theta(2,3)-YZ    
        Theta(3,2)=Theta(3,2)-YZ    
      ENDDO
!
      END SUBROUTINE PrincMomTens
!
!------------------------------------------------------------------
!
      SUBROUTINE RotateMol(XYZ,U33)
!
! Rotate Cartesian set by U33 Unitary matrix
! XYZ(I)=XYZ(I)*U33
!
        REAL(DOUBLE),DIMENSION(:,:) :: XYZ,U33
        INTEGER                     :: I,J,K,NatmsLoc
        TYPE(DBL_VECT)              :: AuxVect
!
        NatmsLoc=SIZE(XYZ,2)
!
        CALL New(AuxVect,3)
        DO I=1,NatmsLoc
          CALL DGEMM_NNc(1,3,3,One,Zero,XYZ(1:3,I),U33,AuxVect%D)
          XYZ(1:3,I)=AuxVect%D 
        ENDDO
        CALL Delete(AuxVect)
!
      END SUBROUTINE RotateMol
!
!-------------------------------------------------------------------
!
      SUBROUTINE CROSS_PRODUCT(V1,V2,CrossProd)
      REAL(DOUBLE),DIMENSION(:) :: V1,V2,CrossProd
!
        CrossProd(1)=V1(2)*V2(3)-V1(3)*V2(2)
        CrossProd(2)=V1(3)*V2(1)-V1(1)*V2(3)
        CrossProd(3)=V1(1)*V2(2)-V1(2)*V2(1)
!
      END SUBROUTINE CROSS_PRODUCT
!
!------------------------------------------------------------------
!
      SUBROUTINE GDIISSelect(EigVals,Selection)
        REAL(DOUBLE),DIMENSION(:) :: EigVals
        INTEGER,     DIMENSION(:) :: Selection
        REAL(DOUBLE)              :: GapWidth,MeanWidth,Fluct
        REAL(DOUBLE)              :: MaxEig,MinEig
        INTEGER                   :: I,J,K,L,NDim,LastDomain
        REAL(DOUBLE)              :: GDIISBandWidth
        REAL(DOUBLE)              :: Sum,SumX,Percent
        INTEGER                   :: MinDomCount    
        TYPE(DBL_VECT)            :: AuxVect
        TYPE(INT_VECT)            :: Ordered,Domain,DomainCount
        REAL(DOUBLE)              :: GDIISMaxMem
!
        GDIISMaxMem=GeOpCtrl%GDIISMaxMem
        GDIISBandWidth=GeOpCtrl%GDIISBandWidth 
        MinDomCount=GeOpCtrl%GDIISMinDomCount
!
        NDim=SIZE(EigVals)
        I=SIZE(Selection)
        IF(NDim/=I.OR.NDim==0) THEN
          CALL Halt('Dimensionality error in GDIISSelect')
        ENDIF
!
        CALL New(AuxVect,NDim)
        CALL New(Ordered,NDim)
        CALL New(Domain,NDim)
        CALL New(DomainCount,NDim)
!
! Put eigenvalues onto a logarithmic scale
!
          DO I=1,NDim
            IF(Selection(I)==1) THEN  
              AuxVect%D(I)=DLOG10(ABS(EigVals(I)))
            ELSE
              AuxVect%D(I)=-1.D99
            ENDIF
          ENDDO
!
! Form ordered set with lowest log10 first
!
          DO I=1,NDim
            Ordered%I(I)=I 
          ENDDO
          DO I=2,NDim
          DO J=2,NDim-I+2
            K=Ordered%I(J)
            L=Ordered%I(J-1)
            IF(AuxVect%D(K)<AuxVect%D(L)) THEN
              Ordered%I(J)=L
              Ordered%I(J-1)=K
            ENDIF
          ENDDO
          ENDDO
!
! Calculate average spectral distance and fluctuation
!
            MeanWidth=Zero
            J=0
          DO I=2,NDim
            K=Ordered%I(I)
            L=Ordered%I(I-1)
            IF(Selection(K)/=1.OR.Selection(L)/=1) CYCLE 
            J=J+1
            MeanWidth=MeanWidth+AuxVect%D(K)-AuxVect%D(L)
          ENDDO
            IF(J>0) THEN
              MeanWidth=MeanWidth/DBLE(J)
            ELSE
              MeanWidth=1.D99
            ENDIF
!
            Fluct=Zero
          DO I=2,NDim
            K=Ordered%I(I)
            L=Ordered%I(I-1)
            IF(Selection(K)/=1.OR.Selection(L)/=1) CYCLE 
            Fluct=Fluct+ABS(MeanWidth-(AuxVect%D(K)-AuxVect%D(L)))
          ENDDO
            IF(J/=0) THEN
              Fluct=Fluct/J
            ELSE
              Fluct=1.D99
            ENDIF
!
! First, group log10-s into domains
!
            DomainCount%I=0 
            Domain%I(1)=1 
            DomainCount%I(1)=1 
!
          DO I=2,NDim
            K=Ordered%I(I)
            L=Ordered%I(I-1)
            IF(Selection(K)/=1) THEN
              Domain%I(I)=1  
              DomainCount%I(1)=DomainCount%I(1)+1
              CYCLE
            ENDIF
            IF(AuxVect%D(K)-AuxVect%D(L)>MeanWidth+1.01D0*Fluct) THEN
              Domain%I(I)=Domain%I(I-1)+1
            ELSE
              Domain%I(I)=Domain%I(I-1)
            ENDIF
              K=Domain%I(I)
              DomainCount%I(K)=DomainCount%I(K)+1
          ENDDO
              LastDomain=Domain%I(NDim)
!
              K=Ordered%I(NDim)
              L=Ordered%I(NDim-1)
          IF(LastDomain>1.AND.DomainCount%I(LastDomain)==1) THEN
            IF(AuxVect%D(K)-AuxVect%D(L)>0.9D0) THEN
              IF(DomainCount%I(LastDomain)==1) LastDomain=LastDomain-1
            ELSE
              LastDomain=LastDomain-1
              Domain%I(K)=LastDomain
            ENDIF
          ENDIF
!
! Find maximum coefficient in LastDomain
!
          MaxEig=-1.D99
          DO I=1,NDim
            K=Ordered%I(I)
            IF(Selection(K)/=1) CYCLE
!!!!        IF(Domain%I(I)/=LastDomain) CYCLE
!           MaxEig=MAX(MaxEig,AuxVect%D(K))
            MaxEig=MAX(MaxEig,EigVals(K))
          ENDDO
!
             J=0
             SumX=Zero
          DO I=NDim,1,-1
             K=Ordered%I(I)
!!!!         IF(Selection(K)/=1.OR.Domain%I(I)/=LastDomain) CYCLE
             IF(Selection(K)/=1) CYCLE
             Sum=AuxVect%D(K)
             SumX=SumX+EigVals(K)
             Percent=ABS(EigVals(K))/ABS(SumX)*100.D0
!write(*,50) i,k,EigVals(K),Percent
!50 format(2I4,2X,E12.6,' percent= ',F8.2)
            IF(Percent>GDIISBandWidth) THEN
              J=J+1
              Selection(K)=1
            ELSE
              Selection(K)=0
            ENDIF
            IF(J>=GDIISMaxMem) EXIT
          ENDDO
!
1000    CONTINUE
!        
        CALL Delete(DomainCount)
        CALL Delete(Domain)        
        CALL Delete(Ordered)        
        CALL Delete(AuxVect)        
!
      END SUBROUTINE GDIISSelect
!
!----------------------------------------------------------------
!
      SUBROUTINE GDIISSave(DimGDIIS,EtotCurr,EtotPrev)
        INTEGER                :: I,J,IGeom,SRMemory,RefMemory
        INTEGER                :: GDIISMemory,DimGDIIS,CartGradMemory
        TYPE(DBL_VECT)         :: Auxvect
        REAL(DOUBLE)           :: EtotCurr,EtotPrev
!
        IF(GeOpCtrl%ActStep<=GeOpCtrl%GDIISInit) RETURN
!
        CALL New(Auxvect,DimGDIIS)
!
! If step was unsuccesful restore last succesful state of 
! GDIIS reference space
!
        IF(EtotCurr>EtotPrev.AND.GeOpCtrl%ActStep/=GeOpCtrl%GDIISInit) THEN  
!
          CALL Get(SRMemory,'SRMemory'//'SAVE')
          CALL Get(RefMemory,'RefMemory'//'SAVE')
          CALL Get(CartGradMemory,'CartGradMemory'//'SAVE')
          IF(SRMemory/=RefMemory) THEN
            CALL Halt('SRMemory/=RefMemory in GDIISSave')
          ENDIF
          IF(SRMemory/=CartGradMemory) THEN
            CALL Halt('SRMemory/=CartGradMemory in GDIISSave')
          ENDIF
          GDIISMemory=SRMemory
!
          CALL Put(SRMemory,'SRMemory')
          CALL Put(RefMemory,'RefMemory')
          CALL Put(CartGradMemory,'CartGradMemory')
          DO IGeom=1,GDIISMemory
            CALL Get(AuxVect,'SR'//TRIM(IntToChar(IGeom))//'SAVE')
            CALL Put(AuxVect,'SR'//TRIM(IntToChar(IGeom)))
            CALL Get(AuxVect,'Ref'//TRIM(IntToChar(IGeom))//'SAVE')
            CALL Put(AuxVect,'Ref'//TRIM(IntToChar(IGeom)))
            CALL Get(AuxVect,'CartGrad'//TRIM(IntToChar(IGeom))//'SAVE')
            CALL Put(AuxVect,'CartGrad'//TRIM(IntToChar(IGeom)))
          ENDDO
!
        ELSE 
!
! Get recent Ref and SR structures
!
          CALL Get(SRMemory,'SRMemory')
          CALL Get(RefMemory,'RefMemory')
          CALL Get(CartGradMemory,'CartGradMemory')
          IF(SRMemory/=RefMemory) THEN
            CALL Halt('SRMemory/=RefMemory in GDIISSave')
          ENDIF
          IF(SRMemory/=CartGradMemory) THEN
            CALL Halt('SRMemory/=CartGradMemory in GDIISSave')
          ENDIF
          GDIISMemory=SRMemory
!
          CALL Put(SRMemory,'SRMemory'//'SAVE')
          CALL Put(RefMemory,'RefMemory'//'SAVE')
          CALL Put(CartGradMemory,'CartGradMemory'//'SAVE')
          DO IGeom=1,GDIISMemory
            CALL Get(AuxVect,'SR'//TRIM(IntToChar(IGeom)))
            CALL Put(AuxVect,'SR'//TRIM(IntToChar(IGeom))//'SAVE')
            CALL Get(AuxVect,'Ref'//TRIM(IntToChar(IGeom)))
            CALL Put(AuxVect,'Ref'//TRIM(IntToChar(IGeom))//'SAVE')
            CALL Get(AuxVect,'CartGrad'//TRIM(IntToChar(IGeom)))
            CALL Put(AuxVect,'CartGrad'//TRIM(IntToChar(IGeom))//'SAVE')
          ENDDO
!
        ENDIF
!
        CALL Delete(AuxVect)
!
      END SUBROUTINE GDIISSave
!
!-------------------------------------------------------
!
    SUBROUTINE GDIIS2(GMLoc)
    IMPLICIT NONE
    TYPE(CRDS)     :: GMLoc
    INTEGER        :: I,J,K,L,GDIISMemoryIn,DimGDIIS,IGeom,DimOverl
    INTEGER        :: SRMemory,RefMemory,CartGradMemory
    INTEGER        :: LastStruct,ILow,NCart,GDIISMemory
    INTEGER        :: INFO,NewDim 
    TYPE(DBL_RNK2) :: SRDispl,AMat,SRStruct,ActCarts
    TYPE(DBL_RNK2) :: TrfDispl,TrfGrad,EigVect
    TYPE(DBL_RNK2) :: RefGrad
    TYPE(DBL_RNK2) :: AuxStruct,RefStruct
    TYPE(DBL_VECT) :: Coeffs,AuxVect,AuxVectS,Scale,EigVal
    TYPE(DBL_VECT) :: GrdDotGrd,GrdDotDX,DXDotDX,MetricA
    TYPE(INT_VECT) :: Selection
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
    REAL(DOUBLE)                   :: Sum,S1,S2
    REAL(DOUBLE)                   :: Eig,MFact,MaxEig
    REAL(DOUBLE)                   :: SumX,SumY,SumZ
    REAL(DOUBLE)                   :: MinGrdDotGrd,MinGrdDotDX
    REAL(DOUBLE)                   :: MinDXDotDX,ScaleMin
    REAL(DOUBLE)                   :: EDI,ED0,EDT 
!
! Current implementation runs on Cartesian displacements only
! Input GMLoc contains the actual geometry
! CurGeom is set to the last (simple relaxation) structure
!
      CALL Get(SRMemory,'SRMemory')
      CALL Get(RefMemory,'RefMemory')
      CALL Get(CartGradMemory,'CartGradMemory')
      NCart=3*GMLoc%Natms
      DimGDIIS=NCart   !!! later dimension may become NIntC
      IF(SRMemory/=RefMemory) CALL Halt('SRMemory/=RefMemory in GDIIS')
      IF(SRMemory/=CartGradMemory) CALL Halt('SRMemory/=CartGradMemory in GDIIS')
      GDIISMemory=SRMemory
!
      IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
        CALL OpenAscii(OutFile,Out)
        WRITE(Out,*) 'GDIISMemory= ',GDIISMemory
        WRITE(*,*) 'GDIISMemory= ',GDIISMemory
        CLOSE(Unit=Out,STATUS='KEEP')
      ENDIF
!
    GMTag=''
#ifdef MMech
    IF(HasMM()) GMTag='GM_MM'
#endif
!
! Get Displacements from geometries, stored in HDF
! and fill them into SRDispl columns.
!
      CALL New(SRDispl,(/DimGDIIS,GDIISMemory/))
      CALL New(SRStruct,(/DimGDIIS,GDIISMemory/))
      CALL New(RefStruct,(/DimGDIIS,GDIISMemory/))
      CALL New(RefGrad,(/DimGDIIS,GDIISMemory/))
      SRStruct%D=Zero
      CALL New(AuxVect,DimGDIIS)
      CALL New(Scale,GDIISMemory)
      CALL New(AuxVectS,GDIISMemory)
      CALL New(GrdDotDX,GDIISMemory)
      CALL New(DXDotDX,GDIISMemory)
      CALL New(GrdDotGrd,GDIISMemory)
      CALL New(Selection,GDIISMemory)
      Selection%I=1
      CALL New(MetricA,GDIISMemory)
      CALL New(Coeffs,GDIISMemory)
!
! Get recent Ref and SR structures, as well as the Gradients
!
      DO IGeom=1,GDIISMemory
        CALL Get(AuxVect,'SR'//TRIM(IntToChar(IGeom)))
        SRStruct%D(:,IGeom)=AuxVect%D
        CALL Get(AuxVect,'Ref'//TRIM(IntToChar(IGeom)))
        RefStruct%D(:,IGeom)=AuxVect%D
        CALL Get(AuxVect,'CartGrad'//TRIM(IntToChar(IGeom)))
        RefGrad%D(:,IGeom)=AuxVect%D
      ENDDO
!
! Get simple relaxation displacements ('error vectors')
! and normalize them!
!
      SRDispl%D=SRStruct%D-RefStruct%D
!
! Now calculate overlap ('A' matrix). 
! Presently only non-sparse representation is available
!
      CALL New(AMat,(/GDIISMemory,GDIISMemory/))
!
!     CALL DGEMM_TNc(GDIISMemory,DimGDIIS,GDIISMemory,One,Zero,  &
!                  SRDispl%D,SRDispl%D,AMat%D)
      CALL DGEMM_TNc(GDIISMemory,DimGDIIS,GDIISMemory,One,Zero,  &
                   RefGrad%D,RefGrad%D,AMat%D)
!
! Scale AMat before diagonalization
!
         Scale%D=Zero
         ScaleMin=1.D99
       DO I=1,GDIISMemory
         Sum=One/SQRT(AMat%D(I,I))
         Scale%D(I)=Sum
         ScaleMin=MIN(ScaleMin,Sum)
       ENDDO
!scale%d=one
!ScaleMin=One
       DO I=1,GDIISMemory
         SumX=Scale%D(I)
         DO J=1,GDIISMemory
           SumY=Scale%D(J)
           AMat%D(I,J)=SumX*AMat%D(I,J)*SumY 
         ENDDO
       ENDDO
!
! Now, calculate eigenvalues and eigenvectors of Overlap
!
      CALL SetDSYEVWork(GDIISMemory)
!
        BLKVECT%D=AMat%D
        CALL DSYEV('V','U',GDIISMemory,BLKVECT%D,BIGBLOK,BLKVALS%D, &
          BLKWORK%D,BLKLWORK,INFO)
        IF(INFO/=SUCCEED) &
        CALL Halt('DSYEV hosed in GDIIS. INFO='&
                   //TRIM(IntToChar(INFO)))
!
! Now, construct adaptive reference if requested
!
!       CALL TrfGDIIS(RefStruct,SRStruct,SRDispl,RefGrad,AMat,Scale, &
!                    BLKVECT,BLKVALS)
!
! Transform iterative subspace into new metric
!
        IF(GeOpCtrl%GDIISMetricOn) THEN
          CALL TrfMetric(RefGrad,RefStruct,SRStruct,SRDispl,BLKVALS)
        ENDIF
!
! Rescale eigenvalues
!
        Sum=One/BLKVALS%D(1)
        BLKVALS%D=Sum*BLKVALS%D
!
! Rescale scalefactors
!
          Sum=One/ScaleMin
          Scale%D=Sum*Scale%D
!
! Carry out Ut*Scale
!
        CALL DGEMM_TNc(GDIISMemory,GDIISMemory,1,One,Zero,  &
          BLKVECT%D,Scale%D,AuxVectS%D)
!
! Scale Ut*Scale
!
          SumX=1.D99
        DO I=1,GDIISMemory
          Sum=AuxVectS%D(I)
          SumX=MIN(SumX,ABS(Sum))
        ENDDO
          SumX=One/SumX
          AuxVectS%D=SumX*AuxVectS%D
!
! Calculate pre-coeffs in scaled metric
!
          SumX=1.D99
        DO I=1,GDIISMemory
          Coeffs%D(I)=AuxVectS%D(I)/BLKVALS%D(I)
          SumX=MIN(SumX,Coeffs%D(I))
        ENDDO
          SumX=One/SumX
          Coeffs%D=SumX*Coeffs%D
!
! Calculate coeffs in the scaled metric
!
        CALL DGEMM_NNc(GDIISMemory,GDIISMemory,1,One,Zero,  &
          BLKVECT%D,Coeffs%D,AuxVectS%D)
          Coeffs%D=AuxVectS%D
!
! Calculate coeffs in the original metric
!
        DO I=1,GDIISMemory
          Coeffs%D(I)=Scale%D(I)*Coeffs%D(I)
        ENDDO
!
! Select out subspace for final GDIIS
!
        CALL GDIISSelect(Coeffs%D,Selection%I)
!
      CALL UnSetDSYEVWork()
!
! Rescale coeffs to get a sum of One.
!
          Sum=Zero
        DO I=1,GDIISMemory 
          IF(Selection%I(I)==1) THEN 
            Sum=Sum+Coeffs%D(I)  
          ELSE
            Coeffs%D(I)=Zero
          ENDIF
        ENDDO
        Sum=One/Sum
        Coeffs%D=Sum*Coeffs%D
!
! Calculate new geometry, linearcombine previous steps 
!
          AuxVect%D=Zero
        DO I=1,GDIISMemory
          IF(Selection%I(I)/=1) CYCLE
          AuxVect%D=AuxVect%D+Coeffs%D(I)*SRStruct%D(:,I)
        ENDDO
!
! Restore original metric on resulting GDIIS geometry
!
        IF(GeOpCtrl%GDIISMetricOn) THEN
          MFact=One/GeOpCtrl%GDIISMetric
          AuxVect%D=MFact*AuxVect%D
        ENDIF
!
! Fill new geometry into GM array
!
        CALL CartRNK1ToCartRNK2(AuxVect%D,GMLoc%Carts%D)
!
! Save unitary transformed structures and gradients
!
        CALL Put(0,'SrMemory')
        CALL Put(0,'RefMemory')
        CALL Put(0,'CartGradMemory')
!       CALL Put(0,'IntGradMemory')
        J=MAX(1,GDIISMemory-5)
        DO I=J,GDIISMemory
!       DO I=1,GDIISMemory
          IF(Selection%I(I)==1) THEN
            AuxVect%D=SRStruct%D(:,I)
            CALL PutSRStep(Vect_O=AuxVect,Tag_O='SR',Metric_O=.FALSE.)
            AuxVect%D=RefStruct%D(:,I)
            CALL PutSRStep(Vect_O=AuxVect,Tag_O='Ref',Metric_O=.FALSE.)
            AuxVect%D=RefGrad%D(:,I)
            CALL PutSRStep(Vect_O=AuxVect,Tag_O='CartGrad',&
                           Metric_O=.FALSE.)
          ENDIF
        ENDDO
!
!call prtxyz(GMLoc%Carts%D,Title_O='Final Carts.')
!
! Tidy up
!
       CALL Delete(Selection)
       CALL Delete(GrdDotGrd)
       CALL Delete(GrdDotDX)
       CALL Delete(DXDotDX)
       CALL Delete(AuxVectS)
       CALL Delete(Scale)
       CALL Delete(AuxVect)
       CALL Delete(Coeffs)
       CALL Delete(MetricA)
       CALL Delete(AMat)
       CALL Delete(RefGrad)
       CALL Delete(RefStruct)
       CALL Delete(SRStruct)
       CALL Delete(SRDispl)
!
     END SUBROUTINE GDIIS2
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
!----------------------------------------------------------------
!
      SUBROUTINE TrfMetric(RefGrad,RefStruct,SRStruct,SRDispl,BLKVALS)
        TYPE(DBL_RNK2) :: RefGrad,RefStruct,SRStruct,SRDispl
        TYPE(DBL_VECT) :: BLKVALS
        INTEGER        :: I,J,DimGDIIS,GDIISMemory
        REAL(DOUBLE)   :: MFact,Sum
!
        DimGDIIS   =SIZE(RefGrad%D,1)
        GDIISMemory=SIZE(RefGrad%D,2)
!
        Sum=SQRT(BLKVALS%D(1))
        MFact=One/Sum
        GeOpCtrl%GDIISMetric=Sum*GeOpCtrl%GDIISMetric
!
        RefGrad%D=MFact*RefGrad%D
        RefStruct%D=Sum*RefStruct%D
        SRStruct%D=Sum*SRStruct%D
        SRDispl%D=Sum*SRDispl%D
        MFact=MFact*MFact
        BLKVALS%D=MFact*BLKVALS%D
!
      END SUBROUTINE TrfMetric
!
!--------------------------------------------------------------------
!
      SUBROUTINE TrfGDIIS(RefStruct,SRStruct,SRDispl,RefGrad,AMat,&
                     Scale,TrfMatr,EigVals)
        TYPE(DBL_RNK2) :: RefStruct,SRStruct,SRDispl,RefGrad,AMat,TrfMatr
        TYPE(DBL_RNK2) :: TrfGrad,TrfDispl
        TYPE(DBL_VECT) :: Scale,EigVals,ColmnSum,ColmnSum2
        TYPE(DBL_RNK2) :: DiagMatr,AuxStruct
        INTEGER        :: I,J,K,GDIISMemory,DimGDIIS
        REAL(DOUBLE)   :: Sum,Sum2,SumX
        LOGICAL        :: WhichTrf
!
        DimGDIIS=SIZE(RefStruct%D,1)
        GDIISMemory=SIZE(RefStruct%D,2)
        CALL New(ColmnSum,GDIISMemory)
        CALL New(ColmnSum2,GDIISMemory)
        CALL New(AuxStruct,(/DimGDIIS,GDIISMemory/))
        CALL New(TrfGrad,(/GDIISMemory,GDIISMemory/))
        CALL New(TrfDispl,(/GDIISMemory,GDIISMemory/))
call pprint(amat%d,' amat scaled',unit_o=6)
sum=zero
do i=1,gdiismemory
sum=sum+dot_product(TrfGrad%D(:,i),SRDispl%D(:,i))
enddo
write(*,*) 'energy lowering in original system= ',sum
!
! Carry out Ut*Scale matrix matrix multipl. (Rows of U are multiplied.)
! Scale is normalized at this point 
! such that its lowest component is One.
!
call pprint(trfmatr%d,'TrfMatr original ',unit_o=6)
        DO I=1,GDIISMemory
!         Sum=One/Scale%D(I)
          Sum=Scale%D(I)
          TrfMatr%D(I,:) =Sum*TrfMatr%D(I,:)
        ENDDO
call pprint(trfmatr%d,'TrfMatr*scale ',unit_o=6)
!
! Now, calculate sums of column-components
! invert them and normalize the result
!
          ColmnSum%D=Zero
        DO I=1,GDIISMemory
          Sum=Zero
          DO J=1,GDIISMemory
            Sum=Sum+TrfMatr%D(J,I)
          ENDDO
            Sum2=DOT_PRODUCT(TrfMatr%D(:,I),TrfMatr%D(:,I))
write(*,*) i,' Sum= ',Sum
          ColmnSum%D(I)=Sum
          ColmnSum2%D(I)=Sum2
        ENDDO
write(*,*) 'ColmnSum= ',ColmnSum%D
write(*,*) 'ColmnSum2= ',ColmnSum2%D
!
! Transform Cartesian coordinates and displacements
!
         CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero, &
           RefStruct%D,TrfMatr%D,AuxStruct%D)
           RefStruct%D=AuxStruct%D
         CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero, &
           SRStruct%D,TrfMatr%D,AuxStruct%D)
           SRStruct%D=AuxStruct%D
         CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero, &
           SRDispl%D,TrfMatr%D,AuxStruct%D)
           SRDispl%D=AuxStruct%D
         DO I=1,GDIISMemory
           Sum=One/ColmnSum%D(I)
           RefStruct%D(:,I)=Sum*RefStruct%D(:,I)
           SRStruct%D(:,I) =Sum*SRStruct%D(:,I)
           SRDispl%D(:,I)  =Sum*SRDispl%D(:,I)
         ENDDO
call pprint(refstruct%d,'refstruct transformed ',unit_o=6)
call pprint(SRstruct%d,'srstruct transformed ',unit_o=6)
call pprint(srdispl%d,'srdispl transformed ',unit_o=6)
!
! Transform gradients
!
         CALL DGEMM_NNc(DimGDIIS,GDIISMemory,GDIISMemory,One,Zero, &
           RefGrad%D,TrfMatr%D,AuxStruct%D)
           RefGrad%D=AuxStruct%D
         DO I=1,GDIISMemory
!          Sum=ColmnSum%D(I)*ColmnSum%D(I)/ColmnSum2%D(I)
           Sum=One/ColmnSum%D(I)
           RefGrad%D(:,I)=Sum*RefGrad%D(:,I)
         ENDDO
call pprint(refgrad%d,'refgrad transformed ',unit_o=6)
CALL DGEMM_TNc(GDIISMemory,DimGDIIS,GDIISMemory,One,Zero, &
  RefGrad%D,RefGrad%D,AMat%D)
call pprint(amat%d,'grad overlap new ',unit_o=6)
sum=zero
do i=1,gdiismemory
sum=sum+dot_product(TrfGrad%D(:,i),SRDispl%D(:,i))
enddo
write(*,*) 'energy lowering in transformed system= ',sum
stop
!
        CALL Delete(TrfDispl)
        CALL Delete(TrfGrad)
        CALL Delete(AuxStruct)
        CALL Delete(ColmnSum2)
        CALL Delete(ColmnSum)
!
      END SUBROUTINE TrfGDIIS
!
!--------------------------------------------------------------------
!
END MODULE InCoords

