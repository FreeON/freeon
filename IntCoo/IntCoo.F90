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
SUBROUTINE TOPOLOGY_12(NATOMS,NBONDS,BONDI,BONDJ,TOP12,InfFile)
! Set up a table which shows the atom numbers of atoms 
! connected to a certain atom by the input bonds (topology mtr)
! Here, the generation of the topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL,INTENT(INOUT) :: TOP12
TYPE(INT_RNK2) :: TOP12_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBONDS,NATOMS,NMAX12
INTEGER,DIMENSION(1:NBONDS) :: BONDI,BONDJ
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    NMAX12=5
    K=NMAX12+1
    CALL NEW(TOP12,(/NATOMS,K/))
    TOP12%I(1:NATOMS,1:NMAX12+1)=0 
!
    DO I=1,NBONDS
!
      II=BONDI(I)
      JJ=BONDJ(I)
      NI=TOP12%I(II,1)
      NJ=TOP12%I(JJ,1)
!
! check matrix size, increase size if necessary
      IF(NI>NMAX12 .OR. NJ>NMAX12) THEN
        NMAX12=NMAX12+5
        CALL NEW(TOP12_2,(/NATOMS,NMAX12+1/))
        TOP12_2%I(1:NATOMS,1:NMAX12+1)=0 
        TOP12_2%I(1:NATOMS,1:NMAX12+1-5)=TOP12%I(1:NATOMS,1:NMAX12+1-5)
        CALL DELETE(TOP12)
        CALL NEW(TOP12,(/NATOMS,NMAX12+1/))
        TOP12%I(:,:)=TOP12_2%I(:,:)
        CALL DELETE(TOP12_2)
      ENDIF
!
      IF(NI/=0) THEN
      DO M=1,NI         
        IF(TOP12%I(II,M+1)==JJ) THEN
          EXIT
        ELSE
          TOP12%I(II,1)=NI+1
          TOP12%I(II,1+(NI+1))=JJ
          EXIT
        ENDIF
      ENDDO 
      ELSE
          TOP12%I(II,1)=NI+1
          TOP12%I(II,1+(NI+1))=JJ
      ENDIF
!
      IF(NJ/=0) THEN
      DO M=1,NJ         
        IF(TOP12%I(JJ,M+1)==II) THEN
          EXIT
        ELSE
          TOP12%I(JJ,1)=NJ+1
          TOP12%I(JJ,1+(NJ+1))=II
          EXIT
        ENDIF
      ENDDO 
      ELSE
          TOP12%I(JJ,1)=NJ+1
          TOP12%I(JJ,1+(NJ+1))=II
      ENDIF
!
    ENDDO
!
    IF(PRESENT(InfFile)) THEN
      CALL OpenHDF(InfFile)
      CALL Put(NMAX12,'NMAX12')
      CALL Put(TOP12,'TOP12')
      CALL CloseHDF()
    ENDIF
    IF(.NOT.PRESENT(TOP12)) THEN
      CALL DELETE(TOP12)
    ENDIF
!

END SUBROUTINE TOPOLOGY_12 
!--------------------------------------------------------------
!
SUBROUTINE TOPOLOGY_13(NATOMS,TOP12,TOP13,InfFile)
! Set up a table which shows the atom numbers of atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL,INTENT(IN) :: TOP12
TYPE(INT_RNK2),OPTIONAL,INTENT(OUT) :: TOP13
TYPE(INT_RNK2) :: TOP13_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NATOMS,NMAX13,NMAX12,KK,IN12,JN12
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    IF(.NOT.PRESENT(TOP12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL OpenHDF(InfFile)
        CALL Get(NMAX12,'NMAX12')
        K=NMAX12+1
        CALL NEW(TOP12,(/NATOMS,K/))
        CALL Get(TOP12,'TOP12')
        CALL CloseHDF()
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing top_12 matrix')
      ENDIF
    ENDIF
!
    NMAX13=10
    K=NMAX13+1
    CALL NEW(TOP13,(/NATOMS,K/))
    TOP13%I(1:NATOMS,1:NMAX13+1)=0 
!
    DO II=1,NATOMS
      IN12=TOP12%I(II,1)
      DO J=1,IN12
        JJ=TOP12%I(II,J+1)
        JN12=TOP12%I(JJ,1)
          DO K=1,JN12
          KK=TOP12%I(JJ,K+1)
          IF(II/=KK) THEN
!
! continue here
!
      NI=TOP13%I(II,1)
!
! check matrix size, increase size if necessary
      IF(NI>NMAX13) THEN
        NMAX13=NMAX13+10
        CALL NEW(TOP13_2,(/NATOMS,NMAX13+1/))
        TOP13_2%I(1:NATOMS,1:NMAX13+1)=0 
        TOP13_2%I(1:NATOMS,1:NMAX13+1-10)=TOP13%I(1:NATOMS,1:NMAX13+1-10)
        CALL DELETE(TOP13)
        CALL NEW(TOP13,(/NATOMS,NMAX13+1/))
        TOP13%I(1:NATOMS,1:NMAX13+1)=TOP13_2%I(1:NATOMS,1:NMAX13+1)
        CALL DELETE(TOP13_2)
      ENDIF
!
      IF(NI/=0) THEN
      DO M=1,NI         
        IF(TOP13%I(II,M+1)==KK) THEN
          EXIT
        ELSE
          TOP13%I(II,1)=NI+1
          TOP13%I(II,1+(NI+1))=KK
          EXIT
        ENDIF
      ENDDO 
      ELSE
          TOP13%I(II,1)=NI+1
          TOP13%I(II,1+(NI+1))=KK
      ENDIF
!
        ENDIF !!! II/=KK
        ENDDO !!!! KK
      ENDDO !!!! JJ
    ENDDO !!!! II
!
    IF(PRESENT(InfFile)) THEN
      CALL OpenHDF(InfFile)
      CALL Put(NMAX13,'NMAX13')
      CALL Put(TOP13,'TOP13')
      CALL CloseHDF()
    ENDIF
!
    IF(.NOT.PRESENT(TOP12)) CALL DELETE(TOP12)
    IF(.NOT.PRESENT(TOP13)) CALL DELETE(TOP13)
!
END SUBROUTINE TOPOLOGY_13 
!--------------------------------------------------------------
!
SUBROUTINE TOPOLOGY_14(NATOMS,TOP12,TOP14,InfFile)
! Set up a table which shows the atom numbers of atoms 
! beeing second neighbours of a certain atom.
!
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL,INTENT(IN) :: TOP12
TYPE(INT_RNK2),OPTIONAL,INTENT(OUT) :: TOP14
TYPE(INT_RNK2) :: TOP14_2
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
INTEGER :: NATOMS,NMAX14,NMAX12,IN12,JN12,KN12
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
    IF(.NOT.PRESENT(TOP12)) THEN
      IF(PRESENT(InfFile)) THEN
        CALL OpenHDF(InfFile)
        CALL Get(NMAX12,'NMAX12')
        K=NMAX12+1
        CALL NEW(TOP12,(/NATOMS,K/))
        CALL Get(TOP12,'TOP12')
        CALL CloseHDF()
      ELSE
        CALL MondoHalt(INTC_ERROR,'Missing top_12 matrix')
      ENDIF
    ENDIF
!
    NMAX14=10
    K=NMAX14+1
    CALL NEW(TOP14,(/NATOMS,K/))
    TOP14%I(1:NATOMS,1:NMAX14+1)=0 
!
    DO II=1,NATOMS
      IN12=TOP12%I(II,1)
      DO J=1,IN12
        JJ=TOP12%I(II,J+1)
        JN12=TOP12%I(JJ,1)
          DO K=1,JN12
          KK=TOP12%I(JJ,K+1)
          IF(II/=KK) THEN
            KN12=TOP12%I(KK,1)
              DO L=1,KN12
              LL=TOP12%I(KK,L+1)
          IF(JJ/=LL.AND.II/=LL) THEN
!
! continue here
!
      NI=TOP14%I(II,1)
!
! check matrix size, increase size if necessary
      IF(NI>NMAX14) THEN
        NMAX14=NMAX14+10
        CALL NEW(TOP14_2,(/NATOMS,NMAX14+1/))
        TOP14_2%I(1:NATOMS,1:NMAX14+1)=0 
        TOP14_2%I(1:NATOMS,1:NMAX14+1-10)=TOP14%I(1:NATOMS,1:NMAX14+1-10)
        CALL DELETE(TOP14)
        CALL NEW(TOP14,(/NATOMS,NMAX14+1/))
        TOP14%I(1:NATOMS,1:NMAX14+1)=TOP14_2%I(1:NATOMS,1:NMAX14+1)
        CALL DELETE(TOP14_2)
      ENDIF
!
      IF(NI/=0) THEN
      DO M=1,NI         
        IF(TOP14%I(II,M+1)==LL) THEN
          EXIT
        ELSE
          TOP14%I(II,1)=NI+1
          TOP14%I(II,1+(NI+1))=LL
          EXIT
        ENDIF
      ENDDO 
      ELSE
          TOP14%I(II,1)=NI+1
          TOP14%I(II,1+(NI+1))=LL
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
      CALL Put(NMAX14,'NMAX14')
      CALL Put(TOP14,'TOP14')
      CALL CloseHDF()
    ENDIF
!
    IF(.NOT.PRESENT(TOP12)) CALL DELETE(TOP12)
    IF(.NOT.PRESENT(TOP14)) CALL DELETE(TOP14)
!
END SUBROUTINE TOPOLOGY_14 
!--------------------------------------------------------------
!
!--------------------------------------------------------------
!
!--------------------------------------------------------------
SUBROUTINE SORT_INTO_BOX(BOXSIZE,C,NATOMS)
!
! sort the atoms of a molecule into boxes
IMPLICIT NONE
REAL(DOUBLE) :: BOXSIZE,VBIG,C(1:3,NATOMS),BXMIN,BXMAX,BYMIN,BYMAX,BZMIN,BZMAX
INTEGER,ALLOCATABLE,DIMENSION(:) :: BOXI,BOXJ
INTEGER,ALLOCATABLE,DIMENSION(:) :: ISIGN 
INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: BOXCOUNTER 
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
  ALLOCATE(BOXI(1:NBOX+1))
  ALLOCATE(BOXJ(1:NATOMS))
  ALLOCATE(ISIGN(1:NATOMS))
!
  ALLOCATE(BOXCOUNTER(1:NX,1:NY,1:NZ))
  BOXCOUNTER(1:NX,1:NY,1:NZ)=0
!
! COUNT NUMBER OF ATOMS IN THE BOX
!
  DO I=1,NATOMS
!
! identify box
!     
    IX=INT((C(1,I)-BXMIN)/BOXSIZE)+1
    IY=INT((C(1,I)-BYMIN)/BOXSIZE)+1
    IZ=INT((C(1,I)-BZMIN)/BOXSIZE)+1
!     
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY
    BOXCOUNTER(IX,IY,IZ)=BOXCOUNTER(IX,IY,IZ)+1
    ISIGN(I)=IORD*(NATOMS+1)+BOXCOUNTER(IX,IY,IZ) !!! shows both box and index within box
!     
  ENDDO
!     
! Count pointers to individual boxes within array A
!     
  DO IZ=1,NZ
  DO IX=1,NX
  DO IY=1,NY
!
    IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY
    BOXI(IORD)=BOXI(IORD)+BOXCOUNTER(IX,IY,IZ)+1
!     
  ENDDO
  ENDDO
  ENDDO
!
! Set up content of boxes as represented in A, sparse row-wise
!
  DO I=1,NATOMS
!
    IORD=INT(ISIGN(I)/(NATOMS+1))
    IADD=ISIGN(I)-IORD*(NATOMS+1)
    BOXJ(BOXI(IORD)-1+IADD)=I
!     
  ENDDO
!
! Now check ordering algorithm
!
  DO I=1,NBOX  
  DO J=BOXI(I),BOXI(I+1)-1
    JJ=BOXJ(BOXI(J)) 
    write(*,100) i,C(1,JJ)
  ENDDO
  ENDDO
100 format(I8,3F12.8)
!
! DEALLOCATE(BOXI)
! DEALLOCATE(BOXJ)
  DEALLOCATE(BOXCOUNTER)
!
END SUBROUTINE SORT_INTO_BOX
!
!----------------------------------------------------------------
!
SUBROUTINE TOPOLOGIES_MM(NATOMS,NBONDS,BONDI,BONDJ,InfFile,TOP12)
! Set up a table which shows the atom numbers of atoms 
! connected to a certain atom by the input bonds (topology mtr)
! Here, the generation of the topology mtr is based on input list
! of bonds
IMPLICIT NONE
TYPE(INT_RNK2),OPTIONAL :: TOP12
TYPE(INT_RNK2) :: TOP13,TOP14
INTEGER :: I,J,K,L,N,M,II,JJ,NI,NJ,NBONDS,NATOMS
INTEGER,DIMENSION(1:NBONDS) :: BONDI,BONDJ
CHARACTER(LEN=DEFAULT_CHR_LEN),OPTIONAL :: InfFile
!
   CALL TOPOLOGY_12(NATOMS,NBONDS,BONDI,BONDJ,TOP12,InfFile)
   CALL TOPOLOGY_13(NATOMS,TOP12,TOP13,InfFile)
   CALL TOPOLOGY_14(NATOMS,TOP12,TOP14,InfFile)
!
END SUBROUTINE TOPOLOGIES_MM
!----------------------------------------------------------------


END MODULE IntCoo
