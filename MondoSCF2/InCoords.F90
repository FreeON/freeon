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
   USE CholFactor
   USE ControlStructures
!
IMPLICIT NONE
!
CONTAINS
!--------------------------------------------------------------
!
   SUBROUTINE Topology_12(NatmsLoc,NBond,BondI,BondJ,Top12,Tag)
     !
     ! Set up a table which shows the atom numbers of Atoms 
     ! connected to a certain atom by the input bonds (Topology mtr)
     ! Here, the generation of the Topology mtr is based on input list
     ! of bonds
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)             :: Top12
     TYPE(INT_RNK2)             :: Top12_2
     INTEGER                    :: I,J,K,L,N,M,II,JJ,NI,NJ
     INTEGER                    :: NBond,NatmsLoc,NMax12
     INTEGER,DIMENSION(1:NBond) :: BondI,BondJ
     CHARACTER(LEN=*)           :: Tag
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
       !
       IF(NI>=NMax12 .OR. NJ>=NMax12) THEN
         NMax12=NMax12+5
         CALL New(Top12_2,(/NatmsLoc,NMax12+1/))
         Top12_2%I(1:NatmsLoc,1:NMax12+1)=0 
         Top12_2%I(1:NatmsLoc,1:NMax12+1-5)=&
                            Top12%I(1:NatmsLoc,1:NMax12+1-5)
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
     CALL WriteINT_RNK2(Top12%I,TRIM(Tag)//'Top12')
   END SUBROUTINE Topology_12 
!
!--------------------------------------------------------------
!
   SUBROUTINE Topology_13(NatmsLoc,Top12,Top13,Tag)
     ! Set up a table which shows the atom numbers of Atoms 
     ! being second neighbours of a certain atom.
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)  :: Top12
     TYPE(INT_RNK2)  :: Top13
     TYPE(INT_RNK2)  :: Top13_2
     INTEGER         :: I,J,K,L,N,M,II,JJ,NI,NJ,NatmsLoc
     INTEGER         :: NMax13,NMax12,KK,IN12,JN12
     CHARACTER(LEN=*):: Tag
     !
     NMax12=Size(Top12%I,2)-1
     !
     NMax13=30
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
             !   check matrix Size, increase Size if necessary
             !    
             IF(NI>=NMax13) THEN
               NMax13=NMax13+30
               CALL New(Top13_2,(/NatmsLoc,NMax13+1/))
               Top13_2%I(1:NatmsLoc,1:NMax13+1)=0 
               Top13_2%I(1:NatmsLoc,1:NMax13+1-30)=&
                                       Top13%I(1:NatmsLoc,1:NMax13+1-30)
               CALL Delete(Top13)
               CALL New(Top13,(/NatmsLoc,NMax13+1/))
               Top13%I(1:NatmsLoc,1:NMax13+1)=&
                                        Top13_2%I(1:NatmsLoc,1:NMax13+1)
               CALL Delete(Top13_2)
             ENDIF
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
           ENDIF !!! II/=KK
         ENDDO !!!! KK
       ENDDO !!!! JJ
     ENDDO !!!! II
     !
     CALL WriteINT_RNK2(Top13%I,TRIM(Tag)//'Top13')
   END SUBROUTINE Topology_13 
!
!--------------------------------------------------------------
!
   SUBROUTINE Topology_14(NatmsLoc,Top12,Top14,Tag)
     ! Set up a table which shows the atom numbers of Atoms 
     ! being second neighbours of a certain atom.
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)         :: Top12
     TYPE(INT_RNK2)         :: Top14
     TYPE(INT_RNK2)         :: Top14_2
     INTEGER                :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
     INTEGER                :: NatmsLoc,NMax14,NMax12,IN12,JN12,KN12
     CHARACTER(LEN=*)       :: Tag
     !
     NMax12=Size(Top12%I,2)-1
     !
     NMax14=30
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
                   NMax14=NMax14+30
                   CALL New(Top14_2,(/NatmsLoc,NMax14+1/))
                   Top14_2%I(1:NatmsLoc,1:NMax14+1)=0 
                   Top14_2%I(1:NatmsLoc,1:NMax14+1-30)=&
                                 Top14%I(1:NatmsLoc,1:NMax14+1-30)
                   CALL Delete(Top14)
                   CALL New(Top14,(/NatmsLoc,NMax14+1/))
                   Top14%I(1:NatmsLoc,1:NMax14+1)=&
                                 Top14_2%I(1:NatmsLoc,1:NMax14+1)
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
     CALL WriteINT_RNK2(Top14%I,TRIM(Tag)//'Top14')
   END SUBROUTINE Topology_14 
!
!--------------------------------------------------------------
!
   SUBROUTINE SORT_INTO_Box1(BoxSize,XYZ,NatmsLoc, &
                       NX,NY,NZ,BXMIN,BYMIN,BZMIN)
     !
     ! sort the Atoms of a molecule into Boxes
     !
     IMPLICIT NONE
     INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NatmsLoc
     REAL(DOUBLE) :: BoxSize,VBIG,XYZ(1:3,NatmsLoc)
     REAL(DOUBLE) :: BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
     !
     CALL BoxBorders(XYZ,BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax)
     !
     NX=INT((BXMax-BXMIN)/BoxSize)+1
     NY=INT((BYMax-BYMIN)/BoxSize)+1
     NZ=INT((BZMax-BZMIN)/BoxSize)+1
     NBox=NX*NY*NZ
     !
   END SUBROUTINE SORT_INTO_Box1
!
!--------------------------------------------------------------
!
   SUBROUTINE SORT_INTO_Box2(BoxSize,XYZ,NatmsLoc, &
     NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,ISet)
     !
     ! BoxI(I) : contains the ordering (box) number of the 
     !           first atom of the I-th Box (like in sparse row-wise)
     ! BoxJ(J) : gives the original serial number of the atom desribed 
     !           by the J-th ordering number
     ! XYZ: contains Cartesian coordinates of Atoms
     ! BoxSize: linear Box Size
     !
     IMPLICIT NONE
     INTEGER,OPTIONAL :: ISET 
     TYPE(INT_VECT)   :: BoxI,BoxJ
     TYPE(INT_RNK2)   :: ISign
     TYPE(INT_RNK3)   :: BoxCounter 
     REAL(DOUBLE)     :: BoxSize,VBIG
     REAL(DOUBLE)     :: BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
     INTEGER          :: IORD,IADD,NatmsLoc
     INTEGER          :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     TYPE(INT_VECT),OPTIONAL      :: BoxI1,BoxJ1
     !
     VBIG=1.D+90 
     NBox=NX*NY*NZ
     CALL New(BoxI,NBox+1)
     CALL New(BoxJ,NatmsLoc)
     CALL New(ISign,(/2,NatmsLoc/))
     CALL New(BoxCounter,(/NX,NY,NZ/))
     BoxCounter%I(1:NX,1:NY,1:NZ)=0
     BoxI%I(1:NBox+1)=0
     !
     ! Count number of atoms in the box
     !
     DO I=1,NatmsLoc
       IX=INT((XYZ(1,I)-BXMIN)/BoxSize)+1
       IY=INT((XYZ(2,I)-BYMIN)/BoxSize)+1
       IZ=INT((XYZ(3,I)-BZMIN)/BoxSize)+1
       IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY !order parameter: ZXY
       BoxCounter%I(IX,IY,IZ)=BoxCounter%I(IX,IY,IZ)+1
       ISign%I(1,I)=IORD
       ISign%I(2,I)=BoxCounter%I(IX,IY,IZ) 
     ENDDO
     !     
     ! Count pointers to individual Boxes within array A
     !     
     BoxI%I(1)=1
     DO IZ=1,NZ
     DO IX=1,NX
     DO IY=1,NY
       IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY
       BoxI%I(IORD+1)=BoxI%I(IORD)+BoxCounter%I(IX,IY,IZ)
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
     IF(PRESENT(BoxI1)) THEN
       BoxI1%I(:)=BoxI%I(:)
     ENDIF
     IF(PRESENT(BoxJ1)) THEN
       BoxJ1%I(:)=BoxJ%I(:)
     ENDIF
     !
     CALL Delete(ISign)
     CALL Delete(BoxCounter)
     CALL Delete(BoxI)
     CALL Delete(BoxJ)
   END SUBROUTINE SORT_INTO_Box2
!
!----------------------------------------------------------------
!
   SUBROUTINE Topologies_MM(NatmsLoc,NBond,BondI,BondJ,SCR)
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
     TYPE(INT_RNK2)                 :: Top12,Top13,Top14
     TYPE(INT_RNK2)                 :: Top_Excl,Top_Excl14
     INTEGER                        :: NBond,NatmsLoc
     INTEGER,DIMENSION(1:NBond)     :: BondI,BondJ
     CHARACTER(LEN=DCL)             :: PathName
     CHARACTER(LEN=*)               :: SCR
     !
     PathName=TRIM(SCR)//'MM'
     CALL Topology_12(NatmsLoc,NBond,BondI,BondJ,Top12,PathName)
     CALL Topology_13(NatmsLoc,Top12,Top13,PathName)
     CALL Topology_14(NatmsLoc,Top12,Top14,PathName)
     CALL Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl,PathName)
     CALL Delete(Top_Excl)
     CALL Excl_List14(NatmsLoc,Top12,Top13,Top14,Top_Excl14,PathName)
     CALL Delete(Top_Excl14)
     CALL Delete(Top12)
     CALL Delete(Top13)
     CALL Delete(Top14)
   END SUBROUTINE Topologies_MM
!
!----------------------------------------------------------------
!
   SUBROUTINE Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl,Tag)
     !
     ! This subroutine merges Topological information
     ! to get the list for Exclusion energy calculation
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)          :: Top_Excl_Out,Top12,Top13,Top14
     TYPE(INT_RNK2)          :: Top_Excl,Top_New
     CHARACTER(LEN=*)        :: Tag
     INTEGER                 :: NatmsLoc,NMax12,NMax13,NMax14
     INTEGER                 :: NMax_Excl,NNew,NOLD
     INTEGER                 :: NMax_Excl_Out,I,J,K,KK,JJ
     !
     NMax12=Size(Top12%I,2)-1
     NMax13=Size(Top13%I,2)-1
     NMax14=Size(Top14%I,2)-1
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
     CALL WriteINT_RNK2(Top_Excl%I,TRIM(Tag)//'TOP_Excl')
   END SUBROUTINE Excl_LIST 
!
!----------------------------------------------------------------
!
   SUBROUTINE Excl_List14(NatmsLoc,Top12,Top13,Top14,Top_Excl,Tag)
     !
     ! This subroutine merges Topological information
     ! to get the list for Exclusion energy calculation
     ! of Atoms in 14 distance. From the Top14 list
     ! Top13 and Top12 occurences must be filtered out
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)          :: Top_Excl_Out,Top12,Top13,Top14
     TYPE(INT_RNK2)          :: Top_Excl,Top_New
     CHARACTER(LEN=*)        :: Tag
     INTEGER                 :: NatmsLoc,NMax12,NMax13,NMax14
     INTEGER                 :: NMax_Excl,NNew,NOLD
     INTEGER                 :: NMax_Excl_Out,I,J,K,KK,JJ,NExcl,N12,N13
     !
     NMax12=Size(Top12%I,2)-1
     NMax13=Size(Top13%I,2)-1
     NMax14=Size(Top14%I,2)-1
     !
     ! Initialize Top_Excl to the Size of Top14 and to zero 
     !
     NMax_Excl=NMax14
     CALL New(Top_Excl,(/NatmsLoc,NMax_Excl+1/))
     Top_Excl%I=0
     !
     ! Now merge Topologies, in order to avoid double counting in 
     ! Exclusion energies
     !
     NMax_Excl_Out=0
     DO I=1,NatmsLoc
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
     ENDDO 
     CALL WriteINT_RNK2(Top_Excl%I,TRIM(Tag)//'TOP_Excl14')
   END SUBROUTINE Excl_List14 
!
!----------------------------------------------------------------
!
   SUBROUTINE BMatrix(XYZ,NIntC,IntCs,B,LinCrit,DoClssTrf,ThreeAt)
     !
     ! This subroutine calculates the sparse, TYPE(BMATR) type 
     ! representation of the B matrix 
     ! In current version Cartesian coordinates 
     ! must be passed in in AU
     ! Linear bendings _must_ always appear in pairs, 
     ! Defd as LINB1 and LINB2
     !
     IMPLICIT NONE
     INTEGER :: NIntC,NIntC2,I,J,K,L,NatmsLoc,IntCoo,ThreeAt(1:3)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INTC) :: IntCs
     TYPE(BMATR):: B
     LOGICAL       :: DoClssTrf
     REAL(DOUBLE)  :: LinCrit
     !
     NatmsLoc=SIZE(XYZ,2)
     CALL INTCValue(IntCs,XYZ,LinCrit)
     !
     ! allocate B matrix
     !
     CALL NEW(B,NIntC)
     B%IB(:,:)=IntCs%Atoms(:,:)
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
         CALL STRE(XYZ(1:3,I),XYZ(1:3,J),B%B(IntCoo,1:12))
         !   write(*,100) IntCs%Def(IntCoo)(1:5), &
         !   B%IB(IntCoo,1:4),B%B(IntCoo,1:12)
         !
       ELSE IF(IntCs%Def(IntCoo)(1:4)=='BEND') THEN
         ! bend
         I=IntCs%Atoms(IntCoo,1)
         J=IntCs%Atoms(IntCoo,2) !!! central atom
         K=IntCs%Atoms(IntCoo,3) 
         CALL BEND(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),B%B(IntCoo,1:12))
         !         write(*,100) IntCs%Def(IntCoo)(1:5), &
         !           B%IB(IntCoo,1:4),B%B(IntCoo,1:12)
       ELSE IF(IntCs%Def(IntCoo)(1:4)=='OUTP') THEN
         ! out of plane
         I=IntCs%Atoms(IntCoo,1) !!! end atom
         J=IntCs%Atoms(IntCoo,2) !!! central atom
         K=IntCs%Atoms(IntCoo,3) !!! Def plane
         L=IntCs%Atoms(IntCoo,4) !!! Def plane
         CALL OutP(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ(1:3,L),&
           B%B(IntCoo,1:12))
         !write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:4), &
         ! B%B(IntCoo,1:12)*AngstromsToAu
       ELSE IF(IntCs%Def(IntCoo)(1:4)=='TORS') THEN
         ! torsion of i-j-k-l
         I=IntCs%Atoms(IntCoo,1)
         J=IntCs%Atoms(IntCoo,2) 
         K=IntCs%Atoms(IntCoo,3) 
         L=IntCs%Atoms(IntCoo,4) 
         CALL TORS(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
           XYZ(1:3,L),B%B(IntCoo,1:12),IntCs%Active(IntCoo))
         !write(*,*) IntCoo
         !write(*,100) IntCs%Def(IntCoo)(1:5),B%IB(IntCoo,1:4), &
         !B%B(IntCoo,1:12)
       ELSE IF(IntCs%Def(IntCoo)(1:5)=='LINB1') THEN
         ! linear bendig of i-j-k
         IF(IntCs%Def(IntCoo+1)(1:5)/='LINB2') &
            CALL Halt('LINB Definitions are not paired!')
         I=IntCs%Atoms(IntCoo,1)
         J=IntCs%Atoms(IntCoo,2) 
         K=IntCs%Atoms(IntCoo,3) 
         CALL LinB2(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),&
                    BB1=B%B(IntCoo,1:12),BB2=B%B(IntCoo+1,1:12))
        !CALL LINB(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
        !  B%B(IntCoo,1:12),B%B(IntCoo+1,1:12))
         !write(*,100) IntCs%Def(IntCoo)(1:5), &
         ! B%IB(IntCoo,1:4),B%B(IntCoo,1:12)
         !write(*,100) IntCs%Def(IntCoo)(1:5), &
         ! B%IB(IntCoo+1,1:4),B%B(IntCoo+1,1:12)
       ELSE IF(IntCs%Def(IntCoo)(1:5)=='LINB2') THEN
         CYCLE
       ELSE IF(IntCs%Def(IntCoo)(1:5)=='CARTX' ) THEN
         I=IntCs%Atoms(IntCoo,1)
         CALL BCART('X',B%B(IntCoo,1:12),IntCs%Constraint(I))
       ELSE IF(IntCs%Def(IntCoo)(1:5)=='CARTY' ) THEN
         I=IntCs%Atoms(IntCoo,1)
         CALL BCART('Y',B%B(IntCoo,1:12),IntCs%Constraint(I))
       ELSE IF(IntCs%Def(IntCoo)(1:5)=='CARTZ' ) THEN
         I=IntCs%Atoms(IntCoo,1)
         CALL BCART('Z',B%B(IntCoo,1:12),IntCs%Constraint(I))
       ENDIF
       !
     ENDDO !!!! loop over internal coords
     100 format(A5,4I6,/,6F12.6,/,6F12.6)
     !
     ! Clean columns which cause redundancy or are related to Cartesian
     ! constraints
     !
     CALL CleanBConstr(IntCs,B,NatmsLoc)
     IF(.NOT.DoClssTrf) THEN
       CALL CleanB(ThreeAt,B,NIntC)
     ENDIF
   END SUBROUTINE BMatrix
!
!-----------------------------------------------------------------------
!
   SUBROUTINE Stre(XI,XJ,BB)
     !
     !  Calculate the B matrix elements for a bond stretch as Defined
     !  by Wilson.
     !
     IMPLICIT NONE
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,RIJ
     REAL(DOUBLE) :: dijsq,t,dij
     INTEGER :: I,J,K,L,M,N
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
       T=T/Dij
       BB(M)=-T
       BB(3+M)=T
     ENDDO 
   END SUBROUTINE STRE
!
!----------------------------------------------------------------
!
   SUBROUTINE Bend(XI,XJ,XK,BB)
     !
     !  This subroutine computes the b matrix elements of a valence
     !  angle bending coordinate as defined by Wilson.
     !  i and k are the numbers of the end Atoms.
     !  j is the number of the central atom.
     !  noat: number of Atoms
     !  ic: serial number of internal coordinate
     !
     IMPLICIT NONE
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,RJI,RJK,EJI,EJK
     REAL(DOUBLE) :: DIJSq,Tp,T,DJKSq,DxSq,Dx,DOtj,DJK,DJI
     REAL(DOUBLE) :: SinJ,SMI,SMK,SUM
     INTEGER  :: I,J,K,L,M,N,II,JJ,KK,IER
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
     IF(ABS(DotJ) > 0.99995D0) RETURN    !!!! colinearity
     SinJ=SQRT(One-DotJ*DotJ)
     DO 40 M=1,3
       SMI=Dx*(DotJ*Eji(M)-Ejk(M))/(Dji*SinJ)
       BB(M) = SMI
       SMK=DX*(DotJ*Ejk(M)-Eji(M))/(DJK*SinJ)
       BB(6+M) = SMK
       Sum=SMI+SMK
       BB(3+M) = -Sum
     40 CONTINUE
     RETURN
     50 IER=1
     RETURN
     60 IER=-1
     !WRITE(6,1000)
       BB(:)=Zero 
     RETURN
     !
     1000 format(' i-j-k is colinear - use linear bend!'//&
       'For this bending the row of the B matrix will be deleted!')
     !
   END SUBROUTINE BEND
!
!---------------------------------------------------------------------
!
   SUBROUTINE OutP(XI,XJ,XK,XL,BB)
     !
     !  i is the end atom (atom wagged with respect to j-k-l plane).
     !  j is the apex atom (Atoms i, k and l are attached to j).
     !  k and l are the anchor Atoms (Define the j-k-l plane).
     !
     IMPLICIT NONE
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,XL,RJI,RJK,RJL
     REAL(DOUBLE),DIMENSION(1:3) :: EJI,EJK,EJL,C1
     REAL(DOUBLE) :: DIJSq,Tp,t,DJKSq,djlsq,dx,DOtj
     REAL(DOUBLE) :: DJK,DJI,djl,cosi,cosk,cosl,Tpp
     REAL(DOUBLE) :: sinsin,sini,dot,sint,cost,tant,cossin,sml
     REAL(DOUBLE) :: sinj,SMI,SMK,SUM
     INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
     REAL(DOUBLE) :: BB(1:12)
     !
     Dx=One 
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
       T=XL(M)-XJ(M)
       Rjl(M)=T
       DjlSq=DjlSq+T*T
     20 CONTINUE      
     Dji=SQRT(DijSq)
     Djk=SQRT(DjkSq)
     Djl=SQRT(DjlSq)
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
     IF(ABS(CosI) > 0.99995D0) RETURN  !!!! GO TO 70
     SinSin=One-CosI*CosI
     SinI=SQRT(SinSin)
     C1(1)=Ejk(2)*Ejl(3)-Ejk(3)*Ejl(2)
     C1(2)=Ejk(3)*Ejl(1)-Ejk(1)*Ejl(3)
     C1(3)=Ejk(1)*Ejl(2)-Ejk(2)*Ejl(1)
     Dot=Eji(1)*C1(1)+Eji(2)*C1(2)+Eji(3)*C1(3)
     SinT=Dot/SinI
     !     if(abs(sint).gt.0.00005d0) then
     !        write(6,1020)nob
     !        write(6,'(15x,a36)')intch(nob)
     !     end if
     IF(ABS(SinT)>0.99995D0) RETURN !!!! GO TO 80 !!!! sint shows OutP angle
     CosT=SQRT(One-SinT*SinT)
     TanT=SinT/CosT
     CosSin=CosT*SinI
     DO 50 M=1,3
       T=C1(M)/CosSin
       SMI=(T-TanT*Eji(M))/Dji
       BB(M) = Dx*SMI
       SMK=T*(CosI*CosK-CosL)/(SinSin*Djk)
       BB(6+M) = Dx*SMK
       SML=T*(CosI*CosL-CosK)/(SinSin*Djl)
       BB(9+M) = Dx*SML
       Sum=SMI+SMK+SML
       BB(3+M) = -Dx*Sum
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
   SUBROUTINE Tors(XI,XJ,XK,XL,BB,Active)
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
     REAL(DOUBLE) :: SMI,SMK,dij,dijsq,dx,DJK,t,dxsq,DJKSq
     REAL(DOUBLE) :: cosj,sin2j,smj,dlksq,dlk,cosk,sin2k,sml 
     INTEGER :: I,J,K,L,M,N,II,JJ,KK,LL,IER
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
     IF(ABS(CosJ)>0.99995d0) THEN   !!!! GO TO 120
       BB(:)=Zero
       Active=.FALSE.
       RETURN
     ENDIF
     Sin2J=(One-CosJ*CosJ)
     Cr(1)=Eij(2)*Ejk(3)-Eij(3)*Ejk(2)
     Cr(2)=Eij(3)*Ejk(1)-Eij(1)*Ejk(3)
     Cr(3)=Eij(1)*Ejk(2)-Eij(2)*Ejk(1)
     DO 60 M=1,3
       T=Cr(M)/Sin2J
       SMI=T*Dij
       BB(M) = -Dx*SMI
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
     IF(ABS(CosK)>0.99995d0) THEN !!!! GO TO 120
       BB(:)=Zero
       Active=.FALSE.
       RETURN
     ENDIF
     Sin2K=(One-CosK*CosK)
     Cr(1)=Elk(3)*Ejk(2)-Elk(2)*Ejk(3)
     Cr(2)=Elk(1)*Ejk(3)-Elk(3)*Ejk(1)
     Cr(3)=Elk(2)*Ejk(1)-Elk(1)*Ejk(2)
     DO 90 M=1,3
       T=Cr(M)/Sin2K
       SML=T*Dlk
         BB(9+M) = -Dx*SML
       SMJ=T*CosK*DJK
       SJ(M)=SJ(M)+SMJ
       SMK=SML-SMJ
       SK(M)=SK(M)+SMK
     90 CONTINUE
     DO 100 M=1,3
       SMJ=SJ(M)
       BB(3+M) = Dx*SMJ
       SMK=SK(M)
       BB(6+M) = Dx*SMK
     100 CONTINUE
     RETURN
     110 CONTINUE
     IEr=1
     RETURN
     120 CONTINUE
     IEr=-1
     !WRITE(6,1030)
     ! colinear bonds in torsion, delete the 
     !corresponding row of the B matrix!
     BB(:)=Zero
     Active=.FALSE.
     RETURN
     !
     1030 format(' i-j-k or j-k-l is colinear (no torsion'//&
     'possibble), this row of the B matrix will be deleted.')
     !
   END SUBROUTINE TORS
!
!---------------------------------------------------------------------
!
   SUBROUTINE LinB2(XI,XJ,XK,BB1,BB2,Value1,Value2)
     !
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK
     REAL(DOUBLE),DIMENSION(1:12),OPTIONAL:: BB1,BB2
     INTEGER                     :: I,J,K,L,M,N
     REAL(DOUBLE)                :: Sum,ValChk,Prefact
     REAL(DOUBLE),OPTIONAL       :: Value1,Value2
     TYPE(DBL_VECT)              :: VectIJ,VectJK,VectA,VectAA,A,AA 
     LOGICAL                     :: Active
     !
     CALL New(VectIJ,3)
     CALL New(VectJK,3)
     CALL New(VectA ,3)
     CALL New(VectAA,3)
     CALL New(A,3)
     CALL New(AA,3)
     !
     ! Set up coordinates of auxiliary point a
     ! possibly perpendicular to I-J-K plane
     !
     VectIJ%D=XI-XJ
     VectJK%D=XJ-XK
     CALL CROSS_PRODUCT(VectIJ%D,VectJK%D,VectA%D)
     Sum=DOT_PRODUCT(VectA%D,VectA%D)
     IF(Sum<1.D-4) THEN
       !set up auxiliary vector not parallel with I-J-K line
       Sum=DOT_PRODUCT(VectIJ%D,VectIJ%D)+DOT_PRODUCT(VectJK%D,VectJK%D)
       VectIJ%D(1)=VectJK%D(1)+Sum
       VectIJ%D(2)=VectJK%D(2)+Sum
       VectIJ%D(3)=VectJK%D(3)-Sum
       CALL CROSS_PRODUCT(VectIJ%D,VectJK%D,VectA%D)
       Sum=DOT_PRODUCT(VectA%D,VectA%D)
       Sum=One/SQRT(Sum)
       VectA%D=Sum*VectA%D
       !
       CALL CROSS_PRODUCT(VectIJ%D,VectA%D,VectAA%D)
       Sum=DOT_PRODUCT(VectAA%D,VectAA%D)
       Sum=One/SQRT(Sum)
       VectAA%D=Sum*VectAA%D
     ELSE 
       VectAA%D=XI+XK-Two*XJ
       Sum=DOT_PRODUCT(VectAA%D,VectAA%D)
       Sum=One/SQRT(Sum)
       VectAA%D=Sum*VectAA%D
     ENDIF
     !
     ! generate point a, and the corresponding matrix element of B
     !
     A%D=XJ+VectA%D
     CALL TORSValue(XI,A%D,XJ,XK,0.D0,ValChk,Active)
     Prefact=SIGN(One,ValChk)
     IF(PRESENT(BB1)) THEN
       CALL TORS(XI,A%D,XJ,XK,BB1,Active)
       BB1(4:6)=BB1(7:9)
       BB1(7:9)=BB1(10:12)
       BB1(10:12)=Zero
       BB1=Prefact*BB1
     ENDIF
     IF(PRESENT(Value1)) THEN
       Value1=Prefact*ValChk
     ENDIF
     !
     ! generate point aa, and the corresponding matrix element of B
     !
     AA%D=XJ+VectAA%D
     CALL TORSValue(XI,AA%D,XJ,XK,0.D0,ValChk,Active)
     Prefact=SIGN(One,ValChk)
     IF(PRESENT(BB2)) THEN
       CALL TORS(XI,AA%D,XJ,XK,BB2,Active)
       BB2(4:6)=BB2(7:9)
       BB2(7:9)=BB2(10:12)
       BB2(10:12)=Zero
       BB2=Prefact*BB2
     ENDIF
     IF(PRESENT(Value2)) THEN
       Value2=ValChk
     ENDIF
     !
     CALL Delete(VectIJ)
     CALL Delete(VectJK)
     CALL Delete(VectA)
     CALL Delete(VectAA)
     CALL Delete(A)
     CALL Delete(AA)
   END SUBROUTINE LinB2
!
!----------------------------------------------------------------
!
   SUBROUTINE BCART(CHAR,B,Constraint)
     INTEGER :: I,J
     REAL(DOUBLE),DIMENSION(1:12) :: B
     CHARACTER :: CHAR
     LOGICAL   :: Constraint
     !
     B=Zero
     IF(CHAR=='X') THEN
       IF(Constraint) THEN
         B(1)=Zero 
       ELSE
         B(1)=One
       ENDIF
     ENDIF
     !
     IF(CHAR=='Y') THEN
       IF(Constraint) THEN
         B(2)=Zero 
       ELSE
         B(2)=One  
       ENDIF
     ENDIF
     !
     IF(CHAR=='Z') THEN
       IF(Constraint) THEN
         B(3)=Zero 
       ELSE
         B(3)=One 
       ENDIF
     ENDIF
   END SUBROUTINE BCART
!
!----------------------------------------------------------------
!
   SUBROUTINE DefineIntCoos(NatmsLoc,XYZ,AtNum,IntSet,Refresh, &
                            IntCs,NIntC,CtrlCoord,SCRPath)
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
     INTEGER :: NMax_Excl,NMax12,ILast,NOutP
     REAL(DOUBLE),DIMENSION(1:3,1:NatmsLoc) :: XYZ
     TYPE(INT_RNK2) :: BondIJ 
     TYPE(INT_RNK2) :: AngleIJK,TorsionIJKL,OutPIJKL
     INTEGER        :: NX,NY,NZ,NBOX,NAngle,NTorsion,NLinb,Refresh
     REAL(DOUBLE)   :: BXMIN,BYMIN,BZMIN,BoxSize,Fact,Value
     TYPE(DBL_VECT) :: CritRad,BondLength   
     TYPE(INT_RNK2) :: Top12,Top13,Top14,Top_Excl,Top_VDW
     TYPE(INTC)                     :: IntCs
     TYPE(INT_VECT)                 :: BoxI,BoxJ
     TYPE(INT_VECT)                 :: HBondMark,HAngleMark,HTorsMark
     INTEGER,DIMENSION(:)           :: AtNum
     TYPE(CoordCtrl)                :: CtrlCoord
     CHARACTER(LEN=*)               :: SCRPath
     CHARACTER(LEN=1)               :: LongRangeTag
     LOGICAL                        :: HBondOnly
     !
     NIntC=0
     NBond=0
     NAngle=0
     NTorsion=0
     NOutP=0
     !
     ! first sort atoms into boxes      
     !
     BoxSize=2.5D0*AngstromsToAU !in A
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
     !now define bonding scheme, based on Slater or Van der Waals radii
     !
     IF(IntSet==1) THEN
       N=SIZE(SLRadii,1)
       CALL New(CritRad,N)
       Fact=1.3D0 !!! Scaling factor for Slater Radii
       CritRad%D=Fact*SLRadii*AngstromsToAU
       HBondOnly=.FALSE.
     ELSE
       N=SIZE(VDWRadii,1)
       CALL New(CritRad,N)
       Fact=CtrlCoord%VDWFact !!! Scaling factor for VDW Radii
       CritRad%D=Fact*VDWRadii*AngstromsToAU
      !HBondOnly=.TRUE.
       HBondOnly=.FALSE.
     ENDIF
     !
     CALL BondList(NatmsLoc,XYZ,NBond,AtNum,SCRPath,IntSet, &
                   BoxI,BoxJ,NBox,NX,NY,NZ,CritRad,HBondOnly, &
                   HBondMark,BondLength,BondIJ)
     !
     IF(IntSet==1) THEN
       !
       ! Now define covalent topology matrices
       !
       CALL Topology_12(NatmsLoc,NBond,BondIJ%I(1,1:NBond),&
         BondIJ%I(2,1:NBond),Top12,SCRPath)
       CALL Topology_13(NatmsLoc,Top12,Top13,SCRPath)
       CALL Topology_14(NatmsLoc,Top12,Top14,SCRPath)
       CALL Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl,SCRPath)
       CALL Delete(Top_Excl)
       CALL Delete(Top14)
       CALL Delete(Top13)
     ELSE
       !
      !CALL LinearHBond(NatmsLoc,AtNum,NBond, &
      !                 HBondMark,BondIJ,BondLength)
     ENDIF
     !
     ! Now define bond angles and torsions
     !
     IF(IntSet==1) THEN
       !
       CALL AngleList(NatmsLoc,Top12,AngleIJK, &
                      NAngle,BondIJ,NBond,HBondMark,.FALSE.)
       !
       CALL TorsionList(NatmsLoc,Top12,BondIJ, &
                        TorsionIJKL,NTorsion,HBondMark)
     ELSE
       CALL ReadINT_RNK2(Top12,TRIM(SCRPath)//'Top12',I,J)
         NMax12=J-1
       CALL VDWTop(Top12,BondIJ)
       CALL AngleList(NatmsLoc,Top12,AngleIJK, &
                      NAngle,BondIJ,NBond,HBondMark,.TRUE.)
       CALL TorsionList(NatmsLoc,Top12,BondIJ, &
                        TorsionIJKL,NTorsion,HBondMark)
       CALL OutPList(XYZ,Top12,CtrlCoord%OutPCrit,NOutP,OutPIJKL)
     ENDIF
     !
     ! Fill Data into IntCs
     !
     NIntC=NBond+NAngle+NTorsion+NOutP
     IF(NIntC/=0) THEN          
       CALL New(IntCs,NIntC)
       IntCs%Def='     '
       IntCs%FCType=' '
       IntCs%Atoms(:,:)=0   
       IntCs%Value=Zero
       IntCs%Constraint=.FALSE.
       IntCs%ConstrValue=Zero   
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
       ILast=NBond+NAngle+NTorsion
       DO I=1,NOutP
         IntCs%Def(ILast+I)='OUTP '
         IntCs%Atoms(ILast+I,1:4)=OutPIJKL%I(1:4,I)
       ENDDO
!      ! Put FCType
!      IF(IntSet==1) THEN
!        IntCs%FCType(:)='C' ! Covalent bonding scheme
!      ELSE IF(IntSet==2) THEN
!        J=0
!        DO I=1,NIntC
!          IF(IntCs%Def(I)(1:4)=='STRE') THEN
!!!!!!!!!!!! take more care to the marks
!            IF(HBondMark%I(I)==1) THEN
!              IntCs%FCType(I)='H' ! true H-bond  
!            ELSE
!              IntCs%FCType(I)='V' ! Van der Waals
!            ENDIF
!         !ELSE IF(HasAngle(IntCs%Def(I))) THEN
!          ENDIF
!        ENDDO
!      ELSE
!        IntCs%FCType(:)=' '
!      ENDIF
     ENDIF
     !
     ! tidy up
     IF(IntSet==2) CALL Delete(OutPIJKL)
     CALL Delete(HBondMark)
  !  CALL Delete(HAngleMark)
     CALL Delete(TorsionIJKL)
     CALL Delete(AngleIJK)
     CALL Delete(CritRad)
     CALL Delete(BondLength)
     CALL Delete(BondIJ)
     CALL Delete(BoxI)
     CALL Delete(BoxJ)
     !
     ! Check for bending - lin.bending transions
     ! also for long range torsions!
     !
     IF(IntSet==1) THEN
       LongRangeTag='C'
     ELSE
       LongRangeTag='V'
     ENDIF
     IF(Refresh/=5) CALL ChkBendToLinB(IntCs,NIntC,XYZ,Top12, &
                                       CtrlCoord,SCRPath,LongRangeTag)
     CALL Delete(Top12)
     !
   END SUBROUTINE DefineIntCoos
!
!--------------------------------------------------------
!
   SUBROUTINE VDWTop(Top12,BondIJ)
     TYPE(INT_RNK2) :: Top12,BondIJ,Top12New
     TYPE(INT_VECT) :: Addition
     INTEGER        :: NatmsLoc,NBond,I,J,Dim2,AddDim,I1,I2,II1,II2
     !
     NatmsLoc=SIZE(Top12%I,1)
     Dim2=SIZE(Top12%I,2)
     NBond=SIZE(BondIJ%I,2)
     CALL New(Addition,NatmsLoc)
     !
     Addition%I=0
     DO I=1,NBond    
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       Addition%I(I1)=Addition%I(I1)+1
       Addition%I(I2)=Addition%I(I2)+1
     ENDDO
     AddDim=MAXVAL(Addition%I)
     !
     CALL New(Top12New,(/NatmsLoc,Dim2+AddDim/))
     Top12New%I=0
     DO I=1,NatmsLoc
       DO J=1,Dim2 ; Top12New%I(I,J)=Top12%I(I,J) ; ENDDO
     ENDDO
     DO I=1,NBond
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       II1=Top12New%I(I1,1)+1
       Top12New%I(I1,1)=II1
       Top12New%I(I1,II1+1)=I2
       II2=Top12New%I(I2,1)+1
       Top12New%I(I2,1)=II2
       Top12New%I(I2,II2+1)=I1
     ENDDO
     !
     CALL Delete(Top12)
     CALL New(Top12,(/NatmsLoc,Dim2+AddDim/))
     Top12%I=Top12New%I
     !
     CALL Delete(Top12New)
     CALL Delete(Addition)
   END SUBROUTINE VDWTop
!
!--------------------------------------------------------
!
   SUBROUTINE BondList(NatmsLoc,XYZ,NBond,AtNum,SCRPath,IntSet, &
          BoxI,BoxJ,NBox,NX,NY,NZ,CritRad,HBondOnly, &
          HBondMark,BondLength,BondIJ)
     IMPLICIT NONE
     INTEGER                     :: I,J,NatmsLoc,NBond
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(DBL_VECT)              :: CritRad !!! VdW or Slater Radii
     REAL(DOUBLE)                :: R12,R12_2,CritDist,HBondMax
     INTEGER                     :: IZ,IX,IY,I1,I2,JJ1,JJ2
     INTEGER                     :: NBOX,IORD,IORDD
     INTEGER                     :: IZD,IXD,IYD,NX,NY,NZ,NJJ1,NJJ2
     INTEGER                     :: NMax12,JJ,IntSet,IDimExcl,NBondEst
     LOGICAL                     :: HBondOnly
     LOGICAL                     :: FillBondIJ
     TYPE(INT_RNK2)              :: BondIJ,BondIJAux
     TYPE(DBL_VECT)              :: DVect,BondLength,BondLengthAux
     INTEGER,DIMENSION(:)        :: AtNum
     TYPE(INT_VECT)              :: BoxI,BoxJ,HBondMark,HBondMarkAux
     TYPE(INT_RNK2)              :: Top_Excl,TopCharge
     CHARACTER(LEN=*)            :: SCRPath
     LOGICAL                     :: FoundHBond,FoundMetLig
     !     
     NBondEst=NatmsLoc*20 ! estimate
     NBond=NBondEst
     CALL New(BondIJAux,(/2,NBond/))
     CALL New(BondLengthAux,NBond)
     CALL New(HBondMarkAux,NBond)
     HBondMarkAux%I=0
     !
     HBondMax=2.3D0*AngstromsToAu ! in Au 
     !
     IF(IntSet==2) THEN
       CALL ReadINT_RNK2(Top_Excl,TRIM(SCRPath)//'TOP_Excl',I,J)
       IF(I/=NatmsLoc) CALL Halt('Dimension error in BondList.')
       CALL ReadINT_RNK2(TopCharge,TRIM(SCRPath)//'Top12',I,J)
       IF(I/=NatmsLoc) CALL Halt('Dimension error in BondList.')
       !
       NMax12=J-1
       DO I=1,NatmsLoc
         DO J=1,TopCharge%I(I,1)
           JJ=TopCharge%I(I,1+J)
           TopCharge%I(I,1+J)=AtNum(JJ)
         ENDDO
       ENDDO
     ENDIF
     !     
     CALL New(DVect,3)
     !
     !  Go through all boxes and their neighbours
     !
     NBond=0
     DO IZ=1,NZ
       DO IX=1,NX
         DO IY=1,NY
           ! absolute index of a box
           IOrd=NX*NY*(IZ-1)+NY*(IX-1)+IY
           DO I1=BoxI%I(IOrd),BoxI%I(IOrd+1)-1
             JJ1=BoxJ%I(I1) !!! atom in central box
             NJJ1=AtNum(JJ1)
             ! second atom may come from central or neigbouring Boxes
             !and must be an MM atom,LJ is not calculated for QM-QMpairs
             DO IZD=-1,1
               IF(IZ+IZD>0 .AND. IZ+IZD<=NZ) THEN
                 DO IXD=-1,1
                   IF(IX+IXD>0 .AND. IX+IXD<=NX) THEN
                     DO IYD=-1,1
                       IF(IY+IYD>0 .AND. IY+IYD<=NY) THEN
                         IOrdD=NX*NY*(IZ-1+IZD)+NY*(IX-1+IXD)+IY+IYD
                         DO I2=BoxI%I(IOrdD),BoxI%I(IOrdD+1)-1
                           JJ2=BoxJ%I(I2) !!! second atom
                           NJJ2=AtNum(JJ2)
                           FoundHBond=.FALSE.
                           FoundMetLig=.FALSE.
                           IF(JJ2<=JJ1) CYCLE !!!avoid double counting
                           IF(HasMetal(NJJ1).OR.HasMetal(NJJ2)) THEN
                             IF(HasLigand(NJJ1).OR.HasLigand(NJJ2)) THEN
                               FoundMetLig=.TRUE.
                             ELSE
                               CYCLE
                             ENDIF
                           ENDIF
                           IF(IntSet==2) THEN
                             IDimExcl=Top_Excl%I(JJ1,1)
                             IF(ANY(Top_Excl%I(JJ1,2:1+IDimExcl)==JJ2) &
                                .AND.IDimExcl/=0) CYCLE
                             IF(HBondOnly) THEN !!! only true H-bonds survive
                               IF(HasHBond(NJJ1,NJJ2)) THEN
                                 FoundHBond=.TRUE.
                               ELSE
                                 CYCLE
                               ENDIF
                             ELSE !!! VdW and true H-bonds survive
                               IF(NJJ1==1.OR.NJJ2==1) THEN 
                                 IF(HasHBond(NJJ1,NJJ2)) THEN
                                   FoundHBond=.TRUE.
                                 ELSE 
                                   CYCLE
                                 ENDIF
                               ENDIF
                             ENDIF
                             IF(.NOT.(FoundHBond.OR.FoundMetLig)) CYCLE
                           ENDIF
                           IF(FoundHBond) THEN
                             CritDist=HBondMax
                           ELSE
                             CritDist=CritRad%D(NJJ1)+CritRad%D(NJJ2)
                           ENDIF
                           DVect%D(:)=XYZ(:,JJ1)-XYZ(:,JJ2)
                           R12_2=DOT_PRODUCT(DVect%D,DVect%D)
                           R12=SQRT(R12_2)
                           IF(R12<CritDist) THEN
                             NBond=NBond+1
                             IF(NBond>NBondEst) THEN
                               CALL ResizeBonds(BondIJAux,BondLengthAux,&
                                    HBondMarkAux,NBondEst,NatmsLoc)
                             ENDIF 
                             BondIJAux%I(1,NBond)=JJ1
                             BondIJAux%I(2,NBond)=JJ2
                             BondLengthAux%D(NBond)=R12
                             IF(FoundHBond) HBondMarkAux%I(NBond)=1
                           ENDIF
                         ENDDO
                       ENDIF
                     ENDDO
                   ENDIF
                 ENDDO
               ENDIF
             ENDDO
           ENDDO !!! central box atoms
         ENDDO
       ENDDO
     ENDDO !!! ends on central box indices
     !
     CALL New(BondIJ,(/2,NBond/))
     CALL New(BondLength,NBond)
     CALL New(HBondMark,NBond)
     DO I=1,NBond
       BondIJ%I(1:2,I)=BondIJAux%I(1:2,I)
       BondLength%D(I)=BondLengthAux%D(I)
       HBondMark%I(I)=HBondMarkAux%I(I)
     ENDDO
     ! tidy up
     CALL Delete(BondIJAux)      
     CALL Delete(BondLengthAux)      
     CALL Delete(HBondMarkAux)      
     CALL Delete(DVect)      
     IF(IntSet==2) THEN
       CALL Delete(Top_Excl)
       CALL Delete(TopCharge)
     ENDIF
   END SUBROUTINE BondList
!
!----------------------------------------------------------------
!
   SUBROUTINE ResizeBonds(BondIJ,BondLength,HBondMark,NBond,NatmsLoc)
     TYPE(INT_RNK2)     :: BondIJ,BOndIJ2
     TYPE(DBL_VECT)     :: BondLength,BondLength2
     TYPE(INT_VECT)     :: HBondMark,HBondMark2
     INTEGER            :: I,J,NBond,NBondNew,NatmsLoc
     !
     CALL New(BondIJ2,(/2,NBond/))
     CALL New(BondLength2,NBond)
     CALL New(HBondMark2,NBond)
     HBondMark2%I=0
     BondIJ2%I=BondIJ%I
     BondLength2%D=BondLength%D
     HBondMark%I=HBondMark2%I
     CALL Delete(BondIJ)
     CALL Delete(BondLength)
     CALL Delete(HBondMark)
     !
     NBondNew=NBond+NatmsLoc*10 
     CALL New(BondIJ,(/2,NBondNew/))
     CALL New(BondLength,NBondNew)
     CALL New(HBondMark,NBondNew)
     DO I=1,NBond
       BondIJ%I(1:2,I)=BondIJ2%I(1:2,I)
       BondLength%D(I)=BondLength2%D(I)
       HBondMark%I(I)=HBondMark2%I(I)
     ENDDO 
     NBond=NBondNew
     !
     CALL Delete(BondIJ2) 
     CALL Delete(BondLength2) 
     CALL Delete(HBondMark2) 
   END SUBROUTINE ResizeBonds
!
!----------------------------------------------------------------
!
   SUBROUTINE OutPList(XYZ,Top12,OutPCrit,NOutP,OutPIJKL)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INT_RNK2)              :: Top12,OutPIJKL,OutPIJKL2
     INTEGER                     :: NatmsLoc,NOutP,I,J
     INTEGER                     :: I1,I2,I3,I4
     REAL(DOUBLE)                :: A1,A2,A3,A,Conv,OutPCrit,AngleSum
     !
     NatmsLoc=SIZE(XYZ,2)
     Conv=180.D0/Pi
     NOutP=0
     DO I=1,NatmsLoc
       IF(Top12%I(I,1)==3) NOutP=NOutP+3
      !IF(Top12%I(I,1)>=3) NOutP=NOutP+3
     ENDDO
     CALL New(OutPIJKL2,(/4,NOutP/))
     !
     NOutP=0
     DO I=1,NatmsLoc
      !IF(Top12%I(I,1)>=3) THEN
       IF(Top12%I(I,1)==3) THEN
         I1=I
         CALL OutPSelect(I1,Top12,XYZ,I2,I3,I4,AngleSum)
         IF(ABS(AngleSum*Conv-360.D0)<OutPCrit) THEN
           OutPIJKL2%I(1,NOutP+1)=I2
           OutPIJKL2%I(2,NOutP+1)=I1
           OutPIJKL2%I(3,NOutP+1)=I3
           OutPIJKL2%I(4,NOutP+1)=I4
           !
           OutPIJKL2%I(1,NOutP+2)=I3
           OutPIJKL2%I(2,NOutP+2)=I1
           OutPIJKL2%I(3,NOutP+2)=I2
           OutPIJKL2%I(4,NOutP+2)=I4
           !
           OutPIJKL2%I(1,NOutP+3)=I4
           OutPIJKL2%I(2,NOutP+3)=I1
           OutPIJKL2%I(3,NOutP+3)=I2
           OutPIJKL2%I(4,NOutP+3)=I3
           NOutP=NOutP+3
         ENDIF
       ENDIF
     ENDDO
     CALL New(OutPIJKL,(/4,NOutP/))
     OutPIJKL%I(1:4,1:NOutP)=OutPIJKL2%I(1:4,1:NOutP)
     CALL Delete(OutPIJKL2)
   END SUBROUTINE OutPList
!
!-------------------------------------------------
!
   SUBROUTINE VDWFilter(BondIJ,NBond,Top_Excl)
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
     CALL Delete(BondIJ) 
     NBond=NBond2
     CALL New(BondIJ,(/2,NBond/))
     BONDIJ%I(1:2,1:NBond)=BONDIJ2%I(1:2,1:NBond)
     CALL Delete(BondIJ2) 
   END SUBROUTINE VDWFilter
!
!------------------------------------------------------------
!
   SUBROUTINE TorsionList(NatmsLoc,Top12,BondIJ, &
                           TorsionIJKL,NTorsion,HBondMark)
     !
     ! this routine generates bond-angles associated with WDV bonds
     !
     IMPLICIT NONE
     INTEGER              :: NatmsLoc,I,I1,I2,N1,N2,J,J1,J2,NBond
     INTEGER              :: NTorsion
     INTEGER              :: I1Num,I2Num,JJ,I3,N3,JJ1,JJ2,NTIni
     TYPE(INT_VECT)       :: HBondMark
     TYPE(INT_RNK2)       :: Top12   
     TYPE(INT_RNK2)       :: BondIJ  
     TYPE(INT_RNK2)       :: TorsionIJKL,TorsionIJKLAux
     LOGICAL              :: DoSelect
     !    
     DoSelect=.FALSE.
     NBond=SIZE(BondIJ%I,2)
     NTorsion=0
     DO I=1,NBond
      !IF(HBondMark%I(I)/=0) CYCLE
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       N1=Top12%I(I1,1)
       N2=Top12%I(I2,1)
       NTorsion=NTorsion+N1*N2
     ENDDO
     !    
     CALL New(TorsionIJKLAux,(/4,NTorsion/))
     !    
     NTorsion=0
     DO I=1,NBond 
      !IF(HBondMark%I(I)/=0) CYCLE
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       N1=Top12%I(I1,1)
       N2=Top12%I(I2,1)
       NTIni=NTorsion
       DO J1=1,N1
         JJ1=Top12%I(I1,J1+1)
         DO J2=1,N2
           JJ2=Top12%I(I2,J2+1)
           IF(JJ1==JJ2.OR.JJ1==I2.OR.JJ2==I1) CYCLE
           NTorsion=NTorsion+1
           TorsionIJKLAux%I(1,NTorsion)=JJ1
           TorsionIJKLAux%I(2,NTorsion)=I1
           TorsionIJKLAux%I(3,NTorsion)=I2
           TorsionIJKLAux%I(4,NTorsion)=JJ2
           IF(DoSelect) EXIT
         ENDDO
         IF(DoSelect.AND.NTIni/=NTorsion) EXIT
       ENDDO
     ENDDO
     !
     CALL New(TorsionIJKL,(/4,NTorsion/))
     DO I=1,NTorsion
       DO J=1,4 ; TorsionIJKL%I(J,I)=TorsionIJKLAux%I(J,I) ; ENDDO
     ENDDO
     CALL Delete(TorsionIJKLAux)
   END SUBROUTINE TorsionList
!
!-------------------------------------------------
!
   SUBROUTINE AngleList(NatmsLoc,Top12,AngleIJK, &
                        NAngle,BondIJ,NBond,HBondMark,DoVDW)
     !
     ! this routine generates bond-angles associated with WDV bonds
     !
     IMPLICIT NONE
     INTEGER              :: NatmsLoc,I,I1,I2,N1,N2,J,J1,J2,NBond,NAngle
     INTEGER              :: JJ,II1,II2,K,IC,N,IC2,IntSet
     TYPE(INT_RNK2)       :: Top12   
     TYPE(INT_RNK2)       :: AngleIJK,AngleIJKAux
     TYPE(INT_RNK2)       :: BondIJ
     TYPE(INT_VECT)       :: Sign,HBondMark   
     LOGICAL              :: DoVDW
     !    
     CALL New(Sign,NatmsLoc)
     Sign%I=0
     NAngle=0
     DO I=1,NBond
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       N1=Top12%I(I1,1)
       N2=Top12%I(I2,1)
       NAngle=NAngle+N1+N2
     ENDDO
     !    
     CALL New(AngleIJKAux,(/3,NAngle/))
     !    
     NAngle=0
     DO I=1,NBond
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       N1=Top12%I(I1,1)
       N2=Top12%I(I2,1)
     ! IF(HBondMark%I(I)/=0) THEN
     !   IF(N1/=0.OR.N2/=0) CYCLE
     ! ENDIF
       IF(Sign%I(I1)==0) THEN
         DO J1=1,N1
           II1=Top12%I(I1,J1+1)
           DO J2=J1+1,N1
             II2=Top12%I(I1,J2+1)
             IF(DoVDW.AND.(II1/=I2.AND.II2/=I2)) CYCLE
             NAngle=NAngle+1
             AngleIJKAux%I(1,NAngle)=II1
             AngleIJKAux%I(2,NAngle)=I1
             AngleIJKAux%I(3,NAngle)=II2
           ENDDO
         ENDDO
        !IF(N1>=4) NAngle=NAngle-1
         Sign%I(I1)=1
       ENDIF
       !
       IF(Sign%I(I2)==0) THEN
         DO J1=1,N2
           II1=Top12%I(I2,J1+1)
           DO J2=J1+1,N2
             II2=Top12%I(I2,J2+1)
             IF(DoVDW.AND.(II1/=I1.AND.II2/=I1)) CYCLE
             NAngle=NAngle+1
             AngleIJKAux%I(1,NAngle)=II1
             AngleIJKAux%I(2,NAngle)=I2
             AngleIJKAux%I(3,NAngle)=II2
           ENDDO
         ENDDO
        !IF(N2>=4) NAngle=NAngle-1
         Sign%I(I2)=1
       ENDIF
     ENDDO
     !
     CALL New(AngleIJK,(/3,NAngle/))
     DO I=1,NAngle
       DO J=1,3 ; AngleIJK%I(J,I)=AngleIJKAux%I(J,I) ; ENDDO
     ENDDO
     CALL Delete(Sign)
     CALL Delete(AngleIJKAux)
   END SUBROUTINE AngleList
!
!-------------------------------------------------
!
   SUBROUTINE GetIntCs(XYZ,AtNumIn,IntCs,NIntC,Refresh,& 
                       SCRPath,CtrlCoord,CtrlConstr)
     !
     ! This subroutine constructs the IntCs array, which holds
     ! definitions of internal coordinates to be used in the 
     ! forthcoming geometry manipulation procedure.
     ! Refresh=1 : Refresh all definitions
     !        =2 : Refresh only definitions based on VDW interaction
     !        =3 : Do not refresh definitions, use the one from HDF
     !        =4 : Refresh/generate only the covalent coordinates
     !        =5 : only the extra coordinates from input
     ! WARNING! In the present version 
     ! bending -> linear bending transitions are 
     ! always checked and refreshed
     ! Later, check also linear bending -> bending transitions !
     ! 
     IMPLICIT NONE
     TYPE(INTC)                  :: IntCs,IntC_Cov,IntC_VDW
     TYPE(INTC)                  :: IntC_Extra,IntC_New
     TYPE(INTC)                  :: IntC_Cart
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(CoordCtrl)             :: CtrlCoord
     TYPE(Constr)                :: CtrlConstr
     LOGICAL                     :: DoFixMM
     REAL(DOUBLE),DIMENSION(:)   :: AtNumIn
     INTEGER                     :: NIntC,NIntC_Cov,NIntC_VDW
     INTEGER                     :: NIntC_Extra,NNew,Nintc_New
     INTEGER                     :: NIntC_Cart
     INTEGER                     :: I,J,K,Refresh,NatmsLoc,II,ILast,III
     INTEGER                     :: I1,I2,I3,I4,NMax12
     INTEGER                     :: NStreGeOp,NBendGeOp
     INTEGER                     :: NLinBGeOp,NOutPGeOp,NTorsGeOp
     TYPE(INT_VECT)              :: AtNum
     TYPE(INT_RNK2)              :: Top12
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     !
     NatmsLoc=SIZE(XYZ,2)
     !
     NIntC_Cart=0
     !
     CALL New(AtNum,NatmsLoc)
     DO I=1,NatmsLoc
       AtNum%I(I)=INT(AtNumIn(I))
     ENDDO
     !
     IF(Refresh==1) Then !!! Total refresh
       !define covalent bonding scheme
            CALL DefineIntCoos(NatmsLoc,XYZ,AtNum%I,1,Refresh, &
                               IntC_Cov,NIntC_Cov,CtrlCoord,SCRPath)
       !define Van der Waals bonding scheme
            CALL DefineIntCoos(NatmsLoc,XYZ,AtNum%I,2,Refresh, &
                               IntC_VDW,NIntC_VDW,CtrlCoord,SCRPath)
       !
     ELSE IF(Refresh==2) THEN !!! refresh VDW terms
       !
       ! Refresh only the VDW terms
       !
       CALL ReadIntCs(IntC_Cov,TRIM(SCRPath)//'IntC_Cov')
       NIntC_Cov=SIZE(IntC_Cov%Def)
       CALL DefineIntCoos(NatmsLoc,XYZ,AtNum%I,2,Refresh, &
                            IntC_VDW,NIntC_VDW,CtrlCoord,SCRPath)
       !
     ELSE IF(Refresh==3) THEN !!! no refresh, get everything from disk
       CALL ReadIntCs(IntC_Cov,TRIM(SCRPath)//'IntC_Cov')
       CALL ReadIntCs(IntC_VDW,TRIM(SCRPath)//'IntC_VDW')
       NIntC_Cov=SIZE(IntC_Cov%Def)
       NIntC_VDW=SIZE(IntC_VDW%Def)
       !
     ELSE IF(Refresh==4) THEN !!! covalent bonds only
       CALL DefineIntCoos(NatmsLoc,XYZ,AtNum%I,1,Refresh, &
                          IntC_Cov,NIntC_Cov,CtrlCoord,SCRPath)
       NIntC_VDW=0 
       !
     ELSE IF(Refresh==5) THEN !!! use only extra coords from input
       NIntC_Cov=0 
       NIntC_VDW=0 
     ENDIF
     !
     ! Save Covalent and VDW Terms onto disk
     !
     CALL WriteIntCs(IntC_Cov,TRIM(SCRPath)//'IntC_Cov')
     CALL WriteIntCs(IntC_VDW,TRIM(SCRPath)//'IntC_VDW')  
     !
     ! Get extra internal coordinates and constraints
     !
     NIntC_Extra=CtrlCoord%NExtra
     IF(NIntC_Extra/=0) &
       CALL ReadIntCs(IntC_Extra,TRIM(SCRPath)//'IntC_Extra')
     !
     ! Merge INTC arrays
     !
     CtrlCoord%NCov=NIntC_Cov
     NIntC=NIntC_Cov+NIntC_VDW+NIntC_Extra+NIntC_Cart
     !
     CALL New(IntCs,NIntC)
     !
     ILast=0
     IF(NIntC_Cov/=0) THEN
       CALL Set_INTC_EQ_INTC(IntC_Cov,IntCs,1,NIntC_Cov,ILast+1)
       CALL Delete(IntC_Cov)
     ENDIF
       ILast=NIntC_Cov
     IF(NIntC_VDW/=0) THEN 
       CALL Set_INTC_EQ_INTC(IntC_VDW,IntCs,1,NIntC_VDW,ILast+1)
       CALL Delete(IntC_VDW)
     ENDIF
     !
     IF(NIntC_Extra/=0) THEN
       ILast=ILast+NIntC_VDW
       CALL Set_INTC_EQ_INTC(IntC_Extra,IntCs,1,NIntC_Extra,ILast+1)
       CALL Delete(IntC_Extra)
     ENDIF
     !
     IF(NIntC_Cart/=0) THEN
       ILast=MAX(ILast+NIntC_Extra,ILast+NIntC_VDW)
       CALL Set_INTC_EQ_INTC(IntC_Cart,IntCs,1,NIntC_Cart,ILast+1)
       CALL Delete(IntC_Cart)
     ENDIF
     !
     IF(.NOT.(NIntC==0.OR.NIntC_Extra==0.OR.Refresh==5)) THEN
       CALL CleanINTC(IntCs,NIntC,NIntC_Cov,NIntC_VDW,NIntC_Extra)
     ENDIF             
     !
     ! Set active all internal coords defd so far.
     ! 'Linear torsions will be deactivated when 
     ! their value is calculated
     ! by some subroutine (eg. INTCValue).
     ! The set of active coordinates may vary 
     ! during the process of optimization.
     !
     IntCs%Active=.TRUE.
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
     CtrlCoord%NStre=NStreGeOp
     CtrlCoord%NBend=NBendGeOp
     CtrlCoord%NLinB=NLinBGeOp
     CtrlCoord%NOutP=NOutPGeOp
     CtrlCoord%NTors=NTorsGeOp
     !
     ! Save current internals to disk
     !
     CALL WriteIntCs(IntCs,TRIM(SCRPath)//'IntCs')
     CALL Delete(AtNum)
   END SUBROUTINE GetIntCs
!
!-------------------------------------------------------------------
!
   SUBROUTINE CleanINTC(IntCs,NIntC,NIntC_Cov,NIntC_VDW,NIntC_Extra)
     !
     ! Now filter out repeated definitions, 
     ! which may have occured in INTC_Extra
     ! keep those values, which were defd. in extras.
     ! Do not check for Cartesians.
     !
     TYPE(INTC) :: IntCs,IntC_New
     INTEGER    :: NIntC_Cov,NIntC_VDW,NIntC,I,J,NNew,II,III,ILast
     INTEGER    :: NIntC_Extra
     !
     ILast=NIntC_Cov+NIntC_VDW
     II=0
     DO III=1,NIntC
       DO I=ILast+1,ILast+NIntC_Extra
           IF(IntCs%Def(I)(1:5)=='CART ') CYCLE
         DO J=1,ILast
           IF(IntCs%Def(J)(1:5)=='BLANK') CYCLE
           IF(IntCs%Atoms(J,1)==IntCs%Atoms(I,1).AND.&
              IntCs%Atoms(J,2)==IntCs%Atoms(I,2).AND.&
              IntCs%Atoms(J,3)==IntCs%Atoms(I,3).AND.&
              IntCs%Atoms(J,4)==IntCs%Atoms(I,4)) THEN
              II=II+1
              IntCs%Def(J)(1:5)='BLANK'
              IntCs%Atoms(J,1:4)=0      
           ENDIF
         ENDDO
       ENDDO
     ENDDO
     !
     ! Compress IntCs array, get rid of BLANK-s
     !
     IF(II/=0) THEN
       IF(ANY(IntCs%Def(:)(1:5)=='BLANK')) THEN
         NNew=NIntC
         DO I=1,NIntc
           IF(IntCs%Def(I)(1:5)=='BLANK') NNew=NNew-1
         ENDDO
         CALL New(IntC_New,NNew)
         NNew=0    
         DO I=1,NIntc
           IF(IntCs%Def(I)(1:5)/='BLANK') THEN
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
     ENDIF
   END SUBROUTINE CleanINTC
!
!-------------------------------------------------------
!
   SUBROUTINE INTCValue(IntCs,XYZ,LinCrit)
     !
     ! Determine value of internal coordinates.
     ! Input coordintes are now in atomic units!
     ! 
     IMPLICIT NONE
     TYPE(INTC) :: IntCs
     INTEGER :: NIntCs,I,J,K,L,I1,I2,I3,I4,NatmsLoc
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,LinCrit
     !
     IF(AllocQ(IntCs%Alloc)) THEN
       NIntCs=SIZE(IntCs%Def(:))
     ELSE
       NIntCs=0
       RETURN
     ENDIF
     NatmsLoc=SIZE(XYZ,2)
     !
     IntCs%Value=Zero 
     !
     DO I=1,NIntCs
       IF(.NOT.IntCs%Active(I)) THEN
         IntCs%Value(I)=Zero 
         CYCLE
       ENDIF
       I1=IntCs%Atoms(I,1)
       I2=IntCs%Atoms(I,2)
       I3=IntCs%Atoms(I,3)
       I4=IntCs%Atoms(I,4)
       IF(IntCs%Def(I)(1:4)=='STRE') THEN
         CALL STREValue(XYZ(1:3,I1),XYZ(1:3,I2),IntCs%Value(I))
         !
       ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
         CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3), &
           IntCs%Value(I))
         !
       ELSE IF(IntCs%Def(I)(1:5)=='LINB1') THEN
         CALL LinB2(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),&
           Value1=IntCs%Value(I),Value2=IntCs%Value(I+1))
         !
       ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
         CALL TORSValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),&
           XYZ(1:3,I4),LinCrit,IntCs%Value(I),IntCs%Active(I))
         !
       ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
         CALL OUTPValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I4),&
           XYZ(1:3,I3),LinCrit,IntCs%Value(I),IntCs%Active(I))
         !
       ELSE IF(IntCs%Def(I)(1:5)=='CARTX') THEN
         IntCs%Value(I)=XYZ(1,I1)
       ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
         IntCs%Value(I)=XYZ(2,I1)
       ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
         IntCs%Value(I)=XYZ(3,I1)
       ENDIF
       !
     ENDDO 
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
     REAL(DOUBLE) :: sinj,SMI,SMK,SUM,AIJK
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
     IF(ABS(DotJ)>0.999999D0) THEN
       DotJ=SIGN(ABS(DotJ),DotJ)
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
   SUBROUTINE OUTPValue(XIin,XJin,XKin,XLin,LinCrit,VOutP,Active)
     REAL(DOUBLE),DIMENSION(1:3) :: XIin,XJin,XKin,XLin
     REAL(DOUBLE)                :: LinCrit,VOutP,Sum
     LOGICAL                     :: Active
     !
     CALL TORSValue(XIin,XJin,XKin,XLin,LinCrit,VOutP,Active)
     !
     IF(VOutP>Zero) THEN
       VOutP=PI-VOutP
     ELSE
       VOutP=-(PI+VOutP)
     ENDIF
   END SUBROUTINE OUTPValue
!
!-------------------------------------------------------------
!
   SUBROUTINE TORSValue(XIin,XJin,XKin,XLin,LinCrit,VTors,Active)
     !
     REAL(DOUBLE),DIMENSION(1:3) :: XIin,XJin,XKin,XLin
     TYPE(DBL_VECT)              :: XI,XJ,XK,XL,V1,V2
     TYPE(DBL_RNK2)              :: Rot
     REAL(DOUBLE)                :: VTors,CosPhi,Sum,Conv
     REAL(DOUBLE)                :: V1ABS,V2ABS,LinCrit,TwoPi
     INTEGER                     :: I,J,K,L
     LOGICAL                     :: Active,ActiveIn,NonLinear
     !
     Conv=180.D0/PI
     TwoPi=Two*PI
     ActiveIn=Active
     NonLinear=.TRUE.
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
     ! Check for linearity of I-J-K
     !
     V1%D=XI%D-XJ%D
     V2%D=XJ%D-XK%D
     V1ABS=SQRT(DOT_PRODUCT(V1%D,V1%D))
     V2ABS=SQRT(DOT_PRODUCT(V2%D,V2%D))
     IF(ABS(V1ABS)<1.D-3.OR.ABS(V2ABS)<1.D-3) THEN
       Active=.FALSE. 
       NonLinear=.FALSE.
     ENDIF
     SUM=ABS(DOT_PRODUCT(V1%D,V2%D))/V1ABS/V2ABS
     SUM=ACOS(SUM)*Conv
     IF(Sum<LinCrit) THEN
       Active=.FALSE. 
       NonLinear=.FALSE.
     ENDIF
     ! 
     ! Check for linearity of J-K-L
     !
     V1%D=XL%D-XK%D
     V2%D=XK%D-XJ%D
     V1ABS=SQRT(DOT_PRODUCT(V1%D,V1%D))
     V2ABS=SQRT(DOT_PRODUCT(V2%D,V2%D))
     IF(ABS(V1ABS)<1.D-3.OR.ABS(V2ABS)<1.D-3) THEN
       Active=.FALSE. 
       NonLinear=.FALSE.
     ENDIF
     SUM=ABS(DOT_PRODUCT(V1%D,V2%D))/V1ABS/V2ABS
     SUM=ACOS(SUM)*Conv
     IF(Sum<LinCrit) THEN
       Active=.FALSE. 
       NonLinear=.FALSE.
     ENDIF
     !
     IF(NonLinear) THEN
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
       IF(ABS(CosPhi)>0.999999D0) THEN
         CosPhi=SIGN(ABS(CosPhi),CosPhi)
       ENDIF
       VTors=ACOS(CosPhi)
       !
       ! Take orientation of torsional chain into account
       !
       IF(XL%D(2)>Zero) VTors=TwoPi-VTors
       Sum=TwoPi-VTors
       IF(Sum<VTors) VTors=-Sum
       !
     ENDIF
     CALL Delete(V2)
     CALL Delete(V1)
     CALL Delete(Rot)
     CALL Delete(XI)
     CALL Delete(XJ)
     CALL Delete(XK)
     CALL Delete(XL)
   END SUBROUTINE TORSValue
!
!--------------------------------------------------------------------
!
   SUBROUTINE GetSpBMatr(XYZ,B,SpB) 
     !
     ! Generate vibrational B matrix in quasi-sparse 
     ! TYPE(BMATR) representation and transform into sparse one
     !
     TYPE(BCSR)                           :: SpB
     INTEGER                              :: NIntC,NCart,I,NatmsLoc,NZ
     TYPE(BMATR)                          :: B
     REAL(DOUBLE),DIMENSION(:,:)          :: XYZ
     !
     NIntC=SIZE(B%IB,1)
     NatmsLoc=SIZE(XYZ,2)
     !
     NCart=3*NatmsLoc
     !
     ! Now, turn B matrix into sparse blocked representation
     !
     CALL Set_BCSR_EQ_BMATR(SpB,B)
     !     CALL PPrint(SpB,'SpB',Unit_O=6)
     !
     ! Generate sparse transpose B
     !
     !     CALL PPrint(SpB,'SpB2',Unit_O=6)
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
     ! Map back angles into the ranges the iterative 
     ! back-transformation can cope with.
     TwoPi=Two*PI
     !
     DO I=1,NIntC
       SUM=VectAux(I)
       IF(IntCs%Def(I)(1:4)=='BEND'.OR. &
          IntCs%Def(I)(1:4)=='LINB') THEN 
         J=INT(SUM/TwoPI)
         SUM=SUM-J*TwoPi
         IF(SUM<Zero) SUM=TwoPi+SUM
         IF(SUM>PI) SUM=TwoPI-SUM
       ELSE IF(IntCs%Def(I)(1:4)=='TORS'.OR. &
               IntCs%Def(I)(1:4)=='OUTP') THEN
         J=INT(SUM/TwoPi)
         SUM=SUM-J*TwoPi
         IF(SUM<-PI) SUM=TwoPi+SUM
         IF(SUM> PI) SUM=SUM-TwoPi
       ENDIF
       VectAux(I)=SUM
     ENDDO
   END SUBROUTINE MapBackAngle
!
!-------------------------------------------------------
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
     ! so that all coordinates be greater than zero for central cell 
     !
     ! WARNING! Current implementation is optimal for 
     ! cells with rectangular
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
     ! First, calculate maximum Cartesian 
     ! extension of the elementary cell
     !
     DX=GMLoc%PBC%BoxShape(1,1)
     DY=GMLoc%PBC%BoxShape(1,2)
     DZ=GMLoc%PBC%BoxShape(1,3)
     DO I=2,3
       DX=MAX(DX,GMLoc%PBC%BoxShape(I,1))
       DY=MAX(DY,GMLoc%PBC%BoxShape(I,2))
       DZ=MAX(DZ,GMLoc%PBC%BoxShape(I,3))
     ENDDO
     MaxCart=DX
     MaxCart=MAX(MaxCart,DY) 
     MaxCart=MAX(MaxCart,DZ) 
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
   END SUBROUTINE LJCell
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
   END SUBROUTINE SetOneLJCell
!
!---------------------------------------------------------------------
!
   SUBROUTINE GetFullBMatr(CartsLoc,B,FullB,FullBt)
     !
     ! Generate vibrational B matrix in quasi-sparse 
     ! TYPE(BMATR) representation
     ! and transform into sparse one
     !
     TYPE(CRDS) :: GMLoc
     INTEGER    :: NIntC,NCart,I,J,K,NatmsLoc,JJ,KK,LL,II
     TYPE(BMATR):: B
     TYPE(DBL_RNK2) :: FullB,FullBt
     REAL(DOUBLE),DIMENSION(:,:) :: CartsLoc
     REAL(DOUBLE) :: SUM
     !
     NatmsLoc=SIZE(CartsLoc,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(B%IB,1)
     !
     ! Turn B matrix into full representation
     !
     FullB%D=Zero
     FullBt%D=Zero
     DO I=1,NIntC
       DO J=1,4 
         K=B%IB(I,J) 
         IF(K==0) CYCLE
         DO JJ=1,3
           KK=3*(K-1)+JJ
           LL=3*(J-1)+JJ
           SUM=B%B(I,LL)
           FullB%D(I,KK)=SUM
           FullBt%D(KK,I)=SUM
         ENDDO
       ENDDO
       !write(*,200) i,b%ib(i,1:12)
       !write(*,201) i,b%b(i,1:12)
       !write(*,100) i,(fullb%d(i,k),k=1,ncart)
       !
     ENDDO
     !100 format(I4,100F7.3)
     !200 format(I4,100I4)
     !201 format(I4,100F7.3)
     !call pprint(fullb,'fullb= ',unit_o=6) 
     !
   END SUBROUTINE GetFullBMatr
!
!-------------------------------------------------------
!
   SUBROUTINE GetFullGcInv(B,FullGcInv,NIntC,NCart)
     !
     TYPE(DBL_RNK2) :: FullGcInv
     TYPE(BMATR)    :: B
     TYPE(DBL_RNK2) :: FullGc,TestMat
     INTEGER :: NIntC,NCart,I,J,K
     REAL(DOUBLE) :: SUM
     TYPE(INT_VECT) :: ISpB,JSpB,IGc,JGc
     TYPE(DBL_VECT) :: ASpB,AGc
     !
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     CALL GetGc(NCart,ISpB,JSpB,ASpB,IGc,JGc,AGc)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL New(TestMat,(/NCart,NCart/))
     CALL Sp1x1ToFull(IGc%I,JGc%I,AGc%D,NCart,NCart,FullGc)
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
     !
     ! Get SVD inverse
     !
     CALL SetDSYEVWork(NCart)
     CALL FunkOnSqMat(NCart,Inverse,FullGc%D,FullGcInv%D)
     CALL UnSetDSYEVWork()
     !
     ! Test inverse
     !
     !    CALL DGEMM_NNc(NCart,NCart,NCart,One,Zero,&
     !      FullGcInv%D,FullGc%D,TestMat%D)
     !    DO I=1,NCart
     !      write(*,100) (TestMat%D(I,J),J=1,NCart)
     !    ENDDO
     !100 FORMAT(10F8.4)
     !
     CALL Delete(TestMat)
     CALL Delete(FullGc)
   END SUBROUTINE GetFullGcInv
!
!-------------------------------------------------------
!
   SUBROUTINE CartToInternal(XYZ,IntCs,VectCart,VectInt, &
     TrfGrd,CtrlCoord,CtrlTrf,Print,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ   
     REAL(DOUBLE),DIMENSION(:)    :: VectCart,VectInt
     TYPE(DBL_VECT)  :: VectCartAux,VectIntAux
     TYPE(DBL_VECT)  :: VectCartAux2,VectIntAux2
     REAL(DOUBLE)    :: DiffMax,RMSD
     REAL(DOUBLE)    :: Sum
     INTEGER         :: NCart,I,II,J,NIntC
     INTEGER         :: NatmsLoc,Print
     TYPE(INTC)      :: IntCs
     TYPE(Cholesky)  :: CholData
     TYPE(BMATR)     :: B
     TYPE(GrdTrf)    :: TrfGrd
     TYPE(CoordCtrl) :: CtrlCoord
     TYPE(TrfCtrl)   :: CtrlTrf
     CHARACTER(LEN=*):: SCRPath
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc        
     NIntC=SIZE(IntCs%Def)
     !
     CALL New(VectCartAux,NCart)
     CALL New(VectCartAux2,NCart)
     CALL New(VectIntAux,NIntC)
     CALL New(VectIntAux2,NIntC)
     !
     VectInt=Zero
     !
     IF(Print>=DEBUG_GEOP_MIN) THEN
       WRITE(*,111) NIntC 
       WRITE(Out,111) NIntC
       111 FORMAT('Gradient transformation, No. Int. Coords= ',I7)
       IF(.NOT.CtrlTrf%DoClssTrf) THEN
         WRITE(*,112) CtrlTrf%ThreeAt
         WRITE(Out,112) CtrlTrf%ThreeAt
       ENDIF
       112 FORMAT('Three-atoms reference system used, atoms are ',3I4)
     ENDIF
     !
     ! Cartesian --> Internal transformation
     !
     ! Get B matrix and Bt*B inverse
     !
     CALL GetBMatInfo(SCRPath,NIntC,B,CholData)
     IF(.NOT.CtrlTrf%DoClssTrf) THEN
       CALL RotToAt(VectCart,CtrlTrf%RotAt2ToX)
       CALL RotToAt(VectCart,CtrlTrf%RotAt3ToXY)
       CALL ZeroTrRots(VectCart,CtrlTrf%ThreeAt)
     ENDIF
     !
     DO II=1,TrfGrd%MaxIt_GrdTrf
       !
       VectCartAux%D=Zero
       !
       ! gc-Bt*gi
       !
       CALL CALC_BxVect(B,VectInt,VectCartAux%D,Trp_O=.TRUE.)
       VectCartAux%D=VectCart-VectCartAux%D
       !
       ! GcInv*[gc-Bt*gi]
       !
       CALL CALC_GcInvCartV(CholData,VectCartAux%D,VectCartAux2%D)
       !
       ! B*GcInv*[gc-Bt*gi]
       !
       CALL CALC_BxVect(B,VectIntAux%D,VectCartAux2%D)
       !
       ! Check convergence
       !
       DiffMax=Zero
       DO I=1,NIntC ; DiffMax=MAX(DiffMax,ABS(VectIntAux%D(I))) ; ENDDO
       RMSD=DOT_PRODUCT(VectIntAux%D,VectIntAux%D)
       RMSD=SQRT(RMSD/DBLE(NIntC))
       !
       ! IF DiffMax is too large, eg. due to the 'bad' quality 
       ! of the preconditioner, rescale gradients
       !
       IF(DiffMax>TrfGrd%MaxGradDiff) THEN
         IF(Print>=DEBUG_GEOP_MIN) THEN
           WRITE(*,*) 'Rescale Step from ',DiffMax,' to ',TrfGrd%MaxGradDiff
           WRITE(Out,*) 'Rescale Step from ',DiffMax,' to ',TrfGrd%MaxGradDiff
         ENDIF
         SUM=TrfGrd%MaxGradDiff/DiffMax
         VectIntAux%D(:)=SUM*VectIntAux%D(:)
         DiffMax=TrfGrd%MaxGradDiff
       ENDIF
       !
       ! gi+B*GcInv*[gc-Bt*gi]
       !
       VectInt=VectInt+VectIntAux%D
       !
       ! Review iteration
       !
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,110) II,DiffMax,RMSD
         WRITE(Out,110) II,DiffMax,RMSD
       ENDIF
       110  FORMAT('Grad Trf, step= ',I3,' MaxChange= ',F12.6,&
                   ' ChangeNorm= ',F12.6)
       !      
       IF(DiffMax<TrfGrd%GrdTrfCrit) EXIT 
     ENDDO
     !
     IF(II>=TrfGrd%MaxIt_GrdTrf) THEN
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,777) 
         WRITE(*,778) 
         WRITE(Out,777) 
         WRITE(Out,778) 
         777 FORMAT('Stop Gradient Trf, max. number '//&
                      'of Iterations exceeded!')
         778 FORMAT('Use current gradient vector!')
       ENDIF
     ELSE
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,120) II
         WRITE(Out,120) II
       ENDIF
     ENDIF
     120  FORMAT('Gradient transformation converged in ',I3,' steps')
     !
     ! Tidy up
     !
     CALL Delete(VectIntAux2)
     CALL Delete(VectIntAux)
     CALL Delete(VectCartAux2)
     CALL Delete(VectCartAux)
     CALL DeleteBMatInfo(B,CholData)
     !
   END SUBROUTINE CartToInternal
!
!------------------------------------------------------------------
!
   SUBROUTINE InternalToCart(XYZ,IntCs,VectInt,Print, &
       GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:) :: VectInt
     TYPE(DBL_VECT)            :: VectCart
     TYPE(DBL_VECT)            :: VectCartAux,VectIntAux
     TYPE(DBL_VECT)            :: VectCartAux2
     TYPE(DBL_VECT)            :: VectIntReq
     TYPE(DBL_RNK2)            :: ActCarts
     REAL(DOUBLE)              :: DiffMax,RMSD,RMSDOld
     REAL(DOUBLE)              :: Sum
     REAL(DOUBLE)              :: ConstrMax,ConstrRMS
     REAL(DOUBLE)              :: ConstrRMSOld
     REAL(DOUBLE)              :: ConstrMaxCrit,RMSCrit
     INTEGER                   :: NCart,I,IStep,J,NIntC,NConstr
     INTEGER                   :: NatmsLoc
     INTEGER                   :: NCartConstr
     TYPE(INTC)                :: IntCs
     TYPE(BMATR)               :: B
     LOGICAL                   :: RefreshB,RefreshAct
     LOGICAL                   :: DoIterate
     TYPE(Cholesky)            :: CholData
     TYPE(BackTrf)             :: GBackTrf
     TYPE(Constr)              :: GConstr
     TYPE(TrfCtrl)             :: GTrfCtrl
     TYPE(CoordCtrl)           :: GCoordCtrl
     CHARACTER(LEN=*)          :: SCRPath
     LOGICAL                   :: DoClssTrf,Print2
     INTEGER                   :: Print
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc   
     NIntC=SIZE(IntCs%Def)
     Print2= Print>=DEBUG_GEOP_MAX
     !
     ! Refresh B matrix during iterative back-trf?
     !
     RefreshB=.TRUE.
     RefreshAct=.TRUE.
     !
     ! Auxiliary arrays
     !
     CALL New(ActCarts,(/3,NatmsLoc/))
     CALL New(VectCart,NCart)
     CALL New(VectCartAux,NCart)
     CALL New(VectCartAux2,NCart)
     CALL New(VectIntAux,NIntC)
     CALL New(VectIntReq,NIntC)
     !
     ! The required new value of internal coordinates
     !
     VectIntReq%D=VectInt
     CALL MapBackAngle(IntCs,NIntC,VectIntReq%D) 
     !
     CALL INTCValue(IntCs,XYZ,GCoordCtrl%LinCrit)
     !
     !initialization of new Cartesians
     !
     IF(GTrfCtrl%DoClssTrf) THEN
       ActCarts%D=XYZ
       CALL CartRNK2ToCartRNK1(VectCart%D,ActCarts%D)
     ELSE
       CALL CartRNK2ToCartRNK1(VectCart%D,XYZ)
       CALL TranslToAt1(VectCart%D,GTrfCtrl%ThreeAt)
       CALL RotToAt(VectCart%D,GTrfCtrl%RotAt2ToX)
       CALL RotToAt(VectCart%D,GTrfCtrl%RotAt3ToXY)
       !transform Cartesian constraints!
       CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
       CALL ReSetConstr(IntCs,ActCarts%D)
     ENDIF
     !
     ! Internal --> Cartesian transformation
     !
     IF(Print>=DEBUG_GEOP_MIN) THEN
       WRITE(*,450) NIntC
       WRITE(Out,450) NIntC
     ENDIF
     450 FORMAT('Iterative back-transformation, No. Int. Coords=',I7)
     !
     ConstrMax=GConstr%ConstrMaxCrit*10.D0
     ConstrRMS=1.D0
     ConstrRMSOld=2.D0
     RMSD=1.D+9
     !
     DO IStep=1,GBackTrf%MaxIt_CooTrf
       !
       ! Get B and refresh values of internal coords
       !
       IF(IStep==1) THEN
         CALL GetBMatInfo(SCRPath,NIntC,B,CholData)
       ELSE IF(RefreshB.AND.RefreshAct) THEN
         CALL RefreshBMatInfo(IntCs,ActCarts%D,GTrfCtrl, &
                            GCoordCtrl,Print,SCRPath,DoNewRef_O=.FALSE.)
         CALL GetBMatInfo(SCRPath,NIntC,B,CholData)
       ELSE
         CALL INTCValue(IntCs,ActCarts%D,GCoordCtrl%LinCrit)
       ENDIF
       !
       ! Calculate difference between required and actual internals
       !
       VectIntAux%D=VectIntReq%D-IntCs%Value
       CALL MapAngleDispl(IntCs,NIntC,VectIntAux%D) 
       !
       ! Check convergence on constraints
       !
       IF(GConstr%NConstr/=0.AND..NOT.GConstr%DoLagr) THEN
         ConstrRMSOld=ConstrRMS
         CALL ConstrConv(IntCs,VectIntAux%D,ConstrMax,ConstrRMS)
       ENDIF
       !
       ! Do transformation
       ! 
       ! Bt*[phi_r-phi_a]
       !
       CALL CALC_BxVect(B,VectIntAux%D,VectCartAux%D,Trp_O=.TRUE.)
       !
       ! GcInv*Bt*[phi_r-phi_a]
       !
       CALL CALC_GcInvCartV(CholData,VectCartAux%D,VectCartAux2%D)
       !
       ! Project out rotations and translations 
       !
       IF(.NOT.GTrfCtrl%DoClssTrf) THEN 
         IF(GTrfCtrl%DoTranslOff) &
           CALL TranslsOff(VectCartAux2%D,Print2)
         IF(GTrfCtrl%DoRotOff) &
           CALL RotationsOff(VectCartAux2%D,ActCarts%D,Print2)
       ENDIF
       !
       ! Check convergence
       !
       RMSDOld=RMSD
       CALL ScaleDispl(VectCartAux2%D,GBackTrf%MaxCartDiff, &
                       DiffMax,RMSD) 
       !
       ! Refresh B matrix?  
       !
       IF(DiffMax>GBackTrf%DistRefresh) THEN
         RefreshAct=.TRUE.
       ELSE
         RefreshAct=.FALSE.
       ENDIF
       !
       ! Modify Cartesians
       !
       IF(.NOT.GConstr%DoLagr) THEN
         CALL SetCartConstr(VectCart%D,VectCartAux2%D, &
                            IntCs,GConstr%NCartConstr)
       ENDIF
       VectCart%D=VectCart%D+VectCartAux2%D
       CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
       !
       ! Review iteration
       !
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,210) IStep,DiffMax,RMSD
         WRITE(Out,210) IStep,DiffMax,RMSD
       ENDIF
       210  FORMAT('Step= ',I3,'   Max_DX= ',F12.6,'  X_RMSD= ',F12.6)
       !      
       CALL BackTrfConvg(GConstr,GBackTrf, &
         DoIterate,DiffMax,RMSD,RMSDOld,ConstrMax, &
         ConstrRMS,ConstrRMSOld,IStep,RefreshAct)
       !
       IF(DoIterate) THEN
         IF(RefreshB.AND.RefreshAct) CALL DeleteBMatInfo(B,CholData)
       ELSE
         EXIT
       ENDIF
     ENDDO          
     !
     IF(IStep>=GBackTrf%MaxIt_CooTrf) THEN
       IF(RMSD>0.01D0) THEN
         CALL Halt('Iterative backtransformation did not converge')
       ENDIF
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,180) 
         WRITE(Out,180) 
         WRITE(*,190) 
         WRITE(Out,190) 
       ENDIF
     ELSE
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,220) IStep
         WRITE(Out,220) IStep
       ENDIF
     ENDIF
     180  FORMAT('Stop Coord Back-Trf, max. number of Iterations exceeded!')
     190  FORMAT('Use Current Geometry!')
     220  FORMAT('Coordinate back-transformation converged in ',&
                 I3,' steps')
     !
     ! Fill new Cartesians into XYZ  
     !
     IF(.NOT.GTrfCtrl%DoClssTrf) THEN
       CALL CartRNK2ToCartRNK1(VectCart%D,ActCarts%D)
       CALL RotToAt(VectCart%D,GTrfCtrl%RotAt3ToXY,Rev_O=.TRUE.)
       CALL RotToAt(VectCart%D,GTrfCtrl%RotAt2ToX,Rev_O=.TRUE.)
       CALL TranslToAt1(VectCart%D,GTrfCtrl%ThreeAt, &
                        Vect_O=-GTrfCtrl%TranslAt1)
       CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
       CALL ReSetConstr(IntCs,ActCarts%D)
     ENDIF
     XYZ=ActCarts%D
     !
     ! Final internal coordinates
     !
     !CALL INTCValue(IntCs,XYZ,GCoordCtrl%LinCrit)
     !CALL PrtIntCoords(IntCs,IntCs%Value,'Final Internals')
     !
     ! Tidy up
     !
     CALL Delete(VectIntReq)
     CALL Delete(VectIntAux)
     CALL Delete(VectCartAux2)
     CALL Delete(VectCartAux)
     CALL Delete(VectCart)
     CALL Delete(ActCarts)
     CALL DeleteBMatInfo(B,CholData)
   END SUBROUTINE InternalToCart
!
!----------------------------------------------------------
!
   SUBROUTINE PrtIntCoords(IntCs,Value,CHAR)
     !
     TYPE(INTC) :: IntCs
     INTEGER    :: I,NIntC,J
     REAL(DOUBLE),DIMENSION(:) :: Value
     REAL(DOUBLE) :: SUM,SumConstr,Conv
     CHARACTER(LEN=*) :: CHAR   
     !
     Conv=180.D0/PI
     IF(AllocQ(IntCs%Alloc)) THEN
       NIntC=SIZE(IntCs%Def)
     ELSE
       NIntC=0 
       RETURN
     ENDIF
     !
     WRITE(*,*) TRIM(CHAR)
     WRITE(Out,*) TRIM(CHAR)
     WRITE(*,*) '         INTERNAL COORDINATES'
     WRITE(*,*) '       DEFINITION  ATOMS_INVOLVED         VALUE'//&
                '        CONSTRAINT    ACTIVE'
     WRITE(Out,*) 'INTERNAL COORDINATES'
     WRITE(Out,*) '       DEFINITION  ATOMS_INVOLVED      VALUE'//&
                  '        CONSTRAINT    ACTIVE'
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)=='STRE') THEN
         SUM=Value(I)/AngstromsToAu
         SUMConstr=IntCs%ConstrValue(I)/AngstromsToAu
       ELSE IF(IntCs%Def(I)(1:4)=='BEND'.OR. &
               IntCs%Def(I)(1:4)=='TORS'.OR. &
               IntCs%Def(I)(1:4)=='LINB'.OR. &
               IntCs%Def(I)(1:4)=='OUTP') THEN
         SUM=Value(I)*Conv
         SUMConstr=IntCs%ConstrValue(I)*Conv
       ELSE IF(IntCs%Def(I)(1:4)=='CART') THEN
         SUM=Value(I)
         SUMConstr=IntCs%ConstrValue(I)
       ENDIF
       WRITE(*,111) I,IntCs%Def(I),IntCs%Atoms(I,1:4),SUM, &
         IntCs%Constraint(I),SumConstr,IntCs%Active(I)
       WRITE(Out,111) I,IntCs%Def(I),IntCs%Atoms(I,1:4),SUM, &
         IntCs%Constraint(I),SumConstr,IntCs%Active(I)
     ENDDO
     !      
     111 FORMAT(I7,2X,A5,2X,4I5,2X,F12.6,L5,F12.6,L5)
     222 FORMAT(I7,2X,A5,2X,4I5,2X,3F12.6,L5,F12.6,L5)
     !
     !
   END SUBROUTINE PrtIntCoords
!
!----------------------------------------------------------
!
   SUBROUTINE RotationsOff(CartVect,XYZ,Print)
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
     INTEGER                     :: NCart,NatmsLoc,I,J,INFO
     LOGICAL                     :: Print
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
     !
     ! test orthogonalities
     !
     !      CALL TranslsOff(Rot1%D,.TRUE.)
     !      CALL TranslsOff(Rot2%D,.TRUE.)
     !      CALL TranslsOff(Rot3%D,.TRUE.)
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
     IF(SUM>1.D-6) THEN
       SUM1=SUM1*SUM1/SUM*100.D0
       SUM2=SUM2*SUM2/SUM*100.D0
       SUM3=SUM3*SUM3/SUM*100.D0
     ENDIF
     IF(Print) THEN
       WRITE(*,100) SUM1,SUM2,SUM3
     ENDIF
     100  FORMAT('Rot1= ',F7.3,'%    Rot2= ', &
                  F7.3,'%     Rot3= ',F7.3,'% ')
     CALL Delete(Vect)
     CALL Delete(Rot3)
     CALL Delete(Rot2)
     CALL Delete(Rot1)
     CALL Delete(CMCarts)
     CALL Delete(Theta2)
     CALL Delete(Theta)
   END SUBROUTINE RotationsOff
!
!----------------------------------------------------------
!
   SUBROUTINE TranslsOff(CartVect,Print)
     REAL(DOUBLE),DIMENSION(:) :: CartVect
     REAL(DOUBLE)              :: SUM,SUM1,SUM2,SUM3
     TYPE(DBL_VECT)            :: Tr1,Tr2,Tr3
     INTEGER                   :: I,J,NCart,NatmsLoc
     LOGICAL                   :: Print
     NCart=SIZE(CartVect)
     NatmsLoc=NCart/3
     CALL New(Tr1,NCart) 
     CALL New(Tr2,NCart) 
     CALL New(Tr3,NCart) 
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
     IF(Print) THEN
       WRITE(*,100) SUM1,SUM2,SUM3
     ENDIF
     100  FORMAT(' Tr1= ',F7.3,'%     Tr2= ',F7.3,'%      Tr3= ',&
                  F7.3,'% ')
     !
     CALL Delete(Tr3)
     CALL Delete(Tr2)
     CALL Delete(Tr1)
     !
   END SUBROUTINE TranslsOff
!
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
     INTEGER      :: I,J,Step,III
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
     IF(ABS(Sum1) < 1.D-6) THEN  !!! V1 & V2 are parallel
       Rot=Zero
       Sum=DOT_PRODUCT(V1,V2)
       DO I=1,3 ; Rot(I,I)=Sum ; ENDDO
     ELSE      
       Sum=One/SQRT(Sum1)
       CrossProd%D=SUM*CrossProd%D
       CosPhi=DOT_PRODUCT(V1,V2)     
       SinPhi=SQRT(SUM1)
       Sum=One-CosPhi
       DO III=1,2    
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
         IF(SumM>1.D-6) THEN
           SinPhi=-SinPhi
           Step=Step+1
           IF(Step > 2) THEN
             Call Halt('Rotation unsuccesful for torsvalue')
           ENDIF
         ELSE
           ! write(*,*) 'rotation succesful'
           EXIT    
         ENDIF
       ENDDO
     ENDIF
     !
100  CONTINUE
     CALL Delete(CrossProd)
     CALL Delete(Vect)
     !
   END SUBROUTINE Rotate
!
!-------------------------------------------------------
!
   SUBROUTINE MapAngleDispl(IntCs,NIntC,VectInt) 
     !
     ! This routine maps displacements of angles:
     ! all angle displacements must be between -PI and +PI.
     !
     TYPE(INTC)                :: IntCs
     INTEGER                   :: NIntC,I,II
     REAL(DOUBLE),DIMENSION(:) :: VectInt
     REAL(DOUBLE)              :: Sum,ASum,TwoPi
     !
     TwoPi=Two*PI
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)=='BEND'.OR. &
          IntCs%Def(I)(1:4)=='LINB'.OR. &
          IntCs%Def(I)(1:4)=='TORS'.OR. &
          IntCs%Def(I)(1:4)=='OUTP') THEN
         Sum=VectInt(I)
         IF(Sum>PI) Sum=Sum-TwoPi
         IF(Sum<-PI) Sum=TwoPi+Sum
         VectInt(I)=Sum
       ENDIF
     ENDDO
   END SUBROUTINE MapAngleDispl
!
!--------------------------------------------------------------------
!
   SUBROUTINE CartRNK1ToCartRNK2(VectCart,ActCarts,Add_O)
     REAL(DOUBLE),DIMENSION(:)   :: VectCart
     REAL(DOUBLE),DIMENSION(:,:) :: ActCarts(:,:)
     INTEGER :: I,J,NatmsLoc
     LOGICAL,OPTIONAL :: Add_O
     !
     NatmsLoc=SIZE(ActCarts,2)
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
   END SUBROUTINE CartRNK2ToCartRNK1
!
!----------------------------------------------------------------------
!
   SUBROUTINE ChkBendToLinB(IntCs,NIntC,XYZ,Top12, &
                            CtrlCoord,SCRPath,LongRangeTag)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER                     :: NIntC,Nintc_New
     TYPE(INT_VECT)              :: LinAtom,MarkLinb
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,Conv
     INTEGER                     :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
     INTEGER                     :: I,J,K,NatmsLoc
     TYPE(INT_RNK2)              :: LinBBridge,Top12
     TYPE(CoordCtrl)             :: CtrlCoord
     CHARACTER(LEN=*)            :: SCRPath
     CHARACTER(LEN=1)            :: LongRangeTag
     !
     ! Now check for bending -> linear bending transitions
     !
     NatmsLoc=SIZE(XYZ,2)
     CALL New(MarkLinB,NIntC)
     Conv=180.D0/PI
     !
     MarkLinB%I=0
          NLinB=0
     DO I=1,NIntC 
       IF(IntCs%Def(I)(1:4)=='BEND') THEN
         I1=IntCs%Atoms(I,1) 
         I2=IntCs%Atoms(I,2) 
         I3=IntCs%Atoms(I,3) 
         CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value)
         IF(ABS(Value-PI)*Conv < CtrlCoord%LinCrit) THEN  
           IF(Top12%I(I2,1)>2) THEN
             IntCs%Active(I)=.FALSE.
             CYCLE
           ENDIF
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
     CALL LongRangeTors(CtrlCoord,SCRPath,Top12,IntCs,NIntC, &
                        XYZ,MarkLinB,LongRangeTag)
     CALL Delete(MarkLinB)
     !
   END SUBROUTINE ChkBendToLinB    
!
!-------------------------------------------------------
!
   SUBROUTINE SetConstraint(IntCs,XYZ,Displ,LinCrit,NConstr, &
                            DoInternals)
     TYPE(INTC)     :: IntCs
     INTEGER        :: I,J,NIntC,LastIntcGeom,NDim,JJ,NConstr
     TYPE(DBL_VECT) :: Displ
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ  
     REAL(DOUBLE)   :: LinCrit
     LOGICAL        :: DoInternals
     !
     IF(NConstr==0)RETURN
     IF(AllocQ(IntCs%Alloc)) THEN
       NIntC=SIZE(IntCs%Def)
     ELSE
       NIntC=0
     ENDIF
     NDim=SIZE(Displ%D)
     IF(NDim/=NIntC.AND.DoInternals) &
         Call Halt('Dimensionality error in SetConstraint')
     !
     ! Get values of constraints 
     !
     IF(DoInternals) THEN
       CALL INTCValue(IntCs,XYZ,LinCrit)
       DO I=1,NIntC
         IF(IntCs%Constraint(I)) THEN
           Displ%D(I)=IntCs%ConstrValue(I)-IntCs%Value(I)
           !write(*,100) i,IntCs%ConstrValue(I),&
           !IntCs%Value(I),Displ%D(I)
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
           Displ%D(JJ)=IntCs%ConstrValue(I)-IntCs%Value(I)
         ENDIF
       ENDDO	
     ENDIF
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
   END SUBROUTINE ConstrConv
!
!-------------------------------------------------------------
!
   SUBROUTINE SetCartConstr(Carts,CartDispl,IntCs,NConstr)
     REAL(DOUBLE),DIMENSION(:) :: CartDispl,Carts
     INTEGER                   :: I,J,JJ,NIntC,NConstr
     TYPE(INTC)                :: IntCs
     !
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
           Carts(JJ)=IntCs%ConstrValue(I)
           CartDispl(JJ)=Zero
         ENDIF
       ENDDO
     ENDIF      
   END SUBROUTINE SetCartConstr
!
!---------------------------------------------------------------
!
   SUBROUTINE ScaleDispl(CartDispl,MaxCartDiff,DiffMax,RMSD)
     REAL(DOUBLE),DIMENSION(:) :: CartDispl 
     REAL(DOUBLE)              :: MaxCartDiff,DiffMax,RMSD,Sum
     INTEGER                   :: I,NCart
     !
     NCart=SIZE(CartDispl)
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
   END SUBROUTINE ScaleDispl
!
!-------------------------------------------------------
!
   SUBROUTINE BackTrfConvg(GConstr,GBackTrf, &
     DoIterate,DiffMax,RMSD,RMSDOld,ConstrMax, &
     ConstrRMS,ConstrRMSOld,IStep,RefreshAct)
     !
     REAL(DOUBLE) :: DiffMax,RMSD,RMSDOld
     REAL(DOUBLE) :: ConstrMax,ConstrRMS
     REAL(DOUBLE) :: ConstrRMSOld
     INTEGER      :: IStep
     LOGICAL      :: DoIterate,RefreshAct
     LOGICAL      :: ConvConstr
     TYPE(Constr) :: GConstr
     TYPE(BackTrf):: GBackTrf
     !
     DoIterate=.TRUE.
     ConvConstr=.TRUE.
     !
     IF(GConstr%NConstr/=0) THEN
       IF(IStep>1) THEN
         ConvConstr=(GConstr%DoLagr.OR. &
           ConstrMax<GConstr%ConstrMaxCrit.OR. &
           (ConstrRMS>ConstrRMSOld*GBackTrf%RMSCrit.AND.IStep>5)) 
       ELSE
         ConvConstr=.FALSE.
       ENDIF
     ENDIF
     !
     DoIterate=(DiffMax>GBackTrf%CooTrfCrit)
     IF(RMSD>RMSDOld*GBackTrf%RMSCrit) THEN 
       IF(DiffMax<GBackTrf%MaxCartDiff*GBackTrf%RMSCrit &
          .AND.IStep>10) THEN
         DoIterate=.FALSE.
       ELSE 
         RefreshAct=.TRUE.
       ENDIF
     ENDIF
     DoIterate=DoIterate.OR.(.NOT.ConvConstr)
     DoIterate=(DoIterate.AND.IStep<=GBackTrf%MaxIt_CooTrf)
     !
   END SUBROUTINE BackTrfConvg
!
!----------------------------------------------------------------
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
   END SUBROUTINE RotateMol
!
!-------------------------------------------------------------------
!
   SUBROUTINE CROSS_PRODUCT(V1,V2,CrossProd,Add_O)
     REAL(DOUBLE),DIMENSION(:) :: V1,V2,CrossProd
     LOGICAL,OPTIONAL          :: Add_O
     REAL(DOUBLE),DIMENSION(3) :: VectAux
     !
     VectAux(1)=V1(2)*V2(3)-V1(3)*V2(2)
     VectAux(2)=V1(3)*V2(1)-V1(1)*V2(3)
     VectAux(3)=V1(1)*V2(2)-V1(2)*V2(1)
     IF(PRESENT(Add_O)) THEN
       IF(Add_O) CrossProd=CrossProd+VectAux
     ELSE
       CrossProd=VectAux
     ENDIF
   END SUBROUTINE CROSS_PRODUCT
!
!------------------------------------------------------------------
!
   SUBROUTINE CALC_BxVect(B,VectInt,VectCart,Trp_O)
     TYPE(BMATR)               :: B
     REAL(DOUBLE),DIMENSION(:) :: VectInt,VectCart
     INTEGER                   :: NCart,NIntC,I,J,JJ,K,KK,LL
     REAL(DOUBLE)              :: Sum
     LOGICAL,OPTIONAL          :: Trp_O
     !
     NCart=SIZE(VectCart)
     NIntC=SIZE(VectInt)
     !
     IF(PRESENT(Trp_O)) THEN
       IF(Trp_O) THEN !!! Bt*VectInt
         VectCart=Zero
         DO I=1,NIntC
           Sum=VectInt(I)
           DO J=1,4
             K=B%IB(I,J)
             IF(K==0) CYCLE
               LL=3*(J-1)
               KK=3*(K-1)
             DO JJ=1,3
               LL=LL+1
               KK=KK+1
               VectCart(KK)=VectCart(KK)+Sum*B%B(I,LL)
             ENDDO
           ENDDO
         ENDDO 
       ENDIF
     ELSE
       VectInt=Zero
       DO I=1,NIntC
         Sum=Zero
         DO J=1,4
           K=B%IB(I,J)
           IF(K==0) CYCLE
             LL=3*(J-1)
             KK=3*(K-1)
           DO JJ=1,3
             LL=LL+1
             KK=KK+1
             Sum=Sum+B%B(I,LL)*VectCart(KK)
           ENDDO
         ENDDO
         VectInt(I)=Sum
       ENDDO
     ENDIF 
   END SUBROUTINE CALC_BxVect
!
!----------------------------------------------------------------------
!
   SUBROUTINE CALC_GcInvCartV(CholData,VectCartAux,VectCartAux2)
     TYPE(Cholesky)            :: CholData
     REAL(DOUBLE),DIMENSION(:) :: VectCArtAux,VectCartAux2
     TYPE(DBL_VECT)            :: VectCartAux3
     TYPE(DBL_VECT)            :: Vect1,Vect2,Vect3
     INTEGER                   :: NCart,NDim,NIntC,I
     !
     NCart=SIZE(VectCArtAux)
     CALL New(VectCartAux3,NCart)
     VectCartAux2=VectCartAux
     !
     CALL PermVect(VectCartAux2,VectCartAux3%D,CholData%Perm%I)
     CALL ScaleVect(VectCartAux3%D,CholData%GcScale%D)
     CALL CholFactSolve(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,CholData%ChDiag%D, &
       VectCartAux3%D,NCart,VectCartAux2)
     VectCartAux3%D=VectCartAux2
     CALL ScaleVect(VectCartAux3%D,CholData%GcScale%D)
     CALL PermVect(VectCartAux3%D,VectCartAux2,CholData%IPerm%I)
     CALL Delete(VectCartAux3)
   END SUBROUTINE CALC_GcInvCartV
!
!----------------------------------------------------------------------
!
   SUBROUTINE GetBMatInfo(SCRPath,NIntC,B,CholData)
     TYPE(BMATR)    :: B
     TYPE(Cholesky) :: CholData
     INTEGER        :: NIntC
     CHARACTER(LEN=*):: SCRPath
     !
     CALL ReadBMATR(B,TRIM(SCRPath)//'B')
     NIntC=SIZE(B%IB,1)
     CALL ReadChol(CholData,TRIM(SCRPath)//'CholFact')
   END SUBROUTINE GetBMatInfo
!
!-------------------------------------------------------------------
!
   SUBROUTINE RefreshBMatInfo(IntCs,XYZ,GTrfCtrl, &
                              GCoordCtrl,Print,SCRPath,DoNewRef_O)
     TYPE(INTC)                   :: IntCs
     TYPE(Cholesky)               :: CholData
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     INTEGER                      :: NatmsLoc,NCart,NIntC
     TYPE(BMATR)                  :: B
     INTEGER                      :: Print
     LOGICAL                      :: Print2,DoNewRef
     LOGICAL,OPTIONAL             :: DoNewRef_O
     CHARACTER(LEN=*)             :: SCRPath
     CHARACTER(LEN=DCL)           :: SCRPathWork
     TYPE(DBL_RNK2)               :: XYZWork
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def)
     Print2=(Print>=DEBUG_GEOP_MAX)
     IF(PRESENT(DoNewRef_O)) THEN
       DoNewRef=DoNewRef_O
     ELSE
       DoNewRef=.TRUE.
     ENDIF
     !
     CALL New(XYZWork,(/3,NatmsLoc/))
     XYZWork%D=XYZ
     SCRPathWork=SCRPath
     IF(.NOT.GTrfCtrl%DoClssTrf.AND.DoNewRef) THEN
       CALL CALC_XYZRot(XYZWork%D,IntCs,SCRPath,GTrfCtrl%ThreeAt,&
                        GTrfCtrl%Linearity,GTrfCtrl%TranslAt1, &
                        GTrfCtrl%RotAt2ToX,GTrfCtrl%RotAt3ToXY)
     ENDIF
     !
     ! Calculate B matrix in Atomic Units, 
     ! and compute Cholesky factor.
     !
     CALL BMatrix(XYZWork%D,NIntC,IntCs,B,GCoordCtrl%LinCrit, &
                  GTrfCtrl%DoClssTrf,GTrfCtrl%ThreeAt)
     CALL WriteBMATR(B,TRIM(SCRPathWork)//'B')
     !
     CALL CholFact(B,NCart,CholData,GTrfCtrl%DoClssTrf, &
                   Print2,GTrfCtrl%ThreeAt)
     CALL WriteChol(CholData,TRIM(SCRPathWork)//'CholFact')
     !
     CALL Delete(CholData)
     CALL Delete(B)
     CALL Delete(XYZWork)
   END SUBROUTINE RefreshBMatInfo
!
!---------------------------------------------------------------------
!
   SUBROUTINE DeleteBMatInfo(B,CholData)
     TYPE(BMATR)    :: B
     TYPE(Cholesky) :: CholData
     !
     CALL Delete(B)
     CALL Delete(CholData)
   END SUBROUTINE DeleteBMatInfo
!
!---------------------------------------------------------------------
!
   SUBROUTINE CALC_XYZRot(XYZ,IntCs,SCRPath,ThreeAt,Linearity, &
                          TranslAt1,RotAt2ToX,RotAt3ToXY)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     TYPE(INTC)                   :: IntCs
     REAL(DOUBLE),DIMENSION(3)    :: TranslAt1
     REAL(DOUBLE),DIMENSION(3,3)  :: RotAt2ToX,RotAt3ToXY
     TYPE(DBL_RNK2)  :: XYZRot
     TYPE(DBL_VECT)  :: Vect1,Vect2,Vect3
     TYPE(DBL_RNK2)  :: Rot
     TYPE(INT_RNK2)  :: Top12
     INTEGER         :: I,J,NatmsLoc,NMax12
     INTEGER         :: At1,At2,At3,ThreeAt(1:3)
     LOGICAL         :: Linearity
     CHARACTER(LEN=*):: SCRPath
     !
     ! Subroutine to calculate reference geometry for singularityless
     ! coordinate transformation
     !
     CALL ReadINT_RNK2(Top12,TRIM(SCRPath)//'Top12',I,J)
     NMax12=J-1
     NatmsLoc=SIZE(XYZ,2)
     CALL New(XYZRot,(/3,NatmsLoc/))
     CALL New(Vect1,3)
     CALL New(Vect2,3)
     CALL New(Vect3,3)
     CALL New(Rot,(/3,3/))
     XYZRot%D=XYZ
     !
     TranslAt1=Zero 
     RotAt2ToX=Zero
     RotAt3ToXY=Zero
     DO I=1,3 
       RotAt2ToX(I,I)=One
       RotAt3ToXY(I,I)=One
     ENDDO
     !
     ! First, find three atoms, which form a large triangle, 
     ! if molecule is not linear
     !
     CALL ThreeAtoms(XYZ,IntCs,Top12,ThreeAt,Linearity)
     At1=ThreeAt(1)
     At2=ThreeAt(2)
     At3=ThreeAt(3)
     !
     ! Place At1 to origin
     !
     Vect1%D=XYZ(:,At1)
     DO I=1,NatmsLoc
       XYZRot%D(:,I)=XYZ(:,I)-Vect1%D   
     ENDDO
     TranslAt1=Vect1%D
     !!!!CALL Put(Vect1,'Translation')
     !
     ! Place At2 onto X axis
     !
     Vect1%D=Zero
     Vect1%D(1)=One
     Vect2%D=XYZRot%D(:,At2)-XYZRot%D(:,At1) 
     CALL Rotate(Vect1%D,Vect2%D,Rot%D)
     !!!!CALL Put(Rot,'RotAt2')
     RotAt2ToX=Rot%D
     DO I=1,NatmsLoc
       Vect1%D=XYZRot%D(:,I)
       CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,Vect1%D,XYZRot%D(:,I))
     ENDDO
     !
     ! Linearity
     !
     IF(.NOT.Linearity) THEN         
       !
       ! Place At3 onto XY plane
       !
       Vect1%D=XYZRot%D(:,At1)-XYZRot%D(:,At3)
       Vect2%D=XYZRot%D(:,At2)-XYZRot%D(:,At3)
       CALL CROSS_PRODUCT(Vect1%D,Vect2%D,Vect3%D)
       Vect1%D=Zero
       Vect1%D(3)=One 
       CALL Rotate(Vect1%D,Vect3%D,Rot%D)
       RotAt3ToXY=Rot%D
       !!!!CALL Put(Rot,'RotAt3')
       DO I=1,NatmsLoc
         Vect1%D=XYZRot%D(:,I)
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,Vect1%D,XYZRot%D(:,I))
       ENDDO
       !!!!CALL Put(XYZRot,'XYZRot')
     ENDIF
     !
     XYZ=XYZRot%D
     CALL Delete(XYZRot)
     CALL Delete(Top12)
     CALL Delete(Rot)
     CALL Delete(Vect1)
     CALL Delete(Vect2)
     CALL Delete(Vect3)
   END SUBROUTINE CALC_XYZRot
!
!----------------------------------------------------------------------
!
   SUBROUTINE ThreeAtoms(XYZ,IntCs,Top12,ThreeAt,Linearity)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     TYPE(INTC)                   :: IntCs
     TYPE(DBL_VECT)               :: DistVect
     INTEGER                      :: I,J,NatmsLoc,NIntC
     INTEGER                      :: At1,At2,At3,ThreeAt(1:3)
     REAL(DOUBLE)                 :: D12,D13,D23,CX,CY,CZ
     REAL(DOUBLE)                 :: Dist,Sum,Dist12,Area2
     LOGICAL                      :: Linearity
     TYPE(DBL_RNK2)               :: XYZRot
     TYPE(INT_RNK2)               :: Top12 
     !
     NatmsLoc=SIZE(XYZ,2)
     NIntC=SIZE(IntCs%Def)
     CALL New(XYZRot,(/3,NatmsLoc/))
     CALL New(DistVect,NatmsLoc)
     DistVect%D=0
     XYZRot%D=XYZ
     !
     ! First, try to select three atoms from constraints
     !
     At1=0
     At2=0
     At3=0
     CALL ThreeConstr(IntCs,Top12,NatmsLoc,At1,At2,At3)
     !
     ! Find three atoms, which form a large triangle, 
     ! if molecule is not linear
     !
     CALL CenterOfMass(CX,CY,CZ,XYZ_O=XYZRot%D,Move_O=.True.)
     !
     ! Find farthest point from COM
     !
     IF(At1==0) THEN
       Dist=-One
       DO I=1,NatmsLoc
         Sum=XYZRot%D(1,I)**2+XYZRot%D(2,I)**2+XYZRot%D(3,I)**2
         Sum=Sum*Top12%I(I,1)**2
         IF(Sum>Dist) THEN
           Dist=Sum
           At1=I
         ENDIF
       ENDDO
     ENDIF
     !
     ! Find farthest point from At1
     !
     DO I=1,NatmsLoc
       DistVect%D(I)=SQRT((XYZRot%D(1,I)-XYZRot%D(1,At1))**2+&
                          (XYZRot%D(2,I)-XYZRot%D(2,At1))**2+&
                          (XYZRot%D(3,I)-XYZRot%D(3,At1))**2)
     ENDDO
     IF(At2==0) THEN
       Dist=-One
       DO I=1,NatmsLoc
         Sum=Top12%I(I,1)*DistVect%D(I)
         IF(Sum>Dist) THEN
           Dist=Sum
           At2=I
         ENDIF
       ENDDO
     ENDIF
     !
     ! Find point farthest from At1 and At2
     !
     IF(At3==0) THEN
       D12=DistVect%D(At2)
       Dist=-1.D99
       DO I=1,NatmsLoc
         IF(I==At1.OR.I==At2) CYCLE
         D23=SQRT((XYZRot%D(1,I)-XYZRot%D(1,At2))**2+&
                  (XYZRot%D(2,I)-XYZRot%D(2,At2))**2+&
                  (XYZRot%D(3,I)-XYZRot%D(3,At2))**2)
         D13=DistVect%D(I)
         Sum=0.5D0*(D12+D13+D23)
         Area2=Sum*(Sum-D12)*(Sum-D13)*(Sum-D23)
!write(*,*) 'area2= ',at1,at2,i,Area2,Top12%I(I,1)*Area2,dist
         Area2=Top12%I(I,1)*Area2
         IF(Area2>Dist) THEN
           Dist=Area2
           At3=I
!write(*,*) 'found bigger= ',at3
         ENDIF
       ENDDO
     ENDIF
     !
     ! Check linearity
     !
     Linearity=.FALSE.
     IF(NatmsLoc==2) THEN
       At3=0
       Linearity=.TRUE.
     ELSE
       D12=SQRT(DOT_PRODUCT(XYZRot%D(:,At1)-XYZRot%D(:,At2), &
           XYZRot%D(:,At1)-XYZRot%D(:,At2)))
       D13=SQRT(DOT_PRODUCT(XYZRot%D(:,At1)-XYZRot%D(:,At3), &
           XYZRot%D(:,At1)-XYZRot%D(:,At3)))
       D23=SQRT(DOT_PRODUCT(XYZRot%D(:,At2)-XYZRot%D(:,At3), &
           XYZRot%D(:,At2)-XYZRot%D(:,At3)))
       IF(ABS(D12-(D13+D23))<0.001D0) Linearity=.TRUE.
     ENDIF
     ThreeAt(1)=At1
     ThreeAt(2)=At2
     ThreeAt(3)=At3
     ! Tidy up
     CALL Delete(XYZRot)
     CALL Delete(DistVect)
   END SUBROUTINE ThreeAtoms
!
!--------------------------------------------------------------------
!
   SUBROUTINE WriteIntCs(IntCs,FileName)
     TYPE(INTC)                     :: IntCs
     CHARACTER(LEN=*)               :: FileName
     INTEGER                        :: Dim1,Dim2
     !
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     IF(.NOT.AllocQ(IntCs%Alloc)) THEN
       Dim1=0
     ELSE
       Dim1=SIZE(IntCs%Def)
     ENDIF
     WRITE(99) Dim1
     IF(Dim1==0) RETURN
     WRITE(99) IntCs%Def
     WRITE(99) IntCs%FCType
     WRITE(99) IntCs%Atoms
     WRITE(99) IntCs%Value
     WRITE(99) IntCs%Constraint
     WRITE(99) IntCs%ConstrValue
     WRITE(99) IntCs%Active
     CLOSE(99,STATUS='KEEP')
   END SUBROUTINE WriteIntCs
!
!---------------------------------------------------------------------
!
   SUBROUTINE ReadIntCs(IntCs,FileName)
     TYPE(INTC)       :: IntCs
     CHARACTER(LEN=*) :: FileName   
     INTEGER          :: Dim1,Dim2
     LOGICAL          :: Exists
     !
     INQUIRE(File=FileName,EXIST=Exists)
     IF(.NOT.Exists) &
       CALL Halt('File does not exist error in ReadIntCs')
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     READ(99) Dim1
     CALL New(IntCs,Dim1)
     IF(Dim1==0) RETURN
     READ(99) IntCs%Def
     READ(99) IntCs%FcType
     READ(99) IntCs%Atoms
     READ(99) IntCs%Value
     READ(99) IntCs%Constraint
     READ(99) IntCs%ConstrValue
     READ(99) IntCs%Active
     CLOSE(99,STATUS='KEEP')
   END SUBROUTINE ReadIntCs
!
!--------------------------------------------------------------------
!
   SUBROUTINE WriteChol(CholData,FileName)
     TYPE(Cholesky)   :: CholData
     CHARACTER(LEN=*) :: FileName
     INTEGER          :: ChNon0,NCart
     !
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     IF(.NOT.AllocQ(CholData%GcScale%Alloc)) THEN
       RETURN
     ENDIF
     NCart=SIZE(CholData%Perm%I)
     ChNon0=CholData%ChRowPt%I(NCart+1)-1
     WRITE(99) NCart 
     WRITE(99) ChNon0
     WRITE(99) CholData%Perm%I
     WRITE(99) CholData%IPerm%I
     WRITE(99) CholData%ChRowPt%I
     WRITE(99) CholData%ChColPt%I
     WRITE(99) CholData%GcScale%D
     WRITE(99) CholData%ChDiag%D
     WRITE(99) CholData%ChFact%D
     CLOSE(99,STATUS='KEEP')
   END SUBROUTINE WriteChol
!
!---------------------------------------------------------------
!
   SUBROUTINE ReadChol(CholData,FileName)
     TYPE(Cholesky)   :: CholData
     CHARACTER(LEN=*) :: FileName
     INTEGER          :: ChNon0,NCart
     LOGICAL          :: Exists
     !
     INQUIRE(File=FileName,EXIST=Exists)
     IF(.NOT.Exists) &
       CALL Halt('File does not exist error in ReadChol')
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     READ(99) NCart 
     READ(99) ChNon0
     CALL New_Chol(CholData,NCart,ChNon0)
     READ(99) CholData%Perm%I 
     READ(99) CholData%IPerm%I 
     READ(99) CholData%ChRowPt%I
     READ(99) CholData%ChColPt%I
     READ(99) CholData%GcScale%D
     READ(99) CholData%ChDiag%D
     READ(99) CholData%ChFact%D
     CLOSE(99,STATUS='KEEP')
   END SUBROUTINE ReadChol
!
!--------------------------------------------------------------------
!
   SUBROUTINE CleanBConstr(IntCs,B,NatmsLoc)
     TYPE(BMATR) :: B
     TYPE(INTC)  :: IntCs
     INTEGER     :: I,J,K,L,NIntC,JJ,LL,NatmsLoc,NCart
     TYPE(INT_VECT):: IConstr
     !
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def)
     CALL New(IConstr,NCart)
     IConstr%I=1  
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         J=3*(B%IB(I,1)-1)
         IF(IntCs%Def(I)(1:5)=='CARTX') IConstr%I(J+1)=0
         IF(IntCs%Def(I)(1:5)=='CARTY') IConstr%I(J+2)=0
         IF(IntCs%Def(I)(1:5)=='CARTZ') IConstr%I(J+3)=0
       ENDIF
     ENDDO 
     !
     DO I=1,NIntC
       DO J=1,4
         JJ=B%IB(I,J)
         IF(JJ==0) EXIT
         JJ=3*(JJ-1)
         LL=3*(J-1)
         DO L=1,3
           IF(IConstr%I(JJ+L)==0) B%B(I,LL+L)=Zero
         ENDDO
       ENDDO
     ENDDO 
     CALL Delete(IConstr)
   END SUBROUTINE CleanBConstr
!
!--------------------------------------------------------------------
!
   SUBROUTINE CleanB(ThreeAt,B,NIntC)
     TYPE(BMATR) :: B
     INTEGER     :: I,J,K,L,NIntC,J1,J2,ThreeAt(1:3)
     !
     DO I=1,NIntC
       DO J=1,4
         IF(B%IB(I,J)==ThreeAt(1)) THEN
           J2=3*J
           J1=J2-2
           B%B(I,J1:J2)=Zero
         ENDIF
         IF(B%IB(I,J)==ThreeAt(2)) THEN
           J2=3*J
           J1=J2-1
           B%B(I,J1:J2)=Zero
         ENDIF
         IF(B%IB(I,J)==ThreeAt(3)) THEN
           J2=3*J
           B%B(I,J2)=Zero
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE CleanB
!
!--------------------------------------------------------------------
!
   SUBROUTINE WriteINT_RNK2(Top,FileName)
     INTEGER,DIMENSION(:,:) :: Top
     CHARACTER(LEN=*)       :: FileName
     INTEGER                :: NDim1,NDim2
     !
     NDim1=SIZE(Top,1)
     NDim2=SIZE(Top,2)
     !
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     WRITE(99) NDim1,NDim2 
     WRITE(99) Top
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE WriteINT_RNK2
!
!-------------------------------------------------------------------
!
   SUBROUTINE ReadINT_RNK2(Top,FileName,NDim1,NDim2)
     TYPE(INT_RNK2)         :: Top
     CHARACTER(LEN=*)       :: FileName
     INTEGER                :: NDim1,NDim2
     LOGICAL                :: Exists
     !
     INQUIRE(File=FileName,EXIST=Exists)
     IF(.NOT.Exists) &
       CALL Halt('File does not exist error in ReadINT_RNK2 for '// &
                  FileName)
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND(99)
       READ(99) NDim1,NDim2 
       CALL New(Top,(/NDim1,NDim2/))
       READ(99) Top%I
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE ReadINT_RNK2
!
!-------------------------------------------------------------------
!
   SUBROUTINE PrtXYZ(AtNum,XYZ,FileName,Title)
     REAL(DOUBLE),DIMENSION(:)     :: AtNum
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ
     CHARACTER(LEN=*)              :: FileName
     CHARACTER(LEN=*)              :: Title
     CHARACTER(LEN=1)              :: Char 
     INTEGER                       :: I,NatmsLoc
     !
     NatmsLoc=SIZE(XYZ,2)
     OPEN(File=FileName,Unit=99,FORM='FORMATTED',STATUS='UNKNOWN')
     DO 
       READ(99,33,END=1) Char
     ENDDO
     33 format(a1)
     1    CONTINUE
     WRITE(99,*) NatmsLoc 
     WRITE(99,*) Title 
     DO I=1,NatmsLoc
       WRITE(99,100) INT(AtNum(I)),XYZ(1:3,I)/AngstromsToAu
     ENDDO 
     100  FORMAT(I4,2X,3F20.10)
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE PrtXYZ
!
!---------------------------------------------------------------------
!
   SUBROUTINE BoxBorders(XYZ,BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: BXMIN,BXMax
     REAL(DOUBLE)                :: BYMIN,BYMax
     REAL(DOUBLE)                :: BZMIN,BZMax
     REAL(DOUBLE)                :: VBig
     INTEGER                     :: NatmsLoc,I
     !
     !find borders of the global Box
     !
     NatmsLoc=SIZE(XYZ,2)
     VBig=1.D99
     BXMIN= VBIG
     BXMax=-VBIG
     BYMIN= VBIG
     BYMax=-VBIG
     BZMIN= VBIG
     BZMax=-VBIG
     DO I=1,NatmsLoc
       IF(XYZ(1,I)<BXMIN) BXMIN=XYZ(1,I)
       IF(XYZ(1,I)>BXMax) BXMax=XYZ(1,I)
       IF(XYZ(2,I)<BYMIN) BYMIN=XYZ(2,I)
       IF(XYZ(2,I)>BYMax) BYMax=XYZ(2,I)
       IF(XYZ(3,I)<BZMIN) BZMIN=XYZ(3,I)
       IF(XYZ(3,I)>BZMax) BZMax=XYZ(3,I)
     ENDDO
   END SUBROUTINE BoxBorders
! 
!-------------------------------------------------------------------
!
   SUBROUTINE SetMaxBlkS(Mat)
     TYPE(BCSR) :: Mat
     INTEGER    :: I,II
     ! Set MaxNon0-s for overlap matrices
     II=0
     DO I=1,Mat%NAtms
       II=II+(Mat%RowPt%I(I+1)-Mat%RowPt%I(I))**2
     ENDDO
     MaxBlkS=II+1
   END SUBROUTINE SetMaxBlkS
! 
!-------------------------------------------------------------------
!
   SUBROUTINE ReadBMATR(B,FileName)
     TYPE(BMATR) :: B
     INTEGER     :: NIntC
     LOGICAL     :: Exists
     CHARACTER(LEN=*) :: FileName
     !
     INQUIRE(File=FileName,EXIST=Exists)
     IF(.NOT.Exists) &
       CALL Halt('File does not exist error in ReadBMATR for '// &
                  FileName)
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     READ(99) NIntC
     CALL New(B,NIntC)
     READ(99) B%IB
     READ(99) B%B
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE ReadBMATR
! 
!-------------------------------------------------------------------
!
   SUBROUTINE WriteBMATR(B,FileName)
     TYPE(BMATR) :: B
     INTEGER     :: NIntC
     CHARACTER(LEN=*) :: FileName
     !
     NIntC=SIZE(B%IB,1)
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     WRITE(99) NIntC
     WRITE(99) B%IB
     WRITE(99) B%B
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE WriteBMATR
! 
!-------------------------------------------------------------------
!
   SUBROUTINE LongRangeTors(CtrlCoord,SCRPath,Top12, &
                            IntCs,NIntC,XYZ,MarkLinB,LongRangeTag)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER                     :: NIntC,NIntc_New
     TYPE(INT_VECT)              :: LinAtom,MarkLinb
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,Conv,Dist12
     INTEGER                     :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
     INTEGER                     :: I,J,K,NatmsLoc,III
     TYPE(INT_RNK2)              :: LinBBridge,Top12
     TYPE(CoordCtrl)             :: CtrlCoord
     CHARACTER(LEN=*)            :: SCRPath
     CHARACTER(LEN=1)            :: LongRangeTag
     LOGICAL                     :: RepeatChk,DoVdW
     !
     Conv=180.D0/PI
     NIntC=SIZE(IntCs%Def)
     !
     NMax12=SIZE(Top12%I,2)-1
     !
     CALL NEW(LinBBridge,(/2,NIntc/))
     CALL NEW(LinAtom,NIntc)
     LinBBridge%I(1:2,:)=0
     LinAtom%I(:)=0
     NLinB=0
     DO I=1,NIntc
       IF(IntCs%Def(I)(1:5)=='LINB1'.AND. &
         LinAtom%I(IntCs%Atoms(I,1))==0) THEN 
         NLinB=NLinB+1 
         LinAtom%I(IntCs%Atoms(I,1))=1
         LinAtom%I(IntCs%Atoms(I,2))=1
         LinAtom%I(IntCs%Atoms(I,3))=1
         LinBBridge%I(1,NLinB)=IntCs%Atoms(I,1)
         LinBBridge%I(2,NLinB)=IntCs%Atoms(I,2)
         ! now, go on left side, use Top12, then 
         ! go on right, then define torsions
         I1=LinBBridge%I(2,NLinB)
         I2=LinBBridge%I(1,NLinB)
         DO III=1,NIntC
           RepeatChk=.FALSE.
           DO J=1,Top12%I(I2,1)
             I3=Top12%I(I2,1+J) 
             CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value)
             IF(ABS(Value-PI)*Conv<CtrlCoord%LinCrit) THEN 
               LinBBridge%I(1,NLinB)=I3
               LinAtom%I(I3)=1
               I1=I2
               I2=I3
               RepeatChk=.TRUE.
               EXIT
             ENDIF
           ENDDO
           IF(.NOT.RepeatChk) EXIT
         ENDDO
         !
         I1=LinBBridge%I(1,NLinB)
         I2=LinBBridge%I(2,NLinB)
         DO III=1,NIntC
           RepeatChk=.FALSE.
           DO J=1,Top12%I(I2,1)
             I3=Top12%I(I2,1+J) 
             CALL BENDValue(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value)
             IF(ABS(Value-PI)*Conv<CtrlCoord%LinCrit) THEN 
               LinBBridge%I(2,NLinB)=I3
               LinAtom%I(I3)=1
               I1=I2
               I2=I3
               RepeatChk=.TRUE.
               EXIT
             ENDIF
           ENDDO
           IF(.NOT.RepeatChk) EXIT
         ENDDO
       ENDIF
       !
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
                 IntC_New%FCType(NIntC)=LongRangeTag
                 IntC_New%Atoms(NIntC,1)=I1
                 IntC_New%Atoms(NIntC,2)=I2
                 IntC_New%Atoms(NIntC,3)=I3
                 IntC_New%Atoms(NIntC,4)=I4
                 IntC_New%Value(NIntC)=Zero
                 IntC_New%Constraint(NIntC)=.FALSE.
                 IntC_New%ConstrValue(NIntC)=Zero   
                 IntC_New%Active(NIntC)=.TRUE. 
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
     !
   END SUBROUTINE LongRangeTors
!
!----------------------------------------------------------------------
!
   SUBROUTINE RedundancyOff(Displ,SCRPath,Print)
     REAL(DOUBLE),DIMENSION(:) :: Displ
     REAL(DOUBLE)              :: Perc 
     TYPE(BMATR)               :: B
     TYPE(Cholesky)            :: CholData
     INTEGER                   :: NIntC,NCart
     TYPE(DBL_VECT)            :: Vect1,Vect2,Displ2
     TYPE(DBL_RNK2)            :: FullGcInv
     INTEGER                   :: Print
     CHARACTER(LEN=*)          :: SCRPath
     !
     CALL GetBMatInfo(SCRPath,NIntC,B,CholData)
     IF(NIntC/=SIZE(Displ)) &
       CALL Halt('Dimension error in RedundancyOff')
     NCart=SIZE(CholData%IPerm%I)
     CALL New(Vect1,NCart)
     CALL New(Vect2,NCart)
     CALL New(Displ2,NIntC)
     Displ2%D=Displ
     !
     CALL CALC_BxVect(B,Displ,Vect1%D,Trp_O=.TRUE.)
     !
     CALL New(FullGcInv,(/NCart,NCart/))
     CALL GetFullGcInv(B,FullGcInv,NIntC,NCart)
     CALL DGEMM_NNc(NCart,NCart,1,One,Zero,FullGcInv%D,Vect1%D,Vect2%D)
     CALL Delete(FullGcInv)
     !
    !CALL CALC_GcInvCartV(CholData,Vect1%D,Vect2%D)
     !
     CALL CALC_BxVect(B,Displ,Vect2%D)
     !
     Perc=DOT_PRODUCT(Displ,Displ2%D)/DOT_PRODUCT(Displ2%D,Displ2%D) 
     Perc=(One-ABS(Perc))*100.D0
     IF(Print>=DEBUG_GEOP_MAX) THEN
       WRITE(*,100) Perc
       WRITE(Out,100) Perc
     ENDIF
100  FORMAT("Percentage of Redundancy projected out= ",F8.2)
     !
     CALL Delete(Displ2)
     CALL Delete(Vect1)
     CALL Delete(Vect2)
     CALL Delete(B)
     CALL Delete(CholData)
     ! 
   END SUBROUTINE RedundancyOff
! 
!-------------------------------------------------------------------
!
   SUBROUTINE PrimToDeloc(VectInt,VectCart,B,CholData)
     REAL(DOUBLE),DIMENSION(:) :: VectInt,VectCart
     TYPE(BMATR)               :: B
     TYPE(Cholesky)            :: CholData
     INTEGER                   :: I,J,NCart,NIntC
     TYPE(DBL_VECT)            :: VectCartAux
     !
     NCart=SIZE(VectCart)
     NIntC=SIZE(VectInt)
     CALL New(VectCartAux,NCart)
     !
     CALL CALC_BxVect(B,VectInt,VectCartAux%D,Trp_O=.TRUE.)
     !
     CALL PermVect(VectCartAux%D,VectCart,CholData%Perm%I)
     CALL ScaleVect(VectCart,CholData%GcScale%D)
     !
     CALL UtInvX(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,VectCart,NCart,VectCartAux%D)
     !
     DO I=1,NCart
       VectCartAux%D(I)=SQRT(CholData%ChDiag%D(I))*VectCartAux%D(I)
     ENDDO
     !
     VectCart=VectCartAux%D
     !
     CALL Delete(VectCartAux)
   END SUBROUTINE PrimToDeloc
! 
!-------------------------------------------------------------------
!
   SUBROUTINE DelocToPrim(NewDelocs,NewPrims,B,CholData)
     REAL(DOUBLE),DIMENSION(:) :: NewDelocs,NewPrims
     INTEGER                   :: NCart,NIntC
     TYPE(Cholesky)            :: CholData
     TYPE(BMATR)               :: B
     TYPE(DBL_VECT)            :: VectCart,VectInt
     INTEGER                   :: I
     !
     NCart=SIZE(NewDelocs)
     NIntC=SIZE(NewPrims)
     CALL New(VectCart,NCart)
     !
     DO I=1,NCart
       NewDelocs(I)=SQRT(CholData%ChDiag%D(I))*NewDelocs(I)
     ENDDO
     !
     CALL UInvX(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,NewDelocs,NCart,VectCart%D)
     !
     CALL ScaleVect(VectCart%D,CholData%GcScale%D)
     CALL PermVect(VectCart%D,NewDelocs,CholData%IPerm%I)
     CALL CALC_BxVect(B,NewPrims,NewDelocs)
     !
     CALL Delete(VectCart)
   END SUBROUTINE DelocToPrim
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasAngle(Char)
     CHARACTER(LEN=*) :: Char
     HasAngle=(Char(1:4)=='BEND'.OR. &
               Char(1:4)=='LINB'.OR. &
               Char(1:4)=='OUTP'.OR. &
               Char(1:4)=='TORS')
   END FUNCTION HasAngle
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasBendLinB(Char)
     CHARACTER(LEN=*) :: Char
     HasBendLinB=(Char(1:4)=='BEND'.OR. &
                  Char(1:4)=='LINB')
   END FUNCTION HasBendLinB
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasTorsOutP(Char)
     CHARACTER(LEN=*) :: Char
     HasTorsOutP=(Char(1:4)=='OUTP'.OR. &
                  Char(1:4)=='TORS')
   END FUNCTION HasTorsOutP
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrInvHess(IntCs,SCRPath,GHessian, &
                          LagrMult,NCart,NConstr, &
                          IHessM,JHessM,AHessM)
     TYPE(BMATR)      :: B
     TYPE(Hessian)    :: GHessian
     TYPE(INTC)       :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: LagrMult
     INTEGER          :: I,J,NIntC,NCart,NConstr
     CHARACTER(LEN=*) :: SCRPath 
     TYPE(INT_VECT)   :: IHessXX,JHessXX,IHessXL,JHessXL
     TYPE(INT_VECT)   :: IHessM,JHessM
     TYPE(DBL_VECT)   :: AHessXX,AHessXL,AHessM
     !
     NIntC=SIZE(IntCs%Def)
     !
     CALL ReadBMATR(B,TRIM(SCRPath)//'B')
     IF(SIZE(B%IB,1)/=NIntC) &
       CALL Halt('Dimension error in LagrInvHess.')
     CALL GetHessXX(IntCs,B,GHessian,LagrMult, &
                    IHessXX,JHessXX,AHessXX,NCart,SCRPath_O=SCRPath)
     CALL GetHessXL(IntCs,B,IHessXL,JHessXL,AHessXL, &
                    NCart,NConstr,SCRPath_O=SCRPath)
     CALL MergeXLXX(IHessXL,JHessXL,AHessXL, &
                    IHessXX,JHessXX,AHessXX, &
                    IHessM,JHessM,AHessM,SCRPath_O=SCRPath)
     !
     CALL Delete(IHessXX)
     CALL Delete(JHessXX)
     CALL Delete(AHessXX)
     CALL Delete(IHessXL)
     CALL Delete(JHessXL)
     CALL Delete(AHessXL)
     CALL Delete(B)
   END SUBROUTINE LagrInvHess
!
!------------------------------------------------------------------
!
   SUBROUTINE GetHessXX(IntCs,B,GHessian,LagrMult, &
                        IHessXX,JHessXX,AHessXX,NCart,SCRPath_O)
     TYPE(INTC)       :: IntCs
     INTEGER          :: I,J,NIntC,NCart,NZSpB,NConstraint
     REAL(DOUBLE),DIMENSION(:) :: LagrMult
     TYPE(Hessian)    :: GHessian
     TYPE(BMATR)      :: B,BSc
     REAL(DOUBLE)     :: BScale
     TYPE(INT_VECT)   :: ISpB,JSpB,ISpBSc,JSpBSc,IHessXX,JHessXX
     TYPE(INT_VECT)   :: ISpBt,JSpBt
     TYPE(DBL_VECT)   :: ASpB,ASpBSc,AHessXX
     TYPE(DBL_VECT)   :: ASpBt
     CHARACTER(LEN=*) :: SCRPath_O
     !
     CALL Set_BMATR_EQ_BMATR(BSc,B)
     NIntC=SIZE(IntCs%Def)
     !
     NConstraint=0
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)=='STRE') THEN
         BScale=GHessian%Stre
       ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
         BScale=GHessian%Bend
       ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
         BScale=GHessian%Tors
       ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
         BScale=GHessian%OutP
       ELSE IF(IntCs%Def(I)(1:4)=='LINB') THEN
         BScale=GHessian%LinB
       ELSE IF(IntCs%Def(I)(1:4)=='CART') THEN
         BScale=Zero
       ELSE
         BScale=One
       ENDIF
       IF(IntCs%Constraint(I)) THEN
         NConstraint=NConstraint+1
         !IF(IntCs%Def(I)(1:4)/='CART') THEN
         !  BScale=BScale-LagrMult(NConstraint)
         !ENDIF
       ENDIF
       DO J=1,12 ; BSc%B(I,J)=BScale*BSc%B(I,J) ; ENDDO
     ENDDO
     !
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     CALL BtoSpB_1x1(BSc,ISpBSc,JSpBSc,ASpBSc)
     !
     CALL New(ISpBt,NCart+1)
     NZSpB=ISpB%I(NIntC+1)-1
     CALL New(JSpBt,NZSpB)
     CALL New(ASpBt,NZSpB)
     CALL TransPose1x1(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart, &
          ISpBt%I,JSpBt%I,ASpBt%D,'full')
     !
     CALL MatMul_1x1(ISpBt%I,JSpBt%I,ASpBt%D, &
                     ISpBSc%I,JSpBSc%I,ASpBSc%D, &
                     IHessXX,JHessXX,AHessXX,NCart,NIntC,NCart)
     !
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(ISpBt)
     CALL Delete(JSpBt)
     CALL Delete(ASpBt)
     CALL Delete(ISpBSc)
     CALL Delete(JSpBSc)
     CALL Delete(ASpBSc)
     CALL Delete(BSc)
   END SUBROUTINE GetHessXX
!
!------------------------------------------------------------------
!
   SUBROUTINE GetHessXL(IntCs,B,IHessXLs,JHessXLs,AHessXLs, &
                        NCart,NConstr,SCRPath_O)
     TYPE(INTC)       :: IntCs
     TYPE(BMATR)      :: B
     TYPE(INT_VECT)   :: IHessXL,JHessXL,IHessXLt,JHessXLt
     TYPE(INT_VECT)   :: IHessXLs,JHessXLs
     TYPE(DBL_VECT)   :: AHessXL,AHessXLt,AHessXLs
     INTEGER          :: NCart,NConstr,NDim,NNew,NZ,I,J,K,JJ,KK,NIntC
     CHARACTER(LEN=*),OPTIONAL :: SCRPath_O
     !
     NDim=NCart+NConstr 
     NIntC=SIZE(IntCs%Def)
     CALL New(IHessXL,NDim+1)
     IHessXL%I(1:NCart+1)=1
     NNew=NCart
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NNew=NNew+1
         IF(NNew-NCart>NConstr) CALL Halt('Dimension error in GetHessXL') 
         NZ=0
         DO J=1,4 
           JJ=B%IB(I,J)
           IF(JJ==0) EXIT
           NZ=NZ+3
         ENDDO
         IHessXL%I(NNew+1)=IHessXL%I(NNew)+NZ
       ENDIF
     ENDDO
     !
     NZ=IHessXL%I(NDim+1)-1
     CALL New(JHessXL,NZ)
     CALL New(AHessXL,NZ)
     !
     NZ=0
     NNew=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NNew=NNew+1
         IF(NNew>NConstr)CALL Halt('Dimension error #2 in GetHessXL') 
         DO J=1,4 
           JJ=B%IB(I,J)
           IF(JJ==0) EXIT
           JJ=3*(JJ-1)
           DO K=1,3 
             NZ=NZ+1
             KK=3*(J-1)+K
             JHessXL%I(NZ)=JJ+K
             AHessXL%D(NZ)=-B%B(I,KK)
           ENDDO
         ENDDO
       ENDIF
     ENDDO
     !
     ! Build transpose of XL and merge XL with it.
     !
     CALL New(IHessXLt,NDim+1)
     CALL New(JHessXLt,NZ)
     CALL New(AHessXLt,NZ)
     CALL TransPose1x1(IHessXL%I,JHessXL%I,AHessXL%D, &
       NDim,NDim,IHessXLt%I,JHessXLt%I,AHessXLt%D,'full')
     CALL AddMat_1x1(IHessXL%I,JHessXL%I,AHessXL%D, &
                     IHessXLt%I,JHessXLt%I,AHessXLt%D, &
                     IHessXLs,JHessXLs,AHessXLs,NDim,NDim)
     !IF(PRESENT(SCRPath_O)) THEN
     !  CALL Plot_1x1(IHessXL%I,JHessXL%I,TRIM(SCRPath_O)//'XL',NDim)
     !  CALL Plot_1x1(IHessXLt%I,JHessXLt%I,TRIM(SCRPath_O)//'XLt',NDim)
     !  CALL Plot_1x1(IHessXLs%I,JHessXLs%I,TRIM(SCRPath_O)//'XLs',NDim)
     !ENDIF
     CALL Delete(IHessXLt)
     CALL Delete(JHessXLt)
     CALL Delete(AHessXLt)
     CALL Delete(IHessXL)
     CALL Delete(JHessXL)
     CALL Delete(AHessXL)
   END SUBROUTINE GetHessXL
!
!------------------------------------------------------------------
!
   SUBROUTINE MergeXLXX(IHessXL,JHessXL,AHessXL, &
                        IHessXX,JHessXX,AHessXX, &
                        IHessM,JHessM,AHessM,SCRPath_O)
     TYPE(INT_VECT)  :: IHessXL,JHessXL,IHessXX,JHessXX,IHessP
     TYPE(INT_VECT)  :: IHessM,JHessM
     TYPE(DBL_VECT)  :: AHessXL,AHessXX,AHessM
     INTEGER         :: NDimXL,NZXL,NDimXX,NZXX
     CHARACTER(LEN=*),OPTIONAL :: SCRPath_O
     !
     NDimXL=SIZE(IHessXL%I)-1
     NZXL=IHessXL%I(NDimXL+1)-1 
     NDimXX=SIZE(IHessXX%I)-1
     NZXX=IHessXX%I(NDimXX+1)-1 
     IF(NDimXX>NDimXL) CALL Halt('Dimension error in MergeXLXX')
     CALL New(IHessP,NDimXL+1)  
     IHessP%I(1:NDimXX+1)=IHessXX%I(1:NDimXX+1) 
     IHessP%I((NDimXX+2):(NDimXL+1))=IHessP%I(NDimXX+1)
     !
     CALL AddMat_1x1(IHessP%I,JHessXX%I,AHessXX%D, &
                     IHessXL%I,JHessXL%I,AHessXL%D, &
                     IHessM,JHessM,AHessM,NDimXL,NDimXL)
     !IF(PRESENT(SCRPath_O)) THEN
     !  CALL Plot_1x1(IHessP%I,JHessXX%I,TRIM(SCRPath_O)//'XXP',NDimXL)
     !  CALL Plot_1x1(IHessXL%I,JHessXL%I,TRIM(SCRPath_O)//'XL',NDimXL)
     !  CALL Plot_1x1(IHessM%I,JHessM%I,TRIM(SCRPath_O)//'M',NDimXL)
     !ENDIF
     CALL Delete(IHessP)
   END SUBROUTINE MergeXLXX
!
!------------------------------------------------------------------
!
   SUBROUTINE DiagDispl(IHessL,JHessL,AHessL,SCRPath, &
                        Grad,GradMult,Displ,LagrDispl)
     TYPE(INT_VECT)             :: IHessL,JHessL
     TYPE(DBL_VECT)             :: AHessL,Vect0,Vect1,Vect2
     REAL(DOUBLE),DIMENSION(:)  :: Grad,GradMult,Displ,LagrDispl
     TYPE(DBL_RNK2)             :: FullMat,InvMat
     INTEGER                    :: I,J,NLagr,NIntC,NDim,NCart
     CHARACTER(LEN=*)           :: SCRPath
     TYPE(BMATR)                :: B 
     !
     NLagr=SIZE(GradMult)
     IF(NLagr/=SIZE(LagrDispl)) &
       CALL Halt('Dimension error #1 in DiagDispl.')
     NIntC=SIZE(Grad)
     IF(NIntC/=SIZE(Displ)) &
       CALL Halt('Dimension error #2 in DiagDispl.')
     NDim=SIZE(IHessL%I)-1
     NCart=NDim-NLagr
     !
     CALL ReadBMATR(B,TRIM(SCRPath)//'B')
     IF(NIntC/=SIZE(B%IB,1)) &
       CALL Halt('Dimension error #3 in DiagDispl.')
     !
     CALL New(Vect0,NCart)
     CALL CALC_BxVect(B,Grad,Vect0%D,Trp_O=.TRUE.)
     !
     CALL Sp1x1ToFull(IHessL%I,JHessL%I,AHessL%D,NDim,NDim,FullMat)
     !
     ! Calc. Abs.Val. Inverse of the Hessian
     !
     CALL New(InvMat,(/NDim,NDim/))
     CALL SetDSYEVWork(NDim**2)
     CALL FunkOnSqMat(NDim,AbsInv,FullMat%D,InvMat%D,Unit_O=6)
     CALL UnSetDSYEVWork()
     !
     CALL New(Vect1,NDim)
     Vect1%D(1:NCart)=-Vect0%D(1:NCart)
     IF(NDim>NCart) Vect1%D(NCart+1:NDim)=-GradMult(1:NLagr)
     !
     CALL New(Vect2,NDim)
       CALL DGEMM_NNc(NDim,NDim,1,One,Zero,InvMat%D,Vect1%D,Vect2%D)
     CALL Delete(Vect1)
     !
     Vect0%D(1:NCart)=Vect2%D(1:NCart)
     CALL CALC_BxVect(B,Displ,Vect0%D)
     IF(NDim>NCart) LagrDispl(1:NLagr)=Vect2%D(NCart+1:NDim)
     CALL Delete(Vect2)
     !
     CALL Delete(Vect0)
     CALL Delete(InvMat)
     CALL Delete(FullMat)
     CALL Delete(B)
   END SUBROUTINE DiagDispl
!
!------------------------------------------------------------------
!
   SUBROUTINE BMatrConstr(IBc,JBc,ABc,XYZ,IntCs,GCoordCtrl, &
                          GTrfCtrl,SCRPath,NLagr,NCart)
     TYPE(INT_VECT)        :: IBc,JBc,JBc2
     TYPE(DBL_VECT)        :: ABc,ABc2
     INTEGER               :: I,J,L,LL,JJ,NLagr,NCart,NIntC,NConstr,NZ
     CHARACTER(LEN=*)      :: SCRPath
     TYPE(INTC)            :: IntCs
     TYPE(BMATR)           :: B     
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(TrfCtrl)         :: GTrfCtrl
     TYPE(CoordCtrl)       :: GCoordCtrl
     !
     NIntC=SIZE(IntCs%Def)
     !!! calling the whole B matrix is not necessary, 
     !!! effort can be saved by calling it only for constraints !!!
     CALL BMatrix(XYZ,NIntC,IntCs,B,GCoordCtrl%LinCrit, &
                  GTrfCtrl%DoClssTrf,GTrfCtrl%ThreeAt)
     !
     CALL New(IBc,NLagr+1)
     CALL New(JBc2,NLagr*12)
     CALL New(ABc2,NLagr*12)
     IBc%I(1)=1
     NZ=0
     NConstr=0 
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NLagr)CALL Halt('Dimension error #2 in BMatrConstr')
         DO J=1,4
           LL=IntCs%Atoms(I,J)
           IF(LL==0) EXIT
           JJ=(J-1)*3
           LL=(LL-1)*3
           DO L=1,3
             NZ=NZ+1
             JBc2%I(NZ)=LL+L 
             ABc2%D(NZ)=B%B(I,JJ+L)
           ENDDO
         ENDDO
         IBc%I(NConstr+1)=NZ+1
       ENDIF
     ENDDO
     !
     CALL New(JBc,NZ)
     CALL New(ABc,NZ)
     JBc%I(1:NZ)=JBc2%I(1:NZ)
     ABc%D(1:NZ)=ABc2%D(1:NZ)
     CALL Delete(JBc2)
     CALL Delete(ABc2)
     !
     CALL Delete(B) 
   END SUBROUTINE BMatrConstr
!
!------------------------------------------------------------------
!
   SUBROUTINE DistMatr(ITop,JTop,ATop,IntCs,NatmsLoc,NStre)
     TYPE(INT_VECT)       :: ITop,JTop,NAux
     TYPE(DBL_VECT)       :: ATop
     TYPE(INTC)           :: IntCs
     INTEGER              :: NatmsLoc,I,J,NBond,NStre,I1,I2,NIntC
     INTEGER              :: I1off,I2Off
     !
     NIntC=SIZE(IntCs%Def)
     CALL New(ITop,NatmsLoc+1)
     CALL New(JTop,2*NStre)
     CALL New(ATop,2*NStre)
     CALL New(NAux,NatmsLoc)
     NBond=0
     NAux%I=0
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)=='STRE') THEN
         I1=IntCs%Atoms(I,1)
         I2=IntCs%Atoms(I,2)
         NAux%I(I1)=NAux%I(I1)+1
         NAux%I(I2)=NAux%I(I2)+1
       ENDIF 
     ENDDO
     !
     ITop%I=0
     ITop%I(1)=1
     DO I=1,NatmsLoc
       ITop%I(I+1)=ITop%I(I)+NAux%I(I)
     ENDDO
     !
     NAux%I=0
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)=='STRE') THEN
         I1=IntCs%Atoms(I,1)
         I2=IntCs%Atoms(I,2)
         NAux%I(I1)=NAux%I(I1)+1
         NAux%I(I2)=NAux%I(I2)+1
         I1off=ITop%I(I1)+NAux%I(I1)-1
         I2off=ITop%I(I2)+NAux%I(I2)-1
         JTop%I(I1off)=I2
         ATop%D(I1off)=IntCs%Value(I)
         JTop%I(I2off)=I1
         ATop%D(I2off)=IntCs%Value(I)
       ENDIF 
     ENDDO
     !
     CALL Delete(NAux)
   END SUBROUTINE DistMatr
!
!------------------------------------------------------------------
!
   FUNCTION HasHBond(NJJ1,NJJ2) 
     LOGICAL                :: HasHBond
     INTEGER                :: NJJ1,NJJ2
     !
     HasHBond=.FALSE.
     IF((NJJ1==1.AND.HasLigand(NJJ2)).OR.(NJJ2==1.AND.HasLigand(NJJ1))) THEN
       HasHBond=.TRUE.
     ENDIF
   END FUNCTION HasHBond
!
!------------------------------------------------------------------
!
   SUBROUTINE LinearHBond(NatmsLoc,AtNum,NBond, &
                          HBondMark,BondIJ,BondLength)
     ! This subroutine selects out those H-bonds, which are most close
     ! to linearity from the set of all local H-bonds.
     INTEGER,DIMENSION(:)   :: AtNum    
     INTEGER                :: NatmsLoc,NBond,I,J,J1,J2,NI1,NI2,MaxHBond
     INTEGER                :: NBondNew,I1,I2,IB1,IB2,ISelect1,ISelect2
     TYPE(INT_RNK2)         :: BondIJ,BondIJNew,HBondList
     TYPE(INT_VECT)         :: CountHBonds,HBondMark,HBondMarkNew,Sign
     TYPE(DBL_VECT)         :: BondLength,BondLengthNew
     REAL(DOUBLE)           :: DistMin,D1,D2,D12
     !
     MaxHBond=10
     CALL New(CountHBonds,NatmsLoc)
     CALL New(HBondList,(/NatmsLoc,MaxHBond/))
     CALL New(Sign,NatmsLoc)
     CALL New(BondIJNew,(/2,NBond/))
     CALL New(HBondMarkNew,NBond)
     CALL New(BondLengthNew,NBond)
     !
     ! Fill Non-H-bonds into new array
     !
     NBondNew=0
     DO I=1,NBond
       IF(HBondMark%I(I)==0) THEN
         NBondNew=NBondNew+1
         BondIJNew%I(1:2,NBondNew)=BondIJ%I(1:2,I)
         HBondMarkNew%I(NBondNew)=HBondMark%I(I)
         BondLengthNew%D(NBondNew)=BondLength%D(I)
       ENDIF
     ENDDO
     !
     CountHBonds%I=0
     HBondList%I=0
     Sign%I=0
     DO I=1,NBond
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       NI1=AtNum(I1)
       NI2=AtNum(I2)
       IF(HBondMark%I(I)==1) THEN
         IF(NI1==1) THEN
           CountHBonds%I(I1)=CountHBonds%I(I1)+1
           IF(CountHBonds%I(I1)>MaxHBond) &
             CALL Halt('Dimension error in LinearHBond')
           HBondList%I(I1,CountHBonds%I(I1))=I
         ENDIF
         IF(NI2==1) THEN
           CountHBonds%I(I2)=CountHBonds%I(I2)+1
           IF(CountHBonds%I(I2)>MaxHBond) &
             CALL Halt('Dimension error in LinearHBond')
           HBondList%I(I2,CountHBonds%I(I2))=I
         ENDIF
       ENDIF
     ENDDO
     !
     DO I=1,NatmsLoc
       IF(CountHBonds%I(I)/=0.AND.Sign%I(I)==0) THEN
         ISelect1=HBondList%I(I,1)
         ISelect2=0
         DistMin=BondLength%D(Iselect1)
         DO J1=1,CountHBonds%I(I)
           IB1=HBondList%I(I,J1)
           D1=BondLength%D(IB1)
           DO J2=J1+1,CountHBonds%I(I)
             IB2=HBondList%I(I,J2)
             D2=BondLength%D(IB2)
             D12=D1+D2
             IF(DistMin>D12) THEN
               DistMin=D12
               ISelect1=IB1
               ISelect2=IB2
             ENDIF
           ENDDO
         ENDDO
         NBondNew=NBondNew+1
         BondIJNew%I(1:2,NBondNew)=BondIJ%I(1:2,ISelect1)
         HBondMarkNew%I(NBondNew)=HBondMark%I(ISelect1)
         BondLengthNew%D(NBondNew)=BondLength%D(ISelect1)
         IF(ISelect2/=0) THEN
           NBondNew=NBondNew+1
           BondIJNew%I(1:2,NBondNew)=BondIJ%I(1:2,ISelect2)
           HBondMarkNew%I(NBondNew)=HBondMark%I(ISelect2)
           BondLengthNew%D(NBondNew)=BondLength%D(ISelect2)
         ENDIF
         Sign%I(I)=1
       ENDIF  
     ENDDO
     !
     CALL Delete(BondIJ)
     CALL Delete(HBondMark)
     CALL Delete(BondLength)
     NBond=NBondNew
     CALL New(BondIJ,(/2,NBond/))
     CALL New(HBondMark,NBond)
     CALL New(BondLength,NBond)
     DO I=1,NBond
       BondIJ%I(1:2,I)=BondIJNew%I(1:2,I)
       HBondMark%I(I)=HBondMarkNew%I(I)
       BondLength%D(I)=BondLengthNew%D(I)
     ENDDO
     !
     CALL Delete(Sign)
     CALL Delete(HBondList)
     CALL Delete(CountHBonds)
     CALL Delete(BondIJNew)
     CALL Delete(HBondMarkNew)
     CALL Delete(BondLengthNew)
   END SUBROUTINE LinearHBond
!
!------------------------------------------------------------------
!
   SUBROUTINE TranslToAt1(VectCart,ThreeAt,Vect_O)
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: Vect_O
     REAL(DOUBLE),DIMENSION(:)          :: VectCart
     REAL(DOUBLE)                       :: V1,V2,V3
     INTEGER,DIMENSION(3)               :: ThreeAt
     INTEGER                            :: I,J,NCart,NatmsLoc,Istart
     !
     NCart=SIZE(VectCart)
     NatmsLoc=NCart/3
     IF(NCart/=3*NatmsLoc) CALL Halt('Dimension error in TranslToAt1')
     !
     IStart=3*(ThreeAt(1)-1)
     IF(PRESENT(Vect_O)) THEN
       V1=Vect_O(1)
       V2=Vect_O(2)
       V3=Vect_O(3)
     ELSE
       V1=VectCart(IStart+1)
       V2=VectCart(IStart+2)
       V3=VectCart(IStart+3)
     ENDIF
     !
     DO I=1,NatmsLoc
       J=3*(I-1)  
       VectCart(J+1)=VectCart(J+1)-V1
       VectCart(J+2)=VectCart(J+2)-V2
       VectCart(J+3)=VectCart(J+3)-V3
     ENDDO
   END SUBROUTINE TranslToAt1
!
!------------------------------------------------------------------
!
   SUBROUTINE RotToAt(VectCart,Rot,Rev_O)
     REAL(DOUBLE),DIMENSION(:)  :: VectCart
     REAL(DOUBLE),DIMENSION(3,3):: Rot     
     INTEGER                    :: I,J,NCart,NatmsLoc
     TYPE(DBL_VECT)             :: Vect,Vect2
     LOGICAL,OPTIONAL           :: Rev_O
     LOGICAL                    :: Reverse
     !
     NCart=SIZE(VectCart)
     NatmsLoc=NCart/3
     CALL New(Vect,3)
     CALL New(Vect2,3)
     Reverse=.FALSE.
     IF(PRESENT(Rev_O)) THEN
       Reverse=Rev_O
     ENDIF
     !
     DO I=1,NatmsLoc
       J=3*(I-1)
       Vect%D(1:3)=VectCart(J+1:J+3)
       IF(Reverse) THEN
         CALL DGEMM_TNc(3,3,1,One,Zero,Rot,Vect%D,Vect2%D)
       ELSE
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot,Vect%D,Vect2%D)
       ENDIF
       VectCart(J+1:J+3)=Vect2%D
     ENDDO 
     CALL Delete(Vect2)
     CALL Delete(Vect)
   END SUBROUTINE RotToAt
!
!------------------------------------------------------------------
!
   SUBROUTINE ZeroTrRots(VectCart,ThreeAt)
     REAL(DOUBLE),DIMENSION(:) :: VectCart
     INTEGER,DIMENSION(3)      :: ThreeAt
     INTEGER                   :: I,J
     !
     I=3*(ThreeAt(1)-1)
     VectCart(I+1:I+3)=Zero
     I=3*(ThreeAt(2)-1)
     VectCart(I+2:I+3)=Zero
     I=ThreeAt(3)
     VectCart(3*I)=Zero
   END SUBROUTINE ZeroTrRots
!
!------------------------------------------------------------------
!
   SUBROUTINE ReSetConstr(IntCs,XYZ)
     TYPE(INTC)                  ::  IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER                     :: I,J,NIntC
     !
     NIntC=SIZE(IntCs%Def)
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         J=IntCs%Atoms(I,1)
         IF(IntCs%Def(I)(1:5)=='CARTX') IntCs%ConstrValue(I)=XYZ(1,J)
         IF(IntCs%Def(I)(1:5)=='CARTY') IntCs%ConstrValue(I)=XYZ(2,J)
         IF(IntCs%Def(I)(1:5)=='CARTZ') IntCs%ConstrValue(I)=XYZ(3,J)
       ENDIF 
     ENDDO
   END SUBROUTINE ReSetConstr
!
!------------------------------------------------------------------
!
   SUBROUTINE ThreeConstr(IntCs,Top12,NatmsLoc,At1,At2,At3)
     TYPE(INTC)     :: IntCs
     INTEGER        :: At1,At2,At3,I,J,NIntC,II,JMax,NatmsLoc
     TYPE(INT_RNK2) :: CountC,Top12
     TYPE(INT_VECT) :: AtVect,Count2
     !
     NIntC=SIZE(IntCs%Def)
     CALL New(CountC,(/NatmsLoc,3/))
     CALL New(Count2,NatmsLoc)
     CALL New(AtVect,3)
     CountC%I=0
     Count2%I=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         J=IntCs%Atoms(I,1)
         IF(IntCs%Def(I)(1:5)=='CARTX') THEN
           CountC%I(J,1)=1
         ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
           CountC%I(J,2)=1
         ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
           CountC%I(J,3)=1
         ENDIF
       ENDIF 
     ENDDO
     !
     DO I=1,NatmsLoc
       Count2%I(I)=CountC%I(I,1)+CountC%I(I,2)+CountC%I(I,3)
     ENDDO
     !
     At1=0
     At2=0
     At3=0
     AtVect%I=0
     DO II=1,3
       JMax=MAXVAL(Count2%I)
       IF(JMax>0) THEN
         JMax=0
         DO I=1,NatmsLoc
           J=Top12%I(I,1)*Count2%I(I)
           IF(J>JMax) THEN
             AtVect%I(II)=I
             JMax=J
           ENDIF
         ENDDO
         Count2%I(AtVect%I(II))=0
       ENDIF
     ENDDO
     At1=AtVect%I(1)
     At2=AtVect%I(2)
     At3=AtVect%I(3)
     !
     ! fix only the origin
     !
     At2=0
     At3=0
     !
     CALL Delete(AtVect)
     CALL Delete(Count2)
     CALL Delete(CountC)
   END SUBROUTINE ThreeConstr
!
!------------------------------------------------------------------
!
   FUNCTION HasMetLig(JJ1,JJ2,NJJ1,NJJ2) 
     LOGICAL                :: HasMetLig
     INTEGER                :: JJ1,JJ2,NJJ1,NJJ2
     !
     HasMetLig=.FALSE.
     IF((HasMetal(NJJ1).AND.HasLigand(NJJ2)).OR. &
        (HasMetal(NJJ2).AND.HasLigand(NJJ1))) THEN
       HasMetLig=.TRUE.
     ENDIF
   END FUNCTION HasMetLig
!
!------------------------------------------------------------------
!
   FUNCTION HasLigand(ICharge) 
     LOGICAL :: HasLigand
     INTEGER :: ICharge
     !
     HasLigand=.FALSE.
     IF(ANY(HBondList(:)==ICharge)) HasLigand=.TRUE.
   END FUNCTION HasLigand
!
!------------------------------------------------------------------
!
   FUNCTION HasMetal(ICharge)
     LOGICAL    :: HasMetal
     INTEGER    :: ICharge,I,J
     !
     HasMetal=.FALSE.
     IF(( 3<=ICharge.AND.ICharge<= 5).OR. &
        (11<=ICharge.AND.ICharge<=14).OR. &
        (19<=ICharge.AND.ICharge<=33).OR. &
        (37<=ICharge.AND.ICharge<=52).OR. &
        (55<=ICharge.AND.ICharge<=84).OR. &
        (87<=ICharge.AND.ICharge<=113)) THEN
       HasMetal=.TRUE.
     ENDIF
   END FUNCTION HasMetal
!
!------------------------------------------------------------------
!
   SUBROUTINE OutPSelect(I1,Top12,XYZ,II2,II3,II4,AngleSum)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INT_RNK2)              :: Top12
     INTEGER                     :: I1,I2,I3,I4,NDim
     INTEGER                     :: II2,II3,II4
     INTEGER                     :: J2,J3,J4
     REAL(DOUBLE)                :: AngleSum,Sum,Conv,Planar,ASum,TwoPi
     REAL(DOUBLE)                :: A1,A2,A3
     !
     AngleSum=-One
     Planar=1.D99
     TwoPi=Two*Pi
     Conv=180.D0/PI
     NDim=Top12%I(I1,1)
     II2=Top12%I(I1,2)
     II3=Top12%I(I1,3)
     II4=Top12%I(I1,4)
     DO J2=1,NDim
       I2=Top12%I(I1,J2+1)
       DO J3=J2+1,NDim
         I3=Top12%I(I1,J3+1)
         CALL BENDValue(XYZ(1:3,I2),XYZ(1:3,I1),XYZ(1:3,I3),A1)
         DO J4=J3+1,NDim
           I4=Top12%I(I1,J4+1)
           CALL BENDValue(XYZ(1:3,I2),XYZ(1:3,I1),XYZ(1:3,I4),A2)
           CALL BENDValue(XYZ(1:3,I3),XYZ(1:3,I1),XYZ(1:3,I4),A3)
           Sum=A1+A2+A3
           ASum=ABS(TwoPI-Sum)
           IF(Planar>ASum) THEN
             Planar=ASum
             AngleSum=Sum
             II2=I2
             II3=I3
             II4=I4
           ENDIF
         ENDDO
       ENDDO
     ENDDO
   END SUBROUTINE OutPSelect
!
!------------------------------------------------------------------
!
   END MODULE InCoords

