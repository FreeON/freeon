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
MODULE ONX2List
!H=================================================================================
!H MODULE ONX2List
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB
!H
!H  PRIVATE:
!H  o SUB
!H
!H  OPTIONS:
!H  DEBUGING: Use -DONX2_DBUG to print some stuff.
!H  INFO    : Use -DONX2_INFO to print some stuff.
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
!#define ONX2_DBUG
#ifdef ONX2_DBUG
#define ONX2_INFO
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONX2DataType
  USE InvExp
  USE ONXParameters
  !
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  !PRIVATE
  !
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: MakeList
  PUBLIC  :: MakeGList
  PUBLIC  :: MakeHList
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: GetAtomPair_
  PRIVATE :: GetAtomPairG_
  !
CONTAINS
  !
  !

#ifdef ONX2_PARALLEL
  SUBROUTINE MakeList(List,GMc,BSc,GMp,BSp,CS_OUT,Ptr1,Nbr1,Ptr2,Nbr2)
#else
  SUBROUTINE MakeList(List,GMc,BSc,GMp,BSp,CS_OUT)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE MakeList(List,GMc,BSc,GMp,BSp,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CList), DIMENSION(:), POINTER :: List
    TYPE(CRDS)                         :: GMc,GMp
    TYPE(BSET)                         :: BSc,BSp
    TYPE(CellSet)                      :: CS_OUT
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT)              :: Ptr1
    TYPE(INT_VECT)              :: Ptr2
    INTEGER       , INTENT(IN ) :: Nbr1
    INTEGER       , INTENT(OUT) :: Nbr2
#endif
    !-------------------------------------------------------------------
    TYPE(ANode), POINTER        :: NodeA
    TYPE(AtomInfo)              :: ACAtmInfo
    INTEGER                     :: AtA,AtC,KA,KC,CFA,CFC,iCell,CFAC
    INTEGER                     :: NFPair,I,IntType,LocNInt,NBFA,NBFC
    INTEGER                     :: iErr
    REAL(DOUBLE)                :: Dum,AC2,NInts,ThresholdDistance
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT)              :: Tmp2
    INTEGER                     :: iC,Idx2
#endif
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(  CS_OUT%NCells*MaxShelPerAtmBlk**2) :: RInt
    INTEGER      , DIMENSION(3*CS_OUT%NCells*MaxShelPerAtmBlk**2) :: Indx
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(MaxFuncPerAtmBlk**4) :: C
    TYPE(AtomPr) , DIMENSION(:), ALLOCATABLE :: ACAtmPair
    !-------------------------------------------------------------------
    REAL(DOUBLE) , EXTERNAL :: DGetAbsMax
    !-------------------------------------------------------------------
    !
    integer :: oa,lda,ob,ldb,oc,ldc,od,ldd
    integer :: isize,idum
    real(double) :: ChkSum
    !
    ChkSum=0.0d0
    !
    ! Allocate arrays.
    ALLOCATE(ACAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In MakeList: Allocation problem.')
    !
    !Simple check Simple check Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,max(BSc%BfKnd%I(GMc%AtTyp%I(i)),BSp%BfKnd%I(GMp%AtTyp%I(i))))
       IF((max(BSc%BfKnd%I(GMc%AtTyp%I(i)),BSp%BfKnd%I(GMp%AtTyp%I(i))))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(max(BSc%BfKnd%I(GMc%AtTyp%I(i)),BSp%BfKnd%I(GMp%AtTyp%I(i))))**4
          write(*,*) 'MaxShelPerAtmBlk',MaxShelPerAtmBlk
          write(*,*) 'SIZE(ACAtmPair)=',MaxShelPerAtmBlk**2*CS_OUT%NCells
          write(*,*) 'SIZE(C)',MaxFuncPerAtmBlk**4
          STOP 'In MakeList: Incrase the size of C'
       ENDIF
    enddo
    !Simple check Simple check Simple check Simple check Simple check Simple check
    !
    idum=0
    NULLIFY(NodeA)
    NInts=0.0d0
    !
#ifdef ONX2_PARALLEL
    CALL New(Tmp2,NAtoms)
    CALL INT_VECT_EQ_INT_SCLR(NAtoms,Tmp2%I(1),0)
#endif
    !
    ThresholdDistance=-10.d0*DLOG(Thresholds%Dist)
    !
#ifdef ONX2_PARALLEL
    DO iC = 1,Nbr1
       AtC = Ptr1%I(iC)
#else
    DO AtC=1,NAtoms ! Run over AtC
#endif
       !
       KC=GMp%AtTyp%I(AtC)
       NBFC=BSp%BfKnd%I(KC)
       ACAtmInfo%Atm2X=GMp%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GMp%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GMp%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       DO AtA=1,NAtoms ! Run over AtA
          !
          KA=GMc%AtTyp%I(AtA)
          NBFA=BSc%BfKnd%I(KA)
          ACAtmInfo%Atm1X=GMc%Carts%D(1,AtA)
          ACAtmInfo%Atm1Y=GMc%Carts%D(2,AtA)
          ACAtmInfo%Atm1Z=GMc%Carts%D(3,AtA)
          ACAtmInfo%K1=KA
          !
          ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
          ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
          ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
          !
          ! Initialize some cell variables.
          NFPair=0
          !
          DO iCell=1,CS_OUT%NCells ! Run over R
             !
             ! Check the interatomic distance between boxes.
             AC2 =  (ACAtmInfo%Atm12X-CS_OUT%CellCarts%D(1,iCell))**2+ &
                    (ACAtmInfo%Atm12Y-CS_OUT%CellCarts%D(2,iCell))**2+ &
                    (ACAtmInfo%Atm12Z-CS_OUT%CellCarts%D(3,iCell))**2
             !
             ! Set the range for range of exchange.
             DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
             !
             ! Cycle if needed.
#ifdef GTRESH
             IF(AC2.GT.ThresholdDistance) CYCLE
#endif
             !
             ! Get the atom pair.
             CALL GetAtomPair_(ACAtmInfo,ACAtmPair,BSc,BSp,CS_OUT%CellCarts%D(1,iCell))
             !
             CFAC=0
             DO CFA=1,BSc%NCFnc%I(KA) ! Run over blkfunc on A
                DO CFC=1,BSp%NCFnc%I(KC) ! Run over blkfunc on C
                   !
                   CFAC=CFAC+1
                   idum=idum+ACAtmPair(CFAC)%SP%L
                   !
                   ! Compute integral type.
                   IntType=ACAtmPair(CFAC)%SP%IntType
!write(*,*) 'IntType',IntType
                   !
                   LocNInt=(BSc%LStop%I(CFA,KA)-BSc%LStrt%I(CFA,KA)+1)**2* &
                           (BSp%LStop%I(CFC,KC)-BSp%LStrt%I(CFC,KC)+1)**2
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(LocNInt,C(1),0.0D0)
                   !
                   ! The integral interface.
!VW                   INCLUDE 'ERIListInterface.Inc'
                   INCLUDE 'ERIListInterfaceB.Inc'
                   !
                   Dum=DSQRT(DGetAbsMax(LocNInt,C(1)))
                   ChkSum=ChkSum+SUM(ABS(C(1:LocNInt)))
                   !
#ifdef GONX2_DBUG
                   WRITE(*,'(A,E22.15,6(A,I6))') 'Dum',Dum, &
                          ' LocNInt',LocNInt,' IntType',IntType,' CFC',CFC,' CFA',CFA,' AtA',AtA,' AtC',AtC
#endif
                   !
                   ! Keep the cell if needed.
#ifdef GTRESH
                IF(Dum.GT.Thresholds%TwoE) THEN
#endif
                   NFPair=NFPair+1
                   RInt(NFPair)=Dum
                   Indx(3*(NFPair-1)+1)=CFA
                   Indx(3*(NFPair-1)+2)=CFC
                   Indx(3*(NFPair-1)+3)=iCell
                   !
#ifdef GTRESH
                ELSE
                   !write(*,*) 'Int smaller that treshold',Dum,Thresholds%TwoE
                ENDIF
#endif
                   !
                   NInts=NInts+DBLE(LocNInt)
                ENDDO ! End over blkfunc on C
             ENDDO ! End over blkfunc on A
             !
          ENDDO ! End R
          !
          ! Check if no cell.
          IF(NFPair.EQ.0) CYCLE
          !
          ! Allocate a new node.
          ALLOCATE(NodeA)
          NULLIFY(NodeA%AtmNext)
          NodeA%Atom=AtA
          NodeA%NFPair=NFPair
          !
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
#ifdef ONX2_PARALLEL
          Tmp2%I(AtA)=1
#endif
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          !
          ! Order the Cell list.
          CALL QSDis(RInt(1),Indx(1),NFPair,-2)
          !
          ! Allocate the Cell array.
          ALLOCATE(NodeA%Indx(3,NFPair),NodeA%RInt(NFPair),STAT=iErr)
          IF(iErr.NE.0) CALL Halt('In MakeList: Allocation problem.')
          !
          ! Copy the arrays.
          DO I=1,NFPair
             NodeA%RInt(I)=RInt(I)
             NodeA%Indx(1,I)=Indx(3*(I-1)+1)
             NodeA%Indx(2,I)=Indx(3*(I-1)+2)
             NodeA%Indx(3,I)=Indx(3*(I-1)+3)
          ENDDO
          !
          ! Insert the new node (at the right place).
          CALL InsertNode(List(AtC)%GoList,NodeA)
          !
       ENDDO ! End AtA
       !
    ENDDO ! End AtC
    !
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
#ifdef ONX2_PARALLEL
    Nbr2=SUM(Tmp2%I)
    CALL New(Ptr2,Nbr2)
    CALL INT_VECT_EQ_INT_SCLR(Nbr2,Ptr2%I(1),0)
    Idx2=1
    DO AtA=1,NAtoms
       IF(Tmp2%I(AtA).EQ.1) THEN
          Ptr2%I(Idx2)=AtA
          Idx2=Idx2+1
       ENDIF
    ENDDO
    CALL Delete(Tmp2)
#endif
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    !write(*,*) 'Number of primitive pairs ',idum
    !write(*,*) 'ChkSum',ChkSum
    !
    ! DeAllocate arrays.
    DEALLOCATE(ACAtmPair,STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In MakeGList: DeAllocation problem.')
    !
#ifdef ONX2_INFO
#ifdef ONX2_PARALLEL
  IF(MyID.EQ.ROOT) THEN
#endif
    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) 'MakeList Statistic.'
    WRITE(*,'(A,F22.1)') ' Nbr ERI=',NInts
    WRITE(*,*) '-------------------------------------'
#ifdef ONX2_PARALLEL
  ENDIF
#endif
#endif
    !
  END SUBROUTINE MakeList
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE MakeGList(List,GMc,BSc,CS_OUT,Ptr1,Nbr1,Ptr2,Nbr2)
#else
  SUBROUTINE MakeGList(List,GMc,BSc,CS_OUT)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE MakeGList(List,GMc,BSc,CS_OUT)
!H  Does the thresholding based on the ERIs.
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CList) , DIMENSION(:), POINTER :: List
    TYPE(CRDS)                          :: GMc
    TYPE(BSET)                          :: BSc
    TYPE(CellSet)                       :: CS_OUT
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT)              :: Ptr1
    TYPE(INT_VECT)              :: Ptr2
    INTEGER       , INTENT(IN ) :: Nbr1
    INTEGER       , INTENT(OUT) :: Nbr2
#endif
    !-------------------------------------------------------------------
    TYPE(ANode) , POINTER       :: NodeA
    TYPE(AtomInfo)              :: ACAtmInfo
    INTEGER                     :: AtA,AtC,KA,KC,CFA,CFC,iCell,CFAC
    INTEGER                     :: NFPair,I,IntType,LocNInt,NBFA,NBFC
    INTEGER                     :: iErr
    REAL(DOUBLE)                :: Dum,AC2,NInts,ThresholdDistance
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT)              :: Tmp2
    INTEGER                     :: iC,Idx2
#endif
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(  CS_OUT%NCells*MaxShelPerAtmBlk**2) :: RInt
    INTEGER      , DIMENSION(3*CS_OUT%NCells*MaxShelPerAtmBlk**2) :: Indx
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(12*MaxFuncPerAtmBlk**4) :: C
    TYPE(AtomPr) , DIMENSION(:), ALLOCATABLE :: ACAtmPair
    !-------------------------------------------------------------------
    REAL(DOUBLE) , EXTERNAL :: DGetAbsMax
    !-------------------------------------------------------------------
    !
    integer :: oa,lda,ob,ldb,oc,ldc,od,ldd
    integer :: isize,idum
    real(double) :: ChkSum
    !
    ChkSum=0.0d0
    !
    ! Allocate arrays.
    ALLOCATE(ACAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In MakeGList: Allocation problem.')
    !
    !Simple check Simple check Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BSc%BfKnd%I(GMc%AtTyp%I(i)))
       IF((BSc%BfKnd%I(GMc%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(BSc%BfKnd%I(GMc%AtTyp%I(i)))**4
          write(*,*) 'MaxShelPerAtmBlk',MaxShelPerAtmBlk
          write(*,*) 'SIZE(ACAtmPair)=',MaxShelPerAtmBlk**2*CS_OUT%NCells
          write(*,*) 'SIZE(C)',12*MaxFuncPerAtmBlk**4
          STOP 'In MakeGList: Incrase the size of C'
       ENDIF
    enddo
    !Simple check Simple check Simple check Simple check Simple check Simple check
    !
    !
    idum=0
    NULLIFY(NodeA)
    NInts=0.0D0
    !
#ifdef ONX2_PARALLEL
    CALL New(Tmp2,NAtoms)
    CALL INT_VECT_EQ_INT_SCLR(NAtoms,Tmp2%I(1),0)
#endif
    !
    !
    ThresholdDistance=-10.d0*DLOG(Thresholds%Dist)
    !
#ifdef ONX2_PARALLEL
    DO iC = 1,Nbr1
       AtC = Ptr1%I(iC)
#else
    DO AtC=1,NAtoms ! Run over AtC
#endif
       !
       KC=GMc%AtTyp%I(AtC)
       NBFC=BSc%BfKnd%I(KC)
       !
       ACAtmInfo%Atm2X=GMc%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GMc%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GMc%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       DO AtA=1,NAtoms ! Run over AtA
          !
          KA=GMc%AtTyp%I(AtA)
          NBFA=BSc%BfKnd%I(KA)
          !
          ACAtmInfo%Atm1X=GMc%Carts%D(1,AtA)
          ACAtmInfo%Atm1Y=GMc%Carts%D(2,AtA)
          ACAtmInfo%Atm1Z=GMc%Carts%D(3,AtA)
          ACAtmInfo%K1=KA
          !
          ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
          ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
          ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
          !
          ! Initialize some cell variables.
          NFPair=0
          !
          DO iCell=1,CS_OUT%NCells ! Run over R
             !
             ! Check the interatomic distance between boxes.
             AC2 =  (ACAtmInfo%Atm12X-CS_OUT%CellCarts%D(1,iCell))**2+ &
                    (ACAtmInfo%Atm12Y-CS_OUT%CellCarts%D(2,iCell))**2+ &
                    (ACAtmInfo%Atm12Z-CS_OUT%CellCarts%D(3,iCell))**2
             !
             ! Set the range for range of exchange.
             DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
             !
             ! Cycle if needed.
#ifdef GTRESH
             IF(AC2.GT.ThresholdDistance) CYCLE
#endif
             !
             ! Get the atom pair.
             CALL GetAtomPair_(ACAtmInfo,ACAtmPair,BSc,BSc,CS_OUT%CellCarts%D(1,iCell))
             !CALL GetAtomPairG_(ACAtmInfo,ACAtmPair,BSc,CS_OUT%CellCarts%D(1,iCell))
             !
             CFAC=0
             DO CFA=1,BSc%NCFnc%I(KA) ! Run over blkfunc on A
                DO CFC=1,BSc%NCFnc%I(KC) ! Run over blkfunc on C
                   !
                   CFAC=CFAC+1
                   idum=idum+ACAtmPair(CFAC)%SP%L
                   !
                   ! Compute integral type.
                   IntType=ACAtmPair(CFAC)%SP%IntType
                   !
                   LocNInt=(BSc%LStop%I(CFA,KA)-BSc%LStrt%I(CFA,KA)+1)**2* &
                           (BSc%LStop%I(CFC,KC)-BSc%LStrt%I(CFC,KC)+1)**2
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(LocNInt,C(1),0.0D0)
                   !
                   ! The integral interface.
!VW                   INCLUDE 'ERIListInterface.Inc'
                   INCLUDE 'ERIListInterfaceB.Inc'
                   !INCLUDE 'DERIListInterface.Inc'
                   !
                   Dum=DSQRT(DGetAbsMax(LocNInt,C(1)))
                   ChkSum=ChkSum+SUM(ABS(C(1:LocNInt)))
                   !
#ifdef GONX2_DBUG
                   WRITE(*,'(A,E22.15,6(A,I6))') 'Dum',Dum, &
                          ' LocNInt',LocNInt,' IntType',IntType,' CFC',CFC,' CFA',CFA,' AtA',AtA,' AtC',AtC
#endif
                   !
                   ! Keep the cell if needed.
#ifdef GTRESH
                IF(Dum.GT.Thresholds%TwoE) THEN
#endif
                   NFPair=NFPair+1
                   RInt(NFPair)=Dum
                   Indx(3*(NFPair-1)+1)=CFA
                   Indx(3*(NFPair-1)+2)=CFC
                   Indx(3*(NFPair-1)+3)=iCell
                   !
#ifdef GTRESH
                ELSE
                   !write(*,*) 'Int smaller that treshold',Dum,Thresholds%TwoE
                ENDIF
#endif
                   !
                   NInts=NInts+DBLE(LocNInt)
                ENDDO ! End over blkfunc on A
             ENDDO ! End over blkfunc on C
             !
          ENDDO ! End R
          !
          ! Check if no function pair.
          IF(NFPair.EQ.0) CYCLE
          !
          ! Allocate a new node.
          ALLOCATE(NodeA)
          NULLIFY(NodeA%AtmNext)
          NodeA%Atom=AtA
          NodeA%NFPair=NFPair
          !
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
#ifdef ONX2_PARALLEL
          Tmp2%I(AtA)=1
#endif
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          !
          ! Order the Cell list.
          CALL QSDis(RInt(1),Indx(1),NFPair,-2)
          !
          ! Allocate the Cell array.
          ALLOCATE(NodeA%Indx(3,NFPair),NodeA%RInt(NFPair),STAT=iErr)
          IF(iErr.NE.0) CALL Halt('In MakeGList: Allocation problem.')
          !
          ! Copy the arrays.
          DO I=1,NFPair
             NodeA%RInt(I)=RInt(I)
             NodeA%Indx(1,I)=Indx(3*(I-1)+1)
             NodeA%Indx(2,I)=Indx(3*(I-1)+2)
             NodeA%Indx(3,I)=Indx(3*(I-1)+3)
          ENDDO
          !
          ! Insert the new node (at the right place).
          CALL InsertNode(List(AtC)%GoList,NodeA)
          !
       ENDDO ! End AtA
       !
    ENDDO ! End AtC
    !
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
#ifdef ONX2_PARALLEL
    Nbr2=SUM(Tmp2%I)
    CALL New(Ptr2,Nbr2)
    CALL INT_VECT_EQ_INT_SCLR(Nbr2,Ptr2%I(1),0)
    Idx2=1
    DO AtA=1,NAtoms
       IF(Tmp2%I(AtA).EQ.1) THEN
          Ptr2%I(Idx2)=AtA
          Idx2=Idx2+1
       ENDIF
    ENDDO
    CALL Delete(Tmp2)
#endif
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    !write(*,*) 'Number of primitive pairs ',idum
    !!write(*,*) 'ChkSum',ChkSum
    !
    ! DeAllocate arrays.
    DEALLOCATE(ACAtmPair,STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In MakeGList: DeAllocation problem.')
    !
#ifdef GONX2_INFO
#ifdef ONX2_PARALLEL
    IF(MyID.EQ.ROOT) THEN
#endif
       WRITE(*,*) '-------------------------------------'
       WRITE(*,*) 'MakeList Statistic.'
       WRITE(*,'(A,F22.1)') ' Nbr ERI=',NInts
       WRITE(*,*) '-------------------------------------'
#ifdef ONX2_PARALLEL
    ENDIF
#endif
#endif
    !
  END SUBROUTINE MakeGList
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE MakeHList(List,GMc,BSc,CS_OUT,Ptr1,Nbr1,Ptr2,Nbr2)
#else
  SUBROUTINE MakeHList(List,GMc,BSc,CS_OUT)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE MakeHList(List,GMc,BSc,CS_OUT)
!H  Does the thresholding based on the ERIs.
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CList) , DIMENSION(:), POINTER :: List
    TYPE(CRDS)                          :: GMc
    TYPE(BSET)                          :: BSc
    TYPE(CellSet)                       :: CS_OUT
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT)              :: Ptr1
    TYPE(INT_VECT)              :: Ptr2
    INTEGER       , INTENT(IN ) :: Nbr1
    INTEGER       , INTENT(OUT) :: Nbr2
#endif
    !-------------------------------------------------------------------
    TYPE(ANode) , POINTER       :: NodeA
    TYPE(AtomInfo)              :: ACAtmInfo
    INTEGER                     :: AtA,AtC,KA,KC,CFA,CFC,iCell,CFAC
    INTEGER                     :: NFPair,I,IntType,LocNInt,NBFA,NBFC
    INTEGER                     :: iErr
    REAL(DOUBLE)                :: Dum,AC2,NInts,ThresholdDistance
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT)              :: Tmp2
    INTEGER                     :: iC,Idx2
#endif
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(  CS_OUT%NCells*MaxShelPerAtmBlk**2) :: RInt
    INTEGER      , DIMENSION(3*CS_OUT%NCells*MaxShelPerAtmBlk**2) :: Indx
    !-------------------------------------------------------------------
    REAL(DOUBLE) , DIMENSION(12*MaxFuncPerAtmBlk**4) :: C
    TYPE(AtomPr) , DIMENSION(:), ALLOCATABLE :: ACAtmPair
    !-------------------------------------------------------------------
    REAL(DOUBLE) , EXTERNAL :: DGetAbsMax
    !-------------------------------------------------------------------
    !
    integer :: oa,lda,ob,ldb,oc,ldc,od,ldd
    integer :: isize,idum
    real(double) :: ChkSum
    !
    ChkSum=0.0d0
    !
    ! Allocate arrays.
    ALLOCATE(ACAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In MakeGList: Allocation problem.')
    !
    !Simple check Simple check Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BSc%BfKnd%I(GMc%AtTyp%I(i)))
       IF((BSc%BfKnd%I(GMc%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(BSc%BfKnd%I(GMc%AtTyp%I(i)))**4
          write(*,*) 'MaxShelPerAtmBlk',MaxShelPerAtmBlk
          write(*,*) 'SIZE(ACAtmPair)=',MaxShelPerAtmBlk**2*CS_OUT%NCells
          write(*,*) 'SIZE(C)',12*MaxFuncPerAtmBlk**4
          STOP 'In MakeGList: Incrase the size of C'
       ENDIF
    enddo
    !Simple check Simple check Simple check Simple check Simple check Simple check
    !
    !
    idum=0
    NULLIFY(NodeA)
    NInts=0.0D0
    !
#ifdef ONX2_PARALLEL
    CALL New(Tmp2,NAtoms)
    CALL INT_VECT_EQ_INT_SCLR(NAtoms,Tmp2%I(1),0)
#endif
    !
    !
    ThresholdDistance=-10.d0*DLOG(Thresholds%Dist)
    !
#ifdef ONX2_PARALLEL
    DO iC = 1,Nbr1
       AtC = Ptr1%I(iC)
#else
    DO AtC=1,NAtoms ! Run over AtC
#endif
       !
       KC=GMc%AtTyp%I(AtC)
       NBFC=BSc%BfKnd%I(KC)
       !
       ACAtmInfo%Atm2X=GMc%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GMc%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GMc%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       DO AtA=1,NAtoms ! Run over AtA
          !
          KA=GMc%AtTyp%I(AtA)
          NBFA=BSc%BfKnd%I(KA)
          !
          ACAtmInfo%Atm1X=GMc%Carts%D(1,AtA)
          ACAtmInfo%Atm1Y=GMc%Carts%D(2,AtA)
          ACAtmInfo%Atm1Z=GMc%Carts%D(3,AtA)
          ACAtmInfo%K1=KA
          !
          ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
          ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
          ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
          !
          ! Initialize some cell variables.
          NFPair=0
          !
          DO iCell=1,CS_OUT%NCells ! Run over R
             !
             ! Check the interatomic distance between boxes.
             AC2 =  (ACAtmInfo%Atm12X-CS_OUT%CellCarts%D(1,iCell))**2+ &
                    (ACAtmInfo%Atm12Y-CS_OUT%CellCarts%D(2,iCell))**2+ &
                    (ACAtmInfo%Atm12Z-CS_OUT%CellCarts%D(3,iCell))**2
             !
             ! Set the range for range of exchange.
             DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
             !
             ! Cycle if needed.
#ifdef GTRESH
             IF(AC2.GT.ThresholdDistance) CYCLE
#endif
             !
             ! Get the atom pair.
             CALL GetAtomPair_(ACAtmInfo,ACAtmPair,BSc,BSc,CS_OUT%CellCarts%D(1,iCell))
             !CALL GetAtomPairH_(ACAtmInfo,ACAtmPair,BSc,CS_OUT%CellCarts%D(1,iCell))
             !
             CFAC=0
             DO CFA=1,BSc%NCFnc%I(KA) ! Run over blkfunc on A
                DO CFC=1,BSc%NCFnc%I(KC) ! Run over blkfunc on C
                   !
                   CFAC=CFAC+1
                   idum=idum+ACAtmPair(CFAC)%SP%L
                   !
                   ! Compute integral type.
                   IntType=ACAtmPair(CFAC)%SP%IntType
                   !
                   LocNInt=(BSc%LStop%I(CFA,KA)-BSc%LStrt%I(CFA,KA)+1)**2* &
                           (BSc%LStop%I(CFC,KC)-BSc%LStrt%I(CFC,KC)+1)**2
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(LocNInt,C(1),0.0D0)
                   !
                   ! The integral interface.
                   INCLUDE 'ERIListInterfaceB.Inc'
                   !INCLUDE 'd2ERIListInterface.Inc'
                   !
                   Dum=DSQRT(DGetAbsMax(LocNInt,C(1)))
                   ChkSum=ChkSum+SUM(ABS(C(1:LocNInt)))
                   !
#ifdef GONX2_DBUG
                   WRITE(*,'(A,E22.15,6(A,I6))') 'Dum',Dum, &
                          ' LocNInt',LocNInt,' IntType',IntType,' CFC',CFC,' CFA',CFA,' AtA',AtA,' AtC',AtC
#endif
                   !
                   ! Keep the cell if needed.
#ifdef GTRESH
                IF(Dum.GT.Thresholds%TwoE) THEN
#endif
                   NFPair=NFPair+1
                   RInt(NFPair)=Dum
                   Indx(3*(NFPair-1)+1)=CFA
                   Indx(3*(NFPair-1)+2)=CFC
                   Indx(3*(NFPair-1)+3)=iCell
                   !
#ifdef GTRESH
                ELSE
                   !write(*,*) 'Int smaller that treshold',Dum,Thresholds%TwoE
                ENDIF
#endif
                   !
                   NInts=NInts+DBLE(LocNInt)
                ENDDO ! End over blkfunc on A
             ENDDO ! End over blkfunc on C
             !
          ENDDO ! End R
          !
          ! Check if no function pair.
          IF(NFPair.EQ.0) CYCLE
          !
          ! Allocate a new node.
          ALLOCATE(NodeA)
          NULLIFY(NodeA%AtmNext)
          NodeA%Atom=AtA
          NodeA%NFPair=NFPair
          !
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
#ifdef ONX2_PARALLEL
          Tmp2%I(AtA)=1
#endif
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
          !
          ! Order the Cell list.
          CALL QSDis(RInt(1),Indx(1),NFPair,-2)
          !
          ! Allocate the Cell array.
          ALLOCATE(NodeA%Indx(3,NFPair),NodeA%RInt(NFPair),STAT=iErr)
          IF(iErr.NE.0) CALL Halt('In MakeGList: Allocation problem.')
          !
          ! Copy the arrays.
          DO I=1,NFPair
             NodeA%RInt(I)=RInt(I)
             NodeA%Indx(1,I)=Indx(3*(I-1)+1)
             NodeA%Indx(2,I)=Indx(3*(I-1)+2)
             NodeA%Indx(3,I)=Indx(3*(I-1)+3)
          ENDDO
          !
          ! Insert the new node (at the right place).
          CALL InsertNode(List(AtC)%GoList,NodeA)
          !
       ENDDO ! End AtA
       !
    ENDDO ! End AtC
    !
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
#ifdef ONX2_PARALLEL
    Nbr2=SUM(Tmp2%I)
    CALL New(Ptr2,Nbr2)
    CALL INT_VECT_EQ_INT_SCLR(Nbr2,Ptr2%I(1),0)
    Idx2=1
    DO AtA=1,NAtoms
       IF(Tmp2%I(AtA).EQ.1) THEN
          Ptr2%I(Idx2)=AtA
          Idx2=Idx2+1
       ENDIF
    ENDDO
    CALL Delete(Tmp2)
#endif
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    ! PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL PARALLEL
    !write(*,*) 'Number of primitive pairs ',idum
    !!write(*,*) 'ChkSum',ChkSum
    !
    ! DeAllocate arrays.
    DEALLOCATE(ACAtmPair,STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In MakeGList: DeAllocation problem.')
    !
#ifdef GONX2_INFO
#ifdef ONX2_PARALLEL
    IF(MyID.EQ.ROOT) THEN
#endif
       WRITE(*,*) '-------------------------------------'
       WRITE(*,*) 'MakeList Statistic.'
       WRITE(*,'(A,F22.1)') ' Nbr ERI=',NInts
       WRITE(*,*) '-------------------------------------'
#ifdef ONX2_PARALLEL
    ENDIF
#endif
#endif
    !
  END SUBROUTINE MakeHList
  !
  !
  SUBROUTINE GetAtomPair_(AtmInfo,AtmPair,BSc,BSp,PBC)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPair_(AtmInfo,AtmPair,BSc,BSp,PBC)
!H
!H---------------------------------------------------------------------------------
    USE Thresholding
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)                           :: AtmInfo
    TYPE(AtomPr)  , DIMENSION(:)             :: AtmPair
    TYPE(BSET)                               :: BSc,BSp
    REAL(DOUBLE)  , DIMENSION(3), INTENT(IN) :: PBC
    !-------------------------------------------------------------------
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ,MAXII
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StartL1,StopL1,StartL2,StopL2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ,Cnt
    LOGICAL      :: Switch
    !-------------------------------------------------------------------
    !
    RX=PBC(1)
    RY=PBC(2)
    RZ=PBC(3)
    !
    ! AtmInfo must be related to the atoms in the working cell ONLY.
    ! Then we add the PBC's to have the right interatomic distance.
    R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
    !
    MAXII=0
    CF12=0
    DO CF1=1,BSc%NCFnc%I(AtmInfo%K1)
       MinL1=BSc%ASymm%I(1,CF1,AtmInfo%K1)
       MaxL1=BSc%ASymm%I(2,CF1,AtmInfo%K1)
       Type1=MaxL1*(MaxL1+1)/2+MinL1+1
       StartL1=BSc%LStrt%I(CF1,AtmInfo%K1)
       StopL1=BSc%LStop%I(CF1,AtmInfo%K1)
       !
       DO CF2=1,BSp%NCFnc%I(AtmInfo%K2)
          CF12=CF12+1
          !
          MinL2=BSp%ASymm%I(1,CF2,AtmInfo%K2)
          MaxL2=BSp%ASymm%I(2,CF2,AtmInfo%K2)
          Type2=MaxL2*(MaxL2+1)/2+MinL2+1
          StartL2=BSp%LStrt%I(CF2,AtmInfo%K2)
          StopL2=BSp%LStop%I(CF2,AtmInfo%K2)
          !
          !>>>>
          Switch=Type1.LT.Type2
          AtmPair(CF12)%SP%Switch=Switch
          IF(Switch) THEN
             AtmPair(CF12)%SP%IntType=Type2*100+Type1
          ELSE
             AtmPair(CF12)%SP%IntType=Type1*100+Type2
          ENDIF
          !oIF(Type1.LT.Type2) THEN
          !o   AtmPair(CF12)%SP%IntType=Type2*100+Type1
          !oELSE
          !o   AtmPair(CF12)%SP%IntType=Type1*100+Type2
          !oENDIF
          !<<<<
          !
          II=0
          !
          ! We assume the primitives are ordered (exponants in decressing order).
          DO I1=BSc%NPFnc%I(CF1,AtmInfo%K1),1,-1
             Z1=BSc%Expnt%D(I1,CF1,AtmInfo%K1)
             JJ=0
             !
             !We assume the primitives are ordered (exponants in decressing order).
             DO I2=BSp%NPFnc%I(CF2,AtmInfo%K2),1,-1
                Z2=BSp%Expnt%D(I2,CF2,AtmInfo%K2)
                Expt=Z1+Z2
                InvExpt=1.0d0/Expt
                XiR12=Z2*Z1*InvExpt*R12
                !
                IF(XiR12<PrimPairDistanceThreshold) THEN
                   JJ=JJ+1
                   IJ=JJ+II
                   Cnt=BSc%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)*BSp%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                   AtmPair(CF12)%SP%Cst(1,IJ)=Expt
                   !
                   ! AtmInfo must be related to the atoms in the working cell ONLY.
                   ! Then we add the PBC's to have the right atomic position.
                   AtmPair(CF12)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                   AtmPair(CF12)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                   AtmPair(CF12)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                   AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt*Cnt
                   !
                   !>>>>
                   IF((Type1.NE.2.AND.Type2.EQ.2).OR.(Type2.NE.2.AND.Type1.EQ.2))THEN
                      AtmPair(CF12)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                      AtmPair(CF12)%SP%Cst(7,IJ)=BIG_DBL
                      AtmPair(CF12)%SP%Cst(8,IJ)=BIG_DBL
                   ELSEIF(Type1.EQ.2.AND.Type2.EQ.2)THEN
                      AtmPair(CF12)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                      AtmPair(CF12)%SP%Cst(7,IJ)=BSc%CCoef%D(StartL1+1,I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                      AtmPair(CF12)%SP%Cst(8,IJ)=BSc%CCoef%D(StartL1  ,I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2+1,I2,CF2,AtmInfo%K2)/Cnt
                   ELSE
                      AtmPair(CF12)%SP%Cst(6,IJ)=BIG_DBL
                      AtmPair(CF12)%SP%Cst(7,IJ)=BIG_DBL
                      AtmPair(CF12)%SP%Cst(8,IJ)=BIG_DBL
                   ENDIF
                   !
                   !oIF(IntType.EQ.0102.OR.IntType.EQ.0201) THEN
                   !o   AtmPair(CF12)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1  ,I1,CF1,AtmInfo%K1)* &
                   !o                              BSp%CCoef%D(StartL2  ,I2,CF2,AtmInfo%K2)/Cnt
                   !oELSEIF(IntType.EQ.0202) THEN
                   !o   AtmPair(CF12)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1  ,I1,CF1,AtmInfo%K1)* &
                   !o                              BSp%CCoef%D(StartL2  ,I2,CF2,AtmInfo%K2)/Cnt
                   !o   AtmPair(CF12)%SP%Cst(7,IJ)=BSc%CCoef%D(StartL1+1,I1,CF1,AtmInfo%K1)* &
                   !o                              BSp%CCoef%D(StartL2  ,I2,CF2,AtmInfo%K2)/Cnt
                   !o   AtmPair(CF12)%SP%Cst(8,IJ)=BSc%CCoef%D(StartL1  ,I1,CF1,AtmInfo%K1)* &
                   !o                              BSp%CCoef%D(StartL2+1,I2,CF2,AtmInfo%K2)/Cnt
                   !oENDIF
                   !<<<<
                ELSE
                   ! We can skipp out the loop because the primitives are ordered.
                   EXIT
                ENDIF
             ENDDO
             ! We can skipp out the loop because we did not get any significant primitives.
             IF(JJ.EQ.0) EXIT
             II=II+JJ
             !
          ENDDO
          !
          MAXII=MAX(MAXII,II)
          AtmPair(CF12)%SP%L=II
          !
          ! We reorder the atomic positions if Type2 > Type1.
          ! Needed for the integral evaluations.
          IF(Type1.GE.Type2) THEN
             AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
             AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
             AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
             !
             ! AtmInfo must be related to the atoms in the working cell ONLY.
             ! Then we add the PBC's to have the right atomic position.
             AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
             AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
             AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ
          ELSE
             !
             ! AtmInfo must be related to the atoms in the working cell ONLY.
             ! Then we add the PBC's to have the right atomic position.
             AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm2X+RX
             AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm2Y+RY
             AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm2Z+RZ
             AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm1X
             AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm1Y
             AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm1Z
          ENDIF
          !
#ifdef ONX2_DBUG
          write(*,*) 'Printing from GetAtomPair'
          write(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
                                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
          !
       ENDDO
    ENDDO
    !
    IF(CF12.GT.SIZE(AtmPair,DIM=1)) STOP 'Increase the size of -AtmPair-'
    IF(MAXII.GT.SIZE(AtmPair(1)%SP%Cst,DIM=2)) STOP 'Increase the size of -AtmPair(1)%SP%Cst-'
    !
  END SUBROUTINE GetAtomPair_
  !
  !
  SUBROUTINE GetAtomPairG_(AtmInfo,AtmPair,BS,PBC)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPairG(AtmInfo,AtmPair,BS,PBC)
!H
!H---------------------------------------------------------------------------------
    USE Thresholding
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)                           :: AtmInfo
    TYPE(AtomPr)  , DIMENSION(:)             :: AtmPair
    TYPE(BSET)                               :: BS
    REAL(DOUBLE)  , DIMENSION(3), INTENT(IN) :: PBC
    !-------------------------------------------------------------------
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StartL1,StopL1,StartL2,StopL2
    INTEGER      :: ISwitch1,ISwitch2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    LOGICAL      :: Switch
    !-------------------------------------------------------------------
    !
    RX=PBC(1)
    RY=PBC(2)
    RZ=PBC(3)
    !
    ! AtmInfo must be related to the atoms in the working cell ONLY.
    ! Then we add the PBC's to have the right interatomic distance.
    R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
    !
    CF12=0
    DO CF1=1,BS%NCFnc%I(AtmInfo%K1)
       MinL1=BS%ASymm%I(1,CF1,AtmInfo%K1)
       MaxL1=BS%ASymm%I(2,CF1,AtmInfo%K1)
       Type1=MaxL1*(MaxL1+1)/2+MinL1+1
       StartL1=BS%LStrt%I(CF1,AtmInfo%K1)
       StopL1=BS%LStop%I(CF1,AtmInfo%K1)

       if(Type1==2) stop 'SP shell not yet supported.'

       DO CF2=1,BS%NCFnc%I(AtmInfo%K2)
          CF12=CF12+1
          !
          MinL2=BS%ASymm%I(1,CF2,AtmInfo%K2)
          MaxL2=BS%ASymm%I(2,CF2,AtmInfo%K2)
          Type2=MaxL2*(MaxL2+1)/2+MinL2+1
          AtmPair(CF12)%SP%IntType=Type1*100+Type2
          StartL2=BS%LStrt%I(CF2,AtmInfo%K2)
          StopL2=BS%LStop%I(CF2,AtmInfo%K2)
          !
          ! .NOT.(Do we need to switch the centers?)
          Switch=Type1.GE.Type2
          !
          IF(Switch) THEN
             ISwitch1=6
             ISwitch2=7
          ELSE
             ISwitch1=7
             ISwitch2=6
          ENDIF
          !
          II=0
          !
          ! We assume the primitives are ordered (exponants in decressing order).
          DO I1=BS%NPFnc%I(CF1,AtmInfo%K1),1,-1
             Z1=BS%Expnt%D(I1,CF1,AtmInfo%K1)
             JJ=0
             !
             !We assume the primitives are ordered (exponants in decressing order).
             DO I2=BS%NPFnc%I(CF2,AtmInfo%K2),1,-1
                Z2=BS%Expnt%D(I2,CF2,AtmInfo%K2)
                Expt=Z1+Z2
                InvExpt=1.0d0/Expt
                XiR12=Z2*Z1*InvExpt*R12
                IF(XiR12<PrimPairDistanceThreshold) THEN
                   JJ=JJ+1
                   IJ=JJ+II
                   AtmPair(CF12)%SP%Cst(1,IJ)=Expt
                   !
                   ! AtmInfo must be related to the atoms in the working cell ONLY.
                   ! Then we add the PBC's to have the right atomic position.
                   AtmPair(CF12)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                   AtmPair(CF12)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                   AtmPair(CF12)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                   !AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)* &
                   AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXPInv(XiR12)*InvExpt* &
                                              BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)* &
                                              BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                   !
                   AtmPair(CF12)%SP%Cst(ISwitch1,IJ)=2.0d0*Z1!->ISwitch1=6,7
                   AtmPair(CF12)%SP%Cst(ISwitch2,IJ)=2.0d0*Z2!->ISwitch2=7,6
                   !
                   ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                   ! Here I need to add the correction factor for SP shell.
                   ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                ELSE
                   ! We can skipp out the loop because the primitives are ordered.
                   EXIT
                ENDIF
             ENDDO
             ! We can skipp out the loop because we did not get any significant primitives.
             IF(JJ.EQ.0) EXIT
             II=II+JJ
             !
          ENDDO
          !
          AtmPair(CF12)%SP%L=II
          !
          ! We reorder the atomic positions if Type2 > Type1.
          ! Needed for the integral evaluations.
          IF(Switch) THEN
             !old IF(Type1.GE.Type2) THEN
             AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
             AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
             AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
             !
             ! AtmInfo must be related to the atoms in the working cell ONLY.
             ! Then we add the PBC's to have the right atomic position.
             AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
             AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
             AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ
          ELSE
             !
             ! AtmInfo must be related to the atoms in the working cell ONLY.
             ! Then we add the PBC's to have the right atomic position.
             AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm2X+RX
             AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm2Y+RY
             AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm2Z+RZ
             AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm1X
             AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm1Y
             AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm1Z
          ENDIF
          !
#ifdef ONX2_DBUG
          write(*,*) 'Printing from GetAtomPair'
          write(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
                                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
          !
       ENDDO
    ENDDO
    !
    IF(CF12.GT.SIZE(AtmPair)) STOP 'Increase the size of -AtmPair-'
    !
  END SUBROUTINE GetAtomPairG_
  !
  !
  SUBROUTINE InsertNode(BegPtr,NewPtr)
!H---------------------------------------------------------------------------------
!H SUBROUTINE InsertNode(BegPtr,NewPtr)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(ANode), POINTER :: BegPtr,NewPtr
    !-------------------------------------------------------------------
    TYPE(ANode), POINTER :: TmpPtr
    !-------------------------------------------------------------------
    !
    NULLIFY(TmpPtr)
    !
    ! Check if Beg is null, then insert node and exit.
    IF(.NOT.ASSOCIATED(BegPtr)) THEN
       BegPtr=>NewPtr
       NULLIFY(NewPtr)
       RETURN
    ENDIF
    !
    ! Check if the first element in the list is smaller
    ! than the ones in the new node.
    IF(BegPtr%RInt(1).LE.NewPtr%RInt(1)) THEN
       NewPtr%AtmNext=>BegPtr
       BegPtr=>NewPtr
       NULLIFY(NewPtr)
       RETURN
    ENDIF
    !
    ! General search.
    TmpPtr=>BegPtr
    DO
       ! We are at the end.
       IF(.NOT.ASSOCIATED(TmpPtr%AtmNext)) THEN
          TmpPtr%AtmNext=>NewPtr
          NULLIFY(NewPtr)
          EXIT
       ENDIF
       !
       IF(TmpPtr%AtmNext%RInt(1).LE.NewPtr%RInt(1)) THEN
          NewPtr%AtmNext=>TmpPtr%AtmNext
          TmpPtr%AtmNext=>NewPtr
          NULLIFY(NewPtr)
          EXIT
       ENDIF
       TmpPtr=>TmpPtr%AtmNext
       !
    ENDDO
    !
  END SUBROUTINE InsertNode
  !
  !
  SUBROUTINE AllocList(List,IdxMin,IdxMax)
!H---------------------------------------------------------------------------------
!H SUBROUTINE AllocList(List,LSize)
!H
!H---------------------------------------------------------------------------------
    !
    TYPE(CList), DIMENSION(:), POINTER    :: List
    INTEGER                   , INTENT(IN) :: IdxMin,IdxMax
    !-------------------------------------------------------------------
    INTEGER                                :: I
    !-------------------------------------------------------------------
    !
    ALLOCATE(List(IdxMin:IdxMax))
    !
    DO I=IdxMin,IdxMax
       NULLIFY(List(I)%GoList)
    ENDDO
    !
  END SUBROUTINE AllocList
  !
  !
  SUBROUTINE DeAllocList(List)
!H---------------------------------------------------------------------------------
!H SUBROUTINE DeAllocList(List)
!H
!H---------------------------------------------------------------------------------
    !
    TYPE(CList), DIMENSION(:), POINTER :: List
    !-------------------------------------------------------------------
    TYPE(ANode)              , POINTER :: ListA,ListATmp
    INTEGER                             :: I,IdxMin,IdxMax
    !-------------------------------------------------------------------
    !
    NULLIFY(ListA,ListATmp)
    !
    IdxMin=LBOUND(List,DIM=1)
    IdxMax=UBOUND(List,DIM=1)
    !
    DO I=IdxMin,IdxMax
       !
       ListA=>List(I)%GoList
       !
       DO
          IF(.NOT.ASSOCIATED(ListA)) EXIT
          !
#ifdef POINTERS_IN_DERIVED_TYPES
          IF(.NOT.ASSOCIATED(ListA%Indx).OR. &
                 .NOT.ASSOCIATED(ListA%RInt)) STOP 'Array not allocate for the list.'
          !
          IF(ASSOCIATED(ListA%Indx)) DEALLOCATE(ListA%Indx)
          IF(ASSOCIATED(ListA%RInt)) DEALLOCATE(ListA%RInt)
#else
          IF(.NOT.ALLOCATED(ListA%Indx).OR. &
                 .NOT.ALLOCATED(ListA%RInt)) STOP 'Array not allocate for the list.'
          !
          IF(ALLOCATED(ListA%Indx)) DEALLOCATE(ListA%Indx)
          IF(ALLOCATED(ListA%RInt)) DEALLOCATE(ListA%RInt)
#endif
          !
          ListATmp=>ListA
          ListA=>ListA%AtmNext
          !
          DEALLOCATE(ListATmp)
          !
       ENDDO
       NULLIFY(List(I)%GoList)
       !
    ENDDO
    !
    DEALLOCATE(List)
    !
  END SUBROUTINE DeAllocList
  !
  !
  SUBROUTINE PrintList(List)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PrintList(List)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CList), DIMENSION(:), POINTER :: List
    !-------------------------------------------------------------------
    TYPE(ANode)              , POINTER :: ListA
    INTEGER                            :: IdxMin,IdxMax,I,J,iPrc,iErr
    !-------------------------------------------------------------------
    !
    NULLIFY(ListA)
    !
    IdxMin=LBOUND(List,DIM=1)
    IdxMax=UBOUND(List,DIM=1)
    !
#ifdef ONX2_PARALLEL
    DO iPrc=0,NPrc-1
       CALL MPI_Barrier(MONDO_COMM,IErr)
       IF(MyId.NE.iPrc) CYCLE
#endif
       !
       DO I=IdxMin,IdxMax
          !
          ListA=>List(I)%GoList
          !
          DO
             IF(.NOT.ASSOCIATED(ListA)) EXIT
             !
#ifdef ONX2_PARALLEL
             WRITE(*,'(4(A,I4))') 'AtC',I,', AtA',ListA%Atom,', NFPair',ListA%NFPair,' MyID',MyID
#else
             WRITE(*,'(3(A,I4))') 'AtC',I,', AtA',ListA%Atom,', NFPair',ListA%NFPair
#endif
             !
#ifdef POINTERS_IN_DERIVED_TYPES
             IF(.NOT.ASSOCIATED(ListA%Indx).OR. &
                    .NOT.ASSOCIATED(ListA%RInt)) STOP 'Array not allocate for the list.'
#else
             IF(.NOT.ALLOCATED(ListA%Indx).OR. &
                    .NOT.ALLOCATED(ListA%RInt)) STOP 'Array not allocate for the list.'
#endif
             !
             DO J=1,SIZE(ListA%Indx,DIM=2)
                WRITE(*,'(3(1X,A,I4),A,E22.15)') 'CFA=',ListA%Indx(1,J),'CFC=',ListA%Indx(2,J), &
                                           'Cell=',ListA%Indx(3,J),' RInt=',ListA%RInt(J)
             ENDDO
             !
             ListA=>ListA%AtmNext
             !
          ENDDO
          !
       ENDDO
       !
#ifdef ONX2_PARALLEL
       WRITE(*,*)
    ENDDO
#endif
    !
  END SUBROUTINE PrintList
  !
  !
  SUBROUTINE PrintMatrix(V,M,N,IOpt,IOut_O,SHFTM_O,SHFTN_O,TEXT_O)
    IMPLICIT NONE
!    ++++ PRINT OUT A RECTANGULAR MATRIX +++++
!
!    V    : MATRIX MxN
!    M    : NUMBER OF ROW
!    N    : NUMBER OF COLUMN
!
!    IOpt = 0; MAX COLUMN = 10
!    IOpt = 1; MAX COLUMN =  7
!    IOpt = 2; MAX COLUMN =  5
!    IOpt = 3; MAX COLUMN =  3
!
    INTEGER                                  :: IMax,IMin,NMAX,I,J
    INTEGER                                  :: IOut,SHFTM,SHFTN
    INTEGER                     , INTENT(IN) :: M,N,IOpt
    INTEGER         , OPTIONAL  , INTENT(IN) :: IOut_O,SHFTM_O,SHFTN_O
    REAL(DOUBLE)                             :: V(M,N)
    CHARACTER(LEN=*), OPTIONAL  , INTENT(IN) :: TEXT_O
!
    IOut = 6
    IF(PRESENT(IOut_O)) IOut = IOut_O
!
    SHFTM = 0
    SHFTN = 0
    IF(PRESENT(SHFTM_O)) SHFTM = SHFTM_O
    IF(PRESENT(SHFTN_O)) SHFTN = SHFTN_O
!
    WRITE(IOut,100)
!
    SELECT CASE(IOPT)
    CASE(0); NMAX = 10
    CASE(1); NMAX =  7
    CASE(2); NMAX =  5
    CASE(3); NMAX =  3
    CASE DEFAULT
       ! WRITE ERROR MESSAGE OR PUT SOMETHING THERE
       RETURN
    END SELECT
!
    IF(N.EQ.0) RETURN
!
    IF(PRESENT(TEXT_O)) WRITE(IOUT,*) TEXT_O
    IMAX = 0
!
    DO WHILE (IMAX.LT.M)
       IMIN = IMAX+1
       IMAX = IMAX+NMAX
       IF(IMAX .GT. M) IMAX = M
       SELECT CASE(IOPT)
       CASE(0); WRITE(IOUT,1000) (I+SHFTN,I=IMIN,IMAX)
       CASE(1); WRITE(IOUT,2000) (I+SHFTN,I=IMIN,IMAX)
       CASE(2); WRITE(IOUT,3000) (I+SHFTN,I=IMIN,IMAX)
       CASE(3); WRITE(IOUT,4000) (I+SHFTN,I=IMIN,IMAX)
       END SELECT
       DO J = 1,M
          SELECT CASE(IOPT)
          CASE(0); WRITE(IOUT,1100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(1); WRITE(IOUT,2100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(2); WRITE(IOUT,3100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(3); WRITE(IOUT,4100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          END SELECT
       ENDDO
       WRITE (IOUT,100)
    ENDDO
    !
    WRITE (IOUT,100)
    !
100 FORMAT('')
1000 FORMAT(6X,10(4X,I3,4X))
1100 FORMAT(I5,1X,10F11.5  )
2000 FORMAT(6X,7(6X,I3,6X ))
2100 FORMAT(I5,1X,7F15.10  )
3000 FORMAT(6X,7(7X,I3,6X ))
3100 FORMAT(I5,1X,7E16.8   )
4000 FORMAT(6X,7(7X,I3,6X ))
4100 FORMAT(I5,1X,7E16.8   )
  END SUBROUTINE PrintMatrix
  !
END MODULE ONX2List
