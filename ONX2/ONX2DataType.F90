MODULE ONX2DataType
  !
  USE DerivedTypes
  USE ShellPairStruct
  !
  TYPE RNode
     INTEGER              :: iCell
     REAL(DOUBLE)         :: MaxInt
     TYPE(RNode), POINTER :: RNext
  END TYPE RNode
  !
  TYPE ANode
     INTEGER              :: Atom
     REAL(DOUBLE)         :: MaxInt
     TYPE(ANode), POINTER :: AtmNext
     TYPE(RNode), POINTER :: RNext
  END TYPE ANode
  !
  TYPE ANode2
     INTEGER :: Atom
     INTEGER :: NCell
#ifdef POINTERS_IN_DERIVED_TYPES
     INTEGER     , DIMENSION(:), POINTER :: CellIdx
     REAL(DOUBLE), DIMENSION(:), POINTER :: SqrtInt
#else
     INTEGER     , DIMENSION(:), ALLOCATABLE :: CellIdx
     REAL(DOUBLE), DIMENSION(:), ALLOCATABLE :: SqrtInt
#endif
     TYPE(ANode2), POINTER :: AtmNext
  END TYPE ANode2
  !
  TYPE CList2
     TYPE(ANode2), POINTER :: GoList
  END TYPE CList2
  !
  TYPE CList
     TYPE(ANode), POINTER :: GoList
  END TYPE CList
  !
  TYPE AtomPrG
     TYPE(ShellPairG) :: SP
  END TYPE AtomPrG
  !
  INTEGER, PUBLIC :: MaxFuncPerAtmBlk=0
  INTEGER, PUBLIC :: MaxShelPerAtmBlk=0
  !
CONTAINS
  !
  SUBROUTINE GetBufferSize(GM,BS)
    TYPE(CRDS), INTENT(IN) :: GM
    TYPE(BSET), INTENT(IN) :: BS
    INTEGER :: ISize,I
    !
    MaxFuncPerAtmBlk=0
    MaxShelPerAtmBlk=0
    !
    DO I=1,MAXVAL(GM%AtTyp%I(1:NAtoms))
       MaxFuncPerAtmBlk=MAX(MaxFuncPerAtmBlk,BS%BfKnd%I(I))
       MaxShelPerAtmBlk=MAX(MaxShelPerAtmBlk,BS%NCFnc%I(I))
    ENDDO
    !
    !write(*,*) 'max',MAXVAL(GM%AtTyp%I(1:NAtoms)), &
    !     & ' MaxFuncPerAtmBlk',MaxFuncPerAtmBlk, &
    !     & ' MaxShelPerAtmBlk',MaxShelPerAtmBlk
  END SUBROUTINE GetBufferSize
  !
END MODULE ONX2DataType



