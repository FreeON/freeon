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
     INTEGER     , DIMENSION(:), POINTER :: CellIdx
     REAL(DOUBLE), DIMENSION(:), POINTER :: SqrtInt
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
  TYPE ShellPair
     INTEGER :: N1
     INTEGER :: N2
     REAL(DOUBLE) :: At1x,At1y,At1z
     REAL(DOUBLE) :: At2x,At2y,At2z
     REAL(DOUBLE) :: R12
     REAL(DOUBLE), DIMENSION(10) :: Expt1
     REAL(DOUBLE), DIMENSION(10) :: Expt2
     REAL(DOUBLE), DIMENSION(10) :: Coef1
     REAL(DOUBLE), DIMENSION(10) :: Coef2
  END TYPE ShellPair
#endif

  INTEGER, PUBLIC :: MaxFuncPerAtmBlk=0
  INTEGER, PUBLIC :: MaxShelPerAtmBlk=0

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



