MODULE ONX2DataType
  !
  USE DerivedTypes
  USE ShellPairStruct
  !
  TYPE ANode
     INTEGER :: Atom
     INTEGER :: NFPair
#ifdef POINTERS_IN_DERIVED_TYPES
     INTEGER     , DIMENSION(:,:), POINTER :: Indx
     REAL(DOUBLE), DIMENSION(:  ), POINTER :: RInt
#else
     INTEGER     , DIMENSION(:,:), ALLOCATABLE :: Indx
     REAL(DOUBLE), DIMENSION(:  ), ALLOCATABLE :: RInt
#endif
     TYPE(ANode), POINTER :: AtmNext
  END TYPE ANode
  !
  TYPE CList
     TYPE(ANode), POINTER :: GoList
  END TYPE CList
  !
  INTEGER, PUBLIC :: MaxFuncPerAtmBlk=0
  INTEGER, PUBLIC :: MaxShelPerAtmBlk=0
  !
CONTAINS
  !
!!$  SUBROUTINE GetBufferSize(GM,BS)
!!$    TYPE(CRDS) :: GM
!!$    TYPE(BSET) :: BS
!!$    INTEGER :: ISize,I
!!$    !
!!$    MaxFuncPerAtmBlk=0
!!$    MaxShelPerAtmBlk=0
!!$    !
!!$    !write(*,*) 'NKind=',BS%NKind,' MAXVAL(GM%AtTyp%I(1:NAtoms))=',MAXVAL(GM%AtTyp%I(1:NAtoms))
!!$    if(BS%NKind.ne.MAXVAL(GM%AtTyp%I(1:NAtoms))) stop 'GetBufferSize: problem'
!!$    DO I=1,BS%NKind !MAXVAL(GM%AtTyp%I(1:NAtoms))
!!$       MaxFuncPerAtmBlk=MAX(MaxFuncPerAtmBlk,BS%BfKnd%I(I))
!!$       MaxShelPerAtmBlk=MAX(MaxShelPerAtmBlk,BS%NCFnc%I(I))
!!$    ENDDO
!!$    !
!!$    !write(*,*) 'max',MAXVAL(GM%AtTyp%I(1:NAtoms)), &
!!$    !     & ' MaxFuncPerAtmBlk',MaxFuncPerAtmBlk, &
!!$    !     & ' MaxShelPerAtmBlk',MaxShelPerAtmBlk
!!$  END SUBROUTINE GetBufferSize
  !
  SUBROUTINE GetBufferSize(GMc,BSc,GMp,BSp)
    TYPE(CRDS) :: GMc,GMp
    TYPE(BSET) :: BSc,BSp
    INTEGER :: I
    !
    MaxFuncPerAtmBlk=0
    MaxShelPerAtmBlk=0
    !
    if(BSc%NKind.ne.MAXVAL(GMc%AtTyp%I(1:NAtoms))) stop 'GetBufferSize: problem'
    DO I=1,BSc%NKind !MAXVAL(GMc%AtTyp%I(1:NAtoms))
       MaxFuncPerAtmBlk=MAX(MaxFuncPerAtmBlk,BSc%BfKnd%I(I))
       MaxShelPerAtmBlk=MAX(MaxShelPerAtmBlk,BSc%NCFnc%I(I))
    ENDDO
    !
    if(BSp%NKind.ne.MAXVAL(GMp%AtTyp%I(1:NAtoms))) stop 'GetBufferSize: problem'
    DO I=1,BSp%NKind !MAXVAL(GMp%AtTyp%I(1:NAtoms))
       MaxFuncPerAtmBlk=MAX(MaxFuncPerAtmBlk,BSp%BfKnd%I(I))
       MaxShelPerAtmBlk=MAX(MaxShelPerAtmBlk,BSp%NCFnc%I(I))
    ENDDO
    !
    !write(*,*) &
    !     & ' MaxFuncPerAtmBlk',MaxFuncPerAtmBlk, &
    !     & ' MaxShelPerAtmBlk',MaxShelPerAtmBlk
  END SUBROUTINE GetBufferSize
  !
END MODULE ONX2DataType



