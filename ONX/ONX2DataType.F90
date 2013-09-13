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
!!$    !       ' MaxFuncPerAtmBlk',MaxFuncPerAtmBlk, &
!!$    !       ' MaxShelPerAtmBlk',MaxShelPerAtmBlk
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
    !       ' MaxFuncPerAtmBlk',MaxFuncPerAtmBlk, &
    !       ' MaxShelPerAtmBlk',MaxShelPerAtmBlk
  END SUBROUTINE GetBufferSize
  !
END MODULE ONX2DataType
