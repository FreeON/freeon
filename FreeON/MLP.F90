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
MODULE MLP
!H=================================================================================
!H
!H
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Parse
  USE LinAlg
  IMPLICIT NONE
  PRIVATE
  !
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: MLPDriver
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: MLPDpol
  !
CONTAINS
  !
  SUBROUTINE MLPDriver(D,GM,Args)
!H---------------------------------------------------------------------------------
!H SUBROUTINE MLPDriver()
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL_PRP
    TYPE(DBCSR), INTENT(IN) :: D
#else
    TYPE(BCSR ), INTENT(IN) :: D
#endif
    TYPE(CRDS ), INTENT(IN) :: GM
    TYPE(ARGMT), INTENT(IN) :: Args
    !-------------------------------------------------------------------
    TYPE(DBL_VECT)          :: COrig
    !-------------------------------------------------------------------
    !
    ! Get coordinate origine.
    !TODO
    !
    ! Compute dipole moment.
    CALL MLPDpol(D,GM,COrig,Args)
    !
  END SUBROUTINE MLPDriver
  !
  SUBROUTINE MLPDpol(D,GM,COrig,Args)
!H---------------------------------------------------------------------------------
!H SUBROUTINE MLPDpol()
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL_PRP
    TYPE(DBCSR)   , INTENT(IN)     :: D
#else
    TYPE(BCSR )   , INTENT(IN)     :: D
#endif
    TYPE(CRDS )   , INTENT(IN)     :: GM
    TYPE(DBL_VECT), INTENT(IN)     :: COrig
    TYPE(ARGMT)   , INTENT(IN)     :: Args
    !-------------------------------------------------------------------
#ifdef PARALLEL_PRP
    TYPE(DBCSR)                    :: Tmp1,Tmp2
#else
    TYPE(BCSR )                    :: Tmp1
#endif
    INTEGER                        :: iXYZ,iAtom
    REAL(DOUBLE)                   :: ZNuc
    REAL(DOUBLE), DIMENSION(3)     :: MltE,MltN
    CHARACTER(LEN=1), DIMENSION(3), PARAMETER :: Cart=(/'x','y','z'/)
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: MltMessage
    !-------------------------------------------------------------------
    !
    CALL New(Tmp1)
#ifdef PARALLEL_PRP
    CALL New(Tmp2)
#endif
    !
    ! Compute electronic dipole.
    MltE(:)=Zero
    DO iXYZ=1,3
       CALL Get(Tmp1,TrixFile('D'//Cart(iXYZ),Args))
#ifdef PARALLEL_PRP
       CALL Multiply(D,Tmp1,Tmp2)
       MltE(iXYZ)=Two*Trace(Tmp2)
#else
       MltE(iXYZ)=Two*Trace(D,Tmp1)
#endif
       !WRITE(*,'(A,A,A,3E20.12)') ' D',Cart(iXYZ),' = ',Mlt(iXYZ)
    ENDDO
    !
    ! Compute nuclear dipole.
    MltN(:)=Zero
    DO iAtom=1,NAtoms
       ZNuc=DBLE(GM%AtNum%D(iAtom))
       DO iXYZ=1,3
          MltN(iXYZ)=MltN(iXYZ)+ZNuc*(GM%Carts%D(iXYZ,iAtom)-COrig%D(iXYZ))
       ENDDO
    ENDDO
    !
    ! Print dipole.
    MltMessage=""
    MltMessage=RTRN//' Electronic Dipole Moment: X ='//TRIM(DblToMedmChar(MltE(1)))// &
                                              ', Y ='//TRIM(DblToMedmChar(MltE(2)))// &
                                              ', Z ='//TRIM(DblToMedmChar(MltE(3)))// &
               RTRN//' Nuclear Dipole Moment   : X ='//TRIM(DblToMedmChar(MltN(1)))// &
                                              ', Y ='//TRIM(DblToMedmChar(MltN(2)))// &
                                              ', Z ='//TRIM(DblToMedmChar(MltN(3)))// &
               RTRN//' Total Dipole Moment     : X ='//TRIM(DblToMedmChar(MltE(1)+MltN(1)))// &
                                              ', Y ='//TRIM(DblToMedmChar(MltE(2)+MltN(2)))// &
                                              ', Z ='//TRIM(DblToMedmChar(MltE(3)+MltN(3)))
    !
#ifdef PARALLEL_PRP
    IF(MyId==ROOT)THEN
#endif
       CALL OpenASCII(OutFile,Out)
       WRITE(*  ,* )TRIM(MltMessage)
       WRITE(Out,* )TRIM(MltMessage)
       CLOSE(Out)
#ifdef PARALLEL_PRP
    ENDIF
#endif
    !
    ! Delete some stuff.
    CALL Delete(Tmp1)
#ifdef PARALLEL_PRP
    CALL Delete(Tmp2)
#endif
    !
  END SUBROUTINE MLPDpol
  !
END MODULE MLP
