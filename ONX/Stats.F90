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
MODULE Stats
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   IMPLICIT NONE

   INTERFACE GetAbsMax
     MODULE PROCEDURE GetAbsMax_DBL_VECT, GetAbsMax_INT_VECT
   END INTERFACE

   INTERFACE GetMax
     MODULE PROCEDURE GetMax_INT_VECT, GetMax_DBL_VECT
   END INTERFACE

   INTERFACE GetMin
     MODULE PROCEDURE GetMin_INT_VECT, GetMin_DBL_VECT
   END INTERFACE

   CONTAINS

   FUNCTION GetAbsMax_DBL_VECT(N,A)
     IMPLICIT NONE
     INTEGER,INTENT(IN)         :: N
     TYPE(DBL_VECT), INTENT(IN) :: A
     REAL(DOUBLE)               :: GetAbsMax_DBL_VECT
     INTEGER                    :: I
     REAL(DOUBLE)               :: X
     X=0.0D0
     DO I=1,N
       X=MAX(X,ABS(A%D(I)))
     ENDDO
     GetAbsMax_DBL_VECT=X
   END FUNCTION GetAbsMax_DBL_VECT

   FUNCTION GetAbsMax_INT_VECT(N,A)
     IMPLICIT NONE
     INTEGER,INTENT(IN)         :: N
     TYPE(INT_VECT), INTENT(IN) :: A
     REAL(DOUBLE)               :: GetAbsMax_INT_VECT
     INTEGER                    :: I,J
     J=0
     DO I=1,N
       J=MAX(J,ABS(A%I(I)))
     ENDDO
     GetAbsMax_INT_VECT=J
   END FUNCTION GetAbsMax_INT_VECT

   FUNCTION GetMax_INT_VECT(N,A,M_O)
     INTEGER, INTENT(IN)           :: N
     TYPE(INT_VECT), INTENT(IN)    :: A
     INTEGER, OPTIONAL, INTENT(IN) :: M_O
     INTEGER                       :: GetMax_INT_VECT,I,J,M
     M=0; J=1; IF(PRESENT(M_O)) J=M_O
     DO I=J,N
       M=MAX(M,A%I(I))
      END DO
     GetMax_INT_VECT=M
   END FUNCTION GetMax_INT_VECT

   FUNCTION GetMax_DBL_VECT(N,A,M_O)
     INTEGER, INTENT(IN)           :: N
     TYPE(DBL_VECT), INTENT(IN)    :: A
     INTEGER, OPTIONAL, INTENT(IN) :: M_O
     REAL(DOUBLE)                  :: GetMax_DBL_VECT,M
     INTEGER                       :: I,J
     M=0.0D0; J=1; IF(PRESENT(M_O)) J=M_O
     DO I=J,N
       M=MAX(M,A%D(I))
     END DO
     GetMax_DBL_VECT=M
   END FUNCTION GetMax_DBL_VECT

   FUNCTION GetMin_INT_VECT(N,A,M_O)
     INTEGER, INTENT(IN)           :: N
     TYPE(INT_VECT), INTENT(IN)    :: A
     INTEGER, OPTIONAL, INTENT(IN) :: M_O
     INTEGER                       :: GetMin_INT_VECT,I,J,M
     M=0; J=1; IF(PRESENT(M_O)) J=M_O
      DO I=J,N
       M=MIN(M,A%I(I))
     END DO
     GetMin_INT_VECT=M
   END FUNCTION GetMin_INT_VECT

   FUNCTION GetMin_DBL_VECT(N,A,M_O)
     INTEGER, INTENT(IN)           :: N
     TYPE(DBL_VECT), INTENT(IN)    :: A
     INTEGER, OPTIONAL, INTENT(IN) :: M_O
     REAL(DOUBLE)                  :: GetMin_DBL_VECT,M
     INTEGER                       :: I,J
     M=0.0D0; J=1; IF(PRESENT(M_O)) J=M_O
      DO I=J,N
       M=MIN(M,A%D(I))
     END DO
     GetMin_DBL_VECT=M
   END FUNCTION GetMin_DBL_VECT

END MODULE
