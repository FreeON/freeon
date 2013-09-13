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
MODULE ONXRng
!H=================================================================================
!H MODULE ONXRng
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB RangeOfExchange
!H
!H  PRIVATE:
!H  o SUB RangeOfExchange_BCSR
!H  o SUB RangeOfExchange_FASTMAT
!H
!H Comments:
!H
!H---------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE Macros
  USE ONXParameters
  USE MemMan
  USE Stats
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  PRIVATE
  !
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
#ifdef ONX2_PARALLEL
  PUBLIC  :: RangeOfExchangeFASTMAT
#else
  PUBLIC  :: RangeOfExchangeBCSR
#endif
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  !
!---------------------------------------------------------------------------------
! MOD PROC DECLARATIONS
!---------------------------------------------------------------------------------
  !
CONTAINS
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE RangeOfExchangeFASTMAT(BSc,GMc,BSp,GMp,D)
!H---------------------------------------------------------------------------------
!H SUBROUTINE RangeOfExchangeFASTMAT(BSc,GMc,BSp,GMp,D)
!H
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET   ), INTENT(IN) :: BSc,BSp
    TYPE(CRDS   ), INTENT(IN) :: GMc,GMp
    TYPE(FASTMAT), POINTER    :: D
    !-------------------------------------------------------------------
    TYPE(FASTMAT), POINTER    :: P
    TYPE(SRST   ), POINTER    :: U
    TYPE(DBL_VECT)            :: MR
    INTEGER                   :: AtA,AtB,AtC,AtD
    INTEGER                   :: Row,Col
    INTEGER                   :: KA,KB,NBFA,NBFB
    REAL(DOUBLE)              :: Ax,Ay,Az,Bx,By,Bz
    REAL(DOUBLE)              :: Cx,Cy,Cz,Dx,Dy,Dz
    REAL(DOUBLE)              :: AB2,AC2,CD2
    !-------------------------------------------------------------------
    !
    ! Set some variables.
    NRows=0
    NCols=0
    NElem=0
    DenRange=0.0D0
    !
    IF(.NOT.ASSOCIATED(D)) THEN
       WRITE(*,*) 'ERR: D is null in RangeOfExchange_FastMat!'
       STOP
    ENDIF
    !
    ! Build the internal linked list.
    CALL FlattenAllRows(D)
    !
    P => D%Next
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       Row = P%Row                     !AtC
       !
       Cx=GMp%Carts%D(1,Row)
       Cy=GMp%Carts%D(2,Row)
       Cz=GMp%Carts%D(3,Row)
       !
       ! Go into the linked list.
       U => P%RowRoot
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          Col = U%L                    !AtD
          IF(Col.EQ.U%R) THEN
             Dx = GMp%Carts%D(1,Col)
             Dy = GMp%Carts%D(2,Col)
             Dz = GMp%Carts%D(3,Col)
             CD2= (Cx-Dx)*(Cx-Dx) + &
                  (Cy-Dy)*(Cy-Dy) + &
                  (Cz-Dz)*(Cz-Dz)
             DenRange = MAX(DenRange,SQRT(CD2)*1.1D0)
          ENDIF
          U => U%Next
       ENDDO
       P => P%Next
    ENDDO
    !
    ! Compute exchange range.
    ONXRange=DenRange+2.0D0*DisRange
    !
    ! Collect all the values from slaves.
    IF(MyID==ROOT) THEN
       CALL New(MR,NPrc-1,0)
       MR%D(:)=Zero
    ENDIF
    CALL Gather(ONXRange,MR)
    IF(MyID==ROOT) THEN
       MatRange=GetMax(NPrc-1,MR,0)
       CALL Delete(MR)
    ENDIF
    !
  END SUBROUTINE RangeOfExchangeFASTMAT
  !
#else
  !
  SUBROUTINE RangeOfExchangeBCSR(BSc,GMc,BSp,GMp,D)
!H---------------------------------------------------------------------------------
!H SUBROUTINE RangeOfExchangeBCSR(BSc,GMc,BSp,GMp,D)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(BSET),INTENT(IN)    :: BSc,BSp
    TYPE(CRDS),INTENT(IN)    :: GMc,GMp
    TYPE(BCSR)               :: D
    !-------------------------------------------------------------------
    INTEGER                  :: AtA,AtB,AtC,AtD
    INTEGER                  :: ri,ci
    INTEGER                  :: KA,KB,NBFA,NBFB
    REAL(DOUBLE)             :: Ax,Ay,Az,Bx,By,Bz
    REAL(DOUBLE)             :: Cx,Cy,Cz,Dx,Dy,Dz
    REAL(DOUBLE)             :: AB2,AC2,CD2
    !-------------------------------------------------------------------
    !
    NRows=0
    NCols=0
    NElem=0
    DenRange=0.0D0
    DO AtC=1,NAtoms
       ri=AtC
       Cx=GMp%Carts%D(1,AtC)
       Cy=GMp%Carts%D(2,AtC)
       Cz=GMp%Carts%D(3,AtC)
       DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1
          AtD=D%ColPt%I(ci)
          Dx=GMp%Carts%D(1,AtD)
          Dy=GMp%Carts%D(2,AtD)
          Dz=GMp%Carts%D(3,AtD)
          CD2= (Cx-Dx)*(Cx-Dx) + &
               (Cy-Dy)*(Cy-Dy) + &
               (Cz-Dz)*(Cz-Dz)
          DenRange=MAX(DenRange,SQRT(CD2)*1.1D0)
       END DO
    END DO
!
    ONXRange=DenRange+2.0D0*DisRange
!
    DO AtA=1,NAtoms
       Ax=GMc%Carts%D(1,AtA)
       Ay=GMc%Carts%D(2,AtA)
       Az=GMc%Carts%D(3,AtA)
       DO ATC=1,NAtoms
          Cx=GMp%Carts%D(1,AtC)
          Cy=GMp%Carts%D(2,AtC)
          Cz=GMp%Carts%D(3,AtC)
          AC2= (Ax-Cx)*(Ax-Cx) + &
               (Ay-Cy)*(Ay-Cy) + &
               (Az-Cz)*(Az-Cz)
          IF(SQRT(AC2).LE.DisRange) THEN
             NRows=NRows+1
             EXIT
          ENDIF
       ENDDO ! AtC
    ENDDO ! AtA
!
! Loop over elements of the exchange matrix
!
    DO AtA=1,NAtoms          !New
       KA=GMc%AtTyp%I(AtA)
       NBFA=BSc%BfKnd%I(KA)
       Ax=GMc%Carts%D(1,AtA)
       Ay=GMc%Carts%D(2,AtA)
       Az=GMc%Carts%D(3,AtA)
       DO AtB=1,NAtoms       !New
          KB=GMc%AtTyp%I(AtB)
          NBFB=BSc%BfKnd%I(KB)
          Bx=GMc%Carts%D(1,AtB)
          By=GMc%Carts%D(2,AtB)
          Bz=GMc%Carts%D(3,AtB)
          AB2= (Ax-Bx)*(Ax-Bx) + &
               (Ay-By)*(Ay-By) + &
               (Az-Bz)*(Az-Bz)
          IF (SQRT(AB2).LE.ONXRange) THEN
             NCols=NCols+1
             NElem=NElem+NBFA*NBFB
          END IF
       END DO ! ci
    END DO ! ri
    NElem=NElem*D%NSMat
    MatRange=ONXRange
  END SUBROUTINE RangeOfExchangeBCSR
  !
#endif
!
END MODULE ONXRng
