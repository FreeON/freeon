C------------------------------------------------------------------------------
C    This code is part of the FreeON suite of programs for linear scaling
C    electronic structure theory and ab initio molecular dynamics.
C
C    Copyright (2004). The Regents of the University of California. This
C    material was produced under U.S. Government contract W-7405-ENG-36
C    for Los Alamos National Laboratory, which is operated by the University
C    of California for the U.S. Department of Energy. The U.S. Government has
C    rights to use, reproduce, and distribute this software.  NEITHER THE
C    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
C    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by the
C    Free Software Foundation; either version 2 of the License, or (at your
C    option) any later version. Accordingly, this program is distributed in
C    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
C    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C    PURPOSE. See the GNU General Public License at www.gnu.org for details.
C
C    While you may do as you like with this software, the GNU license requires
C    that you clearly mark derivative software.  In addition, you are encouraged
C    to return derivative works to the FreeON group for review, and possible
C    dissemination in future releases.
C------------------------------------------------------------------------------

C===================================================================
C     Transform a BCSR matrix to a dense matrix
C===================================================================
      SUBROUTINE BCSR_TO_DENS(NRow,NCol,NBasF,NSMat,NAtoms,NBlks,NNon0,
     >                        MSiz,OffS,A,MTrix,RowPt,ColPt,BlkPt)
      IMPLICIT NONE
      INTEGER NRow,NCol,NBasF,NSMat,NAtoms,NBlks,NNon0
      REAL*8  A(NRow,NCol),MTrix (NNon0)
      INTEGER RowPt(NAtoms+1),ColPt(NBlks),BlkPt(NBlks*NSMat), 
     >        MSiz(NAtoms),OffS(NAtoms),
     >        I,J,K,L,P,S,II,JJ,II0,JJ0,JP,MA,NA,iS,iSMat

      !write(*,*) 'NRow',NRow
      !write(*,*) 'NCol',NCol
      !write(*,*) 'NBasF',NBasF
      !write(*,*) 'NSMat',NSMat
      !write(*,*) 'NAtoms',NAtoms
      !write(*,*) 'NBlks',NBlks
      !write(*,*) 'NNon0',NNon0
      !write(*,*) 'MSiz',MSiz
      !write(*,*) 'OffS',OffS
      !write(*,*) 'MTrix',MTrix
      !write(*,*) 'RowPt',RowPt
      !write(*,*) 'ColPt',ColPt
      !write(*,*) 'BlkPt',BlkPt

      CALL DBL_VECT_EQ_DBL_SCLR(NRow*NCol,A(1,1),0.0D0)

      DO 100 I=1,NAtoms
         MA=MSiz(I)
         II0=OffS(I)
         DO 200 JP=RowPt(I),RowPt(I+1)-1
            J=ColPt(JP)
            NA=MSiz(J)
            JJ0=OffS(J)
            DO 300 iSMat=1,NSMat
ccc               iS=NSMat*(JP-1)+iSMat
c               write(*,*) iS
c               IF(BlkPt(iS+1)-BlkPt(iS).LT.1) GOTO 300
ccc               P=BlkPt(iS)
               P=BlkPt(JP)+(iSMat-1)*MA*NA
               IF(iSMat.EQ.1)THEN
                  II=II0
                  JJ=JJ0
               ELSEIF(iSMat.EQ.2)THEN
                  II=II0
                  JJ=JJ0+NBasF
               ELSEIF(iSMat.EQ.3)THEN
                  II=II0+NBasF
                  JJ=JJ0
               ELSEIF(iSMat.EQ.4)THEN
                  II=II0+NBasF
                  JJ=JJ0+NBasF
               ELSE
                  STOP 'BCSR_TO_DENS'
               ENDIF
               DO 400 K=0,NA-1
                  DO 500 L=0,MA-1
                     A(II+L,JJ+K)=MTrix(P+K*MA+L)
 500              CONTINUE
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
      RETURN
      END
C==========================================================
C     Add a double vector to a double vector
C==========================================================
      SUBROUTINE DBL_VECT_PLS_DBL_VECT(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=A(I)+B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=A(I )+B(I )
         A(I1)=A(I1)+B(I1)
         A(I2)=A(I2)+B(I2)
         A(I3)=A(I3)+B(I3)
         A(I4)=A(I4)+B(I4)
         A(I5)=A(I5)+B(I5)
         A(I6)=A(I6)+B(I6)
      ENDDO
      RETURN
      END
C==========================================================
C     Copy a double vector to a double vector
C==========================================================
      SUBROUTINE DBL_VECT_EQ_DBL_VECT(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B(I )
         A(I1)=B(I1)
         A(I2)=B(I2)
         A(I3)=B(I3)
         A(I4)=B(I4)
         A(I5)=B(I5)
         A(I6)=B(I6)
      ENDDO
      RETURN
      END
C==========================================================
C     Copy an integer vector to an integer vector
C==========================================================
      SUBROUTINE INT_VECT_EQ_INT_VECT(N,A,B)
      IMPLICIT NONE
      INTEGER N
      INTEGER A(N),B(N)
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=B(I)
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B(I )
         A(I1)=B(I1)
         A(I2)=B(I2)
         A(I3)=B(I3)
         A(I4)=B(I4)
         A(I5)=B(I5)
         A(I6)=B(I6)
      ENDDO
      RETURN
      END
C==========================================================
C     Initialize a double vector with a double scalar
C==========================================================
      SUBROUTINE DBL_VECT_EQ_DBL_SCLR(N,A,B)
      IMPLICIT NONE
      INTEGER N
      REAL*8  A(N),B
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=B
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B
         A(I1)=B
         A(I2)=B
         A(I3)=B
         A(I4)=B
         A(I5)=B
         A(I6)=B
      ENDDO
      RETURN
      END
C==========================================================
C     Initialize an integer vector with an integer scalar
C==========================================================
      SUBROUTINE INT_VECT_EQ_INT_SCLR(N,A,B)
      IMPLICIT NONE
      INTEGER N
      INTEGER A(N),B
      INTEGER M,I,I1,I2,I3,I4,I5,I6
C----------------------------------------------
      M=MOD(N,7)
      DO I=1,M
         A(I)=B
      ENDDO
      IF(N.LT.7)RETURN
      DO I=M+1,N,7
         I1=I+1
         I2=I+2
         I3=I+3
         I4=I+4
         I5=I+5
         I6=I+6
         A(I )=B
         A(I1)=B
         A(I2)=B
         A(I3)=B
         A(I4)=B
         A(I5)=B
         A(I6)=B
      ENDDO
      RETURN
      END


