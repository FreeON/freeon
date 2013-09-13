C------------------------------------------------------------------------------
C    This code is part of the MondoSCF suite of programs for linear scaling
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
C    to return derivative works to the MondoSCF group for review, and possible
C    disemination in future releases.
C------------------------------------------------------------------------------
C
      FUNCTION DGetAbsMax(N,A)
      IMPLICIT NONE
      INTEGER N
      REAL*8 A(N)
      REAL*8 DGetAbsMax
      INTEGER I
      DGetAbsMax=0.0D0
      DO I=1,N
         DGetAbsMax=MAX(DGetAbsMax,ABS(A(I)))
      ENDDO
      RETURN
      END
C
      SUBROUTINE GetAbsDenBlk2(D,NR,NC,DM,NFR,NFC,IRB,IRE,ICB,ICE)
      IMPLICIT NONE
      INTEGER NR,NC,NFR,NFC
      INTEGER IRB(NFR),IRE(NFR),ICB(NFC),ICE(NFC)
      REAL*8 D(NR,NC),DM(NFR,NFC)
      INTEGER INFC,INFR,I,J,IRB1,ICB1,STDR,STDC
      REAL*8 TMP
      ICB1=0
      DO INFC=1,NFC
         IRB1=0
         STDC=ICE(INFC)-ICB(INFC)+1
         DO INFR=1,NFR
            STDR=IRE(INFR)-IRB(INFR)+1
            TMP=0.0D0
            DO J=ICB1+1,ICB1+STDC
               DO I=IRB1+1,IRB1+STDR
                  TMP=MAX(TMP,ABS(D(I,J)))
               ENDDO
            ENDDO
            DM(INFR,INFC)=TMP
            IRB1=IRB1+STDR
         ENDDO
         ICB1=ICB1+STDC
      ENDDO
      RETURN
      END
C
      SUBROUTINE GetAbsDenBlk(D,NR,NC,NS,DM,NFR,NFC,IRB,IRE,ICB,ICE)
      IMPLICIT NONE
      INTEGER NR,NC,NS,NFR,NFC
      INTEGER IRB(NFR),IRE(NFR),ICB(NFC),ICE(NFC)
      REAL*8 D(NR,NC*NS),DM(NFC,NFR)
      INTEGER INFC,INFR,I,J,IS,N,IRB1,ICB1,STDR,STDC
      REAL*8 TMP
      ICB1=0
      DO INFC=1,NFC
         STDC=ICE(INFC)-ICB(INFC)+1
         IRB1=0
         DO INFR=1,NFR
            STDR=IRE(INFR)-IRB(INFR)+1
            TMP=0.0D0
            DO IS=1,NS
               N=(IS-1)*NC
               DO J=ICB1+1,ICB1+STDC
               DO I=IRB1+1,IRB1+STDR
                  TMP=MAX(TMP,ABS(D(I,J+N)))
               ENDDO
               ENDDO
            ENDDO
            DM(INFC,INFR)=TMP
            IRB1=IRB1+STDR
         ENDDO
         ICB1=ICB1+STDC
      ENDDO
      RETURN
      END
C
      FUNCTION IBinSrch(IVec,IVal,NDim)
      IMPLICIT NONE
      INTEGER IBinSrch
      INTEGER IVal,NDim
      INTEGER IVec(NDim)
      INTEGER IMin,IMax,IMid
      IMin=0
      IMax=NDim+1
      IBinSrch=-1
      DO WHILE(IMax-IMin.GT.1)
         IMid=(IMax+IMin)/2
         IF(IVal.GT.IVec(IMid)) THEN
            IMin=IMid
         ELSEIF(IVal.LT.IVec(IMid)) THEN
            IMax=IMid
         ELSE
            IBinSrch=IMid
            RETURN
         ENDIF
      ENDDO
      RETURN
      END
C
      SUBROUTINE XPose1C(M,N,A,B)
      IMPLICIT NONE
      INTEGER I,J,M,N
      REAL*8 A(M,N)
      REAL*8 B(N,M)
      DO I=1,N
         DO J=I+1,M
            B(I,J)=A(J,I)
         ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE XPose2C(M,N,A,B)
      IMPLICIT NONE
      INTEGER I,J,M,N
      REAL*8 A(M,N)
      REAL*8 B(N,M)
      DO I=1,N
         DO J=1,M
            B(I,J)=A(J,I)
         ENDDO
      ENDDO
      RETURN
      END
C
