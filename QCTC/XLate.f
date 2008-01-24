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
C    TRANSLATION OF SP MULTIPOLE TENSORS
C    Author: Matt Challacombe
C==============================================================================
      SUBROUTINE XLate77(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
C
      REAL*8 Cp(0:*),Sp(0:*),Cpq(0:*),Spq(0:*),Cq(0:*),Sq(0:*)
      REAL*8 CMN,SMN,CN,SN
C
      INTEGER L,LP,LL,M,LDX,K,LK,KK,LKLK,NSTART,NSTOP,
     >        N,NABS,MNABS,KDX,LKDX,ID,LQ
C
      ID(L)=L*(L+1)/2
C
      DO l=0,LP
         ll=ID(l)
         DO m=0,l
            ldx=ll+m
            DO k=0,MIN(l,LQ)
               lk=l-k
               kk=ID(k)
               lklk=ID(lk)
               nStart=MAX(-k, k-l+m)
               nStop =MIN( k,-k+l+m)
               DO n=nStart,nStop
                  nabs=ABS(n)
                  mnabs=ABS(m-n)
                  kdx=kk+nabs                          
                  lkdx=lklk+mnabs
                  IF(m-n.LT.0)THEN
                     cmn=(-1.0D0)**mnabs
                     ! smn=(-1.0D0)**(mnabs+1)
                     smn = -cmn
                  ELSE
                     cmn=1.0D0
                     smn=1.0D0
                  ENDIF
                  IF(n.LT.0)THEN
                     cn=(-1.0D0)**nabs
                     ! sn=(-1.0D0)**(nabs+1)
                     sn = -cn
                  ELSE
                     cn=1.0D0
                     sn=1.0D0
                  ENDIF
                  Cp(ldx)=Cp(ldx)+cn*Cq(kdx)*cmn*Cpq(lkdx)
     >                           -sn*Sq(kdx)*smn*Spq(lkdx)        
                  Sp(ldx)=Sp(ldx)+cn*Cq(kdx)*smn*Spq(lkdx)
     >                           +sn*Sq(kdx)*cmn*Cpq(lkdx)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
