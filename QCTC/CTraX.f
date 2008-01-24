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
C    ECONOMIZED (FULLY FACTORED) CONTRACTION OF SP MULTIPOLE TENSORS
C    Author: Matt Challacombe
C==============================================================================
      SUBROUTINE CTraX77(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
C
      REAL*8 Cp(0:*),Sp(0:*),Cpq(0:*),Spq(0:*),Cq(0:*),Sq(0:*)
      REAL*8 TMP,SGN,TMPC,TWO,TMPS,SGNL,SGNLM
C
      INTEGER L,LP,LQ,LL,M,LDX,K,KK,LLKK,LLKM,N,KDX,LKDX,ID
C
      ID(L)=L*(L+1)/2
C
C     ZERO
C
      IF(LP.EQ.0)THEN
         Tmp=0.0D0
         DO l=0,LQ
            ll=ID(l)
            Cp(0)=Cp(0)+Cpq(ll)*Cq(ll)+Spq(ll)*Sq(ll)            
            DO m=1,l
              ldx=ll+m
              Tmp=Tmp+Cpq(ldx)*Cq(ldx)+Spq(ldx)*Sq(ldx)
           ENDDO
        ENDDO
        Cp(0)=Cp(0)+2.0D0*Tmp
        RETURN
      ENDIF
C
C     ONE
C
      Sgn=1.0D0
      DO l=0,LP
         ll=ID(l)
         TmpC=0.0D0
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            TmpC=TmpC+Cpq(llkk)*Cq(kk)
         ENDDO
         Cp(ll)=Cp(ll)+Sgn*TmpC
         Sgn=-Sgn
      ENDDO
C
C     TWO
C
      Two=-2.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO m=1,l                   
               ldx=ll+m
               llkm=llkk+m
               Cp(ldx)=Cp(ldx)+Two*Cpq(llkm)*Cq(kk)
               Sp(ldx)=Sp(ldx)+Two*Spq(llkm)*Cq(kk)
            ENDDO
            
         ENDDO
         Two=-Two
      ENDDO 
C
C     THREE
C
      Two=2.0D0
      DO l=0,LP
         ll=ID(l)
         TmpC=0.0D0
         DO k=0,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO n=1,k
               kdx=kk+n
               lkdx=llkk+n
               TmpC=TmpC+Cpq(lkdx)*Cq(kdx)+Spq(lkdx)*Sq(kdx)
            ENDDO
         ENDDO
         Cp(ll)=Cp(ll)+Two*TmpC
         Two=-Two
      ENDDO 
C
C     FOUR
C
      Two=-2.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=1,LQ
            kk=ID(k)
            llkk=ID(l+k)
            DO m=1,l                   
               ldx=ll+m
               llkm=llkk+m
               TmpC =0.0D0
               TmpS =0.0D0
               DO n=1,k
                  kdx=kk+n
                  lkdx=llkm+n
                  TmpC=TmpC+Spq(lkdx)*Sq(kdx)+Cpq(lkdx)*Cq(kdx)
                  TmpS=TmpS-Cpq(lkdx)*Sq(kdx)+Spq(lkdx)*Cq(kdx)
               ENDDO
               Cp(ldx)=Cp(ldx)+Two*TmpC
               Sp(ldx)=Sp(ldx)+Two*TmpS
            ENDDO
         ENDDO
         Two=-Two
      ENDDO 
C
C     FIVE
C
      SgnL=1.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=1,LQ
            kk=ID(k)
            llkk=ID(l+k)
            SgnLM=SgnL
            DO m=1,l      
               ldx=ll+m
               llkm=llkk+m
               Two=2.0D0*SgnL   
               DO n=1,MIN(m,k)
                  kdx=kk+n
                  lkdx=llkm-n
                  Cp(ldx)=Cp(ldx)+Two*Cpq(lkdx)*Cq(kdx)
                  Sp(ldx)=Sp(ldx)+Two*Cpq(lkdx)*Sq(kdx)
                  Two=-Two
               ENDDO
               llkm=llkk-m
               TmpC=0.0D0
               TmpS=0.0D0
               DO n=m+1,k
                  kdx=kk+n
                  lkdx=llkm+n
                  TmpC=TmpC+Cpq(lkdx)*Cq(kdx)
                  TmpS=TmpS+Cpq(lkdx)*Sq(kdx)
               ENDDO
               Cp(ldx)=Cp(ldx)+SgnLM*2.0D0*TmpC
               Sp(ldx)=Sp(ldx)+SgnLM*2.0D0*TmpS
               SgnLM=-SgnLM
            ENDDO
         ENDDO
         SgnL=-SgnL
      ENDDO 
C
      SgnL=1.0D0
      DO l=1,LP
         ll=ID(l)
         DO k=1,LQ
            kk=ID(k)
            llkk=ID(l+k)
            SgnLM=SgnL
            DO m=1,l      
               ldx=ll+m
               llkm=llkk+m
               Two=-2.0D0*SgnL               
               DO n=1,MIN(m-1,k)
                  kdx=kk+n
                  lkdx=llkm-n            
                  Cp(ldx)=Cp(ldx)+Two*Spq(lkdx)*Sq(kdx)
                  Sp(ldx)=Sp(ldx)-Two*Spq(lkdx)*Cq(kdx)
                  Two=-Two
               ENDDO
               llkm=llkk-m
               TmpC =0.0D0
               TmpS =0.0D0
               DO n=m+1,k
                  kdx=kk+n
                  lkdx=llkm+n
                  TmpC=TmpC+Spq(lkdx)*Sq(kdx)
                  TmpS=TmpS+Spq(lkdx)*Cq(kdx)
               ENDDO
               Cp(ldx)=Cp(ldx)+SgnLM*2.0D0*TmpC
               Sp(ldx)=Sp(ldx)-SgnLM*2.0D0*TmpS          
               SgnLM=-SgnLM
            ENDDO
         ENDDO
         SgnL=-SgnL
      ENDDO 
C
      RETURN
      END
