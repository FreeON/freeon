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
