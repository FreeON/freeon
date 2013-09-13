      SUBROUTINE CTraX77P(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
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
      WRITE(*,*)' CP = ',CP(0)
      IF(LP.EQ.0)THEN
         Tmp=0.0D0
         DO l=0,LQ
            ll=ID(l)
            Cp(0)=Cp(0)+Cpq(ll)*Cq(ll)+Spq(ll)*Sq(ll)
            WRITE(*,22)l,0,Cp(0),Cpq(ll),Cq(ll),Spq(ll),Sq(ll)
 22         FORMAT(I2,", ",I2,5(D12.6,", "))
            DO m=1,l
              ldx=ll+m
              Tmp=Tmp+Cpq(ldx)*Cq(ldx)+Spq(ldx)*Sq(ldx)
           ENDDO
            WRITE(*,22)l,m,Cp(0)+2D0*Tmp
        ENDDO
        Cp(0)=Cp(0)+2.0D0*Tmp
        RETURN
      ENDIF
      RETURN
      END

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
