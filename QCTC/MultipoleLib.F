C_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_C
C                                                       C
C     Matts Multipole Manipulation Library (MMML)       C 
C                                                       C
C     Matt Challacombe, 1997                            C
C                                                       C
C_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_C
C
C
C---------------------------------------------------------------
C     Translate the IrRegular Momments M_lm[Q-P] 
C     and Contract with Cq
C---------------------------------------------------------------
      SUBROUTINE CTraX(MaxL,LP,LQ,LPQ,LenP,LenQ,LenPQ,
     $                 PQx,PQy,PQz,Cp,Sp,Cpq,Spq,Cq,Sq,
     $                 LegendreP,Cosine,Sine,RToTh,
     $                 FactMlm0,FactMlm1,FactMlm2)
C
      IMPLICIT INTEGER(i-n)
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 Cp(0:LenP)
      REAL*8 Sp(0:LenP)      
      REAL*8 Cq(0:LenQ)
      REAL*8 Sq(0:LenQ)
      REAL*8 Cpq(0:LenPQ)
      REAL*8 Spq(0:LenPQ)
C
      REAL*8 FactMlm0( 0:MaxL)
      REAL*8 FactMlm1( 0:MaxL)
      REAL*8 FactMlm2( 0:MaxL,0:MaxL)
      REAL*8 LegendreP(0:MaxL,0:MaxL)
      REAL*8 Cosine(   0:MaxL)
      REAL*8 Sine(     0:MaxL)
      REAL*8 RToTh(    0:MaxL)
C
      ID(L)=L*(L+1)/2
C
      CALL IrRegular(MaxL,LPQ,LenPQ,PQx,PQy,PQz,Cpq,Spq,
     $               LegendreP,Cosine,Sine,RToTh,
     $               FactMlm0,FactMlm1,FactMlm2)
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
c
c      CALL PrntSH(LP,LenP,Cp,Sp)
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
C     SIX
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
C     LATER ON 
C
      RETURN
      END
C---------------------------------------------------------------
C     Translate the Regular Momments O_lm[r-P]
C---------------------------------------------------------------
      SUBROUTINE XLate(MaxL,LP,LQ,LPQ,LenP,LenQ,LenPQ,
     $                 PQx,PQy,PQz,Cp,Sp,Cpq,Spq,Cq,Sq,
     $                 LegendreP,Cosine,Sine,RToTh,
     $                 FactOlm0,FactOlm1,FactOlm2)
C
      IMPLICIT INTEGER(i-n)
      IMPLICIT REAL*8 (a-h,o-z)
C
      REAL*8 Cp( 0:LenP)
      REAL*8 Sp( 0:LenP)      
      REAL*8 Cpq(0:LenPQ)
      REAL*8 Spq(0:LenPQ)
      REAL*8 Cq( 0:LenQ)
      REAL*8 Sq( 0:LenQ)
C
      REAL*8 LegendreP(0:MaxL,0:MaxL)
      REAL*8 Cosine(   0:MaxL)
      REAL*8 Sine(     0:MaxL)
      REAL*8 RToTh(    0:MaxL)
      REAL*8 FactOlm0( 0:MaxL)
      REAL*8 FactOlm1( 0:MaxL)
      REAL*8 FactOlm2( 0:MaxL,0:MaxL)
C
      ID(L)=L*(L+1)/2
C
      CALL Regular(MaxL,LPQ,LenPQ,PQx,PQy,PQz,Cpq,Spq,
     $          LegendreP,Cosine,Sine,RToTh,
     $          FactOlm0,FactOlm1,FactOlm2) 
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
                     smn=(-1.0D0)**(mnabs+1)
                  ELSE
                     cmn=1.0D0
                     smn=1.0D0
                  ENDIF
                  IF(n.LT.0)THEN
                     cn=(-1.0D0)**nabs
                     sn=(-1.0D0)**(nabs+1)
                  ELSE
                     cn=1.0D0
                     sn=1.0D0
                  ENDIF
                  Cp(ldx)=Cp(ldx)+cn*Cq(kdx)*cmn*Cpq(lkdx)
     $                           -sn*Sq(kdx)*smn*Spq(lkdx)        
                  Sp(ldx)=Sp(ldx)+cn*Cq(kdx)*smn*Spq(lkdx)
     $                           +sn*Sq(kdx)*cmn*Cpq(lkdx)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      RETURN
      END
C---------------------------------------------------------------
C     Set Up the Factorials for the Multipoles
C---------------------------------------------------------------
      SUBROUTINE MultipoleSetUp(MaxEll,Factorial,
     $                          FactOlm0,FactOlm1,FactOlm2,
     $                          FactMlm0,FactMlm1,FactMlm2)
C
      IMPLICIT INTEGER(i-n)
      IMPLICIT REAL*8 (a-h,o-z)
C
      REAL*8 Factorial(0:2*MaxEll)
      REAL*8 FactOlm0( 0:MaxEll)
      REAL*8 FactOlm1( 0:MaxEll)
      REAL*8 FactOlm2( 0:MaxEll,0:MaxEll)
      REAL*8 FactMlm0( 0:MaxEll)
      REAL*8 FactMlm1( 0:MaxEll)
      REAL*8 FactMlm2( 0:MaxEll,0:MaxEll)
C      
      Factorial(0)=1.0d0
      DO l=1,2*MaxEll
         Factorial(l)=Factorial(l-1)*FLOAT(l)         
      ENDDO
C
      DO l=0,MaxEll
         FactOlm1(l)=FLOAT(2*l-1)
         FactMlm1(l)=FLOAT(2*l-1)
      ENDDO
C
      DO l=2,MaxEll
         DO m=0,l-2
            FactOlm2(m,l)=1.0D0/FLOAT((l+m)*(l-m))
         ENDDO
      ENDDO
C
      DO l=0,MaxEll
         DO m=0,MaxEll
            FactMlm2(m,l)=FLOAT((l+m-1)*(l-m-1))
         ENDDO
      ENDDO
C
      Sgn=1.0D0
      DblFact=1.D0
      TwoTimes=1.D0
      DO m=0,MaxEll
         FactOlm0(m)=Sgn*DblFact/Factorial(2*m)
         FactMlm0(m)=Sgn*DblFact
         DblFact=DblFact*TwoTimes
         TwoTimes=TwoTimes+2.0D0
         Sgn=-Sgn
      ENDDO
C
      RETURN
      END
C---------------------------------------------------------------
C     The Regular Multipoles
C---------------------------------------------------------------
      SUBROUTINE Regular(MaxEll,Ell,EllNdx,x,y,z,C,S,
     $                   LegendreP,Cosine,Sine,RToTh,
     $                   FactOlm0,FactOlm1,FactOlm2)
C
      IMPLICIT INTEGER(i-n)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Ell,EllNdx
C
      REAL*8 C(0:EllNdx)
      REAL*8 S(0:EllNdx)
C
      REAL*8 LegendreP(0:MaxEll,0:MaxEll)
      REAL*8 Cosine(   0:MaxEll)
      REAL*8 Sine(     0:MaxEll)
      REAL*8 RToTh(    0:MaxEll)
      REAL*8 FactOlm0( 0:MaxEll)
      REAL*8 FactOlm1( 0:MaxEll)
      REAL*8 FactOlm2( 0:MaxEll,0:MaxEll)
C
      x2=x*x
      y2=y*y      
      R=DSQRT(x2+y2+z*z)
C
      IF(R.EQ.0.0D0)THEN
         C(0)=1.0D0
         S(0)=0.0D0
         DO lm=1,EllNdx
            C(lm)=0.0D0
            S(lm)=0.0D0
          ENDDO
          RETURN
      ENDIF
C 
      CoTan=z/R
      Cosine(0)=1.0D0
      Sine(  0)=0.0D0
C
      Rxy=DSQRT(x2+y2)
      IF(Rxy.NE.0.0D0)THEN
         OneOvRxy=1.0D0/DSQRT(x2+y2)       
         Cosine(1)=y*OneOvRxy
         Sine(  1)=x*OneOvRxy
      ELSE
         Cosine(1)=0.70710678118654752D0
         Sine(  1)=0.70710678118654752D0         
      ENDIF
C
      TwoC=2.0D0*Cosine(1)
      DO m=2,Ell
         m1=m-1
         m2=m-2
         Cosine(m)=TwoC*Cosine(m1)-Cosine(m2)
         Sine(  m)=TwoC*Sine(  m1)-Sine(  m2)
      ENDDO
C
      Rl=1.0D0
      RS=1.0D0
      Sq=DSQRT(1.0D0-CoTan*CoTan)
      DO m=0,Ell
         LegendreP(m,m)=FactOlm0(m)*RS
         RToTh(m)=Rl
         Rl=Rl*R
         RS=RS*Sq
      ENDDO
C
      DO m=0,Ell-1
         LegendreP(m,m+1)=CoTan*LegendreP(m,m)
      ENDDO
C
      DO l=2,Ell         
         Fact3=CoTan*FLOAT(2*l-1)
         DO m=0,l-2
            LegendreP(m,l)=(LegendreP(m,l-1)*Fact3
     $                    - LegendreP(m,l-2))*FactOlm2(m,l)
         ENDDO   
      ENDDO
C
      DO l=0,ell
         o=RToTh(l)
         ldx=(l+1)*l/2
         DO m=0,l
            lm=ldx+m
            C(lm)=o*LegendreP(m,l)*Cosine(m)
            S(lm)=o*LegendreP(m,l)*Sine(  m)
         ENDDO
      ENDDO
C
      RETURN
      END
C---------------------------------------------------------------
C     The IrRegular Multipoles
C---------------------------------------------------------------
      SUBROUTINE IrRegular(MaxEll,Ell,EllNdx,x,y,z,C,S,
     $                     LegendreP,Cosine,Sine,RToTh,
     $                     FactMlm0,FactMlm1,FactMlm2)
C
      IMPLICIT INTEGER(i-n)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Ell,EllNdx
C
      REAL*8 C(0:EllNdx)
      REAL*8 S(0:EllNdx)
C
      REAL*8 LegendreP(0:MaxEll,0:MaxEll)
      REAL*8 Cosine(   0:MaxEll)
      REAL*8 Sine(     0:MaxEll)
      REAL*8 RToTh(    0:MaxEll)
      REAL*8 FactMlm0( 0:MaxEll)
      REAL*8 FactMlm1( 0:MaxEll)
      REAL*8 FactMlm2( 0:MaxEll,0:MaxEll)
C
      x2=x*x
      y2=y*y      
      R=DSQRT(x2+y2+z*z)
      OneOverR=1.0D0/R
C
      Cosine(0)=1.0D0
      Sine(  0)=0.0D0
C
      Rxy=DSQRT(x2+y2)
      IF(Rxy.NE.0.0D0)THEN
         OneOvRxy=1.0D0/DSQRT(x2+y2)       
         Cosine(1)=y*OneOvRxy
         Sine(  1)=x*OneOvRxy
      ELSE
         Cosine(1)=0.70710678118654752D0
         Sine(  1)=0.70710678118654752D0         
      ENDIF
c
      Rl=OneOverR
      CoTan=z*Rl
C
      TwoC=2.0D0*Cosine(1)
      DO m=2,Ell
         m1=m-1
         m2=m-2
         Cosine(m)=TwoC*Cosine(m1)-Cosine(m2)
         Sine(  m)=TwoC*Sine(  m1)-Sine(  m2)
      ENDDO
C             
      DO m=0,Ell
         RToTh(m)=Rl
         Rl=Rl*OneOverR
      ENDDO
C
      Sq=DSQRT(1.0D0-CoTan*CoTan)
      rs=1.D0
      DO m=0,Ell
         LegendreP(m,m)=FactMlm0(m)*rs
         rs=rs*sq
      ENDDO
C 
      DO m=0,Ell-1
         m1=m+1         
         LegendreP(m,m1)=CoTan*FLOAT(2*m+1)*LegendreP(m,m)
      ENDDO
C
      DO l=2,Ell         
         Fact=CoTan*FLOAT(2*l-1)
         DO m=0,l-2
            LegendreP(m,l)=LegendreP(m,l-1)*Fact
     $                    -LegendreP(m,l-2)*FactMlm2(m,l)
         ENDDO
      ENDDO
C
      DO l=0,Ell
         o=RToTh(l)
         ldx=(l+1)*l/2
         DO m=0,l
            lm=ldx+m
            C(lm)=o*LegendreP(m,l)*Cosine(m)
            S(lm)=o*LegendreP(m,l)*Sine(  m)
         ENDDO
      ENDDO
C    
      RETURN
      END
C
C     Print the Multipoles
C
      SUBROUTINE PrntSH(Ell,EllNdx,C,S)
C
      IMPLICIT INTEGER(i-n)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Ell,EllNdx
C
      REAL*8 C(0:EllNdx)
      REAL*8 S(0:EllNdx)
C     
      WRITE(*,*)'====================================='
      DO l=0,Ell
         DO m=0,l
            lmdx=l*(l+1)/2+m
            WRITE(*,99)l,m,lmdx,C(lmdx),S(lmdx)
 99         FORMAT(3(i3),2(2x,E12.6))
         ENDDO
      ENDDO
      RETURN
      END
C---------------------------------------------------------------
C     Given a Gaussian basis Function , calculate the Multipole 
C     Momments
C---------------------------------------------------------------
      SUBROUTINE HGTFToSP(lt,zz,rc,rs,rhgtf)                          
      IMPLICIT INTEGER(i-n)
c                                                                       
      implicit real*8 (a-h,o-z)                                         
c                                                                       
      parameter (re0=0.d0)                                              
c                                                                       
      REAL*8 rc(0:*),rs(0:*),rhgtf(*)                                
c                                                                       
      if(lt.eq.0)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
      elseif(lt.eq.1)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
      elseif(lt.eq.2)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
        rc(3)=zz*(-5.d-1*rhgtf(5)-5.d-1*rhgtf(8)+rhgtf(10))
        rs(3)=re0
        rc(4)=-5.d-1*zz*rhgtf(9)
        rs(4)=-5.d-1*zz*rhgtf(7)
        rc(5)=zz*(-2.5d-1*rhgtf(5)+2.5d-1*rhgtf(8))
        rs(5)=2.5d-1*zz*rhgtf(6)
      elseif(lt.eq.3)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
        rc(3)=zz*(-5.d-1*rhgtf(5)-5.d-1*rhgtf(8)+rhgtf(10))
        rs(3)=re0
        rc(4)=-5.d-1*zz*rhgtf(9)
        rs(4)=-5.d-1*zz*rhgtf(7)
        rc(5)=zz*(-2.5d-1*rhgtf(5)+2.5d-1*rhgtf(8))
        rs(5)=2.5d-1*zz*rhgtf(6)
        rc(6)=zz*(-5.d-1*rhgtf(13)-5.d-1*rhgtf(18)+rhgtf(20))
        rs(6)=re0
        rc(7)=zz*(1.25d-1*rhgtf(12)+3.75d-1*rhgtf(17)-5.d-1*rhgtf(19)
     >   )
        rs(7)=zz*(3.75d-1*rhgtf(11)+1.25d-1*rhgtf(14)-5.d-1*rhgtf(16)
     >   )
        rc(8)=zz*(-2.5d-1*rhgtf(13)+2.5d-1*rhgtf(18))
        rs(8)=2.5d-1*zz*rhgtf(15)
        rc(9)=zz*(1.25d-1*rhgtf(12)-1.25d-1*rhgtf(17))
        rs(9)=zz*(1.25d-1*rhgtf(11)-1.25d-1*rhgtf(14))
      elseif(lt.eq.4)then 
        rc(0)=zz*rhgtf(1)
        rs(0)=re0
        rc(1)=zz*rhgtf(4)
        rs(1)=re0
        rc(2)=-5.d-1*zz*rhgtf(3)
        rs(2)=-5.d-1*zz*rhgtf(2)
        rc(3)=zz*(-5.d-1*rhgtf(5)-5.d-1*rhgtf(8)+rhgtf(10))
        rs(3)=re0
        rc(4)=-5.d-1*zz*rhgtf(9)
        rs(4)=-5.d-1*zz*rhgtf(7)
        rc(5)=zz*(-2.5d-1*rhgtf(5)+2.5d-1*rhgtf(8))
        rs(5)=2.5d-1*zz*rhgtf(6)
        rc(6)=zz*(-5.d-1*rhgtf(13)-5.d-1*rhgtf(18)+rhgtf(20))
        rs(6)=re0
        rc(7)=zz*(1.25d-1*rhgtf(12)+3.75d-1*rhgtf(17)-5.d-1*rhgtf(19)
     >   )
        rs(7)=zz*(3.75d-1*rhgtf(11)+1.25d-1*rhgtf(14)-5.d-1*rhgtf(16)
     >   )
        rc(8)=zz*(-2.5d-1*rhgtf(13)+2.5d-1*rhgtf(18))
        rs(8)=2.5d-1*zz*rhgtf(15)
        rc(9)=zz*(1.25d-1*rhgtf(12)-1.25d-1*rhgtf(17))
        rs(9)=zz*(1.25d-1*rhgtf(11)-1.25d-1*rhgtf(14))
        rc(10)=zz*(3.75d-1*rhgtf(21)+1.25d-1*rhgtf(24)-5.d-1*rhgtf(26
     >   )+3.75d-1*rhgtf(31)-5.d-1*rhgtf(33)+rhgtf(35))
        rs(10)=re0
        rc(11)=zz*(1.25d-1*rhgtf(25)+3.75d-1*rhgtf(32)-5.d-1*rhgtf(34
     >   ))
        rs(11)=zz*(3.75d-1*rhgtf(23)+1.25d-1*rhgtf(28)-5.d-1*rhgtf(30
     >   ))
        rc(12)=zz*(2.5d-1*rhgtf(21)-2.5d-1*rhgtf(26)-2.5d-1*rhgtf(31)
     >   +2.5d-1*rhgtf(33))
        rs(12)=zz*(-1.25d-1*rhgtf(22)-1.25d-1*rhgtf(27)+2.5d-1*rhgtf(
     >   29))
        rc(13)=zz*(1.25d-1*rhgtf(25)-1.25d-1*rhgtf(32))
        rs(13)=zz*(1.25d-1*rhgtf(23)-1.25d-1*rhgtf(28))
        rc(14)=zz*(6.25d-2*rhgtf(21)-6.25d-2*rhgtf(24)+6.25d-2*rhgtf(
     >   31))
        rs(14)=zz*(-6.25d-2*rhgtf(22)+6.25d-2*rhgtf(27))
      else 
        write(*,*)' lt = ',lt 
        stop ' bad indexing in QCrc ' 
      endif 
c                                                                       
      RETURN
      END                                                            
c                                                 


