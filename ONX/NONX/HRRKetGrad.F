      SUBROUTINE HRRKetGrad(W,BD,N,SLDis,LenB,IType)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER N,LenB,IType,NN
      INTEGER I,Ioff
      INTEGER SLDis(N)
      REAL*8 BD(*)
      REAL*8 W(N*LenB,*)
c
c IType = KType*100 + LType
c  
      IF (N.EQ.0) RETURN
      NN=N*LenB
      IF(    IType.EQ.0101) THEN
         RETURN
      ELSEIF(IType.EQ.0301) THEN
         RETURN
      ELSEIF(IType.EQ.0103) THEN
         DO 0103 I=1,NN
            Ioff = SLDis(MOD(I-1,N)+1)
            BDx  = BD(Ioff  )
            BDy  = BD(Ioff+1)
            BDz  = BD(Ioff+2)
            s = W(I,1)
            x = W(I,2)
            y = W(I,3)
            z = W(I,4)
            W(I,1) = x+BDx*s
            W(I,2) = y+BDy*s
            W(I,3) = z+BDz*s
 0103    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0102) THEN
         DO 0102 I=1,NN
            Ioff = SLDis(MOD(I-1,N)+1)
            BDx  = BD(Ioff  )
            BDy  = BD(Ioff+1)
            BDz  = BD(Ioff+2)
            s1 = W(I,1)
            s2 = W(I,2)
            x2 = W(I,3)
            y2 = W(I,4)
            z2 = W(I,5)
            W(I,1)=s1
            W(I,2)=x2+BDx*s2
            W(I,3)=y2+BDy*s2
            W(I,4)=z2+BDz*s2
 0102    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0501) THEN
         RETURN
      ELSEIF(IType.EQ.0203) THEN
         DO 0203 I=1,NN
            Ioff = SLDis(MOD(I-1,N)+1)
            BDx  = BD(Ioff  )
            BDy  = BD(Ioff+1)
            BDz  = BD(Ioff+2)
            s1 =W(I, 1)
            x1 =W(I, 2)
            y1 =W(I, 3)
            z1 =W(I, 4)
            x2 =W(I, 5)
            y2 =W(I, 6)
            z2 =W(I, 7)
            xx2=W(I, 8)
            xy2=W(I, 9)
            yy2=W(I,10)
            xz2=W(I,11)
            yz2=W(I,12)
            zz2=W(I,13)
            W(I, 1)=x1+BDx*s1
            W(I, 2)=y1+BDy*s1
            W(I, 3)=z1+BDz*s1
            W(I, 4)=xx2+BDx*x2
            W(I, 5)=xy2+BDy*x2
            W(I, 6)=xz2+BDz*x2
            W(I, 7)=xy2+BDx*y2
            W(I, 8)=yy2+BDy*y2
            W(I, 9)=yz2+BDz*y2
            W(I,10)=xz2+BDx*z2
            W(I,11)=yz2+BDy*z2
            W(I,12)=zz2+BDz*z2
 0203    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0502) THEN
         DO 0502 I=1,NN
            Ioff = SLDis(MOD(I-1,N)+1)
            BDx  = BD(Ioff  )
            BDy  = BD(Ioff+1)
            BDz  = BD(Ioff+2)
            x1  =W(I, 1)
            y1  =W(I, 2)
            z1  =W(I, 3)
            xx2 =W(I, 4)
            xy2 =W(I, 5)
            yy2 =W(I, 6)
            xz2 =W(I, 7)
            yz2 =W(I, 8)
            zz2 =W(I, 9)
            x3  =W(I,10)
            y3  =W(I,11)
            z3  =W(I,12)
            xx3 =W(I,13)
            xy3 =W(I,14)
            yy3 =W(I,15)
            xz3 =W(I,16)
            yz3 =W(I,17)
            zz3 =W(I,18)
            xx4 =W(I,19)
            xy4 =W(I,20)
            yy4 =W(I,21)
            xz4 =W(I,22)
            yz4 =W(I,23)
            zz4 =W(I,24)
            xxx4=W(I,25)
            xxy4=W(I,26)
            xyy4=W(I,27)
            yyy4=W(I,28)
            xxz4=W(I,29)
            xyz4=W(I,30)
            yyz4=W(I,31)
            xzz4=W(I,32)
            yzz4=W(I,33)
            zzz4=W(I,34)
            W(I, 1)=x1
            W(I, 2)=xx3+BDx*x3
            W(I, 3)=xy3+BDy*x3
            W(I, 4)=xz3+BDz*x3
            W(I, 5)=y1
            W(I, 6)=xy3+BDx*y3
            W(I, 7)=yy3+BDy*y3
            W(I, 8)=yz3+BDz*y3
            W(I, 9)=z1
            W(I,10)=xz3+BDx*z3
            W(I,11)=yz3+BDy*z3
            W(I,12)=zz3+BDz*z3
            W(I,13)=xx2
            W(I,14)=xxx4+BDx*xx4
            W(I,15)=xxy4+BDy*xx4
            W(I,16)=xxz4+BDz*xx4
            W(I,17)=xy2
            W(I,18)=xxy4+BDx*xy4
            W(I,19)=xyy4+BDy*xy4
            W(I,20)=xyz4+BDz*xy4
            W(I,21)=yy2
            W(I,22)=xyy4+BDx*yy4
            W(I,23)=yyy4+BDy*yy4
            W(I,24)=yyz4+BDz*yy4
            W(I,25)=xz2
            W(I,26)=xxz4+BDx*xz4
            W(I,27)=xyz4+BDy*xz4
            W(I,28)=xzz4+BDz*xz4
            W(I,29)=yz2
            W(I,30)=xyz4+BDx*yz4
            W(I,31)=yyz4+BDy*yz4
            W(I,32)=yzz4+BDz*yz4
            W(I,33)=zz2
            W(I,34)=xzz4+BDx*zz4
            W(I,35)=yzz4+BDy*zz4
            W(I,36)=zzz4+BDz*zz4
 0502    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0205) THEN
         DO 0205 I=1,NN
            Ioff = SLDis(MOD(I-1,N)+1)
            BDx  = BD(Ioff  )
            BDy  = BD(Ioff+2)
            BDz  = BD(Ioff+2)
            s1  =W(I, 1)
            x1  =W(I, 2)
            y1  =W(I, 3)
            z1  =W(I, 4)
            x2  =W(I, 5)
            y2  =W(I, 6)
            z2  =W(I, 7)
            xx2 =W(I, 8)
            xy2 =W(I, 9)
            yy2 =W(I,10)
            xz2 =W(I,11)
            yz2 =W(I,12)
            zz2 =W(I,13)
            s3  =W(I,14)
            x3  =W(I,15)
            y3  =W(I,16)
            z3  =W(I,17)
            xx3 =W(I,18)
            xy3 =W(I,19)
            yy3 =W(I,20)
            xz3 =W(I,21)
            yz3 =W(I,22)
            zz3 =W(I,23)
            x4  =W(I,24)
            y4  =W(I,25)
            z4  =W(I,26)
            xx4 =W(I,27)
            xy4 =W(I,28)
            yy4 =W(I,29)
            xz4 =W(I,30)
            yz4 =W(I,31)
            zz4 =W(I,32)
            xxx4=W(I,33)
            xxy4=W(I,34)
            xyy4=W(I,35)
            yyy4=W(I,36)
            xxz4=W(I,37)
            xyz4=W(I,38)
            yyz4=W(I,39)
            xzz4=W(I,40)
            yzz4=W(I,41)
            zzz4=W(I,42)
            W(I, 1)=x1+BDx*s1
            W(I, 2)=y1+BDy*s1
            W(I, 3)=z1+BDz*s1
            W(I, 4)=xx3+BDx*(2.0D0*x3+BDx*s3)
            W(I, 5)=xy3+BDx*(y3+BDy*s3)+BDy*x3
            W(I, 6)=yy3+BDy*(2.0D0*y3+BDy*s3)
            W(I, 7)=xz3+BDx*(z3+BDz*s3)+BDz*x3
            W(I, 8)=yz3+BDz*(y3+BDy*s3)+BDy*z3
            W(I, 9)=zz3+BDz*(2.0D0*z3+BDz*s3)
            W(I,10)=xx2+BDx*x2
            W(I,11)=xy2+BDy*x2
            W(I,12)=xz2+BDz*x2
            W(I,13)=xxx4+BDx*(2.0D0*xx4+BDx*x4)
            W(I,14)=xxy4+BDy*(xx4+BDx*x4)+BDx*xy4
            W(I,15)=xyy4+BDy*(2.0D0*xy4+BDy*x4)
            W(I,16)=xxz4+BDz*(xx4+BDx*x4)+BDx*xz4
            W(I,17)=xyz4+BDz*xy4+BDy*(xz4+BDz*x4)
            W(I,18)=xzz4+BDz*(2.0D0*xz4+BDz*x4)
            W(I,19)=xy2+BDx*y2
            W(I,20)=yy2+BDy*y2
            W(I,21)=yz2+BDz*y2
            W(I,22)=xxy4+BDx*(2.0D0*xy4+BDx*y4)
            W(I,23)=xyy4+BDy*(xy4+BDx*y4)+BDx*yy4
            W(I,24)=yyy4+BDy*(2.0D0*yy4+BDy*y4)
            W(I,25)=xyz4+BDz*xy4+BDx*(yz4+BDz*y4)
            W(I,26)=yyz4+BDz*yy4+BDy*(yz4+BDz*y4)
            W(I,27)=yzz4+BDz*(2.0D0*yz4+BDz*y4)
            W(I,28)=xz2+BDx*z2
            W(I,29)=yz2+BDy*z2
            W(I,30)=zz2+BDz*z2
            W(I,31)=xxz4+BDx*(2.0D0*xz4+BDx*z4)
            W(I,32)=xyz4+BDy*xz4+BDx*(yz4+BDy*z4)
            W(I,33)=yyz4+BDy*(2.0D0*yz4+BDy*z4)
            W(I,34)=xzz4+BDz*xz4+BDx*(zz4+BDz*z4)
            W(I,35)=yzz4+BDz*yz4+BDy*(zz4+BDz*z4)
            W(I,36)=zzz4+BDz*(2.0D0*zz4+BDz*z4)
 0205    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0201) THEN
         RETURN
      ELSEIF(IType.EQ.0202) THEN
         DO 0202 I=1,NN
            Ioff = SLDis(MOD(I-1,N)+1)
            BDx  = BD(Ioff  )
            BDy  = BD(Ioff+1)
            BDz  = BD(Ioff+2)
            s1  =W(I, 1)
            s3  =W(I, 2)
            x3  =W(I, 3)
            y3  =W(I, 4)
            z3  =W(I, 5)
            x2  =W(I, 6)
            y2  =W(I, 7)
            z2  =W(I, 8)
            x4  =W(I, 9)
            y4  =W(I,10)
            z4  =W(I,11)
            xx4 =W(I,12)
            xy4 =W(I,13)
            yy4 =W(I,14)
            xz4 =W(I,15)
            yz4 =W(I,16)
            zz4 =W(I,17)
            W(I, 1)=s1
            W(I, 2)=x3+BDx*s3
            W(I, 3)=y3+BDy*s3
            W(I, 4)=z3+BDz*s3
            W(I, 5)=x2
            W(I, 6)=xx4+BDx*x4
            W(I, 7)=xy4+BDy*x4
            W(I, 8)=xz4+BDz*x4
            W(I, 9)=y2
            W(I,10)=xy4+BDx*y4
            W(I,11)=yy4+BDy*y4
            W(I,12)=yz4+BDz*y4
            W(I,13)=z2
            W(I,14)=xz4+BDx*z4
            W(I,15)=yz4+BDy*z4
            W(I,16)=zz4+BDz*z4
 0202    CONTINUE
      ELSE
         WRITE(*,*) "IType = ",IType
         STOP 'Illegal IType in HRRKetGrad'
      ENDIF
      RETURN
      END
