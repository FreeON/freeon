      SUBROUTINE HRRBraGrad(W,ACx,ACy,ACz,N,NLOCB2,NLOCK2,IType)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER N,NLOCB2,NLOCK2,IType
      REAL*8 W(N,NLOCB2,NLOCK2)
      REAL*8 ACx,ACy,ACz
c
c IType = IType*100 + JType
c  
      IF (N.EQ.0) RETURN
      IF    (IType.EQ.0101) THEN
         RETURN
      ELSEIF(IType.EQ.0301) THEN
         RETURN
      ELSEIF(IType.EQ.0103) THEN
         DO 0103 J=1,NLOCK2
         DO 0103 I=1,N
            s = W(I,1,J)
            x = W(I,2,J)
            y = W(I,3,J)
            z = W(I,4,J)
            W(I,1,J) = x+ACx*s
            W(I,2,J) = y+ACy*s
            W(I,3,J) = z+ACz*s
 0103    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0102) THEN
         DO 0102 J=1,NLOCK2
         DO 0102 I=1,N
            s1 = W(I,1,J)
            s2 = W(I,2,J)
            x2 = W(I,3,J)
            y2 = W(I,4,J)
            z2 = W(I,5,J)
            W(I,1,J)=s1
            W(I,2,J)=x2+ACx*s2
            W(I,3,J)=y2+ACy*s2
            W(I,4,J)=z2+ACz*s2
 0102    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0501) THEN
         RETURN
      ELSEIF(IType.EQ.0203) THEN
         DO 0203 J=1,NLOCK2
         DO 0203 I=1,N
            s1 =W(I, 1,J)
            x1 =W(I, 2,J)
            y1 =W(I, 3,J)
            z1 =W(I, 4,J)
            x2 =W(I, 5,J)
            y2 =W(I, 6,J)
            z2 =W(I, 7,J)
            xx2=W(I, 8,J)
            xy2=W(I, 9,J)
            yy2=W(I,10,J)
            xz2=W(I,11,J)
            yz2=W(I,12,J)
            zz2=W(I,13,J)
            W(I, 1,J)=x1+ACx*s1
            W(I, 2,J)=y1+ACy*s1
            W(I, 3,J)=z1+ACz*s1
            W(I, 4,J)=xx2+ACx*x2
            W(I, 5,J)=xy2+ACy*x2
            W(I, 6,J)=xz2+ACz*x2
            W(I, 7,J)=xy2+ACx*y2
            W(I, 8,J)=yy2+ACy*y2
            W(I, 9,J)=yz2+ACz*y2
            W(I,10,J)=xz2+ACx*z2
            W(I,11,J)=yz2+ACy*z2
            W(I,12,J)=zz2+ACz*z2
 0203    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0502) THEN
         DO 0502 J=1,NLOCK2
         DO 0502 I=1,N
            x1  =W(I, 1,J)
            y1  =W(I, 2,J)
            z1  =W(I, 3,J)
            xx2 =W(I, 4,J)
            xy2 =W(I, 5,J)
            yy2 =W(I, 6,J)
            xz2 =W(I, 7,J)
            yz2 =W(I, 8,J)
            zz2 =W(I, 9,J)
            x3  =W(I,10,J)
            y3  =W(I,11,J)
            z3  =W(I,12,J)
            xx3 =W(I,13,J)
            xy3 =W(I,14,J)
            yy3 =W(I,15,J)
            xz3 =W(I,16,J)
            yz3 =W(I,17,J)
            zz3 =W(I,18,J)
            xx4 =W(I,19,J)
            xy4 =W(I,20,J)
            yy4 =W(I,21,J)
            xz4 =W(I,22,J)
            yz4 =W(I,23,J)
            zz4 =W(I,24,J)
            xxx4=W(I,25,J)
            xxy4=W(I,26,J)
            xyy4=W(I,27,J)
            yyy4=W(I,28,J)
            xxz4=W(I,29,J)
            xyz4=W(I,30,J)
            yyz4=W(I,31,J)
            xzz4=W(I,32,J)
            yzz4=W(I,33,J)
            zzz4=W(I,34,J)
            W(I, 1,J)=x1
            W(I, 2,J)=xx3+ACx*x3
            W(I, 3,J)=xy3+ACy*x3
            W(I, 4,J)=xz3+ACz*x3
            W(I, 5,J)=y1
            W(I, 6,J)=xy3+ACx*y3
            W(I, 7,J)=yy3+ACy*y3
            W(I, 8,J)=yz3+ACz*y3
            W(I, 9,J)=z1
            W(I,10,J)=xz3+ACx*z3
            W(I,11,J)=yz3+ACy*z3
            W(I,12,J)=zz3+ACz*z3
            W(I,13,J)=xx2
            W(I,14,J)=xxx4+ACx*xx4
            W(I,15,J)=xxy4+ACy*xx4
            W(I,16,J)=xxz4+ACz*xx4
            W(I,17,J)=xy2
            W(I,18,J)=xxy4+ACx*xy4
            W(I,19,J)=xyy4+ACy*xy4
            W(I,20,J)=xyz4+ACz*xy4
            W(I,21,J)=yy2
            W(I,22,J)=xyy4+ACx*yy4
            W(I,23,J)=yyy4+ACy*yy4
            W(I,24,J)=yyz4+ACz*yy4
            W(I,25,J)=xz2
            W(I,26,J)=xxz4+ACx*xz4
            W(I,27,J)=xyz4+ACy*xz4
            W(I,28,J)=xzz4+ACz*xz4
            W(I,29,J)=yz2
            W(I,30,J)=xyz4+ACx*yz4
            W(I,31,J)=yyz4+ACy*yz4
            W(I,32,J)=yzz4+ACz*yz4
            W(I,33,J)=zz2
            W(I,34,J)=xzz4+ACx*zz4
            W(I,35,J)=yzz4+ACy*zz4
            W(I,36,J)=zzz4+ACz*zz4
 0502    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0205) THEN
         DO 0205 J=1,NLOCK2
         DO 0205 I=1,N
            s1  =W(I, 1,J)
            x1  =W(I, 2,J)
            y1  =W(I, 3,J)
            z1  =W(I, 4,J)
            x2  =W(I, 5,J)
            y2  =W(I, 6,J)
            z2  =W(I, 7,J)
            xx2 =W(I, 8,J)
            xy2 =W(I, 9,J)
            yy2 =W(I,10,J)
            xz2 =W(I,11,J)
            yz2 =W(I,12,J)
            zz2 =W(I,13,J)
            s3  =W(I,14,J)
            x3  =W(I,15,J)
            y3  =W(I,16,J)
            z3  =W(I,17,J)
            xx3 =W(I,18,J)
            xy3 =W(I,19,J)
            yy3 =W(I,20,J)
            xz3 =W(I,21,J)
            yz3 =W(I,22,J)
            zz3 =W(I,23,J)
            x4  =W(I,24,J)
            y4  =W(I,25,J)
            z4  =W(I,26,J)
            xx4 =W(I,27,J)
            xy4 =W(I,28,J)
            yy4 =W(I,29,J)
            xz4 =W(I,30,J)
            yz4 =W(I,31,J)
            zz4 =W(I,32,J)
            xxx4=W(I,33,J)
            xxy4=W(I,34,J)
            xyy4=W(I,35,J)
            yyy4=W(I,36,J)
            xxz4=W(I,37,J)
            xyz4=W(I,38,J)
            yyz4=W(I,39,J)
            xzz4=W(I,40,J)
            yzz4=W(I,41,J)
            zzz4=W(I,42,J)
            W(I, 1,J)=x1+ACx*s1
            W(I, 2,J)=y1+ACy*s1
            W(I, 3,J)=z1+ACz*s1
            W(I, 4,J)=xx3+ACx*(2.0D0*x3+ACx*s3)
            W(I, 5,J)=xy3+ACx*(y3+ACy*s3)+ACy*x3
            W(I, 6,J)=yy3+ACy*(2.0D0*y3+ACy*s3)
            W(I, 7,J)=xz3+ACx*(z3+ACz*s3)+ACz*x3
            W(I, 8,J)=yz3+ACz*(y3+ACy*s3)+ACy*z3
            W(I, 9,J)=zz3+ACz*(2.0D0*z3+ACz*s3)
            W(I,10,J)=xx2+ACx*x2
            W(I,11,J)=xy2+ACy*x2
            W(I,12,J)=xz2+ACz*x2
            W(I,13,J)=xxx4+ACx*(2.0D0*xx4+ACx*x4)
            W(I,14,J)=xxy4+ACy*(xx4+ACx*x4)+ACx*xy4
            W(I,15,J)=xyy4+ACy*(2.0D0*xy4+ACy*x4)
            W(I,16,J)=xxz4+ACz*(xx4+ACx*x4)+ACx*xz4
            W(I,17,J)=xyz4+ACz*xy4+ACy*(xz4+ACz*x4)
            W(I,18,J)=xzz4+ACz*(2.0D0*xz4+ACz*x4)
            W(I,19,J)=xy2+ACx*y2
            W(I,20,J)=yy2+ACy*y2
            W(I,21,J)=yz2+ACz*y2
            W(I,22,J)=xxy4+ACx*(2.0D0*xy4+ACx*y4)
            W(I,23,J)=xyy4+ACy*(xy4+ACx*y4)+ACx*yy4
            W(I,24,J)=yyy4+ACy*(2.0D0*yy4+ACy*y4)
            W(I,25,J)=xyz4+ACz*xy4+ACx*(yz4+ACz*y4)
            W(I,26,J)=yyz4+ACz*yy4+ACy*(yz4+ACz*y4)
            W(I,27,J)=yzz4+ACz*(2.0D0*yz4+ACz*y4)
            W(I,28,J)=xz2+ACx*z2
            W(I,29,J)=yz2+ACy*z2
            W(I,30,J)=zz2+ACz*z2
            W(I,31,J)=xxz4+ACx*(2.0D0*xz4+ACx*z4)
            W(I,32,J)=xyz4+ACy*xz4+ACx*(yz4+ACy*z4)
            W(I,33,J)=yyz4+ACy*(2.0D0*yz4+ACy*z4)
            W(I,34,J)=xzz4+ACz*xz4+ACx*(zz4+ACz*z4)
            W(I,35,J)=yzz4+ACz*yz4+ACy*(zz4+ACz*z4)
            W(I,36,J)=zzz4+ACz*(2.0D0*zz4+ACz*z4)
 0205    CONTINUE
         RETURN
      ELSEIF(IType.EQ.0201) THEN
         RETURN
      ELSEIF(IType.EQ.0202) THEN
         DO 0202 J=1,NLOCK2
         DO 0202 I=1,N
            s1  =W(I, 1,J)
            s3  =W(I, 2,J)
            x3  =W(I, 3,J)
            y3  =W(I, 4,J)
            z3  =W(I, 5,J)
            x2  =W(I, 6,J)
            y2  =W(I, 7,J)
            z2  =W(I, 8,J)
            x4  =W(I, 9,J)
            y4  =W(I,10,J)
            z4  =W(I,11,J)
            xx4 =W(I,12,J)
            xy4 =W(I,13,J)
            yy4 =W(I,14,J)
            xz4 =W(I,15,J)
            yz4 =W(I,16,J)
            zz4 =W(I,17,J)
            W(I, 1,J)=s1
            W(I, 2,J)=x3+ACx*s3
            W(I, 3,J)=y3+ACy*s3
            W(I, 4,J)=z3+ACz*s3
            W(I, 5,J)=x2
            W(I, 6,J)=xx4+ACx*x4
            W(I, 7,J)=xy4+ACy*x4
            W(I, 8,J)=xz4+ACz*x4
            W(I, 9,J)=y2
            W(I,10,J)=xy4+ACx*y4
            W(I,11,J)=yy4+ACy*y4
            W(I,12,J)=yz4+ACz*y4
            W(I,13,J)=z2
            W(I,14,J)=xz4+ACx*z4
            W(I,15,J)=yz4+ACy*z4
            W(I,16,J)=zz4+ACz*z4
 0202    CONTINUE
      ELSE
         WRITE(*,*) "IType = ",IType
         STOP 'Illegal IType in HRRBra'
      ENDIF
      RETURN
      END
