      INTEGER MXBUF1,MXBUF2,MXINTS,MXBATCH,MXCNT,MXL,MXL2,MXL4
      INTEGER MXB,MInfo
      INTEGER LenS,len0,len1,len2,len3,len4
      REAL*8 re0,re1,re2,Zero,rpi,PSch,Prev,Prev1,Prev2

      PARAMETER (MXBUF1  =  363000)
      PARAMETER (MXBUF2  =   32000)
      PARAMETER (MXINTS  = 3000000)
      PARAMETER (MXBATCH =    3000)
c
c
c
      PARAMETER (MXCNT=10)

      PARAMETER (MXL=2)
      PARAMETER (MXL2=MXL*2,MXL4=MXL*4)
      PARAMETER (MXB=10)
      PARAMETER (MInfo=3)
      PARAMETER (re0=0.0D0,re1=1.0D0,re2=2.0D0)
      PARAMETER (Zero=0.0D0)
      PARAMETER (rpi=3.1415926535898D0)
      PARAMETER (PSch=4.9738746919472D0)
      PARAMETER (Prev=34.986836655254D0)    ! 2 Pi^(5/2)
      PARAMETER (Prev1=5.914967172796D0)    ! Sqrt[2 Pi^(5/2)]
      PARAMETER (Prev2=Prev1*Prev1)
      PARAMETER (LenS=MXL2*(MXL2+3)/2)
      PARAMETER (len0=(mxl+1)*(mxl+2)*(mxl+3)/6)
      PARAMETER (len1=(mxl2+1)*(mxl2+2)*(mxl2+3)/6)
      PARAMETER (len2=mxl2*(mxl2+1)*(mxl2+2)/6)
      PARAMETER (len3=mxl2*(mxl2+1)/2)
      PARAMETER (len4=len2+len3+mxl2+1)
