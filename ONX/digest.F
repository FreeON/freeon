      SUBROUTINE Digest(Ltot,N,LenA,LenB,LenC,KThresh,
     &                  L1,L2,L3,L4,ISw,pK,jK,iK,
     &                  SubInd,Bfn,W,sK,KAT,DA)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 W(N,L2,L1,L4,*)
      REAL*8 sK(*)
      REAL*8 DA(LenC,*)
      REAL*8 KAT(N,LenA,*)
      INTEGER SubInd(3,*)
      INTEGER pK(*),jK(*),iK(*)
      INTEGER Bfn(2,*)
      INTEGER AtA,AtB
      INTEGER rs,cs
      common /cntr/ intcnt

      CALL Rloader(KAT,N*LenA*LenB,0.0D0)


      IF(ISw.EQ.11) THEN
         DO IB=1,L3
         DO ID=1,L4
         DO IA=1,L1
         DO IC=1,L2
            DAcd=DA(IC,ID)
            DO I=1,N
              KAT(I,IA,IB)=KAT(I,IA,IB)-W(I,IC,IA,ID,IB)*DAcd
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(ISw.EQ.01) THEN
         DO IB=1,L3
         DO ID=1,L4
         DO IC=1,L1
            DAcd=DA(IC,ID)
         DO IA=1,L2
            DO I=1,N
               KAT(I,IA,IB)=KAT(I,IA,IB)-W(I,IA,IC,ID,IB)*DAcd
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(ISw.EQ.10) THEN
         DO ID=1,L3
         DO IB=1,L4
         DO IA=1,L1
         DO IC=1,L2
            DAcd=DA(IC,ID)
            DO I=1,N
              KAT(I,IA,IB)=KAT(I,IA,IB)-W(I,IC,IA,IB,ID)*DAcd
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ELSEIF(ISw.EQ.00) THEN
         NN=N*L2
         DO ID=1,L3
         DO IB=1,L4
         DO IC=1,L1
            DAcd=DA(IC,ID)
         DO IA=1,L2
            DO I=1,N
              KAT(I,IA,IB)=KAT(I,IA,IB)-W(I,IA,IC,IB,ID)*DAcd
            ENDDO
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ELSE
         WRITE(*,*) "ISw = ",ISw
         CALL MondoHalt('Illegal ISw in Digest')
      ENDIF

      DO I=1,N
         IndexA = Bfn(1,I)
         IndexB = Bfn(2,I)
         AtA    = SubInd(1,IndexA)
         NBFA   = SubInd(2,IndexA)
         rs     = SubInd(3,IndexA)
         AtB    = SubInd(1,IndexB)
         NBFB   = SubInd(2,IndexB)
         cs     = SubInd(3,IndexB)
         CALL GetAdrB(AtA,AtB,Ind,KThresh,jK,iK) 
         IF(Ind.GT.0) THEN
            CALL PutSubBlk(I,N,NBFA,NBFB,LenA,LenB,rs,cs,
     >                     sK(pK(Ind)),KAT)
         ENDIF
      ENDDO
      RETURN
      END
