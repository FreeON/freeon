      SUBROUTINE CalcK(NAtoms,NCtrt,NPrim2,KType,LType,
     >                 VLoc,SLoc,LngVRR,LngLoc,
     >                 BTable,BfnBuf,SubInd,NBasF,
     >                 PrmL,DisL,VecL,CosL,
     >                 SinL,IsaL,BfnL,SrtL,
     >                 iSA,iSrt,CDrv,IndA,Bfn,Indx,
     >                 DA,sK,iK,jK,pK,
     >                 PrmBuf,DisBuf,CBuf,SBuf,VecBuf,
     >                 WR,WZ,W1,W2,RAC,RBD,CBra,CKet,xNERIs)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INCLUDE "pmonx.inc"
      INTEGER PrmL,DisL,VecL,CosL
      INTEGER SinL,IsaL,BfnL,SrtL
      REAL*8  DA(*),sK(*)
      REAL*8  PrmBuf(PrmL),DisBuf(DisL),CBuf(CosL),SBuf(SinL)
      REAL*8  VecBuf(3,VecL),WR(12*MXBUF1),WZ(5*MXBUF1)
      REAL*8  W1(MXINTS),W2(MXINTS),RAC(3,MXBATCH),RBD(3,MXBATCH)
      REAL*8  CBra(4*MXBUF2),CKet(4*MXBUF2)
      INTEGER iK(*),jK(*),pK(*)
      INTEGER VLoc(LngVRR),SLoc(3,LngLoc/3)
      INTEGER BTable(NPrim2,NPrim2),SubInd(3,NBasF)
      INTEGER BfnBuf(BfnL),iSA(5,IsaL),iSrt(NCtrt,NAtoms),CDrv(60000)
      INTEGER IndA(0:500),Bfn(2,MXBATCH),Indx(2,MXBATCH)
      REAL*8  Cx,Cy,Cz,Dx,Dy,Dz
      REAL*8  Dcd,Tsh
      INTEGER AtC,AtD,CFC,CFD
      INTEGER iTB,LenBra,iPrmB,iDisB,iMomB,IType,KBra
      INTEGER iTK,LenKet,iPrmK,iDisK,iMomK,JType,KKet
      INTEGER KType,LType,KThresh
      INTEGER KC,KD,NKCase,NLCase,Ktot,Ltot,ISw
      INTEGER LAmin,LBmin,LAmax,LBmax
      INTEGER LCmin,LDmin,LCmax,LDmax
      INTEGER BraCodO,KetCodO,BraCodS,KetCodS
      INTEGER BraLmax,KetLmax
      INTEGER LenBi,LenKi,LenBf,LenKf,NumB
      INTEGER LenAA,LenBB,LenCC
      INTEGER I0,I1,I2,iCP,iCL

      INTEGER Leng(10),Lmax(10),Lmin(10)
      INTEGER IDmn(0:4),ILen(0:4),iT(6,6)

      INCLUDE "bkonx.inc"
      DATA Leng / 1, 4, 3,10, 9, 6,20,19,16,10 /
      DATA Lmax / 0, 1, 1, 2, 2, 2, 3, 3, 3, 3 /
      DATA Lmin / 0, 0, 1, 0, 1, 2, 0, 1, 2, 3 /
      DATA IDmn / 1, 4,17,24,36 /
      DATA ILen / 1, 4,10,20,35 /
      DATA iT   / 1, 0, 0, 0, 0, 0,
     >            2, 3, 0, 0, 0, 0,
     >            4, 5, 6, 0, 0, 0,
     >            0, 0, 0, 0, 0, 0,
     >            0, 0, 0, 0, 0, 0,
     >            7, 8, 9, 0, 0, 10 /


      KThresh=0
      iSB=iSrt(CFC,AtC)
      iSK=iSrt(CFD,AtD)
 1000 CONTINUE

      iTB    = iSA(1,iSB)
      LenBra = iSA(2,iSB)
      iPrmB  = iSA(3,iSB)
      iDisB  = iSA(4,iSB)
      iMomB  = iSA(5,iSB)
      IType  = iTB/100
      KBra   = iTB-IType*100

      iTK    = iSA(1,iSK)
      LenKet = iSA(2,iSK)
      iPrmK  = iSA(3,iSK)
      iDisK  = iSA(4,iSK)
      iMomK  = iSA(5,iSK)
      JType  = iTK/100
      KKet   = iTK-JType*100
      Ktot   = KBra*KKet

c
c Get the distribution info, switch centers if necessary.
c
      LAmax=Lmax(IType)
      LAmin=Lmin(IType)
      LCmax=Lmax(KType)
      LCmin=Lmin(KType)
      LBmax=Lmax(JType)
      LBmin=Lmin(JType)
      LDmax=Lmax(LType)
      LDmin=Lmin(LType)
      BraCodO=LAmax*1000+LAmin*100+LCmax*10+LCmin
      KetCodO=LBmax*1000+LBmin*100+LDmax*10+LDmin

      IF(IType.GE.KType) THEN
         LAmax=Lmax(IType)
         LAmin=Lmin(IType)
         LCmax=Lmax(KType)
         LCmin=Lmin(Ktype)
         SgnB=1.0D0
         SwhB=1.0D0
         L1=Leng(IType)
         L2=Leng(Ktype)
      ELSE
         LAmax=Lmax(KType)
         LAmin=Lmin(KType)
         LCmax=Lmax(IType)
         LCmin=Lmin(Itype)
         SgnB=-1.0D0
         SwhB= 0.0D0
         L1=Leng(KType)
         L2=Leng(Itype)
      ENDIF
      IF(JType.GE.LType) THEN
         LBmax=Lmax(JType)
         LBmin=Lmin(JType)
         LDmax=Lmax(LType)
         LDmin=Lmin(LType)
         SgnK=1.0D0
         SwhK=1.0D0
         L3=Leng(JType)
         L4=Leng(LType)
      ELSE
         LBmax=Lmax(LType)
         LBmin=Lmin(LType)
         LDmax=Lmax(JType)
         LDmin=Lmin(JType)
         SgnK=-1.0D0
         SwhK= 0.0D0
         L3=Leng(LType)
         L4=Leng(JType)
      ENDIF
c
c Get the length (iCL) and pointer (iCP) of the appropriate contraction
c driver that we may need for this batch (needed only if the general
c integral code is called).
c
      I0=iT(MIN(IType,KType),MAX(IType,KType))
      I1=iT(MIN(JType,LType),MAX(JType,LType))
      I2=I0+(I1-1)*10
      iCP=CDrv(I2)
      iCL=CDrv(iCP)

      LenAA=Leng(IType)
      LenBB=Leng(JType)
      LenCC=Leng(KType)

      ISw=10*INT(SwhB)+INT(SwhK)

      BraCodS=LAmax*1000+LAmin*100+LCmax*10+LCmin
      KetCodS=LBmax*1000+LBmin*100+LDmax*10+LDmin
      BraLmax=LAmax*10+LCmax
      KetLmax=LBmax*10+LDmax

      LBra=LAmax+LCmax
      LKet=LBmax+LDmax
      Ltot=LBra+LKet
      LenBi=IDmn(LBra)
      LenKi=IDmn(LKet)
      LenBf=Leng(IType)*Leng(KType)
      LenKf=Leng(JType)*Leng(LType)
      NumB=ILen(LBra)
c
c Find skipout spots. All the integral prescreening based on contracted
c basis functions is in here. Basic idea is to use a combination of a
c binary and a linear search routine to avoid testing all possible
c significant integrals.  Additional primitive screening is in the
c routines called in IntegralSelect.inc
c
      CALL Skipout(LenBra,LenKet,DisBuf(iDisB),DisBuf(iDisK),
     >             IndA,TxDcd)

      IF(IndA(0).EQ.0) GOTO 6000
c
c What type of integral code are we going to use?
c Figure the space requirements for these types of integrals:
c
      IF(Ltot.LE.4) THEN
         I1   = L1*L2*L3*L4
         MX1  = MXINTS/I1
         MXBT = MIN(MXBATCH,MX1)
      ELSE
         I0   = LenBi*LenKi*KBra*KKet*2
         I1   = L1*L2*L3*L4
         I2   = MAX(I0,I1)
         MX1  = MXINTS/I2
         MX2  = BTable(KBra,KKet)
         MXBT = MIN(MX1,MX2)
      ENDIF
c
c Now figure out how many integrals we can do given the buffer
c size we have available. This (4000) is the top of the loop over
c groups of integrals within a single batch. I need a better way
c to split up the batches that takes into account that some of the
c integrals are going to get thrown out of the batch because Ixa
c is less that Ixb (symmetry of the exchange matrix).
c
      iBgn=1
      iEnd=0
      iCnt=0
      IF(IndA(1).GT.MXBT) THEN
         WRITE(*,*) "LenBi=",LenBi
         WRITE(*,*) "LenKi=",LenKi
         WRITE(*,*) "I0=",I0
         WRITE(*,*) "I1=",I1
         WRITE(*,*) "MX1=",MX1
         WRITE(*,*) "MX2=",MX2
         WRITE(*,*) "KBra=",KBra
         WRITE(*,*) "KKet=",KKet
         WRITE(*,*) "MXBT=",MXBT
         WRITE(*,*) "IndA=",IndA(1)
         WRITE(*,*) "Out of memory in ONX"
         STOP
      ENDIF
 4000 CONTINUE
      I0=iCnt+IndA(iEnd+1)
      IF(I0.LE.MXBT.AND.iEnd.LT.IndA(0)) THEN
         iCnt=I0
         iEnd=iEnd+1
         GOTO 4000
      ENDIF
c
c Check for multipoles
c
c      INCLUDE "MultipoleSelect.inc"
c      IF(Mults) THEN
c
c         CALL Digest(LTot,N,LenAA,LenBB,LenCC,KThresh,
c     >               L1,L2,L3,L4,ISw,pK,jK,iK,
c     >               SubInd,Bfn,W2,sK,W1,DA)
c
c         xNMults = xNMults + FLOAT(L1*L2*L3*L4*N)
c
c      ENDIF

c
c Find the true length of the current batch and store some info.
c
      IF(BraCodS.GT.0) THEN
         DO I=iBgn,iEnd
            iBra  = I+iDisB-1
            I0=I-iBgn+1
            RAC(1,I0)=SgnB*VecBuf(1,iBra)
            RAC(2,I0)=SgnB*VecBuf(2,iBra)
            RAC(3,I0)=SgnB*VecBuf(3,iBra)
        ENDDO
      ENDIF

      IF(KetCodS.GT.0) THEN
         DO J=1,ABS(IndA(iBgn))
            iKet=J+iDisK-1
            RBD(1,J)=SgnK*VecBuf(1,iKet)
            RBD(2,J)=SgnK*VecBuf(2,iKet)
            RBD(3,J)=SgnK*VecBuf(3,iKet)
         ENDDO
      ENDIF

      N=0
      DO 1231 I=iBgn,iEnd
         IndexA = BfnBuf(I+iDisB-1)
         DO 1232 J=1,IndA(I)
            IndexB = BfnBuf(J+iDisK-1)
            IF(IndexA.GE.IndexB) THEN
               N=N+1
               Indx(1,N)=I
               Indx(2,N)=J
               Bfn(1,N)=IndexA
               Bfn(2,N)=IndexB
            ENDIF
 1232    CONTINUE
 1231 CONTINUE

            IF(N.EQ.0) GOTO 6000

            xNERIs = xNERIs + FLOAT(L1*L2*L3*L4*N)

            IF(LTot.LE.4) THEN
              INCLUDE "IntegralSelect.inc"
              NumB=IDmn(LBra)
              NumK=IDmn(LKet)
              IF(KetCodS.NE.0000) THEN
                CALL HRRKet(W2,RBD,N,Indx,LenBi,NumB,KetCodS)
              ENDIF
              IF(BraCodS.NE.0000) THEN
                CALL HRRBra(W2,W1,RAC,N,Indx,LenBi,LenBf,
     >                      LenKf,BraCodS)
                CALL Digest(LTot,N,LenAA,LenBB,LenCC,KThresh,
     >                      L1,L2,L3,L4,ISw,pK,jK,iK,
     >                      SubInd,Bfn,W1,sK,W2,DA)
              ELSE
                CALL Digest(LTot,N,LenAA,LenBB,LenCC,KThresh,
     >                      L1,L2,L3,L4,ISw,pK,jK,iK,
     >                      SubInd,Bfn,W2,SK,W1,DA)
              ENDIF
            ELSE
c
c Form the initial R[0,m] integrals and use the general code to do the
c vertical recurrence on the bra and the ket side (vector code).
c
               NBK=N*KBra*KKet
               IF(Ltot.EQ.0) THEN
                  CALL TGen0(Ltot,KBra,KKet,N,Indx,
     >                       PrmBuf(iPrmB),PrmBuf(iPrmK),W1)
               ELSEIF(Ltot.LE.2) THEN
                  CALL TGen1t2(Ltot,KBra,KKet,N,Indx,
     >                         SwhB,SwhK,
     >                         PrmBuf(iPrmB),PrmBuf(iPrmK),
     >                         CBra,CKet,WR,WZ,RAC,RBD,W2)
                  CALL VRRs(NBK,LBra,LKet,W1,W2,WR,WZ)
               ELSEIF(Ltot.LE.4) THEN
                  CALL TGen3t4(Ltot,KBra,KKet,N,Indx,
     >                         SwhB,SwhK,
     >                         PrmBuf(iPrmB),PrmBuf(iPrmK),
     >                         CBra,CKet,WR,WZ,RAC,RBD,W2)
                  CALL VRRs(NBK,LBra,LKet,W1,W2,WR,WZ)
               ELSEIF(Ltot.LE.6) THEN
                  CALL TGen5t6(Ltot,KBra,KKet,N,Indx,
     >                         SwhB,SwhK,
     >                         PrmBuf(iPrmB),PrmBuf(iPrmK),
     >                         CBra,CKet,WR,WZ,RAC,RBD,W2)
                  id = LBra+1+LKet*5
                  is = SLoc(1,id)
                  nr = SLoc(2,id)
                  ns = SLoc(3,id)
                  call VRRl(NBK,NR,NS,VLoc(is),VLoc(is+nr),W1,W2,WR,WZ)

               ELSE
                  CALL TGen7t8(Ltot,KBra,KKet,N,Indx,
     >                       SwhB,SwhK,
     >                       PrmBuf(iPrmB),PrmBuf(iPrmK),
     >                       CBra,CKet,WR,WZ,RAC,RBD,W2)
                  id = LBra+1+LKet*5
                  is = SLoc(1,id)
                  nr = SLoc(2,id)
                  ns = SLoc(3,id)
                  call VRRl(NBK,NR,NS,VLoc(is),VLoc(is+nr),W1,W2,WR,WZ)
               ENDIF
c
c Contract the integrals.
c
               CALL Contract(N,KBra,KKet,iCL,CDrv(iCP+1),
     >                       CBra,CKet,W2,W1)
c
c Do the horizontal recurrence on the ket and the bra. Then digest
c the integrals into the exchange matrix.
c
               NumB=IDmn(LBra)
               NumK=IDmn(LKet)
               IF(KetCodS.NE.0000) THEN
                  CALL HRRKetOld(W2,RBD,N,Indx,LenBi,NumB,KetCodS)
               ENDIF
               IF(BraCodS.NE.0000) THEN
                  CALL HRRBraOld(W2,W1,RAC,N,Indx,LenBi,LenBf,
     >                           LenKf,BraCodS)
                  CALL Digest(LTot,N,LenAA,LenBB,LenCC,KThresh,
     >                        L1,L2,L3,L4,ISw,pK,jK,iK,
     >                        SubInd,Bfn,W1,sK,W2,DA)
               ELSE
                  CALL Digest(LTot,N,LenAA,LenBB,LenCC,KThresh,
     >                        L1,L2,L3,L4,ISw,pK,jK,iK,
     >                        SubInd,Bfn,W2,SK,W1,DA)
               ENDIF
            ENDIF
c
c Check to see if there are any more integrals to do in this
c batch before we grab a new set.
c
 6000       CONTINUE
            IF(iEnd.LT.IndA(0)) THEN
               iBgn=iEnd+1
               iCnt=0
               GOTO 4000
            ENDIF
c
c Get a new set of distribution types:
c
            iSB = iSB + 1
            IF(iSA(1,iSB).NE.-1) THEN
               GOTO 1000
            ELSE
               iSB = ISrt(CFC,AtC)
               iSK = iSK+1
               IF(iSA(1,iSK).NE.-1) GOTO 1000
            ENDIF
      RETURN
      END
