      PROGRAM CGMaker
c
c This program is used to generate the G1 and G4 type of driver files
c for the ONX gradient code. The G2 and G3 files are contained in the 
c files fort72 and fort73. Some of the HRR steps that are generated in 
c these files can be skipped because they dont really do anything (for 
c example any HRR calls with 101 as the integral type).
c
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER LType(8)
      INTEGER IUN(3)
c A -------------------------------------------------- 
      LType(1) = 0
      LType(2) = 1
c C --------------------------------------------------
      LType(3) = 0
      LType(4) = 1
c B --------------------------------------------------
      LType(5) = 0
      LType(6) = 1
c D --------------------------------------------------
      LType(7) = 0
      LType(8) = 1

c      CALL CIndx(LType,IUN,0)      ! generates the G1 type files
      CALL PIndx(LType)            ! generates the G4 type files

      STOP
      END


      SUBROUTINE PIndx(LType)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER LType(8),Facm
      INTEGER ITree(2,3,10,10)
      INTEGER IUN(3)

      CALL CIndx(LType,IUN,1)

      CALL GetIndexT(ITree)

      LBra = LType(2)+LType(4)
      LKet = LType(6)+LType(8)
      LAmin = LType(1)
      LAmax = LType(2)
      LCmin = LType(3)
      LCmax = LType(4)
      LBmin = LType(5)
      LBmax = LType(6)
      LDmin = LType(7)
      LDmax = LType(8)
      ICA   = LAmax*(LAmax+1)/2+LAmin+1
      ICC   = LCmax*(LCmax+1)/2+LCmin+1
      ICB   = LBmax*(LBmax+1)/2+LBmin+1
      ICD   = LDmax*(LDmax+1)/2+LDmin+1
      NA    = NFinal(ICA)
      NC    = NFinal(ICC)
      NB    = NFinal(ICB)
      ND    = NFinal(ICD)

      Ioff  = 0
c---------------------------------------------------------------------- 
c Center A:
c---------------------------------------------------------------------- 
      Lpmin = LAmin+1
      Lpmax = LAmax+1
      Lmmin = MAX(0,LAmin-1)
      Lmmax = MAX(0,LAmax-1)

      ICp = Lpmax*(Lpmax+1)/2+Lpmin+1
      IF(LAmax.GT.0) THEN
         ICm = Lmmax*(Lmmax+1)/2+Lmmin+1
      ELSE
         ICm=0
      ENDIF

      LocB=Location(ICp,ICC)
      LocK=Location(ICB,ICD)
      LenB=Length(LocB)
      LenK=Length(LocK)

      DO IDer = 1,3
      DO IB=1,NB
      DO ID=1,ND
      DO IA=1,NA
      DO IC=1,NC

      CALL GetComp(ICA,IA,Lx,Ly,Lz,1)

      Ind0 = IC+NC*(IA-1)+NC*NA*(ID-1)+
     >       NC*NA*ND*(IB-1)+NC*NA*ND*NB*(IDer-1)

      ILA = ITree(1,IDer,IA,ICA)
      LocBra = IC+NC*(ILA-1)
      LocKet = ID+ND*(IB-1)
      NFAC  = NFinal(ICp)*NC
      NFBD  = NB*ND
      NUAC  = IDmn(LBra+1)
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)
      Indp = LocBra + NAC*(LocKet-1) + Ioff

      IF(ICm.NE.0) THEN
         ILA   = ITree(2,IDer,IA,ICA)
         IF(ILA.GT.0) THEN
            LocBra = IC+NC*(ILA-1)
            LocKet = ID+ND*(IB-1)
            NFAC  = NFinal(ICm)*NC
            NFBD  = NB*ND
c            NUAC  = IDmn(LBra-1)
            NUAC  = Length(Location(ICm,ICC))
            NUBD  = LenK
            NACm   = MAX(NFAC,NUAC)
            NBDm   = MAX(NFBD,NUBD)
            Indm  = LocBra+NACm*(LocKet-1) + IUN(1)
            IF(IDer.EQ.1) Facm  = Lx
            IF(IDer.EQ.2) Facm  = Ly
            IF(IDer.EQ.3) Facm  = Lz
         ELSE
            Indm = 0
            Facm = 0
         ENDIF
      ELSE
         Indm = 0
         Facm = 0
      ENDIF

      WRITE(*,9000) Ind0,Indp,Indm,Facm

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      Ioff=Ioff + NAC*NBD
c----------------------------------------------------------------------
c Center C:
c----------------------------------------------------------------------
      Lpmin = LCmin+1
      Lpmax = LCmax+1
      Lmmin = MAX(0,LCmin-1)
      Lmmax = MAX(0,LCmax-1)

      ICp = Lpmax*(Lpmax+1)/2+Lpmin+1
      IF(LCmax.GT.0) THEN
         ICm = Lmmax*(Lmmax+1)/2+Lmmin+1
      ELSE
         ICm=0
      ENDIF

      LocB=Location(ICA,ICp)
      LocK=Location(ICB,ICD)
      LenB=Length(LocB)
      LenK=Length(LocK)

      DO IDer = 1,3
      DO IB=1,NB
      DO ID=1,ND
      DO IA=1,NA
      DO IC=1,NC

      CALL GetComp(ICC,IC,Lx,Ly,Lz,1)

      Ind0 = IC+NC*(IA-1)+NC*NA*(ID-1)+
     >       NC*NA*ND*(IB-1)+NC*NA*ND*NB*(IDer-1+3)

      ILC = ITree(1,IDer,IC,ICC)
      LocBra = ILC+NFinal(ICp)*(IA-1)
      LocKet = ID+ND*(IB-1)
      NFAC  = NA*NFinal(ICp)
      NFBD  = NB*ND
      NUAC  = IDmn(LBra+1)
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)
      Indp = LocBra + NAC*(LocKet-1) + Ioff

      IF(ICm.NE.0) THEN
         ILC   = ITree(2,IDer,IC,ICC)
         IF(ILC.GT.0) THEN
            LocBra = ILC+NFinal(ICm)*(IA-1)
            LocKet = ID+ND*(IB-1)
            NFAC  = NA*NFinal(ICm)
            NFBD  = NB*ND
c            NUAC  = IDmn(LBra-1)
            NUAC  = Length(Location(ICA,ICm))
            NUBD  = LenK
            NACm   = MAX(NFAC,NUAC)
            NBDm   = MAX(NFBD,NUBD)
            Indm  = LocBra+NACm*(LocKet-1) + IUN(2)
            IF(IDer.EQ.1) Facm  = Lx
            IF(IDer.EQ.2) Facm  = Ly
            IF(IDer.EQ.3) Facm  = Lz
         ELSE
            Indm = 0
            Facm = 0
         ENDIF
      ELSE
         Indm = 0
         Facm = 0
      ENDIF

      WRITE(*,9000) Ind0,Indp,Indm,Facm

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      Ioff=Ioff + NAC*NBD
c----------------------------------------------------------------------
c Center B:
c----------------------------------------------------------------------
      Lpmin = LBmin+1
      Lpmax = LBmax+1
      Lmmin = MAX(0,LBmin-1)
      Lmmax = MAX(0,LBmax-1)

      ICp = Lpmax*(Lpmax+1)/2+Lpmin+1
      IF(LBmax.GT.0) THEN
         ICm = Lmmax*(Lmmax+1)/2+Lmmin+1
      ELSE
         ICm=0
      ENDIF


      DO IDer = 1,3
      DO IB=1,NB
      DO ID=1,ND
      DO IA=1,NA
      DO IC=1,NC

      CALL GetComp(ICB,IB,Lx,Ly,Lz,1)

      Ind0 = IC+NC*(IA-1)+NC*NA*(ID-1)+
     >       NC*NA*ND*(IB-1)+NC*NA*ND*NB*(IDer-1+6)

      ILB = ITree(1,IDer,IB,ICB)
      LocBra = IC+NC*(IA-1)
      LocKet = ID+ND*(ILB-1)
      NFAC  = NA*NC
      NFBD  = NFinal(ICp)*ND
      NUAC  = IDmn(LBra)
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)
      Indp = LocBra + NAC*(LocKet-1) + Ioff

      IF(ICm.NE.0) THEN
         ILB   = ITree(2,IDer,IB,ICB)
         IF(ILB.GT.0) THEN
            LocBra = IC+NC*(IA-1)
            LocKet = ID+ND*(ILB-1)
            NFAC  = NA*NC
            NFBD  = NFinal(ICm)*ND
c            NUAC  = IDmn(LBra)
            NUAC  = Length(Location(ICA,ICC))
            NUBD  = LenK
            NACm   = MAX(NFAC,NUAC)
            NBDm   = MAX(NFBD,NUBD)
            Indm  = LocBra+NACm*(LocKet-1) + IUN(3)
            IF(IDer.EQ.1) Facm  = Lx
            IF(IDer.EQ.2) Facm  = Ly
            IF(IDer.EQ.3) Facm  = Lz
         ELSE
            Indm = 0
            Facm = 0
         ENDIF
      ELSE
         Indm = 0
         Facm = 0
      ENDIF

      WRITE(*,9000) Ind0,Indp,Indm,Facm

      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      RETURN
 9000 FORMAT(I4,1X,I4,1X,I4,1X,I4)
      END

      SUBROUTINE GetIndexT(ITree)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER ITree(2,3,10,10)

      DO L=1,10
      DO K=1,10
      DO J=1,3
      DO I=1,2
         ITree(I,J,K,L)=0
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO Lmax=0,2
      DO Lmin=0,Lmax
         IC=Lmax*(Lmax+1)/2+Lmin+1
         ICp=(Lmax+1)*(Lmax+2)/2+Lmin+2
         IF(Lmax.GT.0) THEN
            IF(Lmin.GT.0) THEN
               ICm=(Lmax-1)*(Lmax+0)/2+Lmin+0
            ELSE
               ICm=(Lmax-1)*(Lmax+0)/2+1
            ENDIF
         ELSE
            ICm=0
         ENDIF

         DO I0=1,NFinal(IC)

            CALL GetComp(IC,I0,Lx,Ly,Lz,1)
 
            IF(Lmax.LT.2) THEN
               CALL GetComp(ICp,Ip,Lx+1,Ly,Lz,0)
               ITree(1,1,I0,IC)=Ip
               CALL GetComp(ICp,Ip,Lx,Ly+1,Lz,0)
               ITree(1,2,I0,IC)=Ip
               CALL GetComp(ICp,Ip,Lx,Ly,Lz+1,0)
               ITree(1,3,I0,IC)=Ip
            ENDIF

            IF(Lx.NE.0) THEN
               CALL GetComp(ICm,Im,Lx-1,Ly,Lz,0)
               ITree(2,1,I0,IC)=Im
            ENDIF
            IF(Ly.NE.0) THEN
               CALL GetComp(ICm,Im,Lx,Ly-1,Lz,0)
               ITree(2,2,I0,IC)=Im
            ENDIF
            IF(Lz.NE.0) THEN
               CALL GetComp(ICm,Im,Lx,Ly,Lz-1,0)
               ITree(2,3,I0,IC)=Im
            ENDIF

         ENDDO
      ENDDO
      ENDDO

 9000 format(3I3,"|",3I3,"|",3I3,"|",3I3,"|",3I3,"|",3I3)

      RETURN
      END


      SUBROUTINE GetComp(IC,I,Lx,Ly,Lz,IType)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER Ldata(3,10,10)

      Ldata(1,1,1)=0
      Ldata(2,1,1)=0
      Ldata(3,1,1)=0

      Ldata(1,1,2)=0
      Ldata(2,1,2)=0
      Ldata(3,1,2)=0
      Ldata(1,2,2)=1
      Ldata(2,2,2)=0
      Ldata(3,2,2)=0
      Ldata(1,3,2)=0
      Ldata(2,3,2)=1
      Ldata(3,3,2)=0
      Ldata(1,4,2)=0
      Ldata(2,4,2)=0
      Ldata(3,4,2)=1

      Ldata(1,1,3)=1
      Ldata(2,1,3)=0
      Ldata(3,1,3)=0
      Ldata(1,2,3)=0
      Ldata(2,2,3)=1
      Ldata(3,2,3)=0
      Ldata(1,3,3)=0
      Ldata(2,3,3)=0
      Ldata(3,3,3)=1

      Ldata(1, 1,4)=0
      Ldata(2, 1,4)=0
      Ldata(3, 1,4)=0
      Ldata(1, 2,4)=1
      Ldata(2, 2,4)=0
      Ldata(3, 2,4)=0
      Ldata(1, 3,4)=0
      Ldata(2, 3,4)=1
      Ldata(3, 3,4)=0
      Ldata(1, 4,4)=0
      Ldata(2, 4,4)=0
      Ldata(3, 4,4)=1
      Ldata(1, 5,4)=2
      Ldata(2, 5,4)=0
      Ldata(3, 5,4)=0
      Ldata(1, 6,4)=1
      Ldata(2, 6,4)=1
      Ldata(3, 6,4)=0
      Ldata(1, 7,4)=0
      Ldata(2, 7,4)=2
      Ldata(3, 7,4)=0
      Ldata(1, 8,4)=1
      Ldata(2, 8,4)=0
      Ldata(3, 8,4)=1
      Ldata(1, 9,4)=0
      Ldata(2, 9,4)=1
      Ldata(3, 9,4)=1
      Ldata(1,10,4)=0
      Ldata(2,10,4)=0
      Ldata(3,10,4)=2

      Ldata(1,1,5)=1
      Ldata(2,1,5)=0
      Ldata(3,1,5)=0
      Ldata(1,2,5)=0
      Ldata(2,2,5)=1
      Ldata(3,2,5)=0
      Ldata(1,3,5)=0
      Ldata(2,3,5)=0
      Ldata(3,3,5)=1
      Ldata(1,4,5)=2
      Ldata(2,4,5)=0
      Ldata(3,4,5)=0
      Ldata(1,5,5)=1
      Ldata(2,5,5)=1
      Ldata(3,5,5)=0
      Ldata(1,6,5)=0
      Ldata(2,6,5)=2
      Ldata(3,6,5)=0
      Ldata(1,7,5)=1
      Ldata(2,7,5)=0
      Ldata(3,7,5)=1
      Ldata(1,8,5)=0
      Ldata(2,8,5)=1
      Ldata(3,8,5)=1
      Ldata(1,9,5)=0
      Ldata(2,9,5)=0
      Ldata(3,9,5)=2

      Ldata(1,1,6)=2
      Ldata(2,1,6)=0
      Ldata(3,1,6)=0
      Ldata(1,2,6)=1
      Ldata(2,2,6)=1
      Ldata(3,2,6)=0
      Ldata(1,3,6)=0
      Ldata(2,3,6)=2
      Ldata(3,3,6)=0
      Ldata(1,4,6)=1
      Ldata(2,4,6)=0
      Ldata(3,4,6)=1
      Ldata(1,5,6)=0
      Ldata(2,5,6)=1
      Ldata(3,5,6)=1
      Ldata(1,6,6)=0
      Ldata(2,6,6)=0
      Ldata(3,6,6)=2

      IF(IType.EQ.1) THEN
         Lx = Ldata(1,I,IC)
         Ly = Ldata(2,I,IC)
         Lz = Ldata(3,I,IC)
      ELSE
         N=NFinal(IC)
         I=0
         DO J=1,N
            IF(Lx.EQ.Ldata(1,J,IC).AND.
     >         Ly.EQ.Ldata(2,J,IC).AND.
     >         Lz.EQ.Ldata(3,J,IC)) I=J
         ENDDO
         IF(I.EQ.0) THEN
            WRITE(*,*) "IC=",IC
            WRITE(*,*) "L=",Lx,Ly,Lz
            WRITE(*,*) "Messed up!"
            STOP
         ENDIF
      ENDIF
      RETURN
      END

      SUBROUTINE CIndx(LType,IUN,IType)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER LType(8)
      INTEGER IUN(3),ISC(3)
      INTEGER BraCod,KetCod

      IUN(1) = 0
      IUN(2) = 0
      IUN(3) = 0
      ISC(1) = 0
      ISC(2) = 0
      ISC(3) = 0

      LineN=0

      LBra = LType(2)+LType(4)
      LKet = LType(6)+LType(8)
      IBra = IDmn(LBra+1)
      IKet = IDmn(LKet+1)
c
c  ZA Scaled:
c
      IBraF = IDmn(LBra+1)
      LAmin = LType(1)+1
      LAmax = LType(2)+1
      LCmin = LType(3)
      LCmax = LType(4)
      LBmin = LType(5)
      LBmax = LType(6)
      LDmin = LType(7)
      LDmax = LType(8)
      ICA   = LAmax*(LAmax+1)/2+LAmin+1
      ICC   = LCmax*(LCmax+1)/2+LCmin+1
      ICB   = LBmax*(LBmax+1)/2+LBmin+1
      ICD   = LDmax*(LDmax+1)/2+LDmin+1
      NA    = NFinal(ICA)
      NC    = NFinal(ICC)
      NB    = NFinal(ICB)
      ND    = NFinal(ICD)

      NFAC  = NA*NC
      NFBD  = NB*ND
      NUAC  = IBraF
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)

      BraCod=ICA*100+ICC
      KetCod=ICB*100+ICD
      WRITE(72,9000) KetCod,NAC,0,1
      WRITE(73,9000) BraCod,NAC,NB*ND,1

      LocB = Location(ICA,ICC)
      LocK = Location(ICB,ICD)
      LenB = Length(LocB)
      LenK = Length(LocK)

      IF(IType.EQ.0) WRITE(*,*) "  ZA Scaled"
      DO J=1,LenK
         DO I=1,LenB
            K   = ILocPS(LocB,I)
            L   = ILocPS(LocK,J)
            LC1 = I+IBraF*(J-1)
            LC2 = K+IBra*(L-1)
            LCB = ILocCC(LocB,I)
            LCK = ILocCC(LocK,J)
 
            LineN=LineN+1
            IF(IType.EQ.0) WRITE(*,9000) LC1,LCB,LCK,LC2
         ENDDO
      ENDDO
      ISC(1) = LineN
c
c  ZC Scaled:
c
      IBraF = IDmn(LBra+1)
      NFAC  = NA*NC
      NFBD  = NB*ND
      NUAC  = IBraF
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)

      Ioff  = NAC*NBD
 
      LAmin = LType(1)
      LAmax = LType(2)
      LCmin = LType(3)+1
      LCmax = LType(4)+1
      LBmin = LType(5)
      LBmax = LType(6)
      LDmin = LType(7)
      LDmax = LType(8)
      ICA   = LAmax*(LAmax+1)/2+LAmin+1
      ICC   = LCmax*(LCmax+1)/2+LCmin+1
      ICB   = LBmax*(LBmax+1)/2+LBmin+1
      ICD   = LDmax*(LDmax+1)/2+LDmin+1
      NA    = NFinal(ICA)
      NC    = NFinal(ICC)
      NB    = NFinal(ICB)
      ND    = NFinal(ICD)

      NFAC  = NA*NC
      NFBD  = NB*ND
      NUAC  = IBraF
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)
      BraCod=ICA*100+ICC
      KetCod=ICB*100+ICD
      WRITE(72,9000) KetCod,NAC,0,Ioff+1
      WRITE(73,9000) BraCod,NAC,NB*ND,Ioff+1

      LocB = Location(ICA,ICC)
      LocK = Location(ICB,ICD)
      LenB = Length(LocB)
      LenK = Length(LocK)

      IF(IType.EQ.0) WRITE(*,*) "  ZC Scaled"
      DO J=1,LenK
         DO I=1,LenB
            K    = ILocPS(LocB,I)
            L    = ILocPS(LocK,J)
            LC1 = I+IBraF*(J-1)
            LC2 = K+IBra*(L-1)
            LCB = ILocCC(LocB,I)
            LCK = ILocCC(LocK,J)
            LineN=LineN+1
            IF(IType.EQ.0) WRITE(*,9000) LC1+Ioff,LCB,LCK,LC2
         ENDDO
      ENDDO
      ISC(2) = LineN
c
c  ZB Scaled:
c
      NFAC  = NA*NC
      NFBD  = NB*ND
      NUAC  = IBraF
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)

      Ioff  = Ioff+NAC*NBD
      IBraF = IDmn(LBra)
      LAmin = LType(1)
      LAmax = LType(2)
      LCmin = LType(3)
      LCmax = LType(4)
      LBmin = LType(5)+1
      LBmax = LType(6)+1
      LDmin = LType(7)
      LDmax = LType(8)
      ICA   = LAmax*(LAmax+1)/2+LAmin+1
      ICC   = LCmax*(LCmax+1)/2+LCmin+1
      ICB   = LBmax*(LBmax+1)/2+LBmin+1
      ICD   = LDmax*(LDmax+1)/2+LDmin+1
      NA    = NFinal(ICA)
      NC    = NFinal(ICC)
      NB    = NFinal(ICB)
      ND    = NFinal(ICD)

      NFAC  = NA*NC
      NFBD  = NB*ND
      NUAC  = IBraF
      NUBD  = LenK
      NAC   = MAX(NFAC,NUAC)
      NBD   = MAX(NFBD,NUBD)
      BraCod=ICA*100+ICC
      KetCod=ICB*100+ICD
      WRITE(72,9000) KetCod,NAC,0,Ioff+1
      WRITE(73,9000) BraCod,NAC,NB*ND,Ioff+1

      LocB = Location(ICA,ICC)
      LocK = Location(ICB,ICD)
      LenB = Length(LocB)
      LenK = Length(LocK)

      IF(IType.EQ.0) WRITE(*,*) "  ZB Scaled"
      DO J=1,LenK
         DO I=1,LenB
            K    = ILocPS(LocB,I)
            L    = ILocPS(LocK,J)
            LC1 = I+IBraF*(J-1)
            LC2 = K+IBra*(L-1)
            LCB = ILocCC(LocB,I)
            LCK = ILocCC(LocK,J)
            LineN=LineN+1
            IF(IType.EQ.0) WRITE(*,9000) LC1+Ioff,LCB,LCK,LC2
         ENDDO
      ENDDO
      ISC(3) = LineN

c---------------------------------------------------------------------- 
c  Unscaled:
c---------------------------------------------------------------------- 
      IF(LType(2).GT.0) THEN
         NFAC  = NA*NC
         NFBD  = NB*ND
         NUAC  = IBraF
         NUBD  = LenK
         NAC   = MAX(NFAC,NUAC)
         NBD   = MAX(NFBD,NUBD)
         Ioff  = Ioff+NAC*NBD
         IUN(1) = Ioff
         LAmin = LType(1)-MAX(0,LType(1)-1)
         LAmax = LType(2)-1
         LCmin = LType(3)
         LCmax = LType(4)
         LBmin = LType(5)
         LBmax = LType(6)
         LDmin = LType(7)
         LDmax = LType(8)
         ICA   = LAmax*(LAmax+1)/2+LAmin+1
         ICC   = LCmax*(LCmax+1)/2+LCmin+1
         ICB   = LBmax*(LBmax+1)/2+LBmin+1
         ICD   = LDmax*(LDmax+1)/2+LDmin+1
         NA    = NFinal(ICA)
         NC    = NFinal(ICC)
         NB    = NFinal(ICB)
         ND    = NFinal(ICD)

         LocB = Location(ICA,ICC)
         LocK = Location(ICB,ICD)
         LenB = Length(LocB)
         LenK = Length(LocK)

         IBraF=MAX(IDmn(LBra-1),Length(LocB))

         NFAC  = NA*NC
         NFBD  = NB*ND
         NUAC  = IBraF
         NUBD  = LenK
         NAC   = MAX(NFAC,NUAC)
         NBD   = MAX(NFBD,NUBD)
         BraCod=ICA*100+ICC
         KetCod=ICB*100+ICD
         WRITE(72,9000) KetCod,NAC,0,Ioff+1
         WRITE(73,9000) BraCod,NAC,NB*ND,Ioff+1

c         IF(IType.EQ.0) WRITE(*,*) "  Center A unscaled"
         DO J=1,LenK
            DO I=1,LenB
               K   = ILocPS(LocB,I)
               L   = ILocPS(LocK,J)
               LC1 = I+IBraF*(J-1)
               LC2 = K+IBra*(L-1)
               LCB = ILocCCU(LocB,I)
               LCK = ILocCC(LocK,J)
               LineN=LineN+1
               IF(IType.EQ.0) WRITE(*,9000) LC1+Ioff,LCB,LCK,LC2
            ENDDO
         ENDDO
      ENDIF

      IF(LType(4).GT.0) THEN
         NFAC  = NA*NC
         NFBD  = NB*ND
         NUAC  = IBraF
         NUBD  = LenK
         NAC   = MAX(NFAC,NUAC)
         NBD   = MAX(NFBD,NUBD)
         Ioff  = Ioff+NAC*NBD
         IUN(2) = Ioff
         LAmin = LType(1)
         LAmax = LType(2)
         LCmin = LType(3)-MAX(0,LType(3)-1)
         LCmax = LType(4)-1
         LBmin = LType(5)
         LBmax = LType(6)
         LDmin = LType(7)
         LDmax = LType(8)
         ICA   = LAmax*(LAmax+1)/2+LAmin+1
         ICC   = LCmax*(LCmax+1)/2+LCmin+1
         ICB   = LBmax*(LBmax+1)/2+LBmin+1
         ICD   = LDmax*(LDmax+1)/2+LDmin+1
         NA    = NFinal(ICA)
         NC    = NFinal(ICC)
         NB    = NFinal(ICB)
         ND    = NFinal(ICD)

         LocB = Location(ICA,ICC)
         LocK = Location(ICB,ICD)
         LenB = Length(LocB)
         LenK = Length(LocK)

         IBraF=MAX(IDmn(LBra-1),Length(LocB))

         NFAC  = NA*NC
         NFBD  = NB*ND
         NUAC  = IBraF
         NUBD  = LenK
         NAC   = MAX(NFAC,NUAC)
         NBD   = MAX(NFBD,NUBD)
         BraCod=ICA*100+ICC
         KetCod=ICB*100+ICD
         WRITE(72,9000) KetCod,NAC,0,Ioff+1
         WRITE(73,9000) BraCod,NAC,NB*ND,Ioff+1

c         IF(IType.EQ.0) WRITE(*,*) "  Center C unscaled"
         DO J=1,LenK
            DO I=1,LenB
               K   = ILocPS(LocB,I)
               L   = ILocPS(LocK,J)
               LC1 = I+IBraF*(J-1)
               LC2 = K+IBra*(L-1)
               LCB = ILocCCU(LocB,I)
               LCK = ILocCC(LocK,J)
               LineN=LineN+1
               IF(IType.EQ.0) WRITE(*,9000) LC1+Ioff,LCB,LCK,LC2
            ENDDO
         ENDDO
      ENDIF

      IF(LType(6).GT.0) THEN
         NFAC  = NA*NC
         NFBD  = NB*ND
         NUAC  = IBraF
         NUBD  = LenK
         NAC   = MAX(NFAC,NUAC)
         NBD   = MAX(NFBD,NUBD)
         Ioff  = Ioff+NAC*NBD
         IUN(3) = Ioff
         IBraF = IDmn(LBra)
         LAmin = LType(1)
         LAmax = LType(2)
         LCmin = LType(3)
         LCmax = LType(4)
         LBmin = LType(5)-MAX(0,LType(5)-1)
         LBmax = LType(6)-1
         LDmin = LType(7)
         LDmax = LType(8)
         ICA   = LAmax*(LAmax+1)/2+LAmin+1
         ICC   = LCmax*(LCmax+1)/2+LCmin+1
         ICB   = LBmax*(LBmax+1)/2+LBmin+1
         ICD   = LDmax*(LDmax+1)/2+LDmin+1
         NA    = NFinal(ICA)
         NC    = NFinal(ICC)
         NB    = NFinal(ICB)
         ND    = NFinal(ICD)

         NFAC  = NA*NC
         NFBD  = NB*ND
         NUAC  = IBraF
         NUBD  = LenK
         NAC   = MAX(NFAC,NUAC)
         NBD   = MAX(NFBD,NUBD)
         BraCod=ICA*100+ICC
         KetCod=ICB*100+ICD
         WRITE(72,9000) KetCod,NAC,0,Ioff+1
         WRITE(73,9000) BraCod,NAC,NB*ND,Ioff+1

         LocB = Location(ICA,ICC)
         LocK = Location(ICB,ICD)
         LenB = Length(LocB)
         LenK = Length(LocK)

c         IF(IType.EQ.0) WRITE(*,*) "  Center B unscaled"
         DO J=1,LenK
            DO I=1,LenB
               K   = ILocPS(LocB,I)
               L   = ILocPS(LocK,J)
               LC1 = I+IBraF*(J-1)
               LC2 = K+IBra*(L-1)
               LCB = ILocCC(LocB,I)
               LCK = ILocCCU(LocK,J)
               LineN=LineN+1
               IF(IType.EQ.0) WRITE(*,9000) LC1+Ioff,LCB,LCK,LC2
            ENDDO
         ENDDO
      ENDIF

      IF(IType.EQ.0) THEN
         WRITE(*,*) LineN*4
         WRITE(*,*) ISC(1)
         WRITE(*,*) ISC(2)
         WRITE(*,*) ISC(3)
         WRITE(*,*) LineN
      ENDIF

      RETURN
 9000 FORMAT(I4,1X,I4,1X,I4,1X,I4)
      END


      FUNCTION ILocPS(IP,I)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER IL(10,50)
c
c (p,s|
c      
      IL(1,1) = 2
      IL(1,2) = 3
      IL(1,3) = 4
c
c (s,p|
c
      IL(2,1) = 1
      IL(2,2) = 2
      IL(2,3) = 3
      IL(2,4) = 4
c
c (pd,s|
c
      IL(3,1) =  2
      IL(3,2) =  3
      IL(3,3) =  4
      IL(3,4) =  5
      IL(3,5) =  6
      IL(3,6) =  7
      IL(3,7) =  8
      IL(3,8) =  9
      IL(3,9) = 10
c
c (sp,p|
c
      IL(4,1)  =  1
      IL(4,2)  =  2
      IL(4,3)  =  3
      IL(4,4)  =  4
      IL(4,5)  =  2
      IL(4,6)  =  3
      IL(4,7)  =  4
      IL(4,8)  =  5
      IL(4,9)  =  6
      IL(4,10) =  7
      IL(4,11) =  8
      IL(4,12) =  9
      IL(4,13) = 10
c
c (pd,sp|
c
      IL(5,1)  =  2
      IL(5,2)  =  3
      IL(5,3)  =  4
      IL(5,4)  =  5
      IL(5,5)  =  6
      IL(5,6)  =  7
      IL(5,7)  =  8
      IL(5,8)  =  9
      IL(5,9)  = 10
      IL(5,10) =  2
      IL(5,11) =  3
      IL(5,12) =  4
      IL(5,13) =  5
      IL(5,14) =  6
      IL(5,15) =  7
      IL(5,16) =  8
      IL(5,17) =  9
      IL(5,18) = 10
      IL(5,19) =  5
      IL(5,20) =  6
      IL(5,21) =  7
      IL(5,22) =  8
      IL(5,23) =  9
      IL(5,24) = 10
      IL(5,25) = 11
      IL(5,26) = 12
      IL(5,27) = 13
      IL(5,28) = 14
      IL(5,29) = 15
      IL(5,30) = 16
      IL(5,31) = 17
      IL(5,32) = 18
      IL(5,33) = 19
      IL(5,34) = 20
c
c (sp,pd|
c
      IL(6, 1) =  1
      IL(6, 2) =  2
      IL(6, 3) =  3
      IL(6, 4) =  4
      IL(6, 5) =  2
      IL(6, 6) =  3
      IL(6, 7) =  4
      IL(6, 8) =  5
      IL(6, 9) =  6
      IL(6,10) =  7
      IL(6,11) =  8
      IL(6,12) =  9
      IL(6,13) = 10
      IL(6,14) =  1
      IL(6,15) =  2
      IL(6,16) =  3
      IL(6,17) =  4
      IL(6,18) =  5
      IL(6,19) =  6
      IL(6,20) =  7
      IL(6,21) =  8
      IL(6,22) =  9
      IL(6,23) = 10
      IL(6,24) =  2
      IL(6,25) =  3
      IL(6,26) =  4
      IL(6,27) =  5
      IL(6,28) =  6
      IL(6,29) =  7
      IL(6,30) =  8
      IL(6,31) =  9
      IL(6,32) = 10
      IL(6,33) = 11
      IL(6,34) = 12
      IL(6,35) = 13
      IL(6,36) = 14
      IL(6,37) = 15
      IL(6,38) = 16
      IL(6,39) = 17
      IL(6,40) = 18
      IL(6,41) = 19
      IL(6,42) = 20
c
c (s,s|
c
      IL(7,1) = 1
c
c (sp,s|
c
      IL(8,1)=1
      IL(8,2)=2
      IL(8,3)=3
      IL(8,4)=4
c
c (sp,sp|
c
      IL(9, 1)=1
      IL(9, 2)=2
      IL(9, 3)=3
      IL(9, 4)=4
      IL(9, 5)=2
      IL(9, 6)=5
      IL(9, 7)=1
      IL(9, 8)=2
      IL(9, 9)=3
      IL(9,10)=6
      IL(9,11)=7
      IL(9,12)=3
      IL(9,13)=4
      IL(9,14)=8
      IL(9,15)=9
      IL(9,16)=10
      IL(9,17)=4
c
c (s,sp|
c
      IL(10,1) = 1
      IL(10,2) = 1
      IL(10,3) = 2
      IL(10,4) = 3
      IL(10,5) = 4


      ILocPS = IL(IP,I)
      RETURN
      END

      FUNCTION ILocCCU(IP,I)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER  IC(10,50)
c
c (s,s|
c
      IC(7,1) = 0
c
c (sp,s|
c
      IC(8,1) = 3
      IC(8,2) = 0
      IC(8,3) = 0
      IC(8,4) = 0
c
c (s,sp|
c
      IC(10,1) = 2
      IC(10,2) = 0
      IC(10,3) = 0
      IC(10,4) = 0
      IC(10,5) = 0


      IF(IP.NE.7.AND.IP.NE.8.AND.IP.NE.10) THEN
         WRITE(*,*) "IP=",IP
         WRITE(*,*) "Messed up in ILocCCU"
         STOP
      ENDIF

      ILocCCU = IC(IP,I)
      RETURN
      END

      

      FUNCTION ILocCC(IP,I)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER  IC(10,50)
c
c (p,s|
c
      IC(1,1) = 0
      IC(1,2) = 0
      IC(1,3) = 0
c
c (s,p|
c
      IC(2,1) = 0
      IC(2,2) = 0
      IC(2,3) = 0
      IC(2,4) = 0

c
c (pd,s|
c
      IC(3,1) = 1
      IC(3,2) = 1
      IC(3,3) = 1
      IC(3,4) = 0
      IC(3,5) = 0
      IC(3,6) = 0
      IC(3,7) = 0
      IC(3,8) = 0
      IC(3,9) = 0

c
c (sp,p|
c
      IC(4,1)  = 1
      IC(4,2)  = 1
      IC(4,3)  = 1
      IC(4,4)  = 1
      IC(4,5)  = 0
      IC(4,6)  = 0
      IC(4,7)  = 0
      IC(4,8)  = 0
      IC(4,9)  = 0
      IC(4,10) = 0
      IC(4,11) = 0
      IC(4,12) = 0
      IC(4,13) = 0
c
c (pd,sp|
c
      IC(5,1)  = 1
      IC(5,2)  = 1
      IC(5,3)  = 1
      IC(5,4)  = 2
      IC(5,5)  = 2
      IC(5,6)  = 2
      IC(5,7)  = 2
      IC(5,8)  = 2
      IC(5,9)  = 2
      IC(5,10) = 3
      IC(5,11) = 3
      IC(5,12) = 3
      IC(5,13) = 3
      IC(5,14) = 3
      IC(5,15) = 3
      IC(5,16) = 3
      IC(5,17) = 3
      IC(5,18) = 3
      IC(5,19) = 0
      IC(5,20) = 0
      IC(5,21) = 0
      IC(5,22) = 0
      IC(5,23) = 0
      IC(5,24) = 0
      IC(5,25) = 0
      IC(5,26) = 0
      IC(5,27) = 0
      IC(5,28) = 0
      IC(5,29) = 0
      IC(5,30) = 0
      IC(5,31) = 0
      IC(5,32) = 0
      IC(5,33) = 0
      IC(5,34) = 0

c
c (sp,pd|
c
      IC(6,1)  =  1
      IC(6,2)  =  1
      IC(6,3)  =  1
      IC(6,4)  =  1
      IC(6,5)  =  2
      IC(6,6)  =  2
      IC(6,7)  =  2
      IC(6,8)  =  2
      IC(6,9)  =  2
      IC(6,10) =  2
      IC(6,11) =  2
      IC(6,12) =  2
      IC(6,13) =  2
      IC(6,14) =  3
      IC(6,15) =  3
      IC(6,16) =  3
      IC(6,17) =  3
      IC(6,18) =  3
      IC(6,19) =  3
      IC(6,20) =  3
      IC(6,21) =  3
      IC(6,22) =  3
      IC(6,23) =  3
      IC(6,24) =  0
      IC(6,25) =  0
      IC(6,26) =  0
      IC(6,27) =  0
      IC(6,28) =  0
      IC(6,29) =  0
      IC(6,30) =  0
      IC(6,31) =  0
      IC(6,32) =  0
      IC(6,33) =  0
      IC(6,34) =  0
      IC(6,35) =  0
      IC(6,36) =  0
      IC(6,37) =  0
      IC(6,38) =  0
      IC(6,39) =  0
      IC(6,40) =  0
      IC(6,41) =  0
      IC(6,42) =  0

c
c (s,s|
c
      IC(7,1) = 0

c
c (sp,s|
c
      IC(8,1)=1
      IC(8,2)=0
      IC(8,3)=0
      IC(8,4)=0

c
c (sp,sp|
c
      IC(9, 1)=1
      IC(9, 2)=3
      IC(9, 3)=3
      IC(9, 4)=3
      IC(9, 5)=2
      IC(9, 6)=0
      IC(9, 7)=3
      IC(9, 8)=0
      IC(9, 9)=2
      IC(9,10)=0
      IC(9,11)=0
      IC(9,12)=0
      IC(9,13)=2
      IC(9,14)=0
      IC(9,15)=0
      IC(9,16)=0 
      IC(9,17)=0

      ILocCC = IC(IP,I)
      RETURN
      END


      FUNCTION Location(I,J)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER LocData(10,10)
      DATA LocData / 7, 8, 1, 0, 3, 0, 0, 0, 0, 0,
     >              10, 9, 0, 0, 5, 0, 0, 0, 0, 0, 
     >               2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
     >               0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      IF(I.LE.10.AND.J.LE.10) THEN
         Location = LocData(I,J)
         IF(Location.EQ.0) THEN
            WRITE(*,*) "I=",I," J=",J
            WRITE(*,*) "Messed up in Location"
            STOP
         ENDIF
      ELSE
         WRITE(*,*) "I=",I," J=",J
         WRITE(*,*) "Messed up in Location"
         STOP
      ENDIF
      RETURN
      END
      
      FUNCTION Length(I)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER LenData(10)
      DATA LenData / 3, 4, 9, 13, 34, 42, 1, 4, 17, 5 /

      IF(I.LE.10) THEN
         Length=LenData(I)
      ELSE
         WRITE(*,*) "I=",I
         WRITE(*,*) "Messed up in Length"
         STOP
      ENDIF
      RETURN
      END

      FUNCTION IDmn(L)
c
c Number of locations needed to contract the VRRs
c
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER DimData(0:4)
      DATA DimData /1,4,17,24,36/
c      DATA DimData  /1,4,17,42 /
      
      IF(L.LE.4) THEN
         IDmn=DimData(L)
      ELSE
         STOP
      ENDIF
      RETURN
      END
      

      FUNCTION NFinal(IType)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER NFinal,IType
      INTEGER T2Size(28)
      SAVE T2Size
      DATA T2Size /  1,
     >               4, 3,
     >              10, 9, 6,
     >              20,19,16,10,
     >              35,34,31,25,15,
     >              56,55,52,46,36,21,
     >              84,83,80,74,64,49,28 /
      IF(IType.LE.28) THEN
         NFinal=T2Size(IType)
      ELSE
         WRITE(*,*) "IType=",IType
         STOP
      ENDIF
      RETURN
      END

