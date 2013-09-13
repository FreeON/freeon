!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE MBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  USE SpecFun
  IMPLICIT NONE
CONTAINS
  FUNCTION DBlok(BS,MD,Pair,IXYZ,COrig) RESULT(DVck)
    TYPE(BSET)                                :: BS
    TYPE(DBL_RNK4)                            :: MD
    TYPE(AtomPair)                            :: Pair
    TYPE(DBL_VECT), INTENT(IN)                :: COrig
    INTEGER       , INTENT(IN)                :: IXYZ
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)   :: DVck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: DBlk
!
    INTEGER                                   :: NBFA,NBFB,KA,KB,KK,AtA,AtB
    REAL(DOUBLE)                              :: POx,POy,POz
    REAL(DOUBLE)                              :: Ax,Ay,Az,Bx,By,Bz,AB2
    REAL(DOUBLE)                              :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
                                                 EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB
    INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                                 StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                                 LA,LB,MA,MB,NA,NB
    !
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
    Ax = Pair%A(1)
    Ay = Pair%A(2)
    Az = Pair%A(3)
    Bx = Pair%B(1)
    By = Pair%B(2)
    Bz = Pair%B(3)
    AB2=Pair%AB2

    DBlk(:,:)=Zero
!
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)

       DO CFB=1,BS%NCFnc%I(KB)                    ! Loop over contracted function B
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)

          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in
             DO PFB=1,BS%NPFnc%I(CFB,KB)          ! contracted functions A and B
                ZetaA=BS%Expnt%D(PFA,CFA,KA)
                ZetaB=BS%Expnt%D(PFB,CFB,KB)
                EtaAB=ZetaA+ZetaB
                EtaIn=One/EtaAB
                ZAB  =ZetaA*ZetaB
                XiAB =ZAB*EtaIn
                ExpAB=EXP(-XiAB*AB2)
                PiE32=(Pi/EtaAB)**(1.5D0)
                Ov=ExpAB*PiE32
                Px=(ZetaA*Ax+ZetaB*Bx)*EtaIn
                Py=(ZetaA*Ay+ZetaB*By)*EtaIn
                Pz=(ZetaA*Az+ZetaB*Bz)*EtaIn
                PAx=Px-Ax
                PAy=Py-Ay
                PAz=Pz-Az
                PBx=Px-Bx
                PBy=Py-By
                PBz=Pz-Bz
                !
                MD%D(1,0,0,1)=0.00D+00
                MD%D(2,0,0,1)=0.00D+00
                MD%D(3,0,0,1)=0.00D+00
                CALL MD2TRR(BS%NASym+2,-1,MaxLA+1,MaxLB+1,EtaAB,MD%D, &
                            PAx,PBx,PAy,PBy,PAz,PBz)
                !vw--->
                POx = Px-COrig%D(1)
                POy = Py-COrig%D(2)
                POz = Pz-COrig%D(3)
                !vw<---
                IB=IndexB
                DO LMNB=StartLB,StopLB
                   LB=BS%LxDex%I(LMNB)
                   MB=BS%LyDex%I(LMNB)
                   NB=BS%LzDex%I(LMNB)
                   CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                   IB=IB+1
                   IA=IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      LA=BS%LxDex%I(LMNA)
                      MA=BS%LyDex%I(LMNA)
                      NA=BS%LzDex%I(LMNA)
                      CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                      !vw--->
                      SELECT CASE(IXYZ)
                      CASE(1)
                         DBlk(IA,IB)=DBlk(IA,IB) &
                                    +CA*CB*Ov*((POx*MD%D(1,LA,LB,0)+MD%D(1,LA,LB,1))*MD%D(2,MA,MB,0)*MD%D(3,NA,NB,0))
                      CASE(2)
                         DBlk(IA,IB)=DBlk(IA,IB) &
                                    +CA*CB*Ov*((POy*MD%D(2,MA,MB,0)+MD%D(2,MA,MB,1))*MD%D(1,LA,LB,0)*MD%D(3,NA,NB,0))
                      CASE(3)
                         DBlk(IA,IB)=DBlk(IA,IB) &
                                    +CA*CB*Ov*((POz*MD%D(3,NA,NB,0)+MD%D(3,NA,NB,1))*MD%D(1,LA,LB,0)*MD%D(2,MA,MB,0))
                      END SELECT
                      !vw<---
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DVck = BlockToVect(NBFA,NBFB,DBlk)
    !
  END FUNCTION DBlok
  !
  !
!!$  FUNCTION DipolCorrect(BS,GM,Pair) RESULT(DVck)
!!$    TYPE(BSET)                                :: BS
!!$    TYPE(AtomPair)                            :: Pair
!!$    TYPE(CRDS), INTENT(IN)                    :: GM
!!$    !
!!$    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)   :: DVck
!!$    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: DBlk
!!$    !
!!$    INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
!!$                                                 StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB, &
!!$                                                 L,M,N,L1,L2,LMN,NBFA,NBFB,KA,KB
!!$    TYPE(PrimPair)                            :: Prim
!!$    REAL(DOUBLE)                              :: Qx,Qy,Qz,POx,POy,POz,Pop,SqZ,UQx,UQy,UQz,LQx,LQy,LQz, &
!!$                                                 CoFact,Z,LQ2,UQ2,RL2,TwoZ,Amp,AB2
!!$    REAL(DOUBLE), DIMENSION(0:HGEll+1)        :: LLambdaX,LLambdaY,LLambdaZ, &
!!$                                                 ULambdaX,ULambdaY,ULambdaZ, &
!!$                                                 LambdaX,LambdaY,LambdaZ
!!$    !
!!$    KA   = Pair%KA
!!$    KB   = Pair%KB
!!$    NBFA = Pair%NA
!!$    NBFB = Pair%NB
!!$    AB2=Pair%AB2
!!$    DBlk(:,:)=Zero
!!$!----------------------------------
!!$    Prim%A=Pair%A
!!$    Prim%B=Pair%B
!!$    Prim%AB2=Pair%AB2
!!$    Prim%KA=Pair%KA
!!$    Prim%KB=Pair%KB
!!$!----------------------------------
!!$    DO CFA=1,BS%NCFnc%I(KA)
!!$       IndexA=CFBlokDex(BS,CFA,KA)
!!$       StartLA=BS%LStrt%I(CFA,KA)
!!$       StopLA =BS%LStop%I(CFA,KA)
!!$       MaxLA=BS%ASymm%I(2,CFA,KA)
!!$       DO CFB=1,BS%NCFnc%I(KB)
!!$          IndexB  = CFBlokDex(BS,CFB,KB)
!!$          StartLB = BS%LStrt%I(CFB,KB)
!!$          StopLB  = BS%LStop%I(CFB,KB)
!!$          MaxLB   = BS%ASymm%I(2,CFB,KB)
!!$!----------------------------------
!!$          Prim%CFA=CFA
!!$          Prim%CFB=CFB
!!$          Prim%Ell=MaxLA+MaxLB
!!$!----------------------------------
!!$          DO PFA=1,BS%NPFnc%I(CFA,KA)
!!$             DO PFB=1,BS%NPFnc%I(CFB,KB)
!!$!----------------------------------
!!$                Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
!!$                Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
!!$                Prim%Zeta=Prim%ZA+Prim%ZB
!!$                Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
!!$                Qx=(Prim%ZA*Prim%A(1)+Prim%ZB*Prim%B(1))/Prim%Zeta
!!$                Qy=(Prim%ZA*Prim%A(2)+Prim%ZB*Prim%B(2))/Prim%Zeta
!!$                Qz=(Prim%ZA*Prim%A(3)+Prim%ZB*Prim%B(3))/Prim%Zeta
!!$                IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
!!$                   Prim%PFA=PFA
!!$                   Prim%PFB=PFB
!!$                   Amp=SetBraBlok(Prim,BS)
!!$!----------------------------------
!!$                   IB=IndexB
!!$                   DO LMNB=StartLB,StopLB
!!$                      IB=IB+1
!!$                      IA=IndexA
!!$                      DO LMNA=StartLA,StopLA
!!$                         IA=IA+1
!!$                         Pop=Zero
!!$                         Z=Prim%Zeta
!!$                         TwoZ=Two*Z
!!$                         SqZ=SQRT(Z)
!!$                         CoFact=SqrtPi/(Two*SqZ)
!!$                         ! Lower bounds
!!$                         LQx=GM%PBC%BoxShape%D(1,2)-Qx
!!$                         LQy=-500.0d0-Qy
!!$                         LQz=-500.0d0-Qz
!!$                         LQ2=LQx*LQx+LQy*LQy+LQz*LQz
!!$                         ! Lower Lambdas
!!$                         LLambdaX(0)=-CoFact*ERF(SqZ*LQx)
!!$                         LLambdaY(0)=-CoFact*ERF(SqZ*LQy)
!!$                         LLambdaZ(0)=-CoFact*ERF(SqZ*LQz)
!!$                         LLambdaX(1)=EXP(-Z*LQx*LQx)
!!$                         LLambdaY(1)=EXP(-Z*LQy*LQy)
!!$                         LLambdaZ(1)=EXP(-Z*LQz*LQz)
!!$                         ! Upper bounds
!!$                         UQx=GM%PBC%BoxShape%D(1,1)-Qx
!!$                         UQy=500.0D0-Qy
!!$                         UQz=500.0D0-Qz
!!$                         UQ2=UQx*UQx+UQy*UQy+UQz*UQz
!!$                         ! Upper Lambdas
!!$                         ULambdaX(0)=-CoFact*ERF(SqZ*UQx)
!!$                         ULambdaY(0)=-CoFact*ERF(SqZ*UQy)
!!$                         ULambdaZ(0)=-CoFact*ERF(SqZ*UQz)
!!$                         ULambdaX(1)=EXP(-Z*UQx*UQx)
!!$                         ULambdaY(1)=EXP(-Z*UQy*UQy)
!!$                         ULambdaZ(1)=EXP(-Z*UQz*UQz)
!!$                         ! Generic Lambdas
!!$                         DO L=2,Prim%Ell+1
!!$                            L1=L-1
!!$                            L2=L-2
!!$                            RL2=DBLE(L2)
!!$                            LLambdaX(L)=TwoZ*(LQx*LLambdaX(L1)-RL2*LLambdaX(L2))
!!$                            LLambdaY(L)=TwoZ*(LQy*LLambdaY(L1)-RL2*LLambdaY(L2))
!!$                            LLambdaZ(L)=TwoZ*(LQz*LLambdaZ(L1)-RL2*LLambdaZ(L2))
!!$                            ULambdaX(L)=TwoZ*(UQx*ULambdaX(L1)-RL2*ULambdaX(L2))
!!$                            ULambdaY(L)=TwoZ*(UQy*ULambdaY(L1)-RL2*ULambdaY(L2))
!!$                            ULambdaZ(L)=TwoZ*(UQz*ULambdaZ(L1)-RL2*ULambdaZ(L2))
!!$                         ENDDO
!!$                         DO L=0,Prim%Ell+1
!!$                            LambdaX(L)=LLambdaX(L)-ULambdaX(L)
!!$                            LambdaY(L)=LLambdaY(L)-ULambdaY(L)
!!$                            LambdaZ(L)=LLambdaZ(L)-ULambdaZ(L)
!!$                         ENDDO
!!$                         !
!!$                         DO L=0,Prim%Ell
!!$                            DO M=0,Prim%Ell-L
!!$                               DO N=0,Prim%Ell-M-L
!!$                                  LMN=LMNDex(L,M,N)
!!$                                  Pop=Pop+HGBra%D(LMN,IA,IB)*LambdaX(L)*LambdaY(M)*LambdaZ(N)
!!$                               ENDDO
!!$                            ENDDO
!!$                         ENDDO
!!$                         DBlk(IA,IB)=DBlk(IA,IB)+Pop
!!$                      ENDDO
!!$                   ENDDO
!!$                ENDIF
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    !
!!$    DVck = BlockToVect(NBFA,NBFB,DBlk)
!!$    !
!!$  END FUNCTION DipolCorrect
  !
END MODULE MBlock
