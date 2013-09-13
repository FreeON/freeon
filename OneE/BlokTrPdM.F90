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
MODULE BlokTrPdM
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
CONTAINS
  FUNCTION TrPdM(BS,Pair,P,COrig) RESULT(Vck)
    TYPE(BSET)                                 :: BS
    TYPE(AtomPair)                             :: Pair
    TYPE(DBL_VECT), INTENT(IN)                 :: COrig
    REAL(DOUBLE),DIMENSION(18)                 :: Vck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)    :: P
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,18) :: dDBlk
    REAL(DOUBLE),DIMENSION(3,-1:BS%NASym+1, &
                             -1:BS%NASym+1, &
                              0:1)             :: dEA,dEB
    REAL(DOUBLE),DIMENSION(3,-1:BS%NASym+2, &
                             -1:BS%NASym+2, &
                             -1:2*BS%NAsym+4)  :: E
!
    INTEGER      :: NBFA,NBFB,KA,KB,KK,AtA,AtB
    REAL(DOUBLE) :: Ax,Ay,Az,Bx,By,Bz,AB2
    REAL(DOUBLE) :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
                    EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB,DLA,DLB
    REAL(DOUBLE) :: POx,POy,POz
    INTEGER      :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                    StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                    LA,LB,MA,MB,NA,NB,K,NBF
    REAL(DOUBLE), EXTERNAL :: DDOT
!--------------------------------------------------------------------------------------------------
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
    dDBlk=Zero
!------------------------------------------------------------------------------------------
    DO CFA=1,BS%NCFnc%I(KA)
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          DO PFA=1,BS%NPFnc%I(CFA,KA)
          DO PFB=1,BS%NPFnc%I(CFB,KB)
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
             !vw--->
             POx = Px-COrig%D(1)
             POy = Py-COrig%D(2)
             POz = Pz-COrig%D(3)
             !vw<---
             E(1,0,0,1)=0.0D0
             E(2,0,0,1)=0.0D0
             E(3,0,0,1)=0.0D0
             CALL MD2TRR(BS%NASym+2,-1,MaxLA+2,MaxLB+2,EtaAB,E,PAx,PBx,PAy,PBy,PAz,PBz)
!------------------------------------------------------------------------------------
!            Naive derivatives of the McMurchie Davidson E and T coefficients WRT A
!            See Helgaker and Taylor, TCA v.83, p177 (1992)
             DO LA=0,MaxLA+1
             DLA=DBLE(LA)
             DO LB=0,MaxLB+1
                DLB=DBLE(LB)
                DO K=1,3
                   dEA(K,LA,LB,0)=Two*ZetaA*E(K,LA+1,LB,0)-DLA*E(K,LA-1,LB,0)
                   dEA(K,LA,LB,1)=Two*ZetaA*E(K,LA+1,LB,1)-DLA*E(K,LA-1,LB,1)
                   dEB(K,LA,LB,0)=Two*ZetaB*E(K,LA,LB+1,0)-DLB*E(K,LA,LB-1,0)
                   dEB(K,LA,LB,1)=Two*ZetaB*E(K,LA,LB+1,1)-DLB*E(K,LA,LB-1,1)
                ENDDO
             ENDDO
             ENDDO
!------------------------------------------------------------------------------------
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
                   ! Ax
                   !mu_x
                   dDBlk(IA,IB,1)=dDBlk(IA,IB,1) &
                              +CA*CB*Ov*((POx*dEA(1,LA,LB,0)+dEA(1,LA,LB,1)) &
                              *  E(2,MA,MB,0)*  E(3,NA,NB,0))
                   !mu_y
                   dDBlk(IA,IB,2)=dDBlk(IA,IB,2) &
                              +CA*CB*Ov*((POy*  E(2,MA,MB,0)+  E(2,MA,MB,1)) &
                              *dEA(1,LA,LB,0)*  E(3,NA,NB,0))
                   !mu_z
                   dDBlk(IA,IB,3)=dDBlk(IA,IB,3) &
                              +CA*CB*Ov*((POz*  E(3,NA,NB,0)+  E(3,NA,NB,1)) &
                              *dEA(1,LA,LB,0)*  E(2,MA,MB,0))
                   ! Ay
                   !mu_x
                   dDBlk(IA,IB,4)=dDBlk(IA,IB,4) &
                              +CA*CB*Ov*((POx*  E(1,LA,LB,0)+  E(1,LA,LB,1)) &
                              *dEA(2,MA,MB,0)*  E(3,NA,NB,0))
                   !mu_y
                   dDBlk(IA,IB,5)=dDBlk(IA,IB,5) &
                              +CA*CB*Ov*((POy*dEA(2,MA,MB,0)+dEA(2,MA,MB,1)) &
                              *  E(1,LA,LB,0)*  E(3,NA,NB,0))
                   !mu_z
                   dDBlk(IA,IB,6)=dDBlk(IA,IB,6) &
                              +CA*CB*Ov*((POz*  E(3,NA,NB,0)+  E(3,NA,NB,1)) &
                              *  E(1,LA,LB,0)*dEA(2,MA,MB,0))
                   ! Az
                   !mu_x
                   dDBlk(IA,IB,7)=dDBlk(IA,IB,7) &
                              +CA*CB*Ov*((POx*  E(1,LA,LB,0)+  E(1,LA,LB,1)) &
                              *  E(2,MA,MB,0)*dEA(3,NA,NB,0))
                   !mu_y
                   dDBlk(IA,IB,8)=dDBlk(IA,IB,8) &
                              +CA*CB*Ov*((POy*  E(2,MA,MB,0)+  E(2,MA,MB,1)) &
                              *  E(1,LA,LB,0)*dEA(3,NA,NB,0))
                   !mu_z
                   dDBlk(IA,IB,9)=dDBlk(IA,IB,9) &
                              +CA*CB*Ov*((POz*dEA(3,NA,NB,0)+dEA(3,NA,NB,1)) &
                              *  E(1,LA,LB,0)*  E(2,MA,MB,0))
                   ! Bx
                   !mu_x
                   dDBlk(IA,IB,10)=dDBlk(IA,IB,10) &
                              +CA*CB*Ov*((POx*dEB(1,LA,LB,0)+dEB(1,LA,LB,1)) &
                              *  E(2,MA,MB,0)*  E(3,NA,NB,0))
                   !mu_y
                   dDBlk(IA,IB,11)=dDBlk(IA,IB,11) &
                              +CA*CB*Ov*((POy*  E(2,MA,MB,0)+  E(2,MA,MB,1)) &
                              *dEB(1,LA,LB,0)*  E(3,NA,NB,0))
                   !mu_z
                   dDBlk(IA,IB,12)=dDBlk(IA,IB,12) &
                              +CA*CB*Ov*((POz*  E(3,NA,NB,0)+  E(3,NA,NB,1)) &
                              *dEB(1,LA,LB,0)*  E(2,MA,MB,0))
                   ! By
                   !mu_x
                   dDBlk(IA,IB,13)=dDBlk(IA,IB,13) &
                              +CA*CB*Ov*((POx*  E(1,LA,LB,0)+  E(1,LA,LB,1)) &
                              *dEB(2,MA,MB,0)*  E(3,NA,NB,0))
                   !mu_y
                   dDBlk(IA,IB,14)=dDBlk(IA,IB,14) &
                              +CA*CB*Ov*((POy*dEB(2,MA,MB,0)+dEB(2,MA,MB,1)) &
                              *  E(1,LA,LB,0)*  E(3,NA,NB,0))
                   !mu_z
                   dDBlk(IA,IB,15)=dDBlk(IA,IB,15) &
                              +CA*CB*Ov*((POz*  E(3,NA,NB,0)+  E(3,NA,NB,1)) &
                              *  E(1,LA,LB,0)*dEB(2,MA,MB,0))
                   ! Bz
                   !mu_x
                   dDBlk(IA,IB,16)=dDBlk(IA,IB,16) &
                              +CA*CB*Ov*((POx*  E(1,LA,LB,0)+  E(1,LA,LB,1)) &
                              *  E(2,MA,MB,0)*dEB(3,NA,NB,0))
                   !mu_y
                   dDBlk(IA,IB,17)=dDBlk(IA,IB,17) &
                              +CA*CB*Ov*((POy*  E(2,MA,MB,0)+  E(2,MA,MB,1)) &
                              *  E(1,LA,LB,0)*dEB(3,NA,NB,0))
                   !mu_z
                   dDBlk(IA,IB,18)=dDBlk(IA,IB,18) &
                              +CA*CB*Ov*((POz*dEB(3,NA,NB,0)+dEB(3,NA,NB,1)) &
                              *  E(1,LA,LB,0)*  E(2,MA,MB,0))
                ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    NBF=NBFA*NBFB
    DO K=1,18
       Vck(K)=DDOT(NBF,P(1,1),1,dDBlk(1,1,K),1)
    ENDDO
    !
  END FUNCTION TrPdM
END MODULE BlokTrPdM
