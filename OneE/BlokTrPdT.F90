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
!    BLOKWISE ACCUMULATION OF Tr{P_(A,B).dT^T_(A,B)/dA}
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE BlokTrPdT
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE

  CONTAINS
     FUNCTION TrPdT(BS,Pair,P) RESULT(Vck)
     TYPE(BSET)                                      :: BS
     TYPE(AtomPair)                                  :: Pair
     REAL(DOUBLE),DIMENSION(3)                       :: Vck
     REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)         :: P
     REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3)       :: dTBlk
     REAL(DOUBLE),DIMENSION(3,0:BS%NASym,0:BS%NASym) :: T,dT
     REAL(DOUBLE),DIMENSION(3,-1:BS%NASym+1, &
                              -1:BS%NASym+1)         :: dE
     REAL(DOUBLE),DIMENSION(3,-1:BS%NASym+2, &
                              -1:BS%NASym+2, &
                              -1:2*BS%NAsym+4)       :: E

!
     INTEGER                                   :: NBFA,NBFB,KA,KB,KK,AtA,AtB
     REAL(DOUBLE)                              :: Ax,Ay,Az,Bx,By,Bz,AB2
     REAL(DOUBLE)                              :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
                                                  EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB,DLA,DLB
     REAL(DOUBLE)                              :: Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz
     INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                                  StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                                  LA,LB,MA,MB,NA,NB,K
     REAL(DOUBLE), EXTERNAL                    :: BlkTrace_2
     LOGICAL                                   :: SameAtom
!--------------------------------------------------------------------------------------------------
     SameAtom=.FALSE.
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
     dTBlk=Zero
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
           CALL MD2TRR(BS%NASym+2,-1,MaxLA+2,MaxLB+2,EtaAB,E,PAx,PBx,PAy,PBy,PAz,PBz)
!------------------------------------------------------------------------------------
!          Naive derivatives of the McMurchie Davidson E and T coefficients WRT A
!          See Helgaker and Taylor, TCA v.83, p177 (1992)
           dE=Zero
           IF(.NOT. SameAtom)THEN
              DO LA=0,MaxLA+1
              DO LB=0,MaxLB+1
                 DO K=1,3
                    dE(K,LA,LB)=Two*ZetaA*E(K,LA+1,LB,0)-DBLE(LA)*E(K,LA-1,LB,0)
                 ENDDO
              ENDDO
              ENDDO
           ENDIF
!------------------------------------------------------------------------------------
!          Build McMurchie Davidson style T coefficients and their derivatives.
!          Note that dT_Ak/dAi=0 if k/=i
           DO LA=0,MaxLA
           DO LB=0,MaxLB
              DLA=DBLE(LA)
              DLB=DBLE(LB)
              DO K=1,3
                 T(K,LA,LB)=        DLA*DLB*E(K,LA-1,LB-1,0) &
                           +       Four*ZAB*E(K,LA+1,LB+1,0) &
                           -  ZetaA*Two*DLB*E(K,LA+1,LB-1,0) &
                           -  ZetaB*Two*DLA*E(K,LA-1,LB+1,0)

                 dT(K,LA,LB)=      DLA*DLB*dE(K,LA-1,LB-1) &
                            +     Four*ZAB*dE(K,LA+1,LB+1) &
                            -ZetaA*Two*DLB*dE(K,LA+1,LB-1) &
                            -ZetaB*Two*DLA*dE(K,LA-1,LB+1)
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
!                dT/dAx
                 Txx=dT(1,LA,LB)*   E(2,MA,MB,0)* E(3,NA,NB,0)
                 Tyx=dE(1,LA,LB)*   T(2,MA,MB)*   E(3,NA,NB,0)
                 Tzx=dE(1,LA,LB)*   E(2,MA,MB,0)* T(3,NA,NB)
                 dTBlk(IA,IB,1)=dTBlk(IA,IB,1)+Half*CA*CB*Ov*(Txx+Tyx+Tzx)
!                dT/dAy
                 Txy= T(1,LA,LB)*  dE(2,MA,MB)*   E(3,NA,NB,0)
                 Tyy= E(1,LA,LB,0)*dT(2,MA,MB)*   E(3,NA,NB,0)
                 Tzy= E(1,LA,LB,0)*dE(2,MA,MB)*   T(3,NA,NB)
                 dTBlk(IA,IB,2)=dTBlk(IA,IB,2)+Half*CA*CB*Ov*(Txy+Tyy+Tzy)
!                dT/dAz
                 Txz= T(1,LA,LB)*   E(2,MA,MB,0)*dE(3,NA,NB)
                 Tyz= E(1,LA,LB,0)* T(2,MA,MB)*  dE(3,NA,NB)
                 Tzz= E(1,LA,LB,0)* E(2,MA,MB,0)*dT(3,NA,NB)
                 dTBlk(IA,IB,3)=dTBlk(IA,IB,3)+Half*CA*CB*Ov*(Txz+Tyz+Tzz)
!
              ENDDO
           ENDDO
        ENDDO
        ENDDO
     ENDDO
     ENDDO
!
     DO K=1,3
        Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,P,TRANSPOSE(dTBlk(:,:,K)))
     ENDDO
!
  END FUNCTION TrPdT
END MODULE BlokTrPdT
