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
MODULE BlokTrPd2T
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
CONTAINS
  FUNCTION TrPd2T(BS,Pair,P) RESULT(Vck)
    TYPE(BSET)                                      :: BS
    TYPE(AtomPair)                                  :: Pair
    REAL(DOUBLE),DIMENSION(21)                      :: Vck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)         :: P
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,21)      :: d2TBlk
    REAL(DOUBLE),DIMENSION(3,0:BS%NASym,0:BS%NASym) :: T,dTA,dTB,dTAB,dTAA,dTBB
    REAL(DOUBLE),DIMENSION(3,-1:BS%NASym+1, &
                             -1:BS%NASym+1)         :: dEA,dEB,dEAB,dEAA,dEBB
    REAL(DOUBLE),DIMENSION(3,-2:BS%NASym+3, &
                             -2:BS%NASym+3, &
                             -2:2*BS%NAsym+6)       :: E
!
    INTEGER                                   :: NBFA,NBFB,KA,KB,KK,AtA,AtB
    REAL(DOUBLE)                              :: Ax,Ay,Az,Bx,By,Bz,AB2
    REAL(DOUBLE)                              :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,ZetaB,  &
                                                 EtaAB,EtaIn,ZAB,XiAB,ExpAB,Ov,PiE32,CA,CB,DLA,DLB
    REAL(DOUBLE)                              :: Dumx,Dumy,Dumz
    INTEGER                                   :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                                 StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                                 LA,LB,MA,MB,NA,NB,K
    REAL(DOUBLE), EXTERNAL                    :: BlkTrace_2
    LOGICAL                                   :: SameAtom
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
    d2TBlk=Zero
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
             CALL MD2TRR(BS%NASym+3,-2,MaxLA+3,MaxLB+3,EtaAB,E,PAx,PBx,PAy,PBy,PAz,PBz)
!------------------------------------------------------------------------------------
!          Naive derivatives of the McMurchie Davidson E and T coefficients WRT A
!          See Helgaker and Taylor, TCA v.83, p177 (1992)
             !
             !dE/dA and dE/dB
             dEA=0d0
             dEB=0d0
             !d^2E/dAdB, d^2E/dAdA and d^2E/dBdB.
             dEAB=0d0
             dEAA=0d0
             dEBB=0d0
             DO LA=0,MaxLA+1
             DO LB=0,MaxLB+1
                DLA=DBLE(LA)
                DLB=DBLE(LB)
                dEA(1,LA,LB)=Two*ZetaA*E(1,LA+1,LB,0)-DLA*E(1,LA-1,LB,0)
                dEA(2,LA,LB)=Two*ZetaA*E(2,LA+1,LB,0)-DLA*E(2,LA-1,LB,0)
                dEA(3,LA,LB)=Two*ZetaA*E(3,LA+1,LB,0)-DLA*E(3,LA-1,LB,0)
                !
                dEB(1,LA,LB)=Two*ZetaB*E(1,LA,LB+1,0)-DLB*E(1,LA,LB-1,0)
                dEB(2,LA,LB)=Two*ZetaB*E(2,LA,LB+1,0)-DLB*E(2,LA,LB-1,0)
                dEB(3,LA,LB)=Two*ZetaB*E(3,LA,LB+1,0)-DLB*E(3,LA,LB-1,0)
                !
                dEAB(1,LA,LB)= 4D0*ZetaA*ZetaB*E(1,LA+1,LB+1,0)-2D0*ZetaA*DLB*E(1,LA+1,LB-1,0) &
                              -2D0*ZetaB*DLA  *E(1,LA-1,LB+1,0)+DLA*DLB      *E(1,LA-1,LB-1,0)
                dEAB(2,LA,LB)= 4D0*ZetaA*ZetaB*E(2,LA+1,LB+1,0)-2D0*ZetaA*DLB*E(2,LA+1,LB-1,0) &
                              -2D0*ZetaB*DLA  *E(2,LA-1,LB+1,0)+DLA*DLB      *E(2,LA-1,LB-1,0)
                dEAB(3,LA,LB)= 4D0*ZetaA*ZetaB*E(3,LA+1,LB+1,0)-2D0*ZetaA*DLB*E(3,LA+1,LB-1,0) &
                              -2D0*ZetaB*DLA  *E(3,LA-1,LB+1,0)+DLA*DLB      *E(3,LA-1,LB-1,0)
                !
                dEAA(1,LA,LB)= 4D0*ZetaA**2 *E(1,LA+2,LB,0)-2D0*ZetaA*(2D0*DLA+1D0)*E(1,LA,LB,0) &
                             + DLA*(DLA-1D0)*E(1,LA-2,LB,0)
                dEAA(2,LA,LB)= 4D0*ZetaA**2 *E(2,LA+2,LB,0)-2D0*ZetaA*(2D0*DLA+1D0)*E(2,LA,LB,0) &
                             + DLA*(DLA-1D0)*E(2,LA-2,LB,0)
                dEAA(3,LA,LB)= 4D0*ZetaA**2 *E(3,LA+2,LB,0)-2D0*ZetaA*(2D0*DLA+1D0)*E(3,LA,LB,0) &
                             + DLA*(DLA-1D0)*E(3,LA-2,LB,0)
                !
                dEBB(1,LA,LB)= 4D0*ZetaB**2 *E(1,LA,LB+2,0)-2D0*ZetaB*(2D0*DLB+1D0)*E(1,LA,LB,0) &
                             + DLB*(DLB-1D0)*E(1,LA,LB-2,0)
                dEBB(2,LA,LB)= 4D0*ZetaB**2 *E(2,LA,LB+2,0)-2D0*ZetaB*(2D0*DLB+1D0)*E(2,LA,LB,0) &
                             + DLB*(DLB-1D0)*E(2,LA,LB-2,0)
                dEBB(3,LA,LB)= 4D0*ZetaB**2 *E(3,LA,LB+2,0)-2D0*ZetaB*(2D0*DLB+1D0)*E(3,LA,LB,0) &
                             + DLB*(DLB-1D0)*E(3,LA,LB-2,0)
             ENDDO
             ENDDO
!------------------------------------------------------------------------------------
!          Build McMurchie Davidson style T coefficients and their derivatives.
!          Note that dT_Ak/dAi=0 if k/=i
             DO LA=0,MaxLA
             DO LB=0,MaxLB
                DLA=DBLE(LA)
                DLB=DBLE(LB)
                T(1,LA,LB)=        DLA*DLB*E(1,LA-1,LB-1,0) &
                          +       Four*ZAB*E(1,LA+1,LB+1,0) &
                          -  ZetaA*Two*DLB*E(1,LA+1,LB-1,0) &
                          -  ZetaB*Two*DLA*E(1,LA-1,LB+1,0)
                T(2,LA,LB)=        DLA*DLB*E(2,LA-1,LB-1,0) &
                          +       Four*ZAB*E(2,LA+1,LB+1,0) &
                          -  ZetaA*Two*DLB*E(2,LA+1,LB-1,0) &
                          -  ZetaB*Two*DLA*E(2,LA-1,LB+1,0)
                T(3,LA,LB)=        DLA*DLB*E(3,LA-1,LB-1,0) &
                          +       Four*ZAB*E(3,LA+1,LB+1,0) &
                          -  ZetaA*Two*DLB*E(3,LA+1,LB-1,0) &
                          -  ZetaB*Two*DLA*E(3,LA-1,LB+1,0)
                !
                dTA(1,LA,LB)=      DLA*DLB*dEA(1,LA-1,LB-1) &
                          +       Four*ZAB*dEA(1,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEA(1,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEA(1,LA-1,LB+1)
                dTA(2,LA,LB)=      DLA*DLB*dEA(2,LA-1,LB-1) &
                          +       Four*ZAB*dEA(2,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEA(2,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEA(2,LA-1,LB+1)
                dTA(3,LA,LB)=      DLA*DLB*dEA(3,LA-1,LB-1) &
                          +       Four*ZAB*dEA(3,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEA(3,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEA(3,LA-1,LB+1)
                !
                dTB(1,LA,LB)=      DLA*DLB*dEB(1,LA-1,LB-1) &
                          +       Four*ZAB*dEB(1,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEB(1,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEB(1,LA-1,LB+1)
                dTB(2,LA,LB)=      DLA*DLB*dEB(2,LA-1,LB-1) &
                          +       Four*ZAB*dEB(2,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEB(2,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEB(2,LA-1,LB+1)
                dTB(3,LA,LB)=      DLA*DLB*dEB(3,LA-1,LB-1) &
                          +       Four*ZAB*dEB(3,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEB(3,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEB(3,LA-1,LB+1)
                !
                dTAB(1,LA,LB)=     DLA*DLB*dEAB(1,LA-1,LB-1) &
                          +       Four*ZAB*dEAB(1,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEAB(1,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEAB(1,LA-1,LB+1)
                dTAB(2,LA,LB)=     DLA*DLB*dEAB(2,LA-1,LB-1) &
                          +       Four*ZAB*dEAB(2,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEAB(2,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEAB(2,LA-1,LB+1)
                dTAB(3,LA,LB)=     DLA*DLB*dEAB(3,LA-1,LB-1) &
                          +       Four*ZAB*dEAB(3,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEAB(3,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEAB(3,LA-1,LB+1)
                !
                dTAA(1,LA,LB)=     DLA*DLB*dEAA(1,LA-1,LB-1) &
                          +       Four*ZAB*dEAA(1,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEAA(1,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEAA(1,LA-1,LB+1)
                dTAA(2,LA,LB)=     DLA*DLB*dEAA(2,LA-1,LB-1) &
                          +       Four*ZAB*dEAA(2,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEAA(2,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEAA(2,LA-1,LB+1)
                dTAA(3,LA,LB)=     DLA*DLB*dEAA(3,LA-1,LB-1) &
                          +       Four*ZAB*dEAA(3,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEAA(3,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEAA(3,LA-1,LB+1)
                !
                dTBB(1,LA,LB)=     DLA*DLB*dEBB(1,LA-1,LB-1) &
                          +       Four*ZAB*dEBB(1,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEBB(1,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEBB(1,LA-1,LB+1)
                dTBB(2,LA,LB)=     DLA*DLB*dEBB(2,LA-1,LB-1) &
                          +       Four*ZAB*dEBB(2,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEBB(2,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEBB(2,LA-1,LB+1)
                dTBB(3,LA,LB)=     DLA*DLB*dEBB(3,LA-1,LB-1) &
                          +       Four*ZAB*dEBB(3,LA+1,LB+1) &
                          -  ZetaA*Two*DLB*dEBB(3,LA+1,LB-1) &
                          -  ZetaB*Two*DLA*dEBB(3,LA-1,LB+1)
                !
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
                   !
                   !****************************************************************
                   ! A-B Terms
                   !****************************************************************
                   !d^2 T/dAxdBx
                   Dumx=dTAB(1,LA,LB  )*   E(2,MA,MB,0)*   E(3,NA,NB,0)
                   Dumy=dEAB(1,LA,LB  )*   T(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz=dEAB(1,LA,LB  )*   E(2,MA,MB,0)*   T(3,NA,NB  )
                   d2TBlk(IA,IB,1)=d2TBlk(IA,IB,1)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAxdBy
                   Dumx= dTA(1,LA,LB  )* dEB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy= dEA(1,LA,LB  )* dTB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz= dEA(1,LA,LB  )* dEB(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,2)=d2TBlk(IA,IB,2)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAxdBz
                   Dumx= dTA(1,LA,LB  )*   E(2,MA,MB,0)* dEB(3,NA,NB  )
                   Dumy= dEA(1,LA,LB  )*   T(2,MA,MB  )* dEB(3,NA,NB  )
                   Dumz= dEA(1,LA,LB  )*   E(2,MA,MB,0)* dTB(3,NA,NB  )
                   d2TBlk(IA,IB,3)=d2TBlk(IA,IB,3)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !d^2 T/dAydBx
                   Dumx= dTB(1,LA,LB  )* dEA(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy= dEB(1,LA,LB  )* dTA(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz= dEB(1,LA,LB  )* dEA(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,4)=d2TBlk(IA,IB,4)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAydBy
                   Dumx=   T(1,LA,LB  )*dEAB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy=   E(1,LA,LB,0)*dTAB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz=   E(1,LA,LB,0)*dEAB(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,5)=d2TBlk(IA,IB,5)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAydBz
                   Dumx=   T(1,LA,LB  )* dEA(2,MA,MB  )* dEB(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)* dTA(2,MA,MB  )* dEB(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)* dEA(2,MA,MB  )* dTB(3,NA,NB  )
                   d2TBlk(IA,IB,6)=d2TBlk(IA,IB,6)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !d^2 T/dAzdBx
                   Dumx= dTB(1,LA,LB  )*   E(2,MA,MB,0)* dEA(3,NA,NB  )
                   Dumy= dEB(1,LA,LB  )*   T(2,MA,MB  )* dEA(3,NA,NB  )
                   Dumz= dEB(1,LA,LB  )*   E(2,MA,MB,0)* dTA(3,NA,NB  )
                   d2TBlk(IA,IB,7)=d2TBlk(IA,IB,7)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAzdBy
                   Dumx=   T(1,LA,LB  )* dEB(2,MA,MB  )* dEA(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)* dTB(2,MA,MB  )* dEA(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)* dEB(2,MA,MB  )* dTA(3,NA,NB  )
                   d2TBlk(IA,IB,8)=d2TBlk(IA,IB,8)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAzdBz
                   Dumx=   T(1,LA,LB  )*   E(2,MA,MB,0)*dEAB(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)*   T(2,MA,MB  )*dEAB(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)*   E(2,MA,MB,0)*dTAB(3,NA,NB  )
                   d2TBlk(IA,IB,9)=d2TBlk(IA,IB,9)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !
                   !****************************************************************
                   ! A-A Terms
                   !****************************************************************
                   !d^2 T/dAxdAx
                   Dumx=dTAA(1,LA,LB  )*   E(2,MA,MB,0)*   E(3,NA,NB,0)
                   Dumy=dEAA(1,LA,LB  )*   T(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz=dEAA(1,LA,LB  )*   E(2,MA,MB,0)*   T(3,NA,NB  )
                   d2TBlk(IA,IB,10)=d2TBlk(IA,IB,10)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAxdAy
                   Dumx= dTA(1,LA,LB  )* dEA(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy= dEA(1,LA,LB  )* dTA(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz= dEA(1,LA,LB  )* dEA(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,11)=d2TBlk(IA,IB,11)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAxdAz
                   Dumx= dTA(1,LA,LB  )*   E(2,MA,MB,0)* dEA(3,NA,NB  )
                   Dumy= dEA(1,LA,LB  )*   T(2,MA,MB  )* dEA(3,NA,NB  )
                   Dumz= dEA(1,LA,LB  )*   E(2,MA,MB,0)* dTA(3,NA,NB  )
                   d2TBlk(IA,IB,12)=d2TBlk(IA,IB,12)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !d^2 T/dAydAy
                   Dumx=   T(1,LA,LB  )*dEAA(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy=   E(1,LA,LB,0)*dTAA(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz=   E(1,LA,LB,0)*dEAA(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,13)=d2TBlk(IA,IB,13)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dAydAz
                   Dumx=   T(1,LA,LB  )* dEA(2,MA,MB  )* dEA(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)* dTA(2,MA,MB  )* dEA(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)* dEA(2,MA,MB  )* dTA(3,NA,NB  )
                   d2TBlk(IA,IB,14)=d2TBlk(IA,IB,14)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !d^2 T/dAzdAz
                   Dumx=   T(1,LA,LB  )*   E(2,MA,MB,0)*dEAA(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)*   T(2,MA,MB  )*dEAA(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)*   E(2,MA,MB,0)*dTAA(3,NA,NB  )
                   d2TBlk(IA,IB,15)=d2TBlk(IA,IB,15)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !
                   !****************************************************************
                   ! B-B Terms
                   !****************************************************************
                   !d^2 T/dBxdBx
                   Dumx=dTBB(1,LA,LB  )*   E(2,MA,MB,0)*   E(3,NA,NB,0)
                   Dumy=dEBB(1,LA,LB  )*   T(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz=dEBB(1,LA,LB  )*   E(2,MA,MB,0)*   T(3,NA,NB  )
                   d2TBlk(IA,IB,16)=d2TBlk(IA,IB,16)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dBxdBy
                   Dumx= dTB(1,LA,LB  )* dEB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy= dEB(1,LA,LB  )* dTB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz= dEB(1,LA,LB  )* dEB(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,17)=d2TBlk(IA,IB,17)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dBxdBz
                   Dumx= dTB(1,LA,LB  )*   E(2,MA,MB,0)* dEB(3,NA,NB  )
                   Dumy= dEB(1,LA,LB  )*   T(2,MA,MB  )* dEB(3,NA,NB  )
                   Dumz= dEB(1,LA,LB  )*   E(2,MA,MB,0)* dTB(3,NA,NB  )
                   d2TBlk(IA,IB,18)=d2TBlk(IA,IB,18)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !d^2 T/dBydBy
                   Dumx=   T(1,LA,LB  )*dEBB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumy=   E(1,LA,LB,0)*dTBB(2,MA,MB  )*   E(3,NA,NB,0)
                   Dumz=   E(1,LA,LB,0)*dEBB(2,MA,MB  )*   T(3,NA,NB  )
                   d2TBlk(IA,IB,19)=d2TBlk(IA,IB,19)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !d^2 T/dBydBz
                   Dumx=   T(1,LA,LB  )* dEB(2,MA,MB  )* dEB(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)* dTB(2,MA,MB  )* dEB(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)* dEB(2,MA,MB  )* dTB(3,NA,NB  )
                   d2TBlk(IA,IB,20)=d2TBlk(IA,IB,20)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !d^2 T/dBzdBz
                   Dumx=   T(1,LA,LB  )*   E(2,MA,MB,0)*dEBB(3,NA,NB  )
                   Dumy=   E(1,LA,LB,0)*   T(2,MA,MB  )*dEBB(3,NA,NB  )
                   Dumz=   E(1,LA,LB,0)*   E(2,MA,MB,0)*dTBB(3,NA,NB  )
                   d2TBlk(IA,IB,21)=d2TBlk(IA,IB,21)+Half*CA*CB*Ov*(Dumx+Dumy+Dumz)
                   !****************************************************************
                   !
                ENDDO
             ENDDO
          ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DO K=1,21
       Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,P,TRANSPOSE(d2TBlk(:,:,K)))
    ENDDO
    !
  END FUNCTION TrPd2T

END MODULE BlokTrPd2T
