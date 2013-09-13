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
MODULE OverlapBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  IMPLICIT NONE
  CONTAINS
  FUNCTION SBlok(BS,Pair) RESULT(SVck)
    TYPE(BSET)                              :: BS
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim

    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: SVck

    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: SBlk
    INTEGER                                 :: NBFA,NBFB,KA,KB,KK
    REAL(DOUBLE)                            :: Ax,Ay,Az,Bx,By,Bz,AB2
    REAL(DOUBLE)                            :: PiE32,Amp
    INTEGER                                 :: CFA,CFB,PFA,PFB,IndexA,IndexB,StartLA,StartLB, &
                                               StopLA,StopLB,MaxLA,MaxLB,IA,IB,LMNA,LMNB,     &
                                               LA,LB,MA,MB,NA,NB
!--------------------------------------------------------------------------------------------------
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
    AB2=Pair%AB2
    SBlk=Zero
!----------------------------------
    Prim%A=Pair%A
    Prim%B=Pair%B
    Prim%AB2=Pair%AB2
    Prim%KA=Pair%KA
    Prim%KB=Pair%KB
!----------------------------------
    DO CFA=1,BS%NCFnc%I(KA)
    IndexA=CFBlokDex(BS,CFA,KA)
    StartLA=BS%LStrt%I(CFA,KA)
    StopLA =BS%LStop%I(CFA,KA)
    MaxLA=BS%ASymm%I(2,CFA,KA)
    DO CFB=1,BS%NCFnc%I(KB)
       IndexB  = CFBlokDex(BS,CFB,KB)
       StartLB = BS%LStrt%I(CFB,KB)
       StopLB  = BS%LStop%I(CFB,KB)
       MaxLB   = BS%ASymm%I(2,CFB,KB)
!----------------------------------
       Prim%CFA=CFA
       Prim%CFB=CFB
       Prim%Ell=MaxLA+MaxLB
!----------------------------------
       DO PFA=1,BS%NPFnc%I(CFA,KA)
       DO PFB=1,BS%NPFnc%I(CFB,KB)
!----------------------------------
          Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
          Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
          Prim%Zeta=Prim%ZA+Prim%ZB
          Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
          IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
             Prim%PFA=PFA
             Prim%PFB=PFB
             Amp=SetBraBlok(Prim,BS)
!----------------------------------
             PiE32=(Pi/Prim%Zeta)**(1.5D0)
             IB=IndexB
             DO LMNB=StartLB,StopLB
                IB=IB+1
                IA=IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   SBlk(IA,IB)=SBlk(IA,IB)+HGBra%D(1,IA,IB)*PiE32
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    SVck = BlockToVect(Pair%NA,Pair%NB,SBlk)
  END FUNCTION SBlok
END MODULE OverlapBlock
