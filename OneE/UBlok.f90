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
MODULE ECPBlock
  USE Parse
  USE Indexing
  USE AtomPairs
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  IMPLICIT NONE
  ! Global temporary arrays for explicit fortran intermediates
  REAL(DOUBLE),DIMENSION(5000) :: W,V
  ! Global solid angle.  Soliiidddd.
  REAL(DOUBLE),PARAMETER    :: FourPi=12.56637061435917295385D0
  ! Global array of binomial coefficients
  REAL(DOUBLE),DIMENSION(0:6,0:6) :: Binomial=RESHAPE((/     &
       1, 0, 0, 0, 0, 0, 0,    &
       1, 1, 0, 0, 0, 0, 0,    &
       1, 2, 1, 0, 0, 0, 0,    &
       1, 3, 3, 1, 0, 0, 0,    &
       1, 4, 6, 4, 1, 0, 0,    &
       1, 5,10,10, 5, 1, 0,    &
       1, 6,15,20,15, 6, 1 /),(/7,7/))
  !---------------------------------------------------------------------
  ! PARAMETERS FOR RADIAL INTEGRATION 
    INTEGER,PARAMETER                         :: NPts=32
!    INTEGER,PARAMETER                         :: NPts=64
!    INTEGER,PARAMETER                         :: NPts=128
    INTEGER,PARAMETER                         :: Infty=20
    INTEGER                                   :: IPts,IWts
    REAL(DOUBLE),PARAMETER                    :: Tau=1D-18
    REAL(DOUBLE)                              :: Discriminator
    REAL(DOUBLE),DIMENSION(1:NPts,1:Infty)    :: Points,Weights
    REAL(DOUBLE),DIMENSION(1:NPts)            :: XA,XB,EX
    INCLUDE "QQuad32.Inc"
!    INCLUDE "QQuad64.Inc"
!    INCLUDE "QQuad128.Inc"
CONTAINS  !
  FUNCTION UBlock(BS,Pair,KC,Cx,Cy,Cz) RESULT(UVck)
    TYPE(BSET)                              :: BS
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: UVck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: ECP
    INTEGER :: NBFA,NBFB,KA,KB,KK,KC,CFA,CFB,PFA, &
         PFB,IndexA,IndexB,StartLA,StartLB, &
         StopLA,StopLB,MaxLA,MaxLB,IA,IB,   &
         LMNA,LMNB,LA,LB,MA,MB,NA,NB,AtC,   &
         PFC,ECPGL,ProjL,EllA,EllB,Lambda,LenAB
    REAL(DOUBLE) :: Ax,Ay,Az,Bx,By,Bz,AB2,Cx,Cy,Cz,   &
         ACx,ACy,ACz,BCx,BCy,BCz,AC2,BC2,AC,BC,&
         ZetaA,ZetaB,ZetaAB,ZetaC,Alpha,       &
         ZzAC2,ZzBC2,Px,Py,Pz,Gauss,CA,CB,CC,  & 
         Kx,Ky,Kz,KAx,KAy,KAz,KBx,KBy,KBz,     &
         Kappa,KappaA,KappaB, tmp
    ! Static work arrays; dimensions defined in Modules/GlobalScalars
    REAL(DOUBLE),DIMENSION(0:HGEll,1:HGLen) :: Omega1
    REAL(DOUBLE),DIMENSION(0:HGEll,0:HGEll+ECPEll) :: Keew1
    REAL(DOUBLE),DIMENSION(0:2*PrjEll,0:BFEll+PrjEll,1:HGLen) :: OmegaA,OmegaB
    REAL(DOUBLE),DIMENSION(0:PrjEll+BFEll,0:PrjEll+BFEll,0:HGEll+ECPEll) :: Keew2
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
    ACx=Ax-Cx
    ACy=Ay-Cy
    ACz=Az-Cz
    BCx=Bx-Cx
    BCy=By-Cy
    BCz=Bz-Cz
    AB2=Pair%AB2    
    AC2=(ACx**2+ACy**2+ACz**2)
    BC2=(BCx**2+BCy**2+BCz**2)
    ! Some testing of overlaps
    ECP=Zero
    IF(AC2+BC2<AtomPairDistanceThreshold)THEN
       AC=SQRT(AC2)
       BC=SQRT(BC2)
       ! Go over basis functions for this atom-atom blok
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
             LenAB   = LHGTF(MaxLA+MaxLB)
             DO PFA=1,BS%NPFnc%I(CFA,KA) 
                DO PFB=1,BS%NPFnc%I(CFB,KB)        
                   ZetaA=BS%Expnt%D(PFA,CFA,KA)
                   ZetaB=BS%Expnt%D(PFB,CFB,KB)
                   ZzAC2=ZetaA*AC2
                   ZzBC2=ZetaB*BC2
                   IF(ZzAC2+ZzBC2<2.D2)THEN !PrimPairDistanceThreshold)THEN
                      !                         IF(.TRUE.)THEN
                      ZetaAB=ZetaA+ZetaB
                      Px=(zetaA*Ax+zetaB*Bx)/ZetaAB
                      Py=(zetaA*Ay+zetaB*By)/ZetaAB
                      Pz=(zetaA*Az+zetaB*Bz)/ZetaAB
                      Kx=Px-Cx
                      Ky=Py-Cy
                      Kz=Pz-Cz
                      Kappa=Two*ZetaAB*SQRT(Kx*Kx+Ky*Ky+Kz*Kz)
                      Gauss=EXP(-ZzAC2-ZzBC2)
                      ! Type one angular integrals, including 4 Pi and with mu terms summed
                      CALL AngularOne(MaxLA+MaxLB,Kx,Ky,Kz,Omega1) 
                      ! Go over primitives in the unprojected, type one integrals
                      DO PFC=1,BS%NTyp1PF%I(KC)
                         ! Here is the type 1 primitive contraction coefficient,
                         CC=BS%Typ1CCo%D(PFC,KC)
                         ! the exponent, 
                         ZetaC=BS%Typ1Exp%D(PFC,KC)
                         ! and the radial "ell" value of the primitive Gaussian ECP
                         ECPGL=BS%Typ1Ell%I(PFC,KC)
                         ! Compute type one radial integrals
                         Alpha=ZetaA+ZetaB+ZetaC
                         CALL RadialOne(MaxLA+MaxLB,MaxLA+MaxLB+ECPGL,Alpha,Kappa,ZzAC2+ZzBC2,Keew1)
                         ! Go over basis function symmetries
                         IA=IndexA
                         DO LMNA=StartLA,StopLA
                            LA=BS%LxDex%I(LMNA)
                            MA=BS%LyDex%I(LMNA) 
                            NA=BS%LzDex%I(LMNA)
                            EllA=LA+MA+NA
                            CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                            IA=IA+1
                            IB=IndexB
                            DO LMNB=StartLB,StopLB
                               LB=BS%LxDex%I(LMNB)
                               MB=BS%LyDex%I(LMNB) 
                               NB=BS%LzDex%I(LMNB)
                               EllB=LB+MB+NB
                               CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                               IB=IB+1
                               ECP(IA,IB)=ECP(IA,IB) &
                                    +CA*CB*CC*ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL,         &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1) 
                            ENDDO
                         ENDDO
                      ENDDO
                      ! Loop on ECP projector symmetries 
                      DO ProjL=0,BS%ProjEll%I(KC)
                         ! Go over primitives in the projected, type two integrals
                         DO PFC=1,BS%NTyp2PF%I(ProjL,KC)                     
                            ! Here is the type 2 primitive contraction coefficient,
                            CC=BS%Typ2CCo%D(PFC,ProjL,KC)
                            ! the corresponding exponent, 
                            ZetaC=BS%Typ2Exp%D(PFC,ProjL,KC)     
                            ! and the radial "ell" value of the ECP primitive 
                            ECPGL=BS%Typ2Ell%I(PFC,ProjL,KC)
                            ! Compute the type two radial integrals                            
                            Alpha=ZetaA+ZetaB+ZetaC
                            ! Gaussian prefactor check
                            KappaA=Two*ZetaA*AC
                            KappaB=Two*ZetaB*BC
                            CALL RadialTwo(CC,ProjL+MaxLA,ProjL+MaxLB,MaxLA+MaxLB+ECPGL, &
                                 Alpha,ZzAC2,KappaA,ZzBC2,KappaB,Keew2) 
                            ! Compute the type two angular integrals
                            CALL AngularTwo(ProjL,MaxLA,ACx,ACy,ACz,OmegaA)
                            CALL AngularTwo(ProjL,MaxLB,BCx,BCy,BCz,OmegaB)
                            ! Go over basis function symmetries
                            IA=IndexA
                            DO LMNA=StartLA,StopLA
                               LA=BS%LxDex%I(LMNA)
                               MA=BS%LyDex%I(LMNA) 
                               NA=BS%LzDex%I(LMNA)
                               EllA=LA+MA+NA
                               CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                               IA=IA+1
                               IB=IndexB
                               DO LMNB=StartLB,StopLB
                                  LB=BS%LxDex%I(LMNB)
                                  MB=BS%LyDex%I(LMNB) 
                                  NB=BS%LzDex%I(LMNB)
                                  EllB=LB+MB+NB
                                  CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                  IB=IB+1
                                  ECP(IA,IB)=ECP(IA,IB)+CA*CB*CC &
                                       *ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL,   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2) 
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    UVck=BlockToVect(Pair%NA,Pair%NB,ECP)
  END FUNCTION UBlock

  FUNCTION TrPdU(BS,Pair,KC,Cx,Cy,Cz,P)  RESULT(Vck)
    TYPE(BSET)                                :: BS
    TYPE(CRDS)                                :: GM
    TYPE(AtomPair)                            :: Pair
    TYPE(PrimPair)                            :: Prim
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: P
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3) :: dUBlk
    REAL(DOUBLE),DIMENSION(3)                 :: Vck
    INTEGER :: NBFA,NBFB,KA,KB,KK,KC,CFA,CFB,PFA, &
         PFB,IndexA,IndexB,StartLA,StartLB, &
         StopLA,StopLB,MaxLA,MaxLB,IA,IB,   &
         LMNA,LMNB,LA,LB,MA,MB,NA,NB,AtC,   &
         PFC,ECPGL,ProjL,EllA,EllB,Lambda,LenAB,K
    REAL(DOUBLE) :: Ax,Ay,Az,Bx,By,Bz,AB2,Cx,Cy,Cz,   &
         ACx,ACy,ACz,BCx,BCy,BCz,AC2,BC2,AC,BC,&
         ZetaA,ZetaB,ZetaAB,ZetaC,Alpha,       &
         ZzAC2,ZzBC2,Px,Py,Pz,Gauss,CA,CB,CC,  & 
         Kx,Ky,Kz,KAx,KAy,KAz,KBx,KBy,KBz,     &
         Kappa,KappaA,KappaB, tmp
    ! Static work arrays; dimensions defined in Modules/GlobalScalars
    REAL(DOUBLE),DIMENSION(0:HGEll,1:HGLen) :: Omega1
    REAL(DOUBLE),DIMENSION(0:HGEll,0:HGEll+ECPEll) :: Keew1
    REAL(DOUBLE),DIMENSION(0:2*PrjEll,0:BFEll+PrjEll,1:HGLen) :: OmegaA,OmegaB
    REAL(DOUBLE),DIMENSION(0:PrjEll+BFEll,0:PrjEll+BFEll,0:HGEll+ECPEll) :: Keew2
    REAL(DOUBLE), EXTERNAL  :: BlkTrace_2
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
    ACx=Ax-Cx
    ACy=Ay-Cy
    ACz=Az-Cz
    BCx=Bx-Cx
    BCy=By-Cy
    BCz=Bz-Cz
    AB2=Pair%AB2
    AC2=(ACx**2+ACy**2+ACz**2)
    BC2=(BCx**2+BCy**2+BCz**2)
    dUBlk=Zero
    ! Some testing of overlaps
    IF(AC2+BC2<AtomPairDistanceThreshold)THEN
       AC=SQRT(AC2)
       BC=SQRT(BC2)
       ! Go over basis functions for this atom-atom blok
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
             LenAB   = LHGTF(MaxLA+MaxLB)
             DO PFA=1,BS%NPFnc%I(CFA,KA) 
                DO PFB=1,BS%NPFnc%I(CFB,KB)        
                   ZetaA=BS%Expnt%D(PFA,CFA,KA)
                   ZetaB=BS%Expnt%D(PFB,CFB,KB)
                   ZzAC2=ZetaA*AC2
                   ZzBC2=ZetaB*BC2
                   IF(ZzAC2+ZzBC2<2.D2)THEN 
                      ZetaAB=ZetaA+ZetaB
                      Px=(zetaA*Ax+zetaB*Bx)/ZetaAB
                      Py=(zetaA*Ay+zetaB*By)/ZetaAB
                      Pz=(zetaA*Az+zetaB*Bz)/ZetaAB
                      Kx=Px-Cx
                      Ky=Py-Cy
                      Kz=Pz-Cz
                      Kappa=Two*ZetaAB*SQRT(Kx*Kx+Ky*Ky+Kz*Kz)
                      Gauss=EXP(-ZzAC2-ZzBC2)
                      ! Type one angular integrals, including 4 Pi and with mu terms summed
                      CALL AngularOne(MaxLA+1+MaxLB,Kx,Ky,Kz,Omega1) 
                      ! Go over primitives in the unprojected, type one integrals
                      DO PFC=1,BS%NTyp1PF%I(KC)
                         ! Here is the type 1 primitive contraction coefficient,
                         CC=BS%Typ1CCo%D(PFC,KC)
                         ! the exponent, 
                         ZetaC=BS%Typ1Exp%D(PFC,KC)
                         ! and the radial "ell" value of the primitive Gaussian ECP
                         ECPGL=BS%Typ1Ell%I(PFC,KC)
                         ! Compute type one radial integrals
                         Alpha=ZetaA+ZetaB+ZetaC
                         CALL RadialOne(MaxLA+1+MaxLB,MaxLA+1+MaxLB+ECPGL,Alpha,Kappa,ZzAC2+ZzBC2,Keew1)
                         ! Go over basis function symmetries
                         IA=IndexA
                         DO LMNA=StartLA,StopLA
                            LA=BS%LxDex%I(LMNA)
                            MA=BS%LyDex%I(LMNA) 
                            NA=BS%LzDex%I(LMNA)
                            EllA=LA+MA+NA
                            CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                            IA=IA+1
                            IB=IndexB
                            DO LMNB=StartLB,StopLB
                               LB=BS%LxDex%I(LMNB)
                               MB=BS%LyDex%I(LMNB) 
                               NB=BS%LzDex%I(LMNB)
                               EllB=LB+MB+NB
                               CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                               IB=IB+1
                               dUBlk(IA,IB,1)=dUBlk(IA,IB,1)+CA*CB*CC*(                             &
                                    Two*ZetaA*ECPOneSum(LA+1,MA,NA,LB,MB,NB,ECPGL,                  &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1) &
                                    -DBLE(LA)*ECPOneSum(LA-1,MA,NA,LB,MB,NB,ECPGL,                  &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1))
                               dUBlk(IA,IB,2)=dUBlk(IA,IB,2)+CA*CB*CC*(                             &
                                    Two*ZetaA*ECPOneSum(LA,MA+1,NA,LB,MB,NB,ECPGL,                  &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1) &
                                    -DBLE(MA)*ECPOneSum(LA,MA-1,NA,LB,MB,NB,ECPGL,                  &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1))
                               dUBlk(IA,IB,3)=dUBlk(IA,IB,3)+CA*CB*CC*(                             &
                                    Two*ZetaA*ECPOneSum(LA,MA,NA+1,LB,MB,NB,ECPGL,                  &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1) &
                                    -DBLE(NA)*ECPOneSum(LA,MA,NA-1,LB,MB,NB,ECPGL,                  &
                                    -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,Omega1,Keew1))
                            ENDDO
                         ENDDO
                      ENDDO

!!$                      ! Loop on ECP projector symmetries 
                      DO ProjL=0,BS%ProjEll%I(KC)
                         ! Go over primitives in the projected, type two integrals
                         DO PFC=1,BS%NTyp2PF%I(ProjL,KC)                     
                            ! Here is the type 2 primitive contraction coefficient,
                            CC=BS%Typ2CCo%D(PFC,ProjL,KC)
                            ! the corresponding exponent, 
                            ZetaC=BS%Typ2Exp%D(PFC,ProjL,KC)     
                            ! and the radial "ell" value of the ECP primitive 
                            ECPGL=BS%Typ2Ell%I(PFC,ProjL,KC)
                            ! Compute the type two radial integrals                            
                            Alpha=ZetaA+ZetaB+ZetaC
                            ! Gaussian prefactor check
                            KappaA=Two*ZetaA*AC
                            KappaB=Two*ZetaB*BC
                            CALL RadialTwo(CC,ProjL+MaxLA+1,ProjL+MaxLB,MaxLA+1+MaxLB+ECPGL, &
                                 Alpha,ZzAC2,KappaA,ZzBC2,KappaB,Keew2) 
                            ! Compute the type two angular integrals
                            CALL AngularTwo(ProjL,MaxLA+1,ACx,ACy,ACz,OmegaA)
                            CALL AngularTwo(ProjL,MaxLB  ,BCx,BCy,BCz,OmegaB)
                            ! Go over basis function symmetries
                            IA=IndexA
                            DO LMNA=StartLA,StopLA
                               LA=BS%LxDex%I(LMNA)
                               MA=BS%LyDex%I(LMNA) 
                               NA=BS%LzDex%I(LMNA)
                               EllA=LA+MA+NA
                               CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                               IA=IA+1
                               IB=IndexB
                               DO LMNB=StartLB,StopLB
                                  LB=BS%LxDex%I(LMNB)
                                  MB=BS%LyDex%I(LMNB) 
                                  NB=BS%LzDex%I(LMNB)
                                  EllB=LB+MB+NB
                                  CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                                  IB=IB+1
                                  dUBlk(IA,IB,1)=dUBlk(IA,IB,1)+CA*CB*CC*(                                    &
                                       Two*ZetaA*ECPTwoSum(LA+1,MA,NA,LB,MB,NB,ProjL,ECPGL,                   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2) &
                                       -DBLE(LA)*ECPTwoSum(LA-1,MA,NA,LB,MB,NB,ProjL,ECPGL,                   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2))
                                  dUBlk(IA,IB,2)=dUBlk(IA,IB,2)+CA*CB*CC*(                                    &
                                       Two*ZetaA*ECPTwoSum(LA,MA+1,NA,LB,MB,NB,ProjL,ECPGL,                   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2) &
                                       -DBLE(MA)*ECPTwoSum(LA,MA-1,NA,LB,MB,NB,ProjL,ECPGL,                   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2))
                                  dUBlk(IA,IB,3)=dUBlk(IA,IB,3)+CA*CB*CC*(                                    &
                                       Two*ZetaA*ECPTwoSum(LA,MA,NA+1,LB,MB,NB,ProjL,ECPGL,                   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2) &
                                       -DBLE(NA)*ECPTwoSum(LA,MA,NA-1,LB,MB,NB,ProjL,ECPGL,                   &
                                       -ACx,-ACy,-ACz,-BCx,-BCy,-BCz,OmegaA,OmegaB,Keew2))
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    DO K=1,3
       Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,P,TRANSPOSE(dUBlk(:,:,K)))
    ENDDO
  END FUNCTION TrPdU
  
  FUNCTION ECPOneSum(LA,MA,NA,LB,MB,NB,ECPGL,ACx,ACy,ACz,BCx,BCy,BCz,O,Q) RESULT(V1)
    INTEGER :: LA,MA,NA,LB,MB,NB,ECPGL,EllAB,  &
         ProjL,Emm,Lambda,LambdaA,LambdaB,AlphabetL, &
         a,b,c,d,e,f,ABCDEF, lll,vvv
    REAL(DOUBLE) :: ACx,ACy,ACz,BCx,BCy,BCz,V1, &
         BinA,BinB,BinC,BinD,BinE,BinF,Prim
    REAL(DOUBLE),DIMENSION(0:,1:) :: O
    REAL(DOUBLE),DIMENSION(0:,0:) :: Q

!!$    WRITE(*,*)' ======================================================================'
!!$    DO LLL=0,LA+MA+NA+MB+NB
!!$       DO vvv=0,LA+MA+NA+MB+NB+ECPGL       
!!$          WRITE(*,*)LLL,VVV,' Q2 = ',Q(LLL,VVV)
!!$       ENDDO
!!$    ENDDO

    ACx=ACx+1D-100
    ACy=ACy+1D-100
    ACz=ACz+1D-100
    BCx=BCx+1D-100
    BCy=BCy+1D-100
    BCz=BCz+1D-100
    !
    V1=Zero
    EllAB=LA+MA+NA+LB+MB+NB
    DO a=0,LA
       BinA=ACx**(LA-a)*Binomial(a,LA)
       DO b=0,MA
          BinB=BinA*ACy**(MA-b)*Binomial(b,MA)
          DO c=0,NA
             BinC=BinB*ACz**(NA-c)*Binomial(c,NA)
             DO d=0,LB
                BinD=BinC*BCx**(LB-d)*Binomial(d,LB)
                DO e=0,MB
                   BinE=BinD*BCy**(MB-e)*Binomial(e,MB)
                   DO f=0,NB
                      BinF=BinE*BCz**(NB-f)*Binomial(f,NB)
                      ! All the funny indexing ...
                      AlphabetL=a+b+c+d+e+f+ECPGL
                      ABCDEF=LMNDex(a+d,b+e,c+f)
                      DO Lambda=0,EllAB
                         V1=V1+BinF*Q(Lambda,AlphabetL)*O(Lambda,ABCDEF)
!!$                            Prim=BinF*Q(Lambda,AlphabetL)*O(Lambda,ABCDEF)
!!$                            IF(LA+MA+NA==0.AND.NB==1.AND.ABS(Prim)>1D-10)THEN
!!$                               WRITE(*,*)' abcdef = ',a,b,c,d,e,f   
!!$                               WRITE(*,*)Lambda,AlphabetL,ABCDEF
!!$                               WRITE(*,*)' BinF = ',BinF 
!!$                               WRITE(*,*)' O = ',O(Lambda,ABCDEF)
!!$                               WRITE(*,*)' Q = ',Q(Lambda,AlphabetL)
!!$                            ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO 
          ENDDO
       ENDDO
    ENDDO
    V1=V1*FourPi
  END FUNCTION ECPOneSum

  FUNCTION ECPTwoSum(LA,MA,NA,LB,MB,NB,ProjL,ECPGL,ACx,ACy,ACz,BCx,BCy,BCz,OA,OB,Q) RESULT(V2)
    INTEGER :: LA,MA,NA,LB,MB,NB,ECPGL,EllAB,  &
         ProjL,Emm,LambdaA,LambdaB,AlphabetL, &
         a,b,c,d,e,f,LMNabc,LMNdef
    REAL(DOUBLE) :: ACx,ACy,ACz,BCx,BCy,BCz,V2, &
         BinA,BinB,BinC,BinD,BinE,BinF,         CC

    REAL(DOUBLE),DIMENSION(0:,0:,1:) :: OA,OB
    REAL(DOUBLE),DIMENSION(0:,0:,0:) :: Q


!WRITE(*,*)'============================================================================='
!WRITE(*,*)' ProjL = ',ProjL
!
    ACx=ACx+1D-100
    ACy=ACy+1D-100
    ACz=ACz+1D-100
    BCx=BCx+1D-100
    BCy=BCy+1D-100
    BCz=BCz+1D-100
    !
    V2=Zero
    EllAB=LA+MA+NA+LB+MB+NB
    DO a=0,LA
       BinA=ACx**(LA-a)*Binomial(a,LA)
       DO b=0,MA
          BinB=BinA*ACy**(MA-b)*Binomial(b,MA)
          DO c=0,NA
             BinC=BinB*ACz**(NA-c)*Binomial(c,NA)
             LMNabc=LMNDex(a,b,c)
             DO d=0,LB
                BinD=BinC*BCx**(LB-d)*Binomial(d,LB)
                DO e=0,MB
                   BinE=BinD*BCy**(MB-e)*Binomial(e,MB)
                   DO f=0,NB
                      LMNdef=LMNDex(d,e,f)
                      AlphabetL=a+b+c+d+e+f+ECPGL
                      BinF=BinE*BCz**(NB-f)*Binomial(f,NB)
                      ! Contract the type two (projected) contribution
                      DO LambdaA=0,ProjL+a+b+c
                         DO LambdaB=0,ProjL+d+e+f
                            DO Emm=0,2*ProjL
                               V2=V2+BinF*OA(Emm,LambdaA,LMNabc)*OB(Emm,LambdaB,LMNdef)*Q(LambdaA,LambdaB,AlphabetL) 

#ifdef MAXDEBUG2
!IF(NA==2.AND.NB==0)THEN
!ENDIF
!IF(NA==2.AND.NB==2.AND.c==2.AND.f==2)THEN
   WRITE(*,30)'A ',Emm,LambdaA,LMNabc,OA(Emm,LambdaA,LMNabc),c,ACz**(NA-c),Binomial(c,NA)
   WRITE(*,30)'B ',Emm,LambdaB,LMNdef,OB(Emm,LambdaB,LMNdef),f,BCz**(NB-f),Binomial(f,NB)
30 format(A2,' m = ',I2,' L = ',I2,' LMN = ',I2,' O = ',F12.6,' c/f = ',I2,' A/BCz = ',F12.6,' Bin = ',F12.6)
  IF(ABS(BinF*OA(Emm,LambdaA,LMNabc)*OB(Emm,LambdaB,LMNdef)*Q(LambdaA,LambdaB,AlphabetL))>1D-5)THEN
   WRITE(*,31)'NA',NA,C,ACz,Binomial(C,NA)
   WRITE(*,31)'NB',NB,F,BCz,Binomial(f,NB)
31 format(A2,I2,I2,' ZZ = ',F20.10,' Bin = ',F20.10)
   WRITE(*,32)Emm,LambdaA,LambdaB,AlphabetL,BinF*OA(Emm,LambdaA,LMNabc)*OB(Emm,LambdaB,LMNdef),Q(LambdaA,LambdaB,AlphabetL),cc*V2*FourPi
32 format(' Emm = ',I2,' LA = ',I2,' LB = ',I2,' N = ',I2,' AA = ',F14.8,' QQ = ',F14.8,', CC*V2 = ',F14.8)
ENDIF
!ENDIF
#endif

                            ENDDO
                         ENDDO
                      ENDDO
                      !<-
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    V2=V2*FourPi
  END FUNCTION ECPTwoSum


  SUBROUTINE AngularOne(Ell,X,Y,Z,O) 
    INTEGER                       :: Ell,  l,lmn
    REAL(DOUBLE)                  :: K,X,Y,Z,Kx,Ky,Kz
    REAL(DOUBLE),DIMENSION(0:,1:) :: O
    K=SQRT(X**2+Y**2+Z**2)
    IF(K/=Zero)THEN
       Kx=X/K
       Ky=Y/K
       Kz=Z/K
    ELSE
       Kx=Zero
       Ky=Zero
       Kz=One
    ENDIF
    SELECT CASE(Ell)
    INCLUDE "Omega1.Inc"
    CASE DEFAULT
       CALL Halt(' No explicit code for case Ell = '//TRIM(IntToChar(Ell))//' in AngularOne.')
    END SELECT
!    DO L=0,Ell
!       DO LMN=1,LHGTF(Ell)
!          WRITE(*,*)L,LMN,' Omega = ',O(L,LMN)
!       ENDDO
!    ENDDO
  END SUBROUTINE AngularOne

  SUBROUTINE AngularTwo(ProjL,Ell,X,Y,Z,O)
    INTEGER                               :: ProjL,Ell,Index
    REAL(DOUBLE)                          :: X,Y,Z,Kx,Ky,Kz,K
    REAL(DOUBLE), DIMENSION(0:,0:,1:)     :: O

    integer :: L,P,LMN,M,N,Lambda
    K=SQRT(X**2+Y**2+Z**2)
    IF(K/=Zero)THEN
       Kx=X/K
       Ky=Y/K
       Kz=Z/K
    ELSE
       Kx=Zero
       Ky=Zero
       Kz=One
    ENDIF
    Index=100*Ell+ProjL
    SELECT CASE(Index)
    INCLUDE "Omega2.Inc"
    CASE DEFAULT
       CALL Halt(' No explicit code for case Ell = '//TRIM(IntToChar(Ell))// &
                 ', ProjL = '//TRIM(IntToChar(ProjL))//' in AngularTwo.')
    END SELECT
!    WRITE(*,*)'------------------------------------------------------'
!    DO Lambda=0,Ell
!       DO P=-ProjL,ProjL
!          DO L=0,Ell
!             DO M=0,Ell-L
!                DO N=0,Ell-L-M
!                   LMN=LMNDex(L,M,N)
!                   WRITE(*,22)lambda,ProjL,P,L,M,N,Kx,Ky,Kz
!22                 format(' Print[omega[',I2,",",I2,",",I3,",",I3,",",I3,",",I3, &
!                          ']/.{Kx->',F20.10,',Ky->',F20.10,',Kz->',F20.10,'}]')
!                   WRITE(*,*)' Print[" omega = "',O(P+ProjL,L,LMN),'];'
!                ENDDO
!             ENDDO
!           ENDDO
!        ENDDO
!     ENDDO
  END SUBROUTINE AngularTwo

  SUBROUTINE RadialOne(Ell,N,Alpha,K,ZAC2pZBC2,Q)
    INTEGER                                :: Ell,Ehn,N,I,Lambda,IGrid
    REAL(DOUBLE)                           :: K,Alpha,ZAC2pZBC2,R,Delta,Left,Rhgt,Cntr,DistToZero,AlphaInv,Ov
    REAL(DOUBLE),DIMENSION(1:NPts,0:HGEll) :: M
    REAL(DOUBLE),DIMENSION(0:,0:)          :: Q
    !
    IF(K<1D-8)THEN
       Q=Zero
       AlphaInv=One/Alpha
       Ov=EXP(-ZAC2pZBC2)
       ! Only values with Lambda=0 survive
       Q(0,0)=Half*Ov*SQRT(Pi*AlphaInv)
       Q(0,1)=Half*Ov*AlphaInv
       DO Ehn=1,N-1
          Q(0,Ehn+1)=Half*AlphaInv*DBLE(Ehn)*Q(0,Ehn-1)
       ENDDO
    ELSE
       Discriminator=4D0*Alpha*(ZAC2pZBC2+LOG(Tau))
       IF(K**2<Discriminator)THEN
          Q=Zero
       ELSE
          Rhgt=(K+SQRT(K**2-Discriminator))/(Two*Alpha)
          IGrid=CEILING(Rhgt)
          IF(IGrid>Infty)THEN
             WRITE(*,*)' Rhgt = ',Rhgt
             STOP
          ENDIF
          DO I=1,NPts
             XA(I)=K*Points(I,IGrid)
             EX(I)=EXP(-Alpha*Points(I,IGrid)**2)
          ENDDO
          CALL BesselM(Ell,NPts,ZAC2pZBC2,XA,M)
          DO Ehn=0,N
             DO Lambda=0,Ell
                Q(Lambda,Ehn)=Zero
                DO I=NPts,1,-1
                   Q(Lambda,Ehn)=Q(Lambda,Ehn)+Weights(I,IGrid)*M(I,Lambda)*EX(I)*Points(I,IGrid)**(Ehn)
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF
!    write(*,*)'=============================================================='
!    write(*,*)' Ell = ',Ell,' N = ',N
!    DO Ehn=0,N
!       DO Lambda=0,Ell
!          WRITE(*,*)Lambda,Ehn,Q(Lambda,Ehn)
!       ENDDO
!    ENDDO
  END SUBROUTINE RadialOne

  SUBROUTINE RadialTwo(CC,EllA,EllB,N,Alpha,ZAC2,Ka,ZBC2,Kb,Q)
    REAL(DOUBLE),DIMENSION(1:NPts,0:BFEll+PrjEll) :: MA,MB
    INTEGER                                       :: EllA,EllB,N,LambdaA,LambdaB,I,Ehn,IGrid
    REAL(DOUBLE),DIMENSION(0:,0:,0:)              :: Q
    REAL(DOUBLE)                                 :: ZAC2,ZBC2,R,Delta,Ka,Kb,Alpha,AlphaInv, &
                                                     Left,Rhgt,Cntr,KaKb,DistToZero,CC,Ov
    !--------------------------------------------------------------------------------------------------------------
    IF(ABS(Ka)<1D-8.AND.ABS(Kb)<1D-8)THEN
       Q=Zero
       AlphaInv=One/Alpha
       Ov=EXP(-ZAC2-ZBC2)
       ! Only values with LambdaA==LambdaB==0 accrue
       Q(0,0,0)=Half*Ov*SQRT(Pi*AlphaInv)
       Q(0,0,1)=Half*Ov*AlphaInv
       DO Ehn=1,N-1
          Q(0,0,Ehn+1)=Half*AlphaInv*DBLE(Ehn)*Q(0,0,Ehn-1)
       ENDDO
    ELSE
       KaKb=(Ka+Kb)
       Discriminator=4D0*Alpha*(ZAC2+ZBC2+LOG(Tau))
       IF(KaKb**2<Discriminator)THEN
          Q=Zero
       ELSE
          Rhgt=(KaKb+SQRT(KaKb**2-Discriminator))/(Two*Alpha)
          IGrid=CEILING(Rhgt)
!          WRITE(*,*)' IGRID = ',IGRID
          IF(IGrid>Infty)THEN
             WRITE(*,*)' Rhgt = ',Rhgt
             STOP
          ENDIF
          DO I=1,NPts
             XA(I)=KA*Points(I,IGrid)
             XB(I)=KB*Points(I,IGrid)
             EX(I)=EXP(-Alpha*Points(I,IGrid)**2)
          ENDDO
          CALL BesselM(EllA,NPts,ZAC2,XA,MA)
          CALL BesselM(EllB,NPts,ZBC2,XB,MB)
          DO Ehn=0,N
             DO LambdaA=0,EllA
                DO LambdaB=0,EllB
                   Q(LambdaA,LambdaB,Ehn)=Zero
                   ! Accrue the integral from small to large
                   DO I=NPts,1,-1
                      Q(LambdaA,LambdaB,Ehn)=Q(LambdaA,LambdaB,Ehn)+Weights(I,IGrid) &
                                            *MA(I,LambdaA)*MB(I,LambdaB)*EX(I)*Points(I,IGrid)**(Ehn)
!                      IF(IGRID==1)THEN
!                      WRITE(*,44)I,Points(I),MA(I,LambdaA),MB(I,LambdaB),EX(I)
!                      ENDIF
!44                    format(I3,5(3x,D20.10))

                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF
#ifdef MAXDEBUG2
    IF(IGRID==1)THEN
    write(*,*)'=============================================================='
    WRITE(*,*)' KA = ',KA,' KB = ',KB,' Alpha = ',Alpha
    WRITE(*,*)' ZZABC2 = ',ZAC2+ZBC2
    write(*,*)' Ell = ',EllA,' EllB = ',EllB,' N = ',N
    DO Ehn=0,N
       DO LambdaA=0,EllA
          DO LambdaB=0,EllB
             WRITE(*,*)Ehn,LambdaA,LambdaB,' Q = ',Q(LambdaA,LambdaB,Ehn)
          ENDDO
       ENDDO
    ENDDO
    ENDIF
#endif
  END SUBROUTINE RadialTwo

  SUBROUTINE BesselM(L,NPts,G2,X,M)
    INTEGER                            :: L,NPts,I,J
    REAL(DOUBLE),DIMENSION(1:NPts)     :: X
    REAL(DOUBLE),DIMENSION(1:NPts,0:L) :: M
    REAL(DOUBLE)                       :: G2,EPlus,EMnus,TwoJMnsOne
    SELECT CASE(L)         
    CASE(0)       
       DO I=1,NPts
          IF(X(I)==Zero)THEN
             M(I,0)=One
          ELSE
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             M(I,0)=Half*(EPlus-EMnus)/X(I)
          ENDIF
       ENDDO
    CASE(1)
       DO I=1,NPts
          IF(X(I)==Zero)THEN
             M(I,0)=One
             M(I,1)=Zero
          ELSE
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             M(I,0)=Half*(EPlus-EMnus)/X(I)
             M(I,1)=Half*(EPlus*(X(I)-One)+EMnus*(X(I)+One))/(X(I)*X(I))
          ENDIF
       ENDDO
    CASE DEFAULT
       DO I=1,NPts
          IF(X(I)==Zero)THEN
             M(I,0)=One
             M(I,1)=Zero
          ELSE
             EPlus=Exp( X(I)-G2)
             EMnus=Exp(-X(I)-G2)
             M(I,0)=Half*(EPlus-EMnus)/X(I)
             M(I,1)=Half*(EPlus*(X(I)-One)+EMnus*(X(I)+One))/(X(I)*X(I))
          ENDIF
       ENDDO
       DO J=2,L
          TwoJMnsOne=DBLE(2*J-1)
          DO I=1,NPts          
             IF(X(I)==Zero)THEN
                M(I,J)=Zero
             ELSE
                M(I,J)=M(I,J-2)-TwoJMnsOne*M(I,J-1)/X(I)
             ENDIF
          ENDDO
       ENDDO
    END SELECT
!    WRITE(*,*)' ======================================================'
!    WRITE(*,*)' Ell = ',l,' Ell = ',l,' Ell = ',l,' Ell = ',l
!    WRITE(*,*)' X = ',X
!    WRITE(*,*)' M = ',M
  END SUBROUTINE BesselM

END MODULE ECPBlock
