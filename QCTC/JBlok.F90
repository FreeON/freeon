!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE BLOCK OF THE COULOMB MATRIX
!==============================================================================
MODULE JBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE Thresholding
  USE AtomPairs
  USE BraKetBloks
  USE Multipoles
  IMPLICIT NONE
CONTAINS
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
  SUBROUTINE VMD(zq,Px,Py,Pz,Omega,Upq,Tol,Rho)  
    INTEGER                           :: zq
    REAL(DOUBLE)                      :: Px,Py,Pz,Omega,Upq,Tol
    TYPE(HGRho)                       :: Rho
!
    HGKet%D = Zero
!
  END SUBROUTINE VMD
!------------------------------------------------------------------------------
! Calculate JBlok
!------------------------------------------------------------------------------
  FUNCTION JBlok(BS,MD,Pair,Rho) RESULT(JVck)
    TYPE(BSET)                              :: BS
    TYPE(DBL_RNK4)                          :: MD
    TYPE(AtomPair)                          :: Pair
    TYPE(HGRho)                             :: Rho
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: JVck
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: JBlk
!
    INTEGER                                 :: NBFA,NBFB,KA,KB
    REAL(DOUBLE)                            :: Ax,Ay,Az,Bx,By,Bz,AB2 
    REAL(DOUBLE)                            :: Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz,ZetaA,       &
                                               ZetaB,EtaAB,EtaABIn,XiA,XiB,XiAB,ExpAB,       &
                                               BraEstimate,RTE,RPE,Omega,Upq,Tol,TwoENeglect
    INTEGER                                 :: CFA,CFB,PFA,PFB,StartIA,StartIB,StopIA,StopIB &
                                              ,MaxLA,MaxLB,IA,IB,LMN
    INTEGER                                 :: MaxLMN,LBra,LenBra,zq,oq,or,NQ,LKet,LCode
#ifdef PERIODIC
    INTEGER                                 :: NC
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: JBlk_MM,JBlk_QQ
#endif
!
    KA   = Pair%KA
    KB   = Pair%KB
    NBFA = Pair%NA
    NBFB = Pair%NB
    Ax   = Pair%A(1)
    Ay   = Pair%A(2)
    Az   = Pair%A(3)
    Bx   = Pair%B(1)
    By   = Pair%B(2)
    Bz   = Pair%B(3)
    AB2  = Pair%AB2
!
    JBlk    = Zero
#ifdef PERIODIC
    JBlk_MM = Zero
    JBlk_QQ = Zero
#endif 
!
    MaxLMN      = BS%LMNLen
    TwoENeglect = Thresholds%TwoE
!   
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A
       MaxLA=BS%ASymm%I(2,CFA,KA) 
       StartIA = CFBlokDex(BS,CFA,KA)+1
       StopIA  = StartIA + BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)                    ! Loop over contracted function B
          MaxLB=BS%ASymm%I(2,CFB,KB)
          StartIB = CFBlokDex(BS,CFB,KB)+1
          StopIB  = StartIB + BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)
!
          LBra   = MaxLA+MaxLB
          LenBra = LHGTF(LBra)
!
          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in 
             DO PFB=1,BS%NPFnc%I(CFB,KB)          ! contracted functions A and B      
                ZetaA   = BS%Expnt%D(PFA,CFA,KA)
                ZetaB   = BS%Expnt%D(PFB,CFB,KB)
                EtaAB   = ZetaA+ZetaB 
                EtaABIn = One/EtaAB
                XiAB    = ZetaA*ZetaB*EtaABIn
                IF(TestPrimPair(XiAB,Pair%AB2))THEN
                   XiA     = ZetaA*EtaABIn
                   XiB     = ZetaB*EtaABIn
                   ExpAB   = EXP(-XiAB*AB2)
                   Px=(XiA*Ax+XiB*Bx)
                   Py=(XiA*Ay+XiB*By)
                   Pz=(XiA*Az+XiB*Bz)
                   PAx=Px-Ax
                   PAy=Py-Ay
                   PAz=Pz-Az
                   PBx=Px-Bx
                   PBy=Py-By
                   PBz=Pz-Bz 
!
!                  McMurchie-Davidson 2 term coefficients
!
                   CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D,PAx,PBx,PAy,PBy,PAz,PBz) 
!
!                  Primitive coefficients in a HG representationj
!
                   CALL SetBraBlok(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Phase_O=.TRUE.)
!
!                  Integral estimates (calculation should definately be moved out of inner loop!) 
!
                   CALL BraEst(StartIA,StopIA,StartIB,StopIB,LBra,XiAB,BraEstimate)
!
!                  Calculate the Tolerance
!
                   Tol = TwoENeglect/MAX(1.D-32,BraEstimate)
#ifdef PERIODIC
!
!                  Sum over Cells
!

                   DO NC=1,CSMM1%NCells
                      Px=(XiA*Ax+XiB*Bx)+CSMM1%CellCarts%D(1,NC)
                      Py=(XiA*Ay+XiB*By)+CSMM1%CellCarts%D(2,NC)
                      Pz=(XiA*Az+XiB*Bz)+CSMM1%CellCarts%D(3,NC)
#endif
!
!                     Direct Contribution
!
                      HGKet%D = Zero
                      DO zq=1,Rho%NExpt
                         NQ = Rho%NQ%I(zq)
                         IF(NQ /= 0) THEN
                            RTE   = Rho%Expt%D(zq)*EtaAB
                            RPE   = Rho%Expt%D(zq)+EtaAB
                            Omega = RTE/RPE
                            Upq   = TwoPi5x2/(RTE*DSQRT(RPE))
                            oq    = Rho%OffQ%I(zq)+1
                            or    = Rho%OffR%I(zq)+1
                            LKet  = Rho%Lndx%I(zq)
                            LCode = 10*LBra+LKet
                            INCLUDE 'VMD/VMDBlok.Inc'
                         ENDIF
                      ENDDO

!
!                     Update the Matrix Block
!
                      DO IA = StartIA,StopIA
                         DO IB = StartIB,StopIB
                            DO LMN=1,LenBra
                               JBlk(IA,IB) = JBlk(IA,IB) + HGBra%D(LMN,IA,IB)*HGKet%D(LMN)
                            ENDDO
                         ENDDO
                      ENDDO
#ifdef PERIODIC
                   ENDDO
!
!                  Calculate the Multipole Contribution to the Matrix Element
!                  Contract the Primative MM  with the density MM
!
                   IF(Dimen > 0) THEN
                      Px=(XiA*Ax+XiB*Bx)
                      Py=(XiA*Ay+XiB*By)
                      Pz=(XiA*Az+XiB*Bz)
                      HGBra%D = Zero
                      CALL SetBraBlok(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Phase_O=.FALSE.)
                      DO IA = StartIA,StopIA
                         DO IB = StartIB,StopIB
                            JBlk_MM(IA,IB) = JBlk_MM(IA,IB)+CTraxBraKet(LBra,EtaAB,Px,Py,Pz,HGBra%D(:,IA,IB))
                         ENDDO
                      ENDDO
!
!                     Calculate the Cartiesian Moment Dipole and Quadipole Corrections
!
                      IF(Dimen > 1) THEN
                         DO IA = StartIA,StopIA
                            DO IB = StartIB,StopIB
                               JBlk_QQ(IA,IB) = JBlk_QQ(IA,IB)+QTraxBraKet(LBra,EtaAB,Px,Py,Pz,HGBra%D(:,IA,IB))
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDIF
#endif
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
#ifdef PERIODIC
    JBlk = JBlk+JBlk_MM+JBlk_QQ
#endif
!
    JVck = BlockToVect(NBFA,NBFB,JBlk)
!
  END FUNCTION JBlok
!------------------------------------------------------------------------------
! Calculate Vblok
!------------------------------------------------------------------------------
  FUNCTION VBlok(AtA,BS,MD,Rho)
    TYPE(BSET)                           :: BS
    TYPE(DBL_RNK4)                       :: MD
    TYPE(HGRho)                          :: Rho
!
    INTEGER                              :: AtA,zqA,oqA,orA,zq,oq,or,NQ,LCode,LBra,LKet
    REAL(DOUBLE)                         :: VBlok,Px,Py,Pz,BraEstimate,RTE,RPE,   &
                                            Omega,Upq,Tol,TwoENeglect,ExpA,CoefA
#ifdef PERIODIC
    INTEGER                              :: NC
    REAL(DOUBLE)                         :: VBlok_MM,VBlok_QQ
#endif
!
    zqA   = Rho%NExpt
    oqA   = Rho%OffQ%I(zqA)
    orA   = Rho%OffR%I(zqA)
    ExpA  = Rho%Expt%D(zqA)
    LBra  = Rho%Lndx%I(zqA)
    CoefA = Rho%Co%D(orA+AtA)
!
!   Initialize
!
    VBlok    = Zero
#ifdef PERIODIC
    VBlok_MM = Zero    
    VBlok_QQ = Zero
#endif
    TwoENeglect = Thresholds%TwoE
    BraEstimate = Rho%Est%D(oqA+AtA)
    Tol         = TwoENeglect/MAX(1.D-32,BraEstimate)
!
#ifdef PERIODIC
!
!   Sum over Cells
!
    DO NC=1,CSMM1%NCells
       Px=Rho%Qx%D(oqA+AtA)+CSMM1%CellCarts%D(1,NC)
       Py=Rho%Qy%D(oqA+AtA)+CSMM1%CellCarts%D(2,NC)
       Pz=Rho%Qz%D(oqA+AtA)+CSMM1%CellCarts%D(3,NC)
       IF(InCell_CellSet(CSMM1,NC,Zero,Zero,Zero)) THEN
          Rho%Co%D(orA+AtA) = Zero
       ELSE
          Rho%Co%D(orA+AtA) = CoefA
       ENDIF
!
!      Calculate the Direct Piece
!
       HGKet%D = Zero
       DO zq=1,Rho%NExpt
          NQ = Rho%NQ%I(zq)
          IF(NQ /= 0) THEN
             RTE   = Rho%Expt%D(zq)*ExpA
             RPE   = Rho%Expt%D(zq)+ExpA
             Omega = RTE/RPE
             Upq   = TwoPi5x2/(RTE*DSQRT(RPE))
             oq    = Rho%OffQ%I(zq)+1
             or    = Rho%OffR%I(zq)+1
             LKet  = Rho%Lndx%I(zq)
             LCode = 10*LBra+LKet
             INCLUDE 'VMD/VMDBlok.Inc'
          ENDIF
       ENDDO
       VBlok = VBlok+CoefA*HGKet%D(1)
    ENDDO
!
!   Do the Multipole Piece
!
    IF(Dimen > 0) THEN
       Px    = Rho%Qx%D(oqA+AtA)
       Py    = Rho%Qy%D(oqA+AtA)
       Pz    = Rho%Qz%D(oqA+AtA)
       HGBra%D(1,1,1) = CoefA
       VBlok_MM  = VBlok_MM+CTraxBraKet(0,ExpA,Px,Py,Pz,HGBra%D(:,1,1))
       IF(Dimen > 1) THEN
          VBlok_QQ = VBlok_QQ+QTraxBraKet(0,ExpA,Px,Py,Pz,HGBra%D(:,1,1))
       ENDIF
       VBlok = VBlok+VBlok_MM+VBlok_QQ
    ENDIF
!
#else
    Px    = Rho%Qx%D(oqA+AtA)
    Py    = Rho%Qy%D(oqA+AtA)
    Pz    = Rho%Qz%D(oqA+AtA)
    Rho%Co%D(orA+AtA) = Zero
!
!   Do the Direct Piece
!
    HGKet%D = Zero
    DO zq=1,Rho%NExpt
       NQ = Rho%NQ%I(zq)
       IF(NQ /= 0) THEN
          RTE   = Rho%Expt%D(zq)*ExpA
          RPE   = Rho%Expt%D(zq)+ExpA
          Omega = RTE/RPE
          Upq   = TwoPi5x2/(RTE*DSQRT(RPE))
          Tol   = TwoENeglect/(ExpA*BraEstimate)
          oq    = Rho%OffQ%I(zq)+1
          or    = Rho%OffR%I(zq)+1
          LKet  = Rho%Lndx%I(zq)
          LCode = 10*LBra+LKet
          INCLUDE 'VMD/VMDBlok.Inc'
       ENDIF
    ENDDO
    VBlok = VBlok+CoefA*HGKet%D(1)
#endif
!
!   Reset
!
    Rho%Co%D(orA+AtA) = CoefA
  END FUNCTION VBlok
!--------------------------------------------------------------
! Estimate the Integrals from the Primative Distributions
!--------------------------------------------------------------
  SUBROUTINE BraEst(StartIA,StopIA,StartIB,StopIB,LBra,Zeta,BraEstimate)
    INTEGER                      :: StartIA,StopIA,StartIB,StopIB,LBra
    INTEGER                      :: IA,IB,LenBra,L
    REAL(DOUBLE)                 :: Zeta,SqUpp,BraEstimate,Est
!
    SqUpp    = Sqrt2Pi5x2*Zeta**(-FiveFourths)
    LenBra = LHGTF(LBra)
!
!   Generate The R Function
!
    CALL GenerateRfun(Zeta,LBra)
!
!   Calculate the BraEst
!    
    BraEstimate = Zero
    DO IA=StartIA,StopIA
       DO IB=StartIB,StopIB
          Est = Estimate(LBra,LenBra,HGBra%D(1:LenBra,IA,IB))
          BraEstimate = MAX(BraEstimate,Est)
       ENDDO
    ENDDO
    BraEstimate =SqUpp*BraEstimate
  END SUBROUTINE BraEst
!
END MODULE JBlock





