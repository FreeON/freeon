MODULE JBlock
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE Thresholding
  USE AtomPairs
  USE BraKetBloks
  USE MMoments
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
!
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
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: JBlk_MM
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
#ifdef PERIODIC
!
!                  Sum over Cells
!
                   DO NC=1,CSMM1%NCells
                      Px=(XiA*Ax+XiB*Bx)+CSMM1%CellCarts%D(1,NC)
                      Py=(XiA*Ay+XiB*By)+CSMM1%CellCarts%D(2,NC)
                      Pz=(XiA*Az+XiB*Bz)+CSMM1%CellCarts%D(3,NC)
#else
                      Px=(XiA*Ax+XiB*Bx)
                      Py=(XiA*Ay+XiB*By)
                      Pz=(XiA*Az+XiB*Bz)
#endif
                      PAx=Px-Ax
                      PAy=Py-Ay
                      PAz=Pz-Az
                      PBx=Px-Bx
                      PBy=Py-By
                      PBz=Pz-Bz 
!                     McMurchie-Davidson 2 term coefficients
                      CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D,PAx,PBx,PAy,PBy,PAz,PBz) 
!                     Primitive coefficients in a HG representationj
                      CALL SetBraBlok(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Phase_O=.TRUE.)
!                     Integral estimates (calculation should definately be moved out of inner loop!) 
                      CALL BraEst(StartIA,StopIA,StartIB,StopIB,LBra,XiAB,BraEstimate)
!                     Direct Contribution
                      HGKet%D = Zero
                      Tol     = TwoENeglect/MAX(1.D-32,BraEstimate)
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

!                      WRITE(*,*)' CFA = ',CFA,' CFB = ',CFB,' PFA = ',PFA,' PFB = ',PFB
!                      WRITE(*,*)' HGKet = ',HGKet%D(1:LenBra)
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
!
!                  >>>>>>>>> WHY WOULD YOU DO THIS TWICE CJ ????????
!
!                   Px=(XiA*Ax+XiB*Bx)
!                   Py=(XiA*Ay+XiB*By)
!                   Pz=(XiA*Az+XiB*Bz)
!                   PAx=Px-Ax
!                   PAy=Py-Ay
!                   PAz=Pz-Az
!                   PBx=Px-Bx
!                   PBy=Py-By
!                   PBz=Pz-Bz 
!
!
!                   CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D,PAx,PBx,PAy,PBy,PAz,PBz) 
!                  Calculate the  the Primative
!
!                   CALL GetPrimative(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD)
!
!                  Contract the Primative MM  with the density MM
!
                   DO IA = StartIA,StopIA
                      DO IB = StartIB,StopIB
                         JBlk_MM(IA,IB) = JBlk_MM(IA,IB) &
                              +CTraxBraKet(LBra,EtaAB,Px,Py,Pz,HGBra%D(:,IA,IB))
                      ENDDO
                   ENDDO
#endif
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
#ifdef PERIODIC
    JBlk = JBlk+JBlk_MM
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
    REAL(DOUBLE)                         :: VBlok_MM
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
#endif
    TwoENeglect = Thresholds%TwoE
    BraEstimate = Rho%Est%D(oqA+AtA)
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
       Tol     = TwoENeglect/MAX(1.D-32,ExpA*BraEstimate)
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
       CALL MD2TRR(0,0,0,0,ExpA,MD%D,Px,Px,Py,Py,Pz,Pz) 
       VBlok = VBlok+CoefA*MD%D(1,0,0,0)*MD%D(2,0,0,0)*MD%D(3,0,0,0)*HGKet%D(1)
    ENDDO
!
!   Do the Multipole Piece
!
    Px    = Rho%Qx%D(oqA+AtA)
    Py    = Rho%Qy%D(oqA+AtA)
    Pz    = Rho%Qz%D(oqA+AtA)
    HGBra%D(1,1,1) = CoefA
    VBlok_MM       = VBlok_MM+CTraxBraKet(0,ExpA,Px,Py,Pz,HGBra%D(:,1,1))
!
    VBlok          = VBlok+VBlok_MM
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
    Tol     = TwoENeglect/MAX(1.D-32,ExpA*BraEstimate)
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
    CALL MD2TRR(0,0,0,0,ExpA,MD%D,Px,Px,Py,Py,Pz,Pz) 
    VBlok = VBlok+CoefA*MD%D(1,0,0,0)*MD%D(2,0,0,0)*MD%D(3,0,0,0)*HGKet%D(1)
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
END MODULE JBlock





