MODULE RhoBlok
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraBloks
  IMPLICIT NONE
CONTAINS
!--------------------------------------------------------------
! Add nuclear centers to the density
!--------------------------------------------------------------
  SUBROUTINE AddNukes(GM,Rho)
    TYPE(CRDS)                 :: GM
    TYPE(HGRho)                :: Rho
    INTEGER                    :: IA,Iq,Ir
    REAL(DOUBLE)               :: DDelta
!---------------------------------------------------
!   Initialize the Nuclear centers and Constants
!
    DDelta = Half*(Rho%Expt%D(Rho%NExpt)/Pi)**(ThreeHalves)
!
    DO IA = 1,NAtoms
       Iq = Rho%OffQ%I(Rho%NExpt)+IA
       Ir = Rho%OffR%I(Rho%NExpt)+IA
       Rho%Qx%D(Iq)    = GM%Carts%D(1,IA)
       Rho%Qy%D(Iq)    = GM%Carts%D(2,IA)
       Rho%Qz%D(Iq)    = GM%Carts%D(3,IA)
       Rho%Co%D(Ir)  = (-GM%AtNum%I(IA))*DDelta
    ENDDO
  END SUBROUTINE AddNukes
!--------------------------------------------------------------
! Count the Number of Primative basis funtions contained in atom A and B
!--------------------------------------------------------------
  SUBROUTINE PrimCount(BS,Pair,Rho)
    TYPE(BSET)                 :: BS
    TYPE(AtomPair)             :: Pair
    TYPE(HGRho)                :: Rho
!
    INTEGER                    :: KA,KB,CFA,CFB,PFA,PFB,IE,Endex
    REAL(DOUBLE)               :: AB2,ZetaA,ZetaB,ZetaAB,XiAB,ExpAB
    LOGICAL                    :: AEQB
!      
    KA   = Pair%KA
    KB   = Pair%KB  
    AB2  = Pair%AB2 
    AEQB = Pair%SameAtom
!
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A 
       DO CFB=1,BS%NCFnc%I(KB)                    ! Loop over contracted function B          
          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in A
             DO PFB=1,BS%NPFnc%I(CFB,KB)          ! Loops over primitives in A
                ZetaA = BS%Expnt%D(PFA,CFA,KA)
                ZetaB = BS%Expnt%D(PFB,CFB,KB)
                ZetaAB= ZetaA+ZetaB 
                XiAB  = ZetaA*ZetaB/ZetaAB
                IF(TestPrimPair(XiAB,Pair%AB2))THEN
                   ExpAB = EXP(-XiAB*AB2)
!
!                  Determine the Exponent Index
!
                   Endex=0
                   DO IE=1,Rho%NExpt
                      IF(ABS(ZetaAB-Rho%Expt%D(IE))<1.0D-8) THEN
                         Endex = IE
                         EXIT
                      ENDIF
                   ENDDO
                   IF(Endex == 0) THEN
                      CALL MondoHalt(-100,' Problem in PrimCount, No Exponent found ') 
                   ELSE
                      Rho%NQ%I(Endex) = Rho%NQ%I(Endex) + 1
                   ENDIF
!
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
  END SUBROUTINE PrimCount
!--------------------------------------------------------------
! Calculate the Primative Distributions Generated from atom A and B
!--------------------------------------------------------------
  SUBROUTINE RhoBlk(AtA,AtB,BS,MD,Dmat,Pair,First,Rho)
    TYPE(BSET)                              :: BS
    TYPE(DBL_RNK4)                          :: MD
    TYPE(AtomPair)                          :: Pair
    TYPE(PrimPair)                          :: Prim
    TYPE(HGRho)                             :: Rho
    REAL(DOUBLE)                            :: FacAtom
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: Dmat
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: DD
    INTEGER                                 :: KA,KB,NBFA,NBFB,CFA,CFB,PFA,PFB, &
                                               IndexA,StartLA,StopLA,MaxLA, &
                                               IndexB,StartLB,StopLB,MaxLB, &
                                               IE,Endex,Qndex,Rndex,LMN,LMNA,LMNB, &
                                               IA,IB,LA,MA,NA,LB,MB,NB,LAB,MAB,NAB, &
                                               LenKet,AtA,AtB
    REAL(DOUBLE)                            :: ZetaA,ZetaB,ZetaAB,ZetaIn,XiAB,ExpAB, &
                                               AB2,Ax,Ay,Az,Bx,By,Bz, &
                                               Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz, &
                                               Ex,Exy,Exyz,CA,CB,MaxAmp
!
    TYPE(INT_VECT),SAVE                     :: OffQ,OffR
    LOGICAL                                 :: First,AEQB
!
    IF(First) THEN
       CALL New(OffQ,Rho%NExpt)
       CALL New(OffR,Rho%NExpt)
       DO IA = 1,Rho%NExpt
          OffQ%I(IA) = Rho%OffQ%I(IA)
          OffR%I(IA) = Rho%OffR%I(IA)
       ENDDO
       First = .FALSE.
    ENDIF
!
    KA   = Pair%KA
    KB   = Pair%KB  
    NBFA = Pair%NA
    NBFB = Pair%NB
    AEQB = Pair%SameAtom
!----------------------------------
    Prim%A=Pair%A
    Prim%B=Pair%B
    Prim%AB2=Pair%AB2
    Prim%KA=Pair%KA
    Prim%KB=Pair%KB
!----------------------------------
!
    DD = VectToBlock(NBFA,NBFB,Dmat)
! 
    DO CFA=1,BS%NCFnc%I(KA)                          ! Loop over contracted function A 
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)        
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       DO CFB=1,BS%NCFnc%I(KB)                       ! Loop over contracted function B
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
!----------------------------------
          Prim%CFA=CFA
          Prim%CFB=CFB
          Prim%Ell=MaxLA+MaxLB
!----------------------------------
          DO PFA=1,BS%NPFnc%I(CFA,KA)                ! Loops over primitives in CFA and CFB
!----------------------------------
             Prim%PFA=PFA 
!----------------------------------
             DO PFB=1,BS%NPFnc%I(CFB,KB)             ! Loops over primitives in CFB
!----------------------------------
                Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                Prim%Zeta=Prim%ZA+Prim%ZB
                Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
                IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                   Prim%PFB=PFB
!----------------------------------
!               Set primitive values
!               Primitive coefficients in a HG representation
!--------------------------------------------------
                   MaxAmp=SetBraBlok(Prim,BS)
                   Endex = 0
                   DO IE=1,Rho%NExpt-1
                      IF(ABS(Prim%Zeta-Rho%Expt%D(IE))<1.0D-8) THEN
                         Endex = IE
                         GOTO 100
                      ENDIF
                   ENDDO
100                CONTINUE
                   IF(Endex .EQ. 0) THEN
                      CALL MondoHalt(-100,' Problem in RhoBlok, No Exponent found ') 
                   ENDIF
!
                   Qndex  = OffQ%I(Endex)+1
                   Rndex  = OffR%I(Endex)
                   LenKet = LHGTF(Rho%Lndx%I(Endex))
!
                   OffQ%I(Endex) = OffQ%I(Endex) + 1 
                   OffR%I(Endex) = OffR%I(Endex) + LenKet
!
!                  Store the position of the distribution
!
                   Rho%Qx%D(Qndex)=Prim%P(1)
                   Rho%Qy%D(Qndex)=Prim%P(2)
                   Rho%Qz%D(Qndex)=Prim%P(3)
!
!                  Calculate and Store the Coefficients of the Distribution
!
                   IA=IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      LA=BS%LxDex%I(LMNA)
                      MA=BS%LyDex%I(LMNA)
                      NA=BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         LB=BS%LxDex%I(LMNB)
                         MB=BS%LyDex%I(LMNB)
                         NB=BS%LzDex%I(LMNB)


!WRITE(66,22)AtA,AtB,BS%CCoef%D(LMNA,PFA,CFA,KA),BS%CCoef%D(LMNB,PFB,CFB,KB),Prim%ZA,Prim%ZB,Prim%A,Prim%B
!22 FORMAT('test3[',I3,',',I3,']={ca->',F10.7,',cb->',F10.7,               &
!                       ',za->',F10.7,',zb->',F10.7,               &
!                       ',ax->',F10.7,',ay->',F10.7,',az->',F10.7, &
!                       ',bx->',F10.7,',by->',F10.7,',bz->',F10.7,'}; \n')


                         DO LMN=1,LHGTF(LA+MA+NA+LB+MB+NB)
                            Rho%Co%D(Rndex+LMN)=Rho%Co%D(Rndex+LMN)+HGBra%D(LMN,IA,IB)*DD(IA,IB)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE RhoBlk
!
  SUBROUTINE Integrate_HGRho(Rho)
    TYPE(HGRho)                      :: Rho
    INTEGER                          :: zq,iq,oq,or,NQ,jadd,Ell,LenKet
    REAL(DOUBLE)                     :: RhoSumE,RhoSumN,Dex
!
    RhoSumE = Zero
    DO zq = 1,Rho%NExpt-1
       Dex    = Rho%Expt%D(zq)
       NQ     = Rho%NQ%I(zq)
       oq     = Rho%OffQ%I(zq)
       or     = Rho%OffR%I(zq)
       Ell    = Rho%Lndx%I(zq) 
       LenKet = LHGTF(Ell)
       DO iq = 1,NQ
          jadd = or+(iq-1)*LenKet
          RhoSumE = RhoSumE + Rho%Co%D(jadd+1)*((Pi/Dex)**ThreeHalves)
       ENDDO
    ENDDO
!
    RhoSumN = Zero
    zq     = Rho%NExpt
    Dex    = Rho%Expt%D(zq)
    NQ     = Rho%NQ%I(zq)
    oq     = Rho%OffQ%I(zq)
    or     = Rho%OffR%I(zq)
    Ell    = Rho%Lndx%I(zq) 
    LenKet = LHGTF(Ell)
    DO iq = 1,NQ
       jadd = or+(iq-1)*LenKet
       RhoSumN = RhoSumN + Rho%Co%D(jadd+1)*((Pi/Dex)**ThreeHalves)
    ENDDO
!
    WRITE(*,*) ' Int[Rho_E] = ',RhoSumE
    WRITE(*,*) ' Int[Rho_N] = ',RhoSumN   
    WRITE(*,*) ' Int[Rho_T] = ',RhoSumN+RhoSumE
  END SUBROUTINE Integrate_HGRho
!--------------------------------------------------------------
! Calculate the Number of Distributions from NQ
!--------------------------------------------------------------
  FUNCTION  CalNDist(Rho)
    TYPE(HGRho)                :: Rho
    INTEGER                    :: I,CalNDist
    CalNDist = 0
    DO I=1,Rho%NExpt
       CalNDist = CalNDist+Rho%NQ%I(I)
    ENDDO
  END FUNCTION CalNDist
!--------------------------------------------------------------
! Calculate the Number of Coefficients  from NQ and Lndx
!--------------------------------------------------------------
  FUNCTION CalNCoef(Rho)
    TYPE(HGRho)                :: Rho
    INTEGER                    :: I,CalNCoef
    CalNCoef = 0
    DO I=1,Rho%NExpt
       CalNCoef = CalNCoef+Rho%NQ%I(I)*LHGTF(Rho%Lndx%I(I))
    ENDDO
  END FUNCTION CalNCoef
!--------------------------------------------------------------
! Calculate OffQ from NQ
!--------------------------------------------------------------
  FUNCTION CalOffQ(Rho)
    TYPE(HGRho)                  :: Rho
    INTEGER,DIMENSION(Rho%NExpt) :: CalOffQ   
    INTEGER                      :: I,J,Iq
!
    DO I = 1,Rho%NExpt
       Iq = 0
       DO J = 1,I-1
          Iq = Iq + Rho%NQ%I(J)
       ENDDO
       CalOffQ(I) = Iq
    ENDDO
  END FUNCTION CalOffQ
!--------------------------------------------------------------
! Calculate OffR from NQ and Lndx
!--------------------------------------------------------------
  FUNCTION CalOffR(Rho)
    TYPE(HGRho)                  :: Rho
    INTEGER,DIMENSION(Rho%NExpt) :: CalOffR   
    INTEGER                      :: I,J,Ir
    DO I = 1,Rho%NExpt
       Ir = 0
       DO J = 1,I-1
          Ir = Ir + Rho%NQ%I(J)*LHGTF(Rho%Lndx%I(J))
       ENDDO
       CalOffR(I) = Ir
    ENDDO
  END FUNCTION CalOffR
!
END MODULE RhoBlok
