MODULE RhoBlok
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE BraKetBloks
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
    INTEGER                    :: KA,KB,CFA,CFB,PFA,PFB,IE,Endex,CFBbeg,PFBbeg
    REAL(DOUBLE)               :: AB2,ZetaA,ZetaB,ZetaAB,XiAB,ExpAB
    LOGICAL                    :: AEQB
!      
    KA   = Pair%KA
    KB   = Pair%KB  
    AB2  = Pair%AB2 
    AEQB = Pair%SameAtom
!
    DO CFA=1,BS%NCFnc%I(KA)                       ! Loop over contracted function A 
       IF(AEQB) THEN 
          CFBbeg = CFA
       ELSE
          CFBbeg = 1
       ENDIF
       DO CFB=CFBbeg,BS%NCFnc%I(KB)               ! Loop over contracted function B          
          DO PFA=1,BS%NPFnc%I(CFA,KA)             ! Loops over primitives in 
             IF(AEQB .AND. CFA==CFB) THEN                        ! contracted functions A and B
                PFBbeg = PFA
             ELSE
                PFBbeg = 1
             ENDIF
             DO PFB=PFBbeg,BS%NPFnc%I(CFB,KB)     
                ZetaA = BS%Expnt%D(PFA,CFA,KA)
                ZetaB = BS%Expnt%D(PFB,CFB,KB)
                ZetaAB= ZetaA+ZetaB 
                XiAB  = ZetaA*ZetaB/ZetaAB
                IF(TestPrimPair(XiAB,Pair%AB2))THEN
!
                   ExpAB = EXP(-XiAB*AB2)
!                  Determine the Exponent Index
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
  SUBROUTINE RhoBlk(BS,MD,Dmat,Pair,First,Rho)
    TYPE(BSET)                              :: BS
    TYPE(DBL_RNK4)                          :: MD
    TYPE(AtomPair)                          :: Pair
    TYPE(HGRho)                             :: Rho
    REAL(DOUBLE)                            :: FacAtom,FacState,FacPrim
!
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB) :: Dmat
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB) :: DD
    INTEGER                                 :: KA,KB,NBFA,NBFB,CFA,CFB,PFA,PFB, &
                                               IndexA,StartLA,StopLA,MaxLA, &
                                               IndexB,StartLB,StopLB,MaxLB, &
                                               IE,Endex,Qndex,Rndex,LMN,LMNA,LMNB, &
                                               IA,IB,LA,MA,NA,LB,MB,NB,LAB,MAB,NAB, &
                                               CFBbeg,PFBbeg,LenKet
    REAL(DOUBLE)                            :: ZetaA,ZetaB,ZetaAB,ZetaIn,XiAB,ExpAB, &
                                               AB2,Ax,Ay,Az,Bx,By,Bz, &
                                               Px,Py,Pz,PAx,PAy,PAz,PBx,PBy,PBz, &
                                               Ex,Exy,Exyz,CA,CB,ComCoef
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
    Ax   = Pair%A(1)
    Ay   = Pair%A(2)
    Az   = Pair%A(3)
    Bx   = Pair%B(1)
    By   = Pair%B(2)
    Bz   = Pair%B(3)
    AB2  = Pair%AB2  
    AEQB = Pair%SameAtom
!
    DD = VectToBlock(NBFA,NBFB,Dmat)
!
    DO CFA=1,BS%NCFnc%I(KA)                         ! Loop over contracted function A 
       IndexA  = CFBlokDex(BS,CFA,KA)
       StartLA = BS%LStrt%I(CFA,KA)        
       StopLA  = BS%LStop%I(CFA,KA)
       MaxLA   = BS%ASymm%I(2,CFA,KA)
       IF(AEQB) THEN 
          CFBbeg = CFA
       ELSE
          CFBbeg = 1
       ENDIF
       DO CFB=CFBbeg,BS%NCFnc%I(KB)                    ! Loop over contracted function B
          IndexB  = CFBlokDex(BS,CFB,KB)
          StartLB = BS%LStrt%I(CFB,KB)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          DO PFA=1,BS%NPFnc%I(CFA,KA)                  ! Loops over primitives in CFA and CFB
             IF(AEQB .AND. CFA==CFB) THEN                     
                PFBbeg = PFA
             ELSE
                PFBbeg = 1
             ENDIF
             DO PFB=PFBbeg,BS%NPFnc%I(CFB,KB)          ! Loops over primitives in CFB
                ZetaA=BS%Expnt%D(PFA,CFA,KA)
                ZetaB=BS%Expnt%D(PFB,CFB,KB)
                ZetaAB=ZetaA+ZetaB 
                ZetaIn=One/ZetaAB
                XiAB =ZetaA*ZetaB*ZetaIn
                IF(TestPrimPair(XiAB,Pair%AB2))THEN
                   ExpAB=EXP(-XiAB*AB2)
!                  Determine the Counting Factors
                   IF(AEQB) THEN
                      FacAtom = one
                      IF(CFA == CFB) THEN
                         FacState = one
                         IF(PFA == PFB) THEN
                            FacPrim = one
                         ELSE
                            FacPrim = two
                         ENDIF
                      ELSE
                         FacState = two
                         FacPrim  = one
                      ENDIF
                   ELSE
                      FacAtom  = two
                      FacState = one
                      FacPrim  = one
                   ENDIF
                   ExpAB = FacAtom*FacState*FacPrim*ExpAB
!
!                  Determine the Indexs
!
                   Endex = 0
                   DO IE=1,Rho%NExpt-1
                      IF(ABS(ZetaAB-Rho%Expt%D(IE))<1.0D-8) THEN
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
!                  Calculate and Store the Position of the Distribution
!
                   Px=(ZetaA*Ax+ZetaB*Bx)*ZetaIn
                   Py=(ZetaA*Ay+ZetaB*By)*ZetaIn
                   Pz=(ZetaA*Az+ZetaB*Bz)*ZetaIn
                   PAx=Px-Ax
                   PAy=Py-Ay
                   PAz=Pz-Az
                   PBx=Px-Bx
                   PBy=Py-By
                   PBz=Pz-Bz
                   Rho%Qx%D(Qndex) = Px
                   Rho%Qy%D(Qndex) = Py
                   Rho%Qz%D(Qndex) = Pz
!
!                  Calculate and Store the Coefficients of the Distribution
!
                   CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,ZetaAB,MD%D, &
                        PAx,PBx,PAy,PBy,PAz,PBz) 
!
                   IA=IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      LA=BS%LxDex%I(LMNA)
                      MA=BS%LyDex%I(LMNA)
                      NA=BS%LzDex%I(LMNA)
                      CA=BS%CCoef%D(LMNA,PFA,CFA,KA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         LB=BS%LxDex%I(LMNB)
                         MB=BS%LyDex%I(LMNB)
                         NB=BS%LzDex%I(LMNB)
                         CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                         ComCoef = CA*CB*ExpAB*DD(IA,IB)
                         DO LAB=0,LA+LB
                            Ex = ComCoef*MD%D(1,LA,LB,LAB)
                            DO MAB=0,MA+MB
                               Exy = Ex*MD%D(2,MA,MB,MAB)
                               DO NAB=0,NA+NB
                                  Exyz = Exy*MD%D(3,NA,NB,NAB)
                                  LMN = LMNdex(LAB,MAB,NAB)
                                  Rho%Co%D(Rndex+LMN)=Rho%Co%D(Rndex+LMN)+Exyz
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
!
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE RhoBlk
!--------------------------------------------------------------
! Estimate the Integrals from the Primative Distributions
!--------------------------------------------------------------
  SUBROUTINE RhoEst(Rho)
    TYPE(HGRho)                  :: Rho
    INTEGER                      :: MaxL,zq,iq,oq,or,LKet,LenKet,qadd,radd,NQ
    REAL(DOUBLE)                 :: Zeta,SqUqq
!
    IF(Rho%AllocRE==ALLOCATED_TRUE) THEN
!
       MaxL = 1
       DO zq = 1,Rho%NExpt
          IF(MaxL < Rho%Lndx%I(zq)) MaxL = Rho%Lndx%I(zq)
       ENDDO
!
       DO zq = 1,Rho%NExpt
          NQ     = Rho%NQ%I(zq)
          IF(NQ/=0) THEN
             oq    =Rho%OffQ%I(zq)
             or    =Rho%OffR%I(zq)
             Zeta  =Rho%Expt%D(zq)
             LKet  =Rho%Lndx%I(zq)
             LenKet=LHGTF(LKet)
             SqUqq =Sqrt2Pi5x2*Zeta**(-FiveFourths)
!            Generate The R and AuxR Function
             CALL GenerateRfun(Zeta,LKet)
!            Calculate the Estimate
             DO iq=1,NQ
                qadd=oq+iq
                radd=or+(iq-1)*LenKet
                Rho%Est%D(qadd)=SqUqq*Estimate(LKet,LenKet,Rho%Co%D(radd+1:))  
             ENDDO
!            Sort the estimates
             CALL SortRhoEst(zq,NQ,LenKet,Rho)
          ENDIF
       ENDDO
       CALL Delete(AuxRfun)
       CALL Delete(Rfun)
    ENDIF
!       DO zq = 1,Rho%NExpt
!          NQ     = Rho%NQ%I(zq)
!          IF(NQ/=0) THEN
!             oq    =Rho%OffQ%I(zq)
!             or    =Rho%OffR%I(zq)
!             Zeta  =Rho%Expt%D(zq)
!             LKet  =Rho%Lndx%I(zq)
!             LenKet=LHGTF(LKet)
!             SqUqq   =Sqrt2Pi5x2*Zeta**(-FiveFourths)
!             WRITE(11,*)' ====================================================='
!             WRITE(11,*)' oq, or, zeta, lket = ',oq,or,zeta,lket
!             WRITE(11,*)' SqUqq = ',SqUqq
!             DO iq=1,NQ
!@                qadd=oq+iq
!!                WRITE(11,*)iq,Rho%Est%D(qadd)
!             ENDDO
!          ENDIF
!       ENDDO

  END SUBROUTINE RhoEst
!--------------------------------------------------------------
! Sort the Estimates per Exponent
!--------------------------------------------------------------
  SUBROUTINE SortRhoEst(zq,NQ,LenKet,Rho)
    INTEGER                           :: zq,NQ,LenKet
    TYPE(HGRho)                       :: Rho
    INTEGER                           :: oq,or,iq,qadd,radd,LMN,add,sadd,isort
    INTEGER,DIMENSION(NQ)             :: IntRhoEst
    REAL(DOUBLE),DIMENSION(NQ)        :: RhoEst,QRx,QRy,QRz
    REAL(DOUBLE),DIMENSION(NQ*LenKet) :: RhoCo
!
    oq     = Rho%OffQ%I(zq)
    or     = Rho%OffR%I(zq)
!
    DO iq = 1,NQ
       qadd = oq+iq
       IntRhoEst(iq)= iq
       RhoEst(iq)   = Rho%Est%D(qadd)
    ENDDO
    CALL DblIntSort77(NQ,RhoEst,IntRhoEst,-2)
!
!   Reorder RhoEst QRx,QRy, and QRz and RhoCo
!
    DO iq = 1,NQ
       add  = (iq-1)*LenKet
       qadd = oq+iq
       radd = or+add
       RhoEst(iq)  = Rho%Est%D(qadd)
       QRx(iq)     = Rho%Qx%D(qadd)
       QRy(iq)     = Rho%Qy%D(qadd)
       QRz(iq)     = Rho%Qz%D(qadd)
       DO LMN = 1,LenKet
          RhoCo(add+LMN) = Rho%Co%D(radd+LMN)
       ENDDO
    ENDDO
    DO iq = 1,NQ
       add   = (iq-1)*LenKet
       qadd  = oq+iq
       radd  = or+add
       isort = IntRhoEst(iq)
       sadd  = (isort-1)*LenKet
       Rho%Est%D(qadd) = RhoEst(isort)
       Rho%Qx%D(qadd)    = QRx(isort)
       Rho%Qy%D(qadd)    = QRy(isort)
       Rho%Qz%D(qadd)    = QRz(isort)
       DO LMN = 1,LenKet
          Rho%Co%D(radd+LMN) = RhoCo(sadd+LMN)
       ENDDO
    ENDDO
!
  END SUBROUTINE SortRhoEst
!
!
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
