MODULE BraKetBloks
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
!--------------------------------------------------------------
! Global Objects
!--------------------------------------------------------------
  TYPE(DBL_RNK3)  :: HGBra,dHGBradX,dHGBradY,dHGBradZ
  TYPE(DBL_VECT)  :: HGKet
  TYPE(DBL_VECT)  :: PhFactor,Phase
  TYPE(DBL_VECT)  :: AuxRfun
  TYPE(DBL_RNK3)  :: Rfun
  CONTAINS
!--------------------------------------------------------------
! Allocate global space for computation of BraKet coefficients
! and their corresponding integral estimates
!
  SUBROUTINE NewBraKetBlok(BS)
    TYPE(BSET) :: BS
    INTEGER    :: K,NBF,MaxL,MaxL2,MaxBF
!   Find the bigest 
    MaxBF=0
    DO K=1,BS%NKind
       NBF=BS%BFKnd%I(K)
       MaxBF=MAX(MaxBF,NBF)
    ENDDO
!   Max L on a product function
    MaxL=2*BS%NAsym
    MaxLMN=LHGTF(MaxL)
!   Max L for a two electron integral
    MaxL2=2*(MaxL+1)
!   Allocate space for computing integral estimates
!   
!   OBSOLEATE WITH QCTC, SHOULD REALLY BE ALLOCATING GLOBAL MD INSTEAD!
!
    CALL New(AuxRfun,MaxL2,-1)
    CALL New(Rfun,(/MaxL2,MaxL2,MaxL2/),(/-1,-1,-1/))
!   Allocate space for BraKet contraction
    CALL New(HGKet,MaxLMN)
    CALL New(HGBra,(/MaxLMN,MaxBF,MaxBF/))
!   Allocate and set intermediate phase factors
    CALL New(PhFactor,MaxL,0)
    DO L=0,MaxL; PhFactor%D(L)=(-One)**L; ENDDO
    CALL New(Phase,MaxLMN)
    DO L=0,MaxL
    DO M=0,MaxL-L
    DO N=0,MaxL-L-M
       LMN=LMNDex(L,M,N)
       Phase%D(LMN)=(-One)**(L+M+N)
    ENDDO; ENDDO; ENDDO
  END SUBROUTINE NewBraKetBlok
!
  SUBROUTINE NewGradBraKetBlok(BS)
    TYPE(BSET) :: BS
    INTEGER    :: K,NBF,MaxL,MaxL2,MaxBF
!   Find the biggest 
    MaxBF=0
    DO K=1,BS%NKind
       NBF=BS%BFKnd%I(K)
       MaxBF=MAX(MaxBF,NBF)
    ENDDO
!   Max L on a product function
    MaxL=2*BS%NAsym+1
    MaxLMN=LHGTF(MaxL+1)
!   Max L for a two electron integral
    MaxL2=2*(MaxL+1)
!   Allocate space for computing integral estimates
    CALL New(AuxRfun,MaxL2,-1)
    CALL New(Rfun,(/MaxL2,MaxL2,MaxL2/),(/-1,-1,-1/))
!   Allocate space for BraKet contraction
    CALL New(HGKet,MaxLMN)
    CALL New(HGBra,(/MaxLMN,MaxBF,MaxBF/))
    CALL New(dHGBradX,(/MaxLMN,MaxBF,MaxBF/))
    CALL New(dHGBradY,(/MaxLMN,MaxBF,MaxBF/))
    CALL New(dHGBradZ,(/MaxLMN,MaxBF,MaxBF/))
!   Allocate and set intermediate phase factor
    CALL New(PhFactor,MaxL,0)
    DO L=0,MaxL; PhFactor%D(L)=(-One)**L; ENDDO
    CALL New(Phase,MaxLMN)
    DO L=0,MaxL
    DO M=0,MaxL-L
    DO N=0,MaxL-L-M
       LMN=LMNDex(L,M,N)
       Phase%D(LMN)=(-One)**(L+M+N)
    ENDDO; ENDDO; ENDDO
  ENDSUBROUTINE NewGradBraKetBlok
!
!--------------------------------------------------------------
! Delete the Global Objects 
!--------------------------------------------------------------
  SUBROUTINE DeleteBraKetBlok()
    CALL Delete(AuxRfun)
    CALL Delete(Rfun)
    CALL Delete(HGKet)
    CALL Delete(HGBra)
    CALL Delete(PhFactor)
  END SUBROUTINE DeleteBraKetBlok
  SUBROUTINE DeleteGradBraKetBlok()
    CALL Delete(AuxRfun)
    CALL Delete(Rfun)
    CALL Delete(HGKet)
    CALL Delete(HGBra)
    CALL Delete(dHGBradX)
    CALL Delete(dHGBradY)
    CALL Delete(dHGBradZ)
    CALL Delete(PhFactor)
  ENDSUBROUTINE DeleteGradBraKetBlok
!
!--------------------------------------------------------------
! Generate The Primative
!--------------------------------------------------------------
  SUBROUTINE SetBraBlok(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Print_O,Phase_O)
    INTEGER                  :: CFA,PFA,KA,CFB,PFB,KB
    INTEGER                  :: LMNA,LA,MA,NA,LMNB,LB,MB,NB
    INTEGER                  :: IA,IB,LAB,MAB,NAB,LMN
    INTEGER                  :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB
    REAL(DOUBLE)             :: ExpAB,Ex,Exy,Exyz,CA,CB
    TYPE(BSET)               :: BS
    TYPE(DBL_RNK4)           :: MD
    LOGICAL, OPTIONAL        :: Phase_O,Print_O
!
    HGBra%D = Zero
    IndexA  = CFBlokDex(BS,CFA,KA)
    IndexB  = CFBlokDex(BS,CFB,KB)
    StartLA = BS%LStrt%I(CFA,KA)        
    StopLA  = BS%LStop%I(CFA,KA)
    StartLB = BS%LStrt%I(CFB,KB)        
    StopLB  = BS%LStop%I(CFB,KB)
!
    IF(PRESENT(Phase_O))THEN
       IA = IndexA
       DO LMNA=StartLA,StopLA
          IA=IA+1
          IB=IndexB
          LA=BS%LxDex%I(LMNA)
          MA=BS%LyDex%I(LMNA) 
          NA=BS%LzDex%I(LMNA)
          CA=ExpAB*BS%CCoef%D(LMNA,PFA,CFA,KA)
          DO LMNB=StartLB,StopLB
             IB=IB+1
             LB=BS%LxDex%I(LMNB)
             MB=BS%LyDex%I(LMNB)
             NB=BS%LzDex%I(LMNB)
             CAB = CA*BS%CCoef%D(LMNB,PFB,CFB,KB)
             DO LAB=0,LA+LB
                Ex = CAB*MD%D(1,LA,LB,LAB)
                DO MAB=0,MA+MB
                   Exy = Ex*MD%D(2,MA,MB,MAB)
                   DO NAB=0,NA+NB
                      LMN = LMNdex(LAB,MAB,NAB)
                      Exyz = Exy*MD%D(3,NA,NB,NAB)*PhFactor%D(LAB+MAB+NAB)
                      HGBra%D(LMN,IA,IB) = HGBra%D(LMN,IA,IB)+Exyz
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ELSE
       IA = IndexA
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
             DO LAB=0,LA+LB
                DO MAB=0,MA+MB
                   DO NAB=0,NA+NB
                      LMN = LMNdex(LAB,MAB,NAB)
                      HGBra%D(LMN,IA,IB)=HGBra%D(LMN,IA,IB)   &
                                        +CA*CB*ExpAB          &
                                        *MD%D(1,LA,LB,LAB)    &
                                        *MD%D(2,MA,MB,MAB)    &
                                        *MD%D(3,NA,NB,NAB)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDIF
!
  END SUBROUTINE SetBraBlok
!--------------------------------------------------------------
! Generate The Estimate
!--------------------------------------------------------------
  FUNCTION Estimate(LKet,LenKet,Coef)
    INTEGER                        :: LKet,LenKet,LKet2
    INTEGER                        :: I,J,K,L,M,N,IJK,LMN,IL,JM,KN
    REAL(DOUBLE)                   :: Coef1,Coef2,Estimate
    REAL(DOUBLE),DIMENSION(LenKet) :: Coef
!
    LKet2    = LKet
    Estimate = Zero
    DO K=0,LKet2
       DO J=0,LKet2-K
          DO I=0,LKet2-K-J
             IJK=LMNdex(I,J,K)
             Coef1 = ((-One)**(I+J+K))*Coef(IJK)
             DO N=0,LKet2
                KN = K+N
                DO M=0,LKet2-N
                   JM = J+M
                   DO L=0,LKet2-N-M
                      IL = I+L
                      LMN=LMNdex(L,M,N)
                      Coef2 = Coef(LMN)
                      Estimate = Estimate + Coef1*Coef2*Rfun%D(IL,JM,KN)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    Estimate = DSQRT(DABS(Estimate))
!    WRITE(*,*) Estimate
!
  END FUNCTION Estimate
!--------------------------------------------------------------
! Generate R and AuxRfun
!--------------------------------------------------------------
  SUBROUTINE GenerateRfun(Zeta,LKet)
    INTEGER            :: I0,I1,I2,L,LKet,LKet2
    REAL(DOUBLE)       :: Zeta,Omega,o1,o2
!
    Rfun%D    = Zero
    AuxRfun%D = Zero
!
    Omega  = Half*Zeta
    LKet2  = 2*LKet
    o1     = One
    o2     = -Two*Omega
    DO L=0,LKet2
       AuxRfun%D(L)=o1/DBLE(2*L+1)
       o1=o1*o2
    ENDDO
    DO L=LKet2,0,-1
       DO I0=LKet2-L,1,-1
          Rfun%D(I0,0,0)=DBLE(I0-1)*Rfun%D(I0-2,0,0)
          DO I1=LKet2-L-I0,1,-1
             Rfun%D(I1,I0,0)=DBLE(I1-1)*Rfun%D(I1-2,I0,0)
             DO I2=LKet2-L-I0-I1,1,-1
                Rfun%D(I2,I1,I0)=DBLE(I2-1)*Rfun%D(I2-2,I1,I0)
             ENDDO
             Rfun%D(0,I1,I0)=DBLE(I1-1)*Rfun%D(0,I1-2,I0)
             Rfun%D(I1,0,I0)=DBLE(I1-1)*Rfun%D(I1-2,0,I0)
          ENDDO
          Rfun%D(0,I0,0)=DBLE(I0-1)*Rfun%D(0,I0-2,0)
          Rfun%D(0,0,I0)=DBLE(I0-1)*Rfun%D(0,0,I0-2)
       ENDDO
       Rfun%D(0,0,0)=AuxRfun%D(L)
    ENDDO
!
  END SUBROUTINE GenerateRfun
!
!
!--------------------------------------------------------------
! Determine if the Atom is Being Differentiated wrt
!--------------------------------------------------------------
  FUNCTION OnAtom(I,DiffAtm)
    INTEGER                 :: I,DiffAtm
    REAL(DOUBLE)            :: OnAtom
!
    IF (DiffAtm==I) THEN
       OnAtom=One
    ELSE
       OnAtom=Zero
    ENDIF
  END FUNCTION OnAtom
!
!--------------------------------------------------------------
! Compute the Gradient of SAB_F
! Eq. 2.47 of McMurchies Thesis, pg. 15.  The notation
! is the derivative of SAB wrt F.
!--------------------------------------------------------------
  FUNCTION GradSAB_F(D,DA,DB,E,EA,EB,F,FA,FB,MD,ZetaA,ZetaB,OnA,OnB)
    INTEGER                 :: D,E,F,DA,DB,EA,EB,FA,FB
    REAL(DOUBLE)            :: GradSAB_F,RFA,RFB,ZetaA,ZetaB,OnA,OnB
    TYPE(DBL_RNK4)          :: MD

    RFA=DBLE(FA)
    RFB=DBLE(FB)

    GradSAB_F=MD%D(D,DA,DB,0)*MD%D(E,EA,EB,0)*( &
         OnA*(Two*ZetaA*MD%D(F,FA+1,FB,0) - RFA*MD%D(F,FA-1,FB,0)) + &
         OnB*(Two*ZetaB*MD%D(F,FA,FB+1,0) - RFB*MD%D(F,FA,FB-1,0)))
  END FUNCTION GradSAB_F
!
!--------------------------------------------------------------
! Compute the Gradient of DelPhiA.DelPhiB
! Derivative of Eq. 2.53 of McMurchies Thesis, pg. 17.  The notation
! is the derivative of I_E wrt F.
!--------------------------------------------------------------
  FUNCTION GradI_EF(D,DA,DB,E,EA,EB,F,FA,FB,MD,ZetaA,ZetaB,OnA,OnB)
    INTEGER                 :: D,E,F,DA,DB,EA,EB,FA,FB
    REAL(DOUBLE)            :: GradI_EF,REA,REB,RFA,RFB,ZetaA,ZetaB,OnA,OnB
    TYPE(DBL_RNK4)          :: MD

    REA=DBLE(EA)
    REB=DBLE(EB)
    RFA=DBLE(FA)
    RFB=DBLE(FB)

    GradI_EF=OnA*(MD%D(D,DA,DB,0)*(REA*REB*MD%D(E,EA-1,EB-1,0) - &
         Two*ZetaB*REA*MD%D(E,EA-1,EB+1,0) - Two*ZetaA*REB*MD%D(E,EA+1,EB-1,0) + &
         Four*ZetaA*ZetaB*MD%D(E,EA+1,EB+1,0))*( &
         Two*ZetaA*MD%D(F,FA+1,FB,0)-RFA*MD%D(F,FA-1,FB,0))) + &
         OnB*(MD%D(D,DA,DB,0)*(REA*REB*MD%D(E,EA-1,EB-1,0) - &
         Two*ZetaB*REA*MD%D(E,EA-1,EB+1,0) - Two*ZetaA*REB*MD%D(E,EA+1,EB-1,0) + &
         Four*ZetaA*ZetaB*MD%D(E,EA+1,EB+1,0))*( &
         Two*ZetaB*MD%D(F,FA,FB+1,0)-RFB*MD%D(F,FA,FB-1,0)))

  END FUNCTION GradI_EF
!
!--------------------------------------------------------------
! Compute the Gradient of DelPhiA.DelPhiB
! Derivative of Eq. 2.53 of McMurchies Thesis, pg. 17.  The notation
! is the derivative of I_F wrt F.
!--------------------------------------------------------------
  FUNCTION GradI_FF(D,DA,DB,E,EA,EB,F,FA,FB,MD,ZetaA,ZetaB,OnA,OnB)
    INTEGER                 :: D,E,F,DA,DB,EA,EB,FA,FB
    REAL(DOUBLE)            :: GradI_FF,RFA,RFB,ZetaA,ZetaB,OnA,OnB
    TYPE(DBL_RNK4)          :: MD

    RFA=DBLE(FA)
    RFB=DBLE(FB)

    GradI_FF=MD%D(D,DA,DB,0)*MD%D(E,EA,EB,0)*( &
         OnA*(Two*ZetaA*(RFA*RFB*MD%D(F,FA,FB-1,0) - &
         Two*ZetaB*RFA*MD%D(F,FA,FB+1,0) - Two*ZetaA*RFB*MD%D(F,FA+2,FB-1,0) + &
         Four*ZetaA*ZetaB*MD%D(F,FA+2,FB+1,0)) - ( &
         RFA*RFB*(RFA-One)*MD%D(F,FA-2,FB-1,0) - &
         Two*ZetaB*RFA*(RFA-One)*MD%D(F,FA-2,FB+1,0) - &
         Two*ZetaA*RFB*(RFA+One)*MD%D(F,FA,FB-1,0) + &
         Four*ZetaA*ZetaB*(RFA+One)*MD%D(F,FA,FB+1,0))) + &
         OnB*(Two*ZetaB*(RFA*RFB*MD%D(F,FA-1,FB,0) - &
         Two*ZetaB*RFA*MD%D(F,FA-1,FB+2,0) - Two*ZetaA*RFB*MD%D(F,FA+1,FB,0) + &
         Four*ZetaA*ZetaB*MD%D(F,FA+1,FB+2,0)) - ( &
         RFA*RFB*(RFB-One)*MD%D(F,FA-1,FB-2,0) - &
         Two*ZetaB*RFA*(RFB+One)*MD%D(F,FA-1,FB,0) - &
         Two*ZetaA*RFB*(RFB-One)*MD%D(F,FA+1,FB-2,0) + &
         Four*ZetaA*ZetaB*(RFB+One)*MD%D(F,FA+1,FB,0))))

  END FUNCTION GradI_FF
!
  FUNCTION GradBra(D,DA,DB,DAB,E,EA,EB,EAB,F,FA,FB,FAB,MD,ZetaA,ZetaB,OnA,OnB,PorM)
    INTEGER                 :: D,E,F,DA,DB,EA,EB,FA,FB,DAB,EAB,FAB
    REAL(DOUBLE)            :: GradBra,RFA,RFB,ZetaA,ZetaB,OnA,OnB
    TYPE(DBL_RNK4)          :: MD
    CHARACTER               :: PorM

    RFA=DBLE(FA)
    RFB=DBLE(FB)

    IF (PorM.EQ.'M') THEN
       GradBra=OnA*RFA*MD%D(F,FA-1,FB,FAB)*MD%D(D,DA,DB,DAB)*MD%D(E,EA,EB,EAB) + &
               OnB*RFB*MD%D(F,FA,FB-1,FAB)*MD%D(D,DA,DB,DAB)*MD%D(E,EA,EB,EAB)
    ELSEIF (PorM.EQ.'P') THEN
       GradBra=OnA*Two*ZetaA*MD%D(F,FA+1,FB,FAB)*MD%D(D,DA,DB,DAB)*MD%D(E,EA,EB,EAB) + &
               OnB*Two*ZetaB*MD%D(F,FA,FB+1,FAB)*MD%D(D,DA,DB,DAB)*MD%D(E,EA,EB,EAB)
    ENDIF
  END FUNCTION GradBra
!--------------------------------------------------------------
! Generate The Gradient of the One Electron Pieces
!--------------------------------------------------------------
  SUBROUTINE dOneE(DiffAtm,AtA,AtB,CFA,PFA,KA,CFB,PFB,KB,Ov,BS,MD,OneE, &
       GXP,GYP,GZP,Pair)
    INTEGER                  :: CFA,PFA,KA,CFB,PFB,KB
    INTEGER                  :: LMNA,LA,MA,NA,LMNB,LB,MB,NB
    INTEGER                  :: IA,IB,LAB,MAB,NAB,LMN,DiffAtm,AtA,AtB
    INTEGER                  :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB
    TYPE(AtomPair)           :: Pair
    REAL(DOUBLE)             :: Ov,CA,CB,CAB,CC,RLA,RLB,RMA,RMB,RNA,RNB, &
                                MDxyzAX,MDxyzBX,MDxyzAY,MDxyzBY,MDxyzAZ,MDxyzBZ,     &
                                OnA,OnB,ZetaA,ZetaB,IYX,IXY,IYZ,IZY,IXZ,IZX,IXX,IYY,IZZ
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)        :: GXP,GYP,GZP
    TYPE(BSET),   INTENT(IN)                   :: BS
    TYPE(DBL_RNK4)           :: MD
    CHARACTER                :: OneE
!
    OnA=OnAtom(AtA,DiffAtm)
    OnB=OnAtom(AtB,DiffAtm)
    IndexA  = CFBlokDex(BS,CFA,KA)
    IndexB  = CFBlokDex(BS,CFB,KB)
    StartLA = BS%LStrt%I(CFA,KA)        
    StopLA  = BS%LStop%I(CFA,KA)
    StartLB = BS%LStrt%I(CFB,KB)        
    StopLB  = BS%LStop%I(CFB,KB)
    ZetaA= BS%Expnt%D(PFA,CFA,KA)
    ZetaB= BS%Expnt%D(PFB,CFB,KB)
!
    IA = IndexA
    DO LMNA=StartLA,StopLA
       IA=IA+1
       IB=IndexB
       LA=BS%LxDex%I(LMNA)
       MA=BS%LyDex%I(LMNA) 
       NA=BS%LzDex%I(LMNA)
       CA=Ov*BS%CCoef%D(LMNA,PFA,CFA,KA)
       DO LMNB=StartLB,StopLB
          IB=IB+1
          LB=BS%LxDex%I(LMNB)
          MB=BS%LyDex%I(LMNB)
          NB=BS%LzDex%I(LMNB)
          CAB = CA*BS%CCoef%D(LMNB,PFB,CFB,KB)

          IF (OneE.EQ.'S') THEN
             GXP(IA,IB)=GXP(IA,IB) + &
                  CAB*GradSAB_F(3,NA,NB,2,MA,MB,1,LA,LB,MD,ZetaA,ZetaB,OnA,OnB)
             GYP(IA,IB)=GYP(IA,IB) + &
                  CAB*GradSAB_F(3,NA,NB,1,LA,LB,2,MA,MB,MD,ZetaA,ZetaB,OnA,OnB)
             GZP(IA,IB)=GZP(IA,IB) + &
                  CAB*GradSAB_F(1,LA,LB,2,MA,MB,3,NA,NB,MD,ZetaA,ZetaB,OnA,OnB)
          ELSEIF (OneE.EQ.'T') THEN
             IYX=GradI_EF(3,NA,NB,2,MA,MB,1,LA,LB,MD,ZetaA,ZetaB,OnA,OnB)
             IXY=GradI_EF(3,NA,NB,1,LA,LB,2,MA,MB,MD,ZetaA,ZetaB,OnA,OnB)
             IYZ=GradI_EF(1,LA,LB,2,MA,MB,3,NA,NB,MD,ZetaA,ZetaB,OnA,OnB)
             IZY=GradI_EF(1,LA,LB,3,NA,NB,2,MA,MB,MD,ZetaA,ZetaB,OnA,OnB)
             IXZ=GradI_EF(2,MA,MB,1,LA,LB,3,NA,NB,MD,ZetaA,ZetaB,OnA,OnB)
             IZX=GradI_EF(2,MA,MB,3,NA,NB,1,LA,LB,MD,ZetaA,ZetaB,OnA,OnB)

             IXX=GradI_FF(2,MA,MB,3,NA,NB,1,LA,LB,MD,ZetaA,ZetaB,OnA,OnB)
             IYY=GradI_FF(1,LA,LB,3,NA,NB,2,MA,MB,MD,ZetaA,ZetaB,OnA,OnB)
             IZZ=GradI_FF(1,LA,LB,2,MA,MB,3,NA,NB,MD,ZetaA,ZetaB,OnA,OnB)

             GXP(IA,IB) = GXP(IA,IB) + &
                  Half*CAB*(IXX+IYX+IZX)
             GYP(IA,IB) = GYP(IA,IB) + &
                  Half*CAB*(IXY+IYY+IZY)
             GZP(IA,IB) = GZP(IA,IB) + &
                  Half*CAB*(IXZ+IYZ+IZZ)
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE dOneE
!
!--------------------------------------------------------------
! Generate The Gradient of the Two Electron Pieces
!--------------------------------------------------------------
  SUBROUTINE dTwoE(DiffAtm,AtA,AtB,CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Pair)
    TYPE(AtomPair)                            :: Pair
    INTEGER                  :: CFA,PFA,KA,CFB,PFB,KB,DiffAtm,AtA,AtB
    INTEGER                  :: LMNA,LA,MA,NA,LMNB,LB,MB,NB
    INTEGER                  :: IA,IB,LAB,MAB,NAB,LMN
    INTEGER                  :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB
    REAL(DOUBLE)             :: ExpAB,Ex,Exy,Exyz,CA,CB,CAB,RLA,RLB,RMA,RMB,RNA,RNB, &
                                MDxyzAX,MDxyzBX,MDxyzAY,MDxyzBY,MDxyzAZ,MDxyzBZ,     &
                                OnA,OnB,ZetaA,ZetaB
    TYPE(BSET)               :: BS
    TYPE(DBL_RNK4)           :: MD
!
    OnA=OnAtom(AtA,DiffAtm)
    OnB=OnAtom(AtB,DiffAtm)
    dHGBradX%D = Zero
    dHGBradY%D = Zero
    dHGBradZ%D = Zero
    IndexA  = CFBlokDex(BS,CFA,KA)
    IndexB  = CFBlokDex(BS,CFB,KB)
    StartLA = BS%LStrt%I(CFA,KA)        
    StopLA  = BS%LStop%I(CFA,KA)
    StartLB = BS%LStrt%I(CFB,KB)        
    StopLB  = BS%LStop%I(CFB,KB)
    ZetaA= BS%Expnt%D(PFA,CFA,KA)
    ZetaB= BS%Expnt%D(PFB,CFB,KB)
!
    IA = IndexA
    DO LMNA=StartLA,StopLA
       IA=IA+1
       IB=IndexB
       LA=BS%LxDex%I(LMNA)
       MA=BS%LyDex%I(LMNA) 
       NA=BS%LzDex%I(LMNA)
       CA=ExpAB*BS%CCoef%D(LMNA,PFA,CFA,KA)
       DO LMNB=StartLB,StopLB
          IB=IB+1
          LB=BS%LxDex%I(LMNB)
          MB=BS%LyDex%I(LMNB)
          NB=BS%LzDex%I(LMNB)
          CAB = CA*BS%CCoef%D(LMNB,PFB,CFB,KB)
          DO LAB=0,LA+LB+1
             DO MAB=0,MA+MB+1
                DO NAB=0,NA+NB+1
                   IF (LAB+MAB+NAB.LE.LA+LB+MA+MB+NA+NB+1) THEN
                      LMN = LMNdex(LAB,MAB,NAB)                         
                      IF (LAB.LE.LA+LB-1.AND.MAB.LE.MA+MB.AND.NAB.LE.NA+NB) THEN
                         dHGBradX%D(LMN,IA,IB)=dHGBradX%D(LMN,IA,IB) - &
                              CAB*GradBra(2,MA,MB,MAB,3,NA,NB,NAB,1,LA,LB,LAB,MD,ZetaA,ZetaB,OnA,OnB,'M')
                      ENDIF
                      IF (MAB.LE.MA+MB-1.AND.LAB.LE.LA+LB.AND.NAB.LE.NA+NB) THEN
                         dHGBradY%D(LMN,IA,IB)=dHGBradY%D(LMN,IA,IB) - &
                              CAB*GradBra(1,LA,LB,LAB,3,NA,NB,NAB,2,MA,MB,MAB,MD,ZetaA,ZetaB,OnA,OnB,'M')
                      ENDIF
                      IF (NAB.LE.NA+NB-1.AND.LAB.LE.LA+LB.AND.MAB.LE.MA+MB) THEN
                         dHGBradZ%D(LMN,IA,IB)=dHGBradZ%D(LMN,IA,IB) - &
                              CAB*GradBra(1,LA,LB,LAB,2,MA,MB,MAB,3,NA,NB,NAB,MD,ZetaA,ZetaB,OnA,OnB,'M')
                      ENDIF
                      IF (MAB.LE.MA+MB.AND.NAB.LE.NA+NB) THEN
                         dHGBradX%D(LMN,IA,IB)=dHGBradX%D(LMN,IA,IB) + &
                              CAB*GradBra(2,MA,MB,MAB,3,NA,NB,NAB,1,LA,LB,LAB,MD,ZetaA,ZetaB,OnA,OnB,'P')
                      ENDIF
                      IF (LAB.LE.LA+LB.AND.NAB.LE.NA+NB) THEN
                         dHGBradY%D(LMN,IA,IB)=dHGBradY%D(LMN,IA,IB) + &
                              CAB*GradBra(1,LA,LB,LAB,3,NA,NB,NAB,2,MA,MB,MAB,MD,ZetaA,ZetaB,OnA,OnB,'P')
                      ENDIF
                      IF (MAB.LE.MA+MB.AND.LAB.LE.LA+LB) THEN
                         dHGBradZ%D(LMN,IA,IB)=dHGBradZ%D(LMN,IA,IB) + &
                              CAB*GradBra(1,LA,LB,LAB,2,MA,MB,MAB,3,NA,NB,NAB,MD,ZetaA,ZetaB,OnA,OnB,'P')
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
  END SUBROUTINE dTwoE
!
END MODULE BraKetBloks
