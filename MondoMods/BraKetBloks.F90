MODULE BraKetBloks
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
!--------------------------------------------------------------
! Global Objects
!--------------------------------------------------------------
  TYPE(DBL_VECT)              :: AuxRfun
  TYPE(DBL_RNK3)              :: Rfun
  TYPE(DBL_VECT)              :: HGKet
  TYPE(DBL_RNK3)              :: HGBra
  CONTAINS
!--------------------------------------------------------------
! Setup and Allocate the Global Objects 
!--------------------------------------------------------------
  SUBROUTINE NewBlock(KA,KB,NBFA,NBFB,BS)
    TYPE(BSET)           :: BS
    INTEGER              :: KA,KB,NBFA,NBFB
    INTEGER              :: MaxLMN,MaxL,MaxL2,MaxLA,MaxLB,CFA,CFB
!
    MaxL = 0
    DO CFA=1,BS%NCFnc%I(KA)                      
       MaxLA=BS%ASymm%I(2,CFA,KA) 
       DO CFB=1,BS%NCFnc%I(KB)                   
          MaxLB=BS%ASymm%I(2,CFB,KB)
          MaxL = MAX(MaxL,MaxLA+MaxLB)
       ENDDO
    ENDDO
    MaxL2  = 2*MaxL+2
    MaxLMN = LHGTF(MaxL)
!
    CALL New(AuxRfun,MaxL2,-1)
    CALL New(Rfun,(/MaxL2,MaxL2,MaxL2/),(/-1,-1,-1/))
    CALL New(HGKet,MaxLMN)
    CALL New(HGBra,(/MaxLMN,NBFA,NBFB/))
!
  END SUBROUTINE NewBlock
!--------------------------------------------------------------
! Delete the Global Objects 
!--------------------------------------------------------------
  SUBROUTINE DeleteBlock()
    CALL Delete(AuxRfun)
    CALL Delete(Rfun)
    CALL Delete(HGKet)
    CALL Delete(HGBra)
  END SUBROUTINE DeleteBlock
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
!--------------------------------------------------------------
! Generate The Primative
!--------------------------------------------------------------
  SUBROUTINE GetPrimative(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD)
    INTEGER                  :: CFA,PFA,KA,CFB,PFB,KB
    INTEGER                  :: LMNA,LA,MA,NA,LMNB,LB,MB,NB
    INTEGER                  :: IA,IB,LAB,MAB,NAB,LMN
    INTEGER                  :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB
    REAL(DOUBLE)             :: ExpAB,Ex,Exy,Exyz,CA,CAB
    TYPE(BSET)               :: BS
    TYPE(DBL_RNK4)           :: MD
!
    HGBra%D = Zero
    IndexA  = IBloDex(BS,CFA,KA)
    IndexB  = IBloDex(BS,CFB,KB)
    StartLA = BS%LStrt%I(CFA,KA)        
    StopLA  = BS%LStop%I(CFA,KA)
    StartLB = BS%LStrt%I(CFB,KB)        
    StopLB  = BS%LStop%I(CFB,KB)
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
          DO LAB=0,LA+LB
             Ex = CAB*MD%D(1,LA,LB,LAB)
             DO MAB=0,MA+MB
                Exy = Ex*MD%D(2,MA,MB,MAB)
                DO NAB=0,NA+NB
                   LMN = LMNdex(LAB,MAB,NAB)
                   Exyz = Exy*MD%D(3,NA,NB,NAB)
                   HGBra%D(LMN,IA,IB) = HGBra%D(LMN,IA,IB)+Exyz
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
  END SUBROUTINE GetPrimative
!
END MODULE BraKetBloks
