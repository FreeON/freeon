MODULE BraBloks
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
!--------------------------------------------------------------
! Globals
!--------------------------------------------------------------
  TYPE(DBL_RNK3)  :: HGBra
  TYPE(DBL_RNK4)  :: dHGBra
  TYPE(DBL_RNK4)  :: E
  TYPE(DBL_VECT)  :: PhFactor,Phase
  INTEGER         :: MaxL,MaxL2,MaxLMN,MaxBF
  CONTAINS
!--------------------------------------------------------------
! Allocate global space for computation of Bra coefficients
! and their corresponding integral estimates
!
  SUBROUTINE NewBraBlok(BS,Gradients_O)
    TYPE(BSET) :: BS
    INTEGER    :: K,NBF,L,M,N,LMN
    LOGICAL, OPTIONAL :: Gradients_O
!----------------------------------------------------------------
!   Find the biggest bf block
    MaxBF=0
    DO K=1,BS%NKind
       NBF=BS%BFKnd%I(K)
       MaxBF=MAX(MaxBF,NBF)
    ENDDO
!   Max L on a product function
    MaxL=2*BS%NAsym+1
    MaxLMN=LHGTF(MaxL+1)
!   Allocate space for Bra contractions coefficients and their gradients
    IF(PRESENT(Gradients_O))THEN
!      Allocate space for McMurchie Davidson E coefficients
       CALL New(E,(/3,BS%NASym+1,BS%NASym+1,2*BS%NASym+2/),(/1,-1,-1,-1/))
       CALL New(dHGBra,(/MaxLMN,MaxBF,MaxBF,3/))
    ELSE
!      Allocate space for McMurchie Davidson E coefficients
       CALL New(E,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
       CALL New(HGBra,(/MaxLMN,MaxBF,MaxBF/))
    ENDIF
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
  ENDSUBROUTINE NewBraBlok
!--------------------------------------------------------------
! Delete the Global Objects 
!--------------------------------------------------------------
  SUBROUTINE DeleteBraBlok(Gradients_O)
    LOGICAL, OPTIONAL :: Gradients_O
    IF(PRESENT(Gradients_O))THEN
      CALL Delete(dHGBra)
    ELSE
      CALL Delete(HGBra)
    ENDIF
    CALL Delete(PhFactor)
  END SUBROUTINE DeleteBraBlok
!--------------------------------------------------------------
! Generate a primative bra blok
!--------------------------------------------------------------
  FUNCTION SetBraBlok(Prim,BS,SameAtom_O,Print_O,Gradients_O,K_O,KK_O) RESULT(MaxAmp)
    TYPE(PrimPair)    :: Prim
    TYPE(BSET)        :: BS
    LOGICAL, OPTIONAL :: Print_O,SameAtom_O,Gradients_O
    INTEGER,OPTIONAL  :: K_O,KK_O
    REAL(DOUBLE),DIMENSION(3) :: PA,PB
    INTEGER           :: CFA,PFA,KA,CFB,PFB,KB
    INTEGER           :: LMNA,LA,MA,NA,LMNB,LB,MB,NB,L
    INTEGER           :: IA,IB,LAB,MAB,NAB,LMN,MaxLA,MaxLB
    INTEGER           :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB,K
    REAL(DOUBLE)      :: ZA,ZB,Zeta,Xi,ExpAB,CA,CB,CAB,Amp2,MaxAmp,Fx,Fy,Fz
!---------------------------------------------------------------------------------
    KA=Prim%KA
    KB=Prim%KB
    CFA=Prim%CFA
    CFB=Prim%CFB
    PFA=Prim%PFA
    PFB=Prim%PFB
    ZA=Prim%ZA
    ZB=Prim%ZB
    Zeta=Prim%Zeta
    Xi=Prim%Xi
    ExpAB=EXP(-Xi*Prim%AB2)
    Prim%P=(ZA*Prim%A+ZB*Prim%B)/Prim%Zeta
    PA=Prim%P-Prim%A
    PB=Prim%P-Prim%B
    IndexA  = CFBlokDex(BS,CFA,KA)
    IndexB  = CFBlokDex(BS,CFB,KB)
    StartLA = BS%LStrt%I(CFA,KA)        
    StopLA  = BS%LStop%I(CFA,KA)
    StartLB = BS%LStrt%I(CFB,KB)        
    StopLB  = BS%LStop%I(CFB,KB)
    MaxLA=BS%ASymm%I(2,CFA,KA)
    MaxLB=BS%ASymm%I(2,CFB,KB)
    IF(PRESENT(Gradients_O))THEN
!=================================================================================================
!      Compute McMurchie Davidson E coefficients for HG Primitives
       CALL MD2TRR(BS%NASym+1,-1,MaxLA+1,MaxLB,Zeta,E%D,  &
                   PA(1),PB(1),PA(2),PB(2),PA(3),PB(3)) 
!
       MaxAmp   = Zero
!
       L=LHGTF(MaxLA+MaxLB+1)
       IA = IndexA
       DO LMNA=StartLA,StopLA
          IA=IA+1
          IB=IndexB
          DO LMNB=StartLB,StopLB
             IB=IB+1
             DO K=1,3
                CALL DBL_VECT_EQ_DBL_SCLR(L,dHGBra%D(1:L,IA,IB,K),Zero)
             ENDDO
          ENDDO
       ENDDO
!
       IA=IndexA
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
             CAB=CA*BS%CCoef%D(LMNB,PFB,CFB,KB)
#ifdef EXTREME_DEBUG
!-----------------DEBUG---DEBUG---DEBUG---DEBUG---DEBUG---DEBUG---DEBUG----------------
 WRITE(77,22)LMNA,LMNB,LA,MA,NA,LB,MB,NB,  &
            BS%CCoef%D(LMNA,PFA,CFA,KA),BS%CCoef%D(LMNB,PFB,CFB,KB), &
            Prim%ZA,Prim%ZB,Prim%A,Prim%B
22 FORMAT('test[',I3,',',I3,']={ la->',I2,',ma->',I2,',na->',I2,           &
                               ',lb->',I2,',mb->',I2,',nb->',I2,           &
                               ',ca->',F10.7,',cb->',F10.7,               &
                               ',za->',F10.7,',zb->',F10.7,               &
                               ',ax->',F10.7,',ay->',F10.7,',az->',F10.7, &
                               ',bx->',F10.7,',by->',F10.7,',bz->',F10.7,'}; \n')
#endif
!---------------------------------------------------------------------------------------
!            Naive gradients formulae:  See Helgaker and Taylor, TCA 83, p177 (1992), 
!            Eqs (14), (15), (18) and comments after Eq.(26): Gradients_O==Pair%SameAtom
!
             DO LAB=0,LA+LB+1
             DO MAB=0,MA+MB
             DO NAB=0,NA+NB
                LMN=LMNdex(LAB,MAB,NAB)
                IF(Gradients_O)THEN
                   Fx=E%D(1,LA,LB,LAB-1)
                ELSE
                   Fx=Two*ZA*E%D(1,LA+1,LB,LAB)-DBLE(LA)*E%D(1,LA-1,LB,LAB)
                ENDIF
                dHGBra%D(LMN,IA,IB,1)=dHGBra%D(LMN,IA,IB,1)+CAB*Fx*E%D(2,MA,MB,MAB)*E%D(3,NA,NB,NAB)
             ENDDO;ENDDO;ENDDO;
             DO LAB=0,LA+LB
             DO MAB=0,MA+MB+1
             DO NAB=0,NA+NB
                LMN=LMNdex(LAB,MAB,NAB)
                IF(Gradients_O)THEN
                   Fy=E%D(2,MA,MB,MAB-1)
                ELSE
                   Fy=Two*ZA*E%D(2,MA+1,MB,MAB)-DBLE(MA)*E%D(2,MA-1,MB,MAB)
                ENDIF
                dHGBra%D(LMN,IA,IB,2)=dHGBra%D(LMN,IA,IB,2)+CAB*Fy*E%D(1,LA,LB,LAB)*E%D(3,NA,NB,NAB)
          ENDDO;ENDDO;ENDDO;
             DO LAB=0,LA+LB
             DO MAB=0,MA+MB
             DO NAB=0,NA+NB+1
                LMN=LMNdex(LAB,MAB,NAB)
                IF(Gradients_O)THEN
                   Fz=E%D(3,NA,NB,NAB-1)
                ELSE
                   Fz=Two*ZA*E%D(3,NA+1,NB,NAB)-DBLE(NA)*E%D(3,NA-1,NB,NAB)
                ENDIF
                dHGBra%D(LMN,IA,IB,3)=dHGBra%D(LMN,IA,IB,3)+CAB*Fz*E%D(1,LA,LB,LAB)*E%D(2,MA,MB,MAB)
             ENDDO;ENDDO;ENDDO;
          ENDDO
       ENDDO
!--------------------------------------------------------------------------------------------     
       Amp2=Zero
       DO K=1,3
          DO LMN=1,LHGTF(MaxLA+MaxLB+1)
             Amp2=Amp2+dHGBra%D(LMN,IA,IB,K)**2
          ENDDO
          MaxAmp=MAX(MaxAmp,SQRT(Amp2))
       ENDDO
    ELSE
!=================================================================================================
!      Compute McMurchie Davidson E coefficients for HG Primitives
       CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,Zeta,E%D,PA(1),PB(1),PA(2),PB(2),PA(3),PB(3)) 
!
       L=LHGTF(MaxLA+MaxLB)
       IA = IndexA
       DO LMNA=StartLA,StopLA
          IA=IA+1
          IB=IndexB
          DO LMNB=StartLB,StopLB
             IB=IB+1
             CALL DBL_VECT_EQ_DBL_SCLR(L,HGBra%D(1:L,IA,IB),Zero)
          ENDDO
       ENDDO
!
       IA = IndexA
       MaxAmp=Zero
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
                DO MAB=0,MA+MB
                   DO NAB=0,NA+NB
                      LMN=LMNdex(LAB,MAB,NAB)
                      HGBra%D(LMN,IA,IB)=HGBra%D(LMN,IA,IB)      &
                                        + CAB*E%D(1,LA,LB,LAB)   &
                                             *E%D(2,MA,MB,MAB)   &
                                             *E%D(3,NA,NB,NAB)
                   ENDDO
                ENDDO
             ENDDO
!--------------------------------------------------------------------------------------------     
             Amp2=Zero
             DO LMN=1,LHGTF(MaxLA+MaxLB)
                 Amp2=Amp2+HGBra%D(LMN,IA,IB)**2
             ENDDO
             MaxAmp=MAX(MaxAmp,SQRT(Amp2))
          ENDDO
       ENDDO
    ENDIF
!
  END FUNCTION SetBraBlok
!
END MODULE BraBloks
