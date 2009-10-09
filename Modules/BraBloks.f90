MODULE BraBloks
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE Thresholding
#ifdef MMech
  USE Mechanics
#endif
  !-------------------------------------------------------------------------------
  ! Globals
  !-------------------------------------------------------------------------------
  TYPE(DBL_RNK3)  :: HGBra
  TYPE(DBL_RNK4)  :: E,dE,dD,dHGBra,EcoD
  TYPE(DBL_VECT)  :: PhFactor,Phase,HGSum
  INTEGER         :: MaxL,MaxL2,MaxLMN,MaxBF
CONTAINS
  !-------------------------------------------------------------------------------
  ! Allocate global space for computation of Bra coefficients
  ! and their corresponding integral estimates
  !
  SUBROUTINE NewBraBlok(BS,Gradients_O)
    TYPE(BSET),OPTIONAL :: BS
    INTEGER             :: K,NBF,L,M,N,LMN,N_NAsym
    LOGICAL, OPTIONAL   :: Gradients_O
    !-------------------------------------------------------------------------------
#ifdef MMech
    CALL InitMMech()
#endif
    !   Find the biggest bf block
    MaxBF=0
    IF(PRESENT(BS)) THEN
       DO K=1,BS%NKind
          NBF=BS%BFKnd%I(K)
          MaxBF=MAX(MaxBF,NBF)
       ENDDO
       ! Max L on a product function
       MaxL=2*BS%NAsym+1
       N_NAsym=BS%NAsym
#ifdef MMech
    ELSE IF(MMOnly()) THEN
       MaxBF=1 !maximum of number of basis funct per atm kind
       N_NAsym=0
       MaxL=2*N_NAsym+1
#else
    ELSE
       CALL Halt('No case definition in NewBraBlok')
#endif
    ENDIF
    MaxLMN=LHGTF(MaxL+1)
    ! Allocate space for Bra contractions coefficients and their gradients
    CALL New(dD,(/3,N_NASym+1,N_NASym+1,2*N_NASym+2/),(/1,-1,-1,-1/))
    CALL New(dHGBra,(/MaxLMN,MaxBF,MaxBF,3/))
    ! Allocate space for McMurchie Davidson E coefficients
    CALL New(E,(/3,N_NASym,N_NASym,2*N_NASym/),(/1,0,0,0/))
    CALL New(dE,(/3,N_NASym,N_NASym,2*N_NASym/),(/1,0,0,0/))
    CALL New(HGBra,(/MaxLMN,MaxBF,MaxBF/))
    CALL New(ECoD,(/MaxLMN,MaxBF,MaxBF,3/))
!   Allocate space for summed contributions
    CALL New(HGSum,MaxLMN)
    ! Allocate and set intermediate phase factor
    CALL New(PhFactor,MaxL,0)
    DO L=0,MaxL; PhFactor%D(L)=(-One)**L; ENDDO
       CALL New(Phase,MaxLMN)
       DO L=0,MaxL
          DO M=0,MaxL-L
             DO N=0,MaxL-L-M
                LMN=LMNDex(L,M,N)
                Phase%D(LMN)=(-One)**(L+M+N)
             ENDDO
          ENDDO
       ENDDO
     END SUBROUTINE NewBraBlok
     !-------------------------------------------------------------------------------
     ! Delete the Global Objects
     !-------------------------------------------------------------------------------
     SUBROUTINE DeleteBraBlok(Gradients_O)
       LOGICAL, OPTIONAL :: Gradients_O
       CALL Delete(E)
       CALL Delete(HGSum)
       CALL Delete(dE)
       CALL Delete(dD)
       CALL Delete(EcoD)
       CALL Delete(Phase)
       CALL Delete(HGBra)
       CALL Delete(dHGBra)
       CALL Delete(PhFactor)
     END SUBROUTINE DeleteBraBlok
     !-------------------------------------------------------------------------------
     ! Generate a primative bra blok returning its extent
     !-------------------------------------------------------------------------------
     SUBROUTINE SetBraBlocks(Prim,BS,CompPrim_O,Gradients_O,DerivativeE_O)
       TYPE(PrimPair)    :: Prim
       TYPE(BSET)        :: BS
       LOGICAL, OPTIONAL :: Gradients_O,CompPrim_O,DerivativeE_O
       REAL(DOUBLE)      :: Ext
       REAL(DOUBLE),DIMENSION(3) :: PA,PB
       INTEGER           :: CFA,PFA,KA,CFB,PFB,KB
       INTEGER           :: LMNA,LA,MA,NA,LMNB,LB,MB,NB,Ell,Len
       INTEGER           :: I,IA,IB,LAB,MAB,NAB,LMN,MaxLA,MaxLB
       INTEGER           :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB,K
       REAL(DOUBLE)      :: ZA,ZB,Zeta,Xi,ExpAB,CA,CB,CAB,Amp2,Fx,Fy,Fz
       INTEGER           :: N_NAsym
       !-------------------------------------------------------------------------------
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
       !       WRITE(*,*)" XiAB2 = ",-Xi*Prim%AB2," P = ",Prim%P(1),' OVERLAP = ',ExpAB
       IF(.NOT.Present(CompPrim_O))  &
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
       N_NAsym=BS%NASym
       !
       IF(.NOT.PRESENT(Gradients_O))THEN
          IF(.NOT.PRESENT(DerivativeE_O))THEN
             ! DEFAULT: Compute McMurchie Davidson E coefficients for HG Primitives
             !
             CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,Zeta,E%D,PA(1),PB(1),PA(2),PB(2),PA(3),PB(3))
             !
             Ell=MaxLA+MaxLB
             Len=LHGTF(Ell)
             IA = IndexA
             DO LMNA=StartLA,StopLA
                IA=IA+1
                IB=IndexB
                DO LMNB=StartLB,StopLB
                   IB=IB+1
                   CALL DBL_VECT_EQ_DBL_SCLR(Len,HGBra%D(1,IA,IB),Zero)
                ENDDO
             ENDDO
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
                      DO MAB=0,MA+MB
                         DO NAB=0,NA+NB
                            LMN=LMNdex(LAB,MAB,NAB)
                            HGBra%D(LMN,IA,IB)=HGBra%D(LMN,IA,IB)   &
                                 + CAB*E%D(1,LA,LB,LAB)   &
                                 *E%D(2,MA,MB,MAB)   &
                                 *E%D(3,NA,NB,NAB)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ELSE
             IF(.NOT.DerivativeE_O)THEN
                CALL Halt(' Bad logic in BraBlocks for DerivativeE_O ')
             ELSE
                !
                ! Compute derivatives of the McMurchie Davidson E coefficients with respect to
                ! nuclear center A.   See Helgaker and Taylor, TCA 83, p177 (1992), Eqs (20).
                ! NOTE: <<NOT>> THE SAME AS THE COEFFICIENTS FOR GRADIENTS OF THE
                ! DISTRIBUTIONS, GIVEN BY EQ. (15).
                !
                CALL dMD2TRR(BS%NASym,0,MaxLA,MaxLB,Za,Zb,Zeta,E%D,dE%D,PA(1),PB(1),PA(2),PB(2),PA(3),PB(3))
                !
                Ell=MaxLA+MaxLB
                Len=LHGTF(Ell)
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      CALL DBL_VECT_EQ_DBL_SCLR(Len,HGBra%D(1,IA,IB),Zero)
                      DO K=1,3
                         CALL DBL_VECT_EQ_DBL_SCLR(Len,ECoD%D(1,IA,IB,K),Zero)
                      ENDDO
                   ENDDO
                ENDDO
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
                      CAB=CA*BS%CCoef%D(LMNB,PFB,CFB,KB)
                      DO LAB=0,LA+LB
                         DO MAB=0,MA+MB
                            DO NAB=0,NA+NB
                               LMN=LMNdex(LAB,MAB,NAB)
                               HGBra%D(LMN,IA,IB)=HGBra%D(LMN,IA,IB)   &
                                    + CAB*E%D(1,LA,LB,LAB)   &
                                    *E%D(2,MA,MB,MAB)   &
                                    *E%D(3,NA,NB,NAB)
                               ECoD%D(LMN,IA,IB,1)=ECoD%D(LMN,IA,IB,1) &
                                    + CAB*dE%D(1,LA,LB,LAB)            &
                                    * E%D(2,MA,MB,MAB)                 &
                                    * E%D(3,NA,NB,NAB)
                               ECoD%D(LMN,IA,IB,2)=ECoD%D(LMN,IA,IB,2) &
                                    + CAB* E%D(1,LA,LB,LAB)            &
                                    *dE%D(2,MA,MB,MAB)            &
                                    * E%D(3,NA,NB,NAB)
                               ECoD%D(LMN,IA,IB,3)=ECoD%D(LMN,IA,IB,3) &
                                    + CAB* E%D(1,LA,LB,LAB)            &
                                    * E%D(2,MA,MB,MAB)            &
                                    *dE%D(3,NA,NB,NAB)
                            ENDDO
                         ENDDO
                      ENDDO
                      !                      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                      !                      IF(Ell>=2)HGBra%D(:,:,:)=0D0
                   ENDDO
                ENDDO
             ENDIF
          ENDIF
       ELSE
          ! COEFFICIENTS FOR GRADIENTS OF CARTESIAN GAUSSIAN DISTRIBUTIONS
          ! See Helgaker and Taylor, TCA 83, p177 (1992), Eqs (14), (15), (18)
          ! and comments after Eq.(26).  Gradients_O==Pair%SameAtom
          ! IF GRADIENTS_O IS .TRUE. THEN, GET PARTIAL WITH RESPECT TO HGTF CENTER P.
          ! OTHERWISE, GET PARTIAL WITH RESPECT TO NUCLEAR CENTER A.
          !
          CALL MD2TRR(N_NASym+1,-1,MaxLA+1,MaxLB,Zeta,dD%D,PA(1),PB(1),PA(2),PB(2),PA(3),PB(3))
          !
          Ell=MaxLA+MaxLB+1
          Len=LHGTF(Ell)
          IA = IndexA
          DO LMNA=StartLA,StopLA
             IA=IA+1
             IB=IndexB
             DO LMNB=StartLB,StopLB
                IB=IB+1
                DO K=1,3
                   CALL DBL_VECT_EQ_DBL_SCLR(Len,dHGBra%D(1,IA,IB,K),Zero)
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
                DO LAB=0,LA+LB+1
                   DO MAB=0,MA+MB
                      DO NAB=0,NA+NB
                         LMN=LMNdex(LAB,MAB,NAB)
                         IF(Gradients_O)THEN
                            Fx=dD%D(1,LA,LB,LAB-1)
                         ELSE
                            Fx=Two*ZA*dD%D(1,LA+1,LB,LAB)-DBLE(LA)*dD%D(1,LA-1,LB,LAB)
                         ENDIF
                         dHGBra%D(LMN,IA,IB,1)=dHGBra%D(LMN,IA,IB,1)+CAB*Fx*dD%D(2,MA,MB,MAB)*dD%D(3,NA,NB,NAB)
                      ENDDO
                   ENDDO
                ENDDO
                DO LAB=0,LA+LB
                   DO MAB=0,MA+MB+1
                      DO NAB=0,NA+NB
                         LMN=LMNdex(LAB,MAB,NAB)
                         IF(Gradients_O)THEN
                            Fy=dD%D(2,MA,MB,MAB-1)
                         ELSE
                            Fy=Two*ZA*dD%D(2,MA+1,MB,MAB)-DBLE(MA)*dD%D(2,MA-1,MB,MAB)
                         ENDIF
                         dHGBra%D(LMN,IA,IB,2)=dHGBra%D(LMN,IA,IB,2)+CAB*Fy*dD%D(1,LA,LB,LAB)*dD%D(3,NA,NB,NAB)
                      ENDDO
                   ENDDO
                ENDDO
                DO LAB=0,LA+LB
                   DO MAB=0,MA+MB
                      DO NAB=0,NA+NB+1
                         LMN=LMNdex(LAB,MAB,NAB)
                         IF(Gradients_O)THEN
                            Fz=dD%D(3,NA,NB,NAB-1)
                         ELSE
                            Fz=Two*ZA*dD%D(3,NA+1,NB,NAB)-DBLE(NA)*dD%D(3,NA-1,NB,NAB)
                         ENDIF
                         dHGBra%D(LMN,IA,IB,3)=dHGBra%D(LMN,IA,IB,3)+CAB*Fz*dD%D(1,LA,LB,LAB)*dD%D(2,MA,MB,MAB)
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

       ENDIF
     END SUBROUTINE SetBraBlocks

     FUNCTION SetBraBlok(Prim,BS,SameAtom_O,Gradients_O,Tau_O,ExtraEll_O) RESULT(Ext)
       TYPE(PrimPair)    :: Prim
       TYPE(BSET)        :: BS
       LOGICAL, OPTIONAL :: SameAtom_O,Gradients_O
       REAL(DOUBLE),OPTIONAL:: Tau_O
       INTEGER, OPTIONAL :: ExtraEll_O
       REAL(DOUBLE)      :: Ext
       REAL(DOUBLE),DIMENSION(3) :: PA,PB
       INTEGER           :: CFA,PFA,KA,CFB,PFB,KB
       INTEGER           :: LMNA,LA,MA,NA,LMNB,LB,MB,NB,Ell,Len
       INTEGER           :: I,IA,IB,LAB,MAB,NAB,LMN,MaxLA,MaxLB
       INTEGER           :: IndexA,IndexB,StartLA,StartLB,StopLA,StopLB,K
       REAL(DOUBLE)      :: ZA,ZB,Zeta,Xi,ExpAB,CA,CB,CAB,Amp2,Fx,Fy,Fz
       INTEGER           :: N_NAsym
       !-------------------------------------------------------------------------------
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
       N_NAsym=BS%NASym
       !
       IF(PRESENT(Gradients_O))THEN
          !-------------------------------------------------------------------------------
          !      Compute McMurchie Davidson E coefficients for HG Primitives

          CALL MD2TRR(N_NASym+1,-1,MaxLA+1,MaxLB,Zeta,dD%D,PA(1),PB(1),PA(2),PB(2),PA(3),PB(3))
!!$          CALL MD2TRR(N_NASym+1,-1,MaxLA+1,MaxLB,Zeta,E%D,  &
!!$               PA(1),PB(1),PA(2),PB(2),PA(3),PB(3))
          !
          Ell=MaxLA+MaxLB+1
          Len=LHGTF(Ell)
          IA = IndexA
          DO LMNA=StartLA,StopLA
             IA=IA+1
             IB=IndexB
             DO LMNB=StartLB,StopLB
                IB=IB+1
                DO K=1,3
                   CALL DBL_VECT_EQ_DBL_SCLR(Len,dHGBra%D(1:Len,IA,IB,K),Zero)
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
                !            Naive gradients formulae:  See Helgaker and Taylor, TCA 83, p177 (1992),
                !            Eqs (14), (15), (18) and comments after Eq.(26): Gradients_O==Pair%SameAtom
                DO LAB=0,LA+LB+1
                   DO MAB=0,MA+MB
                      DO NAB=0,NA+NB
                         LMN=LMNdex(LAB,MAB,NAB)
                         IF(Gradients_O)THEN
                            Fx=dD%D(1,LA,LB,LAB-1)
                         ELSE
                            Fx=Two*ZA*dD%D(1,LA+1,LB,LAB)-DBLE(LA)*dD%D(1,LA-1,LB,LAB)
                         ENDIF
                         dHGBra%D(LMN,IA,IB,1)=dHGBra%D(LMN,IA,IB,1)+CAB*Fx*dD%D(2,MA,MB,MAB)*dD%D(3,NA,NB,NAB)
                      ENDDO;ENDDO;ENDDO;
                      DO LAB=0,LA+LB
                         DO MAB=0,MA+MB+1
                            DO NAB=0,NA+NB
                               LMN=LMNdex(LAB,MAB,NAB)
                               IF(Gradients_O)THEN
                                  Fy=dD%D(2,MA,MB,MAB-1)
                               ELSE
                                  Fy=Two*ZA*dD%D(2,MA+1,MB,MAB)-DBLE(MA)*dD%D(2,MA-1,MB,MAB)
                               ENDIF
                               dHGBra%D(LMN,IA,IB,2)=dHGBra%D(LMN,IA,IB,2)+CAB*Fy*dD%D(1,LA,LB,LAB)*dD%D(3,NA,NB,NAB)
                            ENDDO;ENDDO;ENDDO;
                            DO LAB=0,LA+LB
                               DO MAB=0,MA+MB
                                  DO NAB=0,NA+NB+1
                                     LMN=LMNdex(LAB,MAB,NAB)
                                     IF(Gradients_O)THEN
                                        Fz=dD%D(3,NA,NB,NAB-1)
                                     ELSE
                                        Fz=Two*ZA*dD%D(3,NA+1,NB,NAB)-DBLE(NA)*dD%D(3,NA-1,NB,NAB)
                                     ENDIF
                                     dHGBra%D(LMN,IA,IB,3)=dHGBra%D(LMN,IA,IB,3)+CAB*Fz*dD%D(1,LA,LB,LAB)*dD%D(2,MA,MB,MAB)
                                  ENDDO;ENDDO;ENDDO;
                               ENDDO
                            ENDDO
                            ! Compute the extent of this derivative distribution
                            Ext=Zero
                            DO K=1,3
                               CALL DBL_VECT_EQ_DBL_SCLR(Len,HGSum%D(1),Zero)
                               IA = IndexA
                               DO LMNA=StartLA,StopLA
                                  IA=IA+1
                                  IB=IndexB
                                  DO LMNB=StartLB,StopLB
                                     IB=IB+1
                                     DO I=1,Len
                                        HGSum%D(I)=HGSum%D(I)+ABS(dHGBra%D(I,IA,IB,K))
                                     ENDDO
                                  ENDDO
                               ENDDO
                               Ext=MAX(Ext,Extent(Ell,Zeta,HGSum%D(1:Len),Tau_O=Tau_O,ExtraEll_O=ExtraEll_O))
                            ENDDO
                         ELSE
                            !      Compute McMurchie Davidson E coefficients for HG Primitives
                            CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,Zeta,E%D,PA(1),PB(1),PA(2),PB(2),PA(3),PB(3))
                            !
                            Ell=MaxLA+MaxLB
                            Len=LHGTF(Ell)
                            IA = IndexA
                            DO LMNA=StartLA,StopLA
                               IA=IA+1
                               IB=IndexB
                               DO LMNB=StartLB,StopLB
                                  IB=IB+1
                                  CALL DBL_VECT_EQ_DBL_SCLR(Len,HGBra%D(1:Len,IA,IB),Zero)
                               ENDDO
                            ENDDO
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
                                  DO I=1,Len
                                     HGSum%D(I)=HGSum%D(I)+ABS(HGBra%D(I,IA,IB))
                                  ENDDO
                               ENDDO
                            ENDDO
                            Ext=Extent(Ell,Zeta,HGSum%D(1:Len),Tau_O=Tau_O,ExtraEll_O=ExtraEll_O)
                         ENDIF
                       END FUNCTION SetBraBlok

                     END MODULE BraBloks

