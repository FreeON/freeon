MODULE MondoPoles
   USE Derivedtypes
   USE GlobalScalars   
   USE GlobalObjects
   USE ProcessControl
   USE Indexing
   USE Parse
   USE InOut
   USE Macros
   USE Thresholding 
   USE BoundingBox
   USE McMurchie
   IMPLICIT NONE
!=================================================================================================
!  Hierarchical density node
!=================================================================================================
   TYPE PoleNode
      LOGICAL                               :: Leaf     ! Is this a data containing node?
      INTEGER                               :: Bdex     ! Begining index of ORB list for this node
      INTEGER                               :: Edex     ! ENDign index of ORB list for this node
      INTEGER                               :: NQ       ! Number of centers
      INTEGER                               :: Ell      ! Ell type
      REAL(DOUBLE)                          :: Zeta     ! Minimum exponent in this node
      REAL(DOUBLE)                          :: Strength ! Strength of the Pole
      REAL(DOUBLE)                          :: DMax2    ! (Max distance)^2 from node center to DBox 
      TYPE(BBox)                            :: Box      ! Bounding Box of distribution (for PAC)
      TYPE(BBox)                            :: DBox     ! Bounding Box of distribution centers (for MAC)
      TYPE(PoleNode),POINTER                :: Descend  ! Next node in tree descent
      TYPE(PoleNode),POINTER                :: Travrse  ! Next node in tree traversal
#ifdef POINTERS_IN_DERIVED_TYPES
      REAL(DOUBLE),DIMENSION(:),POINTER     :: S        ! Im component of the multipole tensor
      REAL(DOUBLE),DIMENSION(:),POINTER     :: C        ! Re component of the multipole tensor
      REAL(DOUBLE),DIMENSION(:),POINTER     :: Co       ! Coefficients of the HGTF density
#else
      REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: S        ! Im component of the multipole tensor
      REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: C        ! Re component of the multipole tensor
      REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: Co       ! Coefficients of the HGTF density
#endif
   END TYPE                                      
!=====================================================================================================
!
!=====================================================================================================
  INTERFACE HGToSP
     MODULE PROCEDURE HGToSP_PoleNode,HGToSP_Bra
  END INTERFACE
!====================================================================================================
!  Global Array intermediates for multipole computations
!====================================================================================================
   REAL(DOUBLE), DIMENSION(0:2*FFEll2)          :: Factorial
   REAL(DOUBLE), DIMENSION(0:FFEll2)            :: FactOlm0,FactMlm0,Sine,Cosine,CoFact
   REAL(DOUBLE), DIMENSION(0:FFLen2)            :: FactOlm2,FactMlm2,ALegendreP,Spq,Cpq
   REAL(DOUBLE), DIMENSION(0:SPEll+1,0:FFELL)   :: FudgeFactorial
   CONTAINS
!====================================================================================
!     Q->P
!====================================================================================
      SUBROUTINE XLate(Q,P)
         TYPE(PoleNode)            :: Q,P
         REAL(DOUBLE),DIMENSION(3) :: QP
         INTEGER                   :: LP,LQ,LPQ
!------------------------------------------------------------------------------------
         LP=P%Ell
         LQ=Q%Ell
         LPQ=MAX(LP,LQ)
         QP=Q%Box%Center-P%Box%Center
         CALL Regular(LPQ,QP(1),QP(2),QP(3))
         CALL XLate77(LP,LQ,P%C,P%S,Cpq,Spq,Q%C,Q%S)
       END SUBROUTINE XLate
!====================================================================================
!     Contract P with Q
!====================================================================================
      SUBROUTINE CTrax(P,Q,SPKetC,SPKetS)
         TYPE(PrimPair)            :: P
         TYPE(PoleNode)            :: Q
         REAL(DOUBLE),DIMENSION(3) :: PQ
         INTEGER                   :: LP,LQ,LPQ
         REAL(DOUBLE),DIMENSION(0:):: SPKetC,SPKetS
!------------------------------------------------------------------------------------
         LP=P%Ell
         LQ=Q%Ell
         LPQ=LP+LQ        
         PQ=P%P-Q%Box%Center
         CALL IrRegular(LPQ,PQ(1),PQ(2),PQ(3))
         CALL CTraX77(LP,LQ,SPKetC,SPKetS,Cpq,Spq,Q%C,Q%S)
      END SUBROUTINE CTrax
!====================================================================================
!
!====================================================================================
       SUBROUTINE HGToSP_PoleNode(Node)
         TYPE(PoleNode) :: Node
         INTEGER        :: LenHG,LenSP
         REAL(DOUBLE)   :: PiZ
         PiZ=(Pi/Node%Zeta)**(ThreeHalves)
         LenHG=LHGTF(Node%Ell)
         LenSP=LSP(Node%Ell)
         CALL HGToSP_Gen(Node%Ell,PiZ,Node%Co(1:LenHG), &
                         Node%C(0:LenSP),Node%S(0:LenSP))
       END SUBROUTINE HGToSP_PoleNode
!====================================================================================
!
!====================================================================================
       SUBROUTINE HGToSP_Bra(P,HGBra,SPBraC,SPBraS)
         TYPE(PrimPair)                    :: P
         INTEGER                           :: LenHG,LenSP
         REAL(DOUBLE)                      :: PiZ
         REAL(DOUBLE), DIMENSION(1:)       :: HGBra
         REAL(DOUBLE), DIMENSION(0:)       :: SPBraC,SPBraS
!------------------------------------------------------------------------------------
!        Transform <Bra| coefficients from HG to SP
         PiZ=(Pi/P%Zeta)**(ThreeHalves)
         LenHG=LHGTF(P%Ell)
         LenSP=LSP(P%Ell)
         CALL HGToSP_Gen(P%Ell,PiZ,HGBra(1:LenHG), &
                         SPBraC(0:LenSP),SPBraS(0:LenSP))   
       END SUBROUTINE HGToSP_Bra
!====================================================================================
!
!====================================================================================
       SUBROUTINE HGToSP_Gen(L,PiZ,HGCo,C,S) 
          INTEGER                    :: L
          REAL(DOUBLE)               :: PiZ
          REAL(DOUBLE),DIMENSION(1:) :: HGCo      
          REAL(DOUBLE),DIMENSION(0:) :: C,S
          REAL(DOUBLE),DIMENSION(20) :: W
          SELECT CASE(L)
          INCLUDE 'HGToSP.Inc'
          CASE DEFAULT
             CALL Halt('Bad logic in HGToSP_Gen, time to remake HGToSP.Inc')
          END SELECT
       END SUBROUTINE HGToSP_Gen
!====================================================================================
!
!====================================================================================
       SUBROUTINE MultipoleSetUp(Ell)
          INTEGER                          :: Ell,L,M,LDex,LMDex,   &
                                              LP,LQ,MP,MQ
          REAL(DOUBLE)                     :: Sgn,DblFact,TwoTimes, &
                                              DenomP,DenomQ,NumPQ,  &
                                              DegenP,DegenQ
!------------------------------------------------------------------------------------          
!         Factorial(M)=M!
          Factorial(0)=One
          DO L=1,2*Ell 
             Factorial(L)=Factorial(L-1)*DBLE(L)
          ENDDO
!         FactOlm2=
          DO L=2,Ell 
            LDex=L*(L+1)/2
            DO M=0,L-2
               LMDex=LDex+M
               FactOlm2(LMDex)=One/DBLE((L+M)*(L-M))
             ENDDO
          ENDDO
!         FactMlm2=
          DO L=0,Ell
             LDex=L*(L+1)/2
             DO M=0,Ell
                LMDex=LDex+M
                FactMlm2(LMDex)=DBLE((L+M-1)*(L-M-1))
             ENDDO
          ENDDO
!         FactOlm=
!         FactMlm=
          Sgn=One
          DblFact=One
          TwoTimes=One
          DO M=0,Ell
             FactOlm0(M)=Sgn*DblFact/Factorial(2*M)
             FactMlm0(M)=Sgn*DblFact
             DblFact=DblFact*TwoTimes
             TwoTimes=TwoTimes+Two
             Sgn=-Sgn
          ENDDO
!
!         FudgeFactorial(LP,LQ)=[Sum_(MP,MQ) (LP+LQ-MP-MQ)!/((LP+MP)!(LQ+MQ)!)]/{
!         [Sum_MP (LP-MP)!/(LP+MP)!]*[Sum_MQ (LQ-MQ)!/(LQ+MQ)!]}
!
          DO LP=0,SPEll+1
             DO LQ=0,FFEll
                DenomP=One
                DenomQ=One
                DO MP=1,LP
                   DenomP=DenomP+Two*Factorial(LP-MP)/Factorial(LP+MP)
                ENDDO
                DO MQ=1,LQ
                   DenomQ=DenomQ+Two*Factorial(LQ-MQ)/Factorial(LQ+MQ)
                ENDDO
                NumPQ=Zero
                DO MP=0,LP
                   IF(MP==0)THEN
                      DegenP=One
                   ELSE
                      DegenP=Two
                   ENDIF
                   DO MQ=0,LQ
                      IF(MQ==0)THEN
                         DegenQ=One
                      ELSE
                         DegenQ=Two
                      ENDIF
                      NumPQ=NumPQ+DegenP*DegenQ*Factorial(LP+LQ-MP-MQ)/(Factorial(LP+MP)*Factorial(LQ+MQ))
                   ENDDO
                ENDDO
                FudgeFactorial(LP,LQ)=NumPQ/(DenomP*DenomQ)
             ENDDO
          ENDDO
      END SUBROUTINE MultipoleSetup
!====================================================================================
!     Regular Function
!====================================================================================
      SUBROUTINE Regular(Ell,PQx,PQy,PQz)
         INTEGER                   :: Ell
         INTEGER                   :: L,M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
         REAL(DOUBLE)              :: PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
                                      CoFact,PQToThPlsL,PQx,PQy,PQz
!------------------------------------------------------------------------------------
         Cpq=Zero
         Spq=Zero
         PQx2=PQx*PQx
         PQy2=PQy*PQy      
         PQ=SQRT(PQx2+PQy2+PQz*PQz)
         IF(PQ==Zero)THEN
            Cpq(0) = One
            RETURN
         ENDIF
         OneOvPQ=One/PQ
         CoTan=PQz*OneOvPQ
!        Sine and Cosine by recursion
         Cosine(0)=One
         Sine(  0)=Zero
         PQxy=SQRT(PQx2+PQy2)
         IF(PQxy/=Zero)THEN
            Sine(1)=PQx/PQxy
            Cosine(1)=PQy/PQxy
         ELSE
            Sine(1)  =0.70710678118654752D0         
            Cosine(1)=0.70710678118654752D0
         ENDIF
         TwoC=Two*Cosine(1)
         DO M=2,Ell
            M1=M-1
            M2=M-2
            Sine(M)=TwoC*Sine(M1)-Sine(M2)
            Cosine(M)=TwoC*Cosine(M1)-Cosine(M2)
         ENDDO
!        Associated Legendre Polynomials by recursion
         Sq=SQRT(ABS(One-CoTan**2))
         RS=One
         DO M=0,Ell
            MDex=LTD(M)+M
            ALegendreP(MDex)=FactOlm0(M)*RS
            RS=RS*Sq
         ENDDO
         DO M=0,Ell-1
            MDex=LTD(M)+M
            MDex1=LTD(M+1)+M
            ALegendreP(MDex1)=CoTan*ALegendreP(MDex)
         ENDDO
         DO L=2,Ell         
            CoFact=CoTan*DBLE(2*l-1)
            LDex0=LTD(L)
            LDex1=LTD(L-1)
            LDex2=LTD(L-2)
            DO M=0,L-2
               ALegendreP(LDex0+M)=(CoFact*ALegendreP(LDex1+M)-ALegendreP(LDex2+M))*FactOlm2(LDex0+M)
            ENDDO
         ENDDO
!        Regular Spharical Harmonics
         PQToThPlsL=One
         DO L=0,Ell
            LDex=LTD(L)
            DO M=0,L
               LMDex=LDex+M
               Spq(LMDex)=PQToThPlsL*ALegendreP(LMDex)*Sine(M)
               Cpq(LMDex)=PQToThPlsL*ALegendreP(LMDex)*Cosine(M)
            ENDDO
            PQToThPlsL=PQToThPlsL*PQ
         ENDDO
      END SUBROUTINE Regular
!====================================================================================
!     Irregular Function
!====================================================================================
      SUBROUTINE IrRegular(Ell,PQx,PQy,PQz)
         INTEGER                    :: Ell
         INTEGER                    :: L,M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
         REAL(DOUBLE)               :: PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
                                       CoFact,PQToThMnsL,PQx,PQy,PQz
!------------------------------------------------------------------------------------
         Cpq = Zero
         Spq = Zero
         PQx2=PQx*PQx
         PQy2=PQy*PQy      
         PQ=SQRT(PQx2+PQy2+PQz*PQz)
         OneOvPQ=One/PQ
         CoTan=PQz*OneOvPQ
!        Sine and Cosine by recursion
         Cosine(0)=One
         Sine(  0)=Zero
         PQxy=SQRT(PQx2+PQy2)
         IF(PQxy .GT. 1.D-12)THEN
            Sine(1)=PQx/PQxy
            Cosine(1)=PQy/PQxy
         ELSE
            Sine(1)=0.70710678118654752D0         
            Cosine(1)=0.70710678118654752D0
         ENDIF
!
         TwoC=Two*Cosine(1)
         DO M=2,Ell
            M1=M-1
            M2=M-2
            Sine(M)=TwoC*Sine(M1)-Sine(M2)
            Cosine(M)=TwoC*Cosine(M1)-Cosine(M2)
         ENDDO
!        Associated Legendre Polynomials by recursion
         Sq=SQRT(ABS(One-CoTan*CoTan))
!
         RS=One
         DO M=0,Ell
            MDex=LTD(M)+M
            ALegendreP(MDex)=FactMlm0(M)*RS
            RS=RS*Sq
         ENDDO
         DO M=0,Ell-1
            MDex=LTD(M)+M
            MDex1=LTD(M+1)+M
            ALegendreP(MDex1)=CoTan*DBLE(2*M+1)*ALegendreP(MDex)
         ENDDO
         DO L=2,Ell         
            CoFact=CoTan*DBLE(2*l-1)
            LDex0=LTD(L)
            LDex1=LTD(L-1)
            LDex2=LTD(L-2)
            DO M=0,L-2
               ALegendreP(LDex0+M)=CoFact*ALegendreP(LDex1+M)-FactMlm2(LDex0+M)*ALegendreP(LDex2+M)
            ENDDO
         ENDDO
!        IrRegular Spharical Harmonics
         PQToThMnsL=OneOvPQ
         DO L=0,Ell
            LDex=LTD(L)
            DO M=0,l
               LMDex=LDex+M
               Spq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Sine(M)
               Cpq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Cosine(M)
            ENDDO
            PQToThMnsL=PQToThMnsL*OneOvPQ
         ENDDO
      END SUBROUTINE IrRegular
!====================================================================================
!     Irregular Function
!====================================================================================
      SUBROUTINE IrRegularL(Ell,PQx,PQy,PQz)
         INTEGER                    :: Ell
         INTEGER                    :: M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
         REAL(DOUBLE)               :: PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
                                       CoFact,PQToThMnsL,PQx,PQy,PQz
!------------------------------------------------------------------------------------
         Cpq = Zero
         Spq = Zero
         PQx2=PQx*PQx
         PQy2=PQy*PQy      
         PQ=SQRT(PQx2+PQy2+PQz*PQz)
         OneOvPQ=One/PQ
         CoTan=PQz*OneOvPQ
!        Sine and Cosine by recursion
         Cosine(0)=One
         Sine(  0)=Zero
         PQxy=SQRT(PQx2+PQy2)
         IF(PQxy .GT. 1.D-12)THEN
            Sine(1)=PQx/PQxy
            Cosine(1)=PQy/PQxy
         ELSE
            Sine(1)=0.70710678118654752D0         
            Cosine(1)=0.70710678118654752D0
         ENDIF
!
         TwoC=Two*Cosine(1)
         DO M=2,Ell
            M1=M-1
            M2=M-2
            Sine(M)=TwoC*Sine(M1)-Sine(M2)
            Cosine(M)=TwoC*Cosine(M1)-Cosine(M2)
         ENDDO
!        Associated Legendre Polynomials by recursion
         Sq=SQRT(ABS(One-CoTan*CoTan))
!
         RS=One
         DO M=0,Ell
            MDex=LTD(M)+M
            ALegendreP(MDex)=FactMlm0(M)*RS
            RS=RS*Sq
         ENDDO
         DO M=0,Ell-1
            MDex=LTD(M)+M
            MDex1=LTD(M+1)+M
            ALegendreP(MDex1)=CoTan*DBLE(2*M+1)*ALegendreP(MDex)
         ENDDO
!
         CoFact=CoTan*DBLE(2*Ell-1)
         LDex0=LTD(Ell)
         LDex1=LTD(Ell-1)
         LDex2=LTD(Ell-2)
         DO M=0,Ell-2
            ALegendreP(LDex0+M)=CoFact*ALegendreP(LDex1+M)-FactMlm2(LDex0+M)*ALegendreP(LDex2+M)
         ENDDO
!        IrRegular Spharical Harmonics


         PQToThMnsL=OneOvPQ
         LDex=LTD(Ell)
         DO M=0,Ell
            LMDex=LDex+M
            Spq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Sine(M)
            Cpq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Cosine(M)
         ENDDO
         PQToThMnsL=PQToThMnsL*OneOvPQ
!
      END SUBROUTINE IrRegularL
!====================================================================================
!     Compute a multipole strength O_L based on Unsolds theorem
!====================================================================================
      FUNCTION Unsold2(Llow,Lhig,C,S)
        INTEGER                     :: L,Llow,Lhig
        REAL(DOUBLE)                :: Unsold2,U,Exp
        REAL(DOUBLE), DIMENSION(0:) :: C,S
!
        Unsold2=Zero
        DO L=Llow,Lhig
           IF(L==0) THEN 
              U   = Unsold0(L,C,S)
              Unsold2=MAX(Unsold2,U)
           ELSE
              Exp = DBLE(Lhig)/DBLE(L) 
              U   = Unsold0(L,C,S)**Exp
              Unsold2=MAX(Unsold2,U)
           ENDIF
        ENDDO
!
      END FUNCTION Unsold2
!====================================================================================
!     Compute a multipole strength O_L based on Unsolds theorem
!====================================================================================
      FUNCTION Unsold1(L,C,S)
        INTEGER                     :: L,LL
        REAL(DOUBLE)                :: Unsold1,U,Exp
        REAL(DOUBLE), DIMENSION(0:) :: C,S
!
        Unsold1=Zero
        DO LL=0,L
           IF(LL==0) THEN 
              U   = Unsold0(LL,C,S)
              Unsold1=MAX(Unsold1,U)
           ELSE
              Exp = DBLE(L)/DBLE(LL) 
              U   = Unsold0(LL,C,S)**Exp
              Unsold1=MAX(Unsold1,U)
           ENDIF
        ENDDO
!
      END FUNCTION Unsold1
!====================================================================================
!     Compute a multipole strength O_L based on Unsolds theorem
!====================================================================================
      FUNCTION Unsold0(L,C,S)
        INTEGER                     :: I,K,L,M
        REAL(DOUBLE)                :: Unsold0
        REAL(DOUBLE), DIMENSION(0:) :: C,S
!
        K=LTD(L)
        Unsold0=(C(K)**2+S(K)**2)*Factorial(L)**2
        DO M=1,L
           Unsold0=Unsold0+Two*(C(K+M)**2+S(K+M)**2)*Factorial(L+M)*Factorial(L-M)            
        ENDDO
        Unsold0 = SQRT(ABS(Unsold0))
!
      END FUNCTION Unsold0
!====================================================================================
!
!====================================================================================
       SUBROUTINE Print_PoleNode(Node,Tag_O)
          TYPE(PoleNode)                 :: Node
          CHARACTER(LEN=*),OPTIONAL      :: Tag_O 
          CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
          INTEGER                        :: L,M,LMDx
          REAL(DOUBLE)                   :: Cpp,Spp
!
          WRITE(*,*)'======================================================'
          IF(PRESENT(Tag_O))WRITE(*,*)Tag_O          
          IF(Node%Leaf)THEN
             Line='Q[#'//TRIM(IntToChar(Node%Box%Number))//',Tier' &
                 //TRIM(IntToChar(Node%Box%Tier))//',Leaf] = (' &
                 //TRIM(DblToMedmChar(Node%Box%Center(1)))//', '        &
                 //TRIM(DblToMedmChar(Node%Box%Center(2)))//', '        &
                 //TRIM(DblToMedmChar(Node%Box%Center(3)))//') '
          ELSE
             Line='Q[#'//TRIM(IntToChar(Node%Box%Number))//',Tier' &
                 //TRIM(IntToChar(Node%Box%Tier))//'] = (' &
                 //TRIM(DblToMedmChar(Node%Box%Center(1)))//', '        &
                 //TRIM(DblToMedmChar(Node%Box%Center(2)))//', '        &
                 //TRIM(DblToMedmChar(Node%Box%Center(3)))//') '
          ENDIF
          WRITE(*,*)TRIM(Line)
          DO l=0,Node%Ell
             DO m=0,l
               lmdx=LTD(l)+m
               Cpp = Node%C(lmdx)
               Spp = Node%S(lmdx)
               Line='L = '//TRIM(IntToChar(l))//' M = '//TRIM(IntToChar(m)) &
                    //' Cq = '//TRIM(DblToMedmChar(Cpp)) &
                    //' Sq = '//TRIM(DblToMedmChar(Spp))
               WRITE(*,*)TRIM(Line)
            ENDDO
         ENDDO
       END SUBROUTINE Print_PoleNode
!========================================================================================
! Print Fields
!========================================================================================
  SUBROUTINE Print_SP(Ell,C,S,Tag_O,Pre_O)
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O,Pre_O 
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Line
    INTEGER                          :: Ell,L,M,LMDx,IPre
    REAL(DOUBLE)                     :: Cpp,Spp
    REAL(DOUBLE),DIMENSION(0:)       :: C,S
    WRITE(*,*)'======================================================' 
    IF(PRESENT(Tag_O)) WRITE(*,*) Tag_O
    IF(PRESENT(Pre_O)) THEN
       IF(Pre_O == 'Short') IPre = 0
       IF(Pre_O == 'Med'  ) IPre = 1
       IF(Pre_O == 'Long' ) IPre = 2       
    ELSE
       IPre = 1
    ENDIF
!
    IF(IPre == 0) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero 
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero 
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToShrtChar(Cpp)) &
                  //' Sq = '//TRIM(DblToShrtChar(Spp))
             WRITE(*,*) TRIM(Line)
          ENDDO
       ENDDO
    ELSEIF(IPre == 1) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero 
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero 
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToMedmChar(Cpp)) &
                  //' Sq = '//TRIM(DblToMedmChar(Spp))
             WRITE(*,*)TRIM(Line)
          ENDDO
       ENDDO
    ELSEIF(IPre == 2) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero 
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero 
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToChar(Cpp)) &
                  //' Sq = '//TRIM(DblToChar(Spp))
!             IF(.NOT. (Cpp == Zero .AND. Spp == Zero)) THEN
                WRITE(*,*)TRIM(Line)
!             ENDIF
          ENDDO
       ENDDO
    ENDIF
!
  END SUBROUTINE Print_SP
!====================================================================================
!     Xlate in F90
! 
!        OO_{l1,m1}[Q] = Sum_{l2,m2} OO_{l1-l2,m1-m2}[Q-P] OO_{l2,m2}[P]
!
!        |M1-M2| .LE. |L1-L2|= L3 ==> M2 : MAX(-L2,M1-L3) to MIN(L2,M1+L3)
!      
!====================================================================================
      SUBROUTINE XLate90(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
        INTEGER                         :: LP,LQ
        INTEGER                         :: L1,L2,L3,M1,M2,M3,LDX1,LDX2,LDX3,ABSM2,ABSM3
        REAL(DOUBLE)                    :: CN,SN,CMN,SMN
        REAL(DOUBLE),DIMENSION(0:FFLen) :: Cp,Sp,Cq,Sq
        REAL(DOUBLE),DIMENSION(0:FFLen2):: Cpq,Spq
!
        DO L1 = 0,LP
           DO M1 = 0,L1
              LDX1 = LTD(L1)+M1
              DO L2 = 0,MIN(L1,LQ)
                 L3    = L1-L2
                 DO M2 = -L2,L2
                    ABSM2 = ABS(M2)
                    LDX2  = LTD(L2)+ABSM2
                    M3    = M1-M2
                    ABSM3 = ABS(M3)
                    LDX3  = LTD(L3)+ABSM3
                    IF(ABSM3 .LE. L3) THEN
                       IF(M2 .LT. 0) THEN
                          CN = (-One)**(ABSM2)
                          SN = -CN
                       ELSE
                          CN = One
                          SN = One
                       ENDIF
                       IF(M3 .LT. 0) THEN
                          CMN = (-One)**(ABSM3)
                          SMN = -CMN
                       ELSE
                          CMN = One
                          SMN = One
                       ENDIF
                       Cp(LDX1) = Cp(LDX1)+CN*CMN*Cq(LDX2)*Cpq(LDX3) &
                                          -SN*SMN*Sq(LDX2)*Spq(LDX3)
                       Sp(LDX1) = Sp(LDX1)+CN*SMN*Cq(LDX2)*Spq(LDX3) &
                                          +SN*CMN*Sq(LDX2)*Cpq(LDX3)

                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!
      END SUBROUTINE XLate90

END MODULE 
