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
   USE Globals
   IMPLICIT NONE
!=================================================================================================
!  Hierarchical density node
!=================================================================================================
   TYPE PoleNode
      LOGICAL                               :: Leaf     ! Is this a data containing node?
!     Indexes USEd for tree building
      INTEGER                               :: Bdex     ! Begining index of ORB list for this node
      INTEGER                               :: Edex     ! ENDign index of ORB list for this node
      INTEGER                               :: NQ       ! Number of centers
      TYPE(BBox)                            :: Box      ! Bounding box 
      INTEGER                               :: Ell      ! Ell type
      REAL(DOUBLE)                          :: Zeta     ! Minimum exponent in this node
      REAL(DOUBLE)                          :: D2       ! Max distance from node center to BBox
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
!  Global Array intermediates for multipole opperations
!====================================================================================================
   REAL(DOUBLE), DIMENSION(0:SPEll2) :: FactOlm0,FactMlm0,Sine,Cosine,CoFact
   REAL(DOUBLE), DIMENSION(0:SPLen2) :: FactOlm2,FactMlm2,ALegendreP,Spq,Cpq
   REAL(DOUBLE), DIMENSION(0:SPLen)  :: FarFC,FarFS
   CONTAINS
!====================================================================================
!
!====================================================================================
      SUBROUTINE FarField(Q)
         TYPE(PoleNode) :: Q
#ifdef PERIODIC
         CALL PBCTensor(Cpq,Spq)
         CALL CTraX77(SPEll,SPEll,FarFC,FarFS,Cpq,Spq,Q%C,Q%S)
#else
         FarFC=Zero
         FarFS=Zero
#endif
      END SUBROUTINE FarField
!====================================================================================
!     Q->P
!====================================================================================
      SUBROUTINE XLate(Q,P)
         TYPE(PoleNode)            :: Q,P
         REAL(DOUBLE),DIMENSION(3) :: QP
         INTEGER                   :: LP,LQ,LPQ,LenP,LenQ,LenPQ,SPell4,Status,LMDex,L,M
!------------------------------------------------------------------------------------
         QP=(Q%Box%Center-P%Box%Center)
         LP=P%Ell
         LQ=Q%Ell
         LPQ=LP+LQ
         LenPQ=LSP(LPQ)
         CALL Regular(LPQ,QP(1),QP(2),QP(3))       
         CALL XLate77(LP,LQ,P%C,P%S,Cpq,Spq,Q%C,Q%S)
       END SUBROUTINE XLate
!====================================================================================
!
!====================================================================================
      SUBROUTINE CTrax(P,Q,SPKetC,SPKetS)
         TYPE(PrimPair)            :: P
         TYPE(PoleNode)            :: Q
         REAL(DOUBLE),DIMENSION(3) :: PQ
         INTEGER                   :: LP,LQ,LPQ
         REAL(DOUBLE),DIMENSION(0:):: SPKetC,SPKetS
!------------------------------------------------------------------------------------
         PQ=(P%P-Q%Box%Center)
         LP=P%Ell
         LQ=Q%Ell
         LPQ=LP+LQ
         CALL IrRegular90(LPQ,PQ(1),PQ(2),PQ(3))
         CALL CTraX77(LP,LQ,SPKetC,SPKetS,Cpq,Spq,Q%C,Q%S)
      END SUBROUTINE CTrax
!====================================================================================
!
!====================================================================================
       SUBROUTINE HGToSP_PoleNode(Node)
         TYPE(PoleNode) :: Node
         INTEGER        :: LQ,LenQ
         REAL(DOUBLE)   :: PiZ
         LQ=Node%Ell
         LenQ=LSP(LQ)
         PiZ=(Pi/Node%Zeta)**(ThreeHalves)
         CALL HGToSP_Gen(Node%Ell,PiZ,Node%Co,Node%C,Node%S)
       END SUBROUTINE HGToSP_PoleNode
!====================================================================================
!
!====================================================================================
       SUBROUTINE HGToSP_Bra(P,HGBra,SPBraC,SPBraS)
         TYPE(PrimPair)                   :: P
         INTEGER                          :: L,Len
         REAL(DOUBLE)                     :: PiZ,Zeta
         REAL(DOUBLE), DIMENSION(:)       :: HGBra
         REAL(DOUBLE), DIMENSION(0:SPLen) :: SPBraC,SPBraS
#ifdef PERIODIC
         REAL(DOUBLE), DIMENSION(0:SPLen) :: XLBraC,XLBraS
         REAL(DOUBLE),DIMENSION(3)        :: QP
         INTEGER                          :: LP,LQ,LPQ
#endif
!------------------------------------------------------------------------------------
!        Transform <Bra| coefficients from HG to SP
         PiZ=(Pi/P%Zeta)**(ThreeHalves)
         CALL HGToSP_Gen(P%Ell,PiZ,HGBra,SPBraC,SPBraS)
#ifdef PERIODIC
!        Translate <Bra| SPs to Cell center
         QP=(Cell%Box%Center-P%Box%Center)
         LP=P%Ell
         LQ=SPEll
         LPQ=LP+LQ
         CALL Regular(LPQ,QP(1),QP(2),QP(3))       
         CALL XLate77(LP,LQ,SPBraC,SPBraS,Cpq,Spq,XLBraC,XLBraS)         
         SPBraC=XLBraC
         XLBraS=XLBraS
#endif         
       END SUBROUTINE HGToSP_Bra
!====================================================================================
!
!====================================================================================
       SUBROUTINE HGToSP_Gen(L,PiZ,HGCo,C,S) 
          INTEGER                    :: L
          REAL(DOUBLE)               :: PiZ
          REAL(DOUBLE),DIMENSION(:)  :: HGCo      
          REAL(DOUBLE),DIMENSION(0:) :: C,S
          REAL(DOUBLE),DIMENSION(20) :: W
          SELECT CASE(L)
          INCLUDE 'HGToSP.Inc'
          CASE(HGEll+1:) 
             CALL Halt(' Bad logic in HGToSPX')
          END SELECT
       END SUBROUTINE HGToSP_Gen
!====================================================================================
!
!====================================================================================
       SUBROUTINE Print_PoleNode(Node,Tag_O)
          TYPE(PoleNode)                 :: Node
          CHARACTER(LEN=*),OPTIONAL      :: Tag_O 
          CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
          INTEGER                        :: L,M,LMDx
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
               lmdx=l*(l+1)/2+m
!               IF(ABS(Node%C(lmdx))>1.D-20.AND.ABS(Node%S(lmdx))>1.D-20)THEN
                  Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                   //' Cq = '//TRIM(DblToMedmChar(Node%C(lmdx))) &
                   //' Sq = '//TRIM(DblToMedmChar(Node%S(lmdx)))
                  WRITE(*,*)TRIM(Line)
!               ENDIF
            ENDDO
         ENDDO
       END SUBROUTINE Print_PoleNode
!====================================================================================
!
!====================================================================================
       SUBROUTINE MultipoleSetUp(Ell)
          INTEGER                          :: Ell,L,M,LDex,LMDex
          REAL(DOUBLE), DIMENSION(0:2*Ell) :: Factorial
          REAL(DOUBLE)                     :: Sgn,DblFact,TwoTimes
!------------------------------------------------------------------------------------          
          Factorial(0)=One
          DO L=1,2*Ell 
             Factorial(L)=Factorial(L-1)*DBLE(L)
          ENDDO
          DO L=2,Ell 
            LDex=L*(L+1)/2
            DO M=0,L-2
               LMDex=LDex+M
               FactOlm2(LMDex)=One/DBLE((L+M)*(L-M))
             ENDDO
          ENDDO
          DO L=0,Ell
             LDex=L*(L+1)/2
             DO M=0,Ell
                LMDex=LDex+M
                FactMlm2(LMDex)=DBLE((L+M-1)*(L-M-1))
             ENDDO
          ENDDO
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
      END SUBROUTINE MultipoleSetup
!====================================================================================
!
!====================================================================================
      SUBROUTINE IrRegular90(Ell,PQx,PQy,PQz)
         INTEGER      :: Ell,L,M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
         REAL(DOUBLE) :: PQx,PQy,PQz,PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
                         CoFact,PQToThMnsL
         INTEGER      :: LD
         LD(L)=L*(L+1)/2
!------------------------------------------------------------------------------------
         PQx2=PQx*PQx
         PQy2=PQy*PQy      
         PQ=SQRT(PQx2+PQy2+PQz*PQz)
         OneOvPQ=One/PQ
         IF(Ell==0)THEN
            Cpq(0)=OneOvPQ
            Spq(0)=Zero
            RETURN
         ENDIF
         CoTan=PQz*OneOvPQ
!        Sine and Cosine by recursion
         Cosine(0)=One
         Sine(  0)=Zero
         PQxy=SQRT(PQx2+PQy2)
         IF(PQxy/=Zero)THEN
            Sine(1)=PQx/PQxy
            Cosine(1)=PQy/PQxy
         ELSE
            Sine(1)=0.70710678118654752D0         
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
         Sq=SQRT(One-CoTan*CoTan)
         RS=One
         DO M=0,Ell
            MDex=LD(M)+M
            ALegendreP(MDex)=FactMlm0(M)*RS
            RS=RS*Sq
         ENDDO
         DO M=0,Ell-1
            MDex=LD(M)+M
            MDex1=LD(M+1)+M
            ALegendreP(MDex1)=CoTan*DBLE(2*M+1)*ALegendreP(MDex)
         ENDDO
         DO L=2,Ell         
            CoFact=CoTan*DBLE(2*l-1)
            LDex0=LD(L)
            LDex1=LD(L-1)
            LDex2=LD(L-2)
            DO M=0,L-2
               ALegendreP(LDex0+M)=CoFact*ALegendreP(LDex1+M)-FactMlm2(LDex0+M)*ALegendreP(LDex2+M)
            ENDDO
         ENDDO
!        IrRegular Spharical Harmonics
         PQToThMnsL=OneOvPQ
         DO L=0,Ell
            LDex=LD(L)
            DO M=0,l
               LMDex=LDex+M
               Spq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Sine(M)
               Cpq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Cosine(M)
            ENDDO
            PQToThMnsL=PQToThMnsL*OneOvPQ
         ENDDO
      END SUBROUTINE IrRegular90
!====================================================================================
!
!====================================================================================
      SUBROUTINE Regular(Ell,PQx,PQy,PQz)
         INTEGER      :: Ell,L,M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
         REAL(DOUBLE) :: PQx,PQy,PQz,PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
                         CoFact,PQToThPlsL
         INTEGER      :: LD
         LD(L)=L*(L+1)/2
!------------------------------------------------------------------------------------
         PQx2=PQx*PQx
         PQy2=PQy*PQy      
         PQ=DSQRT(PQx2+PQy2+PQz*PQz)
         IF(PQ==Zero)THEN
           Cpq=Zero
           Spq=Zero
           Cpq(0)=One
           RETURN
         ENDIF
         OneOvPQ=One/PQ
         IF(Ell==0)THEN
            Cpq(0)=OneOvPQ
            RETURN
         ENDIF
         CoTan=PQz*OneOvPQ
!        Sine and Cosine by recursion
         Cosine(0)=One
         Sine(  0)=Zero
         PQxy=SQRT(PQx2+PQy2)
         IF(PQxy/=Zero)THEN
            Sine(1)=PQx/PQxy
            Cosine(1)=PQy/PQxy
         ELSE
            Sine(1)=0.70710678118654752D0         
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
         Sq=SQRT(One-CoTan**2)
         RS=One
         DO M=0,Ell
            MDex=LD(M)+M
            ALegendreP(MDex)=FactOlm0(M)*RS
            RS=RS*Sq
         ENDDO
         DO M=0,Ell-1
            MDex=LD(M)+M
            MDex1=LD(M+1)+M
            ALegendreP(MDex1)=CoTan*ALegendreP(MDex)
         ENDDO
         DO L=2,Ell         
            CoFact=CoTan*DBLE(2*l-1)
            LDex0=LD(L)
            LDex1=LD(L-1)
            LDex2=LD(L-2)
            DO M=0,L-2
               ALegendreP(LDex0+M)=(CoFact*ALegendreP(LDex1+M)-ALegendreP(LDex2+M))*FactOlm2(LDex0+M)
            ENDDO
         ENDDO
!        Regular Spharical Harmonics
         PQToThPlsL=One
         DO L=0,Ell
            LDex=LD(L)
            DO M=0,L
               LMDex=LDex+M
               Spq(LMDex)=PQToThPlsL*ALegendreP(LMDex)*Sine(M)
               Cpq(LMDex)=PQToThPlsL*ALegendreP(LMDex)*Cosine(M)
            ENDDO
            PQToThPlsL=PQToThPlsL*PQ
         ENDDO
      END SUBROUTINE Regular
END MODULE 
