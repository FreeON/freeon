MODULE ParseBasis
  USE Order
  USE InOut
  USE PrettyPrint
  USE ControlStructures
  USE BasisSetParameters
  IMPLICIT NONE
  CHARACTER(LEN=9),  PARAMETER :: BASIS_SETS ='BasisSets' 
CONTAINS
  !============================================================================
  ! 
  ! 
  !============================================================================
  SUBROUTINE LoadBasisSets(N,O,G,B)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(BasisSets)    :: B
    INTEGER            :: I,J
    CHARACTER(LEN=DCL) :: BaseFile
   !-------------------------------------------------------------------------!    
    CALL OpenASCII(N%IFile,Inp)
    ! Find out how many basis sets we are going over, and their names
    CALL ParseBasisNames(B)
    ! Allocate aux stuff...
    ALLOCATE(B%BSiz(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%OffS(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%BSets(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%LnDex(1:G%Clones,1:B%NBSets))
    ALLOCATE(B%DExpt(1:G%Clones,1:B%NBSets))
    DO J=1,B%NBSets
       BaseFile=TRIM(N%M_HOME)//'BasisSets/'//TRIM(B%BName(J))//BasF
       CALL OpenASCII(BaseFile,Bas,OldFileQ_O=.TRUE.)
       DO I=1,G%Clones
          IF(.NOT.ParseBasisSets(G%Clone(I),B%BSets(I,J),B%BSiz(I,J),B%OffS(I,J)))   &
               CALL MondoHalt(PRSE_ERROR,'ParseBasisSets failed for basis set file ' &
                              //RTRN//TRIM(BaseFile))
          B%BSets(I,J)%BName=B%BName(J)
          CALL BCSRDimensions(G%Clone(I),B%BSets(I,J),O%AccuracyLevels(J), &
                              B%MxAts(J),B%MxBlk(J),B%MxN0s(J))
          CALL PrimitiveReDistribution(B%BSets(I,J),B%NExpt(J),B%DExpt(I,J),B%Lndex(I,J))
#ifdef FULL_ON_FRONT_END_DEBUG
          PrintFlags%Set=DEBUG_BASISSET
          CALL Print_BSET(B%BSets(I,J))
#endif
      ENDDO
       CLOSE(Bas,STATUS='KEEP')
    ENDDO
    CLOSE(Inp,STATUS='KEEP')
  END SUBROUTINE LoadBasisSets
  !============================================================================
  ! WEAK ATTEMPT TO ESTIMATE NUMBER OF NONZEROS FOR SPARSE MATRICES, NEADS TWIDDLING
  ! AND SHOULD BECOME OBSOLETE IN NEAR FURTURE AS FASTMAT TAKES OVER...
  !============================================================================
  SUBROUTINE BCSRDimensions(G,B,AccL,MaxAtms,MaxBlks,MaxNon0)
    TYPE(CRDS) :: G
    TYPE(BSET) :: B
    INTEGER    :: AccL,MaxAtms,MaxBlks,MaxNon0
    REAL(DOUBLE) :: BWEstim
    REAL(DOUBLE),PARAMETER,DIMENSION(4) :: BandWidth=(/ 1.D3, 1.D3, 1.3D3,1.6D3/)
    REAL(DOUBLE),PARAMETER,DIMENSION(4) :: BWDecay  =(/ 1.D-4,1.D-4,1.D-3,1.D-2/)
    ! Use assymptotics to set the max matrix dimensions
    MaxAtms=1+G%NAtms
    BWEstim=MIN(G%NAtms,CEILING((DBLE(G%NAtms) &
         +BandWidth(AccL)*BWDecay(AccL)*DBLE(G%NAtms)**2) &
         /(One+BWDecay(AccL)*DBLE(G%NAtms)**2) ) ) 
    MaxBlks=1+G%NAtms*BWEstim
    MaxNon0=1+B%NBasF*(DBLE(B%NBasF)*DBLE(BWEstim)/DBLE(G%NAtms))
  END SUBROUTINE BCSRDimensions
  !============================================================================
  ! PARSE A BASIS SET, SET UP ITS INDECIES, NORMALIZE THE PRIMITIVES AND 
  ! COMPUTE BLOCKING FOR SPARSE BLOCKED LINEAR ALGEBRA
  !============================================================================
  FUNCTION ParseBasisSets(G,B,BlkSiz,OffSet)
    TYPE(BSET)                 :: BS,B
    TYPE(CRDS)                 :: G
    TYPE(INT_VECT)             :: BlkSiz,OffSet
    REAL(DOUBLE),    &
       DIMENSION(1:MaxASymt+2) :: Dbls
    LOGICAL                    :: ParseBasisSets
    CHARACTER(LEN=DCL)         :: Line
    INTEGER                    :: I,J,K,L,N,NC,NK,NP,NS,MinL,MaxL,KFound
    !-------------------------------------------------------------------------!
    ! Allocate temporary set
    BS%LMNLen=MaxAsymt
    BS%NCtrt=MaxCntrx
    BS%NPrim=MaxPrmtv
    BS%NAtms=G%NAtms
    BS%NKind=G%NAtms
    CALL New(BS)
    ! Count kinds and load basis set kind index
    BS%NKind=1
    BS%Kinds%I(1)=G%AtNum%D(1)
    DO I=2,BS%NAtms 
       DO J=1,BS%NKind
          IF(BS%Kinds%I(J)==G%AtNum%D(I))GOTO 10
       ENDDO
       BS%NKind=BS%NKind+1
       BS%Kinds%I(BS%NKind)=G%AtNum%D(I)
10     CONTINUE
    ENDDO
    ! now load geometry atom type (kinds pointer) array.
    G%NKind=BS%NKind
    DO K=1,BS%NKind
       DO I=1,G%NAtms
          IF(BS%Kinds%I(K)==G%AtNum%D(I))G%AtTyp%I(I)=K
       ENDDO
    ENDDO
    ! Zero counters
    BS%NASym=0
    BS%NCtrt=0
    BS%NPrim=0
    DO I=1,BS%NKind
       BS%NCFnc%I(I)=0
       DO J=1,MaxCntrx
          BS%NPFnc%I(J,I)=0
       ENDDO
    ENDDO
    ! Parse basis set 
    KFound=0
    REWIND(Bas)
    DO 
       READ(Bas,DEFAULT_CHR_FMT,END=99)Line                 
       DO NK=1,BS%NKind            
          IF(KeyQ(Line,Ats(BS%Kinds%I(NK))).AND.KeyQ(Line,'0'))THEN                    
             NC=0
             KFound=KFound+1
             DO 
                READ(Bas,DEFAULT_CHR_FMT,END=99)Line         
                IF(KeyQ(Line,Stars))GOTO 100
                NC=NC+1                 
                DO K=1,MaxLTyps
                   IF(KeyQ(Line,CLTyps(K)))THEN
                      BS%ASymm%I(1,NC,NK)=LTyps(1,K)                     
                      BS%ASymm%I(2,NC,NK)=LTyps(2,K)
                      DO L=1,MaxPrmtv
                         IF(KeyQ(Line,TRIM(IntToChar(L))))THEN
                            NP=L
                            GOTO 101
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
                RETURN
101             CONTINUE
                BS%NCFnc%I(NK)=NC
                BS%NPFnc%I(NC,NK)=NP
                BS%NCtrt=Max(BS%NCtrt,NC)
                BS%NPrim=Max(BS%NPrim,NP)
                MinL=BS%ASymm%I(1,NC,NK)
                MaxL=BS%ASymm%I(2,NC,NK)
                BS%NAsym=Max(BS%NAsym,MaxL)
                DO NP=1,BS%NPFnc%I(NC,NK)
                   READ(Bas,DEFAULT_CHR_FMT,END=99)Line
                   N=MaxL-MinL+2
                   CALL LineToDbls(Line,N,Dbls)
                   BS%Expnt%D(NP,NC,NK)=Dbls(1)
                   K=1
                   BS%CCoef%D(1:,NP,NC,NK)=Zero
                   DO NS=MinL,MaxL
                      K=K+1
                      BS%CCoef%D(NS+1,NP,NC,NK)=Dbls(K) 
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
100    CONTINUE
    ENDDO
99  CONTINUE
    CLOSE(UNIT=Bas,STATUS='KEEP')
    IF(KFound/=BS%NKind)RETURN
    ! Computing basis set indexing
    CALL BSetIndx(BS)
    ! Normalize the primitives
    CALL ReNormalizePrimitives(BS)
    ! Set the correct dimensions
    B%NAsym=BS%NAsym
    B%LMNLen=BS%LMNLen
    B%NCtrt=BS%NCtrt
    B%NPrim=BS%NPrim
    B%NAtms=BS%NAtms
    B%NKind=BS%NKind
    ! New basis set
    CALL New(B)
    ! Set old eq to new (should go to seteq sometime soon...)
    B%Kinds%I(1:B%NKind)=BS%Kinds%I(1:B%NKind)
    B%NCFnc%I(1:B%NKind)=BS%NCFnc%I(1:B%NKind)
    B%BFKnd%I(1:B%NAtms)=BS%BFKnd%I(1:B%NAtms)
    B%LxDex%I(1:B%LMNLen)=BS%LxDex%I(1:B%LMNLen)
    B%LyDex%I(1:B%LMNLen)=BS%LyDex%I(1:B%LMNLen)
    B%LzDex%I(1:B%LMNLen)=BS%LzDex%I(1:B%LMNLen)
    B%LStrt%I(1:B%NCtrt,1:B%NKind)=BS%LStrt%I(1:B%NCtrt,1:B%NKind)
    B%LStop%I(1:B%NCtrt,1:B%NKind)=BS%LStop%I(1:B%NCtrt,1:B%NKind)
    B%NPFnc%I(1:B%NCtrt,1:B%NKind)=BS%NPFnc%I(1:B%NCtrt,1:B%NKind)
    B%ASymm%I(:,1:B%NCtrt,1:B%NKind)=BS%ASymm%I(:,1:B%NCtrt,1:B%NKind)
    B%Expnt%D(1:B%NPrim,1:B%NCtrt,1:B%NKind)=BS%Expnt%D(1:B%NPrim,1:B%NCtrt,1:B%NKind)
    B%CCoef%D(1:B%LMNLen,1:B%NPrim,1:B%NCtrt,1:B%NKind)=BS%CCoef%D(1:B%LMNLen,1:B%NPrim,1:B%NCtrt,1:B%NKind)
    ! Done with the temp BS
    CALL Delete(BS)
    ! Compute blocking for sparse matrix methods
    CALL New(BlkSiz,B%NAtms)
    CALL New(OffSet,B%NAtms)
    CALL BlockBuild(G,B,BlkSiz,OffSet)
    ParseBasisSets=.TRUE.
  END FUNCTION ParseBasisSets

  SUBROUTINE ReNormalizePrimitives(A)
    TYPE(BSET)       :: A
    TYPE(DBL_RNK2)   :: AuxCoef
    REAL(DOUBLE)     :: Z,C,Expnt,RNorm,ZA,ZB,CA,CB,SQNrm
    INTEGER          :: K,L,M,N,NC,NK,NP,LMN,MinL,MaxL,PFA,PFB,MaxSym,MaxPrm
    REAL(DOUBLE),PARAMETER, &      ! Fact=Sqrt[Pi](2*L-1)!! 2^(-L)
         DIMENSION(0:10)  :: Fact=(/0.17724538509055160273D1, &
                                              0.8862269254527580136D0,  &
                                              0.13293403881791370205D1, &
                                              0.3323350970447842551D1,  &
                                              0.11631728396567448929D2, &
                                              0.5234277778455352018D2,  &
                                              0.28788527781504436100D3, &
                                              0.1871254305797788347D4,  &
                                              0.14034407293483412599D5, &
                                              0.1192924619946090071D6,  &
                                              0.11332783889487855673D7/) 
    !-------------------------------------------------------------------------------------!
    ! Allocate temp space
    MaxSym=0
    MaxPrm=0
    DO K=1,A%NKind
       DO NC=1,A%NCFnc%I(K)
          MaxSym=MAX(MaxSym,A%ASymm%I(2,NC,K))
          MaxPrm=MAX(MaxPrm,A%NPFnc%I(NC,K))
       ENDDO
    ENDDO
    CALL New(AuxCoef,(/MaxSym+1,MaxPrm/))    
    DO K=1,A%NKind
       DO NC=1,A%NCFnc%I(K)
          MinL=A%ASymm%I(1,NC,K)
          MaxL=A%ASymm%I(2,NC,K) 
          DO NP=1,A%NPFnc%I(NC,K)
             Z=A%Expnt%D(NP,NC,K)
             DO L=MinL,MaxL
                AuxCoef%D(L+1,NP)=A%CCoef%D(L+1,NP,NC,K)*Z**(Half*DBLE(L)+0.75D0)
             ENDDO
          ENDDO
          DO LMN=A%LStrt%I(NC,K),A%LStop%I(NC,K)            
             RNorm=0.0D0
             L=A%LxDex%I(LMN)
             M=A%LyDex%I(LMN)
             N=A%LzDex%I(LMN)
             DO PFA=1,A%NPFnc%I(NC,K)
                ZA=A%Expnt%D(PFA,NC,K)
                CA=AuxCoef%D(L+M+N+1,PFA)
                DO PFB=1,A%NPFnc%I(NC,K)
                   ZB=A%Expnt%D(PFB,NC,K)
                   CB=AuxCoef%D(L+M+N+1,PFB)
                   RNorm=RNorm+CA*CB*(ZA+ZB)**(-DBLE(L+M+N)-1.5D0)
                ENDDO
             ENDDO
             SqNrm=1.0D0/SQRT(RNorm*Fact(L)*Fact(M)*Fact(N))
             DO NP=1,A%NPFnc%I(NC,K)
                A%CCoef%D(LMN,NP,NC,K)=AuxCoef%D(L+M+N+1,NP)*SqNrm
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL Delete(AuxCoef)
  END SUBROUTINE ReNormalizePrimitives

  SUBROUTINE BlockBuild(G,B,BS,OS)
    TYPE(CRDS)      :: G
    TYPE(BSET)      :: B         
    TYPE(INT_VECT)  :: BS,OS
    INTEGER         :: NA,NK,NC,Stride
    !-------------------------------------------------------------------------------------!
    B%NBasF=0
    ! Off set starts at 1
    OS%I(1)=1      
    DO NA=1,G%NAtms
       ! Block size is total number of basis functions per atom 
       BS%I(NA)=0
       NK=G%AtTyp%I(NA)
       ! Go over contracted functions
       DO NC=1,B%NCFnc%I(NK)
          ! Add in size of each contraction
          Stride=B%LStop%I(NC,NK)-B%LStrt%I(NC,NK)+1
          BS%I(NA)=BS%I(NA)+Stride
          ! Oh yeah, accumulate basis function counter too...
          B%NBasF=B%NBasF+Stride
       ENDDO
       ! Off set counter from block sizes
       IF(NA.GE.2)OS%I(NA)=OS%I(NA-1)+BS%I(NA-1)         
    ENDDO
  END SUBROUTINE BlockBuild

  SUBROUTINE ParseBasisNames(B)
    TYPE(BasisSets) :: B
    INTEGER                               :: I,ILoc
    !----------------------------------------------------------------------------
    B%NBSets=0
    DO I=1,NSupSets
       IF(OptKeyLocQ(Inp,BASIS_SETS,TRIM(CSets(1,I)),MaxSets,NLoc,Location))THEN
          B%NBSets=B%NBSets+NLoc
          DO ILoc=1,NLoc 
             B%BName(Location(ILoc))=ADJUSTL(CSets(2,I))
          ENDDO
       ENDIF
    ENDDO
    IF(B%NBSets==0)CALL MondoHalt(PRSE_ERROR,TRIM(BASIS_SETS)//' not found in input.')
  END SUBROUTINE ParseBasisNames
  !=============================================================================
  ! PRECOMPUTE PRIMITIVE DISTRIBUTION INFORMATION
  !=============================================================================
  SUBROUTINE PrimitiveReDistribution(BS,NExpt,DExpt,Lndex)
    TYPE(BSET),     INTENT(INOUT):: BS
    TYPE(DBL_VECT), INTENT(OUT)  :: DExpt
    TYPE(INT_VECT), INTENT(OUT)  :: Lndex
    INTEGER,        INTENT(OUT)  :: NExpt
    TYPE(INT_VECT)               :: ITmp,IPnt
    INTEGER                      :: I,K,KA,KB,CFA,CFB,PFA,PFB,NPB
    !---------------------------------------------------------------------------!
    ! First count the types of primitive distributions(K=1), then compute 
    ! their primitive distribution exponents and associated values of 
    ! angular symmetry (K=2)
    DO K=1,2
       IF(K==1)THEN
          NExpt=1
       ELSE
          CALL New(DExpt,NExpt)
          CALL New(Lndex,NExpt)
          NExpt=1
       ENDIF
       DO KA=1,BS%NKind
          DO KB=1,KA-1
             DO CFA=1,BS%NCFnc%I(KA)
                DO CFB=1,BS%NCFnc%I(KB)
                   DO PFA=1,BS%NPFnc%I(CFA,KA)
                      DO PFB=1,BS%NPFnc%I(CFB,KB)
                         IF(K==2)THEN
                            DExpt%D(NExpt)=BS%Expnt%D(PFA,CFA,KA)  &
                                 +BS%Expnt%D(PFB,CFB,KB)  
                            Lndex%I(NExpt)=BS%ASymm%I(2,CFA,KA)    &
                                 +BS%ASymm%I(2,CFB,KB)
                         ENDIF
                         NExpt=NExpt+1          
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO KA=1,BS%NKind
          DO CFA=1,BS%NCFnc%I(KA)
             DO CFB=1,CFA
                DO PFA=1,BS%NPFnc%I(CFA,KA)
                   NPB=BS%NPFnc%I(CFB,KA)
                   IF(CFA.EQ.CFB)NPB=PFA
                   DO PFB=1,NPB
                      IF(K==2)THEN 
                         DExpt%D(NExpt)=BS%Expnt%D(PFA,CFA,KA)  &
                              +BS%Expnt%D(PFB,CFB,KA)  
                         Lndex%I(NExpt)=BS%ASymm%I(2,CFA,KA)    &
                              +BS%ASymm%I(2,CFB,KA)
                      ENDIF
                      NExpt=NExpt+1 
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! Add exponent for nuclear delta function
    DExpt%D(NExpt)=NuclearExpnt
    Lndex%I(NExpt)=0
    ! Sort the exponents in assending order and carry
    ! allong angular symetry indecies
    CALL New(IPnt,NExpt)
    CALL New(ITmp,NExpt)
    DO I=1,NExpt
       IPnt%I(I)=I
    ENDDO
    CALL Sort(DExpt,IPnt,NExpt,2)   
    DO I=1,NExpt
       ITmp%I(I)=Lndex%I(IPnt%I(I))
    ENDDO
    DO I=1,NExpt
       Lndex%I(I)=ITmp%I(I)
    ENDDO
    ! Tidy up
    CALL Delete(ITmp)
    CALL Delete(IPnt)
  END SUBROUTINE PrimitiveReDistribution    
END MODULE ParseBasis
