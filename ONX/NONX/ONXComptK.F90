MODULE ONXComptK
!H=================================================================================
!H MODULE ONXComptK
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB ComputeK: Explicit and general code.
!H
!H  Comments:
!H  
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  USE ONXRGen   , ONLY: RGen
  USE ONXGet    , ONLY: GetIntCode,GetIntSpace,GetSubBlk
  USE ONXInLoop , ONLY: Digest,Scatter
  USE ONXVrr    , ONLY: VRRs,VRRl
  USE ONXHrr    , ONLY: HrrBra,HrrKet
#ifdef PARALLEL_ONX
  USE MondoMPI
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC :: ComputeK
  !
CONTAINS
  !
!old#ifdef PARALLEL_ONX
  SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
!old#else
!old  SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,IS,Drv,SubInd,BfnInd)
!old#endif
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,IS,Drv,SubInd,BfnInd)
!H SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
!H
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
#ifdef PARALLEL_ONX
    TYPE(FASTMAT), POINTER    :: K,D
    TYPE(FASTMAT), POINTER    :: P
    TYPE(SRST   ), POINTER    :: U
!old   TYPE(DBuf), INTENT(INOUT) :: DBC,DBD   ! ONX distribution buffers
#else
    TYPE(BCSR), INTENT(IN   ) :: D
    TYPE(BCSR), INTENT(INOUT) :: K
#endif
    TYPE(DBuf), INTENT(INOUT) :: DBC,DBD   ! ONX distribution buffers   !per

    TYPE(BSET), INTENT(IN   ) :: BSc,BSp   ! basis set info
    TYPE(CRDS), INTENT(IN   ) :: GMc,GMp   ! geometry info
!old    TYPE(DBuf)                :: DB        ! ONX distribution buffers
    TYPE(IBuf)                :: IB        ! ONX 2-e eval buffers
    TYPE(DSL )                :: SB        ! ONX distribution pointers
    TYPE(ISpc)                :: IS        ! Array dimms for 2-e routines
    TYPE(IDrv)                :: Drv       ! VRR/contraction drivers
    TYPE(INT_RNK2)            :: SubInd    ! Index -> BCSR converter
    TYPE(INT_VECT),INTENT(IN) :: BfnInd
!-------------------------------------------------------------------
! Misc. internal variables
!-------------------------------------------------------------------
    INTEGER                   :: iTBra,iCBra,TBra,LBra,ITypeC,ITypeA,CBra
    INTEGER                   :: iTKet,iCKet,TKet,LKet,ITypeD,ITypeB,CKet
    INTEGER                   :: LenBra,iDBra,iPBra,iBra
    INTEGER                   :: LenKet,iDKet,iPKet,iKet
    INTEGER                   :: LTot,ShellC,ShellD,IndexA,IndexB
    INTEGER                   :: AtC,KC,CFC,StartLC,StopLC,StrideC,IndexC,NBFC
    INTEGER                   :: AtD,KD,CFD,StartLD,StopLD,StrideD,IndexD,NBFD
    INTEGER                   :: ri,ci,iPtr,iCP,iCL
    INTEGER                   :: I,J,I0,I1,I2,ISL,NInts
    INTEGER                   :: IBD,IBP,IKD,IKP
    INTEGER                   :: NA,NB,NC,ND
    INTEGER                   :: IntSpace,IntCodeV,IntCodeC
    INTEGER                   :: BraSwitch,KetSwitch,IntSwitch
    REAL(DOUBLE)              :: Dcd,SchB,SchK,Test
    REAL(DOUBLE)              :: ACx,ACy,ACz
    REAL(DOUBLE)              :: TmBeg,TmEnd
    TYPE(DBL_VECT)            :: DA                                    ! Buffer to store block of density matrix.
    LOGICAL                   :: Explicit
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
    INTEGER                   :: LTotal,MaxBatchSize,NFinal,iT
    !
    !
    CALL New(DA,BSp%LMNLen*BSp%LMNLen)
    xNERIs = 0.0D0

!old#ifdef PARALLEL_ONX
    DO iTBra = 1,DBC%LenTC                      ! Loop over angular symmetry types on the Bra
       TBra = DBC%TCode%I(iTBra)
!old#else
!old    DO iTBra = 1,DB%LenTC                       ! Loop over angular symmetry types on the Bra
!old       TBra = DB%TCode%I(iTBra)
!old#endif
       ITypeA = MOD(TBra,100)
       ITypeC = (TBra-ITypeA)/100
       LBra   = LTotal(ITypeA)+LTotal(ITypeC)
!old#ifdef PARALLEL_ONX
       DO iTKet = 1,DBD%LenTC                   ! Loop over angular symmetry types on the Ket
          TKet = DBD%TCode%I(iTKet)
!old#else
!old       DO iTKet = 1,DB%LenTC                    ! Loop over angular symmetry types on the Ket
!old          TKet = DB%TCode%I(iTKet)
!old#endif
          ITypeB = MOD(TKet,100)
          ITypeD = (TKet-ITypeB)/100
          LKet = LTotal(ITypeB)+LTotal(ITypeD)
          LTot = LBra+LKet

          CALL GetIntCode(LTot,TBra,TKet,IntCodeV,IntCodeC,Explicit)
          CALL GetIntSpace(TBra,TKet,LBra,LKet,IS)

          I0  = iT(MIN(ITypeC,ITypeA),MAX(ITypeC,ITypeA))
          I1  = iT(MIN(ITypeD,ITypeB),MAX(ITypeD,ITypeB))
          I2  = I0+(I1-1)*10
          iCP = Drv%CDrv%I(I2)               ! The pointer to the 2e contraction
          iCL = Drv%CDrv%I(iCP)              ! driver
          CALL VRRs(LBra,LKet,Drv)         ! Get the pointers to the VRR table

!old#ifdef PARALLEL_ONX
          DO iCBra = 1,DBC%LenCC             ! Loop over contraction lengths on the Bra
             CBra = DBC%CCode%I(iCBra)
             IF (DBC%TCPop%I(iTBra,iCBra)==1) THEN
!old#else
!old          DO iCBra = 1,DB%LenCC              ! Loop over contraction lengths on the Bra
!old             CBra = DB%CCode%I(iCBra)
!old             IF (DB%TCPop%I(iTBra,iCBra)==1) THEN
!old#endif

!old#ifdef PARALLEL_ONX
                DO iCKet = 1,DBD%LenCC       ! Loop over contraction lengths on the Ket
                   CKet = DBD%CCode%I(iCKet)
                   IF (DBD%TCPop%I(iTKet,iCKet)==1) THEN
!old#else
!old                DO iCKet = 1,DB%LenCC        ! Loop over contraction lengths on the Ket
!old                   CKet = DB%CCode%I(iCKet)
!old                   IF (DB%TCPop%I(iTKet,iCKet)==1) THEN
!old#endif
                      ShellC = 0
#ifdef PARALLEL_ONX
                      P => D%Next                              ! Loop over atom C
                      DO                               
                         IF(.NOT.ASSOCIATED(P)) EXIT   
                         AtC = P%Row                   
#else
                      DO AtC = 1,NAtoms                        ! Loop over atom C
                         ri = AtC
#endif
                         KC = GMp%AtTyp%I(AtC)
                         NBFC = BSp%BfKnd%I(KC)

                         !     !
                         !     !
                         !     !
                         IndexC = 0                                           !C This should go down
                         DO CFC = 1,BSp%NCFnc%I(KC)                           !C
                            ShellC  = ShellC+1  !BfnInd%I(AtC)+CFC            !C I can use something like shellc=shellc+1,only for C
                            StartLC = BSp%LStrt%I(CFC,KC)                     !C
                            StopLC  = BSp%LStop%I(CFC,KC)                     !C
                            StrideC = StopLC-StartLC+1                        !C
!old#ifdef PARALLEL_ONX                                                           !C     !test_para
                            LenBra  = DBC%DisPtr%I(1,ShellC,iTBra,iCBra)      !C     !test_para
                            iDBra   = DBC%DisPtr%I(2,ShellC,iTBra,iCBra)      !C     !test_para
                            iPBra   = DBC%DisPtr%I(3,ShellC,iTBra,iCBra)      !C     !test_para
!old#else                                                                         !C     !test_para
!old                            LenBra  = DB%DisPtr%I(1,ShellC,iTBra,iCBra)       !C
!old                            iDBra   = DB%DisPtr%I(2,ShellC,iTBra,iCBra)       !C
!old                            iPBra   = DB%DisPtr%I(3,ShellC,iTBra,iCBra)       !C
!old#endif                                                                        !C     !test_para
                            !                                                 !C
                            IF(LenBra>0) THEN                                 !C
                            !   ShellD=0   I do not need                      !C
                               !
                               !
                               !
#ifdef PARALLEL_ONX
                               U => P%RowRoot                               ! Loop over atom D (Dcd)
                               DO                                
                                  IF(.NOT.ASSOCIATED(U)) EXIT    
                                  IF(U%L.NE.U%R) THEN            
                                     U => U%Next                 
                                     CYCLE                       
                                  ENDIF
                                  AtD = U%L                      
                                  ! Set Time.                    
                                  TmBeg = MPI_WTIME()            
                                  !call cpu_time(TmBeg)          
#else
                               DO ci = D%RowPt%I(ri),D%RowPt%I(ri+1)-1        ! Loop over atom D (Dcd)
                                  AtD    = D%ColPt%I(ci)
                                  iPtr   = D%BlkPt%I(ci)
#endif
                                  KD     = GMp%AtTyp%I(AtD)
                                  NBFD   = BSp%BfKnd%I(KD)
                                  !                                                !C Part should come here.
                                  !                                                !C Part should come here.
                                  !                                                !C Part should come here.
                                  !                                                !C Part should come here.
!test                            IndexC=0                                          !C ->>> This should go down
!test                            DO CFC=1,BSp%NCFnc%I(KC)                          !C
!test                               ShellC  = ShellC+1 !ShellC  = BfnInd%I(AtC,CFC)!C!ShellC  = BfnInd%I( AtC)+CFC 
!test                               StartLC = BSp%LStrt%I(CFC,KC)                  !C
!test                               StopLC  = BSp%LStop%I(CFC,KC)                  !C
!test                               StrideC = StopLC-StartLC+1                     !C
!test#ifdef PARALLEL_ONX                                 !C
!test                               LenBra  = DBC%DisPtr%I(1,ShellC,iTBra,iCBra)   !C
!test                               iDBra   = DBC%DisPtr%I(2,ShellC,iTBra,iCBra)   !C
!test                               iPBra   = DBC%DisPtr%I(3,ShellC,iTBra,iCBra)   !C
!test#else                                               !C
!test                               LenBra  = DB%DisPtr%I(1,ShellC,iTBra,iCBra)    !C
!test                               iDBra   = DB%DisPtr%I(2,ShellC,iTBra,iCBra)    !C
!test                               iPBra   = DB%DisPtr%I(3,ShellC,iTBra,iCBra)    !C
!test#endif                                              !C
!test                               !                                              !C
!test                               IF(LenBra>0) THEN                              !C
!test                                  !ShellD=0                                   !C <<<-
                                  !                                                !C Part should come here.
                                  !                                                !C Part should come here.
                                  !                                                !C Part should come here.
                                  !                                                !C Part should come here.

                                  IndexD = 0
                                  DO CFD = 1,BSp%NCFnc%I(KD)
                                     !
                                     ShellD  = BfnInd%I(AtD)+CFD                   !I need to replace that for para
                                     StartLD = BSp%LStrt%I(CFD,KD)
                                     StopLD  = BSp%LStop%I(CFD,KD)
                                     StrideD = StopLD-StartLD+1
!old#ifdef PARALLEL_ONX
                                     LenKet  = DBD%DisPtr%I(1,ShellD,iTKet,iCKet)
                                     iDKet   = DBD%DisPtr%I(2,ShellD,iTKet,iCKet)
                                     iPKet   = DBD%DisPtr%I(3,ShellD,iTKet,iCKet)
!old#else
!old                                     LenKet  = DB%DisPtr%I(1,ShellD,iTKet,iCKet)
!old                                     iDKet   = DB%DisPtr%I(2,ShellD,iTKet,iCKet)
!old                                     iPKet   = DB%DisPtr%I(3,ShellD,iTKet,iCKet)
!old#endif
                                     IF (LenKet>0) THEN
#ifdef PARALLEL_ONX
                                        CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                                             &         IndexD+1,U%MTrix(1,1),DA%D(1))
#else
                                        CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                                             &         IndexD+1,D%MTrix%D(iPtr),DA%D(1))
#endif
                                        Dcd = GetAbsMax(StrideC*StrideD,DA)   ! could test on the density here...
!old#ifdef PARALLEL_ONX
                                        IF (DBC%DisBuf%D(iDBra+1)<0.0D0) THEN
!old#else
!old                                        IF (DB%DisBuf%D(iDBra+1)<0.0D0) THEN
!old#endif
                                           BraSwitch = 1
                                           NA = IS%L1
                                           NC = IS%L2
                                        ELSE
                                           BraSwitch = 0
                                           NA = IS%L2
                                           NC = IS%L1
                                        END IF
!old#ifdef PARALLEL_ONX
                                        IF (DBD%DisBuf%D(iDKet+1)<0.0D0) THEN
!old#else
!old                                        IF (DB%DisBuf%D(iDKet+1)<0.0D0) THEN
!old#endif
                                           KetSwitch = 1
                                           NB = IS%L3
                                           ND = IS%L4
                                        ELSE
                                           KetSwitch = 0
                                           NB = IS%L4
                                           ND = IS%L3
                                        END IF

                                        IntSwitch = 10*BraSwitch+KetSwitch

                                        DO I = 1,LenBra
!old#ifdef PARALLEL_ONX
                                           IBD    = iDBra+(I-1)*DBC%MAXC
                                           IBP    = iPBra+(I-1)*DBC%MAXP*CBra
                                           IndexA = ABS(DBC%DisBuf%D(IBD+1))
                                           ACx    = DBC%DisBuf%D(IBD+4)
                                           ACy    = DBC%DisBuf%D(IBD+5)
                                           ACz    = DBC%DisBuf%D(IBD+6)
                                           SchB   = DBC%DisBuf%D(IBD+10)
!old#else
!old                                           IBD    = iDBra+(I-1)*DB%MAXC
!old                                           IBP    = iPBra+(I-1)*DB%MAXP*CBra
!old                                           IndexA = ABS(DB%DisBuf%D(IBD+1))
!old                                           ACx    = DB%DisBuf%D(IBD+4)
!old                                           ACy    = DB%DisBuf%D(IBD+5)
!old                                           ACz    = DB%DisBuf%D(IBD+6)
!old                                           SchB   = DB%DisBuf%D(IBD+10)
!old#endif
                                           Test = Thresholds%TwoE/(Dcd*SchB)
                                           ISL = 0
                                           DO J = 1,LenKet
!old#ifdef PARALLEL_ONX
                                              IKD = iDKet+(J-1)*DBD%MAXC
                                              IndexB = ABS(DBD%DisBuf%D(IKD+1))
                                              IF (IndexA>=IndexB) THEN      ! Symmetry of the K matrix
                                                 SchK = DBD%DisBuf%D(IKD+10)
                                                 IF (SchK<=Test) EXIT       ! ONX skipout
                                                 IKP = iPKet+(J-1)*DBD%MAXP*CKet
!old#else
!old                                              IKD = iDKet+(J-1)*DB%MAXC
!old                                              IndexB = ABS(DB%DisBuf%D(IKD+1))
!old                                              IF (IndexA>=IndexB) THEN      ! Symmetry of the K matrix
!old                                                 SchK = DB%DisBuf%D(IKD+10)
!old                                                 IF (SchK<=Test) EXIT       ! ONX skipout
!old                                                 IKP = iPKet+(J-1)*DB%MAXP*CKet
!old#endif
                                                 ISL = ISL+1
                                                 SB%SLDis%I(ISL) = IKD+4
                                                 SB%SLPrm%I(ISL) = IKP
                                              END IF
                                           END DO ! J, LenKet

                                           IF (ISL>0) THEN
                                              IntSpace = ISL*IS%NVRR  
                                              IF (IntSpace.GT.IB%MAXI) THEN
                                                 ErrorCode = eMAXI
                                                 GOTO 9000
                                              ENDIF

                                              IF(Explicit) THEN
                                                 SELECT CASE (IntCodeV)
                                                 CASE (01010101)  
!old#ifdef PARALLEL_ONX
                                                    CALL Int1111(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                         &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int1111(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (01010201)  
!old#ifdef PARALLEL_ONX
                                                        CALL Int1121(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int1121(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02010101) 
!old#ifdef PARALLEL_ONX
                                                        CALL Int2111(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2111(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (01010202)  
!old#ifdef PARALLEL_ONX
                                                        CALL Int1122(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int1122(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02010201)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2121(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2121(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02010202)  
!old#ifdef PARALLEL_ONX
                                                        CALL Int2122(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2122(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02020101)  
!old#ifdef PARALLEL_ONX
                                                        CALL Int2211(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2211(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02020201)  
!old#ifdef PARALLEL_ONX
                                                        CALL Int2221(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2221(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02020202)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2222(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2222(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (01010601)
!old#ifdef PARALLEL_ONX
                                                        CALL Int1161(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int1161(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (01010602)
!old#ifdef PARALLEL_ONX
                                                        CALL Int1162(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int1162(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (01010606)
!old#ifdef PARALLEL_ONX
                                                        CALL Int1166(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int1166(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02010601)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2161(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2161(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02010602)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2162(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2162(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02020601)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2261(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2261(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06010101)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6111(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6111(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06010201)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6121(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6121(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06010202)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6122(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6122(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                 
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06010601)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6161(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6161(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06020101)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6211(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6211(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06020201)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6221(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6221(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06060101)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6611(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6611(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02010606)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2166(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2166(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02020602)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2262(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2262(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (02020606)
!old#ifdef PARALLEL_ONX
                                                        CALL Int2266(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int2266(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06010602)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6162(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6162(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06010606)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6166(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6166(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06020202)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6222(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6222(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06020601)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6261(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6261(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06020602)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6262(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6262(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06020606)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6266(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6266(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06060201)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6621(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6621(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06060202)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6622(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6622(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06060601)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6661(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6661(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06060602)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6662(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6662(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE (06060606)
!old#ifdef PARALLEL_ONX
                                                        CALL Int6666(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#else
!old                                                    CALL Int6666(ISL,IntCodeC,CBra,CKet,DB%DisBuf%D(IBD),  &                  
!old                                                         &       DB%PrmBuf%D(IBP) ,DB ,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
!old#endif
                                                 CASE DEFAULT
                                                    WRITE(*,*) "IntCode=",IntCodeV,IntCodeC
                                                    CALL Halt(' Illegal integral type in ONX:ComputeKe')
                                                 END SELECT
                                              ELSE
!old#ifdef PARALLEL_ONX
                                                 CALL RGen(ISL,Ltot,CBra,CKet,IB%CB%D(1,1),IB%CK%D(1,1,1), &
                                                      &    DBC%DisBuf%D(IBD),DBC%PrmBuf%D(IBP),IB%W1%D(1), &
                                                      &    DBD,IB,SB,GMc%PBC)                                                 !per
!old#else
!old                                                 CALL RGen(ISL,Ltot,CBra,CKet,IB%CB%D(1,1),IB%CK%D(1,1,1), &
!old                                                      &    DB%DisBuf%D(IBD) ,DB%PrmBuf%D(IBP) ,IB%W1%D(1), &
!old                                                      &    DB ,IB,SB,GMc%PBC)                                                 !per
!old#endif
                                                 CALL VRRl(ISL*CBra*CKet,IS%NVRR,Drv%nr,Drv%ns,   &
                                                      &    Drv%VLOC%I(Drv%is),                    &
                                                      &    Drv%VLOC%I(Drv%is+Drv%nr),IB,          &
                                                      &    IB%W2%D(1),IB%W1%D(1))
                                                 CALL Contract(ISL,CBra,CKet,IS%NVRR,iCL,Drv%CDrv%I(iCP+1),&
                                                      &        IB%CB%D(1,1),IB%CK%D(1,1,1),IB%W1%D(1),IB%W2%D(1))
                                              ENDIF

                                              xNERIs = xNERIs+FLOAT(IS%L1*IS%L2*IS%L3*IS%L4*ISL)
!old#ifdef PARALLEL_ONX
                                              IF(LKet>0) &
                                                   CALL HRRKet(IB%W1%D(1),DBD%DisBuf%D(1),ISL,  &
                                                   &           SB%SLDis%I(1),IS%NB1,IS%NB2,IS%NK1,TKet)
!old#else
!old                                              IF(LKet>0) & 
!old                                                   CALL HRRKet(IB%W1%D(1),DB%DisBuf%D(1),ISL,   &
!old                                                   &           SB%SLDis%I(1),IS%NB1,IS%NB2,IS%NK1,TKet)
!old#endif
                                              IF(LBra>0) THEN 
                                                 CALL HRRBra(IB%W1%D(1),IB%W2%D(1),ACx,ACy,ACz,ISL,   &
                                                      &      IS%NB1,IS%L1*IS%L2,IS%L3*IS%L4,TBra)
                                                 CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4, &
                                                      &      IntSwitch,IB%W1%D(1),IB%W2%D(1),DA%D(1))
!old#ifdef PARALLEL_ONX
                                                 CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DBD,IB%W1%D(1),K)
!old#else
!old                                                 CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB ,IB%W1%D(1),K)
!old#endif
                                              ELSE
                                                 CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4, &
                                                      &      IntSwitch,IB%W2%D(1),IB%W1%D(1),DA%D(1))
!old#ifdef PARALLEL_ONX
                                                 CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DBD,IB%W2%D(1),K)
!old#else
!old                                                 CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB ,IB%W2%D(1),K)
!old#endif
                                              ENDIF
                                           ENDIF  ! ISL
                                        ENDDO ! I, LenBra
                                     ENDIF ! LenKet
                                     IndexD = IndexD+StrideD
                                  ENDDO ! CFD

                                  !                                         !C Should come here
                                  !                                         !C Should come here
              !test                 ENDIF ! LenBra                          !C This part should go up
              !test                 IndexC = IndexC+StrideC                 !C This part should go up
              !test              ENDDO ! CFC                                !C This part should go up
                                  !                                         !C Should come here
                                  !                                         !C Should come here

#ifdef PARALLEL_ONX
                                  ! Set Time.
                                  TmEnd = MPI_WTIME()
                                  !call cpu_time(TmEnd)
                                  !Add Time.
                                  U%Part = U%Part+TmEnd-TmBeg
                                  U => U%Next
#endif 
                               ENDDO ! ci

                               !
                               !
                            ENDIF ! LenBra                                !C This part should go up
                            IndexC = IndexC+StrideC                       !C This part should go up
                         ENDDO ! CFC                                      !C This part should go up
                         !
                         !
#ifdef PARALLEL_ONX
                         P => P%Next
#endif
                      END DO ! AtC
                   END IF ! TCPop on Ket
                END DO ! iCKet
             END IF ! TCPop on Bra
          END DO ! iCBra
       END DO ! iTKet
    END DO ! iTBra
    !
9000 CONTINUE
    !
    CALL Delete(DA)
    !
  END SUBROUTINE ComputeK
  !
  !
END MODULE ONXComptK

