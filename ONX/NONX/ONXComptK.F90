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
  USE ONXInLoop , ONLY: Digest,Scatter,Contract
  USE ONXVrr    , ONLY: VRRs,VRRl
  USE ONXHrr    , ONLY: HrrBra,HrrKet
#ifdef PARALLEL
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
  SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,IS,Drv,SubInd,BfnInd)
!H SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
!H
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(FASTMAT), POINTER    :: K,D
    TYPE(FASTMAT), POINTER    :: P
    TYPE(SRST   ), POINTER    :: U
#else
    TYPE(BCSR), INTENT(IN   ) :: D
    TYPE(BCSR), INTENT(INOUT) :: K
#endif
    TYPE(DBuf), INTENT(INOUT) :: DBC,DBD   ! ONX distribution buffers   !per

    TYPE(BSET), INTENT(IN   ) :: BSc,BSp   ! basis set info
    TYPE(CRDS), INTENT(IN   ) :: GMc,GMp   ! geometry info
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
    REAL(DOUBLE)              :: MondoTimer
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
    INTEGER                   :: LTotal,MaxBatchSize,NFinal,iT
    !
    !
    CALL New(DA,BSp%LMNLen*BSp%LMNLen)
    xNERIs = 0.0D0
    DO iTBra = 1,DBC%LenTC                      ! Loop over angular symmetry types on the Bra
       TBra = DBC%TCode%I(iTBra)
       ITypeA = MOD(TBra,100)
       ITypeC = (TBra-ITypeA)/100
       LBra   = LTotal(ITypeA)+LTotal(ITypeC)
       DO iTKet = 1,DBD%LenTC                   ! Loop over angular symmetry types on the Ket
          TKet = DBD%TCode%I(iTKet)
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
          CALL VRRs(LBra,LKet,Drv)           ! Get the pointers to the VRR table
          DO iCBra = 1,DBC%LenCC             ! Loop over contraction lengths on the Bra
             CBra = DBC%CCode%I(iCBra)
             IF (DBC%TCPop%I(iTBra,iCBra)==1) THEN
                DO iCKet = 1,DBD%LenCC       ! Loop over contraction lengths on the Ket
                   CKet = DBD%CCode%I(iCKet)
                   IF (DBD%TCPop%I(iTKet,iCKet)==1) THEN
                      ShellC = 0
#ifdef PARALLEL
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
                            LenBra  = DBC%DisPtr%I(1,ShellC,iTBra,iCBra)      !C     !test_para
                            iDBra   = DBC%DisPtr%I(2,ShellC,iTBra,iCBra)      !C     !test_para
                            iPBra   = DBC%DisPtr%I(3,ShellC,iTBra,iCBra)      !C     !test_para
                            !                                                 !C
                            IF(LenBra>0) THEN                                 !C
                            !   ShellD=0   I do not need                      !C
                               !
                               !
                               !
#ifdef PARALLEL
                               U => P%RowRoot                               ! Loop over atom D (Dcd)
                               DO                                
                                  IF(.NOT.ASSOCIATED(U)) EXIT    
                                  IF(U%L.NE.U%R) THEN            
                                     U => U%Next                 
                                     CYCLE                       
                                  ENDIF
                                  AtD = U%L                      
                                  ! Set Time.                    
                                  TmBeg = MondoTimer()            
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
!test#ifdef PARALLEL                                 !C
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
                                     LenKet  = DBD%DisPtr%I(1,ShellD,iTKet,iCKet)
                                     iDKet   = DBD%DisPtr%I(2,ShellD,iTKet,iCKet)
                                     iPKet   = DBD%DisPtr%I(3,ShellD,iTKet,iCKet)
                                     IF (LenKet>0) THEN
#ifdef PARALLEL
                                        CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                                             &         IndexD+1,U%MTrix(1,1),DA%D(1))
#else
                                        CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                                             &         IndexD+1,D%MTrix%D(iPtr),DA%D(1))
#endif
                                        Dcd = GetAbsMax(StrideC*StrideD,DA)   ! could test on the density here...
                                        IF (DBC%DisBuf%D(iDBra+1)<0.0D0) THEN
                                           BraSwitch = 1
                                           NA = IS%L1
                                           NC = IS%L2
                                        ELSE
                                           BraSwitch = 0
                                           NA = IS%L2
                                           NC = IS%L1
                                        END IF
                                        IF (DBD%DisBuf%D(iDKet+1)<0.0D0) THEN
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
                                           IBD    = iDBra+(I-1)*DBC%MAXC
                                           IBP    = iPBra+(I-1)*DBC%MAXP*CBra
                                           IndexA = ABS(DBC%DisBuf%D(IBD+1))
                                           ACx    = DBC%DisBuf%D(IBD+4)
                                           ACy    = DBC%DisBuf%D(IBD+5)
                                           ACz    = DBC%DisBuf%D(IBD+6)
                                           SchB   = DBC%DisBuf%D(IBD+10)
                                           Test = Thresholds%TwoE/(Dcd*SchB)
                                           ISL = 0
                                           DO J = 1,LenKet
                                              IKD = iDKet+(J-1)*DBD%MAXC
                                              IndexB = ABS(DBD%DisBuf%D(IKD+1))
                                              IF (IndexA>=IndexB) THEN      ! Symmetry of the K matrix
                                                 SchK = DBD%DisBuf%D(IKD+10)
                                                 IF (SchK<=Test) EXIT       ! ONX skipout
                                                 IKP = iPKet+(J-1)*DBD%MAXP*CKet
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
                                                    CALL Int1111(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                         &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (01010201)  
                                                        CALL Int1121(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02010101) 
                                                        CALL Int2111(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (01010202)  
                                                        CALL Int1122(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02010201)
                                                        CALL Int2121(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02010202)  
                                                        CALL Int2122(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02020101)  
                                                        CALL Int2211(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02020201)  
                                                        CALL Int2221(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02020202)
                                                        CALL Int2222(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (01010601)
                                                        CALL Int1161(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (01010602)
                                                        CALL Int1162(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (01010606)
                                                        CALL Int1166(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02010601)
                                                        CALL Int2161(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02010602)
                                                        CALL Int2162(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02020601)
                                                        CALL Int2261(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06010101)
                                                        CALL Int6111(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06010201)
                                                        CALL Int6121(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06010202)
                                                        CALL Int6122(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06010601)
                                                        CALL Int6161(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06020101)
                                                        CALL Int6211(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06020201)
                                                        CALL Int6221(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06060101)
                                                        CALL Int6611(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02010606)
                                                        CALL Int2166(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02020602)
                                                        CALL Int2262(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (02020606)
                                                        CALL Int2266(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06010602)
                                                        CALL Int6162(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06010606)
                                                        CALL Int6166(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06020202)
                                                        CALL Int6222(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06020601)
                                                        CALL Int6261(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06020602)
                                                        CALL Int6262(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06020606)
                                                        CALL Int6266(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06060201)
                                                        CALL Int6621(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06060202)
                                                        CALL Int6622(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06060601)
                                                        CALL Int6661(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06060602)
                                                        CALL Int6662(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE (06060606)
                                                        CALL Int6666(ISL,IntCodeC,CBra,CKet,DBC%DisBuf%D(IBD),  &
                                                             &       DBC%PrmBuf%D(IBP),DBD,IB,SB,IB%W1%D(1),IB%W2%D(1),GMc%PBC)  !per
                                                 CASE DEFAULT
                                                    WRITE(*,*) "IntCode=",IntCodeV,IntCodeC
                                                    CALL Halt(' Illegal integral type in ONX:ComputeKe')
                                                 END SELECT
                                              ELSE
                                                 CALL RGen(ISL,Ltot,CBra,CKet,IB%CB%D(1,1),IB%CK%D(1,1,1), &
                                                      &    DBC%DisBuf%D(IBD),DBC%PrmBuf%D(IBP),IB%W1%D(1), &
                                                      &    DBD,IB,SB,GMc%PBC)                                                 !per

                                                 CALL VRRl(ISL*CBra*CKet,IS%NVRR,Drv%nr,Drv%ns,   &
                                                      &    Drv%VLOC%I(Drv%is),                    &
                                                      &    Drv%VLOC%I(Drv%is+Drv%nr),IB,          &
                                                      &    IB%W2%D(1),IB%W1%D(1))

                                                 CALL Contract(ISL,CBra,CKet,IS%NVRR,iCL,Drv%CDrv%I(iCP+1),&
                                                      &        IB%CB%D(1,1),IB%CK%D(1,1,1),IB%W1%D(1),IB%W2%D(1))
                                              ENDIF

                                              xNERIs = xNERIs+FLOAT(IS%L1*IS%L2*IS%L3*IS%L4*ISL)
                                              IF(LKet>0) &
                                                   CALL HRRKet(IB%W1%D(1),DBD%DisBuf%D(1),ISL,  &
                                                   &           SB%SLDis%I(1),IS%NB1,IS%NB2,IS%NK1,TKet)
                                              IF(LBra>0) THEN 
                                                 CALL HRRBra(IB%W1%D(1),IB%W2%D(1),ACx,ACy,ACz,ISL,   &
                                                      &      IS%NB1,IS%L1*IS%L2,IS%L3*IS%L4,TBra)
                                                 CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4, &
                                                      &      IntSwitch,IB%W1%D(1),IB%W2%D(1),DA%D(1))
                                                 CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DBD,IB%W1%D(1),K)
                                              ELSE
                                                 CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4, &
                                                      &      IntSwitch,IB%W2%D(1),IB%W1%D(1),DA%D(1))
                                                 CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DBD,IB%W2%D(1),K)
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

#ifdef PARALLEL
                                  ! Set Time.
                                  TmEnd = MondoTimer()
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
#ifdef PARALLEL
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

