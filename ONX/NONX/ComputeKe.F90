SUBROUTINE ComputeKe(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,Drv,SubInd,BfnInd)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  USE GetTables
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)           :: D
  TYPE(DBCSR)           :: K
  TYPE(BCSR)            :: KTotal
#else
  TYPE(BCSR)            :: D
  TYPE(BCSR)            :: K
  INTEGER               :: MyID=0
#endif
  TYPE(BSET),INTENT(IN)     :: BSc,BSp   ! basis set info
  TYPE(CRDS),INTENT(IN)     :: GMc,GMp   ! geometry info
  TYPE(DBuf)                :: DB        ! ONX distribution buffers
  TYPE(IBuf)                :: IB        ! ONX 2-e eval buffers
  TYPE(DSL)                 :: SB        ! ONX distribution pointers
  TYPE(IDrv)                :: Drv       ! VRR/contraction drivers
  TYPE(INT_RNK2)            :: SubInd    ! Index -> BCSR converter
  TYPE(INT_RNK2),INTENT(IN) :: BfnInd
!-------------------------------------------------------------------
! Misc. internal variables
!-------------------------------------------------------------------
  INTEGER               :: iTBra,iCBra,TBra,LBra,ITypeC,ITypeA,CBra
  INTEGER               :: iTKet,iCKet,TKet,LKet,ITypeD,ITypeB,CKet
  INTEGER               :: LenBra,iDBra,iPBra,iBra
  INTEGER               :: LenKet,iDKet,iPKet,iKet
  INTEGER               :: LTot,ShellC,ShellD,IndexA,IndexB
  INTEGER               :: AtC,KC,CFC,StartLC,StopLC,StrideC,IndexC,NBFC
  INTEGER               :: AtD,KD,CFD,StartLD,StopLD,StrideD,IndexD,NBFD
  INTEGER               :: ri,ci,iPtr,iCP,iCL
  INTEGER               :: I,J,I0,I1,I2,ISL,NVRR,NInts
  INTEGER               :: IBD,IBP,IKD,IKP
  INTEGER               :: NB1,NB2,NK1,NK2
  INTEGER               :: NA,NB,NC,ND
  INTEGER               :: L1,L2,L3,L4,IntSpace,IntCode
  INTEGER               :: BraSwitch,KetSwitch,IntSwitch
  REAL(DOUBLE)          :: Dcd,SchB,SchK,Test
  REAL(DOUBLE)          :: ACx,ACy,ACz
  TYPE(DBL_VECT)        :: DA
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
  INTEGER               :: LTotal,MaxBatchSize,NFinal,iT


  CALL GetExpTable(IB)      ! Read in the Exp table
  CALL New(DA,BSp%LMNLen*BSp%LMNLen)
  xNERIs=0.0D0

  DO iTBra=1,DB%LenTC       ! Loop over angular symmetry types on the Bra
    TBra=DB%TCode%I(iTBra)
    ITypeA=MOD(TBra,100)
    ITypeC=(TBra-ITypeA)/100
    LBra=LTotal(ITypeA)+LTotal(ITypeC)

  DO iTKet=1,DB%LenTC       ! Loop over angular symmetry types on the Ket
    TKet=DB%TCode%I(iTKet)
    ITypeB=MOD(TKet,100)
    ITypeD=(TKet-ITypeB)/100
    LKet=LTotal(ITypeB)+LTotal(ITypeD)
    LTot=LBra+LKet

    IF (LTot.LE.0) THEN

                  IntCode=0

    NB1=IDmn(LBra)
    NK1=IDmn(LKet)
    L1=NFinal(ITypeC)
    L2=NFinal(ITypeA)
    L3=NFinal(ITypeD)
    L4=NFinal(ITypeB)
    NB2=L1*L2
    NK2=L3*L4
    NVRR=NB1*NK1
    NInts=NB2*NK2

    I0=iT(MIN(ITypeC,ITypeA),MAX(ITypeC,ITypeA))
    I1=iT(MIN(ITypeD,ITypeB),MAX(ITypeD,ITypeB))
    I2=I0+(I1-1)*10
    iCP=Drv%CDrv%I(I2)            ! The pointer to the 2e contraction
    iCL=Drv%CDrv%I(iCP)           ! driver

    CALL GetGammaTable(LTot,IB)   ! Get the correct gamma fcn table
    CALL VRRs(LBra,LKet,Drv)      ! Get the pointers to the VRR table

  DO iCBra=1,DB%LenCC       ! Loop over contraction lengths on the Bra
    CBra=DB%CCode%I(iCBra)
    IF (DB%TCPop%I(iTBra,iCBra)==1) THEN
  DO iCKet=1,DB%LenCC       ! Loop over contraction lengths on the Ket
    CKet=DB%CCode%I(iCKet)
    IF (DB%TCPop%I(iTKet,iCKet)==1) THEN

      ShellC=0
      DO AtC=1,NAtoms                              ! Loop over atom C
        ri=AtC
        KC=GMp%AtTyp%I(AtC)
        NBFC=BSp%BfKnd%I(KC)
        IndexC=0
        DO CFC=1,BSp%NCFnc%I(KC)
          ShellC  = BfnInd%I(AtC,CFC)
          StartLC = BSp%LStrt%I(CFC,KC)
          StopLC  = BSp%LStop%I(CFC,KC)
          StrideC = StopLC-StartLC+1
          LenBra  = DB%DisPtr%I(1,ShellC,iTBra,iCBra)
          iDBra   = DB%DisPtr%I(2,ShellC,iTBra,iCBra)
          iPBra   = DB%DisPtr%I(3,ShellC,iTBra,iCBra)
          IF (LenBra>0) THEN

      ShellD=0
      DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1        ! Loop over atom D (Dcd)
        AtD    = D%ColPt%I(ci)
        iPtr   = D%BlkPt%I(ci)
        KD     = GMp%AtTyp%I(AtD)
        NBFD   = BSp%BfKnd%I(KD)
        IndexD = 0
        DO CFD=1,BSp%NCFnc%I(KD)
          ShellD  = BfnInd%I(AtD,CFD)
          StartLD = BSp%LStrt%I(CFD,KD)
          StopLD  = BSp%LStop%I(CFD,KD)
          StrideD = StopLD-StartLD+1
          LenKet  = DB%DisPtr%I(1,ShellD,iTKet,iCKet)
          iDKet   = DB%DisPtr%I(2,ShellD,iTKet,iCKet)
          iPKet   = DB%DisPtr%I(3,ShellD,iTKet,iCKet)
          IF (LenKet>0) THEN

            CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                           IndexD+1,D%MTrix%D(iPtr),DA%D(1))
            Dcd=GetAbsMax(StrideC*StrideD,DA)   ! could test on the density here...

            IF (DB%DisBuf%D(iDBra+1)<0.0D0) THEN
              BraSwitch=1
              NA=L1
              NC=L2
            ELSE
              BraSwitch=0
              NA=L2
              NC=L1
            END IF
            IF (DB%DisBuf%D(iDKet+1)<0.0D0) THEN
              KetSwitch=1
              NB=L3
              ND=L4
            ELSE 
              KetSwitch=0
              NB=L4
              ND=L3
            END IF
            IntSwitch=10*BraSwitch+KetSwitch

            DO I=1,LenBra
              IBD=iDBra+(I-1)*DB%MAXC
              IBP=iPBra+(I-1)*DB%MAXP*CBra
              IndexA=ABS(DB%DisBuf%D(IBD+1))
              ACx=DB%DisBuf%D(IBD+4)
              ACy=DB%DisBuf%D(IBD+5)
              ACz=DB%DisBuf%D(IBD+6)
              SchB=DB%DisBuf%D(IBD+10)
              Test=Thresholds%TwoE/(Dcd*SchB)
              ISL=0

              DO J=1,LenKet
                IKD=iDKet+(J-1)*DB%MAXC
                IndexB=ABS(DB%DisBuf%D(IKD+1))
                IF (IndexA>=IndexB) THEN     ! Symmetry of the K matrix
                  SchK=DB%DisBuf%D(IKD+10)
                  IF (SchK<=Test) EXIT       ! ONX skipout
                  ISL=ISL+1
                  IKP=iPKet+(J-1)*DB%MAXP*CKet
                  SB%SLDis%I(ISL)=IKD+4
                  SB%SLPrm%I(ISL)=IKP
                END IF
              END DO ! J, LenKet
 
              IF (ISL>0) THEN

                IntSpace=ISL*CBra*CKet*NVRR    ! WRONG HERE
                IF (IntSpace.GT.IB%MAXI) THEN
                  ErrorCode=eMAXI
                  GOTO 9000
                ENDIF

                SELECT CASE (IntCode)
                  CASE (00000000)  
                    CALL Int1111(ISL,CBra,CKet,DB%DisBuf%D(IBD),  &
                                 DB%PrmBuf%D(IBP),DB,IB,SB,IB%W1%D)
!                  CASE (00001000)  
!                    CALL Int1121()
!                  CASE (00001010)  
!                    CALL Int1122()
!                  CASE (10000000)  
!                    CALL Int2111()
!                  CASE (10001000)  
!                    CALL Int2121()
!                  CASE (10001010)  
!                    CALL Int2122()
!                  CASE (10100000)  
!                    CALL Int2211()
!                  CASE (10101000)  
!                    CALL Int2221()
!                  CASE (10101010) 
!                    CALL Int2222()
!                  CASE (00002200)
!                    CALL Int1161()
!                  CASE (00002210)
!                    CALL Int1162()
!                  CASE (00002222)
!                    CALL Int1166()
!                  CASE (10002200)
!                    CALL Int2161()
!                  CASE (10002210)
!                    CALL Int2162()
!                  CASE (10102200)
!                    CALL Int2261()
!                  CASE (22000000)
!                    CALL Int6111()
!                  CASE (22001000)
!                    CALL Int6121()
!                  CASE (22001010)
!                    CALL Int6122()
!                  CASE (22002200)
!                    CALL Int6161()
!                  CASE (22100000)
!                    CALL Int6211()
!                  CASE (22101000)
!                    CALL Int6221()
!                  CASE (22220000)
!                    CALL Int6611()
                  CASE DEFAULT
                    CALL Halt(' Illegal integral type in ONX:ComputeKe')
                END SELECT

                xNERIs=xNERIs+FLOAT(L1*L2*L3*L4*ISL)

                IF (LKet>0) CALL HRRKet(IB%W1%D,DB%DisBuf%D,ISL,               &
                                        SB%SLDis%I,NB1,NB1,TKet)
                IF (LBra>0) THEN 
                  CALL HRRBra(IB%W1%D,IB%W2%D,ACx,ACy,ACz,ISL,                 &
                              NB1,NB2,NK2,TBra)
                  CALL Digest(ISL,NA,NB,NC,ND,L1,L2,L3,L4,                     &
                              IntSwitch,IB%W1%D,IB%W2%D,DA%D)
                  CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB,IB%W1%D,K)
                ELSE
                  CALL Digest(ISL,NA,NB,NC,ND,L1,L2,L3,L4,                     &
                              IntSwitch,IB%W2%D,IB%W1%D,DA%D)
                  CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB,IB%W2%D,K)
                END IF
              END IF  ! ISL

            END DO ! I, LenBra
          END IF ! LenKet
          IndexD=IndexD+StrideD
        END DO ! CFD
      END DO ! ci
          END IF ! LenBra
          IndexC=IndexC+StrideC
        END DO ! CFC
      END DO ! AtC
    END IF ! TCPop on Ket
  END DO ! iCKet
    END IF ! TCPop on Bra
  END DO ! iCBra
  ENDIF  ! LTot
  END DO ! iTKet
  END DO ! iTBra
9000 CONTINUE
  CALL Delete(DA)
END SUBROUTINE ComputeKe
