SUBROUTINE ComputeKg(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,IS,Drv,SubInd,BfnInd)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
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
  TYPE(ISpc)                :: IS        ! Array dimms for 2-e routines
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
  INTEGER               :: I,J,I0,I1,I2,ISL
  INTEGER               :: IBD,IBP,IKD,IKP
  INTEGER               :: NA,NB,NC,ND
  INTEGER               :: IntSpace,IntCodeC,IntCodeV
  INTEGER               :: BraSwitch,KetSwitch,IntSwitch
  REAL(DOUBLE)          :: Dcd,SchB,SchK,Test
  REAL(DOUBLE)          :: ACx,ACy,ACz
  TYPE(DBL_VECT)        :: DA
  LOGICAL               :: Explicit
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
  INTEGER               :: LTotal,MaxBatchSize,NFinal,iT

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

    CALL GetIntCode(LTot,TBra,TKet,IntCodeV,IntCodeC,Explicit)
    CALL GetIntSpace(TBra,TKet,LBra,LKet,IS)

    IF (.NOT.Explicit) THEN

    I0=iT(MIN(ITypeC,ITypeA),MAX(ITypeC,ITypeA))
    I1=iT(MIN(ITypeD,ITypeB),MAX(ITypeD,ITypeB))
    I2=I0+(I1-1)*10

    iCP=Drv%CDrv%I(I2)            ! The pointer to the 2e contraction
    iCL=Drv%CDrv%I(iCP)           ! driver

    CALL VRRs(LBra,LKet,Drv)      ! Get the pointers to the VRR table

  DO iCBra=1,DB%LenCC       ! Loop over contraction lengths on the Bra
    CBra=DB%CCode%I(iCBra)
    IF (DB%TCPop%I(iTBra,iCBra)==1) THEN
  DO iCKet=1,DB%LenCC       ! Loop over contraction lengths on the Ket
    CKet=DB%CCode%I(iCKet)
    IF (DB%TCPop%I(iTKet,iCKet)==1) THEN

!     MaxInts=MaxBatchSize(iTBra,iCBra,iTKet,iCBra,DB)

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
              NA=IS%L1
              NC=IS%L2
            ELSE
              BraSwitch=0
              NA=IS%L2
              NC=IS%L1
            END IF
            IF (DB%DisBuf%D(iDKet+1)<0.0D0) THEN
              KetSwitch=1
              NB=IS%L3
              ND=IS%L4
            ELSE 
              KetSwitch=0
              NB=IS%L4
              ND=IS%L3
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
                IntSpace=ISL*CBra*CKet*IS%NVRR
                IF (IntSpace.GT.IB%MAXI) THEN
                  ErrorCode=eMAXI
                  GOTO 9000
                ENDIF
                xNERIs=xNERIs+FLOAT(IS%L1*IS%L2*IS%L3*IS%L4*ISL)
                CALL RGen(ISL,Ltot,CBra,CKet,IB%CB%D(1,1),IB%CK%D(1,1,1),      &
                          DB%DisBuf%D(IBD),DB%PrmBuf%D(IBP),IB%W1%D(1),        &
                          DB,IB,SB)
                CALL VRRl(ISL*CBra*CKet,IS%NVRR,Drv%nr,Drv%ns,                 &
                          Drv%VLOC%I(Drv%is),                                  &
                          Drv%VLOC%I(Drv%is+Drv%nr),IB,                        &
                          IB%W2%D(1),IB%W1%D(1))
                CALL Contract(ISL,CBra,CKet,IS%NVRR,iCL,Drv%CDrv%I(iCP+1),     &
                              IB%CB%D(1,1),IB%CK%D,IB%W1%D(1),IB%W2%D(1))
                IF (LKet>0) CALL HRRKet(IB%W1%D,DB%DisBuf%D,ISL,               &
                                        SB%SLDis%I,IS%NB1,IS%NB2,IS%NK1,TKet)
                IF (LBra>0) THEN
                  CALL HRRBra(IB%W1%D(1),IB%W2%D(1),ACx,ACy,ACz,ISL,           &
                              IS%NB1,IS%L1*IS%L2,IS%L3*IS%L4,TBra)
                  CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4,         &
                              IntSwitch,IB%W1%D(1),IB%W2%D(1),DA%D(1))
                  CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB,IB%W1%D(1),K)
                ELSE
                  CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4,         &
                              IntSwitch,IB%W2%D(1),IB%W1%D(1),DA%D(1))
                  CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB,IB%W2%D(1),K)
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
END SUBROUTINE ComputeKg
