SUBROUTINE ComputeKe(BSc,GMc,BSp,GMp,D,K,DB1,DB2,IB,SB,IS,Drv,SubInd,BfnInd)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
#ifdef PARALLEL_ONX
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL_ONX
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
  TYPE(DBuf)                :: DB1,DB2   ! ONX distribution buffers
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
  INTEGER               :: I,J,I0,I1,I2,ISL,NInts
  INTEGER               :: IBD,IBP,IKD,IKP
  INTEGER               :: NA,NB,NC,ND
  INTEGER               :: IntSpace,IntCodeV,IntCodeC
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

  DO iTBra=1,DB1%LenTC       ! Loop over angular symmetry types on the Bra
    TBra=DB1%TCode%I(iTBra)
    ITypeA=MOD(TBra,100)
    ITypeC=(TBra-ITypeA)/100
    LBra=LTotal(ITypeA)+LTotal(ITypeC)

  DO iTKet=1,DB2%LenTC       ! Loop over angular symmetry types on the Ket
    TKet=DB2%TCode%I(iTKet)
    ITypeB=MOD(TKet,100)
    ITypeD=(TKet-ITypeB)/100
    LKet=LTotal(ITypeB)+LTotal(ITypeD)
    LTot=LBra+LKet

    CALL GetIntCode(LTot,TBra,TKet,IntCodeV,IntCodeC,Explicit)
    CALL GetIntSpace(TBra,TKet,LBra,LKet,IS)

    IF (Explicit) THEN

    I0=iT(MIN(ITypeC,ITypeA),MAX(ITypeC,ITypeA))
    I1=iT(MIN(ITypeD,ITypeB),MAX(ITypeD,ITypeB))
    I2=I0+(I1-1)*10
    iCP=Drv%CDrv%I(I2)            ! The pointer to the 2e contraction
    iCL=Drv%CDrv%I(iCP)           ! driver

    CALL VRRs(LBra,LKet,Drv)      ! Get the pointers to the VRR table

  DO iCBra=1,DB1%LenCC       ! Loop over contraction lengths on the Bra
    CBra=DB1%CCode%I(iCBra)
    IF (DB1%TCPop%I(iTBra,iCBra)==1) THEN
  DO iCKet=1,DB2%LenCC       ! Loop over contraction lengths on the Ket
    CKet=DB2%CCode%I(iCKet)
    IF (DB2%TCPop%I(iTKet,iCKet)==1) THEN

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
          LenBra  = DB1%DisPtr%I(1,ShellC,iTBra,iCBra)
          iDBra   = DB1%DisPtr%I(2,ShellC,iTBra,iCBra)
          iPBra   = DB1%DisPtr%I(3,ShellC,iTBra,iCBra)
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
          LenKet  = DB2%DisPtr%I(1,ShellD,iTKet,iCKet)
          iDKet   = DB2%DisPtr%I(2,ShellD,iTKet,iCKet)
          iPKet   = DB2%DisPtr%I(3,ShellD,iTKet,iCKet)
          IF (LenKet>0) THEN

            CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                           IndexD+1,D%MTrix%D(iPtr),DA%D(1))
            Dcd=GetAbsMax(StrideC*StrideD,DA)   ! could test on the density here...

            IF (DB1%DisBuf%D(iDBra+1)<0.0D0) THEN
              BraSwitch=1
              NA=IS%L1
              NC=IS%L2
            ELSE
              BraSwitch=0
              NA=IS%L2
              NC=IS%L1
            END IF
            IF (DB2%DisBuf%D(iDKet+1)<0.0D0) THEN
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
              IBD=iDBra+(I-1)*DB1%MAXC
              IBP=iPBra+(I-1)*DB1%MAXP*CBra
              IndexA=ABS(DB1%DisBuf%D(IBD+1))
              ACx=DB1%DisBuf%D(IBD+4)
              ACy=DB1%DisBuf%D(IBD+5)
              ACz=DB1%DisBuf%D(IBD+6)
              SchB=DB1%DisBuf%D(IBD+10)
              Test=Thresholds%TwoE/(Dcd*SchB)
              ISL=0

              DO J=1,LenKet
                IKD=iDKet+(J-1)*DB2%MAXC
                IndexB=ABS(DB2%DisBuf%D(IKD+1))
                IF (IndexA>=IndexB) THEN     ! Symmetry of the K matrix
                  SchK=DB2%DisBuf%D(IKD+10)
                  IF (SchK<=Test) EXIT       ! ONX skipout
                  ISL=ISL+1
                  IKP=iPKet+(J-1)*DB2%MAXP*CKet
                  SB%SLDis%I(ISL)=IKD+4
                  SB%SLPrm%I(ISL)=IKP
                END IF
              END DO ! J, LenKet
 
              IF (ISL>0) THEN

                IntSpace=ISL*IS%NVRR  
                IF (IntSpace.GT.IB%MAXI) THEN
                  ErrorCode=eMAXI
                  GOTO 9000
                ENDIF

                SELECT CASE (IntCodeV)
                  CASE (01010101)  
                    CALL Int1111(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (01010201)  
                    CALL Int1121(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02010101) 
                    CALL Int2111(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (01010202)  
                    CALL Int1122(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02010201)
                    CALL Int2121(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02010202)  
                    CALL Int2122(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02020101)  
                    CALL Int2211(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02020201)  
                    CALL Int2221(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02020202)
                    CALL Int2222(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (01010601)
                    CALL Int1161(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (01010602)
                    CALL Int1162(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (01010606)
                    CALL Int1166(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02010601)
                    CALL Int2161(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02010602)
                    CALL Int2162(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02020601)
                    CALL Int2261(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06010101)
                    CALL Int6111(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06010201)
                    CALL Int6121(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06010202)
                    CALL Int6122(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06010601)
                    CALL Int6161(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06020101)
                    CALL Int6211(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06020201)
                    CALL Int6221(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06060101)
                    CALL Int6611(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))


                  CASE (02010606)
                    CALL Int2166(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02020602)
                    CALL Int2262(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (02020606)
                    CALL Int2266(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06010602)
                    CALL Int6162(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06010606)
                    CALL Int6166(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06020202)
                    CALL Int6222(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06020601)
                    CALL Int6261(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06020602)
                    CALL Int6262(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06020606)
                    CALL Int6266(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06060201)
                    CALL Int6621(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06060202)
                    CALL Int6622(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06060601)
                    CALL Int6661(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06060602)
                    CALL Int6662(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE (06060606)
                    CALL Int6666(ISL,IntCodeC,CBra,CKet,DB1%DisBuf%D(IBD),  &
                                 DB1%PrmBuf%D(IBP),DB2,IB,SB,IB%W1%D(1),IB%W2%D(1))
                  CASE DEFAULT
                    WRITE(*,*) "IntCode=",IntCodeV,IntCodeC
                    CALL Halt(' Illegal integral type in ONX:ComputeKe')
                END SELECT

                xNERIs=xNERIs+FLOAT(IS%L1*IS%L2*IS%L3*IS%L4*ISL)

                IF (LKet>0) CALL HRRKet(IB%W1%D,DB2%DisBuf%D,ISL,               &
                                        SB%SLDis%I,IS%NB1,IS%NB2,IS%NK1,TKet)
                IF (LBra>0) THEN 
                  CALL HRRBra(IB%W1%D(1),IB%W2%D(1),ACx,ACy,ACz,ISL,           &
                              IS%NB1,IS%L1*IS%L2,IS%L3*IS%L4,TBra)
                  CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4,         &
                              IntSwitch,IB%W1%D(1),IB%W2%D(1),DA%D(1))
                  CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB2,IB%W1%D(1),K)
                ELSE
                  CALL Digest(ISL,NA,NB,NC,ND,IS%L1,IS%L2,IS%L3,IS%L4,         &
                              IntSwitch,IB%W2%D(1),IB%W1%D(1),DA%D(1))
                  CALL Scatter(ISL,NA,NB,IndexA,SB,SubInd,DB2,IB%W2%D(1),K)
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
