SUBROUTINE ComputeXForce(BS,GM,D,XFrc,DB,IB,SB,IS,Drv,SubInd,BfnInd)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  IMPLICIT NONE
  TYPE(BCSR)                :: D
  TYPE(DBL_RNK2)            :: XFrc
  TYPE(BSET),INTENT(IN)     :: BS        ! basis set info
  TYPE(CRDS),INTENT(IN)     :: GM        ! geometry info
  TYPE(DBuf)                :: DB        ! ONX distribution buffers
  TYPE(IBuf)                :: IB        ! ONX 2-e eval buffers
  TYPE(DSL)                 :: SB        ! ONX distribution pointers
  TYPE(ISpc)                :: IS        ! Array dimms for 2-e routines
  TYPE(IDrv)                :: Drv       ! VRR/contraction drivers
  TYPE(GradD)               :: GD        ! Gradient drivers
  TYPE(INT_RNK2)            :: SubInd    ! Index -> BCSR converter
  TYPE(INT_RNK2),INTENT(IN) :: BfnInd
!-------------------------------------------------------------------
! Misc. internal variables
!-------------------------------------------------------------------
  INTEGER               :: iTBra,iCBra,TBra,ITypeC,ITypeA,CBra
  INTEGER               :: iTKet,iCKet,TKet,ITypeD,ITypeB,CKet
  INTEGER               :: LenBra,iDBra,iPBra,IBra
  INTEGER               :: LenKet,iDKet,iPKet,IKet
  INTEGER               :: Ltot,LtotG,LBra,LBraG,LKet,LKetG
  INTEGER               :: ShellC,ShellD,IndexA,IndexB
  INTEGER               :: AtC,KC,CFC,StartLC,StopLC,StrideC,IndexC,NBFC
  INTEGER               :: AtD,KD,CFD,StartLD,StopLD,StrideD,IndexD,NBFD
  INTEGER               :: AtA,NBFA,RS
  INTEGER               :: ri,ci,iPtr
  INTEGER               :: I,J,I0,I1,I2,ISL,NInts
  INTEGER               :: IBD,IBP,IKD,IKP
  INTEGER               :: NB1,NB2,NK1,NK2
  INTEGER               :: NT,N1,N2,N3,N4
  INTEGER               :: NA,NB,NC,ND
  INTEGER               :: L1,L2,L3,L4
  INTEGER               :: BraSwitch,KetSwitch,IntSwitch
  INTEGER               :: MaxInts
  REAL(DOUBLE)          :: Dmax,SchB,SchK,Test
  REAL(DOUBLE)          :: ACx,ACy,ACz
  TYPE(DBL_VECT)        :: Dcd,Dab
  TYPE(INT_VECT)        :: NTmp
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
  INTEGER               :: LTotal,MaxBatchSize,NFinal,iT

  MaxInts = 500 ! This is a hack...
  CALL New(Dcd,BS%LMNLen*BS%LMNLen)
  CALL New(Dab,MaxInts*BS%LMNLen*BS%LMNLen)  ! THIS NEEDS TO BE FIXED
  CALL New(GD)
  CALL New(NTmp,MaxInts)
!
! Need to call a routine here to sort the density matrix -> for skipout loops
! Could sort the density matrix but this will complicate the binary search 
! rountine later on when the integrals are digested... If we sort into sub-blocks
! corresponding to integral types and contraction lengths then this will work.
! This is not true above. This is actually a very big problem with onx gradients...
! It may be that sonx gradients are the only way to do this...
!
  DO iTBra=1,DB%LenTC       ! Loop over angular symmetry types on the Bra
    TBra=DB%TCode%I(iTBra)
    ITypeA=MOD(TBra,100)
    ITypeC=(TBra-ITypeA)/100
    LBra=LTotal(ITypeA)+LTotal(ITypeC)
    LBraG=LBraG+1

  DO iTKet=1,DB%LenTC       ! Loop over angular symmetry types on the Ket
    TKet=DB%TCode%I(iTKet)
    ITypeB=MOD(TKet,100)
    ITypeD=(TKet-ITypeB)/100
    LKet=LTotal(ITypeB)+LTotal(ITypeD)
    LKetG=LKet+1

    NB1=IDmn(LBra)
    NK1=IDmn(LKet)
    L1=NFinal(ITypeC)  ! Is NFinal equal to the stride????
    L2=NFinal(ITypeA)
    L3=NFinal(ITypeD)
    L4=NFinal(ITypeB)
    NB2=L1*L2
    NK2=L3*L4
    NInts=NB2*NK2
    
    I0=iT(MIN(ITypeC,ITypeA),MAX(ITypeC,ITypeA))
    I1=iT(MIN(ITypeD,ITypeB),MAX(ITypeD,ITypeB))
    I2=I0+(I1-1)*10
    Ltot=LBra+LKet
    LtotG=LBraG+LKetG

    CALL GetIntSpace(TBra,TKet,LBraG,LKetG,IS)

    write(*,*) "Getting the VRR table for LBra=",LBraG," and LKet=",LKetG

    CALL VRRs(LBraG,LKetG,Drv)      ! Get the pointers to the VRR table

    write(*,*) "Getting gradient drivers, TBra=",TBra," TKet=",TKet

    CALL GDrivers(TBra,TKet,GD)     ! Get the gradient driver files

!    GD%NCON   = GD%GDrv1%I(GD%LG1-3)
    GD%NLOCB1 = (LBra+2)*(LBra+3)*(LBra+4)/6
    GD%NLOCB2 = IDmn(LBraG)
    GD%NLOCB3 = L1*L2
    GD%NLOCK2 = IDmn(LKetG)
    GD%NLOCK3 = L3*L4
!
! NVRR is not correct here. It is larger than what you get in 
! rgen.... For now fix this by computing more interals than needed
! in the VRR step.
!
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
        KC=GM%AtTyp%I(AtC)
        NBFC=BS%BfKnd%I(KC)
        IndexC=0
        DO CFC=1,BS%NCFnc%I(KC)
          ShellC  = BfnInd%I(AtC,CFC)
          StartLC = BS%LStrt%I(CFC,KC)
          StopLC  = BS%LStop%I(CFC,KC)
          StrideC = StopLC-StartLC+1
          LenBra  = DB%DisPtr%I(1,ShellC,iTBra,iCBra)
          iDBra   = DB%DisPtr%I(2,ShellC,iTBra,iCBra)
          iPBra   = DB%DisPtr%I(3,ShellC,iTBra,iCBra)
          IF (LenBra>0) THEN

      ShellD=0
      DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1        ! Loop over atom D (Dcd)
        AtD    = D%ColPt%I(ci)
        iPtr   = D%BlkPt%I(ci)
        KD     = GM%AtTyp%I(AtD)
        NBFD   = BS%BfKnd%I(KD)
        IndexD = 0
        DO CFD=1,BS%NCFnc%I(KD)
          ShellD  = BfnInd%I(AtD,CFD)
          StartLD = BS%LStrt%I(CFD,KD)
          StopLD  = BS%LStop%I(CFD,KD)
          StrideD = StopLD-StartLD+1
          LenKet  = DB%DisPtr%I(1,ShellD,iTKet,iCKet)
          iDKet   = DB%DisPtr%I(2,ShellD,iTKet,iCKet)
          iPKet   = DB%DisPtr%I(3,ShellD,iTKet,iCKet)
          IF (LenKet>0) THEN

            CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                           IndexD+1,D%MTrix%D(iPtr),Dcd%D(1))
            Dmax=GetAbsMax(StrideC*StrideD,Dcd)   

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
              AtA  = SubInd%I(1,IndexA)
              NBFA = SubInd%I(2,IndexA)
              RS   = SubInd%I(3,IndexA)
              ACx=DB%DisBuf%D(IBD+4)
              ACy=DB%DisBuf%D(IBD+5)
              ACz=DB%DisBuf%D(IBD+6)
              SchB=DB%DisBuf%D(IBD+10)
              Test=Thresholds%TwoE/(Dmax*SchB)
              ISL=0
              DO J=1,LenKet
                IKD=iDKet+(J-1)*DB%MAXC
                IndexB=ABS(DB%DisBuf%D(IKD+1))
!                IF (IndexA>=IndexB) THEN     ! Symmetry of the K matrix THIS IS WRONG?
                  SchK=DB%DisBuf%D(IKD+10)
!                  IF (SchK<=Test) EXIT       ! ONX skipout
                  ISL=ISL+1
                  IKP=iPKet+(J-1)*DB%MAXP*CKet
                  SB%SLDis%I(ISL)=IKD+4
                  SB%SLPrm%I(ISL)=IKP
!                END IF
              END DO ! J, LenKet

              IF (ISL>0) THEN

       write(*,*) "ISL=",ISL

                CALL RGen(ISL,LtotG,CBra,CKet,IB%CB%D(1,1),IB%CK%D(1,1,1), &
                          DB%DisBuf%D(IBD),DB%PrmBuf%D(IBP),IB%W1%D(1),    &
                          DB,IB,SB)

                CALL VRRl(ISL*CBra*CKet,IS%NVRR,Drv%nr,Drv%ns,                &
                          Drv%VLOC%I(Drv%is),                              &
                          Drv%VLOC%I(Drv%is+Drv%nr),IB,                    &
                          IB%W2%D(1),IB%W1%D(1))


                CALL ContractG(ISL,CBra,CKet,IS%NVRR,IB%CB%D,IB%CK%D,IB%W1%D, &
                               IB%W2%D,DB%PrmBuf%D(IBP),DB,SB,GD)


                DO IKet=1,GD%LG2
                  NT=GD%GDrv2%I(1,IKet)
                  N1=GD%GDrv2%I(2,IKet)
                  N2=GD%GDrv2%I(3,IKet)
                  N3=ISL*(GD%GDrv2%I(4,IKet)-1)+1
                  N4=GD%GDrv2%I(5,IKet)

     write(*,*) "Ket=",nt,n1,n2,n3,n4

                  CALL HrrKet(IB%W1%D(N3),DB%DisBuf%D,ISL,    &
                              SB%SLDis%I,N1,N1,N2,NT)


!
!   CALL HRRKet(IB%W1%D,DB%DisBuf%D,ISL,SB%SLDis%I,IS%NB1,IS%NB2,IS%NK1,TKet)
!   CALL HRRKetGrad(IB%W1%D(N3),DB%DisBuf%D,ISL,SB%SLDis%I,N1,NT)
!
                END DO

                DO IBra=1,GD%LG3
                  NT=GD%GDrv3%I(1,IBra)
                  N1=GD%GDrv3%I(2,IBra)
                  N2=GD%GDrv3%I(3,IBra)
                  N3=ISL*(GD%GDrv3%I(4,IBra)-1)+1
                  N4=GD%GDrv3%I(5,IBra)

     write(*,*) "Bra=",nt,n1,n2,n3,n4

                  CALL HrrBra(IB%W1%D(N3),IB%W2%D(N3),ACx,ACy,ACz,ISL,N1,N4,N2,NT)
!     call halt('enough')
!
!   CALL HRRBraGrad(IB%W1%D(N3),ACx,ACy,ACz,ISL,N1,N2,NT)
!
                END DO

                CALL GetGradient(ISL,GD,IB%W1%D,IB%W2%D)

!                write(*,*) "Calling DigestGradient"
                CALL DigestGradient(ISL,NA,NB,NC,ND,L1,L2,L3,L4,IntSwitch,  &
                                    AtA,AtC,AtD,NBFA,RS,SB,DB,D,SubInd,     &
                                    NTmp,Dcd%D,Dab%D,XFrc,IB%W1%D)

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
  END DO ! iTKet
  END DO ! iTBra

  CALL Delete(Dcd)
  CALL Delete(Dab)
  CALL Delete(GD)
  CALL Delete(NTmp)
END SUBROUTINE ComputeXForce
