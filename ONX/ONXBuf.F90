SUBROUTINE ONXBuf(BSc,GMc,BSp,GMp,LenF,NameBuf)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Basis set, coordinates, ect.
!--------------------------------------------------------------------------------
  TYPE(BSET),INTENT(IN)    :: BSc,BSp
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp
  TYPE(BUFL),INTENT(INOUT) :: LenF
  TYPE(INT_VECT)           :: NameBuf
!--------------------------------------------------------------------------------
! Old buffers that may eventually go away...
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)         :: AuxR
  TYPE(DBL_RNK2)         :: SchT,VecT,EBra
  TYPE(DBL_RNK3)         :: TmpBuf,Ctmp,Stmp,RMD
  TYPE(DBL_RNK4)         :: MD
  TYPE(INT_VECT)         :: BfnT
  TYPE(INT_RNK2)         :: BufT
!--------------------------------------------------------------------------------
! Large buffers to hold the sorted distribution and multipole data. 
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)         :: PrmBuf,DisBuf,CBuf,SBuf
  TYPE(DBL_RNK2)         :: VecBuf
  TYPE(INT_VECT)         :: BfnBuf
  TYPE(INT_RNK2)         :: iSA,iSrt
!--------------------------------------------------------------------------------
! Various internal variables
!--------------------------------------------------------------------------------
  REAL(DOUBLE)           :: Test,VSAC
  REAL(DOUBLE)           :: Ax,Ay,Az,Cx,Cy,Cz,ACx,ACy,ACz,AC2
  REAL(DOUBLE)           :: Px,Py,Pz
  REAL(DOUBLE)           :: Za,Zc,Zeta
  REAL(DOUBLE)           :: Cnt,rInt
  INTEGER                :: AtA,CFA,PFA,StartLA,StopLA,StrideA
  INTEGER                :: AtC,CFC,PFC,StartLC,StopLC,StrideC
  INTEGER                :: KA,NBFA,MinLA,MaxLA,IType,NICase 
  INTEGER                :: KC,NBFC,MinLC,MaxLC,KType,NKCase
  INTEGER                :: IndexA,N
  INTEGER                :: CType,LenC
  INTEGER                :: iPrm=1,iDis=1,iMom=1,iCnt1=1,iLen,iB
  INTEGER                :: iBType,i,Lng,iTB,KLng,Ltot,i0,i1,i2
  INTEGER                :: Max4L,MaxSym
  INTEGER                :: BufSize,MXPrm,MXDis,MXMom,MXDs2
#ifdef PARALLEL
#else
  INTEGER                :: MyID=0
#endif

!--------------------------------------------------------------------------------
! Allocate the small buffers
!--------------------------------------------------------------------------------
  Max4L=8         ! Hack for now...
  MaxSym=MAX(BSc%NASym,BSp%NASym)
  CALL New(AuxR,Max4L,0)
  CALL New(SchT,(/500,80/))
  CALL New(VecT,(/3,500/))
  CALL New(EBra,(/Len1,Len0*Len0/))
  CALL New(TmpBuf,(/MXB,MXCNT*MXCNT+MInfo,500/))
  CALL New(Ctmp,(/LenS,Len0*Len0,500/),(/0,1,1/))
  CALL New(Stmp,(/LenS,Len0*Len0,500/),(/0,1,1/))
  CALL New(RMD,(/Max4L,Max4L,Max4L/),(/0,0,0/))
  CALL New(MD,(/3,MaxSym,MaxSym,2*MaxSym/),(/1,0,0,0/))
  CALL New(BufT,(/502,80/))
  CALL New(BfnT,500)
!--------------------------------------------------------------------------------
! Allocate the large buffers
!--------------------------------------------------------------------------------
  BufSize=500000
  MXPrm=BufSize
  MXDis=BufSize/30
  MXMom=BufSize/3
  MXDs2=BufSize/60
  CALL New(PrmBuf,MXPrm)
  CALL New(DisBuf,MXDis)
  CALL New(CBuf,MXMom)
  CALL New(SBuf,MXMom)
  CALL New(VecBuf,(/3,MXDis/))
  CALL New(BfnBuf,MXDis)
  CALL New(iSA,(/5,MXDs2/))
  CALL New(iSrt,(/BSc%NCtrt,NAtoms/))
!--------------------------------------------------------------------------------

  Test=-10.d0*DLOG(Thresholds%Dist)
!  DO AtC=Beg%I(MyId),End%I(MyId)
  DO AtC=1,NAtoms
#ifdef PARALLEL
    IF (NameBuf%I(AtC)==1) THEN
#endif
    KC=GMp%AtTyp%I(AtC)
    NBFC=BSp%BfKnd%I(KC)
    Cx=GMp%Carts%D(1,AtC) 
    Cy=GMp%Carts%D(2,AtC) 
    Cz=GMp%Carts%D(3,AtC) 

    DO CFC=1,BSp%NCFnc%I(KC)
      StartLC=BSp%LStrt%I(CFC,KC)
      StopLC=BSp%LStop%I(CFC,KC)
      MinLC=BSp%ASymm%I(1,CFC,KC)
      MaxLC=BSp%ASymm%I(2,CFC,KC)
      KType=MaxLC*(MaxLC+1)/2+MinLC+1
      NKCase=MaxLC-MinLC+1
      StrideC=StopLC-StartLC+1

      BufT%I=0
      iB=1
      iLen=0
      IndexA=0

      DO AtA=1,NAtoms
        KA=GMc%AtTyp%I(AtA)
        NBFA=BSc%BfKnd%I(KA)
        Ax=GMc%Carts%D(1,AtA)
        Ay=GMc%Carts%D(2,AtA)
        Az=GMc%Carts%D(3,AtA)
        ACx=Ax-Cx
        ACy=Ay-Cy
        ACz=Az-Cz
        AC2=ACx*ACx+ACy*ACy+ACz*ACz
        IF(AC2.GT.Test) THEN
          IndexA=IndexA+BSc%NCFnc%I(KA)
        ELSE
          DO CFA=1,BSc%NCFnc%I(KA)
            IndexA=IndexA+1
            StartLA=BSc%LStrt%I(CFA,KA)
            StopLA=BSc%LStop%I(CFA,KA)
            MinLA=BSc%ASymm%I(1,CFA,KA)
            MaxLA=BSc%ASymm%I(2,CFA,KA)
            IType=MaxLA*(MaxLA+1)/2+MinLA+1
            NICase=MaxLA-MinLA+1
            StrideA=StopLA-StartLA+1
            CType=NICase*10+NKCase
            VecT%D(1,iB)=ACx
            VecT%D(2,iB)=ACy
            VecT%D(3,iB)=ACz
            BfnT%I(iB)=IndexA
            LenC=BSc%NPFnc%I(CFA,KA)*BSp%NPFnc%I(CFC,KC)
            i0=0
            DO PFC=1,BSp%NPFnc%I(CFC,KC)
            DO PFA=1,BSc%NPFnc%I(CFA,KA)
              Za=BSc%Expnt%D(PFA,CFA,KA)
              Zc=BSp%Expnt%D(PFC,CFC,KC)
              Zeta=Za+Zc
              Cnt=BSc%CCoef%D(StopLA,PFA,CFA,KA)*BSp%CCoef%D(StopLC,PFC,CFC,KC)
              Px=(Za*Ax+Zc*Cx)/Zeta
              Py=(Za*Ay+Zc*Cy)/Zeta
              Pz=(Za*Az+Zc*Cz)/Zeta
              VSAC=EXP(-Za*Zc/Zeta*AC2)/Zeta*Prev1*Cnt
              IF(ABS(VSAC).GT.Thresholds%Dist) THEN
                i0=i0+1
                TmpBuf%D(1,i0,iB)=Zeta
                TmpBuf%D(2,i0,iB)=Px
                TmpBuf%D(3,i0,iB)=Py
                TmpBuf%D(4,i0,iB)=Pz
                TmpBuf%D(5,i0,iB)=VSAC
                IF (CType.EQ.11) THEN
                  TmpBuf%D(6,i0,iB)=1.0D0
                ELSE IF (CType.EQ.21.OR.CType.EQ.12) THEN
                  TmpBuf%D(6,i0,iB)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                                    BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
                  TmpBuf%D(7,i0,iB)=1.0D0
                ELSE IF (CType.EQ.22) THEN
                  TmpBuf%D(6,i0,iB)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                                    BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
                  TmpBuf%D(7,i0,iB)=BSc%CCoef%D(StartLA+1,PFA,CFA,KA) * &
                                    BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
                  TmpBuf%D(8,i0,iB)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                                    BSp%CCoef%D(StartLC+1,PFC,CFC,KC) / Cnt
                  TmpBuf%D(9,i0,iB)=1.0D0
                ELSE
                  WRITE(*,*) "CType = ",CType
                  CALL Halt(' Illegal CType in ONXBuf')
                ENDIF 
                TmpBuf%D(10,i0,iB)=EXP(-Za*Zc/Zeta*AC2)*Cnt
              ENDIF
            END DO ! PFA
            END DO ! PFC
            LenC=i0

            rInt=0.0D0
            IF(LenC.gt.0) THEN
              CALL ONXSch(TmpBuf%D(1,1,iB),MaxSym,Max4L,MinLA,MaxLA,            &
                          MinLC,MaxLC,                                          &
                          BSc%LxDex%I(1),BSc%LyDex%I(1),BSc%LzDex%I(1),         &
                          Ax,Ay,Az,Cx,Cy,Cz,AC2,IType,KType,NICase,NKCase,      &
                          LenC,EBra%D(1,1),Ctmp%D(0,1,iB),Stmp%D(0,1,iB),       &
                          MD%D(1,0,0,0),AuxR%D(0),RMD%D(0,0,0),rInt,            &
                          Thresholds%Dist)
            ENDIF

            IF (rInt.GT.Thresholds%Dist) THEN
              DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
              DO i0=1,LenC
                IF (CType.EQ.11) THEN
                  TmpBuf%D(7,i0,iB)=1.0D0
                  TmpBuf%D(8,i0,iB)=1.0D0
                  TmpBuf%D(9,i0,iB)=1.0D0
                ELSE IF (CType.EQ.12.OR.CType.EQ.21) THEN
                  TmpBuf%D(7,i0,iB)=TmpBuf%D(6,i0,iB)
                  TmpBuf%D(8,i0,iB)=TmpBuf%D(6,i0,iB)
                  TmpBuf%D(9,i0,iB)=1.0D0
                END IF
              END DO
              iB=iB+1
              iBType=100*IType+LenC
              DO i=1,80
                IF (iBType.EQ.BufT%I(1,i)) THEN   ! add to an existing type
                  Lng=BufT%I(2,i)+1
                  IF (Lng.GT.500) CALL Halt(' Blown 1 in ONX')
                  BufT%I(2,i)=Lng
                  BufT%I(Lng+2,i)=iB-1
                  SchT%D(Lng,i)=rInt
                  EXIT
                ELSE IF (BufT%I(1,i).EQ.0) THEN   ! create a new type
                  BufT%I(1,i)=iBType
                  Lng=1
                  BufT%I(2,i)=Lng
                  BufT%I(Lng+2,i)=iB-1
                  SchT%D(Lng,i)=rInt
                  iLen=iLen+1
                  EXIT
                END IF
                IF (i.EQ.80)  CALL Halt(' Blown 2 in ONX')
              END DO  
            END IF ! (rInt.GT.Thresholds%Dist)
          END DO ! CFA
        END IF ! (AC2.GT.Test)
      END DO ! AtA

      DO i=1,iLen
        CALL ONXSort2(SchT%D(1,i),BufT%I(3,i),BufT%I(2,i),-2)
      END DO

      ISrt%I(CFC,AtC)=iCnt1
      DO i=1,iLen
        iTB=BufT%I(1,i)
        IType=iTB/100
        KLng=iTB-IType*100
        iSA%I(1,iCnt1)=iTB
        iSA%I(2,iCnt1)=BufT%I(2,I)
        iSA%I(3,iCnt1)=iPrm
        iSA%I(4,iCnt1)=iDis
        iSA%I(5,iCnt1)=iMom
        Ltot=MaxLC+Lmax(IType)
        I1=Ltot*(Ltot+3)/2
        I2=Leng(IType)*Leng(KType)
        N=BufT%I(2,i)
       
        CALL AddDis(N,KLng,I1,I2,iPrm,iDis,iMom,SchT%D(1,i),           &
                    BfnT%I(1),BufT%I(3,i),VecT%D(1,1),TmpBuf%D(1,1,1), &
                    PrmBuf%D(iPrm),DisBuf%D(iDis),BfnBuf%I(iDis),      &
                    VecBuf%D(1,iDis),Ctmp%D(0,1,1),Stmp%D(0,1,1),      &
                    CBuf%D(iMom),SBuf%D(iMom))

        iCnt1=iCnt1+1
      END DO 
   
      iSA%I(1,iCnt1)=-1
      iCnt1=iCnt1+1

    END DO ! CFC
#ifdef PARALLEL
  END IF
#endif
  END DO ! Atc
!--------------------------------------------------------------------------------
! Set the file buffer lengths
!--------------------------------------------------------------------------------
  LenF%PrmL=iPrm
  LenF%DisL=iDis
  LenF%VecL=iDis
  LenF%CosL=iMom
  LenF%SinL=iMom
  LenF%iSAL=iCnt1
  LenF%BfnL=iDis
  LenF%SrtL=BSc%NCtrt*NAtoms
!--------------------------------------------------------------------------------
! Write the buffer arrays to disk.
!--------------------------------------------------------------------------------
!  write(*,*) "iPrm = ",iPrm
!  write(*,*) "iDis = ",iDis
!  write(*,*) "iMom = ",iMom
!  write(*,*) "iCnt1 = ",iCnt1
!  write(*,*) "---------"
  CALL ONXIO(MyID,FWrite,FilePrm,TypeR,LenF%PrmL,PrmBuf%D(1))
  CALL ONXIO(MyID,FWrite,FileDis,TypeR,LenF%DisL,DisBuf%D(1))
  CALL ONXIO(MyID,FWrite,FileVec,TypeR,3*LenF%VecL,VecBuf%D(1,1))
  CALL ONXIO(MyID,FWrite,FileCos,TypeR,LenF%CosL,CBuf%D(1))
  CALL ONXIO(MyID,FWrite,FileSin,TypeR,LenF%SinL,SBuf%D(1))
  CALL ONXIO(MyID,FWrite,FileiSA,TypeI,5*LenF%iSAL,iSA%I(1,1))
  CALL ONXIO(MyID,FWrite,FileBfn,TypeI,LenF%BfnL,BfnBuf%I(1))
  CALL ONXIO(MyID,FWrite,FileSrt,TypeI,LenF%SrtL,iSrt%I(1,1))
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  CALL Delete(AuxR)
  CALL Delete(SchT)
  CALL Delete(VecT)
  CALL Delete(EBra)
  CALL Delete(TmpBuf)
  CALL Delete(Ctmp)
  CALL Delete(Stmp)
  CALL Delete(RMD)
  CALL Delete(MD)
  CALL Delete(BufT)
  CALL Delete(BfnT)
  CALL Delete(PrmBuf)
  CALL Delete(DisBuf)
  CALL Delete(CBuf)
  CALL Delete(SBuf)
  CALL Delete(VecBuf)
  CALL Delete(BfnBuf)
  CALL Delete(iSA)
  CALL Delete(iSrt)

END SUBROUTINE ONXBuf
  
