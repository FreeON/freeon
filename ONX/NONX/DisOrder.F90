SUBROUTINE DisOrder(BSc,GMc,BSp,GMp,NameBuf)  
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE Thresholding
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Basis set, coordinates, ect.
!--------------------------------------------------------------------------------
  TYPE(BSET),INTENT(IN)    :: BSc,BSp
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp
  TYPE(INT_VECT)           :: NameBuf
!--------------------------------------------------------------------------------
! Temporary buffers to hold the distribution data
!--------------------------------------------------------------------------------
  TYPE(DBL_RNK2)           :: TBufC
  TYPE(DBL_RNK3)           :: TBufP
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)           :: W1,W2
  TYPE(DBL_RNK2)           :: CD,WK
  TYPE(DBL_RNK3)           :: WR,WZ
  TYPE(INT_VECT)           :: CDrv,VLoc,MLDis
  TYPE(INT_RNK2)           :: SLoc
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  REAL(DOUBLE)           :: Test,VSAC
  REAL(DOUBLE)           :: ACx,ACy,ACz,AC2,x,y,z
  REAL(DOUBLE)           :: Zeta,Za,Zc,Cnt,rInt

  INTEGER                :: MXINT,MXD,MXC,MXP
  INTEGER                :: LngVRR,LngLoc,LngDrv
  INTEGER                :: id,is,nr,ns
  INTEGER                :: KonAC,LDis,NLOCD2,NLOCD3,NInts
  INTEGER                :: iB,I0,I1,NPrim,iCP,iCL,IKType
  INTEGER                :: AtA,CFA,PFA,StartLA,StopLA,StrideA
  INTEGER                :: AtC,CFC,PFC,StartLC,StopLC,StrideC
  INTEGER                :: KA,NBFA,MinLA,MaxLA,IType,NICase,IndexA
  INTEGER                :: KC,NBFC,MinLC,MaxLC,KType,NKCase,IndexC
  INTEGER                :: CType,II,JJ,IJ
!--------------------------------------------------------------------------------
! F77 function calls
!--------------------------------------------------------------------------------
  INTEGER                :: IDmn,NFinal,iT

  NPrim=MAX(BSc%NPrim,BSp%NPrim)

  MXINT=10000
  MXD=1000
  MXC=11
  MXP=8

  1000 CONTINUE

  CALL New(WK,(/4,NPrim**2/))
  CALL New(CD,(/NPrim**2,3/))
  CALL New(WR,(/6,NPrim**2,NPrim**2/))
  CALL New(WZ,(/3,NPrim**2,NPrim**2/))
  CALL New(W1,MXINT)
  CALL New(W2,MXINT)
  CALL New(CDrv,60000)
  CALL New(TBufC,(/MXC,MXD/))
  CALL New(TBufP,(/MXP,NPrim*NPrim+MInfo,MXD/))
  CALL New(MLDis,1)

  CALL CCDriver(CDrv%I(1),LngDrv)
  CALL VRRLng(LngVRR,LngLoc)
  CALL New(VLoc,LngVRR)
  CALL New(SLoc,(/3,LngLoc/))
  CALL VRRDriver(VLOC%I(1),SLOC%I(1,1),LngVRR,LngLoc)

  Test=-10.d0*DLOG(Thresholds%Dist) 
  MLDis%I(1)=5
  IndexC=0
  DO AtC=1,NAtoms
#ifdef PARALLEL
    IF (NameBuf%I(AtC)==1) THEN
#endif
    KC=GMp%AtTyp%I(AtC)
    NBFC=BSp%BfKnd%I(KC)
    DO CFC=1,BSp%NCFnc%I(KC)
      IndexC=IndexC+1
      StartLC=BSp%LStrt%I(CFC,KC)
      StopLC=BSp%LStop%I(CFC,KC)
      MinLC=BSp%ASymm%I(1,CFC,KC)
      MaxLC=BSp%ASymm%I(2,CFC,KC)
      KType=MaxLC*(MaxLC+1)/2+MinLC+1
      NKCase=MaxLC-MinLC+1
      StrideC=StopLC-StartLC+1         

      iB=1 
      IndexA=0  

  DO AtA=1,NAtoms
    KA=GMc%AtTyp%I(AtA)
    NBFA=BSc%BfKnd%I(KA)
    ACx=GMc%Carts%D(1,AtA)-GMp%Carts%D(1,AtC)
    ACy=GMc%Carts%D(2,AtA)-GMp%Carts%D(2,AtC) 
    ACz=GMc%Carts%D(3,AtA)-GMp%Carts%D(3,AtC) 
    AC2=ACx*ACx+ACy*ACy+ACz*ACz    
    IF (AC2.GT.Test) THEN
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
        II=MAX(IndexA,IndexC)
        JJ=MIN(IndexA,IndexC)
        IJ=II*(II-1)/2+JJ 
        TBufC%D( 1,iB)=DFLOAT(IJ)+Small
        TBufC%D( 2,iB)=DFLOAT(IndexA)+Small
        IF(IType.GE.KType) THEN
          IKType=IType*100+KType
          TBufC%D( 3,iB)=DFLOAT(IndexA)+Small
          TBufC%D( 4,iB)=DFLOAT(IndexC)+Small
          TBufC%D( 5,iB)=ACx
          TBufC%D( 6,iB)=ACy
          TBufC%D( 7,iB)=ACz
          TBufC%D( 8,iB)=GMc%Carts%D(1,AtA)
          TBufC%D( 9,iB)=GMc%Carts%D(2,AtA)
          TBufC%D(10,iB)=GMc%Carts%D(3,AtA)
          x=ACx
          y=ACy
          z=ACz
        ELSE 
          IKType=KType*100+IType
          TBufC%D( 3,iB)=DFLOAT(IndexC)+Small
          TBufC%D( 4,iB)=DFLOAT(IndexA)+Small
          TBufC%D( 5,iB)=-ACx
          TBufC%D( 6,iB)=-ACy
          TBufC%D( 7,iB)=-ACz
          TBufC%D( 8,iB)=GMp%Carts%D(1,AtC)
          TBufC%D( 9,iB)=GMp%Carts%D(2,AtC)
          TBufC%D(10,iB)=GMp%Carts%D(3,AtC)
          x=-ACx
          y=-ACy
          z=-ACz
        ENDIF

        I0=0
        DO PFC=1,BSp%NPFnc%I(CFC,KC)
        DO PFA=1,BSc%NPFnc%I(CFA,KA)
          Za=BSc%Expnt%D(PFA,CFA,KA)
          Zc=BSp%Expnt%D(PFC,CFC,KC)
          Zeta=Za+Zc
          Cnt=BSc%CCoef%D(StopLA,PFA,CFA,KA)*BSp%CCoef%D(StopLC,PFC,CFC,KC)
          VSAC=EXP(-Za*Zc/Zeta*AC2)/Zeta*Prev1*Cnt
          IF (ABS(VSAC).GT.Thresholds%Dist*10D-2) THEN
            I0=I0+1
            TBufP%D(1,I0,iB)=Zeta
            TBufP%D(2,I0,iB)=(Za*GMc%Carts%D(1,AtA)+Zc*GMp%Carts%D(1,AtC))/Zeta
            TBufP%D(3,I0,iB)=(Za*GMc%Carts%D(2,AtA)+Zc*GMp%Carts%D(2,AtC))/Zeta
            TBufP%D(4,I0,iB)=(Za*GMc%Carts%D(3,AtA)+Zc*GMp%Carts%D(3,AtC))/Zeta
            TBufP%D(5,I0,iB)=VSAC
            IF (CType.EQ.11) THEN
              TBufP%D(6,I0,iB)=1.0D0
              TBufP%D(7,I0,iB)=1.0D0
              TBufP%D(8,I0,iB)=1.0D0
            ELSEIF (CType.EQ.12.OR.CType.EQ.21) THEN
              TBufP%D(6,I0,iB)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              TBufP%D(7,I0,iB)=TBufP%D(6,I0,iB)
              TBufP%D(8,I0,iB)=TBufP%D(6,I0,iB)
            ELSEIF (CType.EQ.22) THEN
              TBufP%D(6,I0,iB)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              TBufP%D(7,I0,iB)=BSc%CCoef%D(StartLA+1,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              TBufP%D(8,I0,iB)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC+1,PFC,CFC,KC) / Cnt
            ELSE
              WRITE(*,*) 'CType=',CType
              CALL Halt(' Illegal CType in ONX:DisOrder')
            END IF ! CType
          END IF ! test on VSAC
        ENDDO ! PFA
        ENDDO ! PFC

        IF (I0.GT.0) THEN
          KonAC=I0
          LDis=MaxLA+MaxLC
          NLOCD2=IDmn(LDis)
          NLOCD3=NFinal(IType)*NFinal(KType)
          NInts=NLOCD3*NLOCD3
          I0=iT(MIN(IType,KType),MAX(IType,KType))
          I1=11*I0-10
          iCP=CDrv%I(I1)
          iCL=CDrv%I(iCP)
          
          CALL RGen1C(2*LDis,KonAC,WK%D(1,1),CD%D(1,1),WR%D(1,1),  &
                      WZ%D(1,1),W1%D(1),TBufP%D(1,1,iB),TBufC%D(1,iB))

          id = LDis+1+LDis*5
          is = SLoc%I(1,id)
          nr = SLoc%I(2,id)
          ns = SLoc%I(3,id)

          CALL VRRl(KonAC,nr,ns,VLoc%I(is),VLoc%I(is+nr),W2%D(1),W1%D(1), &
                    WR%D(1,1),WZ%D(1,1))

          CALL Contract(1,KonAC,KonAC,iCL,CDrv%I(iCP+1),CD,CD,W1,W2)

          IF(LDis.NE.0) THEN
            CALL HRRKetOld(W1%D(1),TBufC%D(1,iB),1,MLDis%I(1),NLOCD2,MXC,IKType)
            CALL HRRBraOld(W1%D(1),W2%D(1),x,y,z,1,NLOCD2,NLOCD3,NLOCD3,IKType)
            CALL GetBigInt(NInts,rInt,W2%D(1))
          ELSE
            CALL GetBigInt(NInts,rInt,W1%D(1))
          ENDIF

          write(*,*) "rint = ",rInt

        END IF ! I0
      ENDDO ! CFA
    END IF ! test on AC2
  ENDDO ! AtA

    END DO ! CFC
#ifdef PARALLEL
  END IF
#endif
  END DO ! Atc


  CALL Delete(SLOC)
  CALL Delete(VLOC)
  CALL Delete(MLDis)
  CALL Delete(TBufP)
  CALL Delete(TBufC)
  CALL Delete(CDrv)
  CALL Delete(W2)
  CALL Delete(W1)
  CALL Delete(WZ)
  CALL Delete(WR)
  CALL Delete(CD)
  CALL Delete(WK)

END SUBROUTINE DisOrder
 
