SUBROUTINE DisOrder(BSc,GMc,BSp,GMp,DB,IB,SB,Drv,NameBuf)  
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Basis set, coordinates, ect...
!--------------------------------------------------------------------------------
  TYPE(BSET),INTENT(IN)    :: BSc,BSp   ! basis set info
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp   ! geometry info
  TYPE(DBuf)               :: DB        ! ONX distribution buffers
  TYPE(IBuf)               :: IB        ! ONX 2-e eval buffers
  TYPE(DSL)                :: SB        ! ONX distribution pointers
  TYPE(IDrv)               :: Drv       ! VRR/contraction drivers
  TYPE(INT_VECT)           :: NameBuf   ! for parallel implementation
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  TYPE(DBL_RNK3)         :: SchT
  TYPE(INT_RNK3)         :: BufT
  TYPE(INT_RNK2)         :: BufN
  REAL(DOUBLE)           :: Test,VSAC
  REAL(DOUBLE)           :: ACx,ACy,ACz,AC2,x,y,z
  REAL(DOUBLE)           :: Zeta,Za,Zc,Cnt,rInt,XiAB
  INTEGER                :: Lng,i,j,n
  INTEGER                :: KonAC,LDis,NLOCD2,NLOCD3,NInts,NVRR
  INTEGER                :: iBf,I0,I1,iCP,iCL,IKType
  INTEGER                :: AtA,CFA,PFA,StartLA,StopLA,StrideA
  INTEGER                :: AtC,CFC,PFC,StartLC,StopLC,StrideC
  INTEGER                :: KA,NBFA,MinLA,MaxLA,IType,NICase,IndexA
  INTEGER                :: KC,NBFC,MinLC,MaxLC,KType,NKCase,IndexC
  INTEGER                :: CType,II,JJ,IJ
  INTEGER                :: iKonAC,iIKType
  INTEGER                :: iDis,iPrm, itmp
  LOGICAL                :: Found
!--------------------------------------------------------------------------------
! Function calls
!--------------------------------------------------------------------------------
  INTEGER                :: NFinal,iT

  CALL New(BufN,(/DB%MAXT,DB%NPrim*DB%NPrim/))
  CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
  CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
  DB%LenCC=0
  DB%LenTC=0
  iDis=1
  iPrm=1
  Test=-10.d0*DLOG(Thresholds%Dist) 
  SB%SLDis%I(1)=5
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

      BufN%I=0
      iBf=1 
      IndexA=0  

  DO AtA=1,NAtoms
    KA=GMc%AtTyp%I(AtA)
    NBFA=BSc%BfKnd%I(KA)
    ACx=GMc%Carts%D(1,AtA)-GMp%Carts%D(1,AtC)
    ACy=GMc%Carts%D(2,AtA)-GMp%Carts%D(2,AtC) 
    ACz=GMc%Carts%D(3,AtA)-GMp%Carts%D(3,AtC) 
    AC2=ACx*ACx+ACy*ACy+ACz*ACz    
!
! Need to call something like SetAtomPair here, the problem 
! is that it assumes a single basis set (as in J) but I need 
! dual basis sets for K.
!
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
        IF(IType.GE.KType) THEN
          IKType=IType*100+KType
          DB%TBufC%D( 1,iBf)=DFLOAT(IJ)+Small
          DB%TBufC%D( 2,iBf)=-DFLOAT(IndexA)-Small
          DB%TBufC%D( 3,iBf)=DFLOAT(IndexC)+Small
          DB%TBufC%D( 4,iBf)=DFLOAT(IndexA)+Small
          DB%TBufC%D( 5,iBf)=ACx
          DB%TBufC%D( 6,iBf)=ACy
          DB%TBufC%D( 7,iBf)=ACz
          DB%TBufC%D( 8,iBf)=GMc%Carts%D(1,AtA)
          DB%TBufC%D( 9,iBf)=GMc%Carts%D(2,AtA)
          DB%TBufC%D(10,iBf)=GMc%Carts%D(3,AtA)
          x=ACx
          y=ACy
          z=ACz
        ELSE 
          IKType=KType*100+IType
          DB%TBufC%D( 1,iBf)=DFLOAT(IJ)+Small
          DB%TBufC%D( 2,iBf)=DFLOAT(IndexA)+Small
          DB%TBufC%D( 3,iBf)=DFLOAT(IndexC)+Small
          DB%TBufC%D( 4,iBf)=DFLOAT(IndexA)+Small
          DB%TBufC%D( 5,iBf)=-ACx
          DB%TBufC%D( 6,iBf)=-ACy
          DB%TBufC%D( 7,iBf)=-ACz
          DB%TBufC%D( 8,iBf)=GMp%Carts%D(1,AtC)
          DB%TBufC%D( 9,iBf)=GMp%Carts%D(2,AtC)
          DB%TBufC%D(10,iBf)=GMp%Carts%D(3,AtC)
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
          XiAB=Za*Zc/Zeta
          IF (TestPrimPair(XiAB,AC2)) THEN
            I0=I0+1
            DB%TBufP%D(1,I0,iBf)=Zeta
            DB%TBufP%D(2,I0,iBf)=(Za*GMc%Carts%D(1,AtA)+Zc*GMp%Carts%D(1,AtC))/Zeta
            DB%TBufP%D(3,I0,iBf)=(Za*GMc%Carts%D(2,AtA)+Zc*GMp%Carts%D(2,AtC))/Zeta
            DB%TBufP%D(4,I0,iBf)=(Za*GMc%Carts%D(3,AtA)+Zc*GMp%Carts%D(3,AtC))/Zeta
            DB%TBufP%D(5,I0,iBf)=VSAC
            IF (CType.EQ.11) THEN
              DB%TBufP%D(6,I0,iBf)=1.0D0
              DB%TBufP%D(7,I0,iBf)=1.0D0
              DB%TBufP%D(8,I0,iBf)=1.0D0
            ELSEIF (CType.EQ.12.OR.CType.EQ.21) THEN
              DB%TBufP%D(6,I0,iBf)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                                   BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              DB%TBufP%D(7,I0,iBf)=DB%TBufP%D(6,I0,iBf)
              DB%TBufP%D(8,I0,iBf)=DB%TBufP%D(6,I0,iBf)
            ELSEIF (CType.EQ.22) THEN
              DB%TBufP%D(6,I0,iBf)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                                   BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              DB%TBufP%D(7,I0,iBf)=BSc%CCoef%D(StartLA+1,PFA,CFA,KA) * &
                                   BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              DB%TBufP%D(8,I0,iBf)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                                   BSp%CCoef%D(StartLC+1,PFC,CFC,KC) / Cnt
            ELSE
              WRITE(*,*) 'CType=',CType
              CALL Halt(' Illegal CType in ONX:DisOrder')
            END IF ! CType
          END IF ! test on primitive distributions
        ENDDO ! PFA
        ENDDO ! PFC

        IF (I0>0) THEN
          KonAC=I0
          LDis=MaxLA+MaxLC
          NLOCD2=IDmn(LDis)
          NLOCD3=NFinal(IType)*NFinal(KType)
          NVRR =NLOCD2*NLOCD2
          NInts=NLOCD3*NLOCD3

          I0=iT(MIN(IType,KType),MAX(IType,KType))
          I1=11*I0-10
          iCP=Drv%CDrv%I(I1)
          iCL=Drv%CDrv%I(iCP)

          IF (KonAC*KonAC*NVRR>IB%MAXI.OR.NInts>IB%MAXI) THEN
            ErrorCode=eMAXI
            GOTO 9000
          ENDIF

          CALL RGen1C(2*LDis,iBf,KonAC,IB%CB%D,IB%WR,IB%WZ, &
                      IB%W1%D,DB%TBufP,DB%TBufC)
          CALL VRRs(LDis,LDis,Drv)

          CALL VRRl(KonAC*KonAC,NVRR,Drv%nr,Drv%ns,Drv%VLOC%I(Drv%is),  &
                    Drv%VLOC%I(Drv%is+Drv%nr),                          &
                    IB%W2%D,IB%W1%D,IB%WR%D,IB%WZ%D)

          CALL Contract(1,KonAC,KonAC,NVRR,iCL,Drv%CDrv%I(iCP+1), &
                        IB%CB%D,IB%CB%D,IB%W1%D,IB%W2%D)

          IF(LDis.NE.0) THEN
            CALL HRRKet(IB%W1%D,DB%TBufC%D(1,iBf),1,SB%SLDis%I,NLOCD2,NLOCD2,IKType)
            CALL HRRBra(IB%W1%D,IB%W2%D,x,y,z,1,NLOCD2,NLOCD3,NLOCD3,IKType)
            rInt = DSQRT(GetAbsMax(NInts,IB%W2))
          ELSE
            rInt = DSQRT(GetAbsMax(NInts,IB%W1))
          ENDIF

          IF(rInt>Thresholds%Dist) THEN 
            DisRange=MAX(DisRange,SQRT(AC2)*1.01D0)
            DB%TBufC%D(11,iBf)=rInt
            Found=.FALSE.
            DO I=1,DB%LenCC       ! look for the contraction type
              IF (KonAC==DB%CCode%I(I)) THEN
                iKonAC=I
                Found=.TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT.Found) THEN  ! if not found make a new CC type
              DB%LenCC=DB%LenCC+1
              iKonAC=DB%LenCC
              IF (DB%LenCC>DB%MAXK) THEN
                ErrorCode=eMAXK
                GOTO 9000
              END IF
              DB%CCode%I(DB%LenCC)=KonAC
            END IF
            Found=.FALSE.
            DO I=1,DB%LenTC       ! look for the angular symmetry type
              IF (IKType==DB%TCode%I(I)) THEN
                iIKType=I
                Found=.TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT.Found) THEN  ! if not found make a new L type
              DB%LenTC=DB%LenTC+1
              iIKType=DB%LenTC
              IF (DB%LenTC.GT.DB%MAXT) THEN
                ErrorCode=eMAXT
                GOTO 9000
              END IF
              DB%TCode%I(DB%LenTC)=IKType
            END IF
            N=BufN%I(iIKType,iKonAC)+1
            IF (N>DB%MAXD) THEN
              ErrorCode=eMAXD
              GOTO 9000
            END IF
            BufN%I(iIKType,iKonAC)=N
            BufT%I(N,iIKType,iKonAC)=iBf
            SchT%D(N,iIKType,iKonAC)=rInt
            iBf=iBf+1
            IF (iBf>DB%MAXD) THEN
              ErrorCode=eMAXD
              GOTO 9000
            ENDIF
          END IF ! rInt
        END IF ! I0
      END DO ! CFA
    END IF ! test on AC2
  END DO ! AtA

  DO I=1,DB%LenCC
    KonAC=DB%CCode%I(I)
    DO J=1,DB%LenTC
      N=BufN%I(J,I)
      IF (N>0) THEN
        DB%TCPop%I(J,I)=1
        IF (iDis+N*DB%MAXC>DB%MAXDis) THEN
          ErrorCode=eMAXDis
          GOTO 9000
        END IF
        IF (iPrm+N*(KonAC+MInfo)*DB%MAXP>DB%MAXPrm) THEN
          ErrorCode=eMAXPrm
          GOTO 9000
        END IF
        CALL QuickSortDis(SchT%D(1,J,I),BufT%I(1,J,I),N,-2)
        CALL PutDis(N,iDis,iPrm,I,J,IndexC,KonAC,BufT,DB)
      END IF
    END DO
  END DO

    END DO ! CFC
#ifdef PARALLEL
  END IF
#endif
  END DO ! Atc
  ErrorCode=eAOK
  9000 CONTINUE
  CALL Delete(BufN)
  CALL Delete(BufT)
  CALL Delete(SchT)

END SUBROUTINE DisOrder
 
