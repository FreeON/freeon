MODULE DOrder
   USE DerivedTypes
   USE GlobalScalars
   USE PrettyPrint
   USE Thresholding
   USE ONXParameters
   USE ONXMemory
   USE Stats
   IMPLICIT NONE
   REAL(DOUBLE),DIMENSION(3) :: PBC
!  Local buffers allocated in MemInit outside of periodic loops
   TYPE(DBL_RNK3)         :: SchT 
   TYPE(INT_RNK3)         :: BufT
   TYPE(INT_RNK2)         :: BufN
   CONTAINS 
      SUBROUTINE DisOrder(BSc,GMc,BSp,GMp,DB,IB,SB,Drv,NameBuf)  
!--------------------------------------------------------------------------------
! Basis set, coordinates, ect...
!--------------------------------------------------------------------------------
  TYPE(BSET),INTENT(IN)    :: BSc,BSp   ! basis set info
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp   ! geometry info
  TYPE(DBuf)               :: DB        ! ONX distribution buffers
  TYPE(IBuf)               :: IB        ! ONX 2-e eval buffers
  TYPE(DSL)                :: SB        ! ONX distribution pointers
  TYPE(ISpc)               :: IS
  TYPE(IDrv)               :: Drv       ! VRR/contraction drivers
  TYPE(INT_VECT)           :: NameBuf   ! for parallel implementation
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  REAL(DOUBLE)           :: Test,VSAC
  REAL(DOUBLE)           :: Ax,Ay,Az,Cx,Cy,Cz,ACx,ACy,ACz,AC2,x,y,z
  REAL(DOUBLE)           :: Zeta,Za,Zc,Cnt,rInt,XiAB
  INTEGER                :: Lng,i,j,n
  INTEGER                :: KonAC,LDis,NInts
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
  EXTERNAL               INT_VECT_EQ_INT_SCLR

  DB%LenCC=0
  DB%LenTC=0

! The following zeroing of DB%DisPtr%I and DB%TcPop%I 
! is is really expensive for periodic systems.  Unfortunate that 
! ComputeK depends critically on these being carefully zeroed.
!
!  DB%TcPop%I=0  
!  DB%DisPtr%I=0 <- This one really costs
   N=SIZE(DB%TcPop%I,1)*SIZE(DB%TcPop%I,2)
   CALL INT_VECT_EQ_INT_SCLR(N,DB%TcPop%I(1,1),0)
   N=SIZE(DB%DisPtr%I,1)*SIZE(DB%DisPtr%I,2)* &
     SIZE(DB%DisPtr%I,3)*SIZE(DB%DisPtr%I,4)
   CALL INT_VECT_EQ_INT_SCLR(N,DB%DisPtr%I(1,1,1,1),0)

  iDis=1
  iPrm=1
  Test=-10.d0*DLOG(Thresholds%Dist) 
  SB%SLDis%I(1)=5
  IndexC=0
  DO AtC=1,NAtoms
#ifdef PARALLEL_ONX
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
!     Another expensive zeroing...
      N=SIZE(BufN%I,1)*SIZE(BufN%I,2)
      CALL INT_VECT_EQ_INT_SCLR(N,BufN%I(1,1),0)
!      BufN%I=0

      iBf=1 
      IndexA=0  

  DO AtA=1,NAtoms
    KA=GMc%AtTyp%I(AtA) 
    NBFA=BSc%BfKnd%I(KA)
#ifdef PERIODIC
    Ax=GMc%Carts%D(1,AtA)
    Ay=GMc%Carts%D(2,AtA)
    Az=GMc%Carts%D(3,AtA)
    Cx=GMp%Carts%D(1,AtC)+PBC(1)
    Cy=GMp%Carts%D(2,AtC)+PBC(2)
    Cz=GMp%Carts%D(3,AtC)+PBC(3)
#else
    Ax=GMc%Carts%D(1,AtA)
    Ay=GMc%Carts%D(2,AtA)
    Az=GMc%Carts%D(3,AtA)
    Cx=GMp%Carts%D(1,AtC)
    Cy=GMp%Carts%D(2,AtC)
    Cz=GMp%Carts%D(3,AtC)
#endif
    ACx=Ax-Cx
    ACy=Ay-Cy
    ACz=Az-Cz
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
          DB%TBufC%D( 1,iBf)=DBLE(IJ)+Small
          DB%TBufC%D( 2,iBf)=-DBLE(IndexA)-Small
          DB%TBufC%D( 3,iBf)=DBLE(IndexC)+Small
          DB%TBufC%D( 4,iBf)=DBLE(IndexA)+Small
          DB%TBufC%D( 5,iBf)=ACx
          DB%TBufC%D( 6,iBf)=ACy
          DB%TBufC%D( 7,iBf)=ACz
          DB%TBufC%D( 8,iBf)=Ax
          DB%TBufC%D( 9,iBf)=Ay
          DB%TBufC%D(10,iBf)=Az
          x=ACx
          y=ACy
          z=ACz
        ELSE 
          IKType=KType*100+IType
          DB%TBufC%D( 1,iBf)=DBLE(IJ)+Small
          DB%TBufC%D( 2,iBf)=DBLE(IndexA)+Small
          DB%TBufC%D( 3,iBf)=DBLE(IndexC)+Small
          DB%TBufC%D( 4,iBf)=DBLE(IndexA)+Small
          DB%TBufC%D( 5,iBf)=-ACx
          DB%TBufC%D( 6,iBf)=-ACy
          DB%TBufC%D( 7,iBf)=-ACz
          DB%TBufC%D( 8,iBf)=Cx
          DB%TBufC%D( 9,iBf)=Cy
          DB%TBufC%D(10,iBf)=Cz
          x=-ACx
          y=-ACy
          z=-ACz
        ENDIF

        LDis=MaxLA+MaxLC
        CALL GetIntSpace(IKType,IKType,LDis,LDis,IS)

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
            DB%TBufP%D(2,I0,iBf)=(Za*Ax+Zc*Cx)/Zeta
            DB%TBufP%D(3,I0,iBf)=(Za*Ay+Zc*Cy)/Zeta
            DB%TBufP%D(4,I0,iBf)=(Za*Az+Zc*Cz)/Zeta
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
          NInts=IS%L1*IS%L2*IS%L3*IS%L4
          I0=iT(MIN(IType,KType),MAX(IType,KType))
          I1=11*I0-10
          iCP=Drv%CDrv%I(I1)
          iCL=Drv%CDrv%I(iCP)

          IF (2*KonAC*KonAC*IS%NVRR>IB%MAXI.OR.NInts>IB%MAXI) THEN
            ErrorCode=eMAXI
            GOTO 9000
          ENDIF

          ! Beware copy in copy out...
          CALL RGen1C(2*LDis,iBf,KonAC,IB%CB%D(1,1),IB%WR,IB%WZ, &
                      IB%W1%D(1),DB%TBufP,DB%TBufC)
          CALL VRRs(LDis,LDis,Drv)
          CALL VRRl(KonAC*KonAC,IS%NVRR,Drv%nr,Drv%ns,                  &
                    Drv%VLOC%I(Drv%is),                                 &
                    Drv%VLOC%I(Drv%is+Drv%nr),IB,                       &
                    IB%W2%D(1),IB%W1%D(1))
          CALL Contract(1,KonAC,KonAC,IS%NVRR,iCL,Drv%CDrv%I(iCP+1),    &
                        IB%CB%D(1,1),IB%CB%D,IB%W1%D(1),IB%W2%D(1))

          IF(LDis.NE.0) THEN
            CALL HRRKet(IB%W1%D,DB%TBufC%D(1,iBf),1,SB%SLDis%I,IS%NB1,  &
                        IS%NB2,IS%NK1,IKType)
            CALL HRRBra(IB%W1%D,IB%W2%D,x,y,z,1,IS%NB1,IS%L1*IS%L2,     &
                        IS%L3*IS%L4,IKType)
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
        IF (iPrm+N*(KonAC+DB%MInfo)*DB%MAXP>DB%MAXPrm) THEN
          ErrorCode=eMAXPrm
          GOTO 9000
        END IF
        CALL QuickSortDis(SchT%D(1,J,I),BufT%I(1,J,I),N,-2)
        CALL PutDis(N,iDis,iPrm,I,J,IndexC,KonAC,BufT,DB)
      END IF
    END DO
  END DO

    END DO ! CFC
#ifdef PARALLEL_ONX
  END IF
#endif
  END DO ! Atc
  ErrorCode=eAOK
  9000 CONTINUE

END SUBROUTINE DisOrder
END MODULE DOrder 
