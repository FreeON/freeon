SUBROUTINE DisOrder(BSc,GMc,BSp,GMp,DB,IB,NameBuf)  
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Basis set, coordinates, ect.
!--------------------------------------------------------------------------------
  TYPE(BSET),INTENT(IN)    :: BSc,BSp
  TYPE(CRDS),INTENT(IN)    :: GMc,GMp
  TYPE(DBuf)               :: DB
  TYPE(IBuf)               :: IB
  TYPE(INT_VECT)           :: NameBuf
!--------------------------------------------------------------------------------
! Temporary buffers to hold the distribution data
!--------------------------------------------------------------------------------
  TYPE(DBL_RNK2)           :: TBufC
  TYPE(DBL_RNK3)           :: TBufP
  TYPE(INT_RNK2)           :: BufN
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  TYPE(INT_VECT)           :: MLDis
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  REAL(DOUBLE)           :: Test,VSAC
  REAL(DOUBLE)           :: ACx,ACy,ACz,AC2,x,y,z
  REAL(DOUBLE)           :: Zeta,Za,Zc,Cnt,rInt,XiAB
  INTEGER                :: LngVRR,LngLoc,LngDrv
  INTEGER                :: id,is,nr,ns,i
  INTEGER                :: KonAC,LDis,NLOCD2,NLOCD3,NInts,NVRR
  INTEGER                :: iBf,I0,I1,NPrim,iCP,iCL,IKType
  INTEGER                :: AtA,CFA,PFA,StartLA,StopLA,StrideA
  INTEGER                :: AtC,CFC,PFC,StartLC,StopLC,StrideC
  INTEGER                :: KA,NBFA,MinLA,MaxLA,IType,NICase,IndexA
  INTEGER                :: KC,NBFC,MinLC,MaxLC,KType,NKCase,IndexC
  INTEGER                :: CType,II,JJ,IJ,LngTmp
  INTEGER                :: LenCC,LenTC,iKonAC,iIKType
  LOGICAL                :: Found
!--------------------------------------------------------------------------------
! Function calls
!--------------------------------------------------------------------------------
  INTEGER                :: NFinal,iT

  NPrim=MAX(BSc%NPrim,BSp%NPrim)
  CALL New(TBufC,(/DB%MAXC,DB%MAXD/))
  CALL New(TBufP,(/DB%MAXP,NPrim*NPrim+MInfo,DB%MAXD/))
  CALL New(MLDis,1)
  CALL New(BufN,(/DB%MAXT,NPrim*NPrim/))


  Test=-10.d0*DLOG(Thresholds%Dist) 
  MLDis%I(1)=5
  BufN%I=0
  IndexC=0
  LenCC=0
  LenTC=0
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
        TBufC%D( 1,iBf)=DFLOAT(IJ)+Small
        TBufC%D( 2,iBf)=DFLOAT(IndexA)+Small
        IF(IType.GE.KType) THEN
          IKType=IType*100+KType
          TBufC%D( 3,iBf)=DFLOAT(IndexA)+Small
          TBufC%D( 4,iBf)=DFLOAT(IndexC)+Small
          TBufC%D( 5,iBf)=ACx
          TBufC%D( 6,iBf)=ACy
          TBufC%D( 7,iBf)=ACz
          TBufC%D( 8,iBf)=GMc%Carts%D(1,AtA)
          TBufC%D( 9,iBf)=GMc%Carts%D(2,AtA)
          TBufC%D(10,iBf)=GMc%Carts%D(3,AtA)
          x=ACx
          y=ACy
          z=ACz
        ELSE 
          IKType=KType*100+IType
          TBufC%D( 3,iBf)=DFLOAT(IndexC)+Small
          TBufC%D( 4,iBf)=DFLOAT(IndexA)+Small
          TBufC%D( 5,iBf)=-ACx
          TBufC%D( 6,iBf)=-ACy
          TBufC%D( 7,iBf)=-ACz
          TBufC%D( 8,iBf)=GMp%Carts%D(1,AtC)
          TBufC%D( 9,iBf)=GMp%Carts%D(2,AtC)
          TBufC%D(10,iBf)=GMp%Carts%D(3,AtC)
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
            TBufP%D(1,I0,iBf)=Zeta
            TBufP%D(2,I0,iBf)=(Za*GMc%Carts%D(1,AtA)+Zc*GMp%Carts%D(1,AtC))/Zeta
            TBufP%D(3,I0,iBf)=(Za*GMc%Carts%D(2,AtA)+Zc*GMp%Carts%D(2,AtC))/Zeta
            TBufP%D(4,I0,iBf)=(Za*GMc%Carts%D(3,AtA)+Zc*GMp%Carts%D(3,AtC))/Zeta
            TBufP%D(5,I0,iBf)=VSAC
            IF (CType.EQ.11) THEN
              TBufP%D(6,I0,iBf)=1.0D0
              TBufP%D(7,I0,iBf)=1.0D0
              TBufP%D(8,I0,iBf)=1.0D0
            ELSEIF (CType.EQ.12.OR.CType.EQ.21) THEN
              TBufP%D(6,I0,iBf)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              TBufP%D(7,I0,iBf)=TBufP%D(6,I0,iBf)
              TBufP%D(8,I0,iBf)=TBufP%D(6,I0,iBf)
            ELSEIF (CType.EQ.22) THEN
              TBufP%D(6,I0,iBf)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              TBufP%D(7,I0,iBf)=BSc%CCoef%D(StartLA+1,PFA,CFA,KA) * &
                               BSp%CCoef%D(StartLC,PFC,CFC,KC) / Cnt
              TBufP%D(8,I0,iBf)=BSc%CCoef%D(StartLA,PFA,CFA,KA) * &
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
          iCP=IB%CDrv%I(I1)
          iCL=IB%CDrv%I(iCP)

          IF (KonAC*KonAC*NVRR>IB%MAXI.OR.NInts>IB%MAXI) THEN
            ErrorCode=eMAXI
            GO TO 1000
          ENDIF

          CALL RGen1C(2*LDis,iBf,KonAC,IB%WK,IB%CD%D(1,1),IB%WR,IB%WZ, &
                      IB%W1%D(1),TBufP,TBufC)
          CALL VRRs(LDis,LDis,id,is,nr,ns,IB%SLOC)
          CALL VRRl(KonAC*KonAC,NVRR,nr,ns,IB%VLOC%I(is),IB%VLOC%I(is+nr), &
                    IB%W2%D(1),IB%W1%D(1),IB%WR%D(1,1,1),IB%WZ%D(1,1,1))

          CALL Contract(1,KonAC,KonAC,NVRR,iCL,IB%CDrv%I(iCP+1), &
                        IB%CD%D(1,1),IB%CD%D(1,1),IB%W1%D(1),IB%W2%D(1))

          IF(LDis.NE.0) THEN
            CALL HRRKetOld(IB%W1%D(1),TBufC%D(1,iBf),1,MLDis%I(1),NLOCD2,DB%MAXC,IKType)
            CALL HRRBraOld(IB%W1%D(1),IB%W2%D(1),x,y,z,1,NLOCD2,NLOCD3,NLOCD3,IKType)
            rInt = AbsMax(NInts,IB%W2)
          ELSE
            rInt = AbsMax(NInts,IB%W1)
          ENDIF

          IF(rInt>Thresholds%Dist) THEN 
            iBf=iBf+1
            IF (iBf>DB%MAXD) THEN
              ErrorCode=eMAXD
              GO TO 1000
            ENDIF
            Found=.FALSE.
            DO I=1,LenCC       ! look for the contraction type
              IF (KonAC==DB%CCode%I(I)) THEN
                iKonAC=I
                Found=.TRUE.
                EXIT
              END IF
            END DO
            IF (.NOT.Found) THEN  ! if not found make a new CC type
              LenCC=LenCC+1
              iKonAC=LenCC
              IF (LenCC>DB%MAXC) THEN
                ErrorCode=eMAXC
                GO TO 1000
              END IF
              DB%CCode%I(LenCC)=KonAC
            END IF
            Found=.FALSE.
            DO I=1,LenTC       ! look for the angular symmetry type
              IF (IKType==DB%TCode%I(I)) THEN
                iIKType=I
                Found=.TRUE.
                EXIT
              ENDIF
            ENDDO
            IF (.NOT.Found) THEN  ! if not found make a new L type
              LenTC=LenTC+1
              iIKType=LenTC
              IF (LenTC.GT.DB%MAXT) THEN
                ErrorCode=eMAXT
                GO TO 1000
              END IF
              DB%TCode%I(LenTC)=IKType
            END IF
            BufN%I(iIKType,iKonAC)=BufN%I(iIKType,iKonAC)+1
          END IF ! rInt

        END IF ! I0
      ENDDO ! CFA
    END IF ! test on AC2
  ENDDO ! AtA

! sort the distributions

! call adddis


    END DO ! CFC
#ifdef PARALLEL
  END IF
#endif
  END DO ! Atc

  ErrorCode=0
1000 CONTINUE
  CALL Delete(MLDis)
  CALL Delete(TBufP)
  CALL Delete(TBufC)
  CALL Delete(BufN)

END SUBROUTINE DisOrder
 
