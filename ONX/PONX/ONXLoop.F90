SUBROUTINE ONXLoop(BSc,GMc,BSp,GMp,D,K,PrmBuf,DisBuf,CBuf,SBuf,VecBuf,BfnBuf, &
                   iSA,iSrt,SubInd,LenF)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE Thresholding
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)      :: BSc,BSp
  TYPE(CRDS),INTENT(IN)      :: GMc,GMp
#ifdef PARALLEL
  TYPE(DBCSR),INTENT(IN)     :: D
  TYPE(DBCSR),INTENT(INOUT)  :: K
#else
  TYPE(BCSR),INTENT(IN)      :: D
  TYPE(BCSR),INTENT(INOUT)   :: K
#endif
  TYPE(DBL_VECT),INTENT(IN)  :: PrmBuf,DisBuf,CBuf,SBuf
  TYPE(DBL_RNK2),INTENT(IN)  :: VecBuf
  TYPE(INT_VECT),INTENT(IN)  :: BfnBuf
  TYPE(INT_RNK2),INTENT(IN)  :: iSA,iSrt,SubInd
  TYPE(BUFL),INTENT(IN)      :: LenF
!--------------------------------------------------------------------------------
! Temporary ERI evaluation buffers. These buffers do not grow with 
! system size, and can be adjusted to balance performace with memory use. 
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)             :: WR,WZ,CBra,CKet,W1,W2
  TYPE(DBL_RNK2)             :: RAC,RBD
  TYPE(INT_VECT)             :: CDrv,IndA
  TYPE(INT_VECT)             :: SLOC,VLOC
  TYPE(INT_RNK2)             :: Bfn,Indx,BTable
!--------------------------------------------------------------------------------
! Misc. internal variables...
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)             :: DA
  REAL(DOUBLE)               :: ErrorFC
  REAL(DOUBLE)               :: ESw=22.0D0,T,TE
  REAL(DOUBLE)               :: Tsh,ThreshM,ThreshP,TxDcd,Dcd
  INTEGER                    :: i,LngDrv,NPrim2,ri,ci,iPtr,IndexC,IndexD
  INTEGER                    :: LngVRR,LngLoc
  INTEGER                    :: AtC,CFC,StartLC,StopLC,StrideC
  INTEGER                    :: AtD,CFD,StartLD,StopLD,StrideD
  INTEGER                    :: KC,NBFC,MinLC,MaxLC,KType,NKCase
  INTEGER                    :: KD,NBFD,MinLD,MaxLD,LType,NLCase
  INTEGER                    :: LCmin,LDmin,LCmax,LDmax
  REAL(DOUBLE)               :: Cx,Cy,Cz,Dx,Dy,Dz,xNERIs,xNMults
!--------------------------------------------------------------------------------
  INCLUDE "bkonx_F90.inc"
!--------------------------------------------------------------------------------
  NPrim2=MAX(BSc%NPrim,BSp%NPrim)**2
!--------------------------------------------------------------------------------
! Allocate temporary ERI evaluation buffers.
!--------------------------------------------------------------------------------
  CALL New(WR,12*MXBUF1)
  CALL New(WZ,5*MXBUF1)
  CALL New(CBra,4*MXBUF2)
  CALL New(CKet,4*MXBUF2)
  CALL New(W1,MXINTS)
  CALL New(W2,MXINTS)
  CALL New(RAC,(/3,MXBATCH/))
  CALL New(RBD,(/3,MXBATCH/))
  CALL New(CDrv,60000)
  CALL New(IndA,500,0)
  CALL New(Bfn,(/2,MXBATCH/))
  CALL New(Indx,(/2,MXBATCH/))
  CALL New(BTable,(/NPrim2,NPrim2/))
  CALL New(DA,BSp%LMNLen*BSp%LMNLen)

  xNERIs=0.0D0
  xNMults=0.0D0

  DO i=10,22
    T=SQRT(FLOAT(i))
    TE=ErrorFC(T)
    IF (TE.LE.Thresholds%TwoE) THEN
      ESw=FLOAT(i)
      EXIT
    ENDIF
  END DO

  Tsh=Thresholds%TwoE*1.0D-2

  CALL CCDriver(CDrv%I(1),LngDrv)
  CALL VRRLng(LngVRR,LngLoc)
  CALL New(VLOC,LngVRR)
  CALL New(SLOC,LngLoc)
  CALL VRRDriver(VLOC%I(1),SLOC%I(1),LngVRR,LngLoc)
  CALL BuffTbl(NPrim2,BTable%I(1,1))

#ifdef PARALLEL
  DO AtC=Beg%I(MyID),End%I(MyID)
    ri=AtC-Beg%I(MyID)+1
#else
  DO AtC=1,NAtoms
    ri=AtC
#endif
    KC=GMp%AtTyp%I(AtC)
    NBFC=BSp%BfKnd%I(KC)
    Cx=GMp%Carts%D(1,AtC)
    Cy=GMp%Carts%D(2,AtC)
    Cz=GMp%Carts%D(3,AtC)

    DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1
      AtD=D%ColPt%I(ci)
      iPtr=D%BlkPt%I(ci)
      KD=GMp%AtTyp%I(AtD)
      NBFD=BSp%BfKnd%I(KD)
      Dx=GMp%Carts%D(1,AtD)
      Dy=GMp%Carts%D(2,AtD)
      Dz=GMp%Carts%D(3,AtD)
 
      IndexC=0
      DO CFC=1,BSp%NCFnc%I(KC)
        StartLC=BSp%LStrt%I(CFC,KC)
        StopLC=BSp%LStop%I(CFC,KC)
        StrideC=StopLC-StartLC+1
        MinLC=BSp%ASymm%I(1,CFC,KC)
        MaxLC=BSp%ASymm%I(2,CFC,KC)
        KType=MaxLC*(MaxLC+1)/2+MinLC+1

        IndexD=0
        DO CFD=1,BSp%NCFnc%I(KD)
          StartLD=BSp%LStrt%I(CFD,KD)
          StopLD=BSp%LStop%I(CFD,KD)
          StrideD=StopLD-StartLD+1
          MinLD=BSp%ASymm%I(1,CFD,KD)
          MaxLD=BSp%ASymm%I(2,CFD,KD)
          LType=MaxLD*(MaxLD+1)/2+MinLD+1

          CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD, &
                         IndexC+1,IndexD+1,D%MTrix%D(iPtr),DA%D(1))

          Dcd=0.0D0
          DO i=1,StrideC*StrideD
            Dcd=MAX(Dcd,ABS(DA%D(i)))
          END DO

          IF(Dcd.LE.Thresholds%TwoE) EXIT
          TxDcd=Thresholds%TwoE/Dcd
          ThreshM=TxDcd
          ThreshP=TxDcd

#ifdef PARALLEL
          CALL CalcK(NAtoms,BSc%NCtrt,NPrim2,KType,LType,            &
                     VLOC%I(1),SLOC%I(1),LngVRR,LngLoc,              &
                     BTable%I(1,1),BfnBuf%I(1),SubInd%I(1,1),NBasF,  &
                     LenF%PrmL,LenF%DisL,LenF%VecL,LenF%CosL,        &
                     LenF%SinL,LenF%IsaL,LenF%BfnL,LenF%SrtL,        &
                     iSA%I(1,1),iSrt%I(1,1),CDrv%I(1),IndA%I(0),     &
                     Bfn%I(1,1),Indx%I(1,1),DA%D(1),K%MTrix%D(1),    &
                     K%GRwPt%I(1),K%ColPt%I(1),K%BlkPt%I(1),         &
                     PrmBuf%D(1),DisBuf%D(1),CBuf%D(1),SBuf%D(1),    &
                     VecBuf%D(1,1),WR%D(1),WZ%D(1),W1%D(1),W2%D(1),  &
                     RAC%D(1,1),RBD%D(1,1),CBra%D(1),CKet%D(1),      &
                     ESw,ThreshM,ThreshP,xNERIs,xNMults)
#else
          CALL CalcK(NAtoms,BSc%NCtrt,NPrim2,KType,LType,            &
                     VLOC%I(1),SLOC%I(1),LngVRR,LngLoc,              &
                     BTable%I(1,1),BfnBuf%I(1),SubInd%I(1,1),NBasF,  &
                     LenF%PrmL,LenF%DisL,LenF%VecL,LenF%CosL,        &
                     LenF%SinL,LenF%IsaL,LenF%BfnL,LenF%SrtL,        &
                     iSA%I(1,1),iSrt%I(1,1),CDrv%I(1),IndA%I(0),     &
                     Bfn%I(1,1),Indx%I(1,1),DA%D(1),K%MTrix%D(1),    &
                     K%RowPt%I(1),K%ColPt%I(1),K%BlkPt%I(1),         &
                     PrmBuf%D(1),DisBuf%D(1),CBuf%D(1),SBuf%D(1),    &
                     VecBuf%D(1,1),WR%D(1),WZ%D(1),W1%D(1),W2%D(1),  &
                     RAC%D(1,1),RBD%D(1,1),CBra%D(1),CKet%D(1),      &
                     ESw,ThreshM,ThreshP,xNERIs,xNMults)
#endif

          IndexD=IndexD+StrideD
        END DO ! CFD
        IndexC=IndexC+StrideC
      END DO ! CFC

    END DO ! AtD
  END DO ! AtC


  CALL Delete(WR)
  CALL Delete(WZ)
  CALL Delete(CBra)
  CALL Delete(CKet)
  CALL Delete(W1)
  CALL Delete(W2)
  CALL Delete(RAC)
  CALL Delete(RBD)
  CALL Delete(CDrv)
  CALL Delete(IndA)
  CALL Delete(Bfn)
  CALL Delete(Indx)
  CALL Delete(BTable)
  CALL Delete(DA)

  write(*,*) ">>>>>> PONX:xNERIs  = ",xNERIs
  write(*,*) ">>>>>> PONX:xNMults = ",xNMults

END SUBROUTINE ONXLoop
