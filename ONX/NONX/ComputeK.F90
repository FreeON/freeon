SUBROUTINE ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,MB,Drv)
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
  TYPE(BSET),INTENT(IN) :: BSc,BSp   ! basis set info
  TYPE(CRDS),INTENT(IN) :: GMc,GMp   ! geometry info
  TYPE(DBuf)            :: DB        ! ONX distribution buffers
  TYPE(IBuf)            :: IB        ! ONX 2-e eval buffers
  TYPE(DML)             :: MB        ! ONX distribution pointers
  TYPE(IDrv)            :: Drv       ! VRR/contraction drivers
!-------------------------------------------------------------------
! Misc. internal variables
!-------------------------------------------------------------------
  INTEGER               :: iTBra,iCBra,TBra,LBra,ITypeC,ITypeA,CBra
  INTEGER               :: iTKet,iCKet,TKet,LKet,ITypeD,ITypeB,CKet
  INTEGER               :: LenBra,iDBra,iPBra,iBra
  INTEGER               :: LenKet,iDKet,iPKet,iKet
  INTEGER               :: LTot,ShellC,ShellD
  INTEGER               :: AtC,KC,CFC,StartLC,StopLC,StrideC,IndexC,NBFC
  INTEGER               :: AtD,KD,CFD,StartLD,StopLD,StrideD,IndexD,NBFD
  INTEGER               :: ri,ci,iPtr
  TYPE(INT_VECT)        :: IndA
  REAL(DOUBLE)          :: Dcd,TxDcd,SchB,SchK
  TYPE(DBL_VECT)        :: DA
!-------------------------------------------------------------------
! Function calls
!-------------------------------------------------------------------
  INTEGER               :: LTotal,MaxBatchSize

  CALL GetExpTable(IB)
  CALL New(DA,BSp%LMNLen*BSp%LMNLen)

  DO iTBra=1,DB%LenTC       ! Loop over angular symmetry types on the Bra
    TBra=DB%TCode%I(iTBra)
    ITypeC=MOD(TBra,100)
    ITypeA=(TBra-ITypeC)/100
    LBra=LTotal(ITypeA)+LTotal(ITypeC)

  DO iTKet=1,DB%LenTC       ! Loop over angular symmetry types on the Ket
    TKet=DB%TCode%I(iTKet)
    ITypeD=MOD(TKet,100)
    ITypeB=(TKet-ITypeD)/100
    LKet=LTotal(ITypeB)+LTotal(ITypeD)
      
    LTot=LBra+LKet
    CALL GetGammaTable(LTot,IB)

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
          ShellC = ShellC+1
          StartLC=BSp%LStrt%I(CFC,KC)
          StopLC=BSp%LStop%I(CFC,KC)
          StrideC=StopLC-StartLC+1
          LenBra = DB%DisPtr%I(1,ShellC,iTBra,iCBra)
          iDBra  = DB%DisPtr%I(2,ShellC,iTBra,iCBra)
          iPBra  = DB%DisPtr%I(3,ShellC,iTBra,iCBra)
          IF (LenBra>0) THEN
          CALL New(IndA,LenBra,0)                  ! Should probably move this out into
                                                   ! DisOrder so that this is only done 
                                                   ! once... 
      ShellD=0
      DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1        ! Loop over atom D (Dcd)
        AtD=D%ColPt%I(ci)
        iPtr=D%BlkPt%I(ci)
        KD=GMp%AtTyp%I(AtD)
        NBFD=BSp%BfKnd%I(KD)

          write(*,*) "AtC=",AtC," NBFC=",NBFC,LBra,CBra
          write(*,*) "AtD=",AtD," NBFD=",NBFD,LKet,CKet,iPtr

        IndexD=0
        DO CFD=1,BSp%NCFnc%I(KD)
          ShellD = ShellD+1
          StartLD=BSp%LStrt%I(CFD,KD)
          StopLD=BSp%LStop%I(CFD,KD)
          StrideD=StopLD-StartLD+1
          LenKet = DB%DisPtr%I(1,ShellD,iTKet,iCKet)
          iDKet  = DB%DisPtr%I(2,ShellD,iTKet,iCKet)
          iPKet  = DB%DisPtr%I(3,ShellD,iTKet,iCKet)
          IF (LenKet>0) THEN

            write(*,*) "Len = ",LenBra,LenKet
            CALL GetSubBlk(NBFC,NBFD,StrideC,StrideD,IndexC+1,  &
                           IndexD+1,D%MTrix%D(iPtr),DA%D(1))
            Dcd=GetAbsMax(NBFC*NBFD,DA)
            TxDcd=Thresholds%TwoE/Dcd

            CALL SkipOut(LenBra,LenKet,DB%SchT%D(1,iTBra,iCBra), &
                         DB%SchT%D(1,iTKet,iCKet),IndA%I(0),TxDcd)

            DO IBra=0,LenBra
               write(*,*) "IndA=",IndA%I(IBra),IBra
!              SchB=DB%SchT%D(IBra,iTBra,iCBra)
!              write(*,*) "SchB = ",SchB
            END DO
!
!            DO IKet=1,LenKet
!              SchK=DB%SchT%D(IKet,iTKet,iCKet)
!              write(*,*) "SchK = ",SchK
!            END DO


          END IF ! LenKet
          IndexD=IndexD+StrideD
        END DO ! CFD
      END DO ! ci
          END IF ! LenBra
          IndexC=IndexC+StrideC
          CALL Delete(IndA)
        END DO ! CFC
      END DO ! AtC

    END IF ! TCPop on Ket
  END DO ! iCKet
    END IF ! TCPop on Bra
  END DO ! iCBra

  END DO ! iTKet
  END DO ! iTBra
  

  CALL Delete(DA)
  CALL Halt('enough for now')

END SUBROUTINE ComputeK
