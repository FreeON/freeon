MODULE GONX2ComputDK
!H=================================================================================
!H MODULE GONX2ComputDK
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB 
!H
!H  PRIVATE:
!H  o SUB 
!H
!H  OPTIONS:
!H  DEBUGING: Use -DGONX2_DBUG to print some stuff.
!H  INFO    : Use -DGONX2_INFO to print some stuff.
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
#ifdef GONX2_DBUG
#define GONX2_INFO
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONX2DataType
  USE InvExp
  USE ONXParameters
  !
#ifdef PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
  !
  USE LinAlg
  USE ONX2List
  !
  IMPLICIT NONE
  !PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
!  PUBLIC  :: ComputK
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetAtomPairG
!  PRIVATE :: 
  !
CONTAINS
  !
  !
  SUBROUTINE ComputDK(D,GradX,ListC,ListD,GM,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE ComputDK(D,GradX,ListC,ListD,GM,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    USE ONXGet, ONLY: GetAdrB
    !
    IMPLICIT NONE
    !
    !
    ! GradX = Dcd*(ac(R)|bd(R'))'*Dab
    !
    !-------------------------------------------------------------------
    TYPE(BCSR)                 , INTENT(INout) :: D
    TYPE(DBL_RNK2)             , INTENT(INOUT) :: GradX
    TYPE(CList2) , DIMENSION(:), POINTER       :: ListC,ListD
    TYPE(CRDS)                 , INTENT(IN   ) :: GM
    TYPE(BSET)                 , INTENT(IN   ) :: BS
    TYPE(CellSet)              , INTENT(IN   ) :: CS_OUT
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER      :: AtAListTmp,AtAList,AtBListTmp,AtBList
    TYPE(AtomInfo)             :: ACAtmInfo,BDAtmInfo
    INTEGER                    :: AtA,AtB,AtC,AtD,KA,KB,KC,KD,CFA,CFB,CFC,CFD
    INTEGER                    :: ci,iPtrD,iPtrD2,iPtrK,NBFC,NBFD,NBFA,NBFB
    INTEGER                    :: CFAC,CFBD
    INTEGER                    :: Off,Ind
    REAL(DOUBLE)               :: Dcd
    !-------------------------------------------------------------------
    REAL(DOUBLE), DIMENSION(12*MaxFuncPerAtmBlk**4) :: C
    REAL(DOUBLE), DIMENSION(   MaxFuncPerAtmBlk**2) :: Work
    TYPE(AtomPrG), DIMENSION(  MaxShelPerAtmBlk**2*CS_OUT%NCells) :: ACAtmPair,BDAtmPair
    !-------------------------------------------------------------------
    INTEGER      :: LocNInt
    REAL(DOUBLE) :: NInts
    !-------------------------------------------------------------------
    TYPE(ONX2OffSt)               :: OffSet
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL :: DGetAbsMax
    REAL(DOUBLE), EXTERNAL :: DDOT
    !-------------------------------------------------------------------
    INTEGER                :: ACR,BDR,IXYZ,NIntBlk,Indx,SumInt
    REAL(DOUBLE)           :: TmpGradA,TmpGradC,TmpGradB,TmpGradD
    !REAL(DOUBLE), PARAMETER :: ThresholdTwoE=-1.0D-15
    REAL(DOUBLE), PARAMETER :: ThresholdTwoE=-1.0D-12
    !-------------------------------------------------------------------
    real(double) :: MaxInt
    integer :: MaxCont,IntType
    integer :: iint ,isize,i,iii

    integer, dimension(NAtoms) :: BCol
    TYPE(INT_VECT) :: BColIdx

    !
    NULLIFY(AtAListTmp,AtAList,AtBListTmp,AtBList)
    !
    !Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,BS%BfKnd%I(GM%AtTyp%I(i)))
       IF(12*(BS%BfKnd%I(GM%AtTyp%I(i)))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',12*BS%BfKnd%I(GM%AtTyp%I(i))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo
    !write(*,*) 'size C=',12*isize**4
    if(CS_OUT%NCells.GT.SIZE(ACAtmPair)) then
       write(*,*) 'size(ACAtmPair)',size(ACAtmPair),'.LT.',CS_OUT%NCells
       STOP 'Incrase the size of ACAtmPair and BDAtmPair'
    endif
    !write(*,*) 'CS_OUT%NCells=',CS_OUT%NCells
    !Simple check Simple check Simple check Simple check
    !
    !CALL multiply(D,2.0d0)
    !CALL Print_BCSR(D,'Density matrix',Unit_O=6)
    !CALL multiply(D,0.5d0)


    !write(*,*) 'atd(14)', D%ColPt%I(D%RowPt%I(14):D%RowPt%I(14+1)-1)
    !stop

    CALL New(BColIdx,Natoms)

    SumInt=0.0d0
    MaxCont=0
    iint=0
    !
    LocNInt=0
    NInts=0.0d0
    !
    DO AtC=1,NAtoms ! Run over AtC.
       KC=GM%AtTyp%I(AtC)
       NBFC=BS%BfKnd%I(KC)
       ACAtmInfo%Atm2X=GM%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GM%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GM%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       ! Get AtA List.
       AtAListTmp=>ListC(AtC)%GoList
       !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       DO ci=D%RowPt%I(AtC),D%RowPt%I(AtC+1)-1 ! Run over AtD
          AtD=D%ColPt%I(ci)
          iPtrD=D%BlkPt%I(ci)

!vw       DO AtD=1,NAtoms
!vw          CALL GetAdrB(AtC,AtD,Ind,D,1)
!vw          write(*,*) 'AtC=',AtC,' AtD=',AtD,' Ind=',Ind
!vw          IF(Ind.GT.0) THEN
!vw          iPtrD=D%BlkPt%I(Ind)
          !write(*,*) 'D',D%MTrix%D(iPtrD:iPtrD+NBFC*NBFD-1)
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

          KD=GM%AtTyp%I(AtD)
          NBFD=BS%BfKnd%I(KD)
          BDAtmInfo%Atm2X=GM%Carts%D(1,AtD)
          BDAtmInfo%Atm2Y=GM%Carts%D(2,AtD)
          BDAtmInfo%Atm2Z=GM%Carts%D(3,AtD)
          BDAtmInfo%K2=KD

          !write(*,*) 'Info',KC,KD,NBFC,NBFD
          !write(*,*) GM%Carts%D(:,AtC),GM%Carts%D(:,AtD)
          !
          ! Get max of the block density matrix.
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
          !
#ifdef GONX2_DBUG
          WRITE(*,*) 'Max(Dcd)=',Dcd
#endif
          !
          ! Get AtB List.
          AtBListTmp=>ListD(AtD)%GoList
          !
          AtAList=>AtAListTmp
          !
          DO ! Run over AtA
             AtA=AtAList%Atom
!old           DO AtA=1,NAtoms
             KA=GM%AtTyp%I(AtA)
             NBFA=BS%BfKnd%I(KA)
             !
             ACAtmInfo%NCell=GetNonNeglCell(AtAList,AtBListTmp%SqrtInt(1),ThresholdTwoE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SKIP OUT SKIP OUT SKIP OUT SKIP OUT SKIP OUT !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             !
             !
             ACAtmInfo%Atm1X=GM%Carts%D(1,AtA)
             ACAtmInfo%Atm1Y=GM%Carts%D(2,AtA)
             ACAtmInfo%Atm1Z=GM%Carts%D(3,AtA)
             ACAtmInfo%K1=KA
             !
             ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
             ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
             ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
             !
             ! Get atom pair for BD.
             CALL GetAtomPairG(ACAtmInfo,AtAList,ACAtmPair,BS,CS_OUT)
             !
             AtBList=>AtBListTmp
             !
             !BCol=0
             !do iii=D%RowPt%I(AtA),D%RowPt%I(AtA+1)-1
             !   BCol(D%ColPt%I(iii))=iii
             !enddo
             
             CALL GetColIdx(AtA,D,BColIdx)
             
             DO ! Run over AtB
                AtB=AtBList%Atom 
!old            DO AtB=1,NAtoms
                !
!vw             CALL GetAdrB(AtA,AtB,Ind,D,1)
                !Ind=BCol(AtB)
                Ind=BColIdx%I(AtB)
                IF(Ind.GT.0) THEN ! Skip out if no density matrix element.
                   iPtrD2 = D%BlkPt%I(Ind)
                   !
                   BDAtmInfo%NCell=GetNonNeglCell(AtBList,AtAList%SqrtInt(1),ThresholdTwoE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SKIP OUT SKIP OUT SKIP OUT SKIP OUT SKIP OUT !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   !
                   KB=GM%AtTyp%I(AtB)
                   NBFB=BS%BfKnd%I(KB)
                   !
                   BDAtmInfo%Atm1X=GM%Carts%D(1,AtB)
                   BDAtmInfo%Atm1Y=GM%Carts%D(2,AtB)
                   BDAtmInfo%Atm1Z=GM%Carts%D(3,AtB)
                   BDAtmInfo%K1=KB
                   !
                   BDAtmInfo%Atm12X=BDAtmInfo%Atm1X-BDAtmInfo%Atm2X
                   BDAtmInfo%Atm12Y=BDAtmInfo%Atm1Y-BDAtmInfo%Atm2Y
                   BDAtmInfo%Atm12Z=BDAtmInfo%Atm1Z-BDAtmInfo%Atm2Z
                   !
                   ! Get atom pair for BD.
                   CALL GetAtomPairG(BDAtmInfo,AtBList,BDAtmPair,BS,CS_OUT)
                   !
                   NIntBlk=NBFA*NBFB*NBFC*NBFD
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(12*NIntBlk,C(1),0.0d0)
                   !C=BIG_DBL
                   !
                   CFAC=0
!!!!!!!!!!!!!!!!!!!!
                   DO ACR=1,ACAtmInfo%NCell
!!!!!!!!!!!!!!!!!!!!
                   OffSet%A=1
                   DO CFA=1,BS%NCFnc%I(KA) ! Run over blkfunc on A
                      OffSet%C=1
                      DO CFC=1,BS%NCFnc%I(KC) ! Run over blkfunc on C
                         CFAC=CFAC+1
                         MaxCont=MAX(MaxCont,ACAtmPair(CFAC)%SP%L)
                         CFBD=0
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         DO BDR=1,BDAtmInfo%NCell
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         OffSet%B=1
                         DO CFB=1,BS%NCFnc%I(KB) ! Run over blkfunc on B
                            OffSet%D=1
                            DO CFD=1,BS%NCFnc%I(KD) ! Run over blkfunc on D
                               CFBD=CFBD+1
                               !
                               ! Compute integral type.
                               IntType=ACAtmPair(CFAC)%SP%IntType*10000+BDAtmPair(CFBD)%SP%IntType
                               !
!if(IntType.ne.1010101) goto 20
!if(IntType.ne.3010101.and.IntType.ne.1030101.and.IntType.ne.1010301.and.IntType.ne.1010103) goto 20
!if(IntType.ne.3030101.and.IntType.ne.3010301.and.IntType.ne.3010103.and. &
!&  IntType.ne.1030301.and.IntType.ne.1030103.and.IntType.ne.1010303) goto 20
!if(IntType.ne.3030101.and.IntType.ne.1010303) goto 20
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!if(IntType.ne.3030101) goto 20
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!if(IntType.ne.1030303.and.IntType.ne.3010303.and.IntType.ne.3030103.and.IntType.ne.3030301) goto 20
!if(IntType.ne.3030303) goto 20
!if(IntType.ne.3010101) goto 20
!if(IntType.ne.1010103) goto 20
!if(IntType.ne.1030101) goto 20
                               !write(*,*) 'IntType',IntType,NIntBlk
                               ! The integral interface.
                               INCLUDE 'DERIInterface.Inc'
                               !INCLUDE 'dERIInclude.Inc'
                               !
                               NInts=NInts+DBLE(LocNInt)
                               OffSet%D=OffSet%D+BS%LStop%I(CFD,KD)-BS%LStrt%I(CFD,KD)+1
                            ENDDO ! End blkfunc on D
                            OffSet%B=OffSet%B+BS%LStop%I(CFB,KB)-BS%LStrt%I(CFB,KB)+1
                         ENDDO ! End blkfunc on B
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!
                         OffSet%C=OffSet%C+BS%LStop%I(CFC,KC)-BS%LStrt%I(CFC,KC)+1
                      ENDDO ! End blkfunc on C
                      OffSet%A=OffSet%A+BS%LStop%I(CFA,KA)-BS%LStrt%I(CFA,KA)+1
                   ENDDO ! End blkfunc on A
!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$                   if(ata==2.and.atb==2.and.atc==2.and.atd==2) then
!!$                      write(*,*) NIntBlk,NBFA*NBFB*NBFC*NBFD
!!$                      do iii=1,NBFA*NBFB*NBFC*NBFD
!!$                         if(abs(C(iii)).gt.1d-10)write(*,'(2X,I4,2X,E25.16)') iii,C(iii)
!!$                         maxint=MAX(maxint,abs(c(iii)))
!!$                      enddo
!!$                      write(*,*) 'trace',sum(abs(C(1:NIntBlk))),' maxint',sqrt(maxint)
!!$                      write(*,*) ''
!!$                   endif
!!$                   if(ata==atb.and.atc==atd) then
!!$                      maxint=0.0d0
!!$                      do iii=1,NBFA*NBFB*NBFC*NBFD*12
!!$                         !if(abs(C(iii)).gt.1d-10)write(*,'(2X,I4,2X,E25.16)') iii,C(iii)
!!$                         maxint=MAX(maxint,abs(c(iii)))
!!$                      enddo
!!$                      write(*,'(A,4I3,2X,E25.16)') ' maxint',AtA,AtC,AtB,AtD,sqrt(maxint)
!!$                   endif
!!$                do iii=3*NIntBlk+1,3*NIntBlk+NIntBlk
!!$                   if(abs(C(iii)).gt.1d-10)write(*,'(2X,I4,2X,E25.16)') iii,C(iii)
!!$                enddo
!!$                write(*,*) 'trace',sum(abs(C(3*NIntBlk+1:3*NIntBlk+NIntBlk)))
!!$                write(*,*) ''
!!$                do iii=6*NIntBlk+1,6*NIntBlk+NIntBlk
!!$                   if(abs(C(iii)).gt.1d-10)write(*,'(2X,I4,2X,E25.16)') iii,C(iii)
!!$                enddo
!!$                write(*,*) 'trace',sum(abs(C(6*NIntBlk+1:6*NIntBlk+NIntBlk)))
!!$                write(*,*) ''
!!$                do iii=9*NIntBlk+1,9*NIntBlk+NIntBlk
!!$                   if(abs(C(iii)).gt.1d-10)write(*,'(2X,I4,2X,E25.16)') iii,C(iii)
!!$                enddo
!!$                write(*,*) 'trace',sum(abs(C(9*NIntBlk+1:9*NIntBlk+NIntBlk)))
!!$                write(*,*) ''
!!$                do iii=9*NIntBlk+1,9*NIntBlk+NIntBlk
!!$                   if(abs(C(iii)).gt.1d-10)write(*,'(2X,I4,2X,E25.16)') iii,C(iii)
!!$                enddo

                   !
                   ! Compute the exchange-forces.

                   !if(.not.(ata==14.or.atb==14.or.atc==14.or.atd==14)) goto 20

                   DO IXYZ=1,3
                      !old work=BIG_DBL
                      !
                      ! AtA
                      !if(ixyz.eq.1)CALL PrintMatrix(D%MTrix%D(iPtrD),NBFC,NBFD,2,TEXT_O='Int Dx matrix')
                      Indx=(IXYZ-1)*NIntBlk+1
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradA=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
                      GradX%D(IXYZ,AtA)=GradX%D(IXYZ,AtA)+TmpGradA
                      !write(*,'(A,2X,E24.16)') 'Grad=',TmpGradA

                      !old work=BIG_DBL
                      !
                      ! AtC
                      Indx=(3+IXYZ-1)*NIntBlk+1
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradC=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
                      GradX%D(IXYZ,AtC)=GradX%D(IXYZ,AtC)+TmpGradC

                      !old work=BIG_DBL
                      !
                      ! AtB
                      Indx=(6+IXYZ-1)*NIntBlk+1
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,1.0d0,C(Indx), &
                           &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0,    &
                           &     Work(1),1)
                      TmpGradB=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
                      GradX%D(IXYZ,AtB)=GradX%D(IXYZ,AtB)+TmpGradB
                      !
                      ! AtD
                      !Indx=(9+IXYZ-1)*NIntBlk+1
                      !CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(Indx), &
                      !     &     NBFA*NBFB,D%MTrix%D(iPtrD),1,0.0d0, &
                      !     &     Work(1),1)
                      !TmpGradD=-DDOT(NBFA*NBFB,D%MTrix%D(iPtrD2),1,Work(1),1)
                      !GradX%D(IXYZ,AtD)=GradX%D(IXYZ,AtD)+TmpGradD
                      GradX%D(IXYZ,AtD)=GradX%D(IXYZ,AtD)-(TmpGradA+TmpGradC+TmpGradB)
         !if(ata==1.and.atb==1.and.atc==1.and.atd==1)write(*,*)'TmpGrad',TmpGradA,TmpGradB,TmpGradC
         !if(ata==6.or.atb==6.or.atc==6.or.atd==6)write(*,'(A,2X,E24.16)')'Grad',GradX%D(1,6)
!                 if(ata==6.or.atb==6.or.atc==6.or.atd==6)then
!                    if(ixyz==1)write(*,'(A,2X,E24.16,4I4)')'Grad',Gradx%D(1,6),AtA,AtB,AtC,AtD
!                 endif

                   ENDDO
                   !
20 continue
                   
                   !
                ENDIF ! Skip out if no density matrix element.
                !
                IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT
                AtBList=>AtBList%AtmNext
             ENDDO ! End AtB
             !
             IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT
             AtAList=>AtAList%AtmNext
          ENDDO ! End AtA
          !
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!vw          ENDIF
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       ENDDO ! End AtD
       !
    ENDDO ! End AtC
    !
    CALL Delete(BColIdx)
    !
  END SUBROUTINE ComputDK
  !
  !
  SUBROUTINE GetAtomPairG(AtmInfo,List,AtmPair,BS,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPairG(AtmInfo,List,AtmPair,BS,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    USE Thresholding
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)              , INTENT(IN   ) :: AtmInfo
    TYPE(ANode2  ), POINTER                     :: List
    TYPE(AtomPrG ), DIMENSION(:), INTENT(INOUT) :: AtmPair
    TYPE(BSET)                  , INTENT(IN   ) :: BS
    TYPE(CellSet)               , INTENT(IN   ) :: CS_OUT
    !-------------------------------------------------------------------
    INTEGER      :: CF12,CF1,CF2,I1,I2,II,JJ,IJ,iCell,Cell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StopL1,StartL1,StopL2,StartL2
    INTEGER      :: ISwitch1,ISwitch2
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ
    LOGICAL      :: Switch
    !-------------------------------------------------------------------
    !
    ! dbuging
#ifdef GONX2_DBUG
    AtmPair(:)%SP%IntType=BIG_INT
#endif
    !
    CF12=0
    DO iCell=1,AtmInfo%NCell
       Cell=List%CellIdx(iCell)
       RX=CS_OUT%CellCarts%D(1,Cell)
       RY=CS_OUT%CellCarts%D(2,Cell)
       RZ=CS_OUT%CellCarts%D(3,Cell)
       !
       ! AtmInfo must be related to the atoms in the working cell ONLY. 
       ! Then we add the PBC's to have the right interatomic distance.
       R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
       !

       !write(*,*) 'R12',R12


       DO CF1=1,BS%NCFnc%I(AtmInfo%K1)
          MinL1=BS%ASymm%I(1,CF1,AtmInfo%K1)
          MaxL1=BS%ASymm%I(2,CF1,AtmInfo%K1)
          Type1=MaxL1*(MaxL1+1)/2+MinL1+1
          StartL1=BS%LStrt%I(CF1,AtmInfo%K1)
          StopL1=BS%LStop%I(CF1,AtmInfo%K1)

          if(Type1==2) stop 'SP shell not yet supported.'

          DO CF2=1,BS%NCFnc%I(AtmInfo%K2)
             CF12=CF12+1
             !
             MinL2=BS%ASymm%I(1,CF2,AtmInfo%K2)
             MaxL2=BS%ASymm%I(2,CF2,AtmInfo%K2)
             Type2=MaxL2*(MaxL2+1)/2+MinL2+1
             StartL2=BS%LStrt%I(CF2,AtmInfo%K2)
             StopL2=BS%LStop%I(CF2,AtmInfo%K2)
             !
             AtmPair(CF12)%SP%IntType=Type1*100+Type2
             !
             ! .NOT.(Do we need to switch the centers?)
             Switch=Type1.GE.Type2
             !
             IF(Switch) THEN
                ISwitch1=6
                ISwitch2=7
             ELSE
                ISwitch1=7
                ISwitch2=6
             ENDIF
             !
             II=0
             !
             ! We assume the primitives are ordered (exponants in decressing order).
             DO I1=BS%NPFnc%I(CF1,AtmInfo%K1),1,-1
                Z1=BS%Expnt%D(I1,CF1,AtmInfo%K1)
                JJ=0
                !
                ! We assume the primitives are ordered (exponants in decressing order).
                DO I2=BS%NPFnc%I(CF2,AtmInfo%K2),1,-1
                   Z2=BS%Expnt%D(I2,CF2,AtmInfo%K2)
                   Expt=Z1+Z2
                   InvExpt=1.0d0/Expt
                   XiR12=Z2*Z1*InvExpt*R12
                   IF(XiR12<PrimPairDistanceThreshold) THEN
                      JJ=JJ+1
                      IJ=JJ+II
                      AtmPair(CF12)%SP%Cst(1,IJ)=Expt
                      !
                      ! AtmInfo must be related to the atoms in the working cell ONLY. 
                      ! Then we add the PBC's to have the right atomic position.
                      AtmPair(CF12)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                      AtmPair(CF12)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                      AtmPair(CF12)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                      AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt* &
                      !AtmPair(CF12)%SP%Cst(5,IJ)=5.914967172796D0*EXPInv(XiR12)*InvExpt* &
                           &                     BS%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)* &
                           &                     BS%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                      !
                      AtmPair(CF12)%SP%Cst(ISwitch1,IJ)=2.0d0*Z1!->ISwitch1=6,7
                      AtmPair(CF12)%SP%Cst(ISwitch2,IJ)=2.0d0*Z2!->ISwitch2=7,6
                      !
                      ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                      ! Here I need to add the correction factor for SP shell.
                      ! TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
                   ELSE
                      ! We can skipp out the loop because the primitives are ordered.
                      EXIT
                   ENDIF
                ENDDO
                ! We can skipp out the loop if we did not get any significant primitives.
                IF(JJ.EQ.0) EXIT
                II=II+JJ
             ENDDO
             !
             AtmPair(CF12)%SP%L=II
             !
             ! We reorder the atomic positions if Type2 > Type1.
             ! Needed for the integral evaluations.
             IF(Switch) THEN
                AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
                AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
                AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
                !
                ! AtmInfo must be related to the atoms in the working cell ONLY. 
                ! Then we add the PBC's to have the right atomic position.
                AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
                AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
                AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ
             ELSE
                !
                ! AtmInfo must be related to the atoms in the working cell ONLY. 
                ! Then we add the PBC's to have the right atomic position.
                AtmPair(CF12)%SP%AtmInfo%Atm1X=AtmInfo%Atm2X+RX
                AtmPair(CF12)%SP%AtmInfo%Atm1Y=AtmInfo%Atm2Y+RY
                AtmPair(CF12)%SP%AtmInfo%Atm1Z=AtmInfo%Atm2Z+RZ
                AtmPair(CF12)%SP%AtmInfo%Atm2X=AtmInfo%Atm1X
                AtmPair(CF12)%SP%AtmInfo%Atm2Y=AtmInfo%Atm1Y
                AtmPair(CF12)%SP%AtmInfo%Atm2Z=AtmInfo%Atm1Z
             ENDIF
             !
#ifdef GONX2_DBUG
             WRITE(*,'(3(A,I3),A,I5)') 'Type1',Type1,' Type2',Type2, &
                  &                    ' IntType',AtmPair(CF12)%SP%IntType
#endif
             !
          ENDDO
       ENDDO
    ENDDO ! NCell
    !
  END SUBROUTINE GetAtomPairG
  !
  !
  INTEGER FUNCTION GetNonNeglCell(List,DFac,Trsh)
!H---------------------------------------------------------------------------------
!H INTEGER FUNCTION GetNonNeglCell(List,DFac,Trsh)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(ANode2), POINTER    :: List
    REAL(DOUBLE), INTENT(IN) :: DFac,Trsh
    !-------------------------------------------------------------------
    INTEGER                  :: I
    !-------------------------------------------------------------------
    !
    GetNonNeglCell=0
    !
    DO I=1,List%NCell
       IF(List%SqrtInt(I)*DFac.LE.Trsh) EXIT
       GetNonNeglCell=GetNonNeglCell+1
    ENDDO
    !
  END FUNCTION GetNonNeglCell
  !
  !
  SUBROUTINE GetColIdx(At,D,ColIdx)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetColIdx(At,D,ColIdx)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER       , INTENT(IN ) :: At
    TYPE(BCSR)    , INTENT(IN ) :: D
    TYPE(INT_VECT), INTENT(OUT) :: ColIdx
    !-------------------------------------------------------------------
    INTEGER                     :: Ci
    !-------------------------------------------------------------------
    !
    CALL INT_VECT_EQ_INT_SCLR(NAtoms,ColIdx%I(1),0)
    !
    DO Ci=D%RowPt%I(At),D%RowPt%I(At+1)-1
       ColIdx%I(D%ColPt%I(Ci))=Ci
    ENDDO
    !
  END SUBROUTINE GetColIdx
  !
  !
END MODULE GONX2ComputDK










