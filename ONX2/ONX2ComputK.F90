MODULE ONX2ComputK
!H=================================================================================
!H MODULE ONX2ComputK
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB 
!H
!H  PRIVATE:
!H  o SUB 
!H
!H  OPTIONS:
!H  DEBUGING: Use -DONX2_DBUG to print some stuff.
!H  INFO    : Use -DONX2_INFO to print some stuff.
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONX2DataType
  USE InvExp
  USE ONXParameters
  !
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
#endif
  !
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: ComputK
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: GetNonNFPair
  PRIVATE :: GetAtomPair
  !
CONTAINS
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE ComputK(DFM,KxFM,ListC,ListD,OffArrC,OffArrP,GMc,BSc,GMp,BSp,CS_OUT)
#else
  SUBROUTINE ComputK(D,Kx,ListC,ListD,OffArrC,OffArrP,GMc,BSc,GMp,BSp,CS_OUT)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE ComputK(D,Kx,ListC,ListD,OffArrC,OffArrP,GMc,BSc,GMp,BSp,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    !
    USE ONXGet, ONLY: GetAdrB
    !
    IMPLICIT NONE
    !
    !
    ! Kab = Dcd*(ac(R)|bd(R'))
    !
    !-------------------------------------------------------------------
#ifdef ONX2_PARALLEL
    TYPE(fastmat)              , POINTER :: DFM,KxFM
    TYPE(FASTMAT)              , POINTER :: P,Q
    TYPE(SRST   )              , POINTER :: U,V
#else
    TYPE(BCSR)                           :: D,Kx
#endif
    TYPE(INT_RNK2)                       :: OffArrC,OffArrP
    TYPE(CList ) , DIMENSION(:), POINTER :: ListC,ListD
    TYPE(CRDS)                           :: GMc,GMp
    TYPE(BSET)                           :: BSc,BSp
    TYPE(CellSet)                        :: CS_OUT
    !-------------------------------------------------------------------
    TYPE(ANode ), POINTER      :: AtAListTmp,AtAList,AtBListTmp,AtBList
    TYPE(AtomInfo)             :: ACAtmInfo,BDAtmInfo
    INTEGER                    :: AtA,AtB,AtC,AtD,KA,KB,KC,KD,CFA,CFB,CFC,CFD
    INTEGER                    :: ci,iPtrD,iPtrK,NBFC,NBFD,NBFA,NBFB
    INTEGER                    :: NIntBlk,iErr,iFAC,iFBD,NCFncD
    INTEGER                    :: Off,Ind
    INTEGER                    :: LocNInt,IntType
#ifdef ONX2_PARALLEL
    REAL(DOUBLE)               :: TmBeg,TmEnd
#endif
    REAL(DOUBLE)               :: Dcd,NInts,NIntsTot
    !-------------------------------------------------------------------
    REAL(DOUBLE), DIMENSION(MaxFuncPerAtmBlk**4) :: C
    REAL(DOUBLE), DIMENSION(MaxShelPerAtmBlk**2) :: DMcd
    TYPE(AtomPr), DIMENSION(:), ALLOCATABLE :: ACAtmPair,BDAtmPair
    !-------------------------------------------------------------------
    TYPE(ONX2OffSt) :: OffSet
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL     :: DGetAbsMax
    !-------------------------------------------------------------------
    integer :: i,isize,aalen,bblen,cclen,ddlen
    real(double) :: t1,t2,tt
    tt=0.0d0
    !
    ! Initialize.
    NULLIFY(AtAListTmp,AtAList,AtBListTmp,AtBList)
#ifdef ONX2_PARALLEL
    NULLIFY(P,Q,U,V)
#endif
    !
    ! Allocate arrays.
    ALLOCATE(ACAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells), &
         &   BDAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells),STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ComputK: Allocation problem.')
    !
    !
    !Simple check Simple check Simple check Simple check
    isize=0
    do i=1,natoms
       isize=MAX(isize,max(BSc%BfKnd%I(GMc%AtTyp%I(i)),BSp%BfKnd%I(GMp%AtTyp%I(i))))
       IF((max(BSc%BfKnd%I(GMc%AtTyp%I(i)),BSp%BfKnd%I(GMp%AtTyp%I(i))))**4.GT.SIZE(C)) THEN
          write(*,*) 'size',(max(BSc%BfKnd%I(GMc%AtTyp%I(i)),BSp%BfKnd%I(GMp%AtTyp%I(i))))**4
          STOP 'Incrase the size of C'
       ENDIF
    enddo
    !write(*,*) 'size C=',isize**4
    if(CS_OUT%NCells.GT.SIZE(ACAtmPair)) then
       write(*,*) 'size(ACAtmPair)',size(ACAtmPair),'.LT.',CS_OUT%NCells
       STOP 'Incrase the size of ACAtmPair and BDAtmPair'
    endif
    !write(*,*) 'CS_OUT%NCells=',CS_OUT%NCells
    !write(*,*) 'CS_OUT%CellCarts=',CS_OUT%CellCarts%D(1,:)
    !write(*,*) 'Thresholds%Dist',Thresholds%Dist,' Thresholds%TwoE',Thresholds%TwoE
    !Simple check Simple check Simple check Simple check
    !
    LocNInt=0
    NInts=0.0d0
    NIntsTot=0.0d0
    !
#ifdef ONX2_PARALLEL
    P => DFM%Next ! Run over atom C
    DO                               
       IF(.NOT.ASSOCIATED(P)) EXIT   
       AtC = P%Row
#else
    DO AtC=1,NAtoms ! Run over AtC.
#endif
       KC=GMp%AtTyp%I(AtC)
       NBFC=BSp%BfKnd%I(KC)
       ACAtmInfo%Atm2X=GMp%Carts%D(1,AtC)
       ACAtmInfo%Atm2Y=GMp%Carts%D(2,AtC)
       ACAtmInfo%Atm2Z=GMp%Carts%D(3,AtC)
       ACAtmInfo%K2=KC
       !
       ! Get AtA List.
       AtAListTmp=>ListC(AtC)%GoList
       !
#ifdef ONX2_PARALLEL
       U => P%RowRoot ! Run over AtD
       DO                                
          IF(.NOT.ASSOCIATED(U)) EXIT    
          IF(U%L.NE.U%R) THEN            
             U => U%Next                 
             CYCLE                       
          ENDIF
          AtD = U%L                      
          ! Set Time.                    
          TmBeg = MPI_WTIME()            
#else
       DO ci=D%RowPt%I(AtC),D%RowPt%I(AtC+1)-1 ! Run over AtD
          AtD=D%ColPt%I(ci)
          iPtrD=D%BlkPt%I(ci)
#endif
          KD=GMp%AtTyp%I(AtD)
          NBFD=BSp%BfKnd%I(KD)
          NCFncD=BSp%NCFnc%I(KD)

          BDAtmInfo%Atm2X=GMp%Carts%D(1,AtD)
          BDAtmInfo%Atm2Y=GMp%Carts%D(2,AtD)
          BDAtmInfo%Atm2Z=GMp%Carts%D(3,AtD)
          BDAtmInfo%K2=KD
          !
          ! Get max of the block density matrix.
#ifdef ONX2_PARALLEL
          Dcd=DGetAbsMax(NBFC*NBFD,U%MTrix(1,1))
          !write(*,*) 'Dcd',Dcd,'AtC=',AtC,'AtD=',AtD,'MyID',MyID
#else
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
#endif
          !
#ifdef ONX2_PARALLEL
          CALL GetAbsDenBlk(U%MTrix(1,1),NBFC,NBFD,DMcd(1),      &
               &            BSp%NCFnc%I(KC),BSp%NCFnc%I(KD),     &
               &            BSp%LStrt%I(1,KC),BSp%LStop%I(1,KC), &
               &            BSp%LStrt%I(1,KD),BSp%LStop%I(1,KD)  )
#else
          CALL GetAbsDenBlk(D%MTrix%D(iPtrD),NBFC,NBFD,DMcd(1),  &
               &            BSp%NCFnc%I(KC),BSp%NCFnc%I(KD),     &
               &            BSp%LStrt%I(1,KC),BSp%LStop%I(1,KC), &
               &            BSp%LStrt%I(1,KD),BSp%LStop%I(1,KD)  )
#endif
          !
#ifdef ONX2_DBUG
          WRITE(*,*) 'Max(Dcd)=',Dcd
#endif
          !
          ! Get AtB List.
          AtBListTmp=>ListD(AtD)%GoList
          !
          AtAList=>AtAListTmp
          !
          RnOvA: DO ! Run over AtA
             AtA=AtAList%Atom
             KA=GMc%AtTyp%I(AtA)
             NBFA=BSc%BfKnd%I(KA)
             !
             ACAtmInfo%NFPair=GetNonNFPair(AtAList,AtBListTmp%RInt(1)*Dcd,Thresholds%TwoE &
#ifdef GTRESH
             & )
#else
             & *(-1d0))
#endif
             IF(ACAtmInfo%NFPair.EQ.0) EXIT RnOvA
             !
             ! Find the row in Kx.
#ifdef ONX2_PARALLEL
             Q => FindFastMatRow_1(KxFM,AtA)
#endif
             !
             ACAtmInfo%Atm1X=GMc%Carts%D(1,AtA)
             ACAtmInfo%Atm1Y=GMc%Carts%D(2,AtA)
             ACAtmInfo%Atm1Z=GMc%Carts%D(3,AtA)
             ACAtmInfo%K1=KA
             !
             ACAtmInfo%Atm12X=ACAtmInfo%Atm1X-ACAtmInfo%Atm2X
             ACAtmInfo%Atm12Y=ACAtmInfo%Atm1Y-ACAtmInfo%Atm2Y
             ACAtmInfo%Atm12Z=ACAtmInfo%Atm1Z-ACAtmInfo%Atm2Z
             !
             ! Get atom pair for BD.
call cpu_time(t1)
             CALL GetAtomPair(ACAtmInfo,AtAList,ACAtmPair,BSc,BSp,CS_OUT)
call cpu_time(t2)
tt=tt+t2-t1
             !
             AtBList=>AtBListTmp
             !
             RnOvB: DO ! Run over AtB
                AtB=AtBList%Atom 
                !
!WRITE(*,*)' AtA = ',Ata, ' AtB = ',AtB
                IF(AtB.LE.AtA) THEN ! Symmetry of the K matrix
                   !
                   BDAtmInfo%NFPair=GetNonNFPair(AtBList,AtAList%RInt(1)*Dcd,Thresholds%TwoE &
#ifdef GTRESH
                   & )
#else
                   & *(-1d0))
#endif
                   IF(BDAtmInfo%NFPair.EQ.0) EXIT RnOvB
                   !
                   KB=GMc%AtTyp%I(AtB)
                   NBFB=BSc%BfKnd%I(KB)
                   !
                   BDAtmInfo%Atm1X=GMc%Carts%D(1,AtB)
                   BDAtmInfo%Atm1Y=GMc%Carts%D(2,AtB)
                   BDAtmInfo%Atm1Z=GMc%Carts%D(3,AtB)
                   BDAtmInfo%K1=KB
                   !
                   BDAtmInfo%Atm12X=BDAtmInfo%Atm1X-BDAtmInfo%Atm2X
                   BDAtmInfo%Atm12Y=BDAtmInfo%Atm1Y-BDAtmInfo%Atm2Y
                   BDAtmInfo%Atm12Z=BDAtmInfo%Atm1Z-BDAtmInfo%Atm2Z
                   !
                   ! Get atom pair for BD.
call cpu_time(t1)
                   CALL GetAtomPair(BDAtmInfo,AtBList,BDAtmPair,BSc,BSp,CS_OUT)
call cpu_time(t2)
tt=tt+t2-t1
                   !
                   NIntBlk=NBFA*NBFB*NBFC*NBFD
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(NIntBlk,C(1),0.0d0)
                   !CALL DBL_VECT_EQ_DBL_SCLR(NIntBlk,C(1),BIG_DBL)
                   !
                   RnOvFAC: DO iFAC=1,ACAtmInfo%NFPair
#ifdef GTRESH
                      IF(Dcd*AtAList%RInt(iFAC)*AtBList%RInt(1).LT.Thresholds%TwoE) EXIT RnOvFAC
#endif
                      CFA=AtAList%Indx(1,iFAC)
                      AALen=BSc%LStop%I(CFA,KA)-BSc%LStrt%I(CFA,KA)+1
                      CFC=AtAList%Indx(2,iFAC)
                      CCLen=BSc%LStop%I(CFC,KC)-BSc%LStrt%I(CFC,KC)+1
                      OffSet%A=OffArrC%I(CFA,KA)
                      OffSet%C=OffArrP%I(CFC,KC)
                      !
                      RnOvFBD: DO iFBD=1,BDAtmInfo%NFPair 
                         CFB=AtBList%Indx(1,iFBD)
                         BBLen=BSp%LStop%I(CFB,KB)-BSp%LStrt%I(CFB,KB)+1
                         CFD=AtBList%Indx(2,iFBD)
                         DDLen=BSp%LStop%I(CFD,KD)-BSp%LStrt%I(CFD,KD)+1
!#ifdef GTRESH
!                         IF(DMcd((CFC-1)*NCFncD+CFD)* &
!                              & AtAList%RInt(iFAC)*AtBList%RInt(iFBD)>Thresholds%TwoE) THEN
!                            IF(Dcd*AtAList%RInt(iFAC)*AtBList%RInt(iFBD)<Thresholds%TwoE) EXIT RnOvFBD
!#endif
                            OffSet%B=OffArrC%I(CFB,KB)
                            OffSet%D=OffArrP%I(CFD,KD)
                            !
                            ! Compute integral type.
                            IntType=ACAtmPair(iFAC)%SP%IntType*10000+BDAtmPair(iFBD)%SP%IntType

                            INCLUDE 'ERIInterfaceB.Inc'

!                            CALL ShellPrint(NBFA,NBFB,NBFC,NBFD,AALen,BBLen,CCLen,DDLen,  &
!                                            OffSet%A,OffSet%B,OffSet%C,OffSet%D,IntType,C(1))

                            NInts=NInts+DBLE(LocNInt)
!#ifdef GTRESH
!                         ENDIF
!#endif
                      ENDDO RnOvFBD
                   ENDDO RnOvFAC
                   !
                   ! Get address for Kx and digest the block of integral.
#ifdef ONX2_PARALLEL
                   V => InsertSRSTNode(Q%RowRoot,AtB)
                   IF(.NOT.ASSOCIATED(V%MTrix)) THEN
                      ALLOCATE(V%MTrix(NBFA,NBFB),STAT=MemStatus)
                      CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB,V%MTrix(1,1),0.0d0)
                   ENDIF
                   !
                   CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(1), &
                        &     NBFA*NBFB,U%MTrix(1,1),1,1.0d0, &
                        &     V%MTrix(1,1),1)
#else
                   CALL GetAdrB(AtA,AtB,Ind,Kx,0)
                   iPtrK = Kx%BlkPt%I(Ind)
                   !
                   CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(1), &
                        &     NBFA*NBFB,D%MTrix%D(iPtrD),1,1.0d0, &
                        &     Kx%MTrix%D(IPtrK),1)
#endif
                   !
                ENDIF
                !
                IF(.NOT.ASSOCIATED(AtBList%AtmNext)) EXIT RnOvB
                AtBList=>AtBList%AtmNext
             ENDDO RnOvB ! End AtB
             !
             IF(.NOT.ASSOCIATED(AtAList%AtmNext)) EXIT RnOvA
             AtAList=>AtAList%AtmNext
          ENDDO RnOvA ! End AtA
          !
#ifdef ONX2_PARALLEL
          ! Set Time.
          TmEnd = MPI_WTIME()
          !Add Time.
          U%Part = U%Part+TmEnd-TmBeg
          U => U%Next
#endif
       ENDDO ! End AtD
       !
#ifdef ONX2_PARALLEL
       P => P%Next
#endif
    ENDDO ! End AtC
    !
    ! DeAllocate arrays.
    DEALLOCATE(ACAtmPair,BDAtmPair,STAT=iErr)
    IF(iErr.NE.0) CALL Halt('In ComputK: Deallocation problem.')
    !
    write(*,*) 'time',tt
    !
!!$#ifdef ONX2_PARALLEL
!!$    NIntsTot=Reduce(NInts)
!!$    IF(MyID.EQ.ROOT) THEN
!!$       WRITE(*,100) NIntsTot,CS_OUT%NCells**2*DBLE(NBasF)**4, &
!!$            &       NIntsTot/(CS_OUT%NCells**2*DBLE(NBasF)**4)*100D0
!!$    ENDIF
!!$#else
    WRITE(*,100) NInts,CS_OUT%NCells**2*DBLE(NBasF)**4, &
         &       NInts/(CS_OUT%NCells**2*DBLE(NBasF)**4)*100D0
!!$#endif
100 FORMAT(' NInts = ',E8.2,' NIntTot = ',E8.2,' Ratio = ',E8.2,'%')
    !
  END SUBROUTINE ComputK
  !
  !
  INTEGER FUNCTION GetNonNFPair(List,DFac,Trsh)
!H---------------------------------------------------------------------------------
!H INTEGER FUNCTION GetNonNFPair(List,DFac,Trsh)
!H
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(ANode ), POINTER    :: List
    REAL(DOUBLE), INTENT(IN) :: DFac,Trsh
    !-------------------------------------------------------------------
    INTEGER                  :: I
    !-------------------------------------------------------------------
    !
    GetNonNFPair=0
    !
    DO I=1,List%NFPair
       IF(List%RInt(I)*DFac.LE.Trsh) EXIT
       GetNonNFPair=GetNonNFPair+1
    ENDDO
    !
  END FUNCTION GetNonNFPair
  !
  SUBROUTINE GetAtomPair(AtmInfo,List,AtmPair,BSc,BSp,CS_OUT)
!H---------------------------------------------------------------------------------
!H SUBROUTINE GetAtomPair(AtmInfo,List,AtmPair,BSc,BSp,CS_OUT)
!H
!H---------------------------------------------------------------------------------
    USE Thresholding
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(AtomInfo)               :: AtmInfo
    TYPE(AtomPr  ), DIMENSION(:) :: AtmPair
    TYPE(BSET)                   :: BSc,BSp
    TYPE(CellSet)                :: CS_OUT
    TYPE(ANode   ), POINTER      :: List
    !-------------------------------------------------------------------
    INTEGER      :: CF1,CF2,I1,I2,II,JJ,IJ,Cell
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2,IntType
    INTEGER      :: StopL1,StartL1,StopL2,StartL2,iNFPair
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ,Cnt
    !-------------------------------------------------------------------
    !
    DO iNFPair=1,AtmInfo%NFPair
       CF1 =List%Indx(1,iNFPair) !A,B
       CF2 =List%Indx(2,iNFPair) !C,D
       Cell=List%Indx(3,iNFPair)
       RX=CS_OUT%CellCarts%D(1,Cell)
       RY=CS_OUT%CellCarts%D(2,Cell)
       RZ=CS_OUT%CellCarts%D(3,Cell)
       !
       ! AtmInfo must be related to the atoms in the working cell ONLY.
       ! Then we add the PBC's to have the right interatomic distance.
       R12=(AtmInfo%Atm12X-RX)**2+(AtmInfo%Atm12Y-RY)**2+(AtmInfo%Atm12Z-RZ)**2
       !
       MinL1=BSc%ASymm%I(1,CF1,AtmInfo%K1)
       MaxL1=BSc%ASymm%I(2,CF1,AtmInfo%K1)
       Type1=MaxL1*(MaxL1+1)/2+MinL1+1
       StartL1=BSc%LStrt%I(CF1,AtmInfo%K1)
       StopL1=BSc%LStop%I(CF1,AtmInfo%K1)
       !
       MinL2=BSp%ASymm%I(1,CF2,AtmInfo%K2)
       MaxL2=BSp%ASymm%I(2,CF2,AtmInfo%K2)
       Type2=MaxL2*(MaxL2+1)/2+MinL2+1
       StartL2=BSp%LStrt%I(CF2,AtmInfo%K2)
       StopL2=BSp%LStop%I(CF2,AtmInfo%K2)
       !>>>>>
       IntType=Type1*100+Type2
       AtmPair(iNFPair)%SP%IntType=IntType !Type1*100+Type2
       !<<<<<
       II=0
       !
       ! We assume the primitives are ordered (exponants in decressing order).
       DO I1=BSc%NPFnc%I(CF1,AtmInfo%K1),1,-1
          Z1=BSc%Expnt%D(I1,CF1,AtmInfo%K1)
          JJ=0
          !
          ! We assume the primitives are ordered (exponants in decressing order).
          DO I2=BSp%NPFnc%I(CF2,AtmInfo%K2),1,-1
             Z2=BSp%Expnt%D(I2,CF2,AtmInfo%K2)
             Expt=Z1+Z2
             InvExpt=1.0d0/Expt
             XiR12=Z2*Z1*InvExpt*R12
             IF(XiR12<PrimPairDistanceThreshold) THEN
                JJ=JJ+1
                IJ=JJ+II
                Cnt=BSc%CCoef%D(StopL1,I1,CF1,AtmInfo%K1)*BSp%CCoef%D(StopL2,I2,CF2,AtmInfo%K2)
                AtmPair(iNFPair)%SP%Cst(1,IJ)=Expt
                !
                ! AtmInfo must be related to the atoms in the working cell ONLY.
                ! Then we add the PBC's to have the right atomic position.
                AtmPair(iNFPair)%SP%Cst(2,IJ)=(Z1*AtmInfo%Atm1X+Z2*(AtmInfo%Atm2X+RX))*InvExpt
                AtmPair(iNFPair)%SP%Cst(3,IJ)=(Z1*AtmInfo%Atm1Y+Z2*(AtmInfo%Atm2Y+RY))*InvExpt
                AtmPair(iNFPair)%SP%Cst(4,IJ)=(Z1*AtmInfo%Atm1Z+Z2*(AtmInfo%Atm2Z+RZ))*InvExpt
                AtmPair(iNFPair)%SP%Cst(5,IJ)=5.914967172796D0*EXP(-XiR12)*InvExpt*Cnt
                IF((Type1.NE.2.AND.Type2==2).OR.(Type2.NE.2.AND.Type1==2))THEN
                   AtmPair(iNFPair)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,I2,CF2,AtmInfo%K2)/Cnt
                   AtmPair(iNFPair)%SP%Cst(7,IJ)=BIG_DBL
                   AtmPair(iNFPair)%SP%Cst(8,IJ)=BIG_DBL
                ELSEIF(Type1==2.AND.Type2==2)THEN
                   AtmPair(iNFPair)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                   AtmPair(iNFPair)%SP%Cst(7,IJ)=BSc%CCoef%D(StartL1+1,I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
                   AtmPair(iNFPair)%SP%Cst(8,IJ)=BSc%CCoef%D(StartL1  ,I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2+1,I2,CF2,AtmInfo%K2)/Cnt
                ELSE
                   AtmPair(iNFPair)%SP%Cst(6,IJ)=BIG_DBL
                   AtmPair(iNFPair)%SP%Cst(7,IJ)=BIG_DBL
                   AtmPair(iNFPair)%SP%Cst(8,IJ)=BIG_DBL
                ENDIF
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
       AtmPair(iNFPair)%SP%L=II
       !
       ! We reorder the atomic positions if Type2 > Type1.
       ! Needed for the integral evaluations.
       IF(Type1.GE.Type2) THEN
          AtmPair(iNFPair)%SP%AtmInfo%Atm1X=AtmInfo%Atm1X
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Y=AtmInfo%Atm1Y
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Z=AtmInfo%Atm1Z
          !
          ! AtmInfo must be related to the atoms in the working cell ONLY.
          ! Then we add the PBC's to have the right atomic position.
          AtmPair(iNFPair)%SP%AtmInfo%Atm2X=AtmInfo%Atm2X+RX
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Y=AtmInfo%Atm2Y+RY
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Z=AtmInfo%Atm2Z+RZ
       ELSE
          !
          ! AtmInfo must be related to the atoms in the working cell ONLY.
          ! Then we add the PBC's to have the right atomic position.
          AtmPair(iNFPair)%SP%AtmInfo%Atm1X=AtmInfo%Atm2X+RX
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Y=AtmInfo%Atm2Y+RY
          AtmPair(iNFPair)%SP%AtmInfo%Atm1Z=AtmInfo%Atm2Z+RZ
          AtmPair(iNFPair)%SP%AtmInfo%Atm2X=AtmInfo%Atm1X
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Y=AtmInfo%Atm1Y
          AtmPair(iNFPair)%SP%AtmInfo%Atm2Z=AtmInfo%Atm1Z
       ENDIF
    ENDDO
  END SUBROUTINE GetAtomPair
  !
  !
END MODULE ONX2ComputK

