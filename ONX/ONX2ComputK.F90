!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
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

!use VScratchB
!USE int_interface
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
    TYPE(FASTMAT)              , POINTER :: DFM,KxFM
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
    INTEGER                    :: LocNInt,IntType,NSMat
    INTEGER                    :: OT,OAL,OBL,LDAL,LDBL,GOAL,GOBL
    INTEGER                    :: OA,OB,OC,OD,LDA,LDB,LDC,LDD,GOA,GOB,GOC,GOD
#ifdef ONX2_PARALLEL
    REAL(DOUBLE)               :: TmBeg,TmEnd
#endif
    REAL(DOUBLE)               :: Dcd,NInts,NIntsTot
    !-------------------------------------------------------------------
    REAL(DOUBLE), DIMENSION(MaxFuncPerAtmBlk**4) :: C!,cc
    REAL(DOUBLE), DIMENSION(MaxShelPerAtmBlk**2) :: DMcd
    TYPE(AtomPr), DIMENSION(:), ALLOCATABLE :: ACAtmPair,BDAtmPair
    !-------------------------------------------------------------------
    REAL(DOUBLE), EXTERNAL :: DGetAbsMax
#ifdef ONX2_PARALLEL
    REAL(DOUBLE), EXTERNAL :: MondoTimer
#endif

    LOGICAL :: UpperT
    !-------------------------------------------------------------------
    integer :: i,isize,aalen,bblen,cclen,ddlen
    real(double) :: t1,t2,tt,fac
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
             BDAtmPair(MaxShelPerAtmBlk**2*CS_OUT%NCells),STAT=iErr)
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
    !Check for wrong InvBoxSh and BoxShape
    IF(    ABS(GMc%PBC%InvBoxSh%D(2,1)).GT.1D-15.OR.ABS(GMc%PBC%InvBoxSh%D(3,1)).GT.1D-15.OR.&
           ABS(GMc%PBC%InvBoxSh%D(3,2)).GT.1D-15.OR.ABS(GMc%PBC%BoxShape%D(2,1)).GT.1D-15.OR.&
           ABS(GMc%PBC%BoxShape%D(3,1)).GT.1D-15.OR.ABS(GMc%PBC%BoxShape%D(3,2)).GT.1D-15) THEN
       WRITE(*,*) 'The following guys MUST be ZERO!'
       WRITE(*,*) 'InvBoxSh:',GMc%PBC%InvBoxSh%D(2,1),GMc%PBC%InvBoxSh%D(3,1),GMc%PBC%InvBoxSh%D(3,2)
       WRITE(*,*) 'BoxShape:',GMc%PBC%BoxShape%D(2,1),GMc%PBC%BoxShape%D(3,1),GMc%PBC%BoxShape%D(3,2)
       STOP 'STOP in ONX2ComputK.F90'
    ENDIF
    !Simple check Simple check Simple check Simple check
    !
    LocNInt=0
    NInts=0.0d0
    NIntsTot=0.0d0
    !
#ifdef ONX2_PARALLEL
    NSMat=DFM%NSMat
    P => DFM%Next ! Run over atom C
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       AtC = P%Row
#else
    NSMat=D%NSMat
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
          TmBeg = MondoTimer()
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
#else
          Dcd=DGetAbsMax(NBFC*NBFD,D%MTrix%D(iPtrD))
#endif
          !
#ifdef ONX2_PARALLEL
          CALL GetAbsDenBlk(U%MTrix(1,1),NBFC,NBFD,NSMat,DMcd(1),&
                            BSp%NCFnc%I(KC),BSp%NCFnc%I(KD),     &
                            BSp%LStrt%I(1,KC),BSp%LStop%I(1,KC), &
                            BSp%LStrt%I(1,KD),BSp%LStop%I(1,KD)  )
#else
          CALL GetAbsDenBlk(D%MTrix%D(iPtrD),NBFC,NBFD,NSMat,DMcd(1), &
                            BSp%NCFnc%I(KC),BSp%NCFnc%I(KD),          &
                            BSp%LStrt%I(1,KC),BSp%LStop%I(1,KC),      &
                            BSp%LStrt%I(1,KD),BSp%LStop%I(1,KD)       )
#endif
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
               )
#else
               *(-1d0))
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
             !call cpu_time(t1)
             CALL GetAtomPair(ACAtmInfo,AtAList,ACAtmPair,BSc,BSp,CS_OUT)
             !call cpu_time(t2)
             !tt=tt+t2-t1
             !
             AtBList=>AtBListTmp
             !
             RnOvB: DO ! Run over AtB
                AtB=AtBList%Atom

                IF(.NOT.NoSym)THEN
                   UpperT=AtB.LE.AtA
                ELSE
                   UpperT=.TRUE.
                ENDIF
                !

                IF(UpperT)THEN

                   BDAtmInfo%NFPair=GetNonNFPair(AtBList,AtAList%RInt(1)*Dcd,Thresholds%TwoE &
#ifdef GTRESH
                     )
#else
                     *(-1d0))
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
                   !call cpu_time(t1)
                   CALL GetAtomPair(BDAtmInfo,AtBList,BDAtmPair,BSc,BSp,CS_OUT)
                   !call cpu_time(t2)
                   !tt=tt+t2-t1
                   !
                   NIntBlk=NBFA*NBFB*NBFC*NBFD
                   !
                   CALL DBL_VECT_EQ_DBL_SCLR(NIntBlk,C(1),0.0d0)

!cc=0d0
!c=0d0
                   !CALL DBL_VECT_EQ_DBL_SCLR(NIntBlk,C(1),BIG_DBL)
                   !
                   RnOvFAC: DO iFAC=1,ACAtmInfo%NFPair
#ifdef GTRESH
                      IF(Dcd*AtAList%RInt(iFAC)*AtBList%RInt(1).LT.Thresholds%TwoE) EXIT RnOvFAC
#endif
                      CFA=AtAList%Indx(1,iFAC)
                      CFC=AtAList%Indx(2,iFAC)
                      !
                      ! BraSwitch
                      IF(ACAtmPair(iFAC)%SP%Switch) THEN
                         OAL=OffArrP%I(CFC,KC)-1;LDAL=NBFA*NBFB; !A
                         OBL=OffArrC%I(CFA,KA)  ;LDBL=1        ; !C
                      ELSE
                         OAL=OffArrC%I(CFA,KA)  ;LDAL=1        ; !A
                         OBL=OffArrP%I(CFC,KC)-1;LDBL=NBFA*NBFB; !C
                      ENDIF
                      !
                      RnOvFBD: DO iFBD=1,BDAtmInfo%NFPair
                         CFB=AtBList%Indx(1,iFBD)
                         CFD=AtBList%Indx(2,iFBD)
#ifdef GTRESH
                         IF(DMcd((CFC-1)*NCFncD+CFD)* &
                                AtAList%RInt(iFAC)*AtBList%RInt(iFBD)>Thresholds%TwoE) THEN
                            IF(Dcd*AtAList%RInt(iFAC)*AtBList%RInt(iFBD)<Thresholds%TwoE) EXIT RnOvFBD
#endif
                            !
                            ! KetSwitch
                            IF(BDAtmPair(iFBD)%SP%Switch) THEN
                               OC=OffArrP%I(CFD,KD)-1;LDC=NBFA*NBFB*NBFC; !B
                               OD=OffArrC%I(CFB,KB)-1;LDD=NBFA          ; !D
                            ELSE
                               OC=OffArrC%I(CFB,KB)-1;LDC=NBFA          ; !B
                               OD=OffArrP%I(CFD,KD)-1;LDD=NBFA*NBFB*NBFC; !D
                            ENDIF
                            !
                            ! BraKetSwitch
                            IF(ACAtmPair(iFAC)%SP%IntType.LT.BDAtmPair(iFBD)%SP%IntType) THEN
                               OA =OC ;OC =OAL ;
                               LDA=LDC;LDC=LDAL;
                               !
                               OB =OD ;OD =OBL ;
                               LDB=LDD;LDD=LDBL;
                            ELSE
                               OA =OAL ;OB =OBL ;
                               LDA=LDAL;LDB=LDBL;
                            ENDIF
                            !
                            ! Compute integral type.
                            IntType=ACAtmPair(iFAC)%SP%IntType*10000+BDAtmPair(iFBD)%SP%IntType
                            !
                            INCLUDE 'ERIInterfaceB.Inc'
!!$                            INCLUDE '/n/srv2/vweber/MONDO/ONX2/ERIInterface.Inc'


!if(maxval(abs(c-cc)).gt.1d-10) goto 1000

                            NInts=NInts+DBLE(LocNInt)

#ifdef GTRESH
                         ENDIF
#endif
                      ENDDO RnOvFBD
                   ENDDO RnOvFAC
                   !


!goto 1001
!1000 continue
!write(*,*) 'IntType',IntType
!                   CALL Print2E(AtA,AtC,AtB,AtD,GMc,BSc,GMp,BSp,C)
!                   CALL Print2E(AtA,AtC,AtB,AtD,GMc,BSc,GMp,BSp,cc)
!stop 9999
!1001 continue


                   !CALL Print2E(AtA,AtC,AtB,AtD,GMc,BSc,GMp,BSp,C)
                   !
                   ! Get address for Kx and digest the block of integral.
#ifdef ONX2_PARALLEL
                   if(atb.le.0) stop 'In ComputK'
                   V => InsertSRSTNode(Q%RowRoot,AtB)
                   IF(.NOT.ASSOCIATED(V%MTrix)) THEN
                      ALLOCATE(V%MTrix(NBFA,NBFB*NSMat),STAT=MemStatus)
                      CALL DBL_VECT_EQ_DBL_SCLR(NBFA*NBFB*NSMat,V%MTrix(1,1),0.0d0)
                   ENDIF
                   !
                   ! The variable -Fac- is needed to get the
                   ! right full K matrix when filling out in ONX.
                   Fac=-1.0d0
                   iF(AtA.EQ.AtB) Fac=-0.5d0
                   IF(NSMat.EQ.1) THEN
                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,Fac,C(1), &
                                 NBFA*NBFB,U%MTrix(1,1),1,1.0d0, &
                                 V%MTrix(1,1),1)
                      !CALL DGEMM_NNC(NBFA*NBFB,NBFC*NBFD,1,Fac,1D0, &
                      !               C(1),U%MTrix(1,1),V%MTrix(1,1))
                   ELSE
                      CALL DGEMM('N','N',NBFA*NBFB,NSMat,NBFC*NBFD,Fac,C(1), &
                                 NBFA*NBFB,U%MTrix(1,1),NBFC*NBFD,1.0d0, &
                                 V%MTrix(1,1),NBFA*NBFB)
                   ENDIF

#else
                   CALL GetAdrB(AtA,AtB,Ind,Kx,0)
                   iPtrK = Kx%BlkPt%I(Ind)
                   !
                   IF(NSMat.EQ.1) THEN


                      CALL DGEMV('N',NBFA*NBFB,NBFC*NBFD,-1.0d0,C(1), &
                                 NBFA*NBFB,D%MTrix%D(iPtrD),1,1.0d0, &
                                 Kx%MTrix%D(IPtrK),1)

!                      WRITE(*,*)AtA,AtB,AtC,AtD,D%MTrix%D(iPTRD),C(1),Kx%MTrix%D(IPtrK)/2D0
!33                    FORMAT(4(I2,", "),3(D12.6,", "))

                   ELSE
                      CALL DGEMM('N','N',NBFA*NBFB,NSMat,NBFC*NBFD,-1.0d0,C(1), &
                                 NBFA*NBFB,D%MTrix%D(iPtrD),NBFC*NBFD,1.0d0, &
                                 Kx%MTrix%D(IPtrK),NBFA*NBFB)
                   ENDIF
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
          TmEnd = MondoTimer()
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
#ifdef ONX2_PARALLEL
    ! We need to remove the empty Row we may have created.
    CALL FASTMAT_RmEmptyRow(KxFM)
#endif
    !
!!$    WRITE(*,100) NInts,CS_OUT%NCells**2*DBLE(NBasF)**4, &
!!$                 NInts/(CS_OUT%NCells**2*DBLE(NBasF)**4)*100D0
!!$100 FORMAT(' NInts = ',E8.2,' NIntTot = ',E8.2,' Ratio = ',E8.2,'%')
    !

  END SUBROUTINE ComputK
  !
  !
  SUBROUTINE Print2E(AtA,AtC,AtB,AtD,GMc,BSc,GMp,BSp,C)
!H---------------------------------------------------------------------------------
!H SUBROUTINE Print2E(AtA,AtC,AtB,AtD,BSc,BSp,C)
!H
!H Print the 2-E integral (ac|bd) matrix strored in the form (ab)x(cd).
!H !! Doesn't work when SCFActn=='BasisSetSwitch' (need OffSc and OffSp) !!
!H---------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(CRDS)   :: GMc,GMp
    TYPE(BSET)   :: BSc,BSp
    REAL(DOUBLE) :: C(*)
    INTEGER      :: AtA,AtB,AtC,AtD
    !-------------------------------------------------------------------
    INTEGER      :: KA,KB,KC,KD,NBFA,NBFB,NBFC,NBFD
    INTEGER      :: FA,FB,FC,FD,IJKL,iiii
    !-------------------------------------------------------------------
    !
    WRITE(*,'(4(A,I3))') 'AtA',AtA,' AtC',AtC,' AtB',AtB,' AtD',AtD
    !
    KA=GMc%AtTyp%I(AtA)
    NBFA=BSc%BfKnd%I(KA)
    KB=GMc%AtTyp%I(AtB)
    NBFB=BSc%BfKnd%I(KB)
    !
    KC=GMp%AtTyp%I(AtC)
    NBFC=BSp%BfKnd%I(KC)
    KD=GMp%AtTyp%I(AtD)
    NBFD=BSp%BfKnd%I(KD)
    !
    IJKL=0
    iiii=0
    DO FD=1,NBFD
    DO FC=1,NBFC
    DO FB=1,NBFB
    DO FA=1,NBFA
       IJKL=IJKL+1
       IF(ABS(C(IJKL)).GT.1D-15) THEN
          iiii=iiii+1
          WRITE(*,'(A,4I3,I8,A,E23.15)') 'Int(', &
                 OffS%I(AtA)+FA-1, &
                 OffS%I(AtC)+FC-1, &
                 OffS%I(AtB)+FB-1, &
               ! OffS%I(AtD)+FD-1,')=',C(IJKL)
                 OffS%I(AtD)+FD-1,IJKL,')=',C(IJKL)
          !
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !
    write(*,'(2(A,I5))') '#Int_tot',NBFA*NBFC*NBFB*NBFD,' #IntN0',iiii
    !
END SUBROUTINE Print2E


#ifdef ONX2_PARALLEL
  SUBROUTINE FASTMAT_RmEmptyRow(A)
    TYPE(FASTMAT), POINTER :: A
    TYPE(FASTMAT), POINTER :: P,Q
    TYPE(SRST   ), POINTER :: U
    NULLIFY(P,Q,U)
    P=>A%Next
    Q=>A
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       U=>P%RowRoot
       IF(.NOT.ASSOCIATED(U%Left).AND..NOT.ASSOCIATED(U%Right)) THEN
          !write(*,*) 'The Row',P%Row,' MyID',MyID
          Q%Next=>P%Next
          CALL Delete_SRST_1(U)
          DEALLOCATE(P)
          P=>Q%Next
       ELSE
          Q=>P
          P=>P%Next
       ENDIF
    ENDDO
  END SUBROUTINE FASTMAT_RmEmptyRow
#endif
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
    INTEGER      :: MinL1,MaxL1,Type1,MinL2,MaxL2,Type2
    INTEGER      :: StopL1,StartL1,StopL2,StartL2,iNFPair
    REAL(DOUBLE) :: Z1,Z2,Expt,InvExpt,R12,XiR12,RX,RY,RZ,Cnt
    LOGICAL      :: Switch
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
       !
       Switch=Type1.LT.Type2
       AtmPair(iNFPair)%SP%Switch=Switch
       IF(Switch) THEN
          AtmPair(iNFPair)%SP%IntType=Type2*100+Type1
       ELSE
          AtmPair(iNFPair)%SP%IntType=Type1*100+Type2
       ENDIF
       !
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
                !
                IF((Type1.NE.2.AND.Type2==2).OR.(Type2.NE.2.AND.Type1==2))THEN
                   AtmPair(iNFPair)%SP%Cst(6,IJ)=BSc%CCoef%D(StartL1,  I1,CF1,AtmInfo%K1) * &
                                                 BSp%CCoef%D(StartL2,  I2,CF2,AtmInfo%K2)/Cnt
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
