!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
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
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------
!    FAST O(N Lg N) COMPUTATION OF GRADIENTS OF THE COULOMB ENERGY 
!    WRT TO NUCLEAR COORDINATES
!    Author: Matt Challacombe
!==============================================================================

#include "MondoConfig.h"

PROGRAM JForce
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree
  USE PBCFarField
  USE BlokTrPdJ
  USE SetXYZ 

  IMPLICIT NONE

  TYPE(BCSR)                     :: P,PT,P2
  TYPE(AtomPair)                 :: Pair
  TYPE(HGLL),POINTER             :: RhoHead
  TYPE(DBL_VECT)                 :: Frc,JFrc
  INTEGER                        :: AtA,AtB,A1,A2,B1,B2,MA,NB,MN1,MN,JP,Q
  REAL(DOUBLE)                   :: JFrcChk,ETot
  REAL(DOUBLE)                   :: JFORCE_TotalTime_Start
  CHARACTER(LEN=6),PARAMETER     :: Prog='JForce'
  INTEGER                        :: NLink,NC,I,J,K
  REAL(DOUBLE),DIMENSION(3)      :: A,B,nlm
  REAL(DOUBLE),DIMENSION(15)     :: F_nlm
  TYPE(DBL_RNK2)                 :: LatFrc_J,LatFrc_J_PFF,LatFrc_J_Dip
  LOGICAL                        :: NoWrap=.FALSE. ! WRAPPING IS ON
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,PMat2Use
#if defined(PARALLEL_CLONES)
  INTEGER                        :: oldClone, rank
#endif

  JFORCE_TotalTime_Start=MTimer()
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)

  PMat2Use=TrixFile('D',Args,1) ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  CALL Get(P,TRIM(PMat2Use))
  !
  ! Set thresholds local to QCTC (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
  ! RhoHead is the start of a linked density list
  ALLOCATE(RhoHead)
  RhoHead%LNum=0

  CALL MakeRhoList(GM,BS,P,NLink,RhoHead,'JForce',NoWrap_O=NoWrap)
  ! Add in the nuclear charges
  CALL AddNukes(GM,RhoHead)
  NLink=NLink+GM%NAtms
  ! Load density into arrays and delete the linked list
  CALL Collate(GM,RhoHead,Rho,'JForce',RhoPoles,NLink)

!!!  CALL DeleteHGLL(RhoHead)

  MaxPFFFEll=GM%PBC%PFFMaxEll
  ! Local expansion order of the multipoles to use in the tree
  MaxPoleEll=MIN(2*(BS%NASym+4),MaxPFFFEll)
  IF(MaxPoleEll<2*(BS%NASym+1)) &
     CALL Halt('Bombed in QCTC. Please set PFFMaxEll larger ')
  ! Allocate some memory for bra HG shenanigans 
  CALL NewBraBlok(BS)
  ! Find the total energy from past calculations
  ! Over-ride initialization to big_dbl in punchhdf, so that
  ! PFF is computed correctly in PBCFarFieldSetUp
  IF(SCFCycl=='0'.AND.CurGeom=='1')THEN
     ETot=1D2
  ELSE
     CALL Get(Etot,'Etot')
  ENDIF
  ! For now, set PFFFEll to be the largest expansion length.  Reset later in FarFieldSetUp
  MaxPFFFEll=MAX(MaxPoleEll,MaxPFFFEll)

  IF(NoWrap)THEN
     Mssg=ProcessName('JForce','No wrap')
  ELSE
     Mssg=ProcessName('JForce','Wrapping on')
  ENDIF

!  Mssg=TRIM(Mssg)//' Cluster Size = '//TRIM(IntToChar(MinCluster))//', MaxPFFEll = '//TRIM(IntToChar(MaxPFFFEll))
!  WRITE(*,*)TRIM(Mssg)

  ! Initialize addressing for tensor contraction loops
  CALL TensorIndexingSetUp()
  ! Setup global arrays for computation of multipole tensors ...
  CALL MultipoleSetUp()
  ! Initialize some counters
  MaxTier=0
  RhoLevel=0
  PoleNodes=0
  ! Initialize the root node   
  CALL NewPoleNode(PoleRoot,0)
  ! Initialize the auxiliary density arrays
  CALL InitRhoAux
  ! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
  ! Set up the crystal field, compute related energies and forces
  CALL New(LatFrc_J,(/3,3/))
  CALL New(LatFrc_DWRAP,(/3,3/))

  CALL PBCFarFieldSetUp(GM,Rho,'JForce',MaxPFFFEll,ETot,LatFrc_J)

  ! Delete the auxiliary density arrays
  CALL DeleteRhoAux
  ! Delete the Density
  CALL Delete(Rho)
  ! Some allocations 
  CALL New(JFrc,3*GM%Natms)
  JFrc%D(:)=Zero
  ALLOCATE(TempHerm%Coef(1:HGLen))
  ALLOCATE(NNearCount(1:CS_IN%NCells))
  ! TreeCode counters
  NInts=0
  NPrim=0
  NFarAv=0
  NNearAv=0
  JWalk_Time=0D0
  Integral_Time=0D0
  Multipole_Time=0D0
  NNearCount=0D0
  ! Rescale the density matrix for U/G theory. ??????? WTF ??????
  IF(P%NSMat.GT.1) CALL DSCAL(P%NNon0,0.5D0,P%MTrix%D(1),1)


  ! Do forces
  DO AtA=1,GM%Natms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     CALL dNuke(GM,AtA,JFrc%D(A1:A2),LatFrc_J%D,NoWrap_O=NoWrap)
     !
     ! Start AtB Loop   
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1 
        AtB=P%ColPt%I(JP)
        B1=3*(AtB-1)+1
        B2=3*AtB
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN 
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
           MN=MN1+1
           A=Pair%A
           B=Pair%B
           !
           !Quick and dirty!
           SELECT CASE(P%NSMat)
           CASE(1)
              !We don't need to do anything!
           CASE(2)
              !We add up the two density martices!
              CALL DAXPY(MN,1D0,P%MTrix%D(Q+MN),1,P%MTrix%D(Q),1)
           CASE(4)
              !We add up the diagonal density martices!
              CALL DAXPY(MN,1D0,P%MTrix%D(Q+3*MN),1,P%MTrix%D(Q),1)
           CASE DEFAULT;CALL Halt(' JForce: P%NSMat doesn''t have an expected value! ')
           END SELECT
           !
           DO NC=1,GM%OvCells%NCells 
              Pair%A=A
              Pair%B=B+GM%OvCells%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                      +(Pair%A(2)-Pair%B(2))**2 &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair))THEN
!!                  WRITE(*,*)' C = ',NC
                  CALL TrPdJ2(Pair,P%MTrix%D(Q:Q+MN1),GM,JFrc%D(A1:A2),LatFrc_J%D,NoWrap_O=NoWrap)
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  ! Zero the lower triange of the lattice forces ...
  DO I=1,3
     DO J=1,I-1
        LatFrc_J%D(I,J)=Zero
     ENDDO
  ENDDO
!!$
!!$  PrintFlags%Key=DEBUG_MAXIMUM	
!!$  PrintFlags%MM=DEBUG_FRC
!!$  CALL Print_LatForce(GM,LatFrc_J%D,'TOTAL Lattice Force')
!!$  CALL Print_LatForce(GM,LatFrc_J%D,'TOTAL Lattice Force',Unit_O=6)

  ! ... and add them into the rest of the lattice gradients.
  GM%PBC%LatFrc%D=GM%PBC%LatFrc%D+LatFrc_J%D
  ! Add in the atomic forces as well 
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     GM%Gradients%D(1:3,AtA) = GM%Gradients%D(1:3,AtA)+JFrc%D(A1:A2)
  ENDDO
  !
!!  WRITE(*,11)' JFORCE Total Time = ',MTimer()-JFORCE_TotalTime_Start
  !
  K=0
  DO I=1,CS_IN%NCells
     !     WRITE(*,*)I,NNearCount(I)
     IF(NNearCount(I)==0D0)K=K+1
  ENDDO
!!$  WRITE(*,*)' % of NoPAC = ',DBLE(K)/DBLE(CS_IN%NCells)
!!$  WRITE(*,11)' Decompos_Time = ',Decompose_Time
!!$  WRITE(*,11)' TreeMake_Time = ',TreeMake_Time
!!$  WRITE(*,11)' JWalking_Time = ',JWalk_Time
!!$  WRITE(*,11)' Integral_Time = ',Integral_Time
!!$  WRITE(*,11)' Multipol_Time = ',Multipole_Time
!!$  WRITE(*,11)' Total J Time  = ',Decompose_Time+TreeMake_Time+JWalk_Time+Multipole_Time+Integral_Time
  !  WRITE(*,11)' Total JWalks  = ',DBLE(NPrim)
  !  WRITE(*,11)' Av  Ints/Prim = ',DBLE(NInts)/DBLE(NPrim)
  !  WRITE(*,11)' Av  # NF/Prim = ',DBLE(NNearAv)/DBLE(NPrim)
  !  WRITE(*,11)' Av  # FF/Prim = ',DBLE(NFarAv)/DBLE(NPrim)
  !  WRITE(*,11)' Time per INode= ',Integral_Time/DBLE(NNearAv)
  !  WRITE(*,11)' Time per MNode= ',Multipole_Time/DBLE(NFarAv)
11 FORMAT(A20,D12.6)

  ! Tidy Up
  DEALLOCATE(TempHerm%Coef)
!!$  ! Do some printing
!!$  CALL Print_Force(GM,JFrc,'J Force')
!!$  CALL Print_Force(GM,JFrc,'J Force',Unit_O=6)
  ! Do some checksumming and IO 
  CALL PChkSum(JFrc,    'dJ/dR',Proc_O=Prog)  
  CALL PChkSum(LatFrc_J,'LFrcJ',Proc_O=Prog)  
  ! Save Forces to Disk
#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL Put(GM, Tag_O = CurGeom)

    oldClone = MyClone
    DO rank = 1, MSize(MPI_COMM_WORLD)-1
      CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(GM, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

      ! Put to correct HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
      CALL Put(GM, Tag_O = CurGeom)
    ENDDO
    MyClone = oldClone

    ! Reopen old HDFGroup.
    CALL CloseHDFGroup(H5GroupID)
    H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID = H5GroupID
  ELSE
    CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(GM, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
  ENDIF
#else
  CALL Put(GM, Tag_O = CurGeom)
#endif
  ! Tidy Up  
  CALL Delete(JFrc)
  CALL Delete(LatFrc_J)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(RhoPoles)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  ! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 

  CALL MondoLog(DEBUG_MAXIMUM, "JForce", 'Total='//TRIM(DblToMedmChar(MTimer()-JFORCE_TotalTime_Start)) &
    //'; Bisect='//TRIM(DblToShrtChar(Decompose_Time))//', Tree='//TRIM(DblToShrtChar(TreeMake_Time)), "Timing")
  CALL MondoLog(DEBUG_MAXIMUM, "JForce", 'Walk='//TRIM(DblToShrtChar(JWalk_Time))//', Ints='//TRIM(DblToShrtChar(Integral_Time)) &
    //', Mults='//TRIM(DblToShrtChar(Multipole_Time)), "Timing")

  CALL ShutDown(Prog)
END PROGRAM JForce
