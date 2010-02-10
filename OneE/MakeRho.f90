! vim: tw=0
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
! COMPUTE THE DENSITY IN A HGTF BASIS FROM THE DENSITY MATRIX
! BASED ON AHMADI AND ALMLOF, CPL 246 p.364 (1995)
! Authors: Matt Challacombe and C.J. Tymczak and C. K. Gan
!----------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM MakeRho
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE BraBloks
  USE RhoBlok
  USE RhoTools
  USE MondoLogger

#ifdef PARALLEL
  USE MondoMPI
#endif

  IMPLICIT NONE

#ifdef PARALLEL
  TYPE(DBCSR)                     :: Dmat,D1,D2
  INTEGER                         :: LocalAtom,NumAtoms
  REAL(DOUBLE)                    :: TotRSumE,TotRSumE2,TotRSumN
  INTEGER                         :: IErr
  TYPE(CMPoles)                   :: SMP
#else
  TYPE(BCSR)                      :: Dmat,D1,D2
#endif
  INTEGER                         :: NC
  REAL(DOUBLE),DIMENSION(3)       :: B
  TYPE(AtomPair)                  :: Pair
  TYPE(BSET)                      :: BS
  TYPE(CRDS)                      :: GM
  TYPE(ARGMT)                     :: Args
  TYPE(HGRho)                     :: Rho
  TYPE(HGRho_new)                 :: RhoA
  TYPE(CMPoles)                   :: MP
  TYPE(INT_VECT)                  :: Stat
  TYPE(DBL_VECT)                  :: PTmp
  INTEGER                         :: P,R,AtA,AtB,NN
  INTEGER                         :: NDist,NCoef,Pbeg,Pend,NDist_old,NDist_new
  INTEGER                         :: PcntDist,OldFileID
  REAL(DOUBLE)                    :: RSumE,RSumE2,RSumN,dNel,RelRhoErr
  CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg1,Mssg2,Prog1,Prog2
  CHARACTER(LEN=*),PARAMETER      :: Prog='MakeRho'

#if defined(PARALLEL_CLONES)
  INTEGER                         :: oldClone, rank
#endif

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif

  ! Chose a density matrix
  IF(SCFActn=='BasisSetSwitch')THEN
    ! Get the previous information
    CALL Get(BS,PrvBase)
    CALL Get(GM,CurGeom)
    CALL SetThresholds(PrvBase)
    CALL Get(BSiz,'atsiz',PrvBase)
    CALL Get(OffS,'atoff',PrvBase)
    CALL Get(NBasF,'nbasf',PrvBase)
    CALL Get(Dmat,TrixFile('D',Args,-1))
  ELSEIF(SCFActn=='Restart'.OR. SCFActn=='RestartBasisSwitch')THEN
    ! Get the current geometry from the current HDF first
    CALL Get(GM,CurGeom)
    ! then close current group and HDF
    CALL CloseHDFGroup(H5GroupID)
    CALL CloseHDF(HDFFileID)
    ! and open the old group and HDF
    HDF_CurrentID=OpenHDF(Restart)
    OldFileID=HDF_CurrentID
    CALL New(Stat,3)
    CALL Get(Stat,'current_state')
    HDF_CurrentID=OpenHDFGroup(HDF_CurrentID,"Clone #"//TRIM(IntToChar(MyClone)))
    ! Get old basis set stuff
    SCFCycl=TRIM(IntToChar(Stat%I(1)))
    CurBase=TRIM(IntToChar(Stat%I(2)))
    CurGeom=TRIM(IntToChar(Stat%I(3)))
    CALL Get(BS,CurBase)
    ! Compute a sparse matrix blocking scheme for the old BS
    CALL BlockBuild(GM,BS,BSiz,OffS)
    NBasF=BS%NBasF
#ifdef PARALLEL
    CALL BCast(BSiz)
    CALL BCast(OffS)
    CALL BCast(NBasF)
#endif
    ! Close the old hdf up
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(OldFileID)
    ! Reopen current group and HDF
    HDFFileID=OpenHDF(H5File)
    H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID=H5GroupID
    CALL Get(Dmat,TrixFile('D',Args,0))
  ELSE
    ! Get the current information
    CALL Get(BS,CurBase)
    CALL Get(GM,CurGeom)
    SELECT CASE(SCFActn)

    CASE('InkFok')
      CALL Get(D1,TrixFile('D',Args,-1))
      CALL Get(D2,TrixFile('D',Args,0))
      CALL Multiply(D1,-One)
      CALL Add(D1,D2,DMat)
      IF(HasHF(ModelChem))THEN
        CALL Filter(D1,DMat)
        CALL Put(D1,TrixFile('DeltaD',Args,0))
      ENDIF
      CALL Delete(D1)
      CALL Delete(D2)

    CASE('ForceEvaluation')
      ! Index check...
      CALL MondoLog(DEBUG_MAXIMUM, TRIM(Prog), "Index Check: getting Dmat from "//TRIM(TrixFile('D',Args,1)))
      CALL Get(Dmat,TrixFile('D',Args,1))

    CASE('StartResponse')
      CALL Halt('MakeRho: SCFActn cannot be equal to <StartResponse>')

    CASE('DensityPrime')
      CALL Get(Dmat,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,0))

    CASE('Core')
      CALL Halt("I do not know what to do....")

    CASE DEFAULT
      ! Default
      CALL Get(Dmat,TrixFile('D',Args,0))

    END SELECT
  ENDIF
  CALL New(PTmp,4*MaxBlkSize**2)
  CALL NewBraBlok(BS)

  ! We are getting NSMat now from hdf via StartUp.
  !NSMat=DMat%NSMat

  !--------------------------------------------------------------
  ! Main loops: First pass calculates the size.
  ! Second pass calculates the density
  !-------------------------------------------------------------
  ! Initailize
  NDist = 0
  NCoef = 0
  ! Loop over atoms and count primatives
#ifdef PARALLEL
  DO LocalAtom = 1, Dmat%NAtms
    AtA  = Beg%I(MyID)+(LocalAtom-1)
    Pbeg = Dmat%RowPt%I(LocalAtom)
    Pend = Dmat%RowPt%I(LocalAtom+1)-1
#else
  DO AtA=1,NAtoms
    Pbeg = Dmat%RowPt%I(AtA)
    Pend = Dmat%RowPt%I(AtA+1)-1
#endif
    DO P = Pbeg,Pend
      AtB = Dmat%ColPt%I(P)
      IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
        B = Pair%B
        DO NC = 1,CS_OUT%NCells
          Pair%B = B+CS_OUT%CellCarts%D(:,NC)
          Pair%AB2 = (Pair%A(1)-Pair%B(1))**2 &
               + (Pair%A(2)-Pair%B(2))**2 &
               + (Pair%A(3)-Pair%B(3))**2
          IF(TestAtomPair(Pair)) THEN
            CALL PrimCount(BS,Pair,NDist,NCoef)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  ! Allocate the Density
  SELECT CASE(NSMat)
  CASE(1)
    CALL New_HGRho_new(RhoA,(/NDist,NCoef,1/))
  CASE(2)
    CALL New_HGRho_new(RhoA,(/NDist,NCoef,3/))!<<< SPIN 3-> rho_tot,rho_a,rho_b
  CASE(4)
    CALL New_HGRho_new(RhoA,(/NDist,NCoef,4/))!<<< SPIN 4-> rho_tot,rho_a,rho_b,rho_ab
  CASE DEFAULT
    CALL Halt('MakeRho: NSMat not valid!')
  END SELECT

  ! Initialize  Counters
  NDist = 0
  NCoef = 0
#ifdef PARALLEL
  DO LocalAtom = 1, Dmat%NAtms
    AtA  = Beg%I(MyID)+(LocalAtom-1)
    Pbeg = Dmat%RowPt%I(LocalAtom)
    Pend = Dmat%RowPt%I(LocalAtom+1)-1
#else
  DO AtA=1,NAtoms
    Pbeg = Dmat%RowPt%I(AtA)
    Pend = Dmat%RowPt%I(AtA+1)-1
#endif
    DO P=Pbeg,Pend
      AtB = Dmat%ColPt%I(P)
      R   = Dmat%BlkPt%I(P)
      NN=BSiz%I(AtA)*BSiz%I(AtB)
      ! Quick and dirty.
      SELECT CASE(NSMat)
      CASE(1)
        !We need to copy the matrix!
        CALL DCOPY(NN,Dmat%MTrix%D(R),1,PTmp%D(1),1)
      CASE(2)
        !We copy the first matrix!
        CALL DCOPY(NN,Dmat%MTrix%D(R),1,PTmp%D(1),1)
        !(Pa+Pb)/2->Ptot
        CALL DAXPY(NN,1D0,Dmat%MTrix%D(R+NN),1,PTmp%D(1),1)
        ! Scale the total density matrix.
        CALL DSCAL(NN,0.5D0,PTmp%D(1),1)
        !Pa
        CALL DCOPY(NN,Dmat%MTrix%D(R),1,PTmp%D(NN+1),1)
        !Pb
        CALL DCOPY(NN,Dmat%MTrix%D(R+NN),1,PTmp%D(2*NN+1),1)
      CASE(4)
        !We copy the first matrix!
        CALL DCOPY(NN,Dmat%MTrix%D(R),1,PTmp%D(1),1)
        !(Pa+Pb)/2->Ptot
        CALL DAXPY(NN,1D0,Dmat%MTrix%D(R+3*NN),1,PTmp%D(1),1)
        ! Scale the total density matrix.
        CALL DSCAL(NN,0.5D0,PTmp%D(1),1)
        ! TODO something here for HiCu densities, i.e rho_a,rho_b,rho_ab
        !...
      CASE DEFAULT
        CALL Halt('[MakeRho] NSMat doesn''t have an expected value! ')
      END SELECT

      IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
        B = Pair%B
        DO NC = 1,CS_OUT%NCells
          Pair%B = B+CS_OUT%CellCarts%D(:,NC)
          Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
               + (Pair%A(2)-Pair%B(2))**2 &
               + (Pair%A(3)-Pair%B(3))**2
          IF(TestAtomPair(Pair)) THEN
            NN = Pair%NA*Pair%NB
            !CALL RhoBlk(BS,Dmat%MTrix%D(R),Pair,NDist,NCoef,RhoA)
            CALL RhoBlk(BS,PTmp%D(1),Pair,NDist,NCoef,RhoA)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDDO

  ! Don't add in nuclear charges if incremental fock builds or CPSCF.
  IF(SCFActn/='InkFok'.AND. SCFActn/='StartResponse'.AND.SCFActn/='DensityPrime') THEN
#ifdef PARALLEL
    CALL AddDist(RhoA,GM,NuclearExpnt,Beg%I(MyID),End%I(MyID))
#else
    CALL AddDist(RhoA,GM,NuclearExpnt,1,GM%NAtms)
#endif
  ENDIF
  ! Prune negligible distributions from the electronic density
  NDist_old = RhoA%NDist
  CALL Prune_Rho_new(Thresholds%Dist,RhoA)
  NDist_new = RhoA%NDist
  !************************************
!!$  EllPrune=4
!!$  CALL PruneEll_Rho_new(EllPrune,RhoA)
!!$  WRITE(*,*) 'Ell = ',EllPrune
!!$  WRITE(*,*) 'Number of Dists = ', RhoA%NDist,NDist_new
  !************************************
  !
  ! Fold distributions back into the box; For ForceEvaluation, rho is not folded
  ! This is turned off for now inorder to make the lattice force calculation
  ! correct
  ! IF(SCFActn .NE. 'ForceEvaluation')
  ! CALL Fold_Rho_new(GM,RhoA)
  !
  ! Compute integrated electron and nuclear densities
#ifdef PARALLEL
  NumAtoms = End%I(MyID)-Beg%I(MyID)+1
  RSumE    = Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-NumAtoms)
  RSumN    = Integrate_HGRho_new(RhoA,1,RhoA%NDist-NumAtoms+1,RhoA%NDist)
  TotRSumE = AllReduce(RSumE)
  TotRSumN = AllReduce(RSumN)
  RSumE    = TotRSumE
  RSumN    = TotRSumN
  IF(NSMat.EQ.1)THEN
    RSumE2=Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-NumAtoms)
    TotRSumE2 = Reduce(RSumE2)
    IF(MyID.EQ.0) THEN
      CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_tot = '//TRIM(DblToChar(TotRSumE2)))
    ENDIF
  ELSEIF(NSMat.EQ.2)THEN
    RSumE2=Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-NumAtoms)
    TotRSumE2 = Reduce(RSumE2)
    IF(MyID.EQ.0) CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_tot = '//TRIM(DblToChar(TotRSumE2)))
    RSumE2=Integrate_HGRho_new(RhoA,2,1,RhoA%NDist-NumAtoms)
    TotRSumE2 = Reduce(RSumE2)
    IF(MyID.EQ.0) CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_alp = '//TRIM(DblToChar(TotRSumE2)))
    RSumE2=Integrate_HGRho_new(RhoA,3,1,RhoA%NDist-NumAtoms)
    TotRSumE2 = Reduce(RSumE2)
    IF(MyID.EQ.0) CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_bet = '//TRIM(DblToChar(TotRSumE2)))
  ELSEIF(NSMat.EQ.4)THEN
    RSumE2=Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-NumAtoms)
    IF(MyID.EQ.0) CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_tot = '//TRIM(DblToChar(TotRSumE2)))
  ENDIF
  IF(MyID.EQ.0) CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_nuc = '//TRIM(DblToChar(RSumN)))
#else
  RSumE  =  Integrate_HGRho_new(RhoA,1,1                    ,RhoA%NDist-GM%NAtms)
  RSumN  =  Integrate_HGRho_new(RhoA,1,RhoA%NDist-GM%NAtms+1,RhoA%NDist         )
  IF(NSMat.EQ.1)THEN
    RSumE2=Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-GM%NAtms)
    CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_tot = '//TRIM(DblToChar(RSumE2)))
  ELSEIF(NSMat.EQ.2)THEN
    RSumE2=Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-GM%NAtms)
    CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_tot = '//TRIM(DblToChar(RSumE2)))
    RSumE2=Integrate_HGRho_new(RhoA,2,1,RhoA%NDist-GM%NAtms)
    CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_alp = '//TRIM(DblToChar(RSumE2)))
    RSumE2=Integrate_HGRho_new(RhoA,3,1,RhoA%NDist-GM%NAtms)
    CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_bet = '//TRIM(DblToChar(RSumE2)))
  ELSEIF(NSMat.EQ.4)THEN
    RSumE2=Integrate_HGRho_new(RhoA,1,1,RhoA%NDist-GM%NAtms)
    CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_tot = '//TRIM(DblToChar(RSumE2)))
  ENDIF
  CALL MondoLog(DEBUG_MAXIMUM, "MakeRho", 'Rho_nuc = '//TRIM(DblToChar(RSumN)))
#endif
  ! Calculate dipole and quadrupole moments
  CALL New(MP)
  CALL CalRhoPoles_new(MP,RhoA,GM)
#ifdef PARALLEL
  CALL New(SMP)
  CALL MPI_Reduce(MP%Dpole%D(1),SMP%Dpole%D(1),3,MPI_DOUBLE_PRECISION,MPI_SUM,ROOT,MONDO_COMM,IErr)
  CALL MPI_Reduce(MP%Qpole%D(1),SMP%Qpole%D(1),6,MPI_DOUBLE_PRECISION,MPI_SUM,ROOT,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    MP%Dpole%D(1:3) = SMP%Dpole%D(1:3)
    MP%Qpole%D(1:6) = SMP%Qpole%D(1:6)
  ENDIF
#endif
  ! Convert to the old format
  CALL ConvertToOldRho(Rho,RhoA)
  ! Do Some Outputing
  IF(SCFActn=='InkFok')THEN
    Prog1='InkFok'
    Prog2=Prog1
  ELSE
    Prog1='Pruned Rho'
    Prog2='Moments'
  ENDIF
  dNel     = Two*(RSumE+RSumN)+TotCh
  RelRhoErr= ABS(dNel)/DBLE(NEl)
  PcntDist=FLOOR(1.D2*DBLE(NDist_new)/DBLE(NDist_old))
  Mssg1='dNel = '//TRIM(DblToShrtChar(dNel))//', kept '  &
       //TRIM(IntToChar(PcntDist))//'% of distributions.'
  Mssg2='<r> = ('//TRIM(DblToShrtChar(MP%DPole%D(1))) &
       //', '//TRIM(DblToShrtChar(MP%DPole%D(2)))       &
       //', '//TRIM(DblToShrtChar(MP%DPole%D(3)))       &
       //'), <r^2> = '//TRIM(DblToShrtChar(             &
       MP%QPole%D(1)+MP%QPole%D(2)+MP%QPole%D(3)))
  ! Output pruning and multipole stats
  IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
#ifdef PARALLEL
    IF(MyID == ROOT) THEN
#endif
      CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg1, Prog1)
      CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg2, Prog1)
#ifdef PARALLEL
    ENDIF
#endif
  ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
#if defined(PARALLEL) || defined(PARALLEL_CLONES)
    IF(MyID == ROOT) THEN
#endif
      CALL MondoLog(DEBUG_MAXIMUM, Prog, Mssg1, Prog1)
#if defined(PARALLEL) || defined(PARALLEL_CLONES)
    ENDIF
#endif
  ENDIF
  ! Check error
  IF(RelRhoErr>Thresholds%Dist*5.D3)THEN
#ifdef PARALLEL
    IF(MyID == ROOT) THEN
#endif
      CALL Warn(ProcessName(Prog)//'relative error in density = '//TRIM(DblToShrtChar(RelRhoErr)) &
           //'. Distribution threshold = '//TRIM(DblToShrtChar(Thresholds%Dist))      &
           //'. Total charge lost = '//TRIM(DblToShrtChar(dNel)))
#ifdef PARALLEL
    ENDIF
#endif
  ENDIF
  !------------------------------------------------------------------------------------
  ! Put Rho and MPs to disk
  SELECT CASE(SCFActn)

  CASE('ForceEvaluation')
#if defined(PARALLEL)
    CALL Put(Rho, 'Rho'//IntToChar(MyID), Args, 1)
#elif defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
      CALL Put(Rho, "Rho", Args, 1)
      CALL Put(MP)

      oldClone = MyClone
      DO rank = 1, MSize(MPI_COMM_WORLD)-1
        CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(Rho, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(MP, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

        ! Put to correct HDFGroup.
        CALL CloseHDFGroup(H5GroupID)
        H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID = H5GroupID
        CALL Put(Rho, "Rho", Args, 1)
        CALL Put(MP)
      ENDDO
      MyClone = oldClone

      ! Reopen old HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
    ELSE
      !CALL MondoLog(DEBUG_NONE, Prog, "sending density and multipoles to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
      CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(Rho, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(MP, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    ENDIF
#else
    CALL Put(Rho, "Rho", Args, 1)
    CALL Put(MP)
#endif

  CASE('InkFok')
#ifdef PARALLEL
    CALL Halt("[FIXME]")
#elif defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
      CALL MondoLog(DEBUG_NONE, Prog, "writing density and multipoles to hdf", "Clone "//TRIM(IntToChar(MyClone)))
      CALL Halt("[FIXME]")
      CALL Put(Rho,'DeltaRho',Args,0)
      CALL Put(MP,'Delta'//TRIM(SCFCycl))
    ELSE
      CALL MondoLog(DEBUG_NONE, Prog, "sending density and multipoles to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
    ENDIF
#else
    CALL Put(Rho,'DeltaRho',Args,0)
    CALL Put(MP,'Delta'//TRIM(SCFCycl))
#endif

  CASE DEFAULT
#if defined(PARALLEL)
    CALL Put(Rho, "Rho"//IntToChar(MyID), Args, 0)
#elif defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
      CALL Put(Rho, "Rho", Args, 0)
      CALL Put(MP)

      oldClone = MyClone
      DO rank = 1, MSize(MPI_COMM_WORLD)-1
        CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(Rho, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(MP, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

        ! Put to correct HDFGroup.
        CALL CloseHDFGroup(H5GroupID)
        H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID = H5GroupID
        CALL Put(Rho, "Rho", Args, 0)
        CALL Put(MP)
      ENDDO
      MyClone = oldClone

      ! Reopen old HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
    ELSE
      !CALL MondoLog(DEBUG_MAXIMUM, Prog, "sending density and multipoles to clone 1", "Clone "//TRIM(IntToChar(MyClone)))
      CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(Rho, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(MP, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    ENDIF
#else
    CALL Put(Rho, "Rho", Args, 0)
    CALL Put(MP)
#endif

  END SELECT

  CALL PChkSum(Rho,'Rho',.TRUE.,Prog)

  ! Tidy up
  CALL Delete(GM)
  CALL Delete(Dmat)
  CALL Delete(BS)
  CALL Delete(PTmp)
  CALL DeleteBraBlok()
  CALL Delete(Rho)
  CALL Delete_HGRho_new(RhoA)
  CALL ShutDown(Prog)

END PROGRAM MakeRho
