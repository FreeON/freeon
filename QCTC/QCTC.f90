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
!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!    Authors:  Matt Challacombe and CJ Tymczak
!===============================================================================

#include "MondoConfig.h"

PROGRAM QCTC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE BraBloks
  USE QCTCThresholds
  USE PoleTree
  USE Globals
  USE PBCFarField
  USE JGen
  USE NukularE
  USE Density
  USE Clock
  USE TreeWalk

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  USE MondoMPI
#endif

  IMPLICIT NONE

  TYPE(BCSR)                     :: J,  P
  TYPE(BCSR)                     :: T1,T2
  TYPE(BCSR)                     :: Dmat,D1,D2
  TYPE(INT_VECT)                 :: Stat
  TYPE(HGLL),POINTER             :: RhoHead
  REAL(DOUBLE)                   :: E_Nuc_Tot,Etot,SdvErrorJ,MaxErrorJ,JMExact,E_PFF,D_PFF

  CHARACTER(LEN=4),PARAMETER     :: Prog='QCTC'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  INTEGER                        :: I,K,NLink,OldFileID, oldMyClone, IErr

  REAL(DOUBLE)                   :: QCTC_TotalTime_Start

  TYPE(TIME)                     :: TimeMakeJ,TimeMakeTree,TimeNukE

  LOGICAL                        :: NoWrap=.FALSE. ! WRAPPING IS ON

#if defined(PARALLEL_CLONES)
  INTEGER :: oldClone, rank
#endif

  ! Start timer.
  QCTC_TotalTime_Start=MTimer()

  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
  NukesOn=.TRUE.
  ElectronsOn=.TRUE.

  ! ---------------------------------------------------------------------------------
  ! Begin building a density
  !
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
     NukesOn=.TRUE.
     ElectronsOn=.TRUE.
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
     ! Close the old hdf up
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     HDF_CurrentID=HDFFileID
     !
     CALL Get(Stat,'current_state')
     SCFCycl=TRIM(IntToChar(Stat%I(1)))
     CurBase=TRIM(IntToChar(Stat%I(2)))
     CurGeom=TRIM(IntToChar(Stat%I(3)))
     !
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
     !
     CALL Get(Dmat,TrixFile('D',Args,0))
     !
     NukesOn=.TRUE.
     ElectronsOn=.TRUE.
  ELSE
     ! Get the current information
     CALL Get(BS,CurBase)
     CALL Get(GM,CurGeom)
     IF(SCFActn=='InkFok')THEN ! Incremental Fock build
        ! Maybe the difference density build should go elswhere?
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
        !
        NukesOn=.FALSE.
        ElectronsOn=.TRUE.
     ELSEIF(SCFActn=='FockPrimeBuild')THEN
        CALL Get(Dmat,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,0))
        NukesOn=.FALSE.
        ElectronsOn=.TRUE.
     ELSEIF(SCFActn=='TD-SCF')THEN
        CALL Get(Dmat,TrixFile(Args%C%C(3),Args,0))
        NukesOn=.FALSE.
        ElectronsOn=.TRUE.
     ELSEIF(SCFActn=='GuessEqCore')THEN
        ! No electrons, just the nukes
        NukesOn=.TRUE.
        ElectronsOn=.FALSE.
     ELSE
        ! Default
        CALL Get(Dmat,TrixFile('D',Args,0))
        NukesOn=.TRUE.
        ElectronsOn=.TRUE.
     ENDIF
  ENDIF
  ! Set thresholds local to the density build
  CALL SetLocalThresholds(Thresholds%TwoE)
  ! -------------------------------------------------------------------------
  ! Build the density
  Density_Time_Start=MTimer()
  ! -------------------------------------------------------------------------
  ! RhoHead is the start of a linked density list
  ALLOCATE(RhoHead)
  RhoHead%LNum=0
  ! Here, the LL is filled out
  IF(ElectronsOn)THEN
     CALL MakeRhoList(GM,BS,DMat,NLink,RhoHead,'QCTC',NoWrap_O=NoWrap)
  ENDIF
  ! Add in the nuclear charges only in certain cases
  IF(NukesOn)THEN
     CALL AddNukes(GM,RhoHead,NoWrap)
     NLink=NLink+GM%NAtms
  ENDIF
  ! Load density into arrays and delete the linked list
  CALL Collate(GM,RhoHead,Rho,'QCTC',RhoPoles,NLink)
  ! Delete is broken on some compilers, avoid for now
  !  CALL DeleteHGLL(RhoHead)
  ! -------------------------------------------------------------------------
  Density_Time=MTimer()-Density_Time_Start
  ! Done building density
  ! -------------------------------------------------------------------------
  ! Take care of some micilaneous book keeping
  ! -------------------------------------------------------------------------
  ! Set some values of the angular symmetry for use in computing the Lorentz field
  MaxPFFFEll=GM%PBC%PFFMaxEll
  ! Local expansion order of the multipoles to use in the tree
  MaxPoleEll=MIN(2*BS%NASym+6,MaxPFFFEll)
  IF(MaxPoleEll<2*(BS%NASym+1)) &
     CALL Halt('Bombed in QCTC. Please set PFFMaxEll larger ')
  ! Find the total energy from past calculations
  ! Over-ride initialization to big_dbl in punchhdf, so that
  ! PFF is computed correctly in PBCFarFieldSetUp
  IF(SCFCycl=='0'.AND.CurGeom=='1')THEN
     ETot=1D2
  ELSE
    CALL Get(Etot,'Etot')
  ENDIF

  ! Now that we are done with the density, lets make sure we have the
  ! current basis set, matrix block sizes etc:
  CALL Get(BS,CurBase)
  CALL Get(BSiz,'atsiz',CurBase)
  CALL Get(OffS,'atoff',CurBase)
  CALL Get(NBasF,'nbasf',CurBase)

  ! Set space for bra blocking
  CALL NewBraBlok(BS)

  ! Thresholds local to J matrix build (may be different from those used to build density)
  CALL SetThresholds(CurBase)

  CALL SetLocalThresholds(Thresholds%TwoE)
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
  ! Here we set the max multipole expansion relative
  ! to the max angular symmetry
  NInts=0
  NPrim=0
  NFarAv=0
  NNearAv=0
  ! -------------------------------------------------------------------------
  ! Done with book keeping.
  ! -------------------------------------------------------------------------
  ! Now build the multipole tree
  TreeMake_Time_Start=MTimer()
  ! -------------------------------------------------------------------------
  ! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
  ! Set up the crystal field
  CALL PBCFarFieldSetUp(GM,Rho,'QCTC',MaxPFFFEll,Etot)
  ! Delete the auxiliary density arrays
  CALL DeleteRhoAux
  ! Delete the Density
  !  CALL Delete(Rho)
  ! -------------------------------------------------------------------------
  TreeMake_Time=MTimer()-TreeMake_Time_Start
  ! Done making the pole-tree
  ! -------------------------------------------------------------------------
  ! Compute the Coulomb matrix J in O(N Lg N)
  CALL Elapsed_Time(TimeMakeJ,'Init')
  JWalk_Time=0D0
  Integral_Time=0D0
  Multipole_Time=0D0
  ALLOCATE(NNearCount(1:CS_IN%NCells))
  NNearCount=0D0
  ! -------------------------------------------------------------------------
  CALL New(J)
  CALL MakeJ(J,NoWrap_O=NoWrap) ! <<<<<<<<<
  ! -------------------------------------------------------------------------
  CALL Elapsed_TIME(TimeMakeJ,'Accum')
  ! Done making Coulomb matrix
  ! -------------------------------------------------------------------------
  IF(SCFActn=='InkFok')THEN
     ! Add in correction if incremental J build
     CALL New(T1)
     CALL New(T2)
     CALL Get(T1,TrixFile('J',Args,-1))
     CALL Add(T1,J,T2)
     CALL Filter(T1,T2)
     CALL Delete(T2)
  ELSE
     CALL Filter(T1,J)
  ENDIF
  ! Put J to disk
  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
     CALL Put(T1,TrixFile('JPrime'//TRIM(Args%C%C(3)),Args,0))
  ELSE
     CALL Put(T1,TrixFile('J',Args,0))
  ENDIF
  ! Nuclear-total density Coulomb interaction
  IF(NukesOn)THEN
     CALL Elapsed_Time(TimeNukE,'Init')
     E_Nuc_Tot=NukE(GM,Rho,NoWrap)
     CALL Elapsed_Time(TimeNukE,'Accum')
  ELSEIF(SCFActn=='InkFok')THEN
     CALL Get(E_Nuc_Tot,'E_NuclearTotal',Stats_O=Previous)
     CALL Elapsed_Time(TimeNukE,'Init')
     E_Nuc_Tot=E_Nuc_Tot+NukE(GM,Rho,NoWrap)
     CALL Elapsed_Time(TimeNukE,'Accum')
  ELSE
     E_Nuc_Tot=Zero
  ENDIF

  CALL MondoLog(DEBUG_MAXIMUM, Prog, "E_NuclearTotal = "//TRIM(DblToChar(E_Nuc_Tot))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))

  PerfMon%FLOP=Zero
  CALL MondoLog(DEBUG_MAXIMUM,Prog,'Coulomb Energy = '//TRIM(DblToChar(E_Nuc_Tot+Trace(DMat,T1)))//" hartree")

  CALL MondoLog(DEBUG_MAXIMUM,Prog,'CPUSec='//TRIM(DblToMedmChar(MTimer()-QCTC_TotalTime_Start)) &
                               //'; RhoBld='//TRIM(DblToShrtChar(Density_Time))                  &
                               //', TreeBld='//TRIM(DblToShrtChar(TreeMake_Time)))
  CALL MondoLog(DEBUG_MAXIMUM,Prog,'Walk='//TRIM(DblToShrtChar(JWalk_Time))                     &
                                //', Ints='//TRIM(DblToShrtChar(Integral_Time))                  &
                                //', Mults='//TRIM(DblToShrtChar(Multipole_Time)))

  IF(CS_IN%NCells>1)THEN
     CALL MondoLog(DEBUG_MAXIMUM,Prog, &
                   '% Ints in UC= '//TRIM(FltToShrtChar(1D2*NNearCount(CS_IN%NCells)/(NLeaf*NPrim))) &
                 //', % Total Ints = '//TRIM(FltToShrtChar(                                &
                     1D2*SUM(NNearCount(1:CS_IN%NCells))/(NLeaf*NPrim*CS_IN%NCells) )))
  ELSE
     CALL MondoLog(DEBUG_MAXIMUM,Prog, &
                   '% NFInts = '//TRIM(FltToShrtChar(1D2*NNearCount(1)/(NLeaf*NPrim))))
  ENDIF

  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
     CALL PChkSum(T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     CALL Plot(   T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
  ELSE
     CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')!,Unit_O=6)
     CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
  ENDIF

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
     CALL Put(E_Nuc_Tot, 'E_NuclearTotal', Stats_O = Current)

     oldClone = MyClone
     DO rank = 1, MSize(MPI_COMM_WORLD)-1
        CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(E_Nuc_Tot, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

        ! Put to correct HDFGroup.
        CALL CloseHDFGroup(H5GroupID)
        H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID = H5GroupID
        CALL Put(E_Nuc_Tot, 'E_NuclearTotal', Stats_O = Current)
     ENDDO
     MyClone = oldClone

     ! Reopen old HDFGroup.
     CALL CloseHDFGroup(H5GroupID)
     H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID = H5GroupID
  ELSE
     CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
     CALL Send(E_Nuc_Tot, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
  ENDIF
#else
  CALL Put(E_Nuc_Tot, 'E_NuclearTotal', Stats_O = Current)
#endif

  !-------------------------------------------------------------------------------
  ! All done
  !-------------------------------------------------------------------------------
  CALL ShutDown(Prog)

END PROGRAM QCTC
