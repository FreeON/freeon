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
!    COMPUTES AN ORTHOGONAL GUESS DENSITY MATRIX EITHER FROM
!    A SUPERPOSITION OF DIAGONAL, ATOMIC LEWIS STRUCTURE BLOCKS
!    OR FROM A ORTHOGONAL DENSITY MATRIX COORESPONDING TO A DIFFERENT
!    (HOPEFULLY CLOSE) GEOMETRY OR LEVEL OF ACCURACY
!    Author: Matt Challacombe
!-------------------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM P2Use
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE PBlokGuess
  USE MondoLogger
#ifdef PARALLEL
  USE MondoMPI
#endif
  USE DenMatMethods

  IMPLICIT NONE

  TYPE(ARGMT)                   :: Args
  TYPE(BSET)                    :: BS,BS_old
  TYPE(CRDS)                    :: GM,GM_old
#ifdef PARALLEL
  TYPE(BCSR)  :: P_BCSR,F_BCSR
  TYPE(DBCSR) :: F,P,P0,PV,X,S,S0,S1,Tmp1,Tmp2,T1,T2
#else
  TYPE(BCSR) :: F,P,P0,PV,X,S,S0,S1,Tmp1,Tmp2,T1,T2
#endif

  TYPE(INT_VECT)                :: Stat
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: MaxDS,NoiseLevel, alpha, beta, v_scale, MDDeltaTime,A1,B1
  INTEGER                       :: MDDampStep, m_step, mm_step, Imin, Imin_divergence_check
  REAL(DOUBLE)                  :: Scale,Fact,ECount,RelNErr, DeltaP,OldDeltaP, &
       DensityDev,dN,MaxGDIff,GDIff,OldN,M,PNon0s,PSMin,PSMax, &
       Ipot_Error,Ipot_Error_start,Norm_Error,Lam,DLam,TError0,SFac,Dum,Fmin,Fmax,Occ3,Occ2,Occ1,Occ0
  INTEGER                       :: I,J,JP,AtA,Q,R,T,KA,NBFA,NPur,PcntPNon0, &
       OldFileID,ICart,N,NStep,iGEO,iBAS,iSCF,DMPOrder,MM,ICycle,Cycle
  CHARACTER(LEN=2)              :: Cycl
  LOGICAL                       :: Present,DoingMD,AOSPExit
  CHARACTER(LEN=DEFAULT_CHR_LEN):: Mssg,BName,FileName,DMFile,logtag
  CHARACTER(LEN=8)              :: MDGeuss
  CHARACTER(LEN=5),PARAMETER    :: Prog='P2Use'

  !-------------------------------------------------------------------------------
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  Cycl=IntToChar(Args%I%I(1))

  iSCF = Args%I%I(1)
  iBAS = Args%I%I(2)
  iGEO = Args%I%I(3)

  ! Select spin factor for R/U/G theory.
  SFac = 2D0
  IF(NSMat > 1) SFac = 1D0
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "CurBase = "//TRIM(IntToChar(Args%I%I(2)))// &
    ", NSMat = "//TRIM(IntToChar(NSMat)))

  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)

  ! Start logging.
  logtag = TRIM(Prog)//":"//TRIM(SCFActn)
  ! CALL MondoLog(DEBUG_NONE, "P2Use", "SCFActn = "//TRIM(SCFActn))
  ! Do what needs to be done CASE by CASE

  SELECT CASE(SCFActn)

    ! P=0
  CASE('GuessEqCore')

    CALL New(P,NSMat_O=NSMat)
    CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
    DO I=1,NAtoms
      BlkP%D(:,I)=Zero
    ENDDO
    CALL SetToI(P,BlkP)
    CALL Delete(BlkP)
    ! IO for the non-orthogonal P
    CALL Put(P,TrixFile('D',Args,0))
    CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
    CALL PPrint( P,'P['//TRIM(Cycl)//']')
    CALL Plot(   P,'P_'//TRIM(Cycl))
    CALL Delete(P)

    ! Density SuperPosition
  CASE('DensitySuperposition','DMDGeuss')
    CALL New(P,NSMat_O=NSMat)
    CALL New(Tmp1)
    CALL New(Tmp2)
    CALL Get(BName,'bsetname',CurBase)

    IF(INDEX(BName,'STO')/=0)THEN
      ! Compute a diagonal guess as the superposition of atomic lewis structure
      ! occupancies--works only for minimal (STO) basis sets
      CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
      DO I=1,NAtoms
        IF(GM%AtNum%D(I) < 105.D0) THEN
          CALL FillPBlok(BSiz%I(I),INT(GM%AtNum%D(I)),BlkP%D(:,I))
          ! Rescale for charged molecules.
          IF(TotCh.NE.Zero)BlkP%D(:,I)=NEl*BlkP%D(:,I)/(NEl+TotCh)
        ENDIF
      ENDDO

      IF(NSMat.EQ.1)THEN
        !Set the P
        CALL SetToI(P,BlkP)
      ELSEIF(NSMat.EQ.2)THEN
         ! Check for the case of an uncharged singlet
         IF(NAlph==NBeta)THEN
            ! We have a perfectly paired singlet.  Alternate the spin density to avoid the restricted solution.
            BlkP%D=Zero
            DO I=1,NAtoms
               IF(MOD(I,2)==0) &
                    CALL FillPBlok(BSiz%I(I),INT(GM%AtNum%D(I)),BlkP%D(:,I))
            ENDDO
            A1=SUM(BlkP%D)
            BlkP%D=BlkP%D*(Dble(NAlph)/A1)
            CALL SetToI(P,BlkP,Expert_O=1)
            BlkP%D=Zero
            DO I=1,NAtoms
               IF(MOD(I,2)==1) &
                    CALL FillPBlok(BSiz%I(I),INT(GM%AtNum%D(I)),BlkP%D(:,I))

            ENDDO
            B1=SUM(BlkP%D)
            BlkP%D=BlkP%D*(Dble(NBeta)/B1)
            CALL SetToI(P,BlkP,Expert_O=2)
         ELSE
            !Set the P_alpha
            CALL DSCAL(MaxBlkSize**2*NAtoms,2D0*DBLE(NAlph)/DBLE(NEl),BlkP%D(1,1),1)
            CALL SetToI(P,BlkP,Expert_O=1)
            !Set the P_beta
            CALL DSCAL(MaxBlkSize**2*NAtoms,DBLE(NBeta)/DBLE(NAlph),BlkP%D(1,1),1)
            CALL SetToI(P,BlkP,Expert_O=2)
         ENDIF
      ELSEIF(NSMat.EQ.4)THEN
        !Set the P_alpha
        CALL DSCAL(MaxBlkSize**2*NAtoms,2D0*DBLE(NAlph)/DBLE(NEl)  ,BlkP%D(1,1),1)
        CALL SetToI(P,BlkP,Expert_O=1)

        !Set the off diag bloks P_alpha,beta P_beta,alpha
        Dum=-0.1D0
        CALL DSCAL(MaxBlkSize**2*NAtoms,Dum/DBLE(NAlph),BlkP%D(1,1),1)
        CALL SetToI(P,BlkP,Expert_O=2)
        CALL SetToI(P,BlkP,Expert_O=3)

        !Set the P_beta
        CALL DSCAL(MaxBlkSize**2*NAtoms,DBLE(NBeta)/Dum,BlkP%D(1,1),1)
        CALL SetToI(P,BlkP,Expert_O=4)
        ! BlkP%D=0d0
        ! CALL SetToI(P,BlkP,Expert_O=2)
        ! CALL SetToI(P,BlkP,Expert_O=3)
      ELSE
        CALL Halt("[P2Use:"//TRIM(SCFActn)//"] Something wrong when computing P_guess!")
      ENDIF

      ! Check for the correct electron count
      TrP=Trace(P)
      IF(ABS(TrP-DBLE(NEl/SFac))>1.D-10) THEN
        CALL Warn("[P2Use:"//TRIM(SCFActn)//"] TrP = "//TRIM(DblToChar(TrP)))
      ENDIF
      CALL Delete(BlkP)
    ELSE
      CALL SetToI(P)
#ifdef PARALLEL
      IF(MyID == ROOT) &
#endif
      CALL Warn("[P2Use:"//TRIM(SCFActn)//"] Attempting to use density superposition with a non STO basis set. Going for scaled I.")

      IF(NSMat.EQ.1)THEN
        !Set the P
        CALL Multiply(P,DBLE(NEl)/(Two*DBLE(NBasF)))
      ELSEIF(NSMat.EQ.2)THEN
         IF(NAlph==NBeta)THEN
            CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
            ! We have a perfectly paired singlet.  Alternate the spin density to avoid the restricted solution.
            A1=Zero
            BlkP%D=Zero
            DO I=1,NAtoms
               IF(MOD(I,2)==0)THEN
                  BlkP%D(:,I)=One
                  A1=A1+BS%BFKnd%I(GM%AtTyp%I(I))
               ENDIF
            ENDDO

            BlkP%D=BlkP%D*(Dble(NAlph)/A1)

!            WRITE(*,*)' Alpha, Beta = ',SUM(BlkP%D)

            CALL SetToI(P,BlkP,Expert_O=1)
!
            B1=Zero
            BlkP%D=Zero
            DO I=1,NAtoms
               IF(MOD(I,2)==1)THEN
                  BlkP%D(:,I)=One
                  B1=B1+BS%BFKnd%I(GM%AtTyp%I(I))
               ENDIF
            ENDDO

            BlkP%D=BlkP%D*(Dble(NBeta)/B1)
            CALL SetToI(P,BlkP,Expert_O=2)

            WRITE(*,*)' Alpha, Beta = ',SUM(BlkP%D)


            CALL Delete(BlkP)
         ELSE
            !Set the P_alpha
            CALL Multiply(P,DBLE(NAlph)/DBLE(NBasF),Expert_O=1)
            !Set the P_beta
            CALL Multiply(P,DBLE(NBeta)/DBLE(NBasF),Expert_O=2)
         ENDIF
      ELSEIF(NSMat.EQ.4)THEN
        !Set the P_alpha
        CALL Multiply(P,DBLE(NAlph)/DBLE(NBasF),Expert_O=1)
        !Set the P_beta
        CALL Multiply(P,DBLE(NBeta)/DBLE(NBasF),Expert_O=4)
        !Set the P_ab and P_ba
        CALL Multiply(P,DBLE(NEl)/DBLE(NBasF),Expert_O=2)
        CALL Multiply(P,DBLE(NEl)/DBLE(NBasF),Expert_O=3)
      ELSE
        CALL Halt("[P2Use:"//TRIM(SCFActn)//"] Something wrong when computing P_guess!")
      ENDIF

      TrP=Trace(P)
      IF(ABS(TrP-DBLE(NEl/SFac))>1.D-10) THEN
        CALL Warn("[P2Use:"//TRIM(SCFActn)//"] TrP = "//TRIM(DblToChar(TrP)))
      ENDIF
    ENDIF

    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(X,TrixFile('X',Args))    ! T1=S_new^(-1/2)
      CALL Multiply(X,P,Tmp1)           ! T2=S_new^(-1/2).P_old
      CALL Multiply(Tmp1,X,Tmp2)        ! P_new_AO=S_new^(-1/2).P_old.S_new^(-1/2)
      CALL Filter(P,Tmp2)               ! T1=Filter[P_new_AO]
    ELSE
      CALL Get(X,TrixFile('Z',Args))    ! T1=Z_new
      CALL Multiply(X,P,Tmp1)           ! T2=Z.P_old
      CALL Get(X,TrixFile('ZT',Args))   ! T1=Z^T
      CALL Multiply(Tmp1,X,Tmp2)        ! P_new_AO=Z.P_old.Z^T
      CALL Filter(P,Tmp2)               ! T1=Filter[P_new_AO]
    ENDIF
    ! IO for the non-orthogonal P
    CALL Put(P,TrixFile('D',Args,0))
    CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
    CALL PPrint( P,'P['//TRIM(Cycl)//']')
    CALL Plot(   P,'P_'//TRIM(Cycl))
    CALL Delete(P)
    CALL Delete(X)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! Restarting without Geometry of BasisSet Change
  CASE('Restart')

    CALL New(P,NSMat_O=NSMat)
    CALL New(S)
    CALL New(Tmp1)
    CALL New(Tmp2)
    ! Close Current Group
    CALL CloseHDFGroup(H5GroupID)
    CALL CloseHDF(HDFFileID)
    ! Open old group and HDF
    OldFileID=OpenHDF(Restart)
    HDF_CurrentID=OpenHDF(Restart)
    ! Get old basis set stuff
    CALL New(Stat,3)
    CALL Get(Stat,'current_state')
    SCFCycl=TRIM(IntToChar(Stat%I(1)))
    CurBase=TRIM(IntToChar(Stat%I(2)))
    CurGeom=TRIM(IntToChar(Stat%I(3)))
    ! Open the old group
    HDF_CurrentID=OpenHDFGroup(OldFileID,"Clone #"//TRIM(IntToChar(MyClone)))
    ! Get the old AO-dM
    CALL Halt(' Bad logic in P2Use.  HDF Checkpoined DM is disabled for now. ')
!!!    CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Get(S,TrixFile('S',Args))
    CALL Multiply(P,S,Tmp1)
    TError0 = ABS(SFac*Trace(Tmp1)-DBLE(NEl))/DBLE(NEl)
    IF(TError0>1D-3) THEN
      !Print out the checksum for P, S and P*S.
      CALL PChkSum(P   ,'CurrentDM',Prog)
      CALL PChkSum(S   ,'Overlap  ',Prog)
      CALL PChkSum(Tmp1,'P*S      ',Prog)
      CALL Halt(' Possible geometry, density matrix mismatch '//Rtrn// &
           ' Relative error in Tr(S.P)-Nel = '//DblToMedmChar(TError0)//Rtrn// &
           ' NEl = '//DblToMedmChar(DBLE(NEl))//Rtrn// &
           ' Try Guess=(Restart,Reparse) with previous geometry ' )
    ENDIF

    !CALL AOSP2(P,S,Tmp1,Tmp2,.TRUE.)
    !CALL AOSP2(P,S,Tmp1,Tmp2,.FALSE.)
    !CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog,Unit_O=6)

    ! IO for the non-orthogonal P
    CALL Put(P,TrixFile('D',Args,0))
    CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
    CALL PPrint( P,'P['//TRIM(Cycl)//']')
    CALL Plot(   P,'P_'//TRIM(Cycl))
    ! Close Old group
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(OldFileID)
    ! Reopen current group and HDF
    HDFFileID=OpenHDF(H5File)
    H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID=H5GroupID
    ! Put the DM into the hdf

!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(P,'CurrentDM',CheckPoint_O=.TRUE.)

    ! Clean Up
    CALL Delete(P)
    !CALL Delete(S)
    !CALL Delete(Tmp1)
    !CALL Delete(Tmp2)

    ! Restarting with BasisSet Change
  CASE('RestartBasisSwitch')

    ! Close Current Group
    CALL CloseHDFGroup(H5GroupID)
    CALL CloseHDF(HDFFileID)
    ! Open old group and HDF
    OldFileID=OpenHDF(Restart)
    HDF_CurrentID=OpenHDF(Restart)
    ! Get old basis set stuff
    CALL New(Stat,3)
    CALL Get(Stat,'current_state')
    SCFCycl=TRIM(IntToChar(Stat%I(1)))
    CurBase=TRIM(IntToChar(Stat%I(2)))
    CurGeom=TRIM(IntToChar(Stat%I(3)))
    ! Open the old group
    HDF_CurrentID=OpenHDFGroup(OldFileID,"Clone #"//TRIM(IntToChar(MyClone)))
    ! Get the Old Basis Set and Geometry
    CALL Get(BS_old,CurBase)
    CALL Get(GM_old,CurGeom)
    ! Compute a new sparse matrix blocking scheme for the old BS
    CALL BlockBuild(GM_old,BS_old,BSiz,OffS)
    NBasF=BS_old%NBasF
#ifdef PARALLEL
    CALL BCast(BSiz)
    CALL BCast(OffS)
    CALL BCast(NBasF)
#endif
    ! Get the old AO-DM
    CALL Halt(' Bad logic in P2Use.  HDF Checkpoined DM is disabled for now. ')
!!    CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
    ! IO for the non-orthogonal P
    CALL Put(P,TrixFile('D',Args,0))
    CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
    CALL PPrint( P,'P['//TRIM(Cycl)//']')
    CALL Plot(   P,'P_'//TRIM(Cycl))
    ! Close Old group
    CALL CloseHDFGroup(HDF_CurrentID)
    CALL CloseHDF(OldFileID)
    ! Reopen current group and HDF
    HDFFileID=OpenHDF(H5File)
    H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID=H5GroupID
    CALL Delete(P)

    !
    ! Density Matrix Verlet's
    !
  CASE("DMLinear")

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO < 3) THEN
      CALL Halt('[P2Use:Linear] No previous density matrix defined')
    ENDIF

    ! P(n) = 2 D(n-1) - D(n-2)

    ! Get D(p-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBas, iGEO-1 /)))
    CALL Multiply(Tmp1,2.0D0)
    CALL SetEq(P,Tmp1)

    ! Get D(p-2)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1,-1.0D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Purify P
!    MM = 0
!    CALL New(P0)
!    CALL SetEq(P0,P)
!    ! Do SP2 iterations
!    Occ0 = 0.D0
!    Occ1 = 0.D0
!    Occ2 = 0.D0
!    Occ3 = 0.D0
!    Imin = 4
!    DO I=1,100
!      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
!      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
!      Occ3 = Occ2
!      Occ2 = Occ1
!      Occ1 = Occ0
!    ENDDO
!    CALL Delete(P0)
!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P)  = "//TRIM(DblToChar(TrP)))
!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

  CASE("DMTRBO")

    ! Get the current time step size.
    CALL Get(MDDeltaTime, "MDDeltaTime", Tag_O = TRIM(IntToChar(iGEO)))
    CALL MondoLog(DEBUG_NONE, logtag, "MDDeltaTime = "//TRIM(FltToChar(MDDeltaTime*InternalTimeToFemtoseconds)))

    IF(iGEO < 4) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL Get(alpha, "MDalpha")
    CALL Get(MDDampStep, "MDDampStep")

    IF(iGEO > MDDampStep) THEN
      ! Initial damping in case of noise.
      alpha = 0.0D0
    ENDIF

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO == 4) THEN
      ! Initial boundary conditions: Save D(p-1) as P(p-1).
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,3
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO
    ENDIF

    ! Debugging: check.... P(n-1)-D(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get P(n-1)
      CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF

    ! Debugging: check.... P(n-1)-D_tilde(n-1)
    !
    ! Get D_tilde(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, -1.0D0)

    ! Get P(n-1)
    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Add(Tmp1,Tmp2,P)

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D_tilde) = "//TRIM(DblToChar(FNorm(P))))

    ! Debugging: check.... D(n-1)-D_tilde(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get D_tilde(n-1)
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(D-D_tilde) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF
    ! End Debugging.

    ! Notation:
    !
    ! P(n) = Auxiliary non-self-consistent dynamical variable. (DOPsave)
    ! D(n) = SCF[P(n)] (DOsave)
    !
    ! P(n) = 2*D(n-1)-P(n-2)-0.5*alpha*(D(n-1)-2*P(n-2)+D(n-3))

    ! Get D(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 2.0D0-0.5D0*alpha)
    CALL SetEq(P,Tmp1)

    ! Get P(n-2)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1,-1.0D0+alpha)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get D(n-3)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-3 /)))
    CALL Multiply(Tmp1,-0.5D0*alpha)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Save P(p,0)
    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P
    MM = 0
    CALL New(P0)
    CALL SetEq(P0,P)
    ! Do SP2 iterations
    Occ0 = 0.D0
    Occ1 = 0.D0
    Occ2 = 0.D0
    Occ3 = 0.D0
    Imin = 4
    DO I=1,100
      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
      Occ3 = Occ2
      Occ2 = Occ1
      Occ1 = Occ0
    ENDDO
    CALL Delete(P0)
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P)  = "//TRIM(DblToChar(TrP)))
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! DMTRBO with solid super-duper high order damping.
  CASE("DMTRBO_Damp_dt3")

    IF(iGEO < 5) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO == 5) THEN
      ! Initial boundary conditions: Save D(p-1) as P(p-1).
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,4
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO
    ENDIF

    ! Debugging: check.... P(n-1)-D(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get P(n-1)
      CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF

    ! Debugging: check.... P(n-1)-D_tilde(n-1)
    !
    ! Get D_tilde(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, -1.0D0)

    ! Get P(n-1)
    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Add(Tmp1,Tmp2,P)

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D_tilde) = "//TRIM(DblToChar(FNorm(P))))

    ! Debugging: check.... D(n-1)-D_tilde(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get D_tilde(n-1)
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(D-D_tilde) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF
    ! End Debugging.

    !  P(n+1) = 1.692*D(n) + 0.008*P(n) - 0.55*P(n-1) + 0*P(n-2) - 0.15*P(n-3)
    !
    ! We reverse the addition order, sorted from smallest prefactor to largest,
    ! for numerical reasons.

    ! Get P(n-3)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-4 /)))
    CALL Multiply(Tmp1, -0.15D0)
    CALL SetEq(P,Tmp1)

    ! Get P(n-1)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1, -0.55D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 0.008D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get D(n)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 1.692D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Save P(p,0)
    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P
    MM = 0
    CALL New(P0)
    CALL SetEq(P0,P)
    ! Do SP2 iterations
    Occ0 = 0.D0
    Occ1 = 0.D0
    Occ2 = 0.D0
    Occ3 = 0.D0
    Imin = 4
    DO I=1,100
      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
      Occ3 = Occ2
      Occ2 = Occ1
      Occ1 = Occ0
    ENDDO
    CALL Delete(P0)
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(DblToChar(TrP)))
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk

!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! DMTRBO with solid super-duper high order damping.
  CASE("DMTRBO_Damp_dt5")

    IF(iGEO < 6) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO == 6) THEN
      ! Initial boundary conditions: Save D(p-1) as P(p-1).
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,5
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO
    ENDIF

    ! Debugging: check.... P(n-1)-D(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get P(n-1)
      CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF

    ! Debugging: check.... P(n-1)-D_tilde(n-1)
    !
    ! Get D_tilde(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, -1.0D0)

    ! Get P(n-1)
    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Add(Tmp1,Tmp2,P)

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D_tilde) = "//TRIM(DblToChar(FNorm(P))))

    ! Debugging: check.... D(n-1)-D_tilde(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get D_tilde(n-1)
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(D-D_tilde) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF
    ! End Debugging.

    !  P(n+1) = 1.75*D(n) + 0.079*P(n) - 0.658*P(n-1) - 0.114*P(n-2)
    !           - 0.114*P(n-3) + 0.057*P(n-4)
    !
    ! We reverse the addition order, sorted from smallest prefactor to largest,
    ! for numerical reasons.

    ! Get P(n-4)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-5 /)))
    CALL Multiply(Tmp1, 0.057D0)
    CALL SetEq(P,Tmp1)

    ! Get P(n-3)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-4 /)))
    CALL Multiply(Tmp1, -0.114D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-2)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-3 /)))
    CALL Multiply(Tmp1, -0.114D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-1)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1, -0.658D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 0.079D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get D(n)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 1.75D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Save P(p,0)
    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P
    MM = 0
    CALL New(P0)
    CALL SetEq(P0,P)
    ! Do SP2 iterations
    Occ0 = 0.D0
    Occ1 = 0.D0
    Occ2 = 0.D0
    Occ3 = 0.D0
    Imin = 4
    DO I=1,100
      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
      Occ3 = Occ2
      Occ2 = Occ1
      Occ1 = Occ0
    ENDDO
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(DblToChar(TrP)))
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! DMTRBO with super-duper high order damping.
  CASE("DMTRBO_Damp_dt7")

    IF(iGEO < 7) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO == 7) THEN
      ! Initial boundary conditions: Save D(p-1) as P(p-1).
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,6
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO
    ENDIF

    ! Debugging: check.... P(n-1)-D(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get P(n-1)
      CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF

    ! Debugging: check.... P(n-1)-D_tilde(n-1)
    !
    ! Get D_tilde(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, -1.0D0)

    ! Get P(n-1)
    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Add(Tmp1,Tmp2,P)

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D_tilde) = "//TRIM(DblToChar(FNorm(P))))

    ! Debugging: check.... D(n-1)-D_tilde(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get D_tilde(n-1)
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(D-D_tilde) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF
    ! End Debugging.

    !  P(n+1) = 1.804*D(n) + 0.088*P(n) - 0.748*P(n-1) - 0.144*P(n-2)
    !           - 0.054*P(n-3) + 0.072*P(n-4) - 0.018*P(n-5)
    !
    ! We reverse the addition order, sorted from smallest prefactor to largest,
    ! for numerical reasons.

    ! Get P(n-5)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-6 /)))
    CALL Multiply(Tmp1, -0.018D0)
    CALL SetEq(P,Tmp1)

    ! Get P(n-4)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-5 /)))
    CALL Multiply(Tmp1, 0.072D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-3)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-4 /)))
    CALL Multiply(Tmp1, -0.054D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-2)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-3 /)))
    CALL Multiply(Tmp1, -0.144D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-1)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1, -0.748D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 0.088D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get D(n)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 1.804D0)
    CALL Add(P,Tmp1,Tmp2)

    ! WE NEED TO CHECK PUTTING A FILTER HERE!!
    CALL SetEq(P,Tmp2)

    ! Save P(p,0)
    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P

    ! For Fermi-Dirac or temperature smearing, do _not_ purify! Else (T=0) do.
!!    CALL New(P0)
!!    CALL SetEq(P0,P)
!!    ! Do SP2 iterations
!!    Occ0 = 0.D0
!!    Occ1 = 0.D0
!!    Occ2 = 0.D0
!!    Occ3 = 0.D0
!!    Imin = 4
!!    DO I=1,100
!!      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
!!      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
!!      Occ3 = Occ2
!!      Occ2 = Occ1
!!      Occ1 = Occ0
!!    ENDDO
!!    CALL Delete(P0)
!!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
!!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(DblToChar(TrP)))
!!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! DMTRBO with even more super-duper high order damping.
  CASE("DMTRBO_Damp_dt9")

    IF(iGEO < 8) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO == 8) THEN
      ! Initial boundary conditions: Save D(p-1) as P(p-1).
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,7
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO
    ENDIF

    ! Debugging: check.... P(n-1)-D(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get P(n-1)
      CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF

    ! Debugging: check.... P(n-1)-D_tilde(n-1)
    !
    ! Get D_tilde(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, -1.0D0)

    ! Get P(n-1)
    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Add(Tmp1,Tmp2,P)

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D_tilde) = "//TRIM(DblToChar(FNorm(P))))

    ! Debugging: check.... D(n-1)-D_tilde(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get D_tilde(n-1)
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(D-D_tilde) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF
    ! End Debugging.

    ! P(n+1) = 1.838*D(n) + 0.085*P(n) - 0.802*P(n-1) - 0.1485*P(n-2)
    !          - 0.011*P(n-3) + 0.066*P(n-4) - 0.033*P(n-5) + 0.0055*P(n-6)
    !
    ! We reverse the addition order, sorted from smallest prefactor to largest,
    ! for numerical reasons.

    ! Get P(n-6)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-7 /)))
    CALL Multiply(Tmp1, 0.0055D0)
    CALL SetEq(P,Tmp1)

    ! Get P(n-5)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-6 /)))
    CALL Multiply(Tmp1, -0.033D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-4)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-5 /)))
    CALL Multiply(Tmp1, 0.066D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-3)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-4 /)))
    CALL Multiply(Tmp1, -0.011D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-2)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-3 /)))
    CALL Multiply(Tmp1, -0.1485D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-1)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1, -0.802D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 0.085D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get D(n)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 1.838D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Save P(p,0)
    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P
!    MM = 0
!    CALL New(P0)
!    CALL SetEq(P0,P)
!    ! Do SP2 iterations
!    Occ0 = 0.D0
!    Occ1 = 0.D0
!    Occ2 = 0.D0
!    Occ3 = 0.D0
!    Imin = 4
!    DO I=1,100
!      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
!      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
!      Occ3 = Occ2
!      Occ2 = Occ1
!      Occ1 = Occ0
!    ENDDO
!    CALL Delete(P0)
!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(DblToChar(TrP)))
!    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! DMTRBO with industrial strength super-duper high order damping.
  CASE("DMTRBO_Damp_dt11")

    IF(iGEO < 9) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO == 9) THEN
      ! Initial boundary conditions: Save D(p-1) as P(p-1).
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,8
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO
    ENDIF

    ! Debugging: check.... P(n-1)-D(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get P(n-1)
      CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF

    ! Debugging: check.... P(n-1)-D_tilde(n-1)
    !
    ! Get D_tilde(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, -1.0D0)

    ! Get P(n-1)
    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Add(Tmp1,Tmp2,P)

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(P-D_tilde) = "//TRIM(DblToChar(FNorm(P))))

    ! Debugging: check.... D(n-1)-D_tilde(n-1)
    !
    ! Get D(n-1)
    INQUIRE(FILE=TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1), EXIST=Present)
    IF(Present) THEN
      CALL Get(Tmp1, TrixFile("OrthoD", Stats_O = Args%I%I(4:6), Offset_O = 1))
      CALL Multiply(Tmp1, -1.0D0)

      ! Get D_tilde(n-1)
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Add(Tmp1,Tmp2,P)

      CALL MondoLog(DEBUG_MAXIMUM, logtag, "FNorm(D-D_tilde) = "//TRIM(DblToChar(FNorm(P))))
    ENDIF
    ! End Debugging.

    !  P(n+1) = 1.861*D(n) + 0.0814*P(n) - 0.8416*P(n-1) - 0.1408*P(n-2)
    !           + 0.0176*P(n-3) + 0.0512*P(n-4) - 0.04*P(n-5) + 0.0128*P(n-6)
    !           - 0.0016*P(n-7)
    !
    ! We reverse the addition order, sorted from smallest prefactor to largest,
    ! for numerical reasons.

    ! Get P(n-7)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-8 /)))
    CALL Multiply(Tmp1, -0.0016D0)
    CALL SetEq(P,Tmp1)

    ! Get P(n-6)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-7 /)))
    CALL Multiply(Tmp1, 0.0128D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-5)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-6 /)))
    CALL Multiply(Tmp1, -0.04D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-4)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-5 /)))
    CALL Multiply(Tmp1, 0.0512D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-3)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-4 /)))
    CALL Multiply(Tmp1, 0.0176D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-2)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-3 /)))
    CALL Multiply(Tmp1, -0.1408D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n-1)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
    CALL Multiply(Tmp1, -0.8416D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get P(n)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 0.0814D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Get D(n)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1, 1.861D0)
    CALL Add(P,Tmp1,Tmp2)
    CALL SetEq(P,Tmp2)

    ! Save P(p,0)
    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P
    MM = 0
    CALL New(P0)
    CALL SetEq(P0,P)
    ! Do SP2 iterations
    Occ0 = 0.D0
    Occ1 = 0.D0
    Occ2 = 0.D0
    Occ3 = 0.D0
    Imin = 4
    DO I=1,100
      CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
      IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) EXIT
      Occ3 = Occ2
      Occ2 = Occ1
      Occ1 = Occ0
    ENDDO
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(DblToChar(TrP)))
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    ! DMSymplectic for 4th order symplecitc integration scheme 4.617
  CASE("DMSymplectic")

    IF(iGEO < 7) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "No previous density matrix defined")
      CALL Halt("["//TRIM(logtag)//"] Fatal error")
    ENDIF

    CALL Get(alpha, "MDalpha")
    CALL Get(MDDampStep, "MDDampStep")

    IF(iGEO > MDDampStep) THEN
      alpha = 0.0D0
    ENDIF

    CALL New(P)
    CALL New(PV)
    CALL New(Tmp1)
    CALL New(Tmp2)

    ! Compute Pnew in ortho space: This is specifically for MD.  Save D(p-1) as
    ! P(p-1), where p < 5.

    ! Calculate symplectic counter.
    m_step = MOD(iGEO-2,4)+1
    IF(iGEO == 7) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Initial boundary condition")
      DO I=1,6
        CALL Get(Tmp1, TrixFile("DOsave",  Stats_O = (/ iSCF, iBAS, iGEO-I /)))
        CALL Put(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-I /)))
      ENDDO

      ! PV(2) = PV(1) + 4.61D0*b(1)*[D(1)-P(1)]
      ! P(2) = P(1) + a(1)*PV(2)

      ! PV(3) = PV(2) + 4.61D0*b(2)*[D(2)-P(2)]
      ! P(3) = P(2) + a(2)*PV(3)

      ! PV(4) = PV(3) + 4.61D0*b(3)*[D(3)-P(3)]
      ! P(4) = P(3) + a(3)*PV(4)

      ! PV(5) = PV(4) + 4.61D0*b(4)*[D(4)-P(4)]
      ! P(5) = P(4) + a(5)*PV(5)

      ! PV(n) = PV(n-1) + 4.61D0*b(m_step)*[D(n-1)-P(n-1)]
      ! P(n) = P(n-1) + a(m_step)*PV(n)

      ! mm = MOD(iGEO-1-2,4)+1
      ! P(n-1) = P(n-2) + a(mm)*PV(n-1)
      ! PV(n-1) = [P(n-1)-P(n-2)]/a(mm)

      ! Calculate initial velocity v(n-1) = [P(n-1)-P(n-2)]/a(mm)
      mm_step = MOD(iGEO-3,4)+1

      CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Multiply(Tmp1,1.D0/Symplectic_4th_Order_a(mm_step))
      CALL Get(Tmp2, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
      CALL Multiply(Tmp2,-1.D0/Symplectic_4th_Order_a(mm_step))
      CALL Add(Tmp1,Tmp2,PV)
      CALL Multiply(PV,1.D0)

      CALL Put(PV, TrixFile("PVsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))

    ENDIF

    ! PV(n) = PV(n-1) + 4.61D0*b(m_step)*[D(n-1)-P(n-1)]
    ! P(n) = P(n-1) + a(m_step)*PV(n)

    ! Get D(p-1)
    CALL Get(PV, TrixFile("PVsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))

    ! Get D(n-1), P(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1,4.61D0*Symplectic_4th_Order_b(m_step))

    CALL Get(Tmp2, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp2,-4.61D0*Symplectic_4th_Order_b(m_step))
    CALL Add(Tmp1,Tmp2,P)

#ifdef ANDERS
    beta = 0.00275D0
    IF(iGEO > MDDampStep) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Anders beta hack active now")
      CALL Multiply(P,1.D0-beta)

      CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
      CALL Multiply(Tmp1,beta*Symplectic_4th_Order_b(m_step)*(-14.D0/5.D0))
      CALL Add(P,Tmp1,Tmp2)
      CALL SetEq(P,Tmp2)

      CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-2 /)))
      CALL Multiply(Tmp1,beta*Symplectic_4th_Order_b(m_step)*(42.D0/5.D0))
      CALL Add(P,Tmp1,Tmp2)
      CALL SetEq(P,Tmp2)

      CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-3 /)))
      CALL Multiply(Tmp1,beta*Symplectic_4th_Order_b(m_step)*(-48.D0/5.D0))
      CALL Add(P,Tmp1,Tmp2)
      CALL SetEq(P,Tmp2)

      CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-4 /)))
      CALL Multiply(Tmp1,beta*Symplectic_4th_Order_b(m_step)*(27.D0/5.D0))
      CALL Add(P,Tmp1,Tmp2)
      CALL SetEq(P,Tmp2)

      CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-5 /)))
      CALL Multiply(Tmp1,beta*Symplectic_4th_Order_b(m_step)*(-8.D0/5.D0))
      CALL Add(P,Tmp1,Tmp2)
      CALL SetEq(P,Tmp2)

      CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-6 /)))
      CALL Multiply(Tmp1,beta*Symplectic_4th_Order_b(m_step)*(1.D0/5.D0))
      CALL Add(P,Tmp1,Tmp2)
      CALL SetEq(P,Tmp2)
    ENDIF
#endif

    CALL Add(PV,P,Tmp1)
    CALL SetEq(PV,Tmp1)

    ! PV(n) = PV(n-1) + 4.61D0*b(m_step)*[D(n-1)-P(n-1)]
    CALL Put(PV, TrixFile("PVsave", Args))

    ! Get P(n-1)
    CALL Get(Tmp1, TrixFile("DOPsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1,1.D0-alpha)

    ! P(n) = P(n-1) + a(m_step)*PV(n)
    CALL Multiply(PV,Symplectic_4th_Order_a(m_step))
    CALL Add(Tmp1,PV,Tmp2)

    ! Get D(n-1)
    CALL Get(Tmp1, TrixFile("DOsave", Stats_O = (/ iSCF, iBAS, iGEO-1 /)))
    CALL Multiply(Tmp1,alpha)
    CALL Add(Tmp1,Tmp2,P)

    CALL Put(P, TrixFile("DOPsave", Args))

    ! Purify P
#ifdef PARALLEL
    CALL SetEq(P_BCSR,P)
    CALL SpectralBounds(P_BCSR,Fmin,Fmax)
    CALL Delete(P_BCSR)
#else
    CALL SpectralBounds(P,Fmin,Fmax)
#endif
    CALL Add(P,-Fmin)
    CALL Multiply(P,One/(Fmax-Fmin))
    MM = 0
    DO I=1,40
      CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
      IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
        IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
      ENDIF
    ENDDO
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "purified after "//TRIM(IntToChar(I))//" iterations")
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(DblToChar(TrP)))
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P2) = "//TRIM(DblToChar(TrP2)))

    ! Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)

    ! Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(PV)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

  CASE('FMVerlet0')

    CALL New(F)
    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO .LE. 2) THEN
      CALL Halt("[P2Use:FMVerlet0] No Previous Fock Matrix Defined")
    ENDIF

    !    Compute Fnew in ortho space: This is specifically for MD
    !
    !    Save F(p-1,n) as F(p-1,0), when p==3
    IF(iGEO==3) THEN
      DO I=1,2
        FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
        CALL Get(Tmp1,FileName)
        FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
        CALL Put(Tmp1,FileName)
      ENDDO
    ENDIF
    !    F(p,0) = 2*F(p-1,n)-F(p-2,0)
    !    Get F(p-1,n)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-1))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
    CALL Get(Tmp1,FileName)
    CALL Multiply(Tmp1, 2.0D0)
    CALL SetEq(F,Tmp1)
    !    Get F(p-2,0)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-2))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
    CALL Get(Tmp1,FileName)
    CALL Multiply(Tmp1,-1.0D0)
    CALL Add(F,Tmp1,Tmp2)
    CALL SetEq(F,Tmp2)
    !    Save F(p,0)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO  ))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
    CALL Put(F,FileName)
    !
    !    Solve for the Density Matrix using SP2
    !
    !    Guess P from F
#ifdef PARALLEL
    CALL SetEq(F_BCSR,F)
    CALL FockGuess(F_BCSR,P,Half*DBLE(NEl),1)
    CALL Delete(F_BCSR)
#else
    CALL FockGuess(F,P,Half*DBLE(NEl),1)
#endif
    !    Solve for P
    DO I=1,40
      CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
      IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
        IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
      ENDIF
    ENDDO
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(FltToChar(TrP)))

    !    Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)
    !    Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)
    !    Clean Up
    CALL Delete(F)
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

  CASE('FMVerlet1')

    CALL New(F)
    CALL New(P)
    CALL New(Tmp1)
    CALL New(Tmp2)

    IF(iGEO .LE. 4) THEN
      CALL Halt("[P2Use:FMVerlet1] No Previous Fock Matrix Defined")
    ENDIF

    !    Compute Fnew in ortho space: This is specifically for MD
    !
    !    Save F(p-1,n) as F(p-1,0), when p==3
    IF(iGEO==5) THEN
      DO I=1,4
        FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
        CALL Get(Tmp1,FileName)
        FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
        CALL Put(Tmp1,FileName)
      ENDDO
    ENDIF
    !    F(p,0) = 4*P(p-1,n)-6*P(p-2,n)+4*P(p-3,n)-P(p-4,0)
    !    Get F(p-1,n)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-1))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
    CALL Get(Tmp1,FileName)
    CALL Multiply(Tmp1, 4.0D0)
    CALL SetEq(F,Tmp1)
    !    Get F(p-2,0)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-2))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
    CALL Get(Tmp1,FileName)
    CALL Multiply(Tmp1,-6.0D0)
    CALL Add(F,Tmp1,Tmp2)
    CALL SetEq(F,Tmp2)
    !    Get F(p-3,0)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-3))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
    CALL Get(Tmp1,FileName)
    CALL Multiply(Tmp1, 4.0D0)
    CALL Add(F,Tmp1,Tmp2)
    CALL SetEq(F,Tmp2)
    !    Get F(p-4,0)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-4))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
    CALL Get(Tmp1,FileName)
    CALL Multiply(Tmp1,-1.0D0)
    CALL Add(F,Tmp1,Tmp2)
    CALL SetEq(F,Tmp2)
    !    Save F(p,0)
    FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO  ))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
    CALL Put(F,FileName)
    !
    !
    !    Solve for the Density Matrix using SP2
    !
    !    Guess P from F
#ifdef PARALLEL
    CALL SetEq(F_BCSR,F)
    CALL FockGuess(F_BCSR,P,Half*DBLE(NEl),1)
    CALL Delete(F_BCSR)
#else
    CALL FockGuess(F,P,Half*DBLE(NEl),1)
#endif
    !    Solve for P
    DO I=1,40
      CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
      IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
        IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
      ENDIF
    ENDDO
    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace(P) = "//TRIM(FltToChar(TrP)))

    !    Convert to AO Rep
    INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
    IF(Present)THEN
      CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Multiply(Tmp2,Tmp1,P)
    ELSE
      CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Tmp1,P,Tmp2)
      CALL Get(Tmp1,TrixFile('ZT',Args))
      CALL Multiply(Tmp2,Tmp1,P)
    ENDIF
    CALL Filter(Tmp1,P)
    !    Put to Disk
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
    CALL Put(Tmp1,TrixFile('D',Args,0))
    CALL PChkSum(Tmp1,'P[0]',Prog)
    !    Clean Up
    CALL Delete(F)
    CALL Delete(P)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

    !
    !    Geometry Change
    !
  CASE('Extrapolate','DMProj0','DMProj1','DMProj2','DMProj3','DMProj4')

    ! Allocate
    CALL New(P)
    CALL New(P0)
    CALL New(S)
    CALL New(S0)
    CALL New(S1)
    CALL New(T1)
    CALL New(T2)
    CALL New(Tmp1)
    CALL New(Tmp2)

    ! Get Matrices
    DO ICycle=1,1000
      DMFile=TRIM(SCRName)//'_Geom#'//TRIM(IntToChar(Current(3)-1)) &
           //'_Base#'//TRIM(IntToChar(Current(2))) &
           //'_Cycl#'//TRIM(IntToChar(ICycle)) &
           //'_Clone#'//TRIM(IntToChar(MyClone)) &
           //'.D'
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Looking for "//TRIM(DMFile))
      INQUIRE(FILE=DMFile,EXIST=Present)
      IF(.NOT.Present)THEN
        Cycle=ICycle-1
        DMFile=TRIM(SCRName)//'_Geom#'//TRIM(IntToChar(Current(3)-1)) &
             //'_Base#'//TRIM(IntToChar(Current(2))) &
             //'_Cycl#'//TRIM(IntToChar(Cycle)) &
             //'_Clone#'//TRIM(IntToChar(MyClone)) &
             //'.D'
        CALL MondoLog(DEBUG_MAXIMUM, logtag, "On extrapolation, P2Use is opening DM "//TRIM(DMFile))
        EXIT
      ENDIF
    ENDDO

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "Cycle = "//TRIM(IntToChar(Cycle)))

    IF(Cycle<0)THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Assuming this is a restart! If there is no restart, its going to die...")
      CALL Get(S0,TrixFile('S',Args,Stats_O=(/Current(1),Current(2),Current(3)-1/)))
      CALL Get(S1,TrixFile('S',Args,Stats_O=Current))
      ! Close Current Group
      CALL CloseHDFGroup(H5GroupID)
      CALL CloseHDF(HDFFileID)
      ! Open old group and HDF
      OldFileID=OpenHDF(Restart)
      HDF_CurrentID=OpenHDF(Restart)
      ! Get old basis set stuff
      CALL New(Stat,3)
      CALL Get(Stat,'current_state')
      SCFCycl=TRIM(IntToChar(Stat%I(1)))
      CurBase=TRIM(IntToChar(Stat%I(2)))
      CurGeom=TRIM(IntToChar(Stat%I(3)))
      ! Open the old group
      HDF_CurrentID=OpenHDFGroup(OldFileID,"Clone #"//TRIM(IntToChar(MyClone)))
      ! Get the old AO-DM
      CALL Halt(' Bad logic in P2Use.  HDF Checkpoined DM is disabled for now. ')
!!!      CALL Get(P0,'CurrentDM',CheckPoint_O=.TRUE.)
      ! Close Old group
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(OldFileID)
      ! Reopen current group and HDF
      HDFFileID=OpenHDF(H5File)
      H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID=H5GroupID
    ELSE
      CALL Get(S0,TrixFile('S',Args,Stats_O=Previous))
      CALL Get(S1,TrixFile('S',Args,Stats_O=Current ))
      DMFile=TRIM(SCRName)//'_Geom#'//TRIM(IntToChar(Current(3)-1)) &
           //'_Base#'//TRIM(IntToChar(Current(2))) &
           //'_Cycl#'//TRIM(IntToChar(Cycle)) &
           //'_Clone#'//TRIM(IntToChar(MyClone)) &
           //'.D'
      CALL Get(P0,DMFile)
      CALL Get(DoingMD ,'DoingMD')
      IF(DoingMD) THEN
        DMPOrder=0
        CALL Get(MDGeuss ,"MDGeuss")
        SELECT CASE(MDGeuss)

        CASE('DMProj0')
          DMPOrder=0

        CASE('DMProj1')
          DMPOrder=1

        CASE('DMProj2')
          DMPOrder=2

        CASE('DMProj3')
          DMPOrder=3

        CASE('DMProj4')
          DMPOrder=4

        CASE DEFAULT
          CALL Halt("illegal DMPOrder ("//TRIM(MDGeuss)//")")
        END SELECT
        iGEO     = Args%I%I(3)
        DMPOrder = MIN(MAX(iGEO-2,0),DMPOrder)
        CALL DMPProj(iGEO,DMPOrder,P0,Tmp1,Tmp2)
      ENDIF
    ENDIF
    ! Initial Trace Error
#ifdef PARALLEL
    CALL Multiply(P0,S0,Tmp1)
    TError0 = ABS(Trace(Tmp1)-DBLE(NEl)/SFac)
    IF(MyID.EQ.ROOT) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace Error: Tr[P0,S0] = "//TRIM(IntToChar(TError0)))
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace Error: Tr[P0,S1] = "//TRIM(IntToChar(TError0)))
    ENDIF
    CALL Multiply(P0,S1,Tmp1)
    TError0 = ABS(Trace(Tmp1)-DBLE(NEl)/SFac)
    IF(MyID.EQ.ROOT)THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace Error: Tr[P0,S1] = "//TRIM(IntToChar(TError0)))
    ENDIF
#else
    TError0 = ABS(Trace(P0,S1)-DBLE(NEl)/SFac)
    IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace Error: Tr[P0,S0] = "//TRIM(FltToChar(ABS(Trace(P0,S0)-DBLE(NEl)/SFac))))
      CALL MondoLog(DEBUG_MAXIMUM, logtag, "Trace Error: Tr[P0,S1] = "//TRIM(FltToChar(ABS(Trace(P0,S1)-DBLE(NEl)/SFac))))
    ENDIF
#endif

    ! Calculate norm error and idempotency error in last geometry.
    CALL Multiply(S0, P0, T1)  ! T1 = S0.P0
    CALL Multiply(P0, T1, T2)  ! T2 = P0.S0.P0
    TrP = Trace(P0, S0)
    TrP2 = Trace(T2, S0)

    Norm_Error = TrP-DBLE(NEl)/SFac
    Ipot_Error = TrP2-TrP

    Ipot_Error_start = Ipot_Error

    CALL MondoLog(DEBUG_MAXIMUM, logtag, "initial norm error = "//TRIM(DblToShrtChar(Norm_Error))// &
      ", initial idempotency error = "//TRIM(DblToShrtChar(Ipot_Error_start)))

    ! Initialize
    NStep       = 0
    Lam         = Zero
    DLam        = One

    ! Purify
    DO
      NStep = NStep + 1
      Lam  = Lam + DLam
      ! Set Up S and P
      CALL SetEq(Tmp1,S0)
      CALL SetEq(Tmp2,S1)
      CALL Multiply(Tmp1,One-Lam)
      CALL Multiply(Tmp2,Lam)
      CALL Add(Tmp1,Tmp2,S)
      CALL SetEq(P,P0)

#ifdef PARALLEL
      IF(MyId==ROOT)THEN
#endif
        CALL MondoLog(DEBUG_MEDIUM,Prog, "Nstep  = "//TRIM(IntToChar(NStep))// &
          ", Lambda = "//TRIM(FltToMedmChar(Lam)), "AO-DMX")

#ifdef PARALLEL
      ENDIF
#endif

#ifdef PARALLEL
      CALL Multiply(P,S,Tmp1)
      TrP = Trace(Tmp1)
#else
      TrP = Trace(P,S)
#endif
      Norm_Error = TrP-DBLE(NEl)/SFac
      Ipot_Error = One

      AOSPExit     = .FALSE.

      Occ0 = 0.D0
      Occ1 = TrP
      Occ2 = 0.D0
      Occ3 = 0.D0
      Imin_divergence_check = 6
      Imin = 10

      DO I = 1, 20
        ! Calculate new norm error.
        CALL Multiply(S, P, T1)   ! T1 = S.P
        CALL Multiply(P, T1, T2)  ! T2 = P.S.P

        TrP = Trace(P, S)
        TrP2 = Trace(T2, S)

        Norm_Error = TrP-DBLE(NEl)/SFac
        Ipot_Error = TrP2-TrP

        IF(ABS(2*TrP-TrP2-DBLE(Nel)/SFac) > ABS(TrP2-DBLE(Nel)/SFac)) THEN
          CALL AOSP2(P,S,Tmp1,Tmp2,.TRUE.)
        ELSE
          CALL AOSP2(P,S,Tmp1,Tmp2,.FALSE.)
        ENDIF

        Occ0 = Trace(P, S)

#ifdef PARALLEL
        IF(MyId==ROOT)THEN
#endif
          PNon0s=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
          CALL MondoLog(DEBUG_MEDIUM, Prog, 'norm error = '//TRIM(DblToMedmChar(Norm_Error))// &
            ", idempotency error = "//TRIM(DblToMedmChar(Ipot_Error))//", %Non0 = "//TRIM(FltToShrtChar(PNon0s)), &
            "AO-DMX "//TRIM(IntToChar(NStep))//", "//TRIM(IntToChar(I)))

#ifdef PARALLEL
        ENDIF
#endif
        ! Logic

        IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) THEN
          CALL MondoLog(DEBUG_MAXIMUM, Prog, "converged in "//TRIM(IntToChar(I))//" iterations")
          CALL MondoLog(DEBUG_MAXIMUM, Prog, "idempotency error = "//TRIM(DblToChar(ABS(Occ0-Occ1))))
          CALL MondoLog(DEBUG_MAXIMUM, Prog, "previous idempotency error = "//TRIM(DblToChar(ABS(Occ2-Occ3))))
          AOSPExit = .TRUE.
        ENDIF
        Occ3 = Occ2
        Occ2 = Occ1
        Occ1 = Occ0

        IF(I >= Imin_divergence_check) THEN
          IF(ABS(Occ0-Occ1)/ABS(Ipot_Error_start) >= 1.0D0 .AND. &
             ABS(Occ0-Occ1)/ABS(Ipot_Error_start) <  2.0D0) THEN
            CALL MondoLog(DEBUG_MAXIMUM, Prog, "converged in "//TRIM(IntToChar(I))//" iterations")
            CALL MondoLog(DEBUG_MAXIMUM, Prog, "Idempotency error between 1 and 2 times initial error")
            AOSPExit = .TRUE.
          ELSEIF(ABS(Occ0-Occ1)/ABS(Ipot_Error_start) > 10.0D0) THEN
            CALL MondoLog(DEBUG_MAXIMUM, Prog, "failed convergence")
            CALL MondoLog(DEBUG_MAXIMUM, Prog, "idempotency error diverging, greater 10 times initial error")
            AOSPExit = .TRUE.
          ENDIF
        ENDIF

        IF(AOSPExit) EXIT
      ENDDO

      ! Check convergence of the last Dlam step.
      IF(.NOT. AOSPExit .AND. ABS(Occ0-Occ1)/Abs(Ipot_Error_start) < 1.0D-2 .AND. ABS(Occ0-Occ1) < 1.0D-3) THEN
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "forced convergence, exceeded maximum number of inner iterations")
        CALL SetEq(P0, P)
      ELSEIF(.NOT. AOSPExit) THEN
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "could not converge in inner loop, repeating with half step")
        Lam  = Lam - DLam
        DLam = Half*DLam

        IF(DLam < 1.0D-2) THEN
          CALL Halt("[P2Use:"//TRIM(SCFActn)//"] DLam < 0.01, density matrix extrapolater failed to converge in "//TRIM(IntToChar(NStep))//" steps")
        ENDIF

      ELSEIF(AOSPExit) THEN
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "converged in this Dlam step")
        CALL SetEq(P0, P)
      ENDIF

      IF(Lam > (One-1.0D-14)) THEN
        CALL MondoLog(DEBUG_MAXIMUM, Prog, "density matrix converged in "//TRIM(IntToChar(NStep))//" steps")
        EXIT
      ENDIF

    ENDDO

    ! Save back to be sure.
    CALL PChkSum(P,'P[0]',Prog)
    CALL Put(P,TrixFile('D',Args,0))
!!!    DM archiveal has been put on hold, pending better
!!!    way of dealing with HDF in parallel.
!!!    CALL Put(P,'CurrentDM',CheckPoint_O=.TRUE.)

    ! Clean Up
    CALL Delete(P)
    CALL Delete(P0)
    CALL Delete(S)
    CALL Delete(S0)
    CALL Delete(S1)
    CALL Delete(T1)
    CALL Delete(T2)
    CALL Delete(Tmp1)
    CALL Delete(Tmp2)

  CASE DEFAULT
    CALL Halt("Unknown option "//TRIM(SCFActn))

  END SELECT

  ! Tidy up ...
  CALL Delete(GM)
  CALL Delete(BS)
  CALL ShutDown(Prog)

END PROGRAM P2Use
