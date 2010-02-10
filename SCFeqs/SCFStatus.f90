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

#include "MondoConfig.h"

PROGRAM SCFStatus
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
  USE Functionals
  USE MondoLogger
#ifdef PARALLEL
  USE MondoMPI
#endif

  IMPLICIT NONE

  TYPE(ARGMT)                     :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                     :: P,Tmp1,Tmp2,Tmp3
#else
  TYPE(BCSR)                      :: P,Tmp1,Tmp2,Tmp3
#endif
  REAL(DOUBLE)                    :: E_el_tot,E_Nuc_Tot,E_es_tot,E_ECPs,KinE,ExchE,Exc,Gap,Etot,DMax,Virial,DIISErr,S2,SFac,Entropy
#ifdef MMech
  REAL(DOUBLE)                    :: EBOND,EANGLE,ETorsion,ELJ,EOutOfPlane,MM_COUL,MM_ENERGY
  REAL(DOUBLE)                    :: E_C_EXCL,E_LJ_EXCL
#endif
  LOGICAL                         :: HasECPs
  CHARACTER(LEN=DEFAULT_CHR_LEN)  :: SCFMessage
  CHARACTER(LEN=DEFAULT_CHR_LEN)  :: SCFTag
  CHARACTER(LEN=9),PARAMETER      :: Prog='SCF'
#if defined(PARALLEL_CLONES)
  INTEGER                         :: oldClone, rank
#endif

  !  Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)

  ! Reset derived types.
  CALL Initialize(P)
  CALL Initialize(Tmp1)
  CALL Initialize(Tmp2)
  CALL Initialize(Tmp3)

  !---------------------------------------------
  !  Get the density matrix
  IF(SCFActn=='BasisSetSwitch' .OR. SCFActn=="RestartBasisSwitch") THEN
    ! If switching the density matrix or using a previous one from
    ! restart use i+1 density matrix--its all that is available
    CALL Get(P,TrixFile('D',Args,1))
  ELSE
    CALL Get(P,TrixFile('D',Args,0))
  ENDIF
  !---------------------------------------------
  ! Rescaling factor for R/U/G theory.
  SFac=1D0
  IF(P%NSMat.GT.1) SFac=0.5D0 !<<< SPIN
  !---------------------------------------------
  !  COMPUTE SOME EXPECTATION VALUES
  !
  !  S**2
  CALL Get(Tmp1,TrixFile('S',Args))
  S2 = GetS2(P,Tmp1,Tmp2)

  !  KinE=<T>=Tr{P.T}
  CALL Get(Tmp1,TrixFile('T',Args))
#ifdef PARALLEL
  CALL Multiply(P,Tmp1,Tmp2)
  KinE=Two*Trace(Tmp2)
#else
  KinE=Two*Trace(P,Tmp1)
#endif
  KinE=KinE*SFac !<<< SPIN

  CALL Get(HasECPs,'hasecps',Tag_O=CurBase)
  IF(HasECPs)THEN
    ! Get the pseudopotential matrix U
    CALL Get(Tmp1,TrixFile('U',Args))                        ! Tmp1=U_{ECP}
#ifdef PARALLEL
    CALL Multiply(P,Tmp1,Tmp2)
    E_ECPs=Two*Trace(Tmp2)
#else
    E_ECPs=Two*Trace(P,Tmp1)
#endif
  ELSE
    E_ECPs=Zero
  ENDIF
  E_ECPs=E_ECPs*SFac !<<< SPIN
  ! E_el_tot=<Vee+Vne>=Tr{P.(Vee+Vne)}
  CALL Get(Tmp1,TrixFile('J',Args,0))
#ifdef PARALLEL
  CALL Multiply(P,Tmp1,Tmp2)
  E_el_tot=Trace(Tmp2)
#else
  E_el_tot=Trace(P,Tmp1)
#endif
  E_el_tot=E_el_tot*SFac !<<< SPIN
  ! Total electrostatic energy icluding ECPs
  E_el_tot=E_el_tot+E_ECPs

  ExchE=Zero
  Exc=Zero
  IF(SCFActn/="GuessEqCore")THEN
    ! ExchE=<Kx>=Tr{P.K}
    IF(HasHF(ModelChem))THEN
      CALL Get(Tmp1,TrixFile('K',Args,0))
#ifdef PARALLEL
      CALL Multiply(P,Tmp1,Tmp2)
      ExchE=ExactXScale(ModelChem)*Trace(Tmp2)
#else
      ExchE=ExactXScale(ModelChem)*Trace(P,Tmp1)
#endif
    ENDIF
    ! Get the exchange correlation energy
    IF(HasDFT(ModelChem)) CALL Get(Exc,'Exc',Stats_O=Current)
  ENDIF
  ExchE=ExchE*SFac !<<< SPIN
  !  Exc = Exc*SFac !<<< SPIN
  !  Get E_Nuc_Tot =<Vnn+Vne>
  CALL Get(E_Nuc_Tot,'E_NuclearTotal',Stats_O=Current)
  ! Total electrostatic energy
  E_es_tot=E_el_tot+E_Nuc_Tot
  ! Total SCF energy
  Etot=KinE+E_es_tot+Exc+ExchE

  CALL MondoLog(DEBUG_MAXIMUM, Prog, "KinE      = "//TRIM(DblToChar(KinE))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "E_el_tot  = "//TRIM(DblToChar(E_el_tot))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "E_Nuc_Tot = "//TRIM(DblToChar(E_Nuc_Tot))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "E_es_tot  = "//TRIM(DblToChar(E_es_tot))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "Exc       = "//TRIM(DblToChar(Exc))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "Exch      = "//TRIM(DblToChar(ExchE))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))
  CALL MondoLog(DEBUG_MAXIMUM, Prog, "Etot      = "//TRIM(DblToChar(Etot))//" hartree", "Clone "//TRIM(IntToChar(MyClone)))

  !  The Virial
  Virial=E_es_tot/KinE
  !--------------------------------------------------------
  !  Find the largest block of the delta density matrix
  !  Allows for checking between extrapolated or projected DMs
  IF(SCFActn=='BasisSetSwitch'.OR. SCFActn=="RestartBasisSwitch" .OR. SCFActn=='NumForceEvaluation')THEN
    DMax=Max(P)
  ELSE
    CALL Get(Tmp1,TrixFile('D',Args,0))
    CALL Get(Tmp2,TrixFile('D',Args,1))
    CALL Multiply(Tmp1,-One)
    CALL Add(Tmp1,Tmp2,P)
    DMax=Max(P)
  ENDIF

#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL Put(E_el_tot, 'E_ElectronicTotal')
    CALL Put(Etot, 'Etot')
    CALL Put(Etot, 'Etot', Stats_O = Current)
    CALL Put(DMax, 'DMax')

    oldClone = MyClone
    DO rank = 1, MSize(MPI_COMM_WORLD)-1
      CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(E_el_tot, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(Etot, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(DMax, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

      ! Put to correct HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
      CALL Put(E_el_tot, 'E_ElectronicTotal')
      CALL Put(Etot, 'Etot')
      CALL Put(Etot, 'Etot', Stats_O = Current)
      CALL Put(DMax, 'DMax')
    ENDDO
    MyClone = oldClone

    ! Reopen old HDFGroup.
    CALL CloseHDFGroup(H5GroupID)
    H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID = H5GroupID
  ELSE
    CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(E_el_tot, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(Etot, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(DMax, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
  ENDIF
#else
  CALL Put(E_el_tot,'E_ElectronicTotal')
  CALL Put(Etot,'Etot')
  CALL Put(Etot,'Etot',Stats_O=Current)
  CALL Put(DMax,'DMax')
#endif

  !  IO for the delta density matrix
  IF(SCFActn=='InkFok')THEN
    CALL Put(P,TrixFile('DeltaD',Args,1))
    CALL PChkSum(P,'DeltaP['//TRIM(NxtCycl)//']',Prog)
    CALL PPrint( P,'DeltaP['//TRIM(NxtCycl)//']')
    CALL Plot(   P,'DeltaP_'//TRIM(NxtCycl))
  ENDIF
  !-----------------------------------------------------------
  !  Get DIIS Err and HOMO-LUMO Gap
  IF(Current(1)>=1.AND.SCFActn/='NumForceEvaluation') THEN
    CALL Get(DIISErr,'diiserr')
  ELSE
    DIISErr=Zero
  ENDIF
  CALL Get(Gap,'HomoLumoGap')
  CALL Get(Entropy, "Entropy")
  !-----------------------------------------------------------
  !  PRINT STATISTICS

  IF(NClones>1)THEN
     SCFTag='['//TRIM(SCFCycl)//','//TRIM(CurBase)//','//TRIM(CurGeom)//','//TRIM(CurClone)//']'
  ELSE
     SCFTag='['//TRIM(SCFCycl)//','//TRIM(CurBase)//','//TRIM(CurGeom)//']'
  ENDIF

  SCFMessage=""
  IF(Gap/=Zero)SCFMessage=TRIM(SCFMessage)//'Gap = '//TRIM(DblToShrtChar(-Gap))//" hartree"
  IF(Entropy > Zero) SCFMessage = TRIM(SCFMessage)//", Entropy = "//TRIM(DblToShrtChar(Entropy))//" hartree"
  IF(P%NSMat/=1.AND.Args%I%I(1)/=0)SCFMessage=TRIM(SCFMessage)//' <S^2> = '//TRIM(FltToShrtChar(S2))
  IF(TRIM(SCFMessage) /= "") &
  CALL MondoLog(DEBUG_MAXIMUM,Prog,SCFMessage,SCFTag)

  SCFMessage='<T> = '//TRIM(DblToMedmChar(KinE))//', <V> = '//TRIM(DblToMedmChar(E_es_tot))//" hartree"
  IF(ExchE/=Zero)SCFMessage=TRIM(SCFMessage)//', <HF> = '//TRIM(DblToMedmChar(ExchE))//" hartree"
  IF(Exc/=Zero)SCFMessage=TRIM(SCFMessage)//', <DFT> = '//TRIM(DblToMedmChar(Exc))//" hartree"
  CALL MondoLog(DEBUG_MAXIMUM,Prog,SCFMessage,SCFTag)

  IF(SCFActn=='BasisSetSwitch')THEN
     SCFMessage='Basis set switch, MxD = '//TRIM(DblToShrtChar(DMax))
  ELSE
     SCFMessage='<SCF> = '//TRIM(DblToChar(Etot))//' hartree, '//TRIM(DblToChar(Etot*au2eV))//' eV, dD = '//TRIM(DblToShrtChar(DMax))
  ENDIF

  ! Add in DIIS error
  IF(DIISErr/=Zero)SCFMessage=TRIM(SCFMessage)//', DIIS = '//TRIM(DblToShrtChar(DIISErr))
  CALL MondoLog(DEBUG_NONE,Prog,SCFMessage,SCFTag)

  CALL Delete(P)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL ShutDown(Prog)

CONTAINS

  REAL(DOUBLE) FUNCTION GetS2(P,S,Tmp1)
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR)    :: P,S,Tmp1
#else
    TYPE(BCSR)     :: P,S,Tmp1
#endif
    TYPE(BCSR)     :: Tmp2
    TYPE(DBL_RNK2) :: D1,D2
    INTEGER        :: I
    REAL(DOUBLE)   :: Sx,Sy,Sz

    CALL Initialize(Tmp2)

    GetS2=0D0
    CALL Multiply(P,S,Tmp1)
    CALL SetEq(Tmp2,Tmp1)
    IF(MyID.EQ.0)THEN
      CALL SetEq(D1,Tmp2)
      CALL New(D2,(/NBasF,NBasF/))
      SELECT CASE(P%NSMat)
      CASE(1)
        !
      CASE(2)
        CALL DGEMM('N','N',NBasf,NBasF,NBasF,1D0,D1%D(1,1),NBasF, &
             D1%D(1,NBasF+1),NBasF,0D0,D2%D(1,1),NBasF)
        DO I=1,NBasF
          GetS2=GetS2+D2%D(I,I)
        ENDDO
        GetS2=0.25D0*(NAlph-NBeta)**2+0.5D0*(NAlph+NBeta)-GetS2
      CASE(4)
        !
      CASE DEFAULT;CALL Halt('GetS2: Something is wrong there!')
      END SELECT
      CALL Delete(D1)
      CALL Delete(D2)
      CALL Delete(Tmp2)
    ENDIF
  END FUNCTION GetS2
END PROGRAM SCFStatus
