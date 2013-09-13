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

#include "MondoConfig.h"

PROGRAM ODA
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE DenMatMethods
  IMPLICIT NONE
  !-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                    ::  &
#else
  TYPE(BCSR)                     ::  &
#endif
  P,PTilde,F,FTilde,T,K0,K1,T1,T2,T3
  REAL(DOUBLE)                   :: e0,e1,e0p,e1p,a3,b3,c3,d3,EMns,EPls,EMin,        &
       LMns,LPls,L,L1,ENucTotTilde,                     &
       DIISErr,Enuc0,Enuc1,Exc0,Exc1
  REAL(DOUBLE)                   :: Tmp1,Tmp2,Tmp3,Tmp4,alph,SFac
  INTEGER                        :: I,iSCF
  LOGICAL                        :: Present,HasECPs
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: ODAMssg,ODATag,MatFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='ODA'
  REAL(DOUBLE)                   :: TrP0T,TrP1T,TrP0F0,TrP1F1,TrP0F1,TrP1F0, &
       TrP0K1,TrP0K0,TrP1K0,TrP1K1,IKPS_Error,IKPS_denom
#if defined(PARALLEL_CLONES)
  INTEGER                          :: oldClone, rank
#endif
  !-------------------------------------------------------------------------
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  iSCF=Args%I%I(1)
  ! Suss for matrix threshold overide
  CALL SussTrix('ODATrix',Prog)
  ! Allocate All Read from Disk Matrices, what a waste
  CALL New(P)
  CALL New(PTilde)
  CALL New(F)
  CALL New(FTilde)
  CALL New(K0)
  CALL New(K1)
  ! Allocate temp Matrices
  CALL New(T1)
  CALL New(T2)
  CALL New(T3)
  ! Get the Matrices
  CALL Get(PTilde,TrixFile('D',Args,-1))
  CALL Get(FTilde,TrixFile('F',Args,-1))
  CALL Get(P,TrixFile('D',Args,0))
  !---------------------------------------------
  ! Rescaling factor for R/U/G theory.
  SFac=1D0
  IF(P%NSMat.GT.1) SFac=0.5D0
  !WRITE(*,*) "ODA, NSMat = ", P%NSMat, ", Sfac = ", SFac
  !---------------------------------------------
  CALL Get(F,TrixFile('F',Args,0))
  ! Get Kinectic Energy and ECPs
  CALL Get(T,TrixFile('T',Args))
  CALL Get(HasECPs,'hasecps',Tag_O=CurBase)
  IF(HasECPs) THEN
    CALL Get(T1,TrixFile('U',Args))
    CALL Add(T,T1,T2)
    CALL SetEq(T,T2)
  ENDIF
  ! Calculate Exchange Asymmetry
  IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
    IF(HasHF(ModelChem)) THEN
      CALL Get(K0,TrixFile('K',Args,-1))
      CALL Get(K1,TrixFile('K',Args,0))
#ifdef PARALLEL
      CALL Multiply(PTilde,K1,T1)
      TrP0K1 =  Trace(T1)
      CALL Multiply(P,K0,T1)
      TrP1K0 =  Trace(T1)
      CALL SetEq(T1,K0)
      CALL Multiply(T1,-One)
      CALL Add(K1,T1,T2)
      CALL Multiply(T2,T2,T3)
      IKPS_denom = SQRT(ABS(Trace(T3)))
      IKPS_Error =  ABS(TrP0K1-TrP1K0)/IKPS_denom
#else
      TrP0K1 =  Trace(PTilde,K1)
      TrP1K0 =  Trace(P,K0)
      CALL SetEq(T1,K0)
      CALL Multiply(T1,-One)
      CALL Add(K1,T1,T2)
      IKPS_denom = SQRT(ABS(Trace(T2,T2)))
      IKPS_Error =  ABS(TrP0K1-TrP1K0)/IKPS_denom
#endif
      CALL MondoLog(DEBUG_NONE, "ODA", "IKPS_Error = "//TRIM(DblToChar(IKPS_Error)))
    ENDIF
  ENDIF
  ! Get Energies: E_nuc and E_xc and K_xc matrices
  Current(1)=Current(1)-1
  CALL Get(Enuc0,'E_NuclearTotal',Stats_O=Current)
  Current(1)=Current(1)+1
  CALL Get(Enuc1,'E_NuclearTotal',Stats_O=Current)
  IF(HasDFT(ModelChem)) THEN
    Current(1)=Current(1)-1
    CALL Get(Exc0 ,'Exc'           ,Stats_O=Current)
    CALL Get(K0,TrixFile('Kxc',Args,-1))
    Current(1)=Current(1)+1
    CALL Get(Exc1 ,'Exc'           ,Stats_O=Current)
    CALL Get(K1,TrixFile('Kxc',Args,0))
  ENDIF
  ! Compute the Endpoint energies and Derivatives
#ifdef PARALLEL
  ! Compute the Traces
  CALL Multiply(PTilde,T,T1)
  TrP0T = Trace(T1)*SFac
  CALL Multiply(P,T,T1)
  TrP1T = Trace(T1)*SFac
  CALL Multiply(PTilde,FTilde,T1)
  TrP0F0 = Trace(T1)*SFac
  CALL Multiply(P,FTilde,T1)
  TrP1F0 = Trace(T1)*SFac
  CALL Multiply(PTilde,F,T1)
  TrP0F1 = Trace(T1)*SFac
  CALL Multiply(P,F,T1)
  TrP1F1 = Trace(T1)*SFac
  IF(HasDFT(ModelChem)) THEN
    CALL Multiply(PTilde,K0,T1)
    TrP0K0 = Trace(T1)
    CALL Multiply(P,K1,T1)
    TrP1K1 = Trace(T1)
    CALL Multiply(PTilde,K1,T1)
    TrP0K1 = Trace(T1)
    CALL Multiply(P,K0,T1)
    TrP1K0 = Trace(T1)
  ENDIF
  e0  = TrP0T+TrP0F0+Enuc0
  e1  = TrP1T+TrP1F1+Enuc1
  e0p = Enuc1-Enuc0+TrP1T-TrP0T+TrP1F0+TrP0F1-Two*TrP0F0
  e1p = Enuc1-Enuc0+TrP1T-TrP0T+Two*TrP1F1-TrP1F0-TrP0F1
  IF(HasDFT(ModelChem)) THEN
    e0  = e0  + Exc0 - TrP0K0
    e1  = e1  + Exc1 - TrP1K1
    e0p = e0p + TrP1K0-TrP0K1
    e1p = e1p + TrP1K0-TrP0K1
  ENDIF
#else
  !!! HEY, WE ARE OVERCOMPUTING HERE!!  CLEAN UP THE REDUNDANT TRACE ETC!!!
  ! Avoid assumption of two electron integral symmetry which may be lost
  ! in the case of small cell PBC HF and also due to excesive thresholding.
  e0  = Trace(PTilde,T)*SFac+Trace(PTilde,FTilde)*SFac + Enuc0
!  CALL PPrint(P," P1 ",Unit_O=6)
  e1  = Trace(P,T)*SFac     +Trace(P,F)*SFac           + Enuc1
  e0p = Enuc1-Enuc0+Trace(P,T)*SFac-Trace(PTilde,T)*SFac+Trace(P,FTilde)*SFac+Trace(PTilde,F)*SFac-Two*Trace(PTilde,FTilde)*SFac
  e1p = Enuc1-Enuc0+Trace(P,T)*SFac-Trace(PTilde,T)*SFac+Two*Trace(P,F)*SFac-Trace(P,FTilde)*SFac-Trace(PTilde,F)*SFac
  IF(HasDFT(ModelChem)) THEN
    e0  = e0  + Exc0 - Trace(PTilde,K0)*SFac
    e1  = e1  + Exc1 - Trace(P,K1)*SFac
    e0p = e0p + (Trace(P,K0)-Trace(PTilde,K1))*SFac
    e1p = e1p + (Trace(P,K0)-Trace(PTilde,K1))*SFac
  ENDIF
#endif
  ! Find the mixing parameter L from the
  ! cubic E3(L)=a3+b3*L+c3*L^2+d3*L^3
  !
!!$    WRITE(*,*) "e0  = ",e0
!!$    WRITE(*,*) "e1  = ",e1
!!$    WRITE(*,*) "e0p = ",e0p
!!$    WRITE(*,*) "e1p = ",e1p
!!$    WRITE(*,*) " "
  !
  a3=e0
  b3=e0p
  c3=-3D0*e0-2D0*e0p+3D0*e1-e1p
  d3=2D0*e0+e0p-2D0*e1+e1p
  !
!!$    WRITE(*,*) "a3 = ",a3
!!$    WRITE(*,*) "b3 = ",b3
!!$    WRITE(*,*) "c3 = ",c3
!!$    WRITE(*,*) "d3 = ",d3
  !
  IF(ABS(d3)<1D-6.OR.c3*c3-3*b3*d3<Zero)THEN
    L=-Half*b3/c3
    L1=One-L
    EMin=a3+b3*L+c3*L**2+d3*L**3
  ELSE
    IF(c3*c3 < 3*b3*d3) THEN
      CALL Halt("ODA: c3*c3-3*b3*d3 = "//TRIM(FltToShrtChar(c3*c3-3*b3*d3)))
    ENDIF

    LMns=MAX(Zero,MIN(One,(-c3-SQRT(c3*c3-3*b3*d3))/(3*d3)))
    LPls=MAX(Zero,MIN(One,(-c3+SQRT(c3*c3-3*b3*d3))/(3*d3)))
    EMns=a3+b3*LMns+c3*LMns**2+d3*LMns**3
    EPls=a3+b3*LPls+c3*LPls**2+d3*LPls**3
    IF(EMns<EPls)THEN
      L=LMns
      L1=One-LMns
    ELSE
      L=LPls
      L1=One-LPls
    ENDIF
    EMin=MIN(EMns,EPls)
  ENDIF
  ! End point and Midpoint checks
  IF((L<=Zero.OR.L>One) .OR. (EMin>e0 .AND. EMin > e1))THEN
    IF(e0<e1)THEN
      L=1D-2
      L1=One-L
      EMin=a3+b3*L+c3*L**2
    ELSE
      L=One
      L1=Zero
      EMin=e1
    ENDIF
  ENDIF
  ! Parse for Mixing OverRide
  CALL OpenASCII(InpFile,Inp)
  IF(OptDblQ(Inp,'ODALambda',LMns)) THEN
    L   = LMns
    L1  = One-L
    EMin=a3+b3*L+c3*L**2
  ENDIF
  CLOSE(Inp)
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
    ODATag='Mix = '//TRIM(FltToShrtChar(L))
    ODAMssg="<SCF> = "//TRIM(DblToChar(EMin))//" hartree, "//TRIM(DblToChar(EMin*au2eV))//' eV, 3 = '//TRIM(DblToShrtChar(d3))
    CALL MondoLog(DEBUG_MEDIUM,Prog,ODAMssg,ODATag)
#ifdef PARALLEL
  ENDIF
#endif
  ! Compute PTilde_N=(1-L)*PTilde_(N-1)+L*P_N, then put bto disk
  CALL Multiply(P,L)
  CALL Multiply(PTilde,L1)
  CALL Add(P,PTilde,T2)
  CALL Put(T2,TrixFile('D',Args,0))
  !CALL Put(T2,'CurrentDM',CheckPoint_O=.TRUE.)
  CALL PChkSum(T2,'Pao['//TRIM(NxtCycl)//']',Prog)
  ! Compute FTilde_N ~ (1-L)*FTilde_(N-1)+L*F_N
  CALL Multiply(F,L)
  CALL Multiply(FTilde,L1)
  CALL Add(F,FTilde,T3)
  CALL Put(T3,TrixFile('F',Args,0))
  CALL PChkSum(T3,'Fao['//TRIM(NxtCycl)//']',Prog)
  !----------------------------------------------------------------------
  ! Convert FTilde to ortho rep
  INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
  IF(Present)THEN
    CALL Get(K0,TrixFile('X',Args))         ! Z=S^(-1/2)
    !    F_ortho
    CALL Multiply(K0,T3,T1)                 ! Z*F
    CALL Multiply(T1,K0,T3)                 ! (Z*F)*Z
    CALL Filter(T1,T3)                      ! Filter
    CALL Put(T1,TrixFile('OrthoF',Args,0))
    CALL PChkSum(T1,'For['//TRIM(NxtCycl)//']',Prog)
    CALL SetEq(FTilde,T1)
  ELSE
    CALL Get(K0,TrixFile('ZT',Args))        ! ZT=S^(L)
    CALL Get(K1,TrixFile('Z',Args))         ! Z=S^(-L)
    !    F_ortho
    CALL Multiply(K0,T3,T1)                 ! ZT*F
    CALL Multiply(T1,K1,T3)                 ! (ZT*F)*Z
    CALL Filter(T1,T3)                      ! Filter
    CALL Put(T1,TrixFile('OrthoF',Args,0))
    CALL PChkSum(T1,'For['//TRIM(NxtCycl)//']',Prog)
    CALL SetEq(FTilde,T1)
  ENDIF
  ! Convert PTilde to ortho rep
  IF(iSCF>1)THEN
    CALL Get(PTilde,TrixFile('OrthoD',Args,-1))
    CALL Get(P,TrixFile('OrthoD',Args,0))
    CALL Multiply(P,L)
    CALL Multiply(PTilde,L1)
    CALL Add(P,PTilde,T1)
    CALL PChkSum(T1,'Por['//TRIM(NxtCycl)//']',Prog)
    CALL Put(T1,TrixFile('OrthoD',Args,0))
    CALL SetEq(PTilde,T1)
    !    Compute the DIIS error
    CALL Multiply(FTilde,PTilde,T1)
    CALL Multiply(PTilde,FTilde,T2)
    CALL Multiply(T2,-One)
    CALL Add(T1,T2,T3)
    DIISErr=SQRT(Dot(T3,T3))/DBLE(NBasF)
  ELSE
    DIISErr=One
  ENDIF


  ! JTilde_N = (1-L)*JTilde_(N-1)+L*J_N
  CALL Get(PTilde,TrixFile('J',Args,-1))
  CALL Get(P,     TrixFile('J',Args,0))
  CALL Multiply(P,L)
  CALL Multiply(PTilde,L1)
  CALL Add(P,PTilde,T1)
  CALL Put(T1,TrixFile('J',Args,0))
  ! KTilde_N = (1-L)*KTilde_(N-1)+L*K_N
  IF(HasHF(ModelChem))THEN
    CALL Get(PTilde,TrixFile('K',Args,-1))
    CALL Get(P     ,TrixFile('K',Args,0))
    CALL Multiply(P,L)
    CALL Multiply(PTilde,L1)
    CALL Add(P,PTilde,T1)
    CALL Put(T1,TrixFile('K',Args,0))
  ENDIF
  ! ENucTotTilde_N = (1-L)*ENucTotTilde_(N-1)+L*ENucTotTilde_N
  ENucTotTilde=L*Enuc1+L1*Enuc0

#if defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
       CALL Put(DIISErr,'diiserr')
       CALL Put(EMin   ,'ODAEnergy')
       CALL Put(ENucTotTilde,'E_NuclearTotal',Stats_O=Current)
       oldClone = MyClone
       DO rank = 1, MSize(MPI_COMM_WORLD)-1
          CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
          CALL Recv(DIISErr, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
          CALL Recv(EMin   , rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
          CALL Recv(ENucTotTilde, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
          ! Put to correct HDFGroup.
          CALL CloseHDFGroup(H5GroupID)
          H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
          HDF_CurrentID = H5GroupID
          CALL Put(DIISErr, 'diiserr')
          CALL Put(EMin   ,'ODAEnergy')
          CALL Put(ENucTotTilde,'E_NuclearTotal',Stats_O=Current)
       ENDDO
       MyClone = oldClone
       ! Reopen old HDFGroup.
       CALL CloseHDFGroup(H5GroupID)
       H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
       HDF_CurrentID = H5GroupID
    ELSE
      CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(DIISErr, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(EMin,    ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(ENucTotTilde,ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    ENDIF
#else
    CALL Put(DIISErr,'diiserr')
    CALL Put(EMin   ,'ODAEnergy')
    CALL Put(ENucTotTilde,'E_NuclearTotal',Stats_O=Current)
#endif

  ! Tidy up
  CALL Delete(P)
  CALL Delete(PTilde)
  CALL Delete(F)
  CALL Delete(FTilde)
  ! Delete Temps
  CALL Delete(T1)
  CALL Delete(T2)
  CALL Delete(T3)
!  STOP
  !
  CALL ShutDown(Prog)
  !
END PROGRAM ODA
