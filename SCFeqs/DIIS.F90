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

PROGRAM DIIS
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE SetXYZ
  USE LinAlg
  USE MatFunk
  USE MondoLogger
#ifdef PARALLEL
  USE MondoMPI
#endif

  IMPLICIT NONE

#ifdef PARALLEL
  TYPE(DBCSR) :: F,P,S,X,Z,ZT,EI,EJ,Tmp1,Tmp2
#else
  TYPE(BCSR)  :: F,P,S,X,Z,ZT,EI,EJ,Tmp1,Tmp2
#endif

  TYPE(ARGMT)                      :: Args
  TYPE(INT_VECT)                   :: IWork,Idx,SCFOff,DIISInfo
  TYPE(DBL_VECT)                   :: V,DIISCo,AbsDIISCo
  TYPE(DBL_RNK2)                   :: B,BB,BInv,BTmp
  TYPE(CMPoles),DIMENSION(2)       :: MP
  REAL(DOUBLE),DIMENSION(3)        :: AvDP,DeltaDPi,DeltaDPj
  REAL(DOUBLE),DIMENSION(6)        :: AvQP,DeltaQPi,DeltaQPj
  REAL(DOUBLE)                     :: DIISErr,C0,C1,Damp,EigThresh,RelDevDP,RelDevQP,Sellers,DMax
  INTEGER                          :: I,J,I0,J0,K,N,iSCF,BMax,DoDIIS,iDIIS,iOffSet,DIISBeg
  INTEGER                          :: DIISDelay, DIISFirstSCF
  REAL(DOUBLE)                     :: noise
  CHARACTER(LEN=2)                 :: Cycl,NxtC
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: Mssg,FFile
  LOGICAL                          :: Present,Sloshed
  INTEGER                          :: IPresent,JPresent
  CHARACTER(LEN=4),PARAMETER       :: Prog='DIIS'
  CHARACTER(LEN=DCL)               :: NAME
#if defined(PARALLEL_CLONES)
  INTEGER                          :: oldClone, rank
#endif

  ! Initial setup
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  iSCF=Args%I%I(1)

  ! Parse for DIIS options
  CALL OpenASCII(InpFile,Inp)

  ! Damping coefficient for first cycle
  !
  ! When we do not DIIS we construct the new Fock matrix by simple linear
  ! mixing. A Damp factor of 1 means that we only keep the old Fock matrix,
  ! i.e. no propagation and the SCF cylce stalls. A Damp factor of 0 means
  ! that we only use the new Fock matrix, i.e. full propagation.
  IF(.NOT.OptDblQ(Inp,'DIISDamp',Damp)) THEN
    Damp = 0.5D0
  ENDIF

  ! Max number of equations to keep in DIIS
  IF(.NOT.OptIntQ(Inp,'DIISDimension',BMax)) THEN
    BMax = 8
  ELSEIF(BMax == 0) THEN
    ! DIISDimension set to 0, which means linear mixing.
  ENDIF

  ! Set the delay, i.e. the first SCF cycle for which we will use DIIS.
  IF(.NOT.OptIntQ(Inp, "DIISDelay", DIISDelay)) THEN
    DIISDelay = 2
  ENDIF

  ! Set the first SCF cylce that will be considered for the B matrix.
  IF(.NOT.OptIntQ(Inp, "DIISFirstSCF", DIISFirstSCF)) THEN
    DIISFirstSCF = DIISDelay
  ENDIF

  IF(DIISDelay < 1) THEN
    CALL Warn("[DIIS] DIISDelay has be at least 1")
    DIISDelay = 1
  ENDIF

  IF(DIISDelay < DIISFirstSCF) THEN
    CALL Warn("[DIIS] DIISDelay has to be greater than DIISFirstSCF!")
    DIISFirstSCF = DIISDelay
  ENDIF

  IF(BMax.GT.DIIS_MAX_MATRIX_SIZE)THEN
    CALL Warn('Requested DIISDimension '//TRIM(IntToChar(BMax))// &
              ' greater than DIIS_MAX_MATRIX_SIZE '//TRIM(IntToChar(DIIS_MAX_MATRIX_SIZE))// &
              '! '//RTRN//'Reseting DIISDimension to DIIS_MAX_MATRIX_SIZE.')
    BMax=DIIS_MAX_MATRIX_SIZE
  ENDIF

  IF(.NOT.OptDblQ(Inp, "DIISRandomNoise", noise)) THEN
    noise = 0.0D0
  ELSE
    IF(noise > 1.0D0) THEN
      noise = 1.0D0
    ENDIF

    IF(noise < 0.0D0) THEN
      noise = 0.0D0
    ENDIF

    CALL MondoLog(DEBUG_NONE, Prog, "DIIS random noise = "//TRIM(FltToChar(noise)))
  ENDIF

  CLOSE(Inp)

  ! Allocations
  CALL New(P)
  CALL New(F)
  CALL New(S)
  CALL New(X)
  CALL New(Z)
  CALL New(ZT)
  CALL New(EI)
  CALL New(EJ)
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(MP(1))
  CALL New(MP(2))

#ifdef USE_AO
  ! The current DIIS error in atomic orbital basis. The DIIS error is given by
  !
  ! FPS - SPF
  !
  ! Journal of Computational Chemistry 3, 556, P. Pulay, Improved SCF
  ! Convergence Acceleration (1982).
  CALL Get(F,TrixFile('F',Args,0))
  CALL Get(P,TrixFile('D',Args,0))

  ! In HF, we might be off here by a factor of 2. This does not change the DIIS
  ! coefficients however, only the normalization factor, DIISCo%D(N).
  CALL Multiply(P,2.0D0)

  CALL Get(S,TrixFile("S",Args))

  INQUIRE(FILE=TrixFile("X", Args), EXIST=Present)
  IF(Present) THEN
    CALL Get(X, TrixFile("X", Args))
  ELSE
    CALL Get(Z, TrixFile("Z", Args))
    CALL Get(ZT, TrixFile("ZT", Args))
  ENDIF

  ! Calculate commutator.
  CALL Multiply(S,P,Tmp1)
  CALL Multiply(Tmp1,F,EI)
  CALL Multiply(F,P,Tmp1)
  CALL Multiply(Tmp1,S,EI,-One)

  CALL PPrint(F, "Fockian")
  CALL PPrint(P, "Density")
  CALL PPrint(S, "S")
  IF(Present) THEN
    CALL PPrint(X, "X")
  ELSE
    CALL PPrint(Z, "Z")
    CALL PPrint(ZT, "ZT")
  ENDIF
  CALL PPrint(EI, "DIIS Error Matrix")

  ! Transform error matrix into orthogonal basis.
  !
  ! e^ortho = Z^t e Z.
  IF(Present) THEN
    CALL Multiply(X, EI, Tmp1)
    CALL Multiply(Tmp1, X, EI)
  ELSE
    CALL Multiply(ZT, EI, Tmp1)
    CALL Multiply(Tmp1, Z, EI)
  ENDIF

#else

  ! The DIIS error using the Fockian and the density in orthogonal basis
  ! representation.
  CALL Get(F, TrixFile("OrthoF", Args, 0))
  CALL Get(P, TrixFile("OrthoD", Args, 0))

  ! Calculate commutator.
  CALL Multiply(F, P, EI)
  CALL Multiply(P, F, EI, -One)

#endif

  CALL Put(EI, TrixFile("OrthoE_DIIS", Args, 0))
  !CALL PPrint(EI, "DIIS Error Matrix in orthogonal basis")

  DIISErr=SQRT(Dot(EI,EI))/DBLE(NBasF)

  ! Consider just damping, certainly on first go through
  IF(iSCF < DIISDelay)THEN
    DoDIIS = 0   ! No DIIS, but damp non-extrapolated Fock matrices
  ELSEIF(BMax /= 0)THEN
    DoDIIS = 1   ! We are doing DIIS, extrapolating non-extrapolated Fock matrices
  ELSE
    DoDIIS = 0   ! We are purely damping, using previously extrapolated Fock matrices
  ENDIF

  ! Build the B matrix.
  IF(DoDIIS == 1)THEN

    N = MIN(iSCF-DIISFirstSCF+1, BMax)+1

    CALL New(B,(/N,N/))
    CALL New(DIISInfo,2)
    CALL Get(DIISInfo,'diisinfo')
    CALL New(BTmp,(/DIIS_MAX_MATRIX_SIZE,DIIS_MAX_MATRIX_SIZE/))
    CALL Get(BTmp,'diismtrix')

    IF(DIISInfo%I(1).EQ.Args%I%I(2).AND.DIISInfo%I(2).EQ.Args%I%I(3)) THEN
      ! We didn't BS switch or new geom or oda, then build a part of the B matrix.
      B%D(1:N-2,1:N-2)=BTmp%D(1:N-2,1:N-2)
      DIISFirstSCF=iSCF
    ELSE
      ! We did BS switch or new geom or oda, then build the full B matrix.
      DIISInfo%I(1)=Args%I%I(2)
      DIISInfo%I(2)=Args%I%I(3)
    ENDIF

#if defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
      CALL Put(DIISInfo, 'diisinfo')

      oldClone = MyClone
      DO rank = 1, MSize(MPI_COMM_WORLD)-1
        CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(DIISInfo, UBOUND(DIISInfo%I, 1), rank, PUT_TAG, M_O = LBOUND(DIISInfo%I, 1), comm_O = MPI_COMM_WORLD)

        ! Put to correct HDFGroup.
        CALL CloseHDFGroup(H5GroupID)
        H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID = H5GroupID
        CALL Put(DIISInfo, 'diisinfo')
      ENDDO
      MyClone = oldClone

      ! Reopen old HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
    ELSE
      CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(DIISInfo, UBOUND(DIISInfo%I, 1), ROOT, PUT_TAG, M_O = LBOUND(DIISInfo%I, 1), comm_O = MPI_COMM_WORLD)
    ENDIF
#else
    CALL Put(DIISInfo,'diisinfo')
#endif
    CALL Delete(DIISInfo)

    DO I = DIISFirstSCF, iSCF
      FFile=TrixFile('OrthoE_DIIS',Args,I-iSCF)
      INQUIRE(FILE=FFile,EXIST=Present)
      IF(Present) THEN
        IPresent = 0
      ELSE
        IPresent = 1
      ENDIF

      IF(IPresent /= 0) THEN

#ifdef USE_AO
        ! Construct missing E_DIIS.
        CALL Get(F,TrixFile('F',Args,I-iSCF))
        CALL Get(P,TrixFile('D',Args,I-iSCF))
        CALL Multiply(P,2.0D0)

        INQUIRE(FILE=TrixFile("X", Args), EXIST=Present)
        IF(Present) THEN
          CALL Get(X, TrixFile("X", Args))
        ELSE
          CALL Get(Z, TrixFile("Z", Args))
          CALL Get(ZT, TrixFile("ZT", Args))
        ENDIF

        CALL Multiply(S,P,Tmp1)
        CALL Multiply(Tmp1,F,EI)
        CALL Multiply(F,P,Tmp1)
        CALL Multiply(Tmp1,S,EI,-One)

        IF(Present) THEN
          CALL Multiply(X, EI, Tmp1)
          CALL Multiply(Tmp1, X, EI)
        ELSE
          CALL Multiply(ZT, EI, Tmp1)
          CALL Multiply(Tmp1, Z, EI)
        ENDIF
#else
        ! Construct missing E_DIIS.
        CALL Get(F, TrixFile("OrthoF", Args, 0))
        CALL Get(P, TrixFile("OrthoD", Args, 0))

        CALL Multiply(F, P, EI)
        CALL Multiply(P, F, EI, -One)
#endif

        CALL Put(EI, TrixFile("OrthoE_DIIS", Args, I-iSCF))

      ELSE
        CALL Get(EI,TrixFile('OrthoE_DIIS',Args,I-iSCF))
      ENDIF

      ! We dont filter E for obvious reasons
      DO J = iSCF-N+2, I
        FFile=TrixFile('OrthoE_DIIS',Args,J-iSCF)
        INQUIRE(FILE=FFile,EXIST=Present)
        JPresent=1
        IF(Present) JPresent=0

        IF(JPresent.EQ.0) then
          CALL Get(EJ,TrixFile('OrthoE_DIIS',Args,J-iSCF))
        ELSE
          CALL MondoLog(DEBUG_NONE, Prog, "no previous information on error matrix E_DIIS["//TRIM(IntToChar(J))//"] found, constructing it")
          CALL Halt("he?")
        ENDIF

        IF(iSCF-DIISFirstSCF > 0) THEN
          I0 = I-DIISFirstSCF+1
        ELSE
          I0 = N-1
        ENDIF
        J0 = J-(iSCF-N+2)+1

        ! B_ij = Trace(e_i e_j^{dagger})
        B%D(I0, J0)=Dot(EI,EJ)
        B%D(J0, I0)=B%D(I0, J0)

      ENDDO
    ENDDO

    B%D(N,1:N)=One
    B%D(1:N,N)=One
    B%D(N,N)=Zero

    IF(N.LT.BMax+1) THEN
      ! We didn't reach the size of the matrix.
      BTmp%D(1:N-1,1:N-1)=B%D(1:N-1,1:N-1)
    ELSE
      ! We reach the size of the matrix, we reduce it.
      BTmp%D(1:N-2,1:N-2)=B%D(2:N-1,2:N-1)
    ENDIF
#if defined(PARALLEL_CLONES)
    IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
      CALL Put(BTmp, 'diismtrix')

      oldClone = MyClone
      DO rank = 1, MSize(MPI_COMM_WORLD)-1
        CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
        CALL Recv(BTmp, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

        ! Put to correct HDFGroup.
        CALL CloseHDFGroup(H5GroupID)
        H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID = H5GroupID
        CALL Put(BTmp, 'diismtrix')
      ENDDO
      MyClone = oldClone

      ! Reopen old HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
    ELSE
      CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Send(BTmp, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    ENDIF
#else
    CALL Put(BTmp, 'diismtrix')
#endif
    CALL Delete(BTmp)

    CALL PPrint(B, "DIIS B matrix")

    ! Solve the least squares problem to obtain new DIIS coeficients.
    CALL New(DIISCo,N)
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      CALL New(BInv,(/N,N/))
      CALL New(V,N)
      V%D=Zero
      V%D(N)=One
      CALL SetDSYEVWork(N)
      BInv%D=Zero
      IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
        CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE.,EigenThresh_O=1D-12, &
             PrintValues_O=.TRUE.,PrintCond_O=.TRUE.,Prog_O=Prog)
      ELSE
        CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE.,EigenThresh_O=1D-12)
      ENDIF
      CALL UnSetDSYEVWork()
      CALL DGEMV('N',N,N,One,BInv%D,N,V%D,1,Zero,DIISCo%D,1)
      CALL Delete(V)
      CALL Delete(BInv)
#ifdef PARALLEL
    ENDIF
    CALL BCast(DIISCo)
#endif

  ELSEIF(iSCF > 0) THEN

    N = 3
    CALL New(DIISCo,2)

    ! Damping on the second cycle
    DIISCo%D(1)=Damp
    DIISCo%D(2)=One-Damp

    ! Add some random noise to the coefficients.
    DIISCo%D(1) = ABS(DIISCo%D(1)+Random((/0.0D0, noise/)))
    DIISCo%D(2) = One-DIISCo%D(1)

  ELSE

    N = 3
    CALL New(DIISCo, 2)

    DIISCo%D(1) = Zero
    DIISCo%D(2) = One

  ENDIF

  ! IO
#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL Put(DIISErr, 'diiserr')

    oldClone = MyClone
    DO rank = 1, MSize(MPI_COMM_WORLD)-1
      CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(DIISErr, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

      ! Put to correct HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
      CALL Put(DIISErr, 'diiserr')
    ENDDO
    MyClone = oldClone

    ! Reopen old HDFGroup.
    CALL CloseHDFGroup(H5GroupID)
    H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID = H5GroupID
  ELSE
    CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(DIISErr, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
  ENDIF
#else
  CALL Put(DIISErr, 'diiserr')
#endif

  ! Allocate some indecies for re-ordering
  CALL New(Idx,N)
  CALL New(SCFOff,N)
  CALL New(AbsDIISCo,N)

  ! Reorder the DIIS, starting with smallest values and summing to the largest
  DO I=1,N-1
    Idx%I(I)=I
    IF(DoDIIS == 1) THEN
      ! DIIS
      SCFOff%I(I) = I-(N-1)
    ELSE
      ! Linear mixing...
      SCFOff%I(I) = I-2
    ENDIF
    AbsDIISCo%D(I)=ABS(DIISCo%D(I))
  ENDDO

  CALL Sort(AbsDIISCo,Idx,N-1,1)

  ! Start with a matrix of diagonal zeros...
  CALL SetToI(F)
  CALL Multiply(F,Zero)

  ! And do the summation
  DO I=1,N-1
    CALL MondoLog(DEBUG_MAXIMUM, Prog, "Adding "//TRIM(DblToChar(DIISCo%D(Idx%I(I))))//" * F["//TRIM(IntToChar(iSCF+SCFOff%I(Idx%I(I))))//"]")
    IF(DIISCo%D(Idx%I(I)) /= Zero) THEN
#ifdef SUM_AO
      CALL Get(Tmp1,TrixFile('F',Args,SCFOff%I(Idx%I(I))))
#else
      CALL Get(Tmp1,TrixFile('OrthoF',Args,SCFOff%I(Idx%I(I))))
#endif
      CALL Multiply(Tmp1,DIISCo%D(Idx%I(I)))
      CALL Add(F,Tmp1,EI)
      CALL SetEq(F,EI)
    ENDIF
  ENDDO

  Mssg=""
  DO I=1,N-2
     IF(MOD(I,6)==0)THEN
        IF(I.NE.1) CALL MondoLog(DEBUG_MEDIUM, Prog, Mssg, 'C-1 coefficients')
        Mssg=TRIM(DblToMedmChar(DIISCo%D(Idx%I(I))))//','
     ELSE
        Mssg=TRIM(Mssg)//'  '//TRIM(DblToMedmChar(DIISCo%D(Idx%I(I))))//','
     ENDIF
  ENDDO
  Mssg=TRIM(Mssg)//'  '//TRIM(DblToMedmChar(DIISCo%D(Idx%I(N-1))))
  CALL MondoLog(DEBUG_MEDIUM, Prog, Mssg, 'C-1 coefficients')

#ifdef SUM_AO
  ! F^ortho = Z^t F Z.
  IF(Present) THEN
    CALL Multiply(X, F, Tmp1)
    CALL Multiply(Tmp1, X, F)
  ELSE
    CALL Multiply(ZT, F, Tmp1)
    CALL Multiply(Tmp1, Z, F)
  ENDIF
#endif

  ! IO for the orthogonal, extrapolated F
  CALL Put(F,TrixFile('F_DIIS',Args,0))
  CALL PChkSum(F,'F_DIIS['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint(F,'F_DIIS['//TRIM(SCFCycl)//'] in orthogonal basis')
  CALL Plot(F,'F_DIIS_'//TRIM(SCFCycl))

  ! Tidy up
  CALL Delete(Idx)
  CALL Delete(SCFOff)
  CALL Delete(AbsDIISCo)
  CALL Delete(F)
  CALL Delete(DIISCo)
  CALL Delete(P)
  CALL Delete(S)
  CALL Delete(X)
  CALL Delete(Z)
  CALL Delete(ZT)
  CALL Delete(EI)
  CALL Delete(EJ)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL ShutDown(Prog)
END PROGRAM DIIS
