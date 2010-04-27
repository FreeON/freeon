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

PROGRAM RHEqs
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
  USE DenMatMethods,ONLY:MDiag_DSYEVX,MDiag_DSYEVD
  USE MondoLogger
#ifdef NAG
  USE F90_UNIX_ENV
#endif

  IMPLICIT NONE

  TYPE(BCSR)                     :: sP,sF,FT,sX,sTmp1,sTmp2
  TYPE(DBL_RNK2)                 :: X,F,MO,P
  TYPE(DBL_VECT)                 :: EigenV
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: CJK,HOMO,LUMO,dt
  REAL(DOUBLE)                   :: Mu,Entropy,Ek,Z,A1,A2,H1,H3,H4,Fk,Sigma,Dum
  REAL(DOUBLE)                   :: kB,Temperature,Occ,OccErr,EPS,dOccdMu,Ne,Beta,Xk
  INTEGER                        :: I,J,K,LgN,NRow,NCol,m
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FMatrix,PMatrix,XFile,smearing
  CHARACTER(LEN=5),PARAMETER     :: Prog='RHEqs'
  LOGICAL                        :: Present,DensityArchive
  EXTERNAL :: DGEMM_NT
#if defined(PARALLEL_CLONES)
  INTEGER                        :: oldClone, rank
  REAL(DOUBLE)                   :: HOMOLUMO
#endif

  CALL StartUp(Args,Prog,Serial_O=.TRUE.)

  !--------------------------------------------------------------------
  ! Initialize and parse some variables.
  !
  CALL OpenASCII(InpFile,Inp)
  Smearing='NoSmearing'
  IF(OptKeyQ(Inp,'Smearing','MP')) Smearing='Methfessel-Paxton'
  IF(OptKeyQ(Inp,'Smearing','Fermi-Dirac')) Smearing='Fermi-Dirac'

  IF(Smearing == "Fermi-Dirac") THEN
    IF(.NOT.OptDblQ(Inp,'SmearingValue',Sigma)) Sigma=0.002D0
    IF(.NOT.OptDblQ(Inp, "SmearingTemperature", Temperature)) THEN
      CALL Halt("Set SmearingTemperature > 0")
    ENDIF
  ENDIF
  CLOSE(Inp)
  !--------------------------------------------------------------------
  !
  !
  CALL New(sF)
  FMatrix=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FMatrix,EXIST=Present)
  IF(Present)THEN
    CALL Get(sF,FMatrix)
  ELSE
    CALL Get(sF,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF

  ! We get NSMat from hdf via StartUp.
  !NSMat=sF%NSMat
  IF(NSMat.GT.1.AND.Smearing.NE.'NoSmearing') CALL Halt('Smearing with unrestricted are not supported!')

  SELECT CASE(NSMat)
  CASE(1)
    NRow=  NBasF
    NCol=  NBasF
  CASE(2)
    NRow=  NBasF
    NCol=2*NBasF
  CASE(4)
    NRow=2*NBasF
    NCol=2*NBasF
  CASE DEFAULT
    CALL Halt(' RHeqs: sF%NSMat doesn''t have an expected value! ')
  END SELECT

!!$  ! Some symmetry checks.
!!$  CALL New(FT)
!!$
!!$  CALL Xpose(sF, FT)
!!$  FT%MTrix%D = -FT%MTrix%D
!!$  CALL Add(sF, FT, sTmp1)
!!$  CALL MondoLog(DEBUG_MINIMUM, "RHeqs", "FNorm(F-FT) = "//TRIM(DblToChar(FNorm(sTmp1))))
!!$  CALL Delete(FT)

  CALL New(F,(/NRow,NCol/))
  CALL SetEq(F,sF)
  CALL Delete(sF)

  !write(*,*) 'RHeqs: NSMat',NSMat
  !write(*,*) 'RHeqs: NAlph',NAlph
  !write(*,*) 'RHeqs: NBeta',NBeta

  CALL New(EigenV,NCol)
  CALL SetEq(EigenV,Zero)

  SELECT CASE(NSMat)
  CASE(1)
    ! We just have one matrix.
    IF(Smearing.EQ.'NoSmearing')THEN
      CALL MDiag_DSYEVX(F,NBasF,MIN(NBasF,Nel/2+1),EigenV,0)
    ELSE
      CALL MDiag_DSYEVD(F,NBasF,EigenV,0)
    ENDIF
    !
    HOMO=EigenV%D(NEl/2)
    LUMO=EigenV%D(MIN(NBasF,Nel/2+1))
  CASE(2)
    ! We have a block diagonal matrix, we diag each blocks separately.
    IF(Smearing.EQ.'NoSmearing')THEN
      CALL MDiag_DSYEVX(F,NBasF,MIN(NBasF,NAlph+1),EigenV,0)
      CALL MDiag_DSYEVX(F,NBasF,MIN(NBasF,NBeta+1),EigenV,NBasF)
    ELSE
      CALL MDiag_DSYEVD(F,NBasF,EigenV,0)
      CALL MDiag_DSYEVD(F,NBasF,EigenV,NBasF)
    ENDIF
    !
    IF(EigenV%D(NAlph)-EigenV%D(NAlph+1).LT.EigenV%D(NBasF+NBeta)-EigenV%D(NBasF+NBeta+1)) THEN
      HOMO=EigenV%D(NAlph  )
      LUMO=EigenV%D(MIN(NBasF,NAlph+1))
    ELSE
      HOMO=EigenV%D(NBasF+NBeta  )
      LUMO=EigenV%D(NBasF+MIN(NBasF,NBeta+1))
    ENDIF
  CASE(4)
    ! We just have one matrix.
    IF(Smearing.EQ.'NoSmearing')THEN
      CALL MDiag_DSYEVX(F,2*NBasF,MIN(NBasF,Nel+1),EigenV,0)
    ELSE
      CALL MDiag_DSYEVD(F,2*NBasF,EigenV,0)
    ENDIF
    !
    HOMO=EigenV%D(NEl)
    LUMO=EigenV%D(MIN(NBasF,Nel+1))
  CASE DEFAULT
    CALL Halt(' RHeqs: NSMat doesn''t have an expected value! ')
  END SELECT
  !
  Mu=(HOMO+LUMO)*0.5D0
  CALL MondoLog(DEBUG_MEDIUM, Prog, "Chemical potential Mu = "//TRIM(DblToChar(Mu)))

  CALL MondoLog(DEBUG_MEDIUM, Prog, 'HOMO = '//TRIM(DblToMedmChar(HOMO))//', LUMO = '//TRIM(DblToMedmChar(LUMO)))

  !--------------------------------------------------------------
  ! Make a new closed shell, orthogonal density matrix
  !
  CALL New(P,(/NRow,NCol/))
  CALL DBL_VECT_EQ_DBL_SCLR(NRow*NCol,P%D(1,1),Zero)

  SELECT CASE(Smearing)
  CASE('Methfessel-Paxton')
    Entropy=0.D0
    DO K=1,NBasF
      ! PRB 40, 3616, 1989.
      ! Second order smearing N=2.
      Ek=EigenV%D(K)
      Z=(Ek-Mu)/Sigma
      A1=-1D0/( 4D0*SqrtPi)
      A2= 1D0/(32D0*SqrtPi)
      H1=2D0*Z
      H3=Z*(8D0*Z**2-12D0)
      H4=Z**2*(16D0*Z**2-48D0)+12D0
      ! Fractional occupation.
      Fk=0.5D0*(1D0-ERF(Z))+EXP(-Z**2)*(A1*H1+A2*H3)
      ! Entropic correction to the energy (Comp. Mat. Sci. 6, 15, 1996).
      Entropy=Entropy+Sigma*0.5D0*A2*H4*EXP(-Z**2)
      !write(*,*) Fk,Entropy
      DO J=1,NBasF
        CJK=F%D(J,K)*Fk
        DO I=1,NBasF
          P%D(I,J)=P%D(I,J)+F%D(I,K)*CJK
        ENDDO
      ENDDO
    ENDDO
    CALL MondoLog(DEBUG_MEDIUM, Prog, "Sigma = "//TRIM(DblToShrtChar(Sigma))// &
      ", Entropic correction per atom = "//TRIM(DblToShrtChar(Entropy/DBLE(NAtoms))))
  CASE('Fermi-Dirac')
    ! kB = 6.33366256E-06 ! Ry/K
    kB = 2.D0*6.33366256E-06 ! au/K
    ! kB = 8.61739E-05  ! eV/K
    CALL MondoLog(DEBUG_MAXIMUM, Prog, "kB = "//TRIM(DblToChar(kB)))

    ! Make this in input parameter.
    CALL MondoLog(DEBUG_MAXIMUM, Prog, "smearing temperature = "//TRIM(FltToChar(Temperature))//" K")

    ! Leave for now as hack.
    m = 6

    Beta = 1.D0/(kB*Temperature)
    Ne=Half*DBLE(NEl)
    EPS = 1E-10
    OccErr = 1.0D0
    DO WHILE (OccErr.GT.EPS)
      Entropy=0.D0
      Occ = 0.D0
      dOccdMu = 0.D0
      DO K=1,NBasF
        Ek=EigenV%D(K)

        Fk = 1.D0/(EXP(Beta*(Ek-Mu))+1.D0)

!        Xk = 0.5D0 - Beta*(Ek - Mu)/(2**(2+m))
!        DO I = 1,m
!          Xk = Xk*Xk/(2*Xk*(Xk-1.D0)+1.D0)
!        ENDDO
!        Fk = Xk

        Occ = Occ + Fk
        IF(Fk.GT.0.D0) THEN
          IF(Fk.LT.1.D0) THEN
            Entropy = Entropy - 2*Temperature*kB*(Fk*LOG(Fk)+(1.D0-Fk)*LOG(1.D0-Fk))
          ENDIF
        ENDIF
        dOccdMu = dOccdMu + Beta*Fk*(1.D0-Fk)
      ENDDO
      CALL MondoLog(DEBUG_MEDIUM, Prog, "Mu_fore = "//TRIM(DblToChar(Mu)))
      Mu = Mu + (Ne - Occ)/dOccdMu
      OccErr = ABS(Occ-Ne)
      CALL MondoLog(DEBUG_MEDIUM, Prog, "Mu = "//TRIM(DblToChar(Mu)))
      CALL MondoLog(DEBUG_MEDIUM, Prog, "Occ = "//TRIM(DblToChar(Occ)))
      CALL MondoLog(DEBUG_MEDIUM, Prog, "Ne = "//TRIM(DblToChar(Ne)))
      CALL MondoLog(DEBUG_MEDIUM, Prog, "Entropy = "//TRIM(DblToChar(Entropy))//" hartree")
    ENDDO
    DO K=1,NBasF
      Ek=EigenV%D(K)
      Fk = 1.D0/(EXP(Beta*(Ek-Mu))+1.D0)
      DO J=1,NBasF
        CJK=F%D(J,K)*Fk
        DO I=1,NBasF
          P%D(I,J)=P%D(I,J)+F%D(I,K)*CJK
        ENDDO
      ENDDO
    ENDDO
  CASE('NoSmearing')
    !
    SELECT CASE(NSMat)
    CASE(1)
      ! We have one density matrix to build.
      CALL DGEMM_NT(NBasF,Nel/2,NBasF,0D0,F%D(1,1),F%D(1,1),P%D(1,1))
      !CALL DGEMM('N','T',NBasF,NBasF,Nel/2,1D0,F%D(1,1), &
      !           NBasF,F%D(1,1),NBasF,0D0,P%D(1,1),NBasF)
    CASE(2)
      ! We have two density matrices to build.
      CALL DGEMM_NT(NBasF,NAlph,NBasF,0D0,F%D(1,      1),F%D(1,      1),P%D(1,      1))
      CALL DGEMM_NT(NBasF,NBeta,NBasF,0D0,F%D(1,NBasF+1),F%D(1,NBasF+1),P%D(1,NBasF+1))
      !CALL DGEMM('N','T',NBasF,NBasF ,NAlph,1D0,F%D(1,      1), &
      !           NBasF,F%D(1,      1),NBasF,0D0,P%D(1,      1),NBasF)
      !CALL DGEMM('N','T',NBasF,NBasF ,NBeta,1D0,F%D(1,NBasF+1), &
      !           NBasF,F%D(1,NBasF+1),NBasF,0D0,P%D(1,NBasF+1),NBasF)
    CASE(4)
      ! We have one density matrix to build.
      CALL DGEMM_NT(2*NBasF,Nel,2*NBasF,0D0,F%D(1,1),F%D(1,1),P%D(1,1))
      !CALL DGEMM('N','T',2*NBasF,2*NBasF,Nel,1D0,F%D(1,1), &
      !           2*NBasF,F%D(1,1),2*NBasF,0D0,P%D(1,1),2*NBasF)
    CASE DEFAULT
      CALL Halt('RHeqs: NSMat doesn''t have an expected value!')
    END SELECT
    !
  CASE DEFAULT
    CALL Halt('RHEqs: Doesn''t regonize this Smearing <'//TRIM(Smearing)//'>')
  END SELECT
  !
  CALL Delete(EigenV)

  ! Put entropy and HOMO-LUMO gap to hdf.
#if defined(PARALLEL_CLONES)
  IF(MRank(MPI_COMM_WORLD) == ROOT) THEN
    CALL Put(HOMO-LUMO, "HomoLumoGap")
    CALL Put(Entropy, "Entropy")

    oldClone = MyClone
    DO rank = 1, MSize(MPI_COMM_WORLD)-1
      CALL Recv(MyClone, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(HOMOLUMO, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)
      CALL Recv(Entropy, rank, PUT_TAG, comm_O = MPI_COMM_WORLD)

      ! Put to correct HDFGroup.
      CALL CloseHDFGroup(H5GroupID)
      H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
      HDF_CurrentID = H5GroupID
      CALL Put(HOMOLUMO, "HomoLumoGap")
      CALL Put(Entropy, "Entropy")
    ENDDO
    MyClone = oldClone

    ! Reopen old HDFGroup.
    CALL CloseHDFGroup(H5GroupID)
    H5GroupID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(MyClone)))
    HDF_CurrentID = H5GroupID
  ELSE
    CALL Send(MyClone, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(HOMO-LUMO, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
    CALL Send(Entropy, ROOT, PUT_TAG, comm_O = MPI_COMM_WORLD)
  ENDIF
#else
  CALL Put(HOMO-LUMO, "HomoLumoGap")
  CALL Put(Entropy, "Entropy")
#endif

  CALL SetEq(sX,P,nsmat_o=nsmat)          !  sX=P
  CALL New(sP,nsmat_o=nsmat)
  CALL Filter(sP,sX)        !  sP=Filter[sX]
  CALL Put(sP,TrixFile('OrthoD',Args,1))
  CALL PChkSum(sP,'OrthoP['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(sP,'OrthoP['//TRIM(NxtCycl)//']')

  CALL Plot(sP,'OrthoP['//TRIM(NxtCycl)//']')

  ! Transform to non-orthogonal rep
  XFile=TrixFile('X',Args)
  INQUIRE(FILE=XFile,EXIST=Present)

  IF(Present)THEN
    CALL Get(sX,XFile)                  ! X =S^{-1/2}
    CALL Multiply(sX,sP,sTmp1)          ! T1=S^{-1/2}.P_Orthog
    CALL Multiply(sTmp1,sX,sP)          ! T1=S^{-1/2}.P_Orthog.S^{-1/2}
  ELSE
    CALL Get(sX,TrixFile('Z',Args))     ! X=Z=L^{-1}
    CALL Multiply(sX,sP,sTmp1)          ! T1=Z.P_Orthog
    CALL Get(sX,TrixFile('ZT',Args))    ! X =Z^t=L^{-t}
    CALL Multiply(sTmp1,sX,sP)          ! F=Z.P_Orthog.Z^t
  ENDIF
  CALL Filter(sTmp1,sP)                 ! T1 =P_AO=Filter[Z.P_Orthog.Z^t]

  ! Archive the AO-DM
  CALL Get(DensityArchive,'ArchiveDensity')
  !IF(DensityArchive) THEN
  !  CALL Put(sTmp1,'CurrentDM',CheckPoint_O=.TRUE.)
  !ENDIF

  CALL Put(sTmp1,TrixFile('D',Args,1))
  CALL PChkSum(sTmp1,'P['//TRIM(NxtCycl)//']',Prog)
  CALL PPrint(sTmp1,'P['//TRIM(NxtCycl)//']')
  CALL Plot(sTmp1,'P['//TRIM(NxtCycl)//']')

  CALL Delete(sX)
  CALL Delete(sP)
  CALL Delete(sTmp1)
  CALL ShutDown(Prog)
END PROGRAM  RHEqs
