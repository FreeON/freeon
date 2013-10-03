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
!-------------------------------------------------------------------------------
!                         May 13th 2002
! Anders M. N. Niklasson: "Expansion Algorithm for the Density Matrix".
! Constructs the density matrix from the Hamiltonian in terms of a
! trace correcting purification expansion with 2nd order purifications.
!-------------------------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM DMP_SP2 ! Density matrix purification, SP2 variation
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
  USE MondoLogger

  IMPLICIT NONE

#ifdef PARALLEL
  TYPE(DBCSR) :: F,P,POld,Tmp1,Tmp2
  TYPE(BCSR)  :: F_BCSR
#else
  TYPE(BCSR)  :: F,FT,P,PT,POld,Tmp1,Tmp2
#endif

  !-------------------------------------------------------------------------------
  ! Trace Setting SP2
  !-------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Ne,Lambda, idempotency_error,Occ0,Occ1,Occ2,Occ3
  INTEGER                        :: I, J, MM, Imin
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='SP2'
  TYPE(DBL_RNK2)                 :: PDense
  !-------------------------------------------------------------------------------
#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif

  ! Check for matrix threshold overide
  CALL SussTrix('TCTwoTrix',Prog)

  CALL New(F)
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
    CALL Get(F,FFile)
  ELSE
    CALL Get(F,TrixFile('OrthoF',Args,0))
  ENDIF

  ! Initialize
  CALL New(P)
  CALL New(POld)
  CALL New(Tmp1)
  CALL New(Tmp2)

  MM=0
  Ne=Half*DBLE(NEl)

  ! Guess P from F
#ifdef PARALLEL
  CALL SetEq(F_BCSR,F)
  CALL FockGuess(F_BCSR,P,Ne,1)
  CALL Delete(F_BCSR)
#else
  CALL FockGuess(F,P,Ne,1)
#endif

  IF(PrintFlags%MM == DEBUG_PRT_SP2) THEN
    CALL SetEq(PDense, F)
    CALL MondoLog(DEBUG_NONE, Prog, "Fockian, " &
      //TRIM(IntToChar(size(PDense%D, 1)))//"x" &
      //TRIM(IntToChar(size(PDense%D, 2)))//" matrix, " &
      //TRIM(IntToChar(size(PDense%D, 1)*size(PDense%D, 2)))//" elements")
    DO I = 1, SIZE(PDense%D, 1)
      DO J = 1, SIZE(PDense%D, 2)
        CALL MondoLog(DEBUG_NONE, Prog, TRIM(IntToChar(I)) &
          //" "//TRIM(IntToChar(J)) &
          //" "//TRIM(DblToChar(PDense%D(I, J))))
      ENDDO
    ENDDO

    CALL SetEq(PDense, P)
    CALL MondoLog(DEBUG_NONE, Prog, "Fock guess, " &
      //TRIM(IntToChar(size(PDense%D, 1)))//"x" &
      //TRIM(IntToChar(size(PDense%D, 2)))//" matrix, " &
      //TRIM(IntToChar(size(PDense%D, 1)*size(PDense%D, 2)))//" elements")
    DO I = 1, SIZE(PDense%D, 1)
      DO J = 1, SIZE(PDense%D, 2)
        CALL MondoLog(DEBUG_NONE, Prog, TRIM(IntToChar(I)) &
          //" "//TRIM(IntToChar(J)) &
          //" "//TRIM(DblToChar(PDense%D(I, J))))
      ENDDO
    ENDDO
  ENDIF

  ! Do SP2 iterations
  Occ0 = 0.D0
  Occ1 = 0.D0
  Occ2 = 0.D0
  Occ3 = 0.D0
  Imin = 20
  DO I=1,100
    CALL TC2(P,Tmp1,Tmp2,Half*DBLE(NEl),Occ0,I)
    IF(PrintFlags%MM == DEBUG_PRT_SP2) THEN
      CALL MondoLog(DEBUG_NONE, Prog, "iteration "//TRIM(IntToChar(I)) &
        //": Tr(P) = "//TRIM(DblToChar(Occ0)) &
        //", Ne/2 = "//TRIM(DblToChar(Half*DBLE(NEl))) &
        //", Tr(P)-Ne/2 = "//TRIM(DblToChar(Occ0-Half*DBLE(NEl))))
    ENDIF
    IF(IdmpCnvrgChck(Occ0,Occ1,Occ2,Occ3,Imin,I)) THEN
      CALL MondoLog(DEBUG_MAXIMUM, Prog, "converged in "//TRIM(IntToChar(I))//" iterations")
      CALL MondoLog(DEBUG_MAXIMUM, Prog, "Idempotency error = "//TRIM(DblToChar(ABS(Occ0-Occ1))))
      CALL MondoLog(DEBUG_MAXIMUM, Prog, "Previous idempotency error = "//TRIM(DblToChar(ABS(Occ2-Occ3))))
      EXIT
    ENDIF
    Occ3 = Occ2
    Occ2 = Occ1
    Occ1 = Occ0
  ENDDO

  IF(PrintFlags%MM == DEBUG_PRT_SP2) THEN
    CALL SetEq(PDense, P)
    CALL MondoLog(DEBUG_NONE, Prog, "Converged density, " &
      //TRIM(IntToChar(size(PDense%D, 1)))//"x" &
      //TRIM(IntToChar(size(PDense%D, 2)))//" matrix, " &
      //TRIM(IntToChar(size(PDense%D, 1)*size(PDense%D, 2)))//" elements")
    DO I = 1, SIZE(PDense%D, 1)
      DO J = 1, SIZE(PDense%D, 2)
        CALL MondoLog(DEBUG_NONE, Prog, TRIM(IntToChar(I)) &
          //" "//TRIM(IntToChar(J)) &
          //" "//TRIM(DblToChar(PDense%D(I, J))))
      ENDDO
    ENDDO
  ENDIF

  ! If we are called without the DIIS Fockian, consider a levelshift
  CALL OpenASCII(InpFile,Inp)
  IF(.NOT.PRESENT.AND.OptDblQ(Inp,'LevelShift',Lambda))THEN
    ! Get the Fock matrix back ...
    CALL Get(F,TrixFile('OrthoF',Args,0))
    ! Construct the virtual projector, Q=I-P
    CALL SetEq(POld,P)
    CALL Multiply(POld,-One)
    CALL Add(POld,One)
    ! The shifted Fockian is F[Lambda] = P.F + (1+Lambda)*Q.F
    Lambda=One+ABS(Lambda)
    CALL Multiply(P,F,Tmp1)
    CALL Multiply(POld,F,Tmp2)
    CALL Multiply(Tmp2,Lambda)
    CALL Add(Tmp1,Tmp2,POld)
    ! Put the Fock matrix back all tidy like
    CALL Put(POld,TrixFile('OrthoF',Args,0))
    Lambda=ABS(Lambda)-One
    CALL MondoLog(DEBUG_NONE, "SP2", TRIM(ProcessName(Prog)) &
      //' LevlShift = '//TRIM(DblToMedmChar(Lambda)))
  ENDIF
  CLOSE(Inp)

  ! Orthogonal put and xform to AO rep and put
  CALL PutXForm(Prog,Args,P,POld,Tmp1)

  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(POld)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)

  CALL ShutDown(Prog)

END PROGRAM DMP_SP2
