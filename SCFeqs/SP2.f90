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
  TYPE(DBCSR)                     :: F,P,Pold,Tmp1,Tmp2
  TYPE(BCSR)                      :: F_BCSR
#else
  TYPE(BCSR)                     :: F,FT,P,PT,Pold,Tmp1,Tmp2
#endif
  !-------------------------------------------------------------------------------
  ! Trace Setting SP2
  !-------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  REAL(DOUBLE)                   :: Ne,Lambda, idempotency_error
  INTEGER                        :: I, I2, MM
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='SP2'
  !-------------------------------------------------------------------------------
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  ! Suss for matrix threshold overide
  CALL SussTrix('SPTwoTrix',Prog)
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
  CALL New(Pold)
  CALL New(Tmp1)
  CALL New(Tmp2)

  ! Some symmetry checks.
  CALL New(PT)
  CALL New(FT)

  CALL Xpose(F, FT)
  FT%MTrix%D = -FT%MTrix%D
  CALL Add(F, FT, Tmp1)
  CALL MondoLog(DEBUG_NONE, Prog, "FNorm(F-FT) = "//TRIM(FltToChar(FNorm(Tmp1))))

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
  CALL SetEq(Pold,P)
  ! Do SP2 iterations
  DO I=1,150
    CALL SP2(P,Tmp1,Tmp2,Ne,MM)
    IF(CnvrgChck(Prog,I,Ne,MM,F,P,POld,Tmp1,Tmp2)) EXIT
    !idempotency_error = ABS(Trace(P) - Trace(Tmp1))
    !WRITE(*,*) "I = ", I, ", idempotency_error = ", idempotency_error
    !IF(idempotency_error < 1d-3) THEN
    !  WRITE(*,*) "I = ", I, ", idempotency_error = ", idempotency_error
    !  DO I2 = 1, 6
    !    CALL SP2(P, Tmp1, Tmp2, Ne, MM)
    !  ENDDO
    !  EXIT
    !ENDIF
  ENDDO
  CALL OpenASCII(InpFile,Inp)
  ! If we are called without the DIIS Fockian, consider a levelshift
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

  ! Symmetry check.
  CALL Xpose(P, PT)
  PT%MTrix%D = -PT%MTrix%D
  CALL Add(P, PT, Tmp2)
  CALL MondoLog(DEBUG_NONE, "SP2", "FNorm(P-PT) = "//TRIM(FltToChar(FNorm(Tmp2))))

  ! Tidy up
  CALL Delete(PT)
  CALL Delete(FT)
  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(Pold)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)

  CALL ShutDown(Prog)

END PROGRAM DMP_SP2
