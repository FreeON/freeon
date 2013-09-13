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
PROGRAM DMP_TS4 ! Density matrix purification, TS4 variation
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
! Trace purserving TS4
!-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                    :: F,P,POld,Tmp1,Tmp2,Tmp3
  TYPE(BCSR)                     :: F_BCSR
#else
  TYPE(BCSR)                     :: F,P,POld,Tmp1,Tmp2,Tmp3
#endif
  REAL(DOUBLE)                   :: Ne
  INTEGER                        :: I,MM
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='TS4'
!-------------------------------------------------------------------------
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  ! Suss for matrix threshold overide
  CALL SussTrix('NTFourTrix',Prog)
  ! Get the Fock matrix
  CALL New(F)
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))
  ENDIF
  ! Allocate some more matrices
  CALL New(P)
  CALL New(Pold)
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(Tmp3)
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
  ! Do TS4 iterations
  DO I=1,100
     CALL TS4(P,Tmp1,Tmp2,Tmp3,Ne,MM)
     IF(CnvrgChck(Prog,I,Ne,MM,F,P,POld,Tmp1,Tmp2))EXIT
  ENDDO
  ! Delete some obsolete matrices
  CALL Delete(F)
  ! Orthogonal put and xform to AO rep and put
  CALL PutXForm(Prog,Args,P,POld,Tmp1)
  ! Tidy up
  CALL Delete(P)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  CALL Delete(Tmp3)
  CALL ShutDown(Prog)
END PROGRAM DMP_TS4
