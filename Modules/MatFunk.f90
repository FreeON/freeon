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
MODULE MatFunk
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE PrettyPrint
  !
  INTEGER,SAVE   :: BIGBLOK,BLKLWORK
  TYPE(DBL_VECT) :: BLKWORK,BLKVALS
  TYPE(DBL_RNK2) :: BLKVECT,BLKMAT1,BLKMAT2
CONTAINS
  !
  SUBROUTINE SetDSYEVWork(BIGBLOK_O)
    INTEGER,OPTIONAL,INTENT(IN) :: BIGBLOK_O
    IF(PRESENT(BIGBLOK_O))THEN
      BIGBLOK=BIGBLOK_O
    ELSE
      BIGBLOK=MaxBlkSize
    ENDIF
    BLKLWORK=MAX(1,3*BIGBLOK)
    CALL New(BLKWORK,BLKLWORK)
    CALL New(BLKVALS,BIGBLOK)
    CALL New(BLKVECT,(/BIGBLOK,BIGBLOK/))
    CALL New(BLKMAT1,(/BIGBLOK,BIGBLOK/))
    CALL New(BLKMAT2,(/BIGBLOK,BIGBLOK/))
  END SUBROUTINE SetDSYEVWork
  !
  SUBROUTINE UnSetDSYEVWork()
    CALL Delete(BLKWORK)
    CALL Delete(BLKVALS)
    CALL Delete(BLKVECT)
    CALL Delete(BLKMAT1)
    CALL Delete(BLKMAT2)
  END SUBROUTINE UnSetDSYEVWork
  !
  SUBROUTINE FunkOnSqMat(N,Funk,A,FOnA,PosDefMat_O,PrintValues_O,CoNo_O,EigenThresh_O,PrintCond_O,Prog_O,FileName_O,Unit_O)
    INTEGER,                    INTENT(IN)    :: N
    REAL(DOUBLE),DIMENSION(N,N),INTENT(IN)    :: A
    REAL(DOUBLE),DIMENSION(N,N),INTENT(OUT)   :: FOnA
    LOGICAL, OPTIONAL                         :: PosDefMat_O,PrintCond_O,PrintValues_O
    REAL(DOUBLE),OPTIONAL                     :: EigenThresh_O
    CHARACTER(Len=*),OPTIONAL                 :: Prog_O
    CHARACTER(LEN=*), OPTIONAL,INTENT(IN)     :: FileName_O
    INTEGER,          OPTIONAL,INTENT(IN)     :: Unit_O
    LOGICAL                                   :: PosDefMat,PrintCond,PrintValues
    INTEGER                                   :: I,INFO,PU
    REAL(DOUBLE),OPTIONAL,INTENT(OUT)         :: CoNo_O
    REAL(DOUBLE)                              :: EigenThreshold
    REAL(DOUBLE)                              :: CondA,Mn,Mx
    REAL(DOUBLE), EXTERNAL                    :: Funk
    CHARACTER(LEN=DEFAULT_CHR_LEN)            :: String
    !----------------------------------------------------------------------------
    !
    IF(PRESENT(EigenThresh_O))THEN
      EigenThreshold=EigenThresh_O
    ELSE
      EigenThreshold=1.D-10
    ENDIF
    !
    IF(PRESENT(PrintCond_O))THEN
      PrintCond=PrintCond_O
    ELSE
      PrintCond=.FALSE.
    ENDIF
    !
    IF(PRESENT(PrintValues_O))THEN
      PrintValues=PrintValues_O
    ELSE
      PrintValues=.FALSE.
    ENDIF
    !
    IF(PRESENT(PosDefMat_O))THEN
      PosDefMat=PosDefMat_O
    ELSE
      PosDefMat=.TRUE.
    ENDIF
    !
    BLKVECT%D(1:N,1:N)=A(1:N,1:N)
    !     CALL PPrint(A,'A',Unit_O=6)

    CALL DSYEV('V','U',N,BLKVECT%D,BIGBLOK,BLKVALS%D,BLKWORK%D,BLKLWORK,INFO)
    IF(INFO/=SUCCEED) THEN
      CALL Halt('DSYEV hosed in FunkOnSqMat. INFO='//TRIM(IntToChar(INFO)))
    ENDIF

    IF(PrintValues)THEN
      CALL PPrint(BLKVALS,'EigenValues',N,FileName_O=FileName_O,Unit_O=Unit_O)
    ENDIF

    Mx=0.0D0
    Mn=1.D30
    DO I=1,N
      Mx=MAX(Mx,ABS(BLKVALS%D(I)))
      Mn=MIN(Mn,ABS(BLKVALS%D(I)))
    ENDDO
    CondA=Mx/Mn

    IF(PRESENT(CoNo_O))CoNo_O=CondA

    IF(PrintCond)THEN
      String = "Cond# = "//TRIM(DblToShrtChar(CondA))//', MIN(E) = '//TRIM(DblToShrtChar(Mn))
      IF(PRESENT(Prog_O))THEN
        CALL MondoLog(DEBUG_NONE, ProcessName(Prog_O), TRIM(String))
      ELSE
        CALL MondoLog(DEBUG_NONE, "", TRIM(String))
      ENDIF
    ENDIF
    !
    !     Apply the function to eigenvalues of the matrix, projecting out "zeros"
    !
    BLKMAT1%D=Zero
    IF(PosDefMat) THEN
      DO I=1,N
        IF(BLKVALS%D(I)<EigenThreshold)THEN
          BLKMAT1%D(I,I)=Zero
        ELSEIF(BLKVALS%D(I)<-1D-14)THEN
          CALL Halt('In MatFunk, apparently a non-positive definate matrix: Eval = ' &
               //TRIM(DblToChar(BLKVALS%D(I)))//', Perhaps increase the accuracy level?')
        ELSE
          BLKMAT1%D(I,I)=Funk(BLKVALS%D(I))
        ENDIF
      ENDDO
    ELSE
      DO I=1,N
        BLKMAT1%D(I,I)=Funk(BLKVALS%D(I))
      ENDDO
    ENDIF
    !
    !     Reconstruct the functionalized matrix
    !
    CALL DGEMM_NNc(N,N,N,One,Zero,BLKVECT%D(1:N,1:N),BLKMAT1%D(1:N,1:N),BLKMAT2%D(1:N,1:N))
    CALL DGEMM_NTc(N,N,N,One,Zero,BLKMAT2%D(1:N,1:N),BLKVECT%D(1:N,1:N),     FOnA(1:N,1:N))

  END SUBROUTINE FunkOnSqMat
  !
  !
  FUNCTION Inverse(A) RESULT(InvA)
    REAL(DOUBLE), INTENT(IN)  :: A
    REAL(DOUBLE) :: InvA
    InvA=One/A
  END FUNCTION Inverse
  !
  FUNCTION SqRoot(A) RESULT(SqRtA)
    REAL(DOUBLE), INTENT(IN)  :: A
    REAL(DOUBLE) :: SqRtA
    SqRtA=SQRT(A)
  END FUNCTION SqRoot

  FUNCTION InvSqRt(A) RESULT(InvSqRtA)
    REAL(DOUBLE), INTENT(IN)  :: A
    REAL(DOUBLE) :: InvSqRtA
    InvSqRtA=One/SQRT(A)
  END FUNCTION InvSqRt

  FUNCTION AbsInv(A) RESULT(AbsInvA)
    REAL(DOUBLE), INTENT(IN)  :: A
    REAL(DOUBLE) :: AbsInvA
    AbsInvA=One/ABS(A)
  END FUNCTION AbsInv
  !
END MODULE MatFunk
