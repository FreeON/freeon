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
!    ROUTINES FOR THE IMPLIMENTATION OF PERIODIC BOUNDARY CONDITIONS
!    Author: C.J. Tymczak

MODULE CellSets
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE PrettyPrint
  USE Thresholding
  USE Order
  USE MemMan
  USE MondoLogger

  IMPLICIT NONE

CONTAINS
  !--------------------------------------------------------------------------
  ! Print the CellSet
  !--------------------------------------------------------------------------
  SUBROUTINE PPrint_CellSet(CS,Name,FileName_O,Unit_O)
    TYPE(CellSet)                    :: CS
    CHARACTER(LEN=*)                 :: Name
    CHARACTER(LEN=*),OPTIONAL        :: FileName_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: OutU
    INTEGER                          :: NC
    REAL(DOUBLE)                     :: RMax,R2

    IF(.NOT. AllocQ(CS%Alloc)) THEN
      CALL Halt(' Cells are  not allocated in PPrint_CellSet')
    ENDIF
    IF(PRESENT(Unit_O)) THEN
      OutU=Unit_O
    ELSE
      OutU=Out
    ENDIF
    IF(PRESENT(FileName_O) .AND. OutU /= 6) THEN
      CALL OpenASCII(FileName_O,OutU)
    ELSEIF(OutU /= 6) THEN
      CALL OpenASCII(OutFile,OutU)
    ENDIF
    !
    RMax = Zero
    DO NC = 1,CS%NCells
      R2 = CS%CellCarts%D(1,NC)**2+CS%CellCarts%D(2,NC)**2+CS%CellCarts%D(3,NC)**2
      IF(R2 > RMax) RMax = R2
    ENDDO
    RMax = SQRT(RMax)

    IF(PrintFlags%Key==DEBUG_MEDIUM) THEN
      WRITE(OutU,30)
      WRITE(OutU,10) Name
      WRITE(OutU,11) CS%NCells,RMax
      WRITE(OutU,30)
    ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
      WRITE(OutU,30)
      WRITE(OutU,10) Name
      WRITE(OutU,11) CS%NCells,RMax
      IF(CS%NCells /= 0) THEN
        WRITE(OutU,20)
        WRITE(OutU,31)
        DO NC = 1,CS%NCells
          WRITE(OutU,21) NC,CS%CellCarts%D(1,NC),CS%CellCarts%D(2,NC),CS%CellCarts%D(3,NC)
        ENDDO
      ENDIF
      WRITE(OutU,30)
    ENDIF
    RETURN
10  FORMAT(1x,A)
11  FORMAT(1x,'Number of Cells    = ',I8,2x,'   Maxium Range   = ',D15.8)
12  FORMAT(1x,'Distance Threshold = ',D10.4,'   Pair Threshold = ',D15.8)
20  FORMAT(1x,'<Coordinates of the Cells>')
21  FORMAT(2x,'Cell #',I3,2x,'Q = (',D15.8,', ',D15.8,', ',D15.8,' )')
30  FORMAT(1x,'=========================================================================')
31  FORMAT(1x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
  END SUBROUTINE PPrint_CellSet
  !--------------------------------------------------------------------------
  ! Test if a Cell in CS has Coordinates (x,y,z)
  !--------------------------------------------------------------------------
  FUNCTION InCell_CellSet(CS,X,Y,Z)
    TYPE(CellSet)      :: CS
    INTEGER            :: NC
    REAL(DOUBLE)       :: X,Y,Z
    LOGICAL            :: InCell_CellSet
    !
    InCell_CellSet = .FALSE.
    DO NC = 1,CS%NCells
      IF(ABS(CS%CellCarts%D(1,NC)-X) < 1.0D-10 .AND. &
           ABS(CS%CellCarts%D(2,NC)-Y) < 1.0D-10 .AND. &
           ABS(CS%CellCarts%D(3,NC)-Z) < 1.0D-10) THEN
        InCell_CellSet = .TRUE.
        RETURN
      ENDIF
    ENDDO
    !
  END FUNCTION InCell_CellSet
  !--------------------------------------------------------------------------
  ! Set up the CellSet from {-N,N} minus the Cells in {-M,M}
  !--------------------------------------------------------------------------
  SUBROUTINE New_CellSet_Cube(CS,AW,MAT,N,MaxCell_O)
    TYPE(CellSet)                        :: CS
    INTEGER,DIMENSION(3)                 :: AW
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    INTEGER,DIMENSION(3)                 :: N
    INTEGER,OPTIONAL                     :: MaxCell_O
    !
    INTEGER                              :: I,J,K,NCELL
    REAL(DOUBLE)                         :: X,Y,Z
    !
    IF(N(1) .LT. 0) CALL Halt('N(1) is Miss Dimensioning in New_CellSet_Box')
    IF(N(2) .LT. 0) CALL Halt('N(2) is Miss Dimensioning in New_CellSet_Box')
    IF(N(3) .LT. 0) CALL Halt('N(3) is Miss Dimensioning in New_CellSet_Box')
    !
    IF(AW(1)==0) N(1)=0
    IF(AW(2)==0) N(2)=0
    IF(AW(3)==0) N(3)=0
    !
    NCELL = 0
    DO I = -N(1),N(1)
      DO J = -N(2),N(2)
        DO K = -N(3),N(3)
          NCELL = NCELL+1
        ENDDO
      ENDDO
    ENDDO
    !
    IF(PRESENT(MaxCell_O)) THEN
      IF(NCELL > MaxCell_O) THEN
        CALL MondoLogPlain('NCELL   = '//TRIM(IntToChar(NCELL)))
        CALL MondoLogPlain('MaxCell = '//TRIM(IntToChar(MaxCell_O)))
        CALL Halt('NCELL is Greater then MaxCell')
      ELSE
        CALL New_CellSet(CS,MaxCell_O)
        CS%NCells=NCELL
      ENDIF
    ELSE
      CALL New_CellSet(CS,NCELL)
      CS%NCells=NCELL
    ENDIF
    !
    NCELL = 0
    DO I = -N(1),N(1)
      DO J = -N(2),N(2)
        DO K = -N(3),N(3)
          NCELL = NCELL+1
          X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
          Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
          Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
          CS%CellCarts%D(1,NCELL)= X
          CS%CellCarts%D(2,NCELL)= Y
          CS%CellCarts%D(3,NCELL)= Z
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE New_CellSet_Cube



  SUBROUTINE New_CellSet_Sphere(CS,AW,MAT,Radius,Rmin_O,MaxCell_O)
    TYPE(CellSet)                        :: CS
    INTEGER,DIMENSION(3)                 :: AW
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    REAL(DOUBLE)                         :: Radius,Radius_min
    REAL(DOUBLE),OPTIONAL                :: Rmin_O
    INTEGER,OPTIONAL                     :: MaxCell_O
    !
    INTEGER                              :: I,J,K
    INTEGER                              :: IXM,IYM,IZM,NCELL
    REAL(DOUBLE)                         :: X,Y,Z,Rad,R
    !
    IF(PRESENT(Rmin_O)) THEN
      Radius_min = Rmin_O
    ELSE
      Radius_min = Zero
    ENDIF
    !
    IXM = 0
    IYM = 0
    IZM = 0
    IF(AW(1)==1) IXM = 2*(1+INT(Radius/SQRT(MAT(1,1)**2+MAT(2,1)**2+MAT(3,1)**2)))
    IF(AW(2)==1) IYM = 2*(1+INT(Radius/SQRT(MAT(1,2)**2+MAT(2,2)**2+MAT(3,2)**2)))
    IF(AW(3)==1) IZM = 2*(1+INT(Radius/SQRT(MAT(1,3)**2+MAT(2,3)**2+MAT(3,3)**2)))
    !
    NCELL = 0
    DO I=-IXM,IXM
      DO J=-IYM,IYM
        DO K=-IZM,IZM
          X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
          Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
          Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
          R = SQRT(X*X+Y*Y+Z*Z)
          IF(R .GE. Radius_min .AND. R .LT. Radius) THEN
            NCELL = NCELL+1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    !
    IF(PRESENT(MaxCell_O)) THEN
      IF(NCELL > MaxCell_O) THEN
        CALL MondoLogPlain('NCELL   = '//TRIM(IntToChar(NCELL)))
        CALL MondoLogPlain('MaxCell = '//TRIM(IntToChar(MaxCell_O)))
        CALL Halt('NCELL is Greater then MaxCell')
      ELSE
        CALL New_CellSet(CS,MaxCell_O)
        CS%NCells=NCELL
      ENDIF
    ELSE
      CALL New_CellSet(CS,NCELL)
      CS%NCells=NCELL
    ENDIF
    !
    NCELL=0
    DO I=-IXM,IXM
      DO J=-IYM,IYM
        DO K=-IZM,IZM
          X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
          Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
          Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
          R = SQRT(X*X+Y*Y+Z*Z)
          IF(R .GE. Radius_min .AND. R .LT. Radius) THEN
            NCELL = NCELL+1
            CS%CellCarts%D(1,NCELL)= X
            CS%CellCarts%D(2,NCELL)= Y
            CS%CellCarts%D(3,NCELL)= Z
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE New_CellSet_Sphere
  !--------------------------------------------------------------------------
  ! Set up the set of cells out to some radius MinExp and Threshold
  !--------------------------------------------------------------------------
  SUBROUTINE New_CellSet_Sphere2(CS,AW,MAT,Option,CellSetsThreshold,MinExpt,NElect,Rmin_O,MaxCell_O)
    TYPE(CellSet)                        :: CS
    INTEGER,DIMENSION(3)                 :: AW
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    REAL(DOUBLE),DIMENSION(3)            :: UpperUC,LowerUC
    REAL(DOUBLE)                         :: CellSetsThreshold,MinExpt,MinSepSq,Sep2,T
    REAL(DOUBLE),OPTIONAL                :: Rmin_O
    REAL(DOUBLE),DIMENSION(3,8)          :: UCCorners
    REAL(DOUBLE),DIMENSION(3)            :: XYZVec
    INTEGER,OPTIONAL                     :: MaxCell_O
    CHARACTER(LEN=*)                     :: Option
    !
    INTEGER                              :: I,J,K,L,M,Corners,NElect, HardP,HardO
    INTEGER                              :: IXM,IYM,IZM,NCELL
    REAL(DOUBLE)                         :: Rad,R,TTemp
    LOGICAL                              :: InCell
    !
    !  Make the Corners
    !

    Corners=0
    DO I=0,1*AW(1)
      DO J=0,1*AW(2)
        DO K=0,1*AW(3)
          Corners=Corners+1
          UCCorners(:,Corners)=(/I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3),  &
               I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3),  &
               I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)/)
        ENDDO
      ENDDO
    ENDDO
    !
    IXM = 0
    IYM = 0
    IZM = 0
    IF(AW(1)==1) IXM=20
    IF(AW(2)==1) IYM=20
    IF(AW(3)==1) IZM=20
    !
    NCELL = 0
    DO I=-IXM,IXM
      DO J=-IYM,IYM
        DO K=-IZM,IZM
          XYZVec(1)=I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
          XYZVec(2)=I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
          XYZVec(3)=I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
          MinSepSq = MiniumSeperation(XYZVec,Corners,UCCorners)+1.0D-10
          IF(Option=='Penetration')THEN
            T=MinExpt*MinSepSq
            InCell=DBLE(NElect)*EXP(-T)/MinSepSq>CellSetsThreshold
          ELSEIF(Option=='Overlap')THEN
            T=Half*MinExpt*MinSepSq
            InCell=EXP(-T)>CellSetsThreshold
          ELSE
            CALL Halt(' Unknown option to New_CellSet_Sphere2 ')
          ENDIF
          IF(InCell)THEN
            NCELL = NCELL+1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    IF(PRESENT(MaxCell_O)) THEN
      IF(NCELL > MaxCell_O) THEN
        CALL MondoLogPlain('NCELL   = '//TRIM(IntToChar(NCELL)))
        CALL MondoLogPlain('MaxCell = '//TRIM(IntToChar(MaxCell_O)))
        CALL Halt('NCELL is Greater then MaxCell')
      ELSE
        CALL New_CellSet(CS,MaxCell_O)
        CS%NCells=NCELL
      ENDIF
    ELSE
      CALL New_CellSet(CS,NCELL)
      CS%NCells=NCELL
    ENDIF
    !
    NCELL = 0
    DO I=-IXM,IXM
      DO J=-IYM,IYM
        DO K=-IZM,IZM
          XYZVec(1)=I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
          XYZVec(2)=I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
          XYZVec(3)=I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
          MinSepSq = MiniumSeperation(XYZVec,Corners,UCCorners)+1.0D-10
          IF(Option=='Penetration')THEN
            T=MinExpt*MinSepSq
            InCell=DBLE(NElect)*EXP(-T)/MinSepSq>CellSetsThreshold
          ELSEIF(Option=='Overlap')THEN
            T=Half*MinExpt*MinSepSq
            InCell=EXP(-T)>CellSetsThreshold
          ELSE
            CALL Halt(' Unknown option to New_CellSet_Sphere2 ')
          ENDIF
          IF(InCell)THEN
            NCELL = NCELL+1
            CS%CellCarts%D(:,NCELL)=XYZVec(:)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CALL MondoLog(DEBUG_NONE, "New_CellSet_Sphere2", 'NEW '//Option//' CELLSET SPHERE, NCELL = '//TRIM(IntToChar(NCELL)))

  END SUBROUTINE New_CellSet_Sphere2
  !--------------------------------------------------------------------------
  !  Minimum Seperation Between Celll Corners
  !--------------------------------------------------------------------------
  FUNCTION MiniumSeperation(XYZVec,Corners,UCCorners) RESULT(MinSepSq)
    REAL(DOUBLE)                         :: MinSepSq
    REAL(DOUBLE),DIMENSION(3,8)          :: UCCorners
    REAL(DOUBLE),DIMENSION(3)            :: XYZVec,Dist
    INTEGER                              :: I,J,Corners
    !
    MinSepSq=1D20
    DO I=1,Corners
      DO J=1,Corners
        Dist(:) = (XYZVec+UCCorners(:,I))-UCCorners(:,J)
        MinSepSq=MIN(MinSepSq,Dist(1)**2+Dist(2)**2+Dist(3)**2)
      ENDDO
    ENDDO
    !
  END FUNCTION MiniumSeperation
  !--------------------------------------------------------------------------
  ! Sort the Cells From Large R to Small R, or small R to Large R
  !--------------------------------------------------------------------------
  SUBROUTINE Sort_CellSet(CS,Order_O)
    TYPE(CellSet)                        :: CS
    TYPE(INT_VECT)                       :: IPnt
    TYPE(DBL_VECT)                       :: Vec
    TYPE(DBL_RNK2)                       :: CCarts
    INTEGER                              :: NC
    INTEGER, OPTIONAL                    :: Order_O
    INTEGER                              :: Order

    IF(PRESENT(Order_O))THEN
      Order=Order_O
    ELSE
      Order=2
    ENDIF

    CALL New(IPnt  ,CS%NCells)
    CALL New(Vec   ,CS%NCells)
    CALL New(CCarts,(/3,CS%NCells/))
    DO NC=1,CS%NCells
      IPnt%I(NC)= NC
      Vec%D(NC) = SQRT(CS%CellCarts%D(1,NC)**2+CS%CellCarts%D(2,NC)**2+CS%CellCarts%D(3,NC)**2)
    ENDDO
    !
    CALL Sort(Vec,IPnt,CS%NCells,Order)
    !
    DO NC=1,CS%NCells
      CCarts%D(1,NC) = CS%CellCarts%D(1,IPnt%I(NC))
      CCarts%D(2,NC) = CS%CellCarts%D(2,IPnt%I(NC))
      CCarts%D(3,NC) = CS%CellCarts%D(3,IPnt%I(NC))
    ENDDO
    !
    DO NC=1,CS%NCells
      CS%CellCarts%D(1,CS%NCells-NC+1) = CCarts%D(1,NC)
      CS%CellCarts%D(2,CS%NCells-NC+1) = CCarts%D(2,NC)
      CS%CellCarts%D(3,CS%NCells-NC+1) = CCarts%D(3,NC)
    ENDDO
    CALL Delete(IPnt)
    CALL Delete(Vec)
    CALL Delete(CCarts)
    !
  END SUBROUTINE Sort_CellSet
  !
END MODULE CellSets
