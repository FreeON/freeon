!    ROUTINES FOR THE IMPLIMENTATION OF PERIODIC BOUNDARY CONDITIONS
!    Author: C.J. Tymczak
MODULE CellSets
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE PrettyPrint
  USE Thresholding
  IMPLICIT NONE
#ifdef PERIODIC
  CONTAINS
!--------------------------------------------------------------------------
! Delete the CellSet
!--------------------------------------------------------------------------
  SUBROUTINE New_CellSet(CS,NCELL,Dimen_O)
    TYPE(CellSet)                    :: CS   
    INTEGER                          :: NCELL,Dimen
    INTEGER,OPTIONAL                 :: Dimen_O
!
    IF(Present(Dimen_O)) THEN
       Dimen = Dimen_O
    ELSE
       Dimen = 3
    ENDIF
    CS%NCells = NCELL
    IF(AllocQ(CS%Alloc)) THEN
       CALL Delete(CS%CellCarts)
       CALL New(CS%CellCarts,(/Dimen,CS%NCells/))    
    ELSE
       CALL New(CS%CellCarts,(/Dimen,CS%NCells/))    
    ENDIF
    CS%Alloc=ALLOCATED_TRUE
!
  END SUBROUTINE New_CellSet
!--------------------------------------------------------------------------
! Delete the CellSet
!--------------------------------------------------------------------------
  SUBROUTINE Delete_CellSet(CS)
    TYPE(CellSet)                  :: CS
!
    IF(AllocQ(CS%Alloc)) THEN
       CS%Alloc  = ALLOCATED_FALSE
       CS%NCells = 0
       CALL Delete(CS%CellCarts)
    ENDIF
!
  END SUBROUTINE Delete_CellSet
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
!
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
!
    IF(PrintFlags%Key==DEBUG_MEDIUM) THEN
       WRITE(OutU,30)
       WRITE(OutU,10) Name
       WRITE(OutU,11) CS%NCells,RMax
!       WRITE(OutU,12) Thresholds%Dist,SQRT(Thresholds%Pair)
       WRITE(OutU,30)
    ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
       WRITE(OutU,30)
       WRITE(OutU,10) Name
       WRITE(OutU,11) CS%NCells,RMax
!       WRITE(OutU,12) Thresholds%Dist,SQRT(Thresholds%Pair)
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
! Test if Cell(N) has Coordinates (x,y,z) 
!--------------------------------------------------------------------------
  FUNCTION CellN_XYZ(CS,NC,X,Y,Z)
    TYPE(CellSet)      :: CS
    INTEGER            :: NC
    REAL(DOUBLE)       :: X,Y,Z
    LOGICAL            :: CellN_XYZ
!
    IF(ABS(CS%CellCarts%D(1,NC)-X) < 1.0D-12 .AND. &
       ABS(CS%CellCarts%D(2,NC)-Y) < 1.0D-12 .AND. &
       ABS(CS%CellCarts%D(3,NC)-Z) < 1.0D-12) THEN
       CellN_XYZ = .TRUE.
    ELSE
       CellN_XYZ = .FALSE.
    ENDIF
!
  END FUNCTION CellN_XYZ
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
       IF(CellN_XYZ(CS,NC,X,Y,Z)) THEN
          InCell_CellSet = .TRUE.
          RETURN
       ENDIF
    ENDDO
!
  END FUNCTION InCell_CellSet
!--------------------------------------------------------------------------
! Set up the CellSet from {-N,N} minus the Cells in {-M,M}
!--------------------------------------------------------------------------
  SUBROUTINE New_CellSet_Cube(CS,MAT,N,M_O)
    TYPE(CellSet)                        :: CS
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    INTEGER,DIMENSION(3)                 :: N,M
    INTEGER,OPTIONAL,DIMENSION(3)        :: M_O
!
    INTEGER                              :: I,J,K,NCELL
    REAL(DOUBLE)                         :: X,Y,Z
!
    IF(PRESENT(M_O)) THEN 
       M = M_O
    ELSE
       M = 0
    ENDIF
!
    IF(N(1) .LT. 0) CALL Halt('N(1) is Miss Dimensioning in New_CellSet_Box')
    IF(N(2) .LT. 0) CALL Halt('N(2) is Miss Dimensioning in New_CellSet_Box')
    IF(N(3) .LT. 0) CALL Halt('N(3) is Miss Dimensioning in New_CellSet_Box')
!
    IF(M(1) .LT. 0) CALL Halt('M(1) is Miss Dimensioning in New_CellSet_Box')
    IF(M(2) .LT. 0) CALL Halt('M(2) is Miss Dimensioning in New_CellSet_Box')
    IF(M(3) .LT. 0) CALL Halt('M(3) is Miss Dimensioning in New_CellSet_Box')
!
    NCELL = 0
    DO I = -N(1),N(1)
       DO J = -N(2),N(2)
          DO K = -N(3),N(3)
             IF(IJKTest(I,J,K,M(1),M(2),M(3))) THEN
                NCELL = NCELL+1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
!
    CALL New_CellSet(CS,NCELL)
!
    NCELL = 0
    DO I = -N(1),N(1)
       DO J = -N(2),N(2)
          DO K = -N(3),N(3)
             IF(IJKTest(I,J,K,M(1),M(2),M(3))) THEN
                NCELL = NCELL+1
                X  = I*MAT(1,1)+J*MAT(1,2)+K*MAT(1,3)
                Y  = I*MAT(2,1)+J*MAT(2,2)+K*MAT(2,3)
                Z  = I*MAT(3,1)+J*MAT(3,2)+K*MAT(3,3)
                CS%CellCarts%D(1,NCELL)= X
                CS%CellCarts%D(2,NCELL)= Y
                CS%CellCarts%D(3,NCELL)= Z
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CS%Alloc=ALLOCATED_TRUE
!
  END SUBROUTINE New_CellSet_Cube
!--------------------------------------------------------------------------
! Set up the set of cells out to some radius R
!--------------------------------------------------------------------------
  SUBROUTINE New_CellSet_Sphere(CS,AW,MAT,Radius,Rmin_O)
    TYPE(CellSet)                        :: CS
    LOGICAL,DIMENSION(3)                 :: AW
    REAL(DOUBLE),DIMENSION(3,3)          :: MAT
    REAL(DOUBLE)                         :: Radius,Radius_min
    REAL(DOUBLE),OPTIONAL                :: Rmin_O
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
    IF(AW(1)) IXM = 2*(1+INT(Radius/SQRT(MAT(1,1)**2+MAT(2,1)**2+MAT(3,1)**2)))
    IF(AW(2)) IYM = 2*(1+INT(Radius/SQRT(MAT(1,2)**2+MAT(2,2)**2+MAT(3,2)**2)))
    IF(AW(3)) IZM = 2*(1+INT(Radius/SQRT(MAT(1,3)**2+MAT(2,3)**2+MAT(3,3)**2)))
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
    CALL New_CellSet(CS,NCELL)
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
    CS%Alloc=ALLOCATED_TRUE
!
  END SUBROUTINE New_CellSet_Sphere
!--------------------------------------------------------------------------
! Set up the set of cells out to some radius R
!--------------------------------------------------------------------------
  FUNCTION IJKTest(I,J,K,MI,MJ,MK)
    INTEGER               :: I,J,K,MI,MJ,MK
    LOGICAL               :: IJKTest
!
    IJKTest = .FALSE.
    IF(ABS(I) .GE. MI) THEN 
       IJKTest = .TRUE.
       RETURN
    ENDIF
    IF(ABS(J) .GE. MJ) THEN
       IJKTest = .TRUE.
       RETURN
    ENDIF
    IF(ABS(K) .GE. MK) THEN
       IJKTest = .TRUE.
       RETURN
    ENDIF
!   
  END FUNCTION IJKTest
!
#endif
!
END MODULE CellSets
