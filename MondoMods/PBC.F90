!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    TYPE CellSet and Assoicated Functions
!
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
  SUBROUTINE New_CellSet(CS,NCELL)
    TYPE(CellSet)                    :: CS   
    INTEGER                          :: NCELL
!
    CS%NCells = NCELL
    IF(AllocQ(CS%Alloc)) THEN
       CALL Delete(CS%CellCarts)
       CALL New(CS%CellCarts,(/3,CS%NCells/))    
    ELSE
       CALL New(CS%CellCarts,(/3,CS%NCells/))    
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
!
!--------------------------------------------------------------------------
! Set up the CellSet from {-N,N} minus the Cells in {-M,M}
!--------------------------------------------------------------------------
  SUBROUTINE New_CellSet_Cube(CS,GM,N,M_O)
    TYPE(CellSet)                  :: CS
    TYPE(CRDS)                     :: GM
    INTEGER                        :: I,J,K,NCELL
    REAL(DOUBLE)                   :: X,Y,Z
    INTEGER,DIMENSION(3)           :: N
    INTEGER,OPTIONAL,DIMENSION(3)  :: M_O
!
    IF(.NOT. PRESENT(M_O)) THEN
       IF(.NOT. GM%AutoW(1)) N(1) = 0
       IF(.NOT. GM%AutoW(2)) N(2) = 0
       IF(.NOT. GM%AutoW(3)) N(3) = 0
       NCELL = (2*N(1)+1)*(2*N(2)+1)*(2*N(3)+1)
       CALL New_CellSet(CS,NCELL)
!
       NCELL = 0
       DO I = -N(1),N(1)
          DO J = -N(2),N(2)
             DO K = -N(3),N(3)
                NCELL = NCELL+1
                X  = I*GM%BoxShape%D(1,1)
                Y  = I*GM%BoxShape%D(1,2)+J*GM%BoxShape%D(2,2)
                Z  = I*GM%BoxShape%D(1,3)+J*GM%BoxShape%D(2,3) &
                     + K*GM%BoxShape%D(3,3)
                CS%CellCarts%D(1,NCELL)= X
                CS%CellCarts%D(2,NCELL)= Y
                CS%CellCarts%D(3,NCELL)= Z
             ENDDO
          ENDDO
       ENDDO
!
    ELSE
       IF(M_O(1) > N(1)) CALL Halt('M_O(1) is Miss Dimensioning in New_CellSet_Box')
       IF(M_O(2) > N(2)) CALL Halt('M_O(2) is Miss Dimensioning in New_CellSet_Box')
       IF(M_O(3) > N(3)) CALL Halt('M_O(3) is Miss Dimensioning in New_CellSet_Box')
       DO I = 1,3
          IF(.NOT. GM%AutoW(I)) THEN
             N(I)   = 0
             M_O(I) = 0
          ENDIF
       ENDDO
       NCELL = (2*N(1)+1)*(2*N(2)+1)*(2*N(3)+1)-(2*M_O(1)+1)*(2*M_O(2)+1)*(2*M_O(3)+1)
       CALL New_CellSet(CS,NCELL)
       IF(NCELL == 0) RETURN
!
       NCELL = 0
       DO I = -N(1),N(1)
          DO J = -N(2),N(2)
             DO K = -N(3),N(3)
                IF(I .LT. -M_O(1) .OR. I .GT. M_O(1) .OR. & 
                   J .LT. -M_O(2) .OR. J .GT. M_O(2) .OR. &
                   K .LT. -M_O(3) .OR. K .GT. M_O(3) ) THEN
                   NCELL = NCELL+1
                   X  = I*GM%BoxShape%D(1,1)
                   Y  = I*GM%BoxShape%D(1,2)+J*GM%BoxShape%D(2,2)
                   Z  = I*GM%BoxShape%D(1,3)+J*GM%BoxShape%D(2,3) &
                        + K*GM%BoxShape%D(3,3)
                   CS%CellCarts%D(1,NCELL)= X
                   CS%CellCarts%D(2,NCELL)= Y
                   CS%CellCarts%D(3,NCELL)= Z
                ENDIF
             ENDDO
          ENDDO
       ENDDO
!
    ENDIF
    CS%Alloc=ALLOCATED_TRUE
!
  END SUBROUTINE New_CellSet_Cube
!--------------------------------------------------------------------------
! Set up the set of cells out to some radius R
!--------------------------------------------------------------------------
  SUBROUTINE New_CellSet_Sphere(CS,GM,Radius)
    TYPE(CellSet)                  :: CS
    TYPE(CRDS)                     :: GM
    REAL(DOUBLE)                   :: Radius
    INTEGER                        :: I,J,K
    INTEGER                        :: IXM,IYM,IZM,NCELL
    REAL(DOUBLE)                   :: X,Y,Z,Rad,R,R2
!
    IXM = 0
    IYM = 0
    IZM = 0
    IF(GM%AutoW(1)) THEN
       IXM = 1+INT(SQRT(Radius/(GM%BoxShape%D(1,1)**2)))
    ENDIF
    IF(GM%AutoW(2)) THEN
       IYM = 1+INT(SQRT(Radius/(GM%BoxShape%D(1,2)**2 &
            +GM%BoxShape%D(2,2)**2)))
    ENDIF
    IF(GM%AutoW(3)) THEN
       IZM = 1+INT(SQRT(Radius/(GM%BoxShape%D(1,3)**2 &
            +GM%BoxShape%D(2,3)**2 &
            +GM%BoxShape%D(3,3)**2)))
    ENDIF
!
    NCELL = 0
    DO I=-IXM,IXM
       DO J=-IYM,IYM
          DO K=-IZM,IZM
             X  = I*GM%BoxShape%D(1,1)
             Y  = I*GM%BoxShape%D(1,2)+J*GM%BoxShape%D(2,2)
             Z  = I*GM%BoxShape%D(1,3)+J*GM%BoxShape%D(2,3)+K*GM%BoxShape%D(3,3)
             R2 = X*X+Y*Y+Z*Z
             IF(ABS(I) .LE. 1 .AND. ABS(J) .LE. 1 .AND. ABS(K) .LE. 1) THEN
                NCELL = NCELL+1
             ELSE
                IF(R2 < Radius) THEN
                   NCELL = NCELL+1
                ENDIF
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
             X  = I*GM%BoxShape%D(1,1)
             Y  = I*GM%BoxShape%D(1,2)+J*GM%BoxShape%D(2,2)
             Z  = I*GM%BoxShape%D(1,3)+J*GM%BoxShape%D(2,3)+K*GM%BoxShape%D(3,3)
             R2 = X*X+Y*Y+Z*Z
             IF(ABS(I) .LE. 1 .AND. ABS(J) .LE. 1 .AND. ABS(K) .LE. 1) THEN
                NCELL = NCELL+1
                CS%CellCarts%D(1,NCELL)= X
                CS%CellCarts%D(2,NCELL)= Y
                CS%CellCarts%D(3,NCELL)= Z
             ELSE
                IF(R2 < Radius) THEN
                   NCELL = NCELL+1
                   CS%CellCarts%D(1,NCELL)= X
                   CS%CellCarts%D(2,NCELL)= Y
                   CS%CellCarts%D(3,NCELL)= Z
                ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CS%Alloc=ALLOCATED_TRUE
!
  END SUBROUTINE New_CellSet_Sphere
#endif
!
END MODULE CellSets
