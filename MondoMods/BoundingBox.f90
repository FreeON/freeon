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
MODULE BoundingBox
   USE DerivedTypes
   IMPLICIT NONE
   CONTAINS
!===============================================================================
!
!===============================================================================
      SUBROUTINE CalCenterAndHalf(A)
         TYPE(BBox) :: A
         A%Center(1:3) = (A%BndBox(1:3,2)+A%BndBox(1:3,1))*Half
         A%Half(1:3) =   (A%BndBox(1:3,2)-A%BndBox(1:3,1))*Half
      END SUBROUTINE CalCenterAndHalf

!===============================================================================
!
!===============================================================================
      SUBROUTINE SetBBox(A,B)
         TYPE(BBox) :: A,B
!--------------------------------------------------------------------------------
         B%BndBox=A%BndBox
         B%Center=A%Center
         B%Half  =A%Half
      END SUBROUTINE SetBBox
!===============================================================================
!
!===============================================================================
      SUBROUTINE SplitBox(Node,Left,Right,ISplit_O,NewWall_O)
         TYPE(BBox)                         :: Node,Left,Right
         REAL(DOUBLE),DIMENSION(2),OPTIONAL :: NewWall_O
         REAL(DOUBLE)                       :: DHalf,LeftWall,RightWall
         INTEGER,                  OPTIONAL :: ISplit_O
         INTEGER                            :: ISplit
!--------------------------------------------------------------------------------
         CALL SetBBox(Node,Left)
         CALL SetBBox(Node,Right)
!        Split the bounding box
         IF(PRESENT(ISplit_O))THEN
            ISplit=ISplit_O
         ELSE
            ISplit=MOD(Node%Tier,3)+1
         ENDIF
         IF(PRESENT(NewWall_O))THEN
            LeftWall=NewWall_O(1)
            RightWall=NewWall_O(2)
         ELSE
            DHalf=Half*(Node%BndBox(ISplit,2)-Node%BndBox(ISplit,1))
            LeftWall=Left%BndBox(ISplit,2)-DHalf
            RightWall=Right%BndBox(ISplit,1)+DHalf
         ENDIF
!        New bounding boxes
         Left%BndBox(ISplit,2)=LeftWall
         Right%BndBox(ISplit,1)=RightWall
!        New sides
         Left%Half(ISplit)=Half*(Left%BndBox(ISplit,2)-Left%BndBox(ISplit,1))
         Right%Half(ISplit)=Half*(Right%BndBox(ISplit,2)-Right%BndBox(ISplit,1))
!        New centers
         Left%Center(ISplit)=Half*(Left%BndBox(ISplit,2)+Left%BndBox(ISplit,1))
         Right%Center(ISplit)=Half*(Right%BndBox(ISplit,2)+Right%BndBox(ISplit,1))
      END SUBROUTINE SplitBox

!===============================================================================
!     Determine if a point with extent is outside a BBox
!===============================================================================
      FUNCTION FurthestPointInBox(Point,Box)
        TYPE(BBox)                  :: Box
        REAL(DOUBLE)                :: PCDist2,FurthestPointInBox
        REAL(DOUBLE),DIMENSION(3)   :: Point
        REAL(DOUBLE),DIMENSION(3,8) :: Corner
        INTEGER                     :: I
        Corner(:,1)=(/Box%BndBox(1,1),Box%BndBox(2,1),Box%BndBox(3,1)/)
        Corner(:,2)=(/Box%BndBox(1,2),Box%BndBox(2,1),Box%BndBox(3,1)/)
        Corner(:,3)=(/Box%BndBox(1,1),Box%BndBox(2,2),Box%BndBox(3,1)/)
        Corner(:,4)=(/Box%BndBox(1,1),Box%BndBox(2,1),Box%BndBox(3,2)/)
        Corner(:,5)=(/Box%BndBox(1,2),Box%BndBox(2,1),Box%BndBox(3,2)/)
        Corner(:,6)=(/Box%BndBox(1,2),Box%BndBox(2,2),Box%BndBox(3,1)/)
        Corner(:,7)=(/Box%BndBox(1,1),Box%BndBox(2,2),Box%BndBox(3,2)/)
        Corner(:,8)=(/Box%BndBox(1,2),Box%BndBox(2,2),Box%BndBox(3,2)/)
        PCDist2=Zero
        DO I=1,8
           PCDist2=MAX(PCDist2,DOT_PRODUCT(Corner(:,I)-Point,Corner(:,I)-Point))
        ENDDO
        FurthestPointInBox=SQRT(PCDist2)
      END FUNCTION FurthestPointInBox


      SUBROUTINE BoxMerge(Left,Right,Box)
        TYPE(BBox)       :: Left,Right,Box
        INTEGER          :: K
        DO K=1,3
           Box%BndBox(K,1)=MIN(Left%BndBox(K,1),Right%BndBox(K,1))
           Box%BndBox(K,2)=MAX(Left%BndBox(K,2),Right%BndBox(K,2))
        ENDDO
        Box%Half(1)  = Half*(Box%BndBox(1,2)-Box%BndBox(1,1))
        Box%Half(2)  = Half*(Box%BndBox(2,2)-Box%BndBox(2,1))
        Box%Half(3)  = Half*(Box%BndBox(3,2)-Box%BndBox(3,1))
        Box%Center(1)= Half*(Box%BndBox(1,2)+Box%BndBox(1,1))
        Box%Center(2)= Half*(Box%BndBox(2,2)+Box%BndBox(2,1))
        Box%Center(3)= Half*(Box%BndBox(3,2)+Box%BndBox(3,1))
      END SUBROUTINE BoxMerge

      SUBROUTINE PointBoxMerge(Left,Point,Box)
        TYPE(BBox)                :: Left,Box
        REAL(DOUBLE),DIMENSION(3) :: Point
        INTEGER                   :: K
        DO K=1,3
           Box%BndBox(K,1)=MIN(Left%BndBox(K,1),Point(K))
           Box%BndBox(K,2)=MAX(Left%BndBox(K,2),Point(K))
        ENDDO
        Box%Half(1)  = Half*(Box%BndBox(1,2)-Box%BndBox(1,1))
        Box%Half(2)  = Half*(Box%BndBox(2,2)-Box%BndBox(2,1))
        Box%Half(3)  = Half*(Box%BndBox(3,2)-Box%BndBox(3,1))
        Box%Center(1)= Half*(Box%BndBox(1,2)+Box%BndBox(1,1))
        Box%Center(2)= Half*(Box%BndBox(2,2)+Box%BndBox(2,1))
        Box%Center(3)= Half*(Box%BndBox(3,2)+Box%BndBox(3,1))
      END SUBROUTINE PointBoxMerge
!===============================================================================
!     Determine if a point with extent is outside a BBox
!===============================================================================
      FUNCTION PointOutSideBox(R,Box)
         LOGICAL                                :: PointOutSideBox
         REAL(DOUBLE),DIMENSION(3),INTENT(IN)   :: R
         TYPE(BBox),  INTENT(IN)                :: Box
         PointOutSideBox=.TRUE.
         IF(R(1)<Box%BndBox(1,1))RETURN
         IF(R(1)>Box%BndBox(1,2))RETURN
         IF(R(2)<Box%BndBox(2,1))RETURN
         IF(R(2)>Box%BndBox(2,2))RETURN
         IF(R(3)<Box%BndBox(3,1))RETURN
         IF(R(3)>Box%BndBox(3,2))RETURN
         PointOutSideBox=.FALSE.
      END FUNCTION PointOutSideBox
!===============================================================================
!     Determine if a BBox A with extent overlaps with BBox B
!===============================================================================
      FUNCTION BoxOutSideBox(A,B)
         LOGICAL                  :: BoxOutSideBox
         TYPE(BBox),INTENT(IN)    :: A,B
         REAL(DOUBLE)             :: Tx,Ty,Tz
         BoxOutSideBox=.TRUE.
         Tx=ABS(A%Center(1)-B%Center(1))
         IF(Tx>A%Half(1)+B%Half(1))RETURN
         Ty=ABS(A%Center(2)-B%Center(2))
         IF(Ty>A%Half(2)+B%Half(2))RETURN
         Tz=ABS(A%Center(3)-B%Center(3))
         IF(Tz>A%Half(3)+B%Half(3))RETURN
         BoxOutSideBox=.FALSE.
     END FUNCTION BoxOutSideBox
!============================================================================
!     Generate an expanded BBox
!============================================================================
      FUNCTION ExpandPoint(Point,Extent,Box_O) RESULT(Expando)
         Type(BBox)                 :: Expando
         Type(BBox),OPTIONAL        :: Box_O
         REAL(DOUBLE),DIMENSION(3)  :: Point
         Real(Double)               :: Extent
         Integer                    :: I

         IF(PRESENT(Box_O))THEN
            Expando%Tier=Box_O%Tier
            Expando%Number=Box_O%Number
         ENDIF
         Expando%BndBox(1:3,1)=Point-(/Extent,Extent,Extent/)
         Expando%BndBox(1:3,2)=Point+(/Extent,Extent,Extent/)
         Expando%Half  (1:3)  =Half*(Expando%BndBox(1:3,2)-Expando%BndBox(1:3,1))
         Expando%Center(1:3)  =Half*(Expando%BndBox(1:3,2)+Expando%BndBox(1:3,1))
       END FUNCTION ExpandPoint

      FUNCTION ExpandBox(Box,Extent) RESULT(Expando)
         Type(BBox)     :: Box,Expando
         Real(Double)   :: Extent
         Integer        :: I
         Expando%Tier=Box%Tier
         Expando%Number=Box%Number
         Expando%BndBox(1:3,1)=Box%BndBox(1:3,1)-Extent
         Expando%BndBox(1:3,2)=Box%BndBox(1:3,2)+Extent
         Expando%Half  (1:3)  =Half*(Expando%BndBox(1:3,2)-Expando%BndBox(1:3,1))
         Expando%Center(1:3)  =Half*(Expando%BndBox(1:3,2)+Expando%BndBox(1:3,1))
      END FUNCTION ExpandBox

      FUNCTION BoxVolume(Box) RESULT(Vol)
         TYPE(BBox)     :: Box
         REAL(Double)   :: Vol
         Vol=(Box%BndBox(1,2)-Box%BndBox(1,1)) &
            *(Box%BndBox(2,2)-Box%BndBox(2,1)) &
            *(Box%BndBox(3,2)-Box%BndBox(3,1))
       END FUNCTION BoxVolume
       !
       FUNCTION BoxMargin(Box) RESULT(Mrg)
         TYPE(BBox)     :: Box
         REAL(Double)   :: Mrg
         Mrg=(Box%BndBox(1,2)-Box%BndBox(1,1)) &
            +(Box%BndBox(2,2)-Box%BndBox(2,1)) &
            +(Box%BndBox(3,2)-Box%BndBox(3,1))
       END FUNCTION BoxMargin
       !
       FUNCTION BoxOverlap(Left,Right) RESULT(Vol)
         TYPE(BBox)     :: Left,Right
         REAL(Double)   :: Vol,Hi,Low
         INTEGER        :: I
         IF(BoxOutSideBox(Left,Right))THEN
            Vol=Zero
         ELSE
            Vol=One
            DO I=1,3
               IF(Left%BndBox(I,1)<Right%BndBox(I,1))THEN
                  Low=Right%BndBox(I,1)
               ELSE
                  Low=Left%BndBox(I,1)
               ENDIF
               IF(Left%BndBox(I,2)<Right%BndBox(I,2))THEN
                  Hi=Right%BndBox(I,2)
               ELSE
                  Hi=Left%BndBox(I,2)
               ENDIF
               Vol=Vol*(Hi-Low)
            ENDDO
         ENDIF
       END FUNCTION BoxOverlap
!==========================================================================
!
!==========================================================================
      SUBROUTINE PrintBBox(Box,U,ZetaMin_O,MaxAmp_O)
         TYPE(BBox)             :: Box
         REAL(DOUBLE), OPTIONAL :: ZetaMin_O,MaxAmp_O
         INTEGER    :: U,Number,Tier,I,J
         WRITE(U,*)'======================================================='
         WRITE(U,33)Box%Number,Box%Tier
         33 FORMAT(' Number = ',I4,' Tier = ',I3)
         IF(PRESENT(ZetaMin_O).AND.PRESENT(MaxAmp_O)) &
            WRITE(22,44)ZetaMin_O,MaxAmp_O
         44 FORMAT(' ZetaMin = ',F12.6,' MaxAmp = ',F12.6)
         DO I=1,3
            WRITE(U,55)I,Box%BndBox(I,1),Box%BndBox(I,2),Box%Center(I),Box%Half(I)
         ENDDO
      55 FORMAT(I2,' Low=',F10.5,', Hi =',F10.5,', Center = ',F10.5,', Half = ',F10.5)
      END SUBROUTINE PrintBBox
END MODULE BoundingBox
