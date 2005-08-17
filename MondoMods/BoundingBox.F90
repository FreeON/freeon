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
!==========================================================================
!
!==========================================================================
      SUBROUTINE PrintBBox(Box,ZetaMin_O,MaxAmp_O)
         TYPE(BBox)             :: Box
         REAL(DOUBLE), OPTIONAL :: ZetaMin_O,MaxAmp_O
         INTEGER    :: Number,Tier,I,J
         WRITE(22,*)'======================================================='
         WRITE(22,33)Box%Number,Box%Tier
         33 FORMAT(' Number = ',I4,' Tier = ',I3)
         IF(PRESENT(ZetaMin_O).AND.PRESENT(MaxAmp_O)) &
            WRITE(22,44)ZetaMin_O,MaxAmp_O
         44 FORMAT(' ZetaMin = ',F12.6,' MaxAmp = ',F12.6)
         DO I=1,3
            WRITE(22,55)I,Box%BndBox(I,1),Box%BndBox(I,2),Box%Center(I),Box%Half(I)
         ENDDO
      55 FORMAT(I2,' Low=',F10.5,', Hi =',F10.5,', Center = ',F10.5,', Half = ',F10.5)
      END SUBROUTINE PrintBBox
END MODULE BoundingBox
