!
!    SORTING AND ORDERING (VIA SPACE FILLING CURVES)
!    Author: Matt Challacombe
MODULE TravelMan
  USE DerivedTypes; USE MemMan
  IMPLICIT NONE
  CONTAINS
  
  SUBROUTINE Anneal(XCoor,YCoor,ZCoor,NCity,TravOrder)
  IMPLICIT NONE
  TYPE(DBL_VECT) :: XCoor,YCoor,ZCoor,XX,YY,ZZ
  INTEGER::LeftCity,RightCity,NCity,I,J,K,NOver,NLimit,NSucc
  TYPE(INT_VECT)::TravOrder,ChosenCity,NewTravOrder
  REAL(DOUBLE) :: TFacR,T,MaxIncCost,TotCost1,TotCost,Dec,IncCost
  REAL(DOUBLE),EXTERNAL::Random
  LOGICAL::Ans

  IF(NCity <= 4) THEN 
    STOP 'Error in Anneal! Too few cities!'
  ENDIF

  CALL GetEndCities1(XCoor,YCoor,ZCoor,NCity,LeftCity,RightCity)
  DO I = 1, NCity
    TravOrder%I(I) = I
  ENDDO  
  TravOrder%I(1) = LeftCity
  TravOrder%I(LeftCity) = 1
  TravOrder%I(NCity) = RightCity
  TravOrder%I(RightCity) = NCity

  ! Calculate the length of the total path
  TotCost = CalTotCost(XCoor,YCoor,ZCoor,TravOrder,NCity)
  WRITE(*,*) 'TotCost = ',TotCost

  CALL NEW(ChosenCity,6)
  CALL NEW(XX,6); CALL NEW(YY,6); CALL NEW(ZZ,6)
  CALL NEW(NewTravOrder,NCity)

  MaxIncCost = -BIG_DBL
  DO I = 1,5000
    Dec = Random()
    WRITE(*,*) 'Dec =', Dec
    IF(Dec < 0.5) THEN
      CALL RevCst(XCoor,YCoor,ZCoor,NCity,TravOrder,&
        ChosenCity,XX,YY,ZZ,IncCost)
      IF(IncCost > MaxIncCost) THEN
        MaxIncCost = IncCost
      ENDIF
      CALL Revers(TravOrder,ChosenCity)
      TotCost = TotCost + IncCost
  
      TotCost1 = CalTotCost(XCoor,YCoor,ZCoor,TravOrder,NCity)
      WRITE(*,*) 'Rev: TotCost,TotCost1 =',TotCost,TotCost1
      IF(abs(TotCost-TotCost1) > 1.0e-7) THEN
        STOP 'TotCost and TotCost1 differ!!'
      ENDIF
    ELSE
      CALL TrnCst(XCoor,YCoor,ZCoor,NCity,TravOrder,ChosenCity,&
        XX,YY,ZZ,IncCost)
      IF(IncCost > MaxIncCost) THEN
        MaxIncCost = IncCost
      ENDIF
      CALL Trnspt(TravOrder,ChosenCity,NCity,NewTravOrder)
      TotCost = TotCost + IncCost
      TotCost1 = CalTotCost(XCoor,YCoor,ZCoor,TravOrder,NCity)
      WRITE(*,*) 'Trn: TotCost,TotCost1 =',TotCost,TotCost1
      IF(abs(TotCost-TotCost1) > 1.0e-7) THEN
        STOP 'TotCost and TotCost1 differ!!'
      ENDIF
    ENDIF
  ENDDO      
  WRITE(*,*) 'NCity,MaxIncCost =',NCity,MaxIncCost

  NOver = 100*NCity
  NLimit = 10*NCity
  TFacR = 0.9
  T = MaxIncCost*0.5
  
  DO J = 1, 100
    NSucc = 0
    DO K = 1,NOver
      Dec = Random() 
      IF(Dec < 0.5) THEN
        CALL RevCst(XCoor,YCoor,ZCoor,NCity,TravOrder,&
          ChosenCity,XX,YY,ZZ,IncCost)
        Ans = (IncCost < 0.0) .OR. (Random().LT.EXP(-IncCost/T))
        IF(Ans) THEN
          NSucc = NSucc + 1
          TotCost = TotCost + IncCost
          CALL Revers(TravOrder,ChosenCity)
        ENDIF
      ELSE
        CALL TrnCst(XCoor,YCoor,ZCoor,NCity,TravOrder,ChosenCity,&
          XX,YY,ZZ,IncCost)
        Ans = (IncCost < 0.0) .OR. (Random().LT.EXP(-IncCost/T))
        IF(Ans) THEN
          NSucc = NSucc + 1
          TotCost = TotCost + IncCost
          CALL Trnspt(TravOrder,ChosenCity,NCity,NewTravOrder)
        ENDIF
      ENDIF
      IF(NSucc .ge. NLimit) Exit  
    ENDDO
    T = T*TFacR
    write(*,*) 'T = ', T
    write(*,*) 'NSucc = ',NSucc
    IF(NSucc == 0) EXIT
  ENDDO
  write(*,*) 'temperature index J = ', J
  write(*,*) 'success number index K = ', K  
  TotCost1 = CalTotCost(XCoor,YCoor,ZCoor,TravOrder,NCity)
  WRITE(*,*) 'Final cost: TotCost,TotCost1 =',TotCost,TotCost1
  IF(ABS(TotCost-TotCost1) > 1.0e-7) THEN
    STOP 'TotCost and TotCost1 differ!!'
  ENDIF
  
  END SUBROUTINE Anneal

!---------------------------------------------------------------------  
  SUBROUTINE TrnCst(XCoor,YCoor,ZCoor,NCity,TravOrder,ChosenCity,&
    XX,YY,ZZ,IncCost)
    TYPE(DBL_VECT)::XCoor,YCoor,ZCoor,XX,YY,ZZ
    INTEGER::Val,IP,NCity,ArrLen,I,J,RanIndex
    TYPE(INT_VECT)::TravOrder,SubscrArr,ChosenCity
    REAL(DOUBLE)::IncCost
    REAL(DOUBLE),EXTERNAL::Random
 
    ArrLen = NCity-2 ! the first and the last cities are fixed
    DO

!   DO I = 1, ArrLen
!     SubscrArr%I(I) = I+1
!   ENDDO
!   DO I = 1, 3
!     RanIndex = 1 + ArrLen*Random()
!     ChosenCity%I(I) = SubscrArr%I(RanIndex)
!     ! remove RanIndex from SubscrArr
!     SubscrArr%I(RanIndex) = SubscrArr%I(ArrLen)
!     ArrLen = ArrLen-1
!   ENDDO
      ChosenCity%I(1) = 2+ArrLen*Random()
      ChosenCity%I(2) = 2+ArrLen*Random()
      ChosenCity%I(3) = 2+ArrLen*Random()

      DO I = 2, 3
        Val = ChosenCity%I(I); IP = I
        DO 
          IF(.NOT.(IP > 1 .AND. Val < ChosenCity%I(IP-1))) EXIT
          ChosenCity%I(IP) = ChosenCity%I(IP-1)
          IP = IP-1
        ENDDO
        ChosenCity%I(IP) = Val
      ENDDO

      IF(ChosenCity%I(1) /= ChosenCity%I(2) .AND. &
         ChosenCity%I(2) /= ChosenCity%I(3)) EXIT
    ENDDO
    ChosenCity%I(4) = ChosenCity%I(3)+1
    ChosenCity%I(5) = ChosenCity%I(1)-1
    ChosenCity%I(6) = ChosenCity%I(2)+1

!   DO I = 1, 6
!     IF(ChosenCity%I(I) < 1 .OR. ChosenCity%I(I) > NCity) THEN
!       WRITE(*,*) 'Wrong indices in TrnCst ! Something is wrong !'
!       STOP
!     ENDIF
!   ENDDO

    DO I = 1, 6
      J = TravOrder%I(ChosenCity%I(I))
      XX%D(I) = XCoor%D(J)
      YY%D(I) = YCoor%D(J)
      ZZ%D(I) = ZCoor%D(J)
    ENDDO
    IncCost = -ALen(xx%D(2),yy%D(2),zz%D(2),xx%D(6),yy%D(6),zz%D(6)) &
              -ALen(xx%D(1),yy%D(1),zz%D(1),xx%D(5),yy%D(5),zz%D(5)) &
              -ALen(xx%D(3),yy%D(3),zz%D(3),xx%D(4),yy%D(4),zz%D(4)) &
              +ALen(xx%D(1),yy%D(1),zz%D(1),xx%D(3),yy%D(3),zz%D(3)) &
              +ALen(xx%D(2),yy%D(2),zz%D(2),xx%D(4),yy%D(4),zz%D(4)) &
              +ALen(xx%D(5),yy%D(5),zz%D(5),xx%D(6),yy%D(6),zz%D(6))

  END SUBROUTINE TrnCst

!---------------------------------------------------------------------  
  SUBROUTINE Trnspt(TravOrder,ChosenCity,NCity,NewTravOrder)
  TYPE(INT_VECT)::TravOrder,ChosenCity,NewTravOrder
  INTEGER::NCity,NN,I,J,INDEX

  NN = ChosenCity%I(5)
  INDEX = 1
  DO I = 1, NN
    NewTravOrder%I(INDEX) = TravOrder%I(I)
    INDEX = INDEX + 1
  ENDDO

  NN = ChosenCity%I(3)-ChosenCity%I(6)+1
  J = ChosenCity%I(6)
  DO I = 1,NN
    NewTravOrder%I(INDEX) = TravOrder%I(J+I-1)
    INDEX = INDEX + 1
  ENDDO

  NN = ChosenCity%I(2)-ChosenCity%I(1)+1
  J = ChosenCity%I(1)
  DO I = 1,NN
    NewTravOrder%I(INDEX) = TravOrder%I(J+I-1)
    INDEX = INDEX + 1
  ENDDO

  NN = NCity-ChosenCity%I(4)+1
  J = ChosenCity%I(4)
  DO I = 1, NN
    NewTravOrder%I(INDEX) = TravOrder%I(J+I-1)
    INDEX = INDEX + 1
  ENDDO

  IF(INDEX /= (NCity+1)) THEN
    STOP 'missing cities !!'
  ENDIF  

  DO I = 1, NCity
    TravOrder%I(I) = NewTravOrder%I(I)
  ENDDO    
  END SUBROUTINE TRNSPT
  
!----------------------------------------------------------------------
  SUBROUTINE RevCst(XCoor,YCoor,ZCoor,NCity,TravOrder,&
    ChosenCity,XX,YY,ZZ,IncCost)
    TYPE(DBL_VECT)::XCoor,YCoor,ZCoor,XX,YY,ZZ
    INTEGER::TmpInt,NCity,ArrLen,I,J,RanIndex1,RanIndex2
    TYPE(INT_VECT)::TravOrder,SubscrArr,ChosenCity
    REAL(DOUBLE)::IncCost
    REAL(DOUBLE),EXTERNAL::Random
 
!   ArrLen = NCity-2 ! the first and the last cities are fixed
!   DO I = 1, ArrLen
!     SubscrArr%I(I) = I+1
!   ENDDO
!   DO I = 1, 2
!     RanIndex = 1 + ArrLen*Random()
!     ChosenCity%I(I) = SubscrArr%I(RanIndex)
!     ! remove RanIndex from SubscrArr
!     SubscrArr%I(RanIndex) = SubscrArr%I(ArrLen)
!     ArrLen = ArrLen-1
!   ENDDO
!   ! WRITE(*,*) 'Chosen Cities : ',ChosenCity%I(1),ChosenCity%I(2)

    ArrLen = NCity-2
    DO 
      RanIndex1 = 2 + ArrLen*Random()
      RanIndex2 = 2 + ArrLen*Random()
      IF(RanIndex1 /= RanIndex2) EXIT
    ENDDO
    IF(RanIndex1 > RanIndex2) THEN
      TmpInt = RanIndex1
      RanIndex1 = RanIndex2
      RanIndex2 = TmpInt
    ENDIF

    ChosenCity%I(1) = RanIndex1
    ChosenCity%I(2) = RanIndex2

    ChosenCity%I(3) = ChosenCity%I(1)-1
    ChosenCity%I(4) = ChosenCity%I(2)+1

    DO I = 1, 4
      IF(ChosenCity%I(I) < 1 .OR. ChosenCity%I(I) > NCity) THEN
        WRITE(*,*) 'Wrong indices in RevCst ! Something is wrong !'
        STOP
      ENDIF
    ENDDO
    
    DO I = 1, 4
      J = TravOrder%I(ChosenCity%I(I))
      XX%D(I) = XCoor%D(J)
      YY%D(I) = YCoor%D(J)
      ZZ%D(I) = ZCoor%D(J)
    ENDDO

    IncCost = -ALen(xx%D(1),yy%D(1),zz%D(1),xx%D(3),yy%D(3),zz%D(3)) &
              -ALen(xx%D(2),yy%D(2),zz%D(2),xx%D(4),yy%D(4),zz%D(4)) &
              +ALen(xx%D(1),yy%D(1),zz%D(1),xx%D(4),yy%D(4),zz%D(4)) &
              +ALen(xx%D(2),yy%D(2),zz%D(2),xx%D(3),yy%D(3),zz%D(3))
  END SUBROUTINE RevCst

!----------------------------------------------------------------------
  SUBROUTINE Revers(TravOrder,ChosenCity)
  INTEGER::NN,I,TmpInt,Left,Right
  TYPE(INT_VECT)::TravOrder,ChosenCity
  
  NN = (ChosenCity%I(2)-ChosenCity%I(1)+1)/2
  DO I = 1,NN
    Left = ChosenCity%I(1)+I-1
    Right = ChosenCity%I(2)-I+1
    TmpInt = TravOrder%I(Left)
    TravOrder%I(Left) = TravOrder%I(Right)
    TravOrder%I(Right) = TmpInt
  ENDDO
  END SUBROUTINE Revers

!----------------------------------------------------------------------
  ! the order is coordinates of 1, and followed by 2
  FUNCTION ALen(x1,y1,z1,x2,y2,z2)
    REAL(DOUBLE)::x1,y1,z1,x2,y2,z2,ALen,dx,dy,dz
    dx = x1-x2; dy = y1-y2; dz = z1-z2
    ALen = SQRT(dx*dx + dy*dy + dz*dz)
  END FUNCTION ALen

!----------------------------------------------------------------------
  FUNCTION CalTotCost(XCoor,YCoor,ZCoor,TravOrder,NCity)
  TYPE(DBL_VECT)::XCoor,YCoor,ZCoor
  TYPE(INT_VECT)::TravOrder
  INTEGER::NCity,I,J,K
  REAL(DOUBLE)::CalTotCost
  
  CalTotCost = 0.0
  DO I = 1,NCity-1
    J = TravOrder%I(I)
    K = TravOrder%I(I+1)
    CalTotCost = CalTotCost + ALen(XCoor%D(J),YCoor%D(J),ZCoor%D(J),&
      XCoor%D(K),YCoor%D(K),ZCoor%D(K))
  ENDDO
  END FUNCTION CalTotCost

!----------------------------------------------------------------------
  ! get the end cities using the tight bounding box
  SUBROUTINE GetEndCities1(XCoor,YCoor,ZCoor,NCity,LeftCity,RightCity)
  TYPE(DBL_VECT)::XCoor,YCoor,ZCoor
  REAL(DOUBLE)::Dis,MinLDis,MinRDis,X,Y,Z,XLEx,YLEx,ZLEx,XREx,YREx,ZREx
  INTEGER::NCity,I,LeftCity,RightCity
  
  XLEx = Big_DBL;  XREx = -Big_DBL
  YLEx = Big_DBL;  YREx = -Big_DBL
  ZLEx = Big_DBL;  ZREx = -Big_DBL
 
  DO I = 1, NCity
    X = XCoor%D(I); Y = YCoor%D(I); Z = ZCoor%D(I)
    IF(X < XLEx) XLEx = X; IF(X > XREx) XREx = X
    IF(Y < YLEx) YLEx = Y; IF(Y > YREx) YREx = Y
    IF(Z < ZLEx) ZLEx = Z; IF(Z > ZREx) ZREx = Z
  ENDDO
  WRITE(*,*) 'XEx : ',XLEx, XREx
  WRITE(*,*) 'YEx : ',YLEx, YREx
  WRITE(*,*) 'ZEx : ',ZLEx, ZREx
    
  MinLDis = Big_DBL; MinRDis = Big_DBL
 
  DO I = 1, NCity
    Dis = ALen(XCoor%D(I),YCoor%D(I),ZCoor%D(I),XLEx,YLEx,ZLEx)
    IF(Dis < MinLDis) THEN
      MinLDis = Dis
      LeftCity = I
    ENDIF
    Dis = ALen(XCoor%D(I),YCoor%D(I),ZCoor%D(I),XREx,YREx,ZREx)
    IF(Dis < MinRDis) THEN
      MinRDis = Dis
      RightCity = I
    ENDIF
  ENDDO

  WRITE(*,*) 'MinLDis = ',MinLDis
  WRITE(*,*) 'MinRDis = ',MinRDis
  WRITE(*,*) 'LeftCity,RightCity = ',LeftCity,RightCity
  END SUBROUTINE GetEndCities1

!----------------------------------------------------------------------
  ! find out the largest end-to-end distance
  SUBROUTINE GetEndCities2(XCOOR,YCOOR,ZCoor,NCity,LeftCity,RightCity)
  TYPE(DBL_VECT) :: XCoor,YCoor,ZCoor
  INTEGER::I,J,Ncity,leftcity,rightcity
  REAL(DOUBLE)::DisIJ,MaxDis

  MaxDis = 0.0
  DO I = 1, NCity
    DO J = I+1, NCity
      DisIJ = ALen(XCoor%D(I),YCoor%D(I),ZCoor%D(I),&
        XCoor%D(J),YCoor%D(J),ZCoor%D(J))
      IF(DisIJ .GE. MaxDis) THEN
        MaxDis = DisIJ
        LeftCity = I
        RightCity = J
      ENDIF
    ENDDO
  ENDDO
  END SUBROUTINE GetEndCities2

END MODULE TravelMan

!---------------------------------------------------------------------------
MODULE AnnealMap
  USE PrettyPrint
  IMPLICIT NONE
  INTEGER,PRIVATE::NCity
CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE TableAnneal(RCoor,NCityV,TravO)
    TYPE(DBL_RNK2)::RCoor
    INTEGER(INT8),ALLOCATABLE::CityRank(:)
    TYPE(INT_VECT)::TravO
    TYPE(INT_RNK3)::AnnealKey
    INTEGER::Imax,IntVect(3),NCityV,I,J,K,LinDim,LM1,ReadNCity,M,Rank,IndexInt
    REAL(DOUBLE)::Ratio,RMin(3),Diff,MaxDiff
    INTEGER,PARAMETER::ReadU=30
    TYPE(INT_VECT)::ReadTravO
   
  
    NCity = NCityV
    CALL OpenASCII('Final_TravO.dat',ReadU,OldFileQ_O=.TRUE.,Rewind_O=.TRUE.)
    Read(ReadU,*) ReadNCity
    WRITE(*,*) 'TableAnneal, ReadNCity = ',ReadNCity
    LinDim = NINT(ReadNCity**(1.0d0/3.0d0))
    WRITE(*,*) 'LinDim = ', LinDim
    IF(LinDim*LinDim*LinDim /= ReadNCity) THEN
      STOP 'Error: Cube root problem in TableAnneal!'
    ENDIF
    LM1 = LinDim-1
    WRITE(*,*) 'LinDim = ', LinDim, ', Lm1=', LM1
    CALL New(AnnealKey,(/LM1,LM1,LM1/),(/0,0,0/))
    CALL New(ReadTravO,ReadNCity)

    Read(ReadU,*) (ReadTravO%I(I),I=1,ReadNCity)
!   WRITE(*,*) (ReadTravO%I(I),I=1,ReadNCity)
    IndexInt = 0
    DO I = 0,LM1
      DO J = 0, LM1
        DO K = 0, LM1
          IndexInt = IndexInt + 1
          Rank = 0
          DO M = 1, ReadNCity
            IF(ReadTravO%I(M) == IndexInt) THEN
              Rank = M
              EXIT
            ENDIF
          ENDDO
          IF(Rank < 1 .OR. Rank > ReadNCity) THEN
            WRITE(*,*) 'IndexInt = ',IndexInt
            WRITE(*,*) 'Rank = ',Rank
            STOP 'Error: Some serious problem in inverse mapping!'
          ENDIF
          AnnealKey%I(K,J,I) = Rank
        ENDDO
      ENDDO
    ENDDO
  
    DO I = 0,LM1
      DO J = 0, LM1
        DO K = 0, LM1
!         WRITE(*,*) 'AnnealKey',K,J,I,AnnealKey%I(K,J,I)
        ENDDO
      ENDDO
    ENDDO
    CLOSE(ReadU)
  
    RMin(1) = Big_DBL
    RMin(2) = Big_DBL
    RMin(3) = Big_DBL
    DO I = 1, NCity
      DO J = 1,3
        RMin(J) = DMin1(RMin(J),RCoor%D(J,I))
      ENDDO
      TravO%I(I) = I
    ENDDO
    ! any -ve value for MaxDiff is okay, since we deal with the 1st quadrant
    MaxDiff = -1.0D0 
    DO I = 1,NCity
      DO J = 1,3
        Diff = RCoor%D(J,I)-RMin(J)
        MaxDiff = DMax1(MaxDiff,Diff)
      ENDDO
    ENDDO
    RATIO = (LinDim*1.D0)/MaxDiff
    IMax = -Big_Int
    ALLOCATE(CityRank(1:NCity))
    DO I = 1,NCity
      DO J = 1, 3
        ! IntVect(J) = DNINT( (RCoor%D(J,I)-RMin(J))*Ratio)
        IntVect(J) =  (RCoor%D(J,I)-RMin(J))*Ratio ! truncate the fraction part.
        IMax = Max(IMax,IntVect(J))
        IF(IntVect(J) == LinDim) THEN
          IntVect(J) = IntVect(J)-1
        ENDIF
        IF(IntVect(J) < 0 .OR. IntVect(J) .GE. LinDim) THEN
          WRITE(*,*) 'IntVect',j, ' =',IntVect(j)
          STOP 'ERROR: An error has occured! array bound problem'
        ENDIF
      ENDDO
      CityRank(I) = AnnealKey%I(IntVect(1),IntVect(2),IntVect(3))
!     WRITE(*,*) 'CityRank',I, ' =', CityRank(I)
    ENDDO
    IF(IMax /= LinDim) THEN
      STOP 'ERROR: Numerical accuracy problem in TableAnneal!'
    ENDIF
    CALL I8Sort(CityRank,TravO%I,NCity,2)
    CALL Delete(AnnealKey)
    CALL Delete(ReadTravO)

  END SUBROUTINE TableAnneal
END MODULE AnnealMap


MODULE Order
   USE DerivedTypes
   USE GlobalCharacters
   USE GlobalScalars
   USE GlobalObjects
   USE ProcessControl
   USE MemMan
   USE TravelMan ; USE AnnealMap
   USE Parse
   IMPLICIT NONE
   INTERFACE Sort
      MODULE PROCEDURE Sort_DBL_INT,Sort_INT_INT, Sort_INT_VECT
   END INTERFACE
   INTERFACE Random
      MODULE PROCEDURE RANDOM_INT,RANDOM_DBL
   END INTERFACE
   INTERFACE 
      SUBROUTINE DblIntSort77(N,X,Y,Ordr)
         USE DerivedTypes
         INTEGER,                  INTENT(IN)    :: N,Ordr
         REAL(DOUBLE),DIMENSION(N),INTENT(INOUT) :: X
         INTEGER,     DIMENSION(N),INTENT(INOUT) :: Y
      END SUBROUTINE
      SUBROUTINE SFCOrder77(N,R,Point,Key,Hilbert)
         USE DerivedTypes
         INTEGER,                     INTENT(IN)    :: N
         LOGICAL,                     INTENT(IN)    :: Hilbert
         REAL(DOUBLE), DIMENSION(3,N),INTENT(INOUT) :: R 
         INTEGER(INT8),DIMENSION(N),  INTENT(INOUT) :: Key 
         INTEGER,      DIMENSION(N),  INTENT(INOUT) :: Point
      END SUBROUTINE

      FUNCTION Interleave(Ix,Iy,Iz)       
         IMPLICIT NONE
         INTEGER,INTENT(IN) :: Ix,Iy,Iz
         INTEGER, PARAMETER :: INT8=SELECTED_INT_KIND(18)    
         INTEGER(INT8)            :: Interleave
       END FUNCTION Interleave
   END INTERFACE
   CONTAINS
!
!      FUNCTION RANDOM_INT(Limits)
!         INTEGER                :: RANDOM_INT
!         INTEGER, DIMENSION(2)  :: Limits
!         REAL(DOUBLE)           :: Delta
!         REAL(DOUBLE), EXTERNAL :: Random
!         Delta=DBLE(Limits(2)-Limits(1)+1)
!         RANDOM_INT=Limits(1)+INT(Delta*Random())
!         RANDOM_INT=Limits(1)+INT(Delta*RAND())
!      END FUNCTION RANDOM_INT
!
      FUNCTION RANDOM_DBL(Limits)
         REAL(DOUBLE)              :: RANDOM_DBL
         REAL(DOUBLE),DIMENSION(2) :: Limits
         REAL(DOUBLE)              :: Delta
         REAL(DOUBLE), EXTERNAL    :: Random
         Delta=Limits(2)-Limits(1)+0.0D0
         RANDOM_DBL=Limits(1)+Delta*Random()
!         RANDOM_DBL=Limits(1)+Delta*Rand()
      END FUNCTION RANDOM_DBL
!
      FUNCTION RANDOM_INT(Limits)
         INTEGER               :: RANDOM_INT,Delta
         INTEGER, SAVE         :: JRan=10408
         INTEGER, DIMENSION(2) :: Limits
         INTEGER, PARAMETER    :: Im=259200,Ia=7141,Ic=54773
         JRan=MOD(JRan*Ia+Ic,Im)
         Delta=Limits(2)-Limits(1)+1
         RANDOM_INT=Limits(1)+(Delta*JRan)/Im
         IF(RANDOM_INT>Limits(2).OR.RANDOM_INT<Limits(1)) &
            CALL Halt(' Limits hosed in RANDOM_INT ')
      END FUNCTION RANDOM_INT
!--------------------------------------------------------------
!    F90 wrapper for SFCOrder77, which circumvents the lack
!    of INTEGER(KIND=8) (INTEGER*8) support for cheazy 
!    F90 compilers (pgf,nag...)
!
     SUBROUTINE SFCOrder(N,R,Point,SFC_KEY)
        INTEGER,        INTENT(IN)    :: N,SFC_KEY
        TYPE(DBL_RNK2), INTENT(INOUT) :: R
        TYPE(INT_VECT), INTENT(INOUT) :: Point
        INTEGER(INT8),ALLOCATABLE, &
                         DIMENSION(:) :: IKey
        TYPE(DBL_VECT)                :: RKey
        INTEGER                       :: I
        TYPE(DBL_VECT)                :: XCoor,YCoor,ZCoor
        
!
        IF(SFC_KEY==SFC_RANDOM)THEN
           CALL New(RKey,N)
           DO I=1,N
              Point%I(I)=I
              CALL RANDOM_NUMBER(RKey%D(I))
           ENDDO            
           CALL Sort_DBL_INT(RKey,Point,N)
           CALL Delete(RKey)           
        ELSEIF(SFC_KEY==SFC_PEANO)THEN
           ALLOCATE(IKey(N))
           CALL SFCOrder77(N,R%D,Point%I,IKey,.FALSE.)
           DEALLOCATE(IKey)
        ELSEIF(SFC_KEY==SFC_HILBERT)THEN
           ALLOCATE(IKey(N))
           CALL SFCOrder77(N,R%D,Point%I,IKey,.TRUE.)
           DEALLOCATE(IKey)
        ELSEIF(SFC_KEY==SFC_TRAVEL)THEN
           CALL NEW(XCoor,N)
           CALL NEW(YCoor,N)
           CALL NEW(ZCoor,N)
           DO I = 1,N
             XCoor%D(I) = R%D(1,I)
             YCoor%D(I) = R%D(2,I)
             ZCoor%D(I) = R%D(3,I)
           ENDDO
           CALL Anneal(XCoor,YCoor,ZCoor,N,Point)
        ELSEIF(SFC_KEY==SFC_TABLETRAV)THEN
          CALL TableAnneal(R,N,Point)
        ELSE
           CALL MondoHalt(-100,'Bad SFC_Key in SFCOrder')
        ENDIF
     END SUBROUTINE SFCOrder

     SUBROUTINE Sort_DBL_INT(X,Y,N_O,Ordr_O)
        TYPE(DBL_VECT), INTENT(INOUT) :: X
        TYPE(INT_VECT), INTENT(INOUT) :: Y
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-2
        N=MIN(SIZE(X%D),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL DblIntSort77(N,X%D,Y%I,Ordr)
    END SUBROUTINE Sort_DBL_INT


     SUBROUTINE Sort_INT_INT(X,Y,N_O,Ordr_O)
        TYPE(INT_VECT), INTENT(INOUT) :: X
        TYPE(INT_VECT), INTENT(INOUT) :: Y
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-1
        N=MIN(SIZE(X%I),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntIntSort77(N,X%I,Y%I,Ordr)
    END SUBROUTINE Sort_INT_INT


     SUBROUTINE Sort_INT_VECT(X,N_O,Ordr_O)
        TYPE(INT_VECT), INTENT(INOUT) :: X
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-1
        N=SIZE(X%I)
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntSort77(N,X%I,Ordr)
    END SUBROUTINE Sort_INT_VECT


END MODULE
