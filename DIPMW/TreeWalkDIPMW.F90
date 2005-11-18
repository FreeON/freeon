MODULE TreeWalkDIPMW
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE McMurchie
  USE InOut
  USE Thresholding
  USE AtomPairs
  USE BraBloks
  USE DIPMWThresholds
  USE DIPMWTree
  IMPLICIT NONE
  LOGICAL PrintFlag
!-------------------------------------------------------------------------------------
! Globals
  REAL(DOUBLE),DIMENSION(HGLen)  :: Ket
  REAL(DOUBLE)                   :: PExtent
  TYPE(BBox)                     :: PBox
  TYPE(PrimPair)                 :: Prim
  CONTAINS 
!-------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE KxcWalkDIPMW(Node)
    TYPE(PIGWNode), POINTER   :: Node
    REAL(DOUBLE)              :: PQx,PQy,PQz
!   On Null box return
!    IF(Node%Null) RETURN
!   Does PBox intersect this PIGW BBox? 
    PQx=PBox%Center(1)-Node%Box%Center(1)
    IF(ABS(PQx)>PBox%Half(1)+Node%Box%Half(1))RETURN
    PQy=PBox%Center(2)-Node%Box%Center(2)
    IF(ABS(PQy)>PBox%Half(2)+Node%Box%Half(2))RETURN
    PQz=PBox%Center(3)-Node%Box%Center(3)
    IF(ABS(PQz)>PBox%Half(3)+Node%Box%Half(3))RETURN
!   Compute the Integral
    SELECT CASE(Prim%Ell)
    CASE(0)
!!$       Ket(1) =Ket(1) +Node%VxcCoef*IntPIGW1(Node%Level,Prim%Zeta,PQx,PQy,PQz)
    CASE(1)
!!$       Ket(1) =Ket(1) +Node%VxcCoef*IntPIGW1(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(2) =Ket(2) +Node%VxcCoef*IntPIGW2(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(3) =Ket(3) +Node%VxcCoef*IntPIGW3(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(4) =Ket(4) +Node%VxcCoef*IntPIGW4(Node%Level,Prim%Zeta,PQx,PQy,PQz)
    CASE(2)       
!!$       Ket(1) =Ket(1) +Node%VxcCoef*IntPIGW1(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(2) =Ket(2) +Node%VxcCoef*IntPIGW2(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(3) =Ket(3) +Node%VxcCoef*IntPIGW3(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(4) =Ket(4) +Node%VxcCoef*IntPIGW4(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(5) =Ket(5) +Node%VxcCoef*IntPIGW5(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(6) =Ket(6) +Node%VxcCoef*IntPIGW6(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(7) =Ket(7) +Node%VxcCoef*IntPIGW7(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(8) =Ket(8) +Node%VxcCoef*IntPIGW8(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(9) =Ket(9) +Node%VxcCoef*IntPIGW9(Node%Level,Prim%Zeta,PQx,PQy,PQz)
!!$       Ket(10)=Ket(10)+Node%VxcCoef*IntPIGW10(Node%Level,Prim%Zeta,PQx,PQy,PQz)
    END SELECT
    IF(Node%LeafType=='EnddLeaf') RETURN
!   Continue Walking
!!$    CALL  KxcWalkPIGW(Node%Direc000)
!!$    CALL  KxcWalkPIGW(Node%Direc100)
!!$    CALL  KxcWalkPIGW(Node%Direc010)
!!$    CALL  KxcWalkPIGW(Node%Direc001)
!!$    CALL  KxcWalkPIGW(Node%Direc011)
!!$    CALL  KxcWalkPIGW(Node%Direc101)
!!$    CALL  KxcWalkPIGW(Node%Direc110)
!!$    CALL  KxcWalkPIGW(Node%Direc111)
!
  END SUBROUTINE KxcWalkDIPMW
!=========================================================================================
! Integral Code: 
!=========================================================================================
! Integrate Int x^n Exp[-x*x],{x,a,b}]
!=========================================================================================
  SUBROUTINE XVector(IPow,A,B,XV)
    INTEGER                     :: IPOW
    REAL(DOUBLE)                :: A,B
    REAL(DOUBLE),PARAMETER      :: SqrtPi2 = 0.88622692545275801365
    REAL(DOUBLE),DIMENSION(0:6) :: XV
!
    SELECT CASE (IPow)
    CASE(0)
       XV(0) = SqrtPi2*(ERF(B)-ERF(A))
    CASE(1)
       Tmp(1) = ERF(A)
       Tmp(2) = ERF(B)
       Tmp(3) = EXP(-A*A)
       Tmp(4) = EXP(-B*B)
!       
       XV(0)  = SqrtPi2*(Tmp(2)-Tmp(1))
       XV(1)  = Half*(Tmp(3)-Tmp(4))
    CASE(2)
       Tmp(1) = ERF(A)
       Tmp(2) = ERF(B)
       Tmp(3) = EXP(-A*A)
       Tmp(4) = EXP(-B*B)
!       
       XV(0)  = SqrtPi2*(Tmp(2)-Tmp(1))
       XV(1)  = Half*(Tmp(3)-Tmp(4))
       XV(2)  = Half*(XV(0)+A*Tmp(3)-B*Tmp(4))
    CASE(3)
       Tmp(1) = A*A
       Tmp(2) = B*B
       Tmp(3) = ERF(A)
       Tmp(4) = ERF(B)
       Tmp(5) = EXP(-Tmp(1))
       Tmp(6) = EXP(-Tmp(2))
!       
       XV(0)  = SqrtPi2*(Tmp(4)-Tmp(3))
       XV(1)  = Half*(Tmp(5)-Tmp(6))
       XV(2)  = Half*(XV(0)+A*Tmp(5)-B*Tmp(6))
       XV(3)  = Half*(Tmp(1)*Tmp(5)-Tmp(2)*Tmp(6))+XV(1)
    CASE(4)
       Tmp(1) = A*A
       Tmp(2) = B*B
       Tmp(3) = ERF(A)
       Tmp(4) = ERF(B)
       Tmp(5) = EXP(-Tmp(1))
       Tmp(6) = EXP(-Tmp(2))
!       
       XV(0)  = SqrtPi2*(Tmp(4)-Tmp(3))
       XV(1)  = Half*(Tmp(5)-Tmp(6))
       XV(2)  = Half*(XV(0)+A*Tmp(5)-B*Tmp(6))
       XV(3)  = Half*(Tmp(1)*Tmp(5)-Tmp(2)*Tmp(6))+XV(1)
    END SELECT

  END SUBROUTINE XVector
!
!
!
END MODULE
