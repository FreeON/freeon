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
! Tmp Matrices
  REAL(DOUBLE),DIMENSION(0:DOrder,0:4,0:100)  :: TmpMatX,TmpMatY,TmpMatZ
  CONTAINS 
!-------------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE KxcWalkDIPMW(Node)
    TYPE(PIGWNode), POINTER                      :: Node
    INTEGER                                      :: Ell,I,J,K,IJK,ILink
    INTEGER                                      :: STypeX,STypeY,STypeZ
    REAL(DOUBLE)                                 :: PQx,PQy,PQz
    REAL(DOUBLE)                                 :: Zeta,SqZeta
    REAL(DOUBLE),DIMENSION(0:DOrder,0:4)         :: WxVec,WyVec,WzVec
!
!   Does PBox intersect this DIPMW BBox? 
!
    IF( ABS(PBox%Center(1)-Node%Box%Center(1))>PBox%Half(1)+Node%Box%Half(1) ) RETURN
    IF( ABS(PBox%Center(2)-Node%Box%Center(2))>PBox%Half(2)+Node%Box%Half(2) ) RETURN
    IF( ABS(PBox%Center(3)-Node%Box%Center(3))>PBox%Half(3)+Node%Box%Half(3) ) RETURN
!   Distances to Wavelet Centers
    PQx = PBox%Center(1)-Node%X(1)
    PQy = PBox%Center(2)-Node%X(2)
    PQz = PBox%Center(3)-Node%X(3)
    Ell     = Prim%Ell
    Zeta    = Prim%Zeta 
    SqZeta  = SQRT(Zeta)
!   Determine Surface Type, Should be moved to the Tree Build {Node%SType(1:3) }
    STypeX=0
    IF(Node%X(1)-Node%Box%BndBox(1,1) < 1.D-12) STypeX=-1
    IF(Node%Box%BndBox(1,2)-Node%X(1) < 1.D-12) STypeX= 1
    STypeY=0
    IF(Node%X(2)-Node%Box%BndBox(2,1) < 1.D-12) STypeY=-1
    IF(Node%Box%BndBox(2,2)-Node%X(2) < 1.D-12) STypeY= 1
    STypeZ=0
    IF(Node%X(3)-Node%Box%BndBox(3,1) < 1.D-12) STypeZ=-1
    IF(Node%Box%BndBox(3,2)-Node%X(3) < 1.D-12) STypeZ= 1
!   Compute Wavelet Integrals
    CALL WVector(Ell,STypeX,Zeta,SqZeta,PQx,Node%DX(1),WxVec)
    CALL WVector(Ell,STypeY,Zeta,SqZeta,PQy,Node%DX(2),WyVec)
    CALL WVector(Ell,STypeZ,Zeta,SqZeta,PQz,Node%DX(3),WzVec)
!   Compute the Ket
    SELECT CASE(Ell)  
    CASE(0)
       DO K=0,Dorder
          DO J=0,Dorder
             DO I=0,Dorder  
                IJK = 1+I+DMax*J+DMax*DMax*K
                Ket(1) = Ket(1)+Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,0)*WzVec(K,0)
             ENDDO
          ENDDO
       ENDDO
    CASE(1)
       DO K=0,Dorder
          DO J=0,Dorder
             DO I=0,Dorder  
                IJK = 1+I+DMax*J+DMax*DMax*K
                Ket(1) = Ket(1)+Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,0)*WzVec(K,0)
                Ket(2) = Ket(2)-Node%WCoef(IJK)*WxVec(I,1)*WyVec(J,0)*WzVec(K,0)
                Ket(3) = Ket(3)-Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,1)*WzVec(K,0)
                Ket(4) = Ket(4)-Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,0)*WzVec(K,1)
             ENDDO
          ENDDO
       ENDDO
    CASE(2)
       DO K=0,Dorder
          DO J=0,Dorder
             DO I=0,Dorder  
                IJK = 1+I+DMax*J+DMax*DMax*K
                Ket(1) = Ket(1) +Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,0)*WzVec(K,0)
                Ket(2) = Ket(2) -Node%WCoef(IJK)*WxVec(I,1)*WyVec(J,0)*WzVec(K,0)
                Ket(3) = Ket(3) -Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,1)*WzVec(K,0)
                Ket(4) = Ket(4) -Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,0)*WzVec(K,1)
!
                Ket(5) = Ket(5) +Node%WCoef(IJK)*WxVec(I,2)*WyVec(J,0)*WzVec(K,0)
                Ket(6) = Ket(6) +Node%WCoef(IJK)*WxVec(I,1)*WyVec(J,1)*WzVec(K,0)
                Ket(7) = Ket(7) +Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,2)*WzVec(K,0)
!
                Ket(8) = Ket(8) +Node%WCoef(IJK)*WxVec(I,1)*WyVec(J,0)*WzVec(K,1)
                Ket(9) = Ket(9) +Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,1)*WzVec(K,1)
                Ket(10)= Ket(10)+Node%WCoef(IJK)*WxVec(I,0)*WyVec(J,0)*WzVec(K,2)
             ENDDO
          ENDDO
       ENDDO
    END SELECT
    IF(Node%LeafType=='EnddLeaf') RETURN
!   Continue Walking
    DO ILink=1,Node%NLinks
       CALL KxcWalkDIPMW(Node%Links(ILink)%Link)
    ENDDO
  END SUBROUTINE KxcWalkDIPMW
!=========================================================================================
! Integral Code: 
!=========================================================================================
! Integrate_{-Half,Half} w[n,x/Half] d^q/dx^q Exp[-bet*(x-Ax)^2]
!=========================================================================================
  SUBROUTINE WVector(Ell,SType,Zeta,SqZeta,RDist,DX0,WV)
    INTEGER                                      :: Ell,SType
    REAL(DOUBLE)                                 :: Zeta,SqZeta,RDist,DX0
    REAL(DOUBLE)                                 :: OneSqZeta,RZeta1,XZeta1,RmXZeta1,RpXZeta1
    REAL(DOUBLE)                                 :: ErfRmX0,Erf0,ErfRpX0,ExpRmX0,Exp0,ExpRpX0
    REAL(DOUBLE),DIMENSION(0:DOrder,0:4)         :: WV    
    REAL(DOUBLE)                                 :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14, &
                                                    t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26, &
                                                    t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38, &
                                                    t39,t40,t41,t42,t43,t44,t45,t46,t47,t48,t49,t50, &
                                                    t51,t52,t53,t54,t55,t56,t57,t58,t59,t60,t61
    REAL(DOUBLE)                                 :: IntL0,IntL1,IntL2,IntL3,IntL4
    REAL(DOUBLE)                                 :: IntD0L0,IntD0L1,IntD0L2,IntD0L3,IntD0L4,IntD0L5,IntD0L6,IntD0L7, &
                                                    IntD1L0,IntD1L1,IntD1L2,IntD1L3,IntD1L4,IntD1L5,IntD1L6,IntD1L7, &
                                                    IntD2L0,IntD2L1,IntD2L2,IntD2L3,IntD2L4,IntD2L5,IntD2L6,IntD2L7, &
                                                    IntD3L0,IntD3L1,IntD3L2,IntD3L3,IntD3L4,IntD3L5,IntD3L6,IntD3L7
    REAL(DOUBLE)                                 :: RZeta2,XZeta2,RZeta3,XZeta3,RZeta4,XZeta4,       &
                                                    RZeta5,XZeta5,RZeta6,XZeta6,RZeta7,XZeta7,       &
                                                    RZeta8,XZeta8,RZeta9,XZeta9,RZeta10,XZeta10,     &
                                                    RZeta11,XZeta11,RZeta12,XZeta12,RZeta13,XZeta13, &
                                                    RZeta14,XZeta14,RZeta15,XZeta15,RZeta16,XZeta16, &
                                                    RZeta17,XZeta17,RZeta18,XZeta18,RZeta19,XZeta19, &
                                                    RZeta20,XZeta20
    
!
    RZeta1   = SqZeta*RDist
    XZeta1   = SqZeta*DX0
    IF(XZeta1 < 0.2D0) THEN
!      Asymtotic form which only uses one Exp function
       Exp0    = EXP(-RZeta1*RZeta1)
!
       SELECT CASE(DOrder)
       CASE(0)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_Asy_0_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_Asy_0_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_Asy_0_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_Asy_0_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_Asy_0_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_Asy_0_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_Asy_0_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_Asy_0_7.Inc'
          CASE(8:)
!             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(1)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_Asy_1_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_Asy_1_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_Asy_1_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_Asy_1_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_Asy_1_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_Asy_1_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_Asy_1_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_Asy_1_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(2)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_Asy_2_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_Asy_2_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_Asy_2_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_Asy_2_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_Asy_2_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_Asy_2_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_Asy_2_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_Asy_2_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(3)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_Asy_3_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_Asy_3_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_Asy_3_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_Asy_3_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_Asy_3_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_Asy_3_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_Asy_3_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_Asy_3_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(4:)
          CALL Halt("Integrals for this DOrder not Implimented")
       END SELECT
    ELSE
!      Standard Form which uses the Erf functions
       RmXZeta1 = RZeta1-XZeta1
       RpXZeta1 = RZeta1+XZeta1
!
       ErfRmX0  = ERF(RmXZeta1)
       Erf0     = ERF(RZeta1)
       ErfRpX0  = ERF(RpXZeta1)
!
       ExpRmX0  = EXP(-RmXZeta1*RmXZeta1)
       Exp0     = EXP(-RZeta1*RZeta1)
       ExpRpX0  = EXP(-RpXZeta1*RpXZeta1)
!
       SELECT CASE(DOrder)
       CASE(0)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_0_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_0_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_0_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_0_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_0_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_0_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_0_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_0_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(1)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_1_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_1_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_1_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_1_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_1_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_1_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_1_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_1_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(2)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_2_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_2_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_2_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_2_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_2_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_2_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_2_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_2_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(3)
          SELECT CASE(Ell)
          CASE(0)
             INCLUDE 'MMA/IntWavelet_3_0.Inc'
          CASE(1)
             INCLUDE 'MMA/IntWavelet_3_1.Inc'
          CASE(2)
             INCLUDE 'MMA/IntWavelet_3_2.Inc'
          CASE(3)
             INCLUDE 'MMA/IntWavelet_3_3.Inc'
          CASE(4)
!             INCLUDE 'MMA/IntWavelet_3_4.Inc'
          CASE(5)
!             INCLUDE 'MMA/IntWavelet_3_5.Inc'
          CASE(6)
!             INCLUDE 'MMA/IntWavelet_3_6.Inc'
          CASE(7)
!             INCLUDE 'MMA/IntWavelet_3_7.Inc'
          CASE(8:)
             CALL Halt("Integrals for this Ell not Implimented")
          END SELECT
       CASE(4:)
          CALL Halt("Integrals for this DOrder not Implimented")
       END SELECT
    ENDIF
!
  END SUBROUTINE WVector
!=========================================================================================
! Integrate_{-Half,Half} w[n,x/Half] d^q/dx^q Exp[-bet*(x-Ax)^2] Asy
!=========================================================================================
!!$  FUNCTION WVector_Asy(Ell,R,DX) RESULT(WV)
!!$!
!!$    R2      = Two*R
!!$    Herm(0) = One
!!$    Herm(1) = R2
!!$    DO I=2,MaxL+Ell
!!$       Herm(I) = R2*Herm(I-1)-(I-2)*Herm(I-2)
!!$    ENDDO
!!$    WV = Zero
!!$    DelP = DX
!!$    DO I=0,MaxL
!!$       DO L=0,Ell
!!$          Tmp = Herm(I+L)*DelP
!!$          DO D=0,DOrder
!!$             WV(D,L) = WV(D,L)+Tmp*Coef(D,L,I)
!!$          ENDDO
!!$       ENDDO
!!$       DelP = DelP*DX
!!$    ENDDO
!!$  END FUNCTION WVector_Asy
!
!
!
END MODULE

