MODULE InCoords
!
   USE ControlStructures
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE GlobalCharacters
   USE InOut
   USE MemMan
   USE SetXYZ
   USE ProcessControl
   USE PrettyPrint
   USE ParsingConstants
   USE LinAlg
   USE AInv   
   USE CholFactor
   USE AtomPairs
   !
   IMPLICIT NONE
!
CONTAINS
!
!----------------------------------------------------------------
!
   SUBROUTINE BMatrix(XYZ,IntCs,B,PBCDim,LinCrit,TorsLinCrit)
     !
     ! This subroutine calculates the sparse, TYPE(BMATR) type 
     ! representation of the B matrix 
     ! In current version Cartesian coordinates 
     ! must be passed in in AU
     ! Linear bendings _must_ always appear in pairs, 
     ! Defd as LINB1 and LINB2
     !
     IMPLICIT NONE
     INTEGER :: I,J,K,L,NatmsLoc,IInt,KA,LA
     INTEGER :: II,IC1,IC2,M,N,NCart,PBCDim,III,ALPHA,BETA,JI,JJ
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value
     TYPE(INTC)                  :: IntCs
     TYPE(BMATR)                 :: B
     REAL(DOUBLE)                :: LinCrit,TorsLinCrit,Frac(12)
     REAL(DOUBLE)                :: Vect(3),Vect1(3)
     REAL(DOUBLE),DIMENSION(3,3) :: BoxShapeT,BoxShape,InvBoxSh
     REAL(DOUBLE),DIMENSION(3,4) :: XYZAux
     REAL(DOUBLE),DIMENSION(12)  :: Aux12
     CHARACTER(LEN=DCL)          :: Char
     TYPE(DBL_VECT)              :: CartAux
     TYPE(INT_VECT)              :: AtmsAux
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     CALL INTCValue(IntCs,XYZ,PBCDim,LinCrit,TorsLinCrit)
     !
     DO J=1,3
       BoxShapeT(J,1:3)=XYZ(1:3,NatmsLoc-3+J)
       BoxShape(1:3,J)=XYZ(1:3,NatmsLoc-3+J)
     ENDDO
     InvBoxSh=InverseMatrix(BoxShape)
     !
     ! allocate B matrix
     !
     CALL NEW(B,IntCs%N)
     B%IB%I=0
     B%B%D=Zero
     !
     DO IInt=1,IntCs%N   
       IF(.NOT.IntCs%Active%L(IInt)) CYCLE
       CALL PBCXYZAux(XYZ,BoxShapeT,XYZAux,IntCs,IInt)
       Char=IntCs%Def%C(IInt)
       IF(PBCDim>0) THEN
         IF(Char(1:5)=='ALPHA') THEN
           Char(1:5)='BEND '
         ELSE IF(Char(1:4)=='BETA') THEN
           Char(1:5)='BEND '
         ELSE IF(Char(1:5)=='GAMMA') THEN
           Char(1:5)='BEND '
         ENDIF
       ENDIF
       !
       IF(Char(1:4)=='STRE') THEN
         ! stre
         B%IB%I(IInt,1:2)=IntCs%Atoms%I(IInt,1:2)
         Aux12=B%B%D(IInt,1:12)
         CALL STRE(XYZAux(1:3,1),XYZAux(1:3,2),BB_O=Aux12)
         B%B%D(IInt,1:12)=Aux12
        !write(*,100) Char(1:5), &
        !B%IB%I(IInt,1:4),B%B%D(IInt,1:12)
         !
       ELSE IF(Char(1:4)=='BEND'.OR. &
               Char(1:5)=='ALPHA'.OR. &
               Char(1:4)=='BETA'.OR. &
               Char(1:5)=='GAMMA') THEN
         ! bend
         B%IB%I(IInt,1:3)=IntCs%Atoms%I(IInt,1:3)
         Aux12=B%B%D(IInt,1:12)
         CALL BEND(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),BB_O=Aux12)
         B%B%D(IInt,1:12)=Aux12
        !write(*,100) Char(1:5), &
        !B%IB%I(IInt,1:4),B%B%D(IInt,1:12)
       ELSE IF(Char(1:4)=='OUTP') THEN
         ! out of plane
         B%IB%I(IInt,1:4)=IntCs%Atoms%I(IInt,1:4)
         Aux12=B%B%D(IInt,1:12)
         CALL OutP(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3), &
                   XYZAux(1:3,4),TorsLinCrit,IntCs%Active%L(IInt), &
                   BB_O=Aux12)
         B%B%D(IInt,1:12)=Aux12
       ! write(*,100) Char(1:5),B%IB%I(IInt,1:4), &
       !  B%B%D(IInt,1:12)*AngstromsToAu
       ELSE IF(Char(1:4)=='TORS') THEN
         ! torsion of i-j-k-l
         B%IB%I(IInt,1:4)=IntCs%Atoms%I(IInt,1:4)
         Aux12=B%B%D(IInt,1:12)
         CALL TORS(XYZAux(1:3,1),XYZAux(1:3,2), &
                   XYZAux(1:3,3),XYZAux(1:3,4), &
                   TorsLinCrit,IntCs%Active%L(IInt), &
                   BB_O=Aux12)
         B%B%D(IInt,1:12)=Aux12
        !write(*,*) IInt
        !write(*,100) Char(1:5),B%IB%I(IInt,1:4), &
        !B%B%D(IInt,1:12)
       ELSE IF(Char(1:6)=='VOLM_L') THEN
         B%IB%I(IInt,1:4)=IntCs%Atoms%I(IInt,1:4)
         Aux12=B%B%D(IInt,1:12)
         CALL VOLUME(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),&
           XYZAux(1:3,4),IntCs%Active%L(IInt),BB_O=Aux12)
         B%B%D(IInt,1:12)=Aux12
       ELSE IF(Char(1:6)=='AREA_L') THEN
         B%IB%I(IInt,1:3)=IntCs%Atoms%I(IInt,1:3)
         Aux12=B%B%D(IInt,1:12)
         CALL AREA(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),&
                   IntCs%Active%L(IInt),BB_O=Aux12)
         B%B%D(IInt,1:12)=Aux12
       ELSE IF(Char(1:5)=='LINB1') THEN
         ! linear bendig of i-j-k
         IF(IntCs%Def%C(IInt+1)(1:5)/='LINB2') &
            CALL Halt('LINB2 Definitions are not paired!')
         B%IB%I(IInt,1:3)=IntCs%Atoms%I(IInt,1:3)
         B%IB%I(IInt+1,1:3)=IntCs%Atoms%I(IInt+1,1:3)
         L=IntCs%Atoms%I(IInt,4)  ! reference atom
         CALL LinB(XYZAux(1:3,1),XYZAux(1:3,2), &
                   XYZAux(1:3,3),XYZAux(1:3,4),L, &
                   BB1=B%B%D(IInt,1:12),BB2=B%B%D(IInt+1,1:12))  
        !write(*,100) Char(1:5), &
        ! B%IB%I(IInt,1:4),B%B%D(IInt,1:12)
        !write(*,100) Char(1:5), &
        ! B%IB%I(IInt+1,1:4),B%B%D(IInt+1,1:12)
       ELSE IF(Char(1:5)=='LINB2') THEN
         CYCLE
       ELSE IF(Char(1:5)=='CARTX' ) THEN
         I=IntCs%Atoms%I(IInt,1)
         B%IB%I(IInt,1)=I
         CALL BCART('X',B%B%D(IInt,1:12),IntCs%Constraint%L(IInt))
       ELSE IF(Char(1:5)=='CARTY' ) THEN
         I=IntCs%Atoms%I(IInt,1)
         B%IB%I(IInt,1)=I
         CALL BCART('Y',B%B%D(IInt,1:12),IntCs%Constraint%L(IInt))
       ELSE IF(Char(1:5)=='CARTZ' ) THEN
         I=IntCs%Atoms%I(IInt,1)
         B%IB%I(IInt,1)=I
         CALL BCART('Z',B%B%D(IInt,1:12),IntCs%Constraint%L(IInt))
       ENDIF
     ENDDO !!!! loop over internal coords
     100 FORMAT(A5,4I6,/,6F12.6,/,6F12.6)
     !
     ! Generate BL, the portion of the B matrix related to lattice
     ! distorsions
     !
     B%BLI%I=0
     B%BL%D=Zero
     IF(PBCDim>0) THEN
       DO IInt=1,IntCs%N
         IF(.NOT.IntCs%Active%L(IInt)) CYCLE
         ! filter out fixed fractional coordinates
         IF(IntCs%Def%C(IInt)(1:4)=='CART') CYCLE
         !
         CALL PBCXYZAux(XYZ,BoxShapeT,XYZAux,IntCs,IInt)
         DO J=1,4
           K=3*(J-1)+1
           L=K+2
           Vect=XYZAux(1:3,J)
           CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect,Frac(K:L))
         ENDDO
         B%BLI%I(IInt)=NCart-9
         II=0
         DO BETA=1,PBCDim
           DO ALPHA=1,3
             II=II+1
             DO J=1,4
               IF(IntCs%Atoms%I(IInt,J)==0) EXIT
               JI=3*(J-1)+ALPHA
               JJ=3*(J-1)+BETA
               B%BL%D(IInt,II)=B%BL%D(IInt,II)+B%B%D(IInt,JI)*Frac(JJ)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     !
     ! Make lattice sum for atomic coordinates
     !
     IF(PBCDim>0) THEN
       CALL New(CartAux,NCart)
       CALL New(AtmsAux,NatmsLoc)
       DO IInt=1,IntCs%N
         IF(.NOT.IntCs%Active%L(IInt)) CYCLE
         !
         ! clean memory
         !
         DO J=1,4
           I=IntCs%Atoms%I(IInt,J)
           IF(I==0) EXIT
            AtmsAux%I(I)=0
           KA=3*(I-1)+1   
           LA=KA+2
            CartAux%D(KA:LA)=Zero
         ENDDO
         !
         ! do lattice sum
         !
         DO J=1,4
           I=IntCs%Atoms%I(IInt,J)
           IF(I==0) EXIT
           K=3*(J-1)+1
           L=K+2
           KA=3*(I-1)+1   
           LA=KA+2
            CartAux%D(KA:LA)=CartAux%D(KA:LA)+B%B%D(IInt,K:L)
            IF(AtmsAux%I(I)==0) THEN
              AtmsAux%I(I)=I
            ELSE
              B%IB%I(IInt,J)=0
            ENDIF
         ENDDO
         !
         ! fill back
         !
         II=0
         B%B%D(IInt,1:12)=Zero
         DO J=1,4
           I=IntCs%Atoms%I(IInt,J)
           IF(I==0) EXIT
           IF(B%IB%I(IInt,J)==0) THEN
             CYCLE
           ELSE
             B%IB%I(IInt,J)=0
           ENDIF
           II=II+1
           B%IB%I(IInt,II)=I
           K=3*(II-1)+1
           L=K+2
           KA=3*(I-1)+1   
           LA=KA+2
           B%B%D(IInt,K:L)=CartAux%D(KA:LA)
         ENDDO
       ENDDO
       CALL Delete(AtmsAux)
       CALL Delete(CartAux)
     ENDIF
     !
     ! Transform molecular B-matrix into fractional coords
     !
     IF(PBCDim>0) THEN
       DO IInt=1,IntCs%N
         IF(.NOT.IntCs%Active%L(IInt)) CYCLE
         DO J=1,4
           IF(IntCs%Atoms%I(IInt,J)==0) EXIT
           K=3*(J-1)+1
           L=K+2
           Vect=B%B%D(IInt,K:L)
           CALL DGEMM_NNc(3,3,1,One,Zero,BoxShapeT,Vect,Vect1)
           B%B%D(IInt,K:L)=Vect1
         ENDDO
       ENDDO
     ENDIF
     !
     ! Clean BL for fixed lattice orientation
     !
     IF(PBCDim>0) THEN
       DO IInt=1,IntCs%N
         B%BL%D(IInt,2)=Zero
         B%BL%D(IInt,3)=Zero
         B%BL%D(IInt,6)=Zero
       ENDDO
     ENDIF
     !
   END SUBROUTINE BMatrix
!-----------------------------------------------------------------------
!
   SUBROUTINE Stre(XI,XJ,BB_O,Value_O)
     !
     REAL(DOUBLE),DIMENSION(1:12),OPTIONAL :: BB_O      
     REAL(DOUBLE),DIMENSION(1:3)           :: XI,XJ,U   
     REAL(DOUBLE),OPTIONAL                 :: Value_O      
     REAL(DOUBLE)                          :: RU           
     !
     U=XI-XJ
     RU=SQRT(DOT_PRODUCT(U,U))
     U=U/RU
     !
     IF(PRESENT(Value_O)) Value_O=RU
     IF(PRESENT(BB_O)) THEN
       BB_O(1:3)=U
       BB_O(4:6)=-U
     ENDIF
   END SUBROUTINE STRE
!
!----------------------------------------------------------------
!
   SUBROUTINE Bend(XI,XJ,XK,BB_O,Value_O,W_O)
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,W,U,V,Ref,Aux
     REAL(DOUBLE),DIMENSION(1:3),OPTIONAL :: W_O
     REAL(DOUBLE),DIMENSION(1:12),OPTIONAL:: BB_O
     REAL(DOUBLE)                :: SumU,Sum2,RU,RV,RW
     REAL(DOUBLE),OPTIONAL       :: Value_O
     INTEGER                     :: I,J
     !
     U=XI-XJ
     V=XK-XJ
     IF(PRESENT(W_O)) THEN
       RW=SQRT(DOT_PRODUCT(W_O,W_O))
       W_O=W_O/RW
       U=U-W_O*DOT_PRODUCT(W_O,U)
       V=V-W_O*DOT_PRODUCT(W_O,V)
     ENDIF
     RU=SQRT(DOT_PRODUCT(U,U))
     RV=SQRT(DOT_PRODUCT(V,V))
     U=U/RU
     V=V/RV
     !
     IF(PRESENT(Value_O)) THEN
       SumU=DOT_PRODUCT(U,V)
       IF(ABS(SumU)>One-1.D-8) THEN
         Value_O=PI
       ELSE
         Value_O=ACOS(SumU)
       ENDIF
       IF(PRESENT(W_O)) THEN
         CALL CROSS_PRODUCT(U,V,Aux)
         Sum2=DOT_PRODUCT(Aux,W_O)
         IF(ABS(Sum2)>1.D-8) THEN
           Value_O=SIGN(Value_O,Sum2)
         ENDIF
       ENDIF
     ENDIF
     !
     IF(PRESENT(BB_O)) THEN
       IF(PRESENT(W_O)) THEN
         W=W_O
       ELSE
         CALL CROSS_PRODUCT(U,V,W)
         RW=SQRT(DOT_PRODUCT(W,W))
         IF(RW<0.01D0) THEN
           Ref=(/One,-One,One/)
           CALL CROSS_PRODUCT(U,Ref,W)
           RW=SQRT(DOT_PRODUCT(W,W))
           IF(RW<0.01D0) THEN
             Ref=(/-One,One,One/)
             CALL CROSS_PRODUCT(U,Ref,W)
             RW=SQRT(DOT_PRODUCT(W,W))
           ENDIF
         ENDIF
         W=W/RW
       ENDIF
       CALL LinBGen(BB_O,W,U,V,RU,RV)
     ENDIF
   END SUBROUTINE Bend
!
!----------------------------------------------------------------
!
   SUBROUTINE OutP(XI,XJ,XK,XL,LinCrit,Active,BB_O,Value_O)
     !
     !  i is the end atom (atom wagged with respect to j-k-l plane).
     !  j is the apex atom (Atoms i, k and l are attached to j).
     !  k and l are the anchor Atoms (Define the j-k-l plane).
     !
     IMPLICIT NONE
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,XL,RJI,RJK,RJL
     REAL(DOUBLE),DIMENSION(1:3) :: C1,SMI,SMK,SML
     REAL(DOUBLE)          :: DJK,DJI,DJL,CosI,CosK,CosL
     REAL(DOUBLE)          :: SinSin,SinI,Dot,SinT,CosT,TanT,CosSin
     REAL(DOUBLE)          :: SinJ,Conv,LinCrit
     INTEGER               :: I,J
     REAL(DOUBLE),OPTIONAL :: BB_O(1:12),Value_O
     LOGICAL               :: Active
     !
     Active=.TRUE.
     Conv=180.D0/PI
     RJI=XI-XJ
     RJK=XK-XJ
     RJL=XL-XJ
     DJI=SQRT(DOT_PRODUCT(RJI,RJI))
     DJK=SQRT(DOT_PRODUCT(RJK,RJK))
     DJL=SQRT(DOT_PRODUCT(RJL,RJL))
     IF(DJI<1.D-4.OR.DJK<1.D-4.OR.DJL<1.D-4) THEN
       Active=.FALSE.
       RETURN 
     ENDIF
     RJI=RJI/DJI
     RJK=RJK/DJK
     RJL=RJL/DJL
     CosI=DOT_PRODUCT(RJK,RJL)
     CosK=DOT_PRODUCT(RJI,RJL)
     CosL=DOT_PRODUCT(RJI,RJK)
     IF(ABS(180.D0-Conv*ACOS(CosI)) < LinCrit) THEN
       Active=.FALSE.
       RETURN  
     ENDIF
     SinSin=One-CosI*CosI
     SinI=SQRT(SinSin)
     CALL CROSS_PRODUCT(RJK,RJL,C1)
     IF(SinI<1.D-4) THEN
       Active=.FALSE.
       RETURN
     ENDIF
     Dot=DOT_PRODUCT(RJI,C1) ! this may be +-
     SinT=Dot/SinI
     IF(ABS(SinT)>0.99995D0) THEN
       Active=.FALSE.
       RETURN !!!! sint shows OutP angle
     ENDIF
     CosT=SQRT(One-SinT*SinT)
     TanT=SinT/CosT
     CosSin=CosT*SinI
     IF(PRESENT(Value_O)) THEN
       Value_O=SIGN(ASIN(SinT),-Dot)
     ENDIF
     IF(PRESENT(BB_O)) THEN
       C1=C1/CosSin
       SMI=(C1-TanT*RJI)/DJI
       SMK=C1*(CosI*CosK-CosL)/(SinSin*Djk)
       SML=C1*(CosI*CosL-CosK)/(SinSin*Djl)
       BB_O(1:3)=SMI
       BB_O(7:9)=SMK
       BB_O(10:12)=SML
       BB_O(4:6)=-SMI-SMK-SML
     ENDIF
   END SUBROUTINE OutP
!
!----------------------------------------------------------------------
!
   SUBROUTINE Tors(XI,XJ,XK,XL,TorsCrit,Active,BB_O,Value_O)
     REAL(DOUBLE),DIMENSION(1:12),OPTIONAL :: BB_O
     REAL(DOUBLE),DIMENSION(1:3)  :: XI,XJ,XK,XL,U,W,V,UxW,VxW,UxV
     REAL(DOUBLE),DIMENSION(1:3)  :: Term1,Term2,Term3
     REAL(DOUBLE),OPTIONAL        :: Value_O
     REAL(DOUBLE)                 :: TorsCrit 
     REAL(DOUBLE)                 :: CosPhiU,CosPhiV
     REAL(DOUBLE)                 :: SinPhiU,SinPhiV
     REAL(DOUBLE)                 :: SinPhiU2,SinPhiV2
     LOGICAL                      :: Active
     REAL(DOUBLE)                 :: RU,RV,RW,SumU,Conv
     INTEGER                      :: I,J,K,L
     !
     Conv=180.D0/PI
     U=XI-XJ
     V=XL-XK
     W=XK-XJ 
     RU=SQRT(DOT_PRODUCT(U,U))
     RV=SQRT(DOT_PRODUCT(V,V))
     RW=SQRT(DOT_PRODUCT(W,W))
     U=U/RU
     V=V/RV
     W=W/RW
     CosPhiU=DOT_PRODUCT(U,W)
     CosPhiV=DOT_PRODUCT(V,-W)
     SinPhiU=SQRT(ABS(One-CosPhiU*CosPhiU))
     SinPhiV=SQRT(ABS(One-CosPhiV*CosPhiV))
     SinPhiU2=SinPhiU*SinPhiU
     SinPhiV2=SinPhiV*SinPhiV
     IF(Conv*ASIN(SinPhiU)<TorsCrit.OR.Conv*ASIN(SinPhiV)<TorsCrit) THEN
       Active=.FALSE.
       IF(PRESENT(Value_O)) Value_O=Zero
       IF(PRESENT(BB_O)) BB_O=Zero
       RETURN
     ENDIF
     !
     CALL CROSS_PRODUCT(U,W,UxW)
     CALL CROSS_PRODUCT(V,W,VxW)
     CALL CROSS_PRODUCT(U,V,UxV)
     !
     IF(PRESENT(Value_O)) THEN
       SumU=DOT_PRODUCT(UxW,VxW)/(SinPhiU*SinPhiV)
       IF(ABS(SumU)>One) THEN
         SumU=SIGN(One,SumU)
       ENDIF
         Value_O=ACOS(SumU)
       SumU=DOT_PRODUCT(W,UxV) ! orientation of torsion
       IF(ABS(SumU)>Zero) THEN
         Value_O=SIGN(Value_O,SumU)
       ENDIF
     ENDIF
     !
     IF(PRESENT(BB_O)) THEN
       Term1=UxW/(RU*SinPhiU2) 
       Term2=VxW/(RV*SinPhiV2) 
       Term3=UxW*CosPhiU/(RW*SinPhiU2)+VxW*CosPhiV/(RW*SinPhiV2)
       BB_O(1:3)=Term1
       BB_O(4:6)=-Term1+Term3
       BB_O(7:9)=Term2-Term3
       BB_O(10:12)=-Term2
     ENDIF
     !
   END SUBROUTINE Tors
!
!----------------------------------------------------------------------
!
   SUBROUTINE LinB(XI,XJ,XK,X4,IRef,BB1,BB2,Value1,Value2)
     !
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK,X4
     REAL(DOUBLE),DIMENSION(1:12),OPTIONAL:: BB1,BB2
     INTEGER                     :: I,J,K,L,M,N,KK,IRef
     REAL(DOUBLE)                :: SumU,CosAlpha,SinAlpha,RJI,RJK,RKI
     REAL(DOUBLE),OPTIONAL       :: Value1,Value2
     REAL(DOUBLE),DIMENSION(1:3) :: VectJI,VectJK,VectKI,Aux1,Aux2
     REAL(DOUBLE),DIMENSION(1:3) :: VectAux,EZ,EY,EX,RIJJK,VectRef
     REAL(DOUBLE),DIMENSION(3,3) :: Rot   
     LOGICAL                     :: Active
     !
     VectJI=XI-XJ
     VectJK=XK-XJ
     VectKI=XK-XI
     RJI=SQRT(DOT_PRODUCT(VectJI,VectJI))
     RJK=SQRT(DOT_PRODUCT(VectJK,VectJK))
     RKI=SQRT(DOT_PRODUCT(VectKI,VectKI))
     VectJI=VectJI/RJI
     VectJK=VectJK/RJK
     VectKI=VectKI/RKI
     EZ=VectKI
     !
     CALL CROSS_PRODUCT(VectJI,VectJK,RIJJK)
     CosAlpha=SQRT(DOT_PRODUCT(RIJJK,RIJJK))
     !
     ! Set up reference coordinate system
     !
    !IF(ABS(ACOS(CosAlpha)*180.D0/PI-180.D0)>1.D0) THEN
    !IF(SelfRerefence) THEN
    !  VectRef=RIJJK
    !  CALL CROSS_PRODUCT(VectRef,EZ,EY)
    !  SumU=SQRT(DOT_PRODUCT(EY,EY))
    !ELSE
         IF(IRef==0) THEN
           VectRef=(/One,Zero,Zero/)
         ELSE
           VectRef=X4-XK
           VectRef=VectRef/SQRT(DOT_PRODUCT(VectRef,VectRef))
         ENDIF
         CALL CROSS_PRODUCT(VectRef,EZ,EY)
         SumU=SQRT(DOT_PRODUCT(EY,EY))
         IF(SumU<0.01D0) THEN
           IF(IRef==0) THEN
             VectRef=(/Zero,One,Zero/)
           ELSE
             VectRef=X4-XI
             VectRef=VectRef/SQRT(DOT_PRODUCT(VectRef,VectRef))
           ENDIF
           CALL CROSS_PRODUCT(VectRef,EZ,EY)
           SumU=SQRT(DOT_PRODUCT(EY,EY))
           IF(SumU<0.01D0) THEN
             IF(IRef==0) THEN
               VectRef=(/Zero,Zero,One/)
             ELSE
               VectRef=X4-XJ
               VectRef=VectRef/SQRT(DOT_PRODUCT(VectRef,VectRef))
             ENDIF
             CALL CROSS_PRODUCT(VectRef,EZ,EY)
             SumU=SQRT(DOT_PRODUCT(EY,EY))
           ENDIF           
         ENDIF           
    !ENDIF           
     EY=EY/SumU
     CALL CROSS_PRODUCT(EY,EZ,EX)
     SumU=SQRT(DOT_PRODUCT(EX,EX))
     EX=EX/SumU
     !
     IF(PRESENT(Value1)) THEN
      !CALL BEND(XI,XJ,XK,Value_O=Value1,W_O=EX)
       VectAux=XJ+EX
       CALL TORS(XI,XJ,VectAux,XK,1.D0,Active,Value_O=Value1) 
     ENDIF
     IF(PRESENT(Value2)) THEN
      !CALL BEND(XI,XJ,XK,Value_O=Value2,W_O=EY)
       VectAux=XJ+EY
       CALL TORS(XI,XJ,VectAux,XK,1.D0,Active,Value_O=Value2) 
     ENDIF
     !
     IF(PRESENT(BB1)) THEN
      !CALL BEND(XI,XJ,XK,BB_O=BB1,W_O=EX)
       VectAux=XJ+EX
       CALL TORS(XI,XJ,VectAux,XK,1.D0,Active,BB_O=BB1) 
       BB1(4:6)=BB1(4:6)+BB1(7:9)
       BB1(7:9)=BB1(10:12)
       BB1(10:12)=Zero
     ENDIF
     !
     IF(PRESENT(BB2)) THEN
      !CALL BEND(XI,XJ,XK,BB_O=BB2,W_O=EY)
       VectAux=XJ+EY
       CALL TORS(XI,XJ,VectAux,XK,1.D0,Active,BB_O=BB2) 
       BB2(4:6)=BB2(4:6)+BB2(7:9)
       BB2(7:9)=BB2(10:12)
       BB2(10:12)=Zero
     ENDIF
     !
   END SUBROUTINE LinB
!
!--------------------------------------------------------------------
!
   SUBROUTINE BCART(CharU,B,Constraint)
     INTEGER :: I,J
     REAL(DOUBLE),DIMENSION(1:12) :: B
     CHARACTER :: CharU
     LOGICAL   :: Constraint
     !
     B=Zero
     IF(CharU=='X') THEN
       IF(Constraint) THEN
         B(1)=Zero 
       ELSE
         B(1)=One
       ENDIF
     ENDIF
     !
     IF(CharU=='Y') THEN
       IF(Constraint) THEN
         B(2)=Zero 
       ELSE
         B(2)=One  
       ENDIF
     ENDIF
     !
     IF(CharU=='Z') THEN
       IF(Constraint) THEN
         B(3)=Zero 
       ELSE
         B(3)=One 
       ENDIF
     ENDIF
   END SUBROUTINE BCART
!
!----------------------------------------------------------------
!
   SUBROUTINE VOLUME(XYZO,XYZA,XYZB,XYZC,Active,Value_O,BB_O)
     REAL(DOUBLE),DIMENSION(3) :: XYZO,XYZA,XYZB,XYZC
     LOGICAL                            :: Active
     REAL(DOUBLE),OPTIONAL              :: Value_O
     REAL(DOUBLE)                       :: Value
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: BB_O                   
     REAL(DOUBLE),DIMENSION(3)          :: VA,VB,VC
     REAL(DOUBLE),DIMENSION(3)          :: AxB,BxC,CxA            
     !
     VA=XYZA-XYZO
     VB=XYZB-XYZO
     VC=XYZC-XYZO
     CALL CROSS_PRODUCT(VB,VC,BxC) 
     IF(PRESENT(Value_O)) THEN
       Value_O=DOT_PRODUCT(VA,BxC)
     ENDIF
     IF(PRESENT(BB_O)) THEN
       CALL CROSS_PRODUCT(VA,VB,AxB) 
       CALL CROSS_PRODUCT(VC,VA,CxA) 
       BB_O(1:3)=-(BxC+CxA+AxB)
       BB_O(4:6)=BxC
       BB_O(7:9)=CxA
       BB_O(10:12)=AxB
     ENDIF
   END SUBROUTINE VOLUME
!
!----------------------------------------------------------------
!
   SUBROUTINE AREA(XYZO,XYZA,XYZB,Active,Value_O,BB_O)
     REAL(DOUBLE),DIMENSION(3) :: XYZO,XYZA,XYZB,XYZC
     LOGICAL                            :: Active
     REAL(DOUBLE),OPTIONAL              :: Value_O
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: BB_O                   
     REAL(DOUBLE)                       :: AB
     REAL(DOUBLE),DIMENSION(3)          :: VA,VB,BE
     REAL(DOUBLE),DIMENSION(3)          :: AxB,Eab,EA,EB
     !
     VA=XYZA-XYZO
     VB=XYZB-XYZO
     CALL CROSS_PRODUCT(VA,VB,AxB) 
     AB=SQRT(DOT_PRODUCT(AxB,AxB))
     IF(AB<1.D-6) THEN
       Active=.FALSE.
       IF(PRESENT(Value_O)) Value_O=Zero
       IF(PRESENT(BB_O)) BB_O=Zero
       RETURN
     ENDIF
     !
     XYZC=AxB/AB
     IF(PRESENT(Value_O)) THEN
       CALL VOLUME(XYZO,XYZA,XYZB,XYZC,Active,Value_O=Value_O)
     ELSE IF(PRESENT(BB_O)) THEN
       CALL VOLUME(XYZO,XYZA,XYZB,XYZC,Active,BB_O=BB_O)
       BB_O(1:3)=-(BB_O(4:6)+BB_O(7:9))
       BB_O(10:12)=Zero
     ENDIF
     !
   END SUBROUTINE AREA
!
!----------------------------------------------------------------
!
   SUBROUTINE DefineIntCoos(XYZ,AtNum,IntCs,NIntC, &
                            Cells,IEq,GCoordCtrl, &
                            Bond,AtmB,TOPS,GConvCr, &
                            ArchMem_O,HFileIn_O,iCLONE_O,iGEO_O)
     !
     ! This routine defines internal coordinates
     ! being used in geometry manipulations.
     ! They may be STRE, BEND, TORS, OutP and LINB1&LINB2 -s.
     !
     IMPLICIT NONE
     INTEGER                        :: I,J,N,NatmsLoc,NIntC
     INTEGER                        :: ILast
     TYPE(INT_RNK2)                 :: TorsionIJKL
     INTEGER                        :: NTorsion
     REAL(DOUBLE),DIMENSION(:,:)    :: XYZ
     TYPE(INTC)                     :: IntCs,IntC_Frag
     TYPE(ATOMBONDS)                :: AtmBTot,AtmB
     TYPE(BONDDATA)                 :: BondTot,Bond
     CHARACTER(LEN=*),OPTIONAL      :: HFileIn_O
     INTEGER,OPTIONAL               :: iCLONE_O,iGEO_O,ArchMem_O
     INTEGER                        :: MaxBonds
     INTEGER,DIMENSION(:,:)         :: Cells
     INTEGER,DIMENSION(:)           :: AtNum,IEq
     TYPE(CoordCtrl)                :: GCoordCtrl
     TYPE(GConvCrit)                :: GConvCr
     TYPE(TOPOLOGY)                 :: TOPS
     TYPE(ANGLEDATA)                :: Angle
     TYPE(OUTPDATA)                 :: OutP 
     !
     NIntC=0
     NTorsion=0
     NatmsLoc=SIZE(XYZ,2)
     !
     CALL BondingScheme(XYZ,AtNum,AtmBTot,BondTot,TOPS, &
                        GCoordCtrl,GConvCr,Cells,IEq)
  !  CALL SortBonds(NatmsLoc,AtmBTot,BondTot)
  !  CALL Set_BONDDATA_EQ_BONDDATA(Bond,BondTot)
  !  CALL Set_AtmB_EQ_AtmB(AtmB,AtmBTot)
     !
   ! IF(PRESENT(HFileIn_O).AND.PRESENT(iCLONE_O).AND.&
   !    PRESENT(iGEO_O).AND.PRESENT(ArchMem_O)) THEN
   !   CALL ArchiveTop(TOPS,ArchMem_O, &
   !                  BondTot,AtmBTot,HFileIn_O,iCLONE_O,iGEO_O)
   ! ENDIF
     !
     ! Now define bond angles and torsions
     !
     CALL AngleList(AtmBTot,BondTot,TOPS,XYZ,GConvCr%NonCovBend, &
                    Angle,OutP,Cells,IEq)
     CALL TorsionList(NatmsLoc,AtmBTot,BondTot,XYZ, &
                      GConvCr%NonCovTors,TOPS, &
                      AtNum,TorsionIJKL,NTorsion,IEq)
     !
     ! Fill Data into IntCs
     !
     NIntC=BondTot%N+Angle%N+NTorsion+OutP%N
     IF(NIntC/=0) THEN          
       CALL New(IntCs,NIntC)
       !
         ILast=0
       DO I=1,BondTot%N
         IntCs%Def%C(ILast+I)(1:10)='STRE      '
         IntCs%Atoms%I(ILast+I,1:2)=BondTot%IJ%I(1:2,I)
       ENDDO
         ILast=BondTot%N
       DO I=1,Angle%N
         IntCs%Def%C(ILast+I)(1:10)='BEND      '
         IntCs%Atoms%I(ILast+I,1:3)=Angle%IJK%I(1:3,I)
       ENDDO
         ILast=ILast+Angle%N
       DO I=1,NTorsion
         IntCs%Def%C(ILast+I)(1:10)='TORS      '
         IntCs%Atoms%I(ILast+I,1:4)=TorsionIJKL%I(1:4,I)
       ENDDO
       ILast=ILast+NTorsion
       DO I=1,OutP%N
         IntCs%Def%C(ILast+I)(1:10)='OUTP      '
         IntCs%Atoms%I(ILast+I,1:4)=OutP%IJKL%I(1:4,I)
       ENDDO
     ENDIF
     !
     ! tidy up
     !
     CALL Delete(TorsionIJKL)
     !
     CALL Delete(AtmBTot)
     CALL Delete(BondTot)
     CALL Delete(Angle)
     CALL Delete(OutP)
   END SUBROUTINE DefineIntCoos
!
!----------------------------------------------------------------
!
   SUBROUTINE ScoreBond(BondScore,Bond)
     TYPE(BONDDATA)       :: Bond
     INTEGER,DIMENSION(:) :: BondScore
     INTEGER              :: I
     !
     IF(SIZE(BondScore)/=Bond%N) CALL Halt('Dimension error in ScoreBond')
     BondScore=1
     DO I=1,Bond%N
       IF(Bond%Type%C(I)(1:5)=='HBond') BondScore(I)=BondScore(I)+1
       IF(Bond%Type%C(I)(1:6)=='MetLig') BondScore(I)=BondScore(I)+1
     ENDDO
   END SUBROUTINE ScoreBond
!
!----------------------------------------------------------------
!
   SUBROUTINE BondExcl(JJ1,JJ2,NJJ1,NJJ2,TOPS, &
                       FoundHBond,FoundMetLig,LonelyAtom,DoExclude)
     INTEGER             :: JJ1,JJ2,NJJ1,NJJ2
     TYPE(TOPOLOGY)      :: TOPS
     INTEGER             :: IDimExcl,IDim12,IDim13,IDim14
     LOGICAL             :: FoundHBond,FoundMetLig,LonelyAtom,DoExclude
     !
     IDimExcl=TOPS%CovExcl%I(JJ1,1)
     DoExclude=.FALSE.
     IF(ANY(TOPS%CovExcl%I(JJ1,2:1+IDimExcl)==JJ2) &
       .AND.IDimExcl/=0) THEN
       DoExclude=.TRUE.
       RETURN
     ENDIF
   ! IF(ANY(TOPS%CovExcl%I(JJ1,2:1+IDimExcl)==JJ2) &
   !   .AND.IDimExcl/=0) THEN
   !   IF(.NOT.(FoundHBond.OR.FoundMetLig)) THEN
   !     DoExclude=.TRUE.
   !     RETURN
   !   ELSE
   !     IDim12=TOPS%Cov12%I(JJ1,1)
   !     IDim13=TOPS%Cov13%I(JJ1,1)
   !     IF(ANY(TOPS%Cov12%I(JJ1,2:1+IDim12)==JJ2) &
   !       .OR.ANY(TOPS%Cov13%I(JJ1,2:1+IDim13)==JJ2)) THEN
   !       DoExclude=.TRUE.
   !       RETURN
   !     ENDIF
   !     IDim14=TOPS%Cov14%I(JJ1,1)
   !     IF(.NOT. &
   !       (ANY(TOPS%Cov14%I(JJ1,2:1+IDim14)==JJ2) &
   !       .AND.IDim14/=0)) THEN
   !       DoExclude=.TRUE.
   !       RETURN
   !     ENDIF
   !   ENDIF
   ! ENDIF
    !IF(.NOT.(FoundHBond.OR.FoundMetLig.OR.&
    !   LonelyAtom)) THEN
    !  DoExclude=.TRUE.
    !  RETURN
    !ENDIF
   END SUBROUTINE BondExcl
!
!----------------------------------------------------------------
!
   SUBROUTINE MoreBondArray(Bond,Incr,DimOld)
     TYPE(BondDATA)     :: Bond,Bond2
     INTEGER            :: DimNew,DimOld,Incr
     !
     DimNew=DimOld+Incr
     CALL Set_BONDDATA_EQ_BONDDATA(Bond2,Bond, &
                                   NewDim_O=DimNew,OldDim_O=DimOld)
     CALL Delete(Bond)
     CALL New(Bond,DimNew)
     CALL Set_BONDDATA_EQ_BONDDATA(Bond,Bond2)
     CALL Delete(Bond2)
   END SUBROUTINE MoreBondArray
!
!----------------------------------------------------------------
!
   SUBROUTINE TorsionList(NatmsLoc,AtmB,Bond,XYZ,NonCovTors, &
                          TOPS,AtNum,TorsionIJKL,NTorsion,IEq)
     IMPLICIT NONE
     INTEGER                     :: NatmsLoc,I,I1,I2
     INTEGER                     :: N1,N2,J,J1,J2,NBond
     INTEGER                     :: NTorsion
     INTEGER                     :: I1Num,I2Num,JJ,I3,N3,JJ1,JJ2,NTIni
     TYPE(ATOMBONDS)             :: AtmB
     TYPE(BONDDATA)              :: Bond
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(INT_RNK2)              :: TorsionIJKL,TorsionIJKLAux
     LOGICAL                     :: SelectTors
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Crit,Alph1,Alph2,SumU,PIHalf
     REAL(DOUBLE)                :: Sum1,Sum2
     INTEGER,DIMENSION(1:4)      :: Atoms
     INTEGER,DIMENSION(:)        :: AtNum,IEq
     INTEGER                     :: I1JJ1,I2JJ2
     LOGICAL                     :: ExcludeVDW,NonCovTors
     !    
     ExcludeVDW=.NOT.NonCovTors
     PIHalf=PI*Half
    !SelectTors=.FALSE.
     SelectTors=.TRUE.
     NBond=SIZE(Bond%IJ%I,2)
     NTorsion=0
     DO I=1,NBond
       I1=Bond%IJ%I(1,I)
       I2=Bond%IJ%I(2,I)
       N1=AtmB%Count%I(I1)
       N2=AtmB%Count%I(I2)
       NTorsion=NTorsion+N1*N2
     ENDDO
     !    
     CALL New(TorsionIJKLAux,(/4,NTorsion/))
     !    
     NTorsion=0
     DO I=1,NBond 
       I1=Bond%IJ%I(1,I)
       I2=Bond%IJ%I(2,I)
       IF((ExcludeVDW.AND.Bond%Type%C(I)(1:3)/='COV').AND. &
          (TOPS%Cov12%I(I1,1)/=0.AND.TOPS%Cov12%I(I2,1)/=0)) CYCLE 
       N1=AtmB%Count%I(I1)
       N2=AtmB%Count%I(I2)
       NTIni=NTorsion
       Crit=Zero
       Atoms=0 
       DO J1=1,N1
         JJ1=AtmB%Atoms%I(I1,J1)
         I1JJ1=AtmB%Bonds%I(I1,J1)
         IF((ExcludeVDW.AND.Bond%Type%C(I1JJ1)(1:3)/='COV').AND. &
            (TOPS%Cov12%I(JJ1,1)/=0)) CYCLE 
         IF(SelectTors) THEN
           CALL BEND(XYZ(1:3,JJ1),XYZ(1:3,I1),XYZ(1:3,I2), &
                     Value_O=Alph1)
           Sum1=DBLE((AtmB%Count%I(JJ1)+1)*AtNum(JJ1))/ABS(Alph1-PIHalf+1.D-10)
         ENDIF
         DO J2=1,N2
           JJ2=AtmB%Atoms%I(I2,J2)
           I2JJ2=AtmB%Bonds%I(I2,J2)
           IF((ExcludeVDW.AND.Bond%Type%C(I2JJ2)(1:3)/='COV').AND. & 
              (TOPS%Cov12%I(JJ2,1)/=0)) CYCLE 
           IF(JJ1==JJ2.OR.JJ1==I2.OR.JJ2==I1) CYCLE
           IF(SelectTors) THEN
             CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,JJ2), &
                       Value_O=Alph2)
             Sum2=DBLE((AtmB%Count%I(JJ2)+1)*AtNum(JJ2))/ABS(Alph2-PIHalf+1.D-10)
             SumU=Sum1+Sum2
             IF(SumU>Crit) THEN
               Crit=SumU
               Atoms=(/JJ1,I1,I2,JJ2/)
             ENDIF
           ELSE
             NTorsion=NTorsion+1
             TorsionIJKLAux%I(1:4,NTorsion)=(/JJ1,I1,I2,JJ2/)
           ENDIF
         ENDDO
       ENDDO
       IF(SelectTors.AND.Atoms(1)/=0) THEN
         NTorsion=NTorsion+1
         TorsionIJKLAux%I(1:4,NTorsion)=Atoms(1:4)
       ENDIF
     ENDDO
     !
     CALL New(TorsionIJKL,(/4,NTorsion/))
     DO I=1,NTorsion
       IF(IEq(TorsionIJKLAux%I(1,I))<IEq(TorsionIJKLAux%I(4,I))) THEN
         DO J=1,4 ; TorsionIJKL%I(J,I)=TorsionIJKLAux%I(J,I) ; ENDDO
       ELSE
         DO J=1,4 ; TorsionIJKL%I(5-J,I)=TorsionIJKLAux%I(J,I) ; ENDDO
       ENDIF
     ENDDO
     CALL Delete(TorsionIJKLAux)
   END SUBROUTINE TorsionList
!
!-----------------------------------------------------------------------
!
   SUBROUTINE GetIntCs(XYZ,AtNumIn,IntCs,Refresh,SCRPath, &
                       CtrlCoord,CtrlConstr,GConvCr, &
                       ArchMem,IntC_Extra, &
                       TOPS,Bond,AtmB,PBCDim,iGEO,HFileIn_O,iCLONE_O)
     !
     ! This subroutine constructs the IntCs array, which holds
     ! definitions of internal coordinates to be used in the 
     ! forthcoming geometry manipulation procedure.
     ! Refresh=1 : Refresh all definitions
     !        =2 : Refresh only definitions based on VDW interaction
     !        =3 : Do not refresh definitions, use the one from HDF
     !        =4 : Refresh/generate only the covalent coordinates
     !        =5 : only the extra coordinates from input
     ! 
     IMPLICIT NONE
     TYPE(INTC)                  :: IntCs,IntC_Bas,IntC_VDW
     TYPE(INTC)                  :: IntC_Extra,IntC_New
     TYPE(INTC)                  :: IntC_L
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(CoordCtrl)             :: CtrlCoord
     TYPE(Constr)                :: CtrlConstr
     TYPE(GConvCrit)             :: GConvCr    
     LOGICAL                     :: DoFixMM
     REAL(DOUBLE),DIMENSION(:)   :: AtNumIn
     CHARACTER(LEN=*),OPTIONAL   :: HFileIn_O
     INTEGER,OPTIONAL            :: iCLONE_O
     INTEGER                     :: PBCDim,iGEO
     INTEGER                     :: NIntC,NIntC_Bas,NIntC_VDW
     INTEGER                     :: NIntC_Extra,NNew,Nintc_New
     INTEGER                     :: MaxBonds,ArchMem
     INTEGER                     :: I,J,K,Refresh,NatmsLoc,II,ILast,III
     INTEGER                     :: I1,I2,I3,I4,NMax12
     INTEGER                     :: NStreGeOp,NBendGeOp
     INTEGER                     :: NLinBGeOp,NOutPGeOp,NTorsGeOp
     TYPE(INT_VECT)              :: AtNum,AtNumRepl
     TYPE(INT_RNK2)              :: Top12
     TYPE(INT_RNK3)              :: CellEq
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(IntCBox)               :: Box
     TYPE(DBL_RNK2)              :: XYZRepl
     TYPE(INT_VECT)              :: IEq
     TYPE(INT_RNK2)              :: Cells
     !
     NatmsLoc=SIZE(XYZ,2)
     NIntC_Bas=0
     NIntC_VDW=0
!!! restart may cause filled up IntCs; 
!!! this here is a temporary solution, until restart gets reprganized
     IF(AllocQ(IntCs%Alloc)) CALL Delete(IntCs) 
     !
     CALL New(AtNum,NatmsLoc)
     DO I=1,NatmsLoc
       AtNum%I(I)=INT(AtNumIn(I))
     ENDDO
     !
     CALL PrepCells(XYZ,AtNum%I,PBCDim,XYZRepl,AtNumRepl, &
                    Cells,CellEq,IEq)
!do i=1,size(XYZRepl%D,2)
!write(*,100) AtNumRepl%I(i),XYZRepl%D(1:3,i)/angstromstoau
!write(out,100) AtNumRepl%I(i),XYZRepl%D(1:3,i)/angstromstoau
!enddo
!100 format(I4,3F20.8)
!stop
     !
     IF(Refresh==1) THEN !!! Total refresh
       IF(PRESENT(HFileIn_O).AND.PRESENT(iCLONE_O)) THEN
         CALL DefineIntCoos(XYZRepl%D,AtNumRepl%I,IntC_Bas,NIntC_Bas, &
                            Cells%I,IEq%I,CtrlCoord,Bond,AtmB,TOPS, &
                            GConvCr,ArchMem_O=ArchMem, &
                            HFileIn_O=HFileIn_O,iCLONE_O=iCLONE_O, &
                            iGEO_O=iGEO)  
       ELSE
         CALL DefineIntCoos(XYZRepl%D,AtNumRepl%I,IntC_Bas,NIntC_Bas, &
                            Cells%I,IEq%I,CtrlCoord,Bond,AtmB,TOPS,& 
                            GConvCr)  
       ENDIF
       !
       ! Check for bending - lin.bending transions
       ! also for long range torsions!
       !
       CALL ChkBendToLinB(IntC_Bas,XYZRepl%D,AtNumRepl%I, &
                          CtrlCoord,TOPS,IEq%I)
       CALL CleanPBCIntCs(IntC_Bas,Cells%I,IEq%I)
       NIntC_Bas=IntC_Bas%N
       IF(GConvCr%ExplLatt) THEN
         CALL LatticeINTC(IntC_L,PBCDim)
        !CALL LatticeINTC(IntC_L,PBCDim,VolumeOnly_O=.TRUE.)
        !CALL ExtraLattice(IntC_L,PBCDim)
       ELSE
         IntC_L%N=0
       ENDIF
     ELSE IF(Refresh==5) THEN !!! use only extra coords from input
       NIntC_Bas=0 
       NIntC_VDW=0 
     ELSE
       CALL Halt('Unknown Refresh option in GetIntCs') 
     ENDIF
     !
     CALL Delete(Cells)
     CALL Delete(CellEq)
     CALL Delete(IEq)
     CALL Delete(XYZRepl)
     CALL Delete(AtNumRepl)
     !
     ! Merge INTC arrays
     !
     CtrlCoord%NCov=NIntC_Bas
     NIntC=NIntC_Bas+NIntC_VDW+IntC_Extra%N+IntC_L%N
     !
     CALL New(IntCs,NIntC)
     !
       ILast=0
     IF(NIntC_Bas/=0) THEN
       CALL Set_INTC_EQ_INTC(IntC_Bas,IntCs,1,NIntC_Bas,ILast+1)
       CALL Delete(IntC_Bas)
       ILast=ILast+NIntC_Bas
     ENDIF
     IF(NIntC_VDW/=0) THEN 
       CALL Set_INTC_EQ_INTC(IntC_VDW,IntCs,1,NIntC_VDW,ILast+1)
       CALL Delete(IntC_VDW)
       ILast=ILast+NIntC_VDW
     ENDIF
     IF(IntC_L%N/=0) THEN
       CALL Set_INTC_EQ_INTC(IntC_L,IntCs,1,IntC_L%N,ILast+1)
       CALL Delete(IntC_L)
       ILast=ILast+IntC_L%N
     ENDIF
     !
     ! Set active all internal coords defd so far.
     ! 'Linear torsions will be deactivated later when 
     ! their value is calculated
     ! by some subroutines (eg. INTCValue).
     ! The set of active coordinates may vary 
     ! during the process of optimization.
     !
   ! IntCs%Active%L=.TRUE.
     !
     IF(IntC_Extra%N/=0) THEN
       CALL Set_INTC_EQ_INTC(IntC_Extra,IntCs,1,IntC_Extra%N,ILast+1)
       ILast=ILast+IntC_Extra%N
     ENDIF
     !
     IF(.NOT.(IntCs%N==0.OR.Refresh==5)) THEN
       CALL CleanINTC(IntCs,NIntC_Bas,NIntC_VDW,IntC_Extra%N)
     ENDIF             
     !
     ! Count number of different internal coord types
     !
     NStreGeOp=0;NBendGeOp=0;NLinBGeOp=0;NOutPGeOp=0;NTorsGeOp=0
     DO I=1,IntCs%N
        IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
          NStreGeOp=NStreGeOp+1
        ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
          NBendGeOp=NBendGeOp+1
        ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
          NLinBGeOp=NLinBGeOp+1
        ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
          NOutPGeOp=NOutPGeOp+1
        ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
          NTorsGeOp=NTorsGeOp+1
        ENDIF
     ENDDO
     CtrlCoord%NStre=NStreGeOp
     CtrlCoord%NBend=NBendGeOp
     CtrlCoord%NLinB=NLinBGeOp
     CtrlCoord%NOutP=NOutPGeOp
     CtrlCoord%NTors=NTorsGeOp
     !
     CALL Delete(AtNum)
   ! CALL PrtIntCoords(IntCs,IntCs%Value%D,'GetIntC Internals',PBCDim_O=PBCDim)
   END SUBROUTINE GetIntCs
!
!---------------------------------------------------------------------
!
   SUBROUTINE ExtraLattice(IntC_L,PBCDim)
     TYPE(INTC)                  :: IntC_L
     INTEGER                     :: I,J,K,L,LL,LRef,II,KK,PBCDim
     !
     IF(PBCDim==0) THEN
       IntC_L%N=0
       RETURN
     ENDIF
     !
     LL=0
     LRef=1
     II=0
     !
     IF(PBCDim==3) THEN
       CALL New(IntC_L,16)
       IntC_L%Atoms%I(1:16,1:2)=LRef
       !
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_A1 '
       IntC_L%Cells%I(II,1:6)=(/0,0,0,1,0,0/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_A2 '
       IntC_L%Cells%I(II,1:6)=(/0,1,0,1,1,0/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_A3 '
       IntC_L%Cells%I(II,1:6)=(/0,0,1,1,0,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_A4 '
       IntC_L%Cells%I(II,1:6)=(/0,1,1,1,1,1/)
       !
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_B1 '
       IntC_L%Cells%I(II,1:6)=(/0,0,0,0,1,0/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_B2 '
       IntC_L%Cells%I(II,1:6)=(/0,0,1,0,1,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_B3 '
       IntC_L%Cells%I(II,1:6)=(/1,0,1,1,1,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_B4 '
       IntC_L%Cells%I(II,1:6)=(/1,0,0,1,1,0/)
       !
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_C1 '
       IntC_L%Cells%I(II,1:6)=(/0,0,0,0,0,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_C2 '
       IntC_L%Cells%I(II,1:6)=(/1,0,0,1,0,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_C3 '
       IntC_L%Cells%I(II,1:6)=(/1,1,0,1,1,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_C4 '
       IntC_L%Cells%I(II,1:6)=(/0,1,0,0,1,1/)
       !
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_X1 '
       IntC_L%Cells%I(II,1:6)=(/0,0,0,1,1,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_X2 '
       IntC_L%Cells%I(II,1:6)=(/1,0,0,0,1,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_X3 '
       IntC_L%Cells%I(II,1:6)=(/1,1,0,0,0,1/)
         II=II+1
         IntC_L%Def%C(II)(1:8)='STRE_X4 '
       IntC_L%Cells%I(II,1:6)=(/0,1,0,1,0,1/)
       !
      ! !IntC_L%Def%C(7)(1:8)='VOLM_L  '
      !  IntC_L%Def%C(7)(1:8)='BLANK   '
      !IntC_L%Atoms%I(7,1:4)=LRef
      !IntC_L%Cells%I(7,1:12)=(/0,0,0,1,0,0,0,1,0,0,0,1/)
     ELSE IF(PBCDim==2) THEN
       CALL New(IntC_L,3)
       IntC_L%Atoms%I(1:2,1:2)=LRef
         IntC_L%Def%C(1)(1:8)='STRE_A  '
       IntC_L%Cells%I(1,1:6)=(/0,0,0,1,0,0/)
         IntC_L%Def%C(2)(1:8)='STRE_B  '
       IntC_L%Cells%I(2,1:6)=(/0,0,0,0,1,0/)
       !
         IntC_L%Def%C(3)(1:8)='GAMMA   '
       IntC_L%Atoms%I(3,1:3)=LRef
       IntC_L%Cells%I(3,1:9)=(/1,0,0,0,0,0,0,1,0/)
       !
      ! !IntC_L%Def%C(4)(1:8)='AREA_L  '
      !  IntC_L%Def%C(4)(1:8)='BLANK   '
      !IntC_L%Atoms%I(4,1:3)=LRef
      !IntC_L%Cells%I(4,1:9)=(/0,0,0,1,0,0,0,1,0/)
     ELSE IF(PBCDim==1) THEN
       CALL New(IntC_L,1)
         IntC_L%Def%C(1)(1:8)='STRE_A  '
       IntC_L%Atoms%I(1,1:2)=LRef
       IntC_L%Cells%I(1,1:6)=(/0,0,0,1,0,0/)
     ENDIF
     IntC_L%Constraint%L(:)=.FALSE.
     IntC_L%Active%L(:)=.TRUE.
   END SUBROUTINE ExtraLattice
!
!-------------------------------------------------------------------
!
   SUBROUTINE LatticeINTC(IntC_L,PBCDim,DoVolume_O,VolumeOnly_O)
     TYPE(INTC)                  :: IntC_L
     INTEGER                     :: I,J,K,L,LL,LRef,II,KK,PBCDim,IVol
     INTEGER                     :: NDim
     LOGICAL,OPTIONAL            :: DoVolume_O,VolumeOnly_O
     LOGICAL                     :: DoVolume,VolumeOnly
     !
     IF(PBCDim==0) THEN
       IntC_L%N=0
       RETURN
     ENDIF
     DoVolume=.FALSE.
     IF(PRESENT(DoVolume_O)) DoVolume=DoVolume_O
     VolumeOnly=.FALSE.
     IF(PRESENT(VolumeOnly_O)) VolumeOnly=VolumeOnly_O
     IVol=0
     IF(DoVolume) IVol=1
     !
     LL=0
     LRef=1
     !
     IF(PBCDim==3) THEN
       IF(.NOT.VolumeOnly) THEN
         CALL New(IntC_L,6+IVol)
         IntC_L%Constraint%L(:)=.FALSE.
         IntC_L%Active%L(:)=.TRUE.
         IntC_L%Atoms%I(1:3,1:2)=LRef
           IntC_L%Def%C(1)(1:8)='STRE_A  '
         IntC_L%Cells%I(1,1:6)=(/0,0,0,1,0,0/)
           IntC_L%Def%C(2)(1:8)='STRE_B  '
         IntC_L%Cells%I(2,1:6)=(/0,0,0,0,1,0/)
           IntC_L%Def%C(3)(1:8)='STRE_C  '
         IntC_L%Cells%I(3,1:6)=(/0,0,0,0,0,1/)
         !
         IntC_L%Atoms%I(4:6,1:3)=LRef
           IntC_L%Def%C(4)(1:8)='ALPHA   '
         IntC_L%Cells%I(4,1:9)=(/0,1,0,0,0,0,0,0,1/)
           IntC_L%Def%C(5)(1:8)='BETA    '
         IntC_L%Cells%I(5,1:9)=(/1,0,0,0,0,0,0,0,1/)
           IntC_L%Def%C(6)(1:8)='GAMMA   '
         IntC_L%Cells%I(6,1:9)=(/1,0,0,0,0,0,0,1,0/)
         !
         IF(DoVolume) THEN
           IntC_L%Def%C(7)(1:8)='VOLM_L  '
           IntC_L%Atoms%I(7,1:4)=LRef
           IntC_L%Cells%I(7,1:12)=(/0,0,0,1,0,0,0,1,0,0,0,1/)
         ENDIF
       ELSE
         CALL New(IntC_L,1)
         IntC_L%Constraint%L(:)=.FALSE.
         IntC_L%Active%L(:)=.TRUE.
         IntC_L%Def%C(1)(1:8)='VOLM_L  '
         IntC_L%Atoms%I(1,1:4)=LRef
         IntC_L%Cells%I(1,1:12)=(/0,0,0,1,0,0,0,1,0,0,0,1/)
       ENDIF
     ELSE IF(PBCDim==2) THEN
       IF(.NOT.VolumeOnly) THEN
         CALL New(IntC_L,3+IVol)
         IntC_L%Constraint%L(:)=.FALSE.
         IntC_L%Active%L(:)=.TRUE.
         IntC_L%Atoms%I(1:2,1:2)=LRef
           IntC_L%Def%C(1)(1:8)='STRE_A  '
         IntC_L%Cells%I(1,1:6)=(/0,0,0,1,0,0/)
           IntC_L%Def%C(2)(1:8)='STRE_B  '
         IntC_L%Cells%I(2,1:6)=(/0,0,0,0,1,0/)
         !
           IntC_L%Def%C(3)(1:8)='GAMMA   '
         IntC_L%Atoms%I(3,1:3)=LRef
         IntC_L%Cells%I(3,1:9)=(/1,0,0,0,0,0,0,1,0/)
         !
         IF(DoVolume) THEN
             IntC_L%Def%C(4)(1:8)='AREA_L  '
            !IntC_L%Def%C(4)(1:8)='BLANK   '
           IntC_L%Atoms%I(4,1:3)=LRef
           IntC_L%Cells%I(4,1:9)=(/0,0,0,1,0,0,0,1,0/)
         ENDIF
       ELSE
         CALL New(IntC_L,1)
         IntC_L%Constraint%L(:)=.FALSE.
         IntC_L%Active%L(:)=.TRUE.
           IntC_L%Def%C(1)(1:8)='AREA_L  '
         IntC_L%Atoms%I(1,1:3)=LRef
         IntC_L%Cells%I(1,1:9)=(/0,0,0,1,0,0,0,1,0/)
       ENDIF
     ELSE IF(PBCDim==1) THEN
       CALL New(IntC_L,1)
         IntC_L%Def%C(1)(1:8)='STRE_A  '
       IntC_L%Constraint%L(:)=.FALSE.
       IntC_L%Active%L(:)=.TRUE.
       IntC_L%Atoms%I(1,1:2)=LRef
       IntC_L%Cells%I(1,1:6)=(/0,0,0,1,0,0/)
     ENDIF
   END SUBROUTINE LatticeINTC
!
!-------------------------------------------------------------------
!
   SUBROUTINE CleanPBCIntCs(IntCs,Cells,IEq)
     TYPE(INTC)                  :: IntCs,IntCs2
     INTEGER,DIMENSION(:,:)      :: Cells
     INTEGER,DIMENSION(:)        :: IEq  
     TYPE(INT_VECT)              :: Sort
     INTEGER                     :: I,I1,I2,J,K,L,N,II,M
     INTEGER                     :: K1,K2,J1,J2
     INTEGER,DIMENSION(3)        :: Cell1,Tr
     INTEGER,DIMENSION(4)        :: DotProds
     LOGICAL                     :: DOAllow,AllCentral
     !
     ! Later, to avoid choices on tors and angles and filters
     ! a sparse overlap matrix based filtering would be advantageous
     ! here
     !
     CALL New(Sort,IntCs%N)
     Sort%I=0
     !
     II=0
     DO I=1,IntCs%N
       DoAllow=.FALSE.
       K1=1
       K2=4
       IF(IntCs%Def%C(I)(1:4)=='LINB') K2=3
     ! IF(IntCs%Def%C(I)(1:4)=='TORS') THEN 
     !   K1=2
     !   K2=3
     ! ENDIF
     ! IF(IntCs%Def%C(I)(1:4)=='BEND'.OR. &
     !    IntCs%Def%C(I)(1:4)=='OUTP'.OR. &
     !    IntCs%Def%C(I)(1:4)=='LINB') THEN
     !   K1=2
     !   K2=2
     ! ENDIF
       Tr=0
       DO J=K1,K2
         K=IntCs%Atoms%I(I,J)
         IF(K==0) EXIT
         IF(ALL(Cells(K,1:3)==0)) THEN
           DoAllow=.TRUE.
         ELSE
        !  IF(Cells(K,1)<Tr(1)) Tr(1)=Cells(K,1)
        !  IF(Cells(K,2)<Tr(2)) Tr(2)=Cells(K,2)
        !  IF(Cells(K,3)<Tr(3)) Tr(3)=Cells(K,3)
         ENDIF
       ENDDO
       IF(.NOT.DoAllow) CYCLE
       !
     ! IF(ANY(Tr<0)) THEN
     !   IntCs%Cells%I(I,1:12)=0    
     !   DO J=K1,K2
     !     K=IntCs%Atoms%I(I,J)
     !     IF(K==0) EXIT
     !     J1=3*(J-1)+1
     !     J2=J1+2
     !     ! translate coordinate to positive space eighth
     !     IntCs%Cells%I(I,J1:J2)=Cells(K,1:3)-Tr
     !   ENDDO
     !   !
     !   DoAllow=.TRUE.
     !   DO J=K1,K2
     !     K=IntCs%Atoms%I(I,J)
     !     IF(K==0) EXIT
     !     J1=3*(J-1)+1
     !     J2=J1+2
     !     IF(ALL(IntCs%Cells%I(I,J1:J2)==0)) THEN
     !       DoAllow=.FALSE.
     !       EXIT
     !     ENDIF
     !   ENDDO
     !   IF(DoAllow.AND.(Tr(2)/=0.OR.Tr(3)/=0)) DoAllow=.FALSE.
     !   IF(.NOT.DoAllow) CYCLE
     ! ENDIF
       !
       II=II+1
       Sort%I(I)=1
       !
       IntCs%Cells%I(I,1:12)=0    
       DO J=1,4
         I1=IntCs%Atoms%I(I,J)
         IF(I1==0) Exit
         IntCs%Atoms%I(I,J)=IEq(I1)
         K=(J-1)*3+1
         L=K+2
         IntCs%Cells%I(I,K:L)=Cells(I1,1:3)
       ENDDO
     ENDDO
     !
     CALL New(IntCs2,II)
     II=0
     DO I=1,IntCs%N
       IF(Sort%I(I)==0) CYCLE
       II=II+1
       CALL Set_INTC_EQ_INTC(IntCs,IntCs2,I,I,II)
     ENDDO
     !
     CALL Delete(IntCs)
     CALL New(IntCs,II)
     CALL Set_INTC_EQ_INTC(IntCs2,IntCs,1,II,1)
     CALL Delete(IntCs2)
     CALL Delete(Sort)
   END SUBROUTINE CleanPBCIntCs
!
!-------------------------------------------------------------------
!
   SUBROUTINE PrepCells(XYZ,AtNum,PBCDim,XYZBig,AtNumRepl, &
                        Cells,CellEq,IEq,Dir_O,NDim_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INT_VECT)              :: IEq,AtNumRepl
     INTEGER,DIMENSION(:)        :: AtNum
     TYPE(INT_RNK2)              :: Cells
     TYPE(INT_RNK3)              :: CellEq
     INTEGER                     :: NDim,NA1,NB1,NC1,NA2,NB2,NC2
     REAL(DOUBLE)                :: XTrans,YTrans,ZTrans
     INTEGER                     :: NatmsLoc,I,J,II,NA,NB,NC
     INTEGER                     :: NCells,IA,IB,IC,PBCDim
     TYPE(DBL_RNK2)              :: XYZBig
     CHARACTER(LEN=1),OPTIONAL   :: Dir_O
     INTEGER,OPTIONAL            :: NDim_O
     !
     NatmsLoc=SIZE(XYZ,2)-3 ! last 3 row are PBC data
     !
     ! Grow a layer of 16 a.u. (appr. 8 A) around central cell
     ! To be able to model far reaching internal coordinates
     !
    !IF(NatmsLoc==0) THEN
    !  CALL Halt('No atoms in prepcell.')
    !ELSE IF(NatmsLoc==1) THEN
    !   NDim=2
    !ELSE IF(NatmsLoc>1) THEN
    !   NDim=1
    !ENDIF
     NDim=2
     !
     CALL New(CellEq,(/NDim,NDim,NDim/),M_O=(/-NDim,-NDim,-NDim/))
     CellEq%I=0
     !
     NA=0 ; NB=0 ; NC=0
     IF(PBCDim>0)  NA=NDim
     IF(PBCDim>1)  NB=NDim
     IF(PBCDim>2)  NC=NDim
     IF(PRESENT(Dir_O)) THEN
       IF(Dir_O=='A') THEN
         NA=1 ; NB=0 ; NC=0
       ENDIF
       IF(Dir_O=='B') THEN
         NA=0 ; NB=1 ; NC=0
       ENDIF
       IF(Dir_O=='C') THEN
         NA=0 ; NB=0 ; NC=1
       ENDIF
     ENDIF
     NA1=-NA ; NA2=NA
     NB1=-NB ; NB2=NB
     NC1=-NC ; NC2=NC
     IF(PRESENT(NDim_O)) THEN
       IF(NDim_O<0) THEN
         NA1=0 ; NA2=1
         NB1=0 ; NB2=1
         NC1=0 ; NC2=1
       ELSE
         NA1=-NDim_O
         NB1=-NDim_O
         NC1=-NDim_O
       ENDIF
     ENDIF
!NA1=0 ; NA2=3
!NB1=0 ; NB2=3
!NC1=0 ; NC2=3
     !
     NCells=(NA2-NA1+1)*(NB2-NB1+1)*(NC2-NC1+1)
     CALL New(Cells,(/NCells*NatmsLoc,3/))
     Cells%I=0
     CALL New(IEq,NCells*NatmsLoc)
     IEq%I=0
     CALL New(AtNumRepl,NCells*NatmsLoc)
     CALL New(XYZBig,(/3,NCells*NatmsLoc/))
     !
     II=0
     CellEq%I(0,0,0)=II+1
     DO I=1,NatmsLoc 
       II=II+1
       XYZBig%D(1:3,II)=XYZ(1:3,I) 
       Cells%I(II,1)=0 
       Cells%I(II,2)=0 
       Cells%I(II,3)=0 
       IEq%I(II)=I ! equivalent atom from central cell
       AtNumRepl%I(II)=AtNum(I) 
     ENDDO
     !
     DO IA=NA1,NA2
       DO IB=NB1,NB2
         DO IC=NC1,NC2
           IF(IA==0.AND.IB==0.AND.IC==0) CYCLE
           CellEq%I(IA,IB,IC)=II+1
           DO I=1,NatmsLoc
             XTrans=IA*XYZ(1,NatmsLoc+1)+&
                   IB*XYZ(1,NatmsLoc+2)+IC*XYZ(1,NatmsLoc+3)
             YTrans=IA*XYZ(2,NatmsLoc+1)+&
                   IB*XYZ(2,NatmsLoc+2)+IC*XYZ(2,NatmsLoc+3)
             ZTrans=IA*XYZ(3,NatmsLoc+1)+&
                   IB*XYZ(3,NatmsLoc+2)+IC*XYZ(3,NatmsLoc+3)
             II=II+1
             XYZBig%D(1,II)=XYZ(1,I)+XTrans
             XYZBig%D(2,II)=XYZ(2,I)+YTrans
             XYZBig%D(3,II)=XYZ(3,I)+ZTrans
             Cells%I(II,1)=IA
             Cells%I(II,2)=IB
             Cells%I(II,3)=IC
             IEq%I(II)=I ! equivalent atom from central cell
             AtNumRepl%I(II)=AtNum(I) 
           ENDDO
         ENDDO
       ENDDO
     ENDDO
   END SUBROUTINE PrepCells
!
!-------------------------------------------------------------------
!
   SUBROUTINE CleanINTC(IntCs,NIntC_Cov,NIntC_VDW,NIntC_Extra)
     !
     ! Now filter out repeated definitions, 
     ! which may have occured in INTC_Extra
     ! keep those values, which were defd. in extras.
     ! Do not check for Cartesians.
     !
     TYPE(INTC) :: IntCs,IntC_New
     INTEGER    :: NIntC_Cov,NIntC_VDW,NIntC,I,J,NNew,II,III,ILast
     INTEGER    :: NIntC_Extra
     INTEGER    :: Cell(12),KK,Atoms(4),LL
     !
     NIntC=IntCs%N
     ILast=NIntC_Cov+NIntC_VDW
     II=0
   ! DO III=1,NIntC
       DO J=1,ILast
         IF(IntCs%Def%C(J)(1:5)=='BLANK') THEN
           II=II+1
           CYCLE
         ENDIF
         DO I=ILast+1,NIntC
           IF(IntCs%Def%C(I)(1:5)=='CART ') CYCLE
           Cell=IntCs%Cells%I(J,:)-IntCs%Cells%I(I,:)
           KK=DOT_PRODUCT(Cell,Cell)
           IF(KK/=0) CYCLE
           Atoms=IntCs%Atoms%I(J,:)-IntCs%Atoms%I(I,:)
           LL=DOT_PRODUCT(Atoms,Atoms)
           IF(LL/=0) CYCLE 
           !
           II=II+1
           IntCs%Def%C(J)(1:5)='BLANK'
           EXIT
         ENDDO
       ENDDO
   ! ENDDO
     !
     ! Compress IntCs array, get rid of BLANK-s
     !
     IF(II/=0) THEN
       IF(ANY(IntCs%Def%C(:)(1:5)=='BLANK')) THEN
         NNew=NIntC
         DO I=1,NIntc
           IF(IntCs%Def%C(I)(1:5)=='BLANK') NNew=NNew-1
         ENDDO
         CALL New(IntC_New,NNew)
         NNew=0    
         DO I=1,NIntc
           IF(IntCs%Def%C(I)(1:5)/='BLANK') THEN
             NNew=NNew+1
             CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NNew)
           ENDIF
         ENDDO
         CALL Delete(IntCs)
         NIntC=NNew
         CALL New(IntCs,NIntC)
           CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntC,1)
         CALL Delete(IntC_New)
       ENDIF
     ENDIF
   END SUBROUTINE CleanINTC
!
!-------------------------------------------------------
!
   SUBROUTINE INTCValue(IntCs,XYZ,PBCDim,LinCrit,TorsLinCrit)
     !
     ! Determine value of internal coordinates.
     ! Input coordintes are now in atomic units!
     ! 
     IMPLICIT NONE
     TYPE(INTC) :: IntCs
     INTEGER :: NIntCs,I,J,K,L,I1,I2,I3,I4,NatmsLoc,PBCDim
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,LinCrit,TorsLinCrit
     REAL(DOUBLE),DIMENSION(3,3) :: BoxShapeT,BoxShape,InvBoxSh
     REAL(DOUBLE),DIMENSION(3)   :: Vect1,Vect2
     REAL(DOUBLE),DIMENSION(3,4) :: XYZAux
     !
     IF(IntCs%N<=0) RETURN
     NIntCs=SIZE(IntCs%Def%C)
     NatmsLoc=SIZE(XYZ,2)-3
     DO J=1,3
       BoxShapeT(J,1:3)=XYZ(1:3,NatmsLoc+J)
       BoxShape(1:3,J)=XYZ(1:3,NatmsLoc+J)
     ENDDO
     IF(PBCDim>0) THEN
       InvBoxSh=InverseMatrix(BoxShape)
     ENDIF
     !
     IntCs%Value%D=Zero 
     !
     DO I=1,NIntCs
       IF(.NOT.IntCs%Active%L(I)) THEN
         IntCs%Value%D(I)=Zero 
         CYCLE
       ENDIF
       CALL PBCXYZAux(XYZ,BoxShapeT,XYZAux,IntCs,I)
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         CALL STRE(XYZAux(1:3,1),XYZAux(1:3,2),Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND'.OR. &
               IntCs%Def%C(I)(1:5)=='ALPHA'.OR. &
               IntCs%Def%C(I)(1:4)=='BETA'.OR. &
               IntCs%Def%C(I)(1:5)=='GAMMA') THEN
         CALL BEND(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3), &
                   Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:5)=='LINB1') THEN
         I4=IntCs%Atoms%I(I,4)
         CALL LinB(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),XYZAux(1:3,4),I4,&
           Value1=IntCs%Value%D(I),Value2=IntCs%Value%D(I+1))  
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         CALL TORS(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),&
           XYZAux(1:3,4),TorsLinCrit,IntCs%Active%L(I), &
           Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
         CALL OUTP(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,4),&
           XYZAux(1:3,3),TorsLinCrit,IntCs%Active%L(I),Value_O=IntCs%Value%D(I))
       ELSE IF(IntCs%Def%C(I)(1:6)=='VOLM_L') THEN
         CALL VOLUME(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),&
           XYZAux(1:3,4),IntCs%Active%L(I),Value_O=IntCs%Value%D(I))
       ELSE IF(IntCs%Def%C(I)(1:6)=='AREA_L') THEN
         CALL AREA(XYZAux(1:3,1),XYZAux(1:3,2),XYZAux(1:3,3),&
                   IntCs%Active%L(I),Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='CART') THEN
         Vect1=XYZAux(1:3,1)
         IF(PBCDim>0) THEN
           CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
           Vect1=Vect2
         ENDIF
         IF(IntCs%Def%C(I)(1:5)=='CARTX') THEN
           IntCs%Value%D(I)=Vect1(1)    
         ELSE IF(IntCs%Def%C(I)(1:5)=='CARTY') THEN
           IntCs%Value%D(I)=Vect1(2)    
         ELSE IF(IntCs%Def%C(I)(1:5)=='CARTZ') THEN
           IntCs%Value%D(I)=Vect1(3)    
         ENDIF
       ENDIF
       !
     ENDDO 
   END SUBROUTINE INTCValue
!
!----------------------------------------------------------------
!
   SUBROUTINE PBCXYZAux(XYZ,BoxShapeT,XYZAux,IntCs,I)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,XYZAux,BoxShapeT
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,II,K
     REAL(DOUBLE)                :: XTrans,YTrans,ZTrans
     INTEGER                     :: IA,IB,IC
     !
     XYZAux=Zero
     DO J=1,4
       IF(IntCs%Atoms%I(I,J)==0) EXIT 
       K=(J-1)*3+1
       IA=IntCs%Cells%I(I,K)
       IB=IntCs%Cells%I(I,K+1)
       IC=IntCs%Cells%I(I,K+2)
       XTrans=IA*BoxShapeT(1,1)+&
             IB*BoxShapeT(2,1)+IC*BoxShapeT(3,1)
       YTrans=IA*BoxShapeT(1,2)+&
             IB*BoxShapeT(2,2)+IC*BoxShapeT(3,2)
       ZTrans=IA*BoxShapeT(1,3)+&
             IB*BoxShapeT(2,3)+IC*BoxShapeT(3,3)
       II=IntCs%Atoms%I(I,J)
       XYZAux(1:3,J)=XYZ(1:3,II)+(/XTrans,YTrans,ZTrans/)
     ENDDO
   END SUBROUTINE PBCXYZAux
!
!------------------------------------------------------------
!
   SUBROUTINE CartToInternal(IntCs,VectCart,VectInt,XYZ,PBCDim, &
                             TrfGrd,CtrlCoord,CtrlTrf,Print,SCRPath)
     REAL(DOUBLE),DIMENSION(:) :: VectCart,VectInt
     TYPE(DBL_VECT)  :: VectCartAux,VectIntAux
     TYPE(DBL_VECT)  :: VectCartAux2,VectCartAux3
     REAL(DOUBLE)    :: DiffMax,RMSD
     REAL(DOUBLE)    :: SumU
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)    :: BoxShapeT(3,3),Vect(3),Vect1(3)
     INTEGER         :: NCart,NatmsLoc,I,II,J,K,L,NIntC,PBCDim
     INTEGER         :: Print
     TYPE(INTC)      :: IntCs
     TYPE(Cholesky)  :: CholData
     TYPE(INT_VECT)  :: ISpB,JSpB,IPerm1,IPerm2
     TYPE(DBL_VECT)  :: ASpB
     TYPE(GrdTrf)    :: TrfGrd
     TYPE(CoordCtrl) :: CtrlCoord
     TYPE(TrfCtrl)   :: CtrlTrf
     CHARACTER(LEN=*):: SCRPath
     LOGICAL         :: Print2
     !
     NCart=SIZE(VectCart)
     NatmsLoc=NCart/3-3
     NIntC=SIZE(IntCs%Def%C)
     Print2=(Print>=DEBUG_GEOP_MIN)
     !
     CALL New(VectCartAux,NCart)
     CALL New(VectCartAux2,NCart)
     CALL New(VectCartAux3,NCart)
     CALL New(VectIntAux,NIntC)
     VectCartAux3%D=VectCart
     !
     ! Convert atomic gradients into fractionals
     !
     IF(PBCDim>0) THEN
       DO J=1,3
         BoxShapeT(J,1:3)=XYZ(1:3,NatmsLoc+J)
       ENDDO
       DO I=1,NatmsLoc
         K=3*(I-1)+1
         L=K+2
         Vect=VectCartAux3%D(K:L)
         CALL DGEMM_NNc(3,3,1,One,Zero,BoxShapeT,Vect,Vect1)
         VectCartAux3%D(K:L)=Vect1
       ENDDO
     ENDIF
     !
     ! Get B matrix and Bt*B inverse
     !
     CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
     !
     IF(Print2) THEN
       WRITE(*,111) NIntC 
       WRITE(Out,111) NIntC
       111 FORMAT('Gradient transformation, No. Int. Coords= ',I7)
       IF(.NOT.CtrlTrf%DoClssTrf) THEN
         WRITE(*,112) CtrlTrf%ThreeAt
         WRITE(*,112) CtrlTrf%ThreeAt_2
         WRITE(Out,112) CtrlTrf%ThreeAt
         WRITE(Out,112) CtrlTrf%ThreeAt_2
       ENDIF
       112 FORMAT('Three-atoms reference system used, atoms are ',3I4)
     ENDIF
     !
     ! Cartesian --> Internal transformation
     !
     VectInt=Zero
     DO II=1,TrfGrd%MaxIt_GrdTrf
       !
       VectCartAux%D=Zero
       !
       ! gc-Bt*gi
       !
       CALL CALC_BxVect(ISpB,JSpB,ASpB,VectInt,VectCartAux%D,Trp_O=.TRUE.)
       VectCartAux%D=VectCartAux3%D-VectCartAux%D
       !
       ! GcInv*[gc-Bt*gi]
       !
       CALL CALC_GcInvCartV(CholData,VectCartAux%D,VectCartAux2%D)
       !
       ! B*GcInv*[gc-Bt*gi]
       !
       CALL CALC_BxVect(ISpB,JSpB,ASpB,VectIntAux%D,VectCartAux2%D)
       !
       ! Check convergence
       !
       DiffMax=Zero
       DO I=1,NIntC ; DiffMax=MAX(DiffMax,ABS(VectIntAux%D(I))) ; ENDDO
       RMSD=DOT_PRODUCT(VectIntAux%D,VectIntAux%D)
       RMSD=SQRT(RMSD/DBLE(NIntC))
       !
       ! IF DiffMax is too large, eg. due to the 'bad' quality 
       ! of the preconditioner, rescale gradients
       !
       IF(DiffMax>TrfGrd%MaxGradDiff) THEN
         IF(Print2) THEN
           WRITE(*,*) 'Rescale Step from ',DiffMax,' to ',TrfGrd%MaxGradDiff
           WRITE(Out,*) 'Rescale Step from ',DiffMax,' to ',TrfGrd%MaxGradDiff
         ENDIF
         SumU=TrfGrd%MaxGradDiff/DiffMax
         VectIntAux%D(:)=SumU*VectIntAux%D(:)
         DiffMax=TrfGrd%MaxGradDiff
       ENDIF
       !
       ! gi+B*GcInv*[gc-Bt*gi]
       !
       VectInt=VectInt+VectIntAux%D
       !
       ! Review iteration
       !
       IF(Print2) THEN
         WRITE(*,110) II,DiffMax,RMSD
         WRITE(Out,110) II,DiffMax,RMSD
       ENDIF
       110  FORMAT('Grad Trf, step= ',I3,' MaxChange= ',F12.6,&
                   ' ChangeNorm= ',F12.6)
       !      
       IF(DiffMax<TrfGrd%GrdTrfCrit) EXIT 
     ENDDO
     !
     IF(II>=TrfGrd%MaxIt_GrdTrf) THEN
       IF(Print2) THEN
         WRITE(*,777) 
         WRITE(*,778) 
         WRITE(Out,777) 
         WRITE(Out,778) 
         777 FORMAT('Stop Gradient Trf, max. number of Iterations exceeded!')
         778 FORMAT('Use current gradient vector!')
       ENDIF
     ELSE
       IF(Print2) THEN
         WRITE(*,120) II
         WRITE(Out,120) II
       ENDIF
     ENDIF
     120  FORMAT('Gradient transformation converged in ',I3,' steps')
     !
     ! Tidy up
     !
     CALL Delete(VectIntAux)
     CALL Delete(VectCartAux3)
     CALL Delete(VectCartAux2)
     CALL Delete(VectCartAux)
     CALL DeleteBMatInfo(ISpB,JSpB,ASpB,CholData)
     !
     ! Project out hard constraints
     !
  !  CALL ProjectBCol(SCRPath,IntCs,XYZ,VectInt,PBCDim,Print2)
!CALL PrtIntCoords(IntCs,VectInt,'aft hard constr filt',PBCDim_O=PBCDim)
   END SUBROUTINE CartToInternal
!
!------------------------------------------------------------------
!
   SUBROUTINE InternalToCart(XYZ,AtNum,IntCs,PredVals,RefPoints,Print, &
                           GBackTrf,GTrfCtrl,GCoordCtrl,GConvCr, &
                           GConstr,PBCDim,&
                           SCRPath,PWDPath,IntCsE,MixMat_O,iGEO_O)
     REAL(DOUBLE),DIMENSION(:,:)          :: XYZ
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: MixMat_O
     REAL(DOUBLE),DIMENSION(:)            :: PredVals,RefPoints,AtNum
     INTEGER,OPTIONAL                     :: iGEO_O
     TYPE(DBL_VECT)             :: VectCart
     TYPE(DBL_VECT)             :: VectCartAux,VectIntAux
     TYPE(DBL_VECT)             :: VectCartAux2,VectIntAux2,IntCDispl
     TYPE(DBL_VECT)             :: VectIntReq
     TYPE(DBL_VECT)             :: ValSt
     TYPE(DBL_RNK2)             :: ActCarts,RotCarts
     REAL(DOUBLE)               :: BoxShape(3,3),Aux9(9)
     REAL(DOUBLE)               :: Vect(3),Frac(3),NewFrac(3)
     REAL(DOUBLE)               :: BoxShapeOld(3,3),InvBoxShOld(3,3)
     REAL(DOUBLE)               :: DiffMax,RMSD,RMSDOld
     REAL(DOUBLE)               :: SumU,ConstrMax,ConstrRMS,Fact,Crit
     REAL(DOUBLE)               :: ConstrRMSOld,ConstrMaxCrit,RMSCrit
     REAL(DOUBLE)               :: BackLinCrit,BackTLinCrit
     INTEGER                    :: NCart,I,IStep,J,NT,K,L
     INTEGER                    :: NIntC,NConstr,IRep,RepMax
     INTEGER                    :: NatmsLoc,NCartConstr,PBCDim
     TYPE(INTC)                 :: IntCs,IntCsE
     TYPE(INT_VECT)             :: ISpB,JSpB
     TYPE(DBL_VECT)             :: ASpB
     LOGICAL                    :: RefreshB,RefreshAct
     LOGICAL                    :: DoIterate
     TYPE(Cholesky)             :: CholData
     TYPE(BackTrf)              :: GBackTrf
     TYPE(Constr)               :: GConstr
     TYPE(TrfCtrl)              :: GTrfCtrl
     TYPE(CoordCtrl)            :: GCoordCtrl
     TYPE(GConvCrit)            :: GConvCr
     CHARACTER(LEN=*)           :: SCRPath,PWDPath
     LOGICAL                    :: Print2,DoRepeat
     INTEGER                    :: Print
     !
     BackLinCrit=1.D-8
     BackTLinCrit=1.D-8
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc   
     NIntC=SIZE(IntCs%Def%C)
     NT=SIZE(PredVals)
     Print2=(Print>=DEBUG_GEOP_MAX)
     IF(.NOT.PRESENT(MixMat_O).AND.NT/=NIntC) THEN
       CALL Halt('MixMat missing from the CALL of InternalToCart.')
     ENDIF
     DoRepeat=.FALSE.
     RepMax=5 
     ! 
     ! Auxiliary arrays
     !
     CALL New(ActCarts,(/3,NatmsLoc/))
     CALL New(VectCart,NCart)
     CALL New(VectCartAux,NCart)
     CALL New(VectCartAux2,NCart)
     CALL New(IntCDispl,NIntC)
     CALL New(VectIntAux,NT)
     CALL New(VectIntAux2,NT)
     CALL New(VectIntReq,NT)
     CALL New(ValSt,NT)
     !
     VectIntReq%D=PredVals
     CALL INTCValue(IntCs,XYZ,PBCDim,BackLinCrit,BackTLinCrit)
     CALL SetBackToRefs(IntCs%Value%D,IntCs,RefPoints)
     CALL SetBackToRefs(VectIntReq%D,IntCs,RefPoints)
     PredVals=VectIntReq%D
     !
     IF(PRESENT(MixMat_O)) THEN
       CALL DGEMM_TNc(NT,NIntC,1,One,Zero, &
                      MixMat_O,IntCs%Value%D,ValSt%D)
     ELSE
       ValSt%D=IntCs%Value%D
     ENDIF
     !
     ! Repeat until convergence
     !
     DO IRep=1,RepMax
       RefreshB=.TRUE.
       RefreshAct=.TRUE.
       DoRepeat=.FALSE.
       !
       ! initialization of new Cartesians
       !
       ActCarts%D=XYZ
       CALL CartRNK2ToCartRNK1(VectCart%D,ActCarts%D)
       !
       ! Internal --> Cartesian transformation
       !
       IF(Print2) THEN
         WRITE(*,450) NIntC
         WRITE(Out,450) NIntC
       ENDIF
       450 FORMAT('Iterative back-transformation, No. Int. Coords=',I7)
       !
       ConstrMax=GConstr%ConstrMaxCrit*10.D0
       ConstrRMS=1.D0
       ConstrRMSOld=2.D0
       RMSD=1.D+9
       !
       DO IStep=1,GBackTrf%MaxIt_CooTrf
         IF(GTrfCtrl%PrtBackTr.AND.PRESENT(iGEO_O)) THEN
           CALL PrtBackTrf(AtNum,ActCarts%D,PBCDim,PWDPath, &
                           IRep,IStep,iGEO_O)
         ENDIF
         !
         CALL INTCValue(IntCs,ActCarts%D,PBCDim, &
                        BackLinCrit,BackTLinCrit)
         CALL SetBackToRefs(IntCs%Value%D,IntCs,RefPoints)
         !
         ! Get B and refresh values of internal coords
         !
         IF(RefreshB.AND.RefreshAct) THEN
           CALL RefreshBMatInfo(IntCs,ActCarts%D,GTrfCtrl,GConvCr, &
                               BackLinCrit,BackTLinCrit,PBCDim, &
                               Print,SCRPath,DoCleanCol_O=.TRUE.)
           CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
         ENDIF
         !
         IF(PRESENT(MixMat_O)) THEN
           CALL DGEMM_TNc(NT,NIntC,1,One,Zero, &
                          MixMat_O,IntCs%Value%D,VectIntAux2%D)
         ELSE
           VectIntAux2%D=IntCs%Value%D
         ENDIF
         !
         ! Calculate difference between required and actual internals
         !
         VectIntAux%D=VectIntReq%D-VectIntAux2%D
         IF(PRESENT(MixMat_O)) THEN
           CALL DGEMM_NNc(NIntC,NT,1,One,Zero, &
                          MixMat_O,VectIntAux%D,IntCDispl%D)
         ELSE
           IntCDispl%D=VectIntAux%D
         ENDIF
         !
         CALL MapAngleDispl(IntCs,IntCDispl%D) 
         !
         ! Check convergence on constraints
         !
         IF(GConstr%NConstr/=0) THEN
           ConstrRMSOld=ConstrRMS
           CALL ConstrConv(IntCs,VectIntAux%D,ConstrMax,ConstrRMS)
         ENDIF
         !
         ! Do transformation
         ! 
         ! Bt*[phi_r-phi_a]
         !
         CALL CALC_BxVect(ISpB,JSpB,ASpB,IntCDispl%D, &
                          VectCartAux%D,Trp_O=.TRUE.)
         !
         ! GcInv*Bt*[phi_r-phi_a]
         !
         CALL CALC_GcInvCartV(CholData,VectCartAux%D,VectCartAux2%D)
         !
         ! Update lattice
         !
         IF(PBCDim>0) THEN
           Aux9=VectCart%D(NCart-8:NCart)+VectCartAux2%D(NCart-8:NCart)
           CALL SetFixedLattice(Aux9,IntCsE,GConstr)
           VectCartAux2%D(NCart-8:NCart)=Aux9-VectCart%D(NCart-8:NCart)
           DO J=1,3
             K=3*(J-1)+1
             L=K+2
             BoxShapeOld(1:3,J)=VectCart%D(NCart-9+K:NCart-9+L)
             BoxShape(1:3,J)=Aux9(K:L)
           ENDDO
           InvBoxShOld=InverseMatrix(BoxShapeOld)
           !
           ! Create Cartesian displacements from fractional ones
           !
           DO J=1,NatmsLoc-3
             K=3*(J-1)+1
             L=K+2
             CALL DGEMM_NNc(3,3,1,One,Zero, &
                            InvBoxShOld,VectCart%D(K:L),Frac)
             NewFrac=Frac+VectCartAux2%D(K:L)
             CALL DGEMM_NNc(3,3,1,One,Zero, &
                            BoxShape,NewFrac,Vect)
             VectCartAux2%D(K:L)=Vect-VectCart%D(K:L)
           ENDDO
         ENDIF
         !
         ! Project out translations. 
         ! Rotations are treated in fractional coordinates.
         !
         IF(GTrfCtrl%DoTranslOff) THEN
           CALL TranslsOff(VectCartAux2%D(1:NCart-9),Print2)
         ENDIF
         IF(GTrfCtrl%DoRotOff) THEN
           CALL RotationsOff(VectCartAux2%D,VectCart%D,Print2,PBCDim)
         ENDIF
         !
         ! Scale Cartesian displacements
         ! (preserves constrained fractionals)
         !
         RMSDOld=RMSD
         CALL ScaleDispl(VectCartAux2%D,GBackTrf%MaxCartDiff, &
                         DiffMax,RMSD)
         !
         ! Refresh B matrix?  
         !
         IF(DiffMax>GBackTrf%DistRefresh) THEN
           RefreshAct=.TRUE.
         ELSE
           RefreshAct=.FALSE.
         ENDIF
         !
         ! Modify Cartesians
         !
         VectCart%D=VectCart%D+VectCartAux2%D
         CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
         !
         ! Review iteration
         !
         IF(Print2) THEN
           WRITE(*,210) IStep,DiffMax,RMSD
           WRITE(Out,210) IStep,DiffMax,RMSD
         ENDIF
         210  FORMAT('Step= ',I3,'   Max_DX= ',F12.6,'  X_RMSD= ',F12.6)
         !      
         CALL BackTrfConvg(GConstr,GBackTrf, &
           DoIterate,DiffMax,RMSD,RMSDOld,ConstrMax, &
           ConstrRMS,ConstrRMSOld,IStep,RefreshAct)
         !
         IF(DoIterate.AND.IStep/=GBackTrf%MaxIt_CooTrf) THEN
           IF(RefreshB.AND.RefreshAct) THEN
             CALL DeleteBMatInfo(ISpB,JSpB,ASpB,CholData)
           ENDIF
         ELSE
           CALL DeleteBMatInfo(ISpB,JSpB,ASpB,CholData)
           EXIT
         ENDIF
       ENDDO          
       !
       IF(IStep>=GBackTrf%MaxIt_CooTrf) THEN
         IF(RMSD>0.001D0) THEN
           IF(Print2) THEN
             WRITE(*,*) IRep,' Rescaling and Repeating back-transformation'
             WRITE(Out,*) IRep,' Rescaling and Repeating back-transformation'
           ENDIF
           IF(IRep>RepMax) THEN
             WRITE(*,*)'Warning! Iterative backtransformation has not converged accurately.'
             WRITE(Out,*)'Warning! Iterative backtransformation has not converged accurately.'
             EXIT            
           ENDIF
           !
           VectIntAux%D=VectIntReq%D-ValSt%D  
           CALL MapAngleDispl(IntCs,VectIntAux%D)
           VectIntReq%D=ValSt%D+0.50D0*VectIntAux%D
           DoRepeat=.TRUE.
           CYCLE
         ENDIF
         IF(Print2) THEN
           WRITE(*,180) 
           WRITE(Out,180) 
           WRITE(*,190) 
           WRITE(Out,190) 
         ENDIF
         EXIT 
       ELSE
         IF(IRep<RepMax) THEN
           CALL INTCValue(IntCs,ActCarts%D,PBCDim, &
                          BackLinCrit,BackTLinCrit)
           CALL SetBackToRefs(IntCs%Value%D,IntCs,RefPoints)
           IntCDispl%D=IntCs%Value%D-ValSt%D
           CALL MapAngleDispl(IntCs,IntCDispl%D) 
           !
           ! check for the size of the displacement
           CALL CheckBigStep(IRep,IntCs,DoRepeat,GCoordCtrl%MaxStre, &
                   GCoordCtrl%MaxAngle,PredVals,VectIntReq%D, &
                   ValSt%D,IntCDispl%D,Print2)
           IF(DoRepeat) THEN
             VectIntAux%D=VectIntReq%D-ValSt%D  
             CALL MapAngleDispl(IntCs,VectIntAux%D)
             VectIntReq%D=ValSt%D+0.50D0*VectIntAux%D
             CYCLE
           ENDIF
         ENDIF 
         IF(Print2) THEN
           WRITE(*,220) IStep
           WRITE(Out,220) IStep
         ENDIF
         EXIT
       ENDIF
     ENDDO !!! Repeat
     !
     180  FORMAT('Stop Coord Back-Trf, max. number of Iterations exceeded!')
     190  FORMAT('Use Current Geometry!')
     220  FORMAT('Coordinate back-transformation converged in ',&
                 I3,' steps')
     !
     ! Fill new Cartesians into XYZ  
     !
     XYZ=ActCarts%D
     !
     ! rotate lattice back to standard orientation
     !
     IF(PBCDim>0) THEN
       CALL New(RotCarts,(/3,NatmsLoc+1/))
       DO J=1,NatmsLoc ; RotCarts%D(1:3,J)=XYZ(1:3,J) ; ENDDO
       RotCarts%D(1:3,NatmsLoc+1)=Zero
       GTrfCtrl%ThreeAt(1)=NatmsLoc+1
       GTrfCtrl%ThreeAt(2)=NatmsLoc-3+1
       GTrfCtrl%ThreeAt(3)=NatmsLoc-3+2
       CALL CALC_XYZRot(RotCarts%D,GTrfCtrl%ThreeAt,&
                        .FALSE.,GTrfCtrl%TranslAt1, &
                        GTrfCtrl%RotAt2ToX,GTrfCtrl%RotAt3ToXY, &
                        DoCopy_O=.TRUE.)
       DO J=1,NatmsLoc ; XYZ(1:3,J)=RotCarts%D(1:3,J) ; ENDDO
       XYZ(2:3,NatmsLoc-3+1)=Zero
       XYZ(3,NatmsLoc-3+2)=Zero
       CALL Delete(RotCarts)
     ENDIF
     !
     ! Tidy up
     !
     CALL Delete(IntCDispl)
     CALL Delete(ValSt)
     CALL Delete(VectIntReq)
     CALL Delete(VectIntAux)
     CALL Delete(VectIntAux2)
     CALL Delete(VectCartAux2)
     CALL Delete(VectCartAux)
     CALL Delete(VectCart)
     CALL Delete(ActCarts)
   END SUBROUTINE InternalToCart
!
!---------------------------------------------------------------------
!
   SUBROUTINE ProjectBCol(SCRPath,IntCs,XYZ,VectInt,PBCDim,Print2)
     TYPE(BMATR)                :: B,BS
     TYPE(INTC)                 :: IntCs
     REAL(DOUBLE),DIMENSION(:)  :: VectInt
     INTEGER                    :: NCart,NatmsLoc,PBCDim
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     CHARACTER(LEN=*)           :: SCRPath
     TYPE(INT_VECT)             :: ISpB,JSpB
     TYPE(DBL_VECT)             :: ASpB
     TYPE(Cholesky)             :: CholData
     LOGICAL                    :: Print2
     REAL(DOUBLE)               :: Fact
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B',B_O=B)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     CALL SetEq(BS,B)
     CALL CleanBConstr(IntCs,BS,NatmsLoc)
     CALL CleanBLConstr(XYZ,IntCs,BS,PBCDim)
     CALL CleanBLRot(XYZ,IntCs,BS,PBCDim)
     B%B%D=B%B%D-BS%B%D
     B%BL%D=B%BL%D-BS%BL%D
     ! B now contains only colums of hard constraints
     ! and lattice rotations
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     CALL CholFact(ISpB,JSpB,ASpB,NCart,IntCs%N, &
                   CholData,Print2,Shift_O=1.D-4)
     CALL POffHardGc(VectInt,ISpB,JSpB,ASpB,CholData, &
                     NCart,IntCs%N,Fact)
     IF(Print2) THEN
       WRITE(*,100) Fact
       WRITE(Out,100) Fact
       100 FORMAT('Percentage of Hard Constraints Projected Out= ',F6.2)
     ENDIF
     !
     CALL Delete(B)
     CALL Delete(BS)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(CholData)
   END SUBROUTINE ProjectBCol
!
!---------------------------------------------------------------------
!
   SUBROUTINE SYMMSTRESS(Grad,XYZ,PBCDim)
     REAL(DOUBLE),DIMENSION(9)   :: Grad,Grad1,Grad2,Delta
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(9,9) :: P
     REAL(DOUBLE)                :: Fact,Fact2,Crit
     INTEGER                     :: PBCDim,I,J,MaxStep
     !
     MaxStep=100
     Crit=1.D-6
     CALL GetPBCProj(PBCDim,P,XYZ_O=XYZ)
     Grad1=Grad
     DO I=1,MaxStep
       CALL DGEMM_NNc(9,9,1,One,Zero,P,Grad1,Grad2)
       Delta=Grad1-Grad2
       Grad2(1)=Grad(1)
       Grad2(4:5)=Grad(4:5)
       Grad2(7:9)=Grad(7:9)
       Fact=Zero 
       DO J=1,9 
         Fact2=ABS(Delta(J))
         IF(Fact<Fact2) Fact=Fact2
       ENDDO
       Grad1=Grad2
       IF(Fact<Crit) EXIT
     ENDDO
     IF(Fact>Crit) THEN
       CALL Halt('Unsuccessful symmetrization of stress in SYMMSTRESS')
     ENDIF
     Grad=Grad1
   END SUBROUTINE SYMMSTRESS
!
!---------------------------------------------------------------------
!
   SUBROUTINE CheckActive(Active,ReqInt)
     LOGICAL,DIMENSION(:)      :: Active
     REAL(DOUBLE),DIMENSION(:) :: ReqInt
     INTEGER                   :: I,J,NIntC
     !
     NIntC=SIZE(Active)
     DO I=1,NIntC
       IF(.NOT.Active(I)) ReqInt(I)=Zero
     ENDDO
   END SUBROUTINE CheckActive
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanConstrIntc(IntCDispl,XYZ,IntCsE,SCRPath,&
                              GTrfCtrl,GCoordCtrl,GConvCr,PBCDim,Print)
     TYPE(TrfCtrl)                :: GTrfCtrl
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(GConvCrit)              :: GConvCr
     INTEGER                      :: PBCDim,Print,NCart,NatmsLoc
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     REAL(DOUBLE),DIMENSION(:)    :: IntCDispl
     TYPE(INTC)                   :: IntCsE,IntCsAux
     CHARACTER(LEN=*)             :: SCRPAth
     TYPE(INT_VECT)               :: ISpB,JSpB
     TYPE(DBL_VECT)               :: ASpB
     TYPE(CHOLESKY)               :: CholData
     TYPE(DBL_VECT)               :: IntAux,CartAux
     REAL(DOUBLE)                 :: Fact
     INTEGER                      :: I,J
     !
     IF(IntCsE%N==0) RETURN
     IF(SIZE(IntCDispl)/=IntCsE%N) CALL Halt('Dimension error in CleanConstrIntc')
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     !
     ! set all constraints active
     !
     CALL New(IntCsAux,IntCsE%N)
     CALL SetEq(IntCsE,IntCsAux,1,IntCsE%N,1)
     DO I=1,IntCsAux%N
       IF(IntCsAux%Constraint%L(I)) IntCsAux%Active%L(I)=.TRUE.
     ENDDO
     !
     CALL New(IntAux,IntCsAux%N)
     CALL New(CartAux,NCart)
     !
     CALL RefreshBMatInfo(IntCsAux,XYZ,GTrfCtrl,GConvCr, &
                          GCoordCtrl%LinCrit, &
                          GCoordCtrl%TorsLinCrit, &
                          PBCDim,Print,SCRPath)
     CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)    
     !
     CALL CALC_BxVect(ISpB,JSpB,ASpB,IntCDispl,CartAux%D,Trp_O=.TRUE.)
     CALL GcInvIter(CartAux%D,ISpB,JSpB,ASpB,CholData,NCart)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,IntAux%D,CartAux%D)
     Fact=DOT_PRODUCT(IntAux%D,IntCDispl)/ &
          (DOT_PRODUCT(IntCDispl,IntCDispl)+1.D-10)*100.D0
     WRITE(*,100) Fact
     IntCDispl=IntAux%D
     100  FORMAT('Constraints projected out: ',F7.3,'%')
     !
     CALL Delete(IntAux)
     CALL Delete(CartAux)
   END SUBROUTINE CleanConstrIntc
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanBoxPars(Vec,PBCDim)
     REAL(DOUBLE),DIMENSION(6) :: Vec
     INTEGER                   :: PBCDim
     IF(PBCDim==0) RETURN
     IF(PBCDim<3) THEN
       Vec(3)=AngstromsToAu
       Vec(6)=Half*PI
     ENDIF
     IF(PBCDim<2) THEN
       Vec(2)=AngstromsToAu
       Vec(4)=Half*PI
       Vec(5)=Half*PI
     ENDIF
   END SUBROUTINE CleanBoxPars
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrtBackTrf(AtNumIn,XYZ,PBCDim,PWDPath,IRep,IStep,iGEO)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNumIn
     CHARACTER(LEN=*)            :: PWDPath
     CHARACTER(LEN=DCL)          :: Title
     INTEGER                     :: I,J,NatmsLoc,IRep,IStep,PBCDim,iGEO
     TYPE(DBL_RNK2)              :: XYZRepl
     TYPE(INT_VECT)              :: IEq,AtNum,AtNumRepl
     TYPE(INT_RNK2)              :: Cells
     TYPE(INT_RNK3)              :: CellEq
     LOGICAL                     :: PrtLVect
     CHARACTER(LEN=1)            :: Dir
     !
     PrtLVect=.TRUE.
     NatmsLoc=SIZE(XYZ,2)
     Title='Geom= '//TRIM(IntToChar(iGEO))//'IRep= '//TRIM(IntToChar(IRep))//' IStep= '//TRIM(IntToChar(IStep))
     !
     CALL New(AtNum,NatmsLoc)
     DO I=1,NatmsLoc
       AtNum%I(I)=INT(AtNumIn(I))
     ENDDO
     !
     IF(PBCDim>0) THEN
       DO I=1,4
         IF(I==1) Dir='A'
         IF(I==2) Dir='B'
         IF(I==3) Dir='C'
         IF(I==4) Dir='T'
         CALL PrepCells(XYZ,AtNum%I,PBCDim,XYZRepl, &
                        AtNumRepl,Cells,CellEq,IEq,Dir_O=Dir,NDim_O=-1)
         IF(PrtLVect) THEN
           CALL PrtXYZ(AtNumRepl%I,XYZRepl%D,TRIM(PWDPath)//Dir//'Back.xyz',&
                       Title,XYZL_O=XYZ)
         ELSE 
           CALL PrtXYZ(AtNumRepl%I,XYZRepl%D,TRIM(PWDPath)//'Back.xyz',&
                       Title)
         ENDIF
         CALL Delete(XYZRepl) 
         CALL Delete(AtNumRepl) 
         CALL Delete(Cells) 
         CALL Delete(CellEq) 
         CALL Delete(IEq) 
       ENDDO
     ELSE
       CALL PrepCells(XYZ,AtNum%I,0,XYZRepl,AtNumRepl, &
                      Cells,CellEq,IEq)
       CALL PrtXYZ(AtNumRepl%I,XYZRepl%D,TRIM(PWDPath)//'Back.xyz',&
                   Title)
       CALL Delete(XYZRepl) 
       CALL Delete(AtNumRepl) 
       CALL Delete(Cells) 
       CALL Delete(CellEq) 
       CALL Delete(IEq) 
     ENDIF
     !
     CALL Delete(AtNum) 
   END SUBROUTINE PrtBackTrf
!
!---------------------------------------------------------------------
!
   SUBROUTINE CheckBigStep(IRep,IntCs,DoRepeat,MaxStre,MaxAngle,PredVals, &
                           VectIntReq,ValSt,IntCDispl,Print2)
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: VectIntReq,ValSt,IntCDispl,PredVals
     REAL(DOUBLE)              :: MaxStre,MaxAngle,Crit,Conv
     REAL(DOUBLE)              :: MaxConv,MaxDispl,Fact,Rigid
     LOGICAL                   :: DoRepeat,Print2
     INTEGER                   :: I,IRep,IMax
     TYPE(DBL_VECT)            :: DReq,DReqAct
     !
     Fact=1.01D0
     IMax=1
     MaxDispl=IntCDispl(1)
     DoRepeat=.FALSE.
     DO I=1,IntCs%N
       IF(.NOT.IntCs%Active%L(I)) CYCLE
       IF(IntCs%Def%C(I)(1:4)=='CART') CYCLE
       !
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         Crit=Fact*MaxStre
       ELSE IF(IntCs%Def%C(I)(1:6)=='VOLM_L') THEN
         Crit=1.D99        
       ELSE IF(IntCs%Def%C(I)(1:6)=='AREA_L') THEN
         Crit=1.D99        
       ELSE
         Crit=Fact*MaxAngle
       ENDIF
       IF(ABS(IntCDispl(I))>Crit) THEN
         IMax=I
         MaxDispl=IntCDispl(I)
         EXIT
       ENDIF
       IF(ABS(IntCDispl(I))>ABS(MaxDispl)) THEN
         IMax=I
         MaxDispl=IntCDispl(I)
       ENDIF
     ENDDO
     !
     IF(IntCs%Def%C(IMax)(1:4)=='STRE') THEN
       Crit=Fact*MaxStre
     ELSE IF(IntCs%Def%C(IMax)(1:6)=='VOLM_L') THEN
       Crit=1.D99        
     ELSE IF(IntCs%Def%C(IMax)(1:6)=='AREA_L') THEN
       Crit=1.D99        
     ELSE 
       Crit=Fact*MaxAngle
     ENDIF
     IF(IntCs%Def%C(IMax)(1:4)=='STRE') THEN
       MaxConv=One/AngstromsToAu
     ELSE IF(IntCs%Def%C(IMax)(1:6)=='VOLM_L') THEN
       MaxConv=One/AngstromsToAu**3
     ELSE IF(IntCs%Def%C(IMax)(1:6)=='AREA_L') THEN
       MaxConv=One/AngstromsToAu**2
     ELSE
       MaxConv=180.D0/PI
     ENDIF
     !
     ! Rigidity
     !
     CALL New(DReq,IntCs%N) 
     CALL New(DReqAct,IntCs%N) 
    !DReq%D=VectIntReq-ValSt
     DReq%D=PredVals-ValSt
     DReqAct%D=VectIntReq-IntCs%Value%D
     CALL MapAngleDispl(IntCs,DReq%D) 
     CALL MapAngleDispl(IntCs,DReqAct%D) 
     DO I=1,IntCs%N
       IF(.NOT.IntCs%Active%L(I)) THEN
         DReq%D(I)=Zero
         DReqAct%D(I)=Zero
       ENDIF
     ENDDO
     Rigid=DOT_PRODUCT(DReqAct%D,DReqAct%D)/DOT_PRODUCT(DReq%D,DReq%D)
     CALL Delete(DReq) 
     CALL Delete(DReqAct) 
     !
     IF(ABS(MaxDispl)>Crit) DoRepeat=.TRUE.
   ! IF(Rigid>0.20D0) DoRepeat=.TRUE.
     IF(Print2) THEN
       WRITE(*,*) 'Rigidity= ',Rigid
       WRITE(Out,*) 'Rigidity= ',Rigid
       IF(DoRepeat) THEN
         WRITE(*,*) IRep,' Repeat from CheckBigStep MaxDispl= ', &
                    MaxConv*MaxDispl,' on ',IMax,IntCs%Def%C(IMax)(1:8)
         WRITE(Out,*) IRep,' Repeat from CheckBigStep MaxDispl= ', &
                    MaxConv*MaxDispl,' on ',IMax,IntCs%Def%C(IMax)(1:8)
       ELSE 
         WRITE(*,*) IRep,'Maximum Displacement from Backtransform.= ', &
                    MaxDispl*MaxConv,' on ',IMax,IntCs%Def%C(IMax)(1:8)
         WRITE(Out,*) IRep,'Maximum Displacement from Backtransform.= ', &
                    MaxDispl*MaxConv,' on ',IMax,IntCs%Def%C(IMax)(1:8)
       ENDIF
     ENDIF
   END SUBROUTINE CheckBigStep
!
!-----------------------------------------------------------------------
!
   SUBROUTINE CheckRedundancy(IntCs,IntCDispl,VectIntReq, &
                              ValSt,Value,XYZ,ActCarts)
     TYPE(INTC)                 :: IntCs
     REAL(DOUBLE),DIMENSION(:)  :: IntCDispl,VectIntReq,ValSt,Value
     REAL(DOUBLE),DIMENSION(:,:):: XYZ,ActCarts
     INTEGER                    :: I,J,NIntC
     REAL(DOUBLE)               :: Fact,RNorm,FNorm,X,Ratio,Crit
     TYPE(DBL_VECT)             :: VectR,VectD
     !
     Crit=0.1D0
     CALL New(VectR,IntCs%N)
     CALL New(VectD,IntCs%N)
     VectR%D=VectIntReq-ValSt ! required displ
     VectD%D=VectR%D-IntCDispl ! displ. diff  
     !
     Ratio=Zero
     DO I=1,IntCs%N
       IF(ABS(IntCDispl(I))>1.D-4) THEN
         X=ABS(VectR%D(I)/IntCDispl(I))
         X=MAX(X,One/X)
         Ratio=MAX(Ratio,X)
       ENDIF
     ENDDO
     Fact=MIN(Ratio,1.2D0)/Ratio
     Fact=MAX(0.5D0,Fact)
     ActCarts=XYZ+Fact*(ActCarts-XYZ)
     CALL Delete(VectD)
     CALL Delete(VectR)
   END SUBROUTINE CheckRedundancy
!
!----------------------------------------------------------
!
   SUBROUTINE SetConstraints(IntCs,PredValue)
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: PredValue
     INTEGER                   :: I,J
     !
     DO I=1,IntCs%N
       IF(IntCs%Constraint%L(I)) THEN
         PredValue(I)=IntCs%ConstrValue%D(I)
       ENDIF
     ENDDO
   END SUBROUTINE SetConstraints
!
!----------------------------------------------------------
!
   SUBROUTINE PrtIntCoords(IntCs,Value,CharU,PBCDim_O)
     !
     TYPE(INTC)       :: IntCs
     INTEGER          :: I,NIntC,J
     REAL(DOUBLE),DIMENSION(:) :: Value
     REAL(DOUBLE)     :: SumU,SumConstr,Conv,ConvC
     CHARACTER(LEN=*) :: CharU   
     INTEGER,OPTIONAL :: PBCDim_O
     LOGICAL          :: DoPrtCells
     !
     Conv=180.D0/PI
     ConvC=One/AngstromsToAu
     NIntC=SIZE(IntCs%Def%C)
     DoPrtCells=.FALSE.
     IF(PRESENT(PBCDim_O)) DoPrtCells=(PBCDim_O>0)
     !
     WRITE(*,*) TRIM(CharU)
     WRITE(Out,*) TRIM(CharU)
     WRITE(*,123) 
     WRITE(*,124) 
     WRITE(Out,123) 
     WRITE(Out,124) 
     123  FORMAT('INTERNAL COORDINATES')
     124  FORMAT('       DEFINITION       ATOMS_INVOLVED      VALUE        CONSTRAINT    ACTIVE')
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         SumU=Value(I)*ConvC
         SUMConstr=IntCs%ConstrValue%D(I)/AngstromsToAu
       ELSE IF(HasAngle(IntCs%Def%C(I))) THEN
         SumU=Value(I)*Conv
         SUMConstr=IntCs%ConstrValue%D(I)*Conv
       ELSE IF(IntCs%Def%C(I)(1:6)=='VOLM_L') THEN 
         SumU=Value(I)*ConvC**3
         SUMConstr=IntCs%ConstrValue%D(I)*ConvC**3
       ELSE IF(IntCs%Def%C(I)(1:6)=='AREA_L') THEN 
         SumU=Value(I)*ConvC**2
         SUMConstr=IntCs%ConstrValue%D(I)*ConvC**2
       ELSE IF(IntCs%Def%C(I)(1:4)=='CART') THEN 
         IF(DoPrtCells) THEN
           SumU=Value(I)
           SUMConstr=IntCs%ConstrValue%D(I)
         ELSE
           SumU=Value(I)*ConvC
           SUMConstr=IntCs%ConstrValue%D(I)*ConvC
         ENDIF
       ENDIF
       WRITE(*,111) I,IntCs%Def%C(I)(1:8),IntCs%Atoms%I(I,1:4),SumU, &
         IntCs%Constraint%L(I),SumConstr,IntCs%Active%L(I)
       WRITE(Out,111) I,IntCs%Def%C(I)(1:8),IntCs%Atoms%I(I,1:4),SumU, &
         IntCs%Constraint%L(I),SumConstr,IntCs%Active%L(I)
       IF(DoPrtCells) THEN
         IF(ANY(IntCs%Cells%I(I,1:12)/=0)) THEN
           WRITE(*,112) IntCs%Cells%I(I,1:12)
           WRITE(Out,112) IntCs%Cells%I(I,1:12)
         ENDIF
       ENDIF
     ENDDO
     !      
     111 FORMAT(I7,2X,A8,2X,4I5,2X,F12.6,L5,F12.6,L5)
     112 FORMAT(15X,4(3I2,2X))
     222 FORMAT(I7,2X,A8,2X,4I5,2X,3F12.6,L5,F12.6,L5)
     !
   END SUBROUTINE PrtIntCoords
!
!----------------------------------------------------------
!
   SUBROUTINE PrtPBCs(XYZ)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(6)   :: Vec
     REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
     INTEGER                     :: I,J,NatmsLoc
     !
     NatmsLoc=SIZE(XYZ,2)-3
     DO I=1,3
       BoxShape(1:3,I)=XYZ(1:3,NatmsLoc+I)
     ENDDO
     CALL CalcBoxPars(Vec,BoxShape)
     !
     WRITE(*,111) 'STRE_A  ',Vec(1)/AngstromsToAu
     WRITE(*,111) 'STRE_B  ',Vec(2)/AngstromsToAu
     WRITE(*,111) 'STRE_C  ',Vec(3)/AngstromsToAu
     WRITE(*,111) 'ALPHA   ',Vec(4)*180.D0/PI
     WRITE(*,111) 'BETA    ',Vec(5)*180.D0/PI
     WRITE(*,111) 'GAMMA   ',Vec(6)*180.D0/PI
     !
     WRITE(Out,111) 'STRE_A  ',Vec(1)/AngstromsToAu
     WRITE(Out,111) 'STRE_B  ',Vec(2)/AngstromsToAu
     WRITE(Out,111) 'STRE_C  ',Vec(3)/AngstromsToAu
     WRITE(Out,111) 'ALPHA   ',Vec(4)*180.D0/PI
     WRITE(Out,111) 'BETA    ',Vec(5)*180.D0/PI
     WRITE(Out,111) 'GAMMA   ',Vec(6)*180.D0/PI
     !
     111 FORMAT('LATTICE  ',A8,2X,20X,2X,F12.6,L5,F12.6,L5)
   END SUBROUTINE PrtPBCs
!
!----------------------------------------------------------
!
   SUBROUTINE RotationsOff(DCarts,Carts,Print,PBCDim)
     REAL(DOUBLE),DIMENSION(:) :: DCarts,Carts
     INTEGER                   :: PBCDim,NCart
     LOGICAL                   :: Print
     !
     NCart=SIZE(DCarts)-9
     IF(PBCDim==0) THEN
       CALL MolRotOff(DCarts(1:NCart),Carts(1:NCart),Print)
     ELSE
       DCarts(NCart+2)=Zero
       DCarts(NCart+3)=Zero
       DCarts(NCart+6)=Zero
      !CALL PBCRotOff(DCarts(NCart+1:NCart+9),Carts(NCart+1:NCart+9), &
      !               Print,PBCDim)
       IF(PBCDim==1) CALL MolRotOff(DCarts(1:NCart),Carts(1:NCart),Print,X_O=.TRUE.)
       IF(PBCDim==2) CALL MolRotOff(DCarts(1:NCart),Carts(1:NCart),Print,Z_O=.TRUE.)
     ENDIF
   END SUBROUTINE RotationsOff
!
!----------------------------------------------------------
!
   SUBROUTINE PBCRotOff(DCarts,Carts,Print,PBCDim)
     REAL(DOUBLE),DIMENSION(:)   :: Carts
     REAL(DOUBLE),DIMENSION(9)   :: DCarts,DCarts2
     INTEGER                     :: PBCDim,NCart,I,J,K,L
     LOGICAL                     :: Print
     REAL(DOUBLE),DIMENSION(9,9) :: P
     REAL(DOUBLE)                :: Fact,Norm
     !
     ! In the standard version A is along X, B in XY, C general
     ! thus projection of rotation is very simple
     !
     DCarts(2:3)=Zero
     DCarts(6)=Zero
     RETURN
     !
     !The lines below refer to a situation, where 
     ! all coordinates of the lattice vectors are allowed to change
     ! while still maintaining a no-rotation framework 
     !
     ! only counter space of rotations is in P, 
     ! space of constraints is not
     CALL GetPBCProj(PBCDim,P,Carts_O=Carts)
     CALL DGEMM_NNc(9,9,1,One,Zero,P,DCarts,DCarts2)
     !
     Norm=DOT_PRODUCT(DCarts,DCarts)
     Fact=DOT_PRODUCT(DCarts2,DCarts2)
     IF(Norm>1.D-6) THEN
       Fact=(Norm-Fact)/Norm*100.D0
     ELSE
       Fact=Zero
     ENDIF
     IF(Print) THEN
       WRITE(*,100) Fact
       WRITE(Out,100) Fact
     ENDIF
     !
     DCarts=DCarts2 
100  FORMAT('PRot= ',F7.3,'%')
   END SUBROUTINE PBCRotOff
!
!----------------------------------------------------------
!
   SUBROUTINE GetPBCProj(PBCDim,P,Carts_O,XYZ_O)
     INTEGER                     :: PBCDim,NCart,I,J,K,L
     INTEGER                     :: NCoinc,Info,DimU,NatmsLoc
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(3,4) :: XYZ
     REAL(DOUBLE),DIMENSION(9,9) :: S,P,P1
     REAL(DOUBLE),DIMENSION(9)   :: Eigs
     TYPE(BMATR)                 :: B
     REAL(DOUBLE),DIMENSION(:,:), OPTIONAL :: XYZ_O
     REAL(DOUBLE),DIMENSION(:)  , OPTIONAL :: Carts_O
     !
     ! This subroutine prepares the projection matrix
     ! that refers to the counterspace of the lattice rotations
     ! in Cartesian representation
     !
     XYZ(1:3,1)=Zero
     IF(PRESENT(Carts_O)) THEN
       NCart=SIZE(Carts_O)-9
       DO I=1,3 
         K=NCart+3*(I-1)+1
         L=K+2
         XYZ(1:3,I+1)=Carts_O(K:L)
       ENDDO 
     ELSE IF(PRESENT(XYZ_O)) THEN
       NatmsLoc=SIZE(XYZ_O,2)-3
       IF(NatmsLoc<0) Call Halt('Dimension error in GetPBCProj')
       DO I=1,3 
         XYZ(1:3,I+1)=XYZ_O(1:3,NatmsLoc+I)
       ENDDO 
     ELSE
       CALL Halt('XYZ or Carts missing from GetPBCProj') 
     ENDIF
     !
     CALL LatticeINTC(IntCs,PBCDim)
     CALL BMatrix(XYZ,IntCs,B,PBCDim,1.D0,1.D0)
     CALL DGEMM_TNc(9,IntCs%N,9,One,Zero,B%BL%D,B%BL%D,S)
     !
     CALL SetDSYEVWork(9)
       BLKVECT%D=S
       CALL DSYEV('V','U',9,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in RotationsOff. INFO='&
                  //TRIM(IntToChar(INFO)))
       S=BLKVECT%D
       DO I=1,9 ; Eigs(I)=BLKVALS%D(I) ; ENDDO
     CALL UnSetDSYEVWork()
     !
     ! Dim is the number of lattice parameters
     !
     IF(PBCDim==1) DimU=1
     IF(PBCDim==2) DimU=3
     IF(PBCDim==3) DimU=6
     IF(DimU<0) CALL Halt('DimU<0 in GetPBCProj.') 
     K=0
     DO I=1,9 
       IF(Eigs(I)>1.D-6) K=K+1
     ENDDO
     IF(K<DimU) CALL Halt('K<DimU in GetPBCProj.')
     Eigs(1:9-DimU)=Zero
     Eigs(9-DimU+1:9)=One
     DO I=1,9 ; P1(1:9,I)=Eigs(I)*S(1:9,I) ; ENDDO
     CALL DGEMM_NTc(9,9,9,One,Zero,P1,S,P)
     !
     CALL Delete(IntCs)
     CALL Delete(B)
   END SUBROUTINE GetPBCProj
!
!----------------------------------------------------------
!
   SUBROUTINE MolRotOff(CartVect,Carts,Print,X_O,Y_O,Z_O)
     REAL(DOUBLE),DIMENSION(:) :: CartVect,Carts
     REAL(DOUBLE)                :: X,Y,Z,XX,YY,ZZ,XY,YZ,ZX
     REAL(DOUBLE)                :: CMX,CMY,CMZ
     REAL(DOUBLE)                :: SumU,SUM1,SUM2,SUM3
     REAL(DOUBLE)                :: V1X,V1Y,V1Z
     REAL(DOUBLE)                :: V2X,V2Y,V2Z
     REAL(DOUBLE)                :: V3X,V3Y,V3Z
     TYPE(DBL_RNK2)              :: Theta,Theta2,CMCarts
     TYPE(DBL_Vect)              :: Rot1,Rot2,Rot3,Vect
     INTEGER                     :: NCart,NatmsLoc,I,J,INFO
     LOGICAL                     :: Print,XOff,YOff,ZOff
     LOGICAL,OPTIONAL            :: X_O,Y_O,Z_O
     !
     XOff=.FALSE.
     YOff=.FALSE.
     ZOff=.FALSE.
     IF(PRESENT(X_O)) XOff=X_O
     IF(PRESENT(Y_O)) YOff=Y_O
     IF(PRESENT(Z_O)) ZOff=Z_O
     !
     NCart=SIZE(CartVect)
     NatmsLoc=NCart/3
     CALL New(Theta,(/3,3/))
     CALL New(Theta2,(/3,3/))
     CALL New(CMCarts,(/3,NatmsLoc/))
     CALL New(Rot1,NCart)
     CALL New(Rot2,NCart)
     CALL New(Rot3,NCart)
     CALL New(Vect,3)
     !
     ! Calculate center of Mass (Masses are unit here)
     !
     CALL CenterOfMass(CMX,CMY,CMZ,Vect_O=Carts)
     !
     ! Mass centered Cartesians
     !
     DO I=1,NatmsLoc
       J=3*(I-1)
       CMCarts%D(1,I)=Carts(J+1)-CMX
       CMCarts%D(2,I)=Carts(J+2)-CMY
       CMCarts%D(3,I)=Carts(J+3)-CMZ
     ENDDO
     !
     ! Build inertial momentum tensor
     !
     IF(.NOT.(XOff.OR.YOff.OR.ZOff)) THEN
       CALL PrincMomTens(CMCarts%D,Theta%D)
       !
       ! Get eigenvectors of inertial momentum tensor
       !
       CALL SetDSYEVWork(3)
         BLKVECT%D(1:3,1:3)=Theta%D(1:3,1:3)
         CALL DSYEV('V','U',3,BLKVECT%D,BIGBLOK,BLKVALS%D, &
         BLKWORK%D,BLKLWORK,INFO)
         IF(INFO/=SUCCEED) &
         CALL Halt('DSYEV hosed in RotationsOff. INFO='&
                    //TRIM(IntToChar(INFO)))
         Theta2%D(1:3,1:3)=BLKVECT%D(1:3,1:3)
       CALL UnSetDSYEVWork()
     ELSE
       Theta2%D=Zero
       IF(XOff) Theta2%D(1,1)=One
       IF(YOff) Theta2%D(2,2)=One
       IF(ZOff) Theta2%D(3,3)=One
     ENDIF
     !
     ! Calculate displacements resulting from rotations
     ! around principal axis by a unit rotation vector pointing along
     ! principal axis. 
     !
     DO I=1,NatmsLoc
       J=(I-1)*3
       ! r1
       CALL CROSS_PRODUCT(Theta2%D(:,1),CMCarts%D(:,I),Rot1%D(J+1:J+3))
       ! r2
       CALL CROSS_PRODUCT(Theta2%D(:,2),CMCarts%D(:,I),Rot2%D(J+1:J+3))
       ! r3
       CALL CROSS_PRODUCT(Theta2%D(:,3),CMCarts%D(:,I),Rot3%D(J+1:J+3))
     ENDDO
     !
     ! Normalize Rot vectors!
     !
     SumU=SQRT(DOT_PRODUCT(Rot1%D,Rot1%D))
     IF(SumU>1.D-6) THEN
       SumU=One/SumU
       Rot1%D=SumU*Rot1%D
     ENDIF
     SumU=SQRT(DOT_PRODUCT(Rot2%D,Rot2%D))
     IF(SumU>1.D-6) THEN
       SumU=One/SumU
       Rot2%D=SumU*Rot2%D
     ENDIF
     SumU=SQRT(DOT_PRODUCT(Rot3%D,Rot3%D))
     IF(SumU>1.D-6) THEN
       SumU=One/SumU
       Rot3%D=SumU*Rot3%D
     ENDIF
     !
     ! test orthogonalities
     !
     !      CALL TranslsOff(Rot1%D,.TRUE.)
     !      CALL TranslsOff(Rot2%D,.TRUE.)
     !      CALL TranslsOff(Rot3%D,.TRUE.)
     !      SumU=DOT_PRODUCT(Rot1%D,Rot1%D)
     !      write(*,*) 'rot test 1 1 ',SumU
     !      SumU=DOT_PRODUCT(Rot1%D,Rot2%D)
     !      write(*,*) 'rot test 1 2 ',SumU
     !      SumU=DOT_PRODUCT(Rot1%D,Rot3%D)
     !      write(*,*) 'rot test 1 3 ',SumU
     !      SumU=DOT_PRODUCT(Rot2%D,Rot2%D)
     !      write(*,*) 'rot test 2 2 ',SumU
     !      SumU=DOT_PRODUCT(Rot2%D,Rot3%D)
     !      write(*,*) 'rot test 2 3 ',SumU
     !      SumU=DOT_PRODUCT(Rot3%D,Rot3%D)
     !      write(*,*) 'rot test 3 3 ',SumU
     !
     ! Now, project out rotations from Cartesian displacement vector
     !
     SumU=DOT_PRODUCT(CartVect(1:NCart),CartVect(1:NCart))
     SUM1=DOT_PRODUCT(Rot1%D(1:NCart),CartVect(1:NCart))
     SUM2=DOT_PRODUCT(Rot2%D(1:NCart),CartVect(1:NCart))
     SUM3=DOT_PRODUCT(Rot3%D(1:NCart),CartVect(1:NCart))
     CartVect(1:NCart)=CartVect(1:NCart)-SUM1*Rot1%D(1:NCart)
     CartVect(1:NCart)=CartVect(1:NCart)-SUM2*Rot2%D(1:NCart)
     CartVect(1:NCart)=CartVect(1:NCart)-SUM3*Rot3%D(1:NCart)
     !
     ! Percentage of rotations
     !
     IF(SumU>1.D-6) THEN
       SUM1=SUM1*SUM1/SumU*100.D0
       SUM2=SUM2*SUM2/SumU*100.D0
       SUM3=SUM3*SUM3/SumU*100.D0
     ENDIF
     IF(Print) THEN
       WRITE(*,100) SUM1,SUM2,SUM3
     ENDIF
     100  FORMAT('Rot1= ',F7.3,'%    Rot2= ', &
                  F7.3,'%     Rot3= ',F7.3,'% ')
     CALL Delete(Vect)
     CALL Delete(Rot3)
     CALL Delete(Rot2)
     CALL Delete(Rot1)
     CALL Delete(CMCarts)
     CALL Delete(Theta2)
     CALL Delete(Theta)
   END SUBROUTINE MolRotOff
!
!----------------------------------------------------------
!
   SUBROUTINE TranslsOff(CartVect,Print)
     REAL(DOUBLE),DIMENSION(:) :: CartVect
     REAL(DOUBLE)              :: SumU,SUM1,SUM2,SUM3
     TYPE(DBL_VECT)            :: Tr1,Tr2,Tr3
     INTEGER                   :: I,J,NCart,NatmsLoc
     LOGICAL                   :: Print
     !
     NCart=SIZE(CartVect)
     NatmsLoc=NCart/3
     CALL New(Tr1,NCart) 
     CALL New(Tr2,NCart) 
     CALL New(Tr3,NCart) 
     Tr1%D=Zero
     Tr2%D=Zero
     Tr3%D=Zero
     SumU=One/SQRT(DBLE(NatmsLoc))
     DO I=1,NatmsLoc
       J=(I-1)*3
       Tr1%D(J+1)=SumU
       Tr2%D(J+2)=SumU
       Tr3%D(J+3)=SumU
     ENDDO
     !
     ! Now, project out translations from CartVect
     !
     SumU =DOT_PRODUCT(CartVect,CartVect)
     IF(SQRT(SumU) > 1.D-12) THEN
       SUM1=DOT_PRODUCT(Tr1%D,CartVect)
         CartVect=CartVect-SUM1*Tr1%D
       SUM2=DOT_PRODUCT(Tr2%D,CartVect)
         CartVect=CartVect-SUM2*Tr2%D
       SUM3=DOT_PRODUCT(Tr3%D,CartVect)
         CartVect=CartVect-SUM3*Tr3%D
       SUM1=SUM1*SUM1/SumU*100.D0
       SUM2=SUM2*SUM2/SumU*100.D0
       SUM3=SUM3*SUM3/SumU*100.D0
     ELSE
       SUM1=Zero
       SUM2=Zero
       SUM3=Zero
     ENDIF
     !
     ! Percentage of translations
     !
     IF(Print) THEN
       WRITE(*,100) SUM1,SUM2,SUM3
     ENDIF
     100  FORMAT(' Tr1= ',F7.3,'%     Tr2= ',F7.3,'%      Tr3= ',&
                  F7.3,'% ')
     !
     CALL Delete(Tr3)
     CALL Delete(Tr2)
     CALL Delete(Tr1)
     !
   END SUBROUTINE TranslsOff
!
!-------------------------------------------------------
!
   SUBROUTINE Rotate(V1,V2,Rot)
     !
     ! Compute matrix Rot, which rotates vector V2 into vector V1
     !
     REAL(DOUBLE),DIMENSION(1:3)      :: V1,V2
     REAL(DOUBLE),DIMENSION(1:3,1:3)  :: Rot  
     REAL(DOUBLE)                 :: SumU,SumM,Sum1,Sum2,Sum3,V1N,V2N
     REAL(DOUBLE) :: PnX,PnY,PnZ
     REAL(DOUBLE) :: CosPhi,SinPhi
     INTEGER      :: I,J,Step,III
     TYPE(DBL_VECT) :: Vect,CrossProd
     !          
     CALL New(Vect,3)
     CALL New(CrossProd,3)
     !
     V1N=DSQRT(DOT_PRODUCT(V1,V1))
     V2N=DSQRT(DOT_PRODUCT(V2,V2))
     V1=V1/V1N
     V2=V2/V2N
     Step=1
     !
     CALL CROSS_PRODUCT(V1,V2,CrossProd%D)
     !
     Sum1=DOT_PRODUCT(CrossProd%D(1:3),CrossProd%D(1:3))
     !
     IF(SQRT(Sum1) < 1.D-12) THEN  !!! V1 & V2 are parallel
       Rot=Zero
       SumU=DOT_PRODUCT(V1,V2)
       DO I=1,3 ; Rot(I,I)=SumU ; ENDDO
     ELSE      
       SumU=One/SQRT(Sum1)
       CrossProd%D=SumU*CrossProd%D
       CosPhi=DOT_PRODUCT(V1,V2)     
       SinPhi=SQRT(SUM1)
       SumU=One-CosPhi
       DO III=1,2    
         Rot(1,1)=CosPhi + SumU*CrossProd%D(1)**2
         Rot(2,2)=CosPhi + SumU*CrossProd%D(2)**2
         Rot(3,3)=CosPhi + SumU*CrossProd%D(3)**2
         !
         Sum1=SumU*CrossProd%D(1)*CrossProd%D(2)
         Sum2=SinPhi*CrossProd%D(3)
         Rot(1,2)=Sum1-Sum2
         Rot(2,1)=Sum1+Sum2
         !
         Sum1=SumU*CrossProd%D(1)*CrossProd%D(3)
         Sum2=SinPhi*CrossProd%D(2)
         Rot(1,3)=Sum1+Sum2
         Rot(3,1)=Sum1-Sum2
         !
         Sum1=SumU*CrossProd%D(2)*CrossProd%D(3)
         Sum2=SinPhi*CrossProd%D(1)
         Rot(2,3)=Sum1-Sum2
         Rot(3,2)=Sum1+Sum2
         !
         ! Test rotation
         !
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot,V2,Vect%D(1:3))
         !
         SumM=DOT_PRODUCT((V1-Vect%D),(V1-Vect%D))
         IF(SumM>1.D-6) THEN
           SinPhi=-SinPhi
           Step=Step+1
           IF(Step > 2) THEN
             Call Halt('Rotation unsuccesful')
           ENDIF
         ELSE
           ! write(*,*) 'rotation succesful'
           EXIT    
         ENDIF
       ENDDO
     ENDIF
     !
     CALL Delete(CrossProd)
     CALL Delete(Vect)
   END SUBROUTINE Rotate
!
!--------------------------------------------------------------------
!
   SUBROUTINE ChkBendToLinB(IntCs,XYZ,AtNum,CtrlCoord,TOPS,IEq)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER                     :: NIntC,Nintc_New
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(INT_VECT)              :: LinAtom,MarkLinb
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER,DIMENSION(:)        :: AtNum
     REAL(DOUBLE)                :: Value,Value2,Conv
     INTEGER                     :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
     INTEGER                     :: I,J,K,L,NatmsLoc
     TYPE(INT_RNK2)              :: LinBBridge
     TYPE(CoordCtrl)             :: CtrlCoord
     REAL(DOUBLE),DIMENSION(3)   :: AuxXYZ
     INTEGER,DIMENSION(:)        :: IEq
     !
     ! Now check for bending -> linear bending transitions
     !
     NIntC=IntCs%N
     IF(NIntC==0) RETURN
     NatmsLoc=SIZE(XYZ,2)
     CALL New(MarkLinB,NIntC)
     Conv=180.D0/PI
     !
     MarkLinB%I=0
          NLinB=0
     DO I=1,NIntC 
       IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         I1=IntCs%Atoms%I(I,1) 
         I2=IntCs%Atoms%I(I,2) 
         I3=IntCs%Atoms%I(I,3) 
         CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value_O=Value)
         IF((ABS(Value-PI)*Conv < CtrlCoord%LinCrit)) THEN  
          !
          ! important for rod-like molecules like H-F polymer 
          ! to have the following lines commented out
          !
          !IF(TOPS%Tot12%I(I2,1)>2) THEN
          ! !IntCs%Def%C(I)(1:5)='BLANK'
          !  IntCs%Active%L(I)=.FALSE.
          !  CYCLE
          !ENDIF
           NLinB=NLinB+1
           MarkLinB%I(I)=1
         ENDIF
       ENDIF
     ENDDO
     !
     IF(NLinB/=0) THEN
       !
       ! Modify IntCs by adding new coords (LinB-s).
       !
       NIntc_New=NIntc+NLinB
     ! NIntc_New=NIntc+2*NLinB
       CALL NEW(IntC_New,NIntc_New)
       !
       NLinB=0
       DO I=1,NIntC
         IF(MarkLinB%I(I)==0) THEN
           CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NLinB+I)
         ELSE
           CALL FindLinBRef(XYZ,IntCs%Atoms%I(I,1:4),TOPS%Tot12%I)
           CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NLinB+I)
           IntC_New%Def%C(NLinB+I)='LINB1'
           NLinB=NLinB+1
           CALL Set_INTC_EQ_INTC(IntCs,IntC_New,I,I,NLinB+I)
           IntC_New%Def%C(NLinB+I)='LINB2'
         ENDIF
       ENDDO
       !
       CALL Delete(IntCs)
       NIntC=NIntC_New
       CALL New(IntCs,NIntC)
         CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntC,1)
       CALL Delete(IntC_New)
     ENDIF
     !
     ! Check for possible cases of long-range torsion
     !
     CALL Delete(MarkLinB)
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:5)=='LINB1') THEN
         I1=IntCs%Atoms%I(I,1)
         I2=IntCs%Atoms%I(I,2)
         I3=IntCs%Atoms%I(I,3)
         L =IntCs%Atoms%I(I,4)
         IF(L==0) THEN
           AuxXYZ=Zero
         ELSE
           AuxXYZ=XYZ(1:3,L)
         ENDIF
         CALL LinB(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),AuxXYZ,L,&
           Value1=IntCs%Value%D(I),Value2=IntCs%Value%D(I+1))  
       ENDIF
     ENDDO 
     !
     ! Now recognize colinear atoms of the molecule and 
     ! introduce long-range torsions. 
     !
     CALL LongRangeIntC(CtrlCoord,TOPS%Tot12,IntCs,XYZ,AtNum,IEq)
     !
   END SUBROUTINE ChkBendToLinB    
!
!-------------------------------------------------------
!
   SUBROUTINE FindLinBRef(XYZ,Atoms,Top12)
     INTEGER,DIMENSION(1:4)  :: Atoms
     INTEGER,DIMENSION(:,:)  :: Top12
     INTEGER                 :: I,J,K,L,I1,I2,I3,N1,N3,IC,N,II
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)            :: Val,ValM,Conv,SumU,RefLinCrit
     !
     Conv=180.D0/PI
     RefLinCrit=20.D0 
     I1=Atoms(1)
     I2=Atoms(2)
     I3=Atoms(3)
     N1=Top12(I1,1)
     N3=Top12(I3,1)
     IF(N1<=1.AND.N3<=1) RETURN
     ValM=Zero
     IF(N1>1) THEN
       II=Top12(I1,2)
       IF(II==I2) II=Top12(I1,3)
       CALL BEND(XYZ(1:3,I2),XYZ(1:3,I1),XYZ(1:3,II),Value_O=ValM)
     ELSE 
       II=Top12(I3,2)
       IF(II==I2) II=Top12(I3,3)
       CALL BEND(XYZ(1:3,I2),XYZ(1:3,I3),XYZ(1:3,II),Value_O=ValM)
     ENDIF
     !
     IF(ABS(180.D0-ValM*Conv)<RefLinCrit) THEN
       Atoms(4)=0 
       RETURN
     ENDIF
     !
     ValM=ABS(120.D0-ValM*Conv)
     DO K=1,2 
       IF(K==1) THEN
         IC=I1
         N=N1
       ELSE
         IC=I3
         N=N3
       ENDIF
       DO I=1,N
         J=Top12(IC,I+1)
         IF(J==I2) CYCLE
         CALL BEND(XYZ(1:3,I2),XYZ(1:3,IC),XYZ(1:3,J),Value_O=Val)
         SumU=ABS(120.D0-Val*Conv)
         IF(SumU<ValM) THEN
           II=J
           ValM=SumU 
         ENDIF
       ENDDO
     ENDDO
     Atoms(4)=II 
   END SUBROUTINE FindLinBRef
!
!-------------------------------------------------------
!
   SUBROUTINE ConstrConv(IntCs,IntDiff,ConstrMax,ConstrRMS)
     TYPE(INTC) :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: IntDiff
     REAL(DOUBLE) :: ConstrMax,ConstrRMS,SumU
     INTEGER      :: I,J,NIntC,NConstr
     !
     NIntC=SIZE(IntCs%Def%C)
     !
     ConstrMax=Zero
     ConstrRMS=Zero
     NConstr=0
     DO I=1,NIntC 
       IF(IntCs%Constraint%L(I)) THEN
         NConstr=NConstr+1
         SumU=ABS(IntDiff(I))
         IF(SumU>ConstrMax) ConstrMax=SumU
         ConstrRMS=ConstrRMS+SumU*SumU
       ENDIF
     ENDDO
     ConstrRMS=SQRT(ConstrRMS/DBLE(NConstr))
   END SUBROUTINE ConstrConv
!
!-------------------------------------------------------------
!
   SUBROUTINE SetFixedLattice(VectCart,IntCs,GConstr)
     REAL(DOUBLE),DIMENSION(:)   :: VectCart
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,K,L,NCart
     REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
     REAL(DOUBLE),DIMENSION(6)   :: Vec
     TYPE(Constr)                :: GConstr
     !
     NCart=SIZE(VectCart)-9
     DO I=1,3
       K=NCart+3*(I-1)+1
       L=K+2
       BoxShape(1:3,I)=VectCart(K:L)
     ENDDO
     CALL CalcBoxPars(Vec,BoxShape)
     CALL SetLattValues(Vec,IntCs,GConstr)
     CALL BoxParsToCart(Vec,BoxShape)
     DO I=1,3
       K=NCart+3*(I-1)+1
       L=K+2
       VectCart(K:L)=BoxShape(1:3,I)
     ENDDO
   END SUBROUTINE SetFixedLattice
!
!--------------------------------------------------------------------
!
   SUBROUTINE SetLattValues(Vec,IntCs,GConstr)
     REAL(DOUBLE),DIMENSION(6) :: Vec
     REAL(DOUBLE),DIMENSION(3) :: Ratio
     TYPE(INTC)                :: IntCs
     INTEGER                   :: I,IRef1,IRef2,IConstr(6)
     TYPE(Constr)              :: GConstr
     REAL(DOUBLE)              :: SumIRef1,SumIRef2,SumR1,SumR2
     !
     IConstr=0
     IRef1=0
     IRef2=0
     DO I=1,IntCs%N
       IF(.NOT.IntCs%Constraint%L(I)) CYCLE
       IF(IntCs%Def%C(I)(1:6)=='STRE_A') THEN
         Vec(1)=IntCs%ConstrValue%D(I)
         IF(GConstr%RatioABC(1)>Zero) IRef1=1
       ENDIF
       IF(IntCs%Def%C(I)(1:6)=='STRE_B') THEN
         Vec(2)=IntCs%ConstrValue%D(I)
         IF(GConstr%RatioABC(2)>Zero) IRef1=2
       ENDIF
       IF(IntCs%Def%C(I)(1:6)=='STRE_C') THEN
         Vec(3)=IntCs%ConstrValue%D(I)
         IF(GConstr%RatioABC(3)>Zero) IRef1=3
       ENDIF
       IF(IntCs%Def%C(I)(1:5)=='ALPHA') THEN
         Vec(4)=IntCs%ConstrValue%D(I)
         IF(GConstr%RatioABC(1)>Zero) IRef2=1
       ENDIF
       IF(IntCs%Def%C(I)(1:4)=='BETA') THEN
         Vec(5)=IntCs%ConstrValue%D(I)
         IF(GConstr%RatioABC(2)>Zero) IRef2=2
       ENDIF
       IF(IntCs%Def%C(I)(1:5)=='GAMMA') THEN
         Vec(6)=IntCs%ConstrValue%D(I)
         IF(GConstr%RatioABC(3)>Zero) IRef2=3
       ENDIF
     ENDDO
     SumIRef1=Zero ; SumR1=Zero
     IF(IRef1==0) THEN
       DO I=1,3
         IF(GConstr%RatioABC(I)>Zero) THEN
           SumIRef1=SumIRef1+Vec(I)
           SumR1=SumR1+GConstr%RatioABC(I)
         ENDIF
       ENDDO 
     ENDIF
     SumIRef2=Zero ; SumR2=Zero
     IF(IRef2==0) THEN
       DO I=1,3
         IF(GConstr%RatioAlpBetGam(I)>Zero) THEN
           SumIRef2=SumIRef2+Vec(3+I)
           SumR2=SumR2+GConstr%RatioAlpBetGam(I)
         ENDIF
       ENDDO 
     ENDIF
     DO I=1,3
       IF(GConstr%RatioABC(I)>Zero) THEN
         IF(IRef1/=0) THEN
           Vec(I)=Vec(IRef1)*GConstr%RatioABC(I)/GConstr%RatioABC(IRef1)
         ELSE
           Vec(I)=GConstr%RatioABC(I)*SumIRef1/SumR1
         ENDIF
       ENDIF
       IF(GConstr%RatioAlpBetGam(I)>Zero) THEN
         IF(IRef2/=0) THEN
           Vec(3+I)=Vec(3+IRef2)*GConstr%RatioAlpBetGam(I)/GConstr%RatioAlpBetGam(IRef2)
         ELSE
           Vec(3+I)=GConstr%RatioAlpBetGam(I)*SumIRef2/SumR2
         ENDIF
       ENDIF
     ENDDO 
   END SUBROUTINE SetLattValues
!
!----------------------------------------------------------------------
!
   SUBROUTINE SetFixedCartesians(Carts,CartDispl,IntCs,NConstr)
     REAL(DOUBLE),DIMENSION(:) :: CartDispl,Carts
     INTEGER                   :: I,J,JJ,NIntC,NConstr
     TYPE(INTC)                :: IntCs
     !
     NIntC=SIZE(IntCs%Def%C)
     !
     ! Make constraints on Cartesians 'hard'
     !
     IF(NConstr>0) THEN
       DO I=1,NIntC
         IF(IntCs%Def%C(I)(1:4)=='CART'.AND.IntCs%Constraint%L(I)) THEN
           J=IntCs%Atoms%I(I,1) 
           JJ=(J-1)*3
           IF(IntCs%Def%C(I)(5:5)=='X') THEN
             JJ=JJ+1
           ELSE IF(IntCs%Def%C(I)(5:5)=='Y') THEN
             JJ=JJ+2
           ELSE IF(IntCs%Def%C(I)(5:5)=='Z') THEN
             JJ=JJ+3
           ENDIF
           Carts(JJ)=IntCs%ConstrValue%D(I)
           CartDispl(JJ)=Zero
         ENDIF
       ENDDO
     ENDIF      
   END SUBROUTINE SetFixedCartesians
!
!---------------------------------------------------------------
!
   SUBROUTINE ScaleDispl(CartDispl,MaxCartDiff,DiffMax,RMSD)
     REAL(DOUBLE),DIMENSION(:) :: CartDispl 
     REAL(DOUBLE)              :: MaxCartDiff,DiffMax,RMSD,SumU
     INTEGER                   :: I,NCart
     !
     NCart=SIZE(CartDispl)
     DiffMax=Zero
     DO I=1,NCart     
       DiffMax=MAX(DiffMax,ABS(CartDispl(I)))
     ENDDO
     !
     ! Scale Displacements
     !
     IF(DiffMax>MaxCartDiff) THEN
       SumU=MaxCartDiff/DiffMax
       CartDispl=SumU*CartDispl
       DiffMax=MaxCartDiff
     ENDIF
     !
     RMSD=DOT_PRODUCT(CartDispl,CartDispl)
     RMSD=SQRT(RMSD/DBLE(NCart))
   END SUBROUTINE ScaleDispl
!
!-------------------------------------------------------
!
   SUBROUTINE BackTrfConvg(GConstr,GBackTrf, &
     DoIterate,DiffMax,RMSD,RMSDOld,ConstrMax, &
     ConstrRMS,ConstrRMSOld,IStep,RefreshAct)
     !
     REAL(DOUBLE) :: DiffMax,RMSD,RMSDOld
     REAL(DOUBLE) :: ConstrMax,ConstrRMS
     REAL(DOUBLE) :: ConstrRMSOld
     INTEGER      :: IStep
     LOGICAL      :: DoIterate,RefreshAct
     LOGICAL      :: ConvConstr
     TYPE(Constr) :: GConstr
     TYPE(BackTrf):: GBackTrf
     !
     DoIterate=.TRUE.
     ConvConstr=.TRUE.
     !
     IF(GConstr%NConstr/=0) THEN
       IF(IStep>1) THEN
         ConvConstr=(ConstrMax<GConstr%ConstrMaxCrit.OR. &
           (ConstrRMS>ConstrRMSOld*GBackTrf%RMSCrit.AND.IStep>5)) 
       ELSE
         ConvConstr=.FALSE.
       ENDIF
     ENDIF
     !
     DoIterate=(DiffMax>GBackTrf%CooTrfCrit)
     IF(RMSD>RMSDOld*GBackTrf%RMSCrit) THEN 
       IF(DiffMax<GBackTrf%MaxCartDiff*GBackTrf%RMSCrit &
          .AND.IStep>GBackTrf%MaxIt_CooTrf) THEN
         !.AND.IStep>20) THEN
         DoIterate=.FALSE.
       ELSE 
         RefreshAct=.TRUE.
       ENDIF
     ENDIF
     DoIterate=DoIterate.OR.(.NOT.ConvConstr)
     DoIterate=(DoIterate.AND.IStep<=GBackTrf%MaxIt_CooTrf)
     !
   END SUBROUTINE BackTrfConvg
!
!----------------------------------------------------------------
!
   SUBROUTINE CenterOfMass(CX,CY,CZ,XYZ_O,Vect_O,Move_O)
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: XYZ_O
     REAL(DOUBLE),DIMENSION(:)  ,OPTIONAL :: Vect_O
     REAL(DOUBLE)                :: CX,CY,CZ
     INTEGER                     :: I,J,NCart,NatmsLoc
     LOGICAL,OPTIONAL            :: Move_O
     !
     IF((PRESENT(XYZ_O).AND.PRESENT(Vect_O)) .OR. &
        (.NOT.PRESENT(XYZ_O).AND..NOT.PRESENT(Vect_O))) THEN
       CALL Halt('Inappropriate input in CenterOfMass')
     ENDIF
     !
     CX=Zero
     CY=Zero
     CZ=Zero
     IF(PRESENT(XYZ_O)) THEN
       NatmsLoc=SIZE(XYZ_O,2)
       DO I=1,NatmsLoc
         CX=CX+XYZ_O(1,I) 
         CY=CY+XYZ_O(2,I) 
         CZ=CZ+XYZ_O(3,I) 
       ENDDO
     ELSE IF(PRESENT(Vect_O)) THEN
       NatmsLoc=SIZE(Vect_O)/3
       NCart=3*NatmsLoc
       DO I=1,NatmsLoc
         J=(I-1)*3
         CX=CX+Vect_O(J+1) 
         CY=CY+Vect_O(J+2) 
         CZ=CZ+Vect_O(J+3) 
       ENDDO
     ENDIF
     CX=CX/DBLE(NatmsLoc)
     CY=CY/DBLE(NatmsLoc)
     CZ=CZ/DBLE(NatmsLoc)
     !
     ! Move coordinates into Center of Mass?
     !
     IF(PRESENT(Move_O)) THEN
       IF(Move_O) THEN
         IF(PRESENT(XYZ_O)) THEN
           XYZ_O(1,:)=XYZ_O(1,:)-CX
           XYZ_O(2,:)=XYZ_O(2,:)-CY
           XYZ_O(3,:)=XYZ_O(3,:)-CZ
         ELSE IF(PRESENT(Vect_O)) THEN
           DO I=1,NatmsLoc
             J=(I-1)*3
             Vect_O(J+1)=Vect_O(J+1)-CX
             Vect_O(J+2)=Vect_O(J+2)-CY
             Vect_O(J+3)=Vect_O(J+3)-CZ
           ENDDO
         ENDIF
       ENDIF
     ENDIF
   END SUBROUTINE CenterOfMass
!
!-------------------------------------------------------------------
!
   SUBROUTINE PrincMomTens(CMCarts,Theta)
     !
     ! Calculates inertial momentum tensor
     ! Pass in mass-centered Cartesians
     !
     REAL(DOUBLE),DIMENSION(:,:) :: Theta,CMCarts
     INTEGER                     :: I,NatmsLoc
     REAL(DOUBLE)                :: X,Y,Z,XX,YY,ZZ,XY,YZ,ZX
     !                
     NatmsLoc=SIZE(CMCarts,2)
     !                
     Theta=Zero
     DO I=1,NatmsLoc
       X=CMCarts(1,I)
       Y=CMCarts(2,I)
       Z=CMCarts(3,I)
       XX=X*X
       YY=Y*Y
       ZZ=Z*Z
       XY=X*Y
       YZ=Y*Z
       ZX=Z*X
       Theta(1,1)=Theta(1,1)+YY+ZZ
       Theta(2,2)=Theta(2,2)+ZZ+XX
       Theta(3,3)=Theta(3,3)+XX+YY
       Theta(1,2)=Theta(1,2)-XY    
       Theta(2,1)=Theta(2,1)-XY    
       Theta(1,3)=Theta(1,3)-ZX    
       Theta(3,1)=Theta(3,1)-ZX    
       Theta(2,3)=Theta(2,3)-YZ    
       Theta(3,2)=Theta(3,2)-YZ    
     ENDDO
   END SUBROUTINE PrincMomTens
!
!------------------------------------------------------------------
!
   SUBROUTINE RotateMol(XYZ,U33)
     !
     ! Rotate Cartesian set by U33 Unitary matrix
     ! XYZ(I)=XYZ(I)*U33
     !
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,U33
     INTEGER                     :: I,J,K,NatmsLoc
     TYPE(DBL_VECT)              :: AuxVect
     !
     NatmsLoc=SIZE(XYZ,2)
     !
     CALL New(AuxVect,3)
     DO I=1,NatmsLoc
       CALL DGEMM_NNc(1,3,3,One,Zero,XYZ(1:3,I),U33,AuxVect%D)
       XYZ(1:3,I)=AuxVect%D 
     ENDDO
     CALL Delete(AuxVect)
   END SUBROUTINE RotateMol
!
!-------------------------------------------------------------------
!
   SUBROUTINE CALC_BxVect(ISpB,JSpB,ASpB,VectInt,VectCart,Trp_O)
     TYPE(INT_VECT)            :: ISpB,JSpB,ISpBt,JSpBt
     TYPE(DBL_VECT)            :: ASpB,ASpBt
     REAL(DOUBLE),DIMENSION(:) :: VectInt,VectCart
     INTEGER                   :: NCart,NIntC,I,J,JJ,K,KK,LL,NZSpB
     REAL(DOUBLE)              :: SumU
     LOGICAL,OPTIONAL          :: Trp_O
     !
     NCart=SIZE(VectCart)
     NIntC=SIZE(VectInt)
     NZSpB=SIZE(JSpB%I)
     !
     IF(PRESENT(Trp_O)) THEN
       IF(Trp_O) THEN !!! Bt*VectInt
         CALL New(ISpBt,NCart+1)
         CALL New(JSpBt,NZSpB)
         CALL New(ASpBt,NZSpB)
         CALL TransPose1x1(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart, &
            ISpBt%I,JSpBt%I,ASpBt%D,'full')
         CALL ProdMatrVect(ISpBt%I,JSpBt%I,ASpBt%D, &
                           VectInt,NCart,VectCart) 
         CALL Delete(ISpBt)
         CALL Delete(JSpBt)
         CALL Delete(ASpBt)
       ENDIF
     ELSE
         CALL ProdMatrVect(ISpB%I,JSpB%I,ASpB%D, &
                           VectCart,NIntC,VectInt) 
     ENDIF 
   END SUBROUTINE CALC_BxVect
!
!----------------------------------------------------------------------
!
   SUBROUTINE CALC_GcInvCartV(CholData,VectCartAux,VectCartAux2)
     TYPE(Cholesky)            :: CholData
     REAL(DOUBLE),DIMENSION(:) :: VectCArtAux,VectCartAux2
     TYPE(DBL_VECT)            :: VectCartAux3
     TYPE(DBL_VECT)            :: Vect1,Vect2,Vect3
     INTEGER                   :: NCart,NDim,NIntC,I
     !
     NCart=SIZE(VectCArtAux)
     CALL New(VectCartAux3,NCart)
     VectCartAux2=VectCartAux
     !
     CALL PermVect(VectCartAux2,VectCartAux3%D,CholData%Perm%I)
     CALL ScaleVect(VectCartAux3%D,CholData%GcScale%D)
     CALL CholFactSolve(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,CholData%ChDiag%D, &
       VectCartAux3%D,NCart,VectCartAux2)
     VectCartAux3%D=VectCartAux2
     CALL ScaleVect(VectCartAux3%D,CholData%GcScale%D)
     CALL PermVect(VectCartAux3%D,VectCartAux2,CholData%IPerm%I)
     CALL Delete(VectCartAux3)
   END SUBROUTINE CALC_GcInvCartV
!
!----------------------------------------------------------------------
!
   SUBROUTINE GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData, &
                          IPerm1_O,IPerm2_O,B_O)
     TYPE(Cholesky)          :: CholData
     CHARACTER(LEN=*)        :: SCRPath
     TYPE(INT_VECT)          :: ISpB,JSpB
     TYPE(DBL_VECT)          :: ASpB
     TYPE(INT_VECT),OPTIONAL :: IPerm1_O,IPerm2_O
     TYPE(BMATR),OPTIONAL    :: B_O
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B',&
                    IPerm1_O=IPerm1_O,IPerm2_O=IPerm2_O,B_O=B_O)
     CALL ReadChol(CholData,TRIM(SCRPath)//'CholFact')
   END SUBROUTINE GetBMatInfo
!
!-------------------------------------------------------------------
!
   SUBROUTINE GetMixedBMat(IntCs,XYZ,PBC,GTrfCtrl,GCoordCtrl,TOPS, &
                           NCartConstr,Print,SCRPath,DoConstr)
     TYPE(INTC)                   :: IntCs
     TYPE(TOPOLOGY)               :: TOPS
     TYPE(Cholesky)               :: CholData
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     INTEGER                      :: NatmsLoc,NCart,NIntC,At1,At2,At3
     INTEGER                      :: I,J,NCartConstr
     TYPE(BMATR)                  :: B,B1,B2
     INTEGER                      :: Print,ILowDim
     LOGICAL                      :: Print2,DoConstr
     CHARACTER(LEN=*)             :: SCRPath
     TYPE(INT_VECT)               :: IPerm1,IPerm2,InvP,Perm
     TYPE(INT_VECT)               :: ISpB1,ISpB2,JSpB1,JSpB2,IGc1,JGc1
     TYPE(INT_VECT)               :: ISpBMix,JSpBMix
     TYPE(DBL_VECT)               :: ASpB1,ASpB2,ASpBMix,AGc1
     TYPE(DBL_RNK2)               :: XYZWork
     TYPE(PBCInfo)                :: PBC
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def%C)
     Print2=(Print>=DEBUG_GEOP_MAX)
     !
     CALL New(XYZWork,(/3,NatmsLoc/))
     XYZWork%D=XYZ
     !
     ! Find two reference systems
     CALL ThreeAtoms(XYZWork%D,TOPS%Tot12, &
                     GTrfCtrl%ThreeAt,GTrfCtrl%Linearity,IntCs_O=IntCs)
       At1=GTrfCtrl%ThreeAt(1)
       At2=GTrfCtrl%ThreeAt(2)
       At3=GTrfCtrl%ThreeAt(3)
     IF(NatmsLoc<=6.OR.(NCartConstr/=0.AND.NCartConstr<18)) THEN
       GTrfCtrl%ThreeAt_2(1)=At1
       GTrfCtrl%ThreeAt_2(2)=At2
       GTrfCtrl%ThreeAt_2(3)=At3
     ELSE
       TOPS%Tot12%I(At1,:)=0
       TOPS%Tot12%I(At2,:)=0
       TOPS%Tot12%I(At3,:)=0
       CALL ThreeAtoms(XYZWork%D,TOPS%Tot12,GTrfCtrl%ThreeAt_2, &
                       GTrfCtrl%Linearity,IntCs_O=IntCs)
     ENDIF
     CALL CALC_XYZRot(XYZWork%D,GTrfCtrl%ThreeAt,&
                      GTrfCtrl%Linearity,GTrfCtrl%TranslAt1, &
                      GTrfCtrl%RotAt2ToX,GTrfCtrl%RotAt3ToXY)
     CALL CALC_XYZRot(XYZWork%D,GTrfCtrl%ThreeAt_2,&
                      GTrfCtrl%Linearity,GTrfCtrl%TranslAt1_2, &
                      GTrfCtrl%RotAt2ToX_2,GTrfCtrl%RotAt3ToXY_2)
     !
     ! Calculate B matrix in Atomic Units in absolute coordinate system
     !
     CALL BMatrix(XYZ,IntCs,B,PBC%Dimen,GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
     !
     ! Calculate rotated B matrices
     !
     CALL Set_BMATR_EQ_BMATR(B1,B)
     CALL Set_BMATR_EQ_BMATR(B2,B)
     CALL RotB(B1,GTrfCtrl%RotAt2ToX,GTrfCtrl%RotAt3ToXY)
     CALL RotB(B2,GTrfCtrl%RotAt2ToX_2,GTrfCtrl%RotAt3ToXY_2)
     CALL BtoSpB_1x1(B1,ISpB1,JSpB1,ASpB1)
     CALL BtoSpB_1x1(B2,ISpB2,JSpB2,ASpB2)
     !
     ! Calculate RCM reordering
     !
     CALL New(InvP,NCart)
     CALL New(Perm,NCart)
     CALL GetGc(NCart,ISpB1,JSpB1,ASpB1,IGc1,JGc1,AGc1,SymbOnly_O=.TRUE.)
    !CALL GetGc(NCart,ISpB1,JSpB1,ASpB1,IGc1,JGc1,AGc1)
     CALL RCMOrder(InvP%I,Perm%I,NCart,IGc1%I,JGc1%I)
     CALL Delete(IGc1)
     CALL Delete(JGc1)
     !
     ! Clean columns of constraints and references
     !
     CALL New(IPerm1,NCart)
     CALL New(IPerm2,NCart)
     IPerm1%I=0
     IPerm2%I=0
     !
     CALL PermCleans(IPerm1%I,Perm%I,IntCs,GTrfCtrl%ThreeAt,DoConstr)
     ILowDim=MAXVAL(IPerm1%I)
     CALL PermCleans(IPerm2%I,Perm%I,IntCs,GTrfCtrl%ThreeAt_2,DoConstr)
     IF(ILowDim/=MAXVAL(IPerm2%I)) &
       CALL Halt('Dimension error in GetMixedBMat')
    !CALL Show1x1(ISpB1%I,JSpB1%I,ASpB1%D,'original B1',NIntC,NCart)
   ! CALL Plot_1x1(ISpB1%I,JSpB1%I,'SpB1',NatmsLoc)
     CALL PermCleanB(ISpB1,JSpB1,ASpB1,IPerm1%I)
    !CALL Show1x1(ISpB1%I,JSpB1%I,ASpB1%D,'permuted B1',NIntC,NCart)
   ! CALL Plot_1x1(ISpB2%I,JSpB2%I,'SpB2',NatmsLoc)
     CALL PermCleanB(ISpB2,JSpB2,ASpB2,IPerm2%I)
     !
     CALL AddMat_1x1(ISpB1%I,JSpB1%I,ASpB1%D, &
       ISpB2%I,JSpB2%I,ASpB2%D,ISpBMix,JSpBMix,ASpBMix,NIntC,NCart)
     !
     CALL Delete(B1)
     CALL Delete(ISpB1)
     CALL Delete(JSpB1)
     CALL Delete(ASpB1)
     CALL Delete(B2)
     CALL Delete(ISpB2)
     CALL Delete(JSpB2)
     CALL Delete(ASpB2)
     CALL Delete(B)
     CALL Delete(XYZWork)
     CALL Delete(InvP)
     CALL Delete(Perm)
     !
     CALL WriteBMATR(ISpBMix,JSpBMix,ASpBMix,TRIM(SCRPath)//'MixB',&
                     IPerm1_O=IPerm1,IPerm2_O=IPerm2)
     CALL Delete(IPerm1)
     CALL Delete(IPerm2)
     !
     CALL CholFact(ISpBMix,JSpBMix,ASpBMix,NCart,NIntC, &
                   CholData,Print2,ILow_O=ILowDim)
     CALL WriteChol(CholData,TRIM(SCRPath)//'MixCholFact')
     !
     CALL DeleteBMatInfo(ISpBMix,JSpBMix,ASpBMix,CholData)
   END SUBROUTINE GetMixedBMat
!
!---------------------------------------------------------------------
!
   SUBROUTINE RefreshBMatInfo(IntCs,XYZ,GTrfCtrl,GConvCr, &
                              LinCrit,TorsLinCrit,PBCDim,Print, &
                              SCRPath,Gi_O,Shift_O,DoCleanCol_O)
     TYPE(INTC)                   :: IntCs
     TYPE(Cholesky)               :: CholData
     TYPE(TrfCtrl)                :: GTrfCtrl
     TYPE(GConvCrit)              :: GConvCr
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     REAL(DOUBLE)                 :: LinCrit,TorsLinCrit
     REAL(DOUBLE),OPTIONAL        :: Shift_O
     REAL(DOUBLE)                 :: ShiftU
     INTEGER                      :: PBCDim,I,K,L,J
     INTEGER                      :: NatmsLoc,NCart,NIntC,NZ
     TYPE(BMATR)                  :: B,BS
     INTEGER                      :: Print
     LOGICAL                      :: Print2,DoCleanB,DoGi
     CHARACTER(LEN=*)             :: SCRPath
     TYPE(INT_VECT)               :: ISpB,JSpB,ISpBt,JSpBt,ISpBC,JSpBC
     TYPE(DBL_VECT)               :: ASpB,ASpBt,ASpBC
     LOGICAL,OPTIONAL             :: Gi_O,DoCleanCol_O
     LOGICAL                      :: DoCleanCol
     REAL(DOUBLE),DIMENSION(9)    :: AuxBL
     REAL(DOUBLE),DIMENSION(9,9)  :: P    
     !
     DoGi=.FALSE.
     IF(PRESENT(Gi_O)) DoGi=Gi_O
     DoCleanCol=.TRUE.
     IF(PRESENT(DoCleanCol_O)) DoCleanCol=DoCleanCol_O
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def%C)
     Print2=(Print>=DEBUG_GEOP_MAX)
     ShiftU=1.D-6
     IF(PRESENT(Shift_O)) THEN
       ShiftU=Shift_O
     ENDIF
     !
     CALL BMatrix(XYZ,IntCs,B,PBCDim,LinCrit,TorsLinCrit)
     CALL SetEq(BS,B)
     !
     ! Clean Cartesian constraints
     ! Clean Lattice constraints
     !
     CALL CleanBConstr(IntCs,B,NatmsLoc)
     CALL CleanBLConstr(XYZ,IntCs,B,PBCDim)
     CALL CleanBLRot(XYZ,IntCs,B,PBCDim)
     !
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     !
     CALL WriteBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B',B_O=BS)
     CALL Delete(BS)
     !
     IF(DoGi) THEN
       CALL CholFactGi(ISpB,JSpB,ASpB,NCart,NIntC, &
                     CholData,Print2,Shift_O=ShiftU)
     ELSE
       CALL CholFact(ISpB,JSpB,ASpB,NCart,NIntC, &
                     CholData,Print2,Shift_O=ShiftU)
     ENDIF
     CALL WriteChol(CholData,TRIM(SCRPath)//'CholFact')
     !
     CALL Delete(CholData)
     CALL Delete(B)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
   END SUBROUTINE RefreshBMatInfo
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanBLConstr(XYZ,IntCs,B,PBCDim)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(9,9) :: P   
     REAL(DOUBLE),DIMENSION(9)   :: Vect1,Vect2
     REAL(DOUBLE)                :: Trace
     TYPE(INTC)                  :: IntCs
     TYPE(BMATR)                 :: B
     INTEGER                     :: I,J,PBCDim
     !
     ! Get projector to constraints
     !
     CALL GetLattProj(XYZ,IntCs,PBCDim,P) 
     P=-P
     DO I=1,9 ; P(I,I)=One+P(I,I) ; ENDDO
     DO I=1,IntCs%N
       IF(B%BLI%I(I)/=0) THEN
         Vect1(1:9)=B%BL%D(I,1:9)
         CALL DGEMM_NNc(9,9,1,One,Zero,P,Vect1,Vect2)
         B%BL%D(I,1:9)=Vect2(1:9)
       ENDIF 
     ENDDO
   END SUBROUTINE CleanBLConstr
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanBLRot(XYZ,IntCs,B,PBCDim)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(9,9) :: P   
     REAL(DOUBLE),DIMENSION(9)   :: Vect1,Vect2
     REAL(DOUBLE)                :: Trace
     TYPE(INTC)                  :: IntCs
     TYPE(BMATR)                 :: B
     INTEGER                     :: I,J,PBCDim
     !
     ! Get projector to counterspace of lattice rotations
     !
   ! CALL GetPBCProj(PBCDim,P,XYZ_O=XYZ)
     DO I=1,IntCs%N
       IF(B%BLI%I(I)/=0) THEN
        !Vect1(1:9)=B%BL%D(I,1:9)
        !CALL DGEMM_NNc(9,9,1,One,Zero,P,Vect1,Vect2)
        !B%BL%D(I,1:9)=Vect2(1:9)
         B%BL%D(I,2)=Zero
         B%BL%D(I,3)=Zero
         B%BL%D(I,6)=Zero
       ENDIF 
     ENDDO
   END SUBROUTINE CleanBLRot
!
!---------------------------------------------------------------------
!
   SUBROUTINE GetLattProj(XYZ,IntCsE,PBCDim,P)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INTC)                  :: IntCsE
     INTEGER                     :: PBCDim,NCoinc,NatmsLoc,I,J,K,L
     INTEGER                     :: NZ,N,NZP
     TYPE(DBL_RNK2)              :: AL,BL
     REAL(DOUBLE)                :: Trace   
     REAL(DOUBLE),DIMENSION(9,9) :: P
     REAL(DOUBLE),DIMENSION(9)   :: AuxBL
     !
     ! Construct Projector of lattice constraints
     ! without lattice rotations PRot
     !
     NatmsLoc=SIZE(XYZ,2)
     !
     IF(PBCDim>0) THEN
       CALL LatticeConstrAB(XYZ,IntCsE,PBCDim,AL,BL,NCoinc)
         DO I=1,3
           K=3*(I-1)+1
           L=K+2
           AuxBL(K:L)=XYZ(1:3,NatmsLoc-3+I)
         ENDDO
       IF(NCoinc>0) THEN
         CALL DGEMM_TNc(9,NCoinc,9,One,Zero,AL%D,BL%D,P)
         CALL Delete(AL)
         CALL Delete(BL)
       ELSE
         P=Zero
       ENDIF
     ELSE
       P=Zero
     ENDIF
   END SUBROUTINE GetLattProj
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanHardConstr(ISpBC,JSpBC,ASpBC,NZ,NCart,NZP,P,IntCs)
     TYPE(INT_VECT)              :: ISpBC,JSpBC,ISpP,JSpP,IConstr
     TYPE(INT_VECT)              :: ISpB2,JSpB2
     TYPE(DBL_VECT)              :: ASpBC,ASpP,ASpB2
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(9,9) :: P
     INTEGER                     :: NCart,NZP,NZ,NConstr,I,J,K
     !
     CALL New(IConstr,NCart)
     IConstr%I=0  
     NConstr=0
     DO I=1,IntCs%N
       IF(IntCs%Constraint%L(I)) THEN
         J=3*(IntCs%Atoms%I(I,1)-1)
         IF(IntCs%Def%C(I)(1:5)=='CARTX') THEN
           NConstr=NConstr+1
           IConstr%I(J+1)=1 
         ENDIF
         IF(IntCs%Def%C(I)(1:5)=='CARTY') THEN
           NConstr=NConstr+1
           IConstr%I(J+2)=1
         ENDIF
         IF(IntCs%Def%C(I)(1:5)=='CARTZ') THEN
           NConstr=NConstr+1
           IConstr%I(J+3)=1
         ENDIF
       ENDIF
     ENDDO 
     !
     NZ=NConstr+NZP
     IF(NZ==0) THEN
       CALL Delete(IConstr)
       RETURN
     ENDIF
     ! 
     CALL New(ISpP,NCart+1)
     CALL New(JSpP,NZ)
     CALL New(ASpP,NZ)
     !
     NZ=0
     ISpP%I(1)=1
     DO I=1,NCart-9       
       IF(IConstr%I(I)==1) THEN
         NZ=NZ+1
         JSpP%I(NZ)=I  
         ASpP%D(NZ)=One
       ENDIF
       ISpP%I(I+1)=NZ+1
     ENDDO
     CALL Delete(IConstr)
     !
     K=NCart-9
     DO I=1,9
       IF(NZP>0) THEN
         DO J=1,9
           NZ=NZ+1
           JSpP%I(NZ)=K+J
           ASpP%D(NZ)=P(I,J)
         ENDDO 
       ENDIF
       ISpP%I(K+I+1)=NZ+1
     ENDDO
     !
     CALL MatMul_1x1(ISpBC%I,JSpBC%I,ASpBC%D,ISpP%I,JSpP%I,ASpP%D, &
                    ISpB2,JSpB2,ASpB2,IntCs%N,NCart,NCart)
     !
     CALL Delete(ISpBC) ; CALL Delete(JSpBC) ; CALL Delete(ASpBC)
     NZ=ISpB2%I(IntCs%N+1)-1
     CALL New(ISpBC,IntCs%N+1)
     CALL New(JSpBC,NZ)
     CALL New(ASpBC,NZ)
     ISpBC%I(1:IntCs%N+1)=ISpB2%I(1:IntCs%N+1)
     JSpBC%I(1:NZ)=JSpB2%I(1:NZ)
     ASpBC%D(1:NZ)=ASpB2%D(1:NZ)
     ! 
     CALL Delete(ISpP) ; CALL Delete(JSpP) ; CALL Delete(ASpP)
     CALL Delete(ISpB2) ; CALL Delete(JSpB2) ; CALL Delete(ASpB2)
   END SUBROUTINE CleanHardConstr
!
!---------------------------------------------------------------------
!
   SUBROUTINE ProjectFromMatr(ISpBC,JSpBC,ASpBC,ISpB,JSpB,ASpB, &
                              NCart,NIntC)
     TYPE(INT_VECT)              :: ISpB,JSpB,IConstr
     TYPE(INT_VECT)              :: ISpBC,JSpBC
     TYPE(INT_VECT)              :: ISpBCT,JSpBCT
     TYPE(INT_VECT)              :: ISpB1,JSpB1
     TYPE(INT_VECT)              :: ISpB2,JSpB2
     TYPE(DBL_VECT)              :: ASpB,ASpBC,ASpBCT,ASpB1,ASpB2
     INTEGER                     :: NCart,NIntC,NZSpB,I,J,NZ
     TYPE(Cholesky)              :: CholData
     LOGICAL                     :: Print2
     !
     Print2=.TRUE.
     CALL CholFact(ISpBC,JSpBC,ASpBC,NCart,NIntC, &
                   CholData,Print2,Shift_O=1.D-6)
     !
     NZSpB=SIZE(JSpBC%I)
     CALL New(ISpBCT,NCart+1)
     CALL New(JSpBCT,NZSpB)
     CALL New(ASpBCT,NZSpB)
     CALL TransPose1x1(ISpBC%I,JSpBC%I,ASpBC%D,NIntC,NCart, &
        ISpBCT%I,JSpBCT%I,ASpBCT%D,'full')
     !
     CALL MatMul_1x1(ISpBCT%I,JSpBCT%I,ASpBCT%D, &
                     ISpB%I,JSpB%I,ASpB%D, &
                     ISpB1,JSpB1,ASpB1,NCart,NIntC,NCart)
     CALL Delete(ISpBCT) ; CALL Delete(JSpBCT) ; CALL Delete(ASpBCT)
     !
     CALL GcInvItMatr(ISpB1,JSpB1,ASpB1, &
                      ISpBC,JSpBC,ASpBC,CholData,NCart)
     !
     CALL MatMul_1x1(ISpBC%I,JSpBC%I,ASpBC%D, &
                     ISpB1%I,JSpB1%I,ASpB1%D, &
                     ISpB2,JSpB2,ASpB2,NIntC,NCart,NCart)
     !
     CALL Delete(ISpB1) ; CALL Delete(JSpB1) ; CALL Delete(ASpB1)
     ASpB2%D=-ASpB2%D
     CALL AddMat_1x1(ISpB%I,JSpB%I,ASpB%D, &
       ISpB2%I,JSpB2%I,ASpB2%D,ISpB1,JSpB1,ASpB1,NIntC,NCart)
     CALL Delete(ISpB) ; CALL Delete(JSpB) ; CALL Delete(ASpB)
     NZSpB=ISpB1%I(NIntC+1)-1
     CALL New(ISpB,NIntC+1)
     CALL New(JSpB,NZSpB)
     CALL New(ASpB,NZSpB)
     ISpB%I(1:NIntC+1)=ISpB1%I(1:NIntC+1) 
     JSpB%I(1:NZSpB)=JSpB1%I(1:NZSpB) 
     ASpB%D(1:NZSpB)=ASpB1%D(1:NZSpB)
     !
     CALL Delete(ISpB1) ; CALL Delete(JSpB1) ; CALL Delete(ASpB1)
     CALL Delete(ISpB2) ; CALL Delete(JSpB2) ; CALL Delete(ASpB2)
     CALL Delete(CholData)
     !
   END SUBROUTINE ProjectFromMatr
!
!---------------------------------------------------------------------
!
   SUBROUTINE GcInvItMatr(IM,JM,AM,ISpBC,JSpBC,ASpBC, &
                          CholData,NCart)
     TYPE(INT_VECT)              :: ISpBC,JSpBC
     TYPE(INT_VECT)              :: IM,JM,IMP,JMP,IMP2,JMP2,ICorr,JCorr
     TYPE(INT_VECT)              :: ICorr2,JCorr2
     TYPE(DBL_VECT)              :: AM,AMP,AMP2,ACorr
     TYPE(DBL_VECT)              :: ASpBC
     TYPE(DBL_VECT)              :: ACorr2
     TYPE(INT_VECT)              :: ISpB,JSpB
     TYPE(DBL_VECT)              :: ASpB
     TYPE(INT_VECT)              :: IGc,JGc  
     TYPE(DBL_VECT)              :: AGc 
     TYPE(Cholesky)              :: CholData
     INTEGER                     :: NCart,MaxIt,NZ,I,J
     REAL(DOUBLE)                :: ConvCrit,MaxCorr
     !
     ConvCrit=1.D-6
     CALL GetGc(NCart,ISpBC,JSpBC,ASpBC,IGc,JGc,AGc)
     MaxIt=10
     !
     CALL New(IMP,NCart+1)
     CALL New(JMP,1)
     CALL New(AMP,1)
     IMP%I=1
     DO I=1,MaxIt
       ! MP2=Gc*MP
       CALL MatMul_1x1(IGc%I,JGc%I,AGc%D,IMP%I,JMP%I,AMP%D, &
                       IMP2,JMP2,AMP2,NCart,NCart,NCart)
       ! Corr=MP2-M
       AMP2%D=-AMP2%D
       CALL AddMat_1x1(IMP2%I,JMP2%I,AMP2%D, &
            IM%I,JM%I,AM%D,ICorr,JCorr,ACorr,NCart,NCart)
         CALL Delete(IMP2) ; CALL Delete(JMP2) ; CALL Delete(AMP2)
       ! Corr2=Chol*(MP2-M)*Chol
       CALL InvMatXMatr(CholData,ICorr%I,JCorr%I,ACorr%D, &
                        ICorr2,JCorr2,ACorr2,NCart,NCart,1.D-8)
         CALL Delete(ICorr) ; CALL Delete(JCorr) ; CALL Delete(ACorr)
       ! MPnew=MP+Corr2
       CALL AddMat_1x1(ICorr2%I,JCorr2%I,ACorr2%D, &
            IMP%I,JMP%I,AMP%D,IMP2,JMP2,AMP2,NCart,NCart)
         CALL Delete(IMP) ; CALL Delete(JMP) ; CALL Delete(AMP)
         NZ=SIZE(JMP2%I)
         CALL New(IMP,NCart+1)
         CALL New(JMP,NZ)
         CALL New(AMP,NZ)
         IMP%I=IMP2%I
         JMP%I=JMP2%I
         AMP%D=AMP2%D
         MaxCorr=MAXVAL(ACorr2%D)
         CALL Delete(ICorr2) ; CALL Delete(JCorr2) ; CALL Delete(ACorr2)
         CALL Delete(IMP2) ; CALL Delete(JMP2) ; CALL Delete(AMP2)
       WRITE(*,*) I,' CorrMax= ',MaxCorr
       IF(MaxCorr<ConvCrit) EXIT
     ENDDO
     !
     CALL Delete(IGc) ; CALL Delete(JGc) ; CALL Delete(AGc)
     CALL Delete(IM)  ; CALL Delete(JM)  ; CALL Delete(AM)
     NZ=SIZE(JMP%I)
     CALL New(IM,NCart+1)
     CALL New(JM,NZ)
     CALL New(AM,NZ)
     IM%I=IMP%I
     JM%I=JMP%I
     AM%D=AMP%D
     CALL Delete(IMP) ; CALL Delete(JMP) ; CALL Delete(AMP)
   END SUBROUTINE GcInvItMatr
!
!---------------------------------------------------------------------
!
   SUBROUTINE LatticeConstrAB(XYZ,IntCsE,PBCDim,AL,BL,NCoinc)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INTC)                  :: IntCsE,IntCs
     TYPE(DBL_RNK2)              :: AL,BL,Aux,Inv
     INTEGER                     :: PBCDim,NCoinc,I,J,II
     TYPE(BMATR)                 :: B
     !
     IF(PBCDim==0) THEN
       NCoinc=0
       RETURN
     ENDIF
     CALL LatticeINTC(IntCs,PBCDim,DoVolume_O=.TRUE.)
     IntCs%Active%L=.FALSE.
     NCoinc=0
     DO I=1,IntCs%N
       IF(IntCsE%Def%C(I)(1:5)=="BLANK") CYCLE
       DO J=1,IntCsE%N
         IF(IntCsE%Constraint%L(J)) THEN
           IF(IntCsE%Def%C(J)(1:8)==IntCs%Def%C(I)(1:8)) THEN
             IntCs%Active%L(I)=.TRUE.
             NCoinc=NCoinc+1
             EXIT
           ENDIF
         ENDIF
       ENDDO
     ENDDO
     !
     IF(NCoinc/=0) THEN
       CALL New(BL,(/NCoinc,9/))
       CALL New(AL,(/NCoinc,9/))
       CALL New(Aux,(/NCoinc,NCoinc/))
       CALL New(Inv,(/NCoinc,NCoinc/))
       CALL BMatrix(XYZ,IntCs,B,PBCDim,1.D0,1.D0)
       II=0
       DO I=1,IntCs%N
         IF(IntCs%Active%L(I)) THEN
           II=II+1
           BL%D(II,1:9)=B%BL%D(I,1:9)
         ENDIF
       ENDDO
       CALL DGEMM_NTc(NCoinc,9,NCoinc,One,Zero,BL%D,BL%D,Aux%D)
       CALL SetDSYEVWork(NCoinc*NCoinc)
       CALL FunkOnSqMat(NCoinc,Inverse,Aux%D,Inv%D,EigenThresh_O=1.D-10)
       CALL UnSetDSYEVWork()
       CALL DGEMM_NNc(NCoinc,NCoinc,9,One,Zero,Inv%D,BL%D,AL%D)
       CALL Delete(Aux)
       CALL Delete(Inv)
       CALL Delete(B)
     ENDIF
     CALL Delete(IntCs)
   END SUBROUTINE LatticeConstrAB
!
!---------------------------------------------------------------------
!
   SUBROUTINE DeleteBMatInfo(ISpB,JSpB,ASpB,CholData)
     TYPE(Cholesky) :: CholData
     TYPE(INT_VECT) :: ISpB,JSpB
     TYPE(DBL_VECT) :: ASpB
     !
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(CholData)
   END SUBROUTINE DeleteBMatInfo
!
!---------------------------------------------------------------------
!
   SUBROUTINE CALC_XYZRot(XYZ,ThreeAt,Linearity, &
                          TranslAt1,RotAt2ToX,RotAt3ToXY,DoCopy_O)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     REAL(DOUBLE),DIMENSION(3)    :: TranslAt1
     REAL(DOUBLE),DIMENSION(3,3)  :: RotAt2ToX,RotAt3ToXY
     TYPE(DBL_RNK2)   :: XYZRot
     TYPE(DBL_VECT)   :: Vect1,Vect2,Vect3
     TYPE(DBL_RNK2)   :: Rot
     INTEGER          :: I,J,NatmsLoc,NMax12
     INTEGER          :: At1,At2,At3,ThreeAt(1:3)
     LOGICAL          :: Linearity
     LOGICAL,OPTIONAL :: DoCopy_O
     !
     ! Subroutine to calculate rotation matrices of reference systems
     !
     NatmsLoc=SIZE(XYZ,2)
     CALL New(XYZRot,(/3,NatmsLoc/))
     CALL New(Vect1,3)
     CALL New(Vect2,3)
     CALL New(Vect3,3)
     CALL New(Rot,(/3,3/))
     XYZRot%D=XYZ
     !
     TranslAt1=Zero 
     RotAt2ToX=Zero
     RotAt3ToXY=Zero
     DO I=1,3 
       RotAt2ToX(I,I)=One
       RotAt3ToXY(I,I)=One
     ENDDO
     !
     At1=ThreeAt(1)
     At2=ThreeAt(2)
     At3=ThreeAt(3)
     !
     ! Place At1 to origin
     !
     Vect1%D=XYZ(:,At1)
     DO I=1,NatmsLoc
       XYZRot%D(:,I)=XYZ(:,I)-Vect1%D   
     ENDDO
     TranslAt1=Vect1%D
     !
     ! Place At2 onto X axis
     !
     Vect1%D=Zero
     Vect1%D(1)=One
     Vect2%D=XYZRot%D(:,At2)-XYZRot%D(:,At1) 
     CALL Rotate(Vect1%D,Vect2%D,Rot%D)
     RotAt2ToX=Rot%D
     DO I=1,NatmsLoc
       Vect1%D=XYZRot%D(:,I)
       CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,Vect1%D,XYZRot%D(:,I))
     ENDDO
     !
     ! Linearity
     !
     IF(.NOT.Linearity) THEN         
       !
       ! Place At3 onto XY plane
       !
       Vect1%D=XYZRot%D(:,At1)-XYZRot%D(:,At3)
       Vect2%D=XYZRot%D(:,At2)-XYZRot%D(:,At3)
       CALL CROSS_PRODUCT(Vect1%D,Vect2%D,Vect3%D)
       Vect1%D=Zero
       Vect1%D(3)=One 
       CALL Rotate(Vect1%D,Vect3%D,Rot%D)
       RotAt3ToXY=Rot%D
       DO I=1,NatmsLoc
         Vect1%D=XYZRot%D(:,I)
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot%D,Vect1%D,XYZRot%D(:,I))
       ENDDO
     ENDIF
     !
     IF(PRESENT(DoCopy_O)) THEN
       IF(DoCopy_O) THEN
         XYZ=XYZRot%D 
       ENDIF
     ENDIF
     !
     CALL Delete(XYZRot)
     CALL Delete(Rot)
     CALL Delete(Vect1)
     CALL Delete(Vect2)
     CALL Delete(Vect3)
   END SUBROUTINE CALC_XYZRot
!
!----------------------------------------------------------------------
!
   SUBROUTINE ThreeAtoms(XYZ,Top12,ThreeAt,Linearity,IntCs_O)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     TYPE(INTC),OPTIONAL          :: IntCs_O
     TYPE(DBL_VECT)               :: DistVect
     INTEGER                      :: I,J,NatmsLoc
     INTEGER                      :: At1,At2,At3,ThreeAt(1:3),ITop
     REAL(DOUBLE)                 :: D12,D13,D23,CX,CY,CZ
     REAL(DOUBLE)                 :: Dist,SumU,Dist12,Area2
     REAL(DOUBLE)                 :: Vect1(1:3),Vect2(1:3),Vect3(1:3)
     LOGICAL                      :: Linearity
     TYPE(DBL_RNK2)               :: XYZRot
     TYPE(INT_RNK2)               :: Top12 
     !
     NatmsLoc=SIZE(XYZ,2)
     CALL New(XYZRot,(/3,NatmsLoc/))
     CALL New(DistVect,NatmsLoc)
     DistVect%D=0
     XYZRot%D=XYZ
     !
     ! First, try to select three atoms from constraints
     !
     At1=0
     At2=0
     At3=0
     IF(PRESENT(IntCs_O)) THEN
       CALL ThreeConstr(IntCs_O,Top12,NatmsLoc,At1,At2,At3)
     ENDIF
     !
     ! Find three atoms, which form a large triangle, 
     ! if molecule is not linear
     !
     CALL CenterOfMass(CX,CY,CZ,XYZ_O=XYZRot%D,Move_O=.True.)
     !
     ! Find farthest point from COM
     !
     IF(At1==0) THEN
       Dist=-One
       DO I=1,NatmsLoc
         SumU=XYZRot%D(1,I)**2+XYZRot%D(2,I)**2+XYZRot%D(3,I)**2
         ITop=MAX(Top12%I(I,1),1)
         SumU=SumU*DBLE(ITop**2)
         IF(SumU>Dist) THEN
           Dist=SumU
           At1=I
         ENDIF
       ENDDO
     ENDIF
     !
     ! Find farthest point from At1
     !
     DO I=1,NatmsLoc
       DistVect%D(I)=SQRT((XYZRot%D(1,I)-XYZRot%D(1,At1))**2+&
                          (XYZRot%D(2,I)-XYZRot%D(2,At1))**2+&
                          (XYZRot%D(3,I)-XYZRot%D(3,At1))**2)
     ENDDO
     IF(At2==0) THEN
       Dist=-One
       DO I=1,NatmsLoc
         ITop=MAX(Top12%I(I,1),1)
         SumU=DBLE(ITop)*DistVect%D(I)
         IF(SumU>Dist) THEN
           Dist=SumU
           At2=I
         ENDIF
       ENDDO
     ENDIF
     !
     ! Find point farthest from At1 and At2
     !
     IF(At3==0) THEN
       D12=DistVect%D(At2)
       Dist=-1.D99
       DO I=1,NatmsLoc
         IF(I==At1.OR.I==At2) CYCLE
         D23=SQRT((XYZRot%D(1,I)-XYZRot%D(1,At2))**2+&
                  (XYZRot%D(2,I)-XYZRot%D(2,At2))**2+&
                  (XYZRot%D(3,I)-XYZRot%D(3,At2))**2)
         D13=DistVect%D(I)
         SumU=0.5D0*(D12+D13+D23)
         Area2=SumU*(SumU-D12)*(SumU-D13)*(SumU-D23)
         ITop=MAX(Top12%I(I,1),1)
         Area2=DBLE(ITop)*Area2
         IF(Area2>Dist) THEN
           Dist=Area2
           At3=I
         ENDIF
       ENDDO
     ENDIF
     !
     ! Check linearity
     !
     Linearity=.FALSE.
     IF(NatmsLoc==2) THEN
       At3=0
       Linearity=.TRUE.
     ELSE
     Vect1=XYZRot%D(:,At1)-XYZRot%D(:,At2)
     Vect2=XYZRot%D(:,At1)-XYZRot%D(:,At3)
     CALL CROSS_PRODUCT(Vect1,Vect2,Vect3)
     D12=SQRT(DOT_PRODUCT(Vect3,Vect3))
     IF(D12<0.001D0) Linearity=.TRUE.
     ENDIF
     ThreeAt(1)=At1
     ThreeAt(2)=At2
     ThreeAt(3)=At3
     ! Tidy up
     CALL Delete(XYZRot)
     CALL Delete(DistVect)
   END SUBROUTINE ThreeAtoms
!
!--------------------------------------------------------------------
!
   SUBROUTINE WriteChol(CholData,FileName)
     TYPE(Cholesky)   :: CholData
     CHARACTER(LEN=*) :: FileName
     INTEGER          :: ChNon0,NCart
     !
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     IF(.NOT.AllocQ(CholData%GcScale%Alloc)) THEN
       RETURN
     ENDIF
     NCart=SIZE(CholData%Perm%I)
     ChNon0=CholData%ChRowPt%I(NCart+1)-1
     WRITE(99) NCart 
     WRITE(99) ChNon0
     WRITE(99) CholData%Perm%I
     WRITE(99) CholData%IPerm%I
     WRITE(99) CholData%ChRowPt%I
     WRITE(99) CholData%ChColPt%I
     WRITE(99) CholData%GcScale%D
     WRITE(99) CholData%ChDiag%D
     WRITE(99) CholData%ChFact%D
     CLOSE(99,STATUS='KEEP')
   END SUBROUTINE WriteChol
!
!---------------------------------------------------------------
!
   SUBROUTINE ReadChol(CholData,FileName)
     TYPE(Cholesky)   :: CholData
     CHARACTER(LEN=*) :: FileName
     INTEGER          :: ChNon0,NCart
     LOGICAL          :: Exists
     !
     INQUIRE(File=FileName,EXIST=Exists)
     IF(.NOT.Exists) &
       CALL Halt('File does not exist error in ReadChol')
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     READ(99) NCart 
     READ(99) ChNon0
     CALL New_Chol(CholData,NCart,ChNon0)
     READ(99) CholData%Perm%I 
     READ(99) CholData%IPerm%I 
     READ(99) CholData%ChRowPt%I
     READ(99) CholData%ChColPt%I
     READ(99) CholData%GcScale%D
     READ(99) CholData%ChDiag%D
     READ(99) CholData%ChFact%D
     CLOSE(99,STATUS='KEEP')
   END SUBROUTINE ReadChol
!
!--------------------------------------------------------------------
!
   SUBROUTINE CleanBConstr(IntCs,B,NatmsLoc,Inv_O,NConstr_O)
     TYPE(BMATR)      :: B
     TYPE(INTC)       :: IntCs
     INTEGER          :: I,J,K,L,NIntC,JJ,LL,NatmsLoc,NCart
     TYPE(INT_VECT)   :: IConstr
     LOGICAL,OPTIONAL :: Inv_O
     LOGICAL          :: Inv
     INTEGER,OPTIONAL :: NConstr_O
     INTEGER          :: NConstr
     !
     IF(IntCs%N<=0) RETURN
     Inv=.FALSE.
     IF(PRESENT(Inv_O)) Inv=Inv_O
     !
     !Zero columns of B matrix that are related to Cartesian constraints
     !
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def%C)
     CALL New(IConstr,NCart)
     IConstr%I=1  
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint%L(I)) THEN
         J=3*(IntCs%Atoms%I(I,1)-1)
         IF(IntCs%Def%C(I)(1:5)=='CARTX') THEN
           NConstr=NConstr+1
           IConstr%I(J+1)=0 
         ENDIF
         IF(IntCs%Def%C(I)(1:5)=='CARTY') THEN
           NConstr=NConstr+1
           IConstr%I(J+2)=0
         ENDIF
         IF(IntCs%Def%C(I)(1:5)=='CARTZ') THEN
           NConstr=NConstr+1
           IConstr%I(J+3)=0
         ENDIF
       ENDIF
     ENDDO 
     IF(PRESENT(NConstr_O)) NConstr_O=NConstr
     !
     DO I=1,NIntC
       DO J=1,4
         JJ=B%IB%I(I,J)
         IF(JJ==0) EXIT
         JJ=3*(JJ-1)
         LL=3*(J-1)
         IF(Inv) THEN
           DO L=1,3 ; IF(IConstr%I(JJ+L)/=0) B%B%D(I,LL+L)=Zero ; ENDDO
         ELSE
           DO L=1,3 ; IF(IConstr%I(JJ+L)==0) B%B%D(I,LL+L)=Zero ; ENDDO
         ENDIF
       ENDDO
     ENDDO 
     CALL Delete(IConstr)
   END SUBROUTINE CleanBConstr
!
!--------------------------------------------------------------------
!
   SUBROUTINE CleanB(ThreeAt,B,NIntC)
     TYPE(BMATR) :: B
     INTEGER     :: I,J,K,L,NIntC,J1,J2,ThreeAt(1:3)
     !
     DO I=1,NIntC
       DO J=1,4
         IF(B%IB%I(I,J)==ThreeAt(1)) THEN
           J2=3*J
           J1=J2-2
           B%B%D(I,J1:J2)=Zero
         ENDIF
         IF(B%IB%I(I,J)==ThreeAt(2)) THEN
           J2=3*J
           J1=J2-1
           B%B%D(I,J1:J2)=Zero
         ENDIF
         IF(B%IB%I(I,J)==ThreeAt(3)) THEN
           J2=3*J
           B%B%D(I,J2)=Zero
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE CleanB
!
!-------------------------------------------------------------------
!
   SUBROUTINE PrtXYZ(AtNum,XYZ,FileName,Title,Vects_O,XYZL_O)
     INTEGER,DIMENSION(:)          :: AtNum
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL   :: Vects_O,XYZL_O
     CHARACTER(LEN=*)              :: FileName
     CHARACTER(LEN=*)              :: Title
     CHARACTER(LEN=1)              :: CharU 
     INTEGER                       :: I,NatmsLoc,II
     !
     NatmsLoc=SIZE(XYZ,2)
     OPEN(File=FileName,Unit=99,FORM='FORMATTED',STATUS='UNKNOWN')
     REWIND(99)
     II=0    
     DO 
       READ(99,33,END=1) CharU
       II=II+1
     ENDDO
     1    CONTINUE
     33 format(a1)
     !
     REWIND(99)
     DO I=1,II
       READ(99,33) CharU 
     ENDDO
     IF(PRESENT(XYZL_O)) THEN
       WRITE(99,*) NatmsLoc+3+1
     ELSE
       WRITE(99,*) NatmsLoc 
     ENDIF
     WRITE(99,*) TRIM(Title)
     IF(PRESENT(Vects_O)) THEN
       DO I=1,NatmsLoc
         WRITE(99,100) AtNum(I),XYZ(1:3,I)/AngstromsToAu,100.D0*Vects_O(1:3,I)
       ENDDO 
     ELSE
       DO I=1,NatmsLoc
         WRITE(99,100) AtNum(I),XYZ(1:3,I)/AngstromsToAu
       ENDDO 
     ENDIF
     !
     IF(PRESENT(XYZL_O)) THEN
       NatmsLoc=SIZE(XYZL_O,2)-3
         WRITE(99,200) Zero,Zero,Zero
       DO I=1,3         
         WRITE(99,200) XYZL_O(1:3,NatmsLoc+I)/AngstromsToAu
       ENDDO 
     ENDIF
     200  FORMAT('  I ',2X,3F20.10,2X,3F20.10)
     100  FORMAT(I4,2X,3F20.10,2X,3F20.10)
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE PrtXYZ
!
!---------------------------------------------------------------------
!
   SUBROUTINE BoxBorders(XYZ,BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: BXMIN,BXMax
     REAL(DOUBLE)                :: BYMIN,BYMax
     REAL(DOUBLE)                :: BZMIN,BZMax
     REAL(DOUBLE)                :: VBig
     INTEGER                     :: NatmsLoc,I
     !
     !find borders of the global Box
     !
     NatmsLoc=SIZE(XYZ,2)
     VBig=1.D99
     BXMIN= VBIG
     BXMax=-VBIG
     BYMIN= VBIG
     BYMax=-VBIG
     BZMIN= VBIG
     BZMax=-VBIG
     DO I=1,NatmsLoc
       IF(XYZ(1,I)<BXMIN) BXMIN=XYZ(1,I)
       IF(XYZ(1,I)>BXMax) BXMax=XYZ(1,I)
       IF(XYZ(2,I)<BYMIN) BYMIN=XYZ(2,I)
       IF(XYZ(2,I)>BYMax) BYMax=XYZ(2,I)
       IF(XYZ(3,I)<BZMIN) BZMIN=XYZ(3,I)
       IF(XYZ(3,I)>BZMax) BZMax=XYZ(3,I)
     ENDDO
   END SUBROUTINE BoxBorders
! 
!-------------------------------------------------------------------
!
   SUBROUTINE ReadBMATR(ISpB,JSpB,ASpB,FileName, &
                        IPerm1_O,IPerm2_O,UMatr_O,B_O)
     INTEGER          :: NDim,NZ,NCart,NDim1,NDim2
     TYPE(INT_VECT)   :: ISpB,JSpB
     TYPE(DBL_VECT)   :: ASpB
     TYPE(INT_VECT),OPTIONAL   :: IPerm1_O,IPerm2_O
     TYPE(DBL_RNK2),OPTIONAL   :: UMatr_O
     TYPE(BMATR),OPTIONAL      :: B_O
     LOGICAL          :: Exists
     CHARACTER(LEN=*) :: FileName
     CHARACTER(LEN=8) :: Aux8
     !
     INQUIRE(File=FileName,EXIST=Exists)
     IF(.NOT.Exists) &
       CALL Halt('File does not exist error in ReadBMATR for '// &
                  FileName)
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     READ(99) NDim,NZ
     CALL New(ISpB,NDim+1)
     CALL New(JSpB,NZ)
     CALL New(ASpB,NZ)
     READ(99) ISpB%I
     READ(99) JSpB%I
     READ(99) ASpB%D
     IF(PRESENT(IPerm1_O)) THEN
       READ(99) NCart
       CALL New(IPerm1_O,NCart)
       READ(99) IPerm1_O%I
     ENDIF
     IF(PRESENT(IPerm2_O)) THEN
       READ(99) NCart
       CALL New(IPerm2_O,NCart)
       READ(99) IPerm2_O%I
     ENDIF
     IF(PRESENT(UMatr_O)) THEN
       READ(99) NDim1
       READ(99) NDim2
       CALL New(UMatr_O,(/NDim1,NDim2/))
       READ(99) UMatr_O%D
     ENDIF
     IF(PRESENT(B_O)) THEN
       READ(99) Aux8
       IF(Aux8=='HAVEBMAT') THEN
         READ(99) NDim1
         CALL New(B_O,NDim1)
         READ(99) B_O%IB%I
         READ(99) B_O%B%D
         READ(99) B_O%BLI%I
         READ(99) B_O%BL%D
       ENDIF
     ENDIF
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE ReadBMATR
! 
!-------------------------------------------------------------------
!
   SUBROUTINE WriteBMATR(ISpB,JSpB,ASpB,FileName, &
                         IPerm1_O,IPerm2_O,UMatr_O,B_O)
     INTEGER          :: NDim,NZ
     TYPE(INT_VECT)   :: ISpB,JSpB
     TYPE(INT_VECT),OPTIONAL   :: IPerm1_O,IPerm2_O
     TYPE(DBL_RNK2),OPTIONAL   :: UMatr_O
     TYPE(BMATR),OPTIONAL      :: B_O
     TYPE(DBL_VECT)   :: ASpB
     CHARACTER(LEN=*) :: FileName
     !
     NDim=SIZE(ISpB%I,1)-1
     NZ=SIZE(JSpB%I)
     OPEN(File=FileName,Unit=99,FORM='UNFORMATTED',STATUS='UNKNOWN')
     WRITE(99) NDim,NZ
     WRITE(99) ISpB%I
     WRITE(99) JSpB%I
     WRITE(99) ASpB%D
     IF(PRESENT(IPerm1_O)) THEN
       WRITE(99) SIZE(IPerm1_O%I)
       WRITE(99) IPerm1_O%I
     ENDIF
     IF(PRESENT(IPerm2_O)) THEN
       WRITE(99) SIZE(IPerm2_O%I)
       WRITE(99) IPerm2_O%I
     ENDIF
     IF(PRESENT(UMatr_O)) THEN
       WRITE(99) SIZE(UMatr_O%D,1)
       WRITE(99) SIZE(UMatr_O%D,2)
       WRITE(99) UMatr_O%D 
     ENDIF
     IF(PRESENT(B_O)) THEN
       WRITE(99) 'HAVEBMAT'
       WRITE(99) SIZE(B_O%IB%I,1)
       WRITE(99) B_O%IB%I
       WRITE(99) B_O%B%D
       WRITE(99) B_O%BLI%I
       WRITE(99) B_O%BL%D
     ELSE
       WRITE(99) 'XXNOBMAT'
     ENDIF
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE WriteBMATR
! 
!-------------------------------------------------------------------
!
   SUBROUTINE LongRangeIntC(CtrlCoord,Top12,IntCs,XYZ,AtNum,IEq)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER                     :: NIntC,NIntc_New
     TYPE(INT_VECT)              :: LinAtom,LinCenter
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,Conv,Dist12,Angle123,Angle234
     REAL(DOUBLE)                :: Angle123S,Angle234S,PiHalf
     INTEGER                     :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
     INTEGER                     :: II1,II4
     INTEGER                     :: I,J,K,NatmsLoc,III
     INTEGER                     :: NStreLinb,NBendLinb,JJ,IC,IEx
     TYPE(INT_RNK2)              :: LinBBridge,Top12
     TYPE(CoordCtrl)             :: CtrlCoord
     LOGICAL                     :: RepeatChk,DoVdW,Found,AllTors
     INTEGER,DIMENSION(:)        :: AtNum,IEq
     !
     AllTors=.FALSE.
     Conv=180.D0/PI
     PiHalf=Half*PI
     NatmsLoc=SIZE(Top12%I,1)
     !
     NMax12=SIZE(Top12%I,2)-1
     !
     CALL NEW(LinBBridge,(/2,IntCs%N/))
     CALL NEW(LinAtom,NatmsLoc)
     CALL NEW(LinCenter,IntCs%N)
     LinBBridge%I=0
     LinAtom%I=0
     LinCenter%I=0
     NLinB=0
     DO I=1,IntCs%N
       IF(IntCs%Def%C(I)(1:5)=='LINB1') THEN   
        !IF(LinAtom%I(IntCs%Atoms%I(I,1))/=0) CYCLE
         NLinB=NLinB+1 
         I1=IntCs%Atoms%I(I,1)
         I2=IntCs%Atoms%I(I,2)
         I3=IntCs%Atoms%I(I,3)
         LinAtom%I(I1)=1
         LinAtom%I(I2)=1
         LinAtom%I(I3)=1
         LinCenter%I(NLinB)=I2
         LinBBridge%I(1,NLinB)=I1
         LinBBridge%I(2,NLinB)=I3
         ! now, go on left side, use Top12, then 
         ! go on right, then define torsions
         I1=LinBBridge%I(2,NLinB)
         I2=LinBBridge%I(1,NLinB)
        !
        ! WARNING! The few lines below that have been commented out
        ! refere to the treatment of molecules like H2C=C=C=C=CH2
        ! they are supposed to find very long range torsions
        !
        !DO III=1,IntCs%N
        !  RepeatChk=.FALSE.
        !  DO J=1,Top12%I(I2,1)
        !    I3=Top12%I(I2,1+J) 
        !    CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3), &
        !              Value_O=Value)
        !    IF(ABS(Value-PI)*Conv<CtrlCoord%LinCrit) THEN 
        !      LinBBridge%I(1,NLinB)=I3
        !      LinAtom%I(I3)=1
        !      I1=I2
        !      I2=I3
        !      RepeatChk=.TRUE.
        !      EXIT
        !    ENDIF
        !  ENDDO
        !  IF(.NOT.RepeatChk) EXIT
        !ENDDO
        !!
        !I1=LinBBridge%I(1,NLinB)
        !I2=LinBBridge%I(2,NLinB)
        !DO III=1,IntCs%N
        !  RepeatChk=.FALSE.
        !  DO J=1,Top12%I(I2,1)
        !    I3=Top12%I(I2,1+J) 
        !    CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3), &
        !              Value_O=Value)
        !    IF(ABS(Value-PI)*Conv<CtrlCoord%LinCrit) THEN 
        !      LinBBridge%I(2,NLinB)=I3
        !      LinAtom%I(I3)=1
        !      I1=I2
        !      I2=I3
        !      RepeatChk=.TRUE.
        !      EXIT
        !    ENDIF
        !  ENDDO
        !  IF(.NOT.RepeatChk) EXIT
        !ENDDO
       ENDIF
     ENDDO 
     !
     ! bridges are set now, add torsions, stretches and angles
     ! first, count number of torsions, to be added
     !
     NTorsLinB=0
     NStreLinB=0
     NBendLinB=0
     DO I=1,NLinB
       I1=LinBBridge%I(1,I)
       I2=LinBBridge%I(2,I)
       NTorsLinB=NTorsLinB+Top12%I(I1,1)*Top12%I(I2,1)
       NStreLinB=NStreLinB+1
       NBendLinB=NBendLinB+Top12%I(I1,1)+Top12%I(I2,1)
     ENDDO
     !
     ! Now generate the INTCs for the new torsions
     !
     NIntc_New=IntCs%N+NTorsLinB+NStreLinB+NBendLinB
     CALL New(IntC_New,NIntc_New)
     NIntC=IntCs%N
     CALL Set_INTC_EQ_INTC(IntCs,IntC_New,1,NIntC,1)
     !
     ! torsions
     !
     DO I=1,NLinB
       Found=.FALSE.
       II1=0
       II4=0
       I2=LinBBridge%I(1,I)
       I3=LinBBridge%I(2,I)
       Angle123S=Zero
       Angle234S=Zero
       DO J=1,Top12%I(I2,1)
         I1=Top12%I(I2,1+J)
         IF(LinCenter%I(I)==I1) CYCLE
         CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),Value_O=Angle123)
         DO K=1,Top12%I(I3,1)
           I4=Top12%I(I3,1+K)
           IF(LinCenter%I(I)==I4) CYCLE
           CALL BEND(XYZ(1:3,I2),XYZ(1:3,I3),XYZ(1:3,I4),Value_O=Angle234)
           IF(I1/=I4.AND.I1/=I3.AND.I2/=I4) THEN
             IF(.NOT.Found) THEN
               II1=I1
               II4=I4
               Found=.TRUE.
             ELSE
               IF(ABS(PiHalf-Angle123)<ABS(PiHalf-Angle123S)) THEN
                 Angle123S=Angle123
                 II1=I1
               ENDIF
               IF(ABS(PiHalf-Angle234)<ABS(PiHalf-Angle234S)) THEN
                 Angle234S=Angle234
                 II4=I4
               ENDIF
              !IF(AtNum(I1)>AtNum(II1)) II1=I1
              !IF(AtNum(I4)>AtNum(II4)) II4=I4
             ENDIF
             !
             IF(AllTors) THEN
               NIntC=NIntC+1
               IntC_New%Def%C(NIntC)(1:10)='TORSL     '
               IntC_New%Atoms%I(NIntC,1)=I1
               IntC_New%Atoms%I(NIntC,2)=I2
               IntC_New%Atoms%I(NIntC,3)=I3
               IntC_New%Atoms%I(NIntC,4)=I4
               IntC_New%Value%D(NIntC)=Zero
               IntC_New%Constraint%L(NIntC)=.FALSE.
               IntC_New%ConstrValue%D(NIntC)=Zero   
               IntC_New%Active%L(NIntC)=.TRUE. 
             ENDIF
           ENDIF
         ENDDO
       ENDDO
       IF(.NOT.AllTors) THEN
         I1=II1
         I4=II4
         IF(Found) THEN
           NIntC=NIntC+1
           IntC_New%Def%C(NIntC)(1:10)='TORSL     '
           IntC_New%Atoms%I(NIntC,1)=I1
           IntC_New%Atoms%I(NIntC,2)=I2
           IntC_New%Atoms%I(NIntC,3)=I3
           IntC_New%Atoms%I(NIntC,4)=I4
           IntC_New%Value%D(NIntC)=Zero
           IntC_New%Constraint%L(NIntC)=.FALSE.
           IntC_New%ConstrValue%D(NIntC)=Zero   
           IntC_New%Active%L(NIntC)=.TRUE. 
         ENDIF
       ENDIF
     ENDDO
     !
     CALL Delete(Intcs)
     CALL New(Intcs,NIntC)
       CALL Set_INTC_EQ_INTC(IntC_New,IntCs,1,NIntC,1)
     CALL Delete(IntC_New)
     !
     CALL Delete(LinCenter)
     CALL Delete(LinAtom)
     CALL Delete(LinBBridge)
   END SUBROUTINE LongRangeIntC
!
!----------------------------------------------------------------------
!
   SUBROUTINE RedundancyOff(Displ,SCRPath,Print,Messg_O)
     REAL(DOUBLE),DIMENSION(:) :: Displ
     REAL(DOUBLE)              :: Perc 
     TYPE(INT_VECT)            :: ISpB,JSpB,IGc,JGc
     TYPE(DBL_VECT)            :: ASpB,AGc
     TYPE(Cholesky)            :: CholData
     INTEGER                   :: NIntC,NCart,I
     TYPE(DBL_VECT)            :: Vect1,Displ2
     INTEGER                   :: Print
     CHARACTER(LEN=*)          :: SCRPath
     CHARACTER(LEN=*),OPTIONAL :: Messg_O
     CHARACTER(LEN=DCL)        :: Messg
     !
write(*,*) 'RedundancyOff hardwired to return'
return
     CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)  
     NIntC=SIZE(ISpB%I)-1
     IF(NIntC/=SIZE(Displ)) &
       CALL Halt('Dimension error in RedundancyOff')
     NCart=SIZE(CholData%IPerm%I)
     CALL New(Vect1,NCart)
     CALL New(Displ2,NIntC)
     Displ2%D=Displ
     !
     CALL CALC_BxVect(ISpB,JSpB,ASpB,Displ,Vect1%D,Trp_O=.TRUE.)
     CALL GcInvIter(Vect1%D,ISpB,JSpB,ASpB,CholData,NIntC)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,Displ,Vect1%D)
     !
     Perc=DOT_PRODUCT(Displ,Displ2%D)/DOT_PRODUCT(Displ2%D,Displ2%D) 
     Perc=(One-ABS(Perc))*100.D0
     IF(PRESENT(Messg_O)) THEN
       Messg=TRIM(Messg_O)// &
         " Percentage of Redundancy projected out= " &
         //TRIM(IntToChar(INT(Perc)))
     ELSE
       Messg= &
         " Percentage of Redundancy projected out= " &
         //TRIM(IntToChar(INT(Perc)))
     ENDIF
     IF(Print>=DEBUG_GEOP_MAX) THEN
       WRITE(*,*) TRIM(Messg)
       WRITE(Out,*) TRIM(Messg)
     ENDIF
100  FORMAT(F8.2)
     !
     CALL Delete(Displ2)
     CALL Delete(Vect1)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(CholData)
     ! 
   END SUBROUTINE RedundancyOff
! 
!-------------------------------------------------------------------
!
   SUBROUTINE PrimToDelocOld(VectInt,VectCart,ISpB,JSpB,ASpB,CholData)
     REAL(DOUBLE),DIMENSION(:) :: VectInt,VectCart
     TYPE(INT_VECT)            :: ISpB,JSpB
     TYPE(DBL_VECT)            :: ASpB
     TYPE(Cholesky)            :: CholData
     INTEGER                   :: I,J,NCart,NIntC
     TYPE(DBL_VECT)            :: VectCartAux
     !
     NCart=SIZE(VectCart)
     NIntC=SIZE(VectInt)
     CALL New(VectCartAux,NCart)
     !
     CALL CALC_BxVect(ISpB,JSpB,ASpB,VectInt,VectCartAux%D,Trp_O=.TRUE.)
     !
     CALL PermVect(VectCartAux%D,VectCart,CholData%Perm%I)
     CALL ScaleVect(VectCart,CholData%GcScale%D)
     !
     CALL UtInvX(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,VectCart,NCart,VectCartAux%D)
     !
     DO I=1,NCart
       VectCartAux%D(I)=SQRT(CholData%ChDiag%D(I))*VectCartAux%D(I)
     ENDDO
     !
     VectCart=VectCartAux%D
     !
     CALL Delete(VectCartAux)
   END SUBROUTINE PrimToDelocOld
! 
!-------------------------------------------------------------------
!
   SUBROUTINE DelocToPrimOld(NewDelocs,NewPrims,ISpB,JSpB,ASpB,CholData)
     REAL(DOUBLE),DIMENSION(:) :: NewDelocs,NewPrims
     INTEGER                   :: NCart,NIntC
     TYPE(Cholesky)            :: CholData
     TYPE(DBL_VECT)            :: VectCart,VectInt
     INTEGER                   :: I
     TYPE(INT_VECT)            :: ISpB,JSpB
     TYPE(DBL_VECT)            :: ASpB
     !
     NCart=SIZE(NewDelocs)
     NIntC=SIZE(NewPrims)
     CALL New(VectCart,NCart)
     !
     DO I=1,NCart
       NewDelocs(I)=SQRT(CholData%ChDiag%D(I))*NewDelocs(I)
     ENDDO
     !
     CALL UInvX(CholData%ChRowPt%I,CholData%ChColPt%I, &
       CholData%ChFact%D,NewDelocs,NCart,VectCart%D)
     !
     CALL ScaleVect(VectCart%D,CholData%GcScale%D)
     CALL PermVect(VectCart%D,NewDelocs,CholData%IPerm%I)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,NewPrims,NewDelocs)
     !
     CALL Delete(VectCart)
   END SUBROUTINE DelocToPrimOld
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasAngle(CharU)
     CHARACTER(LEN=*) :: CharU
     HasAngle=(CharU(1:4)=='BEND'.OR. &
               CharU(1:4)=='LINB'.OR. &
               CharU(1:4)=='OUTP'.OR. &
               CharU(1:4)=='TORS'.OR. &
               CharU(1:5)=='ALPHA'.OR. &
               CharU(1:4)=='BETA'.OR. &
               CharU(1:5)=='GAMMA')
   END FUNCTION HasAngle
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasBendLinB(CharU)
     CHARACTER(LEN=*) :: CharU
     HasBendLinB=(CharU(1:4)=='BEND'.OR. &
                  CharU(1:4)=='LINB')
   END FUNCTION HasBendLinB
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasTorsOutP(CharU)
     CHARACTER(LEN=*) :: CharU
     HasTorsOutP=(CharU(1:4)=='OUTP'.OR. &
                  CharU(1:4)=='TORS')
   END FUNCTION HasTorsOutP
!
!-------------------------------------------------------------------
!
   SUBROUTINE DistMatr(ITop,JTop,ATop,IntCs,NatmsLoc,NStre)
     TYPE(INT_VECT)       :: ITop,JTop,NAux
     TYPE(DBL_VECT)       :: ATop
     TYPE(INTC)           :: IntCs
     INTEGER              :: NatmsLoc,I,J,NBond,NStre,I1,I2,NIntC
     INTEGER              :: I1off,I2Off
     !
     NIntC=SIZE(IntCs%Def%C)
     CALL New(ITop,NatmsLoc+1)
     CALL New(JTop,2*NStre)
     CALL New(ATop,2*NStre)
     CALL New(NAux,NatmsLoc)
     NBond=0
     NAux%I=0
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         I1=IntCs%Atoms%I(I,1)
         I2=IntCs%Atoms%I(I,2)
         NAux%I(I1)=NAux%I(I1)+1
         NAux%I(I2)=NAux%I(I2)+1
       ENDIF 
     ENDDO
     !
     ITop%I=0
     ITop%I(1)=1
     DO I=1,NatmsLoc
       ITop%I(I+1)=ITop%I(I)+NAux%I(I)
     ENDDO
     !
     NAux%I=0
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         I1=IntCs%Atoms%I(I,1)
         I2=IntCs%Atoms%I(I,2)
         NAux%I(I1)=NAux%I(I1)+1
         NAux%I(I2)=NAux%I(I2)+1
         I1off=ITop%I(I1)+NAux%I(I1)-1
         I2off=ITop%I(I2)+NAux%I(I2)-1
         JTop%I(I1off)=I2
         ATop%D(I1off)=IntCs%Value%D(I)
         JTop%I(I2off)=I1
         ATop%D(I2off)=IntCs%Value%D(I)
       ENDIF 
     ENDDO
     !
     CALL Delete(NAux)
   END SUBROUTINE DistMatr
!
!------------------------------------------------------------------
!
   FUNCTION HasLonelyAtm(Top12,JJ1,JJ2,LAtm)
     LOGICAL HasLonelyAtm
     TYPE(INT_RNK2) :: Top12
     INTEGER        :: LAtm,JJ1,JJ2
     !
     IF(Top12%I(JJ1,1)==0) THEN
       HasLonelyAtm=.TRUE.
       LAtm=JJ1
     ELSE IF(Top12%I(JJ2,1)==0) THEN
       HasLonelyAtm=.TRUE.
       LAtm=JJ2
     ELSE
       HasLonelyAtm=.FALSE.
       LAtm=0  
     ENDIF
   END FUNCTION HasLonelyAtm
!
!------------------------------------------------------------------
!
   FUNCTION HasHBond(Top12,AtNum,NJJ1,NJJ2,JJ1,JJ2,HAtm,XYZ,LinCrit) 
     LOGICAL                     :: HasHBond,HasAtt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: LinCrit
     TYPE(INT_RNK2)              :: Top12
     INTEGER,DIMENSION(:)        :: AtNum
     INTEGER                     :: NJJ1,NJJ2,JJ1,JJ2,HAtm,JJE
     !
     HasHBond=.FALSE.
     HAtm=0
     IF(NJJ1/=1.AND.NJJ2/=1) RETURN
     IF((NJJ1==1.AND.HasLigand(NJJ2))) THEN
       HasHBond=HasAttached(AtNum,Top12%I,JJ1,JJ2,JJE,XYZ,LinCrit)
      !HasHBond=.TRUE.
       HAtm=JJ1
     ELSE IF((NJJ2==1.AND.HasLigand(NJJ1))) THEN
       HasHBond=HasAttached(AtNum,Top12%I,JJ2,JJ1,JJE,XYZ,LinCrit)
      !HasHBond=.TRUE.
       HAtm=JJ2
     ENDIF
   END FUNCTION HasHBond
!
!---------------------------------------------------------------------
!
   FUNCTION HasAttached(AtNum,Top12,JJ1,JJ2,JJE,XYZ,LinCrit)
     LOGICAL                      :: HasAttached
     INTEGER,DIMENSION(:,:)       :: Top12
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     INTEGER,DIMENSION(:)         :: AtNum
     INTEGER                      :: JJ1,J,K,JJ2,JJE
     REAL(DOUBLE)                 :: Value,Conv,LinCrit
     !
     ! For H-bonds: A-H...B = JJE-JJ1...JJ2
     !
     HasAttached=.FALSE.
     Conv=180.D0/PI
     JJE=0
     IF(Top12(JJ1,1)==0) THEN
       HasAttached=.TRUE.
       RETURN
     ENDIF
     DO J=1,Top12(JJ1,1)
       JJE=Top12(JJ1,J+1)
       K=AtNum(JJE)
       IF(HasLigand(K)) THEN
           HasAttached=.TRUE.
       ELSE
         CALL BEND(XYZ(1:3,JJE),XYZ(1:3,JJ1),XYZ(1:3,JJ2),Value_O=Value)
         IF((ABS(Value-PI)*Conv < LinCrit)) THEN
           HasAttached=.TRUE.
           EXIT  
         ENDIF
       ENDIF
     ENDDO 
   END FUNCTION HasAttached
!
!---------------------------------------------------------------------
!
   SUBROUTINE SortNonCov2(AtNum,XYZ,Bond,AtmB)
     TYPE(BONDDATA)              :: Bond,BondNew
     TYPE(ATOMBONDS)             :: AtmB
     INTEGER                     :: NatmsLoc
     TYPE(INT_VECT)              :: AllowBond,BondScore
     TYPE(DBL_RNK2)              :: BondVect2
     TYPE(DBL_VECT)              :: BondVect3,Ordered
     INTEGER,DIMENSION(:)        :: AtNum
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ  
     REAL(DOUBLE)                :: Dist,DistK
     REAL(DOUBLE)                :: CondNumb1,CondNumb2
     INTEGER                     :: NBondNew,MaxBonds
     INTEGER                     :: I,J,K,L,I1,I2,NDim,M1,M2,M
     !
     MaxBonds=30
     CondNumb1=0.5D0
     CondNumb2=0.8D0
     NatmsLoc=SIZE(XYZ,2)
     !
     CALL SetDSYEVWork(3)
     CALL New(AllowBond,Bond%N)
     CALL New(BondScore,Bond%N)
     CALL ScoreBond(BondScore%I,Bond)
     !
     AllowBond%I=1
     DO I1=1,NatmsLoc
       NDim=AtmB%Count%I(I1)
       IF(NDim<2) CYCLE
       CALL New(BondVect2,(/NDim,3/))
       CALL New(BondVect3,NDim)
       CALL New(Ordered,NDim)
       CALL D3Bonds(I1,XYZ,BondVect2%D,BondScore%I,AtmB,Bond)
       BondVect3%D=Zero
       DO J=1,3
         DO K=1,NDim 
           L=AtmB%Bonds%I(I1,K)
           IF(AllowBond%I(L)/=0) THEN
             BondVect3%D(K)=BondVect3%D(K)+BondVect2%D(K,J)**2
           ENDIF
         ENDDO
       ENDDO
       Dist=One/SUM(BondVect3%D)
       BondVect3%D=Dist*BondVect3%D
       DO K=1,NDim 
         Ordered%D(K)=DBLE(K)+0.1D0
       ENDDO
       CALL Reorder(BondVect3%D,Ordered%D)
       Dist=Zero
       DO J=NDim,1,-1
         K=INT(Ordered%D(J))
         L=AtmB%Bonds%I(I1,K)
         IF(AllowBond%I(L)==0) CYCLE
         M1=Bond%IJ%I(1,L)
         M2=Bond%IJ%I(2,L)
         M=MIN(AtmB%Count%I(M1),AtmB%Count%I(M2))
         IF(Dist>CondNumb2.AND.M>1) THEN 
       ! IF(Dist>CondNumb2.AND.M>1.AND. &
       !    Bond%Type%C(L)(1:4)/='Frag') THEN
           AllowBond%I(L)=0
           AtmB%Count%I(M1)=AtmB%Count%I(M1)-1
           AtmB%Count%I(M2)=AtmB%Count%I(M2)-1
         ENDIF
         Dist=Dist+BondVect3%D(J)
       ENDDO
       CALL Delete(Ordered)
       CALL Delete(BondVect3)
       CALL Delete(BondVect2)
     ENDDO
     !
     NBondNew=SUM(AllowBond%I)
     CALL New(BondNew,NBondNew)
     NBondNew=0
     DO I=1,Bond%N
       IF(AllowBond%I(I)==1) THEN
         NBondNew=NBondNew+1
         CALL Set_Bond_EQ_Bond(BondNew,NBondNew,Bond,I)
       ENDIF
     ENDDO
     !
     CALL Delete(Bond)
     CALL Set_BONDDATA_EQ_BONDDATA(Bond,BondNew)
     CALL Delete(BondNew)
     CALL Delete(AllowBond)
     CALL Delete(BondScore)
     CALL UnSetDSYEVWork()
   END SUBROUTINE SortNonCov2
!
!----------------------------------------------------------------------
!
   SUBROUTINE D3Bonds(I1,XYZ,BondVect2,BondScore,AtmB,Bond,Eig_O)
     TYPE(BONDDATA)                     :: Bond
     TYPE(ATOMBONDS)                    :: AtmB
     REAL(DOUBLE),DIMENSION(:,:)        :: XYZ,BondVect2
     INTEGER,DIMENSION(:)               :: BondScore     
     TYPE(DBL_RNK2)                     :: BondVect
     REAL(DOUBLE),DIMENSION(3,3)        :: Matr,EigVects
     REAL(DOUBLE),DIMENSION(3)          :: EigVals    
     REAL(DOUBLE)                       :: Dist
     REAL(DOUBLE),DIMENSION(3),OPTIONAL :: Eig_O  
     INTEGER                            :: NDim,I,J,K,I1,II1,II2,Info
     !
     NDim=AtmB%Count%I(I1)
     CALL New(BondVect,(/NDim,3/))
     DO J=1,NDim
       K=AtmB%Bonds%I(I1,J)
       II1=Bond%IJ%I(1,K)
       II2=Bond%IJ%I(2,K)
       Dist=DBLE(BondScore(K))/Bond%Length%D(K)**3
       BondVect%D(J,1:3)=Dist*(XYZ(1:3,II1)-XYZ(1:3,II2))
     ENDDO
     !
     CALL DGEMM_TNc(3,NDim,3,One,Zero,BondVect%D,BondVect%D,Matr)
     BLKVECT%D=Matr
     CALL DSYEV('V','U',3,BLKVECT%D,BIGBLOK,BLKVALS%D, &
     BLKWORK%D,BLKLWORK,INFO)
     IF(INFO/=SUCCEED) &
     CALL Halt('DSYEV hosed in D3Bonds. INFO='&
                //TRIM(IntToChar(INFO)))
     EigVects=BLKVECT%D
     EigVals=BLKVALS%D
     !
     CALL DGEMM_NNc(NDim,3,3,One,Zero,BondVect%D,EigVects,BondVect2)
     !
     DO J=1,3 ; EigVals(J)=ABS(EigVals(J)) ; ENDDO
     EigVals=EigVals/MAXVAL(EigVals) 
     IF(PRESENT(Eig_O)) THEN
       Eig_O=EigVals
     ENDIF
     ! 
     ! Weight directions
     ! 
     DO J=1,3
       Dist=One/SQRT(ABS(EigVals(J)))
       DO K=1,NDim 
         BondVect2(K,J)=Dist*BondVect2(K,J) 
       ENDDO
     ENDDO
     !
     CALL Delete(BondVect)
   END SUBROUTINE D3Bonds
!
!---------------------------------------------------------------------
!
   SUBROUTINE TranslToAt1(VectCart,ThreeAt,Vect_O)
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: Vect_O
     REAL(DOUBLE),DIMENSION(:)          :: VectCart
     REAL(DOUBLE)                       :: V1,V2,V3
     INTEGER,DIMENSION(3)               :: ThreeAt
     INTEGER                            :: I,J,NCart,NatmsLoc,Istart
     !
     NCart=SIZE(VectCart)
     NatmsLoc=NCart/3
     IF(NCart/=3*NatmsLoc) CALL Halt('Dimension error in TranslToAt1')
     !
     IStart=3*(ThreeAt(1)-1)
     IF(PRESENT(Vect_O)) THEN
       V1=Vect_O(1)
       V2=Vect_O(2)
       V3=Vect_O(3)
     ELSE
       V1=VectCart(IStart+1)
       V2=VectCart(IStart+2)
       V3=VectCart(IStart+3)
     ENDIF
     !
     DO I=1,NatmsLoc
       J=3*(I-1)  
       VectCart(J+1)=VectCart(J+1)-V1
       VectCart(J+2)=VectCart(J+2)-V2
       VectCart(J+3)=VectCart(J+3)-V3
     ENDDO
   END SUBROUTINE TranslToAt1
!
!------------------------------------------------------------------
!
   SUBROUTINE RotToAt(VectCart,Rot,Rev_O)
     REAL(DOUBLE),DIMENSION(:)  :: VectCart
     REAL(DOUBLE),DIMENSION(3,3):: Rot     
     INTEGER                    :: I,J,NCart,NatmsLoc
     TYPE(DBL_VECT)             :: Vect,Vect2
     LOGICAL,OPTIONAL           :: Rev_O
     LOGICAL                    :: Reverse
     !
     NCart=SIZE(VectCart)
     NatmsLoc=NCart/3
     CALL New(Vect,3)
     CALL New(Vect2,3)
     Reverse=.FALSE.
     IF(PRESENT(Rev_O)) THEN
       Reverse=Rev_O
     ENDIF
     !
     DO I=1,NatmsLoc
       J=3*(I-1)
       Vect%D(1:3)=VectCart(J+1:J+3)
       IF(Reverse) THEN
         CALL DGEMM_TNc(3,3,1,One,Zero,Rot,Vect%D,Vect2%D)
       ELSE
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot,Vect%D,Vect2%D)
       ENDIF
       VectCart(J+1:J+3)=Vect2%D
     ENDDO 
     CALL Delete(Vect2)
     CALL Delete(Vect)
   END SUBROUTINE RotToAt
!
!------------------------------------------------------------------
!
   SUBROUTINE ZeroTrRots(VectCart,ThreeAt)
     REAL(DOUBLE),DIMENSION(:) :: VectCart
     INTEGER,DIMENSION(3)      :: ThreeAt
     INTEGER                   :: I,J
     !
     I=3*(ThreeAt(1)-1)
     VectCart(I+1:I+3)=Zero
     I=3*(ThreeAt(2)-1)
     VectCart(I+2:I+3)=Zero
     I=ThreeAt(3)
     VectCart(3*I)=Zero
   END SUBROUTINE ZeroTrRots
!
!------------------------------------------------------------------
!
   SUBROUTINE ThreeConstr(IntCs,Top12,NatmsLoc,At1,At2,At3)
     TYPE(INTC)     :: IntCs
     INTEGER        :: At1,At2,At3,I,J,NIntC,II,JMax,NatmsLoc
     TYPE(INT_RNK2) :: CountC,Top12
     TYPE(INT_VECT) :: AtVect,Count2
     !
     NIntC=SIZE(IntCs%Def%C)
     CALL New(CountC,(/NatmsLoc,3/))
     CALL New(Count2,NatmsLoc)
     CALL New(AtVect,3)
     CountC%I=0
     Count2%I=0
     DO I=1,NIntC
       IF(IntCs%Constraint%L(I)) THEN
         J=IntCs%Atoms%I(I,1)
         IF(IntCs%Def%C(I)(1:5)=='CARTX') THEN
           CountC%I(J,1)=1
         ELSE IF(IntCs%Def%C(I)(1:5)=='CARTY') THEN
           CountC%I(J,2)=1
         ELSE IF(IntCs%Def%C(I)(1:5)=='CARTZ') THEN
           CountC%I(J,3)=1
         ENDIF
       ENDIF 
     ENDDO
     !
     DO I=1,NatmsLoc
       Count2%I(I)=CountC%I(I,1)+CountC%I(I,2)+CountC%I(I,3)
     ENDDO
     !
     At1=0
     At2=0
     At3=0
     AtVect%I=0
     DO II=1,3
       JMax=MAXVAL(Count2%I)
       IF(JMax>0) THEN
         JMax=0
         DO I=1,NatmsLoc
           J=Top12%I(I,1)*Count2%I(I)
           IF(J>JMax) THEN
             AtVect%I(II)=I
             JMax=J
           ENDIF
         ENDDO
         IF(AtVect%I(II)/=0) Count2%I(AtVect%I(II))=0
       ENDIF
     ENDDO
     At1=AtVect%I(1)
     At2=AtVect%I(2)
     At3=AtVect%I(3)
     !
     ! fix only the origin
     !
    !At2=0
    !At3=0
     !
     CALL Delete(AtVect)
     CALL Delete(Count2)
     CALL Delete(CountC)
   END SUBROUTINE ThreeConstr
!
!------------------------------------------------------------------
!
   FUNCTION HasMetLig(JJ1,JJ2,NJJ1,NJJ2) 
     LOGICAL                :: HasMetLig
     INTEGER                :: JJ1,JJ2,NJJ1,NJJ2
     !
     HasMetLig=.FALSE.
     IF((HasMetal(NJJ1).AND.HasLigand(NJJ2)).OR. &
        (HasMetal(NJJ2).AND.HasLigand(NJJ1))) THEN
       HasMetLig=.TRUE.
     ENDIF
   END FUNCTION HasMetLig
!
!------------------------------------------------------------------
!
   FUNCTION HasLigand(ICharge) 
     LOGICAL :: HasLigand
     INTEGER :: ICharge
     !
     HasLigand=.FALSE.
     IF(ANY(HBondList(:)==ICharge)) HasLigand=.TRUE.
   END FUNCTION HasLigand
!
!------------------------------------------------------------------
!
   FUNCTION HasMetal(ICharge)
     LOGICAL    :: HasMetal
     INTEGER    :: ICharge,I,J
     !
     HasMetal=.FALSE.
     IF(( 3<=ICharge.AND.ICharge<= 5).OR. &
        (11<=ICharge.AND.ICharge<=14).OR. & ! P included
        (19<=ICharge.AND.ICharge<=34).OR. &
        (37<=ICharge.AND.ICharge<=52).OR. &
        (55<=ICharge.AND.ICharge<=84).OR. &
        (87<=ICharge.AND.ICharge<=113)) THEN
       HasMetal=.TRUE.
     ENDIF
   END FUNCTION HasMetal
!
!------------------------------------------------------------------
!
   SUBROUTINE OutPSelect(I1,Top12,XYZ,II2,II3,II4,AngleSum)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INT_RNK2)              :: Top12
     INTEGER                     :: I1,I2,I3,I4,NDim
     INTEGER                     :: II2,II3,II4
     INTEGER                     :: J2,J3,J4
     REAL(DOUBLE)                :: AngleSum,SumU,Conv,Planar,ASum,TwoPi
     REAL(DOUBLE)                :: A1,A2,A3
     !
     AngleSum=-One
     Planar=1.D99
     TwoPi=Two*Pi
     Conv=180.D0/PI
     NDim=Top12%I(I1,1)
     II2=Top12%I(I1,2)
     II3=Top12%I(I1,3)
     II4=Top12%I(I1,4)
     DO J2=1,NDim
       I2=Top12%I(I1,J2+1)
       DO J3=J2+1,NDim
         I3=Top12%I(I1,J3+1)
         CALL BEND(XYZ(1:3,I2),XYZ(1:3,I1),XYZ(1:3,I3),Value_O=A1)
         DO J4=J3+1,NDim
           I4=Top12%I(I1,J4+1)
           CALL BEND(XYZ(1:3,I2),XYZ(1:3,I1),XYZ(1:3,I4), &
                          Value_O=A2)
           CALL BEND(XYZ(1:3,I3),XYZ(1:3,I1),XYZ(1:3,I4), &
                          Value_O=A3)
           SumU=A1+A2+A3
           ASum=ABS(TwoPI-SumU)
           IF(Planar>ASum) THEN
             Planar=ASum
             AngleSum=SumU
             II2=I2
             II3=I3
             II4=I4
           ENDIF
         ENDDO
       ENDDO
     ENDDO
   END SUBROUTINE OutPSelect
!
!------------------------------------------------------------------
!
   SUBROUTINE RotB(B,RotX,RotXY)
     TYPE(BMATR)                 :: B
     REAL(DOUBLE),DIMENSION(3,3) :: RotX,RotXY,Rot
     REAL(DOUBLE),DIMENSION(3)   :: Vect1,Vect2
     INTEGER                     :: I,J,NIntC,JJ
     !
     CALL DGEMM_NNc(3,3,3,One,Zero,RotXY,RotX,Rot)
     NIntC=SIZE(B%IB%I,1)
     DO I=1,NIntC
       DO J=1,4
         IF(B%IB%I(I,J)==0) CYCLE
         JJ=3*(J-1)
         Vect1=B%B%D(I,JJ+1:JJ+3)
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot,Vect1,Vect2)
         B%B%D(I,JJ+1:JJ+3)=Vect2
       ENDDO
     ENDDO
   END SUBROUTINE RotB
!
!------------------------------------------------------------------
!
   SUBROUTINE PermCleans(IPerm,PermRef,IntCs,ThreeAt,DoConstr)
     INTEGER,DIMENSION(:) :: IPerm,PermRef
     INTEGER,DIMENSION(3) :: ThreeAt
     TYPE(INTC)           :: IntCs
     INTEGER              :: I,J,NCart,NIntC,III,K,II
     TYPE(INT_VECT)       :: CleanList
     LOGICAL              :: DoConstr
     !
     NCart=SIZE(IPerm)
     NIntC=SIZE(IntCs%Def%C)
     CALL New(CleanList,NCart)
     CleanList%I=1
     IF(DoConstr) THEN
       DO I=1,NIntC
         IF(IntCs%Constraint%L(I)) THEN
           J=3*(IntCs%Atoms%I(I,1)-1)
           IF(IntCs%Def%C(I)=='CARTX') THEN
             CleanList%I(PermRef(J+1))=0
           ELSE IF(IntCs%Def%C(I)=='CARTY') THEN
             CleanList%I(PermRef(J+2))=0
           ELSE IF(IntCs%Def%C(I)=='CARTZ') THEN
             CleanList%I(PermRef(J+3))=0
           ENDIF
         ENDIF
       ENDDO 
     ENDIF
     J=3*(ThreeAt(1)-1)
     IF(J>=0) THEN
       DO K=1,3 ; CleanList%I(PermRef(J+K))=0 ; ENDDO
     ENDIF
     J=3*(ThreeAt(2)-1)
     IF(J>=0) THEN
       DO K=2,3 ; CleanList%I(PermRef(J+K))=0 ; ENDDO
     ENDIF
     J=3*(ThreeAt(3)-1)
     IF(J>=0) CleanList%I(PermRef(J+3))=0
     !
     III=0
     IPerm=0
     DO I=1,NCart
       II=PermRef(I)
       IF(CleanList%I(II)/=0) THEN
         III=III+1
         IPerm(I)=III
       ENDIF
     ENDDO
     !
     CALL Delete(CleanList)
   END SUBROUTINE PermCleans
!
!------------------------------------------------------------------
!
   SUBROUTINE PermCleanB(ISpB,JSpB,ASpB,IPerm)
     TYPE(INT_VECT)       :: ISpB,JSpB,ISpB2
     TYPE(DBL_VECT)       :: ASpB         
     INTEGER,DIMENSION(:) :: IPerm
     INTEGER              :: I,J,NIntC,NZ,JJ
     !
     NIntC=SIZE(ISpB%I)-1
     CALL New(ISpB2,NIntC+1)
     ISpB2%I=ISpB%I
     NZ=0
     ISpB%I(1)=1
     DO I=1,NIntC
       DO J=ISpB2%I(I),ISpB2%I(I+1)-1
         JJ=IPerm(JSpB%I(J))
         IF(JJ/=0) THEN
           NZ=NZ+1
           JSpB%I(NZ)=JJ
           ASpB%D(NZ)=ASpB%D(J)
         ENDIF
       ENDDO
       ISpB%I(I+1)=NZ+1
     ENDDO
     !
     CALL Delete(ISpB2)
   END SUBROUTINE PermCleanB
!
!------------------------------------------------------------------
!
   SUBROUTINE PermCleanVect(Vect,Perm)
     REAL(DOUBLE),DIMENSION(:)  :: Vect
     INTEGER,DIMENSION(:)       :: Perm
     INTEGER                    :: I,J,Ncart
     TYPE(DBL_VECT)             :: Vect2
     !
     NCart=SIZE(Vect) 
     CALL New(Vect2,NCart)
     Vect2%D=Zero
     DO I=1,NCart
       IF(Perm(I)==0) CYCLE
       Vect2%D(Perm(I))=Vect(I) 
     ENDDO
     Vect=Vect2%D
     CALL Delete(Vect2)
   END SUBROUTINE PermCleanVect
!
!------------------------------------------------------------------
!
   SUBROUTINE ModifyGrad(VectCart,IPerm1,IPerm2,CtrlTrf)
     REAL(DOUBLE),DIMENSION(:) :: VectCart
     INTEGER,DIMENSION(:)      :: IPerm1,IPerm2
     INTEGER                   :: NCart
     TYPE(DBL_VECT)            :: Vect1,Vect2
     TYPE(TrfCtrl)             :: CtrlTrf
     !
     NCart=SIZE(VectCart)
     CALL New(Vect1,NCart)
     CALL New(Vect2,NCart)
     Vect1%D=VectCart
     Vect2%D=VectCart
     !
     CALL RotToAt(Vect1%D,CtrlTrf%RotAt2ToX)
     CALL RotToAt(Vect1%D,CtrlTrf%RotAt3ToXY)
     CALL PermCleanVect(Vect1%D,IPerm1)
     CALL RotToAt(Vect2%D,CtrlTrf%RotAt2ToX_2)
     CALL RotToAt(Vect2%D,CtrlTrf%RotAt3ToXY_2)
     CALL PermCleanVect(Vect2%D,IPerm2)
     VectCart=Vect1%D+Vect2%D
     !
     CALL Delete(Vect1)
     CALL Delete(Vect2)
   END SUBROUTINE ModifyGrad
!
!------------------------------------------------------------------
!
   SUBROUTINE LinBGen(BB1,W,U,V,RU,RV)
     REAL(DOUBLE),DIMENSION(1:12) ::  BB1
     REAL(DOUBLE),DIMENSION(1:3)  ::  W,V,Aux1,Aux2,U
     REAL(DOUBLE)                 ::  RU,RV
     ! 
     CALL CROSS_PRODUCT(U,W,Aux1)
     Aux1=Aux1/RU
     CALL CROSS_PRODUCT(W,V,Aux2)
     Aux2=Aux2/RV
     BB1(1:3)=Aux1
     BB1(4:6)=-Aux1-Aux2
     BB1(7:9)=Aux2
   END SUBROUTINE LinBGen
!
!------------------------------------------------------------------
!
   SUBROUTINE FilterAngles(VectX,VectY)
     TYPE(DBL_VECT) :: VectX,VectY,VectX2,VectY2
     INTEGER        :: I,J,NDim,NDim2
     TYPE(INT_VECT) :: Mark
     REAL(DOUBLE)   :: Center,FilterWidth
     !
     FilterWidth=PI*0.8D0
     NDim=SIZE(VectX%D)
     CALL New(Mark,NDim)
     !
     Center=Sum(VectX%D)/DBLE(NDim)
     CALL Loose2PIs(Center)
     CALL AngleTo180(Center)
     !
     CALL Delete(Mark) 
   END SUBROUTINE FilterAngles
!
!------------------------------------------------------------------
!
   SUBROUTINE PeriodicAngle(VectX,PredVal)
     REAL(DOUBLE),DIMENSION(:) :: VectX
     REAL(DOUBLE)              :: Center,PredVal
     INTEGER                   :: I,NDim
     !
     NDim=SIZE(VectX)
     Center=VectX(1)
     DO I=1,NDim
       CALL PAngle1(Center,VectX(I))
     ENDDO
     CALL PAngle1(Center,PredVal)
   END SUBROUTINE PeriodicAngle
!
!---------------------------------------------------------------------
!
   SUBROUTINE PAngle1(Center,A)
     REAL(DOUBLE)              :: Center,ToCenter,Dist,Dist2,Vect(3),A
     INTEGER                   :: I,J,II
     Vect(1)=A
     Vect(2)=A+TwoPI
     Vect(3)=A-TwoPI
     II=1
     Dist=ABS(Vect(1)-Center)
     DO J=2,3
       Dist2=ABS(Vect(J)-Center)
       IF(Dist2<Dist) THEN
         II=J
         Dist=Dist2 
       ENDIF 
     ENDDO
     A=Vect(II)
   END SUBROUTINE PAngle1
!
!------------------------------------------------------------------
!
   SUBROUTINE AngleOffCenter(VectX,FitVal,Center)
     REAL(DOUBLE),DIMENSION(:) :: VectX
     REAL(DOUBLE)              :: FitVal,Center
     INTEGER                   :: I,J,NDim
     !
     NDim=SIZE(VectX)
     FitVal=Center+FitVal
     DO I=1,NDim
       VectX(I)=Center+VectX(I)
     ENDDO
   END SUBROUTINE AngleOffCenter
!
!------------------------------------------------------------------
!
   SUBROUTINE AngleTo360(Angle)
     REAL(DOUBLE)              :: Angle
     INTEGER                   :: I,J
     ! Input Angles are supposed to be from [-2*PI:2*PI]
     ! Output angles are in                 [0:2*PI]
     IF(Angle<Zero) Angle=TwoPi+Angle
   END SUBROUTINE AngleTo360
!
!------------------------------------------------------------------
!
   SUBROUTINE AngleTo180(Angle)
     REAL(DOUBLE)              :: Angle
     INTEGER                   :: I,J
     ! Input Angles are supposed to be from [0:2*PI]
     ! Output angles are in                 [-PI:PI]
     IF(Angle>PI) Angle=Angle-TwoPI
   END SUBROUTINE AngleTo180
!
!------------------------------------------------------------------
!
   SUBROUTINE BendTo180(Angle) 
     REAL(DOUBLE) :: Angle
     ! Input Angles are supposed to be from [-2*PI:2*PI]
     ! Output angles are in                 [0:PI]
     IF(Angle<Zero) Angle=TwoPI+Angle
     IF(Angle>PI) Angle=TwoPI-Angle
   END SUBROUTINE BendTo180
!
!------------------------------------------------------------------
!
   SUBROUTINE LinBTo180(Angle) 
     REAL(DOUBLE) :: Angle
     ! Input Angles are supposed to be from [-2*PI:2*PI]
     ! Output angles are in                 [-PI:PI]
     ! also for torsions and out-of-planes
     IF(Angle>PI) Angle=Angle-TwoPi
     IF(Angle<-PI) Angle=TwoPI+Angle
   END SUBROUTINE LinBTo180
!
!------------------------------------------------------------------
!
   SUBROUTINE Loose2PIs(Angle)
     REAL(DOUBLE) :: Angle
     Angle=Angle-TwoPi*INT(Angle/TwoPi)
   END SUBROUTINE Loose2PIs
!
!------------------------------------------------------------------
!
   SUBROUTINE MapBendDispl(Angle,AngleDispl)
     REAL(DOUBLE)  :: Angle,AngleDispl,SumU
     SumU=Angle+AngleDispl
     IF(SumU>PI) AngleDispl=(TwoPi-SumU)-Angle
   END SUBROUTINE MapBendDispl
!
!------------------------------------------------------------------
!
   SUBROUTINE MapLinBDispl(AngleDispl)
     REAL(DOUBLE)  :: AngleDispl
     CALL LinBTo180(AngleDispl)         
   END SUBROUTINE MapLinBDispl
!
!------------------------------------------------------------------
!
   SUBROUTINE MapDAngle(Def,Angle,DAngle)
     REAL(DOUBLE)     :: Angle,DAngle 
     CHARACTER(LEN=*) :: Def
     IF(Def(1:4)=='BEND') THEN
       CALL MapBendDispl(Angle,DAngle)
     ELSE IF(Def(1:4)=='LINB'.OR. &
             Def(1:4)=='OUTP'.OR. &
             Def(1:4)=='TORS') THEN
       CALL MapLinBDispl(DAngle)
     ENDIF
   END SUBROUTINE MapDAngle
!
!---------------------------------------------------------
!
   SUBROUTINE MapBackAngle(IntCs,Values)
     IMPLICIT NONE
     TYPE(INTC) :: IntCs
     INTEGER :: I,J,NIntC
     REAL(DOUBLE),DIMENSION(:) :: Values
     REAL(DOUBLE) :: Angle,TwoPi
     !
     ! Map back angles into the ranges the iterative 
     ! back-transformation can cope with.
     NIntC=SIZE(IntCs%Def%C)
     !
     DO I=1,NIntC
       Angle=Values(I)
       IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         CALL Loose2PIs(Angle)
         CALL BendTo180(Angle) 
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS'.OR. &
               IntCs%Def%C(I)(1:4)=='LINB'.OR. &
               IntCs%Def%C(I)(1:4)=='OUTP') THEN
         CALL Loose2PIs(Angle)
         CALL LinBTo180(Angle) 
       ENDIF
       Values(I)=Angle
     ENDDO
   END SUBROUTINE MapBackAngle
!
!-------------------------------------------------------
!
   SUBROUTINE MapAngleDispl(IntCs,Displ) 
     TYPE(INTC)                :: IntCs
     INTEGER                   :: I,NIntC
     REAL(DOUBLE),DIMENSION(:) :: Displ
     !
     NIntC=SIZE(IntCs%Def%C)
     DO I=1,NIntC
       CALL MapDAngle(IntCs%Def%C(I),IntCs%Value%D(I),Displ(I))
     ENDDO
   END SUBROUTINE MapAngleDispl
!
!------------------------------------------------------------------
!
   SUBROUTINE BuildUMatr(SCRPath,NCart)
     CHARACTER(LEN=*) :: SCRPath
     TYPE(INT_VECT)   :: ISpB,JSpB
     TYPE(DBL_VECT)   :: ASpB
     TYPE(INT_VECT)   :: IGc,JGc
     TYPE(DBL_VECT)   :: AGc
     INTEGER          :: NCart,INFO,I,J,NIntC,NZ
     TYPE(DBL_RNK2)   :: FullB,SQFullGc,FullGc,UMatr,Aux1,Aux2
     REAL(DOUBLE)     :: Fact
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     NIntC=Size(ISpB%I)-1
     !
     CALL New(UMatr,(/NIntC,NCart/))
     CALL Sp1x1ToFull(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart,FullB)
     call pprint(FullB,'FullB',Unit_O=6)
     CALL New(FullGc,(/NCart,NCart/))
     CALL New(SQFullGc,(/NCart,NCart/))
     CALL DGEMM_TNc(NCart,NIntC,NCart,One,Zero,FullB%D,&
                    FullB%D,FullGc%D)
     call pprint(FullGc,'FullGc',Unit_O=6)
     !
     CALL SetDSYEVWork(NCart)
       BLKVECT%D=FullGc%D
       CALL DSYEV('V','U',NCart,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in BuildUMatr. INFO='&
                  //TRIM(IntToChar(INFO)))
       NZ=0
       SQFullGc%D=Zero
       DO I=1,NCart
         IF(BLKVALS%D(I)>1.D-7) THEN
           Fact=SQRT(BLKVALS%D(I))
           Fact=One/SQRT(BLKVALS%D(I))
           NZ=NZ+1
           DO J=1,NCart ; SQFullGc%D(J,NZ)=BLKVECT%D(J,I)*Fact ; ENDDO
         ELSE
           Fact=Zero
         ENDIF
       ENDDO     
     call pprint(SQFullGc,'SQFullGc',Unit_O=6)
     CALL UnSetDSYEVWork()
     CALL DGEMM_NNc(NIntC,NCart,NCart,One,Zero,FullB%D,&
                    SQFullGc%D,UMatr%D)
     !
     call pprint(UMatr,'UMatr',Unit_O=6)
     CALL New(Aux1,(/NCart,NCart/))
     CALL DGEMM_TNc(NCart,NIntC,NCart,One,Zero,UMatr%D,&
                    Umatr%D,Aux1%D)
     call pprint(Aux1,'Unit Matr?',Unit_O=6)
     CALL Delete(Aux1)
     !
     CALL WriteBMATR(ISpB,JSpB,ASpB, &
                     TRIM(SCRPath)//'UMatr',UMatr_O=UMatr)
     !
     CALL Delete(SQFullGc)
     CALL Delete(FullGc)
     CALL Delete(FullB)
     CALL Delete(UMatr)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
   END SUBROUTINE BuildUMatr
!
!------------------------------------------------------------------
!
   SUBROUTINE PrimToDeloc(VectInt,VectDeloc,ISpB,JSpB,ASpB,UMatr,Char_O)
     REAL(DOUBLE),DIMENSION(:) :: VectInt,VectDeloc
     CHARACTER(LEN=*),OPTIONAL :: Char_O
     INTEGER                   :: NCart,NIntC
     TYPE(DBL_RNK2)            :: UMatr
     TYPE(INT_VECT)            :: ISpB,JSpB
     TYPE(DBL_VECT)            :: ASpB 
     TYPE(DBL_VECT)            :: AuxCart
     !
     NCart=SIZE(UMatr%D,1)
     CALL New(AuxCart,NCart)
     !
     IF(PRESENT(Char_O)) THEN
       IF(Char_O(1:4)=='Back') THEN
         CALL DGEMM_NNc(NCart,NCart-6,1,One,Zero,UMatr%D,&
                        VectDeloc,AuxCart%D)
         CALL CALC_BxVect(ISpB,JSpB,ASpB,VectInt,AuxCart%D)
       ELSE
         CALL Halt('Erroneous input in PrimToDeloc.')
       ENDIF
     ELSE
       CALL CALC_BxVect(ISpB,JSpB,ASpB,VectInt,AuxCart%D,Trp_O=.TRUE.)
       CALL DGEMM_TNc(NCart-6,NCart,1,One,Zero,UMatr%D,&
                      AuxCart%D,VectDeloc)
     ENDIF
     !
     CALL Delete(AuxCart)
   END SUBROUTINE PrimToDeloc
!
!---------------------------------------------------------------------
!
   SUBROUTINE DelocP(VectInt,ISpB,JSpB,ASpB,UMatr)
     REAL(DOUBLE),DIMENSION(:)  :: VectInt
     TYPE(INT_VECT)             :: ISpB,JSpB
     TYPE(DBL_VECT)             :: ASpB,VectDeloc
     TYPE(DBL_RNK2)             :: UMatr
     INTEGER                    :: NCart
     !
     NCart=SIZE(UMatr%D,1)
     CALL New(VectDeloc,NCart-6)
     CALL PrimToDeloc(VectInt,VectDeloc%D,ISpB,JSpB,ASpB,UMatr)
     CALL PrimToDeloc(VectInt,VectDeloc%D,ISpB,JSpB,ASpB,UMatr, &
                      Char_O='Back')
     CALL Delete(VectDeloc)
   END SUBROUTINE DelocP
!
!---------------------------------------------------------------------
!
   SUBROUTINE Resize_INT_RNK2(Matrix,Dim1New_O,Dim2New_O)
     TYPE(INT_RNK2)   :: Matrix,Matrix2
     INTEGER,OPTIONAL :: Dim1New_O,Dim2New_O
     INTEGER          :: Dim1Old,Dim2Old,Dim1New,Dim2New,I,J
     !
     Dim1Old=SIZE(Matrix%I,1)
     Dim2Old=SIZE(Matrix%I,2)
     Dim1New=Dim1Old
     Dim2New=Dim2Old
     IF(PRESENT(Dim1New_O)) Dim1New=Dim1New_O
     IF(PRESENT(Dim2New_O)) Dim2New=Dim2New_O
     !
     CALL New(Matrix2,(/Dim1New,Dim2New/))
     Matrix2%I=0   
     DO J=1,Dim2Old
       DO I=1,Dim1Old
         Matrix2%I(I,J)=Matrix%I(I,J) 
       ENDDO
     ENDDO
     CALL Delete(Matrix)
     CALL New(Matrix,(/Dim1New,Dim2New/))
     Matrix%I=Matrix2%I
     CALL Delete(Matrix2)
   END SUBROUTINE Resize_INT_RNK2
!
!---------------------------------------------------------------------
!
   SUBROUTINE ReorderN(VectX,IWork,NDim)
     REAL(DOUBLE),DIMENSION(:) :: VectX
     INTEGER     ,DIMENSION(:) :: IWork
     REAL(DOUBLE)              :: X1,X2
     INTEGER                   :: I,J,NDim,I1,I2
     ! order set by decreasing x
     DO I=1,NDim ; IWork(I)=I ; ENDDO
     DO I=1,NDim-1
       DO J=NDim-1,1,-1
         I1=IWork(J)
         I2=IWork(J+1)
         X1=VectX(I1)
         X2=VectX(I2) 
         IF(X2>X1) THEN
           IWork(J)=I2
           IWork(J+1)=I1
          !VectX(I1)=X2
          !VectX(I2)=X1 
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE ReorderN 
!
!--------------------------------------------------------------------
!
   SUBROUTINE ReorderI(VectX,IWork,NDim)
     REAL(DOUBLE),DIMENSION(:) :: VectX
     INTEGER     ,DIMENSION(:) :: IWork
     REAL(DOUBLE)              :: X1,X2
     INTEGER                   :: I,J,NDim,I1,I2
     ! order set by increasing x
     DO I=1,NDim ; IWork(I)=I ; ENDDO
     DO I=1,NDim-1
       DO J=1,NDim-I
         I1=IWork(J)
         I2=IWork(J+1)
         X1=VectX(I1)
         X2=VectX(I2) 
         IF(X1>X2) THEN
           IWork(J)=I2
           IWork(J+1)=I1
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE ReorderI 
!
!---------------------------------------------------------------------
!
   SUBROUTINE Reorder1D(VectX)
     REAL(DOUBLE),DIMENSION(:) :: VectX
     REAL(DOUBLE)              :: X1,X2
     INTEGER                   :: I,J,NDim
     ! order set by increasing x
     NDim=SIZE(VectX)
     DO I=1,NDim-1
       DO J=1,NDim-I
         X1=VectX(J)
         X2=VectX(J+1)
         IF(X1>X2) THEN
           VectX(J)=X2
           VectX(J+1)=X1
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE Reorder1D
!
!---------------------------------------------------------------------
!
   SUBROUTINE Reorder(VectX,VectY)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: X1,X2,Y1,Y2
     INTEGER                   :: I,J,NDim
     ! order set by increasing x
     NDim=SIZE(VectX)
     DO I=1,NDim-1
       DO J=1,NDim-I
         X1=VectX(J)
         X2=VectX(J+1)
         Y1=VectY(J)
         Y2=VectY(J+1)
         IF(X1>X2) THEN
           VectX(J)=X2
           VectX(J+1)=X1
           VectY(J)=Y2
           VectY(J+1)=Y1
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE Reorder
!
!---------------------------------------------------------------------
!
   SUBROUTINE SortBonds(NatmsLoc,AtmB,Bond)
     TYPE(BondDATA) :: Bond
     TYPE(ATOMBONDS):: AtmB
     INTEGER        :: NatmsLoc,I,I1,I2,MaxBonds
     !
     MaxBonds=20
     IF(AllocQ(AtmB%Alloc)) THEN
       CALL Delete(AtmB)
     ENDIF
     CALL New(AtmB,NatmsLoc,MaxBonds)
     AtmB%N1=NatmsLoc
     AtmB%Count%I=0
     AtmB%Bonds%I=0
     AtmB%Atoms%I=0
     DO I=1,Bond%N
       I1=Bond%IJ%I(1,I)
       I2=Bond%IJ%I(2,I)
       AtmB%Count%I(I1)=AtmB%Count%I(I1)+1
       AtmB%Count%I(I2)=AtmB%Count%I(I2)+1
       IF(AtmB%Count%I(I1)>MaxBonds.OR.AtmB%Count%I(I2)>MaxBonds) THEN
         CALL Resize_INT_RNK2(AtmB%Bonds,Dim2New_O=MaxBonds+10)
         CALL Resize_INT_RNK2(AtmB%Atoms,Dim2New_O=MaxBonds+10)
         MaxBonds=MaxBonds+10
       ENDIF
       AtmB%Bonds%I(I1,AtmB%Count%I(I1))=I
       AtmB%Atoms%I(I1,AtmB%Count%I(I1))=I2
       !
       AtmB%Bonds%I(I2,AtmB%Count%I(I2))=I
       AtmB%Atoms%I(I2,AtmB%Count%I(I2))=I1
     ENDDO
     AtmB%N2=MaxBonds
   END SUBROUTINE SortBonds
!
!---------------------------------------------------------------------
!
   SUBROUTINE IntCBoxes(XYZ,Box,BoxSize_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(IntCBox)               :: Box
     INTEGER                     :: NatmsLoc
     REAL(DOUBLE)                :: BXMIN,BYMIN,BZMIN
     REAL(DOUBLE)                :: BoxSize
     REAL(DOUBLE),OPTIONAL       :: BoxSize_O
     !
     NatmsLoc=SIZE(XYZ,2)
     BoxSize=3.0D0*AngstromsToAU !in A
     IF(PRESENT(BoxSize_O)) BoxSize=BoxSize_O*AngstromsToAU
     CALL SORT_INTO_Box1(BoxSize,XYZ,NatmsLoc,&
                         Box%NX,Box%NY,Box%NZ,BXMIN,BYMIN,BZMIN)
     !
     Box%N=Box%NX*Box%NY*Box%NZ
     CALL New(Box%I,Box%N+1)
     CALL New(Box%J,NatmsLoc)
     !
     CALL SORT_INTO_Box2(BoxSize,XYZ,NatmsLoc,Box%NX,Box%NY,Box%NZ,&
                         BXMIN,BYMIN,BZMIN,Box%I,Box%J)
   END SUBROUTINE IntCBoxes
!
!---------------------------------------------------------------------
!
   SUBROUTINE BondingScheme(XYZ,AtNum,AtmB,Bond,TOPS, &
                            GCoordCtrl,GConvCr,Cells,IEq)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER,DIMENSION(:)        :: AtNum,IEq
     TYPE(BONDDATA)              :: Bond,BondF,Bond1,BondCov,BondVDW
     TYPE(BONDDATA)              :: BondVDWPure
     TYPE(ATOMBONDS)             :: AtmB,AtmBF
     TYPE(TOPOLOGY)              :: TOPS
     INTEGER,DIMENSION(:,:)      :: Cells
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(GConvCrit)             :: GConvCr   
     TYPE(DBL_VECT)              :: CritRad,StRad
     REAL(DOUBLE)                :: Fact,HBondMax,BoxSize
     INTEGER                     :: I,J,N,NatmsLoc,IFrags,NFrag
     INTEGER                     :: MaxBonds,NBondEst,IntSet
     LOGICAL                     :: DoRepeat
     TYPE(INT_VECT)              :: FragID
     TYPE(IntCBox)               :: Box
     !
     !now define bonding scheme, based on Slater or Van der Waals radii
     !
     NatmsLoc=SIZE(XYZ,2)
     HBondMax=2.50D0*AngstromsToAu !*MIN(GCoordCtrl%VDWFact,One) ! in Au
     !
     CALL New(CritRad,NatmsLoc)
     CALL New(StRad,NatmsLoc)
   ! DO IntSet=1,3
     DO IntSet=1,1
       IF(IntSet==1) THEN
         N=SIZE(SLRadii,1)
         Fact=1.3D0 !!! Scaling factor for Slater Radii
        !Fact=1.0D0 !!! Scaling factor for Slater Radii
         DO I=1,NatmsLoc
           CritRad%D(I)=Fact*SLRadii(AtNum(I))*AngstromsToAU
         ENDDO
       ELSE IF(IntSet==2) THEN
         N=SIZE(VDWRadii,1)
         Fact=GCoordCtrl%VDWFact !!! Scaling factor for VDW Radii
         DO I=1,NatmsLoc
           CritRad%D(I)=Fact*VDWRadii(AtNum(I))*AngstromsToAU
         ENDDO
       ELSE IF(IntSet==3) THEN
         N=SIZE(SLRadii,1)
         Fact=1.6D0 !!! Scaling factor for Slater Radii
        !Fact=1.0D0 !!! Scaling factor for Slater Radii
         DO I=1,NatmsLoc
           CritRad%D(I)=Fact*SLRadii(AtNum(I))*AngstromsToAU
         ENDDO
       ENDIF
       StRad%D=CritRad%D
       !
       IF(IntSet==1) THEN
         CALL IntCBoxes(XYZ,Box,BoxSize_O=3.0D0)
         CALL New(FragID,NatmsLoc)
         DO J=1,NatmsLoc ; FragID%I(J)=J ; ENDDO
         NFrag=NatmsLoc
         !
         DO IFrags=1,10000
         ! CALL BondList2(XYZ,AtNum,Box,BondCov,IEq,CritRad%D)  
           CALL BondList(XYZ,AtNum,IntSet,Box,BondCov,TOPS, &
                         CritRad,HbondMax,GCoordCtrl%LinCrit, &
                         GConvCr,IEq,FragID%I,NFrag)
           CALL SortBonds(NatmsLoc,AtmB,BondCov)
           CALL Topology_12(AtmB,TOPS%Cov12)
           CALL SortFragments(TOPS%Cov12%I,FragID%I,NFrag)
          !Fact=(One+0.05D0*(IFrags-1))
           Fact=1.05D0**(IFrags-1)
           IF(NFrag/=1.AND.Fact<1.60D0) THEN
           ! CritRad%D=1.05D0*CritRad%D 
             CritRad%D=Fact*StRad%D
             CALL Delete(TOPS%Cov12)
           ELSE
             EXIT
           ENDIF
           CALL Delete(AtmB)
         ENDDO
       ! IF(NFrag/=1) CALL Halt('Fragmented system after primary recognition.')
         CALL Delete(FragID) 
         CALL Delete(Box)
         CALL Topology_13(NatmsLoc,TOPS%Cov12,TOPS%Cov13)
         CALL Topology_14(NatmsLoc,TOPS%Cov12,TOPS%Cov14)
         CALL Excl_List(NatmsLoc,TOPS%Cov12,TOPS%Cov13,TOPS%Cov14, &
                        TOPS%CovExcl)
       ELSE IF(IntSet==2) THEN
         CALL IntCBoxes(XYZ,Box,BoxSize_O=3.0D0)
         CALL New(FragID,NatmsLoc)
         DO J=1,NatmsLoc ; FragID%I(J)=J ; ENDDO
         NFrag=NatmsLoc
         !
         CALL BondList(XYZ,AtNum,IntSet,Box,Bond1,TOPS, &
                       CritRad,HBondMax,GCoordCtrl%LinCrit, &
                       GConvCr,IEq,FragID%I,NFrag)
         !
         CALL Delete(FragID) 
         CALL Delete(Box)
       ELSE IF(IntSet==3) THEN
         !
         IF(GConvCr%NoFragmConnect) THEN
           BondF%N=0
           CYCLE
         ENDIF
        !CALL ConnectFragments(XYZ,AtNum,BondF,TOPS,Cells,IEq)
         !
         CALL New(FragID,NatmsLoc)
         BondF%N=0
         CALL SetEq(BondVDWPure,Bond1)
         DO IFrags=1,10000
          !CALL SortBonds(NatmsLoc,AtmBF,BondF)
          !CALL SortNonCov2(AtNum,XYZ,BondF,AtmBF)
          !CALL Delete(AtmBF)
           CALL MergeBonds(BondVDWPure,BondF,BondVDW)
           CALL Delete(BondVDWPure)
           CALL SetEq(BondVDWPure,BondVDW)
          !CALL SortBonds(NatmsLoc,AtmB,BondVDW)
          !CALL SortNonCov2(AtNum,XYZ,BondVDW,AtmB)
          !CALL Delete(AtmB)
           CALL VDWTop(TOPS%Tot12,TOPS%Cov12,BondVDW%IJ,BondVDW%N)
           CALL SortFragments(TOPS%Tot12%I,FragID%I,NFrag)
           CALL Delete(BondVDW)
           CALL Delete(TOPS%Tot12)
           CALL Delete(BondF)
          !Fact=1.05D0**(IFrags-1)
           Fact=(One+0.05D0*(IFrags-1))
           IF(NFrag/=1) THEN
             CritRad%D=Fact*StRad%D
           ELSE
             CALL Delete(BondVDW)
             EXIT
           ENDIF
           !
           BoxSize=3.D0*Fact
          !BoxSize=3.D0*DBLE(IFrags)
           CALL IntCBoxes(XYZ,Box,BoxSize_O=BoxSize)
           !
           CALL BondList(XYZ,AtNum,IntSet,Box,BondF,TOPS, &
                         CritRad,HbondMax,GCoordCtrl%LinCrit, &
                         GConvCr,IEq,FragID%I,NFrag)
           !
           CALL Delete(Box)
         ENDDO
         IF(NFrag/=1) CALL Halt('Fragmented system after VDW recognition.')
         CALL Delete(FragID) 
       ENDIF
     ENDDO
     !
     CALL Delete(Bond1)
   ! CALL SortBonds(NatmsLoc,AtmB,BondVDW)
   ! CALL SortNonCov2(AtNum,XYZ,BondVDW,AtmB)
   ! CALL Delete(AtmB)
     !
     CALL VDWTop(TOPS%Tot12,TOPS%Cov12,BondVDWPure%IJ,BondVDWPure%N)
     CALL Topology_13(NatmsLoc,TOPS%Tot12,TOPS%Tot13)
     CALL Topology_14(NatmsLoc,TOPS%Tot12,TOPS%Tot14)
     CALL Excl_List(NatmsLoc,TOPS%Tot12,TOPS%Tot13,TOPS%Tot14, &
                    TOPS%TotExcl)
     !
     CALL MergeBonds(BondCov,BondVDWPure,Bond)
     CALL Delete(BondCov)
     CALL Delete(BondVDWPure)
     CALL SortBonds(NatmsLoc,AtmB,Bond)
     CALL Delete(StRad) 
     CALL Delete(CritRad)
   END SUBROUTINE BondingScheme
!
!--------------------------------------------------------
!
   SUBROUTINE BondList(XYZ,AtNum,IntSet,Box,Bond,TOPS, &
                       CritRad,HBondMax,LinCrit,GConvCr,IEq, &
                       FragID,NFrag)
     IMPLICIT NONE
     INTEGER                     :: I,J,NatmsLoc,NBond
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(BONDDATA)              :: Bond
     TYPE(IntCBox)               :: Box
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(GConvCrit)             :: GConvCr
     TYPE(INT_RNK2)              :: FTop
     TYPE(INT_VECT)              :: IFrag,JFrag
     INTEGER,DIMENSION(:)        :: AtNum,IEq,FragID
     TYPE(DBL_VECT)              :: CritRad
     TYPE(INT_VECT)              :: Neighbors
     REAL(DOUBLE)                :: R12,R12_2,CritDist,HBondMax,LinCrit
     REAL(DOUBLE)                :: OriginalRad,DeltaRep
     INTEGER                     :: IZ,IX,IY,I1,I2,JJ1,JJ2,F1,F2,NFrag
     INTEGER                     :: IORD,IORDD
     INTEGER                     :: IZD,IXD,IYD,NJJ1,NJJ2
     INTEGER                     :: NMax12,JJ,IntSet,IDimExcl,NBondEst
     INTEGER                     :: IDim14,IDim12,IDim13
     INTEGER                     :: HAtm,LAtm,DDimU,NBondOld
     INTEGER                     :: IChk,MaxRepeat
     REAL(DOUBLE),DIMENSION(3)   :: DVect
     LOGICAL                     :: FoundHBond,FoundMetLig
     LOGICAL                     :: LonelyAtom,DoExclude
     LOGICAL                     :: NearestOnly,AtomRepeat
     !     
     NatmsLoc=SIZE(XYZ,2)
     HAtm=0
     NBondEst=Bond%N
     NBond=NBondEst
     CALL New(Neighbors,NatmsLoc)
     Neighbors%I=0
     MaxRepeat=1
     AtomRepeat=.FALSE.
     IF(AtomRepeat) THEN
       DeltaRep=0.05D0 ! 5% increment
       MaxRepeat=INT(10.D0*(1.D0/DeltaRep))
     ENDIF
     NearestOnly=.FALSE.
     IF(NearestOnly) THEN
       CALL FragIJs(FragID,NFrag,IFrag,JFrag)
       CALL New(FTop,(/NFrag,11/))
       FTop%I=0
     ENDIF
     !
     !  Go through all boxes and their neighbours
     !
     DO IZ=1,Box%NZ
       DO IX=1,Box%NX
         DO IY=1,Box%NY
           ! absolute index of a box
           IOrd=Box%NX*Box%NY*(IZ-1)+Box%NY*(IX-1)+IY
           DO I1=Box%I%I(IOrd),Box%I%I(IOrd+1)-1
             JJ1=Box%J%I(I1) !!! atom in central box
             NJJ1=AtNum(JJ1)
             OriginalRad=CritRad%D(JJ1)
             DO IChk=1,MaxRepeat
               ! second atom may come from central or neigbouring Boxes
               !and must be an MM atom,LJ is not calculated for QM-QMpairs
               DO IZD=-1,1
                 IF(IZ+IZD>0 .AND. IZ+IZD<=Box%NZ) THEN
                   DO IXD=-1,1
                     IF(IX+IXD>0 .AND. IX+IXD<=Box%NX) THEN
                       DO IYD=-1,1
                         IF(IY+IYD>0 .AND. IY+IYD<=Box%NY) THEN
                           IOrdD=Box%NX*Box%NY*(IZ-1+IZD)+&
                                 Box%NY*(IX-1+IXD)+IY+IYD
                           DO I2=Box%I%I(IOrdD),Box%I%I(IOrdD+1)-1
                             JJ2=Box%J%I(I2) !!! second atom
                             NJJ2=AtNum(JJ2)
                             IF(IntSet/=2) THEN 
                               F1=FragID(JJ1)
                               F2=FragID(JJ2)
                               IF(F1==F2) CYCLE
                               IF(NearestOnly) THEN
                                 IF(ConnectedF(F1,F2,FTop%I)) CYCLE
                               ENDIF
                             ENDIF
                             FoundHBond=.FALSE.
                             FoundMetLig=.FALSE.
                             LonelyAtom=.FALSE.
                             IF(Neighbors%I(JJ1)==0.AND.IntSet==1 &
                               .AND.AtomRepeat) THEN
                               IF(JJ2==JJ1) CYCLE 
                             ELSE
                               IF(JJ2<=JJ1) CYCLE 
                             ENDIF
                             IF(IntSet==2) THEN
                               FoundMetLig=HasMetLig(JJ1,JJ2,NJJ1,NJJ2)
                               LonelyAtom=HasLonelyAtm(TOPS%Cov12,JJ1,JJ2,LAtm)
                               FoundHBond=HasHBond(TOPS%Cov12,AtNum, &
                                          NJJ1,NJJ2,JJ1,JJ2,HAtm,XYZ,&
                                          LinCrit)
                               CALL BondExcl(JJ1,JJ2,NJJ1,NJJ2,TOPS, &
                                             FoundHBond,FoundMetLig,&
                                             LonelyAtom,DoExclude)
                               IF(GConvCr%HBondOnly.AND. &
                                  .NOT.FoundHBond) DoExclude=.TRUE.
                               IF(DoExclude) CYCLE
                             ENDIF
                             CritDist=CritRad%D(JJ1)+CritRad%D(JJ2)
                             DVect(:)=XYZ(:,JJ1)-XYZ(:,JJ2)
                             R12_2=DOT_PRODUCT(DVect,DVect)
                             R12=SQRT(R12_2)
                             IF(R12<CritDist) THEN
                               IF(IntSet/=2) THEN
                                 ! this will overwrite JJ1,JJ2 and R12
                                 IF(NearestOnly) THEN
                                   CALL MergeFrag(JJ1,JJ2,R12, &
                                        XYZ,FragID,Ftop,IFrag,JFrag)
                                 ENDIF
                               ENDIF
                               !
                               IF(FoundHBond.AND..NOT.LonelyAtom) THEN
                                 IF(R12>HBondMax) CYCLE
                               ENDIF
                               IF(NBond+1>NBondEst) THEN
                                 DDimU=NatmsLoc*10
                                 CALL MoreBondArray(Bond,DDimU,NBondEst)
                                 NBondEst=NBondEst+DDimU
                               ENDIF 
                               NBond=NBond+1
                               IF(IEq(JJ1)<IEq(JJ2)) THEN
                                 Bond%IJ%I(1:2,NBond)=(/JJ1,JJ2/)
                               ELSE
                                 Bond%IJ%I(1:2,NBond)=(/JJ2,JJ1/)
                               ENDIF
                               Bond%Length%D(NBond)=R12
                               Neighbors%I(JJ1)=Neighbors%I(JJ1)+1
                               Neighbors%I(JJ2)=Neighbors%I(JJ2)+1
                               !
                               IF(IntSet==1) THEN
                                 Bond%Type%C(NBond)(1:3)='COV'
                               ELSE IF(IntSet==2) THEN
                                 Bond%Type%C(NBond)(1:3)='VDW'
                                 IF(FoundHBond.AND..NOT.LonelyAtom) THEN
                                   Bond%Type%C(NBond)(1:5)='HBond'
                                 ELSE IF(FoundMetLig) THEN
                                   Bond%Type%C(NBond)(1:6)='MetLig'
                                 ENDIF
                               ELSE IF(IntSet==3) THEN
                                 Bond%Type%C(NBond)(1:4)='Frac'
                               ENDIF
                             ENDIF
                           ENDDO
                         ENDIF
                       ENDDO
                     ENDIF
                   ENDDO
                 ENDIF
               ENDDO
               IF(AtomRepeat) THEN
                 IF(IntSet==1) THEN
                 ! Increase radius for lonely atoms in covalent scheme
                   IF(Neighbors%I(JJ1)==0) THEN 
                     CritRad%D(JJ1)=CritRad%D(JJ1)+0.02D0*OriginalRad
                   ELSE
                    !CritRad%D(JJ1)=OriginalRad
                     EXIT  
                   ENDIF
                 ELSE
                   EXIT
                 ENDIF
               ENDIF
             ENDDO !!! IChk for checking LonelyAtoms
           ENDDO !!! central box atoms
         ENDDO
       ENDDO
     ENDDO !!! ends on central box indices
     !
     ! Compress Bond
     !
     IF(NearestOnly) THEN
       CALL Delete(IFrag)
       CALL Delete(JFrag)
       CALL Delete(FTop)
     ENDIF
     CALL MoreBondArray(Bond,0,NBond)
     CALL Delete(Neighbors)
   END SUBROUTINE BondList
!
!--------------------------------------------------------------
!
   SUBROUTINE BondList2(XYZ,AtNum,Box,Bond,IEq,CritRad)  
     IMPLICIT NONE
     INTEGER                     :: I,J,NatmsLoc,NBond,K,L
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: CritRad
     TYPE(BONDDATA)              :: Bond
     TYPE(IntCBox)               :: Box
     INTEGER,DIMENSION(:)        :: AtNum,IEq
     INTEGER                     :: IZ,IX,IY,I1,I2,JJ1,JJ2
     INTEGER                     :: IORD,IORDD
     INTEGER                     :: IZD,IXD,IYD,NJJ1,NJJ2
     INTEGER                     :: JJ,NBondEst
     INTEGER                     :: NBondOld,DDimU,AuxSize,ICount
     INTEGER                     :: Info,MaxAtomBonds,NDim
     REAL(DOUBLE),DIMENSION(3)   :: DVect
     REAL(DOUBLE),DIMENSION(3,3) :: Theta
     REAL(DOUBLE)                :: R12_2,R12,Fact,Q
     TYPE(INT_VECT)              :: Atoms,IOrder
     TYPE(DBL_VECT)              :: VectB,QVals
     TYPE(DBL_RNK2)              :: Vects,Vects2
     TYPE(ATOMBONDS)             :: AtmB,AtmB2
     !     
     NatmsLoc=SIZE(XYZ,2)
     NBondEst=Bond%N
     NBond=NBondEst
     AuxSize=(1+(NatmsLoc-3)/(Box%NZ+Box%NX+Box%NY))*27*5
     CALL New(Atoms,AuxSize)
     CALL New(IOrder,AuxSize)
     CALL New(Vects,(/AuxSize,3/))
     CALL New(Vects2,(/AuxSize,3/))
     CALL New(VectB,AuxSize)
     CALL New(QVals,AuxSize)
     CALL SetDSYEVWork(3)
     MaxAtomBonds=10
     CALL New(AtmB,NatmsLoc,MaxAtomBonds)
     !
     !  Go through all boxes and their neighbours
     !
     DO IZ=1,Box%NZ
       DO IX=1,Box%NX
         DO IY=1,Box%NY
           ! absolute index of a box
           IOrd=Box%NX*Box%NY*(IZ-1)+Box%NY*(IX-1)+IY
           DO I1=Box%I%I(IOrd),Box%I%I(IOrd+1)-1
             JJ1=Box%J%I(I1) !!! atom in central box
             NJJ1=AtNum(JJ1)
             Theta=Zero
             ICount=0
               DO IZD=-1,1
                 IF(IZ+IZD>0 .AND. IZ+IZD<=Box%NZ) THEN
                   DO IXD=-1,1
                     IF(IX+IXD>0 .AND. IX+IXD<=Box%NX) THEN
                       DO IYD=-1,1
                         IF(IY+IYD>0 .AND. IY+IYD<=Box%NY) THEN
                           IOrdD=Box%NX*Box%NY*(IZ-1+IZD)+&
                                 Box%NY*(IX-1+IXD)+IY+IYD
                           DO I2=Box%I%I(IOrdD),Box%I%I(IOrdD+1)-1
                             JJ2=Box%J%I(I2) !!! second atom
                             NJJ2=AtNum(JJ2)
                             IF(JJ1==JJ2) CYCLE
                             ICount=ICount+1
                             DVect(:)=XYZ(:,JJ1)-XYZ(:,JJ2)
                             R12_2=DOT_PRODUCT(DVect,DVect)
                             R12=SQRT(R12_2)
                             DVect=One/(CritRad(JJ2)/R12)**2/R12*DVect
                             Atoms%I(ICount)=JJ2
                             DO J=1,3 ; Vects%D(ICount,J)=DVect(J) ; ENDDO
                             DO I=1,3
                               DO J=1,3
                                 Theta(I,J)=Theta(I,J)+DVect(I)*DVect(J)
                               ENDDO
                             ENDDO
                           ENDDO
                         ENDIF
                       ENDDO
                     ENDIF
                   ENDDO
                 ENDIF
               ENDDO
               ! get principal directions
               BLKVECT%D=Theta
               CALL DSYEV('V','U',3,BLKVECT%D,BIGBLOK,BLKVALS%D, &
                          BLKWORK%D,BLKLWORK,INFO)
               !
               DO J=1,3 
                 IF(BLKVALS%D(J)/BLKVALS%D(3)>0.1D0) THEN
                   BLKVALS%D(J)=One/SQRT(BLKVALS%D(J))
                 ELSE
                   BLKVALS%D(J)=Zero
                 ENDIF 
               ENDDO 
               CALL DGEMM_NNc(ICount,3,3,One,Zero,Vects%D(1:icount,1:3),&
                              BLKVECT%D,Vects2%D(1:icount,1:3))
               DO I=1,ICount
                 DO J=1,3
                   Vects2%D(I,J)=Vects2%D(I,J)*BLKVALS%D(J)
                 ENDDO
               ENDDO
               !
               ! Recognize bonds
               !
               VectB%D=Zero
               DO I=1,ICount
                 DO J=1,3
                   VectB%D(I)=VectB%D(I)+Vects2%D(I,J)**2 
                 ENDDO
               ENDDO
               Fact=SUM(VectB%D(1:ICount))
               VectB%D=VectB%D/Fact
               !
               CALL ReorderN(VectB%D,IOrder%I,ICount) 
               J=IOrder%I(1)
               QVals%D(1)=Zero
               DO I=ICount-1,1,-1
                 K=IOrder%I(I)
                 L=IOrder%I(I+1)
                 QVals%D(I+1)=ABS(LOG10(VectB%D(L))-LOG10(VectB%D(K)))/ &
                             (ABS(LOG10(VectB%D(J))-LOG10(VectB%D(L)))+1.D-10)
               ENDDO
               ! terminate series of bonds if next point is 
               ! droppable by Q-test
               ! and more than 50% of bonding is already described
               NDim=ICount
               Q=Zero                  
               DO I=1,ICount-1
                 J=IOrder%I(I)
                 Q=Q+VectB%D(J)
                 IF(Q>0.67D0) THEN
                   IF(I/=1.AND.QVals%D(I+1)<QTest90(MIN(10,I+1))) CYCLE
                   NDim=I
                   EXIT
                 ENDIF
               ENDDO
               !
               Q=Zero
               DO I=1,NDim
                 J=IOrder%I(I)
                 JJ2=Atoms%I(J)
                 IF(.NOT.ANY(AtmB%Atoms%I(JJ2,:)==JJ1)) THEN ! avoid double counting of bonds
                   NBond=NBond+1
                   IF(NBond>NBondEst) THEN
                     DDimU=NatmsLoc*10
                     CALL MoreBondArray(Bond,DDimU,NBondEst)
                     NBondEst=NBondEst+DDimU
                   ENDIF 
                   DDimU=MAX(AtmB%Count%I(JJ1),AtmB%Count%I(JJ2))
                   IF(MaxAtomBonds>DDimU) THEN
                     MaxAtomBonds=MaxAtomBonds+10
                     CALL SetEq(AtmB2,AtmB,N2_O=MaxAtomBonds)
                     CALL Delete(AtmB)
                     CALL SetEq(AtmB,AtmB2)
                     CALL Delete(AtmB2)
                   ENDIF 
                   !
                   IF(IEq(JJ1)<IEq(JJ2)) THEN
                     Bond%IJ%I(1:2,NBond)=(/JJ1,JJ2/)
                   ELSE
                     Bond%IJ%I(1:2,NBond)=(/JJ2,JJ1/)
                   ENDIF
                   AtmB%Count%I(JJ1)=AtmB%Count%I(JJ1)+1
                   AtmB%Count%I(JJ2)=AtmB%Count%I(JJ2)+1
                   AtmB%Atoms%I(JJ1,AtmB%Count%I(JJ1))=JJ2
                   AtmB%Atoms%I(JJ2,AtmB%Count%I(JJ2))=JJ1
                   AtmB%Bonds%I(JJ1,AtmB%Count%I(JJ1))=NBond
                   AtmB%Bonds%I(JJ2,AtmB%Count%I(JJ2))=NBond
                   !
                   DVect(:)=XYZ(:,JJ1)-XYZ(:,JJ2)
                   R12_2=DOT_PRODUCT(DVect,DVect)
                   R12=SQRT(R12_2)
                   Bond%Length%D(NBond)=R12
                   Bond%Type%C(NBond)(1:3)='COV'
                 ENDIF
               ENDDO
               !
           ENDDO !!! central box atoms
         ENDDO
       ENDDO
     ENDDO !!! ends on central box indices
     CALL Delete(Atoms)
     CALL Delete(IOrder)
     CALL Delete(Vects)
     CALL Delete(Vects2)
     CALL Delete(VectB)
     CALL Delete(QVals)
     CALL Delete(AtmB)
     !
     ! Compress Bond
     !
     CALL UnSetDSYEVWork()
     CALL MoreBondArray(Bond,0,NBond)
   END SUBROUTINE BondList2
!
!--------------------------------------------------------------
!
   SUBROUTINE Topology_12(AtmB,Top12)
     TYPE(ATOMBONDS):: AtmB
     TYPE(INT_RNK2) :: Top12
     INTEGER        :: I,J,NatmsLoc,MaxDim
     !
     NatmsLoc=SIZE(AtmB%Count%I)
     MaxDim=MAXVAL(AtmB%Count%I)
     CALL New(Top12,(/NatmsLoc,MaxDim+1/))
     Top12%I=0
     DO I=1,NatmsLoc
       Top12%I(I,1)=AtmB%Count%I(I)
       DO J=1,AtmB%Count%I(I)
         Top12%I(I,J+1)=AtmB%Atoms%I(I,J)
       ENDDO 
     ENDDO 
   END SUBROUTINE Topology_12
!
!--------------------------------------------------------------
!
   SUBROUTINE Topology_13(NatmsLoc,Top12,Top13)
     ! Set up a table which shows the atom numbers of Atoms 
     ! being second neighbours of a certain atom.
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)  :: Top12
     TYPE(INT_RNK2)  :: Top13
     TYPE(INT_RNK2)  :: Top13_2
     INTEGER         :: I,J,K,L,N,M,II,JJ,NI,NJ,NatmsLoc
     INTEGER         :: NMax13,NMax12,KK,IN12,JN12
     !
     NMax12=Size(Top12%I,2)-1
     !
     NMax13=30
     K=NMax13+1
     CALL New(Top13,(/NatmsLoc,K/))
     Top13%I(1:NatmsLoc,1:NMax13+1)=0 
     !
     DO II=1,NatmsLoc
       IN12=Top12%I(II,1)
       IF(IN12==0) CYCLE
       DO J=1,IN12
         JJ=Top12%I(II,J+1)
         JN12=Top12%I(JJ,1)
         DO K=1,JN12
         KK=Top12%I(JJ,K+1)
           IF(II/=KK) THEN
             !    
             NI=Top13%I(II,1)
             !    
             !   check matrix Size, increase Size if necessary
             !    
             IF(NI>=NMax13) THEN
               NMax13=NMax13+30
               CALL New(Top13_2,(/NatmsLoc,NMax13+1/))
               Top13_2%I(1:NatmsLoc,1:NMax13+1)=0 
               Top13_2%I(1:NatmsLoc,1:NMax13+1-30)=&
                                       Top13%I(1:NatmsLoc,1:NMax13+1-30)
               CALL Delete(Top13)
               CALL New(Top13,(/NatmsLoc,NMax13+1/))
               Top13%I(1:NatmsLoc,1:NMax13+1)=&
                                        Top13_2%I(1:NatmsLoc,1:NMax13+1)
               CALL Delete(Top13_2)
             ENDIF
             IF(NI/=0) THEN
               IF(ANY(Top13%I(II,2:NI+1)==KK)) THEN
                 CYCLE
               ELSE
                 Top13%I(II,1)=NI+1
                 Top13%I(II,1+(NI+1))=KK
               ENDIF
             ELSE
                 Top13%I(II,1)=NI+1
                 Top13%I(II,1+(NI+1))=KK
             ENDIF
           ENDIF !!! II/=KK
         ENDDO !!!! KK
       ENDDO !!!! JJ
     ENDDO !!!! II
   END SUBROUTINE Topology_13 
!
!--------------------------------------------------------------
!
   SUBROUTINE Topology_14(NatmsLoc,Top12,Top14)
     ! Set up a table which shows the atom numbers of Atoms 
     ! being second neighbours of a certain atom.
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)         :: Top12
     TYPE(INT_RNK2)         :: Top14
     TYPE(INT_RNK2)         :: Top14_2
     INTEGER                :: I,J,K,L,N,M,II,JJ,NI,NJ,KK,LL
     INTEGER                :: NatmsLoc,NMax14,NMax12,IN12,JN12,KN12
     !
     NMax12=Size(Top12%I,2)-1
     !
     NMax14=30
     K=NMax14+1
     CALL New(Top14,(/NatmsLoc,K/))
     Top14%I(1:NatmsLoc,1:NMax14+1)=0 
     !
     DO II=1,NatmsLoc
       IN12=Top12%I(II,1)
       DO J=1,IN12
         JJ=Top12%I(II,J+1)
         JN12=Top12%I(JJ,1)
         DO K=1,JN12
           KK=Top12%I(JJ,K+1)
           IF(II/=KK) THEN
             KN12=Top12%I(KK,1)
             DO L=1,KN12
               LL=Top12%I(KK,L+1)
               IF(JJ/=LL.AND.II/=LL) THEN
                 !
                 NI=Top14%I(II,1)
                 !
                 ! check matrix Size, increase Size if necessary
                 !
                 IF(NI>=NMax14) THEN
                   NMax14=NMax14+30
                   CALL New(Top14_2,(/NatmsLoc,NMax14+1/))
                   Top14_2%I(1:NatmsLoc,1:NMax14+1)=0 
                   Top14_2%I(1:NatmsLoc,1:NMax14+1-30)=&
                                 Top14%I(1:NatmsLoc,1:NMax14+1-30)
                   CALL Delete(Top14)
                   CALL New(Top14,(/NatmsLoc,NMax14+1/))
                   Top14%I(1:NatmsLoc,1:NMax14+1)=&
                                 Top14_2%I(1:NatmsLoc,1:NMax14+1)
                   CALL Delete(Top14_2)
                 ENDIF
                 !        
                 IF(NI/=0) THEN
                   IF(ANY(Top14%I(II,2:NI+1)==LL)) THEN
                     CYCLE
                   ELSE
                     Top14%I(II,1)=NI+1
                     Top14%I(II,1+(NI+1))=LL
                   ENDIF
                 ELSE
                     Top14%I(II,1)=NI+1
                     Top14%I(II,1+(NI+1))=LL
                 ENDIF
                 !
               ENDIF !!! II/=LL and JJ/=LL
             ENDDO !!! LL
           ENDIF !!! II/=KK
         ENDDO !!!! KK
       ENDDO !!!! JJ
     ENDDO !!!! II
   END SUBROUTINE Topology_14 
!
!--------------------------------------------------------------
!
   SUBROUTINE SORT_INTO_Box1(BoxSize,XYZ,NatmsLoc, &
                       NX,NY,NZ,BXMIN,BYMIN,BZMIN)
     !
     ! sort the Atoms of a molecule into Boxes
     !
     IMPLICIT NONE
     INTEGER :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ,IORD,IADD,NatmsLoc
     REAL(DOUBLE) :: BoxSize,VBIG,XYZ(1:3,NatmsLoc)
     REAL(DOUBLE) :: BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
     !
     CALL BoxBorders(XYZ,BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax)
     !
     NX=INT((BXMax-BXMIN)/BoxSize)+1
     NY=INT((BYMax-BYMIN)/BoxSize)+1
     NZ=INT((BZMax-BZMIN)/BoxSize)+1
     NBox=NX*NY*NZ
     !
   END SUBROUTINE SORT_INTO_Box1
!
!--------------------------------------------------------------
!
   SUBROUTINE SORT_INTO_Box2(BoxSize,XYZ,NatmsLoc, &
     NX,NY,NZ,BXMIN,BYMIN,BZMIN,BoxI1,BoxJ1,ISet)
     !
     ! BoxI(I) : contains the ordering (box) number of the 
     !           first atom of the I-th Box (like in sparse row-wise)
     ! BoxJ(J) : gives the original serial number of the atom desribed 
     !           by the J-th ordering number
     ! XYZ: contains Cartesian coordinates of Atoms
     ! BoxSize: linear Box Size
     !
     IMPLICIT NONE
     INTEGER,OPTIONAL :: ISET 
     TYPE(INT_VECT)   :: BoxI,BoxJ
     TYPE(INT_RNK2)   :: ISignU
     TYPE(INT_RNK3)   :: BoxCounter 
     REAL(DOUBLE)     :: BoxSize,VBIG
     REAL(DOUBLE)     :: BXMIN,BXMax,BYMIN,BYMax,BZMIN,BZMax
     INTEGER          :: IORD,IADD,NatmsLoc
     INTEGER          :: I,J,JJ,NX,NY,NZ,NBox,IX,IY,IZ
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     TYPE(INT_VECT),OPTIONAL      :: BoxI1,BoxJ1
     !
     VBIG=1.D+90 
     NBox=NX*NY*NZ
     CALL New(BoxI,NBox+1)
     CALL New(BoxJ,NatmsLoc)
     CALL New(ISignU,(/2,NatmsLoc/))
     CALL New(BoxCounter,(/NX,NY,NZ/))
     BoxCounter%I(1:NX,1:NY,1:NZ)=0
     BoxI%I(1:NBox+1)=0
     !
     ! Count number of atoms in the box
     !
     DO I=1,NatmsLoc
       IX=INT((XYZ(1,I)-BXMIN)/BoxSize)+1
       IY=INT((XYZ(2,I)-BYMIN)/BoxSize)+1
       IZ=INT((XYZ(3,I)-BZMIN)/BoxSize)+1
       IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY !order parameter: ZXY
       BoxCounter%I(IX,IY,IZ)=BoxCounter%I(IX,IY,IZ)+1
       ISignU%I(1,I)=IORD
       ISignU%I(2,I)=BoxCounter%I(IX,IY,IZ) 
     ENDDO
     !     
     ! Count pointers to individual Boxes within array A
     !     
     BoxI%I(1)=1
     DO IZ=1,NZ
     DO IX=1,NX
     DO IY=1,NY
       IORD=NX*NY*(IZ-1)+NY*(IX-1)+IY
       BoxI%I(IORD+1)=BoxI%I(IORD)+BoxCounter%I(IX,IY,IZ)
     ENDDO
     ENDDO
     ENDDO
     !
     ! Set up contents of Boxes as represented in A, sparse row-wise
     !
     DO I=1,NatmsLoc
       IORD=ISignU%I(1,I)
       IADD=ISignU%I(2,I)
       BoxJ%I(BoxI%I(IORD)-1+IADD)=I
     ENDDO
     !
     IF(PRESENT(BoxI1)) THEN
       BoxI1%I(:)=BoxI%I(:)
     ENDIF
     IF(PRESENT(BoxJ1)) THEN
       BoxJ1%I(:)=BoxJ%I(:)
     ENDIF
     !
     CALL Delete(ISignU)
     CALL Delete(BoxCounter)
     CALL Delete(BoxI)
     CALL Delete(BoxJ)
   END SUBROUTINE SORT_INTO_Box2
!
!----------------------------------------------------------------
!
   SUBROUTINE Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl)
     !
     ! This subroutine merges Topological information
     ! to get the list for Exclusion energy calculation
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)          :: Top_Excl_Out,Top12,Top13,Top14
     TYPE(INT_RNK2)          :: Top_Excl,Top_New
     INTEGER                 :: NatmsLoc,NMax12,NMax13,NMax14
     INTEGER                 :: NMax_Excl,NNew,NOLD
     INTEGER                 :: NMax_Excl_Out,I,J,K,KK,JJ
     !
     NMax12=Size(Top12%I,2)-1
     NMax13=Size(Top13%I,2)-1
     NMax14=Size(Top14%I,2)-1
     !
     ! Initialize Top_Excl
     !
     NMax_Excl=NMax12+NMax13+NMax14
     CALL New(Top_Excl,(/NatmsLoc,NMax_Excl+1/))
     Top_Excl%I(:,:)=0
     Top_Excl%I(1:NatmsLoc,1:NMax12+1)=Top12%I(1:NatmsLoc,1:NMax12+1)
     !
     ! Now merge Topologies, in order to avoid double counting in 
     ! Exclusion energies
     !
     NMax_Excl_Out=0
     DO I=1,NatmsLoc
       !
       NNew=0
       NOLD=Top_Excl%I(I,1)
       DO J=1,Top13%I(I,1)
         JJ=Top13%I(I,J+1)
         IF(ANY(Top_Excl%I(I,2:NOLD+1)==JJ)) THEN
           CYCLE
         ELSE
           NNew=NNew+1
           Top_Excl%I(I,NOLD+1+NNew)=JJ
         ENDIF
       ENDDO 
       NNew=NOLD+NNew
       Top_Excl%I(I,1)=NNew
       IF(NMax_Excl_Out<NNew) NMax_Excl_Out=NNew
       !
       NNew=0
       NOLD=Top_Excl%I(I,1)
       DO J=1,Top14%I(I,1)
         JJ=Top14%I(I,J+1)
         IF(ANY(Top_Excl%I(I,2:NOLD+1)==JJ)) THEN
           CYCLE
         ELSE
           NNew=NNew+1
           Top_Excl%I(I,NOLD+1+NNew)=JJ
         ENDIF
       ENDDO 
       NNew=NOLD+NNew
       Top_Excl%I(I,1)=NNew
       IF(NMax_Excl_Out<NNew) NMax_Excl_Out=NNew
     ENDDO 
   END SUBROUTINE Excl_LIST 
!
!----------------------------------------------------------------
!
   SUBROUTINE Excl_List14(NatmsLoc,Top12,Top13,Top14,Top_Excl)
     !
     ! This subroutine merges Topological information
     ! to get the list for Exclusion energy calculation
     ! of Atoms in 14 distance. From the Top14 list
     ! Top13 and Top12 occurences must be filtered out
     !
     IMPLICIT NONE
     TYPE(INT_RNK2)          :: Top_Excl_Out,Top12,Top13,Top14
     TYPE(INT_RNK2)          :: Top_Excl,Top_New
     INTEGER                 :: NatmsLoc,NMax12,NMax13,NMax14
     INTEGER                 :: NMax_Excl,NNew,NOLD
     INTEGER                 :: NMax_Excl_Out,I,J,K,KK,JJ,NExcl,N12,N13
     !
     NMax12=Size(Top12%I,2)-1
     NMax13=Size(Top13%I,2)-1
     NMax14=Size(Top14%I,2)-1
     !
     ! Initialize Top_Excl to the Size of Top14 and to zero 
     !
     NMax_Excl=NMax14
     CALL New(Top_Excl,(/NatmsLoc,NMax_Excl+1/))
     Top_Excl%I=0
     !
     ! Now merge Topologies, in order to avoid double counting in 
     ! Exclusion energies
     !
     NMax_Excl_Out=0
     DO I=1,NatmsLoc
       DO J=1,Top14%I(I,1)
         NExcl=Top_Excl%I(I,1)
         N12=Top12%I(I,1)
         N13=Top13%I(I,1)
         JJ=Top14%I(I,J+1)
         IF(ANY(Top_Excl%I(I,2:NExcl+1)==JJ).OR. &
            ANY(Top12%I(I,2:N12+1)==JJ).OR. & 
            ANY(Top13%I(I,2:N13+1)==JJ)  ) THEN
           CYCLE
         ELSE
           Top_Excl%I(I,1)=NExcl+1
           Top_Excl%I(I,NExcl+2)=JJ
         ENDIF
       ENDDO 
         NExcl=Top_Excl%I(I,1)
       IF(NMax_Excl_Out<NExcl) NMax_Excl_Out=NExcl
     ENDDO 
   END SUBROUTINE Excl_List14 
!
!--------------------------------------------------------
!
   SUBROUTINE VDWTop(TopVDW,Top12,BondIJ,NBond)
     TYPE(INT_RNK2) :: Top12,TopVDW,BondIJ,TopVDWNew
     TYPE(INT_VECT) :: Addition
     INTEGER        :: NatmsLoc,NBond,I,J,Dim2,AddDim,I1,I2,II1,II2
     !
     NatmsLoc=SIZE(Top12%I,1)
     Dim2=SIZE(Top12%I,2)
     CALL New(TopVDW,(/NatmsLoc,Dim2/))
     TopVDW%I=Top12%I
     CALL New(Addition,NatmsLoc)
     !
     Addition%I=0
     DO I=1,NBond    
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       Addition%I(I1)=Addition%I(I1)+1
       Addition%I(I2)=Addition%I(I2)+1
     ENDDO
     AddDim=MAXVAL(Addition%I)
     !
     CALL New(TopVDWNew,(/NatmsLoc,Dim2+AddDim/))
     TopVDWNew%I=0
     DO I=1,NatmsLoc
       DO J=1,Dim2 ; TopVDWNew%I(I,J)=TopVDW%I(I,J) ; ENDDO
     ENDDO
     DO I=1,NBond
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       II1=TopVDWNew%I(I1,1)+1
       TopVDWNew%I(I1,1)=II1
       TopVDWNew%I(I1,II1+1)=I2
       II2=TopVDWNew%I(I2,1)+1
       TopVDWNew%I(I2,1)=II2
       TopVDWNew%I(I2,II2+1)=I1
     ENDDO
     !
     CALL Delete(TopVDW)
     CALL New(TopVDW,(/NatmsLoc,Dim2+AddDim/))
     TopVDW%I=TopVDWNew%I
     !
     CALL Delete(TopVDWNew)
     CALL Delete(Addition)
   END SUBROUTINE VDWTop
!
!---------------------------------------------------------------------
!
   SUBROUTINE AngleList(AtmB,Bond,TOPS,XYZ,NonCovBend, &
                        Angle,OutP,Cells,IEq)
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(ANGLEDATA)             :: Angle
     TYPE(OUTPDATA)              :: OutP 
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER     ,DIMENSION(:,:) :: Cells
     TYPE(INT_VECT)              :: BondScore
     INTEGER,DIMENSION(3)        :: RefBonds
     INTEGER                     :: I,J,K,L,NAngle,NatmsLoc,AllBond
     INTEGER                     :: I1,I2,I3,NDimens,NOutP
     INTEGER,DIMENSION(:)        :: IEq
     LOGICAL                     :: NonCovBend
     !
     CALL SetDSYEVWork(3)
     NatmsLoc=SIZE(XYZ,2)
     NAngle=0
     DO I=1,NatmsLoc
       AllBond=AtmB%Count%I(I)
       NAngle=NAngle+AllBond*(AllBond-1)/2
     ENDDO
     CALL New(Angle,NAngle)
     CALL New(OutP,3*NatmsLoc)
     !
     NAngle=0
     NOutP=0
     !
     CALL New(BondScore,Bond%N)
     CALL ScoreBond(BondScore%I,Bond)
     DO I=1,NatmsLoc
       IF(TOPS%Tot12%I(I,1)==0) THEN
         IF(ANY(Cells(I,1:3)>1).OR.ANY(Cells(I,1:3)<-1)) CYCLE
         CALL Halt('Atom not bound at '//TRIM(IntToChar(I))//' .')
       ENDIF
       CALL AnglesRef(RefBonds,NDimens,BondScore%I,I,XYZ,AtmB,Bond)
      !CALL AngleGen(I,Angle,NAngle,RefBonds,NDimens,XYZ,AtmB,Bond,TOPS)
       CALL FullAngleGen(I,Angle,NAngle,AtmB,Bond,TOPS,NonCovBend,IEq)
       CALL OutPGen(I,OutP,NOutP,RefBonds,NDimens, &
                    XYZ,AtmB,Bond,TOPS,NonCovBend,IEq)
     ENDDO
     CALL Delete(BondScore)
     CALL UnSetDSYEVWork()
     !
     Angle%N=NAngle
     OutP%N=NOutP
   END SUBROUTINE AngleList
!
!----------------------------------------------------------------------
!
   SUBROUTINE OutPGen(IAt,OutP,NOutP,RefBonds,NDimens, &
                      XYZ,AtmB,Bond,TOPS,NonCovBend,IEq)
     INTEGER                       :: IAt,I,J,K,L,NOutP,NDimens
     TYPE(OUTPDATA)                :: OutP
     INTEGER,DIMENSION(:)          :: RefBonds,IEq
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ
     REAL(DOUBLE)                  :: D,DMax
     TYPE(ATOMBONDS)               :: AtmB
     TYPE(BONDDATA)                :: Bond
     TYPE(TOPOLOGY)                :: TOPS
     TYPE(INT_VECT)                :: Mark
     INTEGER                       :: NDim,IBonds(3),I1,I2,IM
     LOGICAL                       :: ExcludeVDW,NonCovBend
     !
     ExcludeVDW=.NOT.NonCovBend
     IF(TOPS%Cov12%I(IAt,1)==0) ExcludeVDW=.FALSE. !lonelyatom-angles
     IF(NDimens==2.AND.TOPS%Tot12%I(IAt,1)>2) THEN
       NDim=AtmB%Count%I(IAt)
       CALL New(Mark,NDim)
       Mark%I=0
       IBonds=0
       DO L=1,3   
         DMax=1.D99
         DO J=1,NDim
           K=AtmB%Bonds%I(IAt,J)
           IF(ExcludeVDW.AND.Bond%Type%C(K)(1:3)/='COV') CYCLE
           D=Bond%Length%D(K) 
           IF(D<DMax.AND.Mark%I(J)==0) THEN
             IBonds(L)=K
             DMax=D
             IM=J
           ENDIF 
         ENDDO 
         Mark%I(IM)=1
       ENDDO 
       CALL Delete(Mark)
       IF(ANY(IBonds==0)) RETURN !!! no outp
       !
       DO L=1,3
         I1=Bond%IJ%I(1,IBonds(L))
         I2=Bond%IJ%I(2,IBonds(L))
         IF(I1==IAt) THEN
           IBonds(L)=I2
         ELSE
           IBonds(L)=I1
         ENDIF
       ENDDO
       !
       NOutP=NOutP+1
       OutP%IJKL%I(2,NOutP)=IAt
       OutP%IJKL%I(1,NOutP)=IBonds(1)
       OutP%IJKL%I(3,NOutP)=IBonds(2)
       OutP%IJKL%I(4,NOutP)=IBonds(3)
       NOutP=NOutP+1
       OutP%IJKL%I(2,NOutP)=IAt
       OutP%IJKL%I(1,NOutP)=IBonds(3)
       OutP%IJKL%I(3,NOutP)=IBonds(1)
       OutP%IJKL%I(4,NOutP)=IBonds(2)
       NOutP=NOutP+1
       OutP%IJKL%I(2,NOutP)=IAt
       OutP%IJKL%I(1,NOutP)=IBonds(2)
       OutP%IJKL%I(3,NOutP)=IBonds(3)
       OutP%IJKL%I(4,NOutP)=IBonds(1)
     ENDIF
     DO I=1,NOutP
       I1=OutP%IJKL%I(3,I)
       I2=OutP%IJKL%I(4,I)
       IF(IEq(I1)>IEq(I2)) THEN
         OutP%IJKL%I(3,I)=I2
         OutP%IJKL%I(4,I)=I1
       ENDIF
     ENDDO
   END SUBROUTINE OutPGen
!
!---------------------------------------------------------------------
!
   SUBROUTINE AngleGen(I,Angle,NAngle,RefBonds,NDimens, &
                       XYZ,AtmB,Bond,TOPS)
     TYPE(ANGLEDATA)             :: Angle
     TYPE(ATOMBONDS)             :: AtmB
     TYPE(BONDDATA)              :: Bond
     TYPE(TOPOLOGY)              :: TOPS
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER,DIMENSION(3)        :: RefBonds,RefAtms1,RefAtms2
     INTEGER                     :: I,J,K,L,NAngle,NDimens
     INTEGER                     :: B1,B2,I1B1,I2B1,I1B2,I2B2
     INTEGER                     :: TopDim,II1,II2
     !                         
     CALL AtmsRef1(RefAtms1,I,RefBonds,NDimens,Bond,AtmB)
     DO B1=1,TOPS%Tot12%I(I,1)
       II1=TOPS%Tot12%I(I,B1+1)
       CALL AtmsRef2(RefAtms2,RefAtms1,XYZ,I,II1)
       DO B2=1,3
         II2=RefAtms2(B2)
         !
         IF(II2==0.OR.II2==II1) CYCLE
         TopDim=TOPS%Tot12%I(II1,1)
         IF(ANY(TOPS%Tot12%I(II1,2:TopDim+1)==II2)) CYCLE
         IF(ANY(RefAtms1(1:3)==II1).AND.II2>II1) CYCLE
         !
         NAngle=NAngle+1
         Angle%IJK%I(2,NAngle)=I
         IF(II1<II2) THEN
           Angle%IJK%I(1,NAngle)=II1
           Angle%IJK%I(3,NAngle)=II2
         ELSE
           Angle%IJK%I(1,NAngle)=II2
           Angle%IJK%I(3,NAngle)=II1
         ENDIF
         !
       ENDDO
     ENDDO
   END SUBROUTINE AngleGen
!
!-------------------------------------------------------------------
!
   SUBROUTINE FullAngleGen(IAt,Angle,NAngle,AtmB,Bond,TOPS, &
                           NonCovBend,IEq)
     INTEGER                :: IAt,I,J,K,L,NAngle
     TYPE(ANGLEDATA)        :: Angle
     TYPE(ATOMBONDS)        :: AtmB
     TYPE(BONDDATA)         :: Bond 
     TYPE(TOPOLOGY)         :: TOPS
     INTEGER                :: B1,B2,I1B1,I2B1,I1B2,I2B2,II1,II2,TopDim
     INTEGER,DIMENSION(:)   :: IEq
     LOGICAL                :: ExcludeVDW,NonCovBend
     !
     ExcludeVDW=.NOT.NonCovBend
     IF(TOPS%Cov12%I(IAt,1)==0) ExcludeVDW=.FALSE. !lonelyatom-angles
     DO J=1,AtmB%Count%I(IAt)
       B1=AtmB%Bonds%I(IAt,J)
       IF(ExcludeVDW.AND.Bond%Type%C(B1)(1:3)/='COV') CYCLE
       I1B1=Bond%IJ%I(1,B1)
       I2B1=Bond%IJ%I(2,B1)
       IF(I1B1==IAt) THEN
         II1=I2B1
       ELSE
         II1=I1B1
       ENDIF
       DO K=J+1,AtmB%Count%I(IAt)
         B2=AtmB%Bonds%I(IAt,K) 
         IF(ExcludeVDW.AND.Bond%Type%C(B2)(1:3)/='COV') CYCLE
         I1B2=Bond%IJ%I(1,B2)
         I2B2=Bond%IJ%I(2,B2)
         IF(I1B2==IAt) THEN
           II2=I2B2
         ELSE
           II2=I1B2
         ENDIF
         TopDim=TOPS%Tot12%I(II1,1)
         IF(ANY(TOPS%Tot12%I(II1,2:TopDim+1)==II2)) CYCLE
         IF(II1==II2) CYCLE
         NAngle=NAngle+1
         Angle%IJK%I(2,NAngle)=IAt
         IF(IEq(II1)<IEq(II2)) THEN
           Angle%IJK%I(1,NAngle)=II1
           Angle%IJK%I(3,NAngle)=II2
         ELSE
           Angle%IJK%I(1,NAngle)=II2
           Angle%IJK%I(3,NAngle)=II1
         ENDIF
       ENDDO 
     ENDDO 
   END SUBROUTINE FullAngleGen
!
!-------------------------------------------------------------------
!
   SUBROUTINE AtmsRef1(RefAtms1,I,RefBonds,NDimens,Bond,AtmB)
     TYPE(ATOMBONDS)      :: AtmB
     TYPE(BONDDATA)       :: Bond
     INTEGER,DIMENSION(3) :: RefAtms1,RefBonds
     INTEGER              :: B2,L,I1B2,I2B2,I,II2,NDimens,NRef
     RefAtms1=0
     NRef=0
     DO B2=1,3
       L=RefBonds(B2)
       IF(L==0) CYCLE
       L=AtmB%Bonds%I(I,L)
       I1B2=Bond%IJ%I(1,L)
       I2B2=Bond%IJ%I(2,L)
       IF(I1B2==I) THEN
         II2=I2B2
       ELSE
         II2=I1B2
       ENDIF
       NRef=NRef+1
       RefAtms1(B2)=II2
       IF(NDimens<=2.AND.NRef==1) EXIT
     ENDDO
   END SUBROUTINE AtmsRef1
!
!--------------------------------------------------------------------
!
   SUBROUTINE AtmsRef2(RefAtms2,RefAtms1,XYZ,I,II1)
     INTEGER,DIMENSION(3)        :: RefAtms2,RefAtms1
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(3)   :: Vect1,Vect2,Vect3,CrossProds
     INTEGER                     :: I,II1,II2,B2,IMin
     !
     RefAtms2=RefAtms1
     CrossProds=Zero
     Vect1=XYZ(:,II1)-XYZ(:,I)
     DO B2=1,3
       II2=RefAtms2(B2)
       IF(II2==0) CYCLE
       Vect2=XYZ(:,II2)-XYZ(:,I)
       Vect2=vect2/SQRT(DOT_PRODUCT(Vect2,Vect2))
       CALL CROSS_PRODUCT(Vect1,Vect2,Vect3) 
       CrossProds(B2)=DOT_PRODUCT(Vect3,Vect3)
     ENDDO
     IMin=1
     DO B2=2,3
       IF(CrossProds(B2)<CrossProds(IMin)) IMin=B2
     ENDDO
     RefAtms2(IMin)=0
   END SUBROUTINE AtmsRef2
!
!---------------------------------------------------------------------
!
   SUBROUTINE AnglesRef(RefBonds,NDimens,BondScore,I1,XYZ,AtmB,Bond)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER     ,DIMENSION(:)   :: BondScore
     INTEGER     ,DIMENSION(3)   :: RefBonds 
     REAL(DOUBLE),DIMENSION(3)   :: EigVals
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB 
     TYPE(DBL_RNK2)              :: BondVect2
     TYPE(INT_VECT)              :: Mark       
     REAL(DOUBLE)                :: AMax,D,CondNum1
     INTEGER                     :: I1,I,J,NDim,IMax,NDimens
     ! 
     RefBonds=0     
     CondNum1=0.01D0 ! This parameter controls recognition of flatness of molecule in the vicinity of atom I1
     NDim=AtmB%Count%I(I1)
     CALL New(BondVect2,(/NDim,3/))
     CALL New(Mark,NDim)
     Mark%I=0
     CALL D3Bonds(I1,XYZ,BondVect2%D,BondScore,AtmB,Bond,Eig_O=EigVals)
     NDimens=3
     DO I=1,3
       IF(EigVals(I)<CondNum1) THEN
         NDimens=NDimens-1
         CYCLE
       ENDIF
       AMax=Zero
       IMax=0
       DO J=1,NDim
         D=ABS(BondVect2%D(J,I))
         IF(D>AMax.AND.Mark%I(J)==0) THEN
           AMax=D 
           IMax=J
         ENDIF
       ENDDO
       RefBonds(I)=IMax  
       Mark%I(IMax)=1
     ENDDO
     CALL Delete(Mark)
     CALL Delete(BondVect2)
   END SUBROUTINE AnglesRef
!
!---------------------------------------------------------------------
!
   SUBROUTINE MergeBonds(BondCov,BondVDW,BondTot)
     TYPE(BONDDATA)  :: BondCov,BondVDW,BondTot
     INTEGER         :: NTot,NCov
     INTEGER         :: I,J,K,L
     !
     NTot=BondCov%N+BondVDW%N
     CALL New(BondTot,NTot)
     BondTot%N=NTot
     DO I=1,BondCov%N
       CALL Set_Bond_EQ_Bond(BondTot,I,BondCov,I)
     ENDDO
     K=BondCov%N
     DO I=1,BondVDW%N
       CALL Set_Bond_EQ_Bond(BondTot,K+I,BondVDW,I)
     ENDDO
   END SUBROUTINE MergeBonds
!
!---------------------------------------------------------------------
!
   SUBROUTINE MergeAtmB(AtmBCov,AtmBVDW,NCov,AtmBTot)   
     TYPE(ATOMBONDS) :: AtmBCov,AtmBVDW,AtmBTot 
     INTEGER         :: NatmsLoc,I,J,K,MaxB,NCov
     !
     NatmsLoc=SIZE(AtmBCov%Count%I)
     AtmBTot%N1=NatmsLoc
     CALL New(AtmBTot%Count,NatmsLoc)
     DO I=1,NatmsLoc
       AtmBTot%Count%I(I)=AtmBCov%Count%I(I)+AtmBVDW%Count%I(I)
     ENDDO
     MaxB=MAXVAL(AtmBTot%Count%I)
     CALL New(AtmBTot%Bonds,(/NatmsLoc,MaxB/))
     CALL New(AtmBTot%Atoms,(/NatmsLoc,MaxB/))
     AtmBTot%Bonds%I=0
     DO I=1,NatmsLoc
       K=AtmBCov%Count%I(I)
       DO J=1,K
         AtmBTot%Bonds%I(I,J)=AtmBCov%Bonds%I(I,J)
         AtmBTot%Atoms%I(I,J)=AtmBCov%Atoms%I(I,J)
       ENDDO
       DO J=1,AtmBVDW%Count%I(I)
         AtmBTot%Bonds%I(I,K+J)=NCov+AtmBVDW%Bonds%I(I,J)
         AtmBTot%Atoms%I(I,K+J)=NCov+AtmBVDW%Atoms%I(I,J)
       ENDDO
     ENDDO
   END SUBROUTINE MergeAtmB
!
!---------------------------------------------------------------------
!
   SUBROUTINE ArchiveTop(TOPS,ArchMem,Bond,AtmB,HFileIn,iCLONE,iGEO)
     TYPE(TOPOLOGY)      :: TOPS
     TYPE(BONDDATA)      :: Bond,Bond2
     TYPE(ATOMBONDS)     :: AtmB,AtmB2
     CHARACTER(LEN=*)    :: HFileIn
     INTEGER             :: iCLONE,iGEO,HDFFileID,ArchMem
     INTEGER             :: I,J,NZ,IBack,NatmsLoc
     CHARACTER(LEN=DCL)  :: Tag
     !
     IBack=MIN(ArchMem-1,iGEO-1)
     NatmsLoc=SIZE(AtmB%Count%I)
     HDFFileID=OpenHDF(HFileIn)
     HDF_CurrentID= &
       OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
     !
     DO I=1,IBack
       Tag=TRIM(IntToChar(iGEO-I))
       CALL Get(Bond2,'Bond',Tag)
       CALL Get(AtmB2,'AtmB',Tag)
       IF(Bond2%N>0) THEN
         CALL MergeBondSets(Bond,Bond2,AtmB,AtmB2)
         CALL Delete(AtmB)
         CALL SortBonds(NatmsLoc,AtmB,Bond)
       ENDIF
       !
       CALL Delete(Bond2)
       CALL Delete(AtmB2)
     ENDDO
     !
     IF(IBack/=0) THEN   
       CALL Delete(TOPS%Tot12)
       CALL Delete(TOPS%Tot13)
       CALL Delete(TOPS%Tot14)
       CALL Delete(TOPS%TotExcl)
       CALL Topology_12(AtmB,TOPS%Tot12)
       CALL Topology_13(NatmsLoc,TOPS%Tot12,TOPS%Tot13)
       CALL Topology_14(NatmsLoc,TOPS%Tot12,TOPS%Tot14)
       CALL Excl_List(NatmsLoc,TOPS%Tot12,TOPS%Tot13,TOPS%Tot14, &
                      TOPS%TotExcl)
     ENDIF
     !
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
   END SUBROUTINE ArchiveTop
! 
!---------------------------------------------------------------------
!
   SUBROUTINE MergeBondSets(Bond,Bond2,AtmB,AtmB2)
     TYPE(BONDDATA) :: Bond,Bond2,BondNew,BondM
     TYPE(ATOMBONDS):: AtmB,AtmB2
     TYPE(INT_VECT) :: AllAtms,BondMark
     INTEGER        :: NatmsLoc,NNew,K,L,I,J
     ! 
     NatmsLoc=AtmB2%N1
     CALL New(AllAtms,NatmsLoc)
     AllAtms%I=0
     CALL New(BondNew,Bond2%N)
     CALL New(BondMark,Bond2%N)
     BondMark%I=0
     !
     NNew=0
     DO I=1,NatmsLoc
       DO J=1,AtmB%Count%I(I)
         L=AtmB%Atoms%I(I,J)
         AllAtms%I(L)=1 
       ENDDO 
       !
       DO J=1,AtmB2%Count%I(I)
         L=AtmB2%Atoms%I(I,J)
         IF(AllAtms%I(L)==0) THEN
           K=AtmB2%Bonds%I(I,J)
           IF(BondMark%I(K)==0) THEN
             NNew=NNew+1
             BondMark%I(K)=1
             CALL Set_Bond_EQ_Bond(BondNew,NNew,Bond2,K)
           ENDIF
         ENDIF
       ENDDO 
       !
       DO J=1,AtmB%Count%I(I)
         L=AtmB%Atoms%I(I,J)
         AllAtms%I(L)=0 
       ENDDO 
     ENDDO
     BondNew%N=NNew
     IF(NNew/=0) THEN
       CALL MergeBonds(Bond,BondNew,BondM)
       CALL Delete(Bond)
       CALL Set_BONDDATA_EQ_BONDDATA(Bond,BondM)
       CALL Delete(BondM)      
     ENDIF
     CALL Delete(BondMark)      
     CALL Delete(BondNew)      
     CALL Delete(AllAtms)      
   END SUBROUTINE MergeBondSets
!
!---------------------------------------------------------------------
!
   SUBROUTINE CutOffDispl(Displ,IntCs,StreCritIn,AngleCritIn,Fact_O)
     REAL(DOUBLE),DIMENSION(:) :: Displ
     REAL(DOUBLE)              :: StreCritIn,AngleCritIn
     TYPE(INTC)                :: IntCs
     INTEGER                   :: I
     REAL(DOUBLE)              :: Fact  
     REAL(DOUBLE),OPTIONAL     :: Fact_O
     !
     Fact=One
     IF(PRESENT(Fact_O)) Fact=Fact_O
     DO I=1,IntCs%N
       IF(HasLattice(IntCs%Def%C(I))) CYCLE
       CALL CtrlDispl(IntCs%Def%C(I),Displ(I),Fact, &
                      StreCritIn,AngleCritIn)
     ENDDO
   END SUBROUTINE CutOffDispl
!
!-------------------------------------------------------------------
!
   SUBROUTINE CtrlDispl(Def,Displ,Fact,StreCritIn,AngleCritIn)
     REAL(DOUBLE)     :: Displ,StreCrit,AngleCrit
     REAL(DOUBLE)     :: StreCritIn,AngleCritIn
     REAL(DOUBLE)     :: Fact
     CHARACTER(LEN=*) :: Def
     !
     IF(Def(1:4)=='STRE') THEN
       StreCrit=StreCritIn*Fact
       IF(ABS(Displ)>StreCrit) Displ=SIGN(StreCrit,Displ)
     ELSE
       AngleCrit=AngleCritIn*Fact
       IF(ABS(Displ)>AngleCrit) Displ=SIGN(AngleCrit,Displ)
     ENDIF
   END SUBROUTINE CtrlDispl
!
!-------------------------------------------------------------------
!
   SUBROUTINE CtrlRange(Displ,Range,NDim)
     REAL(DOUBLE) :: Displ,Range(:),StepMax
     INTEGER      :: NDim
     !
   ! IF(NDim<=3) THEN
   !   StepMax=2.D0*Range
   !  !StepMax=Range
   ! ELSE 
       StepMax=Range(2)
   ! ENDIF
     IF(ABS(Displ)>StepMax) Displ=SIGN(StepMax,Displ)
   END SUBROUTINE CtrlRange
!
!-------------------------------------------------------------------
! 
   SUBROUTINE SetBackToRefs(IntCValues,IntCs,RefPoints)
     REAL(DOUBLE),DIMENSION(:)   :: IntCValues
     REAL(DOUBLE),DIMENSION(:)   :: RefPoints
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,I,J
     REAL(DOUBLE)                :: Center
     !
     NIntC=SIZE(IntCValues,1)
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='LINB'.OR. &
          IntCs%Def%C(I)(1:4)=='OUTP'.OR. &
          IntCs%Def%C(I)(1:4)=='TORS') THEN
          Center=RefPoints(I)
          CALL PAngle1(Center,IntCValues(I))
       ENDIF
     ENDDO
     !
   END SUBROUTINE SetBackToRefs
!
!-------------------------------------------------------------------
! 
   SUBROUTINE IntCsConstr(IntCs,IntCsX,DoReturn)
     TYPE(INTC) :: IntCs,IntCsX
     INTEGER    :: I,J,NC
     LOGICAL    :: DoReturn
     !
     DoReturn=.FALSE.
     NC=0
     DO I=1,IntCs%N
       IF(IntCs%Constraint%L(I)) NC=NC+1
     ENDDO
     IF(NC==0) DoReturn=.TRUE.
     IF(DoReturn) RETURN
     CALL New(IntCsX,NC)
     NC=0
     DO I=1,IntCs%N
       IF(IntCs%Constraint%L(I)) THEN
         NC=NC+1
         CALL Set_INTC_EQ_INTC(Intcs,IntcsX,I,I,NC)
       ENDIF
     ENDDO
   END SUBROUTINE IntCsConstr
!
!-------------------------------------------------------------------
! 
   SUBROUTINE ConnectFragments(XYZ,AtNum,BondF,TOPS,Cells,IEq)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(INT_VECT)              :: ITop,JTop,Perm,IPerm
     TYPE(INT_VECT)              :: Center,IFrag,Mark,FragID
     TYPE(INT_VECT)              :: LastFrag,JFrag,FragMark
     TYPE(DBL_VECT)              :: ATop
     TYPE(DBL_RNK2)              :: CenterFrag
     INTEGER                     :: NZ,NAtmsLoc,I,J,II,K,L,M,N,MM,NN,JJ
     INTEGER                     :: ColMax,ISt,IStN
     INTEGER,DIMENSION(:)        :: AtNum,IEq
     INTEGER,DIMENSION(:,:)      :: Cells
     TYPE(BONDDATA)              :: BondF
     REAL(DOUBLE),DIMENSION(3)   :: Vect,XYZII,XYZJJ
     REAL(DOUBLE)                :: Length,Dist,Dist0
     INTEGER,DIMENSION(3)        :: CellAux2
     LOGICAL                     :: DoC,Longer
     ! 
     ! Calculate sparse topology matrix
     ! 
     NatmsLoc=SIZE(XYZ,2)
     CALL TopToSp1x1(TOPS%Tot12%I,ITop,JTop)
     NZ=SIZE(JTop%I)
     CALL New(FragID,NatmsLoc)
     CALL New(ATop,NZ)
     !
     ! Calculate permutations, which make the tightest ordering
     !
     CALL New(Perm,NatmsLoc)
     CALL New(IPerm,NatmsLoc)
     CALL RCMOrder(IPerm%I,Perm%I,NatmsLoc,ITop%I,JTop%I)
     CALL Perm1x1(Perm%I,ITop%I,JTop%I,ATop%D,Symb_O=.TRUE.)
     !
     ! Now, identify separate blocks in Top
     ! It is easy, based on an ordered representation of the topology
     !
     FragID%I=0 
     II=1
     M=IPerm%I(1)
     FragID%I(M)=II
     K=ITop%I(1)
     L=ITop%I(2)
     ColMax=MAX(1,MAXVAL(JTop%I(K:L-1)))
     DO I=2,NatmsLoc
       K=ITop%I(I)
       L=ITop%I(I+1)
       M=IPerm%I(I)
       IF(ColMax<I) THEN
         II=II+1
       ENDIF
       FragID%I(M)=II
       ColMax=MAX(ColMax,MAXVAL(JTop%I(K:L-1)))
       ColMax=MAX(ColMax,I)
     ENDDO 
     !
     ! Sort fragment atoms
     !
     CALL New(IFrag,II+1)
     CALL New(CenterFrag,(/3,II/))
     !
     IFrag%I=0
     CenterFrag%D=Zero
     DO I=1,NatmsLoc
       M=FragID%I(I)
       IFrag%I(M+1)=IFrag%I(M+1)+1
       CenterFrag%D(1:3,M)=CenterFrag%D(1:3,M)+XYZ(1:3,I)
     ENDDO
     !
     DO I=1,II
       CenterFrag%D(1:3,I)=CenterFrag%D(1:3,I)/DBLE(IFrag%I(I+1))
     ENDDO
     !
     IFrag%I(1)=1
     DO I=1,II
       IFrag%I(I+1)=IFrag%I(I)+IFrag%I(I+1) 
     ENDDO
     !
     CALL New(LastFrag,II)
     CALL New(FragMark,II)
     FragMark%I=0
     CALL New(JFrag,NatmsLoc)
     LastFrag%I=0
     DO I=1,NatmsLoc
       M=FragID%I(I)
       K=IFrag%I(M)+LastFrag%I(M)
       JFrag%I(K)=I
       LastFrag%I(M)=LastFrag%I(M)+1
     ENDDO
     !
     ! Connect fragments
     !
     IF(II>1) THEN
       WRITE(*,*) 'NUMBER OF FRAGMENTS OBSERVED= ',II
       WRITE(Out,*) 'NUMBER OF FRAGMENTS OBSERVED= ',II
     ENDIF
     K=II*(II-1)
     CALL New(BondF,K)
     K=0
     FragMark%I(1)=1
     DO I=1,II-2
       IF(FragMark%I(I)/=0) CYCLE ! only unconnected fragments pass
       !
       MM=0
       Dist=1.D99
       DO J=1,II
         IF(FragMark%I(J)==0) CYCLE ! only already connected frags pass
         XYZII=CenterFrag%D(1:3,I)
         XYZJJ=CenterFrag%D(1:3,J)
         CALL ClosestAtms(M,N,I,J,IFrag%I,JFrag%I, &
                          XYZII,XYZJJ,XYZ)
         Vect=XYZ(1:3,M)-XYZ(1:3,N)
         Dist0=SQRT(DOT_PRODUCT(Vect,Vect))
         IF(Dist0<Dist) THEN
           Dist=Dist0
           MM=M
           NN=N
           JJ=J
         ENDIF
       ENDDO
       IF(MM==0) CALL Halt('Error in fragment recognition')
       M=MM
       N=NN
       J=JJ
       FragMark%I(I)=1
       !
       ! connect only to central cell
       IF(HasCentralCell(I,IFrag%I,JFrag%I,IEq)) THEN
         IF(ANY(Cells(N,1:3)>1).OR.ANY(Cells(N,1:3)<-1)) CYCLE 
       ELSE IF(HasCentralCell(J,IFrag%I,JFrag%I,IEq)) THEN
         IF(ANY(Cells(M,1:3)>1).OR.ANY(Cells(M,1:3)<-1)) CYCLE 
       ELSE
         CYCLE
       ENDIF
       K=K+1
       Vect=XYZ(1:3,M)-XYZ(1:3,N)
       Length=SQRT(DOT_PRODUCT(Vect,Vect))
       IF(IEq(M)<IEq(N)) THEN
         BondF%IJ%I(1:2,K)=(/M,N/)
       ELSE
         BondF%IJ%I(1:2,K)=(/N,M/)
       ENDIF
       BondF%Length%D(K)=Length
       BondF%Type%C(K)(1:4)='Frag'
     ENDDO
     IF(K==0) THEN
       CALL Delete(BondF)
     ELSE
       BondF%N=K
     ENDIF
     !
     CALL Delete(CenterFrag)
     CALL Delete(FragID)
     CALL Delete(FragMark)
     CALL Delete(IFrag)
     CALL Delete(LastFrag)
     CALL Delete(JFrag)
     CALL Delete(ITop)
     CALL Delete(JTop)
     CALL Delete(ATop)
     CALL Delete(Perm)
     CALL Delete(IPerm)
   END SUBROUTINE ConnectFragments
!
!-------------------------------------------------------------------
! 
   SUBROUTINE SortFragments(Top12,FragID,NFrag)
     INTEGER,DIMENSION(:,:)      :: Top12
     TYPE(INT_VECT)              :: ITop,JTop,Perm,IPerm
     INTEGER,DIMENSION(:)        :: FragID
     TYPE(INT_VECT)              :: LastFrag
     TYPE(DBL_VECT)              :: ATop
     INTEGER                     :: NZ,NAtmsLoc,I,J,NFrag,K,L,M,N
     INTEGER                     :: ColMax
     ! 
     ! Calculate sparse topology matrix
     ! 
     NatmsLoc=SIZE(Top12,1)
     CALL TopToSp1x1(Top12,ITop,JTop)
     NZ=SIZE(JTop%I)
     CALL New(ATop,NZ)
     !
     ! Calculate permutations, which make the tightest ordering
     !
     CALL New(Perm,NatmsLoc)
     CALL New(IPerm,NatmsLoc)
     CALL RCMOrder(IPerm%I,Perm%I,NatmsLoc,ITop%I,JTop%I)
     CALL Perm1x1(Perm%I,ITop%I,JTop%I,ATop%D,Symb_O=.TRUE.)
     !
     ! Now, identify separate blocks in Top
     ! It is easy, based on an ordered representation of the topology
     !
     FragID=0 
     NFrag=1
     M=IPerm%I(1)
     FragID(M)=NFrag
     K=ITop%I(1)
     L=ITop%I(2)
     ColMax=MAX(1,MAXVAL(JTop%I(K:L-1)))
     DO I=2,NatmsLoc
       K=ITop%I(I)
       L=ITop%I(I+1)
       M=IPerm%I(I)
       IF(ColMax<I) THEN
         NFrag=NFrag+1
       ENDIF
       FragID(M)=NFrag
       ColMax=MAX(ColMax,MAXVAL(JTop%I(K:L-1)))
       ColMax=MAX(ColMax,I)
     ENDDO 
     CALL Delete(ITop)
     CALL Delete(JTop)
     CALL Delete(ATop)
     CALL Delete(Perm)
     CALL Delete(IPerm)
   END SUBROUTINE SortFragments
!
!-------------------------------------------------------------------
! 
   SUBROUTINE FragIJs(FragID,NFrag,IFrag,JFrag)
     TYPE(INT_VECT)       :: IFrag,JFrag,LastFrag
     INTEGER              :: NFrag,NatmsLoc,I,M,K
     INTEGER,DIMENSION(:) :: FragID
     !
     ! Sort fragment atoms
     !
     NatmsLoc=SIZE(FragID)
     CALL New(IFrag,NFrag+1)
     CALL New(JFrag,NatmsLoc)
     !
     IFrag%I=0
     DO I=1,NatmsLoc
       M=FragID(I)
       IFrag%I(M+1)=IFrag%I(M+1)+1
     ENDDO
     IFrag%I(1)=1
     DO I=1,NFrag
       IFrag%I(I+1)=IFrag%I(I)+IFrag%I(I+1) 
     ENDDO
     !
     CALL New(LastFrag,NFrag)
     LastFrag%I=0
     DO I=1,NatmsLoc
       M=FragID(I)
       K=IFrag%I(M)+LastFrag%I(M)
       JFrag%I(K)=I
       LastFrag%I(M)=LastFrag%I(M)+1
     ENDDO
     CALL Delete(LastFrag)
   END SUBROUTINE FragIJs
!
!-------------------------------------------------------------------
! 
   FUNCTION HasCentralCell(II,IFrag,JFrag,IEq)
     LOGICAL :: HasCentralCell
     INTEGER :: II,I,J
     INTEGER,DIMENSION(:) :: IFrag,JFrag,IEq
     !
     HasCentralCell=.FALSE.
     DO I=IFrag(II),IFrag(II+1)-1
       J=JFrag(I)
       IF(IEq(J)==J) THEN
         HasCentralCell=.TRUE.
         EXIT
       ENDIF
     ENDDO
   END FUNCTION HasCentralCell
!
!-------------------------------------------------------------------
! 
   SUBROUTINE ClosestAtms(MM,NN,II,JJ,IFrag,JFrag, &
                          XYZII,XYZJJ,XYZ)
     INTEGER                     :: II,JJ,M,N,I,J,MaxIt,It,K,L
     INTEGER                     :: MMOld,NNOld,MM,NN
     INTEGER,DIMENSION(:)        :: IFrag,JFrag
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(3)   :: XYZII,XYZJJ
     REAL(DOUBLE)                :: DistM,DistN,Dist,Length
     !
     ! Find the closest atoms from fragments I and J
     ! MM is the bridge-head from fragm II, NN is from JJ.
     MMOld=0 
     NNOld=0 
     MaxIt=MAX(IFrag(II+1)-IFrag(II),IFrag(JJ+1)-IFrag(JJ))
     DO It=1,Maxit
       CALL FragFrag(MM,XYZJJ,II,IFrag,JFrag,XYZ)
       XYZII=XYZ(1:3,MM) 
       CALL FragFrag(NN,XYZII,JJ,IFrag,JFrag,XYZ)
       XYZJJ=XYZ(1:3,NN) 
       IF(MMOld==MM.AND.NNOld==NN) EXIT
       MMOld=MM
       NNOld=NN
     ENDDO 
   END SUBROUTINE ClosestAtms
!
!----------------------------------------------------------------
!
   SUBROUTINE FragFrag(MM,XYZJJ,II,IFrag,JFrag,XYZ)
     INTEGER                     :: M,N,II,MM,K
     INTEGER,DIMENSION(:)        :: IFrag,JFrag
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(3)   :: DistV,XYZJJ
     REAL(DOUBLE)                :: DistM,Dist
     !
     DistM=1.D99
     DO K=IFrag(II),IFrag(II+1)-1 
       M=JFrag(K)
       DistV=XYZ(1:3,M)-XYZJJ
       Dist=SQRT(DOT_PRODUCT(DistV,DistV))
       IF(Dist<DistM) THEN
         MM=M
         DistM=Dist
       ENDIF
     ENDDO 
   END SUBROUTINE FragFrag
!
!-------------------------------------------------------------------
! 
   SUBROUTINE GcInvIter(CartVect,ISpB,JSpB,ASpB,CholData,NIntC) 
     REAL(DOUBLE),DIMENSION(:) :: CartVect
     TYPE(INT_VECT)            :: ISpB,JSpB
     TYPE(DBL_VECT)            :: ASpB,CartVect2,cartVect3,CartVect4,IntVect2
     TYPE(Cholesky)            :: CholData
     INTEGER                   :: I,J,NIntC,NCart,MaxIter
     REAL(DOUBLE)              :: Fact,Tol 
     CHARACTER(LEN=DCL)        :: Messg
     !
     MaxIter=50
     Tol=1.D-8
     NCart=SIZE(CartVect)
     CALL New(CartVect2,NCart)
     CALL New(CartVect3,NCart)
     CALL New(CartVect4,NCart)
     CALL New(IntVect2,NIntC)
     !
     Cartvect2%D=Zero
     DO I=1,MaxIter
       CALL CALC_BxVect(ISpB,JSpB,ASpB,IntVect2%D,CartVect2%D)
       CALL CALC_BxVect(ISpB,JSpB,ASpB,IntVect2%D,CartVect3%D,Trp_O=.TRUE.)
       CartVect3%D=CartVect-CartVect3%D
       CALL CALC_GcInvCartV(CholData,CartVect3%D,CartVect4%D) 
       CartVect2%D=Cartvect2%D+CartVect4%D 
       Fact=DOT_PRODUCT(CartVect4%D,CartVect4%D)
       Fact=SQRT(Fact/DBLE(NCart))
       IF(Fact<Tol) EXIT
     ENDDO 
     IF(Fact>Tol) THEN
       Messg='WARNING! '// &
             'Gc V = V0 has not been solved exactly , Fact= '// &
              TRIM(DblToChar(Fact))//' Tol= '//TRIM(DblToChar(Tol))
       WRITE(*,*) Messg
       WRITE(Out,*) Messg
     ENDIF
     CartVect=CartVect2%D
     !
     CALL Delete(CartVect4)
     CALL Delete(CartVect3)
     CALL Delete(CartVect2)
     CALL Delete(Intvect2)
   END SUBROUTINE GcInvIter
!
!-------------------------------------------------------------------
! 
   SUBROUTINE GiInvIter(ISpB,JSpB,ASpB,CholData,IntA1,IntA2, &
                        NCart,NIntC)
     TYPE(INT_VECT)            :: ISpB,JSpB
     TYPE(DBL_VECT)            :: ASpB
     TYPE(Cholesky)            :: CholData
     REAL(DOUBLE),DIMENSION(:) :: IntA1,IntA2
     TYPE(DBL_VECT)            :: IntA3,IntA4,CartA1
     REAL(DOUBLE)              :: Tol,Fact
     CHARACTER(LEN=DCL)        :: Messg
     INTEGER                   :: II,I,JJ,J,MaxStep,NIntC,NCart
     !
     MaxStep=50 
     Tol=1.D-8
     CALL New(IntA3,NIntC)
     CALL New(IntA4,NIntC)
     CALL New(CartA1,NCart)
     !
     IntA2=Zero
     DO II=1,MaxStep 
       CALL CALC_BxVect(ISpB,JSpB,ASpB,IntA2,CartA1%D,Trp_O=.TRUE.)
       CALL CALC_BxVect(ISpB,JSpB,ASpB,IntA3%D,CartA1%D)
       IntA3%D=IntA1-IntA3%D
       CALL CALC_GcInvCartV(CholData,IntA3%D,IntA4%D)
       IntA2=IntA2+IntA4%D 
       Fact=DOT_PRODUCT(IntA4%D,IntA4%D)
       Fact=SQRT(Fact/DBLE(NIntC))
       IF(Fact<Tol) EXIT
     ENDDO 
     IF(Fact>Tol) THEN
       Messg='WARNING! '// &
             'Unsuccesful projection of constraint forces, Fact= '// &
              TRIM(DblToChar(Fact))//' Tol= '//TRIM(DblToChar(Tol))
       WRITE(*,*) Messg
       WRITE(Out,*) Messg
     ENDIF
     !
     CALL Delete(IntA4)
     CALL Delete(IntA3)
     CALL Delete(CartA1)
   END SUBROUTINE GiInvIter
!
!---------------------------------------------------------------------
!
   SUBROUTINE POffHardGc(Displ,ISpB,JSpB,ASpB,CholData, &
                         NCart,NIntC,Fact)
     REAL(DOUBLE),DIMENSION(:)  :: Displ
     TYPE(Cholesky)             :: CholData
     TYPE(INT_VECT)             :: ISpB,JSpB
     TYPE(DBL_VECT)             :: ASpB
     TYPE(DBL_VECT)             :: IntA1,CartA1,Carts,IntA2
     INTEGER                    :: I,J,NCart,NIntC
     REAL(DOUBLE)               :: Fact
     !
     CALL New(IntA1,NIntC)
     CALL New(CartA1,NCart)
     IntA1%D=Displ
     CALL CALC_BxVect(ISpB,JSpB,ASpB,IntA1%D,CartA1%D,Trp_O=.TRUE.)
     CALL GcInvIter(CartA1%D,ISpB,JSpB,ASpB,CholData,NIntC)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,IntA1%D,CartA1%D)
     !
     Fact=DOT_PRODUCT(IntA1%D,Displ)/DOT_PRODUCT(Displ,Displ)
     Fact=Fact*100.D0
     Displ=Displ-IntA1%D
     CALL Delete(IntA1)
     CALL Delete(CartA1)
   END SUBROUTINE POffHardGc
!
!-------------------------------------------------------------------
!
   SUBROUTINE CleanConstrCart(XYZ,PBCDim,CartGrad,GOpt,SCRPath)
     TYPE(INTC)                 :: IntCs,IntCsX
     REAL(DOUBLE),DIMENSION(:)  :: CartGrad
     REAL(DOUBLE)               :: BoxShapeT(3,3),InvBoxSh(3,3)
     REAL(DOUBLE)               :: BoxShape(3,3),Vect1(3),Vect2(3)
     TYPE(GeomOpt)              :: GOpt
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     INTEGER                    :: Print
     CHARACTER(LEN=*)           :: SCRPath
     TYPE(Cholesky)             :: CholData
     TYPE(INT_VECT)             :: ISpB,JSpB
     TYPE(DBL_VECT)             :: ASpB
     TYPE(DBL_VECT)             :: CartA1,IntA1,IntA2
     INTEGER                    :: II,I,JJ,J,NCart,PBCDim,NatmsLoc,K,L
     LOGICAL                    :: DoReturn
     REAL(DOUBLE)               :: Fact
     !
     ! Generate INTC of Constraints
     !
     CALL CleanLattGradRatio(PBCDim,XYZ,GOpt%Constr,CartGrad)
     CALL IntCsConstr(GOpt%ExtIntCs,IntCsX,DoReturn)
     IF(DoReturn) RETURN
     !
     IntCsX%Constraint%L=.FALSE.
     IntCsX%Active%L=.TRUE.
     NCart=SIZE(CartGrad)
     NatmsLoc=SIZE(XYZ,2)
     IF(PBCDim>0) THEN
       DO I=1,3
         BoxShapeT(I,1:3)=XYZ(1:3,NatmsLoc-3+I)
         BoxShape(1:3,I)=XYZ(1:3,NatmsLoc-3+I)
       ENDDO
       InvBoxSh=InverseMatrix(BoxShapeT)
       DO I=1,NatmsLoc-3
         K=3*(I-1)+1
         L=K+2
         Vect1=CartGrad(K:L)
         CALL DGEMM_NNc(3,3,1,One,Zero,BoxShapeT,Vect1,Vect2)
         CartGrad(K:L)=Vect2
       ENDDO
       ! hardwire constrained lattice gradients to zero
       CartGrad(NCart-7)=Zero
       CartGrad(NCart-6)=Zero
       CartGrad(NCart-3)=Zero
     ENDIF
     !
     CALL New(CartA1,NCart)
     CALL New(IntA1,IntCsX%N)
     CALL New(IntA2,IntCsX%N)
     !
     CALL RefreshBMatInfo(IntCsX,XYZ,GOpt%TrfCtrl, &
                          GOpt%GConvCrit, &
                          GOpt%CoordCtrl%LinCrit, &
                          GOpt%CoordCtrl%TorsLinCrit, &
                          PBCDim,Print,SCRPath,Gi_O=.TRUE., &
                          Shift_O=1.0D-4,DoCleanCol_O=.FALSE.)
     CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
     !
     CALL CALC_BxVect(ISpB,JSpB,ASpB,IntA1%D,CartGrad)
     CALL GiInvIter(ISpB,JSpB,ASpB,CholData,IntA1%D,IntA2%D, &
                    NCart,IntCsX%N)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,IntA2%D,CartA1%D,Trp_O=.TRUE.)
     !
     Fact=DOT_PRODUCT(CartGrad,CartA1%D)/DOT_PRODUCT(CartGrad,CartGrad)
     Fact=Fact*100.D0
     WRITE(*,100) Fact
     WRITE(Out,100) Fact
     100 FORMAT('Percentage of Constraint Force That is Projected Out= ',F8.4)
     CartGrad=CartGrad-CartA1%D
     IF(PBCDim>0) THEN
       DO I=1,NatmsLoc-3
         K=3*(I-1)+1
         L=K+2
         Vect1=CartGrad(K:L)
         CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
         CartGrad(K:L)=Vect2
       ENDDO
     ENDIF
     !
     CALL Delete(IntCsX)
     CALL Delete(IntA2)
     CALL Delete(IntA1)
     CALL Delete(CartA1)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(CholData)
   END SUBROUTINE CleanConstrCart
!
!-------------------------------------------------------------------
!
   SUBROUTINE GetLattGrads(IntCL,CartGrad,XYZ,LattGrad,PBCDim)
     REAL(DOUBLE),DIMENSION(:)   :: CartGrad,LattGrad
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Vect(9)
     INTEGER                     :: NCart,NatmsLoc,NCoinc,NDim,PBCDim
     INTEGER                     :: I,J
     TYPE(INTC)                  :: IntCL,IntCs
     TYPE(DBL_RNK2)              :: AL,BL
     !
     IF(PBCDim==0) THEN
       LattGrad=Zero
       RETURN
     ENDIF
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     ! 
     CALL New(IntCs,IntCL%N)
     CALL SetEq(IntCL,IntCs,1,IntCL%N,1)
     IntCs%Constraint%L=.TRUE.
     CALL LatticeConstrAB(XYZ,IntCs,PBCDim,AL,BL,NCoinc)
     CALL Delete(BL)
     LattGrad=Zero
     Vect=CartGrad(NCart-8:NCart)
     CALL DGEMM_NNc(NCoinc,9,1,One,Zero,AL%D,Vect,LattGrad(1:NCoinc))
     CALL Delete(IntCs)
     CALL Delete(AL)
   END SUBROUTINE GetLattGrads
!
!----------------------------------------------------------------------
!
   FUNCTION HasLattice(Def)
     LOGICAL          :: HasLattice
     CHARACTER(LEN=*) :: Def
     !
     HasLattice= (Def(1:5)=='STRE_'.OR. &
                  Def(1:5)=='ALPHA'.OR. &
                  Def(1:4)=='BETA'.OR. &
                  Def(1:5)=='GAMMA'.OR. &
                  Def(1:6)=='AREA_L'.OR. &
                  Def(1:6)=='VOLM_L')
   END FUNCTION HasLattice
!
!----------------------------------------------------------------------
!
   SUBROUTINE CleanBLRatio(XYZ,Aux9,PBCDim,GConstr,ProjOut_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:,:) :: Aux9
     TYPE(BMATR)                 :: BL
     INTEGER                     :: PBCDim,I,J
     TYPE(Constr)                :: Gconstr
     TYPE(INTC)                  :: INTC_L
     REAL(DOUBLE)                :: Vect1(9),Vect2(9),Proj(9,9)
     REAL(DOUBLE)                :: V(9,2),U(9,2),S(2,2),SI(2,2)
     LOGICAL,OPTIONAL            :: ProjOut_O
     LOGICAL                     :: ProjOut
     !
     ProjOut=.FALSE.
     IF(PRESENT(ProjOut_O)) ProjOut=ProjOut_O
     IF(PBCDim==0) RETURN
     IF(ALL(GConstr%RatioABC<Zero).AND. &
        ALL(GConstr%RatioAlpBetGam<Zero)) RETURN
     !
     CALL LatticeINTC(IntC_L,PBCDim)
     CALL BMatrix(XYZ,IntC_L,BL,PBCDim,1.D0,1.D0)
     !
     Vect1=Zero
     Vect2=Zero
     DO I=1,3
       IF(GConstr%RatioABC(I)>Zero) &
         Vect1=Vect1+GConstr%RatioABC(I)*BL%BL%D(I,:)
       IF(GConstr%RatioAlpBetGam(I)>Zero) &
         Vect2=Vect2+GConstr%RatioAlpBetGam(I)*BL%BL%D(I+3,:)
     ENDDO
     !
     DO J=1,9
       V(J,1)=Vect1(J)
       V(J,2)=Vect2(J)
     ENDDO
     !
     CALL DGEMM_TNc(2,9,2,One,Zero,V,V,S)
     CALL SetDSYEVWork(2*2)
     CALL FunkOnSqMat(2,Inverse,S,SI)
     CALL UnSetDSYEVWork()
     !
     CALL DGEMM_NNc(9,2,2,One,Zero,V,SI,U)
     CALL DGEMM_NTc(9,2,9,One,Zero,U,V,Proj)
   ! Proj=-Proj
   ! DO I=1,9 ; Proj(I,I)=Proj(I,I)+One ; ENDDO
     !
     DO I=1,SIZE(Aux9,1)
       DO J=1,9 ; Vect1(J)=Aux9(I,J) ; ENDDO
       CALL DGEMM_NNc(9,9,1,One,One,Proj,Vect1,Vect2)
       DO J=1,9 ; Aux9(I,J)=Vect2(J) ; ENDDO
     ENDDO
     !
     CALL Delete(BL) 
     CALL Delete(IntC_L)
   END SUBROUTINE CleanBLRatio
!
!-------------------------------------------------------------------
!
   SUBROUTINE MergeFrag(JJ1,JJ2,R12,XYZ,FragID,Ftop,IFrag,JFrag)
     INTEGER                       ::  JJ1,JJ2,II,JJ,MM,NN
     REAL(DOUBLE),DIMENSION(:,:)   ::  XYZ
     REAL(DOUBLE),DIMENSION(3)     ::  XYZII,XYZJJ,Vect
     INTEGER,DIMENSION(:)          ::  FragID
     TYPE(INT_RNK2)                ::  FTop,FTop2
     TYPE(INT_VECT)                ::  IFrag,JFrag
     REAL(DOUBLE)                  ::  R12
     INTEGER                       ::  NFrag,I,J,NewDim,NDim2,MaxDim
     INTEGER                       ::  F1,F2
     !
     NFrag=SIZE(IFrag%I)-1
     II=FragID(JJ1)
     JJ=FragID(JJ2)
     XYZII=XYZ(1:3,II)
     XYZJJ=XYZ(1:3,JJ)
     CALL ClosestAtms(MM,NN,II,JJ,IFrag%I,JFrag%I, &
                          XYZII,XYZJJ,XYZ)
     JJ1=MM
     JJ2=NN
     F1=FragID(JJ1)
     F2=FragID(JJ2)
     Vect=XYZ(1:3,JJ1)-XYZ(1:3,JJ2)
     R12=SQRT(DOT_PRODUCT(Vect,Vect))
     !
     NDim2=SIZE(FTop%I,2)
     MaxDim=MAX(FTop%I(F1,1)+2,FTop%I(F2,1)+2)
     IF(MaxDim>NDim2) THEN
       NewDim=NDim2+10
       CALL New(FTop2,(/NFrag,NewDim/))
       FTop2%I=0
       DO I=1,NFrag
         DO J=1,NDim2
           FTop2%I(I,J)=FTop%I(I,J)
         ENDDO
       ENDDO 
       CALL Delete(FTop)
       CALL New(FTop,(/NFrag,NewDim/))
       DO I=1,NFrag
         DO J=1,NewDim
           FTop%I(I,J)=FTop2%I(I,J)
         ENDDO
       ENDDO 
       CALL Delete(FTop2)
     ENDIF
     !
     FTop%I(F1,1)=FTop%I(F1,1)+1
     FTop%I(F2,1)=FTop%I(F2,1)+1
     FTop%I(F1,FTop%I(F1,1)+1)=F2
     FTop%I(F2,FTop%I(F2,1)+1)=F1
   END SUBROUTINE MergeFrag
!
!-------------------------------------------------------------------
!
   FUNCTION ConnectedF(F1,F2,FTop)
     INTEGER                :: F1,F2,NDim2,J
     INTEGER,DIMENSION(:,:) :: FTop
     LOGICAL                :: ConnectedF
     !
     ConnectedF=.FALSE.
     DO J=2,FTop(F1,1)
       IF(FTop(F1,J)==F2) THEN
         ConnectedF=.TRUE.
         EXIT
       ENDIF
     ENDDO 
   END FUNCTION ConnectedF
!
!-------------------------------------------------------------------
!
   SUBROUTINE CleanLattGradRatio(PBCDim,XYZ,GConstr,CartGrad)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: CartGrad
     REAL(DOUBLE)                :: LattGrad(6),Vect(9)
     REAL(DOUBLE)                :: SumABC,SumAlpBetGam
     REAL(DOUBLE)                :: GSumABC,GSumAlpBetGam
     TYPE(Constr)                :: GConstr
     INTEGER                     :: PBCDim,NatmsLoc,NCoinc,NCart,I,J
     TYPE(INTC)                  :: IntC_L
     TYPE(DBL_RNK2)              :: AL,BL
     INTEGER                     :: IRefABC,IRefAlpBetGam
     !
     IF(.NOT.ANY(GConstr%RatioABC>Zero).AND. &
        .NOT.ANY(GConstr%RatioAlpBetGam>Zero)) RETURN
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     CALL LatticeINTC(IntC_L,PBCDim)
     IntC_L%Constraint%L=.TRUE.
     CALL LatticeConstrAB(XYZ,IntC_L,PBCDim,AL,BL,NCoinc)
     LattGrad=Zero
     DO J=1,9 ; Vect(J)=CartGrad(NCart-9+J) ; ENDDO
     CALL DGEMM_NNc(NCoinc,9,1,One,Zero,AL%D,Vect,LattGrad(1:NCoinc))       
     !
     ! Impose Ratios of lattice gradients
     !
     IRefABC=0
     SumABC=Zero
     GSumABC=Zero
     IRefAlpBetGam=0
     SumAlpBetGam=Zero
     GSumAlpBetGam=Zero
     DO I=1,3
       IF(GConstr%RatioABC(I)>Zero) THEN
         SumABC=SumABC+GConstr%RatioABC(I)
         GSumABC=GSumABC+LattGrad(I)
         IRefABC=I
       ENDIF
       IF(GConstr%RatioAlpBetGam(I)>Zero) THEN
         SumAlpBetGam=SumAlpBetGam+GConstr%RatioAlpBetGam(I)
         GSumAlpBetGam=GSumAlpBetGam+LattGrad(3+I)
         IRefAlpBetGam=I
       ENDIF
     ENDDO
     DO I=1,3
       IF(GConstr%RatioABC(I)>Zero) THEN
         LattGrad(I)=GSumABC/SumABC*GConstr%RatioABC(I)
       ENDIF
       IF(GConstr%RatioAlpBetGam(I)>Zero) THEN
         LattGrad(3+I)=GSumAlpBetGam/SumAlpBetGam*GConstr%RatioAlpBetGam(I)
       ENDIF
     ENDDO
     !
     CALL DGEMM_NNc(1,NCoinc,9,One,Zero,LattGrad(1:NCoinc),BL%D,Vect)       
     DO J=1,9 ; CartGrad(NCart-9+J)=Vect(J) ; ENDDO
     !
     CALL Delete(IntC_L)
     CALL Delete(BL)
     CALL Delete(AL)
   END SUBROUTINE CleanLattGradRatio
!
!-------------------------------------------------------------------
!
   END MODULE InCoords
