MODULE InCoords
!
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
   USE ControlStructures
   !
   IMPLICIT NONE
!
CONTAINS
!
!----------------------------------------------------------------
!
   SUBROUTINE BMatrix(XYZ,NIntC,IntCs,B,LinCrit,TorsLinCrit)
     !
     ! This subroutine calculates the sparse, TYPE(BMATR) type 
     ! representation of the B matrix 
     ! In current version Cartesian coordinates 
     ! must be passed in in AU
     ! Linear bendings _must_ always appear in pairs, 
     ! Defd as LINB1 and LINB2
     !
     IMPLICIT NONE
     INTEGER :: NIntC,NIntC2,I,J,K,L,NatmsLoc,IntCoo
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value
     TYPE(INTC) :: IntCs
     TYPE(BMATR):: B
     REAL(DOUBLE)  :: LinCrit,TorsLinCrit
     !
     NatmsLoc=SIZE(XYZ,2)
     CALL INTCValue(IntCs,XYZ,LinCrit,TorsLinCrit)
     !
     ! allocate B matrix
     !
     CALL NEW(B,NIntC)
     B%IB=0
     B%B=Zero
     !
     DO IntCoo=1,NIntC     
       !
       IF(.NOT.IntCs%Active%L(IntCoo)) CYCLE
       !
       IF(IntCs%Def%C(IntCoo)(1:4)=='STRE') THEN
         ! stre
         B%IB(IntCoo,1:2)=IntCs%Atoms%I(IntCoo,1:2)
         I=IntCs%Atoms%I(IntCoo,1)
         J=IntCs%Atoms%I(IntCoo,2)
         CALL STRE(XYZ(1:3,I),XYZ(1:3,J),BB_O=B%B(IntCoo,1:12))
        !write(*,100) IntCs%Def%C(IntCoo)(1:5), &
        !B%IB(IntCoo,1:4),B%B(IntCoo,1:12)
         !
       ELSE IF(IntCs%Def%C(IntCoo)(1:4)=='BEND') THEN
         ! bend
         B%IB(IntCoo,1:3)=IntCs%Atoms%I(IntCoo,1:3)
         I=IntCs%Atoms%I(IntCoo,1)
         J=IntCs%Atoms%I(IntCoo,2) !!! central atom
         K=IntCs%Atoms%I(IntCoo,3) 
         CALL BEND(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
                   BB_O=B%B(IntCoo,1:12))
        !write(*,100) IntCs%Def%C(IntCoo)(1:5), &
        !B%IB(IntCoo,1:4),B%B(IntCoo,1:12)
       ELSE IF(IntCs%Def%C(IntCoo)(1:4)=='OUTP') THEN
         ! out of plane
         B%IB(IntCoo,1:4)=IntCs%Atoms%I(IntCoo,1:4)
         I=IntCs%Atoms%I(IntCoo,1) !!! end atom
         J=IntCs%Atoms%I(IntCoo,2) !!! central atom
         K=IntCs%Atoms%I(IntCoo,3) !!! Def plane
         L=IntCs%Atoms%I(IntCoo,4) !!! Def plane
         CALL OutP(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ(1:3,L),&
           TorsLinCrit,IntCs%Active%L(IntCoo),BB_O=B%B(IntCoo,1:12))
       ! write(*,100) IntCs%Def%C(IntCoo)(1:5),B%IB(IntCoo,1:4), &
       !  B%B(IntCoo,1:12)*AngstromsToAu
       ELSE IF(IntCs%Def%C(IntCoo)(1:4)=='TORS') THEN
         ! torsion of i-j-k-l
         B%IB(IntCoo,1:4)=IntCs%Atoms%I(IntCoo,1:4)
         I=IntCs%Atoms%I(IntCoo,1)
         J=IntCs%Atoms%I(IntCoo,2) 
         K=IntCs%Atoms%I(IntCoo,3) 
         L=IntCs%Atoms%I(IntCoo,4) 
         CALL TORS(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K), &
           XYZ(1:3,L),TorsLinCrit,IntCs%Active%L(IntCoo), &
           BB_O=B%B(IntCoo,1:12))
        !write(*,*) IntCoo
        !write(*,100) IntCs%Def%C(IntCoo)(1:5),B%IB(IntCoo,1:4), &
        !B%B(IntCoo,1:12)
       ELSE IF(IntCs%Def%C(IntCoo)(1:5)=='LINB1') THEN
         ! linear bendig of i-j-k
         IF(IntCs%Def%C(IntCoo+1)(1:5)/='LINB2') &
            CALL Halt('LINB2 Definitions are not paired!')
         B%IB(IntCoo,1:3)=IntCs%Atoms%I(IntCoo,1:3)
         B%IB(IntCoo+1,1:3)=IntCs%Atoms%I(IntCoo+1,1:3)
         I=IntCs%Atoms%I(IntCoo,1)
         J=IntCs%Atoms%I(IntCoo,2) 
         K=IntCs%Atoms%I(IntCoo,3) 
         L=IntCs%Atoms%I(IntCoo,4)  ! reference atom
         CALL LinB(XYZ(1:3,I),XYZ(1:3,J),XYZ(1:3,K),XYZ,L, &
                   BB1=B%B(IntCoo,1:12),BB2=B%B(IntCoo+1,1:12))  
        !write(*,100) IntCs%Def%C(IntCoo)(1:5), &
        ! B%IB(IntCoo,1:4),B%B(IntCoo,1:12)
        !write(*,100) IntCs%Def%C(IntCoo)(1:5), &
        ! B%IB(IntCoo+1,1:4),B%B(IntCoo+1,1:12)
       ELSE IF(IntCs%Def%C(IntCoo)(1:5)=='LINB2') THEN
         CYCLE
       ELSE IF(IntCs%Def%C(IntCoo)(1:5)=='CARTX' ) THEN
         I=IntCs%Atoms%I(IntCoo,1)
         B%IB(IntCoo,1)=I
         CALL BCART('X',B%B(IntCoo,1:12),IntCs%Constraint%L(I))
       ELSE IF(IntCs%Def%C(IntCoo)(1:5)=='CARTY' ) THEN
         I=IntCs%Atoms%I(IntCoo,1)
         B%IB(IntCoo,1)=I
         CALL BCART('Y',B%B(IntCoo,1:12),IntCs%Constraint%L(I))
       ELSE IF(IntCs%Def%C(IntCoo)(1:5)=='CARTZ' ) THEN
         I=IntCs%Atoms%I(IntCoo,1)
         B%IB(IntCoo,1)=I
         CALL BCART('Z',B%B(IntCoo,1:12),IntCs%Constraint%L(I))
       ENDIF
       !
     ENDDO !!!! loop over internal coords
     100 FORMAT(A5,4I6,/,6F12.6,/,6F12.6)
     !
   END SUBROUTINE BMatrix
!
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
     REAL(DOUBLE)                :: Sum,Sum2,RU,RV,RW
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
       Sum=DOT_PRODUCT(U,V)
       IF(ABS(Sum)>One-1.D-8) THEN
         Value_O=PI
       ELSE
         Value_O=ACOS(Sum)
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
     REAL(DOUBLE)                 :: RU,RV,RW,Sum,Conv
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
     SinPhiU=SQRT(One-CosPhiU*CosPhiU)
     SinPhiV=SQRT(One-CosPhiV*CosPhiV)
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
       Sum=DOT_PRODUCT(UxW,VxW)/(SinPhiU*SinPhiV)
       IF(ABS(Sum)>One) THEN
         Sum=SIGN(One,Sum)
       ENDIF
         Value_O=ACOS(Sum)
       Sum=DOT_PRODUCT(W,UxV) ! orientation of torsion
       IF(ABS(Sum)>Zero) THEN
         Value_O=SIGN(Value_O,Sum)
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
   SUBROUTINE LinB(XI,XJ,XK,XYZ,IRef,BB1,BB2,Value1,Value2)
     !
     REAL(DOUBLE),DIMENSION(1:3) :: XI,XJ,XK
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(1:12),OPTIONAL:: BB1,BB2
     INTEGER                     :: I,J,K,L,M,N,KK,IRef
     REAL(DOUBLE)                :: Sum,CosAlpha,SinAlpha,RJI,RJK,RKI
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
     !
     CALL CROSS_PRODUCT(VectJI,VectJK,RIJJK)
     CosAlpha=DOT_PRODUCT(VectJI,VectJK)
     SinAlpha=SQRT(One-CosAlpha*CosAlpha)
     !
     ! Set up reference coordinate system
     !
   ! IF(SinAlpha>0.001D0) THEN
   !   VectRef=RIJJK
   ! ELSE
       IF(IRef==0) THEN
         VectAux=(/One,Zero,Zero/)
         CALL CROSS_PRODUCT(VectKI,VectAux,VectRef)
         Sum=SQRT(DOT_PRODUCT(VectRef,VectRef))
         IF(Sum<0.01D0) THEN
           VectAux=(/Zero,One,Zero/)
           CALL CROSS_PRODUCT(VectKI,VectAux,VectRef)
         ENDIF           
       ELSE
         VectRef=XYZ(1:3,IRef)-XK
       ENDIF
   ! ENDIF
     EZ=VectKI
     CALL CROSS_PRODUCT(VectRef,EZ,EY)
     Sum=SQRT(DOT_PRODUCT(EY,EY))
     EY=EY/Sum
     CALL CROSS_PRODUCT(EY,EZ,EX)
     Sum=SQRT(DOT_PRODUCT(EX,EX))
     EX=EX/Sum
     !
     IF(PRESENT(Value1)) THEN
       CALL BEND(XI,XJ,XK,Value_O=Value1,W_O=-EX)
     ENDIF
     IF(PRESENT(Value2)) THEN
       CALL BEND(XI,XJ,XK,Value_O=Value2,W_O=EY)
     ENDIF
     !
     IF(PRESENT(BB1)) THEN
       CALL BEND(XI,XJ,XK,BB_O=BB1,W_O=-EX)
     ENDIF
     !
     IF(PRESENT(BB2)) THEN
       CALL BEND(XI,XJ,XK,BB_O=BB2,W_O=EY)
     ENDIF
     !
   END SUBROUTINE LinB
!
!--------------------------------------------------------------------
!
   SUBROUTINE BCART(CHAR,B,Constraint)
     INTEGER :: I,J
     REAL(DOUBLE),DIMENSION(1:12) :: B
     CHARACTER :: CHAR
     LOGICAL   :: Constraint
     !
     B=Zero
     IF(CHAR=='X') THEN
       IF(Constraint) THEN
         B(1)=Zero 
       ELSE
         B(1)=One
       ENDIF
     ENDIF
     !
     IF(CHAR=='Y') THEN
       IF(Constraint) THEN
         B(2)=Zero 
       ELSE
         B(2)=One  
       ENDIF
     ENDIF
     !
     IF(CHAR=='Z') THEN
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
   SUBROUTINE DefineIntCoos(XYZ,AtNum,IntCs,NIntC,GCoordCtrl,Bond, &
                          AtmB,TOPS,ArchMem_O,HFileIn_O,iCLONE_O,iGEO_O)
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
     TYPE(INTC)                     :: IntCs
     TYPE(ATOMBONDS)                :: AtmBCov,AtmBVDW,AtmBTot,AtmB
     TYPE(BONDDATA)                 :: BondCov,BondVDW,BondTot,Bond
     CHARACTER(LEN=*),OPTIONAL      :: HFileIn_O
     INTEGER,OPTIONAL               :: iCLONE_O,iGEO_O,ArchMem_O
     INTEGER                        :: MaxBonds
     INTEGER,DIMENSION(:)           :: AtNum
     TYPE(CoordCtrl)                :: GCoordCtrl
     TYPE(TOPOLOGY)                 :: TOPS
     TYPE(IntCBox)                  :: Box
     TYPE(ANGLEDATA)                :: Angle
     TYPE(OUTPDATA)                 :: OutP 
     !
     NIntC=0
     NTorsion=0
     NatmsLoc=SIZE(XYZ,2)
     !
     CALL IntCBoxes(XYZ,Box)
     !
     CALL BondingScheme(XYZ,AtNum,1,AtmBCov,BondCov,TOPS,Box,GCoordCtrl)
     CALL BondingScheme(XYZ,AtNum,2,AtmBVDW,BondVDW,TOPS,Box,GCoordCtrl)
     CALL MergeBonds(BondCov,BondVDW,BondTot)
     CALL SortBonds(NatmsLoc,AtmBTot,BondTot)
     CALL Set_BONDDATA_EQ_BONDDATA(Bond,BondTot)
     CALL Set_AtmB_EQ_AtmB(AtmB,AtmBTot)
     !
     IF(PRESENT(HFileIn_O).AND.PRESENT(iCLONE_O).AND.&
        PRESENT(iGEO_O).AND.PRESENT(ArchMem_O)) THEN
       CALL ArchiveTop(TOPS,ArchMem_O, &
                       BondTot,AtmBTot,HFileIn_O,iCLONE_O,iGEO_O)
     ENDIF
     !
     ! Now define bond angles and torsions
     !
     CALL AngleList(AtmBTot,BondTot,TOPS,XYZ,Angle,OutP)
     CALL TorsionList(NatmsLoc,TOPS%Tot12,BondTot%IJ,XYZ, &
                      AtNum,TorsionIJKL,NTorsion)
     !
     ! Fill Data into IntCs
     !
     NIntC=BondTot%N+Angle%N+NTorsion+OutP%N
     IF(NIntC/=0) THEN          
       CALL New(IntCs,NIntC)
       IntCs%Def%C='     '
       IntCs%Atoms%I(:,:)=0   
       IntCs%Value%D=Zero
       IntCs%Constraint%L=.FALSE.
       IntCs%ConstrValue%D=Zero   
       !
         ILast=0
       DO I=1,BondTot%N
         IntCs%Def%C(ILast+I)='STRE '
         IntCs%Atoms%I(ILast+I,1:2)=BondTot%IJ%I(1:2,I)
       ENDDO
         ILast=BondTot%N
       DO I=1,Angle%N
         IntCs%Def%C(ILast+I)='BEND '
         IntCs%Atoms%I(ILast+I,1:3)=Angle%IJK%I(1:3,I)
       ENDDO
         ILast=ILast+Angle%N
       DO I=1,NTorsion
         IntCs%Def%C(ILast+I)='TORS '
         IntCs%Atoms%I(ILast+I,1:4)=TorsionIJKL%I(1:4,I)
       ENDDO
       ILast=ILast+NTorsion
       DO I=1,OutP%N
         IntCs%Def%C(ILast+I)='OUTP '
         IntCs%Atoms%I(ILast+I,1:4)=OutP%IJKL%I(1:4,I)
       ENDDO
     ENDIF
     !
     ! tidy up
     !
     CALL Delete(TorsionIJKL)
     !
     CALL Delete(AtmBTot)
     CALL Delete(AtmBCov)
     CALL Delete(AtmBVDW)
     CALL Delete(BondTot)
     CALL Delete(BondCov)
     CALL Delete(BondVDW)
     CALL Delete(Box)
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
       IF(Bond%LonelyAtom%I(I)==1) BondScore(I)=BondScore(I)+1
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
 !   IF(.NOT.(FoundHBond.OR.FoundMetLig.OR.&
 !      LonelyAtom)) THEN
 !     DoExclude=.TRUE.
 !     RETURN
 !   ENDIF
   END SUBROUTINE BondExcl
!
!----------------------------------------------------------------
!
   SUBROUTINE MoreBondArray(Bond,DimNew,DimOld)
     TYPE(BondDATA)     :: Bond,Bond2
     INTEGER            :: DimNew,DimOld   
     !
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
   SUBROUTINE TorsionList(NatmsLoc,Top12,BondIJ,XYZ, &
                          AtNum,TorsionIJKL,NTorsion)
     IMPLICIT NONE
     INTEGER                     :: NatmsLoc,I,I1,I2
     INTEGER                     :: N1,N2,J,J1,J2,NBond
     INTEGER                     :: NTorsion
     INTEGER                     :: I1Num,I2Num,JJ,I3,N3,JJ1,JJ2,NTIni
     TYPE(INT_RNK2)              :: Top12   
     TYPE(INT_RNK2)              :: BondIJ  
     TYPE(INT_RNK2)              :: TorsionIJKL,TorsionIJKLAux
     LOGICAL                     :: SelectTors
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Crit,Alph1,Alph2,Sum,PIHalf
     REAL(DOUBLE)                :: Sum1,Sum2
     INTEGER,DIMENSION(1:4)      :: Atoms
     INTEGER,DIMENSION(:)        :: AtNum
     !    
     PIHalf=PI*Half
     SelectTors=.TRUE.
     NBond=SIZE(BondIJ%I,2)
     NTorsion=0
     DO I=1,NBond
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       N1=Top12%I(I1,1)
       N2=Top12%I(I2,1)
       NTorsion=NTorsion+N1*N2
     ENDDO
     !    
     CALL New(TorsionIJKLAux,(/4,NTorsion/))
     !    
     NTorsion=0
     DO I=1,NBond 
       I1=BondIJ%I(1,I)
       I2=BondIJ%I(2,I)
       N1=Top12%I(I1,1)
       N2=Top12%I(I2,1)
       NTIni=NTorsion
       Crit=Zero
       Atoms=0 
       DO J1=1,N1
         JJ1=Top12%I(I1,J1+1)
         IF(SelectTors) THEN
           CALL BEND(XYZ(1:3,JJ1),XYZ(1:3,I1),XYZ(1:3,I2), &
                     Value_O=Alph1)
           Sum1=DBLE((Top12%I(JJ1,1)+1)*AtNum(JJ1))/ABS(Alph1-PIHalf+1.D-10)
         ENDIF
         DO J2=1,N2
           JJ2=Top12%I(I2,J2+1)
           IF(JJ1==JJ2.OR.JJ1==I2.OR.JJ2==I1) CYCLE
           IF(SelectTors) THEN
             CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,JJ2), &
                       Value_O=Alph2)
             Sum2=DBLE((Top12%I(JJ2,1)+1)*AtNum(JJ2))/ABS(Alph2-PIHalf+1.D-10)
             Sum=Sum1+Sum2
             IF(Sum>Crit) THEN
               Crit=Sum
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
       DO J=1,4 ; TorsionIJKL%I(J,I)=TorsionIJKLAux%I(J,I) ; ENDDO
     ENDDO
     CALL Delete(TorsionIJKLAux)
   END SUBROUTINE TorsionList
!
!-----------------------------------------------------------------------
!
   SUBROUTINE GetIntCs(XYZ,AtNumIn,IntCs,NIntC,Refresh,SCRPath, &
                       CtrlCoord,CtrlConstr,ArchMem,IntC_Extra, &
                       TOPS,Bond,AtmB,HFileIn_O,iCLONE_O,iGEO_O)
     !
     ! This subroutine constructs the IntCs array, which holds
     ! definitions of internal coordinates to be used in the 
     ! forthcoming geometry manipulation procedure.
     ! Refresh=1 : Refresh all definitions
     !        =2 : Refresh only definitions based on VDW interaction
     !        =3 : Do not refresh definitions, use the one from HDF
     !        =4 : Refresh/generate only the covalent coordinates
     !        =5 : only the extra coordinates from input
     ! WARNING! In the present version 
     ! bending -> linear bending transitions are 
     ! always checked and refreshed
     ! Later, check also linear bending -> bending transitions !
     ! 
     IMPLICIT NONE
     TYPE(INTC)                  :: IntCs,IntC_Bas,IntC_VDW
     TYPE(INTC)                  :: IntC_Extra,IntC_New
     TYPE(INTC)                  :: IntC_Cart
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(CoordCtrl)             :: CtrlCoord
     TYPE(Constr)                :: CtrlConstr
     LOGICAL                     :: DoFixMM
     REAL(DOUBLE),DIMENSION(:)   :: AtNumIn
     CHARACTER(LEN=*),OPTIONAL   :: HFileIn_O
     INTEGER,OPTIONAL            :: iCLONE_O,iGEO_O
     INTEGER                     :: NIntC,NIntC_Bas,NIntC_VDW
     INTEGER                     :: NIntC_Extra,NNew,Nintc_New
     INTEGER                     :: NIntC_Cart,MaxBonds,ArchMem
     INTEGER                     :: I,J,K,Refresh,NatmsLoc,II,ILast,III
     INTEGER                     :: I1,I2,I3,I4,NMax12
     INTEGER                     :: NStreGeOp,NBendGeOp
     INTEGER                     :: NLinBGeOp,NOutPGeOp,NTorsGeOp
     TYPE(INT_VECT)              :: AtNum
     TYPE(INT_RNK2)              :: Top12
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(IntCBox)               :: Box
     !
     NatmsLoc=SIZE(XYZ,2)
     IF(AllocQ(IntCs%Alloc)) CALL Delete(IntCs)
     !
     NIntC_Cart=0
     !
     CALL New(AtNum,NatmsLoc)
     DO I=1,NatmsLoc
       AtNum%I(I)=INT(AtNumIn(I))
     ENDDO
     !
     IF(Refresh==1) THEN !!! Total refresh
       IF(PRESENT(HFileIn_O).AND.PRESENT(iCLONE_O)) THEN
         CALL DefineIntCoos(XYZ,AtNum%I,IntC_Bas,NIntC_Bas,CtrlCoord,&
                            Bond,AtmB,TOPS,ArchMem_O=ArchMem, &
                            HFileIn_O=HFileIn_O, &
                            iCLONE_O=iCLONE_O,iGEO_O=iGEO_O)
       ELSE
         CALL DefineIntCoos(XYZ,AtNum%I,IntC_Bas,NIntC_Bas,CtrlCoord,&
                            Bond,AtmB,TOPS)
       ENDIF
     ELSE IF(Refresh==5) THEN !!! use only extra coords from input
       NIntC_Bas=0 
       NIntC_VDW=0 
     ELSE
       CALL Halt('Unknown Refresh option in GetIntCs') 
     ENDIF
     !
     ! Merge INTC arrays
     !
     CtrlCoord%NCov=NIntC_Bas
     NIntC=NIntC_Bas+NIntC_VDW+IntC_Extra%N+NIntC_Cart
     !
     CALL New(IntCs,NIntC)
     !
     ILast=0
     IF(NIntC_Bas/=0) THEN
       CALL Set_INTC_EQ_INTC(IntC_Bas,IntCs,1,NIntC_Bas,ILast+1)
       CALL Delete(IntC_Bas)
     ENDIF
       ILast=NIntC_Bas
     IF(NIntC_VDW/=0) THEN 
       CALL Set_INTC_EQ_INTC(IntC_VDW,IntCs,1,NIntC_VDW,ILast+1)
       CALL Delete(IntC_VDW)
     ENDIF
     !
     IF(IntC_Extra%N/=0) THEN
       ILast=ILast+NIntC_VDW
       CALL Set_INTC_EQ_INTC(IntC_Extra,IntCs,1,IntC_Extra%N,ILast+1)
     ENDIF
     !
     IF(NIntC_Cart/=0) THEN
       ILast=MAX(ILast+IntC_Extra%N,ILast+NIntC_VDW)
       CALL Set_INTC_EQ_INTC(IntC_Cart,IntCs,1,NIntC_Cart,ILast+1)
       CALL Delete(IntC_Cart)
     ENDIF
     !
     IF(.NOT.(NIntC==0.OR.IntC_Extra%N==0.OR.Refresh==5)) THEN
       CALL CleanINTC(IntCs,NIntC,NIntC_Bas,NIntC_VDW,IntC_Extra%N)
     ENDIF             
     !
     ! Set active all internal coords defd so far.
     ! 'Linear torsions will be deactivated when 
     ! their value is calculated
     ! by some subroutine (eg. INTCValue).
     ! The set of active coordinates may vary 
     ! during the process of optimization.
     !
     IntCs%Active%L=.TRUE.
     !
     ! Check for bending - lin.bending transions
     ! also for long range torsions!
     !
     CALL ChkBendToLinB(IntCs,NIntC,XYZ,AtNum%I,CtrlCoord,TOPS)
     !
     ! Count number of different internal coord types
     !
     NStreGeOp=0;NBendGeOp=0;NLinBGeOp=0;NOutPGeOp=0;NTorsGeOp=0
     DO I=1,NIntC
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
     IntCs%N=NIntC
     CALL Delete(AtNum)
    !CALL PrtIntCoords(IntCs,IntCs%Value%D,'GetIntC Internals')
   END SUBROUTINE GetIntCs
!
!-------------------------------------------------------------------
!
   SUBROUTINE CleanINTC(IntCs,NIntC,NIntC_Cov,NIntC_VDW,NIntC_Extra)
     !
     ! Now filter out repeated definitions, 
     ! which may have occured in INTC_Extra
     ! keep those values, which were defd. in extras.
     ! Do not check for Cartesians.
     !
     TYPE(INTC) :: IntCs,IntC_New
     INTEGER    :: NIntC_Cov,NIntC_VDW,NIntC,I,J,NNew,II,III,ILast
     INTEGER    :: NIntC_Extra
     !
     ILast=NIntC_Cov+NIntC_VDW
     II=0
     DO III=1,NIntC
       DO I=ILast+1,ILast+NIntC_Extra
           IF(IntCs%Def%C(I)(1:5)=='CART ') CYCLE
         DO J=1,ILast
           IF(IntCs%Def%C(J)(1:5)=='BLANK') CYCLE
           IF(IntCs%Atoms%I(J,1)==IntCs%Atoms%I(I,1).AND.&
              IntCs%Atoms%I(J,2)==IntCs%Atoms%I(I,2).AND.&
              IntCs%Atoms%I(J,3)==IntCs%Atoms%I(I,3).AND.&
              IntCs%Atoms%I(J,4)==IntCs%Atoms%I(I,4)) THEN
              II=II+1
              IntCs%Def%C(J)(1:5)='BLANK'
              IntCs%Atoms%I(J,1:4)=0      
           ENDIF
         ENDDO
       ENDDO
     ENDDO
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
   SUBROUTINE INTCValue(IntCs,XYZ,LinCrit,TorsLinCrit)
     !
     ! Determine value of internal coordinates.
     ! Input coordintes are now in atomic units!
     ! 
     IMPLICIT NONE
     TYPE(INTC) :: IntCs
     INTEGER :: NIntCs,I,J,K,L,I1,I2,I3,I4,NatmsLoc
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,LinCrit,TorsLinCrit
     !
     NIntCs=SIZE(IntCs%Def%C)
     NatmsLoc=SIZE(XYZ,2)
     !
     IntCs%Value%D=Zero 
     !
     DO I=1,NIntCs
       IF(.NOT.IntCs%Active%L(I)) THEN
         IntCs%Value%D(I)=Zero 
         CYCLE
       ENDIF
       I1=IntCs%Atoms%I(I,1)
       I2=IntCs%Atoms%I(I,2)
       I3=IntCs%Atoms%I(I,3)
       I4=IntCs%Atoms%I(I,4)
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         CALL STRE(XYZ(1:3,I1),XYZ(1:3,I2),Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         CALL BEND(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3), &
                   Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:5)=='LINB1') THEN
         CALL LinB(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),XYZ,I4,&
           Value1=IntCs%Value%D(I),Value2=IntCs%Value%D(I+1))  
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         CALL TORS(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),&
           XYZ(1:3,I4),TorsLinCrit,IntCs%Active%L(I), &
           Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
         CALL OUTP(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I4),&
           XYZ(1:3,I3),TorsLinCrit,IntCs%Active%L(I),Value_O=IntCs%Value%D(I))
         !
       ELSE IF(IntCs%Def%C(I)(1:5)=='CARTX') THEN
         IntCs%Value%D(I)=XYZ(1,I1)
       ELSE IF(IntCs%Def%C(I)(1:5)=='CARTY') THEN
         IntCs%Value%D(I)=XYZ(2,I1)
       ELSE IF(IntCs%Def%C(I)(1:5)=='CARTZ') THEN
         IntCs%Value%D(I)=XYZ(3,I1)
       ENDIF
       !
     ENDDO 
   END SUBROUTINE INTCValue
!
!------------------------------------------------------------
!
   SUBROUTINE GetSpBMatr(XYZ,B,SpB) 
     !
     ! Generate vibrational B matrix in quasi-sparse 
     ! TYPE(BMATR) representation and transform into sparse one
     !
     TYPE(BCSR)                           :: SpB
     INTEGER                              :: NIntC,NCart,I,NatmsLoc,NZ
     TYPE(BMATR)                          :: B
     REAL(DOUBLE),DIMENSION(:,:)          :: XYZ
     !
     NIntC=SIZE(B%IB,1)
     NatmsLoc=SIZE(XYZ,2)
     !
     NCart=3*NatmsLoc
     !
     ! Now, turn B matrix into sparse blocked representation
     !
     CALL Set_BCSR_EQ_BMATR(SpB,B)
     !     CALL PPrint(SpB,'SpB',Unit_O=6)
     !
     ! Generate sparse transpose B
     !
     !     CALL PPrint(SpB,'SpB2',Unit_O=6)
   END SUBROUTINE GetSpBMatr
!
!-------------------------------------------------------
!
   SUBROUTINE LJCell(GMLoc,LJCutOff,AtmMark,LJEps,LJRad, &
     XYZLJCell,AtmMarkLJCell,LJEpsLJCell,LJRadLJCell,NAtomsLJCell)
     !
     ! WARNING! Pass in coordinates in Angstroms!
     !
     ! Calculate coordinates of multiples of the elementary cell
     ! including the central cell, and store them in XYZLJCell.
     ! Storage is saved by filtering out replica atoms, which are
     ! far from the LJ sphere of the central cell. It's worth
     ! for the if statements, since later we can save at least 
     ! the same amount of if-s.
     !
     ! WARNING! Pass in coordinates in wrapped representation,
     ! so that all coordinates be greater than zero for central cell 
     !
     ! WARNING! Current implementation is optimal for 
     ! cells with rectangular
     ! lattice vectors! This is the case for large MD simulations.
     ! For small elementary cells the LJ sum should be fast enough,
     ! even with this technique.
     !
     TYPE(CRDS) :: GMLoc
     REAL(DOUBLE) :: LJCutOff,MaxCart,DX,DY,DZ
     REAL(DOUBLE) :: XTrans,YTrans,ZTrans
     INTEGER :: I,J,K,L,NX,NZ,NY,NAtomsLJCell
     INTEGER :: II,IA,IB,IC
     TYPE(DBL_VECT) :: LJEps,LJRad
     TYPE(INT_VECT) :: AtmMark
     TYPE(DBL_RNK2) :: XYZLJCell
     TYPE(DBL_RNK2) :: XYZLJCell2
     TYPE(DBL_VECT) :: LJEpsLJCell,LJRadLJCell
     TYPE(INT_VECT) :: AtmMarkLJCell
     TYPE(DBL_VECT) :: LJEpsLJCell2,LJRadLJCell2
     TYPE(INT_VECT) :: AtmMarkLJCell2
     !
     ! First, calculate maximum Cartesian 
     ! extension of the elementary cell
     !
     DX=GMLoc%PBC%BoxShape%D(1,1)
     DY=GMLoc%PBC%BoxShape%D(1,2)
     DZ=GMLoc%PBC%BoxShape%D(1,3)
     DO I=2,3
       DX=MAX(DX,GMLoc%PBC%BoxShape%D(I,1))
       DY=MAX(DY,GMLoc%PBC%BoxShape%D(I,2))
       DZ=MAX(DZ,GMLoc%PBC%BoxShape%D(I,3))
     ENDDO
     MaxCart=DX
     MaxCart=MAX(MaxCart,DY) 
     MaxCart=MAX(MaxCart,DZ) 
     !
     IF(LJCutOff<DX) THEN
       NX=1
     ELSE
       NX=INT(LJCutOff/DX)+1
     ENDIF
     !
     IF(LJCutOff<DY) THEN
       NY=1
     ELSE
       NY=INT(LJCutOff/DY)+1
     ENDIF
     !
     IF(LJCutOff<DZ) THEN
       NZ=1
     ELSE
       NZ=INT(LJCutOff/DZ)+1
     ENDIF
     !
     ! Calc. Number of atoms in the additional system
     ! and allocate XYZLJCell
     !
     NAtomsLJCell=(2*NX+1)*(2*NY+1)*(2*NZ+1)*GMLoc%Natms
     CALL New(XYZLJCell2,(/3,NAtomsLJCell/))
     CALL New(AtmMarkLJCell2,NAtomsLJCell)
     CALL New(LJEpsLJCell2,NAtomsLJCell)
     CALL New(LJRadLJCell2,NAtomsLJCell)
     !
     ! First, copy central cell coordinates into new array
     !
     DO I=1,GMLoc%Natms
       XYZLJCell2%D(1:3,I)=GMLoc%Carts%D(1:3,I)
       AtmMarkLJCell2%I(I)=AtmMark%I(I)
       LJEpsLJCell2%D(I)=LJEps%D(I)
       LJRadLJCell2%D(I)=LJRad%D(I)
     ENDDO
     !
     ! Now, copy atoms from LJ region of central cell
     !
     II=GMLoc%Natms
     DO IA=-NX,NX
     DO IB=-NY,NY
     DO IC=-NZ,NZ
       IF(IA==0.AND.IB==0.AND.IC==0) CYCLE
       DO I=1,GMLoc%Natms
         XTrans=GMLoc%Carts%D(1,I)+IA*GMLoc%PBC%BoxShape%D(1,1)+&
               IB*GMLoc%PBC%BoxShape%D(2,1)+IC*GMLoc%PBC%BoxShape%D(3,1)
                IF(XTrans<-LJCutOff .OR. XTrans>DX+LJCutOff) CYCLE
         YTrans=GMLoc%Carts%D(2,I)+IA*GMLoc%PBC%BoxShape%D(1,2)+&
               IB*GMLoc%PBC%BoxShape%D(2,2)+IC*GMLoc%PBC%BoxShape%D(3,2)
                IF(YTrans<-LJCutOff .OR. YTrans>DY+LJCutOff) CYCLE
         ZTrans=GMLoc%Carts%D(3,I)+IA*GMLoc%PBC%BoxShape%D(1,3)+&
               IB*GMLoc%PBC%BoxShape%D(2,3)+IC*GMLoc%PBC%BoxShape%D(3,3)
                IF(ZTrans<-LJCutOff .OR. ZTrans>DZ+LJCutOff) CYCLE
         II=II+1
         XYZLJCell2%D(1,II)=XTrans
         XYZLJCell2%D(2,II)=YTrans
         XYZLJCell2%D(3,II)=ZTrans
         AtmMarkLJCell2%I(II)=AtmMark%I(I)
         LJEpsLJCell2%D(II)=LJEps%D(I)
         LJRadLJCell2%D(II)=LJRad%D(I)
       ENDDO
     ENDDO
     ENDDO
     ENDDO
     !
     ! Compress
     !
     NAtomsLJCell=II 
     CALL New(XYZLJCell,(/3,NAtomsLJCell/))
     CALL New(AtmMarkLJCell,NAtomsLJCell)
     CALL New(LJEpsLJCell,NAtomsLJCell)
     CALL New(LJRadLJCell,NAtomsLJCell)
     XYZLJCell%D(1:3,1:NAtomsLJCell)=XYZLJCell2%D(1:3,1:NAtomsLJCell)
     AtmMarkLJCell%I(1:NAtomsLJCell)=AtmMarkLJCell2%I(1:NAtomsLJCell)
     LJEpsLJCell%D(1:NAtomsLJCell)=LJEpsLJCell2%D(1:NAtomsLJCell)
     LJRadLJCell%D(1:NAtomsLJCell)=LJRadLJCell2%D(1:NAtomsLJCell)
     !
     ! Tidy up
     !
     CALL Delete(XYZLJCell2)
     CALL Delete(AtmMarkLJCell2)
     CALL Delete(LJEpsLJCell2)
     CALL Delete(LJRadLJCell2)
   END SUBROUTINE LJCell
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetOneLJCell(GMLoc,AtmMark,LJEps,LJRad, &
     XYZLJCell,AtmMarkLJCell,LJEpsLJCell,LJRadLJCell,NAtomsLJCell)
     !
     TYPE(CRDS)     :: GMLoc
     INTEGER        :: NAtomsLJCell
     TYPE(DBL_VECT) :: LJEps,LJRad,LJEpsLJCell,LJRadLJCell
     TYPE(INT_VECT) :: AtmMark,AtmMarkLJCell
     TYPE(DBL_RNK2) :: XYZLJCell
     REAL(DOUBLE)   :: LJCutOff
     !
     NAtomsLJCell=GMLoc%Natms
     CALL New(XYZLJCell,(/3,NAtomsLJCell/))
     CALL New(AtmMarkLJCell,NAtomsLJCell)
     CALL New(LJEpsLJCell,NAtomsLJCell)
     CALL New(LJRadLJCell,NAtomsLJCell)
     XYZLJCell%D=GMLoc%Carts%D
     AtmMarkLJCell%I=AtmMark%I
     LJEpsLJCell%D=LJEps%D
     LJRadLJCell%D=LJRad%D
   END SUBROUTINE SetOneLJCell
!
!---------------------------------------------------------------------
!
   SUBROUTINE CartToInternal(XYZ,IntCs,VectCart,VectInt, &
                             TrfGrd,CtrlCoord,CtrlTrf,Print,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ   
     REAL(DOUBLE),DIMENSION(:)    :: VectCart,VectInt
     TYPE(DBL_VECT)  :: VectCartAux,VectIntAux
     TYPE(DBL_VECT)  :: VectCartAux2,VectIntAux2
     REAL(DOUBLE)    :: DiffMax,RMSD
     REAL(DOUBLE)    :: Sum
     INTEGER         :: NCart,I,II,J,NIntC
     INTEGER         :: NatmsLoc,Print
     TYPE(INTC)      :: IntCs
     TYPE(Cholesky)  :: CholData
     TYPE(INT_VECT)  :: ISpB,JSpB,IPerm1,IPerm2
     TYPE(DBL_VECT)  :: ASpB
     TYPE(GrdTrf)    :: TrfGrd
     TYPE(CoordCtrl) :: CtrlCoord
     TYPE(TrfCtrl)   :: CtrlTrf
     CHARACTER(LEN=*):: SCRPath
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc        
     NIntC=SIZE(IntCs%Def%C)
     !
     ! Get B matrix and Bt*B inverse
     !
     IF(CtrlTrf%DoClssTrf) THEN
       CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
     ELSE
       CALL GetBMatInfo(TRIM(SCRPath)//'Mix',ISpB,JSpB,ASpB,CholData, &
                        IPerm1_O=IPerm1,IPerm2_O=IPerm2)
       CALL ModifyGrad(VectCart,IPerm1%I,IPerm2%I,CtrlTrf)
       CALL Delete(IPerm1)
       CALL Delete(IPerm2)
     ENDIF
     !
     CALL New(VectCartAux,NCart)
     CALL New(VectCartAux2,NCart)
     CALL New(VectIntAux,NIntC)
     CALL New(VectIntAux2,NIntC)
     !
     VectInt=Zero
     !
     IF(Print>=DEBUG_GEOP_MIN) THEN
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
     DO II=1,TrfGrd%MaxIt_GrdTrf
       !
       VectCartAux%D=Zero
       !
       ! gc-Bt*gi
       !
       CALL CALC_BxVect(ISpB,JSpB,ASpB,VectInt,VectCartAux%D,Trp_O=.TRUE.)
       VectCartAux%D=VectCart-VectCartAux%D
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
         IF(Print>=DEBUG_GEOP_MIN) THEN
           WRITE(*,*) 'Rescale Step from ',DiffMax,' to ',TrfGrd%MaxGradDiff
           WRITE(Out,*) 'Rescale Step from ',DiffMax,' to ',TrfGrd%MaxGradDiff
         ENDIF
         SUM=TrfGrd%MaxGradDiff/DiffMax
         VectIntAux%D(:)=SUM*VectIntAux%D(:)
         DiffMax=TrfGrd%MaxGradDiff
       ENDIF
       !
       ! gi+B*GcInv*[gc-Bt*gi]
       !
       VectInt=VectInt+VectIntAux%D
       !
       ! Review iteration
       !
       IF(Print>=DEBUG_GEOP_MIN) THEN
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
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,777) 
         WRITE(*,778) 
         WRITE(Out,777) 
         WRITE(Out,778) 
         777 FORMAT('Stop Gradient Trf, max. number '//&
                      'of Iterations exceeded!')
         778 FORMAT('Use current gradient vector!')
       ENDIF
     ELSE
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,120) II
         WRITE(Out,120) II
       ENDIF
     ENDIF
     120  FORMAT('Gradient transformation converged in ',I3,' steps')
     !
     ! Tidy up
     !
     CALL Delete(VectIntAux2)
     CALL Delete(VectIntAux)
     CALL Delete(VectCartAux2)
     CALL Delete(VectCartAux)
     CALL DeleteBMatInfo(ISpB,JSpB,ASpB,CholData)
     !
   END SUBROUTINE CartToInternal
!
!------------------------------------------------------------------
!
   SUBROUTINE InternalToCart(XYZ,IntCs,VectInt,Print, &
       GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath,AtNum_O,DoDeloc_O)
     REAL(DOUBLE),DIMENSION(:,:)        :: XYZ
     REAL(DOUBLE),DIMENSION(:)          :: VectInt
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: AtNum_O
     LOGICAL,OPTIONAL                   :: DoDeloc_O
     LOGICAL                            :: DoDeloc
     TYPE(DBL_VECT)            :: VectCart
     TYPE(DBL_VECT)            :: VectCartAux,VectIntAux
     TYPE(DBL_VECT)            :: VectCartAux2
     TYPE(DBL_VECT)            :: VectIntReq,IntCValSt
     TYPE(DBL_RNK2)            :: ActCarts
     REAL(DOUBLE)              :: DiffMax,RMSD,RMSDOld
     REAL(DOUBLE)              :: Sum,ConstrMax,ConstrRMS
     REAL(DOUBLE)              :: ConstrRMSOld,ConstrMaxCrit,RMSCrit
     INTEGER                   :: NCart,I,IStep,J
     INTEGER                   :: NIntC,NConstr,IRep,RepMax
     INTEGER                   :: NatmsLoc,NCartConstr
     TYPE(INTC)                :: IntCs
     TYPE(INT_VECT)            :: ISpB,JSpB
     TYPE(DBL_VECT)            :: ASpB
     LOGICAL                   :: RefreshB,RefreshAct
     LOGICAL                   :: DoIterate
     TYPE(Cholesky)            :: CholData
     TYPE(BackTrf)             :: GBackTrf
     TYPE(Constr)              :: GConstr
     TYPE(TrfCtrl)             :: GTrfCtrl
     TYPE(CoordCtrl)           :: GCoordCtrl
     CHARACTER(LEN=*)          :: SCRPath
     LOGICAL                   :: DoClssTrf,Print2,DoRepeat
     INTEGER                   :: Print
     TYPE(INT_VECT)            :: ISpBD,JSpBD
     TYPE(DBL_VECT)            :: ASpBD
     TYPE(DBL_RNK2)            :: UMatr
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc   
     NIntC=SIZE(IntCs%Def%C)
     Print2=(Print>=DEBUG_GEOP_MAX)
     DoClssTrf=.TRUE.
     DoRepeat=.FALSE.
     RepMax=5
     IF(PRESENT(DoDeloc_O)) THEN
       DoDeloc=DoDeloc_O
     ELSE
       DoDeloc=.FALSE.
     ENDIF
     !
     ! Auxiliary arrays
     !
     CALL New(ActCarts,(/3,NatmsLoc/))
     CALL New(VectCart,NCart)
     CALL New(VectCartAux,NCart)
     CALL New(VectCartAux2,NCart)
     CALL New(VectIntAux,NIntC)
     CALL New(VectIntReq,NIntC)
     CALL New(IntCValSt,NIntC)
     VectIntReq%D=VectInt
     CALL INTCValue(IntCs,XYZ, &
                    GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
     IntCValSt%D=IntCs%Value%D
     !
     ! Repeat until convergence
     !
     DO IRep=1,RepMax
       RefreshB=.TRUE.
       RefreshAct=.TRUE.
       !
       IF(DoRepeat) THEN
         VectIntAux%D=VectIntReq%D-IntCValSt%D  
         CALL CutOffDispl(VectIntAux%D,IntCs, &
                     GCoordCtrl%MaxStre,GCoordCtrl%MaxAngle,Fact_O=Half)
         VectIntReq%D=IntCValSt%D+VectIntAux%D
       ENDIF
       !
       ! The required new value of internal coordinates
       !
       IF(DoDeloc) THEN
         CALL ReadBMATR(ISpBD,JSpBD,ASpBD,TRIM(SCRPath)//'UMatr', &
                        UMatr_O=UMatr)
         CALL DelocP(VectIntReq%D,ISpBD,JSpBD,ASpBD,UMatr)
       ELSE
         CALL MapBackAngle(IntCs,VectIntReq%D) 
       ENDIF
       !
       !initialization of new Cartesians
       !
      !IF(DoClssTrf) THEN
         ActCarts%D=XYZ
         CALL CartRNK2ToCartRNK1(VectCart%D,ActCarts%D)
      !ELSE
      !  CALL CartRNK2ToCartRNK1(VectCart%D,XYZ)
      !  CALL TranslToAt1(VectCart%D,GTrfCtrl%ThreeAt)
      !  CALL RotToAt(VectCart%D,GTrfCtrl%RotAt2ToX)
      !  CALL RotToAt(VectCart%D,GTrfCtrl%RotAt3ToXY)
      !  !transform Cartesian constraints!
      !  CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
      !  CALL ReSetConstr(IntCs,ActCarts%D)
      !ENDIF
       !
       ! Internal --> Cartesian transformation
       !
       IF(Print>=DEBUG_GEOP_MIN) THEN
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
        !IF(PRESENT(AtNum_O)) THEN
        !  CALL PrtXYZ(Atnum_O,ActCarts%D,'backtrf',&
        !              'Step='//TRIM(IntToChar(IStep)))
        !ENDIF
         !
         ! Get B and refresh values of internal coords
         !
         CALL INTCValue(IntCs,ActCarts%D, &
                        GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
         IF(RefreshB.AND.RefreshAct) THEN
           CALL RefreshBMatInfo(IntCs,ActCarts%D,GTrfCtrl, &
                                GCoordCtrl,Print,SCRPath,.TRUE.)
           CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
         ENDIF
         !
         ! Calculate difference between required and actual internals
         !
         VectIntAux%D=VectIntReq%D-IntCs%Value%D
         IF(DoDeloc) CALL DelocP(VectIntAux%D,ISpBD,JSpBD,ASpBD,UMatr)
         CALL MapAngleDispl(IntCs,VectIntAux%D) 
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
         CALL CALC_BxVect(ISpB,JSpB,ASpB,VectIntAux%D, &
                          VectCartAux%D,Trp_O=.TRUE.)
         !
         ! GcInv*Bt*[phi_r-phi_a]
         !
         CALL CALC_GcInvCartV(CholData,VectCartAux%D,VectCartAux2%D)
         !
         ! Project out rotations and translations 
         !
         IF(DoClssTrf) THEN 
           IF(GTrfCtrl%DoTranslOff) &
             CALL TranslsOff(VectCartAux2%D,Print2)
           IF(GTrfCtrl%DoRotOff) &
             CALL RotationsOff(VectCartAux2%D,ActCarts%D,Print2)
         ENDIF
         !
         ! Check convergence
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
         CALL SetCartConstr(VectCart%D,VectCartAux2%D, &
                            IntCs,GConstr%NCartConstr)
         VectCart%D=VectCart%D+VectCartAux2%D
         CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
         !
         ! Review iteration
         !
         IF(Print>=DEBUG_GEOP_MIN) THEN
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
           IF(Print>=DEBUG_GEOP_MIN) THEN
             WRITE(*,*) 'Rescaling and Repeating back-transformation'
             WRITE(Out,*) 'Rescaling and Repeating back-transformation'
           ENDIF
           DoRepeat=.TRUE.
           IF(IRep==RepMax) THEN
             CALL Halt('Iterative backtransformation did not converge')
           ENDIF
           CYCLE
         ENDIF
         IF(Print>=DEBUG_GEOP_MIN) THEN
           WRITE(*,180) 
           WRITE(Out,180) 
           WRITE(*,190) 
           WRITE(Out,190) 
         ENDIF
         EXIT 
       ELSE
         CALL INTCValue(IntCs,ActCarts%D, &
                        GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
         VectIntAux%D=IntCValSt%D-IntCs%Value%D
         CALL MapAngleDispl(IntCs,VectIntAux%D) 
         IF(IRep<RepMax) THEN
           IF(ABS(MAXVAL(VectIntAux%D))>1.01D0*GCoordCtrl%MaxAngle.OR. &
              ABS(MINVAL(VectIntAux%D))>1.01D0*GCoordCtrl%MaxStre) THEN
             DoRepeat=.TRUE.
             CYCLE
           ENDIF
         ENDIF
         WRITE(*,220) IStep
         WRITE(Out,220) IStep
         EXIT
       ENDIF
       !
     ENDDO !!! Repeat
     !
     180  FORMAT('Stop Coord Back-Trf, max. number of Iterations exceeded!')
     190  FORMAT('Use Current Geometry!')
     220  FORMAT('Coordinate back-transformation converged in ',&
                 I3,' steps')
     !
     ! Fill new Cartesians into XYZ  
     !
    !IF(.NOT.DoClssTrf) THEN
    !  CALL CartRNK2ToCartRNK1(VectCart%D,ActCarts%D)
    !  CALL RotToAt(VectCart%D,GTrfCtrl%RotAt3ToXY,Rev_O=.TRUE.)
    !  CALL RotToAt(VectCart%D,GTrfCtrl%RotAt2ToX,Rev_O=.TRUE.)
    !  CALL TranslToAt1(VectCart%D,GTrfCtrl%ThreeAt, &
    !                   Vect_O=-GTrfCtrl%TranslAt1)
    !  CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
    !  CALL ReSetConstr(IntCs,ActCarts%D)
    !ENDIF
     XYZ=ActCarts%D
     !
     ! Tidy up
     !
     CALL Delete(IntCValSt)
     CALL Delete(VectIntReq)
     CALL Delete(VectIntAux)
     CALL Delete(VectCartAux2)
     CALL Delete(VectCartAux)
     CALL Delete(VectCart)
     CALL Delete(ActCarts)
     IF(DoDeloc) THEN
       CALL Delete(ISpBD)
       CALL Delete(JSpBD)
       CALL Delete(ASpBD)
       CALL Delete(UMatr)
     ENDIF
   END SUBROUTINE InternalToCart
!
!----------------------------------------------------------
!
   SUBROUTINE PrtIntCoords(IntCs,Value,CHAR)
     !
     TYPE(INTC) :: IntCs
     INTEGER    :: I,NIntC,J
     REAL(DOUBLE),DIMENSION(:) :: Value
     REAL(DOUBLE) :: SUM,SumConstr,Conv
     CHARACTER(LEN=*) :: CHAR   
     !
     Conv=180.D0/PI
     NIntC=SIZE(IntCs%Def%C)
     !
     WRITE(*,*) TRIM(CHAR)
     WRITE(Out,*) TRIM(CHAR)
     WRITE(*,*) '         INTERNAL COORDINATES'
     WRITE(*,*) '       DEFINITION  ATOMS_INVOLVED         VALUE'//&
                '        CONSTRAINT    ACTIVE'
     WRITE(Out,*) 'INTERNAL COORDINATES'
     WRITE(Out,*) '       DEFINITION  ATOMS_INVOLVED      VALUE'//&
                  '        CONSTRAINT    ACTIVE'
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         SUM=Value(I)/AngstromsToAu
         SUMConstr=IntCs%ConstrValue%D(I)/AngstromsToAu
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND'.OR. &
               IntCs%Def%C(I)(1:4)=='LINB'.OR. &
               IntCs%Def%C(I)(1:4)=='TORS'.OR. &
               IntCs%Def%C(I)(1:4)=='OUTP') THEN
         SUM=Value(I)*Conv
         SUMConstr=IntCs%ConstrValue%D(I)*Conv
       ELSE IF(IntCs%Def%C(I)(1:4)=='CART') THEN
         SUM=Value(I)
         SUMConstr=IntCs%ConstrValue%D(I)
       ENDIF
       WRITE(*,111) I,IntCs%Def%C(I),IntCs%Atoms%I(I,1:4),SUM, &
         IntCs%Constraint%L(I),SumConstr,IntCs%Active%L(I)
       WRITE(Out,111) I,IntCs%Def%C(I),IntCs%Atoms%I(I,1:4),SUM, &
         IntCs%Constraint%L(I),SumConstr,IntCs%Active%L(I)
     ENDDO
     !      
     111 FORMAT(I7,2X,A5,2X,4I5,2X,F12.6,L5,F12.6,L5)
     222 FORMAT(I7,2X,A5,2X,4I5,2X,3F12.6,L5,F12.6,L5)
     !
   END SUBROUTINE PrtIntCoords
!
!----------------------------------------------------------
!
   SUBROUTINE RotationsOff(CartVect,XYZ,Print)
     REAL(DOUBLE),DIMENSION(:) :: CartVect
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: X,Y,Z,XX,YY,ZZ,XY,YZ,ZX
     REAL(DOUBLE)                :: CMX,CMY,CMZ
     REAL(DOUBLE)                :: SUM,SUM1,SUM2,SUM3
     REAL(DOUBLE)                :: V1X,V1Y,V1Z
     REAL(DOUBLE)                :: V2X,V2Y,V2Z
     REAL(DOUBLE)                :: V3X,V3Y,V3Z
     TYPE(DBL_RNK2)              :: Theta,Theta2,CMCarts
     TYPE(DBL_Vect)              :: Rot1,Rot2,Rot3,Vect
     INTEGER                     :: NCart,NatmsLoc,I,J,INFO
     LOGICAL                     :: Print
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
     CALL CenterOfMass(CMX,CMY,CMZ,XYZ_O=XYZ)
     !
     ! Mass centered Cartesians
     !
     DO I=1,NatmsLoc
       CMCarts%D(1,I)=XYZ(1,I)-CMX
       CMCarts%D(2,I)=XYZ(2,I)-CMY
       CMCarts%D(3,I)=XYZ(3,I)-CMZ
     ENDDO
     !
     ! Build inertial momentum tensor
     !
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
     SUM=SQRT(DOT_PRODUCT(Rot1%D,Rot1%D))
     IF(Sum>1.D-6) THEN
       Sum=One/Sum
       Rot1%D=Sum*Rot1%D
     ENDIF
     SUM=SQRT(DOT_PRODUCT(Rot2%D,Rot2%D))
     IF(Sum>1.D-6) THEN
       Sum=One/Sum
       Rot2%D=Sum*Rot2%D
     ENDIF
     SUM=SQRT(DOT_PRODUCT(Rot3%D,Rot3%D))
     IF(Sum>1.D-6) THEN
       Sum=One/Sum
       Rot3%D=Sum*Rot3%D
     ENDIF
     !
     ! test orthogonalities
     !
     !      CALL TranslsOff(Rot1%D,.TRUE.)
     !      CALL TranslsOff(Rot2%D,.TRUE.)
     !      CALL TranslsOff(Rot3%D,.TRUE.)
     !      SUM=DOT_PRODUCT(Rot1%D,Rot1%D)
     !      write(*,*) 'rot test 1 1 ',sum
     !      SUM=DOT_PRODUCT(Rot1%D,Rot2%D)
     !      write(*,*) 'rot test 1 2 ',sum
     !      SUM=DOT_PRODUCT(Rot1%D,Rot3%D)
     !      write(*,*) 'rot test 1 3 ',sum
     !      SUM=DOT_PRODUCT(Rot2%D,Rot2%D)
     !      write(*,*) 'rot test 2 2 ',sum
     !      SUM=DOT_PRODUCT(Rot2%D,Rot3%D)
     !      write(*,*) 'rot test 2 3 ',sum
     !      SUM=DOT_PRODUCT(Rot3%D,Rot3%D)
     !      write(*,*) 'rot test 3 3 ',sum
     !
     ! Now, project out rotations from Cartesian displacement vector
     !
     SUM=DOT_PRODUCT(CartVect(1:NCart),CartVect(1:NCart))
     SUM1=DOT_PRODUCT(Rot1%D(1:NCart),CartVect(1:NCart))
     SUM2=DOT_PRODUCT(Rot2%D(1:NCart),CartVect(1:NCart))
     SUM3=DOT_PRODUCT(Rot3%D(1:NCart),CartVect(1:NCart))
     CartVect(1:NCart)=CartVect(1:NCart)-SUM1*Rot1%D(1:NCart)
     CartVect(1:NCart)=CartVect(1:NCart)-SUM2*Rot2%D(1:NCart)
     CartVect(1:NCart)=CartVect(1:NCart)-SUM3*Rot3%D(1:NCart)
     !
     ! Percentage of rotations
     !
     IF(SUM>1.D-6) THEN
       SUM1=SUM1*SUM1/SUM*100.D0
       SUM2=SUM2*SUM2/SUM*100.D0
       SUM3=SUM3*SUM3/SUM*100.D0
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
   END SUBROUTINE RotationsOff
!
!----------------------------------------------------------
!
   SUBROUTINE TranslsOff(CartVect,Print)
     REAL(DOUBLE),DIMENSION(:) :: CartVect
     REAL(DOUBLE)              :: SUM,SUM1,SUM2,SUM3
     TYPE(DBL_VECT)            :: Tr1,Tr2,Tr3
     INTEGER                   :: I,J,NCart,NatmsLoc
     LOGICAL                   :: Print
     NCart=SIZE(CartVect)
     NatmsLoc=NCart/3
     CALL New(Tr1,NCart) 
     CALL New(Tr2,NCart) 
     CALL New(Tr3,NCart) 
     Tr1%D=Zero
     Tr2%D=Zero
     Tr3%D=Zero
     Sum=One/DBLE(NatmsLoc)
     DO I=1,NatmsLoc
       J=(I-1)*3
       Tr1%D(J+1)=Sum
       Tr2%D(J+2)=Sum
       Tr3%D(J+3)=Sum
     ENDDO
     !
     ! Now, project out translations from CartVect
     !
     SUM =DOT_PRODUCT(CartVect,CartVect)
     IF(SUM > 1.D-6) THEN
       SUM1=DOT_PRODUCT(Tr1%D,CartVect)
         CartVect=CartVect-SUM1*Tr1%D
       SUM2=DOT_PRODUCT(Tr2%D,CartVect)
         CartVect=CartVect-SUM2*Tr2%D
       SUM3=DOT_PRODUCT(Tr3%D,CartVect)
         CartVect=CartVect-SUM3*Tr3%D
       SUM1=SUM1*SUM1/SUM*100.D0
       SUM2=SUM2*SUM2/SUM*100.D0
       SUM3=SUM3*SUM3/SUM*100.D0
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
     REAL(DOUBLE)                 :: Sum,SumM,Sum1,Sum2,Sum3,V1N,V2N
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
     IF(ABS(Sum1) < 1.D-6) THEN  !!! V1 & V2 are parallel
       Rot=Zero
       Sum=DOT_PRODUCT(V1,V2)
       DO I=1,3 ; Rot(I,I)=Sum ; ENDDO
     ELSE      
       Sum=One/SQRT(Sum1)
       CrossProd%D=SUM*CrossProd%D
       CosPhi=DOT_PRODUCT(V1,V2)     
       SinPhi=SQRT(SUM1)
       Sum=One-CosPhi
       DO III=1,2    
         Rot(1,1)=CosPhi + Sum*CrossProd%D(1)**2
         Rot(2,2)=CosPhi + Sum*CrossProd%D(2)**2
         Rot(3,3)=CosPhi + Sum*CrossProd%D(3)**2
         !
         Sum1=Sum*CrossProd%D(1)*CrossProd%D(2)
         Sum2=SinPhi*CrossProd%D(3)
         Rot(1,2)=Sum1-Sum2
         Rot(2,1)=Sum1+Sum2
         !
         Sum1=Sum*CrossProd%D(1)*CrossProd%D(3)
         Sum2=SinPhi*CrossProd%D(2)
         Rot(1,3)=Sum1+Sum2
         Rot(3,1)=Sum1-Sum2
         !
         Sum1=Sum*CrossProd%D(2)*CrossProd%D(3)
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
   SUBROUTINE ChkBendToLinB(IntCs,NIntC,XYZ,AtNum,CtrlCoord,TOPS)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER                     :: NIntC,Nintc_New
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(INT_VECT)              :: LinAtom,MarkLinb,MarkLongR
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER,DIMENSION(:)        :: AtNum
     REAL(DOUBLE)                :: Value,Value2,Conv
     INTEGER                     :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
     INTEGER                     :: I,J,K,L,NatmsLoc
     TYPE(INT_RNK2)              :: LinBBridge
     TYPE(CoordCtrl)             :: CtrlCoord
     !
     ! Now check for bending -> linear bending transitions
     !
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
           IF(TOPS%Tot12%I(I2,1)>2) THEN
             IntCs%Active%L(I)=.FALSE.
             CYCLE
           ENDIF
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
     CALL New(MarkLongR,NIntC)
     MarkLongR%I=0
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:5)=='LINB1') THEN
         I1=IntCs%Atoms%I(I,1)
         I2=IntCs%Atoms%I(I,2)
         I3=IntCs%Atoms%I(I,3)
         L =IntCs%Atoms%I(I,4)
         CALL LinB(XYZ(1:3,I1),XYZ(1:3,I2),XYZ(1:3,I3),XYZ,L,&
           Value1=IntCs%Value%D(I),Value2=IntCs%Value%D(I+1))  
         MarkLongR%I(I)=1 
         MarkLongR%I(I+1)=1 
       ENDIF
     ENDDO 
     !
     ! Now recognize colinear atoms of the molecule and 
     ! introduce long-range torsions. 
     !
     CALL LongRangeIntC(CtrlCoord,TOPS%Tot12,IntCs,NIntC, &
                        XYZ,MarkLongR,AtNum)
     !
     CALL Delete(MarkLongR)
   END SUBROUTINE ChkBendToLinB    
!
!-------------------------------------------------------
!
   SUBROUTINE FindLinBRef(XYZ,Atoms,Top12)
     INTEGER,DIMENSION(1:4)  :: Atoms
     INTEGER,DIMENSION(:,:)  :: Top12
     INTEGER                 :: I,J,K,L,I1,I2,I3,N1,N3,IC,N,II
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)            :: Val,ValM,Conv,Sum,RefLinCrit
     !
     Conv=180.D0/PI
     RefLinCrit=10.D0
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
         Sum=ABS(120.D0-Val*Conv)
         IF(Sum<ValM) THEN
           II=J
           ValM=Sum 
         ENDIF
       ENDDO
     ENDDO
     Atoms(4)=II 
   END SUBROUTINE FindLinBRef
!
!-------------------------------------------------------
!
   SUBROUTINE SetConstraint(IntCs,XYZ,Displ,LinCrit,TorslinCrit, &
                            NConstr,DoInternals)
     TYPE(INTC)     :: IntCs
     INTEGER        :: I,J,NIntC,LastIntcGeom,NDim,JJ,NConstr
     TYPE(DBL_VECT) :: Displ
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ  
     REAL(DOUBLE)   :: LinCrit,TorsLinCrit
     LOGICAL        :: DoInternals
     !
     IF(NConstr==0)RETURN
     IF(AllocQ(IntCs%Alloc)) THEN
       NIntC=SIZE(IntCs%Def%C)
     ELSE
       NIntC=0
     ENDIF
     NDim=SIZE(Displ%D)
     IF(NDim/=NIntC.AND.DoInternals) &
         Call Halt('Dimensionality error in SetConstraint')
     !
     ! Get values of constraints 
     !
     IF(DoInternals) THEN
       CALL INTCValue(IntCs,XYZ,LinCrit,TorsLinCrit)
       DO I=1,NIntC
         IF(IntCs%Constraint%L(I)) THEN
           Displ%D(I)=IntCs%ConstrValue%D(I)-IntCs%Value%D(I)
           !write(*,100) i,IntCs%ConstrValue%D(I),&
           !IntCs%Value%D(I),Displ%D(I)
           !100 format('int= ',i4,3F12.6)
         ENDIF
       ENDDO	
     ELSE
       DO I=1,NIntC
         IF(IntCs%Constraint%L(I)) THEN
           IF(IntCs%Def%C(I)(1:4)/='CART') CALL Halt('Internal '// &
              'coord constraints in Cartesian optimization')
             J=IntCs%Atoms%I(I,1)
           IF(IntCs%Def%C(I)(1:5)/='CARTX') JJ=(J-1)*3+1
           IF(IntCs%Def%C(I)(1:5)/='CARTY') JJ=(J-1)*3+2
           IF(IntCs%Def%C(I)(1:5)/='CARTZ') JJ=(J-1)*3+3
           Displ%D(JJ)=IntCs%ConstrValue%D(I)-IntCs%Value%D(I)
         ENDIF
       ENDDO	
     ENDIF
   END SUBROUTINE SetConstraint
!
!-------------------------------------------------------
!
   SUBROUTINE ConstrConv(IntCs,IntDiff,ConstrMax,ConstrRMS)
     TYPE(INTC) :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: IntDiff
     REAL(DOUBLE) :: ConstrMax,ConstrRMS,Sum
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
         Sum=ABS(IntDiff(I))
         IF(Sum>ConstrMax) ConstrMax=Sum
         ConstrRMS=ConstrRMS+Sum*Sum
       ENDIF
     ENDDO
     ConstrRMS=SQRT(ConstrRMS/DBLE(NConstr))
   END SUBROUTINE ConstrConv
!
!-------------------------------------------------------------
!
   SUBROUTINE SetCartConstr(Carts,CartDispl,IntCs,NConstr)
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
   END SUBROUTINE SetCartConstr
!
!---------------------------------------------------------------
!
   SUBROUTINE ScaleDispl(CartDispl,MaxCartDiff,DiffMax,RMSD)
     REAL(DOUBLE),DIMENSION(:) :: CartDispl 
     REAL(DOUBLE)              :: MaxCartDiff,DiffMax,RMSD,Sum
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
       Sum=MaxCartDiff/DiffMax
       CartDispl=Sum*CartDispl
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
          .AND.IStep>20) THEN
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
   SUBROUTINE CROSS_PRODUCT(V1,V2,CrossProd,Add_O)
     REAL(DOUBLE),DIMENSION(:) :: V1,V2,CrossProd
     LOGICAL,OPTIONAL          :: Add_O
     REAL(DOUBLE),DIMENSION(3) :: VectAux
     !
     VectAux(1)=V1(2)*V2(3)-V1(3)*V2(2)
     VectAux(2)=V1(3)*V2(1)-V1(1)*V2(3)
     VectAux(3)=V1(1)*V2(2)-V1(2)*V2(1)
     IF(PRESENT(Add_O)) THEN
       IF(Add_O) CrossProd=CrossProd+VectAux
     ELSE
       CrossProd=VectAux
     ENDIF
   END SUBROUTINE CROSS_PRODUCT
!
!------------------------------------------------------------------
!
   SUBROUTINE CALC_BxVect(ISpB,JSpB,ASpB,VectInt,VectCart,Trp_O)
     TYPE(INT_VECT)            :: ISpB,JSpB,ISpBt,JSpBt
     TYPE(DBL_VECT)            :: ASpB,ASpBt
     REAL(DOUBLE),DIMENSION(:) :: VectInt,VectCart
     INTEGER                   :: NCart,NIntC,I,J,JJ,K,KK,LL,NZSpB
     REAL(DOUBLE)              :: Sum
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
                          IPerm1_O,IPerm2_O)
     TYPE(Cholesky)          :: CholData
     CHARACTER(LEN=*)        :: SCRPath
     TYPE(INT_VECT)          :: ISpB,JSpB
     TYPE(DBL_VECT)          :: ASpB
     TYPE(INT_VECT),OPTIONAL :: IPerm1_O,IPerm2_O
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B',&
                    IPerm1_O=IPerm1_O,IPerm2_O=IPerm2_O)
     CALL ReadChol(CholData,TRIM(SCRPath)//'CholFact')
   END SUBROUTINE GetBMatInfo
!
!-------------------------------------------------------------------
!
   SUBROUTINE GetMixedBMat(IntCs,XYZ,GTrfCtrl,GCoordCtrl,TOPS, &
                           NCartConstr,Print,SCRPath,DoConstr)
     TYPE(INTC)                   :: IntCs
     TYPE(TOPOLOGY)               :: TOPS
     TYPE(Cholesky)               :: CholData
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     INTEGER                      :: NatmsLoc,NCart,NIntC,At1,At2,At3
     INTEGER                      :: NMax12,I,J,NCartConstr
     TYPE(BMATR)                  :: B,B1,B2
     INTEGER                      :: Print,ILowDim
     LOGICAL                      :: Print2,DoConstr
     CHARACTER(LEN=*)             :: SCRPath
     TYPE(INT_VECT)               :: IPerm1,IPerm2,InvP,Perm
     TYPE(INT_VECT)               :: ISpB1,ISpB2,JSpB1,JSpB2,IGc1,JGc1
     TYPE(INT_VECT)               :: ISpBMix,JSpBMix
     TYPE(DBL_VECT)               :: ASpB1,ASpB2,ASpBMix,AGc1
     TYPE(DBL_RNK2)               :: XYZWork
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def%C)
     Print2=(Print>=DEBUG_GEOP_MAX)
     !
     NMax12=J-1
     CALL New(XYZWork,(/3,NatmsLoc/))
     XYZWork%D=XYZ
     !
     ! Find two reference systems
     CALL ThreeAtoms(XYZWork%D,IntCs,TOPS%Tot12, &
                     GTrfCtrl%ThreeAt,GTrfCtrl%Linearity)
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
       CALL ThreeAtoms(XYZWork%D,IntCs,TOPS%Tot12, &
                       GTrfCtrl%ThreeAt_2,GTrfCtrl%Linearity)
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
     CALL BMatrix(XYZ,NIntC,IntCs,B, &
                  GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
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
   SUBROUTINE RefreshBMatInfo(IntCs,XYZ,GTrfCtrl,GCoordCtrl,&
                              Print,SCRPath,DoCleanB)
     TYPE(INTC)                   :: IntCs
     TYPE(Cholesky)               :: CholData
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     INTEGER                      :: NatmsLoc,NCart,NIntC
     TYPE(BMATR)                  :: B
     INTEGER                      :: Print
     LOGICAL                      :: Print2,DoCleanB
     CHARACTER(LEN=*)             :: SCRPath
     TYPE(INT_VECT)               :: ISpB,JSpB
     TYPE(DBL_VECT)               :: ASpB
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def%C)
     Print2=(Print>=DEBUG_GEOP_MAX)
     !
     CALL BMatrix(XYZ,NIntC,IntCs,B, &
                  GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
     !
     IF(DoCleanB) CALL CleanBConstr(IntCs,B,NatmsLoc)
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     !
     CALL WriteBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     !
     CALL CholFact(ISpB,JSpB,ASpB,NCart,NIntC, &
                   CholData,Print2,Shift_O=1.D-6)
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
                          TranslAt1,RotAt2ToX,RotAt3ToXY)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     REAL(DOUBLE),DIMENSION(3)    :: TranslAt1
     REAL(DOUBLE),DIMENSION(3,3)  :: RotAt2ToX,RotAt3ToXY
     TYPE(DBL_RNK2)  :: XYZRot
     TYPE(DBL_VECT)  :: Vect1,Vect2,Vect3
     TYPE(DBL_RNK2)  :: Rot
     INTEGER         :: I,J,NatmsLoc,NMax12
     INTEGER         :: At1,At2,At3,ThreeAt(1:3)
     LOGICAL         :: Linearity
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
   SUBROUTINE ThreeAtoms(XYZ,IntCs,Top12,ThreeAt,Linearity)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ
     TYPE(INTC)                   :: IntCs
     TYPE(DBL_VECT)               :: DistVect
     INTEGER                      :: I,J,NatmsLoc,NIntC
     INTEGER                      :: At1,At2,At3,ThreeAt(1:3),ITop
     REAL(DOUBLE)                 :: D12,D13,D23,CX,CY,CZ
     REAL(DOUBLE)                 :: Dist,Sum,Dist12,Area2
     REAL(DOUBLE)                 :: Vect1(1:3),Vect2(1:3),Vect3(1:3)
     LOGICAL                      :: Linearity
     TYPE(DBL_RNK2)               :: XYZRot
     TYPE(INT_RNK2)               :: Top12 
     !
     NatmsLoc=SIZE(XYZ,2)
     NIntC=SIZE(IntCs%Def%C)
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
     CALL ThreeConstr(IntCs,Top12,NatmsLoc,At1,At2,At3)
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
         Sum=XYZRot%D(1,I)**2+XYZRot%D(2,I)**2+XYZRot%D(3,I)**2
         ITop=MAX(Top12%I(I,1),1)
         Sum=Sum*DBLE(ITop**2)
         IF(Sum>Dist) THEN
           Dist=Sum
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
         Sum=DBLE(ITop)*DistVect%D(I)
         IF(Sum>Dist) THEN
           Dist=Sum
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
         Sum=0.5D0*(D12+D13+D23)
         Area2=Sum*(Sum-D12)*(Sum-D13)*(Sum-D23)
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
   SUBROUTINE CleanBConstr(IntCs,B,NatmsLoc)
     TYPE(BMATR) :: B
     TYPE(INTC)  :: IntCs
     INTEGER     :: I,J,K,L,NIntC,JJ,LL,NatmsLoc,NCart
     TYPE(INT_VECT):: IConstr
     !
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def%C)
     CALL New(IConstr,NCart)
     IConstr%I=1  
     DO I=1,NIntC
       IF(IntCs%Constraint%L(I)) THEN
         J=3*(B%IB(I,1)-1)
         IF(IntCs%Def%C(I)(1:5)=='CARTX') IConstr%I(J+1)=0
         IF(IntCs%Def%C(I)(1:5)=='CARTY') IConstr%I(J+2)=0
         IF(IntCs%Def%C(I)(1:5)=='CARTZ') IConstr%I(J+3)=0
       ENDIF
     ENDDO 
     !
     DO I=1,NIntC
       DO J=1,4
         JJ=B%IB(I,J)
         IF(JJ==0) EXIT
         JJ=3*(JJ-1)
         LL=3*(J-1)
         DO L=1,3
           IF(IConstr%I(JJ+L)==0) B%B(I,LL+L)=Zero
         ENDDO
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
         IF(B%IB(I,J)==ThreeAt(1)) THEN
           J2=3*J
           J1=J2-2
           B%B(I,J1:J2)=Zero
         ENDIF
         IF(B%IB(I,J)==ThreeAt(2)) THEN
           J2=3*J
           J1=J2-1
           B%B(I,J1:J2)=Zero
         ENDIF
         IF(B%IB(I,J)==ThreeAt(3)) THEN
           J2=3*J
           B%B(I,J2)=Zero
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE CleanB
!
!-------------------------------------------------------------------
!
   SUBROUTINE PrtXYZ(AtNum,XYZ,FileName,Title,Vects_O)
     REAL(DOUBLE),DIMENSION(:)     :: AtNum
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL   :: Vects_O
     CHARACTER(LEN=*)              :: FileName
     CHARACTER(LEN=*)              :: Title
     CHARACTER(LEN=1)              :: Char 
     INTEGER                       :: I,NatmsLoc
     !
     NatmsLoc=SIZE(XYZ,2)
     OPEN(File=FileName,Unit=99,FORM='FORMATTED',STATUS='UNKNOWN')
     DO 
       READ(99,33,END=1) Char
     ENDDO
     33 format(a1)
     1    CONTINUE
     WRITE(99,*) NatmsLoc 
     WRITE(99,*) Title 
     IF(PRESENT(Vects_O)) THEN
       DO I=1,NatmsLoc
         WRITE(99,100) INT(AtNum(I)),XYZ(1:3,I)/AngstromsToAu,100.D0*Vects_O(1:3,I)
       ENDDO 
     ELSE
       DO I=1,NatmsLoc
         WRITE(99,100) INT(AtNum(I)),XYZ(1:3,I)/AngstromsToAu
       ENDDO 
     ENDIF
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
   SUBROUTINE SetMaxBlkS(Mat)
     TYPE(BCSR) :: Mat
     INTEGER    :: I,II
     ! Set MaxNon0-s for overlap matrices
     II=0
     DO I=1,Mat%NAtms
       II=II+(Mat%RowPt%I(I+1)-Mat%RowPt%I(I))**2
     ENDDO
     MaxBlkS=II+1
   END SUBROUTINE SetMaxBlkS
! 
!-------------------------------------------------------------------
!
   SUBROUTINE ReadBMATR(ISpB,JSpB,ASpB,FileName, &
                        IPerm1_O,IPerm2_O,UMatr_O)
     INTEGER          :: NDim,NZ,NCart,NDim1,NDim2
     TYPE(INT_VECT)   :: ISpB,JSpB
     TYPE(DBL_VECT)   :: ASpB
     TYPE(INT_VECT),OPTIONAL   :: IPerm1_O,IPerm2_O
     TYPE(DBL_RNK2),OPTIONAL   :: UMatr_O
     LOGICAL          :: Exists
     CHARACTER(LEN=*) :: FileName
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
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE ReadBMATR
! 
!-------------------------------------------------------------------
!
   SUBROUTINE WriteBMATR(ISpB,JSpB,ASpB,FileName, &
                         IPerm1_O,IPerm2_O,UMatr_O)
     INTEGER          :: NDim,NZ
     TYPE(INT_VECT)   :: ISpB,JSpB
     TYPE(INT_VECT),OPTIONAL   :: IPerm1_O,IPerm2_O
     TYPE(DBL_RNK2),OPTIONAL   :: UMatr_O
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
     CLOSE(Unit=99,STATUS='KEEP')
   END SUBROUTINE WriteBMATR
! 
!-------------------------------------------------------------------
!
   SUBROUTINE LongRangeIntC(CtrlCoord,Top12,IntCs,NIntC, &
                            XYZ,MarkLongR,AtNum)
     TYPE(INTC)                  :: IntCs,IntC_New
     INTEGER                     :: NIntC,NIntc_New
     TYPE(INT_VECT)              :: LinAtom,MarkLongR,LinCenter
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE)                :: Value,Conv,Dist12
     INTEGER                     :: I1,I2,I3,I4,NMax12,NLinB,NtorsLinb
     INTEGER                     :: II1,II4
     INTEGER                     :: I,J,K,NatmsLoc,III
     INTEGER                     :: NStreLinb,NBendLinb,JJ,IC,IEx
     TYPE(INT_RNK2)              :: LinBBridge,Top12
     TYPE(CoordCtrl)             :: CtrlCoord
     LOGICAL                     :: RepeatChk,DoVdW,Found
     INTEGER,DIMENSION(:)        :: AtNum
     !
     Conv=180.D0/PI
     NIntC=SIZE(IntCs%Def%C)
     NatmsLoc=SIZE(Top12%I,1)
     !
     NMax12=SIZE(Top12%I,2)-1
     !
     CALL NEW(LinBBridge,(/2,NIntc/))
     CALL NEW(LinAtom,NatmsLoc)
     CALL NEW(LinCenter,NatmsLoc)
     LinBBridge%I=0
     LinAtom%I=0
     LinCenter%I=0
     NLinB=0
     DO I=1,NIntc
       IF(IntCs%Def%C(I)(1:5)=='LINB1'.OR.MarkLongR%I(I)/=0) THEN   
         IF(LinAtom%I(IntCs%Atoms%I(I,1))/=0) CYCLE
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
        !DO III=1,NIntC
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
        !DO III=1,NIntC
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
     NIntc_New=NIntC+NTorsLinB+NStreLinB+NBendLinB
     CALL New(IntC_New,NIntc_New)
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
       DO J=1,Top12%I(I2,1)
         I1=Top12%I(I2,1+J)
         IF(LinCenter%I(I)==I1) CYCLE
         DO K=1,Top12%I(I3,1)
           I4=Top12%I(I3,1+K)
           IF(LinCenter%I(I)==I4) CYCLE
           IF(I1/=I4.AND.I1/=I3.AND.I2/=I4) THEN
             IF(.NOT.Found) THEN
               II1=I1
               II4=I4
               Found=.TRUE.
             ELSE
               IF(AtNum(I1)>AtNum(II1)) II1=I1
               IF(AtNum(I4)>AtNum(II4)) II4=I4
             ENDIF
           ENDIF
         ENDDO
       ENDDO
       I1=II1
       I4=II4
       IF(Found) THEN
         NIntC=NIntC+1
         IntC_New%Def%C(NIntC)(1:5)='TORS '
         IntC_New%Atoms%I(NIntC,1)=I1
         IntC_New%Atoms%I(NIntC,2)=I2
         IntC_New%Atoms%I(NIntC,3)=I3
         IntC_New%Atoms%I(NIntC,4)=I4
         IntC_New%Value%D(NIntC)=Zero
         IntC_New%Constraint%L(NIntC)=.FALSE.
         IntC_New%ConstrValue%D(NIntC)=Zero   
         IntC_New%Active%L(NIntC)=.TRUE. 
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
   SUBROUTINE RedundancyOff(Displ,SCRPath,Print)
     REAL(DOUBLE),DIMENSION(:) :: Displ
     REAL(DOUBLE)              :: Perc 
     TYPE(INT_VECT)            :: ISpB,JSpB,IGc,JGc
     TYPE(DBL_VECT)            :: ASpB,AGc
     TYPE(Cholesky)            :: CholData
     INTEGER                   :: NIntC,NCart
     TYPE(DBL_VECT)            :: Vect1,Vect2,Displ2
     TYPE(DBL_RNK2)            :: FullGcInv,FullGc
     INTEGER                   :: Print
     CHARACTER(LEN=*)          :: SCRPath
     !
     CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)  
     NIntC=SIZE(ISpB%I)-1
     IF(NIntC/=SIZE(Displ)) &
       CALL Halt('Dimension error in RedundancyOff')
     NCart=SIZE(CholData%IPerm%I)
     CALL New(Vect1,NCart)
     CALL New(Vect2,NCart)
     CALL New(Displ2,NIntC)
     Displ2%D=Displ
     !
     CALL CALC_BxVect(ISpB,JSpB,ASpB,Displ,Vect1%D,Trp_O=.TRUE.)
     !
     CALL GetGc(NCart,ISpB,JSpB,ASpB,IGc,JGc,AGc)
     CALL New(FullGcInv,(/NCart,NCart/))
     CALL Sp1x1ToFull(IGc%I,JGc%I,AGc%D,NCart,NCart,FullGc)
     !
     CALL SetDSYEVWork(NCart)
     CALL FunkOnSqMat(NCart,Inverse,FullGc%D,FullGcInv%D, &
                      PosDefMat_O=.FALSE.)
     CALL UnSetDSYEVWork()
     !
     CALL DGEMM_NNc(NCart,NCart,1,One,Zero,FullGcInv%D,Vect1%D,Vect2%D)
     CALL Delete(FullGc)
     CALL Delete(FullGcInv)
     !
    !CALL CALC_GcInvCartV(CholData,Vect1%D,Vect2%D)
     !
     CALL CALC_BxVect(ISpB,JSpB,ASpB,Displ,Vect2%D)
     !
     Perc=DOT_PRODUCT(Displ,Displ2%D)/DOT_PRODUCT(Displ2%D,Displ2%D) 
     Perc=(One-ABS(Perc))*100.D0
     IF(Print>=DEBUG_GEOP_MAX) THEN
       WRITE(*,100) Perc
       WRITE(Out,100) Perc
     ENDIF
100  FORMAT("Percentage of Redundancy projected out= ",F8.2)
     !
     CALL Delete(Displ2)
     CALL Delete(Vect1)
     CALL Delete(Vect2)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
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
   LOGICAL FUNCTION HasAngle(Char)
     CHARACTER(LEN=*) :: Char
     HasAngle=(Char(1:4)=='BEND'.OR. &
               Char(1:4)=='LINB'.OR. &
               Char(1:4)=='OUTP'.OR. &
               Char(1:4)=='TORS')
   END FUNCTION HasAngle
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasBendLinB(Char)
     CHARACTER(LEN=*) :: Char
     HasBendLinB=(Char(1:4)=='BEND'.OR. &
                  Char(1:4)=='LINB')
   END FUNCTION HasBendLinB
!
!----------------------------------------------------------------------
!
   LOGICAL FUNCTION HasTorsOutP(Char)
     CHARACTER(LEN=*) :: Char
     HasTorsOutP=(Char(1:4)=='OUTP'.OR. &
                  Char(1:4)=='TORS')
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
   FUNCTION HasHBond(Top12,AtNum,NJJ1,NJJ2,JJ1,JJ2,HAtm) 
     LOGICAL              :: HasHBond
     TYPE(INT_RNK2)       :: Top12
     INTEGER,DIMENSION(:) :: AtNum
     INTEGER              :: NJJ1,NJJ2,JJ1,JJ2,HAtm
     !
     HasHBond=.FALSE.
     HAtm=0
     IF(NJJ1/=1.AND.NJJ2/=1) RETURN
     IF((NJJ1==1.AND.HasLigand(NJJ2))) THEN
      !HasHBond=HasAttached(AtNum,Top12%I,JJ1)
       HasHBond=.TRUE.
       HAtm=JJ1
     ELSE IF((NJJ2==1.AND.HasLigand(NJJ1))) THEN
      !HasHBond=HasAttached(AtNum,Top12%I,JJ2)
       HasHBond=.TRUE.
       HAtm=JJ2
     ENDIF
   END FUNCTION HasHBond
!
!---------------------------------------------------------------------
!
   FUNCTION HasAttached(AtNum,Top12,JJ1)
     LOGICAL                 :: HasAttached
     INTEGER,DIMENSION(:,:)  :: Top12
     INTEGER,DIMENSION(:)    :: AtNum
     INTEGER                 :: JJ1,J,K
     !
     HasAttached=.FALSE.
     DO J=1,Top12(JJ1,1)
       K=AtNum(Top12(JJ1,J+1))
       IF(HasLigand(K)) THEN
         HasAttached=.TRUE.
         EXIT  
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
     INTEGER                     :: I,J,K,L,I1,I2,NDim
     !
     MaxBonds=30
     CondNumb1=0.5D0
     CondNumb2=0.8D0
     NatmsLoc=SIZE(XYZ,2)
     !
     CALL New(AllowBond,Bond%N)
     CALL New(BondScore,Bond%N)
     CALL ScoreBond(BondScore%I,Bond)
     !
     AllowBond%I=1
     DO I1=1,NatmsLoc
       NDim=AtmB%Count%I(I1)
       IF(NDim==0) CYCLE
       CALL New(BondVect2,(/NDim,3/))
       CALL New(BondVect3,NDim)
       CALL New(Ordered,NDim)
       CALL D3Bonds(I1,XYZ,BondVect2%D,BondScore%I,AtmB,Bond)
       BondVect3%D=Zero
       DO J=1,3
         DO K=1,NDim 
           BondVect3%D(K)=BondVect3%D(K)+BondVect2%D(K,J)**2
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
         IF(Dist>CondNumb2) THEN
           AllowBond%I(L)=0
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
     CALL SetDSYEVWork(3)
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
     CALL Halt('DSYEV hosed in SortNonCov2. INFO='&
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
       Dist=SQRT(ABS(EigVals(J)))
       DO K=1,NDim 
         BondVect2(K,J)=Dist*BondVect2(K,J) 
       ENDDO
     ENDDO
     !
     CALL Delete(BondVect)
     CALL UnSetDSYEVWork()
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
   SUBROUTINE ReSetConstr(IntCs,XYZ)
     TYPE(INTC)                  ::  IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER                     :: I,J,NIntC
     !
     NIntC=SIZE(IntCs%Def%C)
     DO I=1,NIntC
       IF(IntCs%Constraint%L(I)) THEN
         J=IntCs%Atoms%I(I,1)
         IF(IntCs%Def%C(I)(1:5)=='CARTX') IntCs%ConstrValue%D(I)=XYZ(1,J)
         IF(IntCs%Def%C(I)(1:5)=='CARTY') IntCs%ConstrValue%D(I)=XYZ(2,J)
         IF(IntCs%Def%C(I)(1:5)=='CARTZ') IntCs%ConstrValue%D(I)=XYZ(3,J)
       ENDIF 
     ENDDO
   END SUBROUTINE ReSetConstr
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
     REAL(DOUBLE)                :: AngleSum,Sum,Conv,Planar,ASum,TwoPi
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
           Sum=A1+A2+A3
           ASum=ABS(TwoPI-Sum)
           IF(Planar>ASum) THEN
             Planar=ASum
             AngleSum=Sum
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
     NIntC=SIZE(B%IB,1)
     DO I=1,NIntC
       DO J=1,4
         IF(B%IB(I,J)==0) CYCLE
         JJ=3*(J-1)
         Vect1=B%B(I,JJ+1:JJ+3)
         CALL DGEMM_NNc(3,3,1,One,Zero,Rot,Vect1,Vect2)
         B%B(I,JJ+1:JJ+3)=Vect2
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
     REAL(DOUBLE)  :: Angle,AngleDispl,Sum
     Sum=Angle+AngleDispl
     IF(Sum>PI) AngleDispl=(TwoPi-Sum)-Angle
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
     INTEGER          :: NCart,INFO,I,J
     TYPE(DBL_RNK2)   :: FullGc,UMatr,Aux1,Aux2
     REAL(DOUBLE)     :: Fact
     !
     CALL New(UMatr,(/NCart,NCart-6/))
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     CALL GetGc(NCart,ISpB,JSpB,ASpB,IGc,JGc,AGc)
     CALL Sp1x1ToFull(IGc%I,JGc%I,AGc%D,NCart,NCart,FullGc)
     !
     CALL SetDSYEVWork(NCart)
       BLKVECT%D=FullGc%D
       CALL DSYEV('V','U',NCart,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in BuildUMatr. INFO='&
                  //TRIM(IntToChar(INFO)))
       DO I=7,NCart
         Fact=One/SQRT(BLKVALS%D(I))
         DO J=1,NCart ; UMatr%D(J,I-6)=BLKVECT%D(J,I)*Fact ; ENDDO
       ENDDO     
     CALL UnSetDSYEVWork()
     !
   ! call pprint(UMatr,'UMatr',Unit_O=6)
   ! CALL New(Aux1,(/NCart-6,NCart/))
   ! CALL New(Aux2,(/NCart-6,NCart-6/))
   ! !     CALL PPrint(SpB,'SpB2',Unit_O=6)
   ! CALL DGEMM_TNc(NCart-6,NCart,NCart,One,Zero,UMatr%D,&
   !                FullGc%D,Aux1%D)
   ! CALL DGEMM_NNc(NCart-6,NCart,NCart-6,One,Zero,Aux1%D, &
   !                UMatr%D,Aux2%D)
   ! call pprint(Aux2,'Unit Matr?',Unit_O=6)
   ! CALL Delete(Aux1)
   ! CALL Delete(Aux2)
     !
     CALL WriteBMATR(ISpB,JSpB,ASpB, &
                     TRIM(SCRPath)//'UMatr',UMatr_O=UMatr)
     !
     CALL Delete(FullGc)
     CALL Delete(UMatr)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(IGc)
     CALL Delete(JGc)
     CALL Delete(AGc)
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
       DO J=I,NDim-1
         I1=IWork(J)
         I2=IWork(J+1)
         X1=VectX(I1)
         X2=VectX(I2) 
         IF(X1<X2) THEN
           IWork(J)=I2
           IWork(J+1)=I1
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE ReorderN 
!
!--------------------------------------------------------------------
!
   SUBROUTINE ReorderI(VectX,IWork)
     REAL(DOUBLE),DIMENSION(:) :: VectX
     INTEGER     ,DIMENSION(:) :: IWork
     REAL(DOUBLE)              :: X1,X2
     INTEGER                   :: I,J,NDim,I1,I2
     ! order set by increasing x
     NDim=SIZE(VectX)
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
   SUBROUTINE IntCBoxes(XYZ,Box)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(IntCBox)               :: Box
     INTEGER                     :: NatmsLoc
     REAL(DOUBLE)                :: BXMIN,BYMIN,BZMIN
     REAL(DOUBLE)                :: BoxSize
     !
     NatmsLoc=SIZE(XYZ,2)
     BoxSize=3.0D0*AngstromsToAU !in A
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
   SUBROUTINE BondingScheme(XYZ,AtNum,IntSet,AtmB,Bond,TOPS,Box,&
                            GCoordCtrl)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     INTEGER,DIMENSION(:)        :: AtNum
     INTEGER                     :: IntSet 
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(IntCBox)               :: Box 
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(DBL_VECT)              :: CritRad
     REAL(DOUBLE)                :: Fact,HBondMax
     INTEGER                     :: I,J,N,NatmsLoc
     LOGICAL                     :: DoRepeat
     !
     !now define bonding scheme, based on Slater or Van der Waals radii
     !
     NatmsLoc=SIZE(XYZ,2)
     HBondMax=2.40D0*AngstromsToAu*MAX(GCoordCtrl%VDWFact,One) ! in Au 
     !
     IF(IntSet==1) THEN
       N=SIZE(SLRadii,1)
       CALL New(CritRad,NatmsLoc)
       Fact=1.3D0 !!! Scaling factor for Slater Radii
       DO I=1,NatmsLoc
         CritRad%D(I)=Fact*SLRadii(AtNum(I))*AngstromsToAU
       ENDDO
     ELSE
       N=SIZE(VDWRadii,1)
       CALL New(CritRad,NatmsLoc)
       Fact=GCoordCtrl%VDWFact !!! Scaling factor for VDW Radii
       DO I=1,NatmsLoc
         CritRad%D(I)=Fact*VDWRadii(AtNum(I))*AngstromsToAU
       ENDDO
     ENDIF
     !
     IF(IntSet==1) THEN
       DoRepeat=.FALSE.
       CALL BondList(XYZ,AtNum,IntSet,Box,Bond,TOPS, &
                     CritRad,HbondMax,DoRepeat)
       CALL SortBonds(NatmsLoc,AtmB,Bond)
       CALL Topology_12(AtmB,TOPS%Cov12)
       CALL Topology_13(NatmsLoc,TOPS%Cov12,TOPS%Cov13)
       CALL Topology_14(NatmsLoc,TOPS%Cov12,TOPS%Cov14)
       CALL Excl_List(NatmsLoc,TOPS%Cov12,TOPS%Cov13,TOPS%Cov14, &
                      TOPS%CovExcl)
     ELSE
       DO I=1,6
         DoRepeat=.FALSE.
         CALL BondList(XYZ,AtNum,IntSet,Box,Bond,TOPS, &
                       CritRad,HBondMax,DoRepeat)
         IF(DoRepeat) THEN
           CALL Delete(Bond)
         ELSE
           EXIT
         ENDIF
         IF(I==6.AND.DoRepeat) CALL Halt('The 6th attempt of recognizing bonds to lonely atoms failed.')
       ENDDO
       CALL SortBonds(NatmsLoc,AtmB,Bond)
       CALL SortNonCov2(AtNum,XYZ,Bond,AtmB)
       CALL Delete(AtmB)
       CALL SortBonds(NatmsLoc,AtmB,Bond)
       CALL VDWTop(TOPS%Tot12,TOPS%Cov12,Bond%IJ)
       CALL Topology_13(NatmsLoc,TOPS%Tot12,TOPS%Tot13)
       CALL Topology_14(NatmsLoc,TOPS%Tot12,TOPS%Tot14)
       CALL Excl_List(NatmsLoc,TOPS%Tot12,TOPS%Tot13,TOPS%Tot14, &
                      TOPS%TotExcl)
     ENDIF
     CALL Delete(CritRad)
   END SUBROUTINE BondingScheme
!
!--------------------------------------------------------
!
   SUBROUTINE BondList(XYZ,AtNum,IntSet,Box,Bond,TOPS, &
                       CritRad,HBondMax,DoRepeat)
     IMPLICIT NONE
     INTEGER                     :: I,J,NatmsLoc,NBond
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(BONDDATA)              :: Bond
     TYPE(IntCBox)               :: Box
     TYPE(TOPOLOGY)              :: TOPS
     INTEGER,DIMENSION(:)        :: AtNum
     TYPE(DBL_VECT)              :: CritRad,CritRadMod 
     REAL(DOUBLE)                :: R12,R12_2,CritDist,HBondMax
     INTEGER                     :: IZ,IX,IY,I1,I2,JJ1,JJ2
     INTEGER                     :: IORD,IORDD
     INTEGER                     :: IZD,IXD,IYD,NJJ1,NJJ2
     INTEGER                     :: NMax12,JJ,IntSet,IDimExcl,NBondEst
     INTEGER                     :: IDim14,IDim12,IDim13
     INTEGER                     :: HAtm,LAtm,DDim,NBondOld
     INTEGER                     :: IChk,MaxBonds
     REAL(DOUBLE),DIMENSION(3)   :: DVect
     LOGICAL                     :: FoundHBond,FoundMetLig
     LOGICAL                     :: LonelyAtom,DoExclude,DoRepeat
     !     
     NatmsLoc=SIZE(XYZ,2)
     MaxBonds=20 ! initial value
     NBondEst=NatmsLoc*MaxBonds
     NBond=NBondEst
     HAtm=0
     CALL New(Bond,NBondEst)
     CALL New(CritRadMod,NatmsLoc)
     CritRadMod%D=CritRad%D
     !
     !  Go through all boxes and their neighbours
     !
     NBond=0
     DO IZ=1,Box%NZ
       DO IX=1,Box%NX
         DO IY=1,Box%NY
           ! absolute index of a box
           IOrd=Box%NX*Box%NY*(IZ-1)+Box%NY*(IX-1)+IY
           NBondOld=NBond
           DO I1=Box%I%I(IOrd),Box%I%I(IOrd+1)-1
             JJ1=Box%J%I(I1) !!! atom in central box
             NJJ1=AtNum(JJ1)
             DO IChk=1,NatmsLoc
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
                             FoundHBond=.FALSE.
                             FoundMetLig=.FALSE.
                             LonelyAtom=.FALSE.
                             IF(JJ2<=JJ1) CYCLE !!!avoid double counting
                             IF(IntSet==2) THEN
                               FoundMetLig=HasMetLig(JJ1,JJ2,NJJ1,NJJ2) 
                               LonelyAtom=HasLonelyAtm(TOPS%Cov12,JJ1,JJ2,LAtm)
                               FoundHBond=HasHBond(TOPS%Cov12,AtNum, &
                                            NJJ1,NJJ2,JJ1,JJ2,HAtm)
                               CALL BondExcl(JJ1,JJ2,NJJ1,NJJ2,TOPS, &
                                             FoundHBond,FoundMetLig,&
                                             LonelyAtom,DoExclude)
                               IF(DoExclude) CYCLE
                             ENDIF
                             IF(FoundHBond.AND..NOT.LonelyAtom) THEN
                               CritDist=HBondMax
                             ELSE
                               CritDist=CritRadMod%D(JJ1)+CritRad%D(JJ2)
                             ENDIF
                             DVect(:)=XYZ(:,JJ1)-XYZ(:,JJ2)
                             R12_2=DOT_PRODUCT(DVect,DVect)
                             R12=SQRT(R12_2)
                             IF(R12<CritDist) THEN
                               IF(NBond+1>NBondEst) THEN
                                 DDim=NatmsLoc*10
                                 CALL MoreBondArray(Bond,NBondEst+DDim, &
                                                    NBondEst)
                                 NBondEst=NBondEst+DDim
                               ENDIF 
                               NBond=NBond+1
                               Bond%IJ%I(1:2,NBond)=(/JJ1,JJ2/)
                               Bond%Length%D(NBond)=R12
                               !
                               IF(IntSet==1) THEN
                                 Bond%Type%C(NBond)(1:3)='COV'
                               ELSE
                                 Bond%Type%C(NBond)(1:3)='VDW'
                               ENDIF
                               Bond%HBExtraSN%I(NBond)=0    
                               Bond%HBExtraNC%I(NBond)=0    
                               Bond%LonelyAtom%I(NBond)=0    
                               IF(FoundHBond.AND..NOT.LonelyAtom) THEN
                                 Bond%Type%C(NBond)(1:5)='HBond'
                                 IF(TOPS%Cov12%I(HAtm,1)==1) THEN
                                   Bond%HBExtraSN%I(NBond)=HAtm 
                                   Bond%HBExtraNC%I(NBond)=AtNum(HAtm)
                                 ENDIF
                               ENDIF
                               IF(LonelyAtom) THEN
                                 Bond%LonelyAtom%I(NBond)=1    
                               ENDIF
                               IF(FoundMetLig) THEN
                                 Bond%Type%C(NBond)(1:6)='MetLig'
                               ENDIF
                             ENDIF
                           ENDDO
                         ENDIF
                       ENDDO
                     ENDIF
                   ENDDO
                 ENDIF
               ENDDO
               IF(IntSet==2) THEN
                 IF(NBondOld==NBond.AND.TOPS%Cov12%I(JJ1,1)==0) THEN
                   DoRepeat=.TRUE.
                   CritRadMod%D(JJ1)=1.2D0*CritRadMod%D(JJ1)
                 ELSE
                   EXIT  
                 ENDIF
               ELSE
                 EXIT
               ENDIF
             ENDDO !!! IChk for checking LonelyAtoms
           ENDDO !!! central box atoms
         ENDDO
       ENDDO
     ENDDO !!! ends on central box indices
     !
     ! Compress Bond
     !
     CALL MoreBondArray(Bond,NBond,NBond)
     IF(DoRepeat) THEN
       CritRad%D=CritRadMod%D
     ENDIF
     CALL Delete(CritRadMod)
   END SUBROUTINE BondList
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
     TYPE(INT_RNK2)   :: ISign
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
     CALL New(ISign,(/2,NatmsLoc/))
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
       ISign%I(1,I)=IORD
       ISign%I(2,I)=BoxCounter%I(IX,IY,IZ) 
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
       IORD=ISign%I(1,I)
       IADD=ISign%I(2,I)
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
     CALL Delete(ISign)
     CALL Delete(BoxCounter)
     CALL Delete(BoxI)
     CALL Delete(BoxJ)
   END SUBROUTINE SORT_INTO_Box2
!
!----------------------------------------------------------------
!
  !SUBROUTINE Topologies_MM(NatmsLoc,NBond,BondI,BondJ,SCR)
  !  !
  !  ! Set up a table which shows the atom numbers of Atoms 
  !  ! connected to a certain atom by the input bonds (Topology mtr)
  !  ! Here, the generation of the Topology mtr is based 
  !  ! EXCLUSIVELY on input list of bonds
  !  ! WARNING! This subroutine is mainly used to set up the
  !  ! topology files necessary for the calculation of 
  !  ! MM exclusion energies!!!
  !  !
  !  IMPLICIT NONE
  !  TYPE(INT_RNK2)                 :: Top12,Top13,Top14
  !  TYPE(INT_RNK2)                 :: Top_Excl,Top_Excl14
  !  INTEGER                        :: NBond,NatmsLoc
  !  INTEGER,DIMENSION(1:NBond)     :: BondI,BondJ
  !  CHARACTER(LEN=DCL)             :: PathName
  !  CHARACTER(LEN=*)               :: SCR
  !  !
  !  PathName=TRIM(SCR)//'MM'
  !  CALL Topology_12(NatmsLoc,NBond,BondI,BondJ,Top12,PathName)
  !  CALL Topology_13(NatmsLoc,Top12,Top13,PathName)
  !  CALL Topology_14(NatmsLoc,Top12,Top14,PathName)
  !  CALL Excl_List(NatmsLoc,Top12,Top13,Top14,Top_Excl,PathName)
  !  CALL Delete(Top_Excl)
  !  CALL Excl_List14(NatmsLoc,Top12,Top13,Top14,Top_Excl14)
  !  CALL Delete(Top_Excl14)
  !  CALL Delete(Top12)
  !  CALL Delete(Top13)
  !  CALL Delete(Top14)
  !END SUBROUTINE Topologies_MM
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
   SUBROUTINE VDWTop(TopVDW,Top12,BondIJ)
     TYPE(INT_RNK2) :: Top12,TopVDW,BondIJ,TopVDWNew
     TYPE(INT_VECT) :: Addition
     INTEGER        :: NatmsLoc,NBond,I,J,Dim2,AddDim,I1,I2,II1,II2
     !
     NatmsLoc=SIZE(Top12%I,1)
     Dim2=SIZE(Top12%I,2)
     CALL New(TopVDW,(/NatmsLoc,Dim2/))
     TopVDW%I=Top12%I
     NBond=SIZE(BondIJ%I,2)
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
   SUBROUTINE AngleList(AtmB,Bond,TOPS,XYZ,Angle,OutP)
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     TYPE(TOPOLOGY)              :: TOPS
     TYPE(ANGLEDATA)             :: Angle
     TYPE(OUTPDATA)              :: OutP 
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INT_VECT)              :: BondScore
     INTEGER,DIMENSION(3)        :: RefBonds
     INTEGER                     :: I,J,K,L,NAngle,NatmsLoc,AllBond
     INTEGER                     :: I1,I2,I3,NDimens,NOutP
     !
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
         CALL Halt('Atom not bound at '//TRIM(IntToChar(I))//' .')
       ENDIF
       CALL AnglesRef(RefBonds,NDimens,BondScore%I,I,XYZ,AtmB,Bond)
      !CALL AngleGen(I,Angle,NAngle,RefBonds,NDimens,XYZ,AtmB,Bond,TOPS)
       CALL FullAngleGen(I,Angle,NAngle,AtmB,Bond,TOPS)
       CALL OutPGen(I,OutP,NOutP,RefBonds,NDimens,XYZ,AtmB,Bond,TOPS)
     ENDDO
     CALL Delete(BondScore)
     !
     Angle%N=NAngle
     OutP%N=NOutP
   END SUBROUTINE AngleList
!
!----------------------------------------------------------------------
!
   SUBROUTINE OutPGen(IAt,OutP,NOutP,RefBonds,NDimens,XYZ,AtmB,Bond,TOPS)
     INTEGER                       :: IAt,I,J,K,L,NOutP,NDimens
     TYPE(OUTPDATA)                :: OutP
     INTEGER,DIMENSION(:)          :: RefBonds
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ
     REAL(DOUBLE)                  :: D,DMax
     TYPE(ATOMBONDS)               :: AtmB
     TYPE(BONDDATA)                :: Bond
     TYPE(TOPOLOGY)                :: TOPS
     TYPE(INT_VECT)                :: Mark
     INTEGER                       :: NDim,IBonds(3),I1,I2,IM
     !
     IF(NDimens==2.AND.TOPS%Tot12%I(IAt,1)>2) THEN
       NDim=AtmB%Count%I(IAt)
       CALL New(Mark,NDim)
       Mark%I=0
       IBonds=0
       DO L=1,3   
         DMax=1.D99
         DO J=1,NDim
           K=AtmB%Bonds%I(IAt,J)
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
       IF(I1>I2) THEN
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
   SUBROUTINE FullAngleGen(IAt,Angle,NAngle,AtmB,Bond,TOPS)
     INTEGER                :: IAt,I,J,K,L,NAngle
     TYPE(ANGLEDATA)        :: Angle
     TYPE(ATOMBONDS)        :: AtmB
     TYPE(BONDDATA)         :: Bond 
     TYPE(TOPOLOGY)         :: TOPS
     INTEGER                :: B1,B2,I1B1,I2B1,I1B2,I2B2,II1,II2,TopDim
     !
     DO J=1,AtmB%Count%I(IAt)
       B1=AtmB%Bonds%I(IAt,J)
       I1B1=Bond%IJ%I(1,B1)
       I2B1=Bond%IJ%I(2,B1)
       IF(I1B1==IAt) THEN
         II1=I2B1
       ELSE
         II1=I1B1
       ENDIF
       DO K=J+1,AtmB%Count%I(IAt)
         B2=AtmB%Bonds%I(IAt,K) 
         I1B2=Bond%IJ%I(1,B2)
         I2B2=Bond%IJ%I(2,B2)
         IF(I1B2==IAt) THEN
           II2=I2B2
         ELSE
           II2=I1B2
         ENDIF
         TopDim=TOPS%Tot12%I(II1,1)
         IF(ANY(TOPS%Tot12%I(II1,2:TopDim+1)==II2)) CYCLE
         NAngle=NAngle+1
         Angle%IJK%I(2,NAngle)=IAt
         IF(II1<II2) THEN
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
     CondNum1=0.01D0
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
       CALL MergeBondSets(Bond,Bond2,AtmB,AtmB2)
       CALL Delete(AtmB)
       CALL SortBonds(NatmsLoc,AtmB,Bond)
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
     IF(Def=='STRE') THEN
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
   END MODULE InCoords

