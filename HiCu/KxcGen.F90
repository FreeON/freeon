MODULE KxcGen
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
  USE CubeTree 
  USE TreeWalk
  IMPLICIT NONE
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    SUBROUTINE MakeKxc(Kxc,CubeRoot)
      TYPE(BCSR)             :: Kxc
      TYPE(CubeNode),POINTER :: CubeRoot
!
      TYPE(AtomPair)         :: Pair
      INTEGER                :: AtA,AtB
      INTEGER                :: JP,K,NA,NB,NAB,P,Q,R,I,J
#ifdef PERIODIC        
      INTEGER                :: NCA,NCB
      REAL(DOUBLE)           :: Ax,Ay,Az,Bx,By,Bz
      LOGICAL                :: NotMakeBlock
#endif     
!---------------------------------------------- 
!     Initialize the matrix and associated indecies
      P=1; 
      R=1; 
      Kxc%RowPt%I(1)=1
      CALL SetEq(Kxc%MTrix,Zero)
!     Loop over atom pairs
      Kxc%NAtms=NAtoms
      DO AtA=1,NAtoms
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
               IF(AtB<=AtA)THEN
!              Compute only the lower triangle of symmetric Kxc
#ifdef PERIODIC
                  Ax = Pair%A(1)
                  Ay = Pair%A(2)
                  Az = Pair%A(3)
                  Bx = Pair%B(1)
                  By = Pair%B(2)           
                  Bz = Pair%B(3)
                  NotMakeBlock = .TRUE.
                  DO NCA = 1,CS%NCells
                     Pair%A(1) = Ax+CS%CellCarts%D(1,NCA)
                     Pair%A(2) = Ay+CS%CellCarts%D(2,NCA)
                     Pair%A(3) = Az+CS%CellCarts%D(3,NCA)
                     DO NCB = 1,CS%NCells
                        Pair%B(1) = Bx+CS%CellCarts%D(1,NCB)
                        Pair%B(2) = By+CS%CellCarts%D(2,NCB)
                        Pair%B(3) = Bz+CS%CellCarts%D(3,NCB)
                        Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                                  + (Pair%A(2)-Pair%B(2))**2 &
                                  + (Pair%A(3)-Pair%B(3))**2
                        IF(TestAtomPair(Pair)) THEN
                           Kxc%MTrix%D(R:R+NAB-1)=Kxc%MTrix%D(R:R+NAB-1)+KxcBlock(Pair,CubeRoot)
                           NotMakeBlock = .FALSE.
                        ENDIF
                     ENDDO
                  ENDDO
#else
                  Kxc%MTrix%D(R:R+NAB-1)=KxcBlock(Pair,CubeRoot)
#endif
               ENDIF
               Kxc%ColPt%I(P)=AtB
               Kxc%BlkPt%I(P)=R
               R=R+NAB
               P=P+1 
               Kxc%RowPt%I(AtA+1)=P        
               IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in Kxc ')
#ifdef PERIODIC
               IF(NotMakeBlock) &
                  CALL Halt(' Making a Zero Block in KxcGen ')
#endif
            ENDIF
         ENDDO
      ENDDO
      Kxc%NBlks=P-1
      Kxc%NNon0=R-1
!     Fill the upper triangle of Kxc
      DO I=1,NAtoms
         DO JP=Kxc%RowPt%I(I),Kxc%RowPt%I(I+1)-1
            J=Kxc%ColPt%I(JP)
            IF(I>J)THEN
               DO K=Kxc%RowPt%I(J),Kxc%RowPt%I(J+1)-1
                  IF(Kxc%ColPt%I(K)==I)THEN
                     Q=Kxc%BlkPt%I(K)                     
                     EXIT
                  ENDIF
               ENDDO               
               P=Kxc%BlkPt%I(JP)
               NA=BS%BFKnd%I(GM%AtTyp%I(I))
               NB=BS%BFKnd%I(GM%AtTyp%I(J))
               NAB=NA*NB
               CALL XPose(NA,NB,Kxc%MTrix%D(P:P+NAB-1),Kxc%MTrix%D(Q:Q+NAB-1))
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE MakeKxc
!=============================================================================================
!
!=============================================================================================
     SUBROUTINE XPose(M,N,A,AT)
       INTEGER                      :: I,J,M,N,IDex,JDex
       REAL(DOUBLE), DIMENSION(M*N) :: A,AT
       DO I=1,N
          DO J=1,M
             IDex=(I-1)*M+J
             JDex=(J-1)*N+I
             AT(JDex)=A(IDex)
          ENDDO
       ENDDO
     END SUBROUTINE XPose
!=======================================================================================
!
!=======================================================================================
     FUNCTION KxcBlock(Pair,CubeRoot) RESULT(Kvct)
       TYPE(AtomPair)                           :: Pair
       TYPE(CubeNode), POINTER                  :: CubeRoot
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: Kblk
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Kvct
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LMN,EllA,EllB
!-------------------------------------------------------------------------------------- 
       KBlk=Zero
       KA=Pair%KA
       KB=Pair%KB
!----------------------------------
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%AB2=Pair%AB2
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
!----------------------------------
       IndexA=0                  
       DO CFA=1,BS%NCFnc%I(KA)    
       DO CFB=1,BS%NCFnc%I(KB) 
          IndexA  = CFBlokDex(BS,CFA,KA)                
          IndexB  = CFBlokDex(BS,CFB,KB)                  
          StartLA = BS%LStrt%I(CFA,KA)        
          StartLB = BS%LStrt%I(CFB,KB)
          StopLA  = BS%LStop%I(CFA,KA)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLA   = BS%ASymm%I(2,CFA,KA)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
!----------------------------------
          Prim%CFA=CFA
          Prim%CFB=CFB
          Prim%Ell=MaxLA+MaxLB
!----------------------------------
          DO PFA=1,BS%NPFnc%I(CFA,KA)          
          DO PFB=1,BS%NPFnc%I(CFB,KB)
!----------------------------------
             Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
             Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
             Prim%Zeta=Prim%ZA+Prim%ZB
             Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
             IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                Prim%PFA=PFA 
                Prim%PFB=PFB
!               Set primitive values
                MaxAmp=SetBraBlok(Prim,BS)
!---------------------------------------------------------------------------
                MaxAmp=LOG(MaxAmp)
!               Evaluate this primitives Ket contribution to Kxc_ab
                Ket=Zero
                PBox%BndBox(1,:)=Prim%P(1)
                PBox%BndBox(2,:)=Prim%P(2)
                PBox%BndBox(3,:)=Prim%P(3)
                PExtent=GaussianExtent(Prim%Zeta,MaxAmp)
                PBox=ExpandBox(PBox,PExtent)
!               Walk the walk
                CALL KxcWalk(CubeRoot)
!               Contract <Bra|Ket> bloks to compute matrix elements of Kxc
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                      DO LMN=1,LHGTF(EllA+EllB)
                         Kblk(IA,IB)=Kblk(IA,IB)+HGBra%D(LMN,IA,IB)*Ket(LMN)
                      ENDDO
                   ENDDO
                ENDDO
!---------------------------------------------------------------------------
             ENDIF 
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
       Kvct=BlockToVect(Pair%NA,Pair%NB,Kblk)
     END FUNCTION KxcBlock
!====================================================================================================
!
!====================================================================================================
END MODULE KxcGen
