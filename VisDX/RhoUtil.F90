!=================================================================================
!  
!=================================================================================
MODULE RhoUtil
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE Macros
   USE AtomPairs
!
   CONTAINS
!---------------------------------------------------------------------------
     SUBROUTINE RhoCubed(Args,Del,Origin,Nx,Ny,Nz)
!
       USE RhoTree
       USE HiCuThresholds
       USE CubeTree
!
       TYPE(ARGMT)                :: Args
       INTEGER                    :: I,J,K,Nx,Ny,Nz,NP
       REAL(DOUBLE),DIMENSION(3)  :: Origin,R
       REAL(DOUBLE)               :: Del,Den
       REAL(DOUBLE)               :: AA=One/AngstromsToAU
#ifdef PERIODIC
       INTEGER                    :: NC
#endif
      CHARACTER(LEN=DCL) :: Name	
!      SET DENSITY THRESHOLDS (LOOSE)
       Del=3.D-1
       TauRho=1.D-3
!      BUILD RHO-TREE
       CALL RhoToTree(Args)
       CALL SetBBox(RhoRoot%Box,Box)
#ifdef PERIODIC
       CALL MakeBoxPeriodic(Box)
#endif
!      DETERMINE CUBE GRID (SMART, UNSTRUCTURED GRIDS DONT WORK WITH DX) 
       Origin=Box%BndBox(:,1)
       Nx=CEILING((Box%BndBox(1,2)-Box%BndBox(1,1))/Del)
       Ny=CEILING((Box%BndBox(2,2)-Box%BndBox(2,1))/Del)
       Nz=CEILING((Box%BndBox(3,2)-Box%BndBox(3,1))/Del)
!      WRITE DENSITY TO FILE
       NAME=TrixFile('RhoCubes',Args)
       WRITE(*,*)' DENSITY WRITTEN TO '//TRIM(NAME)	
       CALL OpenASCII(NAME,66)
       WRITE(66,*)' object 1 class array items ',NX*NY*NZ,' data follows '
       DO I=1,Nx    
          DO J=1,Ny
             DO K=1,Nz
#ifdef PERIODIC
                Den=Zero
                DO NC=1,CS_OUT%NCells
!                    Set Atomic Coordinates
                   R=(/I,J,K/)*Del+Origin+CS_OUT%CellCarts%D(:,NC)
                   Den=Den+RhoAtPoint(RhoRoot,R)
                ENDDO
#else
                R=(/I,J,K/)*Del+Origin
                Den=RhoAtPoint(RhoRoot,R)
#endif
                WRITE(66,7) Den
7               FORMAT(8(2x,F12.8)) 
             ENDDO
          ENDDO
       ENDDO
       WRITE(66,*)' '
       WRITE(66,*)' attribute "dep" string "positions" '
       WRITE(66,*)'object 2 class gridpositions counts',NX,NY,NZ
       WRITE(66,*)' origin ',origin*AA
       WRITE(66,*)' delta ',del*AA,' 0.0  0.0 '
       WRITE(66,*)' delta 0.0 ',del*AA,' 0.0 '
       WRITE(66,*)' delta 0.0  0.0 ',del*AA
       WRITE(66,*)'object 3 class gridconnections counts ',NX,NY,NZ
       WRITE(66,*)' attribute "element type" string "cubes" '
       WRITE(66,*)' attribute "ref" string "positions" '
       WRITE(66,*)'object "electron density" class field '
       WRITE(66,*)' component "data" 1 '
       WRITE(66,*)' component "positions" 2 '
       WRITE(66,*)' component "connections" 3 '
       CLOSE(UNIT=66)
!      Delete the RhoTree
       CALL DeleteRhoTree(RhoRoot)
!
     END SUBROUTINE RhoCubed
!=================================================================================
!     Sums the significant contributions (leaves) to the density at a point
!=================================================================================
     RECURSIVE FUNCTION RhoAtPoint(Node,R)
       USE RhoTree
       USE BoundingBox
       REAL(DOUBLE)                               :: RhoAtPoint
       REAL(DOUBLE),DIMENSION(3)                  :: R
       TYPE(RhoNode), POINTER                     :: Node
       REAL(DOUBLE), DIMENSION(0:HGEll+1)         :: LLambdaX,LLambdaY,LLambdaZ, &
                                                     ULambdaX,ULambdaY,ULambdaZ, &
                                                     LambdaX,LambdaY,LambdaZ
       REAL(DOUBLE)                               :: RQx,RQy,RQz,RQ2,Z,X,W,Sgn,Xpt,Co,        &
                                                     LQx,LQy,LQz,LQ2,UQx,UQy,UQz,UQ2,         &
                                                     LXpt,UXpt,TwoZ,SqZ,CoFact,RL1,TmpX,TmpY, &
                                                     LXptX,LXptY,LXptZ,UXptX,UXptY,UXptZ,RL2
       INTEGER                                    :: I,J,IQ,IC,JQ,JC,KQ,KC,L,Ell,L1,L2,M,N,LMN,IGrid
!-------------------------------------------------------------------------------------------------------
       IF(PointOutSideBox(R,Node%Box))THEN
          RhoAtPoint=Zero
       ELSEIF(Node%Leaf)THEN            
          RhoAtPoint=Zero
          RQx=R(1)-Node%Qx            
          RQy=R(2)-Node%Qy
          RQz=R(3)-Node%Qz
          RQ2=RQx*RQx+RQy*RQy+RQz*RQz
          Xpt=EXP(-Node%Zeta*RQ2)
          IF(Xpt>1.D-16)THEN
             TwoZ=Two*Node%Zeta
             LambdaX(0)=Xpt
             LambdaY(0)=One
             LambdaZ(0)=One
             LambdaX(1)=TwoZ*RQx*Xpt
             Lambday(1)=TwoZ*RQy
             LambdaZ(1)=TwoZ*RQz
             DO L=2,Node%Ell
                L1=L-1
                L2=L-2
                RL1=DBLE(L1)
                LambdaX(L)=TwoZ*(RQx*LambdaX(L1)-RL1*LambdaX(L2))
                LambdaY(L)=TwoZ*(RQy*LambdaY(L1)-RL1*LambdaY(L2))
                LambdaZ(L)=TwoZ*(RQz*LambdaZ(L1)-RL1*LambdaZ(L2))
             ENDDO
             DO L=0,Node%Ell
                DO M=0,Node%Ell-L
                   DO N=0,Node%Ell-M-L
                      LMN=LMNDex(L,M,N)
                      RhoAtPoint=RhoAtPoint+LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Node%Co(LMN)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ELSE
          RhoAtPoint=RhoAtPoint(Node%Descend,R)+RhoAtPoint(Node%Descend%Travrse,R)
       ENDIF
     END FUNCTION RhoAtPoint
!
END MODULE RhoUtil
