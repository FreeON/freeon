MODULE RhoUtil
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE AtomPairs
   USE Macros
!
   USE RhoTree
   USE HiCuThresholds
   USE CubeTree
!  Global variables
   INTEGER                         :: NP
   REAL(DOUBLE),DIMENSION(3)       :: R
   REAL(DOUBLE)                    :: AA=One/AngstromsToAU
   CONTAINS

      SUBROUTINE RhoCubed(Arg,Del,Origin,Nx,Ny,Nz)
        INTEGER                    :: I,J,K,Nx,Ny,Nz
        REAL(DOUBLE),DIMENSION(3)  :: Origin
        REAL(DOUBLE)               :: Del,Den
        TYPE(ARGMT)                :: Arg
        TYPE(CRDS)                 :: GMLoc
#ifdef PERIODIC
        INTEGER                    :: NC
#endif
!---------------------------------------------------------------------------
!     Get the geometry
#ifdef MMech
       IF(HasMM()) THEN
       CALL Get(GMLoc,Tag_O=CurGeom)
       ELSE
#endif
       CALL Get(GMLoc,Tag_O=CurGeom)
#ifdef MMech
       ENDIF
#endif
       CALL PPrint(GMLoc,TrixFile('xyz',Arg,PWD_O=.TRUE.),Geo,'XYZ')
#ifdef PERIODIC
!      Calculate the Number of Cells
       CALL SetCellNumber(GMLoc)
#endif
!      SET DENSITY THRESHOLDS (LOOSE)
       Del=3.D-1
       TauRho=1.D-3
!      BUILD RHO-TREE
       CALL RhoToTree(Arg)
       CALL SetBBox(RhoRoot%Box,Box)
#ifdef PERIODIC
       CALL MakeBoxPeriodic(Box)
#endif
!       DETERMINE CUBE GRID (SMART, UNSTRUCTURED GRIDS DONT WORK WITH DX) 
        Origin=Box%BndBox(:,1)
        Nx=CEILING((Box%BndBox(1,2)-Box%BndBox(1,1))/Del)
        Ny=CEILING((Box%BndBox(2,2)-Box%BndBox(2,1))/Del)
        Nz=CEILING((Box%BndBox(3,2)-Box%BndBox(3,1))/Del)
!       WRITE DENSITY TO FILE
        CALL OpenASCII(TrixFile('RhoCubes',Arg),66)
        WRITE(66,*)' object 1 class array items ',NX*NY*NZ,' data follows '
        DO I=1,Nx    
           DO J=1,Ny
              DO K=1,Nz
#ifdef PERIODIC
                  Den=Zero
                  DO NC=1,CS_OUT%NCells
!                    Set Atomic Coordinates
                     R=(/I,J,K/)*Del+Origin+CS_OUT%CellCarts%D(:,NC)
                     Den=Den+RhoAtPoint(RhoRoot)
                  ENDDO
#else
                  R=(/I,J,K/)*Del+Origin
                  Den=RhoAtPoint(RhoRoot)
#endif
                 WRITE(66,7) Den
               7 FORMAT(8(2x,F12.8)) 
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
!       Delete the RhoTree
        CALL DeleteRhoTree(RhoRoot)
        CALL Delete(GMLoc)
      END SUBROUTINE RhoCubed
!=================================================================================
!     Sums the significant contributions (leaves) to the density at a point
!=================================================================================
      RECURSIVE FUNCTION RhoAtPoint(Node)
         USE RhoTree
         USE BoundingBox
         REAL(DOUBLE)                               :: RhoAtPoint
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
            RhoAtPoint=RhoAtPoint(Node%Descend)+RhoAtPoint(Node%Descend%Travrse)
         ENDIF 
      END FUNCTION RhoAtPoint
END MODULE 
!=================================================================================
!     
!=================================================================================
MODULE PotUtil
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE AtomPairs
   USE Macros
!  Global variables
   REAL(DOUBLE) :: AA=One/AngstromsToAU
   CONTAINS
      SUBROUTINE PotCubed(Arg,Del,Origin,Nx,Ny,Nz)
         USE QCTCThresholds
         USE PoleTree 
         USE MondoPoles
         USE TreeWalk
         USE Globals
#ifdef PERIODIC
         USE PBCFarField
         USE PFFT
#endif
         TYPE(ARGMT)                     :: Arg
         INTEGER                         :: I,J,K,Nx,Ny,Nz
         REAL(DOUBLE),DIMENSION(3)       :: Origin
         REAL(DOUBLE)                    :: Del
         REAL(DOUBLE)                    :: NukE,NukeCo,NukePole,PExtent
         REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
         REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
         TYPE(CRDS)                      :: GMLoc
#ifdef PERIODIC
         INTEGER                         :: NC
         REAL(DOUBLE),DIMENSION(3)       :: PTmp
#endif
!--------------------------------------------------------------
!        Get the geometry
#ifdef MMech
         IF(HasMM()) THEN
           CALL Get(GMLoc,Tag_O='GM_MM'//CurGeom)
         ELSE
#endif
           CALL Get(GMLoc,Tag_O=CurGeom)
#ifdef MMech
         ENDIF
#endif
!        SET THE THRESHOLDS (LOOSE)
         TauPAC=1.D-2
         TauMAC=1.D-2 
!        Get multipoles and density
         CALL Get(RhoPoles,SCFCycl)
         CALL Get(Rho,'Rho',Arg,0)
!        Initialize the auxiliary density arrays
         CALL InitRhoAux
!        Setup global arrays for computation of multipole tensors
         CALL MultipoleSetUp(FFEll2)
!        Build the global PoleTree representation of the total density
         CALL RhoToPoleTree
#ifdef PERIODIC
!        Calculate the Number of Cells
         CALL SetCellNumber(GMLoc)
!        Set the electrostatic background 
         CALL PBCFarFieldSetUp(PoleRoot,GMLoc)
#endif
!        Delete the auxiliary density arrays
         CALL DeleteRhoAux
!        Delete the Density
         CALL Delete(Rho)
!        WRITE POTENTIAL TO FILE
!        Cubes file to scratch directory (use PWD_O=.TRUE. to go to PWD)
         CALL OpenASCII(TrixFile('PotCubes',Arg),77)
         WRITE(77,*)' object 1 class array items ',NX*NY*NZ,' data follows '
         DO I=1,Nx    
            DO J=1,Ny
               DO K=1,Nz
                  NukE=Zero
                  Prim%P=(/I,J,K/)*Del+Origin
!                 Test delta function is equiv to an H nucleus
                  HGBra(1) =-(NuclearExpnt/Pi)**(ThreeHalves)
                  SPBraC(0)=-One
                  Prim%Ell=0
                  Prim%Zeta=NuclearExpnt
!                 Set the MAC
                  DP2=(One/TauMAC)**(Two/DBLE(SPEll+2))
!                 Set the PAC
                  PExtent=Extent(0,NuclearExpnt,HGBra,TauPAC)
!                 Initialize <KET|
                  CALL SetKet(Prim,PExtent)
#ifdef PERIODIC
                  PTmp=Prim%P
                  DO NC=1,CS_IN%NCells
!                    Set Atomic Coordinates
                     Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
                     PBox%Center=Prim%P
!                    Walk the walk
                     CALL VWalk(PoleRoot)
                  ENDDO
!                 Reset the Atomic Coordinates
                  Prim%P=PTmp
!                 Accumulate the atomic contribution
                  NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!                 Add in the Far Field, Dipole and Quadripole  Correction
                  IF(GMLoc%PBC%Dimen>0) THEN
                     NukE=NukE+CTraxFF(Prim,HGBra,GMLoc)
                  ENDIF
#else
!                 Walk the walk
                  CALL VWalk(PoleRoot)
!                 Accumulate the atomic contribution
                  NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0) 
#endif
                  WRITE(77,7)NukE
                7 FORMAT(F22.8) 
              ENDDO
            ENDDO
         ENDDO
         WRITE(77,*)' '
         WRITE(77,*)' attribute "dep" string "positions" '
         WRITE(77,*)'object 2 class gridpositions counts',NX,NY,NZ
         WRITE(77,*)' origin ',origin*AA
         WRITE(77,*)' delta ',del*AA,' 0.0  0.0 '
         WRITE(77,*)' delta 0.0 ',del*AA,' 0.0 '
         WRITE(77,*)' delta 0.0  0.0 ',del*AA
         WRITE(77,*)'object 3 class gridconnections counts ',NX,NY,NZ
         WRITE(77,*)' attribute "element type" string "cubes" '
         WRITE(77,*)' attribute "ref" string "positions" '
         WRITE(77,*)'object "electron density" class field '
         WRITE(77,*)' component "data" 1 '
         WRITE(77,*)' component "positions" 2 '
         WRITE(77,*)' component "connections" 3 '
         CLOSE(UNIT=77)
         CALL Delete(GMLoc)
      END SUBROUTINE PotCubed     
END MODULE 
!=================================================================================
!  
!=================================================================================
PROGRAM VisDX
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE Globals  
  USE RhoUtil
  USE PotUtil 
  REAL(DOUBLE)                :: Del1
  REAL(DOUBLE),DIMENSION(3)   :: Origin1
  INTEGER :: Nx1,Ny1,Nz1
  CHARACTER(LEN=5),PARAMETER  :: Prog='VisDX'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
  CALL RhoCubed(Args,Del1,Origin1,Nx1,Ny1,Nz1)
  CALL PotCubed(Args,Del1,Origin1,Nx1,Ny1,Nz1)
  CALL Delete(Args)
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM VisDX
