MODULE BaderRhoUtil
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

      SUBROUTINE GetArgs(nx,ny,nz,nrho,vasp,BaderTol)
        INTEGER :: nx,ny,nz
        LOGICAL :: vasp
        REAL(DOUBLE) :: BaderTol

        nx=128
        ny=128
        nz=128
        nrho=nx*ny*nz
        vasp=.false.
        BaderTol=1e-4

!        CALL OpenASCII(InpFile,Inp)         
!        IF(OptDblQ(Inp,'TauPAC',TauPAC))THEN
!          Mssg=TRIM(ProcessName('QCTC'))//' TauPAC = '//TRIM(DblToShrtChar(TauPAC))
!
!          CALL OpenASCII(OutFile,Out)         
!          WRITE(Out,*)TRIM(Mssg)
!          CLOSE(Out)
!        ENDIF

      END SUBROUTINE GetArgs


      SUBROUTINE MakeCubeRho(Arg,Origin,Lattice,Ndim,Nx,Ny,Nz,BaderRho,Rcar,Rdir,BaderProj)
        INTEGER                    :: I,J,K,Nx,Ny,Nz,Ndim,tenths_done
        REAL(DOUBLE),DIMENSION(3)  :: Origin,Del,Half
        REAL(DOUBLE),DIMENSION(3,3):: Lattice
        INTEGER,POINTER,DIMENSION(:,:,:)       :: BaderProj
        REAL(DOUBLE),POINTER,DIMENSION(:,:,:)  :: BaderRho
        REAL(DOUBLE),POINTER,DIMENSION(:,:)    :: Rcar,Rdir
        REAL(DOUBLE)               :: Den,Zero,Vol,RhoSum
        TYPE(ARGMT)                :: Arg
#ifdef PERIODIC
        INTEGER                    :: NC
#endif
        ALLOCATE(BaderRho(Nx,Ny,Nz))
        ALLOCATE(BaderProj(Nx,Ny,Nz))
!---------------------------------------------------------------------------
!     Get the geometry
        CALL Get(GM,Tag_O=CurGeom)
!       CALL PPrint(GM,TrixFile('xyz',Arg,PWD_O=.TRUE.),Geo,'XYZ')
#ifdef PERIODIC
!       Calculate the Number of Cells
        CALL SetCellNumber(GM)
#endif
!      SET DENSITY THRESHOLDS (LOOSE)
!       Del=3.D-1
        TauRho=1.D-6
!      BUILD RHO-TREE
        CALL RhoToTree(Arg)
        CALL SetBBox(RhoRoot%Box,Box)
#ifdef PERIODIC
        CALL MakeBoxPeriodic(Box)
#endif
!       DETERMINE CUBE GRID (SMART, UNSTRUCTURED GRIDS DONT WORK WITH DX) 
        Origin=Box%BndBox(:,1)
        Lattice=0
        Lattice(1,1)=(Box%BndBox(1,2)-Box%BndBox(1,1))*AA
        Lattice(2,2)=(Box%BndBox(2,2)-Box%BndBox(2,1))*AA
        Lattice(3,3)=(Box%BndBox(3,2)-Box%BndBox(3,1))*AA
        Del(1)=(Box%BndBox(1,2)-Box%BndBox(1,1))/Nx
        Del(2)=(Box%BndBox(2,2)-Box%BndBox(2,1))/Ny
        Del(3)=(Box%BndBox(3,2)-Box%BndBox(3,1))/Nz

! Get Cartesian coordinates
        Ndim=GM%NAtms
        ALLOCATE(Rdir(Ndim,3),Rcar(Ndim,3))
        DO I=1,Ndim
          Rcar(I,1:3)=GM%Carts%D(1:3,I)*AA
!          write(*,*) I,Rcar(I,:)
!          Rdir(I,1)=Rcar(I,1)/Lattice(1,1)+0.5
!          Rdir(I,2)=Rcar(I,2)/Lattice(2,2)+0.5
!          Rdir(I,3)=Rcar(I,3)/Lattice(3,3)+0.5
          Rdir(I,1)=(Rcar(I,1)-Origin(1)*AA)/Lattice(1,1)
          Rdir(I,2)=(Rcar(I,2)-Origin(2)*AA)/Lattice(2,2)
          Rdir(I,3)=(Rcar(I,3)-Origin(3)*AA)/Lattice(3,3)
!          write(*,*) I,Rdir(I,:)
        ENDDO

! Building the density grid
!       write(*,*)'Building a density grid'
        write(*,'(1A26)') 'CALCULATING DENSITY GRID'
        Zero=0
        Half=(/0.5,0.5,0.5/)
        RhoSum=0
        Vol=Del(1)*Del(2)*Del(3)
        DO I=1,Nx
           IF ((i*10/nx) > tenths_done) THEN
              tenths_done=(i*10/nx)
              WRITE(*,'(1X,1I4,1A6)') (tenths_done*10),'% done'
           END IF
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
!                  R=((/I,J,K/)-Half)*Del+Origin
                  R=(/I,J,K/)*Del+Origin
                  Den=RhoAtPoint(RhoRoot)
#endif

                  BaderRho(i,j,k)=Den*Vol
                  RhoSum=RhoSum+Den*Vol
              ENDDO
           ENDDO
        ENDDO
!        write(*,*)'Done building the density.'
!        WRITE(*,'(1A23,11X,1F12.5,/)') 'NUMBER OF ELECTRONS: ',RhoSum
!       Delete the RhoTree
        CALL DeleteRhoTree(RhoRoot)
        CALL Delete(GM)
      END SUBROUTINE MakeCubeRho
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
PROGRAM Bader
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE Globals  
  USE BaderRhoUtil
  USE ChargeUtil
  USE MatrixUtil

  IMPLICIT NONE
  CHARACTER(LEN=5),PARAMETER  :: Prog='Bader'
  TYPE(ARGMT)                 :: Arg
  INTEGER,POINTER,DIMENSION(:,:,:) :: BaderProj
  REAL(DOUBLE),POINTER,DIMENSION(:,:,:) :: BaderRho
  REAL(DOUBLE),POINTER,DIMENSION(:,:) :: BaderCharge,WSCharge,Rcar,Rdir
  REAL(DOUBLE),POINTER,DIMENSION(:) :: BaderDist,BaderAChg
  INTEGER,POINTER,DIMENSION(:) :: BaderAtom
  REAL(DOUBLE),DIMENSION(3,3) :: Lattice
  REAL(DOUBLE),DIMENSION(3)   :: Origin
  REAL(DOUBLE) :: BaderTol
  LOGICAL :: vasp
  INTEGER :: Ndim,nx,ny,nz,bdim,nrho,wdim

!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
  CALL GetArgs(nx,ny,nz,nrho,vasp,BaderTol)
  CALL MakeCubeRho(Args,Origin,Lattice,Ndim,nx,ny,nz,BaderRho,Rcar,Rdir,BaderProj)
  CALL bader(BaderCharge,bdim,BaderRho,nx,ny,nz,nrho,Rdir,ndim,lattice,BaderAchg,&
             BaderAtom,BaderDist,BaderProj)
  CALL wigner_seitz(WSCharge,wdim,BaderRho,nx,ny,nz,Rcar,Rdir,ndim,nrho,vasp)
  CALL output(BaderCharge,bdim,WSCharge,wdim,ndim,lattice,nx,ny,nz,         &
  &           BaderAchg,BaderDist,BaderAtom,BaderTol)
  CALL write_max_rho(BaderProj,nx,ny,nz)
  CALL write_atm_rho(BaderProj,nx,ny,nz,lattice,Rdir,Ndim,BaderCharge,Bdim,BaderTol)

  CALL Delete(Args)

  CALL ShutDown(Prog)
END PROGRAM Bader
