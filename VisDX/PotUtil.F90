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
!
  USE QCTCThresholds
  USE PoleTree 
  USE MondoPoles
  USE TreeWalk
  USE PBCFarField
  USE PFFTen
!
  CONTAINS
!=================================================================================
!     
!=================================================================================
    SUBROUTINE PotCubed(Args,Del,Origin,Nx,Ny,Nz)
      TYPE(ARGMT)                     :: Args
      INTEGER                         :: I,J,K,Nx,Ny,Nz
      REAL(DOUBLE),DIMENSION(3)       :: Origin
      REAL(DOUBLE)                    :: Del
      REAL(DOUBLE)                    :: NukE,NukeCo,NukePole,PExtent
      REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
      REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
      REAL(DOUBLE)                    :: AA=One/AngstromsToAU 
      INTEGER                         :: NC
      REAL(DOUBLE),DIMENSION(3)       :: PTmp
      CHARACTER(LEN=DCL) :: Name	
!     SET THE THRESHOLDS (LOOSE)
      TauPAC=1.D-3
      TauMAC=1.D-3 
!     Get multipoles and density
      CALL Get(RhoPoles)
      CALL Get(Rho,'Rho',Args,0)
!     Initialize the auxiliary density arrays
      CALL InitRhoAux
!     Setup global arrays for computation of multipole tensors
      CALL MultipoleSetUp()
!     Build the global PoleTree representation of the total density
      CALL RhoToPoleTree
!     Set the electrostatic background 
      CALL PBCFarFieldSetUp(PoleRoot,GM)
!     Delete the auxiliary density arrays
      CALL DeleteRhoAux
!     Delete the Density
      CALL Delete(Rho)
!     WRITE POTENTIAL TO FILE
!     Cubes file to scratch directory (use PWD_O=.TRUE. to go to PWD)
      NAME=TrixFile('PotCubes',Args)
      WRITE(*,*)' POTENTIAL WRITTEN TO '//TRIM(NAME)	
      CALL OpenASCII(NAME,77)
      WRITE(77,*)' object 1 class array items ',NX*NY*NZ,' data follows '
      DO I=1,Nx    
         DO J=1,Ny
            DO K=1,Nz
               NukE=Zero
               Prim%P=(/I,J,K/)*Del+Origin
!              Test delta function is equiv to an H nucleus
               HGBra(1) =-(NuclearExpnt/Pi)**(ThreeHalves)
               SPBraC(0)=-One
               Prim%Ell=0
               Prim%Zeta=NuclearExpnt
!              Set the MAC
               DP2=(One/TauMAC)**(Two/DBLE(SPEll+2))
!              Set the PAC
               PExtent=Extent(0,NuclearExpnt,HGBra,TauPAC)
!              Initialize <KET|
               CALL SetKet(Prim,PExtent)
               PTmp=Prim%P
               DO NC=1,CS_IN%NCells
!                 Set Atomic Coordinates
                  Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
                  PBox%Center=Prim%P
!                 Walk the walk
                  CALL VWalk(PoleRoot)
               ENDDO
!              Reset the Atomic Coordinates
               Prim%P=PTmp
!              Accumulate the atomic contribution
               NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!              Add in the Far Field, Dipole and Quadripole  Correction
               IF(GM%PBC%Dimen>0) THEN
                  NukE=NukE+CTraxFF(Prim,HGBra,GM)
               ENDIF
               WRITE(77,7)NukE
7              FORMAT(F22.8) 
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
    END SUBROUTINE PotCubed
!
END MODULE 

