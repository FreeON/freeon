MODULE ThresholdsMod
  USE DerivedTypes
  IMPLICIT NONE
  !
  ! Thresholds for linear scaling routines:
  !
  ! Sparse blocked matrix thresholds
  REAL(DOUBLE),DIMENSION(4) :: TrixNeglect=(/1.D-5, 1.D-6, 1.D-7,  1.D-8/)
!  REAL(DOUBLE),DIMENSION(4) :: TrixNeglect=(/1.D-4, 1.D-5, 1.D-6,  1.D-7 /)
  ! HiCu threshold
  REAL(DOUBLE),DIMENSION(4) :: CubeNeglect=(/1.D-3, 1.D-5, 1.D-7,  1.D-9 /)
  ! QCTC and ONX threshold
  REAL(DOUBLE),DIMENSION(4) :: TwoENeglect=(/1.D-6, 1.D-8, 1.D-10, 1.D-12/)
  ! Distribution threshold
  REAL(DOUBLE),DIMENSION(4) :: DistNeglect=(/1.D-8, 1.D-10,1.D-12, 1.D-14/)
  !
  ! Convergence criteria for SCF and force routines:
  !
  ! Max error in the relative total energy
  REAL(DOUBLE),DIMENSION(4) :: ETol       =(/ 1.D-5, 1.D-7, 1.D-9,  1.D-11 /)
  ! Max element of the density matrix
  REAL(DOUBLE),DIMENSION(4) :: DTol       =(/ 1.D-2, 1.D-3, 1.D-4,  1.D-5  /)
  ! Max Cartesian gradient
  REAL(DOUBLE),DIMENSION(4) :: GTol       =(/ 5.D-2, 5.D-3, 5.D-4,  5.D-5  /)
  ! Max Cartesian displacement
  REAL(DOUBLE),DIMENSION(4) :: XTol       =(/ 1.D-1, 1.D-2, 1.D-3,  1.D-4  /)
END MODULE ThresholdsMod
