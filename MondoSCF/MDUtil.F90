MODULE MDUtil
!
  USE DerivedTypes
  USE GlobalScalars
  USE SCFLocals
  USE MDLocals
  USE AtomPairs
!
  IMPLICIT NONE
  CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION CrossProduct(A,B)
!
      REAL(DOUBLE),DIMENSION(3)        :: A,B
      REAL(DOUBLE),DIMENSION(3)        :: CrossProduct
!
      CrossProduct(1) = A(2)*B(3)-A(3)*B(2)
      CrossProduct(2) = A(3)*B(1)-A(1)*B(3)
      CrossProduct(3) = A(1)*B(2)-A(2)*B(1)
!
    END FUNCTION CrossProduct
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION threeBYthreeDeterminant(A)
!
!   Author: Hugh Nymeyer
!   Does:   returns determinant of 3x3 real matrix
!   Date:   10-21-02      
!
      REAL(DOUBLE),DIMENSION(3,3)     :: A
      REAL(DOUBLE)                    :: threeBYthreeDeterminant
!
      threeBYthreeDeterminant = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) - &
                                A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3)) + &
                                A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
!
    END FUNCTION threeBYthreeDeterminant
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION threeBYthreeInverse(A)
!
!   Author: Hugh Nymeyer
!   Does:   3x3 matrix inversion from hardcoded expression -- first checks for det!=0
!   Date:   10-19-02      
!
      REAL(DOUBLE),DIMENSION(3,3)                 :: A
      REAL(DOUBLE)                                :: det
      REAL(DOUBLE),DIMENSION(3,3)                 :: threeBYthreeInverse
!
      det = threeBYthreeDeterminant(A)
      if (det.EQ.0.0D0) then
         call MondoHalt(666,' Attempt to invert singular matrix!')
      end if
!
      threeBYthreeInverse(1,1) =  (A(2,2)*A(3,3)-A(2,3)*A(3,2))/det
      threeBYthreeInverse(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))/det
      threeBYthreeInverse(1,3) =  (A(1,2)*A(2,3)-A(1,3)*A(2,2))/det
      threeBYthreeInverse(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))/det
      threeBYthreeInverse(2,2) =  (A(1,1)*A(3,3)-A(1,3)*A(3,1))/det
      threeBYthreeInverse(2,3) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))/det
      threeBYthreeInverse(3,1) =  (A(2,1)*A(3,2)-A(2,2)*A(3,1))/det
      threeBYthreeInverse(3,2) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))/det
      threeBYthreeInverse(3,3) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))/det
!
    END FUNCTION threeBYthreeInverse
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeCM(GM)
!
!     Author: Hugh Nymeyer
!     Does:   Computes center of mass coordinates
!     Date:   10-18-02
!
      TYPE(CRDS)                   :: GM
      REAL(DOUBLE),DIMENSION(3)    :: computeCM
      INTEGER                      :: i
      REAL(DOUBLE)                 :: Mass,MassTotal
!
      computeCM = 0.0D0
      MassTotal=0
      Do i=1,GM%NAtms
         Mass      = GM%AtMss%D(i)
         MassTotal = MassTotal  + Mass
#ifdef PERIODIC
         computeCM = computeCM  + Mass*GM%Carts%D(:,i)
#else
         computeCM = computeCM  + Mass*GM%Carts%D(:,i)
#endif
      End Do
      computeCM = computeCM / MassTotal
!
    END FUNCTION computeCM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    SUBROUTINE removeCM(GM)
!
!     Author: Hugh Nymeyer
!     Does:   Computes center of mass coordinates, translates the center of mass to the 
!             origin, and updates the box wrapped coordinates from the absolute 
!             coordinates if necessary.
!
      TYPE(CRDS)                       :: GM
      REAL(DOUBLE),DIMENSION(3)        :: Xcm
      INTEGER                          :: i
!
      Xcm = computeCM(GM)
!
#ifdef PERIODIC
      Do i=1,GM%Natms
         GM%Carts%D(:,i) = GM%Carts%D(:,i) - Xcm
      End Do
      Call WrapAtoms(GM)
#else
      Do i=1,GM%Natms
         GM%Carts%D(:,i) = GM%Carts%D(:,i) - Xcm
      End Do
#endif
!
    END SUBROUTINE removeCM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeVCM(GM)
!
!     Author: Hugh Nymeyer
!     Does:   Computes the center of mass velocity
!
      TYPE(CRDS)                :: GM
      REAL(DOUBLE),DIMENSION(3) :: computeVCM
      INTEGER                   :: i
      REAL(DOUBLE)              :: Mass,MassTotal
!
      computeVCM=0.0D0
      MassTotal=0
      do i=1,GM%NAtms
         Mass       = GM%AtMss%D(i)
         MassTotal  = MassTotal  + Mass
         computeVCM = computeVCM + Mass*GM%Vects%D(:,i)
      end do
      computeVCM = computeVCM / MassTotal
!
    END FUNCTION COMPUTEVCM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    SUBROUTINE removeVCM(GM)
!
!     Author: Hugh Nymeyer
!     Does:   Computes center of mass velocity and removes it.
!     Date:   10-19-02
!
      TYPE(CRDS)                          :: GM
      REAL(DOUBLE),DIMENSION(3)           :: VCM
      INTEGER                             :: i
!
      VCM = computeVCM(GM)
      do i=1,GM%NAtms
         GM%Vects%D(:,i) = GM%Vects%D(:,i) - VCM(:)
      end do
!       
    END SUBROUTINE removeVCM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeCM_KIN(GM)
!
      TYPE(CRDS)                          :: GM
      REAL(DOUBLE),DIMENSION(3)           :: VCM
      REAL(DOUBLE)                        :: DKE
      REAL(DOUBLE)                        :: TotalMass
      REAL(DOUBLE)                        :: computeCM_KIN
      INTEGER                             :: i
!
      TotalMass=0.0D0
      VCM = computeVCM(GM)
      do i=1,GM%NAtms
         TotalMass = TotalMass + GM%AtMss%D(i)
      end do
      DKE = 0.5D0 * TotalMass * (VCM(1)**2 + VCM(2)**2 + VCM(3)**2)
!
    END FUNCTION computeCM_KIN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeAngMoment(GM)
!
!     Author: Hugh Nymeyer
!     Does:   computes angular momentum pseudo-vector from absolute coordinates
!     Date:   10-19-02 
!
      TYPE(CRDS)                              :: GM
      REAL(DOUBLE),DIMENSION(3)               :: computeAngMoment
      INTEGER                                 :: i
      REAL(DOUBLE),DIMENSION(3)               :: Xcm
!      
      Xcm = computeCM(GM)
!
      computeAngMoment = 0.0D0
      do i=1,GM%NAtms
         computeAngMoment = computeAngMoment + GM%AtMss%D(i) * &
              CrossProduct(GM%Carts%D(:,i)-Xcm,GM%Vects%D(:,i))
      end do
!
    END FUNCTION computeAngMoment
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeInertiaTensor(GM)
!
!     Author: Hugh Nymeyer
!     Does:   computes the inertia tensor from absolute coordinates
!     Date:   10-19-02
!
      TYPE(CRDS)                              :: GM
      REAL(DOUBLE),DIMENSION(3,3)             :: computeInertiaTensor
      REAL(DOUBLE)                            :: Mass,X,Y,Z,XX,XY,XZ,YY,YZ,ZZ
      REAL(DOUBLE),DIMENSION(3)               :: Xcm
      INTEGER                                 :: i
!
      XX=0.0D0
      XY=0.0D0
      XZ=0.0D0
      YY=0.0D0
      YZ=0.0D0
      ZZ=0.0D0
!
      Xcm = computeCM(GM)
!
      do i=1,GM%NAtms
         Mass = GM%AtMss%D(i)
         X    = GM%Carts%D(i,1)-Xcm(1)
         Y    = GM%Carts%D(i,2)-Xcm(2)
         Z    = GM%Carts%D(i,3)-Xcm(3)
         XX = Mass * X * X
         XY = Mass * X * Y
         XZ = Mass * X * Z
         YY = Mass * Y * Y
         YZ = Mass * Y * Z
         ZZ = Mass * Z * Z
      end do
!
      computeInertiaTensor(1,1) = YY + ZZ
      computeInertiaTensor(2,1) = -XY
      computeInertiaTensor(3,1) = -XZ
      computeInertiaTensor(1,2) = -XY
      computeInertiaTensor(2,2) = XX + ZZ
      computeInertiaTensor(3,2) = -YZ
      computeInertiaTensor(1,3) = -XZ
      computeInertiaTensor(2,3) = -YZ
      computeInertiaTensor(3,3) = XX + YY
!
    END FUNCTION computeInertiaTensor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    SUBROUTINE removeAngMoment(GM)
!
!     Author: Hugh Nymeyer
!     Does:   removes angular momentum and puts the amount of kinetic energy in the 
!             angular motion into DKE.  
!     Date:   10-19-02
!
      TYPE(CRDS)                                 :: GM
      REAL(DOUBLE),DIMENSION(3)                  :: omega           ! body rotation vector
      REAL(DOUBLE),DIMENSION(3)                  :: AngMoment       ! angular momentum
      REAL(DOUBLE),DIMENSION(3,3)                :: IT,ITinverse    ! inertia tensor
      INTEGER                                    :: i,j
!
      AngMoment = computeAngMoment(GM)
      IT        = computeInertiaTensor(GM)
      ITinverse = threeBYthreeInverse(IT)
      omega     = MATMUL(ITinverse,AngMoment)
      do i=1,GM%NAtms
         GM%Vects%D(:,i) = GM%Vects%D(:,i) - CrossProduct(omega,GM%Carts%D(:,i))
      end do
!
    END SUBROUTINE removeAngMoment
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeKE(GM)
!
      TYPE(CRDS)                             :: GM
      REAL(DOUBLE)                           :: computeKE
      INTEGER                                :: i
!
      computeKE=0.0D0
      do i=1,GM%NAtms
         computeKE = computeKE + GM%AtMss%D(i) * &
              (GM%Vects%D(1,i)**2 + GM%Vects%D(2,i)**2 + GM%Vects%D(3,i)**2)
      end do
      computeKE = computeKE * 0.5D0
!
    END FUNCTION computeKE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    FUNCTION computeTinstant(Ctrl,GM)
!
!   Author: Hugh Nymeyer
!   Does:   Returns the instantaneous temperature computed from the kinetic energy.
!           Note: Ctrl is passed here, b/c I am assuming that it contains information
!           which might affect the number of degrees of freedom, e.g., constraints.
!
      TYPE(SCFControls)                     :: Ctrl
      TYPE(CRDS)                            :: GM
      REAL(DOUBLE)                          :: computeTinstant
      INTEGER                               :: DegOfFreedom
!
      DegOfFreedom = computeDegOfFreedom(Ctrl,GM)
      computeTinstant = 0.2D1 * computeKE(GM) / DegOfFreedom * HartreesToKelvin
!
    END FUNCTION computeTinstant
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      FUNCTION computeDegOfFreedom(Ctrl,GM)
!
!     Author: Hugh Nymeyer
!     Does:   What it says --- this routine still needs some work to check for other
!             constraints.
!
	IMPLICIT NONE
!
        REAL(DOUBLE)               :: computeDegOfFreedom
        TYPE(SCFControls)          :: Ctrl
	TYPE(CRDS)                 :: GM
!
!       first guess
!
        computeDegOfFreedom = 3 * GM%NAtms
!
!       now subtract for constraints
!
        if (Ctrl%MDC%REM_TRANS) &
             computeDegOfFreedom = computeDegOfFreedom - 3
        if (Ctrl%MDC%REM_ROTAT) &
             computeDegOfFreedom = computeDegOfFreedom - 3
!
      END FUNCTION computeDegOfFreedom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SUBROUTINE AssignBoltzmannVelocity(GM,TEMP0)
!
! Author: Hugh Nymeyer
! Date: 10-29-02
! Does: Assigns velocity components from a Boltzmann distribution at TEMP0
!
        TYPE(CRDS)                       :: GM
        REAL(DOUBLE)                     :: TEMP0
        INTEGER                          :: i,j
!
!   Assign velocity from Gaussian with zero mean and variance kT/mass
!
        do i=1,GM%NAtms
           do j=1,3
              GM%Vects%D(j,i) = Gaussian(0.0D0,TEMP0*KelvinToHartrees/GM%AtMss%D(i))
           end do
        end do
!
      END SUBROUTINE AssignBoltzmannVelocity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      FUNCTION Gaussian(mean,var)
!
! Author: Hugh Nymeyer
! Date: 10-29-02
! Does: Polar form of the Box-Muller algorithm for generating gaussian deviates
!
        REAL(DOUBLE)              :: Gaussian
        REAL(DOUBLE)              :: mean,var
        REAL(DOUBLE),SAVE         :: x1,x2,y1,y2,w
        LOGICAL,SAVE              :: even=.FALSE.
!
        if (.NOT. even) then
           do
              call random_number(x1)
              call random_number(x2)
              x1 = 2.0D0 * x1 - 1.0D0;
              x2 = 2.0D0 * x2 - 1.0D0;
              w = x1 * x1 + x2 * x2;
              if (w .LT. 1.0D0) exit
           end do
!
           w = sqrt( -2.0D0 * log( w ) / w );
           y1 = x1 * w;
           y2 = x2 * w;
!
           y1 = y1 * sqrt(var) + mean
           y2 = y2 * sqrt(var) + mean
!
           even = .TRUE.
           Gaussian=y1
!
        else 
           even = .FALSE.
           Gaussian=y2
        end if
!
      END FUNCTION Gaussian
!--------------------------------------------------------
!
!
  END MODULE MDUtil






