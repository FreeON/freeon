!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE PoleNodeType
  USE Derivedtypes
  USE GlobalScalars
  IMPLICIT NONE
  !=================================================================================================
  !  Hierarchical density node
  !=================================================================================================
  TYPE PoleNode
    LOGICAL                               :: Leaf     ! Is this a data containing node?
    INTEGER                               :: Bdex     ! Begining index of ORB list for this node
    INTEGER                               :: Edex     ! ENDign index of ORB list for this node
    INTEGER                               :: NQ       ! Number of centers
    INTEGER                               :: Ell      ! Ell type
    INTEGER                               :: EllCD    ! Maximium Ell of the Distributuions in the box
    REAL(DOUBLE)                          :: Zeta     ! Minimum exponent in this node
    REAL(DOUBLE)                          :: Strength ! Strength of the Pole
    REAL(DOUBLE)                          :: DMax2    ! (Max distance)^2 from node center to dist
    REAL(DOUBLE)                          :: WCoef    ! Weight for the Gaussian in the Box
    TYPE(BBox)                            :: Box      ! Bounding Box of distribution (for PAC)
    TYPE(PoleNode),POINTER                :: Descend  ! Next node in tree descent
    TYPE(PoleNode),POINTER                :: Travrse  ! Next node in tree traversal
#ifdef POINTERS_IN_DERIVED_TYPES
    REAL(DOUBLE),DIMENSION(:),POINTER     :: S        ! Im component of the multipole tensor
    REAL(DOUBLE),DIMENSION(:),POINTER     :: C        ! Re component of the multipole tensor
    REAL(DOUBLE),DIMENSION(:),POINTER     :: Co       ! Coefficients of the HGTF density
#else
    REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: S        ! Im component of the multipole tensor
    REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: C        ! Re component of the multipole tensor
    REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: Co       ! Coefficients of the HGTF density
#endif
  END TYPE PoleNode
END MODULE PoleNodeType
