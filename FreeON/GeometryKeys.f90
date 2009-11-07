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
MODULE GeometryKeys
  ! Coordinate input delemeters
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY_BEGIN    ='<BeginGeometry>'
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY_NEXT     ='<NextGeometry>'
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY_END      ='<EndGeometry>'

  CHARACTER(LEN=*), PARAMETER :: PRODUCTS_BEGIN    ='<BeginProducts>'
  CHARACTER(LEN=*), PARAMETER :: PRODUCTS_END      ='<EndProducts>'

  CHARACTER(LEN=*), PARAMETER :: REACTANTS_BEGIN   ='<BeginReactants>'
  CHARACTER(LEN=*), PARAMETER :: REACTANTS_END     ='<EndReactants>'

  CHARACTER(LEN=*), PARAMETER :: CLONE_BEGIN       = '<BeginClone'
  CHARACTER(LEN=*), PARAMETER :: CLONE_END         = '<EndClone'

  ! Geometry
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY          ='Geometry'
  CHARACTER(LEN=*), PARAMETER :: IN_AU             ='InAU'
  CHARACTER(LEN=*), PARAMETER :: H_ORDER           ='HOrder'
  CHARACTER(LEN=*), PARAMETER :: Z_ORDER           ='ZOrder'
  CHARACTER(LEN=*), PARAMETER :: RANDOM_ORDER      ='RandomOrder'
  CHARACTER(LEN=*), PARAMETER :: TRAVEL_ORDER      ='TravelOrder'
  CHARACTER(LEN=*), PARAMETER :: TABLETRAV_ORDER   ='TableTravOrder'
  CHARACTER(LEN=*), PARAMETER :: NO_ORDER          ='NoOrder'

  ! Clones
  CHARACTER(LEN=*), PARAMETER :: CLONES            ='Clones'

  ! Charge
  CHARACTER(LEN=*), PARAMETER :: TOTAL_CHARGE      ='Charge'

  ! Multiplicity
  CHARACTER(LEN=*), PARAMETER :: MULTIPLICITY      ='Multiplicity'
END MODULE GeometryKeys
