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


!============================================================================
!                 ** THE MYSTERIES OF MONDO MASSAGE **
!
! WHY DO WE NEED TO TRANSLATE TO GET THE RIGHT NUMBERS??????
! SHOULDNT WRAPPING TAKE CARE OF IT ALL ?????
! ALSO, WHY DOES WRAPPING APPARENTLY NOT BRING THE ATOMS BACK INTO THE UC????
! AND WHY IS ALL THE PBC STUFF LIKE WRAPPING AND TRANSLATION IN ATOM PAIRS INSTEAD
! OF PBC?  AND DO WE REALLY NEED FOUR COORDINATE ARRAYS IN CRDS?
!============================================================================
MODULE Massage
  USE Parse
  USE InOut
  USE PBC
  USE AtomPairs
  USE OptionKeys
  USE DynamicsKeys
  USE GeometryKeys
  USE ControlStructures
  IMPLICIT NONE
CONTAINS
  !============================================================================
  ! ALL REORDERING, RESCALING, WRAPPING AND TRANSLATING OF COORDINATES OCCURS
  ! HERE AND NO WHERE ELSE!
  !============================================================================
  SUBROUTINE MassageCoordinates(O,G,P)
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
    INTEGER          :: I,J,GBeg,GEnd
    !-------------------------------------------------------------------------!
    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
      GBeg=0
      GEnd=G%Clones+1
    ELSE
      GBeg=1
      GEnd=G%Clones
    ENDIF

    IF(G%Clone(GBeg)%InAU) THEN
      CALL MondoLog(DEBUG_NONE, "FreeON", "Input coordinates read in atomic units","MassageCoordinates")
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", "Input coordinates read in Angstroms","MassageCoordinates")
    ENDIF

    DO I=GBeg,GEnd
      CALL ToAtomicUnits(G%Clone(I))
      CALL PeriodicXLate(G%Clone(I))
      CALL SuperCellMe(G%Clone(I))
      CALL PeriodicXLate(G%Clone(I))
    ENDDO
  END SUBROUTINE MassageCoordinates
  !=========================================================================
  !
  !=========================================================================
  SUBROUTINE SuperCellMe(G)
    TYPE(CRDS)                :: G
    REAL(DOUBLE)              :: RSCX,RSCY,RSCZ
    INTEGER                   :: I,J,K,ISCX,ISCY,ISCZ,AT,NC
    REAL(DOUBLE),DIMENSION(3) :: RVec,IVec
    !
    TYPE(DBL_VECT)            :: AtNum,AtMss
    TYPE(DBL_RNK2)            :: Velocity,Carts,Displ,Gradients,BoxCarts,Fext
    TYPE(INT_VECT)            :: AtTyp,CConstrain
    TYPE(CHR10_VECT)          :: AtNam
    !
    IF(G%PBC%SuperCell%I(1)<=1.AND.G%PBC%SuperCell%I(2)<=1.AND.G%PBC%SuperCell%I(3)<=1)RETURN
    !
    CALL New(AtNum,G%NAtms)
    CALL New(AtTyp,G%NAtms)
    CALL New(AtNam,G%NAtms)
    CALL New(AtMss,G%NAtms)
    CALL New(CConstrain,G%NAtms)
    CALL New(Velocity,(/3,G%NAtms/))
    CALL New(Fext,(/3,G%NAtms/))
    CALL New(Carts,(/3,G%NAtms/))
    CALL New(Displ,(/3,G%NAtms/))
    CALL New(Gradients,(/3,G%NAtms/))
    CALL New(BoxCarts,(/3,G%NAtms/))
    !
    AtNum%D=G%AtNum%D
    AtTyp%I=G%AtTyp%I
    AtNam%C=G%AtNam%C
    AtMss%D=G%AtMss%D
    CConstrain%I=G%CConstrain%I
    Velocity%D=G%Velocity%D
    Fext%D=G%Fext%D
    Carts%D=G%Carts%D
    Gradients%D=G%Gradients%D
    BoxCarts%D=G%BoxCarts%D
    !
    ISCX=G%PBC%SuperCell%I(1)
    ISCY=G%PBC%SuperCell%I(2)
    ISCZ=G%PBC%SuperCell%I(3)
    RSCX=DBLE(ISCX)
    RSCY=DBLE(ISCY)
    RSCZ=DBLE(ISCZ)
    NAtoms=G%NAtms*ISCX*ISCY*ISCZ
    !
    CALL Delete(G%AtNum)
    CALL Delete(G%AtTyp)
    CALL Delete(G%AtNam)
    CALL Delete(G%AtMss)
    CALL Delete(G%CConstrain)
    CALL Delete(G%Velocity)
    CALL Delete(G%Fext)
    CALL Delete(G%Carts)
    CALL Delete(G%Displ)
    CALL Delete(G%Gradients)
    CALL Delete(G%BoxCarts)
    !
    CALL New(G%AtNum,NAtoms)
    CALL New(G%AtTyp,NAtoms)
    CALL New(G%AtNam,NAtoms)
    CALL New(G%AtMss,NAtoms)
    CALL New(G%CConstrain,NAtoms)
    CALL New(G%Velocity,(/3,NAtoms/))
    CALL New(G%Fext,(/3,NAtoms/))
    CALL New(G%Carts,(/3,NAtoms/))
    CALL New(G%Displ,(/3,NAtoms/))
    CALL New(G%Gradients,(/3,NAtoms/))
    CALL New(G%BoxCarts,(/3,NAtoms/))
    !
    NC=0
    DO I=1,ISCX
      DO J=1,ISCY
        DO K=1,ISCZ
          IVec=(/DBLE(I-1),DBLE(J-1),DBLE(K-1)/)
          RVec=(/G%PBC%BoxShape%D(1,1)*DBLE(I-1)+G%PBC%BoxShape%D(1,2)*DBLE(J-1)+G%PBC%BoxShape%D(1,3)*DBLE(K-1), &
                 G%PBC%BoxShape%D(2,1)*DBLE(I-1)+G%PBC%BoxShape%D(2,2)*DBLE(J-1)+G%PBC%BoxShape%D(2,3)*DBLE(K-1), &
                 G%PBC%BoxShape%D(3,1)*DBLE(I-1)+G%PBC%BoxShape%D(3,2)*DBLE(J-1)+G%PBC%BoxShape%D(3,3)*DBLE(K-1)/)
          DO AT=1,G%NAtms
            NC=NC+1
            G%AtNum%D(NC)=AtNum%D(AT)
            G%AtTyp%I(NC)=AtTyp%I(AT)
            G%AtNam%C(NC)=AtNam%C(AT)
            G%AtMss%D(NC)=AtMss%D(AT)
            G%CConstrain%I(NC)=CConstrain%I(AT)
            G%Velocity%D(1:3,NC)=Velocity%D(1:3,AT)
            G%Fext%D(1:3,NC)=Fext%D(1:3,AT)
            G%Carts%D(:,NC)=Carts%D(:,AT)+RVec(:)
            G%Displ%D(:,NC)=Displ%D(:,AT)
            G%Gradients%D(:,NC)=Gradients%D(:,AT)
            G%BoxCarts%D(:,NC)=BoxCarts%D(:,AT)+IVec(:)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    CALL Delete(AtNum)
    CALL Delete(AtTyp)
    CALL Delete(AtNam)
    CALL Delete(AtMss)
    CALL Delete(CConstrain)
    CALL Delete(Velocity)
    CALL Delete(Fext)
    CALL Delete(Carts)
    CALL Delete(Displ)
    CALL Delete(Gradients)
    CALL Delete(BoxCarts)
    !
    G%NAtms=NAtoms
    G%NElec=G%NElec*ISCX*ISCY*ISCZ
    G%TotCh=G%TotCh*ISCX*ISCY*ISCZ
    G%NAlph=G%NAlph*ISCX*ISCY*ISCZ
    G%NBeta=G%NBeta*ISCX*ISCY*ISCZ
    !
    G%BndBox%D(1,2)=RSCX*G%BndBox%D(1,2)
    G%BndBox%D(2,2)=RSCY*G%BndBox%D(2,2)
    G%BndBox%D(3,2)=RSCZ*G%BndBox%D(3,2)
    G%PBC%CellCenter%D(1)=RSCX*G%PBC%CellCenter%D(1)
    G%PBC%CellCenter%D(2)=RSCY*G%PBC%CellCenter%D(2)
    G%PBC%CellCenter%D(3)=RSCZ*G%PBC%CellCenter%D(3)
    !
    G%PBC%BoxShape%D(1:3,1)=RSCX*G%PBC%BoxShape%D(1:3,1)
    G%PBC%BoxShape%D(1:3,2)=RSCY*G%PBC%BoxShape%D(1:3,2)
    G%PBC%BoxShape%D(1:3,3)=RSCZ*G%PBC%BoxShape%D(1:3,3)
    !
    G%PBC%InvBoxSh%D=InverseBoxShape(G%PBC%BoxShape%D,G%PBC%Dimen)
    !
    G%PBC%CellVolume=CellVolume(G%PBC%BoxShape%D,G%PBC%AutoW%I)
    G%PBC%CellCenter%D=CellCenter(G%PBC)

    CALL WrapAtoms(G)
    !
    WRITE(*,*)'================================ 1 ================================================='
    CALL Print_PBCInfo(G%PBC,Unit_O=6)
    CALL PPrint(G,Unit_O=6,PrintGeom_O='XYZ')

!!$    Gtmp%PBC%CellVolume= RSCX*RSCY*RSCZ*G%PBC%CellVolume
!!$    Gtmp%PBC%DipoleFAC = G%PBC%DipoleFAC/RSCX*RSCY*RSCZ
!!$    Gtmp%PBC%QupoleFAC = G%PBC%QupoleFAC/RSCX*RSCY*RSCZ
!!$    Gtmp%PBC%TransVec%D(1)  =RSCX*G%PBC%TransVec%D(1)
!!$    Gtmp%PBC%TransVec%D(2)  =RSCY*G%PBC%TransVec%D(2)
!!$    Gtmp%PBC%TransVec%D(3)  =RSCZ*G%PBC%TransVec%D(3)
  END SUBROUTINE SuperCellMe

!!$  SUBROUTINE ReorderCoordinates(GM)
!!$    TYPE(CRDS)     :: GM
!!$    TYPE(DBL_VECT) :: DTemp
!!$    TYPE(INT_VECT) :: ITemp,Kinds,Point
!!$    TYPE(CHR_VECT) :: CHTemp
!!$    INTEGER        :: J
!!$    !----------------------------------------------------------------------------
!!$!    IF(G%Ordrd==SFC_NONE)RETURN
!!$    !
!!$    CALL New(Point,G%NAtms)
!!$    CALL New(DTemp,G%NAtms)
!!$    CALL New(ITemp,G%NAtms)
!!$    CALL New(CHTemp,G%NAtms)
!!$!    CALL SFCOrder(G%NAtms,G%Carts,Point,G%Ordrd)
!!$    ! Hard wired for now...
!!$    CALL SFCOrder(G%NAtms,G%Carts,Point,SFC_HILBERT)
!!$    !        Reorder Coordinates
!!$    DO I=1,3
!!$       DO J=1,G%NAtms
!!$          DTemp%D(J)=G%Carts%D(I,J)
!!$       ENDDO
!!$       DO J=1,G%NAtms
!!$          G%Carts%D(I,J)=DTemp%D(Point%I(J))
!!$       ENDDO
!!$    ENDDO
!!$    DO I=1,3
!!$       DO J=1,G%NAtms
!!$          DTemp%D(J)=G%BoxCarts%D(I,J)
!!$       ENDDO
!!$       DO J=1,G%NAtms
!!$          G%BoxCarts%D(I,J)=DTemp%D(Point%I(J))
!!$       ENDDO
!!$    ENDDO
!!$    DO J=1,G%NAtms
!!$       ITemp%I(J)=G%AtNum%D(J)
!!$    ENDDO
!!$    DO J=1,G%NAtms
!!$       G%AtNum%D(J)=ITemp%I(Point%I(J))
!!$    ENDDO
!!$    !
!!$    DO J=1,G%NAtms
!!$       DTemp%D(J)=G%AtMss%D(J)
!!$    ENDDO
!!$    DO J=1,G%NAtms
!!$       G%AtMss%D(J)=DTemp%D(Point%I(J))
!!$    ENDDO
!!$    !
!!$    DO J=1,G%NAtms
!!$       CHTemp%C(J)=G%AtNam%C(J)
!!$    ENDDO
!!$    DO J=1,G%NAtms
!!$       G%AtNam%C(J)=CHTemp%C(Point%I(J))
!!$    ENDDO
!!$    !
!!$    CALL Delete(CHTemp)
!!$    CALL Delete(Point)
!!$    CALL Delete(DTemp)
!!$    CALL Delete(ITemp)
!!$  END SUBROUTINE ReorderCoordinates
  !============================================================================
  ! RESCALE VALUES IF ORIGINALLY IN ANGSTROMS
  !============================================================================
  SUBROUTINE ToAtomicUnits(G)
    TYPE(CRDS) :: G

    IF(G%InAU) RETURN
    G%InAU=.TRUE.
    G%Carts%D    = AngstromsToAU*G%Carts%D
    G%Velocity%D =(AngstromsToAU/FemtosecondsToInternalTime)*G%Velocity%D
    G%Fext%D     =(AngstromsToAU/au2eV)*G%Fext%D
    !
    G%PBC%CellCenter%D = G%PBC%CellCenter%D*AngstromsToAU
    G%PBC%BoxShape%D   = AngstromsToAU*G%PBC%BoxShape%D
    G%PBC%InvBoxSh%D   = G%PBC%InvBoxSh%D/AngstromsToAU
    G%PBC%CellVolume   = G%PBC%CellVolume*AngstromsToAU**G%PBC%Dimen
    G%PBC%DipoleFAC    = G%PBC%DipoleFAC/(AngstromsToAU**G%PBC%Dimen)
    G%PBC%QupoleFAC    = G%PBC%QupoleFAC/(AngstromsToAU**G%PBC%Dimen)
    !
  END SUBROUTINE ToAtomicUnits
  !============================================================================
  !
  !============================================================================
  SUBROUTINE PeriodicXLate(G)
    TYPE(CRDS)                  :: G
    REAL(DOUBLE),DIMENSION(1:3) :: CMVec
    INTEGER                     :: I
    !-------------------------------------------------------------------------!
    IF(G%PBC%Translate) THEN
      CMVec(:)=Zero
      DO I=1,G%NAtms
        CMVec(:)=CMVec(:)+G%BoxCarts%D(:,I)
      ENDDO
      CMVec(:)=Half-CMVec(:)/DBLE(G%NAtms)
      G%PBC%TransVec%D(:)=FracToAtom(G,CMVec(:))
      DO I=1,3
        IF(G%PBC%AutoW%I(I)==0) G%PBC%TransVec%D(I)=Zero
      ENDDO
      CALL Translate(G,G%PBC%TransVec%D)
    ELSE
      G%PBC%TransVec%D(:)=Zero
    ENDIF
    !
  END SUBROUTINE PeriodicXLate
END MODULE Massage

