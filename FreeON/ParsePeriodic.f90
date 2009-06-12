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

MODULE ParsePeriodic
  USE Parse
  USE InOut
  USE PBC
  USE AtomPairs
  USE PrettyPrint
  USE OptionKeys
  USE PeriodicKeys
  USE ControlStructures

CONTAINS
  !=========================================================================
  !
  !=========================================================================
  SUBROUTINE LoadPeriodic(N,O,G,P)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Periodics)  :: P
    TYPE(PBCInfo)    :: PBC
    INTEGER          :: I,atomID, GBeg,GEnd
    REAL(DOUBLE)     :: totalMass, volume
    !-----------------------------------------------------------------------!
    IF(O%Grad==GRAD_TS_SEARCH_NEB)THEN
      GBeg=0
      GEnd=G%Clones+1
    ELSE
      GBeg=1
      GEnd=G%Clones
    ENDIF
    ! If we are restarting or reguessing, just use values read from HDF ...
    IF(O%Guess==GUESS_EQ_RESTART.OR.O%Guess==GUESS_EQ_NUGUESS)THEN
      ! but check first to see if we are building a supercell
      CALL OpenASCII(N%IFile,Inp)
      DO I=GBeg,GEnd
        CALL SuperCellMe(G%Clone(I)%PBC)
      ENDDO
      CLOSE(UNIT=Inp,STATUS='KEEP')
      RETURN
    ENDIF
    ! ... otherwise, we are reading in new PBC info
    CALL OpenASCII(N%IFile,Inp)
    CALL New(PBC)
    CALL LoadPeriodicOptions(PBC)
    CALL LoadLattice(PBC)
    !
    CLOSE(UNIT=Inp,STATUS='KEEP')
    !
    DO I=GBeg,GEnd
      G%Clone(I)%PBC=PBC
      CALL CalculateCoordArrays(G%Clone(I))

      ! Print out some information regarding the density of the system.
      totalMass = 0.0D0
      DO atomID = 1, G%Clone(I)%NAtms
        totalMass = totalMass+G%Clone(I)%AtMss%D(atomID)
      ENDDO

      CALL MondoLog(DEBUG_NONE, "LoadPeriodic", "mass = " &
        //TRIM(DblToMedmChar(totalMass*amuToKg))//" kg = " &
        //TRIM(DblToMedmChar(totalMass))//" u")

      IF(G%Clone(I)%PBC%Dimen /= 0) THEN
        ! Get the volume.
        volume = CellVolume(G%Clone(I)%PBC%BoxShape%D,G%Clone(I)%PBC%AutoW%I)

        CALL MondoLog(DEBUG_NONE, "LoadPeriodic", "volume = "//TRIM(DblToMedmChar(volume))//" A^3")
        CALL MondoLog(DEBUG_NONE, "LoadPeriodic", "volume per atom = " &
          //TRIM(FltToShrtChar(volume/G%Clone(I)%NAtms))//" A^3 = " &
          //TRIM(FltToShrtChar(volume/G%Clone(I)%NAtms*AngstromsToAU*AngstromsToAU*AngstromsToAU))//" a_B^3")
        CALL MondoLog(DEBUG_NONE, "LoadPeriodic", "density = "//TRIM(DblToMedmChar(totalMass*amuToKg/(volume*1D-30)))//" kg/m^3")
      ENDIF
    ENDDO

  END SUBROUTINE LoadPeriodic
  !=========================================================================
  !
  !=========================================================================
  SUBROUTINE LoadPeriodicOptions(PBC)
    TYPE(PBCInfo)        :: PBC
    INTEGER,DIMENSION(3) :: SC
    INTEGER              :: I,J,NTot,MaxEll
    !-----------------------------------------------------------------------!
    !   Parse the coordinate type
    !
    IF(OptKeyQ(Inp,PBOUNDRY,CRT_FRAC))THEN
      PBC%InAtomCrd=.FALSE.
    ELSE
      PBC%InAtomCrd=.TRUE.
    ENDIF
    !   Parse Periodic Directions
    Ntot = 0
    PBC%AutoW%I=BIG_INT
    IF(FindKey(PBCWRAP,Inp))THEN
      IF(OptKeyLocQ(Inp,PBCWRAP,PBC_TRUE,MaxSets,NLoc,Location)) THEN
        Ntot = NLoc
        DO I=1,NLoc
          PBC%AutoW%I(Location(I)) = 1
        ENDDO
      ENDIF
      PBC%Dimen=NLoc
      IF(OptKeyLocQ(Inp,PBCWRAP,PBC_FALSE,MaxSets,NLoc,Location)) THEN
        Ntot = NTot+NLoc
        DO I=1,NLoc
          PBC%AutoW%I(Location(I)) = 0
        ENDDO
      ENDIF
    ELSE
      PBC%AutoW%I(1:3)=0
      PBC%Dimen=0
    ENDIF
    IF(PBC%Dimen==1.AND.PBC%AutoW%I(1)/=1)THEN
      CALL MondoHalt(PRSE_ERROR,' Wire calculations must be along X, ie PBC=(T,F,F) ')
    ELSEIF(PBC%Dimen==2.AND.PBC%AutoW%I(3)/=0)THEN
      CALL MondoHalt(PRSE_ERROR,' Slab calculations must be along X,Y, ie PBC=(T,T,F) ')
    ENDIF
    ! Parse Translate
    IF(OptKeyQ(Inp,PBOUNDRY,CENTERATOMS))THEN
      PBC%Translate=.TRUE.
    ELSE
      PBC%Translate=.FALSE.
    ENDIF
    ! Check for override of PFF Ell and WS criteria.
    ! Note if WS criteria>1, then strict FMM WS criteria are used.
    ! This can lead to errors in the case of non-cubic cells
    IF(.NOT.OptIntQ(Inp,PFFMXELL,PBC%PFFMaxEll)) THEN
      PBC%PFFMaxEll=16
    ENDIF
    IF(.NOT.OptIntQ(Inp,PFFMXSEP,PBC%PFFWelSep))THEN
      PBC%PFFWelSep=0
    ENDIF
    ! Parse permeability
    IF(.NOT.OptDblQ(Inp,EPSILON,PBC%Epsilon))THEN
      PBC%Epsilon=BIG_DBL ! Numerical infinity
    ENDIF

    ! AtomW feature completely depricated
    ! Atoms are ALWAYS wrapped!
!!$    ! Parse periodic dimensions (Atom Wrap)
!!$    IF(OptKeyQ(Inp,PBOUNDRY,ATOMW_OFF))THEN
!!$       PBC%AtomW=.FALSE.
!!$    ELSE
!!$       PBC%AtomW=.TRUE.
!!$    ENDIF

    ! Generate supercells if asked for
    CALL SuperCellMe(PBC)
  END SUBROUTINE LoadPeriodicOptions
  !
  SUBROUTINE SuperCellMe(PBC)
    TYPE(PBCInfo)        :: PBC
    INTEGER,DIMENSION(3) :: SC
    INTEGER              :: I,K,NTot
    ! Parse supercell options
    SC(:)=0
    IF(FindKey(SUPERC,Inp))THEN
      DO K=1,10
        IF(OptKeyLocQ(Inp,SUPERC,IntToChar(K),MaxSets,NLoc,Location)) THEN
          Ntot = NLoc
          DO I=1,NLoc
            SC(Location(I)) = K
          ENDDO
        ENDIF
      ENDDO
    ELSE
      SC(:)=1
    ENDIF
    IF(SC(1)==0.OR.SC(2)==0.OR.SC(3)==0)THEN
      CALL MondoHalt(PRSE_ERROR,'SuperCell = ('//TRIM(IntToChar(SC(1)))//','// &
           TRIM(IntToChar(SC(2)))//','// &
           TRIM(IntToChar(SC(3)))//') given on input.')
    ENDIF
    PBC%SuperCell%I=1
    IF(PBC%AutoW%I(1)==1)PBC%SuperCell%I(1)=SC(1)
    IF(PBC%AutoW%I(2)==1)PBC%SuperCell%I(2)=SC(2)
    IF(PBC%AutoW%I(3)==1)PBC%SuperCell%I(3)=SC(3)
  END SUBROUTINE SuperCellMe
  !============================================================================
  !
  !============================================================================
  SUBROUTINE LoadLattice(PBC)
    TYPE(PBCInfo)                   :: PBC
    INTEGER                         :: NLvec,NTvec,Dimen,I,J
    CHARACTER(LEN=2)                :: At
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line,LowLine
    REAL(DOUBLE),PARAMETER          :: DegToRad = 1.745329251994329576923D-2
    REAL(DOUBLE),DIMENSION(6)       :: Vec
    REAL(DOUBLE)                    :: AngAB,AngAC,AngBC,Error
    !
    ! Initialize some items that may remain zero
    PBC%TransVec%D=Zero
    PBC%InvBoxSh%D=Zero
    PBC%BoxShape%D=Zero
    PBC%CellCenter%D=Zero
    DO I=1,3
      PBC%BoxShape%D(I,I)=One
      PBC%InvBoxSh%D(I,I)=One
    ENDDO
    !
    IF(PBC%Dimen==0) RETURN
    !
    IF(.NOT. FindKey(BEGIN_PERIODIC,Inp)) THEN
      CALL MondoHalt(PRSE_ERROR,'Lattice Vectors Must Be Suppied for Periodic Systems')
    ENDIF
    !
    NLvec=0
    CALL Align(BEGIN_PERIODIC,Inp)
    DO
      READ(Inp,DEFAULT_CHR_FMT,END=1) Line
      CALL RemoveComments(Line)
      IF(INDEX(Line,END_PERIODIC)==0) THEN
        LowLine=Line
        CALL LowCase(LowLine)
        J=SCAN(LowLine,Lower)
        IF(J/=0)THEN
          CALL LineToGeom(Line,At,Vec)
          IF(At==ALAT_VEC) THEN
            NLvec = NLvec+1
            PBC%BoxShape%D(1:3,1)=Vec(1:3)
          ELSEIF(At==BLAT_VEC) THEN
            NLvec = NLvec+1
            PBC%BoxShape%D(1:3,2)=Vec(1:3)
          ELSEIF(At==CLAT_VEC) THEN
            NLvec = NLvec+1
            PBC%BoxShape%D(1:3,3)=Vec(1:3)
          ENDIF
          PBC%BoxShape%D(2,1)=Zero
          PBC%BoxShape%D(3,1)=Zero
          PBC%BoxShape%D(3,2)=Zero
        ELSE
          PBC%BoxShape%D=0D0
          NLvec=3
          CALL LineToDbls(Line,6,Vec)
          IF(PBC%Dimen==1)THEN
            Vec(4)=90D0
            Vec(5)=90D0
            Vec(6)=90D0
            IF(Vec(2)==0D0)THEN
              Vec(2)=1D0
              CALL Warn(' On input, found b=0.  Reset b=1, and beta=90. ')
            ENDIF
            IF(Vec(3)==0D0)THEN
              Vec(3)=1D0
              CALL Warn(' On input, found c=0.  Reset c=1, and gamma=90. ')
            ENDIF
          ELSEIF(PBC%Dimen==2)THEN
            IF(Vec(3)==0D0)THEN
              Vec(3)=1D0
              CALL Warn(' On input, found c=0.  Reset c=1, and gamma=90. ')
            ENDIF
          ENDIF
          PBC%BoxShape%D(1,1)=Vec(1)
          PBC%BoxShape%D(1,2)=Vec(2)*COS(DegToRad*Vec(6))
          PBC%BoxShape%D(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
          PBC%BoxShape%D(1,3)=Vec(3)*COS(DegToRad*Vec(5))
          PBC%BoxShape%D(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
               -PBC%BoxShape%D(1,2)*PBC%BoxShape%D(1,3))/PBC%BoxShape%D(2,2)
          PBC%BoxShape%D(3,3)=SQRT(Vec(3)**2-PBC%BoxShape%D(1,3)**2-PBC%BoxShape%D(2,3)**2)
        ENDIF
      ELSE
        EXIT
      ENDIF
    ENDDO
    !
    IF(NLVec .NE. 3) THEN
      CALL MondoHalt(PRSE_ERROR,'Lattice Vectors are incorrect')
    ENDIF
    DO I=1,3
      DO J=1,3
        IF(ABS(PBC%BoxShape%D(I,J)).LT. 1.D-12) PBC%BoxShape%D(I,J)=Zero
      ENDDO
    ENDDO
    !
    PBC%InvBoxSh%D=InverseBoxShape(PBC%BoxShape%D,PBC%Dimen)
    !

    RETURN
1   CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
         //TRIM(END_PERIODIC)//'. You may be missing blank '  &
         //'line at the end of the inPut file.')
  END SUBROUTINE LoadLattice
  !=========================================================================
  ! Inflate atomic coordinates if given in fractionals.  If atoms are outside
  ! the unit cell, wrap them back in.  Note that all of this here, and in
  ! MakeGMPeriodic (called by punch HDF), probably behave in a bad way with
  ! respect to velocities etc.   This needs to be double checked for MD restarts
  ! and initial conditions.
  !=========================================================================
  SUBROUTINE CalculateCoordArrays(G)
    TYPE(CRDS) :: G
    INTEGER    :: I
    !-----------------------------------------------------------------------!
    IF(G%PBC%InAtomCrd)THEN
      !-----------------------------------------------------------------------!
      ! The coordinates have been given in atomic coordinates, not fractionals.
      ! Therefore, we need to compute the fractionals (BoxCarts) and the wrapped
      ! atomic positions (Carts).  Note that this is only going to work properly
      ! if the unit cell has been correctly defined at this point.
      !-----------------------------------------------------------------------!
      DO I=1,G%NAtms
        ! Overwrite the atomic positions with fractionals
        G%Carts%D(:,I)=AtomToFrac(G,G%Carts%D(:,I))
        ! wrap the fractionals
        CALL FracCyclic(G,G%Carts%D(:,I))
        ! Boxcarts are the fractional coordinates
        G%BoxCarts%D(:,I)=G%Carts%D(:,I)
        ! Now inflate the wrapped fractionals to get the atomic coordinates back
        G%Carts%D(:,I)=FracToAtom(G,G%Carts%D(:,I))
      ENDDO
    ELSE
      !-----------------------------------------------------------------------!
      ! The coordinates have been given in FRACTIONAL coordinates, not atomic coords.
      ! Following needs to be double checked
      !-----------------------------------------------------------------------!
      G%BoxCarts%D=G%Carts%D
      CALL CalAtomCarts(G)
      CALL WrapAtoms(G)
      !
      !      Convert the Velocities from Fractional to Atomic
      ! Hmm... looks dodgy.  Commented out for now.
!!$       DO I=1,G%NAtms
!!$          G%Velocity%D(:,I)   = FracToAtom(G,G%Velocity%D(:,I))
!!$       ENDDO
    ENDIF
  END SUBROUTINE CalculateCoordArrays
END MODULE ParsePeriodic

