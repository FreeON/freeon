!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines to do Molecular Dyanamics                    !
! Author: Hugh Nymeyer                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!--------------------------------------------------------
!
MODULE MD
!
  USE GlobalCharacters      ! DEFAULT_CHR_LEN
  USE DerivedTypes          ! DBL_VECTOR, etc, etc
  USE GlobalScalars         ! conversion constants
  USE SCFLocals             ! SCFControls, CRDS
  USE MDLocals              ! MDState
  USE MDIO                  ! I/O of special files
  USE MDUtil                ! removeCM, AssignBoltzmannVel....
  USE InOut                 ! Get and Put
  USE AtomPairs             ! WrapAtoms
  USE DrvFrcs               ! Force routines
!
  IMPLICIT NONE
  CONTAINS
!
!--------------------------------------------------------
!********************************************************
!--------------------------------------------------------
!
  SUBROUTINE CALC_MD(Ctrl)
!
! Author: Hugh Nymeyer
! Does: MD setup and MD calls
!
    IMPLICIT NONE
!
    TYPE(SCFControls)               :: Ctrl
    TYPE(CRDS)                      :: GM
    TYPE(MDState)                   :: MDS
    INTEGER                         :: i
!
    IF(Ctrl%MDC%ALGORITHM.EQ.1) THEN   
       Call Velocity_Verlet_Setup(Ctrl,MDS)
       Call Velocity_Verlet(Ctrl,MDS)
    ELSE
       Call HALT('Unknown MD algorithm option')
    END IF
!
  END SUBROUTINE CALC_MD
!
!--------------------------------------------------------
!********************************************************
!--------------------------------------------------------
!
  SUBROUTINE Velocity_Verlet_Setup(Ctrl,MDS)
!
! Author: Hugh Nymeyer
!
    IMPLICIT NONE
!
    TYPE(SCFControls)               :: Ctrl
    TYPE(MDState)                   :: MDS
    TYPE(CRDS)                      :: GM  

    INTEGER :: ISet
!
!------------------------------------
! Miscellaneous setup
!------------------------------------
!
    CtrlVect=SetCtrlVect(Ctrl,'MD')  ! what does this do?
    Call RANDOM_SEED
!
!------------------------------------
! MDS Defaults
!------------------------------------
!
    if (.NOT. Ctrl%MDC%RESTRT) then
       MDS%STEP=1
       MDS%TIME=0.0D0
    end if
!
!------------------------------------
! First set up geometry
!------------------------------------
!
    Call Get(GM,Tag_O=IntToChar(Ctrl%Current(3)))
#ifdef PERIODIC
    Call WrapAtoms(GM)
#endif
!
!------------------------------------
! HACK -- Why aren't AbCarts set?
!------------------------------------
!
!   GM%AbCarts%D = GM%Carts%D
!
!------------------------------------
! then velocities
!------------------------------------
!
!    if (.NOT. Ctrl%MDC%RESTRT) then	
!       if (Ctrl%MDC%TEMP0 .GT. 0.0D0) then
!          Call AssignBoltzmannVelocity(GM,Ctrl%MDC%TEMP0)
!       end if
!    end if
!
!------------------------------------
! zero CM, VCM, and Rotation 
!------------------------------------
!
    if (Ctrl%MDC%REM_TRANS .AND. .NOT. Ctrl%MDC%RESTRT) then 
       Call removeCM(GM)
       Call removeVCM(GM)
    end if
!
    if (Ctrl%MDC%REM_ROTAT) then
       Call removeAngMoment(GM)
    end if
!
!------------------------------------
! read restart if necessary
!------------------------------------
!
    if (Ctrl%MDC%RESTRT) then
       Call OpenRESIN(Ctrl%MDC%RESTRT_IN,MDS%TITLE)
       Call ReadRESIN(Ctrl,MDS,GM)
       Call CloseRESIN
    end if
!
!------------------------------------
! save geometry to current
!------------------------------------
!
    Call Put(GM,Tag_O=IntToChar(Ctrl%Current(3))) 
!
!------------------------------------
! more MDS variables which depend on GM
!------------------------------------
!
    ! Do SCFs possibly over a preliminary number
    ! of basis sets, required untill parallel
    ! makerho works on restart -- MC
    DO ISet=1,Ctrl%NSet-1
       Ctrl%Current=(/0,ISet,1/)
       CALL SetGlobalCtrlIndecies(Ctrl)           
       CALL OneSCF(Ctrl)
    ENDDO
    Ctrl%Current=(/0,Ctrl%NSet,1/)
    CALL SetGlobalCtrlIndecies(Ctrl)           

    Call OneSCF(Ctrl)
    Call Get(MDS%E_POT,'Etot',StatsToChar(Ctrl%Current))
    MDS%E_KIN = computeKE(GM)
    MDS%E_TOT = MDS%E_POT + MDS%E_KIN
    if (.NOT. Ctrl%MDC%RESTRT) &
         MDS%E_TOT_START = MDS%E_TOT
!
!------------------------------------
  END SUBROUTINE Velocity_Verlet_Setup
!
!--------------------------------------------------------
!********************************************************
!--------------------------------------------------------
!
  SUBROUTINE Velocity_Verlet(Ctrl,MDS)
!
! Author: Hugh Nymeyer 10-24-02
! Purpose: Does N steps of verlet integration assuming that
!          all the setup has already been done.
!
! Note: I'm trying to do everything here in units of:
!       energy   : Hartrees
!       distance : Bohrs
!       mass     : Atomic Mass Units
!   --> time     : sqrt(mass*dist**2/energy) 
!
    IMPLICIT NONE
!
    TYPE(SCFControls)             :: Ctrl
    TYPE(MDState)                 :: MDS
    TYPE(CRDS)                    :: GM,GMnew
    INTEGER                       :: ICyc,IBas,IGeo
    INTEGER                       :: Atom       ! atom counter
    TYPE(DBL_VECT)                :: F          ! force array
    INTEGER                       :: i
    INTEGER                       :: STEP
    REAL(DOUBLE)                  :: DT         ! alias for convenience
    REAL(DOUBLE),PARAMETER        :: HALF=0.5D0 ! ditto
!
!--------------------------------------------------------
!   Set aliases
!--------------------------------------------------------
!
    DT=Ctrl%MDC%DT
!
!--------------------------------------------------------
!   Setup
!--------------------------------------------------------
!
    write(6,*)'Velocity Verlet setup'            ! debug
    ICyc=Ctrl%Current(1)
    IBas=Ctrl%Current(2)
    IGeo=Ctrl%Current(3)
    Call Get(GM    ,Tag_O=IntToChar(IGeo)) 
    Call Get(GMnew ,Tag_O=IntToChar(IGeo))   ! to copy all the non-MD variables
    Call New(F,3*GM%NAtms)
!
!--------------------------------------------------------
!   Open Files
!--------------------------------------------------------
!




!
!--------------------------------------------------------
!   compute forces
!--------------------------------------------------------
!
    Ctrl%Project=.TRUE.
    write(6,*)'Computing force once at start'    ! debug
!    Call OneSCF(Ctrl)
    Call Forces(Ctrl)
    Call Get(F,'GradE',Tag_O=IntToChar(IGeo))
    F%D = -F%D    ! Hartrees per Bohr
!
!--------------------------------------------------------
!   enter main dynamics loop
!--------------------------------------------------------
!
    write(6,*)'Main dynamics loop'               ! debug
    do STEP = MDS%STEP,Ctrl%MDC%MAX_STEPS        ! main dynamics loop
!
!      position step: t->t+dt
!
       write(6,*)'Updating positions'            ! debug
#ifdef PERIODIC
       do Atom=1,GM%NAtms
          do i=1,3
             GMnew%Carts%D(i,Atom)      = &
                  GM%Carts%D(i,Atom)    + &
                  DT * GM%Vects%D(i,Atom) + &
                  HALF * DT**2 * F%D(3*(Atom-1)+i) / GM%AtMss%D(Atom)
          end do
       end do
!       Call WrapAtoms(GMnew)
#else
       do Atom=1,GM%NAtms
          do i=1,3
             GMnew%Carts%D(i,Atom)        = &
                  GM%Carts%D(i,Atom)      + &
                  DT * GM%Vects%D(i,Atom) + &
                  HALF * DT**2 * F%D(3*(Atom-1)+i) / GM%AtMss%D(Atom)
          end do
       end do
#endif
!
!      velocity half step: t->t+dt/2
!
       write(6,*)'Updating velocities to t+dt/2'            ! debug
       do Atom=1,GM%NAtms
          do i=1,3
             GMnew%Vects%D(i,Atom)        = &
                  GM%Vects%D(i,Atom)      + &
                  HALF * DT *  F%D(3*(Atom-1)+i) / GM%AtMss%D(Atom)
          end do
       end do
!
!      save coordinates at t+dt and velocities at t+dt/2
!
       write(6,*)'Saving intermediate geometry/velocities'  ! debug
       IGeo = Igeo + 1       
       Ctrl%Current(3)=IGeo
       Call Put(GMnew,Tag_O=IntToChar(IGeo))  
       Call SetGlobalCtrlIndecies(Ctrl)

!
!      compute force at t+dt
!
       write(6,*)'Computing forces at t+dt'                 ! debug
       Call OneSCF(Ctrl)
       Call Forces(Ctrl)
       Call Get(F,'GradE',Tag_O=IntToChar(IGeo))
       F%D = -F%D
!
!      velocity half step: t+dt/2->t+dt     
!
       write(6,*)'Updating velocities to t+dt'              ! debug
       do Atom=1,GM%NAtms
          do i=1,3
             GMnew%Vects%D(i,Atom)        = &
                  GMnew%Vects%D(i,Atom)   + &
                  HALF * DT *  F%D(3*(Atom-1)+i) / GM%AtMss%D(Atom)
          end do
       end do
!
!      remove translation/rotation if needed
!
       if (Ctrl%MDC%REM_TRANS) then
          if (MOD(STEP,Ctrl%MDC%TRANSfreq) .EQ. 0) then
             Call removeCM(GM)
             Call removeVCM(GM)
          end if
       end if
!
       if (Ctrl%MDC%REM_ROTAT) then
          if (MOD(STEP,Ctrl%MDC%ROTATfreq) .EQ. 0) then
             Call removeAngMoment(GM)
          end if
       end if
!
!      save geometry again
!
       write(6,*)'Saving final geometry at time t+dt'       ! debug
       Call Put(GMnew,Tag_O=IntToChar(IGeo))  
!
!      update MDS
!
       MDS%STEP = STEP
       MDS%TIME = MDS%TIME + DT 
       MDS%E_KIN = computeKE(GMnew)
       call Get(MDS%E_POT,'Etot',StatsToChar(Ctrl%Current))
       MDS%E_TOT = MDS%E_POT + MDS%E_KIN
       MDS%TEMPERATURE = computeTinstant(Ctrl,GM)
       MDS%DE_TOT = MDS%E_TOT - MDS%E_TOT_START
!
!      output  CHANGED BY MC TO FLUSH AT EACH CYCLE...
!
       if ( MOD(STEP,Ctrl%MDC%CRDfreq) .EQ. 0)THEN
          Call OpenCRD(Ctrl%MDC%CRD_NAME,MDS%TITLE,Ctrl%MDC%CLOBBER)
          Call WriteCRD(GM,Ctrl)
          CALL CloseCRD
       ENDIF
       if ( MOD(STEP,Ctrl%MDC%VELfreq) .EQ. 0)THEN
          Call OpenVEL(Ctrl%MDC%VEL_NAME,MDS%TITLE,Ctrl%MDC%CLOBBER)
          Call WriteVEL(GM)
          Call CloseVEL
       ENDIF
       if ( MOD(STEP,Ctrl%MDC%ENEfreq) .EQ. 0)THEN
          Call OpenENE(Ctrl%MDC%ENE_NAME,MDS%TITLE,Ctrl%MDC%CLOBBER)
          Call WriteENE(GM,Ctrl,MDS)
          Call CloseENE
       ENDIF
       if ( MOD(STEP,Ctrl%MDC%RESfreq) .EQ. 0)THEN
          Call OpenRESOUT(Ctrl%MDC%RESTRT_OUT,MDS%TITLE,Ctrl%MDC%CLOBBER)
          Call WriteRESOUT(Ctrl,MDS,GM)
          Call CloseRESOUT
       ENDIF
       !
!      update GM structures
!
       Call Get(GM,Tag_O=IntToChar(IGeo)) 
!
!--------------------------------------------------------
!   end main dynamics loop
!--------------------------------------------------------
!
    end do
!
!--------------------------------------------------------
!   cleanup
!--------------------------------------------------------
!



    Call Delete(F)
!
!--------------------------------------------------------
!
  END SUBROUTINE VELOCITY_VERLET
!--------------------------------------------------------
!********************************************************
!--------------------------------------------------------
!
  END MODULE
