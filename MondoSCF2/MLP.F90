MODULE MLP
!H=================================================================================
!H
!H
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Parse
  USE LinAlg
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: MLPDriver
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: MLPDpol
  !
CONTAINS
  !
  SUBROUTINE MLPDriver(D,GM,Args)
!H---------------------------------------------------------------------------------
!H SUBROUTINE MLPDriver()
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL_PRP
    TYPE(DBCSR), INTENT(IN) :: D
#else
    TYPE(BCSR ), INTENT(IN) :: D
#endif
    TYPE(CRDS ), INTENT(IN) :: GM
    TYPE(ARGMT), INTENT(IN) :: Args
    !-------------------------------------------------------------------
    TYPE(DBL_VECT)          :: COrig
    !-------------------------------------------------------------------
    !
    ! Get coordinate origine.
    !TODO
    !
    ! Compute dipole moment.
    CALL MLPDpol(D,GM,COrig,Args)
    !
  END SUBROUTINE MLPDriver
  !
  SUBROUTINE MLPDpol(D,GM,COrig,Args)
!H---------------------------------------------------------------------------------
!H SUBROUTINE MLPDpol()
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
#ifdef PARALLEL_PRP
    TYPE(DBCSR)   , INTENT(IN)     :: D
#else
    TYPE(BCSR )   , INTENT(IN)     :: D
#endif
    TYPE(CRDS )   , INTENT(IN)     :: GM
    TYPE(DBL_VECT), INTENT(IN)     :: COrig
    TYPE(ARGMT)   , INTENT(IN)     :: Args
    !-------------------------------------------------------------------
#ifdef PARALLEL_PRP
    TYPE(DBCSR)                    :: Tmp1,Tmp2
#else
    TYPE(BCSR )                    :: Tmp1
#endif
    INTEGER                        :: iXYZ,iAtom
    REAL(DOUBLE)                   :: ZNuc
    REAL(DOUBLE), DIMENSION(3)     :: MltE,MltN
    CHARACTER(LEN=1), DIMENSION(3), PARAMETER :: Cart=(/'x','y','z'/)
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: MltMessage
    !-------------------------------------------------------------------
    !
    CALL New(Tmp1)
#ifdef PARALLEL_PRP
    CALL New(Tmp2)
#endif
    !
    ! Compute electronic dipole.
    MltE(:)=Zero
    DO iXYZ=1,3
       CALL Get(Tmp1,TrixFile('D'//Cart(iXYZ),Args))
#ifdef PARALLEL_PRP
       CALL Multiply(D,Tmp1,Tmp2)
       MltE(iXYZ)=Two*Trace(Tmp2)
#else
       MltE(iXYZ)=Two*Trace(D,Tmp1)    
#endif
       !WRITE(*,'(A,A,A,3E20.12)') ' D',Cart(iXYZ),' = ',Mlt(iXYZ)
    ENDDO
    !
    ! Compute nuclear dipole.
    MltN(:)=Zero
    DO iAtom=1,NAtoms
       ZNuc=DBLE(GM%AtNum%D(iAtom))
       DO iXYZ=1,3
          MltN(iXYZ)=MltN(iXYZ)+ZNuc*(GM%Carts%D(iXYZ,iAtom)-COrig%D(iXYZ))
       ENDDO
    ENDDO
    !
    ! Print dipole.
    MltMessage=""
    MltMessage=RTRN//' Electronic Dipole Moment: X ='//TRIM(DblToMedmChar(MltE(1)))// &
         &                                    ', Y ='//TRIM(DblToMedmChar(MltE(2)))// &
         &                                    ', Z ='//TRIM(DblToMedmChar(MltE(3)))// &
         &     RTRN//' Nuclear Dipole Moment   : X ='//TRIM(DblToMedmChar(MltN(1)))// &
         &                                    ', Y ='//TRIM(DblToMedmChar(MltN(2)))// &
         &                                    ', Z ='//TRIM(DblToMedmChar(MltN(3)))// &
         &     RTRN//' Total Dipole Moment     : X ='//TRIM(DblToMedmChar(MltE(1)+MltN(1)))// &
         &                                    ', Y ='//TRIM(DblToMedmChar(MltE(2)+MltN(2)))// &
         &                                    ', Z ='//TRIM(DblToMedmChar(MltE(3)+MltN(3)))
    !
#ifdef PARALLEL_PRP
    IF(MyId==ROOT)THEN
#endif
       CALL OpenASCII(OutFile,Out)
       WRITE(*  ,* )TRIM(MltMessage)
       WRITE(Out,* )TRIM(MltMessage)
       CLOSE(Out)
#ifdef PARALLEL_PRP
    ENDIF
#endif
    !
    ! Delete some stuff.
    CALL Delete(Tmp1)
#ifdef PARALLEL_PRP
    CALL Delete(Tmp2)
#endif
    !
  END SUBROUTINE MLPDpol
  !
END MODULE MLP


