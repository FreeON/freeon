PROGRAM Make1E
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE Make1X
#ifdef PARALLEL
  USE MondoMPI
#endif
  !-------------------------------------------------------------------
  IMPLICIT NONE
  !-------------------------------------------------------------------
  TYPE(ARGMT)    :: Args
  TYPE(BSET)     :: BS
  TYPE(CRDS)     :: GM
  TYPE(INT_VECT) :: Stat
  INTEGER        :: OldFileID,I
  CHARACTER(LEN=*),PARAMETER :: Prog='Make1E'
  !-------------------------------------------------------------------
  !
  ! Start up macro.
  !
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  !
  ! Get basis set and geometry.
  !
  SELECT CASE(SCFActn)
  CASE('RestartGeomSwitch')
     CALL Get(BS,Tag_O=CurBase)
     ! Close current HDF group 
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! Open the old group and HDF
     HDF_CurrentID=OpenHDF(Restart)
     OldFileID=HDF_CurrentID 
     CALL New(Stat,3)
     CALL Get(Stat,'current_state')
     HDF_CurrentID=OpenHDFGroup(HDF_CurrentID,"Clone #"//TRIM(IntToChar(MyClone)))
     ! Get old basis set stuff
     CurGeom=TRIM(IntToChar(Stat%I(3)))
     CurBase=TRIM(IntToChar(Stat%I(2)))
     CALL Get(GM,Tag_O=CurGeom)
     CALL Delete(CS_OUT%CellCarts)
     CALL SetCellNumber(GM,CS_OUT)
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
  CASE('OneElectronMatrices') 
     CALL Get(BS,Tag_O=CurBase)
     CALL Get(GM,Tag_O=CurGeom)
  CASE DEFAULT
     CALL Halt('Make1E: Do not recognize this action <'//TRIM(SCFActn)//'>.')
  END SELECT
  !
  IF(SIZE(Args%C%C).LT.3) CALL Halt('Make1E: Didn''t find any Integral to compute!')
  !
  DO I=3,SIZE(Args%C%C)
     !
     ! Select 1E integral to compute.
     !
     SELECT CASE(Args%C%C(I))
     CASE('Overlap');CALL MakeS(Args,GM,BS)
     CASE('Kinetic');CALL MakeT(Args,GM,BS)
     CASE('EDipole');CALL MakeD(Args,GM,BS)
     !CASE('EQdpole');CALL MakeQ(Args,GM,BS)
     !CASE('');CALL MakeL(Args,GM,BS)
     !CASE('');CALL Make?(Args,GM,BS)
     CASE('NucAttr');CALL MakeV(Args,GM,BS)
     !
     ! Select 1E gradient to compute.
     !
     !CALL SGrad(Agrs,GM,BS)
     CASE DEFAULT
        CALL Halt('Make1E: Do not recognize this Args <'//TRIM(Args%C%C(I))//'>.')
     END SELECT
  ENDDO
  !
  ! Tidy up.
  !
  CALL Delete(BS)
  CALL Delete(GM)
  !
  ! We are done!
  !
  CALL ShutDown(Prog)
  !
END PROGRAM Make1E
