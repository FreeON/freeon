PROGRAM CPSCFSts
!H=================================================================================
!H PROGRAM CPSCFSts
!H
!H  OPTIONS:
!H  DEBUGING: Use -DCPSCFSTS_DBUG to print some stuff.
!H  INFO    : Use -DCPSCFSTS_INFO to print some stuff.
!H
!H Comment:
!H
!H Ref:
!H
!H=================================================================================
  !
#ifdef CPSCFSTS_DBUG
#define CPSCFSTS_INFO
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  !-------------------------------------------------------------------
  TYPE(ARGMT)                      :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                      :: PPrm,T,P,Tmp1,Tmp2
#else
  TYPE(BCSR )                      :: PPrm,T,P,Tmp1,Tmp2
#endif
  !-------------------------------------------------------------------
  REAL(DOUBLE)                     :: DPrimMax,DIISErr,Prop,TmpP
  REAL(DOUBLE)    , DIMENSION(3)   :: Tensor
  INTEGER                          :: I,iXYZ,CPSCFCycl,LastSCFCycle
  INTEGER                          :: RespOrder
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: CPSCFMessage
  CHARACTER(LEN=  DEFAULT_CHR_LEN) :: PropName
  CHARACTER(LEN=*), PARAMETER      :: Prog='CPSCFSts'  
  CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: Cart=(/'X','Y','Z'/)
  !
  !For nulear dipole moment.
  integer :: iatom,NTotAtoms
  real(double) :: ZNuc
  real(double), dimension(3) :: MltN,MltE,COrig
  TYPE(CRDS)                 :: GM
  !-------------------------------------------------------------------
  !
  ! Macro the start up.
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  !
  CPSCFCycl = Args%I%I(1)
  !
  ! Get Last SCF cycle.
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  CALL Get(GM,Tag_O=CurGeom)
  !
  ! Get the response order.
  RespOrder=LEN(TRIM(Args%C%C(3)))
  !
  SELECT CASE(RespOrder)
  CASE(1); PropName='Alpha'
  CASE(2); PropName='Beta'
  CASE(3); PropName='Gamma'
  CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Allocate some matrices.
  !-------------------------------------------------------------------
  !
  CALL New(T   )
  CALL New(P   )
  CALL New(Tmp1)
  CALL New(Tmp2)
  CALL New(PPrm)
  !
  !-------------------------------------------------------------------
  ! Get the density matrices.
  !-------------------------------------------------------------------
  !
  SELECT CASE(SCFActn)
  CASE('BasisSetSwitch','Restart')
     ! If switching the density matrix or using a previous one from 
     ! restart use i+1 density matrix--its all that is available
     CALL Get(PPrm,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,1))
  CASE('StartResponse','CPSCFSolving')
     CALL Get(PPrm,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,1))
  CASE DEFAULT; CALL Halt('Do not know this argument <'//TRIM(SCFActn)//'>.')
  END SELECT
  !CALL Print_BCSR(PPrm,'PPrm',Unit_O=6)
  !
  ! Get groud state density matrix.
  CALL Get(P,TrixFile('D',Args,LastSCFCycle-Args%I%I(1)))
  !     
  !-------------------------------------------------------------------
  ! Compute expectation values.
  !-------------------------------------------------------------------
  !
  IF(CPSCFCycl.LE.0) THEN
     COrig(:)=Zero! must be changed.
     !
     MltN(:)=Zero
     DO iAtom=1,NAtoms
        ZNuc=GM%AtNum%D(iAtom)
        DO iXYZ=1,3
           MltN(iXYZ)=MltN(iXYZ)+ZNuc*(GM%Carts%D(iXYZ,iAtom)-COrig(iXYZ))
        ENDDO
     ENDDO
  ENDIF
  !
  Prop=BIG_DBL
  !
  DO iXYZ=1,3
     !
     ! Get Dipole Moment.
     CALL Get(T,TrixFile(TRIM(Args%C%C(4))//Cart(iXYZ),Args))      ! T=M_{x} or whatever...
     !
     ! Compute tensor elements.
#ifdef PARALLEL
     !old CALL Multiply(PPrm,T,Tmp1)
     !old Tensor(iXYZ)=-Two*Trace(Tmp1)
     CALL Multiply(PPrm,T,Tmp1)
     TmpP=Trace(Tmp1)
     SELECT CASE(RespOrder)
     CASE(1); Tensor(iXYZ)=-Two*TmpP
     CASE(2); Tensor(iXYZ)=-Four*TmpP
     CASE(3); Tensor(iXYZ)=-12.0d0*TmpP
     CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
     END SELECT
     !
     IF(CPSCFCycl.LE.0) THEN
        CALL Multiply(P,T,Tmp1)
        MltE(iXYZ)=-Two*Trace(Tmp1)
     ENDIF
#else
     !old Tensor(iXYZ)=-Two*Trace(PPrm,T)
     TmpP=Trace(PPrm,T)
     SELECT CASE(RespOrder)
     CASE(1); Tensor(iXYZ)=-Two*TmpP
     CASE(2); Tensor(iXYZ)=-Four*TmpP
     CASE(3); Tensor(iXYZ)=-12.0d0*TmpP
     CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
     END SELECT
     !
     IF(CPSCFCycl.LE.0) MltE(iXYZ)=-Two*Trace(P,T)
#endif
     !
     IF('X'.EQ.TRIM(Args%C%C(5))) Prop=Tensor(1)
     IF('Y'.EQ.TRIM(Args%C%C(5))) Prop=Tensor(2)
     IF('Z'.EQ.TRIM(Args%C%C(5))) Prop=Tensor(3)
     IF(    'X'.NE.TRIM(Args%C%C(5)).AND.'Y'.NE.TRIM(Args%C%C(5)).AND. &
          & 'Z'.NE.TRIM(Args%C%C(5))) THEN
        CALL Halt('The args in CPSCFStatus is not one of X,Y or Z, args=' &
             &    //TRIM(Args%C%C(5)))
     ENDIF
  ENDDO
  !
  ! Save Prop.
  CALL Put(Prop,'Prop')
  !write(*,*) 'Prop',Prop
  !
!#ifdef PARALLEL
!  IF(MyID.EQ.ROOT) THEN
!#endif
!     IF(CPSCFCycl.LE.0) THEN
!        WRITE(*,'(A,3E20.12)') 'MltN in a.u.',MltN(1),MltN(2),MltN(3)
!        WRITE(*,'(A,3E20.12)') 'MltE in a.u.',MltE(1),MltE(2),MltE(3)
!        WRITE(*,'(A,3E20.12)') 'MltT in a.u.',MltN(1)+MltE(1), &
!             &                                MltN(2)+MltE(2), &
!             &                                MltN(3)+MltE(3)
!     ENDIF
!#ifdef PARALLEL
!  ENDIF
!#endif
  !-------------------------------------------------------------------
  ! Get max Density block.
  !-------------------------------------------------------------------
  !
  ! Find the largest block of the delta density matrix
  ! Allows for checking between extrapolated or projected DMs
  IF(CPSCFCycl.GT.0) THEN
     !write(*,*) '<'//trim(TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,0))//'>'
     !write(*,*) '<'//trim(TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,1))//'>'
     CALL Get(Tmp1,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL Get(Tmp2,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,1))
     !
     CALL Multiply(Tmp1,-One)
     CALL Add(Tmp1,Tmp2,PPrm)
     !
     ! Get Max PPrm.
     DPrimMax=Max(PPrm)
     !write(*,*) 'DPrimMax',DPrimMax
  ELSE
     DPrimMax=BIG_DBL
  ENDIF
  !
  ! Save DPrimMax.
  CALL Put(DPrimMax,'dprimmax')
  !
  ! Delete some arrays.
  CALL Delete(P   )
  CALL Delete(GM  )
  CALL Delete(PPrm)
  CALL Delete(Tmp1)
  CALL Delete(Tmp2)
  !
  ! Get DDIIS Err.
  IF(Current(1)>=1)THEN
     CALL Get(DIISErr,'ddiiserr')
  ELSE
     DIISErr=Zero
  ENDIF
  !
  !-------------------------------------------------------------------
  ! Print statistics.
  !-------------------------------------------------------------------
  !
  CPSCFMessage=''
  IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
     !
     ! Fancy output.
     CPSCFMessage=RTRN//'= = = = = CPSCFCycle #'//TRIM(SCFCycl)                    &
          &           //', Basis #'//TRIM(CurBase)                                 &
          &           //', Geometry #'//TRIM(CurGeom)                              &
          &           //' = = = = ='//RTRN//RTRN
     !
     ! Add in DDIIS error.
     IF(DIISErr/=Zero)                                                             &
          &       CPSCFMessage=TRIM(CPSCFMessage)                                  &
          &           //'   DDIISErr    = '//TRIM(DblToShrtChar(DIISErr))//RTRN
     !
     ! Add in Tensor elements.
     CPSCFMessage=TRIM(CPSCFMessage)                                               &
          &           //'   MaxDelDPrim = '//TRIM(DblToShrtChar(DPrimMax)) //RTRN  &
          &     //'  '//TRIM(PropName)//' '//TRIM(Args%C%C(3))//'X    = '          &
          &           //TRIM(DblToMedmChar(Tensor(1)))//RTRN                       &
          &     //'  '//TRIM(PropName)//' '//TRIM(Args%C%C(3))//'Y    = '          & 
          &           //TRIM(DblToMedmChar(Tensor(2)))//RTRN                       &
          &     //'  '//TRIM(PropName)//' '//TRIM(Args%C%C(3))//'Z    = '          & 
          &           //TRIM(DblToMedmChar(Tensor(3)))//RTRN
     !
  ELSEIF(PrintFlags%Key>=DEBUG_NONE) THEN
     !
     ! Fancy output.
     CPSCFMessage=ProcessName(Prog,'['//TRIM(SCFCycl)//','                         &
          &           //TRIM(CurBase)//','                                         &
          &           //TRIM(CurGeom)//']') 
     !
     CPSCFMessage=TRIM(CPSCFMessage)//' '//TRIM(PropName)//' '                     &
          &                              //TRIM(Args%C%C(3))//' '                  &
          &                              //'X,Y,Z = '                              &
          &                              //TRIM(FltToShrtChar(Tensor(1)))          &
          &                        //', '//TRIM(FltToShrtChar(Tensor(2)))          &
          &                        //', '//TRIM(FltToShrtChar(Tensor(3)))
     !
     IF(CPSCFCycl.GT.0) THEN
        SELECT CASE(RespOrder)
        CASE(1); CPSCFMessage=TRIM(CPSCFMessage)//', dD1 = '//TRIM(DblToShrtChar(DPrimMax))
        CASE(2); CPSCFMessage=TRIM(CPSCFMessage)//', dD2 = '//TRIM(DblToShrtChar(DPrimMax))
        CASE(3); CPSCFMessage=TRIM(CPSCFMessage)//', dD3 = '//TRIM(DblToShrtChar(DPrimMax))
        CASE DEFAULT; CALL Halt('Response order unknown! RespOrder='//TRIM(IntToChar(RespOrder)))
        END SELECT
     ENDIF
     !
     ! Add in DIIS error.
     IF(DIISErr/=Zero)                                                             &
          &       CPSCFMessage=TRIM(CPSCFMessage)                                  &
          &           //', DDIIS = '//TRIM(DblToShrtChar(DIISErr)) !//RTRN
     !
  ENDIF
  !
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
  !
  ! Printing.
  CALL OpenASCII(OutFile,Out)
  WRITE(*,* )TRIM(CPSCFMessage)
  WRITE(Out,* )TRIM(CPSCFMessage)
  CLOSE(Out)
#ifdef PARALLEL
  ENDIF
#endif
  !
  CALL ShutDown(Prog)
  !
END PROGRAM CPSCFSts

