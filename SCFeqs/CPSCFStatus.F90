PROGRAM CPSCFSts
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
  USE Functionals
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  !-------------------------------------------------------------------
  TYPE(ARGMT)                      :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                      :: PPrim,T,Tmp1,Tmp2,Tmp3
#else
  TYPE(BCSR )                      :: PPrim,T,Tmp1,Tmp2,Tmp3
#endif
  !-------------------------------------------------------------------
  REAL(DOUBLE)                     :: E_el_tot,E_nuc_tot,E_es_tot,E_ECPs,KinE,ExchE,Exc
  REAL(DOUBLE)                     :: Gap,Etot,DPrimMax,Virial,DIISErr,Prop
  REAL(DOUBLE), DIMENSION(3)       :: Tensor
  LOGICAL                          :: HasECPs
  INTEGER                          :: I,iXYZ,CPSCFCycl
  CHARACTER(LEN=DEFAULT_CHR_LEN)   :: PostFix
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: CPSCFMessage
  CHARACTER(LEN=*), PARAMETER      :: Prog='CPSCFSts'  
  CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: Cart=(/'X','Y','Z'/)
  !-------------------------------------------------------------------
  !
  ! Macro the start up.
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  !
  CPSCFCycl = Args%I%I(1)
  !
  ! Allocate some matrices.
  CALL New(T    )
  CALL New(Tmp1 )
  CALL New(Tmp2 )
  CALL New(PPrim)
  !
  ! Get the density matrix.
  IF(SCFActn=='BasisSetSwitch'.OR.SCFActn=='Restart')THEN 
     ! If switching the density matrix or using a previous one from 
     ! restart use i+1 density matrix--its all that is available
     CALL Get(PPrim,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,1))
  ELSE
     !IF(SCFActn=='StartResponse') THEN
     CALL Get(PPrim,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,1))
     !ELSE
     !   CALL Get(PPrim,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,0))
     !ENDIF
  ENDIF
  !
  ! Compute expectation values.
  Tensor(:)=BIG_DBL
  Prop=BIG_DBL
  !
  DO iXYZ=1,3
     !
     ! Get Dipole Moment.
     PostFix=''
     PostFix=TRIM(Args%C%C(3))//Cart(iXYZ)
     !
     CALL Get(T,TrixFile(PostFix,Args))      ! T=M_{x} or whatever...
     !
     ! Compute tensor elements.
#ifdef PARALLEL
     CALL Multiply(PPrim,T,Tmp1)
     Tensor(iXYZ)=-Two*Trace(Tmp1)
#else
     Tensor(iXYZ)=-Two*Trace(PPrim,T)
#endif
     !      IF(Cart(iXYZ).EQ.TRIM(Args%C%C(4))) Prop=Tensor(iXYZ)  ! Dangerous
     IF('X'.EQ.TRIM(Args%C%C(4))) Prop=Tensor(1) 
     IF('Y'.EQ.TRIM(Args%C%C(4))) Prop=Tensor(2) 
     IF('Z'.EQ.TRIM(Args%C%C(4))) Prop=Tensor(3) 
     IF(    'X'.NE.TRIM(Args%C%C(4)).AND.'Y'.NE.TRIM(Args%C%C(4)).AND. &
          & 'Z'.NE.TRIM(Args%C%C(4))) THEN
        CALL Halt('The args in CPSCFStatus is not one of X,Y or Z, args=' &
             &    //TRIM(Args%C%C(4)))
     ENDIF
  ENDDO
  !
  ! Save Prop.
  CALL Put(Prop,'Prop')
  !
  ! Find the largest block of the delta density matrix
  ! Allows for checking between extrapolated or projected DMs
  IF(CPSCFCycl.GT.0) THEN
     CALL Get(Tmp1,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,0))
     CALL Get(Tmp2,TrixFile('DPrime'//TRIM(Args%C%C(4)),Args,1))
     CALL Multiply(Tmp1,-One)
     CALL Add(Tmp1,Tmp2,PPrim)
     !
     ! Get Max PPrim.
     DPrimMax=Max(PPrim)
  ELSE
     DPrimMax=BIG_DBL
  ENDIF
  !
  ! Save DPrimMax.
  CALL Put(DPrimMax,'dprimmax')
  !
  ! Delete some arrays.
  CALL Delete(PPrim)
  CALL Delete(Tmp1 )
  CALL Delete(Tmp2 )
  !
  ! Get DDIIS Err.
  IF(Current(1)>=1)THEN
     CALL Get(DIISErr,'ddiiserr')
  ELSE
     DIISErr=Zero
  ENDIF
  !
  ! Print statistics.
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
          &           //'      DDIISErr = '//TRIM(DblToShrtChar(DIISErr))//RTRN
     !
     ! Add in Tensor elements.
     CPSCFMessage=TRIM(CPSCFMessage)                                               &
          &           //'   MaxDelDPrim = '//TRIM(DblToShrtChar(DPrimMax)) //RTRN  &
          &           //'       <T>     = '//TRIM(DblToMedmChar(Tensor(1)))//RTRN  &
          &           //'       <T>     = '//TRIM(DblToMedmChar(Tensor(2)))//RTRN  &
          &           //'       <T>     = '//TRIM(DblToMedmChar(Tensor(3)))//RTRN
     !
  ELSEIF(PrintFlags%Key>=DEBUG_NONE) THEN
     !
     ! Fancy output.
     CPSCFMessage=ProcessName(Prog,'['//TRIM(SCFCycl)//','                           &
          &           //TRIM(CurBase)//','                                         &
          &           //TRIM(CurGeom)//']') 
     !IF(SCFActn=='BasisSetSwitch')THEN
     !  SCFMessage=TRIM(SCFMessage)//' Basis set switch ... '       &
     !                             //' MxD = '//TRIM(DblToShrtChar(DMax)) 
     !ELSEIF(SCFActn=='Restart')THEN
     !   SCFMessage=TRIM(SCFMessage)//' Restarting ... '       &
     !                              //' MxD = '//TRIM(DblToShrtChar(DMax)) 
     !ELSE
     CPSCFMessage=TRIM(CPSCFMessage)//' T = '//TRIM(FltToShrtChar(Tensor(1))) &
          &                         //', '   //TRIM(FltToShrtChar(Tensor(2))) &
          &                         //', '   //TRIM(FltToShrtChar(Tensor(3))) 

     IF(CPSCFCycl.GT.0) &
          & CPSCFMessage=TRIM(CPSCFMessage)//', dD = '//TRIM(DblToShrtChar(DPrimMax))
     !CPSCFMessage=TRIM(CPSCFMessage)//' T = '//TRIM(FltToShrtChar(Tensor(1))) &
     !     &                         //', '   //TRIM(FltToShrtChar(Tensor(2))) &
     !     &                         //', '   //TRIM(FltToShrtChar(Tensor(3))) &
     !     &                         //', dD = '//TRIM(DblToShrtChar(DPrimMax))
     !ENDIF
     !
     ! Add in DIIS error.
     IF(DIISErr/=Zero)                                            &
          &       CPSCFMessage=TRIM(CPSCFMessage)                 &
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

