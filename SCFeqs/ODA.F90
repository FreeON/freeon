PROGRAM ODA 
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE DenMatMethods
  IMPLICIT NONE
  !-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
#ifdef PARALLEL
  TYPE(DBCSR)                    ::  &
#else
  TYPE(BCSR)                     ::  &
#endif
                                    P,PTilde,F,FTilde,T,K0,K1,T1,T2,T3
  REAL(DOUBLE)                   :: e0,e1,e0p,e1p,a3,b3,c3,d3,EMns,EPls,EMin,        &
                                    LMns,LPls,L,L1,ENucTotTilde,                     &
                                    DIISErr,Enuc0,Enuc1,Exc0,Exc1
  REAL(DOUBLE)                   :: Tmp1,Tmp2,Tmp3,Tmp4,alph
  INTEGER                        :: I,iSCF
  LOGICAL                        :: Present,HasECPs
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,MatFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='ODA'
  REAL(DOUBLE)                   :: TrP0T,TrP1T,TrP0F0,TrP1F1,TrP0F1,TrP1F0, &
                                    TrP0K1,TrP0K0,TrP1K0,TrP1K1,EASM
  !-------------------------------------------------------------------------
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  iSCF=Args%I%I(1)
! Suss for matrix threshold overide
  CALL SussTrix('ODATrix',Prog)  
! Allocate All Read from Disk Matrices, what a waste
  CALL New(P)
  CALL New(PTilde)
  CALL New(F)
  CALL New(FTilde)
! Allocate temp Matrices
  CALL New(T1)
  CALL New(T2)  
  CALL New(T3)
! Get the Matrices
  CALL Get(PTilde,TrixFile('D',Args,-1))   
  CALL Get(FTilde,TrixFile('F',Args,-1))  
  CALL Get(P,TrixFile('D',Args,0))
  CALL Get(F,TrixFile('F',Args,0))  
! Get Kinectic Energy and ECPs
  CALL Get(T,TrixFile('T',Args))
  CALL Get(HasECPs,'hasecps',Tag_O=CurBase)
  IF(HasECPs) THEN
     CALL Get(T1,TrixFile('U',Args))
     CALL Add(T,T1,T2)
     CALL SetEq(T,T2)
  ENDIF
! Calculate Exchange Asymmetry
  IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
     IF(HasHF(ModelChem)) THEN
        CALL Get(K0,TrixFile('K',Args,-1))
        CALL Get(K1,TrixFile('K',Args,0))
        CALL OpenASCII(OutFile,Out)
#ifdef PARALLEL
        CALL Multiply(PTilde,K1,T1)
        TrP0K1 =  Trace(T1)
        CALL Multiply(P,K0,T1)
        TrP1K0 =  Trace(T1)
        EASM   =  ABS(TrP0K1-TrP1K0)/ABS(TrP0K1+TrP1K0)
#else
        TrP0K1 =  Trace(PTilde,K1)
        TrP1K0 =  Trace(P,K0)
        EASM   =  ABS(TrP0K1-TrP1K0)/ABS(TrP0K1+TrP1K0)
#endif
        WRITE(Out,'(a18,D11.5)') " Tr[P0*K1-P1*K0] = ",EASM
        WRITE(*  ,'(a18,D11.5)') " Tr[P0*K1-P1*K0] = ",EASM
        CLOSE(Out)
     ENDIF
  ENDIF
! Get Energies: E_nuc and E_xc and K_xc matrices
  Current(1)=Current(1)-1
  CALL Get(Enuc0,'E_NuclearTotal',StatsToChar(Current))
  Current(1)=Current(1)+1
  CALL Get(Enuc1,'E_NuclearTotal',StatsToChar(Current))
  IF(HasDFT(ModelChem)) THEN
     Current(1)=Current(1)-1
     CALL Get(Exc0 ,'Exc'           ,StatsToChar(Current))
     CALL Get(K0,TrixFile('Kxc',Args,-1))
     Current(1)=Current(1)+1
     CALL Get(Exc1 ,'Exc'           ,StatsToChar(Current))
     CALL Get(K1,TrixFile('Kxc',Args,0))
  ENDIF
! Compute the Endpoint energies and Derivatives
#ifdef PARALLEL
! Compute the Traces
  CALL Multiply(PTilde,T,T1)
  TrP0T = Trace(T1)
  CALL Multiply(P,T,T1)
  TrP1T = Trace(T1)
  CALL Multiply(PTilde,FTilde,T1)
  TrP0F0 = Trace(T1)
  CALL Multiply(P,FTilde,T1)
  TrP1F0 = Trace(T1)
  CALL Multiply(PTilde,F,T1)
  TrP0F1 = Trace(T1)
  CALL Multiply(P,F,T1)
  TrP1F1 = Trace(T1)
  IF(HasDFT(ModelChem)) THEN
     CALL Multiply(PTilde,K0,T1)
     TrP0K0 = Trace(T1)
     CALL Multiply(P,K1,T1)
     TrP1K1 = Trace(T1)
     CALL Multiply(PTilde,K1,T1)
     TrP0K1 = Trace(T1)
     CALL Multiply(P,K0,T1)
     TrP1K0 = Trace(T1)
  ENDIF  
  e0  = TrP0T+TrP0F0+Enuc0
  e1  = TrP1T+TrP1F1+Enuc1
  e0p = Enuc1-Enuc0+TrP1T-TrP0T+TrP1F0+TrP0F1-Two*TrP0F0
  e1p = Enuc1-Enuc0+TrP1T-TrP0T+Two*TrP1F1-TrP1F0-TrP0F1
  IF(HasDFT(ModelChem)) THEN
     e0  = e0  + Exc0 - TrP0K0
     e1  = e1  + Exc1 - TrP1K1
     e0p = e0p + TrP1K0-TrP0K1
     e1p = e1p + TrP1K0-TrP0K1
  ENDIF
#else
! Avoid assumption of two electron integral symmetry which may be lost 
! in the case of small cell PBC HF and also due to excesive thresholding. 
  e0  = Trace(PTilde,T)+Trace(PTilde,FTilde) + Enuc0
  e1  = Trace(P,T)     +Trace(P,F)           + Enuc1
  e0p = Enuc1-Enuc0+Trace(P,T)-Trace(PTilde,T)+Trace(P,FTilde)+Trace(PTilde,F)-Two*Trace(PTilde,FTilde)
  e1p = Enuc1-Enuc0+Trace(P,T)-Trace(PTilde,T)+Two*Trace(P,F)-Trace(P,FTilde)-Trace(PTilde,F)
  IF(HasDFT(ModelChem)) THEN
     e0  = e0  + Exc0 - Trace(PTilde,K0)
     e1  = e1  + Exc1 - Trace(P,K1)
     e0p = e0p + (Trace(P,K0)-Trace(PTilde,K1))
     e1p = e1p + (Trace(P,K0)-Trace(PTilde,K1))
  ENDIF
#endif
! Find the mixing parameter L from the
! cubic E3(L)=a3+b3*L+c3*L^2+d3*L^3
  a3=e0
  b3=e0p
  c3=-3D0*e0-2D0*e0p+3D0*e1-e1p
  d3=2D0*e0+e0p-2D0*e1+e1p
  IF(ABS(d3)<1D-6.OR.c3*c3-3*b3*d3<Zero)THEN
     L=-Half*b3/c3
     L1=One-L
     EMin=a3+b3*L+c3*L**2+d3*L**3
  ELSE
     LMns=MAX(Zero,MIN(One,(-c3-SQRT(c3*c3-3*b3*d3))/(3*d3)))
     LPls=MAX(Zero,MIN(One,(-c3+SQRT(c3*c3-3*b3*d3))/(3*d3)))
     EMns=a3+b3*LMns+c3*LMns**2+d3*LMns**3
     EPls=a3+b3*LPls+c3*LPls**2+d3*LPls**3
     IF(EMns<EPls)THEN
        L=LMns
        L1=One-LMns
     ELSE
        L=LPls
        L1=One-LPls
     ENDIF
     EMin=MIN(EMns,EPls)
  ENDIF
  ! End point checks
  IF(L<=Zero.OR.L>One)THEN
     IF(e0<e1)THEN
        L=0.001
        L1=One-L
        EMin=a3+b3*L+c3*L**2
     ELSE
        L=One
        L1=Zero
        EMin=e1
     ENDIF
  ENDIF
#ifdef PARALLEL
  IF(MyId==ROOT)THEN
#endif
     Mssg=' Mix = '//TRIM(FltToShrtChar(L))
     Mssg=ProcessName(Prog,TRIM(Mssg))
     Mssg=TRIM(Mssg)//" <SCF> = "//TRIM(FltToMedmChar(EMin))//', d3 = '//TRIM(DblToShrtChar(d3))
     CALL OpenASCII(OutFile,Out)
     WRITE(*,*)TRIM(Mssg)
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
#ifdef PARALLEL
  ENDIF
#endif
! Compute PTilde_N=(1-L)*PTilde_(N-1)+L*P_N, then put bto disk
  CALL Multiply(P,L)
  CALL Multiply(PTilde,L1)
  CALL Add(P,PTilde,T2)
  CALL Put(T2,TrixFile('D',Args,0))
  CALL Put(T2,'CurrentDM',CheckPoint_O=.TRUE.)
  CALL PChkSum(T2,'Pao['//TRIM(NxtCycl)//']',Prog)
! Compute FTilde_N ~ (1-L)*FTilde_(N-1)+L*F_N
  CALL Multiply(F,L)
  CALL Multiply(FTilde,L1)
  CALL Add(F,FTilde,T3)
  CALL Put(T3,TrixFile('F',Args,0)) 
  CALL PChkSum(T3,'Fao['//TRIM(NxtCycl)//']',Prog)
!----------------------------------------------------------------------
! Convert FTilde to ortho rep
  INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
  IF(Present)THEN
     CALL Get(K0,TrixFile('X',Args))         ! Z=S^(-1/2)
!    F_ortho
     CALL Multiply(K0,T3,T1)                 ! Z*F
     CALL Multiply(T1,K0,T3)                 ! (Z*F)*Z
     CALL Filter(T1,T3)                      ! Filter
     CALL Put(T1,TrixFile('OrthoF',Args,0))  
     CALL PChkSum(T1,'For['//TRIM(NxtCycl)//']',Prog)
     CALL SetEq(FTilde,T1)
  ELSE
     CALL Get(K0,TrixFile('ZT',Args))        ! ZT=S^(L)
     CALL Get(K1,TrixFile('Z',Args))         ! Z=S^(-L)
!    F_ortho
     CALL Multiply(K0,T3,T1)                 ! ZT*F
     CALL Multiply(T1,K1,T3)                 ! (ZT*F)*Z
     CALL Filter(T1,T3)                      ! Filter
     CALL Put(T1,TrixFile('OrthoF',Args,0))  
     CALL PChkSum(T1,'For['//TRIM(NxtCycl)//']',Prog)
     CALL SetEq(FTilde,T1)
  ENDIF
! Convert PTilde to ortho rep
  IF(iSCF>1)THEN
     CALL Get(PTilde,TrixFile('OrthoD',Args,-1))   
     CALL Get(P,TrixFile('OrthoD',Args,0))   
     CALL Multiply(P,L)
     CALL Multiply(PTilde,L1)
     CALL Add(P,PTilde,T1)
     CALL PChkSum(T1,'Por['//TRIM(NxtCycl)//']',Prog)
     CALL Put(T1,TrixFile('OrthoD',Args,0))
     CALL SetEq(PTilde,T1)
!    Compute the DIIS error
     CALL Multiply(FTilde,PTilde,T1)
     CALL Multiply(PTilde,FTilde,T2)
     CALL Multiply(T2,-One)
     CALL Add(T1,T2,T3)
     DIISErr=SQRT(Dot(T3,T3))/DBLE(NBasF)
     CALL Put(DIISErr,'diiserr')
  ELSE
     DIISErr=One
     CALL Put(DIISErr,'diiserr')
  ENDIF
! Store Emin
  CALL Put(EMin   ,'ODAEnergy')
! JTilde_N = (1-L)*JTilde_(N-1)+L*J_N
  CALL Get(PTilde,TrixFile('J',Args,-1))
  CALL Get(P,     TrixFile('J',Args,0))
  CALL Multiply(P,L)
  CALL Multiply(PTilde,L1)
  CALL Add(P,PTilde,T1)
  CALL Put(T1,TrixFile('J',Args,0))
! KTilde_N = (1-L)*KTilde_(N-1)+L*K_N
  IF(HasHF(ModelChem))THEN
     CALL Get(PTilde,TrixFile('K',Args,-1))
     CALL Get(P     ,TrixFile('K',Args,0))
     CALL Multiply(P,L)
     CALL Multiply(PTilde,L1)
     CALL Add(P,PTilde,T1)
     CALL Put(T1,TrixFile('K',Args,0))
  ENDIF
! ENucTotTilde_N = (1-L)*ENucTotTilde_(N-1)+L*ENucTotTilde_N
  ENucTotTilde=L*Enuc1+L1*Enuc0
  CALL Put(ENucTotTilde,'E_NuclearTotal',StatsToChar(Current))
! Tidy up
  CALL Delete(P)
  CALL Delete(PTilde)
  CALL Delete(F)
  CALL Delete(FTilde)
! Delete Temps
  CALL Delete(T1) 
  CALL Delete(T2)
  CALL Delete(T3)
  CALL ShutDown(Prog)
!
END PROGRAM ODA
