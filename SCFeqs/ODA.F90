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
  TYPE(DBCSR)                     ::  &
#else
  TYPE(BCSR)                      ::  &
#endif
                                    P,PTilde,F,FTilde,T1,T2,S
  REAL(DOUBLE)                   :: e0,e1,e0p,e1p,a3,b3,c3,d3,EMns,EPls,EMin, &
       LMns,LPls,L,L1,ENucTot,ENucTotTilde,ExcTilde,Exc,DIISErr
  INTEGER                        :: I,iSCF
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,MatFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='ODA'
  INTEGER,PARAMETER              :: M=-1
  !-------------------------------------------------------------------------
#ifdef PARALLEL
  CALL StartUp(Args,Prog,SERIAL_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
  iSCF=Args%I%I(1)
  ! Suss for matrix threshold overide
  CALL SussTrix('ODA',Prog)  
  ! Allocations 
  CALL New(T1)
  CALL New(T2)
  CALL New(S)
  CALL New(P)
  CALL New(F)
  CALL New(PTilde)
  CALL New(FTilde)
#ifdef OrthogonalODA  
  CALL Get(PTilde,TrixFile('OrthoD',Args,M))   
  CALL Get(FTilde,TrixFile('OrthoF',Args,M))   
#else
  CALL Get(PTilde,TrixFile('D',Args,M))   
  CALL Get(FTilde,TrixFile('F',Args,M))   
#endif
  Current(1)=Current(1)+M
  CALL Get(e0,'Etot',StatsToChar(Current))
  Current(1)=Current(1)-M
  ! Get the current (N) non-tilde values
#ifdef OrthogonalODA  
  CALL Get(P,TrixFile('OrthoD',Args,0))   
  CALL Get(F,TrixFile('OrthoF',Args,0))   
#else
  CALL Get(P,TrixFile('D',Args,0))   
  CALL Get(F,TrixFile('F',Args,0))   
#endif
  ! T1 = P_N-PTilde_{N-1}
  CALL Multiply(PTilde,-One)
  CALL Add(PTilde,P,T1)  
  CALL Get(e1,'Etot',StatsToChar(Current))
  ! Get the previous Fock matrix F[PTilde_(N-1)]
#ifdef PARALLEL
  CALL Multiply(FTilde,T1,S)
  e0p=Two*Trace(S)
  CALL Multiply(F,T1,S)
  e1p=Two*Trace(S)
#else
  e0p=Two*Trace(FTilde,T1)
  e1p=Two*Trace(F,T1)
#endif
  ! WRITE(*,*)' e0 = ',e0
  ! WRITE(*,*)' e1 = ',e1
  ! WRITE(*,*)'e0p = ',e0p
  ! WRITE(*,*)'e1p = ',e1p
  ! Find the mixing parameter L from the
  ! cubic E3(L)=a3+b3*L+c3*L^2+d3*L^3
  a3=e0
  b3=e0p
  c3=-3D0*e0-2D0*e0p+3D0*e1-e1p
  d3=2D0*e0+e0p-2D0*e1+e1p
  IF(ABS(d3)<1D-6)THEN
     L=-Half*b3/c3
     L1=One-L
     EMin=a3+b3*L+c3*L**2
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
        ! Mark of the beast
        L=6.66D-3
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
     Mssg=TRIM(Mssg)//" <SCF> ~ "//TRIM(FltToMedmChar(EMin))//', d3 = '//TRIM(DblToShrtChar(d3))
     CALL OpenASCII(OutFile,Out)
     WRITE(*,*)TRIM(Mssg)
     WRITE(Out,*)TRIM(Mssg)
     CLOSE(Out)
#ifdef PARALLEL
  ENDIF
#endif
  ! PTilde_N=(1-L)*PTilde_(N-1)+L*P_N
  CALL Multiply(P,L)
  CALL Multiply(PTilde,-L1)
  CALL Add(P,PTilde,T1)
  CALL SetEq(PTilde,T1)
  ! FTilde_N ~ (1-L)*FTilde_(N-1)+L*F_N
  CALL Multiply(F,L)
  CALL Multiply(FTilde,L1)
  CALL Add(F,FTilde,T1)
  CALL SetEq(FTilde,T1)
#ifdef OrthogonalODA  
  ! Compute the orthogonal DIIS error
  CALL Multiply(FTilde,PTilde,T1)  
  CALL Multiply(PTilde,FTilde,T1,-One)
#else
  CALL Get(S,TrixFile('S',Args))     
  ! Compute the AO-DIIS error
  CALL Multiply(FTilde,PTilde,T1)
  CALL Multiply(T1,S,P)
  ! F.P.S-S.P.F
  CALL Multiply(S,PTilde,T1)
  CALL Multiply(T1,FTilde,F)
  CALL Multiply(F,-One)
  CALL Add(P,F,T1)
#endif
  ! The DIIS error 
  DIISErr=SQRT(Dot(T1,T1))/DBLE(NBasF)
  CALL Put(DIISErr,'diiserr')
  CALL Put(EMin,'ODAEnergy')
  ! Orthogonal put and xform to AO rep and put of PTilde
  ! Step on density from both this cycle and the next
!   CALL Put(PTilde,TrixFile('OrthoD',Args,0))
#ifdef OrthogonalODA  
  ! Convert P to AO representation
  INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
  IF(Present)THEN
     CALL Get(P,TrixFile('X',Args))   ! Z=S^(-1/2)
     CALL Multiply(P,PTilde,T1)
     CALL Multiply(T1,P,PTilde)
  ELSE
     CALL Get(P,TrixFile('Z',Args))   ! Z=S^(-L)
     CALL Multiply(P,PTilde,T1)
     CALL Get(P,TrixFile('ZT',Args))
     CALL Multiply(T1,P,PTilde)
  ENDIF
  CALL Put(FTilde,TrixFile('OrthoF',Args,0)) 
  CALL Put(PTilde,TrixFile('D',Args,0))
  CALL PChkSum(FTilde,'For~['//TRIM(NxtCycl)//']',Prog)
  CALL PChkSum(PTilde,'Pao~['//TRIM(NxtCycl)//']',Prog)
  IF(iSCF>1)THEN
     ! Step on the orthogonal DM as well
     CALL Get(PTilde,TrixFile('OrthoD',Args,M))   
     CALL Get(P,TrixFile('OrthoD',Args,0))   
     CALL Multiply(P,L)
     CALL Multiply(PTilde,L1)
     CALL Add(P,PTilde,T1)
     CALL Put(T1,TrixFile('OrthoD',Args,0))
  ENDIF
#else
  CALL Put(FTilde,TrixFile('F',Args,0)) 
  CALL PChkSum(FTilde,'Fao~['//TRIM(NxtCycl)//']',Prog)
  ! Convert F to orthogonal representation
  INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
  IF(Present)THEN
     CALL Get(P,TrixFile('X',Args))   ! Z=S^(-1/2)
     CALL Multiply(P,FTilde,T1)
     CALL Multiply(T1,P,FTilde)
  ELSE
     CALL Get(P,TrixFile('ZT',Args))
     CALL Multiply(P,FTilde,T1)
     CALL Get(P,TrixFile('Z',Args))  ! Z=S^(-L)
     CALL Multiply(T1,P,FTilde)
  ENDIF
  CALL Put(FTilde,TrixFile('OrthoF',Args,0)) 
  CALL Put(PTilde,TrixFile('D',Args,0))
  CALL PChkSum(FTilde,'For~['//TRIM(NxtCycl)//']',Prog)
  CALL PChkSum(PTilde,'Pao~['//TRIM(NxtCycl)//']',Prog)
  IF(iSCF>1)THEN
     ! Step on the orthogonal DM as well
     CALL Get(PTilde,TrixFile('OrthoD',Args,M))   
     CALL Get(P,TrixFile('OrthoD',Args,0))   
     CALL Multiply(P,L)
     CALL Multiply(PTilde,L1)
     CALL Add(P,PTilde,T1)
     CALL Put(T1,TrixFile('OrthoD',Args,0))
  ENDIF
#endif
  ! JTilde_N ~ (1-L)*JTilde_(N-1)+L*J_N
  CALL Get(PTilde,TrixFile('J',Args,M))
  CALL Get(P,TrixFile('J',Args,0))
  CALL Multiply(P,L)
  CALL Multiply(PTilde,L1)
  CALL Add(P,PTilde,T1)
  CALL Put(T1,TrixFile('J',Args,0))
  ! KTilde_N ~ (1-L)*KTilde_(N-1)+L*K_N
  IF(HasHF(ModelChem))THEN
     CALL Get(PTilde,TrixFile('K',Args,M))
     CALL Get(P,TrixFile('K',Args,0))
     CALL Multiply(P,L)
     CALL Multiply(PTilde,L1)
     CALL Add(P,PTilde,T1)
     CALL Put(T1,TrixFile('K',Args,0))
  ENDIF
  ! ENucTotTilde_N ~ (1-L)*ENucTotTilde_(N-1)+L*K_N
  Current(1)=Current(1)+M
  CALL Get(ENucTotTilde,'E_NuclearTotal',StatsToChar(Current))
  Current(1)=Current(1)-M
  CALL Get(ENucTot,'E_NuclearTotal',StatsToChar(Current))
  ENucTotTilde=L*ENucTot+L1*ENucTotTilde
  CALL Put(ENucTotTilde,'E_NuclearTotal',StatsToChar(Current))
  ! Here is the big punt in Kohn-Sham ODA; we are 
  ! doing linear interpolation of Exc.  This may suck worse than
  ! the cubic approximation, and we should then recompute Kxc, Exc.
  IF(HasDFT(ModelChem))THEN
     Current(1)=Current(1)+M
     CALL Get(ExcTilde,'Exc',StatsToChar(Current))
     Current(1)=Current(1)-M
     CALL Get(Exc,'Exc',StatsToChar(Current))
     ExcTilde=L*Exc+L1*ExcTilde
     CALL Put(ExcTilde,'Exc',StatsToChar(Current))
  ENDIF
  ! Tidy up
  CALL Delete(P)
  CALL Delete(PTilde)
  CALL Delete(F)
  CALL Delete(FTilde)
  CALL Delete(T1)
  CALL ShutDown(Prog)
END PROGRAM ODA
