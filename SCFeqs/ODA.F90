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

  !-------------------------------------------------------------------------------------
  TYPE(ARGMT)                    :: Args
  TYPE(BCSR)                     :: P,PTilde,F,FTilde,T1
  REAL(DOUBLE)                   :: e0,e1,e0p,e1p,a3,b3,c3,d3,EMns,EPls,EMin, &
       LMns,LPls,L,L1,ENucTot,ENucTotTilde,ExcTilde,Exc
  INTEGER                        :: I,M,iSCF
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,MatFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='ODA'
  !-------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  iSCF=Args%I%I(1)
  ! Suss for matrix threshold overide
  CALL SussTrix('ODA',Prog)  
  ! Allocations 
  CALL New(T1)
  CALL New(P)
  CALL New(F)
  CALL New(PTilde)
  CALL New(FTilde)
  ! Get the previous (N-M) tilde values where M is -1
  M=-1
  CALL Get(PTilde,TrixFile('OrthoD',Args,M))   
  CALL Get(FTilde,TrixFile('OrthoF',Args,M))   
  !CALL PChkSum(PTilde,'P0Tilde',Unit_O=6)
  !CALL PChkSum(FTilde,'F0Tilde',Unit_O=6)
  Current(1)=Current(1)+M
  CALL Get(e0,'Etot',StatsToChar(Current))
  Current(1)=Current(1)-M
  ! Get the current (N) non-tilde values
  CALL Get(P,TrixFile('OrthoD',Args,0))   
  CALL Get(F,TrixFile('OrthoF',Args,0))   
  ! T1 = P_N-PTilde_{N-1}
  CALL Multiply(PTilde,-One)
  CALL Add(PTilde,P,T1)  
  CALL Get(e1,'Etot',StatsToChar(Current))
  ! Get the previous Fock matrix F[PTilde_(N-1)]
  e0p=Two*Trace(FTilde,T1)
  e1p=Two*Trace(F,T1)
  !WRITE(*,*)' e0 = ',e0
  !WRITE(*,*)' e1 = ',e1
  !WRITE(*,*)'e0p = ',e0p
  !WRITE(*,*)'e1p = ',e1p
  ! Find the mixing parameter L from the
  ! cubic E3(L)=a3+b3*L+c3*L^2+d3*L^3
  a3=e0
  b3=e0p
  c3=-3D0*e0-2D0*e0p+3D0*e1-e1p
  d3=2D0*e0+e0p-2D0*e1+e1p
  IF(ABS(d3)<1D-4)THEN
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
        ! We must be done ...
        WRITE(*,*)' NEED SOME LOGIC HERE!!! '
        CALL Halt(' Bad logic in ODA ')
     ELSE
        L=One
        L1=Zero
        EMin=e1
     ENDIF
  ENDIF
  Mssg=' Mix = '//TRIM(FltToShrtChar(L))
  Mssg=ProcessName(Prog,TRIM(Mssg))
  Mssg=TRIM(Mssg)//" <SCF> ~ "//TRIM(FltToMedmChar(EMin))//', d3 = '//TRIM(DblToShrtChar(d3))
  WRITE(*,*)TRIM(Mssg)
  !  CALL Put(L,'ODAMixingParameter')
  ! PTilde_N=(1-L)*PTilde_(N-1)+L*P_N
  CALL Multiply(P,L)
  CALL Multiply(PTilde,-L1)
  CALL Add(P,PTilde,T1)
  ! Orthogonal put and xform to AO rep and put
  Args%I%I(1)=Args%I%I(1)-1
  Args%I%I(4)=Args%I%I(4)-1
  CALL PutXForm(Prog,Args,T1,P,PTilde)
  ! FTilde_N ~ (1-L)*FTilde_(N-1)+L*F_N
  Args%I%I(1)=Args%I%I(1)+1
  Args%I%I(4)=Args%I%I(4)+1
  CALL Multiply(F,L)
  CALL Multiply(FTilde,L1)
  CALL Add(F,FTilde,T1)
  CALL Put(T1,TrixFile('OrthoF',Args,0)) 
  !  CALL Put(T1,TrixFile('F_DIIS',Args,0)) 
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
