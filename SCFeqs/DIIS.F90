PROGRAM DIIS
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE Parse
   USE Macros
   USE SetXYZ
   USE LinAlg
   USE MatFunk
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
#ifdef PARALLEL
   TYPE(DBCSR) &
#else
   TYPE(BCSR)  & 
#endif
                                  :: F,P,E,Tmp1,Tmp2
   TYPE(ARGMT)                    :: Args
   TYPE(INT_VECT)                 :: IWork
   TYPE(DBL_VECT)                 :: V,DIISCo
   TYPE(DBL_RNK2)                 :: B,BB,BInv,BTmp
   TYPE(CMPoles),DIMENSION(2)     :: MP
   REAL(DOUBLE),DIMENSION(3)      :: AvDP,DeltaDPi,DeltaDPj
   REAL(DOUBLE),DIMENSION(6)      :: AvQP,DeltaQPi,DeltaQPj
   REAL(DOUBLE)                   :: DIISErr,C0,C1,Damp,EigThresh, &
                                     RelDevDP,RelDevQP,Sellers,DMax
   INTEGER                        :: I,J,I0,J0,K,N,M,ISCF,BMax
   CHARACTER(LEN=2)               :: Cycl,NxtC
   CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: Mssg
   LOGICAL                        :: Present,Sloshed
   CHARACTER(LEN=4),PARAMETER     :: Prog='DIIS'
!-------------------------------------------------------------------------------------
!  Initial setup
   CALL StartUp(Args,Prog)
   ISCF=Args%I%I(1)
!  Parse for DIIS options
   CALL OpenASCII(InpFile,Inp)  
!  Threshold for projection of small eigenvalues
   IF(.NOT.OptDblQ(Inp,'DIISThresh',EigThresh))EigThresh=1.D-10
!  Damping coefficient for first cycle
   IF(.NOT.OptDblQ(Inp,'DIISDamp',Damp))Damp=2D-1
!  Dont allow damping below 0.05, as this can cause false convergence
   Damp=MAX(Damp,5D-2)
!  Max number of equations to keep in DIIS 
   IF(.NOT.OptIntQ(Inp,'DIISDimension',BMax))BMax=15
   CLOSE(Inp)
!-------------------------------------------------------------------------------------
!  Allocations
   CALL New(P)
   CALL New(F)
   CALL New(E)
   CALL New(Tmp1)
   CALL New(MP(1))
   CALL New(MP(2))
!-------------------------------------------------------------------------------------
!  Create a new error vector E=[F_(i+1),P_i]
   CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonal Fock matrix 
   CALL Get(P,TrixFile('OrthoD',Args,0))    ! the orthogonal density matrix
   CALL Multiply(F,P,E)  
   CALL Multiply(P,F,E,-One)
!  We dont filter E for obvious reasons 
   CALL Put(E,TrixFile('E',Args,0))
!  The DIIS Error 
   DIISErr=SQRT(Dot(E,E))/DBLE(NBasF)
!  Build the B matrix if on second SCF cycle (starting from 0)
!  and if pure damping flag is not on (DIISDimension=0)
   IF(ISCF>1.AND.BMax/=0)THEN
      N=MIN(ISCF+1,BMax+1)
      M=MAX(1,ISCF-BMax+1)
      CALL New(B,(/N,N/))
!-------------------------------------------------------------------------------------
!     Check for charge sloshing
      CALL Get(DMax,'DMax',StatsToChar(Previous))
      IF(DMax>1.D-1.AND.ISCF>5)THEN
         Sloshed=.TRUE.
      ELSE
         Sloshed=.FALSE.
      ENDIF
!-------------------------------------------------------------------------------------
!     Build the B matrix, possibly with Sellers multipole modification
      I0=M-ISCF
      DO I=1,N-1
         IF(Sloshed)THEN
            CALL Get(MP(1),IntToChar(M+I-2))
            CALL Get(MP(2),IntToChar(M+I-1))
            DeltaDPi=MP(1)%DPole%D-MP(2)%DPole%D
            DeltaQPi=MP(1)%QPole%D-MP(2)%QPole%D
         ENDIF
         CALL Get(Tmp1,TrixFile('E',Args,I0))
         J0=I0
         DO J=I,N-1
            CALL Get(E,TrixFile('E',Args,J0))
            B%D(I,J)=Dot(Tmp1,E)
            IF(Sloshed)THEN
               CALL Get(MP(1),IntToChar(M+J-2))
               CALL Get(MP(2),IntToChar(M+J-1))
               DeltaDPj=MP(1)%DPole%D-MP(2)%DPole%D
               DeltaQPj=MP(1)%QPole%D-MP(2)%QPole%D
!              Add in Sellers anti charge sloshing terms
               Sellers=1D-1*DOT_PRODUCT(DeltaDPi,DeltaDPj) &
                      +1D-1*DOT_PRODUCT(DeltaQPi,DeltaQPj)
               B%D(I,J)=B%D(I,J)+Sellers
            ENDIF
            B%D(J,I)=B%D(I,J)
            J0=J0+1
         ENDDO
         I0=I0+1
      ENDDO
      B%D(N,1:N)=One
      B%D(1:N,N)=One
      B%D(N,N)=Zero
!-------------------------------------------------------------------------------------
!     Solve the least squares problem to obtain new DIIS coeficients.
      CALL New(DIISCo,N)
#ifdef PARALLEL
      IF(MyId==ROOT)THEN
#endif
         CALL New(BInv,(/N,N/))
         CALL New(V,N)
         V%D=Zero
         V%D(N)=One
         CALL SetDSYEVWork(N)
         BInv%D=Zero
         IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
            CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,EigenThresh_O=EigThresh, &
                             PrintCond_O=.TRUE.,Prog_O=Prog)
         ELSE
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,EigenThresh_O=EigThresh)
         ENDIF
         CALL UnSetDSYEVWork()
         CALL DGEMV('N',N,N,One,BInv%D,N,V%D,1,Zero,DIISCo%D,1)
         CALL Delete(V)
         CALL Delete(BInv)
#ifdef PARALLEL
      ENDIF
      CALL BCast(DIISCo)
#endif
      IF(Sloshed)THEN
         Mssg=ProcessName(Prog,'Multipole C1')//'DIISCo = '
      ELSE
         Mssg=ProcessName(Prog,'Pulay C1')//'DIISCo = '
      ENDIF
   ELSE
      Mssg=ProcessName(Prog,'Damping')//'Co = '
      N=3
      CALL New(DIISCo,2)
!     Damping on the second cycle
      DIISCo%D(1)=One-Damp
      DIISCo%D(2)=Damp
   ENDIF
!-------------------------------------------------------------------------------------
!  IO
   CALL Put(DIISErr,'diiserr',Tag_O='_'//TRIM(CurGeom)//'_'//TRIM(CurBase)//'_'//TRIM(SCFCycl))
   IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
      CALL OpenASCII(OutFile,Out)
      CALL PrintProtectL(Out)
      DO I=1,N-2
         IF(MOD(I,4)==0)THEN
            Mssg=TRIM(Mssg)//RTRN//ProcessName() &
               //'          '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
         ELSE
            Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
         ENDIF
      ENDDO
      Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(DIISCo%D(N-1)))
      IF(PrintFlags%Key==DEBUG_MAXIMUM) &
      WRITE(*,*)TRIM(Mssg)
      WRITE(Out,*)TRIM(Mssg)
      CALL PrintProtectR(Out)
      CLOSE(Out)
  ENDIF
!-------------------------------------------------------------------------------------
! Extrapolation with DIIS coefficients
  CALL Multiply(F,DIISCo%D(N-1))     
  I0=N-2
  DO I=ISCF-1,M,-1
     CALL Get(Tmp1,TrixFile('OrthoF',Args,I-ISCF))
     CALL Multiply(Tmp1,DIISCo%D(I0))
     CALL Add(F,Tmp1,E)
     IF(I==1)THEN
!       Only filter the end product
        CALL Filter(F,E)
     ELSE
        CALL SetEq(F,E)
     ENDIF
     I0=I0-1
  ENDDO
!-------------------------------------------------------------------------------------
!  IO for the orthogonal, extrapolated F 
   CALL Put(F,TrixFile('F_DIIS',Args,0)) 
   CALL PChkSum(F,'F_DIIS['//TRIM(SCFCycl)//']',Prog)
   CALL PPrint(F,'F_DIIS['//TRIM(SCFCycl)//']')
   CALL Plot(F,'F_DIIS_'//TRIM(SCFCycl))
!-------------------------------------------------------------------------------------
!  Tidy up 
   CALL Delete(F)
   CALL Delete(DIISCo)
   CALL Delete(P)
   CALL Delete(E)
   CALL Delete(Tmp1)
   CALL ShutDown(Prog)   
END PROGRAM DIIS




