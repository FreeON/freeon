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
                                 :: F,P,EI,EJ,Tmp1,Tmp2
  TYPE(ARGMT)                    :: Args
  TYPE(INT_VECT)                 :: IWork,Idx,SCFOff,DIISInfo
  TYPE(DBL_VECT)                 :: V,DIISCo,AbsDIISCo
  TYPE(DBL_RNK2)                 :: B,BB,BInv,BTmp
  TYPE(CMPoles),DIMENSION(2)     :: MP
  REAL(DOUBLE),DIMENSION(3)      :: AvDP,DeltaDPi,DeltaDPj
  REAL(DOUBLE),DIMENSION(6)      :: AvQP,DeltaQPi,DeltaQPj
  REAL(DOUBLE)                   :: DIISErr,C0,C1,Damp,EigThresh, &
       RelDevDP,RelDevQP,Sellers,DMax
  INTEGER                        :: I,J,I0,J0,K,N,M,ISCF,BMax,DoDIIS,iDIIS,iOffSet,DIISBeg
  CHARACTER(LEN=2)               :: Cycl,NxtC
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: Mssg,FFile
  LOGICAL                        :: Present,Sloshed
  INTEGER                        :: IPresent,JPresent
  CHARACTER(LEN=4),PARAMETER     :: Prog='DIIS'
  CHARACTER(LEN=DCL) :: NAME
  !-------------------------------------------------------------------------------------
  !  Initial setup
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  ISCF=Args%I%I(1)
  !  Parse for DIIS options
  CALL OpenASCII(InpFile,Inp)  
  !  Threshold for projection of small eigenvalues
  IF(.NOT.OptDblQ(Inp,'DIISThresh',EigThresh))EigThresh=1.D-10
  !  Damping coefficient for first cycle
  IF(.NOT.OptDblQ(Inp,'DIISDamp',Damp))Damp=1D-1
  !  Dont allow damping below 0.001, as this can cause false convergence
  Damp=MAX(Damp,1D-3)
  !  Dont allow damping above 0.5, as this is silly (DIIS would certainly work better)
  Damp=MIN(Damp,5D-1)
  !  Max number of equations to keep in DIIS 
  IF(.NOT.OptIntQ(Inp,'DIISDimension',BMax)) BMax=8
  IF(BMax.GT.DIIS_MAX_MATRIX_SIZE)THEN
     CALL Warn('Requested DIISDimension '//TRIM(IntToChar(BMax))// &
          &    ' greater than DIIS_MAX_MATRIX_SIZE '//TRIM(IntToChar(DIIS_MAX_MATRIX_SIZE))// &
          &    '! '//RTRN//'Reseting DIISDimension to DIIS_MAX_MATRIX_SIZE.')
     BMax=DIIS_MAX_MATRIX_SIZE
  ENDIF
  CLOSE(Inp)
  !  Allocations
  CALL New(P)
  CALL New(F)
  CALL New(EI)
  CALL New(EJ)
  CALL New(Tmp1)
  CALL New(MP(1))
  CALL New(MP(2))
  !  The current DIIS error 
  CALL Get(F,TrixFile('OrthoF',Args,0))  
  CALL Get(P,TrixFile('OrthoD',Args,0))  
  CALL Multiply(F,P,EI)  
  CALL Multiply(P,F,EI,-One)
  DIISErr=SQRT(Dot(EI,EI))/DBLE(NBasF)
  ! Consider just damping, certainly on first go through
  IF(ISCF<=1)THEN
     DoDIIS=-1  ! No DIIS, but damp non-extrapolated Fock matrices
  ELSEIF(BMax/=0)THEN
     DoDIIS=1   ! We are doing DIIS, extrapolating non-extrapolated Fock matrices
  ELSEIF(BMax==0)THEN
     DoDIIS=0   ! We are purely damping, using previously extrapolated Fock matrices
  ENDIF
  !  Build the B matrix if on second SCF cycle (starting from 0)
  !  and if pure damping flag is not on (DIISDimension=0)
  IF(DoDIIS==1)THEN
     N=MIN(ISCF,BMax)+1
     M=MAX(1,ISCF-BMax+1)
     CALL New(B,(/N,N/))
     !write(*,*) 'M',M,' N',N
     CALL New(DIISInfo,2)
     CALL Get(DIISInfo,'diisinfo')
     CALL New(BTmp,(/DIIS_MAX_MATRIX_SIZE,DIIS_MAX_MATRIX_SIZE/))
     CALL Get(BTmp,'diismtrix')
     IF(DIISInfo%I(1).EQ.Args%I%I(2).AND.DIISInfo%I(2).EQ.Args%I%I(3)) THEN
        ! We didn't BS switch or new geom or oda, then build a part of the B matrix.
        B%D(1:N-2,1:N-2)=BTmp%D(1:N-2,1:N-2)
        DIISBeg=N-1
     ELSE
        ! We did BS switch or new geom or oda, then build the full B matrix.
        !write(*,*) 'We did BS switch or new geom.'
        DIISInfo%I(1)=Args%I%I(2)
        DIISInfo%I(2)=Args%I%I(3)
        DIISBeg=1
     ENDIF
     CALL Put(DIISInfo,'diisinfo')
     CALL Delete(DIISInfo)
     !
     ! Pulays most excellent B matrix     
     I0=M-ISCF+DIISBeg-1
     DO I=DIISBeg,N-1
        !write(*,*) 'I0',I0
        !write(*,*) TrixFile('OrthoF',Args,I0)
        IF(MyID.EQ.ROOT) THEN
           FFile=TrixFile('E_DIIS',Args,I0)
           INQUIRE(FILE=FFile,EXIST=Present)
           IPresent=1
           IF(Present)IPresent=0
        ENDIF
#ifdef PARALLEL
        CALL BCast(IPresent)
#endif
        IF(IPresent.EQ.0) THEN
           !write(*,*) 'We load I '//TrixFile('E_DIIS',Args,I0)
           CALL Get(EI,TrixFile('E_DIIS',Args,I0))
        ELSE
           CALL Get(F,TrixFile('OrthoF',Args,I0))    
           CALL Get(P,TrixFile('OrthoD',Args,I0))    
           CALL Multiply(F,P,EI)  
           CALL Multiply(P,F,EI,-One)
           CALL Put(EI,TrixFile('E_DIIS',Args,I0))
           !write(*,*) 'We save I '//TrixFile('E_DIIS',Args,I0)
        ENDIF
!       We dont filter E for obvious reasons 
        J0=I0-I+1
        DO J=1,I
           !write(*,*) 'J0',J0
           !write(*,*) TrixFile('OrthoF',Args,J0)
           IF(MyID.EQ.ROOT) THEN
              FFile=TrixFile('E_DIIS',Args,J0)
              INQUIRE(FILE=FFile,EXIST=Present)
              JPresent=1
              IF(Present)JPresent=0
           ENDIF
#ifdef PARALLEL
           CALL BCast(JPresent)
#endif
           IF(JPresent.EQ.0) then
              !write(*,*) 'We load J '//TrixFile('E_DIIS',Args,J0)
              CALL Get(EJ,TrixFile('E_DIIS',Args,J0))
           ELSE
              CALL Get(F,TrixFile('OrthoF',Args,J0))    
              CALL Get(P,TrixFile('OrthoD',Args,J0))    
              CALL Multiply(F,P,EJ)  
              CALL Multiply(P,F,EJ,-One)
              CALL Put(EJ,TrixFile('E_DIIS',Args,J0))
              !write(*,*) 'We save J '//TrixFile('E_DIIS',Args,J0)
           ENDIF
           B%D(I,J)=Dot(EI,EJ)
           B%D(J,I)=B%D(I,J)
           J0=J0+1
        ENDDO
        I0=I0+1
     ENDDO
!!$     ! Pulays most excellent B matrix
!!$     I0=M-ISCF
!!$     DO I=1,N-1
!!$        CALL Get(F,TrixFile('OrthoF',Args,I0))    
!!$        CALL Get(P,TrixFile('OrthoD',Args,I0))    
!!$        CALL Multiply(F,P,EI)  
!!$        CALL Multiply(P,F,EI,-One)
!!$!       We dont filter E for obvious reasons 
!!$        J0=I0
!!$        DO J=I,N-1
!!$           CALL Get(F,TrixFile('OrthoF',Args,J0))    
!!$           CALL Get(P,TrixFile('OrthoD',Args,J0))    
!!$           CALL Multiply(F,P,EJ)  
!!$           CALL Multiply(P,F,EJ,-One)
!!$           B%D(I,J)=Dot(EI,EJ)
!!$           B%D(J,I)=B%D(I,J)
!!$           J0=J0+1
!!$        ENDDO
!!$        I0=I0+1
!!$     ENDDO
     B%D(N,1:N)=One
     B%D(1:N,N)=One
     B%D(N,N)=Zero
     IF(N.LT.BMax+1) THEN
        ! We didn't reach the size of the matrix.
        BTmp%D(1:N-1,1:N-1)=B%D(1:N-1,1:N-1)
     ELSE
        ! We reach the size of the matrix, we reduce it.
        !write(*,*) 'We reduce it'
        BTmp%D(1:N-2,1:N-2)=B%D(2:N-1,2:N-1)
     ENDIF
     CALL Put(BTmp,'diismtrix')
     CALL Delete(BTmp)
     ! Solve the least squares problem to obtain new DIIS coeficients.
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
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                EigenThresh_O=EigThresh,PrintCond_O=.TRUE.,Prog_O=Prog)
        ELSE
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                EigenThresh_O=EigThresh)
        ENDIF
        CALL UnSetDSYEVWork()
        CALL DGEMV('N',N,N,One,BInv%D,N,V%D,1,Zero,DIISCo%D,1)
        CALL Delete(V)
        CALL Delete(BInv)
#ifdef PARALLEL
     ENDIF
     CALL BCast(DIISCo)
#endif
     Mssg=ProcessName(Prog,'Pulay C1')//'DIISCo = '
  ELSE
     Mssg=ProcessName(Prog,'Damping')//'Co = '
     N=3
     M=ISCF-N+2
     CALL New(DIISCo,2)
     ! Damping on the second cycle
     DIISCo%D(1)=One-Damp
     DIISCo%D(2)=Damp
  ENDIF
  !  IO
  CALL Put(DIISErr,'diiserr')
  IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
     DO I=1,N-2
        IF(MOD(I,4)==0)THEN
           Mssg=TRIM(Mssg)//RTRN//ProcessName() &
                //'          '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
        ELSE
           Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
        ENDIF
     ENDDO
     Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(DIISCo%D(N-1)))
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        IF(PrintFlags%Key==DEBUG_MAXIMUM) &
             WRITE(Out,*)TRIM(Mssg)
        CALL PrintProtectR(Out)
        CLOSE(Out)
#ifdef PARALLEL
     ENDIF
#endif
  ENDIF
  ! Allocate some indecies for re-ordering
  CALL New(Idx,N)
  CALL New(SCFOff,N)
  CALL New(AbsDIISCo,N)
  ! Reorder the DIIS, starting with smallest values and summing to the largest
  IOffSet=M-ISCF
  DO I=1,N-1
     Idx%I(I)=I
     SCFOff%I(I)=IOffSet
     AbsDIISCo%D(I)=ABS(DIISCo%D(I))
     IOffSet=IOffSet+1
  ENDDO
  CALL Sort(AbsDIISCo,Idx,N-1,1)
  ! Start with a matrix of diagonal zeros...
  CALL SetToI(F)
  CALL Multiply(F,Zero)
  ! And do the summation
  DO I=1,N-1
     CALL Get(Tmp1,TrixFile('OrthoF',Args,SCFOff%I(Idx%I(I))))
     CALL Multiply(Tmp1,DIISCo%D(Idx%I(I)))
     CALL Add(F,Tmp1,EI)
     CALL SetEq(F,EI)
  ENDDO
  ! IO for the orthogonal, extrapolated F 
  CALL Put(F,TrixFile('F_DIIS',Args,0)) 
  CALL PChkSum(F,'F_DIIS['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint(F,'F_DIIS['//TRIM(SCFCycl)//']')
  CALL Plot(F,'F_DIIS_'//TRIM(SCFCycl))
  ! Tidy up 
  CALL Delete(Idx)
  CALL Delete(SCFOff)
  CALL Delete(AbsDIISCo)
  CALL Delete(F)
  CALL Delete(DIISCo)
  CALL Delete(P)
  CALL Delete(EI)
  CALL Delete(EJ)
  CALL Delete(Tmp1)
  CALL ShutDown(Prog)   
END PROGRAM DIIS

! SCFStatus : [0,1,1]    :: <SCF> = 8.51823858, dD = 0.24D+01
! SCFStatus : [1,1,1]    :: <SCF> = -71.64537636, dD = 0.15D-01, DIIS = 0.15D+00
! SCFStatus : [2,1,1]    :: <SCF> = -71.69892926, dD = 0.20D+00, DIIS = 0.14D+00
! SCFStatus : [3,1,1]    :: <SCF> = -72.27825158, dD = 0.16D+00, DIIS = 0.70D-01
! SCFStatus : [4,1,1]    :: <SCF> = -72.41620275, dD = 0.27D-01, DIIS = 0.72D-02
! SCFStatus : [5,1,1]    :: <SCF> = -72.41719223, dD = 0.80D-02, DIIS = 0.36D-02
! SCFStatus : [6,1,1]    :: <SCF> = -72.41750179, dD = 0.94D-03, DIIS = 0.36D-03
! SCFStatus : [7,1,1]    :: <SCF> = -72.41750493, dD = 0.45D-04, DIIS = 0.23D-04



