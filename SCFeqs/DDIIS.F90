PROGRAM DDIIS
!H=================================================================================
!H PROGRAM DDIIS
!H
!H  OPTIONS:
!H  DEBUGING: Use -DDDIIS_DBUG to print some stuff.
!H  INFO    : Use -DDDIIS_INFO to print some stuff.
!H
!H Comment:
!H  -MUST CHECK THE ARGS FOR GROUND STATS DENSITY AND FOCK MATRICES.
!H  -MUST CHECK FOR -ISCF- AND -SCFCycl-.
!H
!H Ref:
!H  V. Weber, C. Daul Chem. Phys. Let. 370, 99-105, 2003.
!H
!H=================================================================================
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
  !
  IMPLICIT NONE
  !-------------------------------------------------------------------
#ifdef PARALLEL
  TYPE(DBCSR)                      :: F,P,FPrim,PPrim,EPrim,Tmp1
#else
  TYPE(BCSR )                      :: F,P,FPrim,PPrim,EPrim,Tmp1
#endif
  TYPE(ARGMT)                      :: Args
  TYPE(INT_VECT)                   :: Idx,SCFOff
  TYPE(DBL_VECT)                   :: V,DIISCo,AbsDIISCo
  TYPE(DBL_RNK2)                   :: B,BInv
  !-------------------------------------------------------------------
  REAL(DOUBLE)                     :: DIISErr,Damp,EigThresh
  INTEGER                          :: I,J,I0,J0,N,M,BMax,DoDIIS,iOffSet
  INTEGER                          :: LastSCFCycle,CPSCFCycl,DDIISStart
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=*), PARAMETER      :: Prog='DDIIS'
  LOGICAL                          :: IsPresent
  !-------------------------------------------------------------------
  !
  !  Initial setup
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  CPSCFCycl=Args%I%I(1)
  !
  ! Get Last SCF cycle.
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  ! Parse for DDIIS options.
  CALL OpenASCII(InpFile,Inp)  
  ! Threshold for projection of small eigenvalues.
  IF(.NOT.OptDblQ(Inp,'DDIISThresh',EigThresh))EigThresh=1.D-10
  ! Damping coefficient for first cycle
  IF(.NOT.OptDblQ(Inp,'DDIISDamp',Damp))Damp=1D-1
  ! Dont allow damping below 0.001, as this can cause false convergence.
  Damp=MAX(Damp,1D-3)
  ! Dont allow damping above 0.5, as this is silly (DIIS would certainly work better).
  Damp=MIN(Damp,5D-1)
  ! Max number of equations to keep in DIIS.
  IF(.NOT.OptIntQ(Inp,'DDIISDimension',BMax))BMax=15
  !IF(.NOT.OptIntQ(Inp,'DDIISStart',DDIISStart)) DDIISStart=1
  DDIISStart=1
  CLOSE(Inp)
  !
  ! Crate or zero the DDIISQueue.
  !  CALL New(DDIISQueue,(/0,BMax/))
  !  CALL SetEq(DDIISQueue%I,0)
  !  IF(CPSCFCycl.EQ.0) THEN
  !     DDIISQueue%I(0)=1
  !     CALL Put(DDIISQueue,'DDIISQueue')
  !  ELSE
  !     CALL Get(DDIISQueue,'DDIISQueue')
  !  ENDIF
  !
  ! Allocations.
  CALL New(P    )
  CALL New(F    )
  CALL New(Tmp1 )
  CALL New(EPrim)
  CALL New(PPrim)
  CALL New(FPrim)
  !
  ! Load the orthogonal Derivative Fock matrix.
  CALL Get(FPrim,TrixFile('OrthoFPrime'//TRIM(Args%C%C(4)),Args,0))
  !
  ! Load the orthogonal Derivative Density matrix.
  CALL Get(PPrim,TrixFile('OrthoDPrime'//TRIM(Args%C%C(4)),Args,0))
  !
  ! Load the orthogonal GS Fock matrix.
  CALL Get(F,TrixFile('OrthoF',Args,LastSCFCycle-Args%I%I(1)))
  !
  ! Load the orthogonal GS Density matrix.
  CALL Get(P,TrixFile('OrthoD',Args,LastSCFCycle-Args%I%I(1)))
  !
  ! Create a new error vector E'=[F_(i+1),P_i]'
  CALL Multiply(PPrim,F,EPrim     )   !E'=P'F
  CALL Multiply(P,FPrim,EPrim, One)   !E'=PF'+P'F
  CALL Multiply(FPrim,P,EPrim,-One)   !E'=F'P-(PF'+P'F)
  CALL Multiply(F,PPrim,EPrim, One)   !E'=FP'+F'P-(PF'+P'F)
  !
  ! We dont filter E' for obvious reasons .
  CALL Put(EPrim,TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,0))
  !
  ! Compute the DIIS error.
  DIISErr=SQRT(Dot(EPrim,EPrim))/DBLE(NBasF)
  !
  !
  IF(CPSCFCycl<=DDIISStart)THEN
     DoDIIS=-1  ! No DIIS, but damp non-extrapolated Fock matrices
  ELSEIF(BMax/=0)THEN
     DoDIIS=1   ! We are doing DIIS, extrapolating non-extrapolated Fock matrices
  ELSEIF(BMax==0)THEN
     DoDIIS=0   ! We are purely damping, using previously extrapolated Fock matrices
  ENDIF
  !  Build the B matrix if on second SCF cycle (starting from 0)
  !  and if pure damping flag is not on (DIISDimension=0)
  IF(DoDIIS==1)THEN
     N=MIN(CPSCFCycl+1,BMax+1)
     M=MAX(1,CPSCFCycl-BMax+1)
     CALL New(B,(/N,N/))
     !
     ! Build the B matrix.
     I0=M-CPSCFCycl
     DO I=1,N-1
        CALL Get(Tmp1,TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,I0))
        J0=I0
        DO J=I,N-1
           CALL Get(EPrim,TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,J0))
           B%D(I,J)=Dot(Tmp1,EPrim)
           B%D(J,I)=B%D(I,J)
           J0=J0+1
        ENDDO
        I0=I0+1
     ENDDO
     B%D(N,1:N)=One
     B%D(1:N,N)=One
     B%D(N,N)=Zero
     !
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
                &           EigenThresh_O=EigThresh,PrintCond_O=.TRUE.,Prog_O=Prog)
        ELSE
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                &           EigenThresh_O=EigThresh)
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
     M=CPSCFCycl-N+2
     CALL New(DIISCo,2)
     !
     ! Damping on the second cycle
     DIISCo%D(1)=One-Damp
     DIISCo%D(2)=Damp
  ENDIF
  !
  ! IO.
  CALL Put(DIISErr,'ddiiserr')
  !
  ! Printing.
  IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
     DO I=1,N-2
        IF(MOD(I,4)==0)THEN
           Mssg=TRIM(Mssg)//RTRN//ProcessName() &
                & //'          '//TRIM(DblToShrtChar(DIISCo%D(I)))//','
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
        IF(PrintFlags%Key==DEBUG_MAXIMUM) WRITE(Out,*)TRIM(Mssg)
        CALL PrintProtectR(Out)
        CLOSE(Out)
#ifdef PARALLEL
     ENDIF
#endif
  ENDIF
  !
  ! Allocate some indecies for re-ordering
  CALL New(Idx,N)
  CALL New(SCFOff,N)
  CALL New(AbsDIISCo,N)
  !
  ! Reorder the DIIS, starting with smallest values and summing to the largest
  IOffSet=M-CPSCFCycl
  DO I=1,N-1
     Idx%I(I)=I
     SCFOff%I(I)=IOffSet
     AbsDIISCo%D(I)=ABS(DIISCo%D(I))
     IOffSet=IOffSet+1
  ENDDO
  CALL Sort(AbsDIISCo,Idx,N-1,1)
  !
  ! Start with a matrix of diagonal zeros...
  CALL SetToI(FPrim)
  CALL Multiply(FPrim,Zero)
  !
  ! And do the summation
  DO I=1,N-1
     CALL Get(Tmp1,TrixFile('OrthoFPrime'//TRIM(Args%C%C(4)),Args,SCFOff%I(Idx%I(I))))
     CALL Multiply(Tmp1,DIISCo%D(Idx%I(I)))
     CALL Add(FPrim,Tmp1,EPrim)
     CALL SetEq(FPrim,EPrim)
  ENDDO
  !
  ! IO for the orthogonal, extrapolated FPrim 
  CALL Put(FPrim,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0)) 
  CALL PChkSum(FPrim,'FPrime_DDIIS'//TRIM(Args%C%C(4))//'['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( FPrim,'FPrime_DDIIS'//TRIM(Args%C%C(4))//'['//TRIM(SCFCycl)//']')
  CALL Plot(   FPrim,'FPrime_DDIIS'//TRIM(Args%C%C(4))//'_'//TRIM(SCFCycl))
  !
  ! Tidy up 
  CALL Delete(F        )
  CALL Delete(P        )
  CALL Delete(Idx      )
  CALL Delete(Tmp1     )
  CALL Delete(FPrim    )
  CALL Delete(PPrim    )
  CALL Delete(EPrim    )
  CALL Delete(SCFOff   )
  CALL Delete(DIISCo   )
  CALL Delete(AbsDIISCo)
  !
  CALL ShutDown(Prog)   
END PROGRAM DDIIS




