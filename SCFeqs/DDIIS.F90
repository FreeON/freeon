PROGRAM DDIIS
!H=================================================================================
!H PROGRAM DDIIS
!H
!H  OPTIONS:
!H  DEBUGING: Use -DDDIIS_DBUG to print some stuff.
!H  INFO    : Use -DDDIIS_INFO to print some stuff.
!H
!H Comment:
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
  REAL(DOUBLE)                     :: CondA,DDIISMaxCond
  INTEGER                          :: I,J,I0,J0,N,M,BMax,DoDIIS,iOffSet
  INTEGER                          :: LastSCFCycle,CPSCFCycl,DDIISStart
  INTEGER                          :: DDIISBeg,DDIISEnd,DDIISCurDim
  CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=*), PARAMETER      :: Prog='DDIIS'
  LOGICAL                          :: IsPresent
  !-------------------------------------------------------------------
  !
  ! Initial setup.
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  CPSCFCycl=Args%I%I(1)
  !
  ! Get Last SCF cycle.
  CALL Get(LastSCFCycle,'lastscfcycle')
  !
  !-------------------------------------------------------------------
  ! Parse for DDIIS options.
  !-------------------------------------------------------------------
  !
  ! Open input
  CALL OpenASCII(InpFile,Inp)  
  !
  ! Threshold for projection of small eigenvalues.
  IF(.NOT.OptDblQ(Inp,'DDIISThresh',EigThresh)) EigThresh=1.0D-10
  IF(EigThresh.LT.Zero) &
       & CALL Halt(' DDIISThresh cannot be smaller than zero, DDIISThresh=' &
       &           //TRIM(DblToShrtChar(EigThresh))//'.')
  !
  ! Damping coefficient for first cycle
  IF(.NOT.OptDblQ(Inp,'DDIISDamp',Damp)) Damp=1.0D-1
  IF(Damp.LT.Zero) &
       & CALL Halt(' DDIISDamp cannot be smaller than zero, DDIISDamp=' &
       &           //TRIM(DblToShrtChar(Damp))//'.')
  ! Dont allow damping below 0.001, as this can cause false convergence.
  Damp=MAX(Damp,1D-3)
  ! Dont allow damping above 0.5, as this is silly (DIIS would certainly work better).
  Damp=MIN(Damp,5D-1)
  !
  ! Max number of equations to keep in DIIS.
  IF(.NOT.OptIntQ(Inp,'DDIISDimension',BMax)) BMax=15
  IF(BMax.LT.0) &
       & CALL Halt(' DDIISDimension cannot be smaller than zero, DDIISDimension=' &
       &           //TRIM(IntToChar(BMax))//'.')
  !
  ! Condition number of B.
  IF(.NOT.OptDblQ(Inp,'DDIISMaxCond',DDIISMaxCond)) DDIISMaxCond=1.0D+8
  IF(DDIISMaxCond.LT.Zero) &
       & CALL Halt(' DDIISMaxCond cannot be smaller than zero, DDIISMaxCond=' &
       &           //TRIM(DblToShrtChar(DDIISMaxCond))//'.')
  !
  ! The iteration where we want to start the DDIIS.
  IF(.NOT.OptIntQ(Inp,'DDIISStart',DDIISStart)) DDIISStart=1
  IF(DDIISStart.LT.1) &
       & CALL Halt(' DDIISStart cannot be smaller than one, DDIISStart=' &
       &           //TRIM(IntToChar(DDIISStart))//'.')
  !
  ! Close the input.
  CLOSE(Inp)
  !
  !-------------------------------------------------------------------
  ! Do we need to do the DDIIS?
  !-------------------------------------------------------------------
  !
  IF(CPSCFCycl<=DDIISStart)THEN
     ! No DIIS yet, but damp non-extrapolated Fock matrices
     DoDIIS=-1
  ELSEIF(BMax/=0)THEN
     ! We are doing DIIS, extrapolating non-extrapolated Fock matrices
     DoDIIS=1   
  ELSEIF(BMax==0)THEN
     ! We are purely damping, using previously extrapolated Fock matrices
     DoDIIS=0
  ENDIF
  !
  !-------------------------------------------------------------------
  ! Initialize or get Beg and End DDIIS variables.
  !-------------------------------------------------------------------
  !
  !CALL New(BTmp,(/BMax+1,BMax+1/))
  IF(CPSCFCycl.LE.1) THEN
     DDIISBeg=1
     DDIISEnd=1
     CALL Put(DDIISBeg,'DDIISBeg')
     CALL Put(DDIISEnd,'DDIISEnd')
     !CALL SetEq(BTmp,Zero)
     !CALL Put(BTmp,'DDIISBMtrix')
  ELSE
     CALL Get(DDIISBeg,'DDIISBeg')
     CALL Get(DDIISEnd,'DDIISEnd')
     !CALL Get(BTmp,'DDIISBMtrix')
  ENDIF
  !
  !-------------------------------------------------------------------
  ! Allocations.
  !-------------------------------------------------------------------
  !
  CALL New(P    )
  CALL New(F    )
  CALL New(Tmp1 )
  CALL New(EPrim)
  CALL New(PPrim)
  CALL New(FPrim)
  !
  !-------------------------------------------------------------------
  ! Loading matrices.
  !-------------------------------------------------------------------
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
  !-------------------------------------------------------------------
  ! Build up new DDIIS error.
  !-------------------------------------------------------------------
  !
  ! Create a new error vector E'=[F_(i+1),P_i]'  (Could be done in a better way!)
  CALL Multiply(PPrim,F,EPrim     )  !E'=P'F
  CALL Multiply(P,FPrim,EPrim, One)  !E'=PF'+P'F
  CALL Multiply(FPrim,P,EPrim,-One)  !E'=F'P-(PF'+P'F)
  CALL Multiply(F,PPrim,EPrim, One)  !E'=FP'+F'P-(PF'+P'F)
  !
  ! We dont filter E' for obvious reasons .
  CALL Put(EPrim,TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,0))
  !
#ifdef DDIIS_DBUG
  WRITE(*,*) 'Save E''=<'//TRIM(TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,0))//'>'
#endif
  !
  ! Compute the DDIIS error.
  DIISErr=SQRT(Dot(EPrim,EPrim))/DBLE(NBasF)
  !
  ! IO save DDIIS error.
  CALL Put(DIISErr,'ddiiserr')
  !
  !-------------------------------------------------------------------
  ! Build the B matrix and solve the linear problem..
  !-------------------------------------------------------------------
  !
  ! Build the B matrix if on second SCF cycle (starting from 0)
  ! and if pure damping flag is not on (DIISDimension=0)
  SELECT CASE(DoDIIS)
  CASE(1)
     ! Do DDIIS
     !
     DDIISCurDim=DDIISEnd-DDIISBeg+1
     N=DDIISCurDim+1
     M=DDIISBeg
     CALL New(B,(/N,N/))
     !
     ! Build the B matrix.
     I0=DDIISBeg-CPSCFCycl
     DO I=1,N-1
#ifdef DDIIS_DEBUG
        WRITE(*,*) 'Load <ei|=<'//TRIM(TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,I0))//'>'
#endif
        CALL Get(Tmp1,TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,I0))
        J0=I0
        DO J=I,N-1
#ifdef DDIIS_DEBUG
           WRITE(*,*) 'Load |ej>=<'//TRIM(TrixFile('EPrime'//TRIM(Args%C%C(4)),Args,J0))//'>'
#endif
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
     CALL SetEq(DIISCo,Zero)
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
     !
     CALL New(BInv,(/N,N/))
     CALL New(V,N)
     !
     ! Solve the linear system and remove 
     ! linear dependancies if needed.
     DO
        V%D=Zero
        V%D(N)=One
        CALL SetDSYEVWork(N)
        BInv%D=Zero
        IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                &           EigenThresh_O=EigThresh,PrintCond_O=.TRUE.,Prog_O=Prog, &
                &           CoNo_O=CondA)
        ELSE
           CALL FunkOnSqMat(N,Inverse,B%D,BInv%D,PosDefMat_O=.FALSE., &
                &           EigenThresh_O=EigThresh, &
                &           CoNo_O=CondA)
        ENDIF
        CALL UnSetDSYEVWork()
        CALL DGEMV('N',N,N,One,BInv%D,N,V%D,1,Zero,DIISCo%D,1)
        !
!#ifdef DDIIS_DBUG
!        CALL PrintMatrix(B%D,N,N,2)
        WRITE(*,*) 'CondA=',CondA
!#endif
        !
        ! Do we need to reduce the size?
        IF(.NOT.(CondA.GT.DDIISMaxCond.AND.DDIISCurDim.GT.1)) EXIT
        !
        ! We need to remove the oldest element.
        DDIISBeg=DDIISBeg+1
        !
        ! New current dimension.
        DDIISCurDim=DDIISEnd-DDIISBeg+1
        !
        ! Set the new size.
        N=DDIISCurDim+1
        M=DDIISBeg
        !
        ! Copy temporarly the B matrix.
        BInv%D=B%D
        CALL Delete(B)
        CALL New(B,(/N,N/))
        B%D=Zero
        !
        ! Copy back the B matrix.
        B%D(1:N,1:N)=BInv%D(2:N+1,2:N+1)
        !
        ! Reallocate and set old arrays.
        CALL Delete(DIISCo)
        CALL New(DIISCo,N)
        DIISCo%D=Zero
        CALL Delete(BInv)
        CALL New(BInv,(/N,N/))
        BInv%D=Zero
     ENDDO
     !
     CALL Delete(V)
     CALL Delete(BInv)
     !
#ifdef PARALLEL
     ENDIF
     CALL BCast(DIISCo)
#endif
     Mssg=ProcessName(Prog,'Pulay C1')//'DIISCo = '
  CASE DEFAULT
     ! Do Damping.
     !
     Mssg=ProcessName(Prog,'Damping')//'Co = '
     N=3
     M=CPSCFCycl-N+2
     CALL New(DIISCo,2)
     !
     ! Damping on the second cycle
     DIISCo%D(1)=One-Damp
     DIISCo%D(2)=Damp
  END SELECT
  !
  !-------------------------------------------------------------------
  ! Printing.
  !-------------------------------------------------------------------
  !
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
  !-------------------------------------------------------------------
  ! Reordering.
  !-------------------------------------------------------------------
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
  !-------------------------------------------------------------------
  ! Increment DDIIS.
  !-------------------------------------------------------------------
  !
  DDIISEnd=DDIISEnd+1
  DDIISCurDim=DDIISEnd-DDIISBeg+1
  IF(DDIISCurDim.GT.BMax) DDIISBeg=DDIISBeg+1
  !
  ! Put in HDF Beg and End DDIIS variables.
  CALL Put(DDIISBeg,'DDIISBeg')
  CALL Put(DDIISEnd,'DDIISEnd')
  !
  !-------------------------------------------------------------------
  ! IO for the orthogonal, extrapolated FPrim 
  !-------------------------------------------------------------------
  !
  CALL Put(FPrim,TrixFile('FPrime_DDIIS'//TRIM(Args%C%C(4)),Args,0)) 
  CALL PChkSum(FPrim,'FPrime_DDIIS'//TRIM(Args%C%C(4))//'['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( FPrim,'FPrime_DDIIS'//TRIM(Args%C%C(4))//'['//TRIM(SCFCycl)//']')
  CALL Plot(   FPrim,'FPrime_DDIIS'//TRIM(Args%C%C(4))//'_'//TRIM(SCFCycl))
  !
  !-------------------------------------------------------------------
  ! Tidy up 
  !-------------------------------------------------------------------
  !
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
  !
CONTAINS
  !
  SUBROUTINE PrintMatrix(V,M,N,IOpt,IOut_O,SHFTM_O,SHFTN_O,TEXT_O)
    IMPLICIT NONE
!    ++++ PRINT OUT A RECTANGULAR MATRIX +++++
!    
!    V    : MATRIX MxN
!    M    : NUMBER OF ROW
!    N    : NUMBER OF COLUMN
!
!    IOpt = 0; MAX COLUMN = 10 
!    IOpt = 1; MAX COLUMN =  7
!    IOpt = 2; MAX COLUMN =  5
!    IOpt = 3; MAX COLUMN =  3 
!
    INTEGER                                  :: IMax,IMin,NMAX,I,J
    INTEGER                                  :: IOut,SHFTM,SHFTN
    INTEGER                     , INTENT(IN) :: M,N,IOpt
    INTEGER         , OPTIONAL  , INTENT(IN) :: IOut_O,SHFTM_O,SHFTN_O
    REAL(DOUBLE), DIMENSION(:,:), INTENT(IN) :: V
    CHARACTER(LEN=*), OPTIONAL  , INTENT(IN) :: TEXT_O
!
    IOut = 6
    IF(PRESENT(IOut_O)) IOut = IOut_O
!
    SHFTM = 0
    SHFTN = 0
    IF(PRESENT(SHFTM_O)) SHFTM = SHFTM_O
    IF(PRESENT(SHFTN_O)) SHFTN = SHFTN_O
!
    WRITE(IOut,100)
!
    SELECT CASE(IOPT)
    CASE(0); NMAX = 10
    CASE(1); NMAX =  7
    CASE(2); NMAX =  5
    CASE(3); NMAX =  3
    CASE DEFAULT
       ! WRITE ERROR MESSAGE OR PUT SOMETHING THERE
       RETURN
    END SELECT
!
    IF(N.EQ.0) RETURN
! 
    IF(PRESENT(TEXT_O)) WRITE(IOUT,*) TEXT_O
    IMAX = 0
!    
    DO WHILE (IMAX.LT.M)
       IMIN = IMAX+1
       IMAX = IMAX+NMAX
       IF(IMAX .GT. M) IMAX = M
       SELECT CASE(IOPT)
       CASE(0); WRITE(IOUT,1000) (I+SHFTN,I=IMIN,IMAX)
       CASE(1); WRITE(IOUT,2000) (I+SHFTN,I=IMIN,IMAX)
       CASE(2); WRITE(IOUT,3000) (I+SHFTN,I=IMIN,IMAX)
       CASE(3); WRITE(IOUT,4000) (I+SHFTN,I=IMIN,IMAX)
       END SELECT
       DO J = 1,N
          SELECT CASE(IOPT)
          CASE(0); WRITE(IOUT,1100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(1); WRITE(IOUT,2100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(2); WRITE(IOUT,3100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          CASE(3); WRITE(IOUT,4100) J+SHFTM,(V(J,I),I=IMIN,IMAX)
          END SELECT
       ENDDO
       WRITE (IOUT,100)
    ENDDO
    !
    WRITE (IOUT,100)
    !    
100 FORMAT('')
1000 FORMAT(6X,10(4X,I3,4X))
1100 FORMAT(I5,1X,10F11.5  )
2000 FORMAT(6X,7(6X,I3,6X ))
2100 FORMAT(I5,1X,7F15.10  )
3000 FORMAT(6X,7(7X,I3,6X ))
3100 FORMAT(I5,1X,7E16.8   )
4000 FORMAT(6X,7(7X,I3,6X ))
4100 FORMAT(I5,1X,7E16.8   )
  END SUBROUTINE PrintMatrix


END PROGRAM DDIIS




