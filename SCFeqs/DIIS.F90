PROGRAM PulayDIIS
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
   TYPE(DBL_RNK2)                 :: B,BOld,BInv,BTmp
   REAL(DOUBLE)                   :: DIISErr
   INTEGER                        :: I,J,K,N,ISCF,JSCF,KSCF
   CHARACTER(LEN=2)               :: Cycl,NxtC
   CHARACTER(LEN=9),PARAMETER     :: Prog='PulayDIIS'
   CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
   LOGICAL                        :: Present
!------------------------------------------------------------------ 
!
!
   CALL StartUp(Args,Prog)
   ISCF=Args%I%I(1)
   Cycl=IntToChar(ISCF)
   NxtC=IntToChar(ISCF+1)
!  Allocations
   CALL New(P)
   CALL New(F)
   CALL New(E)
   CALL New(Tmp1)
!------------------------------------------------
!  Create a new error vector E=[F_(i+1),P_i]
!
   CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix 
!
!   Turn off DIIS
!   CALL Put(F,TrixFile('F_DIIS',Args,0))    ! the orthogonalized Fock matrix 
!   CALL Put(One,'diiserr',Tag_O='_'//TRIM(CurGeom)//'_'//TRIM(CurBase)//'_'//TRIM(SCFCycl))
!   CALL ShutDown(Prog)
!
   CALL Get(P,TrixFile('OrthoD',Args,0))    ! the orthogonalized Density matrix
   CALL Multiply(F,P,E)  
   CALL Multiply(P,F,E,-One)
!  We didnt filter E for obvious reasons 
   CALL Put(E,TrixFile('E',Args,0))
!  The DIIS Error 
   DIISErr=SQRT(Dot(E,E))/DBLE(NBasF)
!----------------------------------------------
!  Build a new B matrix
!
   N=ISCF+1
   CALL New(B,(/N,N/))
!  Start with the last cycles B
   IF(ISCF>1)THEN
      CALL New(BOld,(/ISCF,ISCF/))
      CALL Get(BOld,'bmat_diis')
      B%D(1:ISCF,1:ISCF)=BOld%D(1:ISCF,1:ISCF)
      CALL Delete(BOld)
   ENDIF
!  Add new elements 
   DO JSCF=ISCF-1,1,-1
      KSCF=JSCF-ISCF
      CALL Get(Tmp1,TrixFile('E',Args,KSCF))
      B%D(ISCF,JSCF)=Dot(Tmp1,E)
      B%D(JSCF,ISCF)=B%D(ISCF,JSCF)
   ENDDO       
   B%D(ISCF,ISCF)=Dot(E,E)
   B%D(N,1:N)=One
   B%D(1:N,N)=One
   B%D(N,N)=Zero
!  Save for next cycle
   CALL Put(B,'bmat_diis',UnLimit_O=.TRUE.)
!---------------------------------------------------------------
!  Solve the least squares problem via eigensolution to obtain
!  new DIIS (mixing) coeficients.
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
     CALL FunkOnSqMat(N,Inverse,B%D,BInv%D)
     CALL UnSetDSYEVWork()
     CALL DGEMV('N',N,N,One,BInv%D,N,V%D,1,Zero,DIISCo%D,1)
     CALL Delete(V)
     CALL Delete(BInv)
#ifdef PARALLEL
   ENDIF
   CALL BCast(DIISCo)
#endif
!---------------------------------------------------------------
!  IO
! 
   CALL Put(DIISErr,'diiserr',Tag_O='_'//TRIM(CurGeom)//'_'//TRIM(CurBase)//'_'//TRIM(SCFCycl))
   IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
      CALL OpenASCII(OutFile,Out)
      CALL PrintProtectL(Out)
      Mssg=ProcessName(Prog)//' DIISErr = '//TRIM(DblToMedmChar(DIISErr))
      WRITE(*,*)TRIM(Mssg)
      WRITE(Out,*)TRIM(Mssg)
      Mssg=ProcessName(Prog)//' DIISCo = '
      DO I=1,ISCF-1
         Mssg=TRIM(Mssg)//' '//TRIM(DblToMedmChar(DIISCo%D(I)))//','
      ENDDO
      Mssg=TRIM(Mssg)//' '//TRIM(DblToMedmChar(DIISCo%D(ISCF)))
      WRITE(*,*)TRIM(Mssg)
      WRITE(Out,*)TRIM(Mssg)
      CALL PrintProtectR(Out)
      CLOSE(Out)
  ENDIF
!--------------------------------------------------------------------------------
!  Extrapolation with the DIIS coefficients
!
   CALL Multiply(F,DIISCo%D(ISCF))     
   DO JSCF=iSCF-1,1,-1 
      KSCF=JSCF-ISCF
      CALL Get(Tmp1,TrixFile('OrthoF',Args,KSCF))
      CALL Multiply(Tmp1,DIISCo%D(JSCF))
      CALL Add(F,Tmp1,E)
      IF(JSCF==1)THEN
!        Only filter the end product
         CALL Filter(F,E)
      ELSE
         CALL SetEq(F,E)
      ENDIF
   ENDDO
!--------------------------------------------------------------------
!  IO for the orthogonal, extrapolated F 
!
   CALL Put(F,TrixFile('F_DIIS',Args,0)) 
   CALL PChkSum(F,'F_DIIS['//TRIM(Cycl)//']',Prog)
   CALL PPrint(F,'F_DIIS['//TRIM(Cycl)//']')
   CALL Plot(F,'F_DIIS_'//TRIM(Cycl))
!------------------------------
!  Tidy up 
!
   CALL Delete(F)
   CALL Delete(DIISCo)
   CALL Delete(P)
   CALL Delete(E)
   CALL Delete(Tmp1)
!
   CALL ShutDown(Prog)   
!
END PROGRAM !PulayDIIS




