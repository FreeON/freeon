! Authors: C. J. Tymczak and Anders M Niklasson
!---------------------------------------------------------------------
PROGRAM SDMM2
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
  IMPLICIT NONE
  TYPE(BCSR)                     :: P,X1,X2,X3,T1,T2,T3,F,Z
  TYPE(DBL_RNK2)                 :: B
  TYPE(ARGMT)                    :: Args
!
  REAL(DOUBLE)                   :: EigH,EigL,Con1,Con2,TrF,TrFX,TrXF,Error,E_old,E_new, &
                                    Error_Trace,Error_Energy,TrX1,TrX2,TrX3,Num,Den,Mu,  &
                                    Fact,TrP,Max_Error_FX,Tol_Energy,NumP,SUM
  INTEGER                        :: I,J,ITYPE,PNon0
  LOGICAL                        :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,FFile
  CHARACTER(LEN=5),PARAMETER     :: Prog='SDMM2'
!------------------------------------------------------------------ 
! Start
  CALL StartUp(Args,Prog)
!-------------------------------------------------------------------------
  CALL New(F)
!-----------------------------------------------------------------------------
! Get The Fock Matrix
!
  FFile=TrixFile('F_DIIS',Args,0)
  INQUIRE(FILE=FFile,EXIST=Present)
  IF(Present)THEN
     CALL Get(F,FFile)
  ELSE
     CALL Get(F,TrixFile('OrthoF',Args,0))    ! the orthogonalized Fock matrix
  ENDIF
!=============================================================================
!
!  Niklasson-Tymczak Solver
!
!=============================================================================
!
!  Set up the Density Matrix from the Fock Matrix
!
  Tol_Energy = Thresholds%Trix
  NumP = 0.5D0*DBLE(NEl)
!
  CALL New(X1)
  CALL New(X2)
  CALL New(X3)
  CALL New(T1)
  CALL New(T2)
  CALL New(T3)
! Estimate Eigenvalues
! Lowest
  EigL = 1.0D10
  CALL New(B,(/NBasF,NBasF/))
  CALL SetEq(B,F)
  DO I = 1,NBasF
     SUM = Zero
     DO J = 1,NBasF
        SUM = SUM + B%D(I,J)
     ENDDO
     IF(ABS(EigL) > ABS(SUM)) EigL = SUM
  ENDDO
! Higest
  EigH = Zero
  DO I = 1,NBasF
     SUM = Zero
     DO J = 1,NBasF
        SUM = SUM + B%D(I,J)*((-1)**J)
     ENDDO
     EigH = EigH+ABS(SUM)
  ENDDO
  EigH = EigH/DBLE(NBasF)
! X1
  Con1 = Two/(EigH-EigL)
  Con2 = One-Two*EigH/(EigH-EigL)
  CALL SetEq(X1,F)
  CALL Multiply(X1,Con1)
  CALL Add(X1,Con2)
! Energy
  CALL Multiply(F,X1,T1)
  TrF   = Trace(F)
  TrFX  = Trace(T1)
  E_old = Zero
  E_new = 0.5D0*(TrF-TrFX)
!
! Iterate to Convergencve
!
  DO I=1,40
!    X2
     CALL Multiply(X1,X1,T1)
     CALL Filter(X2,T1)
!    X3
     CALL Multiply(X2,X1,T1)
     CALL Filter(X3,T1)
!    Traces
     TrX1 = Trace(X1)
     TrX2 = Trace(X2)
     TrX3 = Trace(X3)
!    Error in Trace(X1) and Energy
     Error_Trace  = ABS(DBLE(NBasF)-TrX1-Two*NumP)
     Error_Energy = ABS(E_old-E_new)
!
     Num = DBLE(NBasF)-Two*NumP-(1.5D0*TrX1-0.5D0*TrX3)
     Den = DBLE(NBasF)-TrX2
     IF(ABS(Den) < 1.D-12) THEN
        IF(Den > Zero) THEN
           Den = 1.D-12
        ELSE
           Den = -1.D-12
        ENDIF
     ENDIF
     Mu  = num/den 
!
     CALL SetEq(T1,X1)
     CALL SetEq(T2,X2)
     CALL SetEq(T3,X3)
     IF(Mu > 0.5) THEN
        ITYPE = 1
        CALL Multiply(T1, 0.75D0)
        CALL Multiply(T2,-0.75D0)
        CALL Multiply(T3, 0.25D0)
        CALL Add(T1,T2,X3)                      ! 0.75*X - 0.75*X2
        CALL Add(X3,T3,T1)                      ! 0.75*X - 0.75*X2 + 0.25*X3
        CALL Add(T1, 0.75D0)                    ! 0.75*X - 0.75*X2 + 0.25*X3 + 0.75*I 
        CALL Filter(X1,T1)
     ELSEIF(Mu < -0.5) THEN
        ITYPE = 2
        CALL Multiply(T1, 0.75D0)
        CALL Multiply(T2, 0.75D0)
        CALL Multiply(T3, 0.25D0)
        CALL Add(T1,T2,X3)                      ! 0.75*X + 0.75*X2
        CALL Add(X3,T3,T1)                      ! 0.75*X + 0.75*X2 + 0.25*X3
        CALL Add(T1,-0.75D0)                    ! 0.75*X + 0.75*X2 + 0.25*X3 - 0.75*I 
        CALL Filter(X1,T1)
     ELSE
        ITYPE = 3
        CALL Multiply(T1, 1.5D0)
        CALL Multiply(T2,-Mu)
        CALL Multiply(T3,-0.50D0)
        CALL Add(T1,T3,X3)                      ! 1.5*X1 - 0.5*X3
        CALL Add(X3,T2,T1)                      ! 1.5*X1 - 0.5*X3 - Mu*X2
        CALL Add(T1,Mu)                         ! 1.5*X1 - 0.5*X3 + Mu*(I-X2)
        CALL Filter(X1,T1)
     ENDIF
!  
     CALL Multiply(F,X1,T1)
     TrFX = Trace(T1)
     E_old = E_new
     E_new = 0.5D0*(TrF-TrFX)
!        
     CALL Multiply(F,X1,T1)
     CALL Multiply(X1,F,T2)
     CALL Multiply(T2,-1.0D0)
     CALL Add (T1,T2,T3)
     Max_Error_FX = 0.5D0*MAX(T3)
!
     PNon0 = 100.D0*DBLE(X1%NNon0)/(DBLE(NBasF)**2)
!
     IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
        IF(ITYPE==1) THEN
           Mssg = TRIM(Prog) // ':   Mu < -1/2  : '  // TRIM(IntToChar(I)) &
                                     // ' Tr(P.F) = '// TRIM(DblToMedmChar(E_New))          &
                                     // '  %Non0  = '// TRIM(IntToChar(PNon0))              &
                                     // ' |[F,P]| = '// TRIM(DblToShrtChar(Max_Error_FX)) 
        ELSEIF(ITYPE==2) THEN
           Mssg = TRIM(Prog) // ':   Mu > 1/2   : ' // TRIM(IntToChar(I)) &
                                     // ' Tr(P.F) = '// TRIM(DblToMedmChar(E_New))          &
                                     //', %Non0   = '// TRIM(IntToChar(PNon0))              &
                                     // ' |[F,P]| = '// TRIM(DblToShrtChar(Max_Error_FX)) 

        ELSEIF(ITYPE==3) THEN            
           Mssg = TRIM(Prog) // ':Trace Constant: '  // TRIM(IntToChar(I)) &
                                     // ' Tr(P.F) = '// TRIM(DblToMedmChar(E_New))          &
                                     // ' %Non0   = '// TRIM(IntToChar(PNon0))              &
                                     // ' |[F,P]| = '// TRIM(DblToShrtChar(Max_Error_FX)) 
        ENDIF
        WRITE(*,*) TRIM(Mssg)
     ENDIF
     IF(Error_Energy < Tol_Energy .AND. ITYPE == 3) EXIT
  ENDDO
!
  CALL SetEq(P,X1)
  CALL Multiply(P,-0.50D0)
  CALL Add(P,0.50D0)
  TrP = Trace(P)
  Mssg = TRIM(Prog) // '               : |Trace(P)-NEl|= '//TRIM(DblToMedmChar(ABS(Two*TrP-NEl)))
  WRITE(*,*) TRIM(Mssg)
!
  CALL Delete(F)
  CALL Delete(T2)
  CALL Delete(T3)  
  CALL Delete(X1)   
  CALL Delete(X2)
  CALL Delete(X3)
!=============================================================================
!
!  TRANSFORMATION TO AN AO REPRESENTATION AND IO
!
!=============================================================================
!-----------------------------------------------------------------------------
!  IO for the orthogonal P 
!
   CALL Put(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)   
   CALL Put(P,TrixFile('OrthoD',Args,1))
   CALL PChkSum(P,'OrthoP['//TRIM(NxtCycl)//']',Prog)
   CALL PPrint( P,'OrthoP['//TRIM(NxtCycl)//']')
   CALL Plot(   P,'OrthoP_'//TRIM(NxtCycl))
!-----------------------------------------------------------------------------
!  Convert to AO representation
!
   INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
   IF(Present)THEN     
      CALL Get(Z,TrixFile('X',Args))   ! Z=S^(-1/2)
      CALL Multiply(Z,P,T1)
      CALL Multiply(T1,Z,P)
   ELSE
      CALL Get(Z,TrixFile('Z',Args))   ! Z=S^(-L)
      CALL Multiply(Z,P,T1)
      CALL Get(Z,TrixFile('ZT',Args))
      CALL Multiply(T1,Z,P)
   ENDIF
   CALL Filter(T1,P) 
!-----------------------------------------------------------------------------
!  IO for the non-orthogonal P 
!
   CALL Put(T1,TrixFile('D',Args,1))
   CALL Put(Zero,'homolumogap')
   CALL PChkSum(T1,'P['//TRIM(NxtCycl)//']',Prog)
   CALL PPrint(T1,'P['//TRIM(NxtCycl)//']')
   CALL Plot(T1,'P_'//TRIM(NxtCycl))
!-----------------------------------------------------------------------------
!  Tidy up 
!
   CALL Delete(P)
   CALL Delete(T1)
   CALL Delete(Z)
   CALL ShutDown(Prog) 
!
 END PROGRAM SDMM2
