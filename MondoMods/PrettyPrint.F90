!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    GENERIC PRETTY PRINTING FOR MONDOSCF TYPES
!
MODULE PrettyPrint
   USE DerivedTypes
   USE GlobalCharacters
   USE GlobalScalars
   USE GlobalObjects
   USE BasisSetParameters
   USE Clock
   USE Parse
   USE SetXYZ
   USE InOut
   USE Order
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
   INTERFACE PPrint
      MODULE PROCEDURE Print_INT_SCLR,   Print_DBL_SCLR,  &
                       Print_CHR_SCLR,   Print_INT_VECT,  &
                       Print_DBL_VECT,   Print_DBL_RNK2,  &
                       Print_DBL_Rank2A, Print_BSET,      &
                       Print_CRDS,       Print_BCSR,      &
#ifdef PARALLEL
                       Print_DBCSR,      &
#endif
                       Print_MEMS,       Print_TIME,      &
                       Print_HGRHO
   END INTERFACE
   INTERFACE PChkSum   
      MODULE PROCEDURE Print_CheckSum_BCSR
      MODULE PROCEDURE Print_CheckSum_HGRho
#ifdef PARALLEL 
      MODULE PROCEDURE Print_CheckSum_DBCSR 
#endif
   END INTERFACE
   CONTAINS 
      SUBROUTINE TimeStamp(Mssg,Enter_O)
         CHARACTER(LEN=*),INTENT(IN) :: Mssg
         LOGICAL,OPTIONAL,INTENT(IN) :: Enter_O
         CHARACTER(LEN=8)            :: DDate,DDay,HMS
         CHARACTER(LEN=10)           :: TTime
         CHARACTER(LEN=5)            :: Zone
         INTEGER, DIMENSION(8)       :: Values
         LOGICAL                     :: Enter 
         Enter=.TRUE.; IF(PRESENT(Enter_O))Enter=Enter_O
         CALL DATE_AND_TIME(DDate,TTime,Zone,Values)
         DDay=DDate(5:6)//'/'//DDate(7:8)//'/'//DDate(3:4)
         HMS=TTime(1:2)//':'//TTime(3:4)//':'//TTime(5:6)
         CALL OpenASCII(OutFile,Out)
         IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            WRITE(Out,'(A)')'(*'//TRIM(Mssg)//' '//DDay//' @ '//HMS//'*)'
         ELSEIF(Enter)THEN
            WRITE(Out,'(A)')'<<'//TRIM(Mssg)//' '//DDay//' @ '//HMS//'>>'
         ELSE
            WRITE(Out,'(A)')' - - - - - - - - done '//DDay//' @ '//HMS//' - - - - - - - - '
!            WRITE(Out,22)
!            22 FORMAT(72('-'))
         ENDIF
         CLOSE(UNIT=Out,STATUS='KEEP')
      END SUBROUTINE TimeStamp

      FUNCTION OpenPU(FileName_O,Unit_O)
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
         INTEGER                              :: OpenPU,PU
         IF(PRESENT(FileName_O).AND.PRESENT(Unit_O))THEN
            PU=Unit_O
            CALL OpenASCII(FileName_O,PU)
         ELSEIF(PRESENT(FileName_O).AND.(.NOT.PRESENT(Unit_O)))THEN
            PU=Tmp
            CALL OpenASCII(FileName_O,PU)
         ELSEIF(PRESENT(Unit_O))THEN
            IF(Unit_O==6)THEN
               PU=Unit_O
            ELSE
               CALL Halt(' Logic Error 1 in OpenPU ')
            ENDIF
         ELSE
            PU=Out
            CALL OpenASCII(OutFile,PU)
         ENDIF
         OpenPU=PU
      END FUNCTION OpenPU

      SUBROUTINE ClosePU(U)
         INTEGER,INTENT(IN) :: U
         IF(U/=6)CLOSE(UNIT=U,STATUS='KEEP')
      END SUBROUTINE ClosePU

      SUBROUTINE Print_CHR_SCLR(X,Name_O,FileName_O,Unit_O)
         CHARACTER(LEN=*),INTENT(IN)          :: X
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Name_O,FileName_O
         INTEGER, OPTIONAL,INTENT(IN)         :: Unit_O
         INTEGER                              :: PU               
         PU=OpenPU(FileName_O,Unit_O)
         IF(PrintFlags%Fmt==DEBUG_MMASTYLE.AND.PRESENT(Name_O))THEN
            WRITE(PU,11)TRIM(Name_O),TRIM(X)
         ELSEIF(PRESENT(Name_O))THEN
            WRITE(PU,12)TRIM(Name_O),TRIM(X)
         ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            WRITE(PU,13)TRIM(X)
         ELSE
            WRITE(PU,14)TRIM(X)
         ENDIF
         CALL ClosePU(PU)
      11 FORMAT(1x,A,' = ',A,';')
      12 FORMAT(1x,A,' = ',A) 
      13 FORMAT('(* ',A,' *)')
      14 FORMAT(1x,A)
      END SUBROUTINE Print_CHR_SCLR
!----------------------------------------------------------------PRINT AN INTEGER
      SUBROUTINE Print_INT_SCLR(X,Name,FileName_O,Unit_O,Protect_O)
         INTEGER,                   INTENT(IN) :: X
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O 
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O
         LOGICAL,          OPTIONAL,INTENT(IN) :: Protect_O
         INTEGER                               :: PU               
         CHARACTER(LEN=INTERNAL_INT_LEN)       :: CTmp,Id
         INTEGER                               :: I,J,M,N
         LOGICAL                               :: Protect
!--------------------------------------------------------
         CTmp=IntToChar(X)
         IF(PRESENT(Protect_O))THEN
            Protect=Protect_O
         ELSE
            Protect=.TRUE.
         ENDIF
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(Protect) &
               CALL PrintProtectL(PU)
            CALL ClosePU(PU)
#ifdef PARALLEL
         ENDIF
         IF(InParallel)THEN
            DO I=0,NPrc-1
               CALL AlignNodes()         
               IF(MyId==I)THEN
                  PU=OpenPU(FileName_O,Unit_O)              
                  Id=IntToChar(MyId)
                  IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
                     WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] = ',TRIM(CTmp),';'
                  ELSE
                     WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] := ',TRIM(CTmp)
                  ENDIF
                  CALL ClosePU(PU)
               ENDIF
            ENDDO
         ELSE
#endif
            PU=OpenPU(FileName_O,Unit_O)              
            IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
               WRITE(PU,*)TRIM(Name),' = ',TRIM(CTmp),';'
            ELSE
               WRITE(PU,*)TRIM(Name),' := ',TRIM(CTmp)
            ENDIF
            CALL ClosePU(PU)
#ifdef PARALLEL
         ENDIF
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(Protect) &
               CALL PrintProtectR(PU)
            CALL ClosePU(PU)
#ifdef PARALLEL
         ENDIF
#endif
      END SUBROUTINE Print_INT_SCLR
!----------------------------------------------------------------
!     PRINT A DOUBLE
!
      SUBROUTINE Print_DBL_SCLR(X,Name,FileName_O,Unit_O)
         REAL(DOUBLE),              INTENT(IN) :: X
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O
         INTEGER                               :: PU               
         PU=OpenPU(FileName_O,Unit_O)
         IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            WRITE(PU,11)Name,FRACTION(X),EXPONENT(X)
         ELSE
            WRITE(PU,*)Name,' = ',TRIM(DblToChar(X))
         ENDIF
         CALL ClosePU(PU)
      11 FORMAT(1x,A,' = ',F19.16,'*2^(',I4,');')
      END SUBROUTINE Print_DBL_SCLR
!----------------------------------------------------------------
!     PRINT AN INT_VECT
!
      SUBROUTINE Print_INT_VECT(A,Name,N_O,M_O,FileName_O,Unit_O)
         TYPE(INT_VECT),            INTENT(IN) :: A 
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O,M_O,N_O
         TYPE(CHR_VECT)                        :: CA
         CHARACTER(LEN=INTERNAL_INT_LEN)       :: Id
         INTEGER                               :: PU               
         INTEGER                               :: I,J,M,N
         M=1;         IF(PRESENT(M_O))M=M_O
         N=SIZE(A%I); IF(PRESENT(N_O))N=N_O
         CALL New(CA,N,M)
         DO J=M,N; CA%C(J)=IntToChar(A%I(J)); ENDDO
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectL(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
         DO I=0,NPrc-1
            CALL AlignNodes()         
            IF(MyId==I)THEN
               PU=OpenPU(FileName_O,Unit_O)              
               Id=IntToChar(MyId)
               WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] := ',(TRIM(CA%C(J)),', ',J=M,N)
               CALL ClosePU(PU)
            ENDIF
         ENDDO
#else
         PU=OpenPU(FileName_O,Unit_O)              
         WRITE(PU,*)TRIM(Name),' := ',(TRIM(CA%C(J)),', ',J=1,N)
         CALL ClosePU(PU)
#endif
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectR(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
#endif
         CALL Delete(CA)
      END SUBROUTINE Print_INT_VECT
!----------------------------------------------------------------
!     PRINT AN INT_VECT
!
      SUBROUTINE Print_DBL_VECT(A,Name,N_O,M_O,FileName_O,Unit_O)
         TYPE(DBL_VECT),            INTENT(IN) :: A 
         CHARACTER(LEN=*),          INTENT(IN) :: Name
         CHARACTER(LEN=*), OPTIONAL,INTENT(IN) :: FileName_O
         INTEGER,          OPTIONAL,INTENT(IN) :: Unit_O,M_O,N_O
         TYPE(CHR_VECT)                        :: CA
         CHARACTER(LEN=INTERNAL_INT_LEN)       :: Id
         INTEGER                               :: PU               
         INTEGER                               :: I,J,M,N
         M=1;         IF(PRESENT(M_O))M=M_O
         N=SIZE(A%D); IF(PRESENT(N_O))N=N_O
         CALL New(CA,N,M)
         DO J=M,N; CA%C(J)=DblToMedmChar(A%D(J)); ENDDO
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectL(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
         DO I=0,NPrc-1
            CALL AlignNodes()         
            IF(MyId==I)THEN
               PU=OpenPU(FileName_O,Unit_O)              
               Id=IntToChar(MyId)
               WRITE(PU,*)TRIM(Name),'[',TRIM(Id),'] := ',(TRIM(CA%C(J)),', ',J=M,N)
               CALL ClosePU(PU)
            ENDIF
         ENDDO
#else
         PU=OpenPU(FileName_O,Unit_O)              
         WRITE(PU,*)TRIM(Name),' := ',(TRIM(CA%C(J)),', ',J=1,N)
         CALL ClosePU(PU)
#endif
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            PU=OpenPU(FileName_O,Unit_O)
            IF(PU/=6)THEN
               CALL PrintProtectR(PU)
               CLOSE(PU)
            ENDIF
#ifdef PARALLEL
         ENDIF
#endif
         CALL Delete(CA)
      END SUBROUTINE Print_DBL_VECT
!----------------------------------------------------------------PRINT BASIS SET 
      SUBROUTINE Print_BSET(BS)
         TYPE(BSET) :: BS
         INTEGER :: NC,NP,MinL,MaxL
         INTEGER :: I,J,K,L,M
         IF(PrintFlags%Set/=DEBUG_BASISSET)RETURN
         CALL OpenASCII(OutFile,Out)
         CALL PrintProtectL(Out)
         WRITE(Out,*)'Internal representation of the ', &
                     TRIM(BS%BName),' basis set: '//Rtrn
         DO I=1,BS%NKind
            NC=BS%NCFnc%I(I)
            WRITE(Out,1002)Ats(BS%Kinds%I(I)),NC
            DO J=1,NC
               NP=BS%NPFnc%I(J,I)
               MinL=LBegin(BS%ASymm%I(1,J,I))
               MaxL=LEnd(BS%ASymm%I(2,J,I))
               IF(NP==1)THEN
                 WRITE(Out,1103)J,NP
                 IF(MaxL==1)THEN
                    WRITE(Out,1104)
                 ELSE
                    WRITE(Out,1004)
                 ENDIF
               ELSE
                 WRITE(Out,1003)J,NP
                 WRITE(Out,1004)
               ENDIF
               WRITE(Out,1005)(ASymmTyps(L),L=MinL,MaxL)
               DO K=1,NP
                  WRITE(Out,1006)K,BS%Expnt%D(K,J,I), &             
                                (BS%CCoef%D(M,K,J,I),M=MinL,MaxL)
              ENDDO
              WRITE(Out,*)Rtrn
            ENDDO
        ENDDO
        CALL PrintProtectR(Out)
        CLOSE(Out)
   1001 FORMAT(72('='))
   1002 FORMAT(1x,A2,' has ',I2,' associated contractions,')
   1003 FORMAT(1x,'Contraction ',I2,' involves', &
               I2,' primitives : ')
   1004 FORMAT(1x,'Primitive  Exponent',7x,'Normalized Coeficients')
   1103 FORMAT(1x,'Contraction ',I2,' involves', &
               I2,' primitive : ')
   1104 FORMAT(1x,'Primitive  Exponent',7x,'Normalized Coeficient')

   1005 FORMAT(18x,20(12x,A3))
   1006 FORMAT(1x,I2,6x,8(1x,D14.8))
   1007 FORMAT(60('='))
     END SUBROUTINE Print_BSET 
!----------------------------------------------------------------PRINT COORDINATES
     SUBROUTINE Print_CRDS(GM,FileName_O,Unit_O,PrintGeom_O)
        TYPE(CRDS) :: GM         
        INTEGER :: K
        LOGICAL :: Opened
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O,PrintGeom_O
        INTEGER                              :: I,PU
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Mssg
        REAL(DOUBLE)                         :: AA
#ifdef PARALLEL
        IF(MyId==ROOT)THEN
#endif
!          
          PU=OpenPU(FileName_O=FileName_O,Unit_O=Unit_O)
          IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
              WRITE(PU,*)LeftParenStar
              WRITE(PU,43)
           ELSE
              WRITE(PU,44)
           ENDIF
           Mssg='Configuration #'//TRIM(IntToChar(GM%Confg))
           IF(GM%ETotal/=BIG_DBL)THEN
              Mssg='Escf['//TRIM(IntToChar(GM%Confg))//']='//TRIM(DblToChar(GM%ETotal)) 
              WRITE(PU,*)TRIM(Mssg)
              Mssg='Grms['//TRIM(IntToChar(GM%Confg))//']='//TRIM(DblToShrtChar(GM%GradRMS)) 
              WRITE(PU,*)TRIM(Mssg)
              Mssg='Gmax['//TRIM(IntToChar(GM%Confg))//']='//TRIM(DblToShrtChar(GM%GradMax)) 
              WRITE(PU,*)TRIM(Mssg)
!              Mssg='Escf['//TRIM(IntToChar(GM%Confg))//']='//TRIM(DblToMMAChar(GM%ETotal))//';' 
!              WRITE(PU,*)TRIM(Mssg)
!              Mssg='Grms['//TRIM(IntToChar(GM%Confg))//']='//TRIM(DblToMMAChar(GM%GradRMS))//';' 
!              WRITE(PU,*)TRIM(Mssg)
!              Mssg='Gmax['//TRIM(IntToChar(GM%Confg))//']='//TRIM(DblToMMAChar(GM%GradMax))//';' 
!              WRITE(PU,*)TRIM(Mssg)
           ENDIF
           IF(PrintGeom_O=='XYZ')THEN
               AA=One/AngstromsToAU
               DO I=1,GM%NAtms
                 Mssg=Ats(GM%AtNum%I(I))                    &
                 //'   '//DblToMedmChar(GM%Carts%D(1,I)*AA) &
                 //'   '//DblToMedmChar(GM%Carts%D(2,I)*AA) &
                 //'   '//DblToMedmChar(GM%Carts%D(3,I)*AA) 
                 WRITE(PU,*)TRIM(Mssg)
              ENDDO
           ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              IF(GM%InAU)THEN
                 WRITE(PU,45)
              ELSE
                 WRITE(PU,46)
              ENDIF
              IF(GM%Ordrd==SFC_HILBERT)THEN
                 WRITE(PU,17)
              ELSEIF(GM%Ordrd==SFC_PEANO)THEN
                 WRITE(PU,27)
              ELSEIF(GM%Ordrd==SFC_RANDOM)THEN
                 WRITE(PU,47)
              ELSE
                 WRITE(PU,37)
              ENDIF
              WRITE(PU,49)
           ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
              WRITE(PU,*)RightParenStar
              DO I=1,GM%NAtms
                 WRITE(PU,54)I,GM%AtNum%I(I)
              ENDDO
              DO I=1,GM%NAtms
                 WRITE(PU,53)I,(FRACTION(GM%Carts%D(K,I)), &
                                 EXPONENT(GM%Carts%D(K,I)),K=1,3)
              ENDDO
           ENDIF
              DO I=1,GM%NAtms
                 WRITE(PU,52)I,Ats(GM%AtNum%I(I)),(GM%Carts%D(K,I),K=1,3)
              ENDDO
!              WRITE(PU,42) GM%ENucN
              WRITE(PU,444)GM%NElec
              WRITE(PU,43) 
#ifdef PERIODIC
              WRITE(PU,55)
              WRITE(PU,56) GM%AutoW(1),GM%AutoW(2),GM%AutoW(3)
              WRITE(PU,57) GM%TransVec%D(1),GM%TransVec%D(2),GM%TransVec%D(3)
              WRITE(PU,58)  
              WRITE(PU,59) GM%BoxShape%D(1,1),GM%BoxShape%D(2,1),GM%BoxShape%D(3,1)
              WRITE(PU,60) GM%BoxShape%D(1,2),GM%BoxShape%D(2,2),GM%BoxShape%D(3,2)
              WRITE(PU,61) GM%BoxShape%D(1,3),GM%BoxShape%D(2,3),GM%BoxShape%D(3,3)
              WRITE(PU,62)  
              DO I=1,GM%NAtms
                 WRITE(PU,63) I, Ats(GM%AtNum%I(I)),(GM%BoxCarts%D(K,I),K=1,3)
              ENDDO
              WRITE(PU,43)
#endif
!           ENDIF
           CALL ClosePU(PU)
#ifdef PARALLEL
        ENDIF
#endif
        IF(PrintFlags%Fmt==DEBUG_MMASTYLE) &
        CALL PPrint(GM%NElec,'NEl',Protect_O=.FALSE.)
!
   43   FORMAT(72('='))
   44   FORMAT(72('='))
   41   FORMAT(' ENucN=',F19.16,'*2^(',I4,');')
   42   FORMAT(2x,' E_{nn} = ',D22.14,' (AU).')
  444   FORMAT(2x,' N_{el} = ',I8)
   45   FORMAT(2x,' Geometry, originally in AU, has not been rescaled.')
   46   FORMAT(2x,' Geometry, originally in Angstroms, has been converted to AU.')

   17   FORMAT(2x,' Geometry has been reordered using the Hilbert curve.')
   27   FORMAT(2x,' Geometry has been reordered using the Peano curve.')
   47   FORMAT(2x,' Geometry has been randomly reordered.')
   37   FORMAT(2x,' Geometry has not been reordered.')

   49   FORMAT(2x,' Internal representation of geometry: ')
   52   FORMAT(1x,I5,1x,A2,3(1x,D16.8))
   53   FORMAT(' R[',I4,']={',F19.16,'*2^(',I4,')', &
                          ',',F19.16,'*2^(',I4,')', &
                          ',',F19.16,'*2^(',I4,')};')
   54   FORMAT(' Z[',I4,']=',I3,';')
   55   FORMAT(2x,' Representation of the Lattice')
   56   FORMAT(2x,' Periodic Boundry Conditions (x,y,z) => {',L1,',',L1,',',L1,'}')
   57   FORMAT(2x,' Translation Vector = (',D14.8,', ',D14.8,', ',D14.8,')')
   58   FORMAT(2x,' Lattice Vectors: ')
   59   FORMAT(2x,'  a  = (',D14.8,', ',D14.8,', ',D14.8,')')
   60   FORMAT(2x,'  b  = (',D14.8,', ',D14.8,', ',D14.8,')')
   61   FORMAT(2x,'  c  = (',D14.8,', ',D14.8,', ',D14.8,')')
   62   FORMAT(2x,' Internal representation of geometry in Fractional Coordinates: ')
   63   FORMAT(1x,I5,1x,A2,3(1x,F16.10))
     END SUBROUTINE Print_CRDS
!-----------------------------------------------------------------------------
!    Print a BCSR matrix
!
     SUBROUTINE Print_BCSR(A,Name,FileName_O,Unit_O)       
        TYPE(BCSR)                           :: A
        TYPE(DBL_RNK2)                       :: B
        CHARACTER(LEN=*),INTENT(IN)          :: Name
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O                
       
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O
       IF(PrintFlags%Mat/=DEBUG_MATRICES.AND.(.NOT.PRESENT(Unit_O)))RETURN        
#ifdef PARALLEL 
        IF(MyId==ROOT)THEN
#endif
           CALL SetEq(B,A)
           CALL Print_DBL_RNK2(B,Name,FileName_O,Unit_O)
           CALL Delete(B)
#ifdef PARALLEL 
        ENDIF
#endif
     END SUBROUTINE Print_BCSR
#ifdef PARALLEL 
!-----------------------------------------------------------------------------
!    Print a DBCSR matrix
!
     SUBROUTINE Print_DBCSR(A,Name,FileName_O,Node_O,Distrib_O)       
        TYPE(DBCSR), INTENT(INOUT)           :: A
        TYPE(DBL_RNK2)                       :: B
        TYPE(BCSR)                           :: C
        CHARACTER(LEN=*),INTENT(IN)          :: Name
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Name2
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: FileName_O
        INTEGER,OPTIONAL                     :: Node_O
        LOGICAL,OPTIONAL                     :: Distrib_O
        INTEGER                              :: I
        IF(PrintFlags%Mat/=DEBUG_MATRICES)RETURN
        IF(PRESENT(Distrib_O))THEN
           IF(Distrib_O)THEN
              CALL SetEq(B,A) 
              DO I=0,NPrc-1
                 IF(InParallel)CALL AlignNodes()
                 IF(MyId==I)THEN
                    Name2=TRIM(Name)//'['//TRIM(IntToChar(I))//']'
                    CALL Print_DBL_RNK2(B,Name2,OutFile,Unit_O)
                 ENDIF
              ENDDO
              CALL Delete(B)
           ELSE
              CALL Halt(' Logic error 1 in Print_DBCSR ')
           ENDIF          
        ELSEIF(PRESENT(Node_O))THEN
           IF(MyId==Node_O)THEN
              CALL SetEq(B,A) 
              Name2=TRIM(Name)//'['//TRIM(IntToChar(Node_O))//']'
              CALL Print_DBL_RNK2(B,Name2,OutFile)
              CALL Delete(B)
           ENDIF
        ELSE
           CALL SetEq(C,A)
           CALL Print_BCSR(C,Name,Filename_O)
           CALL Delete(C)
        ENDIF
     END SUBROUTINE Print_DBCSR
#endif

     SUBROUTINE Print_DBL_RNK2(A,Name,FileName_O,Unit_O)       
        TYPE(DBL_RNK2)             :: A
        CHARACTER(LEN=*)           :: Name
        INTEGER, OPTIONAL          :: Unit_O                
        CHARACTER(LEN=*), OPTIONAL :: FileName_O
        IF(.NOT.AllocQ(A%Alloc)) &
           CALL Halt(' Matrix not allocated in Print_DBL_RNK2.')
        CALL Print_DBL_Rank2A(A%D,Name,FileName_O,Unit_O)
     END SUBROUTINE Print_DBL_RNK2

     SUBROUTINE Print_DBL_Rank2A(A,Name,FileName_O,Unit_O)
        REAL(DOUBLE),DIMENSION(:,:),INTENT(IN) :: A
        INTEGER                                :: I,J,K,L,M,N,Unit
        CHARACTER(LEN=*)                       :: Name
        INTEGER, OPTIONAL                      :: Unit_O
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN)   :: FileName_O
        Unit=Out; IF(PRESENT(Unit_O))Unit=Unit_O
        IF(PRESENT(FileName_O).AND.Unit/=6)THEN
           CALL OpenASCII(FileName_O,Unit)
        ELSEIF(Unit/=6)THEN            
           CALL OpenASCII(OutFile,Unit)
        ENDIF
        M=SIZE(A,1); N=SIZE(A,2)        
        IF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
           IF(M/=N) &
              CALL Halt(' Print_DBL_Rank2A does not currently'//   &
                        ' support rectangular matrices with'//     & 
                        ' PrintFlags%Fmt==DEBUG_MMASTYLE:'//Rtrn// &
                        ' attempting to print '//TRIM(NAME))
           WRITE(Unit,100)
           WRITE(Unit,101)
           WRITE(Unit,102)N
           DO I=1,N,2
              K=MIN(I+1,N)
              IF(K-I+1.EQ.2)THEN
                 WRITE(Unit,202)(J,J=I,K)  
              ELSEIF(K-I+1.EQ.1)THEN
                 WRITE(Unit,201)(J,J=I,K)  
              ENDIF
              DO L=1,N
                  IF(K-I+1.EQ.2)THEN
                     WRITE(Unit,302)L,L, &
                       (FRACTION(A(L,J)),EXPONENT(A(L,J)),J=I,K)
                  ELSEIF(K-I+1.EQ.1)THEN
                     WRITE(Unit,301)L,L, &
                       (FRACTION(A(L,J)),EXPONENT(A(L,J)),J=I,K)
                  ENDIF
              ENDDO
           ENDDO
           WRITE(Unit,103)N
           WRITE(Unit,104)TRIM(Name)
           WRITE(Unit,99)
        ELSE
           WRITE(Unit,401)TRIM(Name)
!
           DO I=1,N,10
              K=MIN(I+9,N)
              WRITE(Unit,501)(J,J=I,K)  
              DO L=1,M
                WRITE(Unit,701) L,(A(L,J),J=I,K)
              ENDDO
           ENDDO
!
        ENDIF
        IF(Out/=6) CLOSE(Out)
        RETURN
!
   100  FORMAT(' (*',65('='),'*)')
    99  FORMAT(' mv[m_List]:=MatrixForm[ N[ Chop[m,0.001],4]]; ')
   101  FORMAT(1x,'Ap[x_List,y_List]:=Append[x,y];')
   102  FORMAT(1x,'Tmp=Table[{},{i,',I4,'}];')
   103  FORMAT(1x,'Do[Tmp[[i]]=Flatten[Tmp[[i]]],{i,1,',I4,'}];')
   104  FORMAT(1x,A,' = Tmp ; ')
   201  FORMAT('(*           ',1(14x,i4),'   *)')
   202  FORMAT('(*           ',2(14x,i4),'   *)')
   301  FORMAT(1x,'Tmp[[',I4,']]=Ap[Tmp[[',I4,']],', &
                         '{',F19.16,'*2^(',I4,')}];')
   302  FORMAT(1x,'Tmp[[',I4,']]=Ap[Tmp[[',I4,']],', &
                       '{',1(F19.16,'*2^(',I4,'),'), &
                             F19.16,'*2^(',I4,')}];')
   401  FORMAT(2x,A,':')
   501  FORMAT(T2,10I16)
   601  FORMAT(I5,10F10.5)
   701  FORMAT(I5,10D16.8)
!
   END SUBROUTINE Print_DBL_Rank2A
!
!==================================================================
!
!    Print Check Sums
!      
!===============================================================
!
     SUBROUTINE Print_CheckSum_BCSR(A,Name,Proc_O,Unit_O)
        TYPE(BCSR), INTENT(IN)               :: A
        REAL(DOUBLE)                         :: Chk
        CHARACTER(LEN=*)                     :: Name
        INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
        CHARACTER(LEN=*),OPTIONAL            :: Proc_O
        INTEGER                              :: I,PU
        CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
        IF(PrintFlags%Key/=DEBUG_MAXIMUM)RETURN
        Chk=Zero
        DO I=1,A%NNon0
           Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
        ENDDO
!        Chk=SQRT(Chk) avoid this to maintain F77 complience
#ifdef PARALLEL
        IF(MyID==ROOT)THEN
#endif
!           CALL OpenASCII(OutFile,Out)
           IF(PRESENT(Proc_O).AND.PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
              ChkStr=ProcessName(Proc_O)//TRIM(Name) &
                   //' matrix check sum = '//TRIM(DblToChar(Chk))
           ELSEIF(PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
              ChkStr=TRIM(Name)//' matrix check sum = '//TRIM(DblToChar(Chk))
           ELSEIF(PRESENT(Proc_O).AND.PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
              ChkStr='(* '//TRIM(Proc_O)//' *)'//'ChkSum'//TRIM(Name)       &
                   //' = '//TRIM(FltToChar(FRACTION(Chk))) &
                   //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
           ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
              ChkStr='ChkSum'//TRIM(Name)//' = '//TRIM(FltToChar(FRACTION(Chk))) &
                   //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
           ELSE 
              CALL Halt(' Logic error in Print_CheckSum_BCSR')
           ENDIF
           PU=OpenPU(Unit_O=Unit_O)
           WRITE(PU,'(1x,A)')TRIM(ChkStr)
           CALL ClosePU(PU)
#ifdef PARALLEL
        ENDIF
#endif
   END SUBROUTINE Print_CheckSum_BCSR
   
#ifdef PARALLEL
   SUBROUTINE Print_CheckSum_DBCSR(A,Name,Proc_O,Unit_O)
      TYPE(DBCSR), INTENT(IN)   :: A
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Proc_O
      INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
      REAL(DOUBLE)              :: Chk
      REAL(DOUBLE), EXTERNAL    :: DDot
      REAL(DOUBLE)              :: DotPrd
      CHARACTER(LEN=*)          :: Name
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: ChkStr
      INTEGER :: I,PU
      IF(PrintFlags%Key/=DEBUG_MAXIMUM)RETURN
      Chk=Zero
      DO I=1,A%NNon0
         Chk=Chk+A%MTrix%D(I)*A%Mtrix%D(I)
      ENDDO
      DotPrd=Reduce(Chk)
      IF(MyID==ROOT)THEN
!         Chk=SQRT(DotPrd) ! maintain F77 complience for now ...
         Chk=DotPrd
         IF(PRESENT(Proc_O).AND.PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
            ChkStr=ProcessName(Proc_O)//TRIM(Name) &
                 //' matrix check sum = '//TRIM(DblToChar(Chk))
         ELSEIF(PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
            ChkStr=TRIM(Name)//' matrix check sum = '//TRIM(DblToChar(Chk))
         ELSEIF(PRESENT(Proc_O).AND.PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            ChkStr='(* '//TRIM(Proc_O)//' *)'//TRIM(Name)       &
                 //'ChkSum = '//TRIM(FltToChar(FRACTION(Chk))) &
                 //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
         ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
            ChkStr=TRIM(Name)//'ChkSum = '//TRIM(FltToChar(FRACTION(Chk))) &
                 //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
         ELSE 
            CALL Halt(' Logic error in Print_CheckSum_DBCSR')
         ENDIF
         PU=OpenPU(Unit_O=Unit_O)
         WRITE(PU,'(1x,A)')ChkStr
         CALL ClosePU(PU)
      ENDIF
   END SUBROUTINE Print_CheckSum_DBCSR
#endif
!
  SUBROUTINE Print_CheckSum_HGRho(A,Name,Proc_O,Unit_O)
    TYPE(HGRho)                      :: A
    CHARACTER(LEN=*)                 :: Name   
    INTEGER,         OPTIONAL,INTENT(IN) :: Unit_O
    CHARACTER(LEN=*),OPTIONAL            :: Proc_O
    INTEGER                              :: PU,I,L,M,N,LMN,jadd,zq,iq,oq,orr,Ell,LenKet
    REAL(DOUBLE)                         :: Chk
    CHARACTER(LEN=2*DEFAULT_CHR_LEN)     :: ChkStr
    IF(PrintFlags%Key/=DEBUG_MAXIMUM)RETURN
    Chk=Zero
!   
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' Density not allocated in Print_CheckSum_HGRho')
    ENDIF
!
    DO zq=1,A%NExpt-1
       oq =A%OffQ%I(zq)
       orr=A%OffR%I(zq)
       Ell    = A%Lndx%I(zq) 
       LenKet = LHGTF(Ell)
       DO iq=1,A%NQ%I(zq)
          jadd=orr+(iq-1)*LenKet
          DO LMN=1,LenKet  
             Chk=Chk+ABS(A%Co%D(jadd+LMN))
          ENDDO
       ENDDO
    ENDDO
!
    IF(PRESENT(Proc_O).AND.PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
       ChkStr=ProcessName(Proc_O)//TRIM(Name)//' matrix check sum = '//TRIM(DblToChar(Chk))
    ELSEIF(PrintFlags%Fmt/=DEBUG_MMASTYLE)THEN
       ChkStr=TRIM(Name)//' matrix check sum = '//TRIM(DblToChar(Chk))
    ELSEIF(PRESENT(Proc_O).AND.PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
       ChkStr='(* '//TRIM(Proc_O)//' *)'//'ChkSum'//TRIM(Name)       &
                   //' = '//TRIM(FltToChar(FRACTION(Chk))) &
                   //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
     ELSEIF(PrintFlags%Fmt==DEBUG_MMASTYLE)THEN
       ChkStr='ChkSum'//TRIM(Name)//' = '//TRIM(FltToChar(FRACTION(Chk))) &
                      //'*2^('//TRIM(IntToChar(EXPONENT(Chk)))//');'
     ELSE 
       CALL Halt(' Logic error in Print_CheckSum_HGRho')
     ENDIF
!
     PU=OpenPU(Unit_O=Unit_O)
     WRITE(PU,'(1x,A)')TRIM(ChkStr)
     CALL ClosePU(PU)
!
  END SUBROUTINE Print_CheckSum_HGRho

!========================================================================================
!     Print Out the Density 
!========================================================================================
  SUBROUTINE Print_HGRho(A,Name,FileName_O,Unit_O)
    TYPE(HGRho)                      :: A
    CHARACTER(LEN=*)                 :: Name   
    CHARACTER(LEN=*),OPTIONAL        :: FileName_O
    INTEGER,OPTIONAL                 :: Unit_O
    INTEGER                          :: PU
    INTEGER                          :: L,M,N,LMN,NPrim,iadd,jadd,zq,iq,oq,orr,Ell,LenKet
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Strng
!
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' Density not allocated in PPrint_HGRho')
    ENDIF
!
    PU=OpenPU(Unit_O=Unit_O)
!
    IF(PrintFlags%Fmt==DEBUG_MMASTYLE) THEN
       NPrim=0
       WRITE(PU,*)'ClearAll[Rho];'
       WRITE(PU,*)'Rho[R_List]:=Module[{zeta,zs,Q,RQ,ExpRQ,Lx,Ly,Lz,RhoSum},'
       WRITE(PU,*)'RhoSum=0;'
       DO zq=1,A%NExpt-1
          oq =A%OffQ%I(zq)
          orr=A%OffR%I(zq)
          Ell=A%Lndx%I(zq) 
          LenKet=LHGTF(Ell)
          Strng=Squish('zeta='//DblToMMAChar(A%Expt%D(zq))//';')
          WRITE(PU,*)TRIM(Strng)
          WRITE(PU,*)'zs=Sqrt[zeta];'
          IF(A%NQ%I(zq).NE.0) THEN
             DO iq=1,A%NQ%I(zq)
                NPrim=NPrim+1
                iadd=oq+iq
                jadd=orr+(iq-1)*LenKet
                Strng=Squish('Q={'//DblToMMAChar(A%Qx%D(iadd))   &
                             //','//DblToMMAChar(A%Qy%D(iadd))   &  
                             //','//DblToMMAChar(A%Qz%D(iadd))//'};')
                WRITE(PU,*)TRIM(Strng)
                WRITE(PU,*)'RQ=R-Q;'
                WRITE(PU,*)'ExpRQ=Exp[-zeta*RQ.RQ];'
                Strng='Do[Lx[l]=zs^l*HermiteH[l,zs*RQ[[1]]];'//Rtrn &
                   //'    Ly[l]=zs^l*HermiteH[l,zs*RQ[[2]]];'//Rtrn &
                   //'    Lz[l]=zs^l*HermiteH[l,zs*RQ[[3]]];'//Rtrn &
                   //',{l,0,'//TRIM(IntToChar(Ell))//'}]; '
                WRITE(PU,*)TRIM(Strng)
                DO L=0,Ell
                DO M=0,Ell-L
                DO N=0,Ell-L-M
                   LMN=LMNDex(L,M,N)+jadd
                   Strng=Squish('RhoSum=RhoSum+Lx['//IntToChar(L)//']*Ly['//IntToChar(M) &
                              //']*Lz['//IntToChar(N)//']'//'*ExpRQ*'//DblToMMAChar(A%Co%D(LMN))//';')
                   WRITE(PU,*)TRIM(Strng)                    
                ENDDO
                ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       WRITE(PU,*)'RhoSum];'
    ELSEIF(PrintFlags%Key==DEBUG_MEDIUM .OR. PrintFlags%Key==DEBUG_MAXIMUM) THEN 
       NPrim=0
       WRITE(PU,30) Name
       WRITE(PU,31)
       WRITE(PU,32) A%NExpt
       WRITE(PU,33) A%NDist
       WRITE(PU,34) A%NCoef
       WRITE(PU,31)
       IF(PrintFlags%Key==DEBUG_MAXIMUM .AND. A%NDist .LT. 100) THEN
          DO zq=1,A%NExpt
             oq =A%OffQ%I(zq)
             orr=A%OffR%I(zq)
             Ell    = A%Lndx%I(zq) 
             LenKet = LHGTF(Ell)
             IF(A%NQ%I(zq).NE.0) THEN
                DO iq=1,A%NQ%I(zq)
                   NPrim=NPrim+1
                   iadd=oq+iq
                   jadd=orr+(iq-1)*LenKet
                   WRITE(PU,20) NPrim,Ell,A%Expt%D(zq)
                   WRITE(PU,21) A%Qx%D(iadd),A%Qy%D(iadd),A%Qz%D(iadd)
                   WRITE(PU,22)
                   DO L=0,Ell
                      DO M=0,Ell-L
                         DO N=0,Ell-L-M
                            LMN=LMNDex(L,M,N)+jadd
                            WRITE(PU,24) L,M,N
                            WRITE(PU,25) A%Co%D(LMN)
                         ENDDO
                      ENDDO
                   ENDDO
                   WRITE(PU,31)
                ENDDO
             ENDIF
          ENDDO
       ENDIF
    ENDIF
!
    CALL ClosePU(PU)
!
    RETURN
!
10  FORMAT(' zeta[',I4,']=SetPrecision[',F19.16,'*2^(' ,I4, ')},50];')
11  FORMAT(' Ell[',I4,']=',I2,';')
12  FORMAT(' Q[',I4,']=SetPrecision[{',F19.16,'*2^(' ,I4, '),', & 
                                       F19.16,'*2^(' ,I4, '),', & 
                                       F19.16,'*2^(' ,I4, ')},50];')
14  FORMAT(' Co[',I4,', ',I2,', ',I2,', ',I2,']=SetPrecision[',F19.16,'*2^(',I4,'),50];')
15  FORMAT(' NPrim=',I6)
!
20  FORMAT(2x,'Primative #',I4,2x,' Max L = ',I4,2x,' EXPONENT = ',D14.8)
21  FORMAT(2x,'QR     = (',D14.8,', ',D14.8,', ',D14.8,')')
22  FORMAT(2x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
24  FORMAT(3x,'L=',I4,' M=',I4,' N=',I4)
25  FORMAT(3x,'RhoCo(L,M,N) = ',D14.8)
30  FORMAT(1x,A,' in a Hermite Gaussian basis: ')
31  FORMAT(1x,'=========================================================================')
32  FORMAT(1x,' Number of Exponents    = ',I8,2x)
33  FORMAT(1x,' Number of Distribution = ',I8,2x)
34  FORMAT(1x,' Number of Coeffients   = ',I8,2x)
  END SUBROUTINE Print_HGRho
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  FUNCTION BlockIndexString(Name,I,K,J,Id)
    INTEGER                         :: I,K
    INTEGER, OPTIONAL               :: J,Id
    CHARACTER(LEN=*)                :: Name
    CHARACTER(LEN=INTERNAL_INT_LEN) :: cI,cJ,cK,cId
    CHARACTER(LEN=64)               :: BlockIndexString
    cI=IntToChar(I)
    cK=IntToChar(K)
    IF(PRESENT(J).AND.PRESENT(Id))THEN
       cJ=IntToChar(J)
       cId=IntToChar(Id)
       BlockIndexString=TRIM(Name)//'(Id='//TRIM(cId)//',I='// &
            TRIM(cI)//',J='//TRIM(cK)//',K='//TRIM(cJ)//')'
    ELSEIF(PRESENT(J))THEN
       cJ=IntToChar(J)
       BlockIndexString=TRIM(Name)//'('// &
            TRIM(cI)//',J='//TRIM(cJ)//','//TRIM(cK)//')'
    ELSEIF(PRESENT(Id))THEN
       cId=IntToChar(Id)
       BlockIndexString=TRIM(Name)//'(Id='//TRIM(cId)//','// &
            TRIM(cI)//','//TRIM(cK)//')'
    ELSE
       BlockIndexString=TRIM(Name)//'('//TRIM(cI)//','//TRIM(cK)//')'
    ENDIF
  END FUNCTION BlockIndexString
!---------------------------------------------------------------------
!     
!---------------------------------------------------------------------
  SUBROUTINE Elapsed_TIME(T,Init_O,Proc_O)
    TYPE(TIME),           INTENT(INOUT)  :: T
    REAL(DOUBLE)                         :: WallTm,CPUSTm,TimeTot,TimeMax, &
         TimeMin,TimeAve,Imb
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Init_O,Proc_O
    CHARACTER(LEN=2*DEFAULT_CHR_LEN)     :: Mssg
    CHARACTER(LEN=10)                    :: Mssg10
    INTEGER                              :: I,PU
    IF(PRESENT(Init_O))THEN
       IF(Init_O=='Init')THEN
          T%FLOP=0
          T%Wall=Zero
          T%CPUS=Zero
          T%WStrt=WallSec()
          T%CStrt=CPUSec()
       ELSEIF(Init_O=='Start')THEN
          T%WStrt=WallSec()
          T%CStrt=CPUSec()
       ELSEIF(Init_O=='Accum')THEN
          WallTm=WallSec()
          CPUSTm=CPUSec()
          T%Wall=T%Wall+(WallTm-T%WStrt)
          T%CPUS=T%CPUS+(CPUSTm-T%CStrt)
          T%WStrt=WallTm
          T%CStrt=CPUSTm
       ELSE
          CALL Halt('Unknown option '//TRIM(Init_O)//' in Print_Elapsed_TIME')
       ENDIF
    ENDIF
#ifdef PARALLEL
    IF(PRESENT(Proc_O))THEN
       IF(InParallel)THEN
          TimeTot=Reduce(T%Wall,MPI_SUM)
          TimeMax=Reduce(T%Wall,MPI_MAX)
          TimeMin=Reduce(T%Wall,MPI_MIN)
       ENDIF
       IF(MyId==ROOT.AND.InParallel.AND. &
            PrintFlags%Key>=DEBUG_MEDIUM)THEN
          !             Compute relative imbalance 
          TimeAve=TimeTot/DBLE(NPrc)
          IF(TimeTot/=Zero)THEN
             Imb=ABS(TimeMax-TimeAve)/TimeAve
          ELSE
             Imb=Zero
          ENDIF
          Mssg=ProcessName(Proc_O)//'Parallel Statistics '
          PU=OpenPU()
          WRITE(PU,*)TRIM(Mssg)
          Mssg=       '                    Imblnce = '//TRIM(DblToShrtChar(Imb))      &
               //Rtrn//'                     MinTime = '//TRIM(DblToShrtChar(TimeMin)) &
               //Rtrn//'                     AvgTime = '//TRIM(DblToShrtChar(TimeAve)) &
               //Rtrn//'                     MaxTime = '//TRIM(DblToShrtChar(TimeMax)) 
          WRITE(PU,*)TRIM(Mssg)
          CLOSE(Out)          
       ENDIF
    ENDIF
#endif
  END SUBROUTINE Elapsed_TIME
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  FUNCTION ProcessName(Proc,Misc_O) RESULT (Tag)
    CHARACTER(LEN=*)           :: Proc
    CHARACTER(LEN=*), OPTIONAL :: Misc_O
    CHARACTER(LEN=20)          :: Tag
    CHARACTER(LEN=16)          :: Name
    CHARACTER(LEN=3 ),PARAMETER:: Colon =","
    CHARACTER(LEN=4 ),PARAMETER:: Colons=" :: "
    IF(PRESENT(Misc_O))THEN
       Name=TRIM(ADJUSTL(Proc))//Colon//TRIM(Misc_O)
       Tag=Name//Colons
    ELSE
       Name=ADJUSTL(TRIM(Proc))
       Tag=Name//Colons
    ENDIF
  END FUNCTION ProcessName
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE Print_TIME(T,Proc_O,FileName_O,Unit_O,BareBones_O)
    TYPE(TIME),               INTENT(INOUT) :: T
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Proc_O,FileName_O
    INTEGER, OPTIONAL                       :: Unit_O
    REAL(DOUBLE)                            :: Elapsed_CPUS,Elapsed_WALL,FLOPS
    LOGICAL, OPTIONAL                       :: BareBones_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)          :: Mssg,Proc
    INTEGER                                 :: PU
!------------------------------------------------------------------------------------
#ifdef PARALLEL
    IF(InParallel)  &
         CALL AlignNodes()
    IF(MyId==ROOT)THEN
#endif
       Elapsed_CPUS=T%CPUS
       Elapsed_Wall=T%Wall
       PU=OpenPU(FileName_O,Unit_O)
       CALL PrintProtectL(PU)
#ifdef PARALLEL
    ENDIF
    IF(InParallel)THEN
       FLOPS=Reduce(T%FLOP)
    ELSE
#endif
       FLOPS=T%FLOP
#ifdef PARALLEL
    ENDIF
#endif
!----------------------------------------------------
    IF(PRESENT(Proc_O))THEN
       Proc=Proc_O
    ELSE
       Proc=Blnk
    ENDIF
!-------------------------------------------------------------------------
    IF(PRESENT(BareBones_O))THEN
       IF(.NOT.BareBones_O)  &
            CALL Halt(' Logic error in Print_TIME')
#ifdef PARALLEL
       IF(MyId==ROOT)THEN
          WRITE(PU,10)NPrc,Elapsed_Wall,MFlops(FLOPS,Elapsed_Wall)
10        FORMAT(I4,' ',D12.6,' ',I10)
       ENDIF
#else
       WRITE(PU,10)Elapsed_Wall,MFlops(FLOPS,Elapsed_Wall)
10     FORMAT(' ',D12.6,' ',I10)
#endif
       CALL PrintProtectR(PU)
       CLOSE(Out)
       RETURN
    ENDIF
!-------------------------------------------------------------------------
#ifdef PARALLEL
    IF(MyID==ROOT)THEN
#endif
       IF(FLOPS>Zero)THEN
          IF(Elapsed_CPUS/=Zero)THEN
#ifdef PARALLEL
             Mssg=ProcessName(Proc)//'CPU (Sec,MFLOPS) = (' &
                  //TRIM(DblToMedmChar(Elapsed_CPUS))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_CPUS))) &
                  //'), WALL (Sec,MFLOPS) = ('                  &
                  //TRIM(DblToMedmChar(Elapsed_Wall))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall))) &
                  //'), NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(Proc)//'CPU (Sec,MFLOPS) = (' &
                  //TRIM(DblToMedmChar(Elapsed_CPUS))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_CPUS))) &
                  //'), WALL (Sec,MFLOPS) = ('                  &
                  //TRIM(DblToMedmChar(Elapsed_Wall))//', '     &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall))) &
                  //')'
#endif
          ELSE
#ifdef PARALLEL
             Mssg=ProcessName(Proc)//'WALL (Sec,MFLOPS) = (' &
                  //TRIM(DblToMedmChar(Elapsed_Wall))//', '      &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall)))  &
                  //'), NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(Proc)//'WALL (Sec,MFLOPS) = (' &
                  //TRIM(DblToMedmChar(Elapsed_Wall))//', '      &
                  //TRIM(IntToChar(MFlops(FLOPS,Elapsed_Wall)))  &
                  //')'
#endif
          ENDIF
       ELSE
          IF(Elapsed_CPUS>Zero)THEN
#ifdef PARALLEL
             Mssg=ProcessName(Proc)//'CPU Sec = '   &
                  //TRIM(DblToMedmChar(Elapsed_CPUS))   &
                  //', WALL (Sec) = '                   &
                  //TRIM(DblToMedmChar(Elapsed_Wall))   &
                  //', NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(Proc)//'CPU Sec = '  &
                  //TRIM(DblToMedmChar(Elapsed_CPUS))  &
                  //', WALL Sec = '                    &
                  //TRIM(DblToMedmChar(Elapsed_Wall))   
#endif
          ELSE
#ifdef PARALLEL
             Mssg=ProcessName(Proc)//'WALL Sec = '  &
                  //TRIM(DblToMedmChar(Elapsed_Wall))      &
                  //', NProc = '//TRIM(IntToChar(NPrc))
#else
             Mssg=ProcessName(Proc)//'WALL (Sec) = '  &
                  //TRIM(DblToMedmChar(Elapsed_Wall))
#endif
          ENDIF
       ENDIF
       WRITE(PU,*)TRIM(Mssg)
       CALL PrintProtectR(PU)
       CLOSE(Out)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Print_TIME
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------     
  FUNCTION MFlops(Flops,Sec)
    REAL(DOUBLE), INTENT(IN) :: Flops,Sec
    INTEGER                  :: MFlops
    IF(Sec<=Zero.OR.Flops<=Zero)THEN
       MFlops=0
    ELSE
       MFlops=INT(Flops*1.0D-6/Sec)         
    ENDIF
  END FUNCTION MFlops
!---------------------------------------------------------------------
!    PRINT MEMORY STATISTICS
!---------------------------------------------------------------------     
  SUBROUTINE Print_MEMS(A,Proc)
    TYPE(MEMS),INTENT(IN)       :: A
    CHARACTER(LEN=*),INTENT(IN) :: Proc
    INTEGER                     :: I,L,PU
    CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
    IF(PrintFlags%Key/=DEBUG_MAXIMUM)RETURN
    L=LEN(TRIM(Proc))
#ifdef PARALLEL
    IF(InParallel)THEN
       IF(MyId==ROOT)THEN
          PU=OpenPU()
          CALL PrintProtectL(PU)
          CLOSE(PU)
       ENDIF
       DO I=0,NPrc-1
          IF(InParallel)CALL AlignNodes()
          IF(MyId==I)THEN
             Mssg=TRIM(Proc)//'#'//TRIM(IntToChar(I))                    &
                  //' :: Allocs='//TRIM(IntToChar(A%Allocs))               &
                  //',  DeAllocs='//TRIM(IntToChar(A%DeAllocs))            &
                  //', '//TRIM(IntToChar(A%MemTab ))                       &
                  //' bytes are presently allocated.'//Rtrn                    &
                  //Blanks(1:L+3)//' A max of '//TRIM(IntToChar(A%MaxMem)) &
                  //' bytes were allocated.'
             PU=OpenPU() 
             WRITE(PU,*)Mssg 
             CLOSE(PU) 
          ENDIF
       ENDDO
       CALL AlignNodes()
       IF(MyId==ROOT)THEN
          PU=OpenPU()
          CALL PrintProtectR(PU)
          CLOSE(PU)
       ENDIF
    ELSE
#endif
       PU=OpenPU()
       CALL PrintProtectL(PU)
       Mssg=ProcessName(Proc)                                      &
            //'Allocs='//TRIM(IntToChar(A%Allocs))                   &
            //',  DeAllocs='//TRIM(IntToChar(A%DeAllocs))            &
            //', '//TRIM(IntToChar(A%MemTab ))                       &
            //' bytes are presently allocated.'//Rtrn                &
            //Blanks(1:20)//' A max of '//TRIM(IntToChar(A%MaxMem))  &
            //' bytes were allocated.'
       WRITE(PU,*)TRIM(Mssg)
       CALL PrintProtectR(PU)
       CLOSE(PU) 
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE Print_MEMS
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE PrintProtectL(Unit)
    INTEGER :: Unit
    IF(PrintFlags%Fmt==DEBUG_MMASTYLE.AND.Unit/=6) &
         WRITE(Unit,*)LeftParenStar
  END SUBROUTINE PrintProtectL
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
  SUBROUTINE PrintProtectR(Unit)
    INTEGER :: Unit
    IF(PrintFlags%Fmt==DEBUG_MMASTYLE.AND.Unit/=6) &
         WRITE(Unit,*)RightParenStar
  END SUBROUTINE PrintProtectR
!
!==================================================================
!
!    PLOT A BCSR MATRIX
!      
!==================================================================
!      FUNCTION Dot(N,V1,V2)
!         REAL(DOUBLE) :: Dot
!         REAL(DOUBLE), DIMENSION(:)  :: V1,V2
!         INTEGER :: I,N
!         Dot=Zero
!         DO I=1,N
!            Dot=Dot+V1(I)*V2(I)
!         ENDDO
!      END FUNCTION Dot
!
END MODULE 
