!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    HEURISTIC ORDERING VIA SPACE FILLING CURVES AND 
!    GRAPH THEORETICAL ENVELOPE REDUCTION VIA REVERSE CUTHILL MCKEE
!
MODULE Order
   USE DerivedTypes
   USE GlobalCharacters
   USE GlobalScalars
   USE GlobalObjects
   USE ProcessControl
   USE MemMan
   IMPLICIT NONE
   INTERFACE Sort
      MODULE PROCEDURE Sort_DBL_INT,Sort_INT_INT, Sort_INT_VECT
   END INTERFACE
   INTERFACE Random
      MODULE PROCEDURE RANDOM_INT,RANDOM_DBL
   END INTERFACE
   INTERFACE 
      SUBROUTINE DblIntSort77(N,X,Y,Ordr)
         USE DerivedTypes
         INTEGER,                  INTENT(IN)    :: N,Ordr
         REAL(DOUBLE),DIMENSION(N),INTENT(INOUT) :: X
         INTEGER,     DIMENSION(N),INTENT(INOUT) :: Y
      END SUBROUTINE
      SUBROUTINE SFCOrder77(N,R,Point,Key,Hilbert)
         USE DerivedTypes
         INTEGER,                     INTENT(IN)    :: N
         LOGICAL,                     INTENT(IN)    :: Hilbert
         REAL(DOUBLE), DIMENSION(3,N),INTENT(INOUT) :: R 
         INTEGER(INT8),DIMENSION(N),  INTENT(INOUT) :: Key 
         INTEGER,      DIMENSION(N),  INTENT(INOUT) :: Point
      END SUBROUTINE

      FUNCTION Interleave(Ix,Iy,Iz)       
         IMPLICIT NONE
         INTEGER,INTENT(IN) :: Ix,Iy,Iz
         INTEGER, PARAMETER :: INT8=SELECTED_INT_KIND(18)    
         INTEGER(INT8)            :: Interleave
       END FUNCTION Interleave
   END INTERFACE
   INTEGER, PARAMETER     :: SFC_NONE   =3480481
   INTEGER, PARAMETER     :: SFC_PEANO  =5308208
   INTEGER, PARAMETER     :: SFC_HILBERT=4808471
   INTEGER, PARAMETER     :: SFC_RANDOM =5058108
   CONTAINS
!
!      FUNCTION RANDOM_INT(Limits)
!         INTEGER                :: RANDOM_INT
!         INTEGER, DIMENSION(2)  :: Limits
!         REAL(DOUBLE)           :: Delta
!         REAL(DOUBLE), EXTERNAL :: Random
!         Delta=DBLE(Limits(2)-Limits(1)+1)
!         RANDOM_INT=Limits(1)+INT(Delta*Random())
!         RANDOM_INT=Limits(1)+INT(Delta*RAND())
!      END FUNCTION RANDOM_INT
!
      FUNCTION RANDOM_DBL(Limits)
         REAL(DOUBLE)              :: RANDOM_DBL
         REAL(DOUBLE),DIMENSION(2) :: Limits
         REAL(DOUBLE)              :: Delta
         REAL(DOUBLE), EXTERNAL    :: Random
         Delta=Limits(2)-Limits(1)+0.0D0
         RANDOM_DBL=Limits(1)+Delta*Random()
!         RANDOM_DBL=Limits(1)+Delta*Rand()
      END FUNCTION RANDOM_DBL
!
      FUNCTION RANDOM_INT(Limits)
         INTEGER               :: RANDOM_INT,Delta
         INTEGER, SAVE         :: JRan=10408
         INTEGER, DIMENSION(2) :: Limits
         INTEGER, PARAMETER    :: Im=259200,Ia=7141,Ic=54773
         JRan=MOD(JRan*Ia+Ic,Im)
         Delta=Limits(2)-Limits(1)+1
         RANDOM_INT=Limits(1)+(Delta*JRan)/Im
         IF(RANDOM_INT>Limits(2).OR.RANDOM_INT<Limits(1)) &
            CALL Halt(' Limits hosed in RANDOM_INT ')
      END FUNCTION RANDOM_INT
!--------------------------------------------------------------
!    F90 wrapper for SFCOrder77, which circumvents the lack
!    of INTEGER(KIND=8) (INTEGER*8) support for cheazy 
!    F90 compilers (pgf,nag...)
!
     SUBROUTINE SFCOrder(N,R,Point,SFC_KEY)
        INTEGER,        INTENT(IN)    :: N,SFC_KEY
        TYPE(DBL_RNK2), INTENT(INOUT) :: R
        TYPE(INT_VECT), INTENT(INOUT) :: Point
        INTEGER(INT8),ALLOCATABLE, &
                         DIMENSION(:) :: IKey
        TYPE(DBL_VECT)                :: RKey
        INTEGER                       :: I        
!
        IF(SFC_KEY==SFC_RANDOM)THEN
           CALL New(RKey,N)
           DO I=1,N
              Point%I(I)=I
              CALL RANDOM_NUMBER(RKey%D(I))
           ENDDO            
           CALL Sort_DBL_INT(RKey,Point,N)
           CALL Delete(RKey)           
        ELSEIF(SFC_KEY==SFC_PEANO)THEN
           ALLOCATE(IKey(N))
           CALL SFCOrder77(N,R%D,Point%I,IKey,.FALSE.)
           DEALLOCATE(IKey)
        ELSEIF(SFC_KEY==SFC_HILBERT)THEN
           ALLOCATE(IKey(N))
           CALL SFCOrder77(N,R%D,Point%I,IKey,.TRUE.)
           DEALLOCATE(IKey)
        ELSE
           CALL MondoHalt(-100,'Bad SFC_Key in SFCOrder')
        ENDIF
     END SUBROUTINE SFCOrder

     SUBROUTINE Sort_DBL_INT(X,Y,N_O,Ordr_O)
        TYPE(DBL_VECT), INTENT(INOUT) :: X
        TYPE(INT_VECT), INTENT(INOUT) :: Y
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-2
        N=MIN(SIZE(X%D),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL DblIntSort77(N,X%D,Y%I,Ordr)
    END SUBROUTINE Sort_DBL_INT


     SUBROUTINE Sort_INT_INT(X,Y,N_O,Ordr_O)
        TYPE(INT_VECT), INTENT(INOUT) :: X
        TYPE(INT_VECT), INTENT(INOUT) :: Y
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-1
        N=MIN(SIZE(X%I),SIZE(Y%I))
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntIntSort77(N,X%I,Y%I,Ordr)
    END SUBROUTINE Sort_INT_INT


     SUBROUTINE Sort_INT_VECT(X,N_O,Ordr_O)
        TYPE(INT_VECT), INTENT(INOUT) :: X
        INTEGER,OPTIONAL,INTENT(IN)   :: N_O,Ordr_O 
        INTEGER                       :: N,Ordr 
        Ordr=-1
        N=SIZE(X%I)
        IF(PRESENT(N_O))THEN
           IF(N_O>N)CALL Halt(' Dimensioning off in DblSort ')
           N=N_O
        ENDIF
        IF(PRESENT(Ordr_O))Ordr=Ordr_O
        CALL IntSort77(N,X%I,Ordr)
    END SUBROUTINE Sort_INT_VECT


END MODULE
!len=20;  dim=3;
!ix=1; iy=6; iz=0;
!l[1]=IntegerDigits[ix,2,len];
!l[2]=IntegerDigits[iy,2,len];
!l[3]=IntegerDigits[iz,2,len];
!n=Table[0,{i,dim*len}];
!Print[" Lx = ",lx];
!Print[" Ly = ",ly];
!Print[" Lz = ",lz];
!k=len;
!Do[
!	     Do[ 	 		
!		             n[[i-j+1]]=l[dim-j+1][[k]];		
!    				(*    
!      Print[" k  = ",k," ndex = ",i-j+1,", dim-j+1 = ",dim-j+1," l = ",
!        l[dim-j+1][[k]]];*)
!      		,{j,1,dim}];
!	    k=k-1;
!	,{i,dim*len,1,-dim}];
!Print[" N = ",n];
!Print[" Key = ",FromDigits[n,2]];
