!    GENERIC IO ROUTINES FOR MONDOSCF TYPES
!    Author: Matt Challacombe and CK Gan
!----------------------------------------------------------
MODULE InOut
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalCharacters
   USE ProcessControl    
   USE MemMan
   USE Indexing
   USE Parse
#ifdef PARALLEL
   USE MondoMPI
#endif   
   USE SetXYZ
   IMPLICIT NONE   
   INTERFACE Get
      MODULE PROCEDURE Get_INT_SCLR, Get_INT_VECT, Get_INT_RNK2, &
                       Get_INT_RNK3, Get_INT_RNK4, Get_DBL_SCLR, &
                       Get_DBL_VECT, Get_DBL_RNK2, Get_DBL_RNK3, &
                       Get_DBL_RNK4, Get_CHR_SCLR,               &
                       Get_LOG_SCLR, Get_BSET,     Get_CRDS,     &
                       Get_TOLS,     Get_BCSR,                   & 
#ifdef PARALLEL
                       Get_DBCSR,                                &
#endif
#ifdef PERIODIC
                       Get_PBCInfo,                              &
#endif  
                       Get_ARGMT,    Get_HGRho,                  &
#ifdef MMech
                       Get_CHR_VECT, Get_LOG_VECT,     &
                       Get_INTC,     Get_BMATR,        &
#endif
                       Get_CMPoles

   END INTERFACE     
   INTERFACE Put
      MODULE PROCEDURE Put_INT_SCLR, Put_INT_VECT, Put_INT_RNK2, &
                       Put_INT_RNK3, Put_INT_RNK4, Put_DBL_SCLR, &
                       Put_DBL_VECT, Put_DBL_RNK2, Put_DBL_RNK3, &
                       Put_DBL_RNK4, Put_DBL_RNK6, Put_CHR_SCLR, &
                       Put_LOG_SCLR, Put_BSET,     Put_CRDS,     &
#ifdef PARALLEL
                       Put_DBCSR,                                &
#endif
#ifdef PERIODIC
                       Put_PBCInfo,                              &
#endif 
                       Put_TOLS,     Put_BCSR,     Put_HGRho,    &
#ifdef MMech
                       Put_CHR_VECT, Put_LOG_VECT,        &
                       Put_INTC, Put_BMATR,               &
#endif
                       Put_CMPoles
   END INTERFACE     
   TYPE META_DATA
      INTEGER            :: Dimension
      INTEGER            :: DataType
      INTEGER            :: DataId
      INTEGER            :: DataSpc
      INTEGER            :: Status
      INTEGER            :: Unlimited
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: VarName 
   END TYPE
   INTEGER, SAVE      :: FileID
   INTEGER, EXTERNAL  :: HDF5CreateFile,HDF5OpenFile,HDF5CloseFile,         &
                         HDF5CloseData,HDF5ReadIntegerVector,               &
                         HDF5ReadDoubleVector,HDF5WriteIntegerVector,       &
                         HDF5WriteDoubleVector,HDF5SizeOfData
   INTEGER, PARAMETER :: NATIVE_DOUBLE=6, NATIVE_INT32=24
   CONTAINS 
!================================================================
!
!  FILE MANIPULATION ROUTINES
!
!================================================================
!
!     Initialize a HDF file
!
      SUBROUTINE InitHDF(FileName)
         CHARACTER(LEN=*), INTENT(IN) :: FileName
         INTEGER                      :: NC,STATUS
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            NC=StringLen(FileName)
            FileID=HDF5CreateFile(NC,Char2Ints(NC,FileName))
            IF(FileID==FAIL) &
            CALL Halt(' Failed to create the HDF file' &
                      //'<'//TRIM(FileName)//'>.') 
#ifdef PARALLEL
         ENDIF
#endif  
      END SUBROUTINE InitHDF



!=================================================================
!    Open a HDF file
!================================================================= 
      SUBROUTINE OpenHDF(FileName)
         CHARACTER(LEN=*),INTENT(IN) :: FileName
         INTEGER                     :: NC,STATUS

#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            NC=StringLen(FileName)



            FileID=HDF5OpenFile(NC,Char2Ints(NC,FileName))
            IF(FileID==FAIL) &
            CALL Halt(' Failed to open the HDF file' &
                      //'<'//TRIM(FileName)//'>.') 
#ifdef PARALLEL
         ENDIF
#endif  
      END SUBROUTINE OpenHDF


!-----------------------------------------------------------------------
!
!----------------------------------------------------------------------- 
      SUBROUTINE CloseHDF()
         INTEGER :: Status
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            STATUS=HDF5CloseFile(FileID)
            IF(FileID==FAIL) &
            CALL Halt(' Failed to close an HDF file with FileID=' &
                      //'<'//TRIM(IntToChar(FileID))//'>.') 
#ifdef PARALLEL
         ENDIF
#endif  
      END SUBROUTINE CloseHDF
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------  
      SUBROUTINE CreateData(Meta)
         TYPE(META_DATA) :: Meta
         INTEGER         :: NC
         NC=StringLen(Meta%VarName)
         CALL HDF5CreateData(FileID,Meta%DataType,Meta%Dimension,          &
                             NC,Char2Ints(NC,Meta%VarName),Meta%Unlimited, &
                             Meta%DataId,Meta%DataSpc)
         IF(Meta%DataId==FAIL) &
         CALL Halt(' Failed in CreateData:'//TRIM(MetaChar(Meta)))
         Meta%Status=SUCCEED
      END SUBROUTINE CreateData
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      SUBROUTINE OpenData(Meta,Put_O)
         TYPE(META_DATA)  :: Meta
         LOGICAL,OPTIONAL :: Put_O
         INTEGER          :: I,NC,SizeOf,SIZE_OF_HDF5_DATA
         Meta%Status=FAIL
         NC=StringLen(Meta%VarName)
         CALL HDF5OpenData(FileID,NC,Char2Ints(NC,Meta%VarName), &
                           Meta%DataId,Meta%DataSpc)
         IF(Meta%DataId==FAIL)THEN
            IF(PRESENT(Put_O))THEN
               CALL CreateData(Meta)
            ELSE
               CALL Halt(' Failed in OpenData:'//TRIM(MetaChar(Meta)))
            ENDIF
         ELSE
            SizeOf=HDF5SizeOfData(Meta%DataSpc)
            IF(PRESENT(Put_O))THEN
               IF(SizeOf<Meta%Dimension.AND.Meta%UnLimited==1)THEN
                  CALL HDF5ExtendData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
               ELSEIF(SizeOf>Meta%Dimension)THEN
                   CALL HDF5SelectData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
               ELSEIF(SizeOf<Meta%Dimension)THEN
                   CALL Halt(' Failed in OpenData, need to redimension fixed form data:' &
                            //TRIM(MetaChar(Meta)))
               ENDIF
            ELSEIF(SizeOf>Meta%Dimension)THEN
               CALL HDF5SelectData(Meta%DataId,Meta%DataSpc,Meta%Dimension)
            ENDIF
         ENDIF
      END SUBROUTINE OpenData
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      SUBROUTINE CloseData(Meta)
         TYPE(META_DATA) :: Meta
         IF(Meta%Status==FAIL) &
         CALL Halt(' HDF R/W error for '//TRIM(MetaChar(Meta)))
         Meta%Status=HDF5CloseData(Meta%DataId,Meta%DataSpc)
         IF(Meta%Status==FAIL) &
         CALL Halt(' HDF5CloseData error for '//TRIM(MetaChar(Meta)))
      END SUBROUTINE CloseData
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      FUNCTION NameTag(VarName,Tag_O) RESULT(FullName)
         CHARACTER(LEN=*),         INTENT(IN)    :: VarName
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
         CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FullName
         IF(PRESENT(Tag_O))THEN
            FullName=TRIM(VarName)//TRIM(Tag_O)
         ELSE
            FullName=TRIM(VarName)
         ENDIF
         CALL LowCase(FullName)
      END FUNCTION NameTag
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      FUNCTION SetMeta(VarName,DataType,SizeOf,UnLimit_O) RESULT(Meta)
         TYPE(META_DATA)              :: Meta
         CHARACTER(LEN=*), INTENT(IN) :: VarName
         INTEGER,          INTENT(IN) :: DataType,SizeOf
         LOGICAL,OPTIONAL, INTENT(IN) :: UnLimit_O
         IF(PRESENT(UnLimit_O))THEN
            IF(UnLimit_O)THEN
               Meta%UnLimited=1
            ELSE
               Meta%UnLimited=0
            ENDIF
         ELSE
            Meta%UnLimited=0
         ENDIF
         Meta%DataType=DataType
         Meta%Dimension=SizeOf
         Meta%VarName=VarName
         Meta%DataId=FAIL
         Meta%Status=FAIL
      END FUNCTION SetMeta
!-----------------------------------------------------------------------
! 
!-----------------------------------------------------------------------
      FUNCTION MetaChar(Meta) RESULT(CharMeta)
         TYPE(META_DATA)                :: Meta
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: CharMeta
         CHARACTER(LEN=1)               :: TF
         IF(Meta%UnLimited==1)THEN
            TF='T'
         ELSE
            TF='F'
         ENDIF
         CharMeta=Rtrn//'Meta%VarName  = <'//TRIM(Meta%VarName)//'>'       //Rtrn &
                      //'Meta%DataType = '//TRIM(IntToChar(Meta%DataType)) //Rtrn &
                      //'Meta%Dimension= '//TRIM(IntToChar(Meta%Dimension))//Rtrn &
                      //'Meta%DataId   = '//TRIM(IntToChar(Meta%DataId))   //Rtrn &
                      //'Meta%Status   = '//TRIM(IntToChar(Meta%Status))
      END FUNCTION MetaChar
!-----------------------------------------------------------------------

!=====================================================================
      FUNCTION StringLen(String)
         CHARACTER(LEN=*),INTENT(IN) :: String
         INTEGER                     :: StringLen
         StringLen=LEN(TRIM(String))
      END FUNCTION StringLen
!---------------------------------------------------------------------

 
!=====================================================================



      FUNCTION Char2Ints(NC,String)
         CHARACTER(LEN=*),INTENT(IN)          :: String
         INTEGER,         INTENT(IN)          :: NC 
         INTEGER,DIMENSION(NC)                :: Char2Ints
         INTEGER                              :: I 
         DO I=1,NC
            Char2Ints(I)=ICHAR(String(I:I))
         ENDDO
      END FUNCTION Char2Ints
!-----------------------------------------------------------------------
!  
!  
!  
      SUBROUTINE WriteIntegerVector(Meta,A)       
         TYPE(META_DATA)                :: Meta
         INTEGER,                   &
            DIMENSION(Meta%Dimension)   :: A
         Meta%Status=HDF5WriteIntegerVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE WriteIntegerVector
!
      SUBROUTINE WriteDoubleVector(Meta,A)       
         TYPE(META_DATA)                :: Meta
         REAL(DOUBLE),              &
            DIMENSION(Meta%Dimension)   :: A
         Meta%Status=HDF5WriteDoubleVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE WriteDoubleVector
!
      SUBROUTINE ReadIntegerVector(Meta,A)       
         TYPE(META_DATA)               :: Meta
         INTEGER,                    &
           DIMENSION(Meta%Dimension)   :: A
         Meta%Status=HDF5ReadIntegerVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE ReadIntegerVector
!
      SUBROUTINE ReadDoubleVector(Meta,A)       
         TYPE(META_DATA)               :: Meta
         REAL(DOUBLE),              &
            DIMENSION(Meta%Dimension)  :: A
         Meta%Status=HDF5ReadDoubleVector(Meta%DataId,Meta%DataSpc,A)
      END SUBROUTINE ReadDoubleVector
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_SCLR(A,VarName,Tag_O)
       INTEGER,                  INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       INTEGER,DIMENSION(1)                    :: B
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,B)
          CALL CloseData(Meta)
          A=B(1)
#ifdef PARALLEL 
       ENDIF       
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_SCLR
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_VECT(A,VarName,Tag_O)
       TYPE(INT_VECT),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       SIZE(A%I,1),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_VECT
!--------------------------------------------------------------------
#ifdef MMech
!
    SUBROUTINE Get_INTC(A,VarName,Tag_O)
       TYPE(INTC),               INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       INTEGER,ALLOCATABLE,DIMENSION(:)   :: B !=ICHAR(' ')
       TYPE(META_DATA)                         :: Meta
       INTEGER :: II,NN,N,I
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
        NN=SIZE(A%DEF,1)
         N=5  !!! LEN(A%DEF(1))
          Meta=SetMeta(NameTag(Trim(VarName)//'Def',Tag_O),NATIVE_INT32,NN*N,.FALSE.)
          ALLOCATE(B(NN*N))
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,B)
          CALL CloseData(Meta)
        DO II = 1, NN
          DO I=1,N; A%DEF(II)(I:I)=CHAR(B((II-1)*N+I)); ENDDO
        ENDDO
          DEALLOCATE(B)
!
          Meta=SetMeta(NameTag(Trim(VarName)//'Atoms',Tag_O),NATIVE_INT32,NN*4,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%ATOMS)
          CALL CloseData(Meta)
!
          Meta=SetMeta(NameTag(Trim(VarName)//'Value',Tag_O),NATIVE_DOUBLE,NN,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%Value)
          CALL CloseData(Meta)
!
          Meta=SetMeta(NameTag(Trim(VarName)//'Constraint',Tag_O),NATIVE_INT32,NN,.FALSE.)
          ALLOCATE(B(NN))
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,B)
          CALL CloseData(Meta)
              A%Constraint=.FALSE.
            DO I=1,NN
              IF(B(I)==1) A%Constraint(I)=.TRUE.
            ENDDO
          DEALLOCATE(B)
!
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INTC
!
!------------------------------------------------------------
!
    SUBROUTINE Put_INTC(A,VarName,Tag_O)
       INTEGER  :: I,N,II,NN
       TYPE(INTC),               INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,ALLOCATABLE,DIMENSION(:)   :: B  !=ICHAR(' ')
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          N=5 !!! LEN(A%DEF(1))
          NN=SIZE(A%DEF,1)
          IF(N>DEFAULT_CHR_LEN) CALL Halt('Static strings overrun in Put_INTC')
          ALLOCATE(B(N*NN))
       DO II = 1, NN
         DO I=1,N; B((II-1)*N+I)=ICHAR(A%DEF(II)(I:I)); ENDDO
       ENDDO
          Meta=SetMeta(NameTag(TRIM(VarName)//'Def',Tag_O),NATIVE_INT32,NN*N,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,B)
          CALL CloseData(Meta)
          DEALLOCATE(B)
!
          Meta=SetMeta(NameTag(TRIM(VarName)//'Atoms',Tag_O),NATIVE_INT32,4*NN,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,A%ATOMS)
          CALL CloseData(Meta)
!
          Meta=SetMeta(NameTag(TRIM(VarName)//'Value',Tag_O),NATIVE_DOUBLE,NN,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteDoubleVector(Meta,A%Value)
          CALL CloseData(Meta)
!
          ALLOCATE(B(NN))
          B=0
          DO I=1,NN
            IF(A%Constraint(I)) B(I)=1
          ENDDO
          Meta=SetMeta(NameTag(TRIM(VarName)//'Constraint',Tag_O),NATIVE_INT32,NN,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,B)
          CALL CloseData(Meta)
          DEALLOCATE(B)
!
#ifdef PARALLEL 
       ENDIF       
#endif 
!
    END SUBROUTINE Put_INTC
!
    SUBROUTINE Get_BMATR(A,VarName,Tag_O)
       TYPE(BMATR),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(TRIM(VarName)//'IB',Tag_O),NATIVE_INT32,SIZE(A%IB,1)*SIZE(A%IB,2),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%IB)
          CALL CloseData(Meta)
!
          Meta=SetMeta(NameTag(TRIM(VarName)//'B',Tag_O),NATIVE_DOUBLE,SIZE(A%B,1)*SIZE(A%B,2),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%B)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_BMATR
!
    SUBROUTINE Put_BMATR(A,VarName,Tag_O)
       TYPE(BMATR),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(TRIM(VarName)//'IB',Tag_O),NATIVE_INT32,SIZE(A%IB,1)*SIZE(A%IB,2),.FALSE.)
          CALL OpenData(Meta)
          CALL WriteIntegerVector(Meta,A%IB)
          CALL CloseData(Meta)
!
          Meta=SetMeta(NameTag(TRIM(VarName)//'B',Tag_O),NATIVE_DOUBLE,SIZE(A%B,1)*SIZE(A%B,2),.FALSE.)
          CALL OpenData(Meta)
          CALL WriteDoubleVector(Meta,A%B)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
#endif 
    END SUBROUTINE Put_BMATR
!
#endif
!
!--------------------------------------------------------------------
!
!
!----------------------------------------------------------
!
    SUBROUTINE Get_INT_RNK2(A,VarName,Tag_O)
       TYPE(INT_RNK2),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       SIZE(A%I,1)*SIZE(A%I,2),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_RNK2
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_RNK3(A,VarName,Tag_O)
       TYPE(INT_RNK3),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3), &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_RNK3
!--------------------------------------------------------------------
!
!
    SUBROUTINE Get_INT_RNK4(A,VarName,Tag_O)
       TYPE(INT_RNK4),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,             &
                       SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)*SIZE(A%I,4), &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_INT_RNK4
!
    SUBROUTINE Get_DBL_SCLR(A,VarName,Tag_O)
       REAL(DOUBLE),             INTENT(INOUT)   :: A
       CHARACTER(LEN=*),         INTENT(IN)      :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)      :: Tag_O
       REAL(DOUBLE),DIMENSION(1)                 :: B
       TYPE(META_DATA)                           :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,B)
          CALL CloseData(Meta)
          A=B(1)
#ifdef PARALLEL 
       ENDIF       
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_SCLR
!
    SUBROUTINE Get_DBL_VECT(A,VarName,Tag_O)
       TYPE(DBL_VECT),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                       SIZE(A%D,1),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_VECT
!
    SUBROUTINE Get_DBL_RNK2(A,VarName,Tag_O)
       TYPE(DBL_RNK2),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                       SIZE(A%D,1)*SIZE(A%D,2),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK2
!
    SUBROUTINE Get_DBL_RNK3(A,VarName,Tag_O)
       TYPE(DBL_RNK3),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE, &
                       SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3),  &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK3
!
    SUBROUTINE Get_DBL_RNK4(A,VarName,Tag_O)
       TYPE(DBL_RNK4),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,            &
                       SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4), &
                       .FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK4
!
    SUBROUTINE Get_DBL_RNK6(A,VarName,Tag_O)
       TYPE(DBL_RNK6),           INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT) THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,            &
                       SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)  &
                      *SIZE(A%D,5)*SIZE(A%D,6),.FALSE.)
          CALL OpenData(Meta)
          CALL ReadDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_DBL_RNK6
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
    SUBROUTINE Get_CHR_SCLR(A,VarName,Tag_O)
       CHARACTER(LEN=*),         INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       INTEGER,DIMENSION(DEFAULT_CHR_LEN)      :: B=ICHAR(' ')
       INTEGER                                 :: I,N
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          N=LEN(A)
          IF(N>DEFAULT_CHR_LEN) &
             CALL Halt('Static strings overrun in Get_CHR_SCLR')
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       DEFAULT_CHR_LEN,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,B)
          DO I=1,N; A(I:I)=CHAR(B(I)); ENDDO
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_CHR_SCLR
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
    SUBROUTINE Get_LOG_SCLR(A,VarName,Tag_O)
       LOGICAL,                  INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       INTEGER,DIMENSION(1)                    :: ILog
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta)
          CALL ReadIntegerVector(Meta,ILog)
          CALL CloseData(Meta)
          IF(ILog(1)==0)THEN
             A=.FALSE.
          ELSEIF(ILog(1)==1)THEN
             A=.TRUE.
          ELSE
             CALL Halt('Error in Get_LOG_SCLR.')
          ENDIF
#ifdef PARALLEL 
       ENDIF       
       IF(InParallel)CALL BCast(A)
#endif 
    END SUBROUTINE Get_LOG_SCLR
!---------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------
    SUBROUTINE Put_INT_SCLR(A,VarName,Tag_O)
       INTEGER,                  INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,DIMENSION(1)                 :: B
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          B(1)=A
          CALL WriteIntegerVector(Meta,B)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_INT_SCLR
!
    SUBROUTINE Put_INT_VECT(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(INT_VECT),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%I,1)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_INT_VECT
!
    SUBROUTINE Put_INT_RNK2(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(INT_RNK2),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%I,1)*SIZE(A%I,2)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_INT_RNK2
!
    SUBROUTINE Put_INT_RNK3(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(INT_RNK3),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_INT_RNK3
!
    SUBROUTINE Put_INT_RNK4(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(INT_RNK4),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%I,1)*SIZE(A%I,2)*SIZE(A%I,3)*SIZE(A%I,4)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,A%I)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_INT_RNK4
!---------------------------------------------------------------------------------------------
!
!---------------------------------------------------------------------------------------------
    SUBROUTINE Put_DBL_SCLR(A,VarName,Tag_O)
       REAL(DOUBLE),             INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       REAL(DOUBLE),DIMENSION(1)            :: B
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,1,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          B(1)=A
          CALL WriteDoubleVector(Meta,B)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_DBL_SCLR
!
    SUBROUTINE Put_DBL_VECT(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(DBL_VECT),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%D,1)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_DBL_VECT
! 
    SUBROUTINE Put_DBL_RNK2(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(DBL_RNK2),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%D,1)*SIZE(A%D,2)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_DBL_RNK2
! 
    SUBROUTINE Put_DBL_RNK3(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(DBL_RNK3),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_DBL_RNK3
! 
    SUBROUTINE Put_DBL_RNK4(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(DBL_RNK4),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_DBL_RNK4
! 
    SUBROUTINE Put_DBL_RNK6(A,VarName,N_O,Tag_O,UnLimit_O)
       TYPE(DBL_RNK6),           INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       INTEGER,         OPTIONAL,INTENT(IN) :: N_O
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       LOGICAL,         OPTIONAL,INTENT(IN) :: UnLimit_O
       INTEGER                              :: N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          IF(PRESENT(N_O))THEN
             N=N_O
          ELSE
             N=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3) &
              *SIZE(A%D,4)*SIZE(A%D,5)*SIZE(A%D,6)
          ENDIF
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_DOUBLE,N,UnLimit_O)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteDoubleVector(Meta,A%D)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_DBL_RNK6
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
    SUBROUTINE Put_CHR_SCLR(A,VarName,Tag_O)
       CHARACTER(LEN=*),         INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,DIMENSION(DEFAULT_CHR_LEN)   :: B=ICHAR(' ')
       INTEGER                              :: I,N
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          N=LEN(A)
          IF(N>DEFAULT_CHR_LEN)THEN
             CALL Halt('Static strings overrun in Put_CHR_SCLR')
          ENDIF
          DO I=1,N; B(I)=ICHAR(A(I:I)); ENDDO
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       DEFAULT_CHR_LEN,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,B)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_CHR_SCLR
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
    SUBROUTINE Put_LOG_SCLR(A,VarName,Tag_O)
       LOGICAL,                  INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,DIMENSION(1)                 :: ILog
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          ILog(1)=0
          IF(A)ILog(1)=1
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,ILog)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_LOG_SCLR
!----------------------------------------------------------------------------------------
!     Get a basis set
!----------------------------------------------------------------------------------------
      SUBROUTINE Get_BSET(BS,Tag_O)
         TYPE(BSET),              INTENT(OUT) :: BS
         CHARACTER(LEN=*),OPTIONAL            :: Tag_O
         IF(AllocQ(BS%Alloc))CALL Delete(BS)         
         CALL Get(BS%NAtms,'natoms')
         CALL Get(BS%BName,'bsetname',Tag_O=Tag_O)
         CALL Get(BS%NBasF,'nbasf',Tag_O=Tag_O)
         CALL Get(BS%NKind,'nkind',Tag_O=Tag_O)
         CALL Get(BS%NCtrt,'nctrt',Tag_O=Tag_O)
         CALL Get(BS%NPrim,'nprim',Tag_O=Tag_O)
         CALL Get(BS%NASym,'nasym',Tag_O=Tag_O)
         CALL Get(BS%LMNLen,'lmnlen',Tag_O=Tag_O)
         CALL New(BS)
         CALL Get(BS%Kinds,'kind',Tag_O=Tag_O)
         CALL Get(BS%NCFnc,'ncfuncs',Tag_O=Tag_O)
         CALL Get(BS%NPFnc,'npfuncs',Tag_O=Tag_O)
         CALL Get(BS%ASymm,'angsym',Tag_O=Tag_O)
         CALL Get(BS%Expnt,'exponent',Tag_O=Tag_O)
         CALL Get(BS%CCoef,'ccoefcnt',Tag_O=Tag_O)
         CALL Get(BS%BFKnd,'basfperkind',Tag_O=Tag_O)
         CALL Get(BS%LStrt,'lstart',Tag_O=Tag_O)
         CALL Get(BS%LStop,'lstop',Tag_O=Tag_O)
         CALL Get(BS%LxDex,'lxdex',Tag_O=Tag_O)
         CALL Get(BS%LyDex,'lydex',Tag_O=Tag_O)
         CALL Get(BS%LzDex,'lzdex',Tag_O=Tag_O)
      END SUBROUTINE Get_BSET
!------------------------------------------------------------------
!     Put a  basis set
!  
      SUBROUTINE Put_BSET(BS,Tag_O)
         TYPE(BSET),               INTENT(IN) :: BS
         CHARACTER(LEN=*),OPTIONAL            :: Tag_O
         CALL Put(BS%NAtms,'natoms')
         CALL Put(BS%BName,'bsetname',Tag_O=Tag_O)
         CALL Put(BS%NBasF,'nbasf',Tag_O=Tag_O)
         CALL Put(BS%NKind,'nkind',Tag_O=Tag_O)
         CALL Put(BS%NCtrt,'nctrt',Tag_O=Tag_O)
         CALL Put(BS%NPrim,'nprim',Tag_O=Tag_O)
         CALL Put(BS%NASym,'nasym',Tag_O=Tag_O)
         CALL Put(BS%LMNLen,'lmnlen',Tag_O=Tag_O)
         CALL Put(BS%Kinds,'kind',Tag_O=Tag_O)
         CALL Put(BS%BFKnd,'basfperkind',Tag_O=Tag_O)
         CALL Put(BS%NCFnc,'ncfuncs',Tag_O=Tag_O)
         CALL Put(BS%NPFnc,'npfuncs',Tag_O=Tag_O)
         CALL Put(BS%ASymm,'angsym',Tag_O=Tag_O)
         CALL Put(BS%Expnt,'exponent',Tag_O=Tag_O)
         CALL Put(BS%CCoef,'ccoefcnt',Tag_O=Tag_O)
         CALL Put(BS%LStrt,'lstart',Tag_O=Tag_O)
         CALL Put(BS%LStop,'lstop',Tag_O=Tag_O)
         CALL Put(BS%LxDex,'lxdex',Tag_O=Tag_O)
         CALL Put(BS%LyDex,'lydex',Tag_O=Tag_O)
         CALL Put(BS%LzDex,'lzdex',Tag_O=Tag_O)
      END SUBROUTINE Put_BSET
#ifdef PERIODIC
!-------------------------------------------------------------------
!     Get the Periodic Info
!
      SUBROUTINE Get_PBCInfo(PBC,Tag_O)
        TYPE(PBCInfo),            INTENT(OUT) :: PBC 
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: Tag_O
!-------------------------------------------------------------------
!       Items that should not change with geometry
!
        CALL Get(PBC%Dimen     ,'Dimension')
        CALL Get(PBC%AtomW     ,'AtomWrap')
        CALL Get(PBC%InVecForm ,'VectorForm')
        CALL Get(PBC%InAtomCrd ,'AtomicCrd')
        CALL Get(PBC%Translate ,'Translate')
        CALL Get(PBC%Trans_COM ,'Trans_COM')
!
        CALL Get(PBC%AutoW(1)  ,'AutoWrap(1)')        
        CALL Get(PBC%AutoW(2)  ,'AutoWrap(2)')
        CALL Get(PBC%AutoW(3)  ,'AutoWrap(3)')
!----------------------------------------------------------------
!       Items that can change with geometry ...   
!
        CALL Get(PBC%CellVolume,   'CellVolume'          ,Tag_O=Tag_O)
        CALL Get(PBC%Epsilon,      'Epsilon'             ,Tag_O=Tag_O)
        CALL Get(PBC%DipoleFAC,    'DPoleFAC'            ,Tag_O=Tag_O)        
        CALL Get(PBC%QupoleFAC,    'QPoleFAC'            ,Tag_O=Tag_O)
!
        CALL Get(PBC%CellCenter(1),'CellCenter(1)'       ,Tag_O=Tag_O)
        CALL Get(PBC%CellCenter(2),'CellCenter(2)'       ,Tag_O=Tag_O)
        CALL Get(PBC%CellCenter(3),'CellCenter(3)'       ,Tag_O=Tag_O)
!
        CALL Get(PBC%TransVec(1),  'Originvector(1)'     ,Tag_O=Tag_O)
        CALL Get(PBC%TransVec(2),  'Originvector(2)'     ,Tag_O=Tag_O)
        CALL Get(PBC%TransVec(3),  'Originvector(3)'     ,Tag_O=Tag_O)
!
        CALL Get(PBC%BoxShape(1,1),'Boxshape(1,1)'       ,Tag_O=Tag_O)
        CALL Get(PBC%BoxShape(1,2),'Boxshape(1,2)'       ,Tag_O=Tag_O)
        CALL Get(PBC%BoxShape(1,3),'Boxshape(1,3)'       ,Tag_O=Tag_O)       
        CALL Get(PBC%BoxShape(2,1),'Boxshape(2,1)'       ,Tag_O=Tag_O)       
        CALL Get(PBC%BoxShape(2,2),'Boxshape(2,2)'       ,Tag_O=Tag_O)
        CALL Get(PBC%BoxShape(2,3),'Boxshape(2,3)'       ,Tag_O=Tag_O)
        CALL Get(PBC%BoxShape(3,1),'Boxshape(3,1)'       ,Tag_O=Tag_O)
        CALL Get(PBC%BoxShape(3,2),'Boxshape(3,2)'       ,Tag_O=Tag_O)
        CALL Get(PBC%BoxShape(3,3),'Boxshape(3,3)'       ,Tag_O=Tag_O)
!
        CALL Get(PBC%InvBoxSh(1,1),'InverseBoxshape(1,1)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(1,2),'InverseBoxshape(1,2)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(1,3),'InverseBoxshape(1,3)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(2,1),'InverseBoxshape(2,1)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(2,2),'InverseBoxshape(2,2)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(2,3),'InverseBoxshape(2,3)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(3,1),'InverseBoxshape(3,1)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(3,2),'InverseBoxshape(3,2)',Tag_O=Tag_O)
        CALL Get(PBC%InvBoxSh(3,3),'InverseBoxshape(3,3)',Tag_O=Tag_O)
!
      END SUBROUTINE Get_PBCInfo
!-------------------------------------------------------------------
!     Put the Periodic Info
!
      SUBROUTINE Put_PBCInfo(PBC,Tag_O)
        TYPE(PBCInfo),            INTENT(IN)  :: PBC 
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN)  :: Tag_O
!-------------------------------------------------------------------
!       Items that should not change with geometry
!
        CALL Put(PBC%Dimen     ,'Dimension')
        CALL Put(PBC%AtomW     ,'AtomWrap')
        CALL Put(PBC%InVecForm ,'VectorForm')
        CALL Put(PBC%InAtomCrd ,'AtomicCrd')
        CALL Put(PBC%Translate ,'Translate')
        CALL Put(PBC%Trans_COM ,'Trans_COM')
!
        CALL Put(PBC%AutoW(1)  ,'AutoWrap(1)')        
        CALL Put(PBC%AutoW(2)  ,'AutoWrap(2)')
        CALL Put(PBC%AutoW(3)  ,'AutoWrap(3)')
!----------------------------------------------------------------
!       Items that can change with geometry ...   
!
        CALL Put(PBC%CellVolume,   'CellVolume'          ,Tag_O=Tag_O)
        CALL Put(PBC%Epsilon,      'Epsilon'             ,Tag_O=Tag_O)
        CALL Put(PBC%DipoleFAC,    'DPoleFAC'            ,Tag_O=Tag_O)        
        CALL Put(PBC%QupoleFAC,    'QPoleFAC'            ,Tag_O=Tag_O)
!
        CALL Put(PBC%CellCenter(1),'CellCenter(1)'       ,Tag_O=Tag_O)
        CALL Put(PBC%CellCenter(2),'CellCenter(2)'       ,Tag_O=Tag_O)
        CALL Put(PBC%CellCenter(3),'CellCenter(3)'       ,Tag_O=Tag_O)
!
        CALL Put(PBC%TransVec(1),  'Originvector(1)'     ,Tag_O=Tag_O)
        CALL Put(PBC%TransVec(2),  'Originvector(2)'     ,Tag_O=Tag_O)
        CALL Put(PBC%TransVec(3),  'Originvector(3)'     ,Tag_O=Tag_O)
!
        CALL Put(PBC%BoxShape(1,1),'Boxshape(1,1)'       ,Tag_O=Tag_O)
        CALL Put(PBC%BoxShape(1,2),'Boxshape(1,2)'       ,Tag_O=Tag_O)
        CALL Put(PBC%BoxShape(1,3),'Boxshape(1,3)'       ,Tag_O=Tag_O)       
        CALL Put(PBC%BoxShape(2,1),'Boxshape(2,1)'       ,Tag_O=Tag_O)       
        CALL Put(PBC%BoxShape(2,2),'Boxshape(2,2)'       ,Tag_O=Tag_O)
        CALL Put(PBC%BoxShape(2,3),'Boxshape(2,3)'       ,Tag_O=Tag_O)
        CALL Put(PBC%BoxShape(3,1),'Boxshape(3,1)'       ,Tag_O=Tag_O)
        CALL Put(PBC%BoxShape(3,2),'Boxshape(3,2)'       ,Tag_O=Tag_O)
        CALL Put(PBC%BoxShape(3,3),'Boxshape(3,3)'       ,Tag_O=Tag_O)
!
        CALL Put(PBC%InvBoxSh(1,1),'InverseBoxshape(1,1)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(1,2),'InverseBoxshape(1,2)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(1,3),'InverseBoxshape(1,3)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(2,1),'InverseBoxshape(2,1)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(2,2),'InverseBoxshape(2,2)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(2,3),'InverseBoxshape(2,3)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(3,1),'InverseBoxshape(3,1)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(3,2),'InverseBoxshape(3,2)',Tag_O=Tag_O)
        CALL Put(PBC%InvBoxSh(3,3),'InverseBoxshape(3,3)',Tag_O=Tag_O)
!
      END SUBROUTINE Put_PBCInfo
#endif
!-------------------------------------------------------------------
!     Get some coordinates
!
      SUBROUTINE Get_CRDS(GM,Tag_O)
         TYPE(CRDS),           INTENT(OUT)    :: GM
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
         IF(AllocQ(GM%Alloc))CALL Delete(GM)
!-------------------------------------------------------------------
!        Items that should not change with geometry...
!
         CALL Get(GM%NAtms,'natoms',Tag_O=Tag_O)
         CALL Get(GM%Confg,'configuration',Tag_O=Tag_O)
         CALL Get(GM%NElec,'nel',Tag_O=Tag_O)
         CALL Get(GM%NAlph,'nelalpha',Tag_O=Tag_O)
         CALL Get(GM%NBeta,'nelbeta',Tag_O=Tag_O)
         CALL Get(GM%TotCh,'charge',Tag_O=Tag_O)
         CALL Get(GM%NKind,'nkind',Tag_O=Tag_O)
         CALL Get(GM%InAu, 'inau',Tag_O=Tag_O)
         CALL New(GM)
!----------------------------------------------------------------
!        Items that can change with geometry ...       
!
         CALL Get(GM%Ordrd,'reordered',Tag_O=Tag_O)
         CALL Get(GM%AtTyp,'atomtype',Tag_O=Tag_O)
         CALL Get(GM%AtNum,'atomicnumbers',Tag_O=Tag_O)
         CALL Get(GM%AtMss,'atomicmass',   Tag_O=Tag_O)
         CALL Get(GM%Carts,'cartesians',Tag_O=Tag_O)
         CALL Get(GM%Vects,'velocities',Tag_O=Tag_O)
         CALL Get(GM%BndBox,'boundingbox',Tag_O=Tag_O)
#ifdef PERIODIC
         CALL Get(GM%PBC,Tag_O=Tag_O)
         CALL Get(GM%BoxCarts,'LatticeCoord',Tag_O=Tag_O)
         CALL Get(GM%BoxVects,'LatticeVeloc',Tag_O=Tag_O)
#endif
      END SUBROUTINE Get_CRDS
!---------------------------------------------------------------------
!     Put a coordinate set
!
      SUBROUTINE Put_CRDS(GM,Tag_O)
         TYPE(CRDS),               INTENT(IN) :: GM
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
!-------------------------------------------------------------------
!        Items that should not change with geometry...
!
         CALL Put(GM%NAtms,'natoms',Tag_O=Tag_O)
         CALL Put(GM%Confg,'configuration',Tag_O=Tag_O)
         CALL Put(GM%NElec,'nel',Tag_O=Tag_O)
         CALL Put(GM%NAlph,'nelalpha',Tag_O=Tag_O)
         CALL Put(GM%NBeta,'nelbeta',Tag_O=Tag_O)
         CALL Put(GM%TotCh,'charge',Tag_O=Tag_O)
         CALL Put(GM%NKind,'nkind',Tag_O=Tag_O)
         CALL Put(GM%InAu, 'inau',Tag_O=Tag_O)
!----------------------------------------------------------------
!        Items that can change with geometry ...       
!
         CALL Put(GM%Ordrd,'reordered',Tag_O=Tag_O)
         CALL Put(GM%AtNum,'atomicnumbers',Tag_O=Tag_O)
         CALL Put(GM%AtMss,'atomicmass',   Tag_O=Tag_O)
         CALL Put(GM%AtTyp,'atomtype',Tag_O=Tag_O)
         CALL Put(GM%Carts,'cartesians',Tag_O=Tag_O)
         CALL Put(GM%Vects,'velocities',Tag_O=Tag_O)
         CALL Put(GM%BndBox,'boundingbox',Tag_O=Tag_O)
#ifdef PERIODIC
         CALL Put(GM%PBC,Tag_O=Tag_O)
         CALL Put(GM%BoxCarts,'LatticeCoord',Tag_O=Tag_O)
         CALL Put(GM%BoxVects,'LatticeVeloc',Tag_O=Tag_O)
#endif
      END SUBROUTINE Put_CRDS
!---------------------------------------------------------------------
!     Get a BCSR matrix
!
      SUBROUTINE Get_BCSR(A,Name,PFix_O,CheckPoint_O)
         TYPE(BCSR),               INTENT(INOUT) :: A     
         CHARACTER(LEN=*),         INTENT(IN)    :: Name
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
         LOGICAL,         OPTIONAL,INTENT(IN)    :: CheckPoint_O
         REAL(DOUBLE)                            :: Dummy
         CHARACTER(LEN=DEFAULT_CHR_LEN)          :: FileName          
         INTEGER                                 :: I,NAtms,NBlks,NNon0,IOS
         LOGICAL                                 :: Exists,LimitsQ
#ifdef PARALLEL
         IF(MyId==0)THEN
#endif
            IF(PRESENT(CheckPoint_O))THEN
               IF(CheckPoint_O)THEN
                  CALL Get(A%NAtms,TRIM(Name)//'%NAtms')
                  CALL Get(A%NBlks,TRIM(Name)//'%NBlks')
                  CALL Get(A%NNon0,TRIM(Name)//'%NNon0')
                  IF(AllocQ(A%Alloc))THEN
                     LimitsQ=.NOT.                   &
                       (A%NAtms<=SIZE(A%RowPt%I)).AND. &
                       (A%NBlks<=SIZE(A%ColPt%I)).AND. &
                       (A%NBlks<=SIZE(A%BlkPt%I)).AND. &
                       (A%NNon0<=SIZE(A%MTrix%D))
                     IF(LimitsQ)THEN
                        CALL Delete(A)
                        CALL New(A,(/A%NAtms,A%NBlks,A%NNon0/))
                     ENDIF
                  ELSE
                     CALL New(A,(/A%NAtms,A%NBlks,A%NNon0/))
                  ENDIF
                  CALL Get(A%RowPt,TRIM(Name)//'%RowPt')
                  CALL Get(A%ColPt,TRIM(Name)//'%ColPt')
                  CALL Get(A%BlkPt,TRIM(Name)//'%BlkPt')
                  CALL Get(A%MTrix,TRIM(Name)//'%MTrix')
                  RETURN
               ENDIF
            ENDIF
!
            IF(PRESENT(PFix_O))THEN
               FileName=TRIM(Name)//TRIM(PFix_O)
            ELSE
               FileName=Name
            ENDIF
            INQUIRE(FILE=TRIM(FileName),EXIST=Exists)
            IF(Exists)THEN
#ifdef FORMATTED              
              OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD', &
                   FORM='FORMATTED',ACCESS='SEQUENTIAL')
#else
              OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD', &
                   FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
            ELSE
              CALL Halt(' Get_BCSR could not find '//TRIM(FileName))
            ENDIF

#ifdef FORMATTED
            READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)NAtms,NNon0,NBlks
            INCLUDE 'Formats.Inc'
#else
            READ(UNIT=Seq,Err=1,IOSTAT=IOS)NAtms,NNon0,NBlks
#endif
            IF(AllocQ(A%Alloc))THEN
               LimitsQ=.NOT.                         &
                       (NAtms<=SIZE(A%RowPt%I)).AND. &
                       (NBlks<=SIZE(A%ColPt%I)).AND. &
                       (NBlks<=SIZE(A%BlkPt%I)).AND. &
                       (NNon0<=SIZE(A%MTrix%D))
               IF(LimitsQ)THEN
                  CALL Delete(A)
                  CALL New(A,(/NAtoms,NBlks,NNon0/))
               ELSE
                  A%NAtms=NAtms
                  A%NBlks=NBlks
                  A%NNon0=NNon0
               ENDIF
            ELSE
               CALL New(A,(/NAtms,NBlks,NNon0/))
            ENDIF
#ifdef FORMATTED
            READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
            READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
            READ(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
            READ(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(Dummy,i=1,A%NBlks)
            READ(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#else
            READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(I),I=1,NAtoms+1)
            READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(I),I=1,A%NBlks)
            READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(I),I=1,A%NBlks)
            READ(UNIT=Seq,Err=1,IOSTAT=IOS)(Dummy,i=1,A%NBlks)
            READ(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(I),I=1,A%NNon0)
#endif

            CLOSE(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
         ENDIF
#endif
         RETURN
       1 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_BCSR.')
      END SUBROUTINE Get_BCSR
!---------------------------------------------------------------------
!     Put a BCSR matrix
!
      SUBROUTINE Put_BCSR(A,Name,PFix_O,CheckPoint_O)
         TYPE(BCSR),               INTENT(IN) :: A             
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O
         CHARACTER(LEN=*),         INTENT(IN) :: Name
         LOGICAL,         OPTIONAL,INTENT(IN) :: CheckPoint_O
         CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FileName
         LOGICAL                              :: Exists
         INTEGER                              :: I,IOS
#ifdef PARALLEL
         IF(MyId==0)THEN
#endif
            IF(PRESENT(PFix_O))THEN
               FileName=TRIM(Name)//TRIM(PFix_O)
            ELSE
               FileName=Name
            ENDIF
!
            IF(PRESENT(CheckPoint_O))THEN
               IF(CheckPoint_O)THEN
                  CALL Put(A%NAtms,TRIM(Name)//'%NAtms')
                  CALL Put(A%NBlks,TRIM(Name)//'%NBlks')
                  CALL Put(A%NNon0,TRIM(Name)//'%NNon0')
                  CALL Put(A%RowPt,TRIM(Name)//'%RowPt',A%NAtms+1)
                  CALL Put(A%ColPt,TRIM(Name)//'%ColPt',A%NBlks,UnLimit_O=.TRUE.)
                  CALL Put(A%BlkPt,TRIM(Name)//'%BlkPt',A%NBlks,UnLimit_O=.TRUE.)
                  CALL Put(A%MTrix,TRIM(Name)//'%MTrix',A%NNon0,UnLimit_O=.TRUE.)
                  RETURN
               ENDIF
            ENDIF

!           WRITE TO BINARY FILE
            INQUIRE(FILE=FileName,EXIST=Exists)
            IF(Exists)THEN
#if FORMATTED
               OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE', &
                    FORM='FORMATTED',ACCESS='SEQUENTIAL')
#else
               OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE', &
                    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
            ELSE
#if FORMATTED
               OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW', &
                    FORM='FORMATTED',ACCESS='SEQUENTIAL')
#else
               OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW', &
                    FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
#endif
            ENDIF
#ifdef FORMATTED
            WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)NAtoms,A%NNon0,A%NBlks
            WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%RowPt%I(i),i=1,NAtoms+1)
            WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%ColPt%I(i),i=1,A%NBlks)
            WRITE(UNIT=Seq,FMT=55,Err=1,IOSTAT=IOS)(A%BlkPt%I(i),i=1,A%NBlks)
            WRITE(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(BIG_DBL,i=1,A%NBlks)
            WRITE(UNIT=Seq,FMT=66,Err=1,IOSTAT=IOS)(A%MTrix%D(i),i=1,A%NNon0)
            INCLUDE 'Formats.Inc'            
#else
            WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)NAtoms,A%NNon0,A%NBlks
            WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%RowPt%I(i),i=1,NAtoms+1)
            WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%ColPt%I(i),i=1,A%NBlks)
            WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%BlkPt%I(i),i=1,A%NBlks)
            WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(BIG_DBL,i=1,A%NBlks)
            WRITE(UNIT=Seq,Err=1,IOSTAT=IOS)(A%MTrix%D(i),i=1,A%NNon0)
#endif
            CLOSE(UNIT=Seq,STATUS='KEEP')
#ifdef PARALLEL
          ENDIF
#endif
          RETURN
        1 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_BCSR.')
      END SUBROUTINE Put_BCSR
#ifdef PARALLEL
!---------------------------------------------------------------------
!     Get a DBCSR matrix
!
      SUBROUTINE Get_DBCSR(A,Name,PFix_O,CheckPoint_O)
         TYPE(DBCSR),              INTENT(INOUT) :: A     
         LOGICAL,         OPTIONAL,INTENT(IN)    :: CheckPoint_O
         CHARACTER(LEN=*),         INTENT(IN)    :: Name
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: PFix_O
         TYPE(BCSR)                              :: B                          
!----------------------------------------------------------------------------
         CALL Get_BCSR(B,Name,PFix_O,CheckPoint_O)
         CALL SetEq(A,B) 
         CALL Delete(B)
         A%Node=MyId
      END SUBROUTINE Get_DBCSR
!---------------------------------------------------------------------
!     Put a DBCSR matrix
!
      SUBROUTINE Put_DBCSR(A,Name,PFix_O,CheckPoint_O)
         TYPE(DBCSR), INTENT(INOUT)           :: A     
         LOGICAL,         OPTIONAL,INTENT(IN) :: CheckPoint_O
         TYPE(BCSR)                           :: B                          
         CHARACTER(LEN=*),         INTENT(IN) :: Name
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: PFix_O
!----------------------------------------------------------------------------
         CALL SetEq(B,A)
         CALL Put_BCSR(B,Name,PFix_O,CheckPoint_O)
         CALL Delete(B)
!--------------------------------------------
      END SUBROUTINE Put_DBCSR
#endif
!------------------------------------------------------------------
!     Put thresholds 
!
      SUBROUTINE Put_TOLS(NGLCT,Tag_O)
         TYPE(TOLS),               INTENT(IN) :: NGLCT
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
         CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Tag
         CALL Put(NGLCT%Cube,'cubeneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%Trix,'trixneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%Dist,'distneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%TwoE,'twoeneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%ETol,'enregyneglect',Tag_O=Tag_O)
         CALL Put(NGLCT%DTol,'densityneglect',Tag_O=Tag_O)
      END SUBROUTINE Put_TOLS

!------------------------------------------------------------------
!     Get thresholds 
!
      SUBROUTINE Get_TOLS(NGLCT,Tag_O)
         IMPLICIT NONE
         TYPE(TOLS),              INTENT(OUT) :: NGLCT
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
         CALL Get(NGLCT%Cube,'cubeneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%Trix,'trixneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%Dist,'distneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%TwoE,'twoeneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%ETol,'enregyneglect',Tag_O=Tag_O)
         CALL Get(NGLCT%DTol,'densityneglect',Tag_O=Tag_O)
      END SUBROUTINE Get_TOLS

!--------------------------------------------------------------------
!     Get arguments from the command line
!
      SUBROUTINE Get_ARGMT(A)
#ifdef NAG
         USE F90_UNIX
         IMPLICIT NONE
#else
         IMPLICIT NONE
         INTEGER,EXTERNAL               :: IARGC
#endif
         TYPE(ARGMT),INTENT(OUT)        :: A
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Tmp1,Tmp2
         INTEGER                        :: I,NArg,NChar,NInts
         INTEGER,PARAMETER              :: MaxArg=10
#ifdef PARALLEL
         IF(MyId==ROOT)THEN
#endif
            NArg=IARGC()
            NChar=0
            NInts=0
            DO I=1,NArg
               CALL GetArg(I,Tmp1)
               CALL LowCase(Tmp1)
               IF(SCAN(Tmp1,Lower)==0)THEN
                  NInts=NInts+1
               ELSE
                  NChar=NChar+1
               ENDIF
            ENDDO
#ifdef PARALLEL
         ENDIF
         IF(InParallel)THEN
            CALL BCast(NInts)
            CALL BCast(NChar)
         ENDIF
#endif
        CALL New(A,(/NChar,NInts/))
#ifdef PARALLEL
        IF(MyId==ROOT)THEN
#endif
            NChar=0
            NInts=0
            DO I=1,NArg
               CALL GetArg(I,Tmp1)
               Tmp2=Tmp1
               CALL LowCase(Tmp2)
               IF(SCAN(Tmp2,Lower)==0)THEN
                  NInts=NInts+1
                  A%I%I(NInts)=CharToInt(TRIM(Tmp1))
               ELSE
                  NChar=NChar+1
                  A%C%C(NChar)=TRIM(Tmp1)
               ENDIF
            ENDDO
#ifdef PARALLEL
         ENDIF
         IF(InParallel)THEN
            CALL BCast(A%I)
            CALL BCast_CHR_VECT(A%C) 
         ENDIF
#endif
      END SUBROUTINE Get_ARGMT
!
!
!
  SUBROUTINE Get_HGRho(A,Name,Args,SCFCycle)
    TYPE(HGRho)                      :: A
    TYPE(ARGMT)                      :: Args
    INTEGER                          :: SCFCycle,IOS,I,NExpt,NDist,NCoef
    CHARACTER(LEN=*)                 :: Name 
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
    LOGICAL                          :: Exists
!   
    FileName=TrixFile(Name,Args,SCFCycle)
    INQUIRE(FILE=FileName,EXIST=Exists)
    IF(Exists) THEN
       OPEN(UNIT=Seq,FILE=FileName,STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ELSE
       CALL Halt(' Get_HGRho could not find '//TRIM(FileName))
    ENDIF
!
!   Allocate Memory
!
    READ(UNIT=Seq,Err=100,IOSTAT=IOS) NExpt,NDist,NCoef
    CALL New_HGRho(A,(/NExpt,NDist,NCoef/))
!
!   Read In the Density
!
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
    READ(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I)    ,I=1,A%NCoef)
!
    Close(UNIT=Seq,STATUS='KEEP')
    RETURN
100 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Get_HGRho.')
  END SUBROUTINE Get_HGRho
!========================================================================================
! Write  the density to disk 
!========================================================================================
  SUBROUTINE Put_HGRho(A,Name,Args,SCFCycle)
    TYPE(HGRho)                      :: A
    TYPE(ARGMT)                      :: Args
    INTEGER                          :: I,SCFCycle,IOS
    REAL(DOUBLE)                     :: Dummy
    CHARACTER(LEN=*)                 :: Name 
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: FileName
    LOGICAL                          :: Exists
!   
    FileName=TrixFile(Name,Args,SCFCycle)
    INQUIRE(FILE=FileName,EXIST=Exists)
    IF(Exists) THEN
       OPEN(UNIT=Seq,FILE=FileName,STATUS='REPLACE',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ELSE
       OPEN(UNIT=Seq,FILE=FileName,STATUS='NEW',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
    ENDIF
!
!   Write density to disk
!
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS) A%NExpt,A%NDist,A%NCoef
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%NQ%I(I)    ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffQ%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%OffR%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Lndx%I(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Expt%D(I)  ,I=1,A%NExpt)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qx%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qy%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Qz%D(I)    ,I=1,A%NDist)
    WRITE(UNIT=Seq,Err=100,IOSTAT=IOS)(A%Co%D(I)    ,I=1,A%NCoef)
!
    CLOSE(UNIT=Seq,STATUS='KEEP')
    RETURN
100 CALL Halt('IO Error '//TRIM(IntToChar(IOS))//' in Put_HGRho.')
  END SUBROUTINE Put_HGRho
!========================================================================================
! Get Cartesian multipoles
!========================================================================================
  SUBROUTINE Get_CMPoles(A,Tag_O)
    TYPE(CMPoles)                    :: A
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O
    IF(.NOT.AllocQ(A%Alloc))CALL New_CMPoles(A)
    CALL Get(A%DPole,'dipole',Tag_O=Tag_O)
    CALL Get(A%QPole,'quadrupole',Tag_O=Tag_O)
  END SUBROUTINE Get_CMPoles
!========================================================================================
! Put Cartesian multipoles
!========================================================================================
  SUBROUTINE Put_CMPoles(A,Tag_O)
    TYPE(CMPoles)                    :: A
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O
    CALL Put(A%DPole,'dipole',Tag_O=Tag_O)
    CALL Put(A%QPole,'quadrupole',Tag_O=Tag_O)
  END SUBROUTINE Put_CMPoles
!------------------------------------------------------------------
!     Open an ASCII file   
! 
      SUBROUTINE OpenASCII(FileName,Unit,NewFile_O,OldFileQ_O,Rewind_O)
         CHARACTER(LEN=*), INTENT(IN) :: FileName
         INTEGER,          INTENT(IN) :: Unit
         LOGICAL, OPTIONAL,INTENT(IN) :: NewFile_O,OldFileQ_O,Rewind_O
         INTEGER                      :: IOS
         LOGICAL                      :: Opened, Exists
!------------------------------------------------------------------
!        Does the file exist?
!
         INQUIRE(FILE=FileName,OPENED=Opened, &
                 EXIST=Exists,ERR=11,IOSTAT=IOS)
         IF(PRESENT(OldFileQ_O))THEN
             IF(OldFileQ_O.AND.(.NOT.Exists)) &
                CALL HALT(' File '//TRIM(FileName)//' does not exist! ')
          ENDIF
!------------------------------------------------------------------
!        Open a new file
!
         IF(PRESENT(NewFile_O))THEN
            IF(NewFile_O.AND.Exists)THEN
!              Open replace if already exists
               OPEN(UNIT=Unit,FILE=FileName, &
                    ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                    ERR=11,IOSTAT=IOS,STATUS='REPLACE')
            ELSE
!              Open a new file
               OPEN(UNIT=Unit,FILE=FileName, &
                    ACCESS='SEQUENTIAL',FORM='FORMATTED', &
                    ERR=11,IOSTAT=IOS,STATUS='NEW')         
            ENDIF
!------------------------------------------------------------------
!        Open existing file and position at the top
!
         ELSEIF(PRESENT(Rewind_O))THEN
            IF(Rewind_O.AND.Exists)THEN
               OPEN(UNIT=Unit,FILE=FileName, &
                    ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                    POSITION='REWIND',ERR=11,IOSTAT=IOS,STATUS='OLD')
            ELSE
               CALL Halt(' Logic error 2 in OpenASCII ')
            ENDIF
!------------------------------------------------------------------
!        Open existing file and position at the bottom (default)
!
         ELSEIF(Exists.AND.(.NOT.Opened))THEN
            OPEN(UNIT=Unit,FILE=FileName, &
                 ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                 POSITION='APPEND',ERR=11,IOSTAT=IOS,STATUS='OLD')
!------------------------------------------------------------------
!        Create a new file and open it
!
         ELSEIF(Exists.AND.Opened)THEN
           CALL Warn(' File '//TRIM(FileName)//' already open')
         ELSE
           OPEN(UNIT=Unit,FILE=FileName, &
                ACCESS='SEQUENTIAL',FORM='FORMATTED', &
                ERR=11,IOSTAT=IOS,STATUS='NEW')         
         ENDIF
         RETURN
      11 CALL Halt(' OpenASCII ERROR: IOS='//TRIM(IntToChar(IOS))// &
                   ' on file '//TRIM(FileName)//'.') 
      END SUBROUTINE OpenASCII

!--------------------------------------------------------------------------
#ifdef MMech
    SUBROUTINE Put_CHR_VECT(A,VarName,NN,Tag_O)
       INTEGER                              :: I,N,II,NN
    CHARACTER(LEN=*),DIMENSION(1:NN),INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,DIMENSION(DEFAULT_CHR_LEN*NN)   :: B  !=ICHAR(' ')
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          N=LEN(A(1))
          IF(N>DEFAULT_CHR_LEN) CALL Halt('Static strings overrun in Put_CHR_VECT')
         DO II = 1, NN
          DO I=1,N; B((II-1)*N+I)=ICHAR(A(II)(I:I)); ENDDO
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
                       N*NN,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,B)
          CALL CloseData(Meta)
         ENDDO
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_CHR_VECT
#endif

!--------------------------------------------------------------------------
#ifdef MMech
    SUBROUTINE Get_CHR_VECT(A,VarName,NN,Tag_O)
       INTEGER                              :: I,N,II,NN
    CHARACTER(LEN=*),DIMENSION(1:NN),INTENT(INOUT):: A
       CHARACTER(LEN=*),         INTENT(IN)    :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: Tag_O
       INTEGER,DIMENSION(DEFAULT_CHR_LEN*NN)   :: B !=ICHAR(' ')
       TYPE(META_DATA)                         :: Meta
#ifdef PARALLEL 
      IF(MyId==ROOT)THEN
#endif 
         NN=SIZE(A,1)
         N=LEN(A)
         IF(N>DEFAULT_CHR_LEN) &
              CALL Halt('Static strings overrun in Get_CHR_VECT')
         Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32, &
              N*NN,.FALSE.)
         CALL OpenData(Meta)
         DO II = 1, NN
            CALL ReadIntegerVector(Meta,B)
            DO I=1,N; A(II)(I:I)=CHAR(B((II-1)*N+I)); ENDDO
            ENDDO
            CALL CloseData(Meta)
#ifdef PARALLEL 
         ENDIF       
         IF(InParallel)CALL BCast(A)
#endif 
       END SUBROUTINE Get_CHR_VECT
#endif 
!--------------------------------------------------------------------------
#ifdef MMech
    SUBROUTINE Put_LOG_VECT(A,VarName,NN,Tag_O)
       INTEGER                              :: I,N,II,NN
       LOGICAL,DIMENSION(1:NN),  INTENT(IN) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,DIMENSION(1:NN)              :: ILog
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
       DO I = 1, NN
         ILog(I)=0
         IF(A(I))ILog(I)=1
       ENDDO
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL WriteIntegerVector(Meta,ILog)
          CALL CloseData(Meta)
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Put_LOG_VECT
#endif
!--------------------------------------------------------------------------
#ifdef MMech
    SUBROUTINE Get_LOG_VECT(A,VarName,NN,Tag_O)
       INTEGER                              :: I,NN
       LOGICAL,DIMENSION(1:NN),INTENT(INOUT) :: A
       CHARACTER(LEN=*),         INTENT(IN) :: VarName
       CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Tag_O
       INTEGER,DIMENSION(1:NN)              :: ILog
       TYPE(META_DATA)                      :: Meta
#ifdef PARALLEL 
       IF(MyId==ROOT)THEN
#endif 
          Meta=SetMeta(NameTag(VarName,Tag_O),NATIVE_INT32,1,.FALSE.)
          CALL OpenData(Meta,.TRUE.)
          CALL ReadIntegerVector(Meta,ILog)
          CALL CloseData(Meta)
       DO I = 1, NN
         IF(ILog(I)==1) Then
           A(I) = .TRUE.
         ELSE
           A(I) = .FALSE.
         ENDIF
       ENDDO
#ifdef PARALLEL 
       ENDIF       
#endif 
    END SUBROUTINE Get_LOG_VECT
#endif
!------------------------------------------------------------------------


END MODULE

