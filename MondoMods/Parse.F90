!    PARSING AND CHARACTER MANIPULATION ROUTINES FOR MONDOSCF
!    Author: Matt Challacombe and CJ Tymczak
!-----------------------------------------------------------------------
MODULE Parse
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalCharacters
   USE ParsingConstants
   USE ProcessControl       
   USE MemMan
   IMPLICIT NONE   
   CONTAINS
!======================================================================
!     Align file pointer one line past a character key
!======================================================================
      SUBROUTINE Align(Key,Unit)     
         CHARACTER(LEN=*)               :: Key
         INTEGER                        :: Unit
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
         REWIND(UNIT=Unit)
         DO 
            READ(Unit,DEFAULT_CHR_FMT,END=99)Line
            IF(INDEX(Line,Key)/=0)RETURN
         ENDDO
 99      CALL Halt(' Key not found in input : '//TRIM(Key))
      END SUBROUTINE Align
!======================================================================
!     Align file pointer one a lowercase-line past a character key
!======================================================================
      SUBROUTINE AlignLowCase(Key,Unit)     
         CHARACTER(LEN=*)               :: Key
         INTEGER                        :: Unit
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
         REWIND(UNIT=Unit)
         DO 
            READ(Unit,DEFAULT_CHR_FMT,END=99)Line
            Call LowCase(Line)
            IF(INDEX(Line,Key)/=0)RETURN
         ENDDO
 99      CALL Halt(' Key not found in input : '//TRIM(Key))
      END SUBROUTINE AlignLowCase
!======================================================================
!     Align file pointer one line past a character key
!======================================================================
      FUNCTION FindKey(Key,Unit)     
        CHARACTER(LEN=*)               :: Key
        INTEGER                        :: Unit
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
        LOGICAL                        :: FindKey
        FindKey = .TRUE.
        REWIND(UNIT=Unit)
        DO 
           READ(Unit,DEFAULT_CHR_FMT,END=99)Line
           IF(INDEX(Line,Key)/=0)RETURN
        ENDDO
99      CONTINUE
        FindKey = .FALSE.
        RETURN
      END FUNCTION FindKey
!==============================================================
      FUNCTION FindMixedCaseKey(Key,Unit)     
        CHARACTER(LEN=*)               :: Key
        INTEGER                        :: Unit,N
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,LineLC,KeyLC
        LOGICAL                        :: FindMixedCaseKey
        FindMixedCaseKey = .TRUE.
!
        REWIND(UNIT=Unit)
        N=LEN(Key)
        KeyLC(1:N)=Key
        CALL LowCase(KeyLC)
        DO 
           READ(Unit,DEFAULT_CHR_FMT,END=99)Line
           LineLC=Line
           CALL LowCase(LineLC)
           IF(INDEX(LineLC,KeyLC(1:N))/=0)RETURN
        ENDDO
99      CONTINUE
        FindMixedCaseKey = .FALSE.
        RETURN
      END FUNCTION FindMixedCaseKey
!======================================================================
!     Convert a string to all lower case
!======================================================================
      SUBROUTINE LowCase(String)
         INTEGER           :: I,J
         CHARACTER(LEN=*)  :: String
         DO I=1,LEN(String)
            J=INDEX(Upper,String(I:I))
            IF(J.NE.0)String(I:I)=Lower(J:J)
         ENDDO
      END SUBROUTINE LowCase
!======================================================================
!     Determine if an option has a defined key set.
!======================================================================
      FUNCTION OptCharQ(Unit,Option,Char)
         IMPLICIT NONE
         INTEGER,         INTENT(IN)    :: Unit
         CHARACTER(LEN=*),INTENT(IN)    :: Option
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,LCLine,Char,OptLC
         INTEGER                        :: J,K1,K2,L,N, &
                                           S1,S2,LD,Pls1,Mns1
         LOGICAL                        :: OptCharQ
!------------------------------------------------------------------
         REWIND(Unit)
         OptLC=Option
!        CALL LowCase(OptLC)
         LD=LEN(Delimiters)-1 ! to avoid the blank!
         DO 
           READ(Unit,DEFAULT_CHR_FMT,END=1)Line
           LCLine=Line
!          Make options case sensitive 
!          CALL LowCase(LCLine)
           IF(INDEX(LCLine,TRIM(OptLC))/=0) THEN
!              WRITE(*,'(A,A)') 'LCLine is ',TRIM(LCLine)
!              WRITE(*,'(A,A)') 'OptLC is ',TRIM(OptLC)
              CALL LowCase(LCLine)
              J=1
              L=LEN(TRIM(LCLine))   
              DO N=1,L
                 S1=SCAN(LCLine(J:L),Delimiters(1:LD))
                 K1=J-1+S1
                 S2=SCAN(LCLine(K1+1:L),Delimiters(1:LD-1))
                 IF(S1==0)RETURN
                 K2=K1+S2
                 IF(K1==K2)K2=L
                 IF(K2==L.AND.SCAN(LCLine(K2:K2),Delimiters(1:LD))==0)THEN
                    Mns1=0
                 ELSE
                    Mns1=-1
                 ENDIF 
                 Pls1=1
!                 WRITE(*,111 )K1,S1,K2,S2,J,LCLine(K1:K2),LCLine(K1+Pls1:K2+Mns1)
                 IF(S1==0.OR.K1>L.OR.K2>L)RETURN
                 IF(K1==K2)K2=L
                 J=K2
                 IF(SCAN(LCLine(K1+Pls1:K2+Mns1),Lower)/=0)THEN
                    Char=ADJUSTL(Line(K1+Pls1:K2+Mns1))
                    OptCharQ=.TRUE.
                    RETURN
                 ENDIF
              ENDDO                
!111           FORMAT(' K1 = ',I3,' S1 = ',I3,' K2 = ',I3, &
!                     ' S2 = ',I3,'J = ',I3,' ALine1 = <',A,'>, ',' ALine2 = <',A,'>')
           ENDIF
         ENDDO
       1 OptCharQ=.FALSE. 
      END FUNCTION OptCharQ
!======================================================================
!     Determine if an option has a defined key set.
!======================================================================
!-----------------------------------------------------------------------
!     Determine if an option has a defined key set.
!
      FUNCTION OptKeyQ(Unit,Option,Key)
         IMPLICIT NONE
         INTEGER,         INTENT(IN)    :: Unit
         CHARACTER(LEN=*),INTENT(IN)    :: Option,Key
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,OptLC
         LOGICAL                        :: OptKeyQ
         REWIND(Unit)


         OptLC=Option
         CALL LowCase(OptLC)
         DO 
           READ(Unit,DEFAULT_CHR_FMT,END=1)Line
           CALL LowCase(Line)
           IF(INDEX(Line,TRIM(OptLC))/=0)THEN
              IF(KeyQ(Line,Key))THEN
                 OptKeyQ=.TRUE.
                 RETURN
              ENDIF
           ENDIF
         ENDDO
       1 OptKeyQ=.FALSE. 
      END FUNCTION OptKeyQ

!-----------------------------------------------------------------------
!     Determine if an option has a defined integer
!
      FUNCTION OptIntQ(Unit,Option,Int)
         INTEGER,         INTENT(IN)    :: Unit
         CHARACTER(LEN=*),INTENT(IN)    :: Option
         INTEGER,         INTENT(Out)   :: Int
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,OptLC
         LOGICAL                        :: OptIntQ
         INTEGER                        :: J,L,K1,K2,S1,S2
         REWIND(Unit)
         OptLC=Option
         CALL LowCase(OptLC)
         DO 
            READ(Unit,DEFAULT_CHR_FMT,END=1)Line
            CALL LowCase(Line)
            IF(INDEX(Line,TRIM(OptLC))/=0)THEN
!               WRITE(*,*)' OptLC = <',TRIM(OptLC),'>'
!               WRITE(*,*)' Line = <',TRIM(Line),'>'
               J=1
               L=LEN(Line)   
               S1=SCAN(Line(J:L),Numbers)
               IF(S1==0)GOTO 1
               K1=J-1+S1
               S2=SCAN(Line(K1:L),Delimiters)            
               K2=K1-2+S2
!               WRITE(*,*)' S1 = ',S1,' S2 = ',S2,' Line = <',TRIM(Line(K1:K2)),'>'
               Int=CharToInt(TRIM(Line(K1:K2)))
               OptIntQ=.TRUE.
               RETURN
            ENDIF
         ENDDO
       1 OptIntQ=.FALSE. 
      END FUNCTION OptIntQ

!-----------------------------------------------------------------------
!     Determine if an option has a defined double
!
      FUNCTION OptDblQ(Unit,Option,Dbl)
         INTEGER,         INTENT(IN)    :: Unit
         CHARACTER(LEN=*),INTENT(IN)    :: Option
         REAL(DOUBLE),    INTENT(InOut) :: Dbl
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,OptLC
         LOGICAL                        :: OptDblQ
         INTEGER                        :: J,L,K1,K2,S1,S2
         REWIND(Unit)
         OptLC=Option
         CALL LowCase(OptLC)
         DO 
            READ(Unit,DEFAULT_CHR_FMT,END=1)Line
            CALL LowCase(Line)
            IF(INDEX(Line,TRIM(OptLC))/=0)THEN
!               WRITE(*,*)' OptLC = <',TRIM(OptLC),'>'
!               WRITE(*,*)' Line = <',TRIM(Line),'>'
               J=1
               L=LEN(Line)   
               S1=SCAN(Line(J:L),Numbers)
               IF(S1==0)GOTO 1
               K1=J-1+S1
               S2=SCAN(Line(K1:L),Delimiters)            
               K2=K1-2+S2
!               WRITE(*,*)' S1 = ',S1,' S2 = ',S2,' Line = <',TRIM(Line(K1:K2)),'>'
               Dbl=CharToDbl(TRIM(Line(K1:K2)))
!               WRITE(*,*)' Dbl = ',Dbl
               OptDblQ=.TRUE.
               RETURN
            ENDIF
         ENDDO
       1 OptDblQ=.FALSE. 
      END FUNCTION OptDblQ
!-----------------------------------------------------------------------
!     Determine if a key marked _ON BOTH SIDES_ 
!     by Delimiters exists in a line
!
      FUNCTION KeyQ(Line,Key)
         CHARACTER(LEN=*),INTENT(IN)    :: Line,Key
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: TmpLine,TmpKey
         INTEGER                        :: J,L,M,N, &
                                           K1,S1,K2,S2,Mns1,Pls1
         LOGICAL                        :: KeyQ
!         WRITE(*,*)' Line = ',TRIM(Line)
!         WRITE(*,*)' Key  = ',TRIM(Key)
         KeyQ=.FALSE.               
         TmpLine=Line
         TmpKey=Key
         CALL LowCase(TmpLine)
         CALL LowCase(TmpKey)
         J=1
         L=LEN(TRIM(TmpLine))   
         M=LEN(TRIM(TmpKey))
!         WRITE(*,*)' TmpKey = <',TmpKey(1:M),'>'
         DO N=1,L
            S1=SCAN(TmpLine(J:L),Delimiters)
            K1=J-1+S1
            S2=SCAN(TmpLine(K1+1:L),Delimiters)
            IF(S1==0)RETURN
            K2=K1+S2
            IF(K1==K2)K2=L
            IF(K2==L.AND.SCAN(TmpLine(K2:K2),Delimiters)==0)THEN
               Mns1=0
            ELSE
               Mns1=-1
            ENDIF 
!            IF(K1==1)THEN
!               Pls1=0
!            ELSE
               Pls1=1
!            ENDIF 
!            WRITE(*,111 )K1,S1,K2,S2,J,TmpLine(K1:K2),TmpLine(K1+Pls1:K2+Mns1)
111         FORMAT(' K1 = ',I3,' S1 = ',I3,' K2 = ',I3, &
                   ' S2 = ',I3,'J = ',I3,' ALine1 = <',A,'>, ',' ALine2 = <',A,'>')
            IF(S1==0.OR.K1>L.OR.K2>L)RETURN
            IF(K1==K2)K2=L
            J=K2
            IF(TmpLine(K1+Pls1:K2+Mns1)==TmpKey(1:M))THEN
               KeyQ=.TRUE.
               RETURN
            ENDIF
         ENDDO                
      END FUNCTION KeyQ
!-----------------------------------------------------------------------
!    
!-----------------------------------------------------------------------
      FUNCTION OptKeyLocQ(Unit,Option,Key,M,N,Loc)
         INTEGER,         INTENT(IN)    :: Unit
         CHARACTER(LEN=*),INTENT(IN)    :: Option,Key
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,LCLine,Char,OptLC,KeyLC
         INTEGER                        :: M,N,NN,J,K1,K2,L,V,S1,S2,LD,Pls1,Mns1
         INTEGER, DIMENSION(M)          :: Loc
         LOGICAL                        :: OptKeyLocQ
         REWIND(Unit)
         OptLC=Option
         KeyLC=Key
         CALL LowCase(OptLC)
         CALL LowCase(KeyLC)
!         WRITE(*,*)' Opt, key = ',OptLC,KeyLC
         N=0
         NN=0
         Loc=0
         OptKeyLocQ=.FALSE. 
         LD=LEN(Delimiters)-1 ! to avoid the blank!
         DO 
           READ(Unit,DEFAULT_CHR_FMT,END=1)Line
           LCLine=Line
           CALL LowCase(LCLine)
           IF(INDEX(LCLine,TRIM(OptLC))/=0)THEN
!              WRITE(*,*)TRIM(LCLine)
              J=1
              L=LEN(TRIM(LCLine))   
              DO V=1,L
                 S1=SCAN(LCLine(J:L),Delimiters(1:LD))
                 K1=J-1+S1
                 S2=SCAN(LCLine(K1+1:L),Delimiters(1:LD))
                 IF(S1==0)RETURN
                 K2=K1+S2
                 IF(K1==K2)K2=L
                 IF(K2==L.AND.SCAN(LCLine(K2:K2),Delimiters(1:LD))==0)THEN
                    Mns1=0
                 ELSE
                    Mns1=-1
                 ENDIF 
                 Pls1=1
!                 WRITE(*,111 )N,NN,LCLine(K1:K2),LCLine(K1+Pls1:K2+Mns1)
                 IF(S1==0.OR.K1>L.OR.K2>L)RETURN
                 IF(K1==K2)K2=L
                 J=K2
                 IF(K1+Pls1<=K2+Mns1)THEN
                    NN=NN+1
!                    WRITE(*,11)TRIM(ADJUSTL(KeyLC)),TRIM(ADJUSTL(LCLine(K1+Pls1:K2+Mns1)))
!11 format(' Test1 = <',A,'>, ',' Test2 = <',A,'>')
                 IF(TRIM(ADJUSTL(KeyLC))==TRIM(ADJUSTL(LCLine(K1+Pls1:K2+Mns1))))THEN
!                 IF(SCAN(LCLine(K1+Pls1:K2+Mns1),Lower)/=0)THEN
!                    Char=ADJUSTL(Line(K1+Pls1:K2+Mns1))                    
                    N=N+1
                    Loc(N)=NN
                 ENDIF
                 ENDIF
              ENDDO      
111           FORMAT('NLoc = ',I1,' Loc = ',I3,' ALine1 = <',A,'>, ',' ALine2 = <',A,'>')
           ENDIF
         ENDDO
       1 CONTINUE
!         WRITE(*,*)' 1 CONTINUE 1 CONTINUE 1 CONTINUE 1 CONTINUE 1 CONTINUE '
         IF(N==0)THEN
            OptKeyLocQ=.FALSE. 
         ELSE
            OptKeyLocQ=.TRUE.
         ENDIF
      END FUNCTION OptKeyLocQ

!------------------------------------------------------------------
!     Resolve a line into geometry components (AtomType,X,Y,Z) 
!
      SUBROUTINE LineToGeom(Line,At,Carts)
         CHARACTER(LEN=*)                    :: Line
         CHARACTER(LEN=DEFAULT_CHR_LEN)      :: TmpLine
         CHARACTER(LEN=2)                    :: At      
         INTEGER                             :: J,L,N,K1,K2
         REAL(DOUBLE), DIMENSION(:)          :: Carts
!
         Carts=Zero
         TmpLine=Line
!
         CALL LowCase(TmpLine) 
         J=SCAN(TmpLine,Lower)
         IF(J==0)Call Halt(' No characters found in line in LineToGeom ')
         At=TmpLine(J:J+1)
         J=J+2
         L=LEN(Line)   
         DO N=1,SIZE(Carts)
            K1=J-1+SCAN(Line(J:L),Numbers)
            K2=K1-2+SCAN(Line(K1:L),' ') 
            J=K2+1
            Carts(N)=CharToDbl(TRIM(Line(K1:K2)))
         ENDDO 
!   
      END SUBROUTINE LineToGeom
!------------------------------------------------------------------
!
      SUBROUTINE LineToDbls(Line,N,Dbls)
         INTEGER,                   INTENT(IN)  :: N
         CHARACTER(LEN=*),          INTENT(IN)  :: Line
         REAL(DOUBLE), DIMENSION(N),INTENT(OUT) :: Dbls
         CHARACTER(LEN=DEFAULT_CHR_LEN)         :: TmpLine
         INTEGER                                :: I,J,K,L,K1,K2,S1A,S1B,S1,S2
         TYPE(CHR_VECT)                         :: Chars
         CALL LineToChars(Line,Chars)
         L=0
         DO I=1,SIZE(Chars%C)
            CALL LowCase(Chars%C(I))
            K=SCAN(Chars%C(I)(1:4),Lower)
            IF(K==0)THEN
               L=L+1
               Dbls(L)=CharToDbl(TRIM(Chars%C(I)))
            ENDIF
            IF(L==N)EXIT
         ENDDO  
         IF(N/=L)CALL Halt(' Parse error in LineToDbls, N/=L')
         CALL Delete(Chars)
         RETURN
!        OLD STYLE ...
         TmpLine=Line
         CALL LowCase(TmpLine)
         J=1
         L=LEN(TmpLine)   
         I=0
         DO 
            S1B=SCAN(TmpLine(J:L),Numbers)
            S1A=SCAN(TmpLine(J:S1B),Lower)
            S1=S1B
            K1=J-1+S1
            S2=SCAN(TmpLine(K1:L),Delimiters)            
            K2=K1-2+S2
            J=K2+1
!            WRITE(*,*)' S1A = ',S1A,' S1B = ',S1B,' J = ',J,' L = ',L
            IF(S1A==0)THEN
               I=I+1
!               WRITE(*,*)' TMPLINE = ',TRIM(TmpLine(K1:K2))
               Dbls(I)=CharToDbl(TRIM(TmpLine(K1:K2)))
!               WRITE(*,*)' I = ',I,' DBLS = ',Dbls(I)
               IF(I==N)EXIT
            ENDIF
         ENDDO    
      END SUBROUTINE LineToDbls


      SUBROUTINE LineToInts(Line,N,Ints)
         INTEGER,INTENT(IN)                :: N
         CHARACTER(LEN=*),     INTENT(IN)  :: Line
         INTEGER, DIMENSION(N),INTENT(OUT) :: Ints
         INTEGER                           :: I,J,L,K1,K2,S1,S2
         J=1
         L=LEN(Line)   
         DO I=1,N
            S1=SCAN(Line(J:L),Numbers)
            K1=J-1+S1
            S2=SCAN(Line(K1:L),Delimiters)            
            K2=K1-2+S2
            Ints(I)=CharToInt(TRIM(Line(K1:K2)))
            J=K2+1
         ENDDO    
      END SUBROUTINE LineToInts
!------------------------------------------------------------------
!     Convert a character string into a vector of character strings
!
      SUBROUTINE LineToChars(Line,Chars,NULL_O)
         CHARACTER(LEN=*),  INTENT(IN)    :: Line
         LOGICAL, OPTIONAL, INTENT(IN)    :: NULL_O
         CHARACTER(LEN=5*DEFAULT_CHR_LEN) :: TmpLine
         TYPE(CHR_VECT),    INTENT(OUT)   :: Chars
         INTEGER                          :: I,J,L,N,K1,K2,S1
         TmpLine=Blanks
         L=LEN(TRIM(Line))   
         TmpLine(1:L)=Line(1:L)         
         CALL LowCase(TmpLine)
         J=1
         N=0
         DO I=1,L
            S1=SCAN(TmpLine(J:L),Characters)
            IF(S1==0)EXIT
            K1=J-1+S1
            J=K1+1
            DO 
              IF(SCAN(TmpLine(J:J),Characters)==0)EXIT
              J=J+1
            ENDDO
            K2=J-1
            N=N+1
            J=K2+1
         ENDDO   
         IF(PRESENT(NULL_O))N=N+1
         IF(Chars%Alloc==ALLOCATED_TRUE) &
         CALL Delete(Chars)
         CALL New(Chars,N)
         J=1
         N=0
         DO I=1,L
            S1=SCAN(TmpLine(J:L),Characters)
            IF(S1==0)EXIT
            K1=J-1+S1
            J=K1+1
            DO 
              IF(SCAN(TmpLine(J:J),Characters)==0)EXIT
              J=J+1
            ENDDO
            K2=J-1
            Chars%C(N+1)=ADJUSTL(Line(K1:K2))
            N=N+1
            J=K2+1
         ENDDO   
         IF(PRESENT(NULL_O))Chars%C(N+1)=Blnk
111      FORMAT(' K1 = ',I3,' S1 = ',I3,' K2 = ',I3, &
                ' S2 = ',I3,'J = ',I3,' ALine1 = <',A,'>')
112      FORMAT(' K1 = ',I3,' S1 = ',I3,' K2 = ',I3, &
                ' S2 = ',I3,'J = ',I3,' Chars = <',A,'>')
      END SUBROUTINE LineToChars
!------------------------------------------------------------------
!     Convert a character string into an integer
!
      FUNCTION CharToInt(C)
         CHARACTER(LEN=*),INTENT(IN) :: C
         INTEGER                     :: CharToInt
         READ(UNIT=C,FMT=INTERNAL_INT_FMT)CharToInt
     END FUNCTION CharToInt
!------------------------------------------------------------------
!     Convert a character string into a double or DBLE(integer)
!
      FUNCTION CharToDbl(C)
         CHARACTER(LEN=*),INTENT(IN) :: C
         REAL(DOUBLE)                :: CharToDbl
         INTEGER                     :: TempInt
         IF(SCAN(C,'.')==0)THEN
           ! Maybe its an integer, lets try that ...
           READ(UNIT=C,FMT=INTERNAL_INT_FMT)TempInt
           ! Now cast to double
           CharToDbl=DBLE(TempInt)
         ELSE
           ! Some sort of decimal, go for the DBL_FMT
            READ(UNIT=C,FMT=INTERNAL_DBL_FMT)CharToDbl
         ENDIF
     END FUNCTION CharToDbl
!------------------------------------------------------------------
!     Convert an integer into a character string
!
      FUNCTION IntToChar(I)
         INTEGER  :: I
         CHARACTER(LEN=INTERNAL_INT_LEN) :: IntToChar
         WRITE(UNIT=IntToChar,FMT=INTERNAL_INT_FMT)I
         IntToChar=ADJUSTL(IntToChar)
     END FUNCTION IntToChar
!------------------------------------------------------------------
!     Convert a double into a character string
!
      FUNCTION DblToChar(D)
         REAL(DOUBLE),INTENT(IN)         :: D
         CHARACTER(LEN=INTERNAL_DBL_LEN) :: DblToChar
         WRITE(UNIT=DblToChar,FMT=INTERNAL_DBL_FMT)D
         DblToChar=ADJUSTL(DblToChar)
     END FUNCTION DblToChar
!------------------------------------------------------------------
!     Convert a float into a character string
!
      FUNCTION FltToChar(D)
         REAL(DOUBLE),INTENT(IN)         :: D
         CHARACTER(LEN=INTERNAL_FLT_LEN) :: FltToChar
         WRITE(UNIT=FltToChar,FMT=INTERNAL_FLT_FMT)D
         FltToChar=ADJUSTL(FltToChar)
     END FUNCTION FltToChar
!------------------------------------------------------------------
!     Convert a float into a character string
!
      FUNCTION FltToMedmChar(D)
         REAL(DOUBLE),INTENT(IN)         :: D
         CHARACTER(LEN=INTERNAL_FLT_LEN) :: FltToMedmChar
         WRITE(UNIT=FltToMedmChar,FMT='(F18.8)')D
         FltToMedmChar=ADJUSTL(FltToMedmChar)
     END FUNCTION FltToMedmChar
!------------------------------------------------------------------
!     Convert a double into a medium character string
!
      FUNCTION DblToMedmChar(D)
         REAL(DOUBLE),INTENT(IN)         :: D
         CHARACTER(LEN=14) :: DblToMedmChar
         WRITE(UNIT=DblToMedmChar,FMT='(D14.8)')D
         DblToMedmChar=ADJUSTL(DblToMedmChar)
     END FUNCTION DblToMedmChar
!------------------------------------------------------------------
!     Convert a double into a short character string
!
      FUNCTION DblToShrtChar(D)
         REAL(DOUBLE),INTENT(IN)         :: D
         CHARACTER(LEN=8) :: DblToShrtChar
         WRITE(UNIT=DblToShrtChar,FMT='(D8.2)')D
         DblToShrtChar=ADJUSTL(DblToShrtChar)
     END FUNCTION DblToShrtChar

     FUNCTION DblToMMAChar(D)
         REAL(DOUBLE),INTENT(IN)         :: D
         CHARACTER(LEN=30) :: DblToMMAChar
         DblToMMAChar=TRIM(FltToChar(FRACTION(D))) &
            //'*2^('//TRIM(IntToChar(EXPONENT(D)))//')'
     END FUNCTION DblToMMAChar
! 
      FUNCTION TrixFile(PostFix,Args_O,OffSet_O,Name_O,Stats_O,NoTags_O,PWD_O)
         CHARACTER(LEN=*),         INTENT(IN) :: PostFix 
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Name_O
         TYPE(ARGMT),     OPTIONAL,INTENT(IN) :: Args_O
         INTEGER,DIMENSION(3), &
                          OPTIONAL,INTENT(IN) :: Stats_O
         INTEGER,         OPTIONAL,INTENT(IN) :: OffSet_O        
         LOGICAL,         OPTIONAL,INTENT(IN) :: NoTags_O,PWD_O
         INTEGER                              :: OffSet
         INTEGER, DIMENSION(3)                :: Stats
         CHARACTER(LEN=DEFAULT_CHR_LEN)       :: Name,TrixFile,Cycl,Base,Geom
         IF(PRESENT(Name_O))THEN
            IF(PRESENT(PWD_O))THEN
               Name=TRIM(MONDO_PWD)//Name_O
            ELSE
               Name=TRIM(MONDO_SCRATCH)//Name_O
            ENDIF
         ELSEIF(PRESENT(Args_O))THEN
            IF(PRESENT(PWD_O))THEN
               Name=PWDName
            ELSE
               Name=SCRName
            ENDIF               
         ELSE
            CALL Halt(' Neither Name_O or Args_O passed to TrixFile ! ')
         ENDIF
         IF(PRESENT(Stats_O))THEN
            Stats(1:3)=Stats_O(1:3)
         ELSEIF(PRESENT(Args_O))THEN
            Stats(1:3)=Args_O%I%I(1:3)
         ELSE
            CALL Halt(' Neither Stats_O or Args_O passed to TrixFile ! ')
         ENDIF
         IF(PRESENT(OffSet_O))THEN
            Stats(1)=Stats(1)+OffSet_O
            IF(Stats(1)<0)THEN
!              Use previous if current < 0
               Stats=Args_O%I%I(4:6)
            ENDIF
            Cycl=IntToChar(Stats(1))
            Base=IntToChar(Stats(2))
            Geom=IntToChar(Stats(3))
#ifdef PARALLEL_CLONES
            TrixFile=TRIM(Name)//'_Geom#'//TRIM(Geom) &
                               //'_Base#'//TRIM(Base) &
                               //'_Cycl#'//TRIM(Cycl) &
                               //'_Clone#'//TRIM(IntToChar(MyClone)) &
                               //'.'//TRIM(PostFix)            
#else

            TrixFile=TRIM(Name)//'_Geom#'//TRIM(Geom) &
                               //'_Base#'//TRIM(Base) &
                               //'_Cycl#'//TRIM(Cycl) &
                               //'.'//TRIM(PostFix)            
#endif
         ELSEIF(PRESENT(NoTags_O))THEN
#ifdef PARALLEL_CLONES
            TrixFile=TRIM(Name)//'_Clone#'//TRIM(IntToChar(MyClone))//'.'//TRIM(PostFix)            
#else
            TrixFile=TRIM(Name)//'.'//TRIM(PostFix)            
#endif

         ELSE
            Base=IntToChar(Stats(2))
            Geom=IntToChar(Stats(3))
#ifdef PARALLEL_CLONES
            TrixFile=TRIM(Name)//'_Geom#'//TRIM(Geom) &
                               //'_Base#'//TRIM(Base) &
                               //'_Clone#'//TRIM(IntToChar(MyClone)) &
                               //'.'//TRIM(PostFix)            
#else
            TrixFile=TRIM(Name)//'_Geom#'//TRIM(Geom) &
                               //'_Base#'//TRIM(Base) &
                               //'.'//TRIM(PostFix)            
#endif
         ENDIF
      END FUNCTION TrixFile
!------------------------------------------------------------------
      FUNCTION StatsToChar(Stats) RESULT(StatString)
         INTEGER,DIMENSION(3) :: Stats
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: StatString
         StatString='_'//TRIM(IntToChar(Stats(3)))  &
                   //'_'//TRIM(IntToChar(Stats(2))) &
                   //'_'//TRIM(IntToChar(Stats(1)))
      END FUNCTION StatsToChar
!------------------------------------------------------------------
!     Remove all blanks from a string
!
      FUNCTION Squish(S) 
         CHARACTER(LEN=*), INTENT(IN) :: S
         CHARACTER(LEN=LEN(S))        :: Squish
         INTEGER                      :: I,J,LenS
         J=0
         LenS=LEN(S)
         Squish=Blnk
         DO I=1,LenS
            IF(S(I:I)/=Blnk)THEN
               J=J+1
               Squish(J:J)=S(I:I)
            ENDIF
         ENDDO
     END FUNCTION Squish
  !
  !
  LOGICAL FUNCTION ChrCkkIfInt(Chr)
!H---------------------------------------------------------------------------------
!H LOGICAL FUNCTION ChrCkkIfInt(Chr)
!H  Check if <Chr> is formed only with integers.
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN) :: Chr
    !-------------------------------------------------------------------
    INTEGER                      :: I
    !-------------------------------------------------------------------
    ChrCkkIfInt=.FALSE.
    DO I=1,LEN(TRIM(Chr))
       IF(IACHAR(Chr(I:I)).LT.IACHAR('0').OR.IACHAR(Chr(I:I)).GT.IACHAR('9')) RETURN
    ENDDO
    ChrCkkIfInt=.TRUE.
  END FUNCTION ChrCkkIfInt
  !
  !
  LOGICAL FUNCTION OptGetKeyArg(Unit,Key,Arg)
!H---------------------------------------------------------------------------------
!H LOGICAL FUNCTION OptGetKeyArg(Unit,Key,Arg)
!H  This function look for a given Key in the unit Unit and gives back an array
!H  of char with the arguments.
!H                 e.g. Key=(aaaa,bbbb,cccc,1111,0101sss,aa_bbb) or
!H                      Key=(1a1a1a)                             or
!H                      Key=abcd_efg
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    INTEGER         , INTENT(IN   ) :: Unit
    CHARACTER(LEN=*), INTENT(IN   ) :: Key
    TYPE(CHR_VECT)  , INTENT(INOUT) :: Arg
    !-------------------------------------------------------------------
    INTEGER                         :: NDim,ISTAT,Indx,IndxL,IndxR,iDim
    CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line,LineTmp
    !-------------------------------------------------------------------
    !
    ! Initialize.
    OptGetKeyArg=.FALSE.
    IF(AllocQ(Arg%Alloc)) CALL Delete(Arg)
    REWIND(Unit)
    !
    ! Let,s go.
    DO
       !
       ! Read a line from inputfile.
       READ(Unit,DEFAULT_CHR_FMT,IOSTAT=ISTAT) Line
       !write(*,*) 'ISTAT',ISTAT
       IF(ISTAT.NE.0) EXIT
       !write(*,*) 'Line<'//trim(line)//'>'
       CALL LowCase(Line)
       IF(INDEX(Line,Key).EQ.0) CYCLE
       !
       ! We have found the Key.
       OptGetKeyArg=.TRUE.
       Line=TRIM(Line)
       !
       ! Find the number of comma.
       LineTmp=Line
       NDim=0
       DO
          Indx=INDEX(LineTmp,',')
          lineTmp=TRIM(lineTmp(Indx+1:))
          IF(Indx.EQ.0) EXIT
          NDim=NDim+1
       ENDDO
       NDim=NDim+1
       !
       ! Allocate the Arg array.
       CALL New(Arg,NDim)
       !
       ! Find the equal.
       Indx=INDEX(Line,'=')
       IF(Indx.EQ.0) CALL Halt('Could not find an = then looking for the Key=<'//TRIM(Key)//'>')
       Line=Line(Indx+1:)
       !
       ! Find parenthesis.
       IndxL=SCAN(Line,'('     )
       IndxR=SCAN(Line,')',BACK=.TRUE.)
       !
       ! Check for unbalanced parenthesis.
       IF(IndxL.EQ.0.AND.NDim .GT.1) &
            & CALL Halt('Unbalenced left parenthesis then looking for the Key=<' //TRIM(Key)//'>')
       IF(IndxR.EQ.0.AND.NDim .GT.1) &
            & CALL Halt('Unbalenced right parenthesis then looking for the Key=<'//TRIM(Key)//'>')
       IF(IndxR.NE.0.AND.IndxL.EQ.0) &
            & CALL Halt('Unbalenced left parenthesis then looking for the Key=<' //TRIM(Key)//'>')
       IF(IndxL.NE.0.AND.IndxR.EQ.0) &
            & CALL Halt('Unbalenced right parenthesis then looking for the Key=<'//TRIM(Key)//'>')
       !
       ! Remove the parenthesis.
       IF(IndxR.EQ.0) THEN
          Line=Line(IndxL+1:)
       ELSE
          Line=Line(IndxL+1:IndxR-1)
       ENDIF
       !
       ! Get the arguments.
       DO iDim=1,NDim
          IF(iDim.LT.NDim) THEN
             IndxR=SCAN(Line,',')
             Arg%C(iDim)=TRIM(Line(1:IndxR-1))
             Line=Line(IndxR+1:)
          ELSE
             Arg%C(iDim)=TRIM(Line(1:))
          ENDIF
          !
          ! Check for missing arguments.
          IF(LEN(TRIM(Arg%C(iDim))).EQ.0) &
               & CALL Halt('Miss an argument at the position '//IntToChar(iDim) &
               &           //' then looking for the Key=<' //TRIM(Key)//'>')
       ENDDO
       !
       ! Look for duplicate Key.
       READ(Unit,DEFAULT_CHR_FMT,IOSTAT=ISTAT) Line
       IF(ISTAT.NE.0) EXIT
       CALL LowCase(Line)
       IF(INDEX(Line,Key).EQ.0) THEN
          CYCLE
       ELSE
          !
          ! There is another Key!
          OptGetKeyArg=.FALSE.
          CALL Halt('There are duplicate Key then looking for the Key=<' //TRIM(Key)//'>')
       ENDIF
    ENDDO
    !
  END FUNCTION OptGetKeyArg
  !
  !
END MODULE

