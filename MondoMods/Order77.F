      SUBROUTINE SFCOrder77(N,R,Point,Key,Hilbert)
         IMPLICIT NONE
         INTEGER I,J,K,L,M,N,Ix,Iy,Iz,K1,K2
         REAL*8 R(3,N),RMin(3),ShiftDbles,Seg
         INTEGER Point(N)
         INTEGER*8 Key(N)
         INTEGER*8 Interleave,HilbertKey
         EXTERNAL  Interleave,HilbertKey
         LOGICAL   Hilbert
         PARAMETER(ShiftDbles=1.D3) 
         RMin(1)=1.D10
         RMin(2)=1.D10
         RMin(3)=1.D10
         DO J=1,N
            Point(J)=J
            RMin(1)=DMIN1(RMin(1),R(1,J))
            RMin(2)=DMIN1(RMin(2),R(2,J))
            RMin(3)=DMIN1(RMin(3),R(3,J))
         ENDDO
         IF(Hilbert)THEN
            DO J=1,N
               Ix=DNINT((R(1,J)-RMin(1))*ShiftDbles)
               Iy=DNINT((R(2,J)-RMin(2))*ShiftDbles)
               Iz=DNINT((R(3,J)-RMin(3))*ShiftDbles)
               Key(J)=Interleave(Ix,Iy,Iz)
               Key(J)=HilbertKey(Key(J))
            END DO
         ELSE
            DO J=1,N
               Ix=DNINT((R(1,J)-RMin(1))*ShiftDbles)
               Iy=DNINT((R(2,J)-RMin(2))*ShiftDbles)
               Iz=DNINT((R(3,J)-RMin(3))*ShiftDbles)
               Key(J)=Interleave(Ix,Iy,Iz)
            END DO
         ENDIF
         CALL I8Sort(Key,Point,N,2)
         RETURN 
      END !SUBROUTINE SFCOrder77
!-----------------------------------------------------------------
!     Bit shuffle for a triple of integers
! 
      FUNCTION Interleave(Ix,Iy,Iz)
         IMPLICIT NONE
         INTEGER I,K,Ix,Iy,Iz
         INTEGER*8 Interleave
         logical print
         K=0
         Interleave=0
         DO I=0,60,3
            IF(BTEST(Ix,K))Interleave=IBSET(Interleave,I+2)
            IF(BTEST(Iy,K))Interleave=IBSET(Interleave,I+1)
            IF(BTEST(Iz,K))Interleave=IBSET(Interleave,I  )
            K=K+1
         ENDDO
         RETURN
      END ! FUNCTION Interleave
!len=20;  dim=3;
!ix=5391; iy=3911; iz=0;
!l[1]=IntegerDigits[ix,2,len];
!l[2]=IntegerDigits[iy,2,len];
!l[3]=IntegerDigits[iz,2,len];
!n=Table[0,{i,dim*len}];
!Print[" IX = ",BaseForm[ix,2]];
!Print[" IY = ",BaseForm[iy,2]];
!Print[" IZ = ",BaseForm[iz,2]];
!k=len;
!Do[
!   Do[n[[i-j+1]]=l[dim-j+1][[k]];		
!      l[dim-j+1][[k]]];*)
!     ,{j,1,3}];
!   k=k-1;
!  ,{i,dim*len,1,-dim}];
!key=FromDigits[n,2];
!Print[" Key = ",key];
!Print[" Key = ",BaseForm[key,2]];
!-----------------------------------------------------------------
!     Hilbert ordering for three dimensional, 64 bit, interleaved key
!     Based on Algortithm H2 of Faloutsos and Rosen, PODS 89.
!
      FUNCTION HilbertKey(Key)
         IMPLICIT NONE
         INTEGER I,J,K,L
         INTEGER*8 Key
         INTEGER*8 HilbertKey
         INTEGER SubKey(21)
         INTEGER ToState(0:7,12)
         INTEGER ToBinry(0:7,12)
!--------------------------------------------------------------        
!        State table, Based on Fig 3, T. Bially, 
!        IEEE Trans. Inf. Theory, 1969
!                                     0, 1, 2, 3, 4, 5, 6, 7    
         DATA (ToState(I,1 ),I=0,7) / 9,12, 2, 1, 2, 3, 9, 1/
         DATA (ToBinry(I,1 ),I=0,7) / 0, 1, 3, 2, 7, 6, 4, 5/
         DATA (ToState(I,2 ),I=0,7) / 5, 3, 2, 2, 3, 5, 1, 4/
         DATA (ToBinry(I,2 ),I=0,7) / 4, 3, 5, 2, 7, 0, 6, 1/
         DATA (ToState(I,3 ),I=0,7) / 2, 3, 6, 3, 1, 7, 7, 1/     
         DATA (ToBinry(I,3 ),I=0,7) / 6, 5, 1, 2, 7, 4, 0, 3/
         DATA (ToState(I,4 ),I=0,7) /10, 9, 4, 2, 5, 2, 4, 9/
         DATA (ToBinry(I,4 ),I=0,7) / 6, 7, 5, 4, 1, 0, 2, 3/
         DATA (ToState(I,5 ),I=0,7) / 5, 2, 5, 6, 8, 4, 4, 8/      
         DATA (ToBinry(I,5 ),I=0,7) / 2, 1, 5, 6, 3, 0, 4, 7/
         DATA (ToState(I,6 ),I=0,7) / 6, 6, 5, 3, 7, 8, 3, 5/
         DATA (ToBinry(I,6 ),I=0,7) / 2, 5, 3, 4, 1, 6, 0, 7/
         DATA (ToState(I,7 ),I=0,7) / 6, 7,11,12,11, 7, 6, 3/
         DATA (ToBinry(I,7 ),I=0,7) / 4, 5, 7, 6, 3, 2, 0, 1/
         DATA (ToState(I,8 ),I=0,7) / 8, 6,10,11, 8,11, 5, 6/
         DATA (ToBinry(I,8 ),I=0,7) / 2, 3, 1, 0, 5, 4, 6, 7/
         DATA (ToState(I,9 ),I=0,7) /12,10, 1, 4,10,12, 9, 9/ 
         DATA (ToBinry(I,9 ),I=0,7) / 0, 7, 1, 6, 3, 4, 2, 5/
         DATA (ToState(I,10),I=0,7) / 8, 4, 4, 8,10, 9,10,11/ 
         DATA (ToBinry(I,10),I=0,7) / 4, 7, 3, 0, 5, 6, 2, 1/
         DATA (ToState(I,11),I=0,7) / 7, 8,12,10,11,11,10,12/
         DATA (ToBinry(I,11),I=0,7) / 6, 1, 7, 0, 5, 2, 4, 3/
         DATA (ToState(I,12),I=0,7) / 1, 7, 7, 1, 9,12,11,12/
         DATA (ToBinry(I,12),I=0,7) / 0, 3, 7, 4, 1, 2, 6, 5/
!
!        Form 3 bit SubKeys from left to right (Step 3 of H2)
!
         DO I=1,21
            SubKey(I)=0
         ENDDO
         K=1
         DO I=60,0,-3
            IF(BTEST(Key,I  ))SubKey(K)=IBSET(SubKey(K),0)
            IF(BTEST(Key,I+1))SubKey(K)=IBSET(SubKey(K),1)
            IF(BTEST(Key,I+2))SubKey(K)=IBSET(SubKey(K),2)
            IF(K.GT.21)STOP ' gt 21 '
            IF(I.LT.0)STOP ' lt 0 '
            K=K+1
         ENDDO
!
!        Change each 3 bit SubKey according to the output key
!        from the current state and move to a new state (Step 4 of H2)
!
         I=1
         DO K=1,21
            L=SubKey(K)
C            WRITE(*,77)I,K,SubKey(K),ToState(L,I),K,ToBinry(L,I)
            SubKey(K)=ToBinry(L,I)
            I        =ToState(L,I)
         ENDDO
!
!        Reassemble key from sub keys, with SubKey(21) 
!        the rightmost 3 bits (Step 5 of H2)
!
         HilbertKey=0
         I=0
         DO K=21,1,-1
            IF(BTEST(SubKey(K),0))HilbertKey=IBSET(HilbertKey,I  )
            IF(BTEST(SubKey(K),1))HilbertKey=IBSET(HilbertKey,I+1)
            IF(BTEST(SubKey(K),2))HilbertKey=IBSET(HilbertKey,I+2)
            I=I+3
         ENDDO
         RETURN
 77      FORMAT(' State = ',I2,' S(',I2,') = ',I2,
     &          ' NewState = ',I2,' S2(',I2,') = ',I2)
      END ! FUNCTION HilbertKey

!      PROGRAM Driver
!         INTEGER ix,iy,iz,I,J,K,L,M,N
!         PARAMETER(N=8**3)
!         INTEGER Point(N)
!         INTEGER IR(3,N),IR2(3,N)
!         INTEGER*8 Interleave, HilbertKey, Key(N)
!         EXTERNAL  Interleave, HilbertKey          
!         K=0
!         DO Ix=0,7
!            DO Iy=0,7
!               DO Iz=0,7
!                  K=K+1      
!                  IR(1,K)=IX
!                  IR(2,K)=IY
!                  IR(3,K)=IZ
!                  Key(K)=Interleave(Ix,Iy,Iz)
!                  Key(K)=HilbertKey(Key(K))
!               ENDDO
!            ENDDO
!         ENDDO
!         DO I=1,N
!            Point(I)=I
!         ENDDO
!         CALL I8Sort(Key,Point,N,2)
!         DO I=1,N
!            IR2(1,I)=IR(1,Point(I))
!            IR2(2,I)=IR(2,Point(I))
!            IR2(3,I)=IR(3,Point(I))
!         ENDDO
!         DO I=1,N
!            WRITE(77,*)(IR2(K,I),K=1,3)
!         ENDDO
!      STOP
!      END
!



!-----------------------------------------------------------------
!     Quick sort of doubles, carrying allong doubles
! 
      SUBROUTINE DblIntSort77(N,DX,DY,KFLAG)                               
      INTEGER KFLAG, N                                                  
      REAL*8 DX(N)
      REAL*8 R, T, TT
      INTEGER DY(N)                                     
      INTEGER I, IJ, J, K, KK, L, M, NN, TY, TTY 
      INTEGER IL(21), IU(21)                                            
      INTRINSIC ABS, INT                                                
C
C      KFLAG=-2
C
      NN = N                                                            
      IF (NN .LT. 1) THEN                                               
         write(*,*)                                                     
     +      'the number of values to be sorted is not positive.',nn        
         STOP 
         RETURN                                                         
      ENDIF                                                             
C                                                                       
      KK = ABS(KFLAG)                                                   
      IF (KK.NE.1 .AND. KK.NE.2) THEN                                   
         write(*,*)                                                     
     +      'The sort control parameter, K, is not 2, 1, -1, or -2.'    
         RETURN                                                         
      ENDIF                                                             
C                                                                       
C     Alter array DX to get decreasing order if needed                  
C                                                                       
      IF (KFLAG .LE. -1) THEN                                           
         DO 10 I=1,NN                                                   
            DX(I) = -DX(I)                                              
   10    CONTINUE                                                       
      ENDIF                                                             
C                                                                       
  100 M = 1                                                             
      I = 1                                                             
      J = NN                                                            
      R = 0.375D0                                                       
C                                                                       
  110 IF (I .EQ. J) GO TO 150                                           
      IF (R .LE. 0.5898437D0) THEN                                      
         R = R+3.90625D-2                                               
      ELSE                                                              
         R = R-0.21875D0                                                
      ENDIF                                                             
C                                                                       
  120 K = I                                                             
C                                                                       
C     Select a central element of the array and save it in location T   
C                                                                       
      IJ = I + INT((J-I)*R)                                             
      T = DX(IJ)                                                        
      TY = DY(IJ)                                                       
C                                                                       
C     If first element of array is greater than T, interchange with T   
C                                                                       
      IF (DX(I) .GT. T) THEN                                            
         DX(IJ) = DX(I)                                                 
         DX(I) = T                                                      
         T = DX(IJ)                                                     
         DY(IJ) = DY(I)                                                 
         DY(I) = TY                                                     
         TY = DY(IJ)                                                    
      ENDIF                                                             
      L = J                                                             
C                                                                       
C     If last element of array is less than T, interchange with T       
C                                                                       
      IF (DX(J) .LT. T) THEN                                            
         DX(IJ) = DX(J)                                                 
         DX(J) = T                                                      
         T = DX(IJ)                                                     
         DY(IJ) = DY(J)                                                 
         DY(J) = TY                                                     
         TY = DY(IJ)                                                    
C                                                                       
C        If first element of array is greater than T, interchange with T
C                                                                       
         IF (DX(I) .GT. T) THEN                                         
            DX(IJ) = DX(I)                                              
            DX(I) = T                                                   
            T = DX(IJ)                                                  
            DY(IJ) = DY(I)                                              
            DY(I) = TY                                                  
            TY = DY(IJ)                                                 
         ENDIF                                                          
      ENDIF                                                             
C                                                                       
C     Find an element in the second half of the array which is smaller  
C     than T                                                            
C                                                                       
  130 L = L-1                                                           
      IF (DX(L) .GT. T) GO TO 130                                       
C                                                                       
C     Find an element in the first half of the array which is greater   
C     than T                                                            
C                                                                       
  140 K = K+1                                                           
      IF (DX(K) .LT. T ) GO TO 140
C                                                                       
C     Interchange these elements                                        
C                                                                       
      IF (K .LE. L) THEN                                                
         TT = DX(L)                                                     
         DX(L) = DX(K)                                                  
         DX(K) = TT                                                     
         TTY = DY(L)                                                    
         DY(L) = DY(K)                                                  
         DY(K) = TTY                                                    
         GO TO 130                                                      
      ENDIF                                                             
C                                                                       
C     Save upper and lower subscripts of the array yet to be sorted     
C                                                                       
      IF (L-I .GT. J-K) THEN                                            
         IL(M) = I                                                      
         IU(M) = L                                                      
         I = K                                                          
         M = M+1                                                        
      ELSE                                                              
         IL(M) = K                                                      
         IU(M) = J                                                      
         J = L                                                          
         M = M+1                                                        
      ENDIF                                                             
      GO TO 160                                                         
C                                                                       
C     Begin again on another portion of the unsorted array              
C                                                                       
  150 M = M-1                                                           
      IF (M .EQ. 0) GO TO 190                                           
      I = IL(M)                                                         
      J = IU(M)                                                         
C                                                                       
  160 IF (J-I .GE. 1) GO TO 120                                         
      IF (I .EQ. 1) GO TO 110                                           
      I = I-1                                                           
C                                                                       
  170 I = I+1                                                           
      IF (I .EQ. J ) GO TO 150
      T = DX(I+1)                                                       
      TY = DY(I+1)                                                      
      IF (DX(I) .LE. T) GO TO 170                                       
      K = I                                                             
C                                                                       
  180 DX(K+1) = DX(K)                                                   
      DY(K+1) = DY(K)                                                   
      K = K-1                                                           
      IF (T .LT. DX(K)) GO TO 180                                       
      DX(K+1) = T                                                       
      DY(K+1) = TY                                                      
      GO TO 170                                                         
C                                                                       
C     Clean up                                                          
C                                                                       
  190 IF (KFLAG .LE. -1) THEN                                           
         DO 200 I=1,NN                                                  
            DX(I) = -DX(I)                                              
  200    CONTINUE                                                       
      ENDIF                                                             
      RETURN                                                            
      END                                                               
                                                                        
!-----------------------------------------------------------------
!     Quick sort modified to handle INTEGER*8 IX
! 
      SUBROUTINE I8SORT (IX,IY,N,KFLAG)
C***BEGIN PROLOGUE  ISORT
C***PURPOSE  Sort an array and optionally make the same interchanges in
C            an auxiliary array.  The array may be sorted in increasing
C            or decreasing order.  A slightly modified QUICKSORT
C            algorithm is used.
C***LIBRARY   SLATEC
C***CATEGORY  N6A2A
C***TYPE      INTEGER (SSORT-S, DSORT-D, ISORT-I)
C***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Wisniewski, J. A., (SNLA)
C***DESCRIPTION
C
C   ISORT sorts array IX and optionally makes the same interchanges in
C   array IY.  The array IX may be sorted in increasing order or
C   decreasing order.  A slightly modified quicksort algorithm is used.
C
C   Description of Parameters
C      IX - integer array of values to be sorted
C      IY - integer array to be (optionally) carried along
C      N  - number of values in integer array IX to be sorted
C      KFLAG - control parameter
C            =  2  means sort IX in increasing order and carry IY along.
C            =  1  means sort IX in increasing order (ignoring IY)
C            = -1  means sort IX in decreasing order (ignoring IY)
C            = -2  means sort IX in decreasing order and carry IY along.
C
C***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
C                 for sorting with minimal storage, Communications of
C                 the ACM, 12, 3 (1969), pp. 185-187.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   761118  DATE WRITTEN
C   810801  Modified by David K. Kahaner.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891009  Removed unreferenced statement labels.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   901012  Declared all variables; changed X,Y to IX,IY. (M. McClain)
C   920501  Reformatted the REFERENCES section.  (DWL, WRB)
C   920519  Clarified error messages.  (DWL)
C   920801  Declarations section rebuilt and code restructured to use
C           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
C***END PROLOGUE  ISORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      INTEGER*8 IX(N)
      INTEGER   IY(N)
C     .. Local Scalars ..
      REAL R
      INTEGER*8 T,TT
      INTEGER I, IJ, J, K, KK, L, M, NN, TTY, TY       
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
C     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
C***FIRST EXECUTABLE STATEMENT  ISORT
      NN = N
      IF (NN .LT. 1) THEN
C        The number of values to be sorted is not positive
         STOP ' Died from error 1 in I8Sort '
      ENDIF
C
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
C        The sort control parameter, K, is not 2, 1, -1, or -2
        STOP ' Died from error 2 in I8Sort '
      ENDIF
C
C     Alter array IX to get decreasing order if needed
C
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            IX(I) = -IX(I)
   10    CONTINUE
      ENDIF
C
      IF (KK .EQ. 2) GO TO 100
C
C     Sort IX only
C
      M = 1
      I = 1
      J = NN
      R = 0.375E0
C
   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
   30 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
   40 L = L-1
      IF (IX(L) .GT. T) GO TO 40
C
C     Find an element in the first half of the array which is greater
C     than T
C
   50 K = K+1
      IF (IX(K) .LT. T) GO TO 50
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         GO TO 40
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70
C
C     Begin again on another portion of the unsorted array
C
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1
C
   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = IX(I+1)
      IF (IX(I) .LE. T) GO TO 80
      K = I
C
   90 IX(K+1) = IX(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 90
      IX(K+1) = T
      GO TO 80
C
C     Sort IX and carry IY along
C
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0
C
  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
C
  120 K = I
C
C     Select a central element of the array and save it in location T
C
      IJ = I + INT((J-I)*R)
      T = IX(IJ)
      TY = IY(IJ)
C
C     If first element of array is greater than T, interchange with T
C
      IF (IX(I) .GT. T) THEN
         IX(IJ) = IX(I)
         IX(I) = T
         T = IX(IJ)
         IY(IJ) = IY(I)
         IY(I) = TY
         TY = IY(IJ)
      ENDIF
      L = J
C
C     If last element of array is less than T, interchange with T
C
      IF (IX(J) .LT. T) THEN
         IX(IJ) = IX(J)
         IX(J) = T
         T = IX(IJ)
         IY(IJ) = IY(J)
         IY(J) = TY
         TY = IY(IJ)
C
C        If first element of array is greater than T, interchange with T
C
         IF (IX(I) .GT. T) THEN
            IX(IJ) = IX(I)
            IX(I) = T
            T = IX(IJ)
            IY(IJ) = IY(I)
            IY(I) = TY
            TY = IY(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is smaller
C     than T
C
  130 L = L-1
      IF (IX(L) .GT. T) GO TO 130
C
C     Find an element in the first half of the array which is greater
C     than T
C
  140 K = K+1
      IF (IX(K) .LT. T) GO TO 140
C
C     Interchange these elements
C
      IF (K .LE. L) THEN
         TT = IX(L)
         IX(L) = IX(K)
         IX(K) = TT
         TTY = IY(L)
         IY(L) = IY(K)
         IY(K) = TTY
         GO TO 130
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted
C
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160
C
C     Begin again on another portion of the unsorted array
C
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)
C
  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1
C
  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = IX(I+1)
      TY = IY(I+1)
      IF (IX(I) .LE. T) GO TO 170
      K = I
C
  180 IX(K+1) = IX(K)
      IY(K+1) = IY(K)
      K = K-1
      IF (T .LT. IX(K)) GO TO 180
      IX(K+1) = T
      IY(K+1) = TY
      GO TO 170
C
C     Clean up
C
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            IX(I) = -IX(I)
  200    CONTINUE
      ENDIF
      RETURN
      END

      SUBROUTINE IntIntSort77(N,IA,JA, KFLAG)

C       modified to work with integer arrays IA and JA. WMC 97
c
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C         array A.
C         The array IA may be sorted in increasing order or decreasing
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
c
      implicit integer (i-m)
c
      real*8 R
c
      INTEGER N,NN
      DIMENSION IA(N), JA(N)
c
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
c
      INTRINSIC ABS, INT

      NN = N
      IF (NN.LT.1) THEN
        write(*,*)' The number of values to be sorted was not positive.'
        stop
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
        write(*,*)'The sort control parameter, K, was not 1 or -1.'
        stop
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location
C     it, jt, at.
C
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
      JT = JA(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
      ENDIF
      L=J
C
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      END



      SUBROUTINE IntSort77(N,IA,KFLAG)

C       modified to work with only the integer array IA. WMC 99
c
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C         array A.
C         The array IA may be sorted in increasing order or decreasing
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
c
      implicit integer (i-m)
c
      real*8 R
c
      INTEGER N,NN
      DIMENSION IA(N)
c
C     .. Local Arrays ..
      INTEGER IL(21), IU(21)
c
      INTRINSIC ABS, INT

      NN = N
      IF (NN.LT.1) THEN
        write(*,*)' The number of values to be sorted was not positive.'
        stop
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
        write(*,*)'The sort control parameter, K, was not 1 or -1.'
        stop
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location
C     it, jt, at.
C
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
      ENDIF
      L=J
C
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
      END
